package main

import (
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"time"
)

type config struct {
	Machine struct {
		BaseOutput string `json:"base_output"`
	} `json:"machine"`
	Vars struct {
		Input struct {
			File string `json:"file"`
			Type string `json:"type"`
		} `json:"input"`
		Output struct {
			Dir string `json:"dir"`
		} `json:"output"`
	} `json:"vars"`
	Selection struct {
		Charges      []int `json:"charges"`
		Centralities []int `json:"centralities"`
	} `json:"selection"`
	Binning struct {
		Names     []string `json:"names"`
		FileNames []string `json:"file_names"`
	} `json:"binning"`
	Fit struct {
		Freeze map[string]*float64 `json:"freeze"`
	} `json:"fit"`
}

// Individual records are additive schema-1 diagnostics. Older manifests have
// only aggregate counters and cannot support an honest good/bad classification.
type FitDetail struct {
	Charge         int      `json:"ch"`
	Centrality     int      `json:"centr"`
	Bin            int      `json:"bin"`
	Usable         bool     `json:"usable"`
	Finite         bool     `json:"finite"`
	Chi2           *float64 `json:"chi2"`
	Ndf            int      `json:"ndf"`
	Status         int      `json:"status"`
	CovStatus      int      `json:"cov_status"`
	AtLimit        bool     `json:"at_limit"`
	Attempts       int      `json:"attempts"`
	CFImage        string   `json:"cf_image"`
	FitOverCFImage string   `json:"fit_over_cf_image"`
	Reason         string   `json:"reason"`
}

type manifest struct {
	SchemaVersion int    `json:"schema_version"`
	Status        string `json:"status"`
	RunID         string `json:"run_id"`
	CompletedAt   string `json:"completed_at"`
	ExitCode      *int   `json:"exit_code"`
	Config        config `json:"config"`
	Fits          struct {
		Requested int         `json:"requested"`
		Usable    int         `json:"usable"`
		Unusable  int         `json:"unusable"`
		Retried   int         `json:"retried"`
		AtLimit   int         `json:"at_limit"`
		Details   []FitDetail `json:"details"`
	} `json:"fits"`
	Images []struct {
		Path string `json:"path"`
		Size int64  `json:"size_bytes"`
	} `json:"images"`
}

type plot struct {
	source   string
	relative string
	asset    string
	group    string
	title    string
	caption  string
}

type CompactSlide struct {
	Title   string
	Kind    string
	Plots   []plot
	Note    string
	Missing string
}

func main() {
	if err := run(os.Args[1:]); err != nil {
		fmt.Fprintln(os.Stderr, "presentation-builder:", err)
		os.Exit(1)
	}
}

func run(args []string) (err error) {
	flags := flag.NewFlagSet("presentation-builder", flag.ContinueOnError)
	flags.SetOutput(os.Stderr)
	configName := flags.String("config", "config/kt.json", "analysis config file")
	resultsName := flags.String("results", "", "use this results directory instead of the config output")
	repoName := flags.String("repo-root", ".", "repository root")
	sourceName := flags.String("source-dir", "docs/thesis", "directory containing presentation/main.tex")
	buildName := flags.String("build-dir", "docs/thesis/build/presentation", "presentation build directory")
	latexmk := flags.String("latexmk", "latexmk", "latexmk executable")
	maxPlots := flags.Int("plots", 12, "maximum plots in the results chapter; 0 means all")
	mode := flags.String("mode", "compact", "report layout: compact (current results only) or full (thesis slides)")
	strictAssets := flags.Bool("strict-assets", false, "fail when source presentation figures are missing")
	executable := flags.String("executable", "", "analysis executable; set to run it before building")
	if err := flags.Parse(args); err != nil {
		if errors.Is(err, flag.ErrHelp) {
			return nil
		}
		return err
	}
	if flags.NArg() != 0 {
		return errors.New("unexpected positional arguments")
	}
	if *maxPlots < 0 {
		return errors.New("-plots must be non-negative (0 includes all plots)")
	}
	if *mode != "compact" && *mode != "full" {
		return errors.New("-mode must be compact or full")
	}

	repoRoot, err := filepath.Abs(*repoName)
	if err != nil {
		return err
	}
	configPath := resolve(repoRoot, *configName)
	var cfg config
	if err := decodeJSON(configPath, &cfg); err != nil {
		return fmt.Errorf("read config %s: %w", configPath, err)
	}
	if cfg.Machine.BaseOutput == "" {
		return errors.New("config is missing machine.base_output")
	}
	if cfg.Vars.Output.Dir == "" {
		cfg.Vars.Output.Dir = "results"
	}
	results := resolve(repoRoot, cfg.Machine.BaseOutput+"/"+cfg.Vars.Output.Dir)
	if *resultsName != "" {
		results = resolve(repoRoot, *resultsName)
	}
	sourceDir := resolve(repoRoot, *sourceName)
	buildDir := resolve(repoRoot, *buildName)

	if _, err := exec.LookPath(*latexmk); err != nil {
		return fmt.Errorf("latexmk not found: %w", err)
	}
	if err := os.MkdirAll(buildDir, 0o755); err != nil {
		return err
	}
	chapterPath := filepath.Join(buildDir, "presentation-results.tex")
	clearChapter := func() error { return os.WriteFile(chapterPath, nil, 0o644) }
	defer func() {
		if err != nil {
			_ = clearChapter()
		}
	}()

	var previousRunID string
	var analysisExit *int
	if *executable != "" {
		if old, _, oldErr := readManifest(results); oldErr == nil {
			previousRunID = old.RunID
		}
		code, runErr := runAnalysis(*executable, configPath, repoRoot)
		if runErr != nil {
			return runErr
		}
		if code != 0 && code != 2 {
			return fmt.Errorf("analysis exited with code %d; presentation was not built", code)
		}
		analysisExit = &code
	}

	data, digest, manifestErr := readManifest(results)
	if manifestErr != nil {
		if *mode == "compact" || *resultsName != "" || analysisExit != nil || !errors.Is(manifestErr, os.ErrNotExist) {
			return manifestErr
		}
		if err := clearChapter(); err != nil {
			return err
		}
		fmt.Println("Расчёт не найден; собираются только исходные слайды.")
	} else {
		if analysisExit != nil {
			if data.RunID == previousRunID {
				return errors.New("analysis did not create a manifest for this run")
			}
			if data.ExitCode == nil || *analysisExit != *data.ExitCode {
				return errors.New("analysis exit code and run manifest disagree")
			}
		}
		if err := writeChapter(buildDir, results, data, digest, *maxPlots, chapterPath, *mode); err != nil {
			return err
		}
	}

	if err := compilePDF(*latexmk, sourceDir, buildDir, *strictAssets, *mode); err != nil {
		return err
	}
	pdfPath := filepath.Join(buildDir, "presentation.pdf")
	info, err := os.Stat(pdfPath)
	if err != nil || info.Size() == 0 {
		return fmt.Errorf("latexmk did not create %s", pdfPath)
	}
	if manifestErr == nil {
		if _, finalDigest, err := readManifest(results); err != nil || finalDigest != digest {
			return errors.New("analysis results changed during PDF creation; PDF was not published")
		}
		if err := copyAtomic(pdfPath, filepath.Join(results, "presentation.pdf")); err != nil {
			return fmt.Errorf("publish PDF: %w", err)
		}
		if data.Status == "incomplete" {
			fmt.Println("Собран диагностический PDF: часть фитов непригодна.")
		}
		fmt.Println("PDF с результатами:", filepath.Join(results, "presentation.pdf"))
	}
	fmt.Println("PDF сборки:", pdfPath)
	return nil
}

func resolve(root, value string) string {
	if filepath.IsAbs(value) {
		return filepath.Clean(value)
	}
	return filepath.Clean(filepath.Join(root, value))
}

func decodeJSON(path string, destination any) error {
	file, err := os.Open(path)
	if err != nil {
		return err
	}
	defer file.Close()
	decoder := json.NewDecoder(file)
	if err := decoder.Decode(destination); err != nil {
		return err
	}
	if err := decoder.Decode(&struct{}{}); err != io.EOF {
		return errors.New("unexpected data after JSON object")
	}
	return nil
}

func runAnalysis(executable, configPath, repoRoot string) (int, error) {
	command := exec.Command(executable, configPath)
	command.Dir = repoRoot
	command.Stdout = os.Stdout
	command.Stderr = os.Stderr
	if err := command.Run(); err != nil {
		var exitErr *exec.ExitError
		if errors.As(err, &exitErr) {
			return exitErr.ExitCode(), nil
		}
		return 1, err
	}
	return 0, nil
}

func readManifest(results string) (manifest, string, error) {
	path := filepath.Join(results, "run_manifest.json")
	raw, err := os.ReadFile(path)
	if err != nil {
		return manifest{}, "", fmt.Errorf("read %s: %w", path, err)
	}
	var data manifest
	if err := json.Unmarshal(raw, &data); err != nil {
		return manifest{}, "", fmt.Errorf("parse %s: %w", path, err)
	}
	if data.SchemaVersion != 1 || data.RunID == "" || data.Status != "completed" && data.Status != "incomplete" {
		return manifest{}, "", errors.New("run manifest is missing or the analysis has not completed")
	}
	expectedCode := 0
	if data.Status == "incomplete" {
		expectedCode = 2
	}
	if data.ExitCode == nil || *data.ExitCode != expectedCode {
		return manifest{}, "", errors.New("run manifest status and exit code disagree")
	}
	if data.Config.Machine.BaseOutput == "" || data.Config.Vars.Output.Dir == "" {
		return manifest{}, "", errors.New("run manifest has no output directory in its config snapshot")
	}
	manifestResults := data.Config.Vars.Output.Dir
	if !filepath.IsAbs(manifestResults) {
		return manifest{}, "", errors.New("run manifest output path is not absolute")
	}
	if filepath.Clean(manifestResults) != filepath.Clean(results) {
		return manifest{}, "", errors.New("run manifest belongs to another results directory")
	}
	for _, count := range []int{data.Fits.Requested, data.Fits.Usable, data.Fits.Unusable, data.Fits.Retried, data.Fits.AtLimit} {
		if count < 0 {
			return manifest{}, "", errors.New("run manifest has a negative fit counter")
		}
	}
	if data.Fits.Requested != data.Fits.Usable+data.Fits.Unusable {
		return manifest{}, "", errors.New("run manifest fit counters disagree")
	}
	if data.Status == "completed" && data.Fits.Unusable != 0 {
		return manifest{}, "", errors.New("completed manifest contains unusable fits")
	}
	digest := sha256.Sum256(raw)
	return data, hex.EncodeToString(digest[:]), nil
}

func writeChapter(buildDir, results string, data manifest, digest string, maxPlots int, chapterPath, mode string) error {
	plots, err := collectPlots(results, data)
	if err != nil {
		return err
	}
	assetDir := filepath.Join(buildDir, "current-results", safeRunID(data.RunID))
	if err := os.MkdirAll(assetDir, 0o755); err != nil {
		return err
	}
	copied := make(map[string]string)
	prepare := func(items []plot) error {
		for i := range items {
			asset, exists := copied[items[i].source]
			if !exists {
				name := fmt.Sprintf("plot-%04d%s", len(copied)+1, strings.ToLower(filepath.Ext(items[i].source)))
				destination := filepath.Join(assetDir, name)
				if err := copyFile(items[i].source, destination); err != nil {
					return err
				}
				relative, err := filepath.Rel(buildDir, destination)
				if err != nil {
					return err
				}
				asset = filepath.ToSlash(relative)
				copied[items[i].source] = asset
			}
			items[i].asset = asset
		}
		return nil
	}
	var text string
	if mode == "compact" {
		slides := select_compact_slides(data, plots)
		for i := range slides {
			if err := prepare(slides[i].Plots); err != nil {
				return err
			}
		}
		text = render_compact(data, slides)
	} else {
		selected := selectPlots(plots, maxPlots)
		if err := prepare(selected); err != nil {
			return err
		}
		text = render(data, selected, len(plots))
	}
	if _, latestDigest, err := readManifest(results); err != nil || latestDigest != digest {
		return errors.New("analysis results changed while preparing slides")
	}
	return writeAtomic(chapterPath, []byte(text))
}

func collectPlots(results string, data manifest) ([]plot, error) {
	resolvedResults, err := filepath.EvalSymlinks(results)
	if err != nil {
		return nil, err
	}
	labels := makeLabels(data.Config)
	plots := make([]plot, 0, len(data.Images))
	seen := make(map[string]bool, len(data.Images))
	for _, image := range data.Images {
		rel := filepath.Clean(filepath.FromSlash(image.Path))
		if filepath.IsAbs(rel) || rel == "." || rel == ".." || strings.HasPrefix(rel, ".."+string(filepath.Separator)) || seen[rel] {
			return nil, fmt.Errorf("invalid or duplicate image path in run manifest: %q", image.Path)
		}
		seen[rel] = true
		source := filepath.Join(results, rel)
		resolved, err := filepath.EvalSymlinks(source)
		if err != nil {
			return nil, fmt.Errorf("image from this run is missing: %s", source)
		}
		if !inside(resolvedResults, resolved) {
			return nil, fmt.Errorf("image path leaves results directory: %s", image.Path)
		}
		info, err := os.Stat(resolved)
		if err != nil || !info.Mode().IsRegular() || info.Size() == 0 || info.Size() != image.Size {
			return nil, fmt.Errorf("image from this run is empty or changed: %s", source)
		}
		ext := strings.ToLower(filepath.Ext(resolved))
		if ext != ".pdf" && ext != ".png" && ext != ".jpg" && ext != ".jpeg" {
			return nil, fmt.Errorf("unsupported image format %q; use pdf, png or jpg", ext)
		}
		item := labels[filepath.ToSlash(strings.TrimSuffix(rel, filepath.Ext(rel)))]
		if item.title == "" {
			item = plot{group: "other", title: "Результат расчёта", caption: texEscape(filepath.Base(rel))}
		}
		item.source = resolved
		item.relative = filepath.ToSlash(rel)
		plots = append(plots, item)
	}
	return plots, nil
}

func makeLabels(cfg config) map[string]plot {
	labels := make(map[string]plot)
	mode := "$y$"
	if cfg.Vars.Input.Type == "kt" {
		mode = "$k_t$"
	}
	charges := map[int]struct{ slug, label string }{0: {"pos", `$\pi^+$`}, 1: {"neg", `$\pi^-$`}}
	centralities := map[int]string{0: "0-10", 1: "10-30", 2: "30-50", 3: "50-80"}
	for _, chargeID := range cfg.Selection.Charges {
		charge, ok := charges[chargeID]
		if !ok {
			continue
		}
		for _, item := range []struct{ file, group, title, caption string }{
			{"c_all_graphs_" + charge.slug, "radii", "Радиусы и сила корреляции", "зависимость от " + mode},
			{"c_fit_quality_" + charge.slug, "quality", "Качество аппроксимации", "зависимость от " + mode},
			{"c_pvalues_" + charge.slug, "pvalues", "Вероятности согласия", "зависимость от " + mode},
			{"c_all_cross_graphs_" + charge.slug, "cross", "Перекрёстные параметры", "зависимость от " + mode},
			{"c_free_cross_graphs_" + charge.slug, "free_cross", "Свободные перекрёстные параметры", "зависимость от " + mode},
			{"c_chi2_root_over_manual_" + charge.slug, "chi2_comparison", "Сравнение ROOT и ручного расчёта", `$\chi^2_{\mathrm{ROOT}}/\chi^2_{\mathrm{manual}}$; зависимость от ` + mode},
		} {
			labels["dependency/"+item.file] = plot{group: item.group, title: item.title, caption: charge.label + ", " + item.caption}
		}
		for _, centralityID := range cfg.Selection.Centralities {
			centrality, ok := centralities[centralityID]
			if !ok {
				continue
			}
			sample := charge.label + ", центральность " + centrality + `\%`
			key := "all_2d_histos/all_out-long_2d_histos_centr_" + centrality + "_" + charge.slug
			labels[key] = plot{group: "2d", title: "Двумерные проекции КФ", caption: sample + ", плоскость out–long"}
			for i, name := range cfg.Binning.Names {
				if i >= len(cfg.Binning.FileNames) {
					break
				}
				file := cfg.Binning.FileNames[i]
				caption := sample + ", " + mode + ": " + texEscape(short(name, 35))
				labels["all_1d_histos/cfs_"+charge.slug+"_"+centrality+"_"+file] = plot{
					group: "1d", title: "Одномерные проекции КФ", caption: caption}
				labels["dependency/fit_over_cf/fit_over_cf_"+charge.slug+"_"+centrality+"_"+file] = plot{
					group: "fit_ratio", title: "Отношение аппроксимации к КФ", caption: caption}
			}
		}
	}
	return labels
}

func selectPlots(plots []plot, maximum int) []plot {
	if maximum == 0 || len(plots) <= maximum {
		return plots
	}
	order := []string{"radii", "quality", "chi2_comparison", "1d", "2d", "pvalues", "fit_ratio", "cross", "other"}
	groups := make(map[string][]plot, len(order))
	for _, item := range plots {
		groups[item.group] = append(groups[item.group], item)
	}
	selected := make([]plot, 0, maximum)
	for len(selected) < maximum {
		added := false
		for _, group := range order {
			items := groups[group]
			if len(items) == 0 || len(selected) == maximum {
				continue
			}
			selected = append(selected, items[0])
			groups[group] = items[1:]
			added = true
		}
		if !added {
			break
		}
	}
	return selected
}

func reduced_chi2(detail FitDetail) (float64, bool) {
	if detail.Chi2 == nil || detail.Ndf <= 0 || math.IsNaN(*detail.Chi2) || math.IsInf(*detail.Chi2, 0) || *detail.Chi2 < 0 {
		return 0, false
	}
	return *detail.Chi2 / float64(detail.Ndf), true
}

func same_fit(left, right FitDetail) bool {
	return left.Charge == right.Charge && left.Centrality == right.Centrality && left.Bin == right.Bin
}

// A representative good fit has chi2/ndf closest to unity. Prefer a genuinely
// unusable bad example; when all fits are usable show the worst and say so.
func select_fit_examples(details []FitDetail) (*FitDetail, *FitDetail) {
	ordered := append([]FitDetail(nil), details...)
	sort.SliceStable(ordered, func(i, j int) bool {
		left, right := ordered[i], ordered[j]
		if left.Charge != right.Charge {
			return left.Charge < right.Charge
		}
		if left.Centrality != right.Centrality {
			return left.Centrality < right.Centrality
		}
		return left.Bin < right.Bin
	})
	var good, bad *FitDetail
	good_distance := math.Inf(1)
	for i := range ordered {
		detail := &ordered[i]
		value, valid := reduced_chi2(*detail)
		if !detail.Usable || !detail.Finite || detail.Attempts <= 0 || !valid {
			continue
		}
		distance := math.Abs(value - 1)
		if good == nil || distance < good_distance {
			good, good_distance = detail, distance
		}
	}
	worst_value := -1.0
	for i := range ordered {
		detail := &ordered[i]
		if detail.Usable {
			continue
		}
		value, valid := reduced_chi2(*detail)
		if !valid {
			value = -1
		}
		if bad == nil || value > worst_value {
			bad, worst_value = detail, value
		}
	}
	if bad == nil {
		for i := range ordered {
			detail := &ordered[i]
			value, valid := reduced_chi2(*detail)
			if !detail.Usable || !detail.Finite || detail.Attempts <= 0 || !valid || good != nil && same_fit(*good, *detail) {
				continue
			}
			if bad == nil || value > worst_value {
				bad, worst_value = detail, value
			}
		}
	}
	return good, bad
}

func plot_for_path(plots []plot, path string) (plot, bool) {
	if path == "" {
		return plot{}, false
	}
	for _, item := range plots {
		if item.relative == filepath.ToSlash(filepath.Clean(path)) {
			return item, true
		}
	}
	return plot{}, false
}

func plots_for_group(plots []plot, group string) []plot {
	var selected []plot
	for _, slug := range []string{"pos", "neg"} {
		for _, item := range plots {
			stem := strings.TrimSuffix(item.relative, filepath.Ext(item.relative))
			if item.group == group && strings.HasSuffix(stem, "_"+slug) {
				selected = append(selected, item)
				break
			}
		}
	}
	return selected
}

func free_cross_terms(cfg config) []string {
	var free []string
	for _, key := range []string{"r_os", "r_ol", "r_sl"} {
		value, exists := cfg.Fit.Freeze[key]
		if !exists || value == nil {
			free = append(free, key)
		}
	}
	return free
}

func fit_caption(cfg config, detail FitDetail) string {
	charge := map[int]string{0: `$\pi^+$`, 1: `$\pi^-$`}[detail.Charge]
	centrality := map[int]string{0: "0--10", 1: "10--30", 2: "30--50", 3: "50--80"}[detail.Centrality]
	bin := fmt.Sprintf("%d", detail.Bin)
	if detail.Bin >= 0 && detail.Bin < len(cfg.Binning.Names) {
		bin = cfg.Binning.Names[detail.Bin]
	}
	mode := "$y$"
	if cfg.Vars.Input.Type == "kt" {
		mode = "$k_T$"
	}
	return charge + ", центральность " + centrality + `\%, ` + mode + ": " + texEscape(bin)
}

func fit_diagnostic(detail FitDetail) string {
	value := "недоступно"
	if reduced, valid := reduced_chi2(detail); valid {
		value = fmt.Sprintf("%.4g", reduced)
	}
	note := fmt.Sprintf(`$\chi^2/\mathrm{ndf}=%s$; status=%d, cov=%d, попыток=%d`, value, detail.Status, detail.CovStatus, detail.Attempts)
	if !detail.Usable {
		reason := detail.Reason
		if reason == "" {
			switch {
			case detail.Attempts == 0:
				reason = "фит отсутствует"
			case !detail.Finite:
				reason = "неконечные параметры"
			case detail.AtLimit:
				reason = "параметр на границе"
			case detail.Ndf <= 0:
				reason = "неположительное ndf"
			default:
				reason = "результат минимизатора не принят"
			}
		}
		note += `\par\alert{Непригодный фит: ` + texEscape(reason) + `.}`
		if detail.Finite && detail.Attempts > 0 {
			note += ` Пунктирная кривая --- диагностическая попытка.`
		}
	}
	return note
}

func select_compact_slides(data manifest, plots []plot) []CompactSlide {
	good, bad := select_fit_examples(data.Fits.Details)
	slides := make([]CompactSlide, 0, 10)
	for index, detail := range []*FitDetail{good, bad} {
		slide := CompactSlide{Title: "Хороший фит", Kind: "good", Missing: "Пригодный пример с сохранёнными метаданными не найден."}
		if index == 1 {
			slide.Title, slide.Kind = "Плохой фит", "bad"
			slide.Missing = "Непригодный или второй пригодный пример не найден. Плохой фит не создаётся искусственно."
		}
		if detail != nil {
			if index == 1 && detail.Usable {
				slide.Title = "Худший из пригодных фитов"
				slide.Note = "Непригодных фитов нет; выбран максимальный $\\chi^2/\\mathrm{ndf}$ среди остальных пригодных.\\par "
			}
			slide.Note += fit_caption(data.Config, *detail) + `\par ` + fit_diagnostic(*detail)
			if item, exists := plot_for_path(plots, detail.CFImage); exists {
				slide.Plots = []plot{item}
			} else {
				slide.Missing = "Проекция выбранного фита не экспортировалась в этом запуске."
			}
		} else if len(data.Fits.Details) == 0 {
			slide.Missing = "В manifest этого расчёта нет индивидуальных метаданных фитов; качество примера не классифицируется."
			if index == 0 {
				slide.Title = "Проекции КФ: качество не классифицировано"
				for _, item := range plots {
					if item.group == "1d" {
						slide.Plots = []plot{item}
						slide.Note = slide.Missing + `\par ` + item.caption
						break
					}
				}
			}
		}
		slides = append(slides, slide)
	}
	for _, item := range plots_for_group(plots, "radii") {
		slides = append(slides, CompactSlide{Title: "Радиусы и сила корреляции", Kind: "radii", Plots: []plot{item}, Note: item.caption})
	}
	if len(free_cross_terms(data.Config)) > 0 {
		var cross_labels []string
		for _, key := range free_cross_terms(data.Config) {
			cross_labels = append(cross_labels, map[string]string{"r_os": `$R_{os}$`, "r_ol": `$R_{ol}$`, "r_sl": `$R_{sl}$`}[key])
		}
		cross := plots_for_group(plots, "free_cross")
		if len(cross) == 0 {
			cross = plots_for_group(plots, "cross")
		}
		slides = append(slides, CompactSlide{Title: "Свободные перекрёстные параметры", Kind: "cross", Plots: cross,
			Missing: "Свободные перекрёстные параметры не экспортировались.",
			Note:    `Свободные параметры: ` + strings.Join(cross_labels, ", ") + `.`})
	}
	for index, detail := range []*FitDetail{good, bad} {
		title := "Отношение фита к КФ: хороший пример"
		if index == 1 {
			title = "Отношение фита к КФ: плохой пример"
		}
		if index == 1 && detail != nil && detail.Usable {
			title = "Отношение фита к КФ: худший пригодный"
		}
		ratio := CompactSlide{Title: title, Kind: "fit_ratio", Missing: "Для выбранного примера отношение не экспортировалось.",
			Note: `$C_2^{\mathrm{fit}}/C_2$; параметры кривой фиксированы.`}
		if detail != nil {
			if item, exists := plot_for_path(plots, detail.FitOverCFImage); exists {
				ratio.Plots = []plot{item}
				ratio.Note += `\par ` + fit_caption(data.Config, *detail)
				if !detail.Usable {
					ratio.Note += `; \alert{диагностическая попытка}.`
				}
			}
		}
		slides = append(slides, ratio)
	}
	slides = append(slides, CompactSlide{Title: "Качество фитов по всем центральностям", Kind: "quality", Plots: plots_for_group(plots, "quality"), Missing: "Сводные графики качества не экспортировались.", Note: `Слева: $\chi^2/\mathrm{ndf}$. Справа: среднее $C_2^{\mathrm{fit}}/C_2$ по области сравнения. Показаны пригодные фиты.`})
	// This frame is mandatory even for old results: a missing comparison is
	// explicitly reported instead of silently substituting the manual statistic.
	slides = append(slides, CompactSlide{Title: "Сравнение ROOT и ручного chi2", Kind: "chi2_comparison", Plots: plots_for_group(plots, "chi2_comparison"), Missing: "Исходное chi2 ROOT или график сравнения отсутствует в этом расчёте. Нужен новый запуск анализа.", Note: `$\chi^2_{\mathrm{ROOT}}/\chi^2_{\mathrm{manual}}$ для всех попыток с определённым отношением; линия 1 означает совпадение. При $\chi^2_{\mathrm{manual}}=0$ отношение не определено.`})
	slides = append(slides, CompactSlide{Title: "Статистика и статус расчёта", Kind: "statistics"})
	return slides
}

func render_compact_stats(out *strings.Builder, data manifest) {
	status := "Все запрошенные фиты пригодны."
	if data.Status == "incomplete" {
		status = `\alert{Диагностический расчёт: есть непригодные фиты.}`
	}
	fmt.Fprintf(out, "\\small\n%s\\par\\smallskip\n", status)
	fmt.Fprintf(out, "Пригодно %d/%d; непригодно %d; повторных %d; на границах %d.\\par\\medskip\n", data.Fits.Usable, data.Fits.Requested, data.Fits.Unusable, data.Fits.Retried, data.Fits.AtLimit)
	if len(data.Fits.Details) > 0 {
		fmt.Fprint(out, "\\centering\\footnotesize\n\\begin{tabular}{lrrrr}\nВыборка & Всего & Пригодно & Повторно & Граница \\\\ \\hline\n")
		for _, ch := range data.Config.Selection.Charges {
			for _, centr := range data.Config.Selection.Centralities {
				var requested, usable, retried, boundary int
				for _, detail := range data.Fits.Details {
					if detail.Charge == ch && detail.Centrality == centr {
						requested++
						if detail.Usable {
							usable++
						}
						if detail.Attempts > 1 {
							retried++
						}
						if detail.AtLimit {
							boundary++
						}
					}
				}
				if requested > 0 {
					charge := map[int]string{0: `$\pi^+$`, 1: `$\pi^-$`}[ch]
					centrality := map[int]string{0: "0--10", 1: "10--30", 2: "30--50", 3: "50--80"}[centr]
					fmt.Fprintf(out, "%s, %s\\%% & %d & %d & %d & %d \\\\\n", charge, centrality, requested, usable, retried, boundary)
				}
			}
		}
		fmt.Fprint(out, "\\end{tabular}\\par\\medskip\n")
	}
	fmt.Fprintf(out, "{\\scriptsize Вход: %s\\par Завершено (UTC): %s\\par Run ID: \\texttt{%s}\\par}\n", texEscape(filepath.Base(data.Config.Vars.Input.File)), texEscape(data.CompletedAt), texEscape(data.RunID))
	fmt.Fprint(out, "{\\scriptsize Хороший пример: пригодный фит с $\\chi^2/\\mathrm{ndf}$ ближе к 1; плохой: непригодный с максимальным определённым $\\chi^2/\\mathrm{ndf}$. Радиусы и обычные сводные графики содержат только пригодные фиты.}\n")
}

func render_compact(data manifest, slides []CompactSlide) string {
	var out strings.Builder
	if data.Status == "incomplete" {
		fmt.Fprint(&out, "\\def\\ReportDiagnostic{1}\n")
	}
	mode := "$y$"
	if data.Config.Vars.Input.Type == "kt" {
		mode = "$k_T$"
	}
	for _, slide := range slides {
		fmt.Fprintf(&out, "\\begin{frame}{%s: %s}\n", mode, texEscape(slide.Title))
		if slide.Kind == "statistics" {
			render_compact_stats(&out, data)
		} else {
			fmt.Fprint(&out, "\\centering\n")
			if len(slide.Plots) == 0 {
				fmt.Fprintf(&out, "\\small %s\\par\\medskip\n", texEscape(slide.Missing))
			} else if len(slide.Plots) > 1 {
				for _, item := range slide.Plots {
					fmt.Fprintf(&out, "\\begin{minipage}{.49\\linewidth}\\centering\\includegraphics[width=\\linewidth,height=.60\\textheight,keepaspectratio]{%s}\\par{\\scriptsize %s}\\end{minipage}\n", item.asset, item.caption)
				}
			} else {
				height := ".68"
				if slide.Kind == "good" || slide.Kind == "bad" {
					height = ".56"
				}
				if len(slide.Plots) > 1 {
					height = ".29"
				}
				for _, item := range slide.Plots {
					fmt.Fprintf(&out, "\\includegraphics[width=.98\\linewidth,height=%s\\textheight,keepaspectratio]{%s}\\par\n", height, item.asset)
					if len(slide.Plots) > 1 {
						fmt.Fprintf(&out, "{\\scriptsize %s}\\par\n", item.caption)
					}
				}
			}
			if slide.Note != "" {
				fmt.Fprintf(&out, "\\smallskip{\\scriptsize %s}\\par\n", slide.Note)
			}
		}
		fmt.Fprint(&out, "\\end{frame}\n")
	}
	return out.String()
}

func render(data manifest, plots []plot, total int) string {
	diagnostic := data.Status == "incomplete"
	status := "Запрошенные стадии завершены."
	if diagnostic {
		status = `\alert{Диагностический расчёт: есть непригодные фиты.}`
	}
	mode := "$y$"
	if data.Config.Vars.Input.Type == "kt" {
		mode = "$k_t$"
	}
	chargeNames := map[int]string{0: `$\pi^+$`, 1: `$\pi^-$`}
	centralityNames := map[int]string{0: "0-10", 1: "10-30", 2: "30-50", 3: "50-80"}
	charges := make([]string, 0, len(data.Config.Selection.Charges))
	for _, id := range data.Config.Selection.Charges {
		if name, ok := chargeNames[id]; ok {
			charges = append(charges, name)
		}
	}
	centralities := make([]string, 0, len(data.Config.Selection.Centralities))
	for _, id := range data.Config.Selection.Centralities {
		if name, ok := centralityNames[id]; ok {
			centralities = append(centralities, name+`\%`)
		}
	}
	fitLine := fmt.Sprintf("Пригодные фиты: %d из %d; непригодные: %d.", data.Fits.Usable, data.Fits.Requested, data.Fits.Unusable)
	if data.Fits.Requested == 0 {
		fitLine = "Аппроксимация в выбранных стадиях не запрашивалась."
	}
	input := filepath.Base(data.Config.Vars.Input.File)
	var out strings.Builder
	fmt.Fprintf(&out, "\\section{Результаты текущего расчёта}\n")
	fmt.Fprintf(&out, "\\begin{frame}{Параметры и статус расчёта}\n\\small\n%s\n\\begin{itemize}\n", status)
	fmt.Fprintf(&out, "\\item Вход: %s\n", texEscape(input))
	fmt.Fprintf(&out, "\\item Завершение (UTC): %s\n", texEscape(data.CompletedAt))
	fmt.Fprintf(&out, "\\item Биннинг: %s; интервалов: %d.\n", mode, len(data.Config.Binning.Names))
	fmt.Fprintf(&out, "\\item Заряды: %s; центральности: %s.\n", strings.Join(charges, ", "), strings.Join(centralities, ", "))
	fmt.Fprintf(&out, "\\item %s\n", fitLine)
	fmt.Fprintf(&out, "\\item Повторные попытки: %d; параметры на границах: %d.\n", data.Fits.Retried, data.Fits.AtLimit)
	fmt.Fprintf(&out, "\\item Рисунков в презентации: %d из %d сохранённых в этом расчёте.\n", len(plots), total)
	fmt.Fprintf(&out, "\\end{itemize}\n\\end{frame}\n")
	if len(plots) == 0 {
		fmt.Fprintf(&out, "\\begin{frame}{Рисунки расчёта}\nВ этом запуске рисунки не экспортировались.\\par\\medskip Для экспорта включите \\texttt{images.need=true} и выберите формат \\texttt{pdf}.\\end{frame}\n")
	}
	for _, item := range plots {
		title := item.title
		if diagnostic {
			title += " (диагностика)"
		}
		fmt.Fprintf(&out, "\\begin{frame}{%s}\n\\centering\n", texEscape(title))
		fmt.Fprintf(&out, "\\includegraphics[width=.96\\linewidth,height=.72\\textheight,keepaspectratio]{%s}\\par\\smallskip\n", filepath.ToSlash(item.asset))
		fmt.Fprintf(&out, "{\\footnotesize %s}\n\\end{frame}\n", item.caption)
	}
	return out.String()
}

func compilePDF(latexmk, sourceDir, buildDir string, strictAssets bool, mode string) error {
	args := []string{"-cd", "-xelatex", "-interaction=nonstopmode", "-halt-on-error", "-file-line-error"}
	if strictAssets {
		args = append(args, "-g", `-usepretex=\def\ThesisStrictAssets{1}`)
	}
	entry := "presentation/compact.tex"
	if mode == "full" {
		entry = "presentation/main.tex"
	}
	args = append(args, "-outdir="+buildDir, "-jobname=presentation", entry)
	command := exec.Command(latexmk, args...)
	command.Dir = sourceDir
	command.Stdout = os.Stdout
	command.Stderr = os.Stderr
	env := replaceEnv(os.Environ(), "TEXINPUTS", buildDir+"//:")
	command.Env = env
	if err := command.Run(); err != nil {
		return fmt.Errorf("latexmk: %w", err)
	}
	return nil
}

func replaceEnv(env []string, key, value string) []string {
	prefix := key + "="
	for i, item := range env {
		if strings.HasPrefix(item, prefix) {
			return append(env[:i], append([]string{prefix + value + strings.TrimPrefix(item, prefix)}, env[i+1:]...)...)
		}
	}
	return append(env, prefix+value)
}

func copyAtomic(source, destination string) error {
	if err := os.MkdirAll(filepath.Dir(destination), 0o755); err != nil {
		return err
	}
	temporary := destination + ".tmp"
	if err := copyFile(source, temporary); err != nil {
		_ = os.Remove(temporary)
		return err
	}
	if err := os.Rename(temporary, destination); err != nil {
		_ = os.Remove(temporary)
		return err
	}
	return nil
}

func copyFile(source, destination string) error {
	input, err := os.Open(source)
	if err != nil {
		return err
	}
	defer input.Close()
	output, err := os.Create(destination)
	if err != nil {
		return err
	}
	if _, err := io.Copy(output, input); err != nil {
		_ = output.Close()
		return err
	}
	return output.Close()
}

func writeAtomic(destination string, content []byte) error {
	if err := os.MkdirAll(filepath.Dir(destination), 0o755); err != nil {
		return err
	}
	temporary := destination + ".tmp"
	if err := os.WriteFile(temporary, content, 0o644); err != nil {
		return err
	}
	if err := os.Rename(temporary, destination); err != nil {
		_ = os.Remove(temporary)
		return err
	}
	return nil
}

func inside(root, path string) bool {
	relative, err := filepath.Rel(root, path)
	return err == nil && relative != ".." && !strings.HasPrefix(relative, ".."+string(filepath.Separator))
}

func safeRunID(value string) string {
	var out strings.Builder
	for _, char := range value {
		if char >= 'a' && char <= 'z' || char >= 'A' && char <= 'Z' || char >= '0' && char <= '9' || char == '-' || char == '_' {
			out.WriteRune(char)
		} else {
			out.WriteByte('_')
		}
	}
	if out.Len() == 0 {
		return time.Now().UTC().Format("20060102T150405Z")
	}
	return out.String()
}

func short(value string, limit int) string {
	runes := []rune(value)
	if len(runes) <= limit {
		return value
	}
	return string(runes[:limit-1]) + "…"
}

func texEscape(value string) string {
	replacements := strings.NewReplacer(
		`\`, `\textbackslash{}`, "{", `\{`, "}", `\}`,
		"$", `\$`, "&", `\&`, "#", `\#`, "%", `\%`,
		"_", `\_`, "~", `\textasciitilde{}`, "^", `\textasciicircum{}`,
		"\n", " ", "\r", " ",
	)
	return replacements.Replace(value)
}
