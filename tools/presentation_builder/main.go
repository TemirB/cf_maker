package main

import (
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
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
}

type manifest struct {
	SchemaVersion int    `json:"schema_version"`
	Status        string `json:"status"`
	RunID         string `json:"run_id"`
	CompletedAt   string `json:"completed_at"`
	ExitCode      *int   `json:"exit_code"`
	Config        config `json:"config"`
	Fits          struct {
		Requested int `json:"requested"`
		Usable    int `json:"usable"`
		Unusable  int `json:"unusable"`
		Retried   int `json:"retried"`
		AtLimit   int `json:"at_limit"`
	} `json:"fits"`
	Images []struct {
		Path string `json:"path"`
		Size int64  `json:"size_bytes"`
	} `json:"images"`
}

type plot struct {
	source  string
	asset   string
	group   string
	title   string
	caption string
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
		if *resultsName != "" || analysisExit != nil || !errors.Is(manifestErr, os.ErrNotExist) {
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
		if err := writeChapter(buildDir, results, data, digest, *maxPlots, chapterPath); err != nil {
			return err
		}
	}

	if err := compilePDF(*latexmk, sourceDir, buildDir, *strictAssets); err != nil {
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

func writeChapter(buildDir, results string, data manifest, digest string, maxPlots int, chapterPath string) error {
	plots, err := collectPlots(results, data)
	if err != nil {
		return err
	}
	selected := selectPlots(plots, maxPlots)
	assetDir := filepath.Join(buildDir, "current-results", safeRunID(data.RunID))
	if err := os.MkdirAll(assetDir, 0o755); err != nil {
		return err
	}
	for i := range selected {
		name := fmt.Sprintf("plot-%04d%s", i+1, strings.ToLower(filepath.Ext(selected[i].source)))
		destination := filepath.Join(assetDir, name)
		if err := copyFile(selected[i].source, destination); err != nil {
			return err
		}
		selected[i].asset, err = filepath.Rel(buildDir, destination)
		if err != nil {
			return err
		}
		selected[i].asset = filepath.ToSlash(selected[i].asset)
	}
	if _, latestDigest, err := readManifest(results); err != nil || latestDigest != digest {
		return errors.New("analysis results changed while preparing slides")
	}
	text := render(data, selected, len(plots))
	return writeAtomic(chapterPath, []byte(text))
}

func collectPlots(results string, data manifest) ([]plot, error) {
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
		if !inside(results, resolved) {
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
	order := []string{"radii", "quality", "1d", "2d", "pvalues", "fit_ratio", "cross", "other"}
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

func compilePDF(latexmk, sourceDir, buildDir string, strictAssets bool) error {
	args := []string{"-cd", "-xelatex", "-interaction=nonstopmode", "-halt-on-error", "-file-line-error"}
	if strictAssets {
		args = append(args, "-g", `-usepretex=\def\ThesisStrictAssets{1}`)
	}
	args = append(args, "-outdir="+buildDir, "-jobname=presentation", "presentation/main.tex")
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
