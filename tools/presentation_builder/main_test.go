package main

import (
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestPlotsInsideSymlinkedResultsDirectory(t *testing.T) {
	base := t.TempDir()
	realDir, alias := filepath.Join(base, "real"), filepath.Join(base, "alias")
	if err := os.Mkdir(realDir, 0700); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(filepath.Join(realDir, "plot.pdf"), []byte("PDF"), 0600); err != nil {
		t.Fatal(err)
	}
	if err := os.Symlink(realDir, alias); err != nil {
		t.Skip(err)
	}
	var data manifest
	if err := json.Unmarshal([]byte(`{"images":[{"path":"plot.pdf","size_bytes":3}]}`), &data); err != nil {
		t.Fatal(err)
	}
	if plots, err := collectPlots(alias, data); err != nil || len(plots) != 1 {
		t.Fatalf("valid results alias rejected: %v", err)
	}
	outside := filepath.Join(base, "outside.pdf")
	if err := os.WriteFile(outside, []byte("PDF"), 0600); err != nil {
		t.Fatal(err)
	}
	if err := os.Symlink(outside, filepath.Join(realDir, "escape.pdf")); err != nil {
		t.Fatal(err)
	}
	data.Images[0].Path = "escape.pdf"
	if _, err := collectPlots(alias, data); err == nil {
		t.Fatal("image escaping results directory accepted")
	}
}

func fit_detail(ch, centr, bin int, chi2 float64, usable bool) FitDetail {
	return FitDetail{Charge: ch, Centrality: centr, Bin: bin, Chi2: &chi2, Ndf: 10,
		Usable: usable, Finite: true, Attempts: 1,
		CFImage: "cf-example.pdf", FitOverCFImage: "ratio-example.pdf"}
}

func compact_fixture() manifest {
	var data manifest
	data.Status = "incomplete"
	data.Config.Vars.Input.Type = "kt"
	data.Config.Selection.Charges = []int{0, 1}
	data.Config.Selection.Centralities = []int{0, 1, 2, 3}
	zero := 0.0
	data.Config.Fit.Freeze = map[string]*float64{"r_os": &zero, "r_sl": &zero}
	data.Fits.Details = []FitDetail{fit_detail(0, 0, 0, 11, true), fit_detail(1, 3, 2, 900, false)}
	data.Fits.Requested, data.Fits.Usable, data.Fits.Unusable = 2, 1, 1
	return data
}

func TestRepresentativeFitsUseQualityAndDeterministicTies(t *testing.T) {
	details := []FitDetail{fit_detail(1, 2, 1, 10, true), fit_detail(0, 1, 1, 10, true),
		fit_detail(0, 2, 0, 400, false), fit_detail(1, 3, 0, 900, false)}
	good, bad := select_fit_examples(details)
	if good == nil || good.Charge != 0 || good.Centrality != 1 || bad == nil || bad.Charge != 1 {
		t.Fatalf("unexpected representatives: good=%+v bad=%+v", good, bad)
	}
	details[0].Chi2 = nil
	details[1].Chi2 = new(float64)
	*details[1].Chi2 = math.NaN()
	good, _ = select_fit_examples(details)
	if good != nil {
		t.Fatal("missing or nonfinite statistic classified as good")
	}
}

func TestAllUsableUsesHonestlyLabelledWorstDistinctExample(t *testing.T) {
	data := compact_fixture()
	data.Status = "completed"
	data.Fits.Details[1].Usable = true
	slides := select_compact_slides(data, nil)
	if slides[1].Title != "Худший из пригодных фитов" || !strings.Contains(slides[1].Note, "Непригодных фитов нет") {
		t.Fatalf("usable fit mislabeled as bad: %+v", slides[1])
	}
	good, bad := select_fit_examples(data.Fits.Details[:1])
	if good == nil || bad != nil {
		t.Fatal("same fit duplicated as good and bad")
	}
}

func TestCompactDeckBoundsPagesAndAlwaysIncludesChi2Ratio(t *testing.T) {
	data := compact_fixture()
	labels := makeLabels(data.Config)
	var plots []plot
	for path, item := range labels {
		item.relative, item.asset = path+".pdf", "test.pdf"
		plots = append(plots, item)
	}
	for i := 0; i < 500; i++ {
		plots = append(plots, plot{group: "other", asset: "test.pdf"})
	}
	for _, mode := range []string{"kt", "rapidity"} {
		data.Config.Vars.Input.Type = mode
		slides := select_compact_slides(data, plots)
		if len(slides) != 10 {
			t.Fatalf("%s got %d pages, want 10", mode, len(slides))
		}
		comparisons := 0
		for _, slide := range slides {
			if slide.Kind == "chi2_comparison" {
				comparisons++
				if len(slide.Plots) != 2 {
					t.Fatal("one of the two charges is missing from chi2 ratio")
				}
			}
		}
		tex := render_compact(data, slides)
		if comparisons != 1 || strings.Count(tex, `\begin{frame}`) != 10 || !strings.Contains(tex, `\chi^2_{\mathrm{ROOT}}/\chi^2_{\mathrm{manual}}`) {
			t.Fatal("compact report lost mandatory chi2 comparison or repeated frames")
		}
	}
}

func TestOldManifestDoesNotInventFitQualityOrRootStatistic(t *testing.T) {
	var data manifest
	slides := select_compact_slides(data, nil)
	if strings.Contains(slides[0].Title, "Хороший") || !strings.Contains(slides[0].Missing, "метаданных") {
		t.Fatal("legacy manifest was classified without fit metadata")
	}
	for _, slide := range slides {
		if slide.Kind == "chi2_comparison" && (len(slide.Plots) != 0 || !strings.Contains(slide.Missing, "отсутствует")) {
			t.Fatal("missing ROOT chi2 ratio was fabricated")
		}
	}
}

func TestCompactTemplateIncludesOnlyCurrentResultsOnce(t *testing.T) {
	raw, err := os.ReadFile("../../docs/thesis/presentation/compact.tex")
	if err != nil {
		t.Fatal(err)
	}
	if strings.Count(string(raw), `\input{presentation-results.tex}`) != 1 || strings.Contains(string(raw), `\input{content}`) {
		t.Fatal("compact entry point includes old slides or duplicates results")
	}
	raw, err = os.ReadFile("../../docs/thesis/presentation/content.tex")
	if err != nil {
		t.Fatal(err)
	}
	if strings.Count(string(raw), `\input{\jobname-results.tex}`) != 1 || strings.Contains(string(raw), `\InputIfFileExists{presentation-results.tex}`) {
		t.Fatal("full entry point duplicates generated results")
	}
}
