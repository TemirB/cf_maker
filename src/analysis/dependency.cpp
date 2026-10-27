#include "analysis/dependency.h"

#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TString.h>

#include "analysis/graphs.h"
#include "core/binning.h"
#include "core/fs.h"
#include "core/log.h"
#include "draw/draw.h"
#include "fit/types.h"
#include "io/output.h"

namespace
{
void draw_summary_legend(const std::vector<std::pair<TObject*, std::string>>& entries)
{
    TLegend legend(0.68, 0.68, 0.89, 0.89);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextSize(0.035f);
    for (const auto& [object, label] : entries) {
        legend.AddEntry(object, label.c_str(), "lp");
    }
    legend.DrawClone();
}
} // namespace

void make_dependency(Config& cfg, TFile* cf3dFile, TFile* outFile)
{
    const char* mode = cfg.input.type.c_str();

    const FitGrid& fit_results = cfg.fit_results;
    const Bin& bin = cfg.binning;

    std::string dir = cfg.output.dir + "/dependency";
    ensure_dir(dir);
    logging::info("dependency: output dir = " + dir);

    const std::string ext = cfg.general.images.format;
    for (const int ch : cfg.selection.charges) {
        logging::debug("dependency: charge = " + std::string(charge::kNames[ch]));
        std::array<std::unique_ptr<TMultiGraph>, lcms::kCount> mg_radii;
        for (auto& mg : mg_radii) {
            mg = std::make_unique<TMultiGraph>();
        }
        auto mg_lambda = std::make_unique<TMultiGraph>();
        auto mg_chi2_ndf = std::make_unique<TMultiGraph>();
        auto mg_chi2_comparison = std::make_unique<TMultiGraph>();
        auto mg_pvalue = std::make_unique<TMultiGraph>();
        auto mg_fit_over_cf = std::make_unique<TMultiGraph>();
        std::vector<std::pair<TObject*, std::string>> legend_entries;
        std::vector<std::pair<TObject*, std::string>> comparison_legend;

        for (const int centr : cfg.selection.centralities) {
            std::array<TGraphErrors*, lcms::kCount> g_radii{};
            for (int lcms = 0; lcms < lcms::kCount; lcms++) {
                g_radii[lcms] =
                    make_styled_graph(Form("g_R_%s_%s_centr_%s", lcms::kNames[lcms],
                                           charge::kNames[ch], centrality::kNames[centr]),
                                      centr);
            }
            TGraphErrors* g_lambda = make_styled_graph(
                Form("g_L_%s_centr_%s", charge::kNames[ch], centrality::kNames[centr]), centr);
            TGraphErrors* g_chi2_ndf = build_chi2_ndf_graph(cfg, ch, centr);
            TGraphErrors* g_chi2_comparison = build_chi2_root_over_manual_graph(cfg, ch, centr);
            TGraphErrors* g_pvalue = build_pvalue_graph(cfg, ch, centr);
            TGraphErrors* g_fit_over_cf = build_fit_over_cf_graph(cfg, cf3dFile, ch, centr);

            int point = 0;
            for (int b = 0; b < bin.count; b++) {
                const FitResult& res = fit_results[ch][centr][b];
                if (!is_usable_fit(res)) {
                    logging::debug(
                        "dependency: excluding unavailable fit ch=" + std::to_string(ch) +
                        " centr=" + std::to_string(centr) + " bin=" + std::to_string(b));
                    continue;
                }

                double x_val = bin_center(bin, b);

                for (int lcms = 0; lcms < lcms::kCount; lcms++) {
                    g_radii[lcms]->SetPoint(point, x_val, res.r[lcms]);
                    g_radii[lcms]->SetPointError(point, 0, res.e_r[lcms]);
                }

                g_lambda->SetPoint(point, x_val, res.lambda);
                g_lambda->SetPointError(point, 0, res.e_lambda);
                ++point;
            }
            if (point > 0) {
                legend_entries.emplace_back(g_radii[0], centrality::kNames[centr]);
            }

            for (int lcms = 0; lcms < lcms::kCount; lcms++) {
                mg_radii[lcms]->Add(g_radii[lcms], "lp");
            }

            mg_lambda->Add(g_lambda, "lp");
            mg_chi2_ndf->Add(g_chi2_ndf, "lp");
            mg_chi2_comparison->Add(g_chi2_comparison, "lp");
            if (g_chi2_comparison->GetN() > 0) {
                comparison_legend.emplace_back(g_chi2_comparison, centrality::kNames[centr]);
            }
            mg_fit_over_cf->Add(g_fit_over_cf, "lp");
            mg_pvalue->Add(g_pvalue, "lp");
        }

        for (int lcms = 0; lcms < lcms::kCount; lcms++) {
            mg_radii[lcms]->SetName(Form("mg_R_%s_%s", lcms::kNames[lcms], charge::kNames[ch]));
            set_range_with_errors(mg_radii[lcms].get(), 0.1);
            draw::GraphKind kind = draw::GraphKind::Radii;
            if (lcms >= 3) {
                kind = draw::GraphKind::Cross;
            }
            write_mg_with_legend(outFile, mg_radii[lcms].get(), mg_radii[lcms]->GetName(), mode,
                                 Form("R_{%s} (fm)", lcms::kNames[lcms]), legend_entries, kind);
        }

        set_range_with_errors(mg_lambda.get(), 0.1);
        mg_lambda->SetName(Form("mg_L_%s", charge::kNames[ch]));
        write_mg_with_legend(outFile, mg_lambda.get(), mg_lambda->GetName(), mode, "lambda",
                             legend_entries, draw::GraphKind::Lambda);

        set_range_with_errors(mg_chi2_ndf.get(), 0.1);
        mg_chi2_ndf->SetName(Form("mg_chi2_ndf_%s", charge::kNames[ch]));
        write_mg_with_legend(outFile, mg_chi2_ndf.get(), mg_chi2_ndf->GetName(), mode,
                             "#chi^{2}/ndf", legend_entries, draw::GraphKind::chi2_ndf);

        set_range_with_errors(mg_chi2_comparison.get(), 0.1);
        mg_chi2_comparison->SetName(Form("mg_chi2_root_over_manual_%s", charge::kNames[ch]));
        write_mg_with_legend(outFile, mg_chi2_comparison.get(), mg_chi2_comparison->GetName(), mode,
                             "#chi^{2}_{ROOT}/#chi^{2}_{manual}", comparison_legend,
                             draw::GraphKind::Chi2RootOverManual);
        if (cfg.general.images.need) {
            auto* comparison_canvas =
                dynamic_cast<TCanvas*>(outFile->Get(mg_chi2_comparison->GetName()));
            if (comparison_canvas) {
                const std::string image =
                    dir + "/c_chi2_root_over_manual_" + charge::kNames[ch] + "." + ext;
                save_canvas_quiet(comparison_canvas, image.c_str());
            }
        }

        set_range_with_errors(mg_fit_over_cf.get(), 0.1);
        mg_fit_over_cf->SetName(Form("mg_FitOverCF_%s", charge::kNames[ch]));
        write_mg_with_legend(outFile, mg_fit_over_cf.get(), mg_fit_over_cf->GetName(), mode,
                             "<C_{2}^{fit}/C_{2}>", legend_entries, draw::GraphKind::FitOverCF);

        mg_pvalue->SetName(Form("mg_Pvalue_%s", charge::kNames[ch]));
        write_mg_with_legend(outFile, mg_pvalue.get(), mg_pvalue->GetName(), mode, "p_{value}",
                             legend_entries, draw::GraphKind::PValue);

        {
            std::string name = Form("c_all_graphs_%s", charge::kNames[ch]);
            std::string title = Form("Graphs for radii and #lambda(%s)", charge::kNames[ch]);
            auto c = std::make_unique<TCanvas>(name.c_str(), title.c_str(), 3200, 1600);
            c->Divide(2, 2);
            for (int lcms = 0; lcms < 3; lcms++) {
                c->cd(lcms + 1);
                draw_mg_or_report(mg_radii[lcms].get());
                draw_summary_legend(legend_entries);
            }

            c->cd(4);
            if (auto* axis = mg_lambda->GetYaxis()) {
                axis->UnZoom();
            }
            set_range_with_errors(mg_lambda.get());
            draw_mg_or_report(mg_lambda.get());
            draw_summary_legend(legend_entries);
            std::string save_name = dir;
            save_name += "/";
            save_name += name;
            save_name += ".";
            save_name += ext;
            if (cfg.general.images.need) {
                save_canvas_quiet(c.get(), save_name.c_str());
            }
        }

        {
            std::string name = Form("c_fit_quality_%s", charge::kNames[ch]);
            std::string title = Form("Fit quality(%s)", charge::kNames[ch]);
            auto c = std::make_unique<TCanvas>(name.c_str(), title.c_str(), 2400, 1000);
            c->Divide(2, 1);

            c->cd(1);
            gPad->SetLogy();
            if (auto* axis = mg_chi2_ndf->GetYaxis()) {
                axis->UnZoom();
            }
            set_range_with_errors(mg_chi2_ndf.get());
            if (draw_mg_or_report(mg_chi2_ndf.get())) {
                mg_chi2_ndf->GetXaxis()->SetTitle(mode);
                mg_chi2_ndf->GetYaxis()->SetTitle("#chi^{2}/ndf");
            }
            draw_summary_legend(legend_entries);

            c->cd(2);
            if (draw_mg_or_report(mg_fit_over_cf.get())) {
                mg_fit_over_cf->GetXaxis()->SetTitle(mode);
                mg_fit_over_cf->GetYaxis()->SetTitle("<C_{2}^{fit}/C_{2}>");
            }
            draw_summary_legend(legend_entries);

            std::string save_name = dir;
            save_name += "/";
            save_name += name;
            save_name += ".";
            save_name += ext;
            if (cfg.general.images.need) {
                save_canvas_quiet(c.get(), save_name.c_str());
            }
        }

        if (cfg.general.images.need) {
            std::vector<int> free_cross;
            for (int parameter = 3; parameter < 6; ++parameter) {
                if (!cfg.fit.freeze[parameter]) {
                    free_cross.push_back(parameter);
                }
            }
            if (!free_cross.empty()) {
                const std::string name = std::string("c_free_cross_graphs_") + charge::kNames[ch];
                TCanvas canvas(name.c_str(), "Free cross parameters",
                               1000 * static_cast<int>(free_cross.size()), 800);
                canvas.Divide(static_cast<int>(free_cross.size()), 1);
                for (std::size_t index = 0; index < free_cross.size(); ++index) {
                    canvas.cd(static_cast<int>(index) + 1);
                    auto* graph = mg_radii[free_cross[index]].get();
                    if (auto* axis = graph->GetYaxis()) {
                        axis->UnZoom();
                    }
                    set_range_with_errors(graph);
                    draw_mg_or_report(graph);
                    draw_summary_legend(legend_entries);
                }
                save_canvas_quiet(&canvas, (dir + "/" + name + "." + ext).c_str());
            }
        }

        {
            std::string name = Form("c_all_cross_graphs_%s", charge::kNames[ch]);
            std::string title = Form("Graphs for cross radii(%s)", charge::kNames[ch]);
            auto c = std::make_unique<TCanvas>(name.c_str(), title.c_str(), 3200, 1600);
            c->Divide(2, 2);
            int idx = 1;
            for (int lcms = 0; lcms < 3; lcms++) {
                c->cd(idx);
                draw_mg_or_report(mg_radii[lcms + 3].get());
                idx++;
            }

            std::string save_name = dir;
            save_name += "/";
            save_name += name;
            save_name += ".";
            save_name += ext;
            if (cfg.general.images.need) {
                save_canvas_quiet(c.get(), save_name.c_str());
            }
        }

        {
            std::string name = Form("c_pvalues_%s", charge::kNames[ch]);
            std::string title = Form("P-value %s", charge::kNames[ch]);

            auto c = std::make_unique<TCanvas>(name.data(), title.data(), 1600, 1600);
            draw_mg_or_report(mg_pvalue.get());

            {
                write_output_object(*outFile, *mg_pvalue);
            }
            std::string save_name = Form("%s/%s.%s", dir.data(), name.data(), ext.data());
            if (cfg.general.images.need) {
                save_canvas_quiet(c.get(), save_name.data());
            }
        }
    }
}
