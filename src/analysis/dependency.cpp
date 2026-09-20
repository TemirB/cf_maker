#include "analysis/dependency.h"

#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TString.h>

#include "analysis/graphs.h"
#include "core/binning.h"
#include "core/fs.h"
#include "core/log.h"
#include "draw/draw.h"
#include "fit/types.h"

void make_dependency(Config& cfg, TFile* cf3dFile, TFile* outFile)
{
    const char* mode = cfg.input.type.c_str();

    const FitGrid& fit_results = cfg.fit_results;
    const Bin& bin = cfg.binning;

    std::string dir = cfg.output.dir + "/dependency";
    ensure_dir(dir);
    log::Info("dependency: output dir = " + dir);

    const std::string ext = cfg.general.images.format;
    for (const int ch : cfg.selection.charges) {
        log::Debug("dependency: charge = " + std::string(charge::kNames[ch]));
        std::array<std::unique_ptr<TMultiGraph>, lcms::kCount> mg_radii;
        for (auto& mg : mg_radii) {
            mg = std::make_unique<TMultiGraph>();
        }
        auto mg_lambda = std::make_unique<TMultiGraph>();
        auto mg_chi2_ndf = std::make_unique<TMultiGraph>();
        auto mg_pvalue = std::make_unique<TMultiGraph>();
        auto mg_fit_over_cf = std::make_unique<TMultiGraph>();
        std::vector<std::pair<TObject*, std::string>> legend_entries;

        for (const int centr : cfg.selection.centralities) {
            std::array<TGraphErrors*, lcms::kCount> g_radii{};
            for (int lcms = 0; lcms < lcms::kCount; lcms++) {
                g_radii[lcms] = make_styled_graph(Form("g_R_%s_%s_centr_%s", lcms::kNames[lcms],
                                                 charge::kNames[ch], centrality::kNames[centr]),
                                            centr);
            }
            TGraphErrors* g_lambda = make_styled_graph(
                Form("g_L_%s_centr_%s", charge::kNames[ch], centrality::kNames[centr]), centr);
            TGraphErrors* g_chi2_ndf = build_chi2_ndf_graph(cfg, ch, centr);
            TGraphErrors* g_pvalue = build_pvalue_graph(cfg, ch, centr);
            TGraphErrors* g_fit_over_cf = build_fit_over_cf_graph(cfg, cf3dFile, ch, centr);

            legend_entries.emplace_back(g_radii[0], centrality::kNames[centr]);

            for (int b = 0; b < bin.count; b++) {
                const FitResult& res = fit_results[ch][centr][b];

                double x_val = bin_center(bin, b);

                for (int lcms = 0; lcms < lcms::kCount; lcms++) {
                    g_radii[lcms]->SetPoint(b, x_val, res.r[lcms]);
                    g_radii[lcms]->SetPointError(b, 0, res.e_r[lcms]);
                }

                g_lambda->SetPoint(b, x_val, res.lambda);
                g_lambda->SetPointError(b, 0, res.e_lambda);
            }

            for (int lcms = 0; lcms < lcms::kCount; lcms++) {
                mg_radii[lcms]->Add(g_radii[lcms], "lp");
            }

            mg_lambda->Add(g_lambda, "lp");
            mg_chi2_ndf->Add(g_chi2_ndf, "lp");
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
        write_mg_with_legend(outFile, mg_lambda.get(), mg_lambda->GetName(), mode, "lambda", legend_entries,
                          draw::GraphKind::Lambda);

        set_range_with_errors(mg_chi2_ndf.get(), 0.1);
        mg_chi2_ndf->SetName(Form("mg_chi2_ndf_%s", charge::kNames[ch]));
        write_mg_with_legend(outFile, mg_chi2_ndf.get(), mg_chi2_ndf->GetName(), mode, "#chi^{2}/ndf",
                          legend_entries, draw::GraphKind::chi2_ndf);

        set_range_with_errors(mg_fit_over_cf.get(), 0.1);
        mg_fit_over_cf->SetName(Form("mg_FitOverCF_%s", charge::kNames[ch]));
        write_mg_with_legend(outFile, mg_fit_over_cf.get(), mg_fit_over_cf->GetName(), mode, "<fit/CF>",
                          legend_entries, draw::GraphKind::FitOverCF);

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
                mg_radii[lcms]->Draw("APL");
            }

            c->cd(4);
            mg_lambda->Draw("APL");
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
            mg_chi2_ndf->Draw("APL");
            mg_chi2_ndf->GetXaxis()->SetTitle(mode);
            mg_chi2_ndf->GetYaxis()->SetTitle("#chi^{2}/ndf");

            c->cd(2);
            mg_fit_over_cf->Draw("APL");
            mg_fit_over_cf->GetXaxis()->SetTitle(mode);
            mg_fit_over_cf->GetYaxis()->SetTitle("<fit/CF>");

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
            std::string name = Form("c_all_cross_graphs_%s", charge::kNames[ch]);
            std::string title = Form("Graphs for cross radii(%s)", charge::kNames[ch]);
            auto c = std::make_unique<TCanvas>(name.c_str(), title.c_str(), 3200, 1600);
            c->Divide(2, 2);
            int idx = 1;
            for (int lcms = 0; lcms < 3; lcms++) {
                c->cd(idx);
                mg_radii[lcms + 3]->Draw("APL");
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
            mg_pvalue->Draw("APL");

            {
                outFile->cd();
                mg_pvalue->Write();
            }
            std::string save_name = Form("%s/%s.%s", dir.data(), name.data(), ext.data());
            if (cfg.general.images.need) {
                save_canvas_quiet(c.get(), save_name.data());
            }
        }
    }
}
