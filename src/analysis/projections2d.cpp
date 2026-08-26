#include "analysis/projections2d.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TH3.h>
#include <TPad.h>
#include <TString.h>

#include "core/binning.h"
#include "core/fs.h"
#include "core/lcms.h"
#include "core/log.h"
#include "io/input.h"

template <class T> using RootPtr = std::unique_ptr<T>;

namespace
{

void crop_2d(TH2D& h, double q = 0.08)
{
    h.GetXaxis()->SetRangeUser(-q, q);
    h.GetYaxis()->SetRangeUser(-q, q);
}

std::size_t canvas_index(int ch, int centr)
{
    return static_cast<std::size_t>(ch) * centrality::kCount + static_cast<std::size_t>(centr);
}

void write_2d_projection(TFile& outFile, const TH3D& den_source, const TH3D& num_source, LCMSAxis ax1,
                       LCMSAxis ax2, const std::string& tag, TCanvas* canvas, const int y,
                       bool draw, double freezeWidth, double cropWidth,
                       std::vector<std::unique_ptr<TH2D>>& keep_alive)
{
    auto den = RootPtr<TH3D>(static_cast<TH3D*>(den_source.Clone()));
    auto num = RootPtr<TH3D>(static_cast<TH3D*>(num_source.Clone()));
    den->SetDirectory(nullptr);
    num->SetDirectory(nullptr);

    auto num = RootPtr<TH2D>(project_2d(*num, ax1, ax2, freezeWidth));
    auto den = RootPtr<TH2D>(project_2d(*den, ax1, ax2, freezeWidth));
    if (!num || !den) {
        std::cerr << "Project3D failed for " << tag << "\n";
        return;
    }

    auto name = tag + " " + axis_name(ax1) + "-" + axis_name(ax2);
    auto cf_hist = RootPtr<TH2D>(static_cast<TH2D*>(num->Clone(name.c_str())));
    cf_hist->Divide(num.get(), den.get());

    crop_2d(*cf_hist, cropWidth);

    cf_hist->SetTitle(name.data());
    cf_hist->GetXaxis()->SetTitle(("q_{" + axis_name(ax1) + "} [GeV/c]").c_str());
    cf_hist->GetYaxis()->SetTitle(("q_{" + axis_name(ax2) + "} [GeV/c]").c_str());
    cf_hist->GetZaxis()->SetRangeUser(0.95, 1.30);
    cf_hist->SetStats(kFALSE);

    TCanvas c((name).c_str(), "", 650, 600);
    cf_hist->Draw("COLZ");

    outFile.cd();
    c.Write();

    if (draw) {
        canvas->cd(y + 1);
        gPad->SetTicks(1, 1);
        gPad->SetLeftMargin(0.12f);
        gPad->SetBottomMargin(0.12f);

        auto cf_clone =
            std::unique_ptr<TH2D>(static_cast<TH2D*>(cf_hist->Clone(Form("%s_clone", cf_hist->GetName()))));
        cf_clone->Draw("COLZ");

        gPad->Modified();
        gPad->Update();
        keep_alive.push_back(std::move(cf_clone));
    }
}

} // namespace

void make_lcms_2d_projections(const Config& cfg, TFile* in, TFile* out)
{
    const Bin& bin = cfg.binning;
    std::string dir = cfg.output.dir + "/all_2d_histos";
    ensure_dir(dir);
    log::Info("projections_2d: " + std::to_string(cfg.selection.charges.size()) + " charges x " +
              std::to_string(cfg.selection.centralities.size()) + " centralities x " +
              std::to_string(bin.count) + " bins");

    const std::string ext = cfg.general.images.format;
    std::vector<std::unique_ptr<TCanvas>> canvases(static_cast<std::size_t>(charge::kCount) *
                                                   centrality::kCount);
    std::vector<std::unique_ptr<TH2D>> keep_alive;
    for (const int ch : cfg.selection.charges) {
        for (const int centr : cfg.selection.centralities) {
            auto name = Form("all_out-long_2d_histos_centr_%s_%s", centrality::kNames[centr],
                             charge::kNames[ch]);
            auto title = Form("CF_{out-long} at ch=%s, centrality=%s", charge::kNames[ch],
                              centrality::kNames[centr]);
            canvases[canvas_index(ch, centr)] = std::make_unique<TCanvas>(name, title, 2000, 800);
            canvases[canvas_index(ch, centr)]->Divide(5, 2);
        }
    }

    for (const int ch_idx : cfg.selection.charges) {
        for (const int cent_idx : cfg.selection.centralities) {
            for (int b = 0; b < bin.count; b++) {
                auto [den, num] = get_hists(in, ch_idx, cent_idx, b);
                if (!den || !num) {
                    continue;
                }

                std::string cf_name = get_cf_name(ch_idx, cent_idx, cfg.input.type, bin.names[b]);

                write_2d_projection(*out, *den, *num, LCMSAxis::Out, LCMSAxis::Side, cf_name,
                                  canvases[canvas_index(ch_idx, cent_idx)].get(), b, false,
                                  cfg.projections.slice_2d, cfg.projections.crop_2d, keep_alive);
                write_2d_projection(*out, *den, *num, LCMSAxis::Out, LCMSAxis::Long, cf_name,
                                  canvases[canvas_index(ch_idx, cent_idx)].get(), b, true,
                                  cfg.projections.slice_2d, cfg.projections.crop_2d, keep_alive);
                write_2d_projection(*out, *den, *num, LCMSAxis::Side, LCMSAxis::Long, cf_name,
                                  canvases[canvas_index(ch_idx, cent_idx)].get(), b, false,
                                  cfg.projections.slice_2d, cfg.projections.crop_2d, keep_alive);
            }
        }
    }

    for (const int ch : cfg.selection.charges) {
        for (const int centr : cfg.selection.centralities) {
            auto name = Form("%s/all_out-long_2d_histos_centr_%s_%s.%s", dir.c_str(),
                             centrality::kNames[centr], charge::kNames[ch], ext.c_str());
            if (cfg.general.images.need) {
                save_canvas_quiet(canvases[canvas_index(ch, centr)].get(), name);
            }
        }
    }
}
