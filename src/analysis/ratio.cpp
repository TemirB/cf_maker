#include "analysis/ratio.h"

#include <memory>

#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TString.h>

#include "core/binning.h"
#include "core/lcms.h"
#include "core/log.h"
#include "io/input.h"

namespace
{

void Style(TH1D* h, const TString& fmt, LCMSAxis axis, int centr, const std::string& binName)
{
    {
        TString axisStr = axis_name(axis);
        TString title =
            TString::Format("%s at centrality=[%s] %%, bin=[%s] and %s axis", fmt.Data(),
                            centrality::kNames[centr], binName.c_str(), axisStr.Data());

        h->SetTitle(title);

        h->SetMarkerStyle(kFullCircle);
        h->SetMarkerSize(1.0);
        h->SetMarkerColor(kBlack);
        h->SetLineColor(kBlack);
        h->SetLineWidth(2);
    }

    {
        TAxis* xAxis = h->GetXaxis();
        TString axisStr = axis_name(axis);
        TString xTitle = "q_{" + axisStr + "} [GeV/c]";

        xAxis->SetTitle(xTitle);
        xAxis->CenterTitle();
        xAxis->SetTitleSize(0.05f);
        xAxis->SetLabelSize(0.045f);
        xAxis->SetRangeUser(-0.4, 0.4);
    }

    {
        TAxis* yAxis = h->GetYaxis();
        yAxis->SetTitle(fmt);
        yAxis->CenterTitle();
        yAxis->SetTitleSize(0.05f);
        yAxis->SetLabelSize(0.045f);
        yAxis->SetTitleOffset(1.2f);
        yAxis->SetRangeUser(0.5, 1.5);
    }
}

TH1D* project_ratio(TH3D& ratio, int centr, int b, LCMSAxis axis, double sliceWidth,
                    const std::string& binName)
{
    TString name = TString::Format("proj_of_ratios_%d_%d_%s", centr, b, axis_name(axis).data());

    std::unique_ptr<TH1D> r(project_1d(ratio, axis, sliceWidth));
    r->SetName(name);

    TString axisStr = axis_name(axis);
    TString fmt = "C^{++}/C^{--}_{" + axisStr + "}";

    Style(r.get(), fmt, axis, centr, binName);
    return r.release();
}

TH1D* ratio_project(TH3D& neg, TH3D& pos, int centr, int b, LCMSAxis axis, double sliceWidth,
                    const std::string& binName)
{
    TString name = TString::Format("ratio_proj_%d_%d_%s", centr, b, axis_name(axis).data());

    std::unique_ptr<TH1D> n(project_1d(neg, axis, sliceWidth));
    std::unique_ptr<TH1D> p(project_1d(pos, axis, sliceWidth));

    std::unique_ptr<TH1D> ratio(static_cast<TH1D*>(n->Clone(name)));
    ratio->Divide(p.get(), n.get());

    TString axisStr = axis_name(axis);
    TString fmt = "C^{++}(q_{" + axisStr + "}) / C^{--}(q_{" + axisStr + "})";

    Style(ratio.get(), fmt, axis, centr, binName);
    return ratio.release();
}

} // namespace

void do_cf_ratios(Config& cfg, TFile* fCF3D, TFile* fRatioProj, TFile* fProjRatio)
{
    const Bin& bin = cfg.binning;
    logging::info("ratios: " + std::to_string(cfg.selection.centralities.size()) + " centralities x " +
              std::to_string(bin.count) + " bins");

    for (const int centr : cfg.selection.centralities) {
        for (int b = 0; b < bin.count; b++) {
            TString name_pos = get_cf_name(0, centr, cfg.input.type, bin.names[b]);
            TString name_neg = get_cf_name(1, centr, cfg.input.type, bin.names[b]);

            TH3D* hist_neg = dynamic_cast<TH3D*>(fCF3D->Get(name_neg));
            TH3D* hist_pos = dynamic_cast<TH3D*>(fCF3D->Get(name_pos));

            if (!hist_neg || !hist_pos) {
                continue;
            }

            auto neg = std::unique_ptr<TH3D>(static_cast<TH3D*>(hist_neg->Clone()));
            neg->SetDirectory(nullptr);
            auto pos = std::unique_ptr<TH3D>(static_cast<TH3D*>(hist_pos->Clone()));
            pos->SetDirectory(nullptr);

            auto ratio = std::unique_ptr<TH3D>(static_cast<TH3D*>(neg->Clone()));
            ratio->SetDirectory(nullptr);
            ratio->Reset();
            ratio->Divide(pos.get(), neg.get());

            for (auto axis : {LCMSAxis::Out, LCMSAxis::Side, LCMSAxis::Long}) {
                {
                    std::unique_ptr<TH1D> h(ratio_project(
                        *neg, *pos, centr, b, axis, cfg.projections.slice_ratio, bin.names[b]));
                    fRatioProj->cd();
                    h->Write();
                }

                {
                    std::unique_ptr<TH1D> h(project_ratio(
                        *ratio, centr, b, axis, cfg.projections.slice_ratio, bin.names[b]));
                    fProjRatio->cd();
                    h->Write();
                }
            }
        }
    }
}
