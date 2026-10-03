#include "analysis/ratio.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TString.h>

#include "core/binning.h"
#include "core/correlation.h"
#include "core/lcms.h"
#include "core/log.h"
#include "io/input.h"

namespace
{

struct ChargeInput
{
    std::unique_ptr<TH3D> count;
    std::unique_ptr<TH3D> sum;
};

ChargeInput read_charge_input(TFile& input, int charge, int centrality, int bin,
                              CorrelationStatistics statistics)
{
    const auto [count, sum] = get_hists(&input, charge, centrality, bin);
    ChargeInput histograms{std::unique_ptr<TH3D>(count), std::unique_ptr<TH3D>(sum)};
    if (!histograms.count || !histograms.sum) {
        throw std::runtime_error("ratios: missing raw input for charge=" + std::to_string(charge) +
                                 ", centrality=" + std::to_string(centrality) +
                                 ", bin=" + std::to_string(bin));
    }
    validate_correlation_inputs(*histograms.sum, *histograms.count, statistics);
    return histograms;
}

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

TH1D* project_ratio(TH3D& ratio, TH3D& support, int centr, int b, LCMSAxis axis, double sliceWidth,
                    const std::string& binName)
{
    TString name = TString::Format("proj_of_ratios_%d_%d_%s", centr, b, axis_name(axis).data());

    std::unique_ptr<TH1D> r(project_1d(ratio, axis, sliceWidth));
    std::unique_ptr<TH1D> counts(project_1d(support, axis, sliceWidth));
    for (int bin = 0; bin < r->GetNcells(); ++bin) {
        const double count = counts->GetBinContent(bin);
        const double mean = count > 0 ? r->GetBinContent(bin) / count : 0;
        const double error = count > 0 ? r->GetBinError(bin) / count : 0;
        r->SetBinContent(bin, mean);
        r->SetBinError(bin, error);
    }
    r->SetName(name);

    TString axisStr = axis_name(axis);
    TString fmt = "#LT C^{++}/C^{--} #GT_{valid cells," + axisStr + "}";

    Style(r.get(), fmt, axis, centr, binName);
    return r.release();
}

TH1D* ratio_project(ChargeInput& neg, ChargeInput& pos, CorrelationStatistics statistics, int centr,
                    int b, LCMSAxis axis, double sliceWidth, const std::string& binName)
{
    TString name = TString::Format("ratio_proj_%d_%d_%s", centr, b, axis_name(axis).data());

    std::unique_ptr<TH1D> n_count(project_1d(*neg.count, axis, sliceWidth));
    std::unique_ptr<TH1D> n_sum(project_1d(*neg.sum, axis, sliceWidth));
    std::unique_ptr<TH1D> p_count(project_1d(*pos.count, axis, sliceWidth));
    std::unique_ptr<TH1D> p_sum(project_1d(*pos.sum, axis, sliceWidth));
    std::unique_ptr<TH1D> n(static_cast<TH1D*>(n_count->Clone()));
    std::unique_ptr<TH1D> p(static_cast<TH1D*>(p_count->Clone()));
    n->SetDirectory(nullptr);
    p->SetDirectory(nullptr);
    fill_correlation(*n, *n_sum, *n_count, statistics);
    fill_correlation(*p, *p_sum, *p_count, statistics);

    std::unique_ptr<TH1D> ratio(static_cast<TH1D*>(n->Clone(name)));
    ratio->SetDirectory(nullptr);
    if (!ratio->Divide(p.get(), n.get())) {
        throw std::runtime_error("Cannot divide projected charge correlation histograms");
    }
    for (int cell = 0; cell < ratio->GetNcells(); ++cell) {
        if (p_count->GetBinContent(cell) <= 0 || n_count->GetBinContent(cell) <= 0 ||
            n->GetBinContent(cell) == 0 || !std::isfinite(ratio->GetBinContent(cell)) ||
            !std::isfinite(ratio->GetBinError(cell))) {
            ratio->SetBinContent(cell, 0);
            ratio->SetBinError(cell, 0);
        }
    }

    TString axisStr = axis_name(axis);
    TString fmt =
        "#LT w^{++} #GT_{pairs,q_{" + axisStr + "}} / #LT w^{--} #GT_{pairs,q_{" + axisStr + "}}";

    Style(ratio.get(), fmt, axis, centr, binName);
    return ratio.release();
}

} // namespace

void do_cf_ratios(Config& cfg, TFile* fCF3D, TFile* fRatioProj, TFile* fProjRatio)
{
    for (const int charge : {0, 1}) {
        if (std::find(cfg.selection.charges.begin(), cfg.selection.charges.end(), charge) ==
            cfg.selection.charges.end()) {
            throw std::invalid_argument("ratios require both selected charges (0 and 1)");
        }
    }
    if (!fCF3D || !fRatioProj || !fProjRatio) {
        throw std::invalid_argument("ratios require CF and output files");
    }
    TDirectory::TContext context;
    std::unique_ptr<TFile> raw_input(TFile::Open(cfg.input.file.c_str(), "READ"));
    if (!raw_input || raw_input->IsZombie()) {
        throw std::runtime_error("ratios: cannot open raw input " + cfg.input.file);
    }
    const auto statistics = correlation_statistics(*raw_input);
    const Bin& bin = cfg.binning;
    logging::info("ratios: " + std::to_string(cfg.selection.centralities.size()) +
                  " centralities x " + std::to_string(bin.count) + " bins");

    for (const int centr : cfg.selection.centralities) {
        for (int b = 0; b < bin.count; b++) {
            TString name_pos = get_cf_name(0, centr, cfg.input.type, bin.names[b]);
            TString name_neg = get_cf_name(1, centr, cfg.input.type, bin.names[b]);

            TH3D* hist_neg = dynamic_cast<TH3D*>(fCF3D->Get(name_neg));
            TH3D* hist_pos = dynamic_cast<TH3D*>(fCF3D->Get(name_pos));

            if (!hist_neg || !hist_pos) {
                throw std::runtime_error("ratios: missing charge CF for centrality=" +
                                         std::to_string(centr) + ", bin=" + std::to_string(b));
            }

            auto raw_pos = read_charge_input(*raw_input, 0, centr, b, statistics);
            auto raw_neg = read_charge_input(*raw_input, 1, centr, b, statistics);
            // CF and original counts must describe the same cells for the support mask.
            validate_correlation_inputs(*hist_pos, *raw_pos.count,
                                        CorrelationStatistics::FixedReference);
            validate_correlation_inputs(*hist_neg, *raw_neg.count,
                                        CorrelationStatistics::FixedReference);

            auto neg = std::unique_ptr<TH3D>(static_cast<TH3D*>(hist_neg->Clone()));
            neg->SetDirectory(nullptr);
            auto pos = std::unique_ptr<TH3D>(static_cast<TH3D*>(hist_pos->Clone()));
            pos->SetDirectory(nullptr);

            auto ratio = std::unique_ptr<TH3D>(static_cast<TH3D*>(neg->Clone()));
            ratio->SetDirectory(nullptr);
            ratio->Reset();
            if (!ratio->Divide(pos.get(), neg.get())) {
                throw std::runtime_error("Cannot divide charge correlation histograms");
            }
            auto support = std::unique_ptr<TH3D>(static_cast<TH3D*>(neg->Clone()));
            support->SetDirectory(nullptr);
            support->Reset("ICES");
            for (int cell = 0; cell < ratio->GetNcells(); ++cell) {
                const double positive = pos->GetBinContent(cell);
                const double negative = neg->GetBinContent(cell);
                if (raw_pos.count->GetBinContent(cell) > 0 &&
                    raw_neg.count->GetBinContent(cell) > 0 && negative != 0 &&
                    std::isfinite(positive) && std::isfinite(negative) &&
                    std::isfinite(ratio->GetBinContent(cell)) &&
                    std::isfinite(ratio->GetBinError(cell))) {
                    support->SetBinContent(cell, 1);
                } else {
                    ratio->SetBinContent(cell, 0);
                    ratio->SetBinError(cell, 0);
                }
            }

            for (auto axis : {LCMSAxis::Out, LCMSAxis::Side, LCMSAxis::Long}) {
                {
                    std::unique_ptr<TH1D> h(ratio_project(raw_neg, raw_pos, statistics, centr, b,
                                                          axis, cfg.projections.slice_ratio,
                                                          bin.names[b]));
                    fRatioProj->cd();
                    h->Write();
                }

                {
                    std::unique_ptr<TH1D> h(project_ratio(*ratio, *support, centr, b, axis,
                                                          cfg.projections.slice_ratio,
                                                          bin.names[b]));
                    fProjRatio->cd();
                    h->Write();
                }
            }
        }
    }
}
