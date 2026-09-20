#include "analysis/graphs.h"

#include <algorithm>
#include <cmath>

#include <TFile.h>
#include <TGraphErrors.h>
#include <TH3D.h>
#include <TString.h>

#include "core/binning.h"
#include "draw/draw.h"
#include "fit/model.h"
#include "io/input.h"

TGraphErrors* build_chi2_ndf_graph(const Config& cfg, int ch, int centr)
{
    const Bin& bin = cfg.binning;

    TGraphErrors* g = make_styled_graph(
        Form("g_chi2_ndf_%s_centr_%s", charge::kNames[ch], centrality::kNames[centr]), centr);

    int point = 0;
    for (int b = 0; b < bin.count; ++b) {
        const FitResult& res = cfg.fit_results[ch][centr][b];
        if (res.ndf <= 0) {
            continue;
        }

        g->SetPoint(point, bin_center(bin, b), res.chi2_ndf());
        g->SetPointError(point, 0, 0);
        ++point;
    }

    return g;
}

TGraphErrors* build_pvalue_graph(const Config& cfg, int ch, int centr)
{
    const Bin& bin = cfg.binning;

    TGraphErrors* g = make_styled_graph(
        Form("g_pvalue_%s_centr_%s", charge::kNames[ch], centrality::kNames[centr]), centr);

    for (int b = 0; b < bin.count; b++) {
        const FitResult& res = cfg.fit_results[ch][centr][b];

        g->SetPoint(b, bin_center(bin, b), res.p_value);
        g->SetPointError(b, 0, 0);
    }

    return g;
}

namespace
{
struct MeanWithError
{
    double mean{};
    double err{};
    bool ok{false};
};

MeanWithError compute_fit_over_cf_mean(const TH3D& cf, const FitResult& r, double fitRange)
{
    double sum = 0.0;
    double sumSq = 0.0;
    int n = 0;

    for (int ix = 1; ix <= cf.GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= cf.GetNbinsY(); ++iy) {
            for (int iz = 1; iz <= cf.GetNbinsZ(); ++iz) {
                const double qOut = cf.GetXaxis()->Getbin_center(ix);
                const double qSide = cf.GetYaxis()->Getbin_center(iy);
                const double qLong = cf.GetZaxis()->Getbin_center(iz);

                if (std::abs(qOut) > fitRange || std::abs(qSide) > fitRange ||
                    std::abs(qLong) > fitRange) {
                    continue;
                }

                const double cfVal = cf.GetBinContent(ix, iy, iz);
                if (!std::isfinite(cfVal) || cfVal <= 0.0) {
                    continue;
                }

                const double fitVal = eval_cf_3d(r, qOut, qSide, qLong);
                if (!std::isfinite(fitVal) || fitVal <= 0.0) {
                    continue;
                }

                const double ratio = fitVal / cfVal;
                if (!std::isfinite(ratio)) {
                    continue;
                }

                sum += ratio;
                sumSq += ratio * ratio;
                ++n;
            }
        }
    }

    if (n == 0) {
        return {};
    }

    const double mean = sum / n;
    const double variance = std::max(0.0, sumSq / n - mean * mean);
    const double err = std::sqrt(variance / n);
    return {mean, err, true};
}
} // namespace

TGraphErrors* build_fit_over_cf_graph(const Config& cfg, TFile* cf3dFile, int ch, int centr)
{
    const Bin& bin = cfg.binning;

    TGraphErrors* g = make_styled_graph(
        Form("g_FitOverCF_%s_centr_%s", charge::kNames[ch], centrality::kNames[centr]), centr);

    if (!cf3dFile) {
        return g;
    }

    int point = 0;
    for (int b = 0; b < bin.count; ++b) {
        const FitResult& res = cfg.fit_results[ch][centr][b];
        const std::string cf_name = get_cf_name(ch, centr, cfg.input.type, bin.names[b]);
        TH3D* cf_hist = dynamic_cast<TH3D*>(cf3dFile->Get(cf_name.c_str()));
        if (!cf_hist) {
            continue;
        }

        const MeanWithError ratio_stats =
            compute_fit_over_cf_mean(*cf_hist, res, cfg.projections.fit_over_cf_range);
        if (!ratio_stats.ok) {
            continue;
        }

        g->SetPoint(point, bin_center(bin, b), ratio_stats.mean);
        g->SetPointError(point, 0, ratio_stats.err);
        ++point;
    }

    return g;
}
