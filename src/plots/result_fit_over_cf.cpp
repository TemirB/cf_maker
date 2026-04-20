#include "plots.h"

#include <algorithm>
#include <cmath>

#include <TH3D.h>

#include "helpers.h"
#include "fit/fit.h"
#include "fit/types.h"

namespace {
constexpr double kFitRange = 0.3;

struct MeanWithError {
    double mean{};
    double err{};
    bool ok{false};
};

MeanWithError ComputeFitOverCFMean(const TH3D& cf, const FitResult& r) {
    double sum = 0.0;
    double sumSq = 0.0;
    int n = 0;

    for (int ix = 1; ix <= cf.GetNbinsX(); ++ix)
    for (int iy = 1; iy <= cf.GetNbinsY(); ++iy)
    for (int iz = 1; iz <= cf.GetNbinsZ(); ++iz) {
        const double qOut = cf.GetXaxis()->GetBinCenter(ix);
        const double qSide = cf.GetYaxis()->GetBinCenter(iy);
        const double qLong = cf.GetZaxis()->GetBinCenter(iz);

        if (std::abs(qOut) > kFitRange || std::abs(qSide) > kFitRange || std::abs(qLong) > kFitRange) {
            continue;
        }

        const double cfVal = cf.GetBinContent(ix, iy, iz);
        if (!std::isfinite(cfVal) || cfVal <= 0.0) {
            continue;
        }

        const double fitVal = EvalCF3D(r, qOut, qSide, qLong);
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

    if (n == 0) {
        return {};
    }

    const double mean = sum / n;
    const double variance = std::max(0.0, sumSq / n - mean * mean);
    const double err = std::sqrt(variance / n);
    return {mean, err, true};
}
} // namespace

TGraphErrors* BuildFitOverCFGraph(const Context& ctx, TFile* cf3dFile, int ch, int centr) {
    Bin bin = ctx.bining;

    TGraphErrors* g = new TGraphErrors();
    g->SetName(Form("g_FitOverCF_%s_centr_%s", Charge::kNames[ch], Centrality::kNames[centr]));
    g->SetLineColor(colors[centr]);
    g->SetMarkerColor(colors[centr]);
    g->SetMarkerStyle(markers[centr]);

    if (!cf3dFile) {
        return g;
    }

    int point = 0;
    for (int b = 0; b < bin.count; ++b) {
        const FitResult& res = ctx.fitRes[ch][centr][b];
        const std::string cfName = getCFName(ch, centr, bin.type, bin.names[b]);
        TH3D* hCF = (TH3D*)cf3dFile->Get(cfName.c_str());
        if (!hCF) {
            continue;
        }

        const MeanWithError ratioStats = ComputeFitOverCFMean(*hCF, res);
        if (!ratioStats.ok) {
            continue;
        }

        const double xVal = (bin.values[b] + bin.values[b + 1]) / 2.0;
        g->SetPoint(point, xVal, ratioStats.mean);
        g->SetPointError(point, 0, ratioStats.err);
        ++point;
    }

    return g;
}
