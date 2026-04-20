#include "plots.h"

#include <TGraphErrors.h>

#include "helpers.h"
#include "fit/types.h"

TGraphErrors* BuildChi2NdfGraph(const Context& ctx, int ch, int centr) {
    Bin bin = ctx.bining;

    TGraphErrors* g = new TGraphErrors();
    g->SetName(Form("g_Chi2Ndf_%s_centr_%s", Charge::kNames[ch], Centrality::kNames[centr]));
    g->SetLineColor(colors[centr]);
    g->SetMarkerColor(colors[centr]);
    g->SetMarkerStyle(markers[centr]);

    int point = 0;
    for (int b = 0; b < bin.count; ++b) {
        const FitResult& res = ctx.fitRes[ch][centr][b];
        if (res.ndf <= 0) {
            continue;
        }

        const double xVal = (bin.values[b] + bin.values[b + 1]) / 2.0;
        g->SetPoint(point, xVal, res.Chi2Ndf());
        g->SetPointError(point, 0, 0);
        ++point;
    }

    return g;
}
