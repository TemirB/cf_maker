#include "plots.h"

#include <string>

#include <TMultiGraph.h>
#include <TGraphErrors.h>

#include "helpers.h"
#include "fit/types.h"
#include "draw.h"

void MakeRapidityDependence(
    TFile* outFile,
    FitGrid& fitRes
) {
    double yStep = (rapidityValues[1] - rapidityValues[0]) / rapiditySize;

    for (int chIdx = 0; chIdx < chargeSize; chIdx++) {
        for (int centIdx = 0; centIdx < centralitySize; centIdx++) {
            TMultiGraph* mg_R[3] = { new TMultiGraph(), new TMultiGraph(), new TMultiGraph() };
            TMultiGraph* mg_L = new TMultiGraph();

            std::vector<std::pair<TObject*, std::string>> legendEntries;

            for (int ktIdx = 0; ktIdx < ktSize; ktIdx++) {
                TGraphErrors* g_R[3] = { new TGraphErrors(), new TGraphErrors(), new TGraphErrors() };
                TGraphErrors* g_L    = new TGraphErrors();
                for (int lcmsIdx = 0; lcmsIdx < lcmsSize; lcmsIdx++) {
                    g_R[lcmsIdx]->SetName(Form("g_R_%s_%s_kt_%s_centr_%s", LCMS[lcmsIdx].c_str(), chargeNames[chIdx].c_str(), ktNames[ktIdx].c_str(), centralityNames[centIdx].c_str()));
                    g_R[lcmsIdx]->SetLineColor(colors[ktIdx]);
                    g_R[lcmsIdx]->SetMarkerColor(colors[ktIdx]);
                    g_R[lcmsIdx]->SetMarkerStyle(markers[ktIdx]);
                }
                g_L->SetName(Form("g_L_%s_kt_%s_centr_%s", chargeNames[chIdx].c_str(), ktNames[ktIdx].c_str(), centralityNames[centIdx].c_str()));
                g_L->SetLineColor(colors[ktIdx]);
                g_L->SetMarkerColor(colors[ktIdx]);
                g_L->SetMarkerStyle(markers[ktIdx]);

                legendEntries.push_back({g_R[0], ktNames[ktIdx]});

                for (int yIdx = 0; yIdx < rapiditySize; yIdx++) {
                    FitResult res = fitRes[chIdx][centIdx][ktIdx][yIdx];
                    if (!res.ok) continue;

                    double left  = rapidityValues[0] + yIdx * yStep;
                    double right = left + yStep;

                    std::string name = Form("[%.2f;%.2f]", left, right);
                    double xVal = rapidityValues[0] + (yIdx + 0.5) * yStep;
                    


                    for (int lcmsIdx = 0; lcmsIdx < lcmsSize; lcmsIdx++) {
                        g_R[lcmsIdx]->SetPoint(yIdx, xVal, res.R[lcmsIdx]);
                        g_R[lcmsIdx]->SetPointError(yIdx, 0, res.eR[lcmsIdx]);
                    }

                    g_L->SetPoint(yIdx, xVal, res.lambda);
                    g_L->SetPointError(yIdx, 0, res.elambda);
                }

                for (int lcmsIdx = 0; lcmsIdx < lcmsSize; lcmsIdx++) {
                    mg_R[lcmsIdx]->Add(g_R[lcmsIdx], "lp");
                }
                mg_L->Add(g_L, "lp");
            }

            for (int lcmsIdx = 0; lcmsIdx < lcmsSize; lcmsIdx++) {
                mg_R[lcmsIdx]->SetName(Form("mg_R_%s_%s_centr_%s", LCMS[lcmsIdx].c_str(), chargeNames[chIdx].c_str(), centralityNames[centIdx].c_str()));
                setRangeWithErrors(mg_R[lcmsIdx], 0.1);
                writeMGWithLegend(outFile, mg_R[lcmsIdx],
                    mg_R[lcmsIdx]->GetName(),
                    "rapidity",
                    Form("R_{%s} (fm)", LCMS[lcmsIdx].c_str()),
                    legendEntries
                );
            }
            setRangeWithErrors(mg_L, 0.1);
            mg_L->SetName(Form("mg_L_%s_centr_%s", chargeNames[chIdx].c_str(), centralityNames[centIdx].c_str()));
            writeMGWithLegend(outFile, mg_L,
                mg_L->GetName(),
                "rapidity",
                "lambda",
                legendEntries
            );
        }
    }
};