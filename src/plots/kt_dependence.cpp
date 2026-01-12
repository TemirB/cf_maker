#include "plots.h"

#include <string>

#include <TMultiGraph.h>
#include <TGraphErrors.h>

#include "helpers.h"
#include "fit/types.h"
#include "draw.h"

void MakeKtDependence(
    TFile* outFile,
    FitGrid& fitRes
) {
    double yStep = (rapidityValues[1] - rapidityValues[0]) / rapiditySize;

    for (int chIdx = 0; chIdx < chargeSize; chIdx++) {
        for (int yIdx = 0; yIdx < rapiditySize; yIdx++) {
            double left  = rapidityValues[0] + yIdx * yStep;
            double right = left + yStep;
            if (left < -0.6 || right > 0.0) continue;
            std::string name = Form("[%.2f;%.2f]", left, right);
            
            // 3 мультиграфа тут
            TMultiGraph* mg_R[3] = { new TMultiGraph(), new TMultiGraph(), new TMultiGraph() };
            TMultiGraph* mg_L = new TMultiGraph();

            std::vector<std::pair<TObject*, std::string>> legendEntries;
            
            for (int centIdx = 0; centIdx < centralitySize; centIdx++) {
                TGraphErrors* g_R[3] = { new TGraphErrors(), new TGraphErrors(), new TGraphErrors() };
                TGraphErrors* g_L    = new TGraphErrors();
                for (int lcmsIdx = 0; lcmsIdx < lcmsSize; lcmsIdx++) {
                    g_R[lcmsIdx]->SetName(Form("g_R_%s_%s_y_%s_centr_%s", LCMS[lcmsIdx].c_str(), chargeNames[chIdx].c_str(), name.c_str(), centralityNames[centIdx].c_str()));
                    g_R[lcmsIdx]->SetLineColor(colors[centIdx]);
                    g_R[lcmsIdx]->SetMarkerColor(colors[centIdx]);
                    g_R[lcmsIdx]->SetMarkerStyle(markers[centIdx]);
                }
                g_L->SetName(Form("g_L_%s_y_%s_centr_%s", chargeNames[chIdx].c_str(), name.c_str(), centralityNames[centIdx].c_str()));
                g_L->SetLineColor(colors[centIdx]);
                g_L->SetMarkerColor(colors[centIdx]);
                g_L->SetMarkerStyle(markers[centIdx]);

                legendEntries.push_back({g_R[0], centralityNames[centIdx]});

                for (int ktIdx = 0; ktIdx < ktSize; ktIdx++) {
                    FitResult res = fitRes[chIdx][centIdx][ktIdx][yIdx];
                    if (!res.ok) continue;

                    double xval = (ktValues[ktIdx+1] + ktValues[ktIdx]) / 2.0;
                    double xerr = ktValues[ktIdx+1] - xval;

                    for (int lcmsIdx = 0; lcmsIdx < lcmsSize; lcmsIdx++) {
                        g_R[lcmsIdx]->SetPoint(ktIdx, xval, res.R[lcmsIdx]);
                        g_R[lcmsIdx]->SetPointError(ktIdx, xerr, res.eR[lcmsIdx]);
                    }

                    g_L->SetPoint(ktIdx, xval, res.lambda);
                    g_L->SetPointError(ktIdx, xerr, res.elambda);
                }

                for (int lcmsIdx = 0; lcmsIdx < lcmsSize; lcmsIdx++) {
                    mg_R[lcmsIdx]->Add(g_R[lcmsIdx], "lp");
                }
                mg_L->Add(g_L, "lp");
            }

            for (int lcmsIdx = 0; lcmsIdx < lcmsSize; lcmsIdx++) {
                mg_R[lcmsIdx]->SetName(Form("mg_R_%s_%s_y_%s", LCMS[lcmsIdx].c_str(), chargeNames[chIdx].c_str(), name.c_str()));
                setRangeWithErrors(mg_R[lcmsIdx], 0.1);
                writeMGWithLegend(outFile, mg_R[lcmsIdx],
                    mg_R[lcmsIdx]->GetName(),
                    "k_{T} (GeV/c)",
                    Form("R_{%s} (fm)", LCMS[lcmsIdx].c_str()),
                    legendEntries
                );
            }
            setRangeWithErrors(mg_L, 0.1);
            mg_L->SetName(Form("mg_L_%s_y_%s", chargeNames[chIdx].c_str(), name.c_str()));
            writeMGWithLegend(outFile, mg_L,
                mg_L->GetName(),
                "k_{T} (GeV/c)",
                "lambda",
                legendEntries
            );
        }
    }
};