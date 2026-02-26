#include "plots.h"

#include <iostream>
#include <string>

#include <TMultiGraph.h>
#include <TGraphErrors.h>

#include "fit/types.h"
#include "fit/badfit.h"
#include "draw.h"

#include "TCanvas.h"

// bool deletePoint(int cent, int kt, LCMSAxis lcms) {

// }

// bool deletePoint(int lcms, int ch, int centr, int kt) {
//     if (kt == 3 && centr == 0) return true;

//     switch (lcms) {
//     case 0: 
//         if (centr == 1 && (kt == 2) || (kt == 3)) return true;

//         break;
//     case 1:
//         if 
//     }
// }
void MakeKtDependence(
    TFile* outFile,
    FitGrid& fitRes
) {
    for (int chIdx = 0; chIdx < chargeSize; chIdx++) {
        // 3 мультиграфа тут
        TMultiGraph* mg_R[6] = { new TMultiGraph(), new TMultiGraph(), new TMultiGraph(), new TMultiGraph(), new TMultiGraph(), new TMultiGraph() };
        TMultiGraph* mg_L = new TMultiGraph();

        std::vector<std::pair<TObject*, std::string>> legendEntries;
        
        for (int centIdx = 0; centIdx < centralitySize; centIdx++) {
            TGraphErrors* g_R[6] = { new TGraphErrors(), new TGraphErrors(), new TGraphErrors(), new TGraphErrors(), new TGraphErrors(), new TGraphErrors() };
            TGraphErrors* g_L    = new TGraphErrors();
            for (int lcmsIdx = 0; lcmsIdx < lcmsSize; lcmsIdx++) {
                g_R[lcmsIdx]->SetName(
                    Form("g_R_%s_%s_centr_%s", LCMS[lcmsIdx], chargeNames[chIdx], centralityNames[centIdx]));
                g_R[lcmsIdx]->SetLineColor(colors[centIdx]);
                g_R[lcmsIdx]->SetMarkerColor(colors[centIdx]);
                g_R[lcmsIdx]->SetMarkerStyle(markers[centIdx]);
            }
            g_L->SetName(Form("g_L_%s_centr_%s", chargeNames[chIdx], centralityNames[centIdx]));
            g_L->SetLineColor(colors[centIdx]);
            g_L->SetMarkerColor(colors[centIdx]);
            g_L->SetMarkerStyle(markers[centIdx]);

            legendEntries.push_back({g_R[0], centralityNames[centIdx]});

            int shift = 0;
            for (int ktIdx = 0; ktIdx < ktSize; ktIdx++) {
                FitResult res = fitRes[chIdx][centIdx][ktIdx];
                if (IsBadFit(res)) {
                    shift++;
                    continue;
                }

                double xval = (ktValues[ktIdx+1] + ktValues[ktIdx]) / 2.0;
                double xerr = ktValues[ktIdx+1] - xval;

                // std::cout << Form("charge: %d, centrality: %d, kt: %d, chi2: %f, ndf %d", chIdx, centIdx, ktIdx, res.chi2, res.ndf) << std::endl;
                for (int lcmsIdx = 0; lcmsIdx < lcmsSize; lcmsIdx++) {

                    // if (deletePoint(lcmsIdx, chIdx, centIdx, ktIdx)) {
                    //     // удалить, но сместить точку
                    // }

                    g_R[lcmsIdx]->SetPoint(ktIdx - shift, xval, res.R[lcmsIdx]);
                    g_R[lcmsIdx]->SetPointError(ktIdx - shift, xerr, res.eR[lcmsIdx]);

                    // std::cout << Form("R(%s)=%f+-%f\t", LCMS[lcmsIdx], res.R[lcmsIdx], res.eR[lcmsIdx]);
                }
                g_L->SetPoint(ktIdx - shift, xval, res.lambda);
                g_L->SetPointError(ktIdx - shift, xerr, res.elambda);

                // std::cout << Form("lambda=%f\t", res.lambda);

                // std::cout << std::endl;
            }

            for (int lcmsIdx = 0; lcmsIdx < lcmsSize; lcmsIdx++) {
                mg_R[lcmsIdx]->Add(g_R[lcmsIdx], "lp");
            }
            mg_L->Add(g_L, "lp");
        }

        for (int lcmsIdx = 0; lcmsIdx < lcmsSize; lcmsIdx++) {
            mg_R[lcmsIdx]->SetName(Form("mg_R_%s_%s", LCMS[lcmsIdx], chargeNames[chIdx]));
            setRangeWithErrors(mg_R[lcmsIdx], 0.1);
            int type = 0;
            if (lcmsIdx >= 3) {
                type = 1;
            }
            writeMGWithLegend(outFile, mg_R[lcmsIdx],
                mg_R[lcmsIdx]->GetName(),
                "k_{T} (GeV/c)",
                Form("R_{%s} (fm)", LCMS[lcmsIdx]),
                legendEntries,
                type
            );
        }
        setRangeWithErrors(mg_L, 0.1);
        mg_L->SetName(Form("mg_L_%s", chargeNames[chIdx]));
        writeMGWithLegend(outFile, mg_L,
            mg_L->GetName(),
            "k_{T} (GeV/c)",
            "lambda",
            legendEntries,
            2
        );

        // Saved for article
        {
            std::string name = Form("c_all_kt_graphs_%s", chargeNames[chIdx]);
            std::string title = Form("Graphs for radii and #lambda(%s)", chargeNames[chIdx]);
            TCanvas* c = new TCanvas(name.c_str(), title.c_str(), 1000, 800);
            c->Divide(2, 2);
            for (int lcms = 0; lcms < 3; lcms++) {
                c->cd(lcms + 1);
                mg_R[lcms]->Draw("APL");
            }

            c->cd(4);
            mg_L->Draw("APL");

            std::string nSave = name + ".pdf";
            c->SaveAs(nSave.c_str());
        }

        // Saved for article
        {
            std::string name = Form("c_all_kt_cross_graphs_%s", chargeNames[chIdx]);
            std::string title = Form("Graphs for cross radii(%s)", chargeNames[chIdx]);
            TCanvas* c = new TCanvas(name.c_str(), title.c_str(), 1000, 800);
            c->Divide(2, 2);
            int idx = 1;
            for (int lcms = 0; lcms < 3; lcms++) {
                c->cd(idx);
                mg_R[lcms+3]->Draw("APL");
                idx++;
            }

            std::string nSave = name + ".pdf";
            c->SaveAs(nSave.c_str());
        }
    }
};