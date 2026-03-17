#include "plots.h"

#include <string>

#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>

#include "draw.h"

void PrepareGs(
    TGraphErrors* R[6], TGraphErrors* L,
    int ch, int centr
) {
    for (int lcms = 0; lcms < lcmsSize; lcms++) {
        R[lcms]->SetName(
            Form(
                "g_R_%s_%s_centr_%s",
                LCMS[lcms], chargeNames[ch], centralityNames[centr]
            )
        );
        R[lcms]->SetLineColor(colors[centr]);
        R[lcms]->SetMarkerColor(colors[centr]);
        R[lcms]->SetMarkerStyle(markers[centr]);
    }

    L->SetName(
        Form("g_L_%s_centr_%s", chargeNames[ch], centralityNames[centr])
    );
    L->SetLineColor(colors[centr]);
    L->SetMarkerColor(colors[centr]);
    L->SetMarkerStyle(markers[centr]);
}

void MakeDependency(
    TFile* f, FitGrid& fitRes, int binSize, std::string aType
) {
    for (int ch = 0; ch < chargeSize; ch++) {
        TMultiGraph* mg_R[6];
        for (int i = 0; i < 6; i++) { mg_R[i] = new TMultiGraph(); }
        TMultiGraph* mg_L = new TMultiGraph();
        std::vector<std::pair<TObject*, std::string>> legendEntries;

        for (int centr = 0; centr < centralitySize; centr++) {
            TGraphErrors* g_R[6];
            for (int i = 0; i < 6; i++) { g_R[i] = new TGraphErrors(); }
            TGraphErrors* g_L = new TGraphErrors();

            PrepareGs(g_R, g_L, ch, centr);

            legendEntries.push_back({g_R[0], centralityNames[centr]});

            int shift = 0;
            for (int bin = 0; bin < binSize; bin++) {
                FitResult res = fitRes[ch][centr][bin];
                if (IsBadFit(res)) {
                    shift++;
                    continue;
                }
                // double left  = rapidityValues[0] + bin * step;
                // double right = left + step;

                // std::string name = Form("[%.2f;%.2f]", left, right);
                double xVal = rapidityValues[0] + (bin + 0.5) * step;

                for (int lcms = 0; lcms < lcmsSize; lcms++) {
                    g_R[lcms]->SetPoint(bin - shift, xVal, res.R[lcms]);
                    g_R[lcms]->SetPointError(bin - shift, 0, res.eR[lcms]);
                }

                g_L->SetPoint(bin - shift, xVal, res.lambda);
                g_L->SetPointError(bin - shift, 0, res.elambda);
            }
            for (int lcms = 0; lcms < lcmsSize; lcms++) {
                mg_R[lcms]->Add(g_R[lcms], "lp");
            }
            mg_L->Add(g_L, "lp");
        }

        for (int lcms = 0; lcms < lcmsSize; lcms++) {
            mg_R[lcms]->SetName(Form("mg_R_%s_%s", LCMS[lcms], chargeNames[ch]));
            setRangeWithErrors(mg_R[lcms], 0.1);
            int type = 0;
            if (lcms >= 3) type = 1;
            writeMGWithLegend(f, mg_R[lcms],
                mg_R[lcms]->GetName(),
                aType.data(),
                Form("R_{%s} (fm)", LCMS[lcms]),
                legendEntries,
                type
            );
        }

        setRangeWithErrors(mg_L, 0.1);
        mg_L->SetName(Form("mg_L_%s", chargeNames[ch]));
        writeMGWithLegend(f, mg_L,
            mg_L->GetName(),
            aType.data(),
            "lambda",
            legendEntries,
            2
        );

        // Saved for article
        {
            std::string name = Form("c_all_graphs_%s", chargeNames[ch]);
            std::string title = Form("Graphs for radii and #lambda(%s)", chargeNames[ch]);
            TCanvas* c = new TCanvas(name.c_str(), title.c_str(), 3200, 1600);
            c->Divide(2, 2);
            for (int lcms = 0; lcms < 3; lcms++) {
                c->cd(lcms + 1);
                mg_R[lcms]->Draw("APL");
            }

            c->cd(4);
            mg_L->Draw("APL");
            std::string nSave = name + ".png";
            c->SaveAs(nSave.c_str());
        }

        // Saved for article
        {
            std::string name = Form("c_all_cross_graphs_%s", chargeNames[ch]);
            std::string title = Form("Graphs for cross radii(%s)", chargeNames[ch]);
            TCanvas* c = new TCanvas(name.c_str(), title.c_str(), 3200, 1600);
            c->Divide(2, 2);
            int idx = 1;
            for (int lcms = 0; lcms < 3; lcms++) {
                c->cd(idx);
                mg_R[lcms+3]->Draw("APL");
                idx++;
            }

            std::string nSave = name + ".png";
            c->SaveAs(nSave.c_str());
        }
    }
}