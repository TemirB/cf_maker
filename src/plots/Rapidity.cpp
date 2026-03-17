#include "plots.h"

#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>

#include "draw.h"

void Rapidity::MakeDependency(TFile* f, FitGrid& fitRes) const {
    for (int chIdx = 0; chIdx < chargeSize; chIdx++) {
        TMultiGraph* mg_R[6] = { new TMultiGraph(), new TMultiGraph(), new TMultiGraph(), new TMultiGraph(), new TMultiGraph(), new TMultiGraph() };
        TMultiGraph* mg_L = new TMultiGraph();
        std::vector<std::pair<TObject*, std::string>> legendEntries;

        for (int centIdx = 0; centIdx < centralitySize; centIdx++) {
            TGraphErrors* g_R[6] = { new TGraphErrors(), new TGraphErrors(), new TGraphErrors(),new TGraphErrors(), new TGraphErrors(), new TGraphErrors() };
            TGraphErrors* g_L    = new TGraphErrors();
            for (int lcmsIdx = 0; lcmsIdx < lcmsSize; lcmsIdx++) {
                g_R[lcmsIdx]->SetName(
                    Form("g_R_%s_%s_centr_%s", LCMS[lcmsIdx], chargeNames[chIdx], centralityNames[centIdx])
                );
                g_R[lcmsIdx]->SetLineColor(colors[centIdx]);
                g_R[lcmsIdx]->SetMarkerColor(colors[centIdx]);
                g_R[lcmsIdx]->SetMarkerStyle(markers[centIdx]);
            }

            g_L->SetName(
                Form("g_L_%s_centr_%s", chargeNames[chIdx], centralityNames[centIdx])
            );
            g_L->SetLineColor(colors[centIdx]);
            g_L->SetMarkerColor(colors[centIdx]);
            g_L->SetMarkerStyle(markers[centIdx]);

            legendEntries.push_back({g_R[0], centralityNames[centIdx]});

            // auto mDot = badDots[centIdx];
            int shift = 0;
            for (int yIdx = 0; yIdx < rapiditySize; yIdx++) {
                FitResult res = fitRes[chIdx][centIdx][yIdx];
                if (IsBadFit(res)) {
                    shift++;
                    continue;
                }

                double left  = rapidityValues[0] + yIdx * step;
                double right = left + step;

                std::string name = Form("[%.2f;%.2f]", left, right);
                double xVal = rapidityValues[0] + (yIdx + 0.5) * step;

                for (int lcmsIdx = 0; lcmsIdx < lcmsSize; lcmsIdx++) {
                    g_R[lcmsIdx]->SetPoint(yIdx - shift, xVal, res.R[lcmsIdx]);
                    g_R[lcmsIdx]->SetPointError(yIdx - shift, 0, res.eR[lcmsIdx]);
                }

                g_L->SetPoint(yIdx - shift, xVal, res.lambda);
                g_L->SetPointError(yIdx - shift, 0, res.elambda);
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
            if (lcmsIdx >= 3) type = 1;
            writeMGWithLegend(f, mg_R[lcmsIdx],
                mg_R[lcmsIdx]->GetName(),
                "rapidity",
                Form("R_{%s} (fm)", LCMS[lcmsIdx]),
                legendEntries,
                type
            );
        }

        setRangeWithErrors(mg_L, 0.1);
        mg_L->SetName(Form("mg_L_%s", chargeNames[chIdx]));
        writeMGWithLegend(f, mg_L,
            mg_L->GetName(),
            "rapidity",
            "lambda",
            legendEntries,
            2
        );

        // Saved for article
        {
            std::string name = Form("c_all_graphs_%s", chargeNames[chIdx]);
            std::string title = Form("Graphs for radii and #lambda(%s)", chargeNames[chIdx]);
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
            std::string name = Form("c_all_cross_graphs_%s", chargeNames[chIdx]);
            std::string title = Form("Graphs for cross radii(%s)", chargeNames[chIdx]);
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