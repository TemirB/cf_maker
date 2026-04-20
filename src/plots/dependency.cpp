#include "plots.h"

#include <iostream>
#include <string>

#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TTree.h>
#include <TCanvas.h>

#include "helpers.h"
#include "fit/types.h"
#include "fit/fit.h"
#include "draw.h"
#include "context.h"

void MakeDependency(
    Context ctx,
    TFile* cf3dFile,
    TFile* outFile
) {
    FitGrid fitRes = ctx.fitRes;
    Bin bin = ctx.bining;

    std::string dir = ctx.outDir + "/dependency";
    EnsureDir(dir);
    std::string ext = "png";
    for (int ch = 0; ch < Charge::kCount; ch++) {
        // Creating: 
        // 1) 6 (LCMS::kCount) multigraphs, for each projection of radius
        // 2) 1, for lambda
        // 3) 1, legend buffer
        TMultiGraph* mg_R[LCMS::kCount] = { new TMultiGraph(), new TMultiGraph(), new TMultiGraph(), new TMultiGraph(), new TMultiGraph(), new TMultiGraph() };
        TMultiGraph* mg_L = new TMultiGraph();
        TMultiGraph* mg_Chi2Ndf = new TMultiGraph();
        TMultiGraph* mg_FitOverCF = new TMultiGraph();
        std::vector<std::pair<TObject*, std::string>> legendEntries;

        for (int centr = 0; centr < Centrality::kCount; centr++) {
            // Creating:
            // 1) 6 graphs, for each multigraph
            // 2) 1 graph, for lambda multigraph
            TGraphErrors* g_R[LCMS::kCount] = { new TGraphErrors(), new TGraphErrors(), new TGraphErrors(), new TGraphErrors(), new TGraphErrors(), new TGraphErrors() };
            TGraphErrors* g_L = new TGraphErrors();
            TGraphErrors* g_Chi2Ndf = BuildChi2NdfGraph(ctx, ch, centr);
            TGraphErrors* g_FitOverCF = BuildFitOverCFGraph(ctx, cf3dFile, ch, centr);

            for (int lcms = 0; lcms < LCMS::kCount; lcms++) {
                g_R[lcms]->SetName(
                    Form(
                        "g_R_%s_%s_centr_%s", 
                        LCMS::kNames[lcms], Charge::kNames[ch], Centrality::kNames[centr]
                    )
                );
                g_R[lcms]->SetLineColor(colors[centr]);
                g_R[lcms]->SetMarkerColor(colors[centr]);
                g_R[lcms]->SetMarkerStyle(markers[centr]);
            }

            g_L->SetName(Form("g_L_%s_centr_%s", Charge::kNames[ch], Centrality::kNames[centr]));
            g_L->SetLineColor(colors[centr]);
            g_L->SetMarkerColor(colors[centr]);
            g_L->SetMarkerStyle(markers[centr]);

            legendEntries.push_back({g_R[0], Centrality::kNames[centr]});

            int shift = 0;
            for (int b = 0; b < bin.count; b++) {
                FitResult res = fitRes[ch][centr][b];
                // if (IsBadFit(res)) {
                //     shift++;
                //     continue;
                // }

                double xVal = (bin.values[b] + bin.values[b+1])/2.;

                for (int lcms = 0; lcms < LCMS::kCount; lcms++) {
                    g_R[lcms]->SetPoint(b - shift, xVal, res.R[lcms]);
                    g_R[lcms]->SetPointError(b - shift, 0, res.eR[lcms]);
                }

                g_L->SetPoint(b - shift, xVal, res.lambda);
                g_L->SetPointError(b - shift, 0, res.elambda);
            }

            for (int lcms = 0; lcms < LCMS::kCount; lcms++) {
                mg_R[lcms]->Add(g_R[lcms], "lp");
            }

            mg_L->Add(g_L, "lp");
            mg_Chi2Ndf->Add(g_Chi2Ndf, "lp");
            mg_FitOverCF->Add(g_FitOverCF, "lp");
        }

        for (int lcms = 0; lcms < LCMS::kCount; lcms++) {
            mg_R[lcms]->SetName(Form("mg_R_%s_%s", LCMS::kNames[lcms], Charge::kNames[ch]));
            setRangeWithErrors(mg_R[lcms], 0.1);
            int type = 0;
            if (lcms >= 3) type = 1;
            writeMGWithLegend(outFile, mg_R[lcms],
                mg_R[lcms]->GetName(),
                ctx.bining.type,
                Form("R_{%s} (fm)", LCMS::kNames[lcms]),
                legendEntries,
                type
            );
        }

        setRangeWithErrors(mg_L, 0.1);
        mg_L->SetName(Form("mg_L_%s", Charge::kNames[ch]));
        writeMGWithLegend(outFile, mg_L,
            mg_L->GetName(),
            ctx.bining.type,
            "lambda",
            legendEntries,
            2
        );

        setRangeWithErrors(mg_Chi2Ndf, 0.1);
        mg_Chi2Ndf->SetName(Form("mg_Chi2Ndf_%s", Charge::kNames[ch]));
        writeMGWithLegend(outFile, mg_Chi2Ndf,
            mg_Chi2Ndf->GetName(),
            ctx.bining.type,
            "#chi^{2}/ndf",
            legendEntries,
            3
        );

        setRangeWithErrors(mg_FitOverCF, 0.1);
        mg_FitOverCF->SetName(Form("mg_FitOverCF_%s", Charge::kNames[ch]));
        writeMGWithLegend(outFile, mg_FitOverCF,
            mg_FitOverCF->GetName(),
            ctx.bining.type,
            "<fit/CF>",
            legendEntries,
            4
        );

        // Saved for article
        {
            std::string name = Form("c_all_graphs_%s", Charge::kNames[ch]);
            std::string title = Form("Graphs for radii and #lambda(%s)", Charge::kNames[ch]);
            TCanvas* c = new TCanvas(name.c_str(), title.c_str(), 3200, 1600);
            c->Divide(2, 2);
            for (int lcms = 0; lcms < 3; lcms++) {
                c->cd(lcms + 1);
                mg_R[lcms]->Draw("APL");
            }

            c->cd(4);
            mg_L->Draw("APL");
            std::string nSave = dir + "/" + name + "." + ext;
            SaveCanvasQuiet(c, nSave.c_str());
        }

        // Saved for article
        {
            std::string name = Form("c_fit_quality_%s", Charge::kNames[ch]);
            std::string title = Form("Fit quality(%s)", Charge::kNames[ch]);
            TCanvas* c = new TCanvas(name.c_str(), title.c_str(), 2400, 1000);
            c->Divide(2, 1);

            c->cd(1);
            gPad->SetLogy();
            mg_Chi2Ndf->Draw("APL");
            mg_Chi2Ndf->GetXaxis()->SetTitle(bin.type);
            mg_Chi2Ndf->GetYaxis()->SetTitle("#chi^{2}/ndf");

            c->cd(2);
            mg_FitOverCF->Draw("APL");
            mg_FitOverCF->GetXaxis()->SetTitle(bin.type);
            mg_FitOverCF->GetYaxis()->SetTitle("<fit/CF>");

            std::string nSave = dir + "/" + name + "." + ext;
            SaveCanvasQuiet(c, nSave.c_str());
        }

        // Saved for article
        {
            std::string name = Form("c_all_cross_graphs_%s", Charge::kNames[ch]);
            std::string title = Form("Graphs for cross radii(%s)", Charge::kNames[ch]);
            TCanvas* c = new TCanvas(name.c_str(), title.c_str(), 3200, 1600);
            c->Divide(2, 2);
            int idx = 1;
            for (int lcms = 0; lcms < 3; lcms++) {
                c->cd(idx);
                mg_R[lcms+3]->Draw("APL");
                idx++;
            }

            std::string nSave = dir + "/" + name + "." + ext;
            SaveCanvasQuiet(c, nSave.c_str());
        }
    }
}