#include "plots.h"

#include <memory>

#include <TH3.h>
#include <TF3.h>
#include <TCanvas.h>
#include "TLegend.h"
#include "TDirectory.h"
#include <TPaveText.h>
#include <TStyle.h>
#include <TSystem.h>

#include "helpers.h"
#include "fit/fit.h"
#include "draw.h"

void SetSlice(TH3D& h, LCMSAxis axis, double w) {
    ResetRanges(h);

    if (axis != LCMSAxis::Out)  h.GetXaxis()->SetRangeUser(-w, w);
    if (axis != LCMSAxis::Side) h.GetYaxis()->SetRangeUser(-w, w);
    if (axis != LCMSAxis::Long) h.GetZaxis()->SetRangeUser(-w, w);
}


TH3D* MakeFitHistogram(const TF3& fit) {
    auto h = new TH3D(
        "fit3d","fit3d",
        80,-0.4,0.4,
        80,-0.4,0.4,
        80,-0.4,0.4
    );
    h->SetDirectory(nullptr);

    for (int i=1;i<=80;i++)
    for (int j=1;j<=80;j++)
    for (int k=1;k<=80;k++) {
        double x = h->GetXaxis()->GetBinCenter(i);
        double y = h->GetYaxis()->GetBinCenter(j);
        double z = h->GetZaxis()->GetBinCenter(k);
        h->SetBinContent(i,j,k, fit.Eval(x,y,z));
    }
    return h;
}

TH1D* CreateFit(
    TH3D* fit3d,
    const LCMSAxis axis, 
    const std::string& baseName,
    int xF = 1, int xL = -1,
    int yF = 1, int yL = -1,
    int zF = 1, int zL = -1
) {
    TH1D* fit = nullptr;
    switch (axis) {
        case LCMSAxis::Out:
            fit = fit3d->ProjectionX((baseName + "_fit").c_str(), yF, yL, zF, zL);
            break;
        case LCMSAxis::Side:
            fit = fit3d->ProjectionY((baseName + "_fit").c_str(), xF, xL, zF, zL);
            break;
        case LCMSAxis::Long:
            fit = fit3d->ProjectionZ((baseName + "_fit").c_str(), xF, xL, yF, yL);
            break;
    }

    fit->SetDirectory(nullptr);
    return fit;
}

TH1D* BuildLCMSFitFrom3D(
    const TH3D& slicedVolume,
    TH3D* fit3D,
    LCMSAxis axis,
    const std::string& tag
) {
    const char* proj =
        axis == LCMSAxis::Out  ? "x" :
        axis == LCMSAxis::Side ? "y" :
                                 "z";

    auto den = (TH3D*)slicedVolume.Clone();
    auto num = (TH3D*)slicedVolume.Clone();
    den->Reset(); num->Reset();
    den->SetDirectory(nullptr);
    num->SetDirectory(nullptr);

    int x1 = slicedVolume.GetXaxis()->GetFirst();
    int x2 = slicedVolume.GetXaxis()->GetLast();
    int y1 = slicedVolume.GetYaxis()->GetFirst();
    int y2 = slicedVolume.GetYaxis()->GetLast();
    int z1 = slicedVolume.GetZaxis()->GetFirst();
    int z2 = slicedVolume.GetZaxis()->GetLast();

    for (int x=x1;x<=x2;x++)
    for (int y=y1;y<=y2;y++)
    for (int z=z1;z<=z2;z++) {
        den->SetBinContent(x,y,z,1.0);
        num->SetBinContent(x,y,z, fit3D->GetBinContent(x,y,z));
    }

    auto hNum = (TH1D*)num->Project3D(proj)->Clone();
    auto hDen = (TH1D*)den->Project3D(proj)->Clone();

    auto fit = (TH1D*)hDen->Clone(("fit_"+tag+"_"+ToString(axis)).c_str());
    fit->Divide(hNum, hDen);

    return fit;
}

// ============================================
// One LCMS projection → canvas in output ROOT
// ============================================

void AddFitStats(const FitResult& r, double textSize,
                 double x1=0.55, double y1=0.70, double x2=0.88, double y2=0.88) {
                    //0.60,0.75,0.88,0.88
    auto* stats = new TPaveText(x1, y1, x2, y2, "NDC");
    stats->SetBorderSize(1);
    stats->SetFillColor(kWhite);
    stats->SetTextAlign(12);
    stats->SetTextFont(42);
    stats->SetTextSize(textSize);
    
    stats->AddText(Form("#chi^{2}/ndf = %.1f / %d = %.3f", r.chi2, r.ndf, r.chi2/r.ndf));
    stats->AddText(Form("R_{out}  = %.3f #pm %.3f fm", r.R[0], r.eR[0]));
    stats->AddText(Form("R_{side} = %.3f #pm %.3f fm", r.R[1], r.eR[1]));
    stats->AddText(Form("R_{long} = %.3f #pm %.3f fm", r.R[2], r.eR[2]));
    stats->AddText(Form("R_{out-long} = %.3f #pm %.3f fm", r.R[4], r.eR[4]));
    stats->AddText(Form("#lambda  = %.3f #pm %.3f", r.lambda, r.elambda));
    stats->Draw();
}

std::pair<TH1D*, TH1D*> Create1D(
    TH3D& A,
    TH3D& Awei,
    TH3D* fit3d,
    const LCMSAxis axis,
    const std::string& name
) {

    // auto workA = (TH3D*) A.Clone();
    // auto workAwei = (TH3D*) Awei.Clone();
    std::string baseName = name + " " + ToString(axis);

    SetSlice(A, axis, 0.05);
    SetSlice(Awei, axis, 0.05);

    const char* proj =
        axis == LCMSAxis::Out  ? "x" :
        axis == LCMSAxis::Side ? "y" :
                                 "z";


    // const auto [xF, xL] = std::pair{1, A.GetNbinsX()};
    // const auto [yF, yL] = std::pair{1, A.GetNbinsY()};
    // const auto [zF, zL] = std::pair{1, A.GetNbinsZ()};

    // TH1D* hDen = nullptr;
    // switch (axis) {
    //     case LCMSAxis::Out:  hDen = A.ProjectionX((baseName + "_num").c_str(), yF, yL, zF, zL); break;
    //     case LCMSAxis::Side: hDen = A.ProjectionY((baseName + "_num").c_str(), xF, xL, zF, zL); break;
    //     case LCMSAxis::Long: hDen = A.ProjectionZ((baseName + "_num").c_str(), xF, xL, yF, yL); break;
    // }
    // hDen->SetDirectory(nullptr);

    // TH1D* hNum = nullptr;
    // switch (axis) {
    //     case LCMSAxis::Out:  hNum = Awei.ProjectionX((baseName + "_den").c_str(), yF, yL, zF, zL); break;
    //     case LCMSAxis::Side: hNum = Awei.ProjectionY((baseName + "_den").c_str(), xF, xL, zF, zL); break;
    //     case LCMSAxis::Long: hNum = Awei.ProjectionZ((baseName + "_den").c_str(), xF, xL, yF, yL); break;
    // }
    // hNum->SetDirectory(nullptr);

    auto hA    = (TH1D*)A.Project3D(proj)->Clone();
    auto hAwei = (TH1D*)Awei.Project3D(proj)->Clone();

    auto CF = (TH1D*)hA->Clone(name.c_str());
    CF->Divide(hAwei, hA, 1., 1., "B");

    std::string n = name + " " + ToString(axis);
    CF->SetTitle(n.c_str());

    auto fit = BuildLCMSFitFrom3D(A, fit3d, axis, baseName);
    // TH1D* fit = CreateFit(
    //     fit3d, axis, n, 
    //     xF, xL,
    //     yF, yL,
    //     zF, zL
    // );

    return {CF, fit};
}

TCanvas* CreateCanvas(
    TH1* CF, TH1* fit, 
    const char* name, 
    const LCMSAxis axis, 
    const FitResult& r
) {
    TCanvas* c = new TCanvas(name, name, 800, 600);

    c->SetTicks(1,1);
    c->SetLeftMargin(0.12);
    c->SetBottomMargin(0.12);

    Style1DCF(CF, name, ToString(axis).data());
    StyleFit(fit);

    CF->GetYaxis()->SetRangeUser(0.9,1.7);
    CF->GetXaxis()->SetRangeUser(-0.2,0.2);
    CF->SetStats(0);

    fit->GetXaxis()->SetRangeUser(-0.2,0.2);
    gStyle->SetOptFit(0102);

    c->cd();
    CF->Draw("P");
    fit->Draw("L SAME");
    AddFitStats(r, 0.032);

    c->Modified();
    c->Update();

    return c;
}

// ============================================
// Full LCMS pipeline
// ============================================

void MakeLCMS1DProjections(Context ctx, TFile* input, TFile* out)
{
    FitGrid fitRes = ctx.fitRes;
    Bin bin = ctx.bining;

    std::string dir = ctx.outDir + "/all_1d_histos";
    EnsureDir(dir);
    std::string format = "png";

    gSystem->mkdir(dir.data());

    for (int chIdx = 0; chIdx < Charge::kCount; chIdx++)
    for (int centIdx = 0; centIdx < Centrality::kCount; centIdx++)
    for (int b = 0; b < bin.count; b++) {
        TF3* fit = CreateCF3DFit(ctx, chIdx, centIdx, b);

        TH3D* h_A = getNum(input, chIdx, centIdx, b);
        TH3D* h_A_wei = getNumWei(input, chIdx, centIdx, b);

        const FitResult r = fitRes[chIdx][centIdx][b];
        fit->SetParameter(0, r.R[0]);
        fit->SetParameter(1, r.R[1]);
        fit->SetParameter(2, r.R[2]);
        fit->SetParameter(3, r.R[3]);
        fit->SetParameter(4, r.R[4]);
        fit->SetParameter(5, r.R[5]);
        fit->SetParameter(6, r.lambda);

        for (int p = 0; p < 7; p++) fit->FixParameter(p, fit->GetParameter(p));

        TH3D* fit3d = MakeFitHistogram(*fit);

        std::string cf_name = getCFName(chIdx, centIdx, bin.type, bin.names[b]);

        TCanvas* can = new TCanvas(cf_name.data(), cf_name.data(), 2400, 800);
        can->SetTitle(cf_name.data());
        can->Divide(3,1);

        for (int lcms = 0; lcms < 3; lcms++) {
            LCMSAxis axis = LCMSAxis(lcms);

            auto [CF, fit] = Create1D(*h_A, *h_A_wei, fit3d, axis, cf_name);

            std::string canvasName = cf_name + "_" + ToString(axis);

            TCanvas* c = CreateCanvas(CF, fit, canvasName.c_str(), axis, r);

            out->cd();
            c->Write();

            {
                can->cd(lcms + 1);

                CF->Draw("P");
                fit->Draw("L SAME");
                
                if (lcms == 2) {
                    AddFitStats(r, 0.032);
                }

                can->Modified();
                can->Update();
            }
        }

        // Rewrite condition bw overview
        // if (b >=4 && b <=6 && centIdx <= 2) 
        {
            auto name = Form(
                    "%s/cfs_%s_%s_%s.%s", 
                    dir.data(), 
                    Charge::kNames[chIdx], Centrality::kNames[centIdx], bin.names[b], 
                    format.data()
                );
            can->SaveAs(name);
        }
    }
}
