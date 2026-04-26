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
TPaveText* GetFitStats(
    const FitResult& r, double textSize,
    double x1=0.55, double y1=0.70, double x2=0.88, double y2=0.88
) {
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

    return stats;
}

std::pair<TH1D*, TH1D*> Create1D(
    TH3D& A,
    TH3D& Awei,
    TH3D* fit3d,
    const LCMSAxis axis,
    const std::string& name
) {
    std::string baseName = name + " " + ToString(axis);

    SetSlice(A, axis, 0.05);
    SetSlice(Awei, axis, 0.05);

    const char* proj =
        axis == LCMSAxis::Out  ? "x" :
        axis == LCMSAxis::Side ? "y" :
                                 "z";

    auto hA    = (TH1D*)A.Project3D(proj)->Clone();
    auto hAwei = (TH1D*)Awei.Project3D(proj)->Clone();

    auto CF = (TH1D*)hA->Clone(name.c_str());
    CF->Divide(hAwei, hA, 1., 1., "B");

    std::string n = name + " " + ToString(axis);
    CF->SetTitle(n.c_str());

    auto fit = BuildLCMSFitFrom3D(A, fit3d, axis, baseName);

    return {CF, fit};
}

void SaveCanvasToFile(
    TFile* out, TH1D* CF, TH1D* fit, 
    TPaveText* stats, Draw::Style style,
    std::string cf_name, LCMSAxis axis
) {
    std::string name = cf_name + "_" + ToString(axis);
    TCanvas* c = new TCanvas(name.data(), name.data(), 800, 600);
    Style1DCF(CF, name, ToString(axis).data(), style);
    StyleFit(fit, style);

    CF->GetYaxis()->SetRangeUser(0.9,1.7);
    CF->GetXaxis()->SetRangeUser(-0.2,0.2);
    CF->SetStats(0);

    fit->GetXaxis()->SetRangeUser(-0.2,0.2);
    gStyle->SetOptFit(0102);

    c->cd();
    CF->Draw("P");
    fit->Draw("L SAME");
    stats->Draw();

    out->cd();
    c->Write();
}

void DrawCFAndFit(
    TCanvas* c, TH1D* CF, TH1D* fit, 
    TPaveText* stats, int lcms
) {
    c->cd(lcms+1);
    CF->Draw("P");
    fit->Draw("L SAME");

    if (lcms == 2) {
        stats->Draw();
    }
}
// ============================================
// Full LCMS pipeline
// ============================================


TH3D* Create3DFitHist(const FitResult r) {
    double fitLim = 0.2;
    TF3* fit = new TF3(
        "fit3d", CF_fit_3d, 
        -fitLim, fitLim, 
        -fitLim, fitLim, 
        -fitLim, fitLim, 
        7
    );
    fit->SetParameter(0, r.R[0] * r.R[0]);
    fit->SetParameter(1, r.R[1] * r.R[1]);
    fit->SetParameter(2, r.R[2] * r.R[2]);
    fit->SetParameter(3, r.R[3]);
    fit->SetParameter(4, r.R[4]);
    fit->SetParameter(5, r.R[5]);
    fit->SetParameter(6, r.lambda);

    for (int p = 0; p < 7; p++) fit->FixParameter(p, fit->GetParameter(p));

    auto fit3d = new TH3D(
        "fit3d","fit3d",
        80,-0.4,0.4,
        80,-0.4,0.4,
        80,-0.4,0.4
    );
    fit3d->SetDirectory(nullptr);

    for (int i=1;i<=80;i++)
    for (int j=1;j<=80;j++)
    for (int k=1;k<=80;k++) {
        double x = fit3d->GetXaxis()->GetBinCenter(i);
        double y = fit3d->GetYaxis()->GetBinCenter(j);
        double z = fit3d->GetZaxis()->GetBinCenter(k);
        fit3d->SetBinContent(i,j,k, fit->Eval(x,y,z));
    }
    return fit3d;
}

void SetStyle(Draw::Style* style) {
    // marker config
    style->marker.style = 20;
    style->marker.size = 0.5;
    style->marker.color = kBlack;

    // line config
    style->line.width = 1;
    style->line.color = kRed+1;
    style->line.style = 1;

    // label config
    style->label.size = 0.045;

    // title config
    style->title.size = 0.05;
    style->title.offset = 1.2;
}

void DrawCFOverFit(
    TCanvas* c,
    TH1D* CF, TH1D* fit, TPaveText* stats,
    std::string name, int lcms, Draw::Style style
) {
    LCMSAxis axis = LCMSAxis(lcms);
    TH1D* h_fit_over_cf = (TH1D*)CF->Clone("h_fit_over_cf");
    h_fit_over_cf->Divide(fit, CF);

    std::string name_fit_over_cf = name + " " + ToString(axis);
    Style1DCF(
        h_fit_over_cf, name_fit_over_cf.data(), ToString(axis).data(), style
    );
    h_fit_over_cf->GetYaxis()->SetRangeUser(0.9, 1.1);
    h_fit_over_cf->GetXaxis()->SetRangeUser(-0.2, 0.2);
    h_fit_over_cf->SetStats(0);

    c->cd(lcms+1);
    h_fit_over_cf->Draw("P");
    if (lcms == 2) {
        stats->Draw();
    }
}

void MakeLCMS1DProjections(
    Context ctx, TFile* in, TFile* out
) {
    Draw::Style style;
    SetStyle(&style);
    
    FitGrid fitRes = ctx.fitRes;
    Bin bin = ctx.bining;

    std::string dep_dir = ctx.outDir + "/dependency/fit_over_cf";
    std::string dir = ctx.outDir + "/all_1d_histos";
    EnsureDir(dir); EnsureDir(dep_dir);
    std::string format = "png";

    for (int chIdx = 0; chIdx < Charge::kCount; chIdx++)
    for (int centIdx = 0; centIdx < Centrality::kCount; centIdx++)
    for (int b = 0; b < bin.count; b++) {
        const FitResult r = fitRes[chIdx][centIdx][b];
        TH3D* fit3d = Create3DFitHist(r);
        TH3D* h_A = getNum(in, chIdx, centIdx, b);
        TH3D* h_A_wei = getNumWei(in, chIdx, centIdx, b);
        TPaveText* stats = GetFitStats(r, 0.032);

        // Save cf + fit in 1 canvas (3 pad for 3 axis)
        std::string cf_name = getCFName(chIdx, centIdx, bin.type, bin.names[b]);
        TCanvas* can = new TCanvas(cf_name.data(), cf_name.data(), 2400, 800);
        {
            can->SetTitle(cf_name.data());
            can->Divide(3,1);
        }

        // Save fit over cf in 1 canvas (3 pad for 3 axis)
        std::string n_fit_over_cf = Form(
            "fit/cf at charge=%s, centrality=%s, %s=%s",
            Charge::kNames[chIdx], 
            Centrality::kNames[centIdx],
            bin.type, bin.names[b]
        );
        TCanvas* c_fit_over_cf = new TCanvas(n_fit_over_cf.data(), n_fit_over_cf.data(), 2400, 800);
        {
            c_fit_over_cf->SetTitle(n_fit_over_cf.data());
            c_fit_over_cf->Divide(3, 1);
        }

        for (int lcms = 0; lcms < 3; lcms++) {
            LCMSAxis axis = LCMSAxis(lcms);

            auto [CF, fit] = Create1D(*h_A, *h_A_wei, fit3d, axis, cf_name);

            DrawCFOverFit(c_fit_over_cf, CF, fit, stats, n_fit_over_cf, lcms, style);
            SaveCanvasToFile(out, CF, fit, stats, style, cf_name, axis);
            DrawCFAndFit(can, CF, fit, stats, lcms);
        }

        {
            auto name = Form(
                    "%s/cfs_%s_%s_%s.%s", 
                    dir.data(), 
                    Charge::kNames[chIdx], Centrality::kNames[centIdx], bin.names[b], 
                    format.data()
                );
            SaveCanvasQuiet(can, name);
        }

        {
            auto name = Form(
                "%s/fit_over_cf_%s_%s_%s.%s",
                dep_dir.data(),
                Charge::kNames[chIdx], Centrality::kNames[centIdx], bin.names[b],
                format.data()
            );
            SaveCanvasQuiet(c_fit_over_cf, name);
        }
    }
}
