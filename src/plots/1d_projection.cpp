#include "plots.h"

#include <memory>

#include <TH3.h>
#include <TF3.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveText.h>
#include "TLegend.h"
#include "TDirectory.h"

#include "fit/fit.h"
#include "fit/badfit.h"
#include "helpers/root_utils.h"
#include "draw.h"

using std::string;

// ============================================
// Build dense 3D histogram of the fitted CF
// ============================================

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

TH1D* BuildLCMSFitFrom3D(
    const TH3D& slicedVolume,
    const TH3D& fit3D,
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
        num->SetBinContent(x,y,z, fit3D.GetBinContent(x,y,z));
    }

    auto hNum = (TH1D*)num->Project3D(proj)->Clone();
    auto hDen = (TH1D*)den->Project3D(proj)->Clone();

    auto fit = (TH1D*)hDen->Clone(("fit_"+tag+"_"+ToString(axis)).c_str());
    fit->Divide(hNum, hDen);

    return fit;
}

// ============================================
// One LCMS projection -> canvas in output ROOT
// ============================================

void AddFitStats(const FitResult& r, 
                 double x1=0.70, double y1=0.75, double x2=0.88, double y2=0.88) {
                    //0.60,0.75,0.88,0.88
    auto* stats = new TPaveText(x1, y1, x2, y2, "NDC");
    stats->SetBorderSize(1);
    stats->SetFillColor(kWhite);
    stats->SetTextAlign(12);
    stats->SetTextFont(42);
    stats->SetTextSize(0.032);
    
    stats->AddText(Form("#chi^{2}/ndf = %.1f / %d = %.3f", r.chi2, r.ndf, r.chi2/r.ndf));
    stats->AddText(Form("R_{out}  = %.3f #pm %.3f fm", r.R[0], r.eR[0]));
    stats->AddText(Form("R_{side} = %.3f #pm %.3f fm", r.R[1], r.eR[1]));
    stats->AddText(Form("R_{long} = %.3f #pm %.3f fm", r.R[2], r.eR[2]));
    stats->AddText(Form("#lambda  = %.3f #pm %.3f", r.lambda, r.elambda));
    stats->Draw();
}

void Write1DProjection(
    int charge, int centr, int kt,
    TFile& outFile,
    const TH3D& A0,
    const TH3D& Awei0,
    const TH3D& Fit3D0,
    const LCMSAxis axis,
    FitResult r,
    TCanvas* canvas
) {

    std::string tag = getCFName(charge, centr, kt);

    auto A    = (TH3D*)A0.Clone();
    auto Awei = (TH3D*)Awei0.Clone();
    A->SetDirectory(nullptr);
    Awei->SetDirectory(nullptr);

    SetSlice1D(*A,axis);
    SetSlice1D(*Awei,axis);

    const char* proj =
        axis == LCMSAxis::Out  ? "x" :
        axis == LCMSAxis::Side ? "y" :
                                 "z";

    auto hA    = (TH1D*)A->Project3D(proj)->Clone();
    auto hAwei = (TH1D*)Awei->Project3D(proj)->Clone();

    // std::string name = tag + "_" + ToString(axis);
    std::string name = Form(
        "CF at charge=%s, centrality=%s, kt=%s, axis=%s", 
        chargeNames[charge], centralityNames[centr], ktNames[kt], ToString(axis)
    );

    auto CF = (TH1D*)hA->Clone(name.c_str());
    CF->Divide(hAwei, hA, 1., 1., "B");

    auto fit = BuildLCMSFitFrom3D(*A, Fit3D0, axis, tag);

    Style1DCF(CF, name, ToString(axis));
    StyleFit(fit);

    CF->GetYaxis()->SetRangeUser(0.9,1.7);
    CF->GetXaxis()->SetRangeUser(-0.2,0.2);
    CF->SetStats(0);


    fit->GetXaxis()->SetRangeUser(-0.2,0.2);
    gStyle->SetOptFit(0102);

    TCanvas c(name.c_str(), name.c_str(), 800, 600);
    c.SetTicks(1,1);
    c.SetLeftMargin(0.12);
    c.SetBottomMargin(0.12);

    CF->Draw("P");
    fit->Draw("L SAME");
    AddFitStats(r);

    canvas->cd(4*centr + kt + 1);
    gPad->SetTicks(1,1);        // Риски со всех сторон
    gPad->SetLeftMargin(0.12);  // Отступы (важно для Divide!)
    gPad->SetBottomMargin(0.12);
    auto* CF_clone = (TH1D*)CF->Clone(Form("%s_clone", CF->GetName()));
    auto* fit_clone = (TH1D*)fit->Clone(Form("%s_clone", fit->GetName()));

    CF_clone->Draw("P");
    fit_clone->Draw("AL SAME");
    gPad->Modified();
    gPad->Update();             // Принудительное обновление пада

    outFile.cd();
    c.Write();
}

// ============================================
// Full LCMS pipeline
// ============================================

void MakeLCMS1DProjections(TFile* input, TFile* out, FitGrid& fitRes)
{
    std::vector<TCanvas*> canvases = {
        new TCanvas("all_out_1d_histos_pos", "All 1d(+) out projects", 800, 600),
        new TCanvas("all_out_1d_histos_neg", "All 1d(-) out projects", 800, 600),
        new TCanvas("all_side_1d_histos_pos", "All 1d(+) side projects", 800, 600),
        new TCanvas("all_side_1d_histos_neg", "All 1d(-) side projects", 800, 600),
        new TCanvas("all_long_1d_histos_pos", "All 1d(+) long projects", 800, 600),
        new TCanvas("all_long_1d_histos_neg", "All 1d(-) long projects", 800, 600)
    };

    for (int chIdx = 0; chIdx < chargeSize; chIdx++)
    for (int lcmsIdx = 0; lcmsIdx < 3; lcmsIdx++) {
        canvases[3 * chIdx + lcmsIdx]->Divide(4, 4);
    }

    for (int chIdx = 0; chIdx < chargeSize; chIdx++)
    for (int centIdx = 0; centIdx < centralitySize; centIdx++)
    for (int ktIdx = 0; ktIdx < ktSize; ktIdx++) {
        TH3D* h_A = getNum(input, chIdx, centIdx, ktIdx);
        TH3D* h_A_wei = getNumWei(input, chIdx, centIdx, ktIdx);

        auto fit = std::unique_ptr<TF3>(CreateCF3DFit(chIdx, centIdx, ktIdx));

        const FitResult& r = fitRes[chIdx][centIdx][ktIdx];

        fit->SetParameter(0, r.R[0]);
        fit->SetParameter(1, r.R[1]);
        fit->SetParameter(2, r.R[2]);
        fit->SetParameter(3, r.R[3]);
        fit->SetParameter(4, r.R[4]);
        fit->SetParameter(5, r.R[5]);

        fit->SetParameter(6, r.lambda);

        for (int p = 0; p < 7; p++)
            fit->FixParameter(p, fit->GetParameter(p));

        auto Fit3D = MakeFitHistogram(*fit);

        for (int lcmsIdx = 0; lcmsIdx < 3; lcmsIdx++)
            Write1DProjection(
                chIdx, centIdx, ktIdx,
                *out, *h_A, *h_A_wei, *Fit3D, LCMSAxis(lcmsIdx), 
                r, 
                canvases[3 * chIdx + lcmsIdx]
            );
    }

    out->cd();
    for (int chIdx = 0; chIdx < chargeSize; chIdx++)
    for (int lcmsIdx = 0; lcmsIdx < 3; lcmsIdx++) {
        canvases[3 * chIdx + lcmsIdx]->SaveAs(
            Form("all_%s_1d_histos_%s.pdf", LCMS[lcmsIdx], chargeNames[chIdx])
        );
    }
}
