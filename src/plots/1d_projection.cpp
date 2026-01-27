#include "plots.h"

#include <memory>

#include <TH3.h>
#include <TF3.h>
#include <TCanvas.h>
#include "TLegend.h"
#include "TDirectory.h"

#include "helpers.h"
#include "fit/fit.h"
#include "draw.h"

template<class T>
using RootPtr = std::unique_ptr<T>;

using std::string;

// ============================================
// Build dense 3D histogram of the fitted CF
// ============================================

RootPtr<TH3D> MakeFitHistogram(const TF3& fit) {
    auto h = RootPtr<TH3D>(new TH3D(
        "fit3d","fit3d",
        80,-0.4,0.4,
        80,-0.4,0.4,
        80,-0.4,0.4
    ));
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

RootPtr<TH1D> BuildLCMSFitFrom3D(
    const TH3D& slicedVolume,
    const TH3D& fit3D,
    LCMSAxis axis,
    const std::string& tag
) {
    const char* proj =
        axis == LCMSAxis::Out  ? "x" :
        axis == LCMSAxis::Side ? "y" :
                                 "z";

    auto den = RootPtr<TH3D>((TH3D*)slicedVolume.Clone());
    auto num = RootPtr<TH3D>((TH3D*)slicedVolume.Clone());
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

    auto hNum = RootPtr<TH1D>((TH1D*)num->Project3D(proj)->Clone());
    auto hDen = RootPtr<TH1D>((TH1D*)den->Project3D(proj)->Clone());

    auto fit = RootPtr<TH1D>((TH1D*)hDen->Clone(("fit_"+tag+"_"+ToString(axis)).c_str()));
    fit->Divide(hNum.get(), hDen.get());

    return fit;
}

// ============================================
// One LCMS projection → canvas in output ROOT
// ============================================

void Write1DProjection(
    TFile& outFile,
    const TH3D& A0,
    const TH3D& Awei0,
    const TH3D& Fit3D0,
    const LCMSAxis axis,
    const std::string tag
) {
    auto A    = RootPtr<TH3D>((TH3D*)A0.Clone());
    auto Awei = RootPtr<TH3D>((TH3D*)Awei0.Clone());
    A->SetDirectory(nullptr);
    Awei->SetDirectory(nullptr);

    SetSlice1D(*A,axis);
    SetSlice1D(*Awei,axis);

    const char* proj =
        axis == LCMSAxis::Out  ? "x" :
        axis == LCMSAxis::Side ? "y" :
                                 "z";

    auto hA    = RootPtr<TH1D>((TH1D*)A->Project3D(proj)->Clone());
    auto hAwei = RootPtr<TH1D>((TH1D*)Awei->Project3D(proj)->Clone());

    std::string name = tag + "_" + ToString(axis);

    auto CF = RootPtr<TH1D>((TH1D*)hA->Clone(name.c_str()));
    CF->Divide(hAwei.get(), hA.get());

    auto fit = BuildLCMSFitFrom3D(*A, Fit3D0, axis, tag);

    Style1DCF(CF.get(), name);
    StyleFit(fit.get());

    CF->GetYaxis()->SetRangeUser(0.9,1.7);
    CF->GetXaxis()->SetRangeUser(-0.2,0.2);
    fit->GetXaxis()->SetRangeUser(-0.2,0.2);

    TCanvas c(name.c_str(), name.c_str(), 800, 600);
    c.SetTicks(1,1);
    c.SetLeftMargin(0.12);
    c.SetBottomMargin(0.12);

    CF->Draw("P");
    fit->Draw("L SAME");

    TLegend leg(0.60,0.75,0.88,0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.AddEntry(CF.get(),"Data","pe");
    leg.AddEntry(fit.get(),"3D Gaussian fit","l");
    leg.Draw();

    outFile.cd();
    c.Write();
}

// ============================================
// Full LCMS pipeline
// ============================================

void MakeLCMS1DProjections(TFile* input, TFile* out, FitGrid& fitRes)
{
    auto fit = std::unique_ptr<TF3>(CreateCF3DFit());

    for (int chIdx = 0; chIdx < chargeSize; chIdx++)
    for (int centIdx = 0; centIdx < centralitySize; centIdx++)
    for (int yIdx = 0; yIdx < rapiditySize; yIdx++) {
        TH3D* h_A = getNum(input, chIdx, centIdx, yIdx);
        TH3D* h_A_wei = getNumWei(input, chIdx, centIdx, yIdx);

        const FitResult& r = fitRes[chIdx][centIdx][yIdx];

        fit->SetParameter(0, r.R[0]);
        fit->SetParameter(1, r.R[1]);
        fit->SetParameter(2, r.R[2]);
        fit->SetParameter(3, r.lambda);

        for (int p = 0; p < 4; p++)
            fit->FixParameter(p, fit->GetParameter(p));

        auto Fit3D = MakeFitHistogram(*fit);

        std::string cf_name = getCFName(chIdx, centIdx, yIdx);
        for (LCMSAxis ax : {LCMSAxis::Out, LCMSAxis::Side, LCMSAxis::Long})
            Write1DProjection(*out, *h_A, *h_A_wei, *Fit3D, ax, cf_name);
    }
}
