#include "plots.h"

#include <TH3.h>
#include <TF3.h>
#include <TCanvas.h>

#include "helpers.h"
#include "fit/fit.h"

using std::string;

// ============================================
// slicing helper (works on a CLONE only)
// ============================================

inline void SetSlice1D(TH3D& h, LCMSAxis axis, double w = 0.05) {
    h.GetXaxis()->SetRange(0,0);
    h.GetYaxis()->SetRange(0,0);
    h.GetZaxis()->SetRange(0,0);

    if (axis != LCMSAxis::Out)  h.GetXaxis()->SetRangeUser(-w, w);
    if (axis != LCMSAxis::Side) h.GetYaxis()->SetRangeUser(-w, w);
    if (axis != LCMSAxis::Long) h.GetZaxis()->SetRangeUser(-w, w);
}

// ============================================
// Build dense 3D histogram of the fitted CF
// ============================================

std::unique_ptr<TH3D> MakeFitHistogram(const TF3& fit) {
    auto h = std::make_unique<TH3D>(
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
        h->SetBinContent(i, j, k, fit.Eval(x, y, z));
    }
    return h;
}

// ============================================
// One LCMS projection → canvas in output ROOT
// ============================================

void Write1DProjection(
    TFile& outFile,
    const TH3D& A0,
    const TH3D& Awei0,
    const TH3D& Fit3D0,
    LCMSAxis axis,
    const string& tag
) {
    // --- clone everything (ROOT ranges are stateful!)
    auto A    = std::unique_ptr<TH3D>((TH3D*)A0.Clone());
    auto Awei = std::unique_ptr<TH3D>((TH3D*)Awei0.Clone());
    auto Fit3D = std::unique_ptr<TH3D>((TH3D*)Fit3D0.Clone());

    A->SetDirectory(nullptr);
    Awei->SetDirectory(nullptr);
    Fit3D->SetDirectory(nullptr);

    // --- slice in LCMS
    SetSlice1D(*A,axis);
    SetSlice1D(*Awei,axis);
    SetSlice1D(*Fit3D,axis);

    // --- project
    auto* num = Awei->Project3D(ToString(axis).c_str());
    auto* den = A->Project3D(ToString(axis).c_str());
    auto* fit = Fit3D->Project3D(ToString(axis).c_str());

    num->SetDirectory(nullptr);
    den->SetDirectory(nullptr);
    fit->SetDirectory(nullptr);

    // --- CF
    auto* CF = (TH1D*)num->Clone(("CF_"+tag+"_"+ToString(axis)).c_str());
    CF->Divide(num,den);

    delete num;
    delete den;

    // --- draw
    TCanvas c(("c_"+tag+"_"+ToString(axis)).c_str(),"",600,500);
    CF->Draw("E");
    fit->Draw("L SAME");

    outFile.cd();
    c.Write();

    delete fit;
    delete CF;
}

// ============================================
// Full LCMS pipeline
// ============================================

void MakeLCMS1DProjections(TFile* input, TFile* out, FitGrid& fitRes) {
    std::unique_ptr<TF3> fit(CreateCF3DFit());

    for (int ch=0;ch<chargeSize;ch++)
    for (int c =0;c <centralitySize;c++)
    for (int k =0;k <ktSize;k++)
    for (int y =0;y <rapiditySize;y++) {

        auto tag = getPrefix(ch,c,y) + std::to_string(k);

        auto* A    = (TH3D*)input->Get(("bp_"+tag).c_str());
        auto* Awei = (TH3D*)input->Get(("bp_"+tag+"wei").c_str());
        if (!A || !Awei) continue;

        const auto& r = fitRes[ch][c][k][y];
        fit->SetParameters(r.R[0],r.R[1],r.R[2],r.lambda);

        for (int p=0;p<4;p++)
            fit->FixParameter(p,fit->GetParameter(p));

        auto Fit3D = MakeFitHistogram(*fit);

        for (LCMSAxis ax : {LCMSAxis::Out,LCMSAxis::Side,LCMSAxis::Long})
            Write1DProjection(*out, *A, *Awei, *Fit3D, ax, tag);
    }
}
