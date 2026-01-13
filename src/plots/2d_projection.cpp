#include "plots.h"

#include <TH3.h>
#include <TH2.h>
#include <TCanvas.h>

#include <helpers.h>

// =====================================
// Slice helper (works on CLONES only)
// =====================================

inline void SetSlice2D(TH3D& h, LCMSAxis freeze, double w = 0.2) {
    ResetRanges(h);
    FreezeAxis(h, freeze, w);
}

// =====================================
// One LCMS 2D projection
// =====================================

void Write2DProjection(
    TFile& outFile,
    const TH3D& A0, const TH3D& Awei0,
    LCMSAxis ax1, LCMSAxis ax2,
    const std::string& tag
) {
    // --- clone (ROOT ranges are global state!)
    TH3D* A     = (TH3D*)A0.Clone();
    TH3D* Awei  = (TH3D*)Awei0.Clone();

    A->SetDirectory(nullptr);
    Awei->SetDirectory(nullptr);

    // --- freeze third axis
    LCMSAxis freeze;
    if      ((ax1==LCMSAxis::Out  && ax2==LCMSAxis::Side) ||
             (ax1==LCMSAxis::Side && ax2==LCMSAxis::Out))  freeze = LCMSAxis::Long;
    else if ((ax1==LCMSAxis::Out  && ax2==LCMSAxis::Long) ||
             (ax1==LCMSAxis::Long && ax2==LCMSAxis::Out))  freeze = LCMSAxis::Side;
    else                                                   freeze = LCMSAxis::Out;

    SetSlice2D(*A,freeze);
    SetSlice2D(*Awei,freeze);

    std::string proj = ToString(ax1) + ToString(ax2);

    // --- project
    auto* num = (TH2D*)Awei->Project3D(proj.c_str());
    auto* den = (TH2D*)A->Project3D(proj.c_str());

    num->SetDirectory(nullptr);
    den->SetDirectory(nullptr);

    auto name = "CF_" + tag + "_" + ToString(ax1) + "_" + ToString(ax2);
    auto* CF = (TH2D*)num->Clone(name.c_str());
    CF->Divide(num,den);

    delete num;
    delete den;

    // --- axes cosmetics
    CF->GetXaxis()->SetTitle(("q_{"+ToString(ax2)+"} [GeV/c]").c_str());
    CF->GetYaxis()->SetTitle(("q_{"+ToString(ax1)+"} [GeV/c]").c_str());
    CF->GetZaxis()->SetRangeUser(0.95,1.30);

    // --- draw
    TCanvas c(("c_"+name).c_str(),"",650,600);
    CF->Draw("COLZ");

    outFile.cd();
    c.Write();

    delete CF;
}

// =====================================
// Full LCMS 2D pipeline
// =====================================

void MakeLCMS2DProjections(TFile* input, TFile* out) {

    for (int ch=0;ch<chargeSize;ch++)
    for (int c =0;c <centralitySize;c++)
    for (int k =0;k <ktSize;k++)
    for (int y =0;y <rapiditySize;y++) {

        auto tag = getPrefix(ch,c,y) + std::to_string(k);

        auto* A    = (TH3D*)input->Get(("bp_"+tag).c_str());
        auto* Awei = (TH3D*)input->Get(("bp_"+tag+"wei").c_str());
        if (!A || !Awei) continue;

        Write2DProjection(*out,*A,*Awei,LCMSAxis::Out ,LCMSAxis::Side,tag);
        Write2DProjection(*out,*A,*Awei,LCMSAxis::Out ,LCMSAxis::Long,tag);
        Write2DProjection(*out,*A,*Awei,LCMSAxis::Side,LCMSAxis::Long,tag);
    }
}
