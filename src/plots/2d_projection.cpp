#include "plots.h"

#include <iostream>
#include <memory>

#include <TH3.h>
#include <TH2.h>
#include <TCanvas.h>

#include <helpers.h>

template<class T>
using RootPtr = std::unique_ptr<T>;

inline char AxisChar(LCMSAxis ax) {
    if (ax == LCMSAxis::Out)  return 'x';
    if (ax == LCMSAxis::Side) return 'y';
    return 'z';
}

inline LCMSAxis ThirdAxis(LCMSAxis a, LCMSAxis b)
{
    if (a != LCMSAxis::Out  && b != LCMSAxis::Out)  return LCMSAxis::Out;
    if (a != LCMSAxis::Side && b != LCMSAxis::Side) return LCMSAxis::Side;
    return LCMSAxis::Long;
}

inline void Crop2D(TH2D& h, double q = 0.08)
{
    h.GetXaxis()->SetRangeUser(-q, q);
    h.GetYaxis()->SetRangeUser(-q, q);
}

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
    // --- clone input (ROOT global state!)
    auto A    = RootPtr<TH3D>((TH3D*)A0.Clone());
    auto Awei = RootPtr<TH3D>((TH3D*)Awei0.Clone());
    A->SetDirectory(nullptr);
    Awei->SetDirectory(nullptr);

    // --- freeze third LCMS axis
    LCMSAxis freeze = ThirdAxis(ax1, ax2);
    SetSlice2D(*A,freeze);
    SetSlice2D(*Awei,freeze);

    // --- build ROOT projection string
    std::string proj;
    proj += AxisChar(ax1);
    proj += AxisChar(ax2);

    // --- project safely
    auto* pNum = Awei->Project3D(proj.c_str());
    auto* pDen = A->Project3D(proj.c_str());
    if (!pNum || !pDen) {
        std::cerr << "Project3D failed for " << proj << " in " << tag << "\n";
        return;
    }

    auto num = RootPtr<TH2D>((TH2D*)pNum->Clone());
    auto den = RootPtr<TH2D>((TH2D*)pDen->Clone());
    num->SetDirectory(nullptr);
    den->SetDirectory(nullptr);

    // --- CF
    auto name = "CF_" + tag + "_" + ToString(ax1) + "_" + ToString(ax2);
    auto CF = RootPtr<TH2D>((TH2D*)num->Clone(name.c_str()));
    CF->Divide(num.get(), den.get());

    // --- crop to |q| < 0.05
    Crop2D(*CF, 0.05);

    // --- cosmetics
    CF->GetXaxis()->SetTitle(("q_{"+ToString(ax2)+"} [GeV/c]").c_str());
    CF->GetYaxis()->SetTitle(("q_{"+ToString(ax1)+"} [GeV/c]").c_str());
    CF->GetZaxis()->SetRangeUser(0.95,1.30);

    // --- draw
    TCanvas c(("c_"+name).c_str(),"",650,600);
    CF->Draw("COLZ");

    outFile.cd();
    c.Write();
}


// =====================================
// Full LCMS 2D pipeline
// =====================================

void MakeLCMS2DProjections(TFile* input, TFile* out) {

    for (int ch=0;ch<chargeSize;ch++)
    for (int c =0;c <centralitySize;c++)
    for (int k =0;k <ktSize;k++)
    for (int y =0;y <rapiditySize;y++) {
        std::string ending = getPrefix(ch, c, y);
        std::string name   = "bp_" + ending;

        TH3D* A     = (TH3D*) input->Get((name + std::to_string(k)).c_str());
        TH3D* Awei = (TH3D*) input->Get((name + "wei_" + std::to_string(k)).c_str());
        if (!A || !Awei) {
            std::cerr << "Warning: missing " << name << " kt=" << k << "\n";
            continue;
        }
        ending = ending + std::to_string(k);

        Write2DProjection(*out, *A, *Awei, LCMSAxis::Out, LCMSAxis::Side, ending);
        Write2DProjection(*out, *A, *Awei, LCMSAxis::Out , LCMSAxis::Long, ending);
        Write2DProjection(*out, *A, *Awei, LCMSAxis::Side, LCMSAxis::Long, ending);
    }
}
