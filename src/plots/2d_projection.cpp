#include "plots.h"

#include <iostream>
#include <memory>

#include <TH3.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TSystem.h>

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
    const std::string& tag,
    TCanvas* canvas,
    const int y,
    bool draw
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
    auto name = tag + " " + ToString(ax1) + "-" + ToString(ax2);
    auto CF = RootPtr<TH2D>((TH2D*)num->Clone(name.c_str()));
    CF->Divide(num.get(), den.get());

    // --- crop to |q| < n
    Crop2D(*CF, 0.1);

    // --- cosmetics
    CF->SetTitle(name.data());
    CF->GetXaxis()->SetTitle(("q_{"+ToString(ax2)+"} [GeV/c]").c_str());
    CF->GetYaxis()->SetTitle(("q_{"+ToString(ax1)+"} [GeV/c]").c_str());
    CF->GetZaxis()->SetRangeUser(0.95,1.30);
    CF->SetStats(0);

    // --- draw
    TCanvas c((name).c_str(),"",650,600);
    CF->Draw("COLZ");

    outFile.cd();
    c.Write();

    if (draw) {
        canvas->cd(y+1);
        gPad->SetTicks(1,1);
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);

        auto* CF_clone = (TH2D*)CF->Clone(Form("%s_clone", CF->GetName()));
        CF_clone->Draw("COLZ");

        gPad->Modified();
        gPad->Update(); 
    }
}


// =====================================
// Full LCMS 2D pipeline
// =====================================

void MakeLCMS2DProjections(TFile* input, TFile* out) {

    std::vector<TCanvas*> canvases(chargeSize * centralitySize, nullptr);
    for (int ch = 0; ch < chargeSize; ch++)
    for (int centr = 0; centr < centralitySize; centr++) {
        auto name = Form(
            "all_out-long_2d_histos_centr_%s_%s",
            centralityNames[centr], chargeNames[ch]
        );
        auto title = Form(
            "CF_{out-long} at ch=%s, centrality=%s", 
            chargeNames[ch], centralityNames[centr]
        );
        canvases[ch*centralitySize + centr] = new TCanvas(name, title, 2000, 800);
        canvases[ch*centralitySize + centr]->Divide(5, 2);
    }

    for (int chIdx=0; chIdx < chargeSize; chIdx++)
    for (int centIdx =0; centIdx < centralitySize; centIdx++)
    for (int yIdx =0; yIdx < rapiditySize;yIdx++) {
        TH3D* A = getNum(input, chIdx, centIdx, yIdx);
        TH3D* Awei = getNumWei(input, chIdx, centIdx, yIdx);

        std::string cf_name = getCFName(chIdx, centIdx, yIdx); 

        Write2DProjection(*out, *A, *Awei, LCMSAxis::Out, LCMSAxis::Side, cf_name, 
            canvases[chIdx * centralitySize + centIdx], yIdx, false);
        Write2DProjection(*out, *A, *Awei, LCMSAxis::Out , LCMSAxis::Long, cf_name,
            canvases[chIdx * centralitySize + centIdx], yIdx, true);
        Write2DProjection(*out, *A, *Awei, LCMSAxis::Side, LCMSAxis::Long, cf_name,
            canvases[chIdx * centralitySize + centIdx], yIdx, false);
    }

    gSystem->mkdir("all_2d_histos"); 
    for (int ch = 0; ch < chargeSize; ch++)
    for (int centr = 0; centr < centralitySize; centr++) {
         auto name = Form(
            "all_2d_histos/all_out-long_2d_histos_centr_%s_%s.png", 
            centralityNames[centr], chargeNames[ch]
        );
        canvases[ch * centralitySize + centr]->SaveAs(name);
    }
}
