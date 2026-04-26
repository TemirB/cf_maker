#include "plots.h"

#include <iostream>

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
    if (freeze == LCMSAxis::Out)  h.GetXaxis()->SetRangeUser(-w, w);
    if (freeze == LCMSAxis::Side) h.GetYaxis()->SetRangeUser(-w, w);
    if (freeze == LCMSAxis::Long) h.GetZaxis()->SetRangeUser(-w, w);
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
    // --- clone in (ROOT global state!)
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

void MakeLCMS2DProjections(Context ctx, TFile* in, TFile* out) {
    Bin bin = ctx.bining;
    std::string dir = ctx.outDir + "/all_2d_histos";
    EnsureDir(dir);
    std::string ext = "png";
    std::vector<TCanvas*> canvases(Charge::kCount * Centrality::kCount, nullptr);
    for (int ch = 0; ch < Charge::kCount; ch++)
    for (int centr = 0; centr < Centrality::kCount; centr++) {
        auto name = Form(
            "all_out-long_2d_histos_centr_%s_%s",
            Centrality::kNames[centr], Charge::kNames[ch]
        );
        auto title = Form(
            "CF_{out-long} at ch=%s, centrality=%s", 
            Charge::kNames[ch], Centrality::kNames[centr]
        );
        canvases[ch*Centrality::kCount + centr] = new TCanvas(name, title, 2000, 800);
        canvases[ch*Centrality::kCount + centr]->Divide(5, 2);
    }

    for (int chIdx=0; chIdx < Charge::kCount; chIdx++)
    for (int centIdx =0; centIdx < Centrality::kCount; centIdx++)
    for (int b =0; b < bin.count; b++) {
        auto [A, Awei] = getHists(in, chIdx, centIdx, b);

        std::string cf_name = getCFName(chIdx, centIdx, bin.type, bin.names[b]); 

        Write2DProjection(*out, *A, *Awei, LCMSAxis::Out, LCMSAxis::Side, cf_name, 
            canvases[chIdx * Centrality::kCount + centIdx], b, false);
        Write2DProjection(*out, *A, *Awei, LCMSAxis::Out , LCMSAxis::Long, cf_name,
            canvases[chIdx * Centrality::kCount + centIdx], b, true);
        Write2DProjection(*out, *A, *Awei, LCMSAxis::Side, LCMSAxis::Long, cf_name,
            canvases[chIdx * Centrality::kCount + centIdx], b, false);
    }

    for (int ch = 0; ch < Charge::kCount; ch++)
    for (int centr = 0; centr < Centrality::kCount; centr++) {
         auto name = Form(
            "%s/all_out-long_2d_histos_centr_%s_%s.%s", 
            dir.c_str(), Centrality::kNames[centr], Charge::kNames[ch], ext.c_str()
        );
        SaveCanvasQuiet(canvases[ch * Centrality::kCount + centr], name);
    }
}
