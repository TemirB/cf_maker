#include "plots.h"

#include <TFile.h>
#include <TDirectory.h>
#include <TH3D.h>
#include <TH1D.h>

#include "fit/types.h"
#include "helpers.h"

template<class T>
using RootPtr = std::unique_ptr<T>;

TH1D* Project1D(TH3D& h, LCMSAxis axis, double w = 0.05) {
    SetSlice1D(h, axis, w);

    const char* proj =
        axis==LCMSAxis::Out  ? "x" :
        axis==LCMSAxis::Side ? "yIdx" : "z";

    TH1D* out = (TH1D*)h.Project3D(proj);
    out->SetDirectory(nullptr);
    return out;
}


TH1D* MakeRatio1D(TH3D& neg, TH3D& pos, LCMSAxis axis, const char* name) {
    TH1D* n = Project1D(neg, axis);
    TH1D* p = Project1D(pos, axis);

    TH1D* r = (TH1D*)n->Clone(name);
    r->SetDirectory(nullptr);
    r->Divide(n, p);

    delete n;
    delete p;

    return r;
}


TH1D* ProjectRatio1D(TH3D& ratio3D, LCMSAxis axis, const char* name) {
    TH1D* r = Project1D(ratio3D, axis);
    r->SetName(name);
    return r;
}

TH1D* MakeProject1D(
    const TH3D& Neg, const TH3D& Pos,
    const LCMSAxis axis, std::string name
) {
    auto neg = RootPtr<TH3D>((TH3D*)Neg.Clone()); neg->SetDirectory(nullptr);
    auto pos = RootPtr<TH3D>((TH3D*)Pos.Clone()); pos->SetDirectory(nullptr);

    SetSlice1D(*neg, axis);
    SetSlice1D(*pos, axis);

    const char* proj =
        axis == LCMSAxis::Out  ? "x" :
        axis == LCMSAxis::Side ? "y" :
                                 "z";

    auto neg1D = RootPtr<TH1D>((TH1D*)neg->Project3D(proj)->Clone());
    auto pos1D = RootPtr<TH1D>((TH1D*)neg->Project3D(proj)->Clone());

    auto ratio = RootPtr<TH1D>((TH1D*)neg1D->Clone(name.c_str()));
    ratio->Divide(pos1D.get(), neg1D.get());

    return ratio.get();
}

void Make_Ratio_of_1D_Projection(
    TFile& output, 
    TH1D& hNeg, TH1D& hPos, 
    std::string name
) {
    auto pos = RootPtr<TH1D>((TH1D*)hPos.Clone()); pos->SetDirectory(nullptr);
    auto neg = RootPtr<TH1D>((TH1D*)hNeg.Clone()); neg->SetDirectory(nullptr);
    
    auto ratio = RootPtr<TH1D>((TH1D*)pos->Clone(name.c_str()));
    ratio->Divide(neg.get());

    /*
        Красоту делаем тут
        Подписи маштабы etc etc etc
    */
    
    output.cd();
    ratio->Write();
}

void do_CF_ratios(TFile* fCF3D, TFile* fRatioProj, TFile* fProjRatio) {
    for (int centIdx = 0; centIdx < centralitySize; centIdx++)
    for (int ktIdx = 0; ktIdx < ktSize; ktIdx++) {
        std::string nPos = getCFName(1, centIdx, ktIdx);
        std::string nNeg = getCFName(0, centIdx, ktIdx);

        TH3D *hNeg = nullptr, *hPos = nullptr;
        fCF3D->GetObject(nNeg.c_str(), hNeg);
        fCF3D->GetObject(nPos.c_str(), hPos);

        if (!hNeg || !hPos) continue;

        TH3D* neg = (TH3D*)hNeg->Clone(); neg->SetDirectory(nullptr);
        TH3D* pos = (TH3D*)hPos->Clone(); pos->SetDirectory(nullptr);

        TH3D* ratio = (TH3D*)neg->Clone(); ratio->SetDirectory(nullptr);
        ratio->Divide(neg, pos);

        for (auto ax : {LCMSAxis::Out, LCMSAxis::Side, LCMSAxis::Long}) {
            std::string proj_name = Form("ratios_of_proj_%d_%d_%s", centIdx, ktIdx, ToString(ax).c_str());
            TH1D* ratio_of_proj = MakeProject1D(*neg, *pos, ax, proj_name);

            fRatioProj->cd();
            ratio_of_proj->Write();

            std::string proj_of_ratio_name = Form("proj_of_ratios_%d_%d_%s", centIdx, ktIdx, ToString(ax).c_str());
            TH1D* proj_of_ratio = ProjectRatio1D(*ratio, ax, proj_of_ratio_name.c_str());

            fProjRatio->cd();
            proj_of_ratio->Write();
        }


        // for (auto ax : {LCMSAxis::Out, LCMSAxis::Side, LCMSAxis::Long}) {
        //     auto name = Form("ratio_%d_%d_%d", centIdx, ktIdx, (int)ax);

        //     // // Ratio of projections
        //     // TH1D* r1 = MakeRatio1D(*N, *P, ax, name);
        //     // r1->SetDirectory(gDirectory);
        //     // fRatioProj->cd();
        //     // r1->Write();
        //     // delete r1;

        //     // Projection of 3D ratio
        //     TH1D* r2 = ProjectRatio1D(*R, ax, name);
        //     r2->SetDirectory(gDirectory);
        //     fProjRatio->cd();
        //     r2->Write();
        //     delete r2;
        // }

        // delete R;
        // delete N;
        // delete P;
    }
}