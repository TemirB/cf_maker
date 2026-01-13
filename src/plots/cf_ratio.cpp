#include "plots.h"

#include <TFile.h>
#include <TDirectory.h>
#include <TH3D.h>
#include <TH1D.h>

#include "fit/types.h"
#include "helpers.h"

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


void do_CF_ratios(TFile* fCF3D, TFile* fRatioProj, TFile* fProjRatio)
{
    for (int chIdx=0; chIdx<centralitySize; chIdx++)
    for (int ktIdx=0; ktIdx<ktSize; ktIdx++)
    for (int yIdx=0; yIdx<rapiditySize; yIdx++)
    {
        std::string neg = "h3d_CF_q_" + getPrefix(0, chIdx, yIdx)
                        + std::to_string(ktIdx) + "_weighted";
        std::string pos = "h3d_CF_q_" + getPrefix(1, chIdx, yIdx)
                        + std::to_string(ktIdx) + "_weighted";

        TH3D *hNeg = nullptr, *hPos = nullptr;
        fCF3D->GetObject(neg.c_str(), hNeg);
        fCF3D->GetObject(pos.c_str(), hPos);
        if (!hNeg || !hPos) continue;

        // Клонируем входные гистограммы
        TH3D* N = (TH3D*)hNeg->Clone();
        TH3D* P = (TH3D*)hPos->Clone();
        N->SetDirectory(nullptr);
        P->SetDirectory(nullptr);

        // 3D ratio
        TH3D* R = (TH3D*)N->Clone();
        R->SetDirectory(nullptr);
        R->Divide(N, P);

        for (auto ax : {LCMSAxis::Out, LCMSAxis::Side, LCMSAxis::Long})
        {
            auto name = Form("ratio_%d_%d_%d_%d",
                             chIdx, ktIdx, yIdx, (int)ax);

            // Ratio of projections
            TH1D* r1 = MakeRatio1D(*N, *P, ax, name);
            r1->SetDirectory(gDirectory);
            fRatioProj->cd();
            r1->Write();
            delete r1;

            // Projection of 3D ratio
            TH1D* r2 = ProjectRatio1D(*R, ax, name);
            r2->SetDirectory(gDirectory);
            fProjRatio->cd();
            r2->Write();
            delete r2;
        }

        delete R;
        delete N;
        delete P;
    }
}