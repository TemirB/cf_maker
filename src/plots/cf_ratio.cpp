#include "plots.h"

#include <TFile.h>
#include <TDirectory.h>
#include <TH3D.h>
#include <TH1D.h>

#include "fit/types.h"
#include "helpers.h"
#include "analysis/context.h"

TH1D* Project1D(TH3D& h, LCMSAxis axis, double w = 0.05) {
    SetSlice1D(h, axis, w);

    const char* proj =
        axis==LCMSAxis::Out  ? "x" :
        axis==LCMSAxis::Side ? "yIdx" : "z";

    auto* out = (TH1D*)h.Project3D(proj);
    out->SetDirectory(nullptr);
    return out;
}

TH1D* MakeRatio1D(TH3D& neg, TH3D& pos, LCMSAxis axis, const char* name) {
    auto n = std::unique_ptr<TH1D>(Project1D(neg, axis));
    auto p = std::unique_ptr<TH1D>(Project1D(pos, axis));

    auto* r = (TH1D*)n->Clone(name);
    r->Divide(n.get(), p.get());
    return r;
}

TH1D* ProjectRatio1D(TH3D& ratio3D, LCMSAxis axis, const char* name) {
    auto r = Project1D(ratio3D, axis);
    r->SetName(name);
    return r;
}

void do_CF_ratios(RunContext& rCtx) {
    TFile& fRatioProj = rCtx.GetOut("ratio_proj");
    TFile& fProjRatio = rCtx.GetOut("ratio_1d");

    for (int chIdx=0;chIdx<centralitySize;chIdx++)
    for (int ktIdx=0;ktIdx<ktSize;ktIdx++)
    for (int yIdx=0;yIdx<rapiditySize;yIdx++) {

        std::string neg = "h3d_CF_q_" + getPrefix(0, chIdx, yIdx) + std::to_string(ktIdx) + "_weighted";
        std::string pos = "h3d_CF_q_" + getPrefix(1, chIdx, yIdx) + std::to_string(ktIdx) + "_weighted";

        TH3D *hNeg=nullptr, *hPos=nullptr;
        rCtx.input->GetObject(neg.c_str(),hNeg);
        rCtx.input->GetObject(pos.c_str(),hPos);
        if(!hNeg || !hPos) continue;

        auto N = std::unique_ptr<TH3D>((TH3D*)hNeg->Clone());
        auto P = std::unique_ptr<TH3D>((TH3D*)hPos->Clone());
        N->SetDirectory(nullptr);
        P->SetDirectory(nullptr);

        auto R = std::unique_ptr<TH3D>((TH3D*)N->Clone());
        R->Divide(N.get(),P.get());

        for(auto ax : {LCMSAxis::Out,LCMSAxis::Side,LCMSAxis::Long}) {
            auto name = Form("ratio_%d_%d_%d_%d", chIdx, ktIdx, yIdx, (int)ax);

            auto r1 = std::unique_ptr<TH1D>(MakeRatio1D(*N, *P, ax, name));
            fRatioProj.cd();
            r1->Write();

            auto r2 = std::unique_ptr<TH1D>(ProjectRatio1D(*R, ax, name));
            fProjRatio.cd();
            r2->Write();
        }
    }
}

REGISTER_STAGE("cf_ratios", do_CF_ratios);