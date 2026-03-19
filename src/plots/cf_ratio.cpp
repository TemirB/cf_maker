#include "plots.h"

#include <TFile.h>
#include <TDirectory.h>
#include <TH3D.h>
#include <TH1D.h>
#include <TCanvas.h>

#include "fit/types.h"
#include "helpers.h"

TH1D* Project1D_internal(TH3D& h, LCMSAxis axis, double w = 0.05) {
    SetSlice1D(h, axis, w);

    const char* proj =
        axis==LCMSAxis::Out  ? "x" :
        axis==LCMSAxis::Side ? "y" : "z";

    TH1D* out = (TH1D*)h.Project3D(proj);
    out->SetDirectory(nullptr);
    return out;
}

// h - target hist
// fmt - formatted title string (ROOT syntax)
void Style(
    TH1D* h, const TString& fmt,
    LCMSAxis axis, int centr, int b
) {
    // General
    {
        TString axisStr = ToString(axis);
        TString title = TString::Format(
            "%s at centrality=[%s] %%, y=[%s] GeV/c and %s axis",
            fmt.Data(),
            Centrality::kNames[centr], 
            Rapidity::kNames[b], 
            axisStr.Data()
        );

        h->SetTitle(title);

        h->SetMarkerStyle(kFullCircle);
        h->SetMarkerSize(1.0);
        h->SetMarkerColor(kBlack);
        h->SetLineColor(kBlack);
        h->SetLineWidth(2);
    }    

    // X axis
    {
        TAxis* xAxis = h->GetXaxis();
        TString axisStr = ToString(axis);
        TString xTitle = "q_{" + axisStr + "} [GeV/c]";

        xAxis->SetTitle(xTitle);
        xAxis->CenterTitle();
        xAxis->SetTitleSize(0.05);
        xAxis->SetLabelSize(0.045);
        xAxis->SetRangeUser(-0.4, 0.4);
    }

    // Y axis
    {
        TAxis* yAxis = h->GetYaxis();
        yAxis->SetTitle(fmt);
        yAxis->CenterTitle();
        yAxis->SetTitleSize(0.05);
        yAxis->SetLabelSize(0.045);
        yAxis->SetTitleOffset(1.2);
        yAxis->SetRangeUser(0.5, 1.5);
    }
}

TH1D* ProjectRatio(
    TH3D& ratio,
    int centr, int b,
    LCMSAxis axis
) {
    TString name = TString::Format(
        "proj_of_ratios_%d_%d_%s",
        centr, b, ToString(axis).data()
    );

    TH1D* r = Project1D_internal(ratio, axis);
    r->SetName(name);

    TString axisStr = ToString(axis);
    TString fmt = "C^{++}/C^{--}_{" + axisStr + "}";
    
    Style(r, fmt, axis, centr, b);
    return r;
}

TH1D* RatioProject(
    TH3D& Neg, TH3D& Pos,
    int centr, int b,
    LCMSAxis axis
) {
    TString name = TString::Format(
        "ratio_proj_%d_%d_%s",
        centr, b, ToString(axis).data()
    );

    TH1D* n = Project1D_internal(Neg, axis);
    TH1D* p = Project1D_internal(Pos, axis);

    TH1D* ratio = (TH1D*)n->Clone(name);
    ratio->Divide(p, n);

    TString axisStr = ToString(axis);
    TString fmt = "C^{++}(q_{" + axisStr + "}) / C^{--}(q_{" + axisStr + "})";
    
    Style(ratio, fmt, axis, centr, b);
    return ratio;
}

void do_CF_ratios(TFile* fCF3D, TFile* fRatioProj, TFile* fProjRatio, Bin bin) {
    for (int centr = 0; centr < Centrality::kCount; centr++)
    for (int b = 0; b < bin.count; b++) {
        TString nPos = getCFName(1, centr, bin.type, bin.names[b]);
        TString nNeg = getCFName(0, centr, bin.type, bin.names[b]);

        TH3D* hNeg = (TH3D*) fCF3D->Get(nNeg);
        TH3D* hPos = (TH3D*) fCF3D->Get(nPos);

        if (!hNeg || !hPos) continue;

        TH3D* neg = (TH3D*)hNeg->Clone(); neg->SetDirectory(nullptr);
        TH3D* pos = (TH3D*)hPos->Clone(); pos->SetDirectory(nullptr);

        TH3D* ratio = (TH3D*)neg->Clone(); ratio->SetDirectory(nullptr);
        ratio->Reset();
        ratio->Divide(pos, neg);

        for (auto axis : {LCMSAxis::Out, LCMSAxis::Side, LCMSAxis::Long}) {
            // Отношение проекций
            {
                TH1D* h = RatioProject(*neg, *pos, centr, b, axis);
                fRatioProj->cd();
                h->Write();
            }

            // Проекция отношения
            {
                TH1D* h = ProjectRatio(*ratio, centr, b, axis);
                fProjRatio->cd();
                h->Write();
            }
        }
    }
}