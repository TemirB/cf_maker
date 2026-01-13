#include "draw.h"

#include <TFile.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TList.h>
#include <TMultiGraph.h>

#include <limits>
#include <cmath>
#include <algorithm>

void setRangeWithErrors(TMultiGraph* mg, double padFrac) {
    if (!mg) return;

    const TList* list = mg->GetListOfGraphs();
    if (!list || list->GetSize() == 0) return;

    double xmin =  std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymin =  std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();

    for (TObject* obj : *list) {
        auto* g = dynamic_cast<TGraphErrors*>(obj);
        if (!g) continue;

        const int n = g->GetN();
        for (int i = 0; i < n; ++i) {
            double x, y;
            g->GetPoint(i, x, y);
            const double ex = g->GetErrorX(i);
            const double ey = g->GetErrorY(i);

            if (!std::isfinite(x) || !std::isfinite(y)) continue;

            xmin = std::min(xmin, x - ex);
            xmax = std::max(xmax, x + ex);
            ymin = std::min(ymin, y - ey);
            ymax = std::max(ymax, y + ey);
        }
    }

    if (!std::isfinite(xmin) || !std::isfinite(xmax) ||
        !std::isfinite(ymin) || !std::isfinite(ymax)) return;

    const double dx = (xmax - xmin);
    const double dy = (ymax - ymin);
    const double padx = (dx > 0 ? dx * padFrac : 1.0);
    const double pady = (dy > 0 ? dy * padFrac : 1.0);

    xmin -= padx; xmax += padx;
    ymin -= pady; ymax += pady;

    // Работает после Draw("A") или Draw("ALP")
    mg->GetXaxis()->SetLimits(xmin, xmax);
    mg->SetMinimum(ymin);
    mg->SetMaximum(ymax);
}

void writeMGWithLegend(
    TFile& out,
    TMultiGraph* mg,
    const char* name,
    const char* xtitle,
    const char* ytitle,
    const std::vector<std::pair<TObject*, std::string>>& legendEntries
) {
    // --- deep clone for ROOT
    auto* mgc = (TMultiGraph*)mg->Clone(name);

    // clone all graphs inside
    for (auto* obj : *mg->GetListOfGraphs()) {
        auto* g = (TGraph*)obj;
        auto* gc = (TGraph*)g->Clone();
        mgc->Add(gc, "lp");
    }

    TCanvas c(("c_"+std::string(name)).c_str(), "", 700, 600);
    mgc->Draw("A");

    mgc->GetXaxis()->SetTitle(xtitle);
    mgc->GetYaxis()->SetTitle(ytitle);

    TLegend leg(0.65,0.65,0.88,0.88);
    for (auto& [obj, label] : legendEntries)
        leg.AddEntry(obj, label.c_str(), "lp");

    leg.Draw();

    out.cd();
    mgc->Write(name);
    c.Write();

    delete mgc;   // canvas already owns its primitives
}


void writeHist(
    TFile* file,
    TH1D* hist,
    const char* canvasName,
    const char* xTitle,
    const char* yTitle
) {
    TCanvas c(canvasName, canvasName, 1000, 800);
    hist->Draw();
    hist->GetXaxis()->SetTitle(xTitle);
    hist->GetYaxis()->SetTitle(yTitle);

    file->cd();
    hist->Write(hist->GetName(), TObject::kOverwrite);
    c.Write(canvasName);
}