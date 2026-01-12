#include "draw.h"

#include <TFile.h>
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
    TFile* file,
    TMultiGraph* mg,
    const char* canvasName,
    const char* xTitle,
    const char* yTitle,
    const std::vector<std::pair<TObject*, std::string>>& legendEntries
) {
    file->cd();

    TCanvas c(canvasName, canvasName, 1000, 800);
    mg->Draw("A");                       // создаёт оси
    mg->GetXaxis()->SetTitle(xTitle);
    mg->GetYaxis()->SetTitle(yTitle);

    // твой helper для range (если надо)
    // setRangeWithErrors(mg, 0.1);

    TLegend leg(0.65, 0.65, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.03);

    for (auto& [obj, label] : legendEntries) {
        leg.AddEntry(obj, label.c_str(), "lp");
    }
    leg.Draw();

    c.Write(canvasName);                 // в root-файл попадёт картинка со всем
    // mg->Write(mg->GetName(), TObject::kOverwrite); // и сам mg отдельно, если нужно
}

void writeHist(
    TFile* file,
    TH1D* hist,
    const char* canvasName,
    const char* xTitle,
    const char* yTitle
) {
    file->cd();

    TCanvas c(canvasName, canvasName, 1000, 800);
    hist->Draw();
    hist->GetXaxis()->SetTitle(xTitle);
    hist->GetYaxis()->SetTitle(yTitle);

    c.Write(canvasName);
    hist->Write(hist->GetName(), TObject::kOverwrite);
}