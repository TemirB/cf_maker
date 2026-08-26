#include "draw/draw.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TList.h>
#include <TMultiGraph.h>
#include <TPad.h>

namespace
{

inline constexpr std::array<Color_t, 4> kCentralityColors = {kRed, kBlue, kMagenta, kGreen};
inline constexpr std::array<Style_t, 4> kCentralityMarkers = {20, 21, 22, 23};

} // namespace

void set_range_with_errors(TMultiGraph* mg, double padFrac)
{
    if (!mg) {
        return;
    }

    const TList* list = mg->GetListOfGraphs();
    if (!list || list->GetSize() == 0) {
        return;
    }

    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymin = std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();

    for (TObject* obj : *list) {
        auto* graph = dynamic_cast<TGraphErrors*>(obj);
        if (!graph) {
            continue;
        }

        const int n = graph->GetN();
        for (int i = 0; i < n; ++i) {
            double x, y;
            graph->GetPoint(i, x, y);
            const double ex = graph->GetErrorX(i);
            const double ey = graph->GetErrorY(i);

            if (!std::isfinite(x) || !std::isfinite(y)) {
                continue;
            }

            xmin = std::min(xmin, x - ex);
            xmax = std::max(xmax, x + ex);
            ymin = std::min(ymin, y - ey);
            ymax = std::max(ymax, y + ey);
        }
    }

    if (!std::isfinite(xmin) || !std::isfinite(xmax) || !std::isfinite(ymin) ||
        !std::isfinite(ymax)) {
        return;
    }

    const double dx = (xmax - xmin);
    const double dy = (ymax - ymin);
    const double padx = (dx > 0 ? dx * padFrac : 1.0);
    const double pady = (dy > 0 ? dy * padFrac : 1.0);

    xmin -= padx;
    xmax += padx;
    ymin -= pady;
    ymax += pady;

    mg->GetXaxis()->SetLimits(xmin, xmax);
    mg->SetMinimum(ymin);
    mg->SetMaximum(ymax);
}

void write_mg_with_legend(TFile* file, TMultiGraph* mg, const char* canvasName, const char* xTitle,
                       const char* yTitle,
                       const std::vector<std::pair<TObject*, std::string>>& legendEntries,
                       draw::GraphKind kind)
{
    file->cd();

    TCanvas c(canvasName, canvasName, 1000, 800);
    mg->Draw("A");
    mg->GetXaxis()->SetTitle(xTitle);
    mg->GetYaxis()->SetTitle(yTitle);

    switch (kind) {
    case draw::GraphKind::Radii:
        mg->GetYaxis()->SetRangeUser(0., 10.);
        break;
    case draw::GraphKind::Cross:
        mg->GetYaxis()->SetRangeUser(-30, 30);
        break;
    case draw::GraphKind::Lambda:
        mg->GetYaxis()->SetRangeUser(0.7, 1.1);
        break;
    case draw::GraphKind::chi2_ndf:
        mg->GetYaxis()->SetRangeUser(0.1, 257500.0);
        c.SetLogy();
        break;
    case draw::GraphKind::FitOverCF:
        mg->GetYaxis()->SetRangeUser(0.8, 1.2);
        break;
    case draw::GraphKind::PValue:
        mg->GetYaxis()->SetRangeUser(0, 1);
        break;
    }

    TLegend leg(0.65, 0.65, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.03f);

    for (auto& [obj, label] : legendEntries) {
        leg.AddEntry(obj, label.c_str(), "lp");
    }
    leg.Draw();

    c.Write(canvasName);
}

void write_hist(TFile* file, TH1D* hist, const char* canvasName, const char* xTitle,
               const char* yTitle)
{
    file->cd();

    TCanvas c(canvasName, canvasName, 1000, 800);
    hist->Draw();
    hist->GetXaxis()->SetTitle(xTitle);
    hist->GetYaxis()->SetTitle(yTitle);

    c.Write(canvasName);
    hist->Write(hist->GetName(), TObject::kOverwrite);
}

void style_1d_cf(TH1* h, const std::string& name, const char* axis, draw::Style style)
{
    h->SetMarkerStyle(style.marker.style);
    h->SetMarkerSize(style.marker.size);
    h->SetMarkerColor(style.marker.color);
    h->SetLineColor(style.marker.color);
    h->SetLineWidth(style.line.width);

    h->SetTitle(name.c_str());

    auto xAxis = h->GetXaxis();
    auto yAxis = h->GetYaxis();

    {
        xAxis->CenterTitle();
        xAxis->SetTitle(Form("q_{%s} (GeV/c)", axis));
        xAxis->SetLabelSize(style.label.size);
        xAxis->SetTitleSize(style.title.size);
    }
    {
        yAxis->CenterTitle();
        yAxis->SetTitle(Form("CF(q_{%s})", axis));
        yAxis->SetLabelSize(style.label.size);
        yAxis->SetTitleSize(style.title.size);
        yAxis->SetTitleOffset(style.title.offset);
    }
}

void style_fit(TH1* h, draw::Style style)
{
    h->SetLineColor(style.line.color);
    h->SetLineWidth(style.line.width + 2);
    h->SetLineStyle(style.line.style);
}

void fix_margin(TCanvas* c, int pads)
{
    for (int i = 1; i <= pads; i++) {
        TPad* pad = dynamic_cast<TPad*>(c->GetPad(i));
        if (!pad) {
            continue;
        }
        pad->SetLeftMargin(0.15f);
        pad->SetBottomMargin(0.04f);
        pad->SetRightMargin(0.01f);
        pad->SetTopMargin(0.05f);
    }
}

TGraphErrors* make_styled_graph(const char* name, int centr)
{
    auto* graph = new TGraphErrors();
    graph->SetName(name);
    graph->SetLineColor(kCentralityColors[centr]);
    graph->SetMarkerColor(kCentralityColors[centr]);
    graph->SetMarkerStyle(kCentralityMarkers[centr]);
    return graph;
}

draw::Style draw::default_style()
{
    draw::Style style;
    style.marker.style = 20;
    style.marker.size = 0.5f;
    style.marker.color = kBlack;

    style.line.width = 1;
    style.line.color = static_cast<Color_t>(kRed + 1);
    style.line.style = 1;

    style.label.size = 0.045f;

    style.title.size = 0.05f;
    style.title.offset = 1.2f;

    return style;
}
