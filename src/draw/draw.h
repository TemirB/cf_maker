#pragma once

#include <string>
#include <utility>
#include <vector>

#include <TH1.h>

class TCanvas;
class TFile;
class TGraphErrors;
class TH1D;
class TMultiGraph;

namespace draw
{
struct Marker
{
    Style_t style;
    Size_t size;
    Color_t color;
};

struct Line
{
    Color_t color;
    Width_t width;
    Style_t style;
};

struct Label
{
    Float_t size;
};

struct Title
{
    std::string main;
    std::string y_axis;
    std::string x_axis;

    Float_t size;
    Float_t offset;
};

struct Style
{
    Marker marker;
    Line line;
    Label label;
    Title title;
};

[[nodiscard]] Style default_style();

enum class GraphKind : char
{
    Radii,
    Cross,
    Lambda,
    chi2_ndf,
    FitOverCF,
    PValue
};
} // namespace draw

void set_range_with_errors(TMultiGraph* mg, double padFrac = 0.10);

void write_mg_with_legend(TFile* file, TMultiGraph* mg, const char* canvasName, const char* xTitle,
                       const char* yTitle,
                       const std::vector<std::pair<TObject*, std::string>>& legendEntries,
                       draw::GraphKind kind);

void write_hist(TFile* file, TH1D* hist, const char* canvasName, const char* xTitle,
               const char* yTitle);

void style_1d_cf(TH1* h, const std::string& name, const char* axis, draw::Style style);

void style_fit(TH1* h, draw::Style style);

void fix_margin(TCanvas* c, int pads);

[[nodiscard]] TGraphErrors* make_styled_graph(const char* name, int centr);
