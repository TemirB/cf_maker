#pragma once

#include <TH1.h>
#include <TMultiGraph.h>
#include <TFile.h>

namespace Draw {
    struct Marker {
        Style_t style;
        Size_t size;
        Color_t color;
    };

    struct Line {
        Color_t color;
        Width_t width;
        Style_t style;
    };

    struct Label {
        Float_t size;
    };

    struct Title {
        std::string main;
        std::string yAxis;
        std::string xAxis;

        Float_t size;
        Float_t offset;
    };

    struct Style {
        Marker marker;
        Line line;
        Label label;
        Title title;
    };

    // static Style drawConfig =  ;
};

// Устанавливает диапазоны осей TMultiGraph
// с учётом ошибок TGraphErrors
void setRangeWithErrors(TMultiGraph* mg, double padFrac = 0.10);

void writeMGWithLegend(
    TFile* file,
    TMultiGraph* mg,
    const char* canvasName,
    const char* xTitle,
    const char* yTitle,
    const std::vector<std::pair<TObject*, std::string>>& legendEntries,
    const int type
);

void writeHist(
    TFile* file,
    TH1D* hist,
    const char* canvasName,
    const char* xTitle,
    const char* yTitle
);

void Style1DCF(
    TH1* h, std::string name, const char* axis,
    Draw::Style style
);

void StyleFit(TH1* h, Draw::Style style);