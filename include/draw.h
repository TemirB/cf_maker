#pragma once

#include <TH1.h>
#include <TMultiGraph.h>
#include <TFile.h>

// Устанавливает диапазоны осей TMultiGraph
// с учётом ошибок TGraphErrors
void setRangeWithErrors(TMultiGraph* mg, double padFrac = 0.10);

void writeMGWithLegend(
    TFile* file,
    TMultiGraph* mg,
    const char* canvasName,
    const char* xTitle,
    const char* yTitle,
    const std::vector<std::pair<TObject*, std::string>>& legendEntries
);

void writeHist(
    TFile* file,
    TH1D* hist,
    const char* canvasName,
    const char* xTitle,
    const char* yTitle
);

void Style1DCF(TH1* h, std::string name);

void StyleFit(TH1* h);