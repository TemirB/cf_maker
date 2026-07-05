#pragma once

#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>

#include "helpers.h"
#include "config.h"

void ResetRanges(TH3D& h);

void FreezeAxis(TH3D& h, LCMSAxis freeze, double w = 0.2);

// cf_3d.cpp
void BuildAndFit3DCorrelationFunctions(Config& cfg, TFile* inputFile, TFile* outFile);

// dependency
void MakeDependency(Config& cfg, TFile* cf3dFile, TFile* outFile);

// 1d_projection.cpp
void MakeLCMS1DProjections(Config& cfg, TFile* input, TFile* out);

// 2d_projection.cpp
void MakeLCMS2DProjections(Config& cfg, TFile* input, TFile* out);

// ratios
void do_CF_ratios(Config& cfg, TFile* fCF3D, TFile* fRatioProj, TFile* fProjRatio);

// result graphs
TGraphErrors* BuildFitOverCFGraph(const Config& cfg, TFile* cf3dFile, int ch, int centr);
TGraphErrors* BuildChi2NdfGraph(const Config& cfg, int ch, int centr);
TGraphErrors* BuildPvalue(const Config& cfg, int ch, int centr);