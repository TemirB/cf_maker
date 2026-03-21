#pragma once

#include <TFile.h>
#include <TTree.h>

#include "helpers.h"
#include "fit/types.h"

void ResetRanges(TH3D& h);

void FreezeAxis(TH3D& h, LCMSAxis freeze, double w = 0.2);

// cf_3d.cpp
void BuildAndFit3DCorrelationFunctions(TFile* inputFile, TFile* outFile, FitGrid& fitRes, Bin bin);

// dependency
void MakeDependency(TFile* outFile, FitGrid& fitRes, Bin bin);

// 1d_projection.cpp
void MakeLCMS1DProjections(TFile* input, TFile* out, FitGrid& fitRes, Bin bin);

// 2d_projection.cpp
void MakeLCMS2DProjections(TFile* input, TFile* out, Bin bin);

// ratios
void do_CF_ratios(TFile* fCF3D, TFile* fRatioProj, TFile* fProjRatio, Bin bin);
