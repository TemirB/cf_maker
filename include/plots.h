#pragma once

#include <TFile.h>
#include <TTree.h>

#include "helpers.h"
#include "fit/types.h"

void ResetRanges(TH3D& h);

void ProjectOnAxis(TH3D& h, LCMSAxis axis, double w = 0.05);

void FreezeAxis(TH3D& h, LCMSAxis freeze, double w = 0.2);

void SetSlice1D(TH3D& h, LCMSAxis axis, double w = 0.05);

// cf_3d.cpp
void BuildAndFit3DCorrelationFunctions(TFile* inputFile, TFile* outFile, FitGrid& fitRes);

// kt_dependence.cpp
void MakeKtDependence(TFile* outFile, FitGrid& fitRes);

// 1d_projection.cpp
void MakeLCMS1DProjections(TFile* input, TFile* out, FitGrid& fitRes);

// 2d_projection.cpp
void MakeLCMS2DProjections(TFile* input, TFile* out);

// ratios
void do_CF_ratios(TFile* fCF3D, TFile* fRatioProj, TFile* fProjRatio);

// void MakeBadFitMaps(TFile* outFile, TTree* badFits);