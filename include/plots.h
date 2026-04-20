#pragma once

#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>

#include "helpers.h"
#include "context.h"
#include "fit/types.h"

void ResetRanges(TH3D& h);

void FreezeAxis(TH3D& h, LCMSAxis freeze, double w = 0.2);

// cf_3d.cpp
void BuildAndFit3DCorrelationFunctions(Context& ctx, TFile* inputFile, TFile* outFile);

// dependency
void MakeDependency(Context ctx, TFile* cf3dFile, TFile* outFile);

// 1d_projection.cpp
void MakeLCMS1DProjections(Context ctx, TFile* input, TFile* out);

// 2d_projection.cpp
void MakeLCMS2DProjections(Context ctx, TFile* input, TFile* out);

// ratios
void do_CF_ratios(Context ctx, TFile* fCF3D, TFile* fRatioProj, TFile* fProjRatio);

// result graphs
TGraphErrors* BuildFitOverCFGraph(const Context& ctx, TFile* cf3dFile, int ch, int centr);
TGraphErrors* BuildChi2NdfGraph(const Context& ctx, int ch, int centr);
