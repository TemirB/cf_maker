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
void BuildAndFit3DCorrelationFunctions(TFile* inputFile, TFile* outFile, FitGrid& fitRes, Bin bin);

// // kt_dependence.cpp
// void MakeKtDependence(TFile* outFile, FitGrid& fitRes);

// // rapidity_dependence.cpp
// void MakeRapidityDependence(TFile* outFile, FitGrid& fitRes);
void MakeDependency(TFile* outFile, FitGrid& fitRes, Bin bin);

// 1d_projection.cpp
void MakeLCMS1DProjections(TFile* input, TFile* out, FitGrid& fitRes, Bin bin);

// 2d_projection.cpp
void MakeLCMS2DProjections(TFile* input, TFile* out, Bin bin);

// ratios
void do_CF_ratios(TFile* fCF3D, TFile* fRatioProj, TFile* fProjRatio, Bin bin);

// void MakeBadFitMaps(TFile* outFile, TTree* badFits);

std::pair<int, int> GetBinRange(const TAxis* axis, double w);

// void MakeDependency(int binSize);

// class Dependency {
// public:
//     virtual ~Dependency() = default;
//     virtual void MakeDependency(TFile* outFile, FitGrid& fitRes) const = 0;
// };

// class Rapidity : public Dependency {
// public:
//     Rapidity() = default;

//     void MakeDependency(TFile* outFile, FitGrid& fitRes) const override;
// };

// class Kt : public Dependency {
// public:
//     Kt() = default;

//     void MakeDependency(TFile* outFile, FitGrid& fitRes) const override;
// };