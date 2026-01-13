#pragma once

#include <TFile.h>
#include <TTree.h>

#include "helpers.h"
#include "fit/types.h"

static void ResetRanges(TH3D& h) {
    h.GetXaxis()->SetRange(0,0);
    h.GetYaxis()->SetRange(0,0);
    h.GetZaxis()->SetRange(0,0);
};

static void ProjectOnAxis(TH3D& h, LCMSAxis axis, double w = 0.05) {
    if (axis != LCMSAxis::Out)  h.GetXaxis()->SetRangeUser(-w, w);
    if (axis != LCMSAxis::Side) h.GetYaxis()->SetRangeUser(-w, w);
    if (axis != LCMSAxis::Long) h.GetZaxis()->SetRangeUser(-w, w);
};

static void FreezeAxis(TH3D& h, LCMSAxis freeze, double w = 0.2) {
    if (freeze == LCMSAxis::Out)  h.GetXaxis()->SetRangeUser(-w, w);
    if (freeze == LCMSAxis::Side) h.GetYaxis()->SetRangeUser(-w, w);
    if (freeze == LCMSAxis::Long) h.GetZaxis()->SetRangeUser(-w, w);
}

static void SetSlice1D(TH3D& h, LCMSAxis axis, double w = 0.05) {
    ResetRanges(h);
    ProjectOnAxis(h, axis, w);
}

void BuildAndFit3DCorrelationFunctions( int chargeIndex, TFile* inputFile, TFile* outFile, FitGrid& fitRes);

void MakeKtDependence(TFile* outFile, FitGrid& fitRes);

TTree* MakeRapidityDependence(TFile* outFile, FitGrid& fitRes);

void MakeLCMS1DProjections(TFile* input, TFile* out, FitGrid& fitRes);
void MakeLCMS2DProjections(TFile* input, TFile* out);

void do_CF_ratios(TFile* fCF3D, TFile* fRatioProj, TFile* fProjRatio);

void MakeBadFitMaps(TFile* outFile, TTree* badFits);