#pragma once

#include <TFile.h>

#include "helpers.h"
#include "fit/types.h"

void BuildAndFit3DCorrelationFunctions( int chargeIndex, TFile* inputFile, TFile* outFile, FitGrid& fitRes);

void MakeKtDependence(TFile* outFile, FitGrid& fitRes);

void MakeRapidityDependence(TFile* outFile, FitGrid& fitRes);

void MakeLCMS1DProjections(TFile* input, TFile* out, FitGrid& fitRes);
