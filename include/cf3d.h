#pragma once

#include <TFile.h>

#include "fit.h"
#include "helpers.h"

void create_and_fit_3d(
    int chIdx,
    TFile* inputFile, TFile* outFile,
    FitResult (&fitRes)[chargeSize][centralitySize][ktSize][rapiditySize]
);