#pragma once

#include <TFile.h>

#include "fit/types.h"
#include "helpers.h"

void do_kt_diff(
    TFile* outFile,
    FitResult (&fitRes)[chargeSize][centralitySize][ktSize][rapiditySize]
);

void do_rapidity_diff(
    TFile* outFile,
    FitResult (&fitRes)[chargeSize][centralitySize][ktSize][rapiditySize]
);

