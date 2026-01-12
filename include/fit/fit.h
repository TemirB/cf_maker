#pragma once

#include <TF3.h>
#include <TH3.h>

#include "fit/result.h"

// 3D Gaussian CF model
Double_t CF_fit_3d(Double_t* q, Double_t* par);

// Factory for ROOT TF3
TF3* CreateCF3DFit();

// Run fit and extract physics result
FitResult FitCF3D(TH3D* hCF, TF3* fit3d);
