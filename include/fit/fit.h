#pragma once

#include <TF3.h>
#include <TH3.h>

#include "fit/types.h"

Double_t CF_fit_3d(Double_t* q, Double_t* par);

TF3* CreateCF3DFit(int charge, int centrality, int kt);

FitResult FitCF3D(TH3D* hCF, TF3* fit3d);