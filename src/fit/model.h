#pragma once

#include <TF3.h>
#include <TH3.h>

#include "config/config.h"
#include "fit/types.h"

Double_t cf_fit_3d(Double_t* q, Double_t* par);

[[nodiscard]] double eval_cf_3d(const FitResult& r, double qOut, double qSide, double qLong);

[[nodiscard]] TF3* create_cf_3d_fit(const Config& cfg, int charge, int centrality, int y);

[[nodiscard]] FitResult fit_cf_3d(TH3D* hCF, TF3* fit3d, const FitConfig& fitCfg);

[[nodiscard]] FitResult fit_cf_3d_with_retry(TH3D* hCF, const Config& cfg, int charge, int centrality,
                                         int y);
