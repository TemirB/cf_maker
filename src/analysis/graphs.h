#pragma once

#include "config/config.h"

class TFile;
class TGraphErrors;

// Boundary and failed fits remain in the saved fit metadata, but are excluded
// from ordinary physical graphs and model overlays.
[[nodiscard]] bool is_usable_fit(const FitResult& result);

[[nodiscard]] TGraphErrors* build_chi2_ndf_graph(const Config& cfg, int ch, int centr);
[[nodiscard]] TGraphErrors* build_pvalue_graph(const Config& cfg, int ch, int centr);
[[nodiscard]] TGraphErrors* build_fit_over_cf_graph(const Config& cfg, TFile* cf3dFile, int ch,
                                                    int centr);
