#pragma once

#include "config/config.h"

class TFile;
class TGraphErrors;

[[nodiscard]] TGraphErrors* build_chi2_ndf_graph(const Config& cfg, int ch, int centr);
[[nodiscard]] TGraphErrors* build_pvalue_graph(const Config& cfg, int ch, int centr);
[[nodiscard]] TGraphErrors* build_fit_over_cf_graph(const Config& cfg, TFile* cf3dFile, int ch,
                                                    int centr);
