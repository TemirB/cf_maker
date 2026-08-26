#pragma once

#include "config/config.h"

class TFile;

void make_lcms_1d_projections(Config& cfg, TFile* in, TFile* out);
