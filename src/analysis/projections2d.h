#pragma once

#include "config/config.h"

class TFile;

void make_lcms_2d_projections(const Config& cfg, TFile* in, TFile* out);
