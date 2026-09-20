#pragma once

#include "config/config.h"

class TFile;

void do_cf_ratios(Config& cfg, TFile* fCF3D, TFile* fRatioProj, TFile* fProjRatio);
