#pragma once

#include "config/config.h"

class TFile;

// Restore compatible saved fits before any output is created when cf3d is disabled.
void prepare_stage_dependencies(Config& cfg);
void stage_cf3d(Config& cfg);
void stage_dependency(Config& cfg);
void stage_1d_projections(Config& cfg, TFile* input);
void stage_2d_projections(Config& cfg, TFile* input);
void stage_ratios(Config& cfg);
