#pragma once

#include <TDirectory.h>
#include <TH3.h>

#include "draw.h"

void ratio_proj(
    TDirectory* dir,
    TH3D* h3d_neg, TH3D* h3d_pos,
    const std::string& axis,
    const char* name
);

void proj_ratio(
    TDirectory* dir,
    TH3D* h3d,
    const std::string& axis,
    const char* name
);

void do_CF_ratios(
    TFile* outFile,
    FitGrid& fitRes
);