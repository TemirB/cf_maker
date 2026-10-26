#pragma once

#include "config/config.h"

class TFile;

// The complete fit grid and its provenance are stored inside cf3d.root.
// Loading rejects legacy files, altered inputs and incompatible fit configurations.
void write_fit_results(TFile& file, const Config& cfg);
void read_fit_results(TFile& file, Config& cfg);
