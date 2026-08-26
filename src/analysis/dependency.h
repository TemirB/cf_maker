#pragma once

#include "config/config.h"

class TFile;

void make_dependency(Config& cfg, TFile* cf3dFile, TFile* outFile);
