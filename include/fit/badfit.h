#pragma once

#include <vector>

#include <TTree.h>
#include <TFile.h>

#include "fit/types.h"

extern std::vector<BadFitPoint> badPoints;

bool IsBadFit(const FitResult& r);
void CollectBadFits(FitGrid& fitRes, std::vector<BadFitPoint>& badPoints);
TTree* WriteBadFitTree(TFile* f, const std::vector<BadFitPoint>& badPoints);