#pragma once

#include <vector>
#include <string>

#include <TH3.h>
#include <TFile.h>

std::string getCFName(int chIdx, int centIdx, int ktIdx);
TH3D* getNum(TFile* f, int charge, int cent, int ktIdx);
TH3D* getNumWei(TFile* f, int charge, int cent, int ktIdx);