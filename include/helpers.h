#pragma once

#include <string>
#include <TH3.h>
#include <TDirectory.h>

inline constexpr int chargeSize = 2;
inline constexpr int centralitySize = 4;
inline constexpr int ktSize = 4;
inline constexpr int rapiditySize = 10;
inline constexpr int lcmsSize = 3;

const std::array<double, 2> rapidityValues = {-1.0, 1.0};
const std::array<double, 5>      ktValues        = {0.15, 0.25, 0.35, 0.45, 0.60};
const std::array<int, 4> colors = {kRed, kBlue, kMagenta, kGreen};
const std::array<int, 4> markers = {20, 21, 22, 23};

const std::array<std::string, 3> axises          = {"x", "y", "z"};
const std::array<std::string, 3> LCMS            = {"out", "side", "long"};
const std::array<std::string, 2> chargeNames     = {"neg", "pos"};
const std::array<std::string, 4> centralityNames = {"0-10", "10-30", "30-50", "50-80"};
const std::array<std::string, 4> ktNames         = {"0.15-0.25", "0.25-0.35", "0.35-0.45", "0.45-0.60"};

int getIdx(int chIdx, int centIdx, int ktIdx, int yIdx);

std::string getPrefix(int chIdx, int centIdx, int yIdx);

bool sameBinning(const TH3D* a, const TH3D* b);

TDirectory* getOrMakeDir(TDirectory* base, const std::string& name);

void EnsureDir(const std::string& dir);

std::string GetExeDir();

std::string Sha256OfFile(const std::string& fname);
int FindCentralityIndex(const std::string& name);
int FindKtIndex(double low, double high);