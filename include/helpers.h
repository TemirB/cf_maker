#pragma once

#include <string>
#include <TH3.h>
#include <TTree.h>
#include <TFile.h>
#include <TDirectory.h>

#include "fit/types.h"

namespace Charge {
    inline constexpr int kCount = 2;
    inline std::vector<const char*> kNames = { "pos", "neg" };
}

namespace Centrality {
    inline constexpr int kCount = 4;
    inline std::vector<const char*> kNames = { "0-10", "10-30", "30-50", "50-80" };
}

namespace Kt {
    inline constexpr int kCount = 4;
    inline std::vector<double> kValues = { 0.15, 0.25, 0.35, 0.45, 0.60 };
    inline std::vector<const char*> kNames = {"0.15-0.25", "0.25-0.35", "0.35-0.45", "0.45-0.60"};
}

namespace Rapidity {
    inline constexpr int kCount = 10;
    inline std::vector<double> kValues = { -1., -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1. };
    inline std::vector<const char*> kNames = {
        "[-1.0; -0.8]", "[-0.8; -0.6]", "[-0.6; -0.4]", "[-0.4; -0.2]", "[-0.2; 0.0]",
        "[0.0; 0.2]", "[0.2; 0.4]", "[0.4; 0.6]", "[0.6; 0.8]", "[0.8; 1.0]"
    };
}

namespace LCMS {
    inline constexpr int kCount = 6;
    inline std::vector<const char*> kNames = { "out", "side", "long", "out-side", "out-long", "side-long" };
}

struct Bin {
    const char* type;
    int count;
    std::vector<const char*> names;
    std::vector<double> values;
};

inline constexpr double kDefaultRange = 0.05;

inline constexpr std::array<int, 4> colors = {kRed, kBlue, kMagenta, kGreen};
inline constexpr std::array<int, 4> markers = {20, 21, 22, 23};

enum class LCMSAxis { Out, Side, Long };

enum class AnalysisType { Rapidity, Kt, Unknown };
AnalysisType getType(const char* type);

inline std::string ToString(LCMSAxis a) {
    switch(a) {
        case LCMSAxis::Out:  return "out";
        case LCMSAxis::Side: return "side";
        case LCMSAxis::Long: return "long";
    }
    return "";
}

using FitGrid = std::vector<std::vector<std::vector<FitResult>>>;

// Геттеры
TH3D* getNum(TFile* f, int charge, int cent, int ktIdx);
TH3D* getNumWei(TFile* f, int charge, int cent, int ktIdx);
std::string getCFName(int ch, int centr, const char* binType, const char* binName);

bool IsBadFit(const FitResult& r);

void EnsureDir(const std::string& dir);
