#pragma once
#include <array>
#include <string>

#include <Rtypes.h>

// Sizes
inline constexpr int chargeSize = 2;
inline constexpr int centralitySize = 4;
inline constexpr int ktSize = 4;
inline constexpr int rapiditySize = 10;
inline constexpr int lcmsSize = 3;

inline constexpr std::array<double, 2> rapidityValues = {-1.0, 1.0};
inline constexpr std::array<double, 5> ktValues = {0.15, 0.25, 0.35, 0.45, 0.60};
inline constexpr std::array<int, 4> colors = {kRed, kBlue, kMagenta, kGreen};
inline constexpr std::array<int, 4> markers = {20, 21, 22, 23};
inline constexpr std::array<const char*, 3> axises = {"x", "y", "z"};
inline constexpr std::array<const char*, 3> LCMS = {"out", "side", "long"};
inline constexpr std::array<const char*, 2> chargeNames = {"neg", "pos"};
inline constexpr std::array<const char*, 4> centralityNames = {"0-10", "10-30", "30-50", "50-80"};
inline constexpr std::array<const char*, 4> ktNames = {"0.15-0.25", "0.25-0.35", "0.35-0.45", "0.45-0.60"};

// Enums
enum class LCMSAxis { Out, Side, Long };
inline const char* ToString(LCMSAxis a) {
    switch(a) {
        case LCMSAxis::Out:  return "out";
        case LCMSAxis::Side: return "side";
        case LCMSAxis::Long: return "long";
    }
    return "";
}