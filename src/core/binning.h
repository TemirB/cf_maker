#pragma once

#include <array>
#include <string>
#include <vector>

namespace charge
{
inline constexpr int kCount = 2;
inline constexpr std::array<const char*, 2> kNames = {"pos", "neg"};
} // namespace charge

namespace centrality
{
inline constexpr int kCount = 4;
inline constexpr std::array<const char*, 4> kNames = {"0-10", "10-30", "30-50", "50-80"};
} // namespace centrality

namespace kt
{
inline constexpr int kCount = 4;
inline constexpr std::array<double, 5> kValues = {0.15, 0.25, 0.35, 0.45, 0.60};
inline constexpr std::array<const char*, 4> kNames = {"0.15-0.25", "0.25-0.35", "0.35-0.45",
                                                      "0.45-0.60"};
inline constexpr std::array<const char*, 4> kFileNames = {"0.15-0.25", "0.25-0.35", "0.35-0.45",
                                                          "0.45-0.60"};
} // namespace kt

namespace rapidity
{
inline constexpr int kCount = 10;
inline constexpr std::array<double, 11> kValues = {-1., -0.8, -0.6, -0.4, -0.2, 0,
                                                   0.2, 0.4,  0.6,  0.8,  1.};
inline constexpr std::array<const char*, 10> kNames = {
    "[-1.0, -0.8]", "[-0.8, -0.6]", "[-0.6, -0.4]", "[-0.4, -0.2]", "[-0.2, 0.0]",
    "[0.0, 0.2]",   "[0.2, 0.4]",   "[0.4, 0.6]",   "[0.6, 0.8]",   "[0.8, 1.0]"};
inline constexpr std::array<const char*, 10> kFileNames = {
    "-1.0_-0.8", "-0.8_-0.6", "-0.6_-0.4", "-0.4_-0.2", "-0.2_0.0",
    "0.0_0.2",   "0.2_0.4",   "0.4_0.6",   "0.6_0.8",   "0.8_1.0"};
} // namespace rapidity

namespace lcms
{
inline constexpr int kCount = 6;
inline constexpr std::array<const char*, 6> kNames = {"out",      "side",     "long",
                                                      "out-side", "out-long", "side-long"};
} // namespace lcms

struct Bin
{
    int count = 0;
    std::vector<std::string> names;
    std::vector<std::string> file_names;
    std::vector<double> values;
};

[[nodiscard]] inline double bin_center(const Bin& bin, int i)
{
    return (bin.values[i] + bin.values[i + 1]) / 2.0;
}
