#pragma once
#include <array>

struct FitResult {
    std::array<double,3> R{};
    std::array<double,3> eR{};
    double lambda{};
    double elambda{};
    double chi2{};
    int ndf{};
    double pvalue{};
    bool ok{false};

    bool IsFinite() const;
    double Chi2Ndf() const;
    bool IsValid() const;
};
