#pragma once

#include <array>

#include "helpers/constants.h"

struct FitResult {
    std::array<double,6> R{};
    std::array<double,6> eR{};
    double lambda{};
    double elambda{};
    double chi2{};
    int ndf{};
    double pvalue{};
    bool ok{false};

    // Методы
    bool IsFinite() const;
    double Chi2Ndf() const;
    bool IsValid() const;
    //
    void returnAllR(int i) const;
};

using FitGrid = FitResult[chargeSize][centralitySize][ktSize];

struct BadFitPoint {
    int charge;
    int cent;
    int kt;
    double ktVal;
    double chi2ndf;
    double pvalue;
    double lambda, elambda;
    double R[6], eR[6];
};