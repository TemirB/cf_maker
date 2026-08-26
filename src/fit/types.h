#pragma once

#include <array>
#include <vector>

struct FitResult
{
    std::array<double, 6> R{};
    std::array<double, 6> eR{};
    double lambda{};
    double e_lambda{};
    double chi2{};
    int ndf{};
    double p_value{};
    bool ok{false};
    int status = -1;
    int cov_status = -1;
    bool at_limit = false;
    int attempts = 0;
    std::array<double, 49> corr{};

    [[nodiscard]] bool is_finite() const;
    [[nodiscard]] double chi2_ndf() const;
    [[nodiscard]] bool is_valid() const;
    [[nodiscard]] double corr(int i, int j) const
    {
        return corr[static_cast<std::size_t>(i) * 7 + static_cast<std::size_t>(j)];
    }
};

[[nodiscard]] bool is_bad_fit(const FitResult& r);

using FitGrid = std::vector<std::vector<std::vector<FitResult>>>;
