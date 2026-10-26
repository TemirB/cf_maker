#include <TH1D.h>
#include <TROOT.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Synthetic check of ROOT Divide("B") for the producer's filling contract:
// each accepted pair contributes w to A_wei and 1 to A in the same cell.
// The sample variance below assumes independent observations; it does not
// estimate correlations between pairs belonging to the same event.
void check_same_sample_errors()
{
    const std::vector<std::vector<double>> cases = {{0.5, 1.5}, {0., 1.}, {1.e-8}};
    std::cout << "ROOT " << gROOT->GetVersion() << '\n';
    std::cout << "case,N,S,Q,C,sigma_sample,sigma_ROOT_B,sigma_independent\n";
    std::cout << std::setprecision(12);

    for (std::size_t case_index = 0; case_index < cases.size(); ++case_index) {
        const auto name = "same_sample_" + std::to_string(case_index);
        TH1D weighted((name + "_weighted").c_str(), "", 1, 0., 1.);
        TH1D count((name + "_count").c_str(), "", 1, 0., 1.);
        TH1D binomial((name + "_B").c_str(), "", 1, 0., 1.);
        TH1D independent((name + "_independent").c_str(), "", 1, 0., 1.);
        weighted.SetDirectory(nullptr);
        count.SetDirectory(nullptr);
        binomial.SetDirectory(nullptr);
        independent.SetDirectory(nullptr);
        weighted.Sumw2();
        count.Sumw2();
        for (const double weight : cases[case_index]) {
            weighted.Fill(0.5, weight);
            count.Fill(0.5, 1.);
        }
        if (!binomial.Divide(&weighted, &count, 1., 1., "B") ||
            !independent.Divide(&weighted, &count)) {
            throw std::runtime_error("ROOT Divide failed");
        }

        const double n = count.GetBinContent(1);
        const double sum = weighted.GetBinContent(1);
        const double sum_squares = weighted.GetSumw2()->At(1);
        const double mean = sum / n;
        const double sample_error = n > 1 ? std::sqrt((sum_squares - sum * sum / n) / n / (n - 1))
                                          : std::numeric_limits<double>::quiet_NaN();
        const double expected_b_error =
            sum == n
                ? 0.
                : std::sqrt(std::abs(((1 - 2 * mean) * sum_squares + mean * mean * n) / n / n));
        if (std::abs(binomial.GetBinError(1) - expected_b_error) > 1.e-12) {
            throw std::runtime_error("ROOT B differs from the documented formula");
        }
        if (case_index == 0 && (sample_error != 0.5 || binomial.GetBinError(1) != 0.)) {
            throw std::runtime_error("continuous-weight example failed");
        }
        std::cout << case_index << ',' << n << ',' << sum << ',' << sum_squares << ',' << mean
                  << ',' << sample_error << ',' << binomial.GetBinError(1) << ','
                  << independent.GetBinError(1) << '\n';
    }
}
