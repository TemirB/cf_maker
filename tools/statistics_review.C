#include <cmath>
#include <iomanip>
#include <iostream>

#include <TH1D.h>
#include <TH1F.h>

// root -l -b -q tools/statistics_review.C
void statistics_review()
{
    std::cout << std::setprecision(12);
    TH1F sum_float("sum_float", "", 1, 0, 1);
    TH1D sum_double("sum_double", "", 1, 0, 1);
    sum_float.Sumw2();
    sum_double.Sumw2();
    constexpr int kCount = 10000;
    for (int pair = 0; pair < kCount; ++pair) {
        const double weight = pair % 2 == 0 ? 1.00100049 : .99900049;
        sum_float.Fill(.5, weight);
        sum_double.Fill(.5, weight);
    }
    for (TH1* histogram : {static_cast<TH1*>(&sum_float), static_cast<TH1*>(&sum_double)}) {
        const double sum = histogram->GetBinContent(1);
        const double q = histogram->GetSumw2()->At(1);
        const double residual = q - sum * sum / kCount;
        std::cout << histogram->ClassName() << " S=" << sum << " Q=" << q << " C=" << sum / kCount
                  << " Q-S*S/N=" << residual << " raw_sigma="
                  << (residual >= 0 ? std::sqrt(residual / kCount / (kCount - 1)) : -1) << '\n';
    }
    // Weighted means, not a pass/total efficiency. ROOT B is compared with
    // the unbiased variance of the mean for the very same pair weights.
    for (int example = 0; example < 3; ++example) {
        TH1D sum("sum", "", 1, 0, 1), count("count", "", 1, 0, 1);
        TH1D binomial("binomial", "", 1, 0, 1);
        sum.Sumw2();
        count.Sumw2();
        for (int pair = 0; pair < 100; ++pair) {
            const double weight =
                example == 2 ? 2. : (pair % 2 == 0 ? 0. : (example == 0 ? 1. : 2.));
            sum.Fill(.5, weight);
            count.Fill(.5);
        }
        binomial.Divide(&sum, &count, 1, 1, "B");
        const double mean = sum.GetBinContent(1) / 100;
        const double variance = (sum.GetSumw2()->At(1) - 100 * mean * mean) / 100 / 99;
        std::cout << "B example=" << example << " C=" << mean
                  << " sigma_B=" << binomial.GetBinError(1) << " sigma_mean=" << std::sqrt(variance)
                  << '\n';
    }
    double chi2 = 0;
    for (int bin = 0; bin < 5; ++bin) {
        const double value = 1. + .1 * bin;
        const double pull = (value - 1.2) / .1;
        chi2 += pull * pull;
        std::cout << "chi2 bin=" << bin + 1 << " C=" << value
                  << " model=1.2 sigma=.1 contribution=" << pull * pull << '\n';
    }
    std::cout << "chi2=" << chi2 << " ndf=5-1=4 chi2/ndf=" << chi2 / 4 << '\n';
    // Rounding a near-constant weighted sum can dominate its small variance.
    constexpr double kExactSum = 10000.0049;
    constexpr double kResidual = .01;
    const double q = kExactSum * kExactSum / kCount + kResidual;
    for (int digits : {2, 3, 4, 6}) {
        const double scale = std::pow(10., digits);
        const double sum = std::round(kExactSum * scale) / scale;
        const double residual = q - sum * sum / kCount;
        std::cout << "digits=" << digits << " S=" << sum << " C=" << sum / kCount
                  << " residual=" << residual
                  << " sigma=" << std::sqrt(residual / kCount / (kCount - 1)) << '\n';
    }
}
