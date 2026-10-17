#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

#include <Math/MinimizerOptions.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH3D.h>
#include <TROOT.h>

#include "config/config.h"
#include "core/correlation.h"
#include "fit/model.h"

namespace
{
constexpr double kHbarC = 0.197;

double gaussian_bin_average(double lower, double upper, double radius)
{
    const double scale = radius / kHbarC;
    return std::sqrt(std::acos(-1.)) *
           (std::erf(scale * upper) - std::erf(scale * lower)) /
           (2 * scale * (upper - lower));
}

void check_float_variance_policy()
{
    TH1F sum_float("sum_float", "", 1, 0., 1.);
    TH1F count_float("count_float", "", 1, 0., 1.);
    TH1D sum_double("sum_double", "", 1, 0., 1.);
    TH1D count_double("count_double", "", 1, 0., 1.);
    TH1D cf_float("cf_float", "", 1, 0., 1.);
    TH1D cf_double("cf_double", "", 1, 0., 1.);
    sum_float.Sumw2();
    count_float.Sumw2();
    sum_double.Sumw2();
    count_double.Sumw2();
    constexpr int kPairs = 100000;
    for (int pair = 0; pair < kPairs; ++pair) {
        const double weight = pair % 2 == 0 ? 0.9 : 1.1;
        sum_float.Fill(0.5, weight);
        count_float.Fill(0.5, 1.);
        sum_double.Fill(0.5, weight);
        count_double.Fill(0.5, 1.);
    }
    fill_correlation(cf_float, sum_float, count_float);
    fill_correlation(cf_double, sum_double, count_double);
    const double n = count_float.GetBinContent(1);
    const double sum = sum_float.GetBinContent(1);
    const double residual = sum_float.GetSumw2()->At(1) - sum * sum / n;
    const double raw_error = std::sqrt(residual / n / (n - 1));
    std::cout << "float_variance_policy N=" << n << " residual=" << residual
              << " raw_sigma=" << raw_error << " program_float_sigma="
              << cf_float.GetBinError(1) << " program_double_sigma="
              << cf_double.GetBinError(1) << '\n';
    if (!(residual > 0 && cf_double.GetBinError(1) > 0)) {
        throw std::runtime_error("controlled variance fixture is invalid");
    }
}

void check_event_dependence()
{
    constexpr int kEvents = 100;
    constexpr int kPairsPerEvent = 20;
    TH1D sum("event_sum", "", 1, 0., 1.);
    TH1D count("event_count", "", 1, 0., 1.);
    TH1D cf("event_cf", "", 1, 0., 1.);
    sum.Sumw2();
    count.Sumw2();
    double event_sum_squares = 0;
    for (int event = 0; event < kEvents; ++event) {
        const double weight = event % 2 == 0 ? 0.7 : 1.3;
        event_sum_squares += (weight - 1) * (weight - 1);
        for (int pair = 0; pair < kPairsPerEvent; ++pair) {
            sum.Fill(0.5, weight);
            count.Fill(0.5, 1.);
        }
    }
    fill_correlation(cf, sum, count);
    const double event_error = std::sqrt(event_sum_squares / kEvents / (kEvents - 1));
    std::cout << "event_dependence events=" << kEvents << " pairs_per_event="
              << kPairsPerEvent << " pair_sigma=" << cf.GetBinError(1)
              << " event_sigma=" << event_error << " ratio="
              << event_error / cf.GetBinError(1) << '\n';
}

void check_matrix_constraint()
{
    FitResult result;
    result.r = {2., 2., 2., 6., 0., 0.};
    result.lambda = 0.8;
    result.ok = true;
    result.status = 0;
    result.cov_status = 3;
    result.ndf = 100;
    result.chi2 = 100;
    result.p_value = 0.5;
    // The coefficient matrix has eigenvalues -2, 4, 10 fm^2, although all
    // coefficients are inside the configured individual parameter bounds.
    std::cout << "matrix_constraint eigenvalues=-2,4,10 usable=" << is_usable_fit(result)
              << " F(0.1,-0.1,0)=" << eval_cf_3d(result, 0.1, -0.1, 0.)
              << " F(0)=" << eval_cf_3d(result, 0., 0., 0.) << '\n';
}

void check_finite_bins()
{
    constexpr double kOut = 5., kSide = 4., kLong = 6., kLambda = 0.8;
    TH3D histogram("bin_average", "", 20, -0.1, 0.1, 20, -0.1, 0.1, 20, -0.1, 0.1);
    for (int x = 1; x <= histogram.GetNbinsX(); ++x) {
        const double out = gaussian_bin_average(histogram.GetXaxis()->GetBinLowEdge(x),
                                               histogram.GetXaxis()->GetBinUpEdge(x), kOut);
        for (int y = 1; y <= histogram.GetNbinsY(); ++y) {
            const double side = gaussian_bin_average(histogram.GetYaxis()->GetBinLowEdge(y),
                                                    histogram.GetYaxis()->GetBinUpEdge(y), kSide);
            for (int z = 1; z <= histogram.GetNbinsZ(); ++z) {
                const double longitudinal = gaussian_bin_average(
                    histogram.GetZaxis()->GetBinLowEdge(z),
                    histogram.GetZaxis()->GetBinUpEdge(z), kLong);
                histogram.SetBinContent(x, y, z, 1 + kLambda * out * side * longitudinal);
                histogram.SetBinError(x, y, z, 0.001);
            }
        }
    }
    histogram.SetEntries(8000);
    Config cfg;
    cfg.input.type = "kt";
    cfg.fit.use_default_ip = true;
    const auto center_fit = fit_cf_3d_with_retry(&histogram, cfg, 0, 0, 0);
    cfg.fit.use_integral = true;
    const auto integral_fit = fit_cf_3d_with_retry(&histogram, cfg, 0, 0, 0);
    std::cout << "finite_bins truth=" << kOut << ',' << kSide << ',' << kLong << ',' << kLambda
              << " centers=" << center_fit.r[0] << ',' << center_fit.r[1] << ','
              << center_fit.r[2] << ',' << center_fit.lambda << " integrals="
              << integral_fit.r[0] << ',' << integral_fit.r[1] << ',' << integral_fit.r[2] << ','
              << integral_fit.lambda << " status=" << center_fit.status << ','
              << integral_fit.status << '\n';
    if (!is_usable_fit(center_fit) || !is_usable_fit(integral_fit) ||
        std::abs(integral_fit.r[0] - kOut) > 0.002 ||
        std::abs(integral_fit.r[1] - kSide) > 0.002 ||
        std::abs(integral_fit.r[2] - kLong) > 0.002 ||
        std::abs(integral_fit.lambda - kLambda) > 0.0001) {
        throw std::runtime_error("bin-integral closure check failed");
    }
}
} // namespace

int main()
{
    std::cout << std::setprecision(12) << "ROOT " << gROOT->GetVersion() << '\n';
    TH1::AddDirectory(false);
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    check_float_variance_policy();
    check_event_dependence();
    check_matrix_constraint();
    check_finite_bins();
    std::cout << "hbar_c_scale fit/producer=" << kHbarC / 0.1973
              << " relative_percent=" << 100 * (kHbarC / 0.1973 - 1) << '\n';
}
