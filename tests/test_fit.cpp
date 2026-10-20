#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>

#include <Math/MinimizerOptions.h>
#include <TH3D.h>
#include <TMath.h>

#include "config/config.h"
#include "fit/model.h"

namespace
{

#define CHECK(cond)                                                                                \
    do {                                                                                           \
        if (!(cond)) {                                                                             \
            std::cerr << "FAILED: " #cond " at " << __FILE__ << ":" << __LINE__ << "\n";           \
            std::exit(1);                                                                          \
        }                                                                                          \
    } while (0)

void check_manual_statistics()
{
    // Five points fitted by their mean 1.2, with errors 0.1:
    // squared pulls = 4, 1, 0, 1, 4 -> chi2=10, ndf=4, p=6*exp(-5).
    TH3D histogram("manual_statistics", "", 7, -0.35, 0.35, 2, -0.1, 0.1, 1, -0.1, 0.1);
    for (int x = 1; x <= 7; ++x) {
        const bool in_range = x >= 2 && x <= 6;
        histogram.SetBinContent(x, 1, 1, in_range ? 1.0 + 0.1 * (x - 2) : 100.0);
        histogram.SetBinError(x, 1, 1, 0.1);
        // Nonempty bins with zero errors must also be excluded.
        histogram.SetBinContent(x, 2, 1, 100.0);
        histogram.SetBinError(x, 2, 1, 0.0);
    }
    histogram.SetBinContent(0, 1, 1, 1000.0);
    histogram.SetBinError(0, 1, 1, 0.1);
    TF3 constant("manual_constant", "[0]+0*x+0*y+0*z", -0.21, 0.21, -0.1, 0.1, -0.1, 0.1);
    constant.SetParameter(0, 1.2);
    FitConfig cfg;
    const auto statistics = calculate_fit_statistics(histogram, constant, cfg);
    CHECK(statistics);
    CHECK(std::abs(statistics->chi2 - 10.0) < 1e-10);
    CHECK(statistics->ndf == 4);
    CHECK(std::abs(statistics->p_value - 6.0 * std::exp(-5.0)) < 1e-12);

    constant.FixParameter(0, 1.2);
    const auto fixed = calculate_fit_statistics(histogram, constant, cfg);
    CHECK(fixed && fixed->ndf == 5);
    histogram.GetXaxis()->SetRange(3, 5);
    const auto restricted = calculate_fit_statistics(histogram, constant, cfg);
    CHECK(restricted && restricted->ndf == 3);
    CHECK(std::abs(restricted->chi2 - 2.0) < 1e-10);
    histogram.GetXaxis()->SetRange(0, 0);
    cfg.options = "QS0";
    const auto full_range = calculate_fit_statistics(histogram, constant, cfg);
    CHECK(full_range && full_range->ndf == 7 && full_range->chi2 > 1e6);
    cfg.options = "RQS0SERIAL";
    const auto serial = calculate_fit_statistics(histogram, constant, cfg);
    CHECK(serial && serial->ndf == 5 && std::abs(serial->chi2 - 10.0) < 1e-10);

    // Unweighted least squares uses unit errors, including zero-error bins
    // with nonzero content; WW additionally includes empty bins.
    histogram.SetBinContent(4, 2, 1, 0.0);
    cfg.options = "RQSW0";
    const auto unit_errors = calculate_fit_statistics(histogram, constant, cfg);
    CHECK(unit_errors && unit_errors->ndf == 9);
    CHECK(std::abs(unit_errors->chi2 - (0.1 + 4.0 * 98.8 * 98.8)) < 1e-8);
    cfg.options = "RQSWW0";
    const auto empty_included = calculate_fit_statistics(histogram, constant, cfg);
    CHECK(empty_included && empty_included->ndf == 10);
    CHECK(std::abs(empty_included->chi2 - unit_errors->chi2 - 1.2 * 1.2) < 1e-8);
    for (const auto* options : {"RQLS0", "RQPS0", "RQUS0", "RQSWIDTH0"}) {
        cfg.options = options;
        CHECK(!calculate_fit_statistics(histogram, constant, cfg));
    }

    TH3D no_errors("no_statistics_errors", "", 2, -0.1, 0.1, 1, -0.1, 0.1, 1, -0.1, 0.1);
    cfg.options = "RQS0";
    const auto empty = calculate_fit_statistics(no_errors, constant, cfg);
    CHECK(empty && empty->chi2 == 0.0 && empty->ndf == 0 && empty->p_value == 0.0);
}

void check_integral_statistics()
{
    const double edges[] = {-0.5, -0.2, 0.1, 0.5};
    const double y_edges[] = {-0.1, 0.2};
    const double z_edges[] = {-0.2, 0.3};
    TH3D histogram("integral_statistics", "", 3, edges, 1, y_edges, 1, z_edges);
    const double deviations[] = {0.1, -0.2, 0.3};
    for (int x = 1; x <= 3; ++x) {
        const double lower = edges[x - 1], upper = edges[x];
        const double average = 1.2 + (lower * lower + lower * upper + upper * upper) / 3.0;
        histogram.SetBinContent(x, 1, 1, average + deviations[x - 1]);
        histogram.SetBinError(x, 1, 1, 0.1);
    }
    TF3 polynomial("manual_integral", "[0]+x*x+0*y+0*z", -0.5, 0.5, -0.1, 0.2, -0.2, 0.3);
    polynomial.FixParameter(0, 1.2);
    FitConfig cfg;
    cfg.use_integral = true;
    const auto integral = calculate_fit_statistics(histogram, polynomial, cfg);
    CHECK(integral && integral->ndf == 3);
    CHECK(std::abs(integral->chi2 - 14.0) < 1e-8);
    cfg.use_integral = false;
    cfg.options = "RQS0I";
    const auto option_integral = calculate_fit_statistics(histogram, polynomial, cfg);
    CHECK(option_integral && std::abs(option_integral->chi2 - 14.0) < 1e-8);
    cfg.options = "RQS0";
    const auto centers = calculate_fit_statistics(histogram, polynomial, cfg);
    CHECK(centers && std::abs(centers->chi2 - 14.0) > 0.1);
}

} // namespace

int main()
{
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    check_manual_statistics();
    check_integral_statistics();

    const double hc2 = 0.197 * 0.197;
    const double R = 4.0;
    const double lambda = 0.8;

    TH3D hCF("hCF", "hCF", 40, -0.4, 0.4, 40, -0.4, 0.4, 40, -0.4, 0.4);
    for (int ix = 1; ix <= hCF.GetNbinsX(); ++ix) {
        const double qx = hCF.GetXaxis()->GetBinCenter(ix);
        for (int iy = 1; iy <= hCF.GetNbinsY(); ++iy) {
            const double qy = hCF.GetYaxis()->GetBinCenter(iy);
            for (int iz = 1; iz <= hCF.GetNbinsZ(); ++iz) {
                const double qz = hCF.GetZaxis()->GetBinCenter(iz);
                const double v =
                    1.0 + lambda * TMath::Exp(-(R * R * (qx * qx + qy * qy + qz * qz)) / hc2);
                hCF.SetBinContent(ix, iy, iz, v);
                hCF.SetBinError(ix, iy, iz, 0.01);
            }
        }
    }
    hCF.SetEntries(10000);

    Config cfg;
    cfg.input.type = "kt";

    const FitResult r = fit_cf_3d_with_retry(&hCF, cfg, 0, 0, 0);

    CHECK(r.ok);
    CHECK(r.status == 0);
    CHECK(r.cov_status == 1 || r.cov_status == 3);
    CHECK(!r.at_limit);
    CHECK(r.attempts == 1);
    CHECK(std::abs(r.r[0] - R) < 0.2);
    CHECK(std::abs(r.r[1] - R) < 0.2);
    CHECK(std::abs(r.r[2] - R) < 0.2);
    CHECK(std::abs(r.lambda - lambda) < 0.05);
    CHECK(std::abs(r.correlation(0, 0) - 1.0) < 1e-4);
    CHECK(std::abs(r.correlation(0, 1)) <= 1.0);

    // Check the production integration against ROOT on a fit with nonzero
    // residuals, and verify that the three fixed cross terms are not subtracted.
    for (int x = 1; x <= hCF.GetNbinsX(); ++x) {
        for (int y = 1; y <= hCF.GetNbinsY(); ++y) {
            for (int z = 1; z <= hCF.GetNbinsZ(); ++z) {
                const int bin = hCF.GetBin(x, y, z);
                hCF.SetBinContent(bin, hCF.GetBinContent(bin) + ((x + y + z) % 2 ? 0.005 : -0.005));
            }
        }
    }
    std::unique_ptr<TF3> fitted(create_cf_3d_fit(cfg, 0, 0, 0));
    const auto noisy = fit_cf_3d(&hCF, fitted.get(), cfg.fit);
    CHECK(noisy.ok);
    CHECK(noisy.ndf == 20 * 20 * 20 - 4);
    CHECK(std::abs(noisy.chi2 - fitted->GetChisquare()) < 1e-6);
    CHECK(noisy.ndf == fitted->GetNDF());
    CHECK(std::abs(noisy.p_value - fitted->GetProb()) < 1e-10);

    TH3D empty("empty", "", 4, -0.2, 0.2, 4, -0.2, 0.2, 4, -0.2, 0.2);
    const FitResult missing = fit_cf_3d_with_retry(&empty, cfg, 0, 0, 0);
    CHECK(!missing.ok);
    CHECK(missing.attempts == 0);

    std::cout << "All tests passed\n";
    return 0;
}
