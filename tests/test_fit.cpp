#include <cmath>
#include <cstdlib>
#include <iostream>

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

} // namespace

int main()
{
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

    const double hc2 = 0.197 * 0.197;
    const double R = 4.0;
    const double lambda = 0.8;

    TH3D hCF("hCF", "hCF", 40, -0.4, 0.4, 40, -0.4, 0.4, 40, -0.4, 0.4);
    for (int ix = 1; ix <= hCF.GetNbinsX(); ++ix) {
        const double qx = hCF.GetXaxis()->Getbin_center(ix);
        for (int iy = 1; iy <= hCF.GetNbinsY(); ++iy) {
            const double qy = hCF.GetYaxis()->Getbin_center(iy);
            for (int iz = 1; iz <= hCF.GetNbinsZ(); ++iz) {
                const double qz = hCF.GetZaxis()->Getbin_center(iz);
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
    CHECK(std::abs(r.r[0] - R) < 0.2);
    CHECK(std::abs(r.r[1] - R) < 0.2);
    CHECK(std::abs(r.r[2] - R) < 0.2);
    CHECK(std::abs(r.lambda - lambda) < 0.05);
    CHECK(std::abs(r.corr(0, 0) - 1.0) < 1e-4);
    CHECK(std::abs(r.corr(0, 1)) <= 1.0);

    std::cout << "All tests passed\n";
    return 0;
}
