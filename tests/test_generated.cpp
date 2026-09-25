#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

#include <TFile.h>
#include <TH3.h>
#include <TList.h>
#include <TNamed.h>
#include <nlohmann/json.hpp>

#include "fit/model.h"
#include "io/input.h"

int main(int argc, char** argv)
{
    if (argc != 2) {
        return 1;
    }
    TH1::AddDirectory(false);
    TFile input(argv[1], "READ");
    auto* metadata = dynamic_cast<TNamed*>(input.Get("generator_truth"));
    if (input.IsZombie() || !metadata || input.GetListOfKeys()->GetSize() != 65) {
        throw std::runtime_error("Missing generated histograms or metadata");
    }
    const auto truth = nlohmann::json::parse(metadata->GetTitle());
    Config cfg;
    cfg.input.type = "kt";
    cfg.fit.use_default_ip = true;
    for (int ch = 0; ch < 2; ++ch) {
        for (int centr = 0; centr < 4; ++centr) {
            for (int bin = 0; bin < 4; ++bin) {
                const auto [den_ptr, num_ptr] = get_hists(&input, ch, centr, bin);
                std::unique_ptr<TH3D> den(den_ptr), num(num_ptr);
                if (!den || !num) {
                    throw std::runtime_error("Missing histogram pair");
                }
                TH3D cf(*num);
                cf.Divide(num.get(), den.get(), 1., 1., "B");
                for (int x = 1; x <= cf.GetNbinsX(); ++x) {
                    const int y = cf.GetNbinsY() / 2;
                    const int z = cf.GetNbinsZ() / 2;
                    const double qx = cf.GetXaxis()->GetBinCenter(x);
                    const double qy = cf.GetYaxis()->GetBinCenter(y);
                    const double qz = cf.GetZaxis()->GetBinCenter(z);
                    const double expected =
                        1 + 0.7 * std::exp(-(16 * qx * qx + 25 * qy * qy + 36 * qz * qz) /
                                           (0.197 * 0.197));
                    if (std::abs(cf.GetBinContent(x, y, z) - expected) > 1e-12) {
                        throw std::runtime_error("Generated CF does not match Gaussian truth");
                    }
                }
                if (ch == 0 && centr == 0 && bin == 0) {
                    const auto fit = fit_cf_3d_with_retry(&cf, cfg, ch, centr, bin);
                    if (!fit.ok || std::abs(fit.r[0] - truth.at("r_out").get<double>()) > 0.02 ||
                        std::abs(fit.r[1] - truth.at("r_side").get<double>()) > 0.02 ||
                        std::abs(fit.r[2] - truth.at("r_long").get<double>()) > 0.02 ||
                        std::abs(fit.lambda - truth.at("lambda").get<double>()) > 0.002) {
                        throw std::runtime_error("Fit failed to recover Gaussian parameters");
                    }
                }
            }
        }
    }
    std::cout << "Generated data and recovered fit parameters are correct\n";
}
