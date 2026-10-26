#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include <Math/MinimizerOptions.h>
#include <TFile.h>
#include <TH3D.h>
#include <TMemFile.h>
#include <TNamed.h>
#include <TROOT.h>

#include "analysis/cf3d.h"
#include "core/log.h"

namespace
{

constexpr std::array<double, 3> kRadii = {4., 5., 6.};
constexpr double kHc2 = 0.197 * 0.197;

void write_input(const std::filesystem::path& path)
{
    TFile input(path.c_str(), "RECREATE");
    if (input.IsZombie()) {
        throw std::runtime_error("cannot create Minuit test input");
    }
    TNamed statistics("correlation_statistics", "fixed_reference");
    statistics.Write();
    for (int centr = 0; centr < 4; ++centr) {
        const auto den_name = "bp_0_" + std::to_string(centr) + "_num_0";
        const auto num_name = "bp_0_" + std::to_string(centr) + "_num_wei_0";
        TH3D denominator(den_name.c_str(), "", 24, -0.3, 0.3, 24, -0.3, 0.3, 24, -0.3, 0.3);
        TH3D numerator(num_name.c_str(), "", 24, -0.3, 0.3, 24, -0.3, 0.3, 24, -0.3, 0.3);
        for (int x = 1; x <= 24; ++x) {
            const double q_out = denominator.GetXaxis()->GetBinCenter(x);
            for (int y = 1; y <= 24; ++y) {
                const double q_side = denominator.GetYaxis()->GetBinCenter(y);
                for (int z = 1; z <= 24; ++z) {
                    const double q_long = denominator.GetZaxis()->GetBinCenter(z);
                    const double exponent = std::pow(kRadii[0] * q_out, 2) +
                                            std::pow(kRadii[1] * q_side, 2) +
                                            std::pow(kRadii[2] * q_long, 2);
                    const double cf = 1 + (0.5 + 0.1 * centr) * std::exp(-exponent / kHc2);
                    denominator.SetBinContent(x, y, z, 10000);
                    denominator.SetBinError(x, y, z, 0);
                    numerator.SetBinContent(x, y, z, 10000 * cf);
                    numerator.SetBinError(x, y, z, std::sqrt(10000 * cf));
                }
            }
        }
        denominator.Write();
        numerator.Write();
    }
}

void check_run(const std::filesystem::path& directory, const std::string& name,
               const std::string& configured_minimizer, const std::string& root_minimizer,
               const std::string& fit_options, int threads, bool expect_serial,
               int expected_status = 0)
{
    Config cfg;
    cfg.input.file = (directory / "input.root").string();
    cfg.input.type = "kt";
    cfg.fit.minimizer = configured_minimizer;
    cfg.fit.options = fit_options;
    cfg.fit.use_default_ip = true;
    cfg.fit.retry_with_defaults = false;
    cfg.selection.charges = {0};
    cfg.selection.centralities = {0, 1, 2, 3};
    cfg.binning.count = 1;
    cfg.binning.names = {"test"};
    cfg.threads = threads;
    cfg.fit_results = FitGrid(1, std::vector<std::vector<FitResult>>(4, std::vector<FitResult>(1)));

    const auto log_path = directory / (name + ".log");
    std::ofstream(log_path, std::ios::trunc).close();
    logging::init(logging::Level::Debug, log_path.string());
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(root_minimizer.c_str());
    TMemFile output((name + ".root").c_str(), "RECREATE");
    build_and_fit_3d_correlation_functions(cfg, &output);
    logging::init(logging::Level::Error, "");
    if (cfg.threads != threads) {
        throw std::runtime_error("effective worker limit changed requested configuration");
    }
    for (int centr = 0; centr < 4; ++centr) {
        const auto& result = cfg.fit_results[0][centr][0];
        if (result.status != expected_status || result.ok != (expected_status == 0) ||
            std::abs(result.lambda - (0.5 + 0.1 * centr)) > 0.002) {
            throw std::runtime_error(name + ": incorrect fit for centrality " +
                                     std::to_string(centr));
        }
        for (std::size_t axis = 0; axis < kRadii.size(); ++axis) {
            if (std::abs(result.r[axis] - kRadii[axis]) > 0.02) {
                throw std::runtime_error(name + ": incorrect fitted radius");
            }
        }
    }
    std::ifstream log(log_path);
    std::ostringstream contents;
    std::set<std::string> fit_threads;
    std::string line;
    int fits = 0;
    while (std::getline(log, line)) {
        contents << line << '\n';
        if (line.find("cf3d: fitting ch=") == std::string::npos) {
            continue;
        }
        ++fits;
        const auto start = line.find("[t=");
        const auto finish = line.find(']', start);
        if (start == std::string::npos || finish == std::string::npos) {
            throw std::runtime_error("cannot read logged fit worker id");
        }
        fit_threads.insert(line.substr(start, finish - start + 1));
    }
    if (fits != 4 || (expect_serial && fit_threads.size() != 1) ||
        (!expect_serial && fit_threads.size() < 2)) {
        throw std::runtime_error(name + ": incorrect fit worker concurrency");
    }
    const bool warned = contents.str().find("forcing threads=1") != std::string::npos;
    if (warned != expect_serial) {
        throw std::runtime_error(name + ": missing or unnecessary Minuit concurrency warning");
    }
}

} // namespace

int main()
{
    gROOT->SetBatch(true);
    TH1::AddDirectory(false);
    ROOT::EnableThreadSafety();
    const auto directory = std::filesystem::absolute("minuit_threads_test_output");
    std::filesystem::create_directories(directory);
    write_input(directory / "input.root");
    check_run(directory, "legacy", "mInUiT", "Minuit", "RQS0", 4, true);
    // IMPROVE reports status 4000 when an exact minimum has no better neighbour.
    // The fitted parameters and serial execution still must be correct; this
    // test must not reclassify ROOT's nonzero status as an accepted fit.
    check_run(directory, "improve", "Minuit2", "Minuit2", "RQS0m", 4, true, 4000);
    check_run(directory, "automatic", "Minuit", "Minuit", "RQS0", 0, true);
    check_run(directory, "root_default", "Minuit2", "Minuit", "RQS0", 4, true);
    check_run(directory, "minuit2", "Minuit2", "Minuit2", "RQS0", 4, false);
    std::filesystem::remove_all(directory);
    std::cout << "Minuit concurrency tests passed\n";
}
