#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <TFile.h>
#include <TH3.h>
#include <TNamed.h>
#include <TRandom3.h>
#include <nlohmann/json.hpp>

#include "core/binning.h"

namespace
{
struct Options
{
    std::string type = "kt";
    std::array<double, 3> radii = {4.0, 5.0, 6.0};
    double lambda = 0.7;
    double q_max = 0.3;
    double counts = 10000;
    int bins = 40;
    unsigned int seed = 42;
    bool noise = false;
};

double number(const std::string& value)
{
    std::size_t end = 0;
    const double result = std::stod(value, &end);
    if (end != value.size() || !std::isfinite(result)) {
        throw std::runtime_error("Invalid number: " + value);
    }
    return result;
}

Options parse_options(int argc, char** argv)
{
    Options options;
    for (int i = 2; i < argc; ++i) {
        const std::string key = argv[i];
        if (key == "--noise") {
            options.noise = true;
            continue;
        }
        if (i + 1 == argc) {
            throw std::runtime_error("Missing value for " + key);
        }
        const std::string value = argv[++i];
        if (key == "--type") {
            options.type = value;
            continue;
        }
        const double parsed = number(value);
        if (key == "--r-out") {
            options.radii[0] = parsed;
        } else if (key == "--r-side") {
            options.radii[1] = parsed;
        } else if (key == "--r-long") {
            options.radii[2] = parsed;
        } else if (key == "--lambda") {
            options.lambda = parsed;
        } else if (key == "--q-max") {
            options.q_max = parsed;
        } else if (key == "--counts") {
            options.counts = parsed;
        } else if (key == "--bins" && parsed >= 8 && parsed <= 200 &&
                   std::floor(parsed) == parsed) {
            options.bins = static_cast<int>(parsed);
        } else if (key == "--seed" && parsed >= 1 && parsed <= 4294967295.0 &&
                   std::floor(parsed) == parsed) {
            options.seed = static_cast<unsigned int>(parsed);
        } else {
            throw std::runtime_error("Unknown option or invalid value: " + key + " " + value);
        }
    }
    if (options.type != "kt" && options.type != "rapidity") {
        throw std::runtime_error("--type must be kt or rapidity");
    }
    for (double radius : options.radii) {
        if (radius <= 1 || radius >= 9) {
            throw std::runtime_error("Radii must be between 1 and 9 fm (analysis validity range)");
        }
    }
    if (options.lambda <= 0.3 || options.lambda >= 1 || options.q_max < 0.2 || options.q_max > 1 ||
        options.counts < 1 || options.counts > 1e9) {
        throw std::runtime_error("Require 0.3 < lambda < 1, 0.2 <= q-max <= 1, 1 <= counts <= 1e9");
    }
    return options;
}

void generate(const std::filesystem::path& directory, const Options& options)
{
    using nlohmann::json;
    std::filesystem::create_directories(directory);
    const auto input_path = directory / "input.root";
    const auto config_path = directory / "config.json";
    if (std::filesystem::exists(input_path) || std::filesystem::exists(config_path)) {
        throw std::runtime_error(
            "input.root or config.json already exists; choose a new directory");
    }
    TH1::AddDirectory(false);
    TFile output(input_path.c_str(), "CREATE");
    if (output.IsZombie()) {
        throw std::runtime_error("Cannot create " + input_path.string());
    }
    TRandom3 random(options.seed);
    const int count = options.type == "kt" ? kt::kCount : rapidity::kCount;
    constexpr double kHc2 = 0.197 * 0.197;
    for (int ch = 0; ch < charge::kCount; ++ch) {
        for (int centr = 0; centr < centrality::kCount; ++centr) {
            for (int b = 0; b < count; ++b) {
                const std::string prefix = "bp_" + std::to_string(ch) + "_" + std::to_string(centr);
                const std::string suffix = "_" + std::to_string(b);
                TH3D den((prefix + "_num" + suffix).c_str(), "Unweighted pairs", options.bins,
                         -options.q_max, options.q_max, options.bins, -options.q_max, options.q_max,
                         options.bins, -options.q_max, options.q_max);
                TH3D num((prefix + "_num_wei" + suffix).c_str(), "Correlated pairs", options.bins,
                         -options.q_max, options.q_max, options.bins, -options.q_max, options.q_max,
                         options.bins, -options.q_max, options.q_max);
                den.Sumw2();
                num.Sumw2();
                for (int x = 1; x <= options.bins; ++x) {
                    for (int y = 1; y <= options.bins; ++y) {
                        for (int z = 1; z <= options.bins; ++z) {
                            const std::array<double, 3> q = {den.GetXaxis()->GetBinCenter(x),
                                                             den.GetYaxis()->GetBinCenter(y),
                                                             den.GetZaxis()->GetBinCenter(z)};
                            double exponent = 0;
                            for (std::size_t axis = 0; axis < q.size(); ++axis) {
                                exponent += std::pow(options.radii[axis] * q[axis], 2) / kHc2;
                            }
                            const double cf = 1 + options.lambda * std::exp(-exponent);
                            // Synthetic fixed reference; optional Poisson noise in the numerator.
                            const double value = options.noise
                                                     ? random.PoissonD(options.counts * cf)
                                                     : options.counts * cf;
                            den.SetBinContent(x, y, z, options.counts);
                            den.SetBinError(x, y, z, 0);
                            num.SetBinContent(x, y, z, value);
                            num.SetBinError(x, y, z, std::sqrt(value));
                        }
                    }
                }
                if (den.Write() <= 0 || num.Write() <= 0) {
                    throw std::runtime_error("Cannot write histograms");
                }
            }
        }
    }
    const json truth = {{"r_out", options.radii[0]},
                        {"r_side", options.radii[1]},
                        {"r_long", options.radii[2]},
                        {"lambda", options.lambda},
                        {"type", options.type},
                        {"bins", options.bins},
                        {"q_max", options.q_max},
                        {"counts", options.counts},
                        {"noise", options.noise},
                        {"seed", options.seed},
                        {"hbar_c", 0.197}};
    TNamed metadata("generator_truth", truth.dump().c_str());
    if (metadata.Write() <= 0) {
        throw std::runtime_error("Cannot write generator metadata");
    }
    output.Close();
    const json config = {
        {"machine", {{"base_input", directory.string()}, {"base_output", directory.string()}}},
        {"vars",
         {{"input", {{"file", "input.root"}, {"type", options.type}}},
          {"output", {{"dir", "results"}}}}},
        {"fit", {{"use_default_ip", true}}},
        {"images", {{"need", false}}},
        {"threads", 1}};
    std::ofstream config_file(config_path);
    config_file << config.dump(4) << '\n';
    config_file.close();
    if (!config_file) {
        throw std::runtime_error("Cannot write " + config_path.string());
    }
    std::cout << "Generated " << input_path << "\nRun: ./build/main " << config_path << '\n';
}
} // namespace

int main(int argc, char** argv)
{
    if (argc < 2 || std::string(argv[1]) == "--help") {
        std::cout << "Usage: generate_gaussian OUTPUT_DIR [options]\n"
                     "  --type kt|rapidity (kt)\n"
                     "  --r-out 4 --r-side 5 --r-long 6  radii in fm\n"
                     "  --lambda 0.7 --bins 40 --q-max 0.3 (GeV/c)\n"
                     "  --counts 10000  reference count per 3D bin\n"
                     "  --noise  enable Poisson noise; --seed 42 (nonzero)\n";
        return argc < 2 ? 1 : 0;
    }
    try {
        const auto options = parse_options(argc, argv);
        generate(std::filesystem::absolute(argv[1]), options);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "generate_gaussian: " << error.what() << '\n';
        return 1;
    }
}
