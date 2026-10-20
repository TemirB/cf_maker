#include <array>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include <Math/MinimizerOptions.h>
#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>

#include "analysis/pipeline.h"
#include "config/config.h"
#include "core/log.h"
#include "fit/types.h"
#include "io/run_manifest.h"

namespace
{
void log_unusable_fit(const Config& cfg, int ch, int centr, int b)
{
    const FitResult& result = cfg.fit_results[ch][centr][b];
    std::ostringstream message;
    message << std::setprecision(10) << "fit (ch=" << ch << ", centr=" << centr << ", b=" << b
            << "): unusable:";
    if (result.attempts == 0) {
        message << " missing fit;";
    } else {
        if (!result.ok) {
            message << " minimizer result rejected;";
        }
        if (!result.is_finite()) {
            message << " nonfinite fit fields;";
        }
        if (result.ndf <= 0) {
            message << " ndf<=0;";
        }
        if (result.at_limit) {
            message << " parameter at limit;";
        }
    }
    message << " status=" << result.status << ", covStatus=" << result.cov_status
            << ", attempts=" << result.attempts << ", chi2=" << result.chi2
            << ", ndf=" << result.ndf << ", R=[" << result.r[0] << ", " << result.r[1] << ", "
            << result.r[2] << "], cross=[" << result.r[3] << ", " << result.r[4] << ", "
            << result.r[5] << "], lambda=" << result.lambda;

    if (result.at_limit) {
        constexpr std::array<const char*, 7> kParameterNames = {
            "R_out_sq",   "R_side_sq",   "R_long_sq", "R_out_side",
            "R_out_long", "R_side_long", "lambda"};
        message << ", reached_limits={";
        const char* separator = "";
        for (std::size_t i = 0; i < kParameterNames.size(); ++i) {
            if (cfg.fit.freeze[i].has_value()) {
                continue;
            }
            const double value =
                i < 3 ? result.r[i] * result.r[i] : (i < 6 ? result.r[i] : result.lambda);
            const double lower =
                i < 3 ? cfg.fit.radius_sq_min : (i < 6 ? cfg.fit.cross_min : cfg.fit.lambda_min);
            const double upper =
                i < 3 ? cfg.fit.radius_sq_max : (i < 6 ? cfg.fit.cross_max : cfg.fit.lambda_max);
            if (std::isfinite(value) &&
                (value <= lower + kFitLimitTolerance || value >= upper - kFitLimitTolerance)) {
                message << separator << kParameterNames[i] << '=' << value << " in [" << lower
                        << ", " << upper << ']';
                separator = ", ";
            }
        }
        message << '}';
    }
    logging::warn(message.str());
}
} // namespace

// NOLINTNEXTLINE(bugprone-exception-escape): all exceptions are caught below
int main(int argc, char** argv) noexcept
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << "<path/to/config.json>\n";
        return 1;
    }

    try {
        gROOT->SetBatch(kTRUE);
        TH1::AddDirectory(kFALSE);
        ROOT::EnableThreadSafety();

        Config cfg = load(argv[1]);
        cfg.input.file = std::filesystem::absolute(cfg.input.file).lexically_normal().string();
        cfg.output.dir = std::filesystem::absolute(cfg.output.dir).lexically_normal().string();

        auto input = std::make_unique<TFile>(cfg.input.file.c_str(), "READ");
        if (!input || input->IsZombie()) {
            throw std::runtime_error("cannot open input file: " + cfg.input.file);
        }
        prepare_stage_dependencies(cfg);
        write_run_config(cfg);
        RunManifest run_manifest(cfg);

        const std::string logFile =
            cfg.logging.file.empty() ? "" : cfg.output.dir + "/" + cfg.logging.file;
        logging::init(logging::parse_level(cfg.logging.level), logFile);

        ROOT::Math::MinimizerOptions::SetDefaultMinimizer(cfg.fit.minimizer.c_str());

        logging::info("cf_maker: config = " + std::string(argv[1]));
        logging::info("cf_maker: output dir = " + cfg.output.dir);
        logging::info("cf_maker: minimizer = " + cfg.fit.minimizer +
                      ", threads = " + std::to_string(cfg.threads) + " (0 = auto)");

        logging::info("cf_maker: input file = " + cfg.input.file);

        if (cfg.stages.cf3d) {
            stage_cf3d(cfg);
        }
        if (cfg.stages.dependency) {
            stage_dependency(cfg);
        }
        if (cfg.stages.projections_1d) {
            stage_1d_projections(cfg, input.get());
        }
        if (cfg.stages.projections_2d) {
            stage_2d_projections(cfg, input.get());
        }
        if (cfg.stages.ratios) {
            stage_ratios(cfg);
        }

        if (cfg.stages.cf3d || cfg.stages.dependency || cfg.stages.projections_1d) {
            std::size_t unusable = 0;
            std::size_t requested = 0;
            for (const int ch : cfg.selection.charges) {
                for (const int centr : cfg.selection.centralities) {
                    for (int b = 0; b < cfg.binning.count; ++b) {
                        ++requested;
                        if (!is_usable_fit(cfg.fit_results[ch][centr][b])) {
                            ++unusable;
                            log_unusable_fit(cfg, ch, centr, b);
                        }
                    }
                }
            }
            if (unusable > 0) {
                logging::error("cf_maker: " + std::to_string(unusable) + "/" +
                               std::to_string(requested) +
                               " requested fits are missing, failed or "
                               "unusable; diagnostic outputs written to " +
                               cfg.output.dir);
                logging::finish();
                run_manifest.finish(cfg, 2);
                return 2;
            }
        }

        logging::info("cf_maker: requested stages completed");
        logging::finish();
        run_manifest.finish(cfg, 0);
        std::cout << "All outputs written to " << cfg.output.dir << "\n";
        logging::info("cf_maker: all outputs written to " + cfg.output.dir);
        return 0;
    } catch (const std::exception& e) {
        try {
            logging::error(std::string("cf_maker: fatal error: ") + e.what());
        } catch (...) {
            std::cerr << "cf_maker: fatal error: " << e.what() << "\n";
        }
        return 1;
    } catch (...) {
        logging::error("cf_maker: unknown fatal error");
        return 1;
    }
}
