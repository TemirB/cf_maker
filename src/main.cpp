#include <exception>
#include <iostream>
#include <memory>
#include <string>

#include <Math/MinimizerOptions.h>
#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>

#include "analysis/pipeline.h"
#include "config/config.h"
#include "core/log.h"

// NOLINTNEXTLINE(bugprone-exception-escape): all exceptions are caught below
int main(int argc, char** argv) noexcept
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << "<path/to/config.json>\n";
        return 1;
    }

    try {
        TH1::AddDirectory(kFALSE);
        ROOT::EnableThreadSafety();

        Config cfg = load(argv[1]);

        const std::string logFile =
            cfg.logging.file.empty() ? "" : cfg.output.dir + "/" + cfg.logging.file;
        log::Init(log::parse_level(cfg.logging.level), logFile);

        ROOT::Math::MinimizerOptions::SetDefaultMinimizer(cfg.fit.minimizer.c_str());

        log::Info("cf_maker: config = " + std::string(argv[1]));
        log::Info("cf_maker: output dir = " + cfg.output.dir);
        log::Info("cf_maker: minimizer = " + cfg.fit.minimizer +
                  ", threads = " + std::to_string(cfg.threads) + " (0 = auto)");

        auto input = std::make_unique<TFile>(cfg.input.file.c_str(), "READ");
        if (!input || input->IsZombie()) {
            log::Error("cannot open input file: " + cfg.input.file);
            return 1;
        }
        log::Info("cf_maker: input file = " + cfg.input.file);

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

        std::cout << "All outputs written to " << cfg.output.dir << "\n";
        log::Info("cf_maker: all outputs written to " + cfg.output.dir);
        return 0;
    } catch (const std::exception& e) {
        try {
            log::Error(std::string("cf_maker: fatal error: ") + e.what());
        } catch (...) {
            std::cerr << "cf_maker: fatal error: " << e.what() << "\n";
        }
        return 1;
    } catch (...) {
        log::Error("cf_maker: unknown fatal error");
        return 1;
    }
}
