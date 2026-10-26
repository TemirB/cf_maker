#include "analysis/cf3d.h"

#include <algorithm>
#include <cctype>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <Math/MinimizerOptions.h>
#include <TFile.h>

#include "core/binning.h"
#include "core/correlation.h"
#include "core/log.h"
#include "core/parallel.h"
#include "fit/model.h"
#include "io/input.h"
#include "io/output.h"

namespace
{

struct Cf3dTask
{
    int ch;
    int centr;
    int b;
};

struct Cf3dResult
{
    FitResult fit;
    std::unique_ptr<TH3D> cf;
};

std::size_t effective_fit_threads(const Config& cfg)
{
    const auto is_legacy_minuit = [](std::string minimizer) {
        std::transform(
            minimizer.begin(), minimizer.end(), minimizer.begin(),
            [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
        return minimizer == "minuit";
    };
    const bool configured_legacy = is_legacy_minuit(cfg.fit.minimizer);
    // Direct callers may have selected a ROOT default different from cfg.
    const bool root_legacy = is_legacy_minuit(ROOT::Math::MinimizerOptions::DefaultMinimizerType());
    // ROOT's M (IMPROVE) option selects the legacy TMinuit implementation even
    // when the configured default minimizer is Minuit2.
    const bool improve =
        std::any_of(cfg.fit.options.begin(), cfg.fit.options.end(),
                    [](unsigned char character) { return std::toupper(character) == 'M'; });
    if (configured_legacy || root_legacy || improve) {
        if (cfg.threads == 0 || cfg.threads > 1) {
            logging::warn("cf3d: legacy Minuit uses shared state; forcing threads=1 "
                          "(requested " +
                          (cfg.threads == 0 ? std::string("auto") : std::to_string(cfg.threads)) +
                          ", " +
                          (improve ? std::string("fit option M")
                                   : (configured_legacy ? cfg.fit.minimizer
                                                        : std::string("ROOT default Minuit"))) +
                          ")");
        }
        return 1;
    }
    return static_cast<std::size_t>(cfg.threads);
}

} // namespace

void build_and_fit_3d_correlation_functions(Config& cfg, TFile* outFile)
{
    std::vector<Cf3dTask> tasks;
    for (const int ch : cfg.selection.charges) {
        for (const int centr : cfg.selection.centralities) {
            for (int b = 0; b < cfg.binning.count; b++) {
                tasks.push_back({ch, centr, b});
            }
        }
    }

    std::vector<Cf3dResult> results(tasks.size());
    const std::string inputPath = cfg.input.file;
    const std::size_t fit_threads = effective_fit_threads(cfg);

    logging::info("cf3d: " + std::to_string(tasks.size()) + " fit tasks, threads = " +
                  (fit_threads == 0 ? "auto" : std::to_string(fit_threads)));

    ParallelFor(tasks.size(), fit_threads, [&](std::size_t idx) {
        const Cf3dTask& task = tasks[idx];
        Cf3dResult& result = results[idx];

        thread_local std::unique_ptr<TFile> tFile;
        if (!tFile) {
            tFile = std::make_unique<TFile>(inputPath.c_str(), "READ");
            if (!tFile || tFile->IsZombie()) {
                tFile.reset();
                throw std::runtime_error("cannot open input file: " + inputPath);
            }
            logging::debug("cf3d worker: opened input file (thread-local)");
        }
        logging::debug("cf3d: fitting ch=" + std::to_string(task.ch) +
                       " centr=" + std::to_string(task.centr) + " b=" + std::to_string(task.b));
        auto [den_raw, num_raw] = get_hists(tFile.get(), task.ch, task.centr, task.b);
        std::unique_ptr<TH3D> den(den_raw), num(num_raw);
        if (!den || !num) {
            return;
        }
        den->SetDirectory(nullptr);
        num->SetDirectory(nullptr);

        auto cf_hist = std::unique_ptr<TH3D>(static_cast<TH3D*>(num->Clone("cf_hist")));
        fill_correlation(*cf_hist, *num, *den, correlation_statistics(*tFile));

        result.fit = fit_cf_3d_with_retry(cf_hist.get(), cfg, task.ch, task.centr, task.b);
        result.cf = std::move(cf_hist);
    });

    int nOk = 0;
    int n_usable = 0;
    int nRetried = 0;
    int nAtLimit = 0;
    int nMissing = 0;
    for (std::size_t i = 0; i < tasks.size(); ++i) {
        const Cf3dTask& task = tasks[i];
        Cf3dResult& result = results[i];
        cfg.fit_results[task.ch][task.centr][task.b] = result.fit;
        if (!result.cf) {
            nMissing++;
            continue;
        }

        if (result.fit.ok) {
            nOk++;
        }
        if (is_usable_fit(result.fit)) {
            n_usable++;
        }
        if (result.fit.attempts > 1) {
            nRetried++;
        }
        if (result.fit.at_limit) {
            nAtLimit++;
        }

        const std::string cfName =
            get_cf_name(task.ch, task.centr, cfg.input.type, cfg.binning.names[task.b]);

        write_output_object(*outFile, *result.cf, cfName.c_str(), TObject::kOverwrite);
    }

    logging::info("cf3d: fits ok=" + std::to_string(nOk) + "/" + std::to_string(tasks.size()) +
                  ", usable=" + std::to_string(n_usable) + "/" + std::to_string(tasks.size()) +
                  ", retried=" + std::to_string(nRetried) + ", atLimit=" +
                  std::to_string(nAtLimit) + ", missing=" + std::to_string(nMissing));
    if (nAtLimit > 0) {
        logging::warn("cf3d: " + std::to_string(nAtLimit) +
                      " fits have parameters at limits — check fit quality");
    }
}
