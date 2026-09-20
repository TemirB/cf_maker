#include "analysis/cf3d.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <TFile.h>

#include "core/binning.h"
#include "core/log.h"
#include "core/parallel.h"
#include "fit/model.h"
#include "io/input.h"

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

    log::Info("cf3d: " + std::to_string(tasks.size()) +
              " fit tasks, threads = " + (cfg.threads == 0 ? "auto" : std::to_string(cfg.threads)));

    ParallelFor(tasks.size(), static_cast<std::size_t>(cfg.threads), [&](std::size_t idx) {
        const Cf3dTask& task = tasks[idx];
        Cf3dResult& result = results[idx];

        thread_local std::unique_ptr<TFile> tFile;
        if (!tFile) {
            tFile = std::make_unique<TFile>(inputPath.c_str(), "READ");
            if (!tFile || tFile->IsZombie()) {
                tFile.reset();
                throw std::runtime_error("cannot open input file: " + inputPath);
            }
            log::Debug("cf3d worker: opened input file (thread-local)");
        }
        log::Debug("cf3d: fitting ch=" + std::to_string(task.ch) +
                   " centr=" + std::to_string(task.centr) + " b=" + std::to_string(task.b));
        auto [den, num] = get_hists(tFile.get(), task.ch, task.centr, task.b);
        if (!den || !num) {
            return;
        }
        den->SetDirectory(nullptr);
        num->SetDirectory(nullptr);

        auto cf_hist = std::unique_ptr<TH3D>(static_cast<TH3D*>(num->Clone("cf_hist")));
        cf_hist->Reset("ICES");
        cf_hist->Divide(num, den, 1., 1., "B");

        delete den;
        delete num;

        result.fit = fit_cf_3d_with_retry(cf_hist.get(), cfg, task.ch, task.centr, task.b);
        result.cf = std::move(cf_hist);
    });

    int nOk = 0;
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
        if (result.fit.attempts > 1) {
            nRetried++;
        }
        if (result.fit.at_limit) {
            nAtLimit++;
        }

        const std::string cfName =
            get_cf_name(task.ch, task.centr, cfg.input.type, cfg.binning.names[task.b]);

        outFile->cd();
        static_cast<void>(result.cf->Write(cfName.c_str(), TObject::kOverwrite));
    }

    log::Info("cf3d: fits ok=" + std::to_string(nOk) + "/" + std::to_string(tasks.size()) +
              ", retried=" + std::to_string(nRetried) + ", atLimit=" + std::to_string(nAtLimit) +
              ", missing=" + std::to_string(nMissing));
    if (nAtLimit > 0) {
        log::Warn("cf3d: " + std::to_string(nAtLimit) +
                  " fits have parameters at limits — check fit quality");
    }
}
