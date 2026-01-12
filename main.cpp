#include <iostream>
#include <memory>

#include <TFile.h>
#include <TH1.h>

#include "helpers.h"
#include "fit/cache.h"
#include "fit/json.h"
#include "plots.h"

#include "analysis/context.h"
#include "analysis/registry.h"
#include "analysis/runner.h"
#include "context.h"   // CLI context

int main(int argc, char** argv) {
    TH1::AddDirectory(kFALSE);

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <input.root> <output_dir>\n";
        return 1;
    }

    // CLI context (пути, имена файлов)
    auto cli = BuildContext(argv);

    std::unique_ptr<TFile> input(TFile::Open(cli.inputFile.c_str(), "READ"));
    if (!input || input->IsZombie()) {
        std::cerr << "ERROR: cannot open " << cli.inputFile << "\n";
        return 1;
    }

    // Fit cache
    static FitGrid* fitRes{};

    std::string reason;
    auto cache = LoadFitCache(cli.inputFile, cli.outDir, *fitRes, reason);

    if (cache == CacheStatus::Ok) {
        std::cout << "Using cached CF and fits\n";
    } else {
        std::cout << "Recomputing (" << reason << ")\n";

        TFile fCF(cli.cf3dFile.c_str(), "RECREATE");
        if (fCF.IsZombie()) {
            std::cerr << "ERROR: cannot create " << cli.cf3dFile << "\n";
            return 1;
        }

        for (int ch = 0; ch < chargeSize; ch++)
            BuildAndFit3DCorrelationFunctions(ch, input.get(), &fCF, *fitRes);

        fCF.Write();
        WriteFitJson(cli.inputFile, cli.outDir, *fitRes);
    }

    // ===== pipeline context =====
    RunContext rCtx;
    rCtx.input = input.get();
    rCtx.fit   = fitRes;

    // ===== run all stages =====
    RunAllStages(rCtx);

    std::cout << "All outputs written to " << cli.outDir << "\n";
    return 0;
}
