#include <iostream>
#include <TFile.h>
#include <TSystem.h>
#include <TH1.h>

#include "helpers.h"
#include "fit/cache.h"
#include "fit/json.h"
#include "cf3d.h"
#include "plots.h"
#include "proj1d.h"
#include "proj2d.h"
#include "ratios.h"
#include "context.h"

int main(int argc, char** argv) {
    TH1::AddDirectory(kFALSE);

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <input.root> <output_dir>\n";
        return 1;
    }

    auto ctx = BuildContext(argv);

    std::unique_ptr<TFile> input(TFile::Open(ctx.inputFile.c_str(), "READ"));
    if (!input || input->IsZombie()) {
        std::cerr << "ERROR: cannot open " << ctx.inputFile << "\n";
        return 1;
    }

    FitResult fitRes[chargeSize][centralitySize][ktSize][rapiditySize]{};

    std::string reason;
    auto cache = LoadFitCache(ctx.inputFile, ctx.outDir, fitRes, reason);

    if (cache == CacheStatus::Ok) {
        std::cout << "Using cached CF and fits\n";
    } else {
        std::cout << "Recomputing (" << reason << ")\n";

        TFile fCF(ctx.cf3dFile.c_str(), "RECREATE");
        if (fCF.IsZombie()) {
            std::cerr << "ERROR: cannot create " << ctx.cf3dFile << "\n";
            return 1;
        }

        for (int ch=0; ch<chargeSize; ch++)
            create_and_fit_3d(ch, input.get(), &fCF, fitRes);

        fCF.Write();
        WriteFitJson(ctx.inputFile, ctx.outDir, fitRes);
    }

    std::cout << "All outputs written to " << ctx.outDir << "\n";
    return 0;
}
