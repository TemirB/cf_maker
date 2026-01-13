#include <iostream>
#include <TFile.h>
#include <TSystem.h>
#include <TH1.h>

#include "helpers.h"
#include "fit/cache.h"
#include "fit/json.h"
#include "plots.h"
#include "context.h"

int main(int argc, char** argv) {
    TH1::AddDirectory(kFALSE);

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <input.root> <output_dir>\n";
        return 1;
    }

    auto ctx = BuildContext(argv);

    TFile* input = new TFile(ctx.inputFile.c_str(), "READ");
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

        for (int ch=0; ch<chargeSize; ch++) BuildAndFit3DCorrelationFunctions(ch, input, &fCF, fitRes);

        fCF.Write();
        WriteFitJson(ctx.inputFile, ctx.outDir, fitRes);
    }
    
    {
        TFile* f = new TFile("kt.root", "RECREATE");
        MakeKtDependence(f, fitRes);

        f->Write();
        f->Close();
    }
    {
        TFile* f = new TFile("rapidity.root", "RECREATE");
        MakeRapidityDependence(f, fitRes);

        f->Write();
        f->Close();
    }
    {
        TFile* f = new TFile("1d.root", "RECREATE");
        MakeLCMS1DProjections(input, f, fitRes);

        f->Write();
        f->Close();
    }
    {
        TFile* f = new TFile("2d.root", "RECREATE");
        MakeLCMS2DProjections(input, f);

        f->Write();
        f->Close();
    }
    // {
    //     TFile* f1 = new TFile("ratio_projs.root", "RECREATE");
    //     TFile* f2 = new TFile("proj_ratios.root", "RECREATE");
    //     TFile* fCF3D = new TFile(ctx.cf3dFile.c_str(), "READ");
    //     do_CF_ratios(fCF3D, f1, f2);

    //     f1->Write();
    //     f1->Close();
    //     f2->Write();
    //     f2->Close();
    // }

    std::cout << "All outputs written to " << ctx.outDir << "\n";
    return 0;
}
