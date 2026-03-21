#include <iostream>
#include <TFile.h>
#include <TSystem.h>
#include <TH1.h>
#include <TTree.h>

#include "helpers.h"
#include "plots.h"
#include "context.h"

int main(int argc, char** argv) {
    TH1::AddDirectory(kFALSE);

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <input.root> <output_dir> <type = rapidity/kt>\n";
        return 1;
    }

    AnalysisType aType = getType(argv[3]);
    if (aType == AnalysisType::Unknown) {
        std::cerr << "Unknown type: " << argv[3] << std::endl;
        return 1;
    }

    auto ctx = BuildContext(argv);

    TFile* input = new TFile(ctx.inputFile.c_str(), "READ");
    if (!input || input->IsZombie()) {
        std::cerr << "ERROR: cannot open " << ctx.inputFile << "\n";
        return 1;
    }

    Bin bin;
    {
        bin.type = "kt";
        bin.count = Kt::kCount;
        bin.names = Kt::kNames;
        bin.values = Kt::kValues;

        if (aType == AnalysisType::Rapidity) {
            bin.type = "rapidity";
            bin.count = Rapidity::kCount;
            bin.names = Rapidity::kNames;
            bin.values = Rapidity::kValues;
        }
    }

    FitGrid fitRes(
        Charge::kCount,
        std::vector<std::vector<FitResult>>(
            Centrality::kCount,
            std::vector<FitResult>(bin.count)
        )
    );

    TFile fCF(ctx.cf3dFile.c_str(), "RECREATE");
    if (fCF.IsZombie()) {
        std::cerr << "ERROR: cannot create " << ctx.cf3dFile << "\n";
        return 1;
    }

    BuildAndFit3DCorrelationFunctions(input, &fCF, fitRes, bin);

    fCF.Write();

    // ===== Analysis =====
    // dependency
    {
        std::string name = ctx.outDir + "/rapidity.root";
        TFile f(name.c_str(), "RECREATE");

        MakeDependency(&f, fitRes, bin);

        f.Close();
    }

    // 1d_proj
    {
        std::string name = ctx.outDir + "/1d.root";
        TFile* f = new TFile(name.c_str(), "RECREATE");

        MakeLCMS1DProjections(input, f, fitRes, bin);

        f->Write();
        f->Close();
    }
    
    // 2d_proj
    {
        std::string name = ctx.outDir + "/2d.root";
        TFile* f = new TFile(name.c_str(), "RECREATE");

        MakeLCMS2DProjections(input, f, bin);

        f->Write();
        f->Close();
    }
    {
        std::string name1 = ctx.outDir + "/ratio_projs.root";
        std::string name2 = ctx.outDir + "/proj_ratios.root";

        TFile* f1 = new TFile(name1.c_str(), "RECREATE");
        TFile* f2 = new TFile(name2.c_str(), "RECREATE");
        TFile* fCF3D = new TFile(ctx.cf3dFile.c_str(), "READ");
        do_CF_ratios(fCF3D, f1, f2, bin);

        f1->Write();
        f1->Close();
        f2->Write();
        f2->Close();
    }

    std::cout << "All outputs written to " << ctx.outDir << "\n";
    return 0;
}
