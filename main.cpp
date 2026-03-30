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

    Context ctx = BuildContext(argv);

    TFile* input = new TFile(ctx.inputFile.c_str(), "READ");
    if (!input || input->IsZombie()) {
        std::cerr << "ERROR: cannot open " << ctx.inputFile << "\n";
        return 1;
    }

    

    // Cf 3d
    {
        std::string name = ctx.outDir + "/cf3d.root";
        TFile f(name.c_str(), "RECREATE");

        BuildAndFit3DCorrelationFunctions(ctx, input, &f);

        f.Write();
    }
    // ===== Analysis =====
    // dependency
    {
        std::string name = Form("%s/%s.root", ctx.outDir.data(), ctx.bining.type);
        TFile f(name.c_str(), "RECREATE");

        MakeDependency(ctx, &f);

        f.Close();
    }

    // 1d_proj
    {
        std::string name = ctx.outDir + "/1d.root";
        TFile* f = new TFile(name.c_str(), "RECREATE");

        MakeLCMS1DProjections(ctx, input, f);

        f->Write();
        f->Close();
    }
    
    // 2d_proj
    {
        std::string name = ctx.outDir + "/2d.root";
        TFile* f = new TFile(name.c_str(), "RECREATE");

        MakeLCMS2DProjections(ctx, input, f);

        f->Write();
        f->Close();
    }
    {
        std::string name1 = ctx.outDir + "/ratio_projs.root";
        std::string name2 = ctx.outDir + "/proj_ratios.root";
        std::string cf3d = ctx.outDir + "/cf3d.root";

        TFile* f1 = new TFile(name1.c_str(), "RECREATE");
        TFile* f2 = new TFile(name2.c_str(), "RECREATE");
        TFile* fCF3D = new TFile(cf3d.c_str(), "READ");
        do_CF_ratios(ctx, fCF3D, f1, f2);

        f1->Write();
        f1->Close();
        f2->Write();
        f2->Close();
    }

    std::cout << "All outputs written to " << ctx.outDir << "\n";
    return 0;
}
