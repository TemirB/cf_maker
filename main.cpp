#include <iostream>
#include <cmath>

#include <TFile.h>
#include <TH3.h>
#include <TF3.h>
#include <TDirectory.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>

#include "helpers.h"
#include "fit.h"
#include "draw.h"
#include "cf3d.h"
#include "plots.h"
#include "ratios.h"
#include "proj1d.h"
#include "proj2d.h"

int main(int argc, char** argv) {
    TH1::AddDirectory(kFALSE);

    // === Аргументы запуска === 
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
        << " <input.root> <output_dir>"
        << std::endl;
        return 1;
    }

    std::string resDir = "results/";
    std::string inputFile = argv[1];
    std::string outDir    = resDir + argv[2]; EnsureDir(outDir);

    TFile* input = TFile::Open(inputFile.c_str(), "READ");

    // ==== PHASE 1 ===
    // === CF + FIT ===
    std::string cf3dFileName = outDir + "/cf3d.root";
    TFile* fCF = new TFile(cf3dFileName.c_str(), "RECREATE");

    FitResult fitRes[chargeSize][centralitySize][ktSize][rapiditySize]{};

    for (int chIdx=0; chIdx<chargeSize; chIdx++) {
        create_and_fit_3d(chIdx, input, fCF, fitRes);
    }

    fCF->Write();
    fCF->Close();
    delete fCF;

    // ====== PHASE 2 ======
    // === Derived plots ===

    // kT dependence 
    {
        std::string fname = outDir + "/kt_dependence.root";
        TFile fout(fname.c_str(), "RECREATE");
        do_kt_diff(&fout, fitRes);
        fout.Write();
        fout.Close();
    }

    // rapidity dependence
    {
        std::string fname = outDir + "/rapidity_dependence.root";
        TFile fout(fname.c_str(), "RECREATE");
        do_rapidity_diff(&fout, fitRes);
        fout.Write();
        fout.Close();
    }
    
    // --- CF ratios (если нужно) ---
    /*
    {
        std::string fname = outDir + "/cf_ratios.root";
        TFile fout(fname.c_str(), "RECREATE");
        do_CF_ratios(&fout, fitRes);
        fout.Write();
        fout.Close();
    }
    */

    // --- 1D projections ---
    {
        std::string fname = outDir + "/proj1d.root";
        TFile fout(fname.c_str(), "RECREATE");
        do_1d_proj(input, &fout, fitRes);
        fout.Write();
        fout.Close();
    }

    // --- 2D projections ---
    /*
    {
        std::string fname = outDir + "/proj2d.root";
        TFile fout(fname.c_str(), "RECREATE");
        do_2d_proj(input, &fout);
        fout.Write();
        fout.Close();
    }
    */

    input->Close();
    delete input;

    std::cout << "All outputs written to: " << outDir << std::endl;
    return 0;
}


