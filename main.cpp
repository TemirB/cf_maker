#include <iostream>
#include <TFile.h>
#include <TSystem.h>
#include <TH1.h>
#include <TTree.h>

#include "plots.h"
#include "config.h"

int main(int argc, char** argv) {
    TH1::AddDirectory(kFALSE);

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0]
                  << "<path/to/config>\n";
        return 1;
    } 
    std::string path = argv[1];

    Config cfg = Load(path);

    TFile* input = new TFile(cfg.input.file.c_str(), "READ");
    if (!input || input->IsZombie()) {
        std::cerr << "ERROR: cannot open " << cfg.input.file << "\n";
        return 1;
    }

    // Cf 3d
    {
        std::string name = cfg.output.dir + "/cf3d.root";
        TFile f(name.c_str(), "RECREATE");

        BuildAndFit3DCorrelationFunctions(cfg, input, &f);

        f.Write();
    }
    // ===== Analysis =====
    // dependency
    {
        std::string name = Form(
            "%s/%s.root", 
            cfg.output.dir.c_str(), 
            cfg.input.type.c_str()
        );

        std::string cf3dName = cfg.output.dir + "/cf3d.root";
        TFile f(name.c_str(), "RECREATE");
        TFile cf3d(cf3dName.c_str(), "READ");

        MakeDependency(cfg, &cf3d, &f);

        f.Close();
    }

    // 1d_proj
    {
        std::string name = cfg.output.dir + "/1d.root";
        TFile* f = new TFile(name.c_str(), "RECREATE");

        MakeLCMS1DProjections(cfg, input, f);

        f->Write();
        f->Close();
    }
    
    // 2d_proj
    {
        std::string name = cfg.output.dir + "/2d.root";
        TFile* f = new TFile(name.c_str(), "RECREATE");

        MakeLCMS2DProjections(cfg, input, f);

        f->Write();
        f->Close();
    }

    // ratios
    {
        // need to check fake data in proj_ratio
        std::string name1 = cfg.output.dir + "/ratio_projs.root";
        std::string name2 = cfg.output.dir + "/proj_ratios.root";
        std::string cf3d = cfg.output.dir + "/cf3d.root";

        TFile* f1 = new TFile(name1.c_str(), "RECREATE");
        TFile* f2 = new TFile(name2.c_str(), "RECREATE");
        TFile* fCF3D = new TFile(cf3d.c_str(), "READ");
        do_CF_ratios(cfg, fCF3D, f1, f2);

        f1->Write();
        f1->Close();
        f2->Write();
        f2->Close();
    }

    std::cout << "All outputs written to " << cfg.output.dir << "\n";
    return 0;
}
