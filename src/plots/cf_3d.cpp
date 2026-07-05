#include "plots.h"

#include <TFile.h>
#include <TF3.h>

#include "fit/types.h"
#include "fit/fit.h"
#include "helpers.h"
#include "config.h"

void BuildAndFit3DCorrelationFunctions(
    Config& cfg, TFile* inputFile, TFile* outFile
) { 
    for (int ch = 0; ch < Charge::kCount; ch++)
    for (int centr = 0; centr < Centrality::kCount; centr++)
    for (int b = 0; b < cfg.bining.count; b++) {
        TF3* fit3d = CreateCF3DFit(cfg, ch, centr, b);

        auto [A, Awei] = getHists(inputFile, ch, centr, b);
        if (!A || !Awei) {
            delete A;
            delete Awei;
            continue;
        }
        A->SetDirectory(nullptr);
        Awei->SetDirectory(nullptr);

        TH3D* h_CF = (TH3D*)Awei->Clone("h_CF");
        h_CF->Reset("ICES");        
        // I	Integral — сбросить интеграл (общее число заполнений)
        // C	Contents — обнулить содержимое бинов (counts)
        // E	Errors — обнулить ошибки (Sumw2)
        // S	Statistics — сбросить статистику (mean, RMS, entries и т.д.)

        h_CF->Divide(Awei, A, 1., 1., "B");

        delete A;
        delete Awei;

        FitResult r = FitCF3D(h_CF, fit3d);
        cfg.fitRes[ch][centr][b] = r;

        std::string cfName = getCFName(ch, centr, cfg.input.type, cfg.bining.names[b]); 

        outFile->cd();
        h_CF->SetDirectory(outFile);
        h_CF->Write(cfName.c_str(), TObject::kOverwrite);
        h_CF->SetDirectory(nullptr);
    }
};
