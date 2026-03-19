#include "plots.h"

#include <TFile.h>
#include <TF3.h>

#include "fit/types.h"
#include "fit/fit.h"
#include "helpers.h"

void BuildAndFit3DCorrelationFunctions(
    TFile* inputFile, TFile* outFile, FitGrid& fitRes, Bin bin
) {
    for (int ch = 0; ch < Charge::kCount; ch++)
    for (int centr = 0; centr < Centrality::kCount; centr++)
    for (int b = 0; b < bin.count; b++) {
        TF3* fit3d = CreateCF3DFit(ch, centr, b);

        TH3D* h_A = getNum(inputFile, ch, centr, b);
        TH3D* h_A_wei = getNumWei(inputFile, ch, centr, b);
        if (!h_A || !h_A_wei) {
            delete h_A;
            delete h_A_wei;
            continue;
        }
        h_A->SetDirectory(nullptr);
        h_A_wei->SetDirectory(nullptr);

        TH3D* h_CF = (TH3D*)h_A_wei->Clone("h_CF");
        h_CF->Reset("ICES");        
        // I	Integral — сбросить интеграл (общее число заполнений)
        // C	Contents — обнулить содержимое бинов (counts)
        // E	Errors — обнулить ошибки (Sumw2)
        // S	Statistics — сбросить статистику (mean, RMS, entries и т.д.)

        h_CF->Divide(h_A_wei, h_A, 1.0, 1.0, "B");

        delete h_A;
        delete h_A_wei;

        FitResult r = FitCF3D(h_CF, fit3d);
        fitRes[ch][centr][b] = r;

        std::string cfName = getCFName(ch, centr, bin.type, bin.names[b]); 

        outFile->cd();
        h_CF->SetDirectory(outFile);
        h_CF->Write(cfName.c_str(), TObject::kOverwrite);
        h_CF->SetDirectory(nullptr);
    }
};
