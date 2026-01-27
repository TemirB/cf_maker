#include "plots.h"

#include <TFile.h>
#include <TF3.h>

#include "fit/types.h"
#include "fit/fit.h"
#include "helpers.h"

void BuildAndFit3DCorrelationFunctions(
    TFile* inputFile,
    TFile* outFile,
    FitGrid& fitRes
) {
    TF3* fit3d = CreateCF3DFit();
    TH3D* h_CF_work = nullptr;

    for (int chIdx = 0; chIdx < chargeSize; chIdx++)
    for (int centIdx = 0; centIdx < centralitySize; centIdx++)
    for (int yIdx = 0; yIdx < rapiditySize; yIdx++) {
        TH3D* h_A = getNum(inputFile, chIdx, centIdx, yIdx);
        TH3D* h_A_wei = getNumWei(inputFile, chIdx, centIdx, yIdx);
        if (!h_A || !h_A_wei) {
            delete h_A;
            delete h_A_wei;
            continue;
        }
        h_A->SetDirectory(nullptr);
        h_A_wei->SetDirectory(nullptr);

        if (!h_CF_work) {
            h_CF_work = (TH3D*)h_A_wei->Clone("h_CF_work");
            h_CF_work->SetDirectory(nullptr); // не привязывать никуда
            if (!h_CF_work->GetSumw2N()) h_CF_work->Sumw2();
        } else if (!sameBinning(h_CF_work, h_A_wei)) {
            delete h_CF_work;
            h_CF_work = (TH3D*)h_A_wei->Clone("h_CF_work");
            h_CF_work->SetDirectory(nullptr);
            if (!h_CF_work->GetSumw2N()) h_CF_work->Sumw2();
        }

        h_CF_work->Reset("ICES");
        // I	Integral — сбросить интеграл (общее число заполнений)
        // C	Contents — обнулить содержимое бинов (counts)
        // E	Errors — обнулить ошибки (Sumw2)
        // S	Statistics — сбросить статистику (mean, RMS, entries и т.д.)

        h_CF_work->Divide(h_A_wei, h_A);

        delete h_A;
        delete h_A_wei;

        FitResult r = FitCF3D(h_CF_work, fit3d);
        fitRes[chIdx][centIdx][yIdx] = r;

        std::string cfName = getCFName(chIdx, centIdx, yIdx); 
        outFile->cd();
        h_CF_work->SetDirectory(outFile);
        h_CF_work->Write(cfName.c_str(), TObject::kOverwrite);
        h_CF_work->SetDirectory(nullptr);
    }

    delete fit3d;
    delete h_CF_work;
};
