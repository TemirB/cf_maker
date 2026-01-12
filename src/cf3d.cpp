#include "cf3d.h"

#include <TFile.h>

#include "fit.h"
#include "helpers.h"

void create_and_fit_3d(
    int chIdx,
    TFile* inputFile, TFile* outFile,
    FitResult (&fitRes)[chargeSize][centralitySize][ktSize][rapiditySize]
) {
    TDirectory* dir = outFile->GetDirectory("cf3d");
    if (!dir) dir = outFile->mkdir("cf3d");
    dir->cd();

    TF3* fit3d = CreateCF3DFit();
    TH3D* h_CF_work = nullptr;

    for (int centIdx = 0; centIdx < centralitySize; centIdx++) {
        for (int ktIdx = 0; ktIdx < ktSize; ktIdx++) {
            for (int yIdx = 0; yIdx < rapiditySize; yIdx++) {

                std::string mid = getPrefix(chIdx, centIdx, yIdx);

                std::string prefix = "bp_" + mid;
                std::string cfName = "h3d_CF_q_" + mid + std::to_string(ktIdx) + "_weighted";

                TH3D* h_A    = (TH3D*) inputFile->Get((prefix + std::to_string(ktIdx)).c_str()); 
                TH3D* h_A_wei= (TH3D*) inputFile->Get((prefix + "wei_" + std::to_string(ktIdx)).c_str());
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

                // Создаем новую 3D-гистограмму для CF = A_wei / A с тем же бинированием, что и у входных
                h_CF_work->Reset("ICES");
                h_CF_work->Divide(h_A_wei, h_A);

                delete h_A;
                delete h_A_wei;

                // if (emptyForFit(h_CF_work)) {
                //     h_CF_work->Write(cfName.c_str(), TObject::kOverwrite);
                //     continue;
                // }
                
                // if (h_CF_work->GetEntries() == 0) {
                //     std::cerr << "EMPTY CF for "
                //             << "ch=" << chIdx
                //             << " cent=" << centIdx
                //             << " kt=" << ktIdx
                //             << " y=" << yIdx
                //             << std::endl;
                //     h_CF_work->Write(cfName.c_str(), TObject::kOverwrite);
                //     continue;
                // }

                FitResult r = FitCF3D(h_CF_work, fit3d);
                fitRes[chIdx][centIdx][ktIdx][yIdx] = r;

                h_CF_work->Write(cfName.c_str(), TObject::kOverwrite);
            }
        }
    }

    delete fit3d;
    delete h_CF_work;
}