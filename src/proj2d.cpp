#include <TFile.h>
#include <TH2.h>
#include <TDirectory.h>
#include <string>

#include "fit.h"
#include "helpers.h"

bool contains(const std::string str, const std::string sub) {
    return str.find(sub) != std::string::npos;
}

void prepare_for_proj_2d(TH3D*  h_A, TH3D* h_A_wei, const std::string axis) {
    Double_t lRange = -0.2;
    Double_t rRange = 0.2;
    if (contains(axis, axises[2])) {
        h_A_wei->GetZaxis()->SetRangeUser(lRange, rRange);
        h_A->GetZaxis()->SetRangeUser(lRange, rRange);
    } else if (contains(axis, axises[1])) {
        h_A_wei->GetYaxis()->SetRangeUser(lRange, rRange);
        h_A->GetYaxis()->SetRangeUser(lRange, rRange);
    } else {
        h_A_wei->GetXaxis()->SetRangeUser(lRange, rRange);
        h_A->GetXaxis()->SetRangeUser(lRange, rRange);
    }
}

void get_2d_proj(
    TDirectory* dir ,TH3D* h_A, TH3D* h_A_wei,
    const std::string& pr_axis,
    const std::string& LCMS1,
    const std::string& LCMS2,
    const std::string& ending
) {
    // Ресетаем оси на дефолт
    auto resetRanges = [&](TH3D* h){
        h->GetXaxis()->SetRange(0,0);
        h->GetYaxis()->SetRange(0,0);
        h->GetZaxis()->SetRange(0,0);
    };

    resetRanges(h_A);
    resetRanges(h_A_wei);

    // Подготваливаем оси
    prepare_for_proj_2d(h_A, h_A_wei, pr_axis);

    // Формируем проекции
    std::string name    = "CF_{q_" + LCMS1 + "_vs_" + LCMS2 + "}_" + ending;
    TH2D* h_A_proj      = (TH2D*) h_A->Project3D(pr_axis.c_str());
    TH2D* h_A_wei_proj  = (TH2D*) h_A_wei->Project3D(pr_axis.c_str());
    h_A_proj->SetDirectory(nullptr);
    h_A_wei_proj->SetDirectory(nullptr);

    // Создаем 2D гистограмму
    TH2D* h_C_proj = (TH2D*) h_A_wei_proj->Clone(name.c_str());
    h_C_proj->Reset();
    h_C_proj->Divide(h_A_wei_proj, h_A_proj);
    h_C_proj->GetYaxis()->SetRangeUser(-0.13, 0.13);
    h_C_proj->GetYaxis()->SetTitle(("q_{" + LCMS1 + "} [GeV/c]").c_str());
    h_C_proj->GetXaxis()->SetRangeUser(-0.13, 0.13);
    h_C_proj->GetXaxis()->SetTitle(("q_{" + LCMS2 + "} [GeV/c]").c_str());
    h_C_proj->GetZaxis()->SetRangeUser(0.95, 1.3);
    h_C_proj->GetYaxis()->SetTitleOffset(0.8);

    dir->cd();
    h_C_proj->Write();

    delete h_A_proj; delete h_A_wei_proj;
}

void do_2d_proj(TFile* inputFile, TFile* outFile) {
    TDirectory* dir = outFile->GetDirectory("2d_projections");
    if (!dir) dir = outFile->mkdir("2d_projections");
    dir->cd();

    for (int chIdx = 0; chIdx < chargeSize; chIdx++) {
        for (int centIdx = 0; centIdx < centralitySize; centIdx++) {
            for (int ktIdx = 0; ktIdx < ktSize; ktIdx++) {
                for (int yIdx = 0; yIdx < rapiditySize; yIdx++) {
                    std::string ending = getPrefix(chIdx, centIdx, yIdx);
                    std::string name = "bp_" + ending;

                    TH3D* h_A = (TH3D*) inputFile->Get((name + std::to_string(ktIdx)).c_str());
                    TH3D* h_A_wei= (TH3D*) inputFile->Get((name + "wei_" + std::to_string(ktIdx)).c_str());

                    if (!h_A || !h_A_wei) {
                        std::cerr << "Warning: 3D histograms for name=" << name << " ktIdx=" << ktIdx << " not found!" << std::endl;
                        continue;;
                    }
                    
                    ending = ending + std::to_string(ktIdx);

                    for (int i = 0; i < lcmsSize; i++) {
                        std::string firstAxis = axises[i%lcmsSize], secondAxis = axises[(i+1)%lcmsSize];
                        std::string firstLCMS = LCMS[i%lcmsSize], secondLCMS = LCMS[(i+1)%lcmsSize];
                        if (i == lcmsSize - 1) {
                            std::swap(firstAxis, secondAxis);
                            std::swap(firstLCMS, secondLCMS);
                        }
                        std::string pr_axis = firstAxis + secondAxis;

                        get_2d_proj(dir, h_A, h_A_wei, pr_axis, firstLCMS, secondLCMS, ending);
                    }

                    delete h_A; delete h_A_wei;
                }
            }
        }
    }
}