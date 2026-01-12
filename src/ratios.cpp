#include <TDirectory.h>
#include <TH3D.h>
#include <TFile.h>

#include "fit.h"
#include "helpers.h"
#include "draw.h"

void ratio_proj(
    TDirectory* dir,
    TH3D* h3d_neg, TH3D* h3d_pos,
    const std::string& axis,
    const char* name
) {
    if (gDirectory) std::cerr << "gDirectory=" << gDirectory->GetName() << std::endl;
    else std::cerr << "gDirectory = nullptr" << std::endl;

    TH3D* hNegTmp = (TH3D*)h3d_neg->Clone("hNegTmp");
    TH3D* hPosTmp = (TH3D*)h3d_pos->Clone("hPosTmp");
    hNegTmp->SetDirectory(nullptr);
    hPosTmp->SetDirectory(nullptr);

    if (axis == "x") {
        hNegTmp->GetYaxis()->SetRangeUser(-0.05, 0.05);
        hNegTmp->GetZaxis()->SetRangeUser(-0.05, 0.05);
        hPosTmp->GetYaxis()->SetRangeUser(-0.05, 0.05);
        hPosTmp->GetZaxis()->SetRangeUser(-0.05, 0.05);
    } else if (axis == "y") {
        hNegTmp->GetXaxis()->SetRangeUser(-0.05, 0.05);
        hNegTmp->GetZaxis()->SetRangeUser(-0.05, 0.05);
        hPosTmp->GetXaxis()->SetRangeUser(-0.05, 0.05);
        hPosTmp->GetZaxis()->SetRangeUser(-0.05, 0.05);
    } else if (axis == "z") {
        hNegTmp->GetXaxis()->SetRangeUser(-0.05, 0.05);
        hNegTmp->GetYaxis()->SetRangeUser(-0.05, 0.05);
        hPosTmp->GetXaxis()->SetRangeUser(-0.05, 0.05);
        hPosTmp->GetYaxis()->SetRangeUser(-0.05, 0.05);
    }

    TH1D* h1d_neg = (TH1D*) hNegTmp->Project3D(axis.c_str());
    TH1D* h1d_pos = (TH1D*) hPosTmp->Project3D(axis.c_str());
    h1d_neg->SetDirectory(nullptr);
    h1d_pos->SetDirectory(nullptr);

    TH1D* hout = (TH1D*) h1d_neg->Clone(name);
    hout->Reset();
    hout->SetDirectory(nullptr);

    // Ошибку КФ-ов считать по биномиальному распределению для отдельных КФ
    // Ошибку отношения считать как 
    // 1 2
    // (e1/v1)^2 + (e2/v2)^2 = e(v1/v2)^2
    // ДЛЯ САМИХ НАЧАЛЬНЫХЪ КФ С УЧЕТОМ РАЗНЫХ ВЫБОРОК
    hout->Divide(h1d_neg, h1d_pos);
    hout->GetXaxis()->SetTitle(Form("q_{%s} (GeV/c)", axis.c_str()));
    hout->GetYaxis()->SetTitle("CF--/CF++");
    hout->SetTitle(name);

    dir->cd();
    hout->Write();

    // writeHist(
    //     dir,
    //     hout,
    //     name,
    //     Form("q_{%s} (GeV/c)", axis.c_str()),
    //     "CF--/CF++"
    // );

    delete h1d_neg;
    delete h1d_pos;
    delete hNegTmp;
    delete hPosTmp;
}

void proj_ratio(
    TDirectory* dir,
    TH3D* h3d,
    const std::string& axis,
    const char* name
) {
    TH3D* htmp = (TH3D*)h3d->Clone("h_tmp");
    htmp->SetDirectory(nullptr);

    if (axis == "x") {
        htmp->GetYaxis()->SetRangeUser(-0.05, 0.05);
        htmp->GetZaxis()->SetRangeUser(-0.05, 0.05);
    } else if (axis == "y") {
        htmp->GetXaxis()->SetRangeUser(-0.05, 0.05);
        htmp->GetZaxis()->SetRangeUser(-0.05, 0.05);
    } else if (axis == "z") {
        htmp->GetXaxis()->SetRangeUser(-0.05, 0.05);
        htmp->GetYaxis()->SetRangeUser(-0.05, 0.05);
    }

    TH1D* h1d = (TH1D*) htmp->Project3D(axis.c_str());
    h1d->SetDirectory(nullptr);   
    h1d->SetDirectory(dir);

    dir->cd();
    h1d->Write(name);
}

void do_CF_ratios(
    TFile* outFile,
    FitResult (&fitRes)[chargeSize][centralitySize][ktSize][rapiditySize]
) {
    // TODO
    // Отношение проекций CF 3D-гистограмм  proj(CF--)/proj(CF++)
    // Проекция отношения (CF--)/(CF++) по q_out, q_side, q_long 
    TDirectory* cf3d = outFile->GetDirectory("cf3d");
    TDirectory* dRatioProj = outFile->GetDirectory("cf_ratio_projs");
    if (!dRatioProj) dRatioProj = outFile->mkdir("cf_ratio_projs");

    TDirectory* dProjRatio = outFile->GetDirectory("cf_proj_ratios");
    if (!dProjRatio) dProjRatio = outFile->mkdir("cf_proj_ratios");

    TH3D* h_CF_ratio = nullptr;

    for (int centIdx = 0; centIdx < centralitySize; centIdx++) {
        for (int ktIdx = 0; ktIdx < ktSize; ktIdx++) {
            for (int yIdx = 0; yIdx < rapiditySize; yIdx++) {

                const std::string midNeg = getPrefix(0, centIdx, yIdx);
                const std::string midPos = getPrefix(1, centIdx, yIdx);

                const std::string nameNeg = "h3d_CF_q_" + midNeg + std::to_string(ktIdx) + "_weighted";
                const std::string namePos = "h3d_CF_q_" + midPos + std::to_string(ktIdx) + "_weighted";

                TH3D* hNegFile = nullptr;
                TH3D* hPosFile = nullptr;

                cf3d->GetObject(nameNeg.c_str(), hNegFile);
                cf3d->GetObject(namePos.c_str(), hPosFile);
                if (!hNegFile || !hPosFile) continue;

                TH3D* hNeg = (TH3D*) hNegFile->Clone();
                hNeg->SetDirectory(nullptr);
                hNeg->Reset("ICES");
                hNeg->Add(hNegFile);

                TH3D* hPos = (TH3D*) hPosFile->Clone();
                hPos->SetDirectory(nullptr);
                hPos->Reset("ICES");
                hPos->Add(hPosFile);


                if (!h_CF_ratio) {
                    h_CF_ratio = (TH3D*)hNeg->Clone("h_CF_ratio_work");
                    h_CF_ratio->SetDirectory(nullptr);
                    if (!h_CF_ratio->GetSumw2N()) h_CF_ratio->Sumw2();
                } else if (!sameBinning(h_CF_ratio, hNeg)) {
                    delete h_CF_ratio;
                    h_CF_ratio = (TH3D*)hNeg->Clone("h_CF_ratio_work");
                    h_CF_ratio->SetDirectory(nullptr);
                    if (!h_CF_ratio->GetSumw2N()) h_CF_ratio->Sumw2();
                }

                h_CF_ratio->Reset("ICES");
                if (hNeg->GetEntries() == 0 || hPos->GetEntries() == 0) {
                    std::cerr << "EMPTY input CF ratio "
                            << nameNeg << " / " << namePos << std::endl;
                    delete hNeg;
                    delete hPos;
                    continue;
                }

                if (!hNeg->GetSumw2N()) hNeg->Sumw2();
                if (!hPos->GetSumw2N()) hPos->Sumw2();

                h_CF_ratio->Divide(hNeg, hPos); // CF-- / CF++

                for (int lcmsIdx = 0; lcmsIdx < lcmsSize; ++lcmsIdx) {
                    ratio_proj(
                        dRatioProj,
                        hNeg, hPos,
                        axises[lcmsIdx],
                        Form("ratioProj_%s_cent%d_y%d_kt%d",
                             LCMS[lcmsIdx].c_str(), centIdx, yIdx, ktIdx)
                    );

                    TH3D* tmp = (TH3D*) h_CF_ratio->Clone("tmp_proj");
                    tmp->SetDirectory(nullptr);

                    proj_ratio(
                        dProjRatio,
                        tmp,
                        axises[lcmsIdx],
                        Form("projRatio_%s_cent%d_y%d_kt%d",
                             LCMS[lcmsIdx].c_str(), centIdx, yIdx, ktIdx)
                    );
                    
                    delete tmp;
                }

                delete hNeg;
                delete hPos;
            }
        }
    }

    delete h_CF_ratio;
}