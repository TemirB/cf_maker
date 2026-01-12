#include <TH3.h>
#include <TF3.h>
#include <TSystem.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TCanvas.h>

#include "helpers.h"
#include "fit/types.h"
#include "fit/fit.h"

void prepare_for_1d_proj(TH3D* h_A_wei, TH3D* h_A, TH3D* denominator, TH3D* numerator, const std::string& LCMS){
    if (LCMS == "out") {
        h_A_wei->GetYaxis()->SetRangeUser(-0.05, 0.05);
        h_A_wei->GetZaxis()->SetRangeUser(-0.05, 0.05);
        h_A->GetYaxis()->SetRangeUser(-0.05, 0.05);
        h_A->GetZaxis()->SetRangeUser(-0.05, 0.05);
        denominator->GetYaxis()->SetRangeUser(-0.05, 0.05);
        denominator->GetZaxis()->SetRangeUser(-0.05, 0.05);
        numerator->GetYaxis()->SetRangeUser(-0.05, 0.05);
        numerator->GetZaxis()->SetRangeUser(-0.05, 0.05); 
    } else if (LCMS == "side") {
        h_A_wei->GetXaxis()->SetRangeUser(-0.05, 0.05);
        h_A_wei->GetZaxis()->SetRangeUser(-0.05, 0.05);
        h_A->GetXaxis()->SetRangeUser(-0.05, 0.05);
        h_A->GetZaxis()->SetRangeUser(-0.05, 0.05);
        denominator->GetXaxis()->SetRangeUser(-0.05, 0.05);
        denominator->GetZaxis()->SetRangeUser(-0.05, 0.05);
        numerator->GetXaxis()->SetRangeUser(-0.05, 0.05);
        numerator->GetZaxis()->SetRangeUser(-0.05, 0.05);
    } else if (LCMS == "long") {
    h_A_wei->GetXaxis()->SetRangeUser(-0.05, 0.05);
    h_A_wei->GetYaxis()->SetRangeUser(-0.05, 0.05);
    h_A->GetXaxis()->SetRangeUser(-0.05, 0.05);
    h_A->GetYaxis()->SetRangeUser(-0.05, 0.05);
    denominator->GetXaxis()->SetRangeUser(-0.05, 0.05);
    denominator->GetYaxis()->SetRangeUser(-0.05, 0.05);
    numerator->GetXaxis()->SetRangeUser(-0.05, 0.05);
    numerator->GetYaxis()->SetRangeUser(-0.05, 0.05);
    }
}

void get_1d_proj(
    TDirectory* dir ,TH3D* h_A, TH3D* h_A_wei, TF3* fit3d,
    const std::string& axis, const std::string& LCMS, const std::string& ending
) {
    TH3D* denominator = new TH3D("denom", "denom", 80, -0.4, 0.4, 80, -0.4, 0.4, 80, -0.4, 0.4); denominator->SetDirectory(nullptr);
    TH3D* numerator = new TH3D("numer", "numer", 80, -0.4, 0.4, 80, -0.4, 0.4, 80, -0.4, 0.4); numerator->SetDirectory(nullptr);
    // Заполняем numerator значениями аналитической функции (fit3d), denominator = 1 во всех ячейках
    for (int xbin = 1; xbin <= 80; ++xbin) {
        for (int ybin = 1; ybin <= 80; ++ybin) {
            for (int zbin = 1; zbin <= 80; ++zbin) {
                denominator->SetBinContent(xbin, ybin, zbin, 1.0);
                denominator->SetBinError(xbin, ybin, zbin, 0.0);
                double xval = numerator->GetXaxis()->GetBinCenter(xbin);
                double yval = numerator->GetYaxis()->GetBinCenter(ybin);
                double zval = numerator->GetZaxis()->GetBinCenter(zbin);
                double value = fit3d->Eval(xval, yval, zval);
                numerator->SetBinContent(xbin, ybin, zbin, value);
                numerator->SetBinError(xbin, ybin, zbin, 0.0); 
            }
        }
    }

    // хелпер
    auto resetRanges = [&](TH3D* h){
        h->GetXaxis()->SetRange(0,0);
        h->GetYaxis()->SetRange(0,0);
        h->GetZaxis()->SetRange(0,0);
    };

    resetRanges(h_A);
    resetRanges(h_A_wei);

    // Ограничение осей
    prepare_for_1d_proj(h_A_wei, h_A, denominator, numerator, LCMS);

    // Проекция + отвязка от директории
    TH1D* h_A_pr        = (TH1D*) h_A->Project3D(axis.c_str());         h_A_pr->SetDirectory(nullptr);
    TH1D* h_A_wei_pr    = (TH1D*) h_A_wei->Project3D(axis.c_str());     h_A_wei_pr->SetDirectory(nullptr);
    TH1D* h_denom_pr    = (TH1D*) denominator->Project3D(axis.c_str()); h_denom_pr->SetDirectory(nullptr);
    TH1D* h_num_pr      = (TH1D*) numerator->Project3D(axis.c_str());   h_num_pr->SetDirectory(nullptr); 
    delete denominator; delete numerator;

    // Гистограмма КФ и фита
    TH1D* h_CF = (TH1D*) h_A_wei_pr->Clone(("h_CF_" + LCMS + "_" + ending).c_str());    h_CF->Reset();  h_CF->SetDirectory(dir);
    TH1D* h_fit = (TH1D*) h_denom_pr->Clone(("h_fit_" + LCMS + "_" + ending).c_str());  h_fit->Reset(); h_fit->SetDirectory(dir);
    
    // Вычисляем КФ и фит
    h_CF->Divide(h_A_wei_pr, h_A_pr);
    h_fit->Divide(h_num_pr, h_denom_pr);
    delete h_A_pr; delete h_A_wei_pr; delete h_denom_pr; delete h_num_pr;

    // Рисование и сохранение
    TCanvas* c = new TCanvas(ending.c_str(), ("CF with fit at " + ending).c_str());

    c->cd();
    h_CF->Draw("ALP");
    h_fit->Draw("L SAME");

    dir->cd();
    c->Write();
    // h_CF->Write();
    // h_fit->Write();

    // c->SaveAs(("1d/cf_1d_" + ending).c_str());
    delete h_CF; delete h_fit; delete c;
}

void do_1d_proj(
    TFile* inputFile, TFile* outFile,
    FitResult (&fitRes)[chargeSize][centralitySize][ktSize][rapiditySize]
) {
    // Берем все 3d CF по очереди, для каждой из них мутим LCMS проекцию
    // Делаем фитирование аналитической функцией,
    // а затем сохрнаяем в директорию
    TDirectory* dir = outFile->GetDirectory("1d_projections");
    if (!dir) dir = outFile->mkdir("1d_projections");
    dir->cd();

    TF3* fit = CreateCF3DFit();

    // gSystem->Exec("mkdir -p 1d");

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
                        return;
                    }
                    
                    ending = ending + std::to_string(ktIdx);
                    FitResult r = fitRes[chIdx][centIdx][ktIdx][yIdx];

                    fit->SetParameter(0, r.R[0]); fit->SetParError(0, r.eR[0]); // R_out
                    fit->SetParameter(1, r.R[1]); fit->SetParError(1, r.eR[1]); // R_side
                    fit->SetParameter(2, r.R[2]); fit->SetParError(2, r.eR[2]); // R_long
                    fit->SetParameter(3, r.lambda); fit->SetParError(3, r.elambda); // lambda
                    for (int p = 0; p < 4; p++) {
                        fit->FixParameter(p, fit->GetParameter(p));
                    } 
                    
                    for (int i = 0; i < lcmsSize; i++) {
                        get_1d_proj(dir, h_A, h_A_wei, fit, axises[i], LCMS[i], ending);
                    }
                }
            }
        }
    }
}