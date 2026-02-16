#include "fit/fit.h"

#include <TSystem.h>
#include <TFitResult.h>
#include <TMath.h>

const Double_t hc2 = 0.197 * 0.197;

Double_t CF_fit_3d(Double_t* q, Double_t* par) {
    Double_t q_out = q[0];
    Double_t q_side = q[1];
    Double_t q_long = q[2];

    Double_t qRq = par[0]*par[0] * q_out*q_out +
                   par[1]*par[1] * q_side*q_side +
                   par[2]*par[2] * q_long*q_long +
                 2.*par[3]*par[3] * q_out*q_side +
                 2.*par[4]*par[4] * q_out*q_long +
                 2.*par[5]*par[5] * q_side*q_long;

    return 1.0 + par[6] * TMath::Exp(-qRq / hc2);
}

TF3* CreateCF3DFit(int centrality, int kt) {
    TF3* fit3d = new TF3(
        "fit3d", CF_fit_3d, 
        -0.05, 0.05,
        -0.05, 0.05,
        -0.05, 0.05,
        7
    );
    
    double p0, p1, p2, p3, p4, p5, p6;
    
    switch (centrality) {
        case 0: // 0-5%
            p0 = 6.0 - 0.30*kt;
            p1 = 4.6 - 0.5*kt;
            p2 = 4.5 - 0.25*kt;
            p6 = 0.88;
            break;
        case 1: // 5-10%
            p0 = 5.3 - 0.35*kt;
            p1 = 4.0 - 0.30*kt;
            p2 = 3.8 - 0.30*kt;
            p6 = 0.89;
            break;
        case 2: // 10-20%
            p0 = 4.5 - 0.50*kt;
            p1 = 3.3 - 0.40*kt;
            p2 = 3.2 - 0.40*kt;
            p6 = 0.90;
            break;
        case 3: // 20-30%
            p0 = 3.8 - 0.70*kt;  // 3.8 → 1.7 при kt=3
            p1 = 2.7 - 0.55*kt;  // 2.7 → 1.05
            p2 = 2.6 - 0.60*kt;  // 2.6 → 0.8
            p6 = 0.92;
            break;
        default:
            p0 = 5.5 - 0.4*kt;
            p1 = 4.0 - 0.3*kt;
            p2 = 3.6 - 0.3*kt;
            p6 = 0.88;
    }
    
    p3 = 0.10 * TMath::Sqrt(p0 * p1);
    p4 = 0.05 * TMath::Sqrt(p0 * p2);
    p5 = 0.20 * TMath::Sqrt(p1 * p2);
    
    fit3d->SetParameters(p0, p1, p2, p3, p4, p5, p6);
    
    switch (centrality) {
        case 0: // 0-5%
            fit3d->SetParLimits(0, 4.0, 8.5);   // R_out
            fit3d->SetParLimits(1, 2.8, 6.8);   // R_side
            fit3d->SetParLimits(2, 2.8, 6.5);   // R_long
            break;
        case 1: // 5-10%
            fit3d->SetParLimits(0, 3.0, 7.5);
            fit3d->SetParLimits(1, 2.2, 5.8);
            fit3d->SetParLimits(2, 2.1, 5.6);
            break;
        case 2: // 10-20%
            fit3d->SetParLimits(0, 2.0, 6.5);
            fit3d->SetParLimits(1, 1.4, 5.0);
            fit3d->SetParLimits(2, 1.3, 4.8);
            break;
        case 3: // 20-30%
            fit3d->SetParLimits(0, 0.8, 4.8);
            fit3d->SetParLimits(1, 0.7, 4.2);
            fit3d->SetParLimits(2, 0.6, 4.0);
            break;
        default:
            fit3d->SetParLimits(0, 3.5, 8.5);
            fit3d->SetParLimits(1, 2.5, 6.8);
            fit3d->SetParLimits(2, 2.4, 6.5);
    }
    
    fit3d->SetParLimits(3, -30., 30.);
    fit3d->SetParLimits(4, -30., 30.);
    fit3d->SetParLimits(5, -30., 30.);
    
    // Lambda — ПОДНЯТЬ нижнюю границу!
    fit3d->SetParLimits(6, 0.7, 1.0);
    
    fit3d->SetParName(0, "R_out");
    fit3d->SetParName(1, "R_side");
    fit3d->SetParName(2, "R_long");
    fit3d->SetParName(3, "R_os");
    fit3d->SetParName(4, "R_ol");
    fit3d->SetParName(5, "R_sl");
    fit3d->SetParName(6, "lambda");
    
    return fit3d;
}

FitResult FitCF3D(TH3D* hCF, TF3* fit3d) {
    FitResult res{};
    if (!hCF || hCF->GetEntries() == 0 || !fit3d) return res;

    // ROOT::Math::MinimizerOptions::SetDefaultStrategy(2); // более агрессивный поиск
    // ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-3); // смягчить критерий сходимости
    // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Simplex");
    // fit3d->FixParameter(6, 0.9);
    auto fitPtr = hCF->Fit(fit3d, "RMQS0");

    gSystem->RedirectOutput("fit_result.txt", "w");
    fitPtr->Print("V");
    gSystem->RedirectOutput(0);

    if (fitPtr.Get()) {
        res.chi2   = fitPtr->Chi2();
        res.ndf    = fitPtr->Ndf();
        res.pvalue = TMath::Prob(res.chi2, res.ndf);
        res.ok     = (res.chi2 >= 0 && res.ndf > 0);
    }

    for (int i = 0; i < 7; i++) {
        Double_t val = fit3d->GetParameter(i);
        Double_t err = fit3d->GetParError(i);

        if (i == 6) {
            res.lambda = val;
            res.elambda = err;
        } else {
            res.R[i] = val;
            res.eR[i] = err;
        } 
    }

    return res;
}