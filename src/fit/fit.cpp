#include "fit/fit.h"

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

TF3* CreateCF3DFit(int centrality, int rapidity) {
    TF3* fit3d = new TF3("fit3d", CF_fit_3d, -0.05, 0.05, -0.05, 0.05, -0.05, 0.05, 7);
    
    double R_out, R_side, R_long, R_os, R_ol, R_sl, lambda;
    
    switch (centrality) {
        case 0: // 0-5% (центральные)
            R_out = 6.0; R_side = 4.5; R_long = 4.8; lambda = 0.87;
            break;
        case 1: // 5-10%
            R_out = 5.5; R_side = 4.1; R_long = 4.4; lambda = 0.88;
            break;
        case 2: // 10-20%
            R_out = 4.8; R_side = 3.6; R_long = 3.9; lambda = 0.89;
            break;
        case 3: // 20-30% (периферия) — КЛЮЧЕВО: низкие значения!
            R_out = 4.2; R_side = 3.1; R_long = 3.4; lambda = 0.91;
            break;
        default:
            R_out = 5.5; R_side = 4.0; R_long = 4.5; lambda = 0.88;
    }
    
    R_os = 0.0;
    R_ol = 0.0;
    R_sl = 0.0;
    
    fit3d->SetParameters(R_out, R_side, R_long, R_os, R_ol, R_sl, lambda);
    
    switch (centrality) {
        case 0: // 0-5%
            fit3d->SetParLimits(0, 4.5, 7.2); // R_out
            fit3d->SetParLimits(1, 3.2, 5.8); // R_side
            fit3d->SetParLimits(2, 3.3, 6.0); // R_long
            break;
        case 1: // 5-10%
            fit3d->SetParLimits(0, 4.0, 6.5);
            fit3d->SetParLimits(1, 2.9, 5.2);
            fit3d->SetParLimits(2, 3.0, 5.3);
            break;
        case 2: // 10-20%
            fit3d->SetParLimits(0, 3.5, 5.8);
            fit3d->SetParLimits(1, 2.5, 4.6);
            fit3d->SetParLimits(2, 2.6, 4.7);
            break;
        case 3: // 20-30%
            fit3d->SetParLimits(0, 2.8, 5.2); // R_out
            fit3d->SetParLimits(1, 2.2, 4.2); // R_side
            fit3d->SetParLimits(2, 2.3, 4.3); // R_long
            break;
        default:
            fit3d->SetParLimits(0, 4.0, 7.0);
            fit3d->SetParLimits(1, 2.8, 5.5);
            fit3d->SetParLimits(2, 2.9, 5.6);
    }
    
    fit3d->SetParLimits(3, -3.0, 3.0); // R_os
    fit3d->SetParLimits(4, -3.0, 3.0); // R_ol
    fit3d->SetParLimits(5, -3.0, 3.0); // R_sl
    
    fit3d->SetParLimits(6, 0.6, 1.0); // lambda
    
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

    auto fitPtr = hCF->Fit(fit3d, "RMQS0");

    if (fitPtr.Get()) {
        res.chi2   = fitPtr->Chi2();
        res.ndf    = fitPtr->Ndf();
        res.pvalue = TMath::Prob(res.chi2, res.ndf);
        res.ok     = (res.chi2 >= 0 && res.ndf > 0);
    }

    for (int i = 0; i < 7; i++) {
        Double_t val = fit3d->GetParameter(i);
        Double_t err = fit3d->GetParError(i);
        if (i < 6) {
            res.R[i] = val;
            res.eR[i] = err;
        } else if (i == 6) {
            res.lambda = val;
            res.elambda = err;
        }
    }

    return res;
}