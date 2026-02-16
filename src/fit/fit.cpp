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

TF3* CreateCF3DFit() {
    TF3* fit3d = new TF3(
        "fit3d",
        CF_fit_3d,
        -0.05, 0.05,
        -0.05, 0.05,
        -0.05, 0.05,
        7
    );
    fit3d->SetParameters(4.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.5);
    
    fit3d->SetParLimits(0, 0., 10.);
    fit3d->SetParLimits(1, 0., 10.);
    fit3d->SetParLimits(2, 0., 10.);
    fit3d->SetParLimits(3, -5., 5.);
    fit3d->SetParLimits(4, -5., 5.);
    fit3d->SetParLimits(5, -5., 5.);
    fit3d->SetParLimits(6, 0., 1.);

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