#include "fit/fit.h"

#include <TFitResult.h>
#include <TMath.h>

Double_t CF_fit_3d(Double_t* q, Double_t* par) {
    const Double_t hc2 = 0.197 * 0.197;

    return 1.0 + par[3] * TMath::Exp(
        (-(q[0] * par[0]) * (q[0] * par[0])
         -(q[1] * par[1]) * (q[1] * par[1])
         -(q[2] * par[2]) * (q[2] * par[2])) / hc2
    );
}

TF3* CreateCF3DFit() {
    TF3* fit3d = new TF3(
        "fit3d",
        CF_fit_3d,
        -0.1, 0.1,
        -0.1, 0.1,
        -0.1, 0.1,
        4
    );

    fit3d->SetParameters(4.0, 4.0, 4.0, 0.5);

    fit3d->SetParLimits(0, 1.0, 10.0);
    fit3d->SetParLimits(1, 1.0, 10.0);
    fit3d->SetParLimits(2, 1.0, 10.0);
    fit3d->SetParLimits(3, 0.1, 1.0);

    return fit3d;
}

FitResult FitCF3D(TH3D* hCF, TF3* fit3d) {
    FitResult res{};
    if (!hCF || hCF->GetEntries()==0 || !fit3d) return res;

    auto fitPtr = hCF->Fit(fit3d, "RMQS0");

    if (fitPtr.Get()) {
        res.chi2   = fitPtr->Chi2();
        res.ndf    = fitPtr->Ndf();
        res.pvalue = TMath::Prob(res.chi2, res.ndf);
        res.ok     = true;
    }

    for (int i=0;i<3;i++) {
        res.R[i]  = fit3d->GetParameter(i);
        res.eR[i] = fit3d->GetParError(i);
    }
    res.lambda  = fit3d->GetParameter(3);
    res.elambda = fit3d->GetParError(3);

    return res;
}
