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

std::map<LCMSAxis, std::map<int, std::vector<double>>> initialFitPoints = {
    {
        LCMS::Out, {
            {0, {6.2, 5.9, 5.7, 5.4}}, // 3?
            {1, {5.35, 5.15, 4.9, 4.75}}, // 3?
            {2, {4.6, 4.4, 4.2, 4.0}}, // 23?
            {3, {4.2, 4.0, 3.8, 3.6}} // 123?
        }
    },
    {
        LCMS::Side, {
            {0, {4.8, 4.1, 3.8, 3.6}}, // 3?
            {1, {4.1, 3.6, 3.15, 2.8}},
            {2, {3.4, 3.0, 2.8, 2.6}}, // 3?
            {3, {2.8, 2.5, 2.3, 2.15}} // 123?
        }
    }
    {
        {
        LCMS::Long, {
            {0, {4.5, 3.5, 3.0, 2.8}}, // 3?
            {1, {4.3, 3.35, 2.85, 2.6}},
            {2, {3.67, 2.86, 2.5, 2.3}}, // 23?
            {3, {2.8, 2.5, 2.3, 2.15}} // 0123?
        }
    }
    },
};

TF3* CreateCF3DFit(int centrality, int kt) {
    double fitLim = 0.2;
    // if (centrality + kt > 4) {
    //     fitLim = 0.15;
    // } else if (kt == 3) {
    //     fitLim = 0.1;
    // }
    TF3* fit3d = new TF3(
        "fit3d", CF_fit_3d, 
        -fitLim, fitLim,
        -fitLim, fitLim,
        -fitLim, fitLim,
        7
    );

    double Roout = initialFitPoint[LCMS::Out][centrality][kt];
    double Rside = initialFitPoint[LCMS::Side][centrality][kt];
    double Rlong = initialFitPoint[LCMS::Long][centrality][kt];
    
    double Routside, Routlong, Rsidelong, Lambda;
    
    switch (centrality) {
        case 0:
            Lambda = 0.88;
            break;
        case 1:
            Lambda = 0.89;
            break;
        case 2:
            Lambda = 0.90;
            break;
        case 3:
            Lambda = 0.92;
            break;
        default:
            Lambda = 0.88;
    }
    
    Routside = 0; // OS
    Routlong = 0; // OL
    Rsidelong = 0; // SL
    
    fit3d->SetParameters(Rout, Rside, Rlong, Routside, Routlong, Rsidelong, Lambda);

    fit3d->FixParameter(3, 0.0); // OS
    fit3d->FixParameter(5, 0.0); // SL
    
    fit3d->SetParLimits(0, 0.0, 10.0);
    fit3d->SetParLimits(1, 0.0, 10.0);
    fit3d->SetParLimits(2, 0.0, 10.0);
    fit3d->SetParLimits(3, -30., 30.);
    fit3d->SetParLimits(4, -30., 30.);
    fit3d->SetParLimits(5, -30., 30.);
    fit3d->SetParLimits(6, 0.0, 1.0);
    
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