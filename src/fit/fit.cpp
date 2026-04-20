#include "fit/fit.h"

#include <cstring>
#include <algorithm>
#include <cmath>

#include <TFitResult.h>
#include <TMath.h>

#include <fit/initial_parameters.h>
#include <helpers.h>
#include <context.h>

namespace {
constexpr double kHc2 = 0.197 * 0.197;

double ClampSqrt(double v) {
    return std::sqrt(std::max(0.0, v));
}
}

double EvalCF3D(const FitResult& r, double qOut, double qSide, double qLong) {
    // FitResult stores physical radii for out/side/long, while model
    // uses coefficients R^2 in front of q^2 terms.
    const double qRq = (r.R[0] * r.R[0]) * qOut * qOut +
                       (r.R[1] * r.R[1]) * qSide * qSide +
                       (r.R[2] * r.R[2]) * qLong * qLong +
                       2.0 * r.R[3] * qOut * qSide +
                       2.0 * r.R[4] * qOut * qLong +
                       2.0 * r.R[5] * qSide * qLong;

    return 1.0 + r.lambda * TMath::Exp(-qRq / kHc2);
}

Double_t CF_fit_3d(Double_t* q, Double_t* par) {
    const Double_t q_out = q[0];
    const Double_t q_side = q[1];
    const Double_t q_long = q[2];

    const Double_t qRq = par[0] * q_out * q_out +
                         par[1] * q_side * q_side +
                         par[2] * q_long * q_long +
                         2.0 * par[3] * q_out * q_side +
                         2.0 * par[4] * q_out * q_long +
                         2.0 * par[5] * q_side * q_long;

    return 1.0 + par[6] * TMath::Exp(-qRq / kHc2);
}

static InitialParameters ktIp = [](){
    InitialParameters tmp;
    tmp.KtInitialParameters();
    return tmp;
}();
static InitialParameters rapIp = [](){
    InitialParameters tmp;
    tmp.RapidityInitialParameters();
    return tmp;
}();

TF3* CreateCF3DFit(Context ctx, int ch, int centr, int b) {
    double fitLim = 0.20;

    TF3* fit3d = new TF3(
        "fit3d", CF_fit_3d, 
        -fitLim, fitLim, 
        -fitLim, fitLim, 
        -fitLim, fitLim, 
        7
    );

    const InitialParameters& ip = (std::strcmp(ctx.bining.type, "kt") == 0) 
                                ? ktIp 
                                : rapIp;
    
    double Rout = ip.get(ch, "out", centr, b);
    double Rside = ip.get(ch, "side", centr, b);
    double Rlong = ip.get(ch, "long", centr, b);
    double Routside = 0;
    double Routlong = 0;
    double Rsidelong = 0;
    double Lambda = ip.get(ch, "lambda", centr, b);
    
    fit3d->SetParameters(Rout * Rout, Rside * Rside, Rlong * Rlong, Routside, Routlong, Rsidelong, Lambda);

    fit3d->SetParLimits(0, 0.0, 50.);
    fit3d->SetParLimits(1, 0.0, 50.);
    fit3d->SetParLimits(2, 0.0, 50.);
    fit3d->SetParLimits(3, -30., 30.);
    fit3d->SetParLimits(4, -30., 30.);
    fit3d->SetParLimits(5, -30., 30.);
    fit3d->SetParLimits(6, 0.7, 1.0);
    
    // Имена параметров
    fit3d->SetParName(0, "R_out");
    fit3d->SetParName(1, "R_side");
    fit3d->SetParName(2, "R_long");
    fit3d->SetParName(3, "R_os");
    fit3d->SetParName(4, "R_ol");
    fit3d->SetParName(5, "R_sl");
    fit3d->SetParName(6, "lambda");

    // fit3d->FixParameter(0, Rout);
    // fit3d->FixParameter(1, Rside);
    // fit3d->FixParameter(2, Rlong);
    // fit3d->FixParameter(6, Lambda);
    
    return fit3d;
}

FitResult FitCF3D(TH3D* hCF, TF3* fit3d) {
    FitResult res{};
    if (!hCF || hCF->GetEntries() == 0 || !fit3d) return res;

    fit3d->FixParameter(3, 0);
    // fit3d->FixParameter(4, 0);
    fit3d->FixParameter(5, 0);
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
        if (i < 3) {
            const double radius = ClampSqrt(val);
            res.R[i] = radius;
            // sigma(sqrt(x)) = sigma(x) / (2*sqrt(x))
            res.eR[i] = (radius > 0.0) ? (std::abs(err) / (2.0 * radius)) : 0.0;
        } else if (i < 6) {
            res.R[i] = val;
            res.eR[i] = err;
        } else if (i == 6) {
            res.lambda = val;
            res.elambda = err;
        }
    }

    return res;
}