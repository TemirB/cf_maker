#include "fit/model.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>

#include <TFitResult.h>
#include <TMath.h>
#include <TMatrixDSym.h>

#include "core/log.h"
#include "fit/initial_parameters.h"

namespace
{
constexpr double kHc2 = 0.197 * 0.197;

double clamp_sqrt(double v)
{
    return std::sqrt(std::max(0.0, v));
}

const InitialParameters& kt_ip()
{
    static const InitialParameters ip(Binning::Kt);
    return ip;
}

const InitialParameters& rapidity_ip()
{
    static const InitialParameters ip(Binning::Rapidity);
    return ip;
}

const InitialParameters& default_kt_ip()
{
    static const InitialParameters ip(Binning::Kt, true);
    return ip;
}

const InitialParameters& default_rapidity_ip()
{
    static const InitialParameters ip(Binning::Rapidity, true);
    return ip;
}

TF3* make_cf_3d_fit(const Config& cfg, int ch, int centr, int b, bool use_defaults)
{
    const double fit_limit = cfg.fit.q_max;

    TF3* fit3d = new TF3("fit3d", cf_fit_3d, -fit_limit, fit_limit, -fit_limit, fit_limit, -fit_limit, fit_limit, 7);
    const bool is_kt = (cfg.input.type == "kt");

    const InitialParameters& ip =
        use_defaults ? (is_kt ? default_kt_ip() : default_rapidity_ip()) : (is_kt ? kt_ip() : rapidity_ip());

    const double r_out = ip.get(ch, "out", centr, b);
    const double r_side = ip.get(ch, "side", centr, b);
    const double r_long = ip.get(ch, "long", centr, b);
    const double r_out_side = 0;
    const double r_out_long = 0;
    const double r_side_long = 0;
    const double lambda = ip.get(ch, "lambda", centr, b);

    fit3d->SetParameters(r_out * r_out, r_side * r_side, r_long * r_long, r_out_side, r_out_long, r_side_long,
                         lambda);

    fit3d->SetParLimits(0, cfg.fit.radius_sq_min, cfg.fit.radius_sq_max);
    fit3d->SetParLimits(1, cfg.fit.radius_sq_min, cfg.fit.radius_sq_max);
    fit3d->SetParLimits(2, cfg.fit.radius_sq_min, cfg.fit.radius_sq_max);
    fit3d->SetParLimits(3, cfg.fit.cross_min, cfg.fit.cross_max);
    fit3d->SetParLimits(4, cfg.fit.cross_min, cfg.fit.cross_max);
    fit3d->SetParLimits(5, cfg.fit.cross_min, cfg.fit.cross_max);
    fit3d->SetParLimits(6, cfg.fit.lambda_min, cfg.fit.lambda_max);

    fit3d->SetParError(0, 1.0);
    fit3d->SetParError(1, 1.0);
    fit3d->SetParError(2, 1.0);
    fit3d->SetParError(3, 1.0);
    fit3d->SetParError(4, 1.0);
    fit3d->SetParError(5, 1.0);
    fit3d->SetParError(6, 0.1);

    fit3d->SetParName(0, "R_out");
    fit3d->SetParName(1, "R_side");
    fit3d->SetParName(2, "R_long");
    fit3d->SetParName(3, "R_os");
    fit3d->SetParName(4, "R_ol");
    fit3d->SetParName(5, "R_sl");
    fit3d->SetParName(6, "lambda");

    return fit3d;
}

void param_limits(const FitConfig& fitCfg, int i, double& min, double& max)
{
    if (i < 3) {
        min = fitCfg.radius_sq_min;
        max = fitCfg.radius_sq_max;
    } else if (i < 6) {
        min = fitCfg.cross_min;
        max = fitCfg.cross_max;
    } else {
        min = fitCfg.lambda_min;
        max = fitCfg.lambda_max;
    }
}

} // namespace

double eval_cf_3d(const FitResult& r, double qOut, double qSide, double qLong)
{
    // FitResult stores physical radii for out/side/long, while model
    // uses coefficients R^2 in front of q^2 terms.
    const double qRq = (r.r[0] * r.r[0]) * qOut * qOut + (r.r[1] * r.r[1]) * qSide * qSide +
                       (r.r[2] * r.r[2]) * qLong * qLong + 2.0 * r.r[3] * qOut * qSide +
                       2.0 * r.r[4] * qOut * qLong + 2.0 * r.r[5] * qSide * qLong;

    return 1.0 + r.lambda * TMath::Exp(-qRq / kHc2);
}

Double_t cf_fit_3d(Double_t* q, Double_t* par)
{
    const Double_t q_out = q[0];
    const Double_t q_side = q[1];
    const Double_t q_long = q[2];

    const Double_t qRq = par[0] * q_out * q_out + par[1] * q_side * q_side +
                         par[2] * q_long * q_long + 2.0 * par[3] * q_out * q_side +
                         2.0 * par[4] * q_out * q_long + 2.0 * par[5] * q_side * q_long;

    return 1.0 + par[6] * TMath::Exp(-qRq / kHc2);
}

TF3* create_cf_3d_fit(const Config& cfg, int ch, int centr, int b)
{
    return make_cf_3d_fit(cfg, ch, centr, b, cfg.fit.use_default_ip);
}

FitResult fit_cf_3d(TH3D* cf_hist, TF3* fit3d, const FitConfig& fitCfg)
{
    FitResult res{};
    if (!cf_hist || cf_hist->GetEntries() == 0 || !fit3d) {
        return res;
    }

    for (int i = 0; i < 7; i++) {
        const auto& frozen = fitCfg.freeze[i];
        if (!frozen.has_value()) {
            continue;
        }
        double value = frozen.value();
        if (i < 3) {
            value = value * value;
        }
        fit3d->FixParameter(i, value);
    }

    std::string options = fitCfg.options;
    if (fitCfg.use_integral) {
        options += "I";
    }
    if (fitCfg.minos_errors) {
        options += "E";
    }

    auto fit_ptr = cf_hist->Fit(fit3d, options.c_str());

    cf_hist->GetListOfFunctions()->Remove(fit3d);

    if (fit_ptr.Get()) {
        res.chi2 = fit_ptr->Chi2();
        res.ndf = fit_ptr->Ndf();
        res.p_value = fit_ptr->Prob();
        res.status = fit_ptr->Status();
        res.cov_status = fit_ptr->CovMatrixStatus();
        res.ok = (res.chi2 >= 0 && res.ndf > 0) && (res.status == 0) && fit_ptr->is_valid() &&
                 (res.cov_status == 1 || res.cov_status == 3);
    }
    res.attempts = 1;

    for (int i = 0; i < 7; i++) {
        const Double_t val = fit3d->GetParameter(i);
        const Double_t err = fit3d->GetParError(i);
        if (i < 3) {
            const double radius = clamp_sqrt(val);
            res.r[i] = radius;
            res.e_r[i] = (radius > 0.0) ? std::abs(err) / (2.0 * radius) : 0.0;
        } else if (i < 6) {
            res.r[i] = val;
            res.e_r[i] = err;
        } else {
            res.lambda = val;
            res.e_lambda = err;
        }
    }

    for (int i = 0; i < 7; i++) {
        if (fitCfg.freeze[i].has_value()) {
            continue;
        }
        double min = 0.0;
        double max = 0.0;
        param_limits(fitCfg, i, min, max);
        const double val = fit3d->GetParameter(i);
        constexpr double kEps = 1e-4;
        if (val <= min + kEps || val >= max - kEps) {
            res.at_limit = true;
        }
    }

    if (fit_ptr.Get() && (res.cov_status == 1 || res.cov_status == 3)) {
        const TMatrixDSym cov = fit_ptr->GetCovarianceMatrix();
        if (cov.GetNrows() == 7) {
            std::array<double, 7> scale{};
            for (int i = 0; i < 7; i++) {
                scale[i] = (i < 3 && res.r[i] > 0.0) ? 1.0 / (2.0 * res.r[i]) : 1.0;
            }
            std::array<double, 7> sigma{};
            for (int i = 0; i < 7; i++) {
                sigma[i] = (i < 3) ? res.e_r[i] : (i == 6 ? res.e_lambda : res.e_r[i]);
            }
            for (int i = 0; i < 7; i++) {
                for (int j = 0; j < 7; j++) {
                    const double c = cov(i, j) * scale[i] * scale[j];
                    res.corr[static_cast<std::size_t>(i) * 7 + static_cast<std::size_t>(j)] =
                        (sigma[i] > 0.0 && sigma[j] > 0.0) ? c / (sigma[i] * sigma[j]) : 0.0;
                }
            }
        }
    }

    return res;
}

FitResult fit_cf_3d_with_retry(TH3D* cf_hist, const Config& cfg, int ch, int centr, int b)
{
    const std::string task = "fit (ch=" + std::to_string(ch) + ", centr=" + std::to_string(centr) +
                             ", b=" + std::to_string(b) + ")";

    auto fit3d = std::unique_ptr<TF3>(create_cf_3d_fit(cfg, ch, centr, b));
    FitResult best = fit_cf_3d(cf_hist, fit3d.get(), cfg.fit);
    log::Debug(task + ": chi2=" + std::to_string(best.chi2) + ", ndf=" + std::to_string(best.ndf) +
               ", status=" + std::to_string(best.status) + ", covStatus=" +
               std::to_string(best.cov_status) + ", atLimit=" + (best.at_limit ? "true" : "false"));

    if ((best.ok && !best.at_limit) || !cfg.fit.retry_with_defaults) {
        return best;
    }

    log::Info(
        task + ": first attempt " +
        (best.ok ? "hit parameter limits" : "failed (status=" + std::to_string(best.status) + ")") +
        " — retrying with default initial parameters");

    auto altFit = std::unique_ptr<TF3>(make_cf_3d_fit(cfg, ch, centr, b, true));
    FitResult alt = fit_cf_3d(cf_hist, altFit.get(), cfg.fit);

    if (alt.ok && (!best.ok || alt.chi2 < best.chi2)) {
        alt.attempts = 2;
        log::Info(task + ": retry improved the fit");
        return alt;
    }
    log::Info(task + ": retry did not improve — keeping first attempt");
    return best;
}
