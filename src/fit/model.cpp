#include "fit/model.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <string>

#include <Fit/BinData.h>
#include <Fit/DataOptions.h>
#include <Fit/DataRange.h>
#include <HFitInterface.h>
#include <Math/ProbFuncMathCore.h>
#include <TFitResult.h>
#include <TList.h>
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

    TF3* fit3d = new TF3("fit3d", cf_fit_3d, -fit_limit, fit_limit, -fit_limit, fit_limit,
                         -fit_limit, fit_limit, 7);
    const bool is_kt = (cfg.input.type == "kt");

    const InitialParameters& ip = use_defaults ? (is_kt ? default_kt_ip() : default_rapidity_ip())
                                               : (is_kt ? kt_ip() : rapidity_ip());

    const double r_out = ip.get(ch, "out", centr, b);
    const double r_side = ip.get(ch, "side", centr, b);
    const double r_long = ip.get(ch, "long", centr, b);
    const double r_out_side = 0;
    const double r_out_long = 0;
    const double r_side_long = 0;
    const double lambda = ip.get(ch, "lambda", centr, b);

    fit3d->SetParameters(r_out * r_out, r_side * r_side, r_long * r_long, r_out_side, r_out_long,
                         r_side_long, lambda);

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

std::optional<FitStatistics> calculate_fit_statistics(const TH3D& cf_hist, TF3& fit3d,
                                                      const FitConfig& fit_cfg)
{
    std::string options = fit_cfg.options;
    std::transform(options.begin(), options.end(), options.begin(), [](unsigned char character) {
        return static_cast<char>(std::toupper(character));
    });
    // Execution modifiers contain letters which are also single-character fit
    // options (e.g. I/L/R in SERIAL); they do not change the objective.
    for (const std::string modifier : {"MULTITHREAD", "MULTIPROCESS", "SERIAL"}) {
        std::size_t position;
        while ((position = options.find(modifier)) != std::string::npos) {
            options.erase(position, modifier.size());
        }
    }
    if (options.find_first_of("LPU") != std::string::npos ||
        options.find("WIDTH") != std::string::npos) {
        return std::nullopt;
    }

    ROOT::Fit::DataOptions data_options;
    data_options.fIntegral = fit_cfg.use_integral || options.find('I') != std::string::npos;
    data_options.fUseRange = options.find('R') != std::string::npos;
    data_options.fErrors1 = options.find('W') != std::string::npos;
    data_options.fUseEmpty = options.find("WW") != std::string::npos;
    ROOT::Fit::DataRange range;
    if (data_options.fUseRange) {
        std::array<double, 3> lower{}, upper{};
        fit3d.GetRange(lower[0], lower[1], lower[2], upper[0], upper[1], upper[2]);
        for (unsigned int axis = 0; axis < 3; ++axis) {
            range.AddRange(axis, lower[axis], upper[axis]);
        }
    }
    // Reuse only ROOT's bin selection, including zero-error exclusions and
    // active TAxis ranges. The chi-square sum itself is evaluated here.
    ROOT::Fit::BinData data(data_options, range);
    ROOT::Fit::FillData(data, &cf_hist, &fit3d);
    long double chi2 = 0.0L;
    for (unsigned int i = 0; i < data.Size(); ++i) {
        double value = 0.0, inverse_error = 0.0;
        const double* coordinates = data.GetPoint(i, value, inverse_error);
        double expected;
        if (data_options.fIntegral) {
            std::array<double, 3> upper{};
            data.GetBinUpEdgeCoordinates(i, upper.data());
            const double volume = (upper[0] - coordinates[0]) * (upper[1] - coordinates[1]) *
                                  (upper[2] - coordinates[2]);
            expected = fit3d.Integral(coordinates[0], upper[0], coordinates[1], upper[1],
                                      coordinates[2], upper[2], 1e-9) /
                       volume;
        } else {
            expected = fit3d.Eval(coordinates[0], coordinates[1], coordinates[2]);
        }
        const long double residual = (static_cast<long double>(value) - expected) * inverse_error;
        chi2 += residual * residual;
    }

    FitStatistics statistics;
    statistics.chi2 = static_cast<double>(chi2);
    statistics.ndf = std::max(0, static_cast<int>(data.Size()) - fit3d.GetNumberFreeParameters());
    if (!std::isfinite(statistics.chi2)) {
        statistics.p_value = std::numeric_limits<double>::quiet_NaN();
    } else if (statistics.ndf > 0) {
        // P(Chi-square_ndf >= chi2) = Gamma(ndf/2, chi2/2) / Gamma(ndf/2).
        statistics.p_value = ROOT::Math::chisquared_cdf_c(statistics.chi2, statistics.ndf);
    }
    return statistics;
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
        if (const auto statistics = calculate_fit_statistics(*cf_hist, *fit3d, fitCfg)) {
            res.chi2 = statistics->chi2;
            res.ndf = statistics->ndf;
            res.p_value = statistics->p_value;
        } else {
            res.chi2 = fit_ptr->Chi2();
            res.ndf = fit_ptr->Ndf();
            res.p_value = fit_ptr->Prob();
        }
        res.status = fit_ptr->Status();
        res.cov_status = fit_ptr->CovMatrixStatus();
        res.ok = (res.chi2 >= 0 && res.ndf > 0) && (res.status == 0) && fit_ptr->IsValid() &&
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
        if (val <= min + kFitLimitTolerance || val >= max - kFitLimitTolerance) {
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
    logging::debug(
        task + ": chi2=" + std::to_string(best.chi2) + ", ndf=" + std::to_string(best.ndf) +
        ", chi2/ndf=" + std::to_string(best.chi2_ndf()) +
        ", p_value=" + std::to_string(best.p_value) + ", status=" + std::to_string(best.status) +
        ", covStatus=" + std::to_string(best.cov_status) +
        ", atLimit=" + (best.at_limit ? "true" : "false"));

    if (is_usable_fit(best) || best.attempts == 0 || !cfg.fit.retry_with_defaults) {
        return best;
    }

    logging::info(task + ": first attempt " +
                  (best.at_limit ? "hit parameter limits"
                                 : "is unusable (status=" + std::to_string(best.status) + ")") +
                  " — retrying with default initial parameters");

    auto altFit = std::unique_ptr<TF3>(make_cf_3d_fit(cfg, ch, centr, b, true));
    FitResult alt = fit_cf_3d(cf_hist, altFit.get(), cfg.fit);
    const int attempts = best.attempts + alt.attempts;
    best.attempts = attempts;
    alt.attempts = attempts;

    if (is_better_fit(alt, best)) {
        logging::info(task + ": retry improved the fit");
        return alt;
    }
    logging::info(task + ": retry did not improve — keeping first attempt");
    return best;
}
