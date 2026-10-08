#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include <Math/MinimizerOptions.h>
#include <TFile.h>
#include <TROOT.h>

#include "config/config.h"
#include "core/correlation.h"
#include "core/lcms.h"
#include "fit/model.h"
#include "io/input.h"

namespace
{
std::ofstream csv(const std::filesystem::path& directory, const std::string& name,
                  const std::string& header)
{
    std::ofstream file(directory / name);
    file.exceptions(std::ios::failbit | std::ios::badbit);
    file << std::setprecision(17) << header << '\n';
    return file;
}

void set_parameters(TF3& model, const FitResult& result)
{
    for (int i = 0; i < 6; ++i) {
        model.SetParameter(i, i < 3 ? result.r[i] * result.r[i] : result.r[i]);
    }
    model.SetParameter(6, result.lambda);
}

void run_config(const std::string& config_path, const std::filesystem::path& directory)
{
    auto cfg = load(config_path);
    build(cfg);
    std::unique_ptr<TFile> input(TFile::Open(cfg.input.file.c_str(), "READ"));
    if (!input || input->IsZombie()) {
        throw std::runtime_error("cannot read input: " + cfg.input.file);
    }
    const auto statistics = correlation_statistics(*input);
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer(cfg.fit.minimizer.c_str());
    // Per-cell hand calculation below deliberately supports the default
    // center-evaluated, diagonal chi-square objective only.
    if (cfg.fit.use_integral || cfg.fit.options.find_first_of("ILPWU") != std::string::npos ||
        cfg.fit.options.find('R') == std::string::npos) {
        throw std::runtime_error(
            "audit requires center-evaluated chi-square with R, without I/L/P/W/U");
    }
    std::filesystem::create_directory(directory);
    std::filesystem::copy_file(config_path, directory / "config.json");
    auto summary = csv(directory, "fits.csv",
                       "charge,centrality,bin,variant,chi2,ndf,chi2_ndf,status,cov_status,at_limit,"
                       "r_out,r_side,r_long,lambda,e_r_out,e_r_side,e_r_long,e_lambda");
    auto cells = csv(directory, "cells.csv",
                     "charge,centrality,bin,cell,q_out,q_side,q_long,N,S,Q,C,sigma_program,sigma_"
                     "raw,sigma_B,sigma_independent,residual");
    auto rounding = csv(directory, "rounding.csv",
                        "charge,centrality,bin,cell,digits,C_original,C_rounded,delta_C,residual_"
                        "original,residual_rounded,sigma_raw_original,sigma_raw_rounded");
    auto pulls = csv(directory, "chi2_cells.csv",
                     "charge,centrality,bin,variant,cell,C,model,sigma,chi2_contribution");
    auto check = csv(directory, "chi2_check.csv",
                     "charge,centrality,bin,variant,manual_chi2,program_chi2,difference,used_cells,"
                     "free_parameters,manual_ndf,program_ndf");
    auto ratios = csv(directory, "fit_over_cf.csv",
                      "charge,centrality,bin,axis,cell,q,CF,sigma_CF,fit,fit_over_CF,sigma_old,"
                      "sigma_new,delta_sigma");
    auto errors = csv(directory, "errors.csv", "charge,centrality,bin,stage");
    std::ofstream notes(directory / "README.txt");
    notes.exceptions(std::ios::failbit | std::ios::badbit);
    notes << "Input: " << cfg.input.file << "\nROOT: " << gROOT->GetVersion() << "\nStatistics: "
          << (statistics == CorrelationStatistics::PairWeights ? "pair_weights" : "fixed_reference")
          << "\nB and independent-error fits are sensitivity comparisons, not endorsed statistical "
             "models.\n"
             "Rounding changes only stored S, keeps N and Q fixed; it cannot reconstruct lost "
             "float accumulation precision.\n"
             "sigma_raw is before the production roundoff guard. Negative residuals have "
             "unavailable sigma (nan).\n"
             "chi2_cells uses bin centers, selected range, positive errors, and counts free "
             "parameters.\n"
             "Fit/CF holds fitted parameters and projection weights fixed; parameter covariance is "
             "not included.\n"
             "Failed/boundary fits are retained in fits.csv; fit/CF uses only usable production "
             "fits.\n";
    for (int charge : cfg.selection.charges) {
        for (int centrality : cfg.selection.centralities) {
            for (int bin = 0; bin < cfg.binning.count; ++bin) {
                const std::string prefix = std::to_string(charge) + "," +
                                           std::to_string(centrality) + "," + std::to_string(bin) +
                                           ",";
                std::cout << cfg.input.type << " " << prefix << std::endl;
                auto [count_ptr, sum_ptr] = get_hists(input.get(), charge, centrality, bin);
                std::unique_ptr<TH3D> count(count_ptr), sum(sum_ptr);
                if (!count || !sum)
                    throw std::runtime_error("missing input histogram: " + prefix);
                TH3D production(*sum), binomial(*sum), independent(*sum);
                production.SetDirectory(nullptr);
                binomial.SetDirectory(nullptr);
                independent.SetDirectory(nullptr);
                binomial.Divide(sum.get(), count.get(), 1, 1, "B");
                independent.Divide(sum.get(), count.get());
                bool valid = true;
                try {
                    fill_correlation(production, *sum, *count, statistics);
                } catch (const std::exception& error) {
                    valid = false;
                    errors << prefix << "production_moments\n";
                    notes << prefix << " " << error.what() << '\n';
                }
                for (int x = 1; x <= sum->GetNbinsX(); ++x) {
                    for (int y = 1; y <= sum->GetNbinsY(); ++y) {
                        for (int z = 1; z <= sum->GetNbinsZ(); ++z) {
                            const int cell = sum->GetBin(x, y, z);
                            const double n = count->GetBinContent(cell),
                                         s = sum->GetBinContent(cell);
                            if (n <= 0)
                                continue;
                            const double q = sum->GetSumw2N() ? sum->GetSumw2()->At(cell) : NAN;
                            const long double residual = q - static_cast<long double>(s) * s / n;
                            const auto raw_error = [n](long double r) {
                                return n > 1 && r >= 0
                                           ? std::sqrt(static_cast<double>(r / n / (n - 1)))
                                           : NAN;
                            };
                            cells << prefix << cell << ',' << sum->GetXaxis()->GetBinCenter(x)
                                  << ',' << sum->GetYaxis()->GetBinCenter(y) << ','
                                  << sum->GetZaxis()->GetBinCenter(z) << ',' << n << ',' << s << ','
                                  << q << ',' << s / n << ','
                                  << (valid ? production.GetBinError(cell) : NAN) << ','
                                  << raw_error(residual) << ',' << binomial.GetBinError(cell) << ','
                                  << independent.GetBinError(cell) << ',' << residual << '\n';
                            for (int digits : {2, 4, 6, 8}) {
                                const double scale = std::pow(10., digits);
                                const double rounded = std::round(s * scale) / scale;
                                const long double changed =
                                    q - static_cast<long double>(rounded) * rounded / n;
                                rounding << prefix << cell << ',' << digits << ',' << s / n << ','
                                         << rounded / n << ',' << (rounded - s) / n << ','
                                         << residual << ',' << changed << ',' << raw_error(residual)
                                         << ',' << raw_error(changed) << '\n';
                            }
                        }
                    }
                }
                FitResult production_result;
                for (const auto& variant : {std::pair<const char*, TH3D*>{"program", &production},
                                            {"B", &binomial},
                                            {"independent", &independent}}) {
                    if (!valid && variant.second == &production)
                        continue;
                    const auto result =
                        fit_cf_3d_with_retry(variant.second, cfg, charge, centrality, bin);
                    if (variant.second == &production)
                        production_result = result;
                    summary << prefix << variant.first << ',' << result.chi2 << ',' << result.ndf
                            << ',' << result.chi2_ndf() << ',' << result.status << ','
                            << result.cov_status << ',' << result.at_limit << ',' << result.r[0]
                            << ',' << result.r[1] << ',' << result.r[2] << ',' << result.lambda
                            << ',' << result.e_r[0] << ',' << result.e_r[1] << ',' << result.e_r[2]
                            << ',' << result.e_lambda << '\n';
                    if (!result.is_finite() || result.attempts == 0)
                        continue;
                    std::unique_ptr<TF3> model(create_cf_3d_fit(cfg, charge, centrality, bin));
                    set_parameters(*model, result);
                    int free_parameters = 7;
                    for (const auto& frozen : cfg.fit.freeze)
                        if (frozen)
                            --free_parameters;
                    long double total = 0;
                    int used = 0;
                    for (int x = 1; x <= variant.second->GetNbinsX(); ++x) {
                        const double qx = variant.second->GetXaxis()->GetBinCenter(x);
                        for (int y = 1; y <= variant.second->GetNbinsY(); ++y) {
                            const double qy = variant.second->GetYaxis()->GetBinCenter(y);
                            for (int z = 1; z <= variant.second->GetNbinsZ(); ++z) {
                                const double qz = variant.second->GetZaxis()->GetBinCenter(z);
                                if (std::abs(qx) > cfg.fit.q_max || std::abs(qy) > cfg.fit.q_max ||
                                    std::abs(qz) > cfg.fit.q_max)
                                    continue;
                                const int cell = variant.second->GetBin(x, y, z);
                                const double error = variant.second->GetBinError(cell);
                                if (error <= 0 || !std::isfinite(error))
                                    continue;
                                const double value = variant.second->GetBinContent(cell);
                                const double expected = model->Eval(qx, qy, qz);
                                const long double pull =
                                    (static_cast<long double>(value) - expected) / error;
                                total += pull * pull;
                                ++used;
                                pulls << prefix << variant.first << ',' << cell << ',' << value
                                      << ',' << expected << ',' << error << ',' << pull * pull
                                      << '\n';
                            }
                        }
                    }
                    check << prefix << variant.first << ',' << total << ',' << result.chi2 << ','
                          << total - result.chi2 << ',' << used << ',' << free_parameters << ','
                          << std::max(0, used - free_parameters) << ',' << result.ndf << '\n';
                }
                if (!valid || !is_usable_fit(production_result))
                    continue;
                // Reconstruct the pre-change weighted model projection, including
                // ROOT's default error handling, then compare with the explicit formula.
                TH3D model_sum(*count), model_count(*count);
                model_sum.SetDirectory(nullptr);
                model_count.SetDirectory(nullptr);
                model_sum.Reset("ICES");
                model_count.Reset("ICES");
                for (int x = 1; x <= count->GetNbinsX(); ++x)
                    for (int y = 1; y <= count->GetNbinsY(); ++y)
                        for (int z = 1; z <= count->GetNbinsZ(); ++z) {
                            const int cell = count->GetBin(x, y, z);
                            const double n = count->GetBinContent(cell);
                            if (n <= 0)
                                continue;
                            const double value =
                                eval_cf_3d(production_result, count->GetXaxis()->GetBinCenter(x),
                                           count->GetYaxis()->GetBinCenter(y),
                                           count->GetZaxis()->GetBinCenter(z));
                            model_count.SetBinContent(cell, n);
                            model_sum.SetBinContent(cell, n * value);
                        }
                for (auto axis : {LCMSAxis::Out, LCMSAxis::Side, LCMSAxis::Long}) {
                    std::unique_ptr<TH1D> ps(project_1d(*sum, axis, cfg.projections.slice_1d)),
                        pn(project_1d(*count, axis, cfg.projections.slice_1d));
                    TH1D cf(*ps);
                    fill_correlation(cf, *ps, *pn, statistics);
                    std::unique_ptr<TH1D> ms(project_1d(model_sum, axis, cfg.projections.slice_1d)),
                        mn(project_1d(model_count, axis, cfg.projections.slice_1d));
                    TH1D model(*ms), old_ratio(cf);
                    model.Divide(ms.get(), mn.get());
                    old_ratio.Divide(&model, &cf);
                    for (int cell = 1; cell <= cf.GetNbinsX(); ++cell) {
                        const double value = cf.GetBinContent(cell);
                        if (value == 0)
                            continue;
                        const double fit = model.GetBinContent(cell),
                                     error = std::abs(fit / value) * cf.GetBinError(cell) /
                                             std::abs(value);
                        ratios << prefix << axis_name(axis) << ',' << cell << ','
                               << cf.GetBinCenter(cell) << ',' << value << ','
                               << cf.GetBinError(cell) << ',' << fit << ',' << fit / value << ','
                               << old_ratio.GetBinError(cell) << ',' << error << ','
                               << error - old_ratio.GetBinError(cell) << '\n';
                    }
                }
            }
        }
    }
    notes << "\nCompleted. Check errors.csv and fit statuses before interpreting comparisons.\n";
}
} // namespace

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "Usage: statistics_diff NEW_OUTPUT_DIR CONFIG ...\n";
        return 1;
    }
    gROOT->SetBatch(true);
    TH1::AddDirectory(false);
    try {
        const std::filesystem::path output(argv[1]);
        if (std::filesystem::exists(output))
            throw std::runtime_error("output already exists: " + output.string());
        std::filesystem::create_directories(output);
        bool failed = false;
        for (int i = 2; i < argc; ++i) {
            try {
                run_config(argv[i], output / (std::to_string(i - 1) + "_" +
                                              std::filesystem::path(argv[i]).stem().string()));
            } catch (const std::exception& error) {
                failed = true;
                std::cerr << argv[i] << ": " << error.what() << '\n';
            }
        }
        std::cout << "Results: " << std::filesystem::absolute(output) << '\n';
        return failed ? 2 : 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
