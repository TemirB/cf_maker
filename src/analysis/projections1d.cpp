#include "analysis/projections1d.h"

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TH3.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TString.h>

#include "analysis/graphs.h"
#include "core/binning.h"
#include "core/correlation.h"
#include "core/fs.h"
#include "core/lcms.h"
#include "core/log.h"
#include "draw/draw.h"
#include "fit/model.h"
#include "io/input.h"
#include "io/output.h"

namespace
{

std::unique_ptr<TPaveText> get_fit_stats(const FitResult& r, Float_t textSize, double x1 = 0.55,
                                         double y1 = 0.70, double x2 = 0.88, double y2 = 0.88)
{
    auto stats = std::make_unique<TPaveText>(x1, y1, x2, y2, "NDC");
    stats->SetBorderSize(1);
    stats->SetFillColor(kWhite);
    stats->SetTextAlign(12);
    stats->SetTextFont(42);
    stats->SetTextSize(textSize);

    if (!is_usable_fit(r)) {
        stats->AddText("Fit unavailable");
        if (!r.is_finite()) {
            stats->AddText("Non-finite fit result");
        } else if (!r.ok) {
            stats->AddText(r.attempts == 0 ? "Missing fit result" : "Fit failed");
        } else if (r.at_limit) {
            stats->AddText("Parameter at limit");
        } else {
            stats->AddText("Nonpositive ndf");
        }
        stats->AddText(Form("status = %d, covariance = %d", r.status, r.cov_status));
        return stats;
    }
    stats->AddText(Form("#chi^{2}/ndf = %.3f", r.chi2_ndf()));
    stats->AddText(Form("R_{out}  = %.3f #pm %.3f fm", r.r[0], r.e_r[0]));
    stats->AddText(Form("R_{side} = %.3f #pm %.3f fm", r.r[1], r.e_r[1]));
    stats->AddText(Form("R_{long} = %.3f #pm %.3f fm", r.r[2], r.e_r[2]));
    stats->AddText(Form("R_{out-long} = %.3f #pm %.3f fm", r.r[4], r.e_r[4]));
    stats->AddText(Form("#lambda  = %.3f #pm %.3f", r.lambda, r.e_lambda));

    return stats;
}

std::unique_ptr<TH1D> build_lcms_fit_from_3d_weighted(const TH3D& denSource, const FitResult& r,
                                                      LCMSAxis axis, double sliceWidth,
                                                      const std::string& tag)
{
    char proj[2] = {projection_char(axis), '\0'};

    auto sliced =
        std::unique_ptr<TH3D>(static_cast<TH3D*>(denSource.Clone(("fit_slice_3d_" + tag).c_str())));
    slice_projection(*sliced, axis, sliceWidth);
    sliced->SetDirectory(nullptr);

    auto den =
        std::unique_ptr<TH3D>(static_cast<TH3D*>(sliced->Clone(("fit_den_3d_" + tag).c_str())));
    auto num =
        std::unique_ptr<TH3D>(static_cast<TH3D*>(sliced->Clone(("fit_num_3d_" + tag).c_str())));

    den->Reset("ICES");
    num->Reset("ICES");

    den->SetDirectory(nullptr);
    num->SetDirectory(nullptr);

    int x1 = sliced->GetXaxis()->GetFirst();
    int x2 = sliced->GetXaxis()->GetLast();
    int y1 = sliced->GetYaxis()->GetFirst();
    int y2 = sliced->GetYaxis()->GetLast();
    int z1 = sliced->GetZaxis()->GetFirst();
    int z2 = sliced->GetZaxis()->GetLast();

    for (int ix = x1; ix <= x2; ++ix) {
        double qOut = sliced->GetXaxis()->GetBinCenter(ix);

        for (int iy = y1; iy <= y2; ++iy) {
            double qSide = sliced->GetYaxis()->GetBinCenter(iy);

            for (int iz = z1; iz <= z2; ++iz) {
                double qLong = sliced->GetZaxis()->GetBinCenter(iz);

                double w = sliced->GetBinContent(ix, iy, iz);
                if (w <= 0.0 || !std::isfinite(w)) {
                    continue;
                }

                double fit_value = eval_cf_3d(r, qOut, qSide, qLong);
                if (!std::isfinite(fit_value)) {
                    continue;
                }

                den->SetBinContent(ix, iy, iz, w);
                num->SetBinContent(ix, iy, iz, w * fit_value);
            }
        }
    }

    auto proj_num = std::unique_ptr<TH1D>(static_cast<TH1D*>(num->Project3D(proj)));
    auto proj_den = std::unique_ptr<TH1D>(static_cast<TH1D*>(den->Project3D(proj)));
    auto hist_num = std::unique_ptr<TH1D>(static_cast<TH1D*>(proj_num->Clone()));
    auto hist_den = std::unique_ptr<TH1D>(static_cast<TH1D*>(proj_den->Clone()));

    hist_num->SetDirectory(nullptr);
    hist_den->SetDirectory(nullptr);

    auto fit = std::unique_ptr<TH1D>(
        static_cast<TH1D*>(hist_num->Clone(("fit_" + tag + "_" + axis_name(axis)).c_str())));
    fit->SetDirectory(nullptr);
    fit->Divide(hist_num.get(), hist_den.get());
    // This curve is evaluated at fixed fitted parameters and fixed projection
    // weights. ROOT must not infer counting errors from its bin contents.
    for (int bin = 0; bin < fit->GetNcells(); ++bin) {
        fit->SetBinError(bin, 0.0);
    }

    return fit;
}

std::pair<std::unique_ptr<TH1D>, std::unique_ptr<TH1D>>
create_1d(TH3D& den, TH3D& num, const FitResult& r, const LCMSAxis axis, const std::string& name,
          double sliceWidth, CorrelationStatistics statistics)
{
    std::string baseName = name + " " + axis_name(axis);

    std::unique_ptr<TH1D> hist_den(project_1d(den, axis, sliceWidth));
    std::unique_ptr<TH1D> hist_num(project_1d(num, axis, sliceWidth));
    auto cf = std::unique_ptr<TH1D>(static_cast<TH1D*>(hist_den->Clone(name.c_str())));
    cf->SetDirectory(nullptr);
    fill_correlation(*cf, *hist_num, *hist_den, statistics);

    std::string n = name + " " + axis_name(axis);
    cf->SetTitle(n.c_str());

    std::unique_ptr<TH1D> fit;
    if (is_usable_fit(r)) {
        fit = build_lcms_fit_from_3d_weighted(den, r, axis, sliceWidth, baseName);
    }

    return {std::move(cf), std::move(fit)};
}

void save_canvas_to_file(TFile* out, TH1D* cf, TH1D* fit, TPaveText* stats, draw::Style style,
                         std::string cf_name, LCMSAxis axis, const ProjectionConfig& projCfg)
{
    std::string name = cf_name + "_" + axis_name(axis);
    auto c = std::make_unique<TCanvas>(name.data(), name.data(), 800, 600);
    style_1d_cf(cf, name, axis_name(axis).data(), style);
    if (fit) {
        style_fit(fit, style);
        fit->GetXaxis()->SetRangeUser(-projCfg.axis_range_1d, projCfg.axis_range_1d);
    }

    cf->GetYaxis()->SetRangeUser(projCfg.cf_y_min_1d, projCfg.cf_y_max_1d);
    cf->GetXaxis()->SetRangeUser(-projCfg.axis_range_1d, projCfg.axis_range_1d);
    cf->SetStats(kFALSE);

    gStyle->SetOptFit(0102);

    c->cd();
    cf->Draw("P");
    if (fit) {
        fit->Draw("L SAME");
    }
    stats->Draw();

    write_output_object(*out, *c);
}

void draw_cf_and_fit(TCanvas* c, TH1D* cf, TH1D* fit, TPaveText* stats, int lcms)
{
    c->cd(lcms + 1);
    cf->Draw("P");
    if (fit) {
        fit->Draw("L SAME");
    }

    if (lcms == 2) {
        stats->Draw();
    }
}

void draw_cf_over_fit(TCanvas* c, TH1D* cf, TH1D* fit, TPaveText* stats, std::string name, int lcms,
                      draw::Style style, const ProjectionConfig& projCfg,
                      std::vector<std::unique_ptr<TH1D>>& keep_alive)
{
    auto axis = LCMSAxis(lcms);
    std::unique_ptr<TH1D> fit_over_cf(static_cast<TH1D*>(cf->Clone("fit_over_cf")));
    fit_over_cf->Divide(fit, cf);
    for (int bin = 0; bin < fit_over_cf->GetNcells(); ++bin) {
        const double value = cf->GetBinContent(bin);
        const double model = fit->GetBinContent(bin);
        // Conditional diagnostic: the fitted curve is held fixed.
        const double error =
            value != 0 ? std::abs(model / value) * (cf->GetBinError(bin) / std::abs(value)) : 0.0;
        fit_over_cf->SetBinError(bin, error);
    }

    std::string name_fit_over_cf = name + " " + axis_name(axis);
    style_1d_cf(fit_over_cf.get(), name_fit_over_cf.data(), axis_name(axis).data(), style);
    fit_over_cf->GetYaxis()->SetRangeUser(projCfg.fit_over_cf_y_min_1d,
                                          projCfg.fit_over_cf_y_max_1d);
    fit_over_cf->GetXaxis()->SetRangeUser(-projCfg.axis_range_1d, projCfg.axis_range_1d);
    fit_over_cf->SetStats(kFALSE);

    c->cd(lcms + 1);
    fit_over_cf->Draw("P");
    if (lcms == 2) {
        stats->Draw();
    }

    keep_alive.push_back(std::move(fit_over_cf));
}

} // namespace

void make_lcms_1d_projections(Config& cfg, TFile* in, TFile* out)
{
    const auto statistics = correlation_statistics(*in);
    draw::Style style = draw::default_style();

    const FitGrid& fitRes = cfg.fit_results;
    const Bin& bin = cfg.binning;

    std::string dep_dir = cfg.output.dir + "/dependency/fit_over_cf";
    std::string dir = cfg.output.dir + "/all_1d_histos";
    ensure_dir(dir);
    ensure_dir(dep_dir);
    const std::string format = cfg.general.images.format;

    logging::info("projections_1d: " + std::to_string(cfg.selection.charges.size()) +
                  " charges x " + std::to_string(cfg.selection.centralities.size()) +
                  " centralities x " + std::to_string(bin.count) + " bins");

    for (const int ch_idx : cfg.selection.charges)
        for (const int cent_idx : cfg.selection.centralities)
            for (int b = 0; b < bin.count; b++) {
                const FitResult r = fitRes[ch_idx][cent_idx][b];
                auto [den_raw, num_raw] = get_hists(in, ch_idx, cent_idx, b);
                std::unique_ptr<TH3D> den(den_raw), num(num_raw);
                if (!den || !num) {
                    throw std::runtime_error(
                        "projections_1d: missing input histograms for charge=" +
                        std::to_string(ch_idx) + ", centrality=" + std::to_string(cent_idx) +
                        ", bin=" + std::to_string(b));
                }
                validate_correlation_inputs(*num, *den, statistics);
                auto stats = get_fit_stats(r, 0.032f);

                std::string cf_name = get_cf_name(ch_idx, cent_idx, cfg.input.type, bin.names[b]);
                auto canvas = std::make_unique<TCanvas>(cf_name.data(), cf_name.data(), 2400, 800);
                canvas->SetTitle(cf_name.data());
                canvas->Divide(3, 1);
                fix_margin(canvas.get(), 3);

                std::string name_fit_over_cf = Form(
                    "fit/cf at charge=%s, centrality=%s, %s=%s", charge::kNames[ch_idx],
                    centrality::kNames[cent_idx], cfg.input.type.c_str(), bin.names[b].c_str());
                std::unique_ptr<TCanvas> c_fit_over_cf;
                if (is_usable_fit(r)) {
                    c_fit_over_cf = std::make_unique<TCanvas>(name_fit_over_cf.data(),
                                                              name_fit_over_cf.data(), 2400, 800);
                    c_fit_over_cf->SetTitle(name_fit_over_cf.data());
                    c_fit_over_cf->Divide(3, 1);
                }

                std::vector<std::unique_ptr<TH1D>> keep_alive;

                for (int lcms = 0; lcms < 3; lcms++) {
                    auto axis = LCMSAxis(lcms);

                    auto [cf, fit] = create_1d(*den, *num, r, axis, cf_name,
                                               cfg.projections.slice_1d, statistics);

                    if (fit) {
                        draw_cf_over_fit(c_fit_over_cf.get(), cf.get(), fit.get(), stats.get(),
                                         name_fit_over_cf, lcms, style, cfg.projections,
                                         keep_alive);
                    }
                    save_canvas_to_file(out, cf.get(), fit.get(), stats.get(), style, cf_name, axis,
                                        cfg.projections);
                    draw_cf_and_fit(canvas.get(), cf.get(), fit.get(), stats.get(), lcms);

                    keep_alive.push_back(std::move(cf));
                    keep_alive.push_back(std::move(fit));
                }

                if (cfg.general.images.need) {
                    auto image_name = Form("%s/cfs_%s_%s_%s.%s", dir.data(), charge::kNames[ch_idx],
                                           centrality::kNames[cent_idx], bin.file_names[b].c_str(),
                                           format.data());
                    save_canvas_quiet(canvas.get(), image_name);

                    auto cf_over_fit_name = Form(
                        "%s/fit_over_cf_%s_%s_%s.%s", dep_dir.data(), charge::kNames[ch_idx],
                        centrality::kNames[cent_idx], bin.file_names[b].c_str(), format.data());
                    if (c_fit_over_cf) {
                        save_canvas_quiet(c_fit_over_cf.get(), cf_over_fit_name);
                    }
                }
            }
}
