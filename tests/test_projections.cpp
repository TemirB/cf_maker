#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

#include <TCanvas.h>
#include <TH1.h>
#include <TH3.h>
#include <TList.h>
#include <TMemFile.h>
#include <TROOT.h>
#include <TPaveText.h>
#include <TText.h>

#include "analysis/projections1d.h"
#include "analysis/projections2d.h"
#include "core/lcms.h"
#include "io/input.h"

namespace
{
void check_slice_boundaries(TH3D& source)
{
    source.Sumw2();
    // Populate flow cells as well, so slices extending outside the axis have
    // a measurable contribution from underflow and overflow.
    for (int x = 0; x <= 9; ++x) {
        for (int y = 0; y <= 9; ++y) {
            for (int z = 0; z <= 9; ++z) {
                source.SetBinContent(x, y, z, 1);
                source.SetBinError(x, y, z, .5);
            }
        }
    }

    struct SliceCase
    {
        double width;
        int sliced_bins;
    };
    // Aligned cuts exclude bins touching only the boundary. Interior cuts
    // include every intersecting bin, even for a slice narrower than one bin.
    const std::array<SliceCase, 4> cases = {{{.05, 2}, {.075, 4}, {1e-20, 2}, {.25, 10}}};
    for (const auto& slice : cases) {
        for (const auto axis : {LCMSAxis::Out, LCMSAxis::Side, LCMSAxis::Long}) {
            std::unique_ptr<TH1D> projection(project_1d(source, axis, slice.width));
            const double cells = slice.sliced_bins * slice.sliced_bins;
            for (int bin = 1; bin <= 8; ++bin) {
                if (std::abs(projection->GetBinContent(bin) - cells) > 1e-12 ||
                    std::abs(projection->GetBinError(bin) - .5 * std::sqrt(cells)) > 1e-12) {
                    throw std::runtime_error("Wrong 1D projection boundary selection for " +
                                             std::string(source.GetName()) + " " + axis_name(axis));
                }
            }
        }
        std::unique_ptr<TH2D> projection(
            project_2d(source, LCMSAxis::Out, LCMSAxis::Long, slice.width));
        for (int x = 1; x <= 8; ++x) {
            for (int y = 1; y <= 8; ++y) {
                if (std::abs(projection->GetBinContent(x, y) - slice.sliced_bins) > 1e-12 ||
                    std::abs(projection->GetBinError(x, y) - .5 * std::sqrt(slice.sliced_bins)) >
                        1e-12) {
                    throw std::runtime_error("Wrong 2D projection boundary selection for " +
                                             std::string(source.GetName()));
                }
            }
        }
    }
}

void check_invalid_slice_width(TH3D& source)
{
    for (const double invalid : {0., -.05, std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::infinity()}) {
        bool rejected_1d = false;
        bool rejected_2d = false;
        try {
            std::unique_ptr<TH1D> projection(project_1d(source, LCMSAxis::Out, invalid));
        } catch (const std::invalid_argument&) {
            rejected_1d = true;
        }
        try {
            std::unique_ptr<TH2D> projection(
                project_2d(source, LCMSAxis::Out, LCMSAxis::Long, invalid));
        } catch (const std::invalid_argument&) {
            rejected_2d = true;
        }
        if (!rejected_1d || !rejected_2d) {
            throw std::runtime_error("Invalid projection slice width was accepted");
        }
    }
}

void check_slice_boundary_roundoff()
{
    TH3D uniform("uniform_slice_boundaries", "", 8, -.2, .2, 8, -.2, .2, 8, -.2, .2);
    check_slice_boundaries(uniform);

    // Force the same one-ULP mismatch that can occur for uniform bin edges
    // when ROOT is built without fused multiply-add on another platform.
    const std::array<double, 9> edges = {
        -.2, -.15, -.1, std::nextafter(-.05, 0.), 0., std::nextafter(.05, 0.), .1, .15, .2};
    TH3D variable("variable_slice_boundaries", "", 8, edges.data(), 8, edges.data(), 8,
                  edges.data());
    check_slice_boundaries(variable);
    check_invalid_slice_width(uniform);
    check_invalid_slice_width(variable);

    // ROOT Project3D misinterprets an underflow-only range as all regular bins.
    TH3D shifted("underflow_only_slice", "", 8, -.2, .2, 8, .3, .8, 8, -.2, .2);
    bool rejected_1d = false;
    bool rejected_2d = false;
    try {
        std::unique_ptr<TH1D> projection(project_1d(shifted, LCMSAxis::Out, .1));
    } catch (const std::invalid_argument&) {
        rejected_1d = true;
    }
    try {
        std::unique_ptr<TH2D> projection(project_2d(shifted, LCMSAxis::Out, LCMSAxis::Long, .1));
    } catch (const std::invalid_argument&) {
        rejected_2d = true;
    }
    if (!rejected_1d || !rejected_2d) {
        throw std::runtime_error("Unsupported underflow-only slice was accepted");
    }
}

void check_2d_axis_orientation()
{
    struct AxisSpec
    {
        LCMSAxis axis;
        int bins;
        double minimum;
        double maximum;
        int slice_first;
        int slice_last;
    };
    const std::array<AxisSpec, 3> axes = {{{LCMSAxis::Out, 6, -0.3, 0.3, 2, 5},
                                           {LCMSAxis::Side, 8, -0.8, 0.8, 4, 5},
                                           {LCMSAxis::Long, 10, -1.5, 1.5, 5, 6}}};
    TH3D source("axis_orientation", "", 6, -0.3, 0.3, 8, -0.8, 0.8, 10, -1.5, 1.5);
    source.Sumw2();
    for (int out = 1; out <= 6; ++out) {
        for (int side = 1; side <= 8; ++side) {
            for (int longitudinal = 1; longitudinal <= 10; ++longitudinal) {
                source.SetBinContent(out, side, longitudinal,
                                     10000.0 * out + 100.0 * side + longitudinal);
                source.SetBinError(out, side, longitudinal, out + 0.1 * side + 0.01 * longitudinal);
            }
        }
    }

    // Check both orders for every pair: X is always the first requested axis.
    for (std::size_t first = 0; first < axes.size(); ++first) {
        for (std::size_t second = 0; second < axes.size(); ++second) {
            if (first == second) {
                continue;
            }
            const std::size_t frozen = 3 - first - second;
            std::unique_ptr<TH2D> projection(
                project_2d(source, axes[first].axis, axes[second].axis, 0.16));
            const auto check_axis = [](const TAxis& actual, const AxisSpec& expected) {
                if (actual.GetNbins() != expected.bins ||
                    std::abs(actual.GetXmin() - expected.minimum) > 1e-12 ||
                    std::abs(actual.GetXmax() - expected.maximum) > 1e-12) {
                    throw std::runtime_error("Wrong 2D projection axis geometry");
                }
            };
            check_axis(*projection->GetXaxis(), axes[first]);
            check_axis(*projection->GetYaxis(), axes[second]);
            for (int x = 1; x <= axes[first].bins; ++x) {
                for (int y = 1; y <= axes[second].bins; ++y) {
                    std::array<int, 3> bins = {};
                    bins[first] = x;
                    bins[second] = y;
                    double content = 0.0;
                    double variance = 0.0;
                    for (int sliced = axes[frozen].slice_first; sliced <= axes[frozen].slice_last;
                         ++sliced) {
                        bins[frozen] = sliced;
                        content += source.GetBinContent(bins[0], bins[1], bins[2]);
                        const double error = source.GetBinError(bins[0], bins[1], bins[2]);
                        variance += error * error;
                    }
                    if (std::abs(projection->GetBinContent(x, y) - content) > 1e-10 ||
                        std::abs(projection->GetBinError(x, y) - std::sqrt(variance)) > 1e-10) {
                        throw std::runtime_error("Wrong 2D projection contents or slice");
                    }
                }
            }
        }
    }
}

void check_canvas(TMemFile& output, const std::string& name, int dimension, int count)
{
    auto* canvas = dynamic_cast<TCanvas*>(output.Get(name.c_str()));
    if (!canvas) {
        throw std::runtime_error("Missing canvas: " + name);
    }
    int histograms = 0;
    for (auto* object : *canvas->GetListOfPrimitives()) {
        auto* histogram = dynamic_cast<TH1*>(object);
        if (!histogram) {
            continue;
        }
        ++histograms;
        if (histogram->GetDimension() != dimension) {
            throw std::runtime_error("Wrong projection dimension");
        }
        for (int x = 1; x <= histogram->GetNbinsX(); ++x) {
            for (int y = 1; y <= histogram->GetNbinsY(); ++y) {
                if (std::abs(histogram->GetBinContent(histogram->GetBin(x, y)) - 1.5) > 1e-10) {
                    throw std::runtime_error("Wrong projected CF or weighted fit value");
                }
            }
        }
    }
    if (histograms != count) {
        throw std::runtime_error("Wrong number of plotted histograms");
    }
}

void check_unavailable_fit(TCanvas& canvas, const std::string& reason)
{
    bool unavailable = false;
    bool has_reason = false;
    bool has_status = false;
    for (auto* object : *canvas.GetListOfPrimitives()) {
        const auto* stats = dynamic_cast<TPaveText*>(object);
        if (!stats) {
            continue;
        }
        for (auto* line : *stats->GetListOfLines()) {
            const auto* text = dynamic_cast<TText*>(line);
            if (!text) {
                continue;
            }
            const std::string title = text->GetTitle();
            unavailable |= title == "Fit unavailable";
            has_reason |= title == reason;
            has_status |= title.find("status = ") != std::string::npos;
            if (title.find("#chi^{2}/ndf") != std::string::npos) {
                throw std::runtime_error("Unavailable fit still has chi2/ndf statistics");
            }
        }
    }
    if (!unavailable || !has_reason || !has_status) {
        throw std::runtime_error("Data-only projection lacks unavailable-fit diagnostics");
    }
}

void check_data_only_projections(Config cfg, TMemFile& input)
{
    const FitResult successful = cfg.fit_results[0][0][0];
    std::array<FitResult, 5> unavailable;
    unavailable[0] = FitResult{};
    unavailable[1] = successful;
    unavailable[1].ok = false;
    unavailable[1].status = 4;
    unavailable[1].attempts = 1;
    unavailable[2] = successful;
    unavailable[2].r[4] = std::numeric_limits<double>::quiet_NaN();
    unavailable[3] = successful;
    unavailable[3].at_limit = true;
    unavailable[4] = successful;
    unavailable[4].ndf = 0;
    const std::array<const char*, 5> reasons = {"Missing fit result", "Fit failed",
                                                "Non-finite fit result", "Parameter at limit",
                                                "Nonpositive ndf"};
    const std::string name = get_cf_name(0, 0, "kt", "test");
    for (std::size_t state = 0; state < unavailable.size(); ++state) {
        cfg.fit_results[0][0][0] = unavailable[state];
        TMemFile output(("data_only_" + std::to_string(state) + ".root").c_str(), "RECREATE");
        make_lcms_1d_projections(cfg, &input, &output);
        for (const char* axis : {"out", "side", "long"}) {
            const std::string canvas_name = name + "_" + axis;
            check_canvas(output, canvas_name, 1, 1);
            auto* canvas = dynamic_cast<TCanvas*>(output.Get(canvas_name.c_str()));
            check_unavailable_fit(*canvas, reasons[state]);
        }
    }
}
} // namespace

int main()
{
    gROOT->SetBatch(true);
    TH1::AddDirectory(false);
    check_slice_boundary_roundoff();
    check_2d_axis_orientation();
    Config cfg;
    cfg.input.type = "kt";
    cfg.output.dir = "projection_test_output";
    cfg.general.images.need = false;
    cfg.selection.charges = {0};
    cfg.selection.centralities = {0};
    cfg.binning.count = 1;
    cfg.binning.names = {"test"};
    cfg.binning.file_names = {"test"};
    cfg.fit_results = FitGrid(1, std::vector<std::vector<FitResult>>(1, std::vector<FitResult>(1)));
    cfg.fit_results[0][0][0].lambda = 0.5;
    cfg.fit_results[0][0][0].ndf = 1;
    cfg.fit_results[0][0][0].ok = true;
    cfg.fit_results[0][0][0].status = 0;
    cfg.fit_results[0][0][0].cov_status = 3;
    cfg.fit_results[0][0][0].attempts = 1;

    TMemFile input("input.root", "RECREATE");
    TH3D den("bp_0_0_num_0", "", 8, -0.2, 0.2, 8, -0.2, 0.2, 8, -0.2, 0.2);
    TH3D num("bp_0_0_num_wei_0", "", 8, -0.2, 0.2, 8, -0.2, 0.2, 8, -0.2, 0.2);
    den.Sumw2();
    num.Sumw2();
    for (int x = 1; x <= 8; ++x) {
        for (int y = 1; y <= 8; ++y) {
            for (int z = 1; z <= 8; ++z) {
                double weight = x + 2.0 * y + 3.0 * z;
                den.SetBinContent(x, y, z, weight);
                den.SetBinError(x, y, z, std::sqrt(weight));
                num.SetBinContent(x, y, z, 1.5 * weight);
                num.SetBinError(x, y, z, std::sqrt(2.45 * weight));
            }
        }
    }
    den.Write();
    num.Write();

    TMemFile output("output.root", "RECREATE");
    make_lcms_1d_projections(cfg, &input, &output);
    make_lcms_2d_projections(cfg, &input, &output);
    const auto name = get_cf_name(0, 0, "kt", "test");
    for (const auto* axis : {"out", "side", "long"}) {
        check_canvas(output, name + "_" + axis, 1, 2);
    }
    for (const auto* axes : {"out-side", "out-long", "side-long"}) {
        check_canvas(output, name + " " + axes, 2, 1);
    }
    check_data_only_projections(cfg, input);
    std::filesystem::remove_all(cfg.output.dir);
    std::cout << "All projection tests passed\n";
}
