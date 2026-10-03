#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>

#include <TCanvas.h>
#include <TH1.h>
#include <TH3.h>
#include <TList.h>
#include <TMemFile.h>
#include <TROOT.h>

#include "analysis/projections1d.h"
#include "analysis/projections2d.h"
#include "core/lcms.h"
#include "io/input.h"

namespace
{
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
} // namespace

int main()
{
    gROOT->SetBatch(true);
    TH1::AddDirectory(false);
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
    std::filesystem::remove_all(cfg.output.dir);
    std::cout << "All projection tests passed\n";
}
