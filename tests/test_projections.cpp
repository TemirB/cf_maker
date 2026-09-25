#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

#include <TCanvas.h>
#include <TH1.h>
#include <TH3.h>
#include <TList.h>
#include <TMemFile.h>
#include <TROOT.h>

#include "analysis/projections1d.h"
#include "analysis/projections2d.h"
#include "io/input.h"

namespace
{
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
    for (int x = 1; x <= 8; ++x) {
        for (int y = 1; y <= 8; ++y) {
            for (int z = 1; z <= 8; ++z) {
                double weight = x + 2.0 * y + 3.0 * z;
                den.SetBinContent(x, y, z, weight);
                num.SetBinContent(x, y, z, 1.5 * weight);
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
