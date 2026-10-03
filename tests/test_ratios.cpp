#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#include <TFile.h>
#include <TH1.h>
#include <TH3.h>
#include <TMemFile.h>
#include <TROOT.h>

#include "analysis/ratio.h"
#include "core/correlation.h"
#include "io/input.h"

namespace
{
struct TempInput
{
    std::filesystem::path directory;
    std::string file;

    TempInput()
    {
        const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        directory =
            std::filesystem::temp_directory_path() / ("cf-maker-ratios-" + std::to_string(stamp));
        std::filesystem::create_directory(directory);
        file = (directory / "raw.root").string();
    }

    ~TempInput()
    {
        std::error_code error;
        std::filesystem::remove_all(directory, error);
    }
};

Config ratio_config(const std::string& input, double slice = .05)
{
    Config cfg;
    cfg.input.file = input;
    cfg.input.type = "kt";
    cfg.selection.charges = {0, 1};
    cfg.selection.centralities = {0};
    cfg.binning.count = 1;
    cfg.binning.names = {"test"};
    cfg.projections.slice_ratio = slice;
    return cfg;
}

void require_close(double actual, double expected)
{
    if (!std::isfinite(actual) || std::abs(actual - expected) > 1e-12) {
        throw std::runtime_error("Unexpected charge-ratio projection");
    }
}

TH1& result(TFile& output, const std::string& name)
{
    auto* histogram = dynamic_cast<TH1*>(output.Get(name.c_str()));
    if (!histogram) {
        throw std::runtime_error("Missing charge-ratio projection: " + name);
    }
    return *histogram;
}

void set_cell(TH3D& count, TH3D& sum, int x, int y, int z, double pairs, double mean,
              double mean_variance = 0)
{
    count.SetBinContent(x, y, z, pairs);
    count.SetBinError(x, y, z, std::sqrt(pairs));
    sum.SetBinContent(x, y, z, pairs * mean);
    sum.SetBinError(x, y, z, std::sqrt(pairs * mean * mean + pairs * (pairs - 1) * mean_variance));
}

void write_charge(TFile& raw, TMemFile& cf, int charge, TH3D& count, TH3D& sum)
{
    raw.cd();
    count.Write(("bp_" + std::to_string(charge) + "_0_num_0").c_str());
    sum.Write(("bp_" + std::to_string(charge) + "_0_num_wei_0").c_str());
    TH3D correlation(count);
    correlation.SetDirectory(nullptr);
    fill_correlation(correlation, sum, count);
    cf.cd();
    correlation.Write(get_cf_name(charge, 0, "kt", "test").c_str());
}

void check_ratio(double slice, bool sparse)
{
    TempInput temporary;
    Config cfg = ratio_config(temporary.file, slice);
    TMemFile input("charge_cf.root", "RECREATE");
    {
        TFile raw(temporary.file.c_str(), "RECREATE");
        for (int charge = 0; charge < 2; ++charge) {
            TH3D count("count", "", 8, -.2, .2, 8, -.2, .2, 8, -.2, .2);
            TH3D sum("sum", "", 8, -.2, .2, 8, -.2, .2, 8, -.2, .2);
            count.Sumw2();
            sum.Sumw2();
            for (int x = 1; x <= 8; ++x) {
                for (int y = 1; y <= 8; ++y) {
                    for (int z = 1; z <= 8; ++z) {
                        if (!sparse || (y == 4 && z == 4)) {
                            set_cell(count, sum, x, y, z, 101, 1, charge == 0 ? .01 : .04);
                        }
                    }
                }
            }
            write_charge(raw, input, charge, count, sum);
        }
    }
    TMemFile ratio_project("ratio_project.root", "RECREATE");
    TMemFile project_ratio("project_ratio.root", "RECREATE");
    do_cf_ratios(cfg, &input, &ratio_project, &project_ratio);
    auto& mean = result(project_ratio, "proj_of_ratios_0_0_out");
    auto& weighted = result(ratio_project, "ratio_proj_0_0_out");
    require_close(mean.GetBinContent(4), 1);
    require_close(weighted.GetBinContent(4), 1);
    // For width .05 the slice has four cells; empty cells carry no weight.
    if (slice == .05) {
        const double cells = sparse ? 1 : 4;
        require_close(mean.GetBinError(4), std::sqrt(.05 / cells));
        require_close(weighted.GetBinError(4), std::sqrt(5 / (101 * cells - 1)));
    }
    if (std::string(mean.GetYaxis()->GetTitle()).find("valid cells") == std::string::npos ||
        std::string(weighted.GetYaxis()->GetTitle()).find("pairs") == std::string::npos) {
        throw std::runtime_error("Charge-ratio labels must distinguish cells and pair weights");
    }
}

void check_unequal_occupancy()
{
    TempInput temporary;
    Config cfg = ratio_config(temporary.file);
    TMemFile input("unequal_cf.root", "RECREATE");
    {
        TFile raw(temporary.file.c_str(), "RECREATE");
        for (int charge = 0; charge < 2; ++charge) {
            TH3D count("count", "", 8, -.2, .2, 8, -.2, .2, 8, -.2, .2);
            TH3D sum("sum", "", 8, -.2, .2, 8, -.2, .2, 8, -.2, .2);
            count.Sumw2();
            sum.Sumw2();
            set_cell(count, sum, 4, 4, 4, charge == 0 ? 9 : 5, 1);
            set_cell(count, sum, 4, 5, 4, charge == 0 ? 1 : 5, charge == 0 ? 2 : 1);
            write_charge(raw, input, charge, count, sum);
        }
    }
    TMemFile ratio_project("unequal_ratio_project.root", "RECREATE");
    TMemFile project_ratio("unequal_project_ratio.root", "RECREATE");
    do_cf_ratios(cfg, &input, &ratio_project, &project_ratio);
    require_close(result(ratio_project, "ratio_proj_0_0_out").GetBinContent(4), 1.1);
    require_close(result(ratio_project, "ratio_proj_0_0_out").GetBinError(4), .1);
    require_close(result(project_ratio, "proj_of_ratios_0_0_out").GetBinContent(4), 1.5);
}

void check_zero_weights()
{
    TempInput temporary;
    Config cfg = ratio_config(temporary.file);
    TMemFile input("zero_cf.root", "RECREATE");
    {
        TFile raw(temporary.file.c_str(), "RECREATE");
        for (int charge = 0; charge < 2; ++charge) {
            TH3D count("count", "", 8, -.2, .2, 8, -.2, .2, 8, -.2, .2);
            TH3D sum("sum", "", 8, -.2, .2, 8, -.2, .2, 8, -.2, .2);
            count.Sumw2();
            sum.Sumw2();
            // Zero positive weight is a valid observation, not an empty cell.
            set_cell(count, sum, 4, 4, 4, 2, charge == 0 ? 0 : 1);
            set_cell(count, sum, 4, 5, 4, 2, charge == 0 ? 2 : 1);
            // A zero negative CF makes this charge ratio undefined.
            set_cell(count, sum, 5, 4, 4, 2, charge == 0 ? 3 : 0);
            write_charge(raw, input, charge, count, sum);
        }
    }
    TMemFile ratio_project("zero_ratio_project.root", "RECREATE");
    TMemFile project_ratio("zero_project_ratio.root", "RECREATE");
    do_cf_ratios(cfg, &input, &ratio_project, &project_ratio);
    for (auto* output : {&ratio_project, &project_ratio}) {
        const std::string name =
            output == &ratio_project ? "ratio_proj_0_0_out" : "proj_of_ratios_0_0_out";
        auto& histogram = result(*output, name);
        require_close(histogram.GetBinContent(4), 1);
        require_close(histogram.GetBinContent(5), 0);
        require_close(histogram.GetBinError(5), 0);
    }
}

template <typename Action> void require_failure(Action action, const std::string& message)
{
    try {
        action();
    } catch (const std::exception& error) {
        if (std::string(error.what()).find(message) != std::string::npos) {
            return;
        }
        throw;
    }
    throw std::runtime_error("Expected charge-ratio failure: " + message);
}

void check_missing_inputs()
{
    TempInput temporary;
    Config cfg = ratio_config(temporary.file);
    {
        TFile raw(temporary.file.c_str(), "RECREATE");
    }
    TMemFile input("missing_raw_cf.root", "RECREATE");
    TH3D histogram("cf", "", 8, -.2, .2, 8, -.2, .2, 8, -.2, .2);
    input.cd();
    for (int charge = 0; charge < 2; ++charge) {
        histogram.Write(get_cf_name(charge, 0, "kt", "test").c_str());
    }
    TMemFile ratio_project("missing_ratio_project.root", "RECREATE");
    TMemFile project_ratio("missing_project_ratio.root", "RECREATE");
    cfg.selection.charges = {0};
    require_failure([&] { do_cf_ratios(cfg, &input, &ratio_project, &project_ratio); },
                    "both selected charges");
    cfg.selection.charges = {0, 1};
    require_failure([&] { do_cf_ratios(cfg, &input, &ratio_project, &project_ratio); },
                    "missing raw input");
    TMemFile missing_cf("missing_cf.root", "RECREATE");
    require_failure([&] { do_cf_ratios(cfg, &missing_cf, &ratio_project, &project_ratio); },
                    "missing charge CF");
    if (ratio_project.GetNkeys() != 0 || project_ratio.GetNkeys() != 0) {
        throw std::runtime_error("Failed charge-ratio inputs produced partial output");
    }
}
} // namespace

int main()
{
    gROOT->SetBatch(true);
    TH1::AddDirectory(false);
    check_ratio(.05, false);
    check_ratio(.15, false);
    check_ratio(.05, true);
    check_unequal_occupancy();
    check_zero_weights();
    check_missing_inputs();
    std::cout << "Charge-ratio normalization and occupancy tests passed\n";
}
