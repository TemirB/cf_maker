#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH3D.h>
#include <TList.h>
#include <TMultiGraph.h>
#include <TNamed.h>
#include <TObjString.h>
#include <TROOT.h>
#include <nlohmann/json.hpp>

#include "config/config.h"
#include "io/fit_results.h"
#include "io/input.h"

namespace
{
void check(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

template <typename Action> void check_throws(Action action, const char* message)
{
    try {
        action();
    } catch (const std::exception&) {
        return;
    }
    throw std::runtime_error(message);
}

std::string read_bytes(const std::filesystem::path& path)
{
    std::ifstream file(path, std::ios::binary);
    check(static_cast<bool>(file), "cannot read protected output");
    return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

void check_double(double left, double right)
{
    check(left == right || (std::isnan(left) && std::isnan(right)), "fit value did not round-trip");
}

void check_grid(const FitGrid& expected, const FitGrid& actual)
{
    check(expected.size() == actual.size(), "charge count changed");
    for (std::size_t ch = 0; ch < expected.size(); ++ch) {
        check(expected[ch].size() == actual[ch].size(), "centrality count changed");
        for (std::size_t centr = 0; centr < expected[ch].size(); ++centr) {
            check(expected[ch][centr].size() == actual[ch][centr].size(), "bin count changed");
            for (std::size_t b = 0; b < expected[ch][centr].size(); ++b) {
                const auto& left = expected[ch][centr][b];
                const auto& right = actual[ch][centr][b];
                for (std::size_t i = 0; i < left.r.size(); ++i) {
                    check_double(left.r[i], right.r[i]);
                    check_double(left.e_r[i], right.e_r[i]);
                }
                check_double(left.lambda, right.lambda);
                check_double(left.e_lambda, right.e_lambda);
                check_double(left.chi2, right.chi2);
                check_double(left.p_value, right.p_value);
                for (std::size_t i = 0; i < left.corr.size(); ++i) {
                    check_double(left.corr[i], right.corr[i]);
                }
                check(left.ndf == right.ndf && left.ok == right.ok && left.status == right.status &&
                          left.cov_status == right.cov_status && left.at_limit == right.at_limit &&
                          left.attempts == right.attempts,
                      "fit diagnostics did not round-trip");
            }
        }
    }
}

Config make_fixture(const std::filesystem::path& dir)
{
    Config cfg;
    cfg.input.type = "kt";
    cfg.input.file = (dir / "input.root").string();
    cfg.output.dir = (dir / "results").string();
    cfg.binning.values = {0.15, 0.25};
    cfg.binning.names = {"test"};
    cfg.selection.charges = {0};
    cfg.selection.centralities = {0};
    cfg.general.images.need = false;
    build(cfg);

    auto& fit = cfg.fit_results[0][0][0];
    fit.r = {4., 5., 6., -2., 3., -4.};
    fit.e_r = {0.11, 0.12, 0.13, 0.14, 0.15, 0.16};
    fit.lambda = 0.7;
    fit.e_lambda = 0.023;
    fit.chi2 = 12.5;
    fit.ndf = 25;
    fit.p_value = 0.987;
    fit.ok = true;
    fit.status = 0;
    fit.cov_status = 3;
    fit.attempts = 2;
    for (std::size_t i = 0; i < fit.corr.size(); ++i) {
        fit.corr[i] = static_cast<double>(i) / 100.;
    }

    // Preserve nonfinite diagnostics and unavailable fits, rather than turn them into zeros.
    auto& failed = cfg.fit_results[1][3][0];
    failed.r[4] = std::numeric_limits<double>::quiet_NaN();
    failed.chi2 = std::numeric_limits<double>::infinity();
    failed.corr[12] = -std::numeric_limits<double>::infinity();
    failed.status = 4;
    failed.cov_status = 1;
    failed.at_limit = true;
    failed.attempts = 1;

    TFile input(cfg.input.file.c_str(), "RECREATE");
    TH3D den("bp_0_0_num_0", "", 8, -0.2, 0.2, 8, -0.2, 0.2, 8, -0.2, 0.2);
    TH3D num("bp_0_0_num_wei_0", "", 8, -0.2, 0.2, 8, -0.2, 0.2, 8, -0.2, 0.2);
    den.Sumw2();
    num.Sumw2();
    for (int x = 1; x <= 8; ++x) {
        for (int y = 1; y <= 8; ++y) {
            for (int z = 1; z <= 8; ++z) {
                const double qx = den.GetXaxis()->GetBinCenter(x);
                const double qy = den.GetYaxis()->GetBinCenter(y);
                const double qz = den.GetZaxis()->GetBinCenter(z);
                const double quadratic = 16 * qx * qx + 25 * qy * qy + 36 * qz * qz - 4 * qx * qy +
                                         6 * qx * qz - 8 * qy * qz;
                const double mean = 1 + 0.7 * std::exp(-quadratic / (0.197 * 0.197));
                den.SetBinContent(x, y, z, 100);
                den.SetBinError(x, y, z, 10);
                num.SetBinContent(x, y, z, 100 * mean);
                num.SetBinError(x, y, z, std::sqrt(100 * (mean * mean + 0.1)));
            }
        }
    }
    input.cd();
    den.Write();
    num.Write();
    input.Close();

    std::filesystem::create_directories(cfg.output.dir);
    TFile output((cfg.output.dir + "/cf3d.root").c_str(), "RECREATE");
    TH3D cf(num);
    cf.Divide(&num, &den);
    output.cd();
    cf.Write(get_cf_name(0, 0, cfg.input.type, cfg.binning.names[0]).c_str());
    write_fit_results(output, cfg);
    output.Close();
    return cfg;
}

void test_metadata(const std::filesystem::path& dir, const Config& cfg)
{
    TFile file((cfg.output.dir + "/cf3d.root").c_str(), "READ");
    auto* pair_metadata = dynamic_cast<TObjString*>(file.Get("cf_maker_fit_results"));
    check(pair_metadata &&
              nlohmann::json::parse(pair_metadata->GetString().Data())["identity"]["statistics"] ==
                  "weighted-mean-source-roundoff-v2",
          "pair-weight statistics version was not recorded");
    Config loaded = cfg;
    build(loaded);
    read_fit_results(file, loaded);
    check_grid(cfg.fit_results, loaded.fit_results);

    loaded = cfg;
    loaded.threads = 7;
    loaded.output.dir = "another/output";
    loaded.projections.slice_1d = 0.04;
    loaded.general.images.format = "png";
    loaded.stages.cf3d = false;
    loaded.logging.level = "debug";
    loaded.binning.file_names = {"changed_image_filename"};
    build(loaded);
    read_fit_results(file, loaded);
    check_grid(cfg.fit_results, loaded.fit_results);

    for (int change = 0; change < 5; ++change) {
        loaded = cfg;
        build(loaded);
        switch (change) {
        case 0:
            loaded.fit.q_max += 0.01;
            break;
        case 1:
            loaded.fit.freeze[4].reset();
            break;
        case 2:
            loaded.selection.centralities = {1};
            break;
        case 3:
            loaded.binning.values[0] += 0.01;
            break;
        default:
            loaded.binning.names[0] = "renamed";
        }
        check_throws([&] { read_fit_results(file, loaded); }, "incompatible saved fits accepted");
        check(loaded.fit_results[0][0][0].attempts == 0, "failed load modified the fit grid");
    }

    const auto changed_input = dir / "changed_input.root";
    std::filesystem::copy_file(cfg.input.file, changed_input);
    std::ofstream(changed_input, std::ios::binary | std::ios::app) << "changed";
    loaded = cfg;
    loaded.input.file = changed_input.string();
    check_throws([&] { read_fit_results(file, loaded); }, "altered input accepted");

    TFile legacy((dir / "legacy.root").c_str(), "RECREATE");
    check_throws([&] { read_fit_results(legacy, loaded); }, "legacy file accepted");
    legacy.Close();

    // Every field is mandatory; loading a partial record must leave the original grid intact.
    const auto corrupted = dir / "corrupted.root";
    std::filesystem::copy_file(cfg.output.dir + "/cf3d.root", corrupted);
    TFile bad_file(corrupted.c_str(), "UPDATE");
    auto* stored = dynamic_cast<TObjString*>(bad_file.Get("cf_maker_fit_results"));
    check(stored != nullptr, "missing test metadata");
    auto metadata = nlohmann::json::parse(stored->GetString().Data());
    metadata["fit_grid"][0][0][0].erase("status");
    bad_file.cd();
    TObjString incomplete(metadata.dump().c_str());
    incomplete.Write("cf_maker_fit_results", TObject::kOverwrite);
    loaded = cfg;
    build(loaded);
    check_throws([&] { read_fit_results(bad_file, loaded); },
                 "incomplete fit diagnostics accepted");
    check(loaded.fit_results[0][0][0].attempts == 0, "incomplete load modified fit grid");

    metadata["schema_version"] = 999;
    TObjString future_schema(metadata.dump().c_str());
    future_schema.Write("cf_maker_fit_results", TObject::kOverwrite);
    check_throws([&] { read_fit_results(bad_file, loaded); }, "unknown saved-fit schema accepted");

    write_fit_results(bad_file, cfg);
    bad_file.Delete((get_cf_name(0, 0, cfg.input.type, cfg.binning.names[0]) + ";*").c_str());
    check_throws([&] { read_fit_results(bad_file, loaded); }, "saved fit without its CF accepted");
    bad_file.Close();

    const auto fixed_input_path = dir / "fixed_input.root";
    std::filesystem::copy_file(cfg.input.file, fixed_input_path);
    {
        TFile fixed_input(fixed_input_path.c_str(), "UPDATE");
        TNamed marker("correlation_statistics", "fixed_reference");
        marker.Write();
    }
    Config fixed = cfg;
    fixed.input.file = fixed_input_path.string();
    const auto fixed_output_path = dir / "fixed_output.root";
    std::filesystem::copy_file(cfg.output.dir + "/cf3d.root", fixed_output_path);
    TFile fixed_output(fixed_output_path.c_str(), "UPDATE");
    write_fit_results(fixed_output, fixed);
    auto* fixed_text = dynamic_cast<TObjString*>(fixed_output.Get("cf_maker_fit_results"));
    check(fixed_text != nullptr, "missing fixed-reference metadata");
    auto fixed_metadata = nlohmann::json::parse(fixed_text->GetString().Data());
    check(fixed_metadata["identity"]["statistics"] == "fixed_reference_v1",
          "fixed-reference statistics version was not recorded");
    loaded = fixed;
    build(loaded);
    read_fit_results(fixed_output, loaded);
    check_grid(fixed.fit_results, loaded.fit_results);
    fixed_metadata["identity"]["statistics"] = "weighted-mean-source-roundoff-v2";
    TObjString wrong_mode(fixed_metadata.dump().c_str());
    fixed_output.cd();
    wrong_mode.Write("cf_maker_fit_results", TObject::kOverwrite);
    check_throws([&] { read_fit_results(fixed_output, loaded); },
                 "metadata with the wrong statistics model was accepted");
    fixed_output.Close();
}

std::string shell_quote(const std::string& text)
{
    std::string result = "'";
    for (const char ch : text) {
        result += ch == '\'' ? "'\\''" : std::string(1, ch);
    }
    return result + "'";
}

void test_cli(const std::filesystem::path& dir, const Config& cfg, const std::string& executable)
{
    nlohmann::json config = {
        {"machine", {{"base_input", dir.string()}, {"base_output", dir.string()}}},
        {"vars",
         {{"input", {{"file", "input.root"}, {"type", "kt"}}}, {"output", {{"dir", "results"}}}}},
        {"binning", {{"values", cfg.binning.values}, {"names", cfg.binning.names}}},
        {"selection", {{"charges", {0}}, {"centralities", {0}}}},
        {"stages",
         {{"cf3d", false},
          {"dependency", true},
          {"projections_1d", true},
          {"projections_2d", false},
          {"ratios", false}}},
        {"images", {{"need", false}}}};
    const auto config_path = dir / "resume_cli.json";
    const auto run = [&] {
        std::ofstream(config_path) << config.dump(2);
        const std::string command = shell_quote(executable) + " " +
                                    shell_quote(config_path.string()) + " > " +
                                    shell_quote((dir / "cli.log").string()) + " 2>&1";
        return std::system(command.c_str());
    };
    check(run() == 0, "valid CLI resume failed");
    {
        TFile dependency((cfg.output.dir + "/kt.root").c_str(), "READ");
        auto* canvas = dynamic_cast<TCanvas*>(dependency.Get("mg_R_out_pos"));
        check(canvas != nullptr, "missing resumed dependency graph");
        auto* multigraph =
            dynamic_cast<TMultiGraph*>(canvas->GetListOfPrimitives()->FindObject("mg_R_out_pos"));
        check(multigraph && multigraph->GetListOfGraphs(), "missing radius multigraph");
        auto* graph = dynamic_cast<TGraphErrors*>(multigraph->GetListOfGraphs()->First());
        check(graph && graph->GetN() == 1 && graph->GetPointY(0) == cfg.fit_results[0][0][0].r[0],
              "resume dependency graph uses zero fit parameters");
    }
    {
        TFile projections((cfg.output.dir + "/1d.root").c_str(), "READ");
        const std::string canvas_name = get_cf_name(0, 0, "kt", cfg.binning.names[0]) + "_out";
        auto* canvas = dynamic_cast<TCanvas*>(projections.Get(canvas_name.c_str()));
        check(canvas != nullptr, "missing resumed 1D projection");
        TH1D* fit = nullptr;
        for (TObject* object : *canvas->GetListOfPrimitives()) {
            auto* histogram = dynamic_cast<TH1D*>(object);
            if (histogram && std::string(histogram->GetName()).find("fit_") == 0) {
                fit = histogram;
            }
        }
        check(fit && fit->GetBinContent(fit->FindBin(0.)) > 1.01,
              "resume 1D model uses zero fit parameters");
    }

    const std::array<std::string, 5> protected_names = {"run_config.json", "run.log", "cf3d.root",
                                                        "kt.root", "1d.root"};
    std::array<std::string, 5> before{};
    for (std::size_t i = 0; i < protected_names.size(); ++i) {
        before[i] = read_bytes(std::filesystem::path(cfg.output.dir) / protected_names[i]);
    }
    config["fit"]["q_max"] = 0.1;
    check(run() != 0, "incompatible CLI resume returned success");
    for (std::size_t i = 0; i < protected_names.size(); ++i) {
        check(before[i] == read_bytes(std::filesystem::path(cfg.output.dir) / protected_names[i]),
              "failed preflight changed an existing output");
    }
    config.erase("fit");
    config["vars"]["output"]["dir"] = "missing_results";
    check(run() != 0, "missing saved fits returned success");
    check(!std::filesystem::exists(dir / "missing_results"),
          "missing resume data created an output directory");
}
} // namespace

int main(int argc, char** argv)
{
    gROOT->SetBatch(true);
    TH1::AddDirectory(false);
    const auto dir = std::filesystem::temp_directory_path() /
                     ("cf_maker_resume_" +
                      std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(dir);
    const Config cfg = make_fixture(dir);
    test_metadata(dir, cfg);
    if (argc == 2) {
        test_cli(dir, cfg, argv[1]);
    }
    std::cout << "Complete saved fits and dependency preflight are correct\n";
}
