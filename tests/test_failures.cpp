#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <sys/wait.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH3F.h>
#include <TList.h>
#include <TMultiGraph.h>
#include <TNamed.h>
#include <TObjString.h>
#include <TObject.h>
#include <TROOT.h>
#include <nlohmann/json.hpp>

#include "io/output.h"

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
    check(static_cast<bool>(file), "cannot read test file");
    return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

std::string shell_quote(const std::string& text)
{
    std::string result = "'";
    for (const char ch : text) {
        result += ch == '\'' ? "'\\''" : std::string(1, ch);
    }
    return result + "'";
}

class UnwritableObject final : public TObject
{
  public:
    using TObject::Write;
    Int_t Write(const char*, Int_t, Int_t) override
    {
        return 0;
    }
};

void test_checked_output(const std::filesystem::path& dir)
{
    const auto path = dir / "checked.root";
    auto output = create_output_file(path.string());
    TNamed marker("checked_marker", "saved");
    write_output_object(*output, marker);
    finish_output_file(*output);
    check(!output->IsOpen(), "finalized output remains open");
    {
        TFile readback(path.c_str(), "READ");
        auto* saved = dynamic_cast<TNamed*>(readback.Get("checked_marker"));
        check(!readback.IsZombie() && saved && std::string(saved->GetTitle()) == "saved",
              "checked output was not persisted");
        check_throws([&] { write_output_object(readback, marker); },
                     "writing to a read-only file was accepted");
    }
    auto rejected = create_output_file((dir / "rejected_object.root").string());
    UnwritableObject object;
    check_throws([&] { write_output_object(*rejected, object); },
                 "a zero-byte TObject write was accepted");
    rejected->SetBit(TFile::kWriteError);
    check_throws([&] { write_output_object(*rejected, marker); },
                 "a file with an outstanding write error was accepted");
    check_throws([&] { finish_output_file(*rejected); },
                 "finalization ignored an outstanding write error");
    rejected->ResetBit(TFile::kWriteError);
    finish_output_file(*rejected);
    check_throws([&] { (void)create_output_file(dir.string()); },
                 "a directory output was accepted");
}

void create_input(const std::filesystem::path& path, bool include_histograms, bool variance)
{
    TFile output(path.c_str(), "RECREATE");
    if (!include_histograms) {
        output.Close();
        return;
    }
    TH3D den("bp_0_0_num_0", "", 12, -0.2, 0.2, 12, -0.2, 0.2, 12, -0.2, 0.2);
    TH3D num("bp_0_0_num_wei_0", "", 12, -0.2, 0.2, 12, -0.2, 0.2, 12, -0.2, 0.2);
    den.Sumw2();
    num.Sumw2();
    for (int x = 1; x <= den.GetNbinsX(); ++x) {
        for (int y = 1; y <= den.GetNbinsY(); ++y) {
            for (int z = 1; z <= den.GetNbinsZ(); ++z) {
                const double qx = den.GetXaxis()->GetBinCenter(x);
                const double qy = den.GetYaxis()->GetBinCenter(y);
                const double qz = den.GetZaxis()->GetBinCenter(z);
                const double mean =
                    1 +
                    0.7 * std::exp(-(16 * qx * qx + 25 * qy * qy + 36 * qz * qz) / (0.197 * 0.197));
                den.SetBinContent(x, y, z, 1000);
                den.SetBinError(x, y, z, std::sqrt(1000.));
                num.SetBinContent(x, y, z, 1000 * mean);
                num.SetBinError(x, y, z, std::sqrt(1000 * (mean * mean + (variance ? 0.1 : 0.))));
            }
        }
    }
    output.cd();
    den.Write();
    num.Write();
    output.Close();
}

void create_float_input(const std::filesystem::path& path)
{
    TFile output(path.c_str(), "RECREATE");
    for (int charge = 0; charge < 2; ++charge) {
        const auto prefix = "bp_" + std::to_string(charge) + "_0_num_";
        TH3F count((prefix + "0").c_str(), "", 12, -0.2, 0.2, 12, -0.2, 0.2, 12, -0.2, 0.2);
        TH3F sum((prefix + "wei_0").c_str(), "", 12, -0.2, 0.2, 12, -0.2, 0.2, 12, -0.2, 0.2);
        count.Sumw2();
        sum.Sumw2();
        for (int x = 1; x <= count.GetNbinsX(); ++x) {
            for (int y = 1; y <= count.GetNbinsY(); ++y) {
                for (int z = 1; z <= count.GetNbinsZ(); ++z) {
                    const double qx = count.GetXaxis()->GetBinCenter(x);
                    const double qy = count.GetYaxis()->GetBinCenter(y);
                    const double qz = count.GetZaxis()->GetBinCenter(z);
                    const double mean =
                        1 + 0.7 * std::exp(-(16 * qx * qx + 25 * qy * qy + 36 * qz * qz) /
                                           (0.197 * 0.197));
                    const auto fill_pair = [&](double weight) {
                        count.Fill(qx, qy, qz, 1.);
                        sum.Fill(qx, qy, qz, weight);
                    };
                    // These cells also form unresolved 1D/2D projected bins. The
                    // fit range excludes them, so the fitted truth stays Gaussian.
                    if (x <= 2 && (y == 6 || y == 7) && (z == 6 || z == 7)) {
                        if (x == 1) {
                            fill_pair(1.2);
                        } else {
                            fill_pair(1.2 - 1e-6);
                            fill_pair(1.2 + 1e-6);
                        }
                    } else {
                        for (int pair = 0; pair < 64; ++pair) {
                            fill_pair(mean + (pair % 2 == 0 ? -0.25 : 0.25));
                        }
                    }
                }
            }
        }
        const int single_pair_cell = sum.GetBin(1, 6, 6);
        const double rounded_sum = sum.GetBinContent(single_pair_cell);
        check(sum.GetSumw2()->At(single_pair_cell) < rounded_sum * rounded_sum,
              "TH3F fixture did not reproduce single-pair negative centered moments");
        output.cd();
        count.Write();
        sum.Write();
    }
    output.Close();
}

nlohmann::json config_for(const std::filesystem::path& dir, const std::string& input,
                          const std::string& output)
{
    return {{"machine", {{"base_input", dir.string()}, {"base_output", dir.string()}}},
            {"vars", {{"input", {{"file", input}, {"type", "kt"}}}, {"output", {{"dir", output}}}}},
            {"binning", {{"values", {0.15, 0.25}}, {"names", {"test"}}}},
            {"selection", {{"charges", {0}}, {"centralities", {0}}}},
            {"images", {{"need", false}}},
            {"threads", 1},
            {"stages",
             {{"cf3d", true},
              {"dependency", false},
              {"projections_1d", false},
              {"projections_2d", false},
              {"ratios", false}}}};
}

void run_case(const std::filesystem::path& dir, const std::string& executable,
              const std::string& name, const nlohmann::json& config, int expected_status)
{
    const auto config_path = dir / (name + ".json");
    const auto log_path = dir / (name + ".log");
    std::ofstream(config_path) << config.dump(2);
    const std::string command = shell_quote(executable) + " " + shell_quote(config_path.string()) +
                                " > " + shell_quote(log_path.string()) + " 2>&1";
    const int process_status = std::system(command.c_str());
    if (process_status == -1 || !WIFEXITED(process_status) ||
        WEXITSTATUS(process_status) != expected_status) {
        std::cerr << name << ": expected exit " << expected_status << ", actual " << process_status
                  << "\n"
                  << read_bytes(log_path) << "\n";
        throw std::runtime_error("incorrect CLI exit status");
    }
    std::string log = read_bytes(log_path);
    std::transform(log.begin(), log.end(), log.begin(), [](unsigned char character) {
        return static_cast<char>(std::tolower(character));
    });
    const bool success_message = log.find("all outputs written") != std::string::npos;
    check(success_message == (expected_status == 0), "CLI reported success for an incomplete run");
}

void test_float_pipeline(const std::filesystem::path& dir, const std::string& executable)
{
    create_float_input(dir / "float_input.root");
    auto config = config_for(dir, "float_input.root", "float_pipeline");
    config["selection"]["charges"] = {0, 1};
    config["images"] = {{"need", true}, {"format", "png"}};
    config["fit"] = {{"q_max", 0.14}, {"use_default_ip", true}};
    config["projections"] = {{"slice_1d", 0.025}, {"slice_2d", 0.025}, {"slice_ratio", 0.025}};
    config["stages"] = {{"cf3d", true},
                        {"dependency", true},
                        {"projections_1d", true},
                        {"projections_2d", true},
                        {"ratios", true}};
    run_case(dir, executable, "float_pipeline", config, 0);

    const auto output_dir = dir / "float_pipeline";
    TFile cf_file((output_dir / "cf3d.root").c_str(), "READ");
    auto* diagnostics = dynamic_cast<TObjString*>(cf_file.Get("cf_maker_fit_results"));
    check(!cf_file.IsZombie() && diagnostics, "TH3F pipeline did not persist fit diagnostics");
    const auto grid = nlohmann::json::parse(diagnostics->GetString().Data())["fit_grid"];
    TFile projections_1d((output_dir / "1d.root").c_str(), "READ");
    TFile projections_2d((output_dir / "2d.root").c_str(), "READ");
    TFile graphs_file((output_dir / "kt.root").c_str(), "READ");
    check(!projections_1d.IsZombie() && !projections_2d.IsZombie() && !graphs_file.IsZombie(),
          "TH3F pipeline did not create projection and graph files");
    for (int charge = 0; charge < 2; ++charge) {
        const auto fit = grid[charge][0][0];
        check(fit["ok"] == true && fit["at_limit"] == false && fit["ndf"].get<int>() > 0 &&
                  std::abs(fit["r"][0].get<double>() - 4.) < 0.02 &&
                  std::abs(fit["r"][1].get<double>() - 5.) < 0.02 &&
                  std::abs(fit["r"][2].get<double>() - 6.) < 0.02 &&
                  std::abs(fit["lambda"].get<double>() - 0.7) < 0.002,
              "TH3F pipeline did not recover Gaussian fit truth");
        const std::string charge_name = charge == 0 ? "pos" : "neg";
        const auto cf_name = "CF (charge=" + charge_name + ", centrality=0-10, kt=test)";
        auto* cf = dynamic_cast<TH3D*>(cf_file.Get(cf_name.c_str()));
        check(cf && cf->GetBinError(1, 6, 6) == 0. && cf->GetBinError(2, 6, 6) == 0. &&
                  std::abs(cf->GetBinContent(1, 6, 6) - 1.2) < 1e-6 &&
                  std::abs(cf->GetBinContent(2, 6, 6) - 1.2) < 1e-6 &&
                  cf->GetBinError(6, 6, 6) > 0.,
              "TH3F pipeline did not mask unresolved cells and retain resolved variance");
        auto* canvas_1d = dynamic_cast<TCanvas*>(projections_1d.Get((cf_name + "_out").c_str()));
        auto* projected_1d =
            canvas_1d ? dynamic_cast<TH1D*>(canvas_1d->GetPrimitive(cf_name.c_str())) : nullptr;
        check(projected_1d && projected_1d->GetBinError(1) == 0. &&
                  projected_1d->GetBinError(2) == 0. &&
                  std::abs(projected_1d->GetBinContent(1) - 1.2) < 1e-6 &&
                  projected_1d->GetBinError(6) > 0.,
              "TH3F provenance was lost in the saved 1D projection");
        const auto projection_name = cf_name + " out-long";
        auto* canvas_2d = dynamic_cast<TCanvas*>(projections_2d.Get(projection_name.c_str()));
        auto* projected_2d =
            canvas_2d ? dynamic_cast<TH2D*>(canvas_2d->GetPrimitive(projection_name.c_str()))
                      : nullptr;
        check(projected_2d && projected_2d->GetBinError(1, 6) == 0. &&
                  projected_2d->GetBinError(2, 6) == 0. && projected_2d->GetBinError(6, 6) > 0.,
              "TH3F provenance was lost in the saved 2D projection");
        const auto graph_name = "mg_R_out_" + charge_name;
        auto* graph_canvas = dynamic_cast<TCanvas*>(graphs_file.Get(graph_name.c_str()));
        auto* multigraph =
            graph_canvas
                ? dynamic_cast<TMultiGraph*>(graph_canvas->GetPrimitive(graph_name.c_str()))
                : nullptr;
        auto* graph = multigraph && multigraph->GetListOfGraphs()
                          ? dynamic_cast<TGraphErrors*>(multigraph->GetListOfGraphs()->At(0))
                          : nullptr;
        check(graph && graph->GetN() == 1 && std::abs(graph->GetPointY(0) - 4.) < 0.02,
              "TH3F pipeline did not persist the fitted radius graph");
        for (const auto& relative :
             {"dependency/c_all_graphs_" + charge_name + ".png",
              "all_1d_histos/cfs_" + charge_name + "_0-10_test.png",
              "all_2d_histos/all_out-long_2d_histos_centr_0-10_" + charge_name + ".png"}) {
            const auto image = output_dir / relative;
            check(std::filesystem::is_regular_file(image) && std::filesystem::file_size(image) > 0,
                  "TH3F pipeline did not publish a requested nonempty image");
        }
    }
    for (const auto* filename : {"ratio_projs.root", "proj_ratios.root"}) {
        TFile ratios((output_dir / filename).c_str(), "READ");
        const auto name = std::string(filename) == "ratio_projs.root" ? "ratio_proj_0_0_out"
                                                                      : "proj_of_ratios_0_0_out";
        auto* ratio = dynamic_cast<TH1D*>(ratios.Get(name));
        check(!ratios.IsZombie() && ratio && std::abs(ratio->GetBinContent(6) - 1.) < 1e-6 &&
                  ratio->GetBinError(1) == 0. && ratio->GetBinError(2) == 0. &&
                  ratio->GetBinError(6) > 0.,
              "TH3F pipeline did not write a valid charge ratio");
    }
}

void test_cli_failures(const std::filesystem::path& dir, const std::string& executable)
{
    create_input(dir / "input.root", true, true);
    create_input(dir / "no_histograms.root", false, false);
    create_input(dir / "zero_variance.root", true, false);

    run_case(dir, executable, "valid", config_for(dir, "input.root", "valid"), 0);
    auto boundary = config_for(dir, "input.root", "boundary_fit");
    boundary["fit"] = {
        {"limits", {{"radius_sq", {0.0, 100.0}}, {"lambda", {0.3, 0.5}}, {"cross", {0.0, 1.0}}}},
        {"freeze",
         {{"r_out", 4.0},
          {"r_side", 5.0},
          {"r_long", 6.0},
          {"r_os", 0.0},
          {"r_ol", 0.0},
          {"r_sl", 0.0}}}};
    run_case(dir, executable, "boundary_fit", boundary, 2);
    const auto boundary_log = read_bytes(dir / "boundary_fit.log");
    check(boundary_log.find("parameter at limit;") != std::string::npos &&
              boundary_log.find("reached_limits={lambda=") != std::string::npos &&
              boundary_log.find("in [0.3, 0.5]") != std::string::npos,
          "boundary diagnostics did not identify lambda or excluded frozen parameters incorrectly");
    {
        TFile diagnostics((dir / "boundary_fit/cf3d.root").c_str(), "READ");
        auto* text = dynamic_cast<TObjString*>(diagnostics.Get("cf_maker_fit_results"));
        check(!diagnostics.IsZombie() && text, "boundary-fit diagnostics were discarded");
        const auto result = nlohmann::json::parse(text->GetString().Data())["fit_grid"][0][0][0];
        check(result["at_limit"] == true && result["attempts"] == 2,
              "boundary-fit status or performed retry count was not retained");
    }
    auto valid_2d = config_for(dir, "input.root", "valid_2d");
    valid_2d["stages"]["cf3d"] = false;
    valid_2d["stages"]["projections_2d"] = true;
    run_case(dir, executable, "valid_2d", valid_2d, 0);

    auto images = valid_2d;
    images["vars"]["output"]["dir"] = "valid_image";
    images["images"] = {{"need", true}, {"format", "png"}};
    run_case(dir, executable, "valid_image", images, 0);
    const auto image_name = "all_2d_histos/all_out-long_2d_histos_centr_0-10_pos.png";
    const auto image = dir / "valid_image" / image_name;
    check(std::filesystem::is_regular_file(image) && std::filesystem::file_size(image) > 0,
          "successful image stage did not publish a nonempty image");
    images["vars"]["output"]["dir"] = "directory_image";
    const auto image_collision = dir / "directory_image" / image_name;
    std::filesystem::create_directories(image_collision);
    run_case(dir, executable, "directory_image", images, 1);
    check(std::filesystem::is_directory(image_collision),
          "image publication replaced an existing directory");

    auto missing_charge = config_for(dir, "input.root", "missing_charge");
    missing_charge["stages"]["ratios"] = true;
    run_case(dir, executable, "missing_charge", missing_charge, 1);
    check(!std::filesystem::exists(dir / "missing_charge"),
          "missing-charge preflight created an output directory");

    run_case(dir, executable, "missing_fits", config_for(dir, "no_histograms.root", "missing_fits"),
             2);
    check(std::filesystem::is_regular_file(dir / "missing_fits/cf3d.root"),
          "missing-fit diagnostics were discarded");
    run_case(dir, executable, "failed_fits", config_for(dir, "zero_variance.root", "failed_fits"),
             2);
    {
        TFile diagnostics((dir / "failed_fits/cf3d.root").c_str(), "READ");
        auto* text = dynamic_cast<TObjString*>(diagnostics.Get("cf_maker_fit_results"));
        check(!diagnostics.IsZombie() && text, "failed-fit diagnostics were discarded");
        const auto result = nlohmann::json::parse(text->GetString().Data())["fit_grid"][0][0][0];
        check(result["ok"] == false && result["attempts"] == 2,
              "failed-fit status or performed retry count was not retained");
    }
    const auto failed_log = read_bytes(dir / "failed_fits.log");
    check(failed_log.find("usable=0/1, retried=1") != std::string::npos,
          "fit summary did not count an unsuccessful retry");
    check(failed_log.find("fit (ch=0, centr=0, b=0): unusable:") != std::string::npos &&
              failed_log.find("minimizer result rejected;") != std::string::npos &&
              failed_log.find("attempts=2") != std::string::npos,
          "unusable fit diagnostics did not identify the failed task and attempts");

    auto projections = config_for(dir, "no_histograms.root", "missing_2d");
    projections["stages"]["cf3d"] = false;
    projections["stages"]["projections_2d"] = true;
    run_case(dir, executable, "missing_2d", projections, 1);

    auto one_dimension = config_for(dir, "no_histograms.root", "missing_1d");
    one_dimension["stages"]["projections_1d"] = true;
    run_case(dir, executable, "missing_1d", one_dimension, 1);

    const auto collision = dir / "parent_file";
    std::ofstream(collision) << "preserve parent file";
    const auto before = read_bytes(collision);
    run_case(dir, executable, "parent_file", config_for(dir, "input.root", "parent_file/child"), 1);
    check(read_bytes(collision) == before, "output initialization changed a parent file");

    std::filesystem::create_directories(dir / "directory_cf3d/cf3d.root");
    run_case(dir, executable, "directory_cf3d", config_for(dir, "input.root", "directory_cf3d"), 1);
    check(std::filesystem::is_directory(dir / "directory_cf3d/cf3d.root"),
          "output initialization replaced a directory");

    std::filesystem::create_directories(dir / "directory_2d/2d.root");
    auto two_dimension = config_for(dir, "input.root", "directory_2d");
    two_dimension["stages"]["cf3d"] = false;
    two_dimension["stages"]["projections_2d"] = true;
    run_case(dir, executable, "directory_2d", two_dimension, 1);

    std::filesystem::create_directories(dir / "directory_snapshot/run_config.json");
    run_case(dir, executable, "directory_snapshot",
             config_for(dir, "input.root", "directory_snapshot"), 1);

    std::filesystem::create_directories(dir / "directory_log/run.log");
    run_case(dir, executable, "directory_log", config_for(dir, "input.root", "directory_log"), 1);
}
} // namespace

int main(int argc, char** argv)
{
    gROOT->SetBatch(true);
    TH1::AddDirectory(false);
    const auto dir = std::filesystem::temp_directory_path() /
                     ("cf_maker_failures_" +
                      std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directories(dir);
    test_checked_output(dir);
    if (argc == 2) {
        test_float_pipeline(dir, argv[1]);
        test_cli_failures(dir, argv[1]);
    }
    std::cout << "Checked output and CLI failure reporting are correct\n";
}
