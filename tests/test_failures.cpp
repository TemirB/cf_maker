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

#include <TFile.h>
#include <TH3D.h>
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

void test_cli_failures(const std::filesystem::path& dir, const std::string& executable)
{
    create_input(dir / "input.root", true, true);
    create_input(dir / "no_histograms.root", false, false);
    create_input(dir / "zero_variance.root", true, false);

    run_case(dir, executable, "valid", config_for(dir, "input.root", "valid"), 0);
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
        check(result["ok"] == false && result["attempts"].get<int>() > 0,
              "failed-fit status was not retained");
    }

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
        test_cli_failures(dir, argv[1]);
    }
    std::cout << "Checked output and CLI failure reporting are correct\n";
}
