#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <string>

#include <TCanvas.h>
#include <TH1D.h>
#include <TROOT.h>
#include <nlohmann/json.hpp>

#include "config/config.h"
#include "core/fs.h"
#include "io/run_manifest.h"

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

nlohmann::json read_json(const std::filesystem::path& path)
{
    std::ifstream input(path);
    check(static_cast<bool>(input), "cannot read run manifest");
    nlohmann::json value;
    input >> value;
    return value;
}

void check_timestamp(const nlohmann::json& value)
{
    const std::regex utc_pattern(R"(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}Z)");
    check(value.is_string() && std::regex_match(value.get<std::string>(), utc_pattern),
          "run timestamp is not an ISO 8601 UTC timestamp");
}

Config make_config(const std::filesystem::path& output_dir)
{
    Config cfg;
    cfg.machine.base_input = output_dir.parent_path().string();
    cfg.machine.base_output = output_dir.parent_path().string();
    cfg.input.type = "kt";
    cfg.input.file = (output_dir.parent_path() / "input.root").string();
    cfg.output.dir = output_dir.string();
    cfg.binning.values = {0.15, 0.25, 0.35, 0.45, 0.60};
    cfg.binning.names = {"first", "second", "third", "fourth"};
    cfg.selection.charges = {0};
    cfg.selection.centralities = {1};
    cfg.general.images.format = "png";
    build(cfg);
    return cfg;
}

FitResult usable_fit()
{
    FitResult fit;
    fit.ok = true;
    fit.ndf = 12;
    fit.r = {4., 5., 6., 0., 0., 0.};
    fit.lambda = 0.7;
    fit.chi2 = 10.;
    fit.attempts = 1;
    return fit;
}

void save_test_canvas(const std::filesystem::path& filename)
{
    TCanvas canvas("manifest_test_canvas", "", 400, 300);
    TH1D histogram("manifest_test_histogram", "", 10, 0., 1.);
    histogram.Fill(0.4);
    histogram.Draw();
    save_canvas_quiet(&canvas, filename.c_str());
}

void test_complete_run(const std::filesystem::path& dir)
{
    auto cfg = make_config(dir / "completed");
    for (auto& fit : cfg.fit_results[0][1]) {
        fit = usable_fit();
    }
    write_run_config(cfg);
    const auto output_dir = std::filesystem::path(cfg.output.dir);
    const auto old_image = output_dir / "old.png";
    save_test_canvas(old_image);
    check(saved_canvas_paths().empty(), "a canvas saved outside a run entered the collector");
    const auto manifest_path = output_dir / "run_manifest.json";
    std::ofstream(manifest_path) << nlohmann::json{
        {"status", "completed"},
        {"run_id", "previous-run"},
        {"images",
         {{{"path", "old.png"}}}}}.dump();

    RunManifest run(cfg);
    const auto running = read_json(manifest_path);
    check(running.at("status") == "running" && running.at("exit_code").is_null() &&
              running.at("completed_at").is_null() && running.at("images").empty(),
          "starting a run did not invalidate the previous successful manifest");
    check(running.at("schema_version") == 1 && running.at("run_id") != "previous-run" &&
              !running.at("run_id").get<std::string>().empty(),
          "new run did not receive a fresh identity");
    check_timestamp(running.at("started_at"));

    const auto current_image = output_dir / "projections" / "current.png";
    save_test_canvas(current_image);
    save_test_canvas(current_image);
    run.finish(cfg, 0);
    const auto completed = read_json(manifest_path);
    check(completed.at("status") == "completed" && completed.at("exit_code") == 0,
          "successful run was not published as completed");
    check(completed.at("run_id") == running.at("run_id") &&
              completed.at("started_at") == running.at("started_at"),
          "finishing changed the run identity or start time");
    check_timestamp(completed.at("completed_at"));
    check(completed.at("completed_at").get<std::string>() >=
              completed.at("started_at").get<std::string>(),
          "completion preceded the run start");
    check(
        completed.at("fits") ==
            nlohmann::json{
                {"requested", 4}, {"usable", 4}, {"unusable", 0}, {"retried", 0}, {"at_limit", 0}},
        "completed fit counts do not reflect the selected tasks");
    const auto& images = completed.at("images");
    check(images.size() == 1 && images.at(0).at("path") == "projections/current.png" &&
              images.at(0).at("size_bytes").get<std::uintmax_t>() ==
                  std::filesystem::file_size(current_image) &&
              std::filesystem::is_regular_file(old_image),
          "manifest did not preserve precisely the images exported by this run");
    check(saved_canvas_paths().empty(), "finished run left image tracking enabled");
    check_throws([&] { run.finish(cfg, 0); }, "manifest accepted a second completion");
    check(read_json(manifest_path) == completed,
          "second completion changed the published manifest");
}

void test_incomplete_fit_summary(const std::filesystem::path& dir)
{
    auto cfg = make_config(dir / "incomplete");
    auto& fits = cfg.fit_results[0][1];
    fits[0] = usable_fit();
    fits[1] = usable_fit();
    fits[1].at_limit = true;
    fits[1].attempts = 2;
    fits[2].ok = false;
    fits[2].attempts = 2;
    fits[3].ok = true;
    fits[3].ndf = 0;
    // An unselected task must not contribute to any diagnostic count.
    cfg.fit_results[1][0][0] = fits[1];
    write_run_config(cfg);
    RunManifest run(cfg);
    run.finish(cfg, 2);
    const auto manifest = read_json(std::filesystem::path(cfg.output.dir) / "run_manifest.json");
    check(manifest.at("status") == "incomplete" && manifest.at("exit_code") == 2,
          "diagnostic run was reported as complete");
    check(
        manifest.at("fits") ==
            nlohmann::json{
                {"requested", 4}, {"usable", 1}, {"unusable", 3}, {"retried", 2}, {"at_limit", 1}},
        "failed, boundary or underdetermined fits were counted as usable");
    check_timestamp(manifest.at("completed_at"));
}

void test_relative_snapshot_and_unfinished_run(const std::filesystem::path& dir)
{
    const auto previous_dir = std::filesystem::current_path();
    std::filesystem::current_path(dir);
    const auto current_dir = std::filesystem::current_path();
    try {
        auto cfg = make_config("relative/results");
        cfg.input.file = "relative/input.root";
        cfg.fit.q_max = 0.17;
        write_run_config(cfg);
        const auto manifest_path = dir / "relative/results/run_manifest.json";
        {
            RunManifest run(cfg);
            const auto manifest = read_json(manifest_path);
            const auto& config = manifest.at("config");
            check(config.at("vars").at("input").at("file") ==
                          (current_dir / "relative/input.root").string() &&
                      config.at("vars").at("output").at("dir") ==
                          (current_dir / "relative/results").string(),
                  "relative input and output paths were not resolved in the snapshot");
            check(config.at("fit").at("q_max") == cfg.fit.q_max &&
                      config.at("selection").at("charges") == cfg.selection.charges &&
                      config.at("binning").at("values") == cfg.binning.values,
                  "manifest lost the run configuration");
        }
        const auto unfinished = read_json(manifest_path);
        check(unfinished.at("status") == "running" && unfinished.at("exit_code").is_null() &&
                  unfinished.at("completed_at").is_null(),
              "an unfinished run was implicitly marked as successful");
        check(saved_canvas_paths().empty(), "destroyed run retained its image collector");
        RunManifest next_run(cfg);
        check(read_json(manifest_path).at("run_id") != unfinished.at("run_id"),
              "unfinished run prevented a fresh run");
    } catch (...) {
        std::filesystem::current_path(previous_dir);
        throw;
    }
    std::filesystem::current_path(previous_dir);
}

void test_missing_snapshot_and_atomic_failure(const std::filesystem::path& dir)
{
    auto cfg = make_config(dir / "missing_snapshot");
    std::filesystem::create_directories(cfg.output.dir);
    check_throws([&] { RunManifest run(cfg); }, "manifest accepted a missing config snapshot");
    check(!std::filesystem::exists(std::filesystem::path(cfg.output.dir) / "run_manifest.json"),
          "missing snapshot produced a misleading run manifest");

    cfg = make_config(dir / "manifest_directory");
    write_run_config(cfg);
    const auto destination = std::filesystem::path(cfg.output.dir) / "run_manifest.json";
    std::filesystem::create_directory(destination);
    const auto protected_file = destination / "preserve.txt";
    std::ofstream(protected_file) << "preserve existing directory";
    check_throws([&] { RunManifest run(cfg); }, "manifest overwrote an existing directory");
    check(std::filesystem::is_directory(destination) &&
              std::filesystem::is_regular_file(protected_file),
          "failed atomic publication changed the existing directory");
    for (const auto& entry : std::filesystem::directory_iterator(cfg.output.dir)) {
        check(entry.path().filename().string().find(".tmp_") == std::string::npos,
              "failed publication left a temporary manifest");
    }
    check(saved_canvas_paths().empty(), "failed constructor left image tracking active");
    // Another manifest must work even after the constructor failed.
    cfg = make_config(dir / "after_failure");
    write_run_config(cfg);
    RunManifest recovered(cfg);
    recovered.finish(cfg, 2);
}
} // namespace

int main()
{
    try {
        gROOT->SetBatch(true);
        TH1::AddDirectory(false);
        const auto dir =
            std::filesystem::temp_directory_path() /
            ("cf_maker_run_manifest_" +
             std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(dir);
        test_complete_run(dir);
        test_incomplete_fit_summary(dir);
        test_relative_snapshot_and_unfinished_run(dir);
        test_missing_snapshot_and_atomic_failure(dir);
        std::cout << "Run manifest lifecycle, diagnostics and current images are correct\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
