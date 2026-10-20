#include "io/run_manifest.h"

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

#include "config/config.h"
#include "core/fs.h"
#include "fit/types.h"

namespace
{
std::string utc_now()
{
    const auto now = std::chrono::system_clock::now();
    const auto seconds = std::chrono::time_point_cast<std::chrono::seconds>(now);
    const auto milliseconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(now - seconds).count();
    const auto timestamp = std::chrono::system_clock::to_time_t(seconds);
    std::tm utc{};
#ifdef _WIN32
    if (gmtime_s(&utc, &timestamp) != 0) {
#else
    if (!gmtime_r(&timestamp, &utc)) {
#endif
        throw std::runtime_error("cannot format run timestamp");
    }
    std::ostringstream text;
    text << std::put_time(&utc, "%Y-%m-%dT%H:%M:%S") << '.' << std::setw(3) << std::setfill('0')
         << milliseconds << 'Z';
    return text.str();
}

std::string make_run_id()
{
    std::random_device random;
    const auto timestamp = std::chrono::duration_cast<std::chrono::nanoseconds>(
                               std::chrono::system_clock::now().time_since_epoch())
                               .count();
    std::ostringstream text;
    text << timestamp << '-' << std::hex << random() << random();
    return text.str();
}

nlohmann::json fit_summary(const Config& cfg)
{
    std::size_t requested = 0;
    std::size_t usable = 0;
    std::size_t retried = 0;
    std::size_t at_limit = 0;
    if (cfg.stages.cf3d || cfg.stages.dependency || cfg.stages.projections_1d) {
        for (const int ch : cfg.selection.charges) {
            for (const int centr : cfg.selection.centralities) {
                for (int b = 0; b < cfg.binning.count; ++b) {
                    const auto& fit = cfg.fit_results.at(static_cast<std::size_t>(ch))
                                          .at(static_cast<std::size_t>(centr))
                                          .at(static_cast<std::size_t>(b));
                    ++requested;
                    if (is_usable_fit(fit)) {
                        ++usable;
                    }
                    if (fit.attempts > 1) {
                        ++retried;
                    }
                    if (fit.at_limit) {
                        ++at_limit;
                    }
                }
            }
        }
    }
    return {{"requested", requested},
            {"usable", usable},
            {"unusable", requested - usable},
            {"retried", retried},
            {"at_limit", at_limit}};
}

void write_atomic(const std::filesystem::path& destination, const nlohmann::json& value,
                  const std::string& run_id)
{
    const auto temporary =
        destination.parent_path() / ("." + destination.filename().string() + ".tmp_" + run_id);
    try {
        std::ofstream output(temporary, std::ios::out | std::ios::trunc);
        if (!output) {
            throw std::runtime_error("cannot create run manifest: " + temporary.string());
        }
        output << value.dump(2) << '\n';
        output.close();
        if (!output) {
            throw std::runtime_error("cannot write run manifest: " + temporary.string());
        }
        std::filesystem::rename(temporary, destination);
    } catch (...) {
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        throw;
    }
}
} // namespace

struct RunManifest::State
{
    std::filesystem::path output_dir;
    nlohmann::json manifest;
    bool tracking = false;
};

RunManifest::RunManifest(const Config& cfg) : state_(std::make_unique<State>())
{
    state_->output_dir = std::filesystem::absolute(cfg.output.dir).lexically_normal();
    const auto snapshot_path = state_->output_dir / "run_config.json";
    std::ifstream input(snapshot_path);
    if (!input) {
        throw std::runtime_error("cannot read config snapshot: " + snapshot_path.string());
    }
    nlohmann::json snapshot;
    input >> snapshot;
    if (!input.eof() && input.fail()) {
        throw std::runtime_error("cannot read config snapshot: " + snapshot_path.string());
    }
    snapshot.at("vars").at("input").at("file") =
        std::filesystem::absolute(cfg.input.file).lexically_normal().string();
    snapshot.at("vars").at("output").at("dir") = state_->output_dir.string();

    const auto run_id = make_run_id();
    state_->manifest = {
        {"schema_version", 1},
        {"status", "running"},
        {"run_id", run_id},
        {"started_at", utc_now()},
        {"completed_at", nullptr},
        {"exit_code", nullptr},
        {"config", std::move(snapshot)},
        {"fits",
         {{"requested", 0}, {"usable", 0}, {"unusable", 0}, {"retried", 0}, {"at_limit", 0}}},
        {"images", nlohmann::json::array()}};
    begin_canvas_tracking();
    state_->tracking = true;
    try {
        write_atomic(state_->output_dir / "run_manifest.json", state_->manifest, run_id);
    } catch (...) {
        end_canvas_tracking();
        throw;
    }
}

RunManifest::~RunManifest()
{
    if (state_->tracking) {
        end_canvas_tracking();
    }
}

void RunManifest::finish(const Config& cfg, int exit_code)
{
    if (exit_code != 0 && exit_code != 2) {
        throw std::invalid_argument("run manifest requires exit code 0 or 2");
    }
    if (!state_->tracking) {
        throw std::logic_error("run manifest is already finished");
    }
    auto completed = state_->manifest;
    completed["status"] = exit_code == 0 ? "completed" : "incomplete";
    completed["completed_at"] = utc_now();
    completed["exit_code"] = exit_code;
    completed["fits"] = fit_summary(cfg);
    for (const auto& filename : saved_canvas_paths()) {
        const std::filesystem::path image(filename);
        const auto relative = image.lexically_relative(state_->output_dir);
        if (relative.empty() || relative.is_absolute() || *relative.begin() == "..") {
            throw std::runtime_error("saved image is outside output directory: " + filename);
        }
        if (!std::filesystem::is_regular_file(image)) {
            throw std::runtime_error("saved image is missing: " + filename);
        }
        const auto size = std::filesystem::file_size(image);
        if (size == 0) {
            throw std::runtime_error("saved image is empty: " + filename);
        }
        completed["images"].push_back({{"path", relative.generic_string()}, {"size_bytes", size}});
    }
    write_atomic(state_->output_dir / "run_manifest.json", completed,
                 completed.at("run_id").get<std::string>());
    end_canvas_tracking();
    state_->tracking = false;
    state_->manifest = std::move(completed);
}
