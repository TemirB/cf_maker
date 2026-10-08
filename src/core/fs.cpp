#include "core/fs.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <mutex>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include <TCanvas.h>
#include <TError.h>

#include "core/log.h"

namespace
{
std::mutex canvas_paths_mutex;
bool canvas_tracking_enabled = false;
std::vector<std::string> canvas_paths;

void track_saved_canvas(const std::filesystem::path& destination)
{
    const std::lock_guard<std::mutex> lock(canvas_paths_mutex);
    if (!canvas_tracking_enabled) {
        return;
    }
    const auto path = std::filesystem::absolute(destination).lexically_normal().string();
    if (std::find(canvas_paths.begin(), canvas_paths.end(), path) == canvas_paths.end()) {
        canvas_paths.push_back(path);
    }
}
} // namespace

void begin_canvas_tracking()
{
    const std::lock_guard<std::mutex> lock(canvas_paths_mutex);
    if (canvas_tracking_enabled) {
        throw std::logic_error("canvas tracking is already active");
    }
    canvas_paths.clear();
    canvas_tracking_enabled = true;
}

void end_canvas_tracking() noexcept
{
    const std::lock_guard<std::mutex> lock(canvas_paths_mutex);
    canvas_tracking_enabled = false;
    canvas_paths.clear();
}

std::vector<std::string> saved_canvas_paths()
{
    const std::lock_guard<std::mutex> lock(canvas_paths_mutex);
    return canvas_paths;
}

void ensure_dir(const std::string& dir)
{
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    if (ec) {
        throw std::runtime_error("cannot create directory " + dir + ": " + ec.message());
    }
}

void save_canvas_quiet(TCanvas* canvas, const char* filename)
{
    if (!canvas || !filename) {
        throw std::invalid_argument("cannot save a null canvas or filename");
    }

    const std::filesystem::path destination(filename);
    if (destination.empty()) {
        throw std::invalid_argument("cannot save a canvas to an empty filename");
    }
    if (!destination.parent_path().empty()) {
        ensure_dir(destination.parent_path().string());
    }
    const auto temporary =
        destination.parent_path() /
        ("." + destination.stem().string() + ".tmp_" +
         std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) +
         destination.extension().string());
    const Int_t prevLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kWarning;
    try {
        // SaveAs has no return status. A new temporary filename ensures an old
        // image cannot conceal a failed write; preserve the extension for ROOT.
        canvas->SaveAs(temporary.c_str());
        gErrorIgnoreLevel = prevLevel;
        if (!std::filesystem::is_regular_file(temporary) ||
            std::filesystem::file_size(temporary) == 0) {
            throw std::runtime_error("cannot save canvas: " + destination.string());
        }
        std::filesystem::rename(temporary, destination);
        track_saved_canvas(destination);
    } catch (...) {
        gErrorIgnoreLevel = prevLevel;
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        throw;
    }
}
