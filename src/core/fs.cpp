#include "core/fs.h"

#include <filesystem>
#include <chrono>
#include <stdexcept>
#include <string>
#include <system_error>

#include <TCanvas.h>
#include <TError.h>

#include "core/log.h"

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
    } catch (...) {
        gErrorIgnoreLevel = prevLevel;
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        throw;
    }
}
