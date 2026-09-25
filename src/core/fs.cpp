#include "core/fs.h"

#include <filesystem>
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
        logging::error("cannot create directory " + dir + ": " + ec.message());
    }
}

void save_canvas_quiet(TCanvas* canvas, const char* filename)
{
    if (!canvas || !filename) {
        return;
    }

    const Int_t prevLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kWarning;
    canvas->SaveAs(filename);
    gErrorIgnoreLevel = prevLevel;
}
