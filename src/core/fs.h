#pragma once

#include <string>
#include <vector>

class TCanvas;

void ensure_dir(const std::string& dir);
void save_canvas_quiet(TCanvas* canvas, const char* filename);

// Track only canvas files saved by the active application run. Calls outside a
// run (including standalone tests) do not retain any paths.
void begin_canvas_tracking();
void end_canvas_tracking() noexcept;
[[nodiscard]] std::vector<std::string> saved_canvas_paths();
