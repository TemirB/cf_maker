#pragma once

#include <string>

class TCanvas;

void ensure_dir(const std::string& dir);
void save_canvas_quiet(TCanvas* canvas, const char* filename);
