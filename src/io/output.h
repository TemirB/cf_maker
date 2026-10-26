#pragma once

#include <memory>
#include <string>

class TFile;
class TObject;

// Files created here also retain errors returned by the operating-system close.
[[nodiscard]] std::unique_ptr<TFile> create_output_file(const std::string& path);
void write_output_object(TFile& file, TObject& object, const char* name = nullptr, int options = 0);
void finish_output_file(TFile& file);
