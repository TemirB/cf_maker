#pragma once

#include <string>
#include <utility>

#include <TH3D.h>

class TFile;

// Returns independent, directory-detached TH3D copies of TH3F/TH3D input.
// The caller owns both pointers. On failure, returns two null pointers.
[[nodiscard]] std::pair<TH3D*, TH3D*> get_hists(TFile* f, int ch, int centr, int bin);
[[nodiscard]] std::string get_cf_name(int ch, int centr, const std::string& binType,
                                      const std::string& binName);
