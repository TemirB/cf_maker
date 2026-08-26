#pragma once

#include <string>
#include <utility>

#include <TH3D.h>

class TFile;

[[nodiscard]] std::pair<TH3D*, TH3D*> get_hists(TFile* f, int ch, int centr, int bin);
[[nodiscard]] std::string get_cf_name(int ch, int centr, const std::string& binType,
                                    const std::string& binName);
