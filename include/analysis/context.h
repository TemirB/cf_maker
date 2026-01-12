#pragma once

#include <TFile.h>

#include "helpers.h"

struct RunContext {
    TFile*   input;
    FitGrid* fit;

    std::map<std::string, std::unique_ptr<TFile>> outputs;

    TFile& GetOut(const std::string& name) {
        auto& ptr = outputs[name];
        if (!ptr)
            ptr = std::make_unique<TFile>((name + ".root").c_str(), "RECREATE");
        return *ptr;
    }
};
