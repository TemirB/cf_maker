#pragma once

#include <TFile.h>

#include "helpers.h"

struct RunContext {
    TFile* input;
    FitGrid* fit;
    std::string outDir;

    std::map<std::string, TFile*> outputs;

    TFile& GetOut(const std::string& name) {
        auto it = outputs.find(name);
        if (it != outputs.end()) {
            it->second->cd();
            return *it->second;
        }

        std::string path = outDir + "/" + name + ".root";
        auto* f = new TFile(path.c_str(), "RECREATE");
        outputs[name] = f;
        f->cd();
        return *f;
    }
};

