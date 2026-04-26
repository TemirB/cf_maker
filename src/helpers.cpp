#include "helpers.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>

#include <TH3.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TError.h>

#include "fit/types.h"

std::pair<TH3D*, TH3D*> getHists(TFile* f, int ch, int centr, int bin) {
    TString numName = Form("bp_%d_%d_num_%d", ch, centr, bin);
    TH3D* num = (TH3D*) f->Get(numName);

    if (!num) {
        std::cerr << "[Num] NOT FOUND: " << numName << std::endl;
        return {nullptr, nullptr};
    }

    TString weiName = Form("bp_%d_%d_num_wei_%d", ch, centr, bin);
    TH3D* wei = (TH3D*) f->Get(weiName);

    if (!wei) {
        std::cerr << "[NumWei] NOT FOUND: " << weiName << std::endl;
        return {nullptr, nullptr};
    }
    return {num, wei};
}

std::string getCFName(int ch, int centr, const char* binType, const char* binName) {
    return Form(
        "CF at charge=%s, centrality=%s, %s=%s", 
        Charge::kNames[ch], Centrality::kNames[centr], binType, binName
    );
}

bool IsBadFit(const FitResult& r) {
    if (!r.ok) return true;
    if (!r.IsFinite()) return true;
    if (!r.IsValid()) return true;
    return false;
}

AnalysisType getType(const char* type) {
    if (std::strcmp(type, "rapidity") == 0) return AnalysisType::Rapidity;
    if (std::strcmp(type, "kt") == 0) return AnalysisType::Kt;
    return AnalysisType::Unknown;
}

void EnsureDir(const std::string& dir) {
    struct stat st;
    if (stat(dir.c_str(), &st) != 0) {
        mkdir(dir.c_str(), 0755);
    }
}

void SaveCanvasQuiet(TCanvas* canvas, const char* filename) {
    if (!canvas || !filename) {
        return;
    }

    const Int_t prevLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kWarning;
    canvas->SaveAs(filename);
    gErrorIgnoreLevel = prevLevel;
}