#include "helpers.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

#include <TH3.h>
#include <TTree.h>
#include <TFile.h>

#include "fit/types.h"

TH3D* getNum(TFile* f, int ch, int centr, int bin) {
    TString name = Form("bp_%d_%d_num_%d", ch, centr, bin);
    TH3D* h = (TH3D*) f->Get(name);

    if (!h) {
        std::cerr << "[getNum] NOT FOUND: " << name << std::endl;
        return nullptr;
    }
    return h;
};

TH3D* getNumWei(TFile* f, int ch, int centr, int bin) {
    TString name = Form("bp_%d_%d_num_wei_%d", ch, centr, bin);
    TH3D* h = (TH3D*) f->Get(name);

    if (!h) {
        std::cerr << "[getNumWei] NOT FOUND: " << name << std::endl;
        return nullptr;
    }
    return h;
};

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