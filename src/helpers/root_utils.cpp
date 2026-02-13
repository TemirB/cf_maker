#include <helpers/root_utils.h>

#include <iostream>

std::string getCFName(int chIdx, int centIdx, int ktIdx) {
    return Form("CF_%d_%d_%d", chIdx, centIdx, ktIdx);
}

TH3D* getNum(TFile* f, int charge, int cent, int ktIdx) {
    TString name = Form("bp_%d_%d_num_%d", charge, cent, ktIdx);
    TH3D* h = (TH3D*) f->Get(name);

    if (!h) {
        std::cerr << "[getNum] NOT FOUND: " << name << std::endl;
        return nullptr;
    }
    return h;
};

TH3D* getNumWei(TFile* f, int charge, int cent, int ktIdx) {
    TString name = Form("bp_%d_%d_num_wei_%d", charge, cent, ktIdx);
    TH3D* h = (TH3D*) f->Get(name);

    if (!h) {
        std::cerr << "[getNumWei] NOT FOUND: " << name << std::endl;
        return nullptr;
    }
    return h;
};
