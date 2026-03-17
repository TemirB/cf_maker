#include "helpers.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include <TH3.h>
#include <TTree.h>
#include <TFile.h>

#include "fit/types.h"

TH3D* getNum(TFile* f, int charge, int cent, int yIdx) {
    TString name = Form("bp_%d_%d_num_%d", charge, cent, yIdx);
    TH3D* h = (TH3D*) f->Get(name);

    if (!h) {
        std::cerr << "[getNum] NOT FOUND: " << name << std::endl;
        return nullptr;
    }
    return h;
};

TH3D* getNumWei(TFile* f, int charge, int cent, int yIdx) {
    TString name = Form("bp_%d_%d_num_wei_%d", charge, cent, yIdx);
    TH3D* h = (TH3D*) f->Get(name);

    if (!h) {
        std::cerr << "[getNumWei] NOT FOUND: " << name << std::endl;
        return nullptr;
    }
    return h;
};

std::string getCFName(int chIdx, int centIdx, int yIdx) {
    double left = rapidityValues[0] + step*yIdx;
    double right = left + step;
    return Form("CF at charge=%s, centrality=%s, y=[%.2f,%.2f]", chargeNames[chIdx], centralityNames[centIdx], left, right);
}

bool IsBadFit(const FitResult& r) {
    if (!r.ok) return true;
    if (!r.IsFinite()) return true;
    if (!r.IsValid()) return true;
    // if (r.ndf <= 0) return true;
    // if (r.Chi2Ndf() > 3.0) return true;
    return false;
}

void CollectBadFits(
    FitGrid& fitRes,
    std::vector<BadFitPoint>& badPoints
) {
    for (int ch = 0; ch < chargeSize; ch++)
    for (int cent = 0; cent < centralitySize; cent++)
    for (int y = 0; y < rapiditySize; y++){
        const FitResult& res = fitRes[ch][cent][y];
        if (!IsBadFit(res)) continue;

        BadFitPoint p{};
        p.charge = ch;
        p.cent   = cent;
        p.y     = y;

        p.chi2ndf = (res.ndf > 0 ? res.Chi2Ndf() : -1);
        p.pvalue  = res.pvalue;
        p.lambda  = res.lambda;
        p.elambda = res.elambda;

        for (int i = 0; i < 3; i++) {
            p.R[i]  = res.R[i];
            p.eR[i] = res.eR[i];
        }

        p.yVal = 0.5 * (rapidityValues[y] + rapidityValues[y + step]);

        badPoints.push_back(p);
    }
}

TTree* WriteBadFitTree(TFile* f, const std::vector<BadFitPoint>& badPoints) {
    f->cd();

    TTree* t = new TTree("badFits", "Bad fit points");
    BadFitPoint p;

    t->Branch("charge", &p.charge);
    t->Branch("cent", &p.cent);
    t->Branch("y", &p.y);

    t->Branch("chi2ndf", &p.chi2ndf);
    t->Branch("pvalue", &p.pvalue);
    t->Branch("lambda", &p.lambda);
    t->Branch("elambda", &p.elambda);

    t->Branch("R", p.R, "R[3]/D");
    t->Branch("eR", p.eR, "eR[3]/D");
    t->Branch("yVal", &p.yVal);

    for (const auto& b : badPoints) {
        p = b;
        t->Fill();
    }

    t->Write();
    return t;
}
