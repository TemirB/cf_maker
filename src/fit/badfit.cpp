#include <fit/badfit.h>

#include "fit/types.h"
#include "helpers/constants.h"

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
    for (int kt = 0; kt < ktSize; kt++){
        const FitResult& res = fitRes[ch][cent][kt];
        if (!IsBadFit(res)) continue;

        BadFitPoint p{};
        p.charge = ch;
        p.cent   = cent;
        p.kt     = kt;

        p.chi2ndf = (res.ndf > 0 ? res.Chi2Ndf() : -1);
        p.pvalue  = res.pvalue;
        p.lambda  = res.lambda;
        p.elambda = res.elambda;

        for (int i = 0; i < 3; i++) {
            p.R[i]  = res.R[i];
            p.eR[i] = res.eR[i];
        }

        p.ktVal = 0.5 * (ktValues[kt] + ktValues[kt + 1]);

        badPoints.push_back(p);
    }
}

TTree* WriteBadFitTree(TFile* f, const std::vector<BadFitPoint>& badPoints) {
    f->cd();

    TTree* t = new TTree("badFits", "Bad fit points");
    BadFitPoint p;

    t->Branch("charge", &p.charge);
    t->Branch("cent", &p.cent);
    t->Branch("kt", &p.kt);

    t->Branch("chi2ndf", &p.chi2ndf);
    t->Branch("pvalue", &p.pvalue);
    t->Branch("lambda", &p.lambda);
    t->Branch("elambda", &p.elambda);

    t->Branch("R", p.R, "R[3]/D");
    t->Branch("eR", p.eR, "eR[3]/D");
    t->Branch("ktVal", &p.ktVal);

    for (const auto& b : badPoints) {
        p = b;
        t->Fill();
    }

    t->Write();
    return t;
}