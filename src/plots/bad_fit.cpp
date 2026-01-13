#include "plots.h"

#include <TTree.h>
#include <TH2D.h>

#include <helpers.h>

void MakeBadFitMaps(TFile* outFile, TTree* t)
{
    outFile->cd();

    TH2D* hCount = new TH2D(
        "hBadFitCount",
        "Bad fits count; k_{T} (GeV/c); y",
        ktSize, 0.1, 1.2,
        rapiditySize, -0.6, 0.0
    );

    TH2D* hChi2 = new TH2D(
        "hBadFitChi2Ndf",
        "<#chi^{2}/ndf> for bad fits; k_{T} (GeV/c); y",
        ktSize, 0.1, 1.2,
        rapiditySize, -0.6, 0.0
    );

    // ROOT сам усреднит chi2ndf в бинах
    t->Draw("yVal:ktVal>>hBadFitCount", "", "goff");
    t->Draw("yVal:ktVal>>hBadFitChi2Ndf", "chi2ndf", "goff");

    hCount->Write();
    hChi2->Write();
}
