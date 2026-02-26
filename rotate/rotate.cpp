#include <string>

#include <TFile.h>
#include <TMath.h>
#include <TH2.h>

void rotate() {

    TFile* input = new TFile("results/kt_run_896_fpc/2d.root", "READ");
    double angle = TMath::Pi()/4;
    double cos = TMath::Cos(angle);
    double sin = TMath::Sin(angle);

    TH2D* data = (TH2D*) input->Get("CF_0_0_2_side_long");
     // 


}