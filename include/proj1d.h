#include <TH3.h>
#include <TF3.h>
#include <TDirectory.h>

void do_1d_proj(
    TFile* inputFile, TFile* outFile, 
    FitResult (&fitRes)[chargeSize][centralitySize][ktSize][rapiditySize]
);