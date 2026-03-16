#include "plots.h"

#include <TH3D.h>

std::pair<int, int> GetBinRange(const TAxis* axis, double w) {
    int first = axis->FindFixBin(-w);
    int second = axis->FindFixBin(w);

    if (first > second) std::swap(first, second);
    
    return {std::max(1, first), std::min(axis->GetNbins(), second)};
}

const char* findRootProjectAxis(LCMSAxis axis) {
    return axis == LCMSAxis::Out  ? "x" :
           axis == LCMSAxis::Side ? "y" : "z";
}

void ResetRanges(TH3D& h) {
    h.GetXaxis()->SetRange(0,0);
    h.GetYaxis()->SetRange(0,0);
    h.GetZaxis()->SetRange(0,0);
}

void ProjectOnAxis(TH3D& h, LCMSAxis axis, double w) {
    if (axis != LCMSAxis::Out)  h.GetXaxis()->SetRangeUser(-w, w);
    if (axis != LCMSAxis::Side) h.GetYaxis()->SetRangeUser(-w, w);
    if (axis != LCMSAxis::Long) h.GetZaxis()->SetRangeUser(-w, w);
}

void FreezeAxis(TH3D& h, LCMSAxis freeze, double w) {
    if (freeze == LCMSAxis::Out)  h.GetXaxis()->SetRangeUser(-w, w);
    if (freeze == LCMSAxis::Side) h.GetYaxis()->SetRangeUser(-w, w);
    if (freeze == LCMSAxis::Long) h.GetZaxis()->SetRangeUser(-w, w);
}

void SetSlice1D(TH3D& h, LCMSAxis axis, double w) {
    ResetRanges(h);
    ProjectOnAxis(h, axis, w);
}