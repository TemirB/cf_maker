#include "plots.h"

#include <TH3D.h>

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