#include "plots.h"

#include <TH3D.h>

void ResetRanges(TH3D& h) {
    h.GetXaxis()->SetRange(0,0);
    h.GetYaxis()->SetRange(0,0);
    h.GetZaxis()->SetRange(0,0);
}