#include "fit/types.h"
#include <cmath>

bool FitResult::IsFinite() const {
    for(int i = 0; i < 6; i++)
        if(!std::isfinite(R[i]) || !std::isfinite(eR[i])) return false;
    return std::isfinite(lambda) && std::isfinite(elambda);
}

double FitResult::Chi2Ndf() const {
    return ndf>0 ? chi2/ndf : 0.0;
}

bool FitResult::IsValid() const {
    bool basic_valid = ok && IsFinite() && lambda > 0. && lambda < 1.;
    bool r_valid = true;
    for (int i = 0; i < 6; i++) {
        if (R[i] > 9.5 || R[i] < 2.) {
            r_valid = false;
            break;
        }
    }
    return basic_valid && r_valid;
}