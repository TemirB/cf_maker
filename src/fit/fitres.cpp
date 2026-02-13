#include "fit/types.h"
#include <cmath>

#include <iostream>

// void FitResult::returnAllR(int i) const {
//     std::cout << Form("R[%d]=%f+-%f\t", i, R[i], eR[i])
// }

bool FitResult::IsFinite() const {
    for(int i = 0; i < 3; i++)
        if(!std::isfinite(R[i]) || !std::isfinite(eR[i])) return false;
    return std::isfinite(lambda) && std::isfinite(elambda);
}

double FitResult::Chi2Ndf() const {
    return ndf>0 ? chi2/ndf : 0.0;
}

bool FitResult::IsValid() const {
    bool basic_valid = ok && IsFinite() && lambda > 0. && lambda < 1.;
    bool r_valid = true;
    for (int i = 0; i < 3; i++) {
        if (R[i] > 6.5 || R[i] < 2.) {
            r_valid = false;
            break;
        }
    }
    return basic_valid && r_valid;
}