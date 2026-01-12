#include "fit/result.h"
#include <cmath>

bool FitResult::IsFinite() const {
    for(int i=0;i<3;i++)
        if(!std::isfinite(R[i]) || !std::isfinite(eR[i])) return false;
    return std::isfinite(lambda) && std::isfinite(elambda);
}

double FitResult::Chi2Ndf() const {
    return ndf>0 ? chi2/ndf : 0.0;
}

bool FitResult::IsValid() const {
    return ok && IsFinite() && ndf>0;
}
