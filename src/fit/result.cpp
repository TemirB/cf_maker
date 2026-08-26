#include "fit/types.h"

#include <cmath>

bool FitResult::is_finite() const
{
    for (int i = 0; i < 3; i++) {
        if (!std::isfinite(R[i]) || !std::isfinite(eR[i])) {
            return false;
        }
    }
    return std::isfinite(lambda) && std::isfinite(elambda);
}

double FitResult::chi2_ndf() const
{
    return ndf > 0 ? chi2 / ndf : 0.0;
}

bool FitResult::is_valid() const
{
    const bool basicValid = ok && is_finite() && lambda > 0. && lambda < 1.;
    bool rValid = true;
    for (int i = 0; i < 3; i++) {
        if (R[i] > 9 || R[i] < 1.) {
            rValid = false;
            break;
        }
    }
    return basicValid && rValid;
}

bool is_bad_fit(const FitResult& r)
{
    if (!r.ok) {
        return true;
    }
    if (!r.is_finite()) {
        return true;
    }
    if (!r.is_valid()) {
        return true;
    }
    return false;
}
