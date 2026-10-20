#include "fit/types.h"

#include <cmath>

bool FitResult::is_finite() const
{
    for (std::size_t i = 0; i < r.size(); ++i) {
        if (!std::isfinite(r[i]) || !std::isfinite(e_r[i])) {
            return false;
        }
    }
    for (const double value : corr) {
        if (!std::isfinite(value)) {
            return false;
        }
    }
    return std::isfinite(lambda) && std::isfinite(e_lambda) && std::isfinite(chi2) &&
           std::isfinite(p_value);
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
        if (r[i] > 9 || r[i] < 1.) {
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

bool is_usable_fit(const FitResult& result)
{
    return result.ok && result.is_finite() && !result.at_limit && result.ndf > 0;
}

bool is_better_fit(const FitResult& candidate, const FitResult& current)
{
    const bool candidate_usable = is_usable_fit(candidate);
    const bool current_usable = is_usable_fit(current);
    if (candidate_usable != current_usable) {
        return candidate_usable;
    }

    const bool candidate_converged = candidate.ok && candidate.is_finite() && candidate.ndf > 0;
    const bool current_converged = current.ok && current.is_finite() && current.ndf > 0;
    return candidate_converged && (!current_converged || candidate.chi2 < current.chi2);
}
