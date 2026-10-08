#pragma once

class TFile;
class TH1;

// PairWeights describes the sthbtmaker input: N unweighted real pairs,
// S = sum(w), and the numerator's Sumw2 = sum(w*w) for those same pairs.
// FixedReference is an explicit alternative for synthetic data with an exact
// denominator and an independently fluctuating numerator.
enum class CorrelationStatistics
{
    PairWeights,
    FixedReference
};

// Copying TH3F into TH3D preserves values, but cannot restore the precision
// of their accumulation. Keep that provenance through subsequent projections.
enum class MomentStorage
{
    Double,
    Float
};

void set_moment_storage(TH1& histogram, MomentStorage storage);
[[nodiscard]] MomentStorage moment_storage(const TH1& histogram);

// An absent marker means PairWeights, the input contract of sthbtmaker.
// Other producers must write TNamed("correlation_statistics", "fixed_reference")
// to request the fixed-reference model; zero denominator errors do not select it.
[[nodiscard]] CorrelationStatistics correlation_statistics(TFile& input);

// Check the producer's axes and stored second moments before projecting them:
// ROOT can create default errors during projection when input Sumw2 is absent.
void validate_correlation_inputs(const TH1& numerator, const TH1& denominator,
                                 CorrelationStatistics statistics);

// Fills C = S/N. For PairWeights and N > 1, the squared error is the
// unbiased estimate of the variance of the mean:
//     (sum(w*w) - sum(w)*sum(w)/N) / (N*(N-1)).
// Empty bins have C = error = 0. Single-pair bins retain their mean but have
// error = 0: a variance cannot be estimated, and ROOT chi-square fits exclude
// them. Variances unresolved at the source's accumulation precision likewise
// have error = 0; for float sources this is reported in the log. Malformed
// moments or incompatible axes raise an exception instead of choosing a model.
// These are diagonal, independent-pair errors; event/track correlations need
// a separate event-level resampling analysis.
void fill_correlation(TH1& cf, const TH1& numerator, const TH1& denominator,
                      CorrelationStatistics statistics = CorrelationStatistics::PairWeights);
