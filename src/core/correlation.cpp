#include "core/correlation.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include <TAxis.h>
#include <TFile.h>
#include <TH1.h>
#include <TNamed.h>

namespace
{

bool same_axis(const TAxis& left, const TAxis& right)
{
    if (left.GetNbins() != right.GetNbins()) {
        return false;
    }
    for (int bin = 1; bin <= left.GetNbins() + 1; ++bin) {
        if (left.GetBinLowEdge(bin) != right.GetBinLowEdge(bin)) {
            return false;
        }
    }
    return true;
}

bool same_binning(const TH1& left, const TH1& right)
{
    return left.GetDimension() == right.GetDimension() &&
           same_axis(*left.GetXaxis(), *right.GetXaxis()) &&
           (left.GetDimension() < 2 || same_axis(*left.GetYaxis(), *right.GetYaxis())) &&
           (left.GetDimension() < 3 || same_axis(*left.GetZaxis(), *right.GetZaxis()));
}

[[noreturn]] void invalid_bin(const TH1& numerator, int bin, const std::string& reason)
{
    throw std::runtime_error("invalid correlation moments in " + std::string(numerator.GetName()) +
                             " at bin " + std::to_string(bin) + ": " + reason);
}

} // namespace

CorrelationStatistics correlation_statistics(TFile& input)
{
    auto* object = input.Get("correlation_statistics");
    if (!object) {
        return CorrelationStatistics::PairWeights;
    }
    const auto* marker = dynamic_cast<TNamed*>(object);
    if (!marker) {
        throw std::runtime_error("correlation_statistics must be a TNamed input marker");
    }
    const std::string value = marker->GetTitle();
    if (value == "pair_weights") {
        return CorrelationStatistics::PairWeights;
    }
    if (value == "fixed_reference") {
        return CorrelationStatistics::FixedReference;
    }
    throw std::runtime_error("unknown correlation_statistics: " + value);
}

void fill_correlation(TH1& cf, const TH1& numerator, const TH1& denominator,
                      CorrelationStatistics statistics)
{
    if (&cf == &numerator || &cf == &denominator) {
        throw std::invalid_argument("correlation output must be distinct from its inputs");
    }
    validate_correlation_inputs(numerator, denominator, statistics);
    if (!same_binning(cf, numerator)) {
        throw std::invalid_argument("correlation histograms have incompatible axes");
    }

    cf.Reset("ICES");
    if (cf.GetSumw2N() == 0) {
        cf.Sumw2();
    }
    for (int bin = 0; bin < cf.GetNcells(); ++bin) {
        const double count = denominator.GetBinContent(bin);
        const double sum = numerator.GetBinContent(bin);
        const double numerator_error = numerator.GetBinError(bin);
        const double denominator_error = denominator.GetBinError(bin);
        if (!std::isfinite(count) || count < 0 || !std::isfinite(sum) ||
            !std::isfinite(numerator_error) || !std::isfinite(denominator_error)) {
            invalid_bin(numerator, bin, "non-finite or negative count/error");
        }
        if (count == 0) {
            if (sum != 0 || numerator_error != 0 || denominator_error != 0) {
                invalid_bin(numerator, bin, "weighted sum without reference pairs");
            }
            continue;
        }

        const double mean = sum / count;
        if (!std::isfinite(mean)) {
            invalid_bin(numerator, bin, "non-finite mean");
        }
        double error = 0;
        if (statistics == CorrelationStatistics::FixedReference) {
            if (denominator_error != 0) {
                invalid_bin(numerator, bin, "fixed reference must have zero error");
            }
            error = numerator_error / count;
        } else {
            const double count_tolerance =
                64 * std::numeric_limits<double>::epsilon() * std::max(1., count);
            if (std::abs(count - std::round(count)) > count_tolerance) {
                invalid_bin(numerator, bin, "reference is not an unweighted pair count");
            }
            const long double reference_variance =
                static_cast<long double>(denominator_error) * denominator_error;
            if (std::abs(reference_variance - count) > count_tolerance) {
                invalid_bin(numerator, bin, "reference variance differs from its pair count");
            }
            const long double sum_squares = numerator.GetSumw2()->At(bin);
            const long double square_of_sum = static_cast<long double>(sum) * sum / count;
            long double residual = sum_squares - square_of_sum;
            const long double cancellation_tolerance =
                64 * std::numeric_limits<double>::epsilon() * std::max(sum_squares, square_of_sum);
            if (residual < -cancellation_tolerance) {
                invalid_bin(numerator, bin,
                            "Sumw2 is smaller than sum(w)^2/N; check the input moments and "
                            "TH3F accumulation precision (use TH3D in the producer)");
            }
            if (std::abs(residual) <= cancellation_tolerance) {
                residual = 0;
            }
            if (count > 1) {
                error = static_cast<double>(std::sqrt(residual / count / (count - 1)));
            }
        }
        if (!std::isfinite(error)) {
            invalid_bin(numerator, bin, "non-finite error");
        }
        cf.SetBinContent(bin, mean);
        cf.SetBinError(bin, error);
    }
    cf.SetEntries(denominator.GetEntries());
}

void validate_correlation_inputs(const TH1& numerator, const TH1& denominator,
                                 CorrelationStatistics statistics)
{
    if (!same_binning(numerator, denominator)) {
        throw std::invalid_argument("correlation histograms have incompatible axes");
    }
    if (statistics == CorrelationStatistics::PairWeights && numerator.GetSumw2N() == 0) {
        throw std::invalid_argument("weighted-pair correlation requires numerator Sumw2");
    }
}
