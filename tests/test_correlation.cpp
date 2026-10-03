#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMemFile.h>
#include <TNamed.h>

#include "core/correlation.h"
#include "core/lcms.h"

namespace
{

void check_close(double actual, double expected, const std::string& message)
{
    if (!std::isfinite(actual) || std::abs(actual - expected) > 1e-12) {
        throw std::runtime_error(message + ": got " + std::to_string(actual) + ", expected " +
                                 std::to_string(expected));
    }
}

template <class Function> void check_throws(Function&& function, const std::string& message)
{
    try {
        function();
    } catch (const std::exception&) {
        return;
    }
    throw std::runtime_error(message);
}

void check_weighted_mean()
{
    TH1D denominator("den", "", 5, 0, 5);
    TH1D numerator("num", "", 5, 0, 5);
    TH1D cf("cf", "", 5, 0, 5);
    numerator.Sumw2();
    for (double weight : {1., 2.}) {
        denominator.Fill(0.5);
        numerator.Fill(0.5, weight);
    }
    for (double weight : {0., 2.}) {
        denominator.Fill(1.5);
        numerator.Fill(1.5, weight);
    }
    denominator.Fill(2.5);
    numerator.Fill(2.5, 0.7);
    for (double weight : {1., 1.}) {
        denominator.Fill(3.5);
        numerator.Fill(3.5, weight);
    }
    fill_correlation(cf, numerator, denominator);
    check_close(cf.GetBinContent(1), 1.5, "weighted mean");
    check_close(cf.GetBinError(1), 0.5, "finite-N error of [1,2]");
    check_close(cf.GetBinContent(2), 1., "unit mean with nonzero weight spread");
    check_close(cf.GetBinError(2), 1., "finite-N error of [0,2]");
    check_close(cf.GetBinContent(3), 0.7, "single-pair mean retained");
    check_close(cf.GetBinError(3), 0., "single-pair uncertainty not estimable");
    check_close(cf.GetBinContent(4), 1., "constant-weight mean");
    check_close(cf.GetBinError(4), 0., "constant-weight sample variance");
    check_close(cf.GetBinContent(5), 0., "empty-bin content");
    check_close(cf.GetBinError(5), 0., "empty-bin error");
}

void check_moment_validation()
{
    TH1D denominator("invalid_den", "", 1, 0, 1);
    TH1D numerator("invalid_num", "", 1, 0, 1);
    TH1D cf("invalid_cf", "", 1, 0, 1);
    denominator.SetBinContent(1, 2);
    numerator.SetBinContent(1, 2);
    check_throws([&] { fill_correlation(cf, numerator, denominator); }, "missing Sumw2 accepted");
    numerator.Sumw2();
    numerator.SetBinError(1, std::sqrt(2. - 1e-15));
    fill_correlation(cf, numerator, denominator);
    check_close(cf.GetBinError(1), 0., "roundoff at zero sample variance");
    numerator.SetBinError(1, 1);
    check_throws([&] { fill_correlation(cf, numerator, denominator); },
                 "impossible second moment accepted");
    numerator.SetBinError(1, std::sqrt(2.));
    denominator.SetBinContent(1, 2.5);
    check_throws([&] { fill_correlation(cf, numerator, denominator); },
                 "non-integral pair count accepted");
    denominator.SetBinContent(1, 2);
    denominator.SetBinError(1, 0);
    check_throws([&] { fill_correlation(cf, numerator, denominator); },
                 "fixed reference implicitly treated as pair data");
    denominator.SetBinError(1, std::sqrt(2.));
    numerator.SetBinContent(1, std::numeric_limits<double>::quiet_NaN());
    check_throws([&] { fill_correlation(cf, numerator, denominator); }, "non-finite mean accepted");
    numerator.SetBinContent(1, 2);
    TH1D wrong_axis("wrong_axis", "", 1, 0, 2);
    check_throws([&] { fill_correlation(wrong_axis, numerator, denominator); },
                 "incompatible axes accepted");
    check_throws([&] { fill_correlation(numerator, numerator, denominator); },
                 "input/output alias accepted");
    denominator.SetBinContent(1, 0);
    denominator.SetBinError(1, 0);
    check_throws([&] { fill_correlation(cf, numerator, denominator); },
                 "weighted sum without pairs accepted");
}

void check_fixed_reference()
{
    TMemFile input("statistics.root", "RECREATE");
    if (correlation_statistics(input) != CorrelationStatistics::PairWeights) {
        throw std::runtime_error("unmarked producer should use weighted pairs");
    }
    TNamed marker("correlation_statistics", "fixed_reference");
    marker.Write();
    const auto statistics = correlation_statistics(input);
    if (statistics != CorrelationStatistics::FixedReference) {
        throw std::runtime_error("fixed-reference marker not read");
    }
    TH1D denominator("fixed_den", "", 1, 0, 1);
    TH1D numerator("fixed_num", "", 1, 0, 1);
    TH1D cf("fixed_cf", "", 1, 0, 1);
    denominator.SetBinContent(1, 100);
    denominator.SetBinError(1, 0);
    numerator.SetBinContent(1, 150);
    numerator.SetBinError(1, std::sqrt(150.));
    fill_correlation(cf, numerator, denominator, statistics);
    check_close(cf.GetBinContent(1), 1.5, "fixed-reference ratio");
    check_close(cf.GetBinError(1), std::sqrt(150.) / 100., "fixed-reference ratio error");
    denominator.SetBinError(1, 10);
    check_throws([&] { fill_correlation(cf, numerator, denominator, statistics); },
                 "fluctuating denominator accepted as fixed");
    TMemFile invalid_input("invalid_statistics.root", "RECREATE");
    TNamed invalid_marker("correlation_statistics", "unsupported");
    invalid_marker.Write();
    check_throws([&] { static_cast<void>(correlation_statistics(invalid_input)); },
                 "unknown statistical marker accepted");
}

void check_projected_moments()
{
    TH3D denominator("project_den", "", 2, -1, 1, 2, -1, 1, 2, -1, 1);
    TH3D numerator("project_num", "", 2, -1, 1, 2, -1, 1, 2, -1, 1);
    numerator.Sumw2();
    auto fill = [&](double side, double longitudinal, double weight) {
        denominator.Fill(-0.5, side, longitudinal);
        numerator.Fill(-0.5, side, longitudinal, weight);
    };
    for (double weight : {0., 2.}) {
        fill(-0.5, -0.5, weight);
    }
    for (double weight : {1., 3.}) {
        fill(0.5, -0.5, weight);
    }
    fill(-0.5, 0.5, 1.);
    for (double weight : {2., 2.}) {
        fill(0.5, 0.5, weight);
    }

    std::unique_ptr<TH1D> den_1d(project_1d(denominator, LCMSAxis::Out, 1));
    std::unique_ptr<TH1D> num_1d(project_1d(numerator, LCMSAxis::Out, 1));
    TH1D cf_1d(*num_1d);
    fill_correlation(cf_1d, *num_1d, *den_1d);
    check_close(den_1d->GetBinContent(1), 7., "projected pair count");
    check_close(num_1d->GetBinContent(1), 11., "projected weighted sum");
    check_close(std::pow(num_1d->GetBinError(1), 2), 23., "projected Sumw2");
    check_close(cf_1d.GetBinContent(1), 11. / 7., "1D weighted mean");
    check_close(cf_1d.GetBinError(1), std::sqrt(20. / 147.), "1D finite-N error");

    std::unique_ptr<TH2D> den_2d(project_2d(denominator, LCMSAxis::Out, LCMSAxis::Side, 1));
    std::unique_ptr<TH2D> num_2d(project_2d(numerator, LCMSAxis::Out, LCMSAxis::Side, 1));
    TH2D cf_2d(*num_2d);
    fill_correlation(cf_2d, *num_2d, *den_2d);
    int populated_bins = 0;
    for (int bin = 0; bin < cf_2d.GetNcells(); ++bin) {
        if (den_2d->GetBinContent(bin) == 3) {
            ++populated_bins;
            check_close(cf_2d.GetBinContent(bin), 1., "2D first weighted mean");
            check_close(cf_2d.GetBinError(bin), std::sqrt(1. / 3.), "2D first finite-N error");
        } else if (den_2d->GetBinContent(bin) == 4) {
            ++populated_bins;
            check_close(cf_2d.GetBinContent(bin), 2., "2D second weighted mean");
            check_close(cf_2d.GetBinError(bin), std::sqrt(1. / 6.), "2D second finite-N error");
        }
    }
    if (populated_bins != 2) {
        throw std::runtime_error("2D projection did not retain both bins");
    }
}

} // namespace

int main()
{
    TH1::AddDirectory(false);
    check_weighted_mean();
    check_moment_validation();
    check_fixed_reference();
    check_projected_moments();
    std::cout << "Correlation statistics tests passed\n";
}
