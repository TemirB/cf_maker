#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TH3D.h>
#include <TH3F.h>
#include <TKey.h>
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

void check_float_weighted_mean()
{
    TH3F denominator("float_den", "", 7, 0, 7, 1, 0, 1, 1, 0, 1);
    TH3F numerator("float_num", "", 7, 0, 7, 1, 0, 1, 1, 0, 1);
    TH3D cf("float_cf", "", 7, 0, 7, 1, 0, 1, 1, 0, 1);
    denominator.Sumw2();
    numerator.Sumw2();
    const auto fill = [&](double x, double weight) {
        denominator.Fill(x, .5, .5, 1.);
        numerator.Fill(x, .5, .5, weight);
    };
    fill(.5, 1.2);
    for (const double weight : {1.2, 1.2}) {
        fill(1.5, weight);
    }
    for (const double weight : {1.3, 1.3}) {
        fill(2.5, weight);
    }
    for (const double weight : {1.1, 1.10001}) {
        fill(3.5, weight);
    }
    for (const double weight : {1.3, 1.30001}) {
        fill(4.5, weight);
    }
    for (const double weight : {0., 2.}) {
        fill(5.5, weight);
    }
    for (const double weight : {1., 1.2}) {
        fill(6.5, weight);
    }
    const auto residual = [&](int x) {
        const int bin = numerator.GetBin(x, 1, 1);
        const long double sum = numerator.GetBinContent(bin);
        return numerator.GetSumw2()->At(bin) - sum * sum / denominator.GetBinContent(bin);
    };
    if (residual(4) >= 0 || residual(5) <= 0) {
        throw std::runtime_error("Float fixture does not exercise both rounding directions");
    }
    fill_correlation(cf, numerator, denominator);
    for (int x = 1; x <= 5; ++x) {
        check_close(cf.GetBinContent(x, 1, 1),
                    numerator.GetBinContent(x, 1, 1) / denominator.GetBinContent(x, 1, 1),
                    "float-resolution masking must retain the measured mean");
        check_close(cf.GetBinError(x, 1, 1), 0., "unresolved float variance must be masked");
    }
    check_close(cf.GetBinContent(6, 1, 1), 1., "resolved float weighted mean");
    check_close(cf.GetBinError(6, 1, 1), 1., "resolved float finite-N variance");
    if (std::abs(cf.GetBinError(7, 1, 1) - .1) > 1e-6) {
        throw std::runtime_error("Unequal float weights lost their usable variance");
    }

    // Promoting storage does not recover the lost precision: the source tag
    // must accompany the converted moments before their variance is computed.
    TH3D converted_denominator;
    TH3D converted_numerator;
    denominator.Copy(converted_denominator);
    numerator.Copy(converted_numerator);
    set_moment_storage(converted_numerator, MomentStorage::Float);
    fill_correlation(cf, converted_numerator, converted_denominator);
    check_close(cf.GetBinError(5, 1, 1), 0., "converted float variance must remain masked");
    check_close(cf.GetBinError(6, 1, 1), 1., "converted resolved float variance");

    converted_numerator.SetBinError(2, 1, 1, 1);
    check_throws([&] { fill_correlation(cf, converted_numerator, converted_denominator); },
                 "grossly invalid float second moment was masked");
    converted_numerator.SetBinError(2, 1, 1, numerator.GetBinError(2, 1, 1));
    converted_numerator.SetBinError(1, 1, 1, 0);
    check_throws([&] { fill_correlation(cf, converted_numerator, converted_denominator); },
                 "grossly invalid single-pair second moment was masked");
    converted_numerator.SetBinError(1, 1, 1, numerator.GetBinError(1, 1, 1));
    converted_denominator.SetBinError(2, 1, 1, 0);
    check_throws([&] { fill_correlation(cf, converted_numerator, converted_denominator); },
                 "float provenance relaxed the unweighted-count variance check");
}

void check_moment_storage()
{
    TH1F one("float_1d_storage", "", 1, 0, 1);
    TH2F two("float_2d_storage", "", 1, 0, 1, 1, 0, 1);
    TH3F three("float_3d_storage", "", 1, 0, 1, 1, 0, 1, 1, 0, 1);
    for (const TH1* histogram :
         {static_cast<TH1*>(&one), static_cast<TH1*>(&two), static_cast<TH1*>(&three)}) {
        if (moment_storage(*histogram) != MomentStorage::Float) {
            throw std::runtime_error("Native float histogram storage was not recognized");
        }
    }
    TH3D converted("converted_storage", "", 1, 0, 1, 1, 0, 1, 1, 0, 1);
    if (moment_storage(converted) != MomentStorage::Double) {
        throw std::runtime_error("Native double histogram received float provenance");
    }
    set_moment_storage(converted, MomentStorage::Float);
    TH3D copy(converted);
    // ROOT deliberately omits attached functions when copying histograms;
    // conversion and projection callers must transfer provenance explicitly.
    set_moment_storage(copy, moment_storage(converted));
    if (moment_storage(copy) != MomentStorage::Float) {
        throw std::runtime_error("Explicit histogram copy lost source storage provenance");
    }
    TMemFile input("storage_roundtrip.root", "RECREATE");
    converted.Write();
    std::unique_ptr<TObject> readback(input.GetKey(converted.GetName())->ReadObj());
    auto* histogram = dynamic_cast<TH1*>(readback.get());
    if (!histogram || moment_storage(*histogram) != MomentStorage::Float) {
        throw std::runtime_error("ROOT serialization lost source storage provenance");
    }
    histogram->SetDirectory(nullptr);
    set_moment_storage(copy, MomentStorage::Double);
    if (moment_storage(copy) != MomentStorage::Double ||
        moment_storage(converted) != MomentStorage::Float) {
        throw std::runtime_error("Changing provenance affected an independently owned copy");
    }
}

void check_double_accumulation()
{
    TH1D denominator("long_double_den", "", 1, 0, 1);
    TH1D numerator("long_double_num", "", 1, 0, 1);
    TH1D cf("long_double_cf", "", 1, 0, 1);
    numerator.Sumw2();
    for (int pair = 0; pair < 10000; ++pair) {
        denominator.Fill(.5);
        numerator.Fill(.5, 1.1);
    }
    fill_correlation(cf, numerator, denominator);
    check_close(cf.GetBinContent(1), 1.1, "mean after long double-storage accumulation");
    check_close(cf.GetBinError(1), 0., "constant double weights acquired a numerical variance");
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
    check_float_weighted_mean();
    check_moment_storage();
    check_double_accumulation();
    check_moment_validation();
    check_fixed_reference();
    check_projected_moments();
    std::cout << "Correlation statistics tests passed\n";
}
