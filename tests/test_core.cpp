#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

#include "core/binning.h"
#include "fit/types.h"

namespace
{

#define CHECK(cond)                                                                                \
    do {                                                                                           \
        if (!(cond)) {                                                                             \
            std::cerr << "FAILED: " #cond " at " << __FILE__ << ":" << __LINE__ << "\n";           \
            std::exit(1);                                                                          \
        }                                                                                          \
    } while (0)

#define CHECK_NEAR(a, b, eps)                                                                      \
    do {                                                                                           \
        if (std::abs((a) - (b)) > (eps)) {                                                         \
            std::cerr << "FAILED: " #a " = " << (a) << ", " #b " = " << (b) << " at " << __FILE__  \
                      << ":" << __LINE__ << "\n";                                                  \
            std::exit(1);                                                                          \
        }                                                                                          \
    } while (0)

void TestBinning()
{
    CHECK(charge::kCount == 2);
    CHECK(centrality::kCount == 4);
    CHECK(kt::kCount == 4);
    CHECK(rapidity::kCount == 10);
    CHECK(lcms::kCount == 6);

    CHECK(std::string(charge::kNames[0]) == "pos");
    CHECK(std::string(charge::kNames[1]) == "neg");
    CHECK(std::string(centrality::kNames[0]) == "0-10");
    CHECK(std::string(kt::kNames[0]) == "0.15-0.25");
    CHECK(std::string(rapidity::kFileNames[9]) == "0.8_1.0");
    CHECK(std::string(lcms::kNames[3]) == "out-side");

    CHECK(kt::kValues.size() == static_cast<std::size_t>(kt::kCount + 1));
    CHECK(rapidity::kValues.size() == static_cast<std::size_t>(rapidity::kCount + 1));
}

void Testbin_center()
{
    Bin bin;
    bin.count = 2;
    bin.values = {0.0, 2.0, 4.0};

    CHECK_NEAR(bin_center(bin, 0), 1.0, 1e-12);
    CHECK_NEAR(bin_center(bin, 1), 3.0, 1e-12);
}

void TestFitResult()
{
    FitResult r;
    CHECK(r.is_finite());
    CHECK(!r.is_valid());
    CHECK(is_bad_fit(r));

    r.ndf = 0;
    CHECK_NEAR(r.chi2_ndf(), 0.0, 1e-12);

    r.chi2 = 8.0;
    r.ndf = 4;
    CHECK_NEAR(r.chi2_ndf(), 2.0, 1e-12);

    r.ok = true;
    r.lambda = 0.5;
    r.r = {4.0, 4.0, 4.0, 0.0, 0.0, 0.0};
    CHECK(r.is_valid());
    CHECK(!is_bad_fit(r));

    r.r[0] = 10.0;
    CHECK(!r.is_valid());
    r.r[0] = 4.0;

    r.lambda = 1.5;
    CHECK(!r.is_valid());
    r.lambda = 0.5;

    r.r[0] = std::numeric_limits<double>::infinity();
    CHECK(!r.is_finite());
    CHECK(!r.is_valid());
    r.r[0] = 4.0;

    r.ok = false;
    CHECK(is_bad_fit(r));
}

void check_non_finite_field(FitResult& result, double& field)
{
    const double original = field;
    for (const double bad_value :
         {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
          -std::numeric_limits<double>::infinity()}) {
        field = bad_value;
        CHECK(!result.is_finite());
        CHECK(!result.is_valid());
        CHECK(is_bad_fit(result));
    }
    field = original;
    CHECK(result.is_finite());
    CHECK(result.is_valid());
}

void test_fit_result_finiteness()
{
    FitResult result;
    result.ok = true;
    result.r = {4.0, 5.0, 6.0, -0.5, 1.5, 0.5};
    result.e_r = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    result.lambda = 0.7;
    result.e_lambda = 0.01;
    result.chi2 = 12.0;
    result.ndf = 10;
    result.p_value = 0.28;
    result.corr.fill(0.1);
    CHECK(result.is_finite());
    CHECK(result.is_valid());

    for (double& value : result.r) {
        check_non_finite_field(result, value);
    }
    for (double& value : result.e_r) {
        check_non_finite_field(result, value);
    }
    check_non_finite_field(result, result.lambda);
    check_non_finite_field(result, result.e_lambda);
    check_non_finite_field(result, result.chi2);
    check_non_finite_field(result, result.p_value);
    for (double& value : result.corr) {
        check_non_finite_field(result, value);
    }
}

FitResult usable_fit(double chi2)
{
    FitResult result;
    result.ok = true;
    result.r = {4.0, 5.0, 6.0, -0.5, 1.5, 0.5};
    result.e_r = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    result.lambda = 0.7;
    result.e_lambda = 0.01;
    result.chi2 = chi2;
    result.ndf = 10;
    result.p_value = 0.28;
    result.status = 0;
    result.cov_status = 3;
    result.attempts = 1;
    result.corr.fill(0.1);
    return result;
}

void test_fit_candidate_ranking()
{
    const FitResult usable = usable_fit(120.0);
    FitResult boundary = usable_fit(80.0);
    boundary.at_limit = true;
    FitResult failed = usable_fit(1.0);
    failed.ok = false;
    failed.status = 3;

    CHECK(is_usable_fit(usable));
    CHECK(!is_usable_fit(boundary));
    CHECK(!is_usable_fit(failed));

    // A usable retry must win even if the fit on a parameter limit has lower chi-square.
    CHECK(is_better_fit(usable, boundary));
    CHECK(!is_better_fit(boundary, usable));
    CHECK(is_better_fit(boundary, failed));
    CHECK(!is_better_fit(failed, boundary));
    CHECK(!is_better_fit(failed, usable));

    FitResult lower_chi2 = usable_fit(100.0);
    CHECK(is_better_fit(lower_chi2, usable));
    CHECK(!is_better_fit(usable, lower_chi2));
    CHECK(!is_better_fit(usable, usable));

    lower_chi2.at_limit = true;
    FitResult higher_boundary = boundary;
    higher_boundary.chi2 = 110.0;
    CHECK(is_better_fit(lower_chi2, higher_boundary));
    CHECK(!is_better_fit(higher_boundary, lower_chi2));
    CHECK(!is_better_fit(boundary, boundary));

    FitResult failed_retry = failed;
    failed_retry.chi2 = 0.0;
    CHECK(!is_better_fit(failed_retry, failed));
    CHECK(!is_better_fit(failed, failed_retry));

    for (const int ndf : {0, -1}) {
        FitResult no_ndf = usable_fit(0.0);
        no_ndf.ndf = ndf;
        CHECK(!is_usable_fit(no_ndf));
        CHECK(!is_better_fit(no_ndf, usable));
        CHECK(!is_better_fit(no_ndf, boundary));
        CHECK(!is_better_fit(no_ndf, failed));
        CHECK(is_better_fit(usable, no_ndf));
        CHECK(is_better_fit(boundary, no_ndf));
    }
}

void check_non_finite_candidate(FitResult& candidate, double& field)
{
    const FitResult usable = usable_fit(120.0);
    FitResult boundary = usable;
    boundary.at_limit = true;
    FitResult failed = usable;
    failed.ok = false;
    const double original = field;
    for (const double bad_value :
         {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
          -std::numeric_limits<double>::infinity()}) {
        field = bad_value;
        CHECK(!is_usable_fit(candidate));
        CHECK(!is_better_fit(candidate, usable));
        CHECK(!is_better_fit(candidate, boundary));
        CHECK(!is_better_fit(candidate, failed));
        CHECK(is_better_fit(usable, candidate));
        CHECK(is_better_fit(boundary, candidate));
    }
    field = original;
    CHECK(is_usable_fit(candidate));
}

void test_fit_candidate_finiteness()
{
    FitResult candidate = usable_fit(1.0);
    for (double& value : candidate.r) {
        check_non_finite_candidate(candidate, value);
    }
    for (double& value : candidate.e_r) {
        check_non_finite_candidate(candidate, value);
    }
    check_non_finite_candidate(candidate, candidate.lambda);
    check_non_finite_candidate(candidate, candidate.e_lambda);
    check_non_finite_candidate(candidate, candidate.chi2);
    check_non_finite_candidate(candidate, candidate.p_value);
    for (double& value : candidate.corr) {
        check_non_finite_candidate(candidate, value);
    }
}

} // namespace

int main()
{
    TestBinning();
    Testbin_center();
    TestFitResult();
    test_fit_result_finiteness();
    test_fit_candidate_ranking();
    test_fit_candidate_finiteness();
    std::cout << "All tests passed\n";
    return 0;
}
