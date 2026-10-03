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

} // namespace

int main()
{
    TestBinning();
    Testbin_center();
    TestFitResult();
    test_fit_result_finiteness();
    std::cout << "All tests passed\n";
    return 0;
}
