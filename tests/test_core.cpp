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

} // namespace

int main()
{
    TestBinning();
    Testbin_center();
    TestFitResult();
    std::cout << "All tests passed\n";
    return 0;
}
