#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include <TFile.h>
#include <TH3.h>

// Read the original storage types and moments without converting the input.
void inspect_pair_moments(const char* input_path, const char* sum_name = "bp_0_2_num_wei_0",
                          int bin = 1386)
{
    std::unique_ptr<TFile> input(TFile::Open(input_path, "READ"));
    if (!input || input->IsZombie()) {
        throw std::runtime_error("cannot open input ROOT file");
    }
    std::string count_name(sum_name);
    const std::string marker = "_num_wei_";
    const auto position = count_name.find(marker);
    if (position == std::string::npos) {
        throw std::invalid_argument("sum histogram name must contain _num_wei_");
    }
    count_name.replace(position, marker.size(), "_num_");
    auto* numerator = dynamic_cast<TH3*>(input->Get(sum_name));
    auto* denominator = dynamic_cast<TH3*>(input->Get(count_name.c_str()));
    if (!numerator || !denominator) {
        throw std::runtime_error("missing numerator or pair-count histogram");
    }
    if (bin < 0 || bin >= numerator->GetNcells() || bin >= denominator->GetNcells()) {
        throw std::out_of_range("global bin is outside the histograms");
    }
    int x = 0, y = 0, z = 0;
    numerator->GetBinXYZ(bin, x, y, z);
    const double count = denominator->GetBinContent(bin);
    const double sum = numerator->GetBinContent(bin);
    const double count_variance = std::pow(denominator->GetBinError(bin), 2);
    std::cout << std::setprecision(17) << "sum_hist=" << numerator->GetName()
              << " class=" << numerator->ClassName() << '\n'
              << "count_hist=" << denominator->GetName() << " class=" << denominator->ClassName()
              << '\n'
              << "global_bin=" << bin << " xyz=(" << x << ',' << y << ',' << z << ')' << " shape=("
              << numerator->GetNbinsX() << ',' << numerator->GetNbinsY() << ','
              << numerator->GetNbinsZ() << ")\n"
              << "N=" << count << " S=" << sum << " count_variance=" << count_variance
              << " Sumw2_size=" << numerator->GetSumw2N() << '\n';
    if (numerator->GetSumw2N() == 0) {
        std::cout << "Q is unavailable: numerator has no stored Sumw2\n";
        return;
    }
    const long double sum_squares = numerator->GetSumw2()->At(bin);
    std::cout << "Q=" << sum_squares << '\n';
    if (count > 0) {
        const long double square_of_sum = static_cast<long double>(sum) * sum / count;
        const long double tolerance =
            64 * std::numeric_limits<double>::epsilon() * std::max(sum_squares, square_of_sum);
        std::cout << "S^2/N=" << square_of_sum << " Q-S^2/N=" << sum_squares - square_of_sum
                  << " double_arithmetic_tolerance=" << tolerance << '\n';
    }
}
