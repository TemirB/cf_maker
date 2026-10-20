#include <algorithm>
#include <array>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <unistd.h>

#include <TFile.h>
#include <TH3D.h>
#include <TObjString.h>

#include "../external/nlohmann_json/single_include/nlohmann/json.hpp"

namespace manual_chi2_detail
{
struct BinContribution
{
    int ix, iy, iz;
    double q_out, q_side, q_long, value, model, sigma;
    long double pull, chi2;
};

std::pair<int, int> selected_bins(const TAxis& axis, bool use_range, double q_max)
{
    int first = std::max(1, axis.GetFirst());
    int last = std::min(axis.GetNbins(), axis.GetLast());
    if (use_range) {
        // ROOT 6.40 HFitInterface::ExamineRange: when only one bin remains,
        // ROOT keeps it even if its center is outside the requested range.
        first = std::min(std::max(axis.FindFixBin(-q_max), first), last + 1);
        last = std::max(std::min(axis.FindFixBin(q_max), last), first - 1);
        if (first < last) {
            if (axis.GetBinCenter(first) < -q_max) {
                ++first;
            }
            if (axis.GetBinCenter(last) > q_max) {
                --last;
            }
        }
    }
    return {first, last};
}

std::string checked_options(const nlohmann::json& fit)
{
    std::string options = fit.at("options").get<std::string>();
    std::transform(options.begin(), options.end(), options.begin(), [](unsigned char character) {
        return static_cast<char>(std::toupper(character));
    });
    for (const std::string modifier : {"MULTITHREAD", "MULTIPROCESS", "SERIAL"}) {
        std::size_t position;
        while ((position = options.find(modifier)) != std::string::npos) {
            options.erase(position, modifier.size());
        }
    }
    if (fit.at("use_integral").get<bool>() || options.find_first_of("IWLPU") != std::string::npos ||
        options.find("WIDTH") != std::string::npos) {
        throw std::runtime_error(
            "manual_chi2 supports center-evaluated chi-square with saved bin errors only; "
            "I/use_integral, W/WW, L, P, U and WIDTH are unsupported");
    }
    return options;
}

void write_csv(const char* path, const std::vector<BinContribution>& rows)
{
    if (!path || !*path) {
        return;
    }
    const int descriptor = ::open(path, O_WRONLY | O_CREAT | O_EXCL, 0666);
    if (descriptor < 0) {
        throw std::runtime_error(std::string("cannot create CSV (existing files are preserved): ") +
                                 path + ": " + std::strerror(errno));
    }
    std::FILE* output = ::fdopen(descriptor, "w");
    if (!output) {
        ::close(descriptor);
        throw std::runtime_error("cannot open newly created CSV stream");
    }
    bool failed = std::fprintf(output, "ix,iy,iz,q_out,q_side,q_long,C,model,sigma,pull,"
                                       "chi2_contribution\n") < 0;
    for (const auto& row : rows) {
        if (std::fprintf(output, "%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.21Lg,%.21Lg\n",
                         row.ix, row.iy, row.iz, row.q_out, row.q_side, row.q_long, row.value,
                         row.model, row.sigma, row.pull, row.chi2) < 0) {
            failed = true;
            break;
        }
    }
    if (std::fclose(output) != 0) {
        failed = true;
    }
    if (failed) {
        throw std::runtime_error(std::string("cannot finish CSV: ") + path);
    }
    std::cout << "CSV=" << path << '\n';
}

void print_summary(long double chi2, std::size_t points, int free_parameters)
{
    const long long ndf = static_cast<long long>(points) - free_parameters;
    std::cout << "npoints=" << points << " nfree=" << free_parameters << " chi2=" << chi2
              << " ndf=" << points << '-' << free_parameters << '=' << ndf << " chi2/ndf=";
    if (ndf > 0) {
        std::cout << chi2 / ndf;
    } else {
        std::cout << "undefined (ndf <= 0)";
    }
    std::cout << '\n';
}

void print_largest_contributions(const TH3D& histogram, const std::vector<BinContribution>& rows,
                                 long double chi2)
{
    std::vector<const BinContribution*> ordered;
    ordered.reserve(rows.size());
    for (const auto& row : rows) {
        ordered.push_back(&row);
    }
    const auto count = std::min<std::size_t>(5, ordered.size());
    std::partial_sort(ordered.begin(), ordered.begin() + count, ordered.end(),
                      [](const BinContribution* left, const BinContribution* right) {
                          return left->chi2 > right->chi2;
                      });
    const auto precision = std::cout.precision();
    std::cout << std::setprecision(12) << "largest_chi2_contributions:\n";
    long double largest_sum = 0;
    for (std::size_t i = 0; i < count; ++i) {
        const auto& row = *ordered[i];
        largest_sum += row.chi2;
        std::cout << "cell=" << histogram.GetBin(row.ix, row.iy, row.iz) << " xyz=(" << row.ix
                  << ',' << row.iy << ',' << row.iz << ") q=(" << row.q_out << ',' << row.q_side
                  << ',' << row.q_long << ") C2=" << row.value << " C2_fit=" << row.model
                  << " sigma_C2=" << row.sigma << " pull=" << row.pull
                  << " contribution=" << row.chi2
                  << " share_percent=" << (chi2 > 0 ? 100 * row.chi2 / chi2 : 0) << '\n';
    }
    std::cout << "top" << count << "_share_percent=" << (chi2 > 0 ? 100 * largest_sum / chi2 : 0)
              << '\n';
    std::cout.precision(precision);
}
} // namespace manual_chi2_detail

// root -l -b -q 'tools/manual_chi2.C()'
// root -l -b -q 'tools/manual_chi2.C("output/cf3d.root",0,0,0,"/tmp/chi2.csv")'
// Uses full saved fit parameters and CF errors; does not refit or reconstruct errors.
void manual_chi2(const char* cf_path = "", int charge = 0, int centrality = 0, int bin = 0,
                 const char* csv_path = "")
{
    using namespace manual_chi2_detail;
    std::cout << std::setprecision(std::numeric_limits<long double>::max_digits10);
    std::vector<BinContribution> rows;
    long double chi2 = 0;
    if (!cf_path || !*cf_path) {
        std::cout << std::setprecision(12);
        for (int point = 0; point < 5; ++point) {
            const double value = 1.0 + 0.1 * point;
            const long double pull = (static_cast<long double>(value) - 1.2) / 0.1;
            rows.push_back(
                {point + 1, 1, 1, double(point), 0, 0, value, 1.2, 0.1, pull, pull * pull});
            chi2 += pull * pull;
            std::cout << "bin=" << point + 1 << " C=" << value
                      << " model=1.2 sigma=0.1 pull=" << pull << " contribution=" << pull * pull
                      << '\n';
        }
        print_summary(chi2, rows.size(), 1);
        write_csv(csv_path, rows);
        return;
    }

    if (charge < 0 || charge >= 2 || centrality < 0 || centrality >= 4 || bin < 0) {
        throw std::invalid_argument("charge must be 0..1, centrality 0..3, and bin >= 0");
    }
    std::unique_ptr<TFile> input(TFile::Open(cf_path, "READ"));
    if (!input || input->IsZombie()) {
        throw std::runtime_error(std::string("cannot open ROOT file: ") + cf_path);
    }
    const auto* text = dynamic_cast<TObjString*>(input->Get("cf_maker_fit_results"));
    if (!text) {
        throw std::runtime_error("missing cf_maker_fit_results metadata; rerun stages.cf3d=true");
    }
    const auto metadata = nlohmann::json::parse(text->GetString().Data());
    const auto& identity = metadata.at("identity");
    if (metadata.at("schema_version") != 1 || identity.at("model") != "lcms_gaussian_v1") {
        throw std::runtime_error("unsupported saved-fit schema or model");
    }
    for (const auto selection :
         {std::make_pair("charges", charge), std::make_pair("centralities", centrality)}) {
        const auto& selected = identity.at("selection").at(selection.first);
        if (std::find(selected.begin(), selected.end(), selection.second) == selected.end()) {
            throw std::invalid_argument("requested charge/centrality was not selected in this run");
        }
    }
    const auto& fit_config = identity.at("fit");
    const std::string options = checked_options(fit_config);
    const double q_max = fit_config.at("q_max").get<double>();
    if (!std::isfinite(q_max) || q_max <= 0) {
        throw std::runtime_error("invalid saved fit range");
    }
    const auto& freeze = fit_config.at("freeze");
    if (!freeze.is_array() || freeze.size() != 7) {
        throw std::runtime_error("invalid saved parameter freeze array");
    }
    int free_parameters = 0;
    for (const auto& parameter : freeze) {
        free_parameters += parameter.is_null();
    }
    const auto& names = identity.at("binning").at("names");
    if (!names.is_array() || static_cast<std::size_t>(bin) >= names.size()) {
        throw std::out_of_range("bin is outside saved binning");
    }
    const auto& saved = metadata.at("fit_grid").at(charge).at(centrality).at(bin);
    if (saved.at("attempts").get<int>() <= 0) {
        throw std::runtime_error("requested fit has no saved fit attempt");
    }
    const auto radii = saved.at("r").get<std::array<double, 6>>();
    const double lambda = saved.at("lambda").get<double>();
    const double saved_chi2 = saved.at("chi2").get<double>();
    const int saved_ndf = saved.at("ndf").get<int>();
    if (!std::isfinite(lambda) || !std::isfinite(saved_chi2) ||
        !std::all_of(radii.begin(), radii.end(),
                     [](double value) { return std::isfinite(value); })) {
        throw std::runtime_error("saved model parameters or chi2 are non-finite");
    }
    constexpr const char* kCharges[] = {"pos", "neg"};
    constexpr const char* kCentralities[] = {"0-10", "10-30", "30-50", "50-80"};
    const std::string histogram_name = std::string("CF (charge=") + kCharges[charge] +
                                       ", centrality=" + kCentralities[centrality] + ", " +
                                       identity.at("input_type").get<std::string>() + "=" +
                                       names.at(bin).get<std::string>() + ")";
    const auto* histogram = dynamic_cast<TH3D*>(input->Get(histogram_name.c_str()));
    if (!histogram) {
        throw std::runtime_error("missing saved TH3D CF: " + histogram_name);
    }
    const bool use_range = options.find('R') != std::string::npos;
    const auto x_range = selected_bins(*histogram->GetXaxis(), use_range, q_max);
    const auto y_range = selected_bins(*histogram->GetYaxis(), use_range, q_max);
    const auto z_range = selected_bins(*histogram->GetZaxis(), use_range, q_max);
    for (int ix = x_range.first; ix <= x_range.second; ++ix) {
        const double q_out = histogram->GetXaxis()->GetBinCenter(ix);
        for (int iy = y_range.first; iy <= y_range.second; ++iy) {
            const double q_side = histogram->GetYaxis()->GetBinCenter(iy);
            for (int iz = z_range.first; iz <= z_range.second; ++iz) {
                const double q_long = histogram->GetZaxis()->GetBinCenter(iz);
                const double value = histogram->GetBinContent(ix, iy, iz);
                const double sigma = histogram->GetBinError(ix, iy, iz);
                if (!std::isfinite(value) || !std::isfinite(sigma)) {
                    throw std::runtime_error("non-finite saved CF content or error");
                }
                if (sigma <= 0) {
                    continue;
                }
                const double q_r_q =
                    radii[0] * radii[0] * q_out * q_out + radii[1] * radii[1] * q_side * q_side +
                    radii[2] * radii[2] * q_long * q_long + 2 * radii[3] * q_out * q_side +
                    2 * radii[4] * q_out * q_long + 2 * radii[5] * q_side * q_long;
                constexpr double kHc2 = 0.197 * 0.197;
                const double model = 1 + lambda * std::exp(-q_r_q / kHc2);
                const long double pull = (static_cast<long double>(value) - model) / sigma;
                if (!std::isfinite(pull)) {
                    throw std::runtime_error("non-finite model or pull");
                }
                rows.push_back(
                    {ix, iy, iz, q_out, q_side, q_long, value, model, sigma, pull, pull * pull});
                chi2 += pull * pull;
            }
        }
    }
    std::cout << "CF=" << histogram_name << '\n'
              << "fit status=" << saved.at("status") << " covariance=" << saved.at("cov_status")
              << " ok=" << saved.at("ok") << " at_limit=" << saved.at("at_limit") << '\n';
    print_summary(chi2, rows.size(), free_parameters);
    std::cout << "saved chi2=" << saved_chi2 << " ndf=" << saved_ndf
              << " delta_chi2=" << chi2 - saved_chi2 << '\n';
    print_largest_contributions(*histogram, rows, chi2);
    write_csv(csv_path, rows);
    const long long ndf = static_cast<long long>(rows.size()) - free_parameters;
    const double tolerance = 1e-8 * std::max(1.0, std::abs(saved_chi2));
    if (std::abs(chi2 - saved_chi2) > tolerance || std::max(0LL, ndf) != saved_ndf) {
        throw std::runtime_error("manual chi2/ndf does not match saved fit statistics");
    }
    std::cout << "comparison=OK\n";
}
