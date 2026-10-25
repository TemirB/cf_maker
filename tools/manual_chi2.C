#include <algorithm>
#include <array>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <functional>
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
#include <TCanvas.h>
#include <TH3D.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TLatex.h>
#include <TH1.h>
#include <TStyle.h>

#include "../external/nlohmann_json/single_include/nlohmann/json.hpp"

namespace manual_chi2_detail
{
struct BinContribution
{
    int ix, iy, iz;
    double q_out, q_side, q_long, value, model, sigma;
    long double pull, chi2;
};

struct FitAudit
{
    int charge{}, centrality{}, bin{};
    std::string histogram_name, error;
    const TH3D* histogram = nullptr;
    std::vector<BinContribution> rows;
    long double chi2 = 0;
    int free_parameters = 0, saved_ndf = 0, status = -1, cov_status = -1, attempts = 0;
    double saved_chi2 = std::numeric_limits<double>::quiet_NaN();
    double root_chi2 = std::numeric_limits<double>::quiet_NaN();
    int root_ndf = -1;
    double lambda = std::numeric_limits<double>::quiet_NaN();
    bool ok = false, at_limit = false, comparison_ok = false;
    bool root_available = false, root_comparison_ok = false;

    long long ndf() const
    {
        return static_cast<long long>(rows.size()) - free_parameters;
    }

    long double root_over_manual() const
    {
        // The saved primary chi-square is not a substitute for ROOT's original
        // value in older files. An exact zero manual sum also has no ratio.
        if (!root_available || !std::isfinite(chi2) || chi2 <= 0 || root_chi2 < 0) {
            return std::numeric_limits<long double>::quiet_NaN();
        }
        return root_chi2 / chi2;
    }
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

nlohmann::json read_metadata(TFile& input)
{
    const auto* text = dynamic_cast<TObjString*>(input.Get("cf_maker_fit_results"));
    if (!text) {
        throw std::runtime_error("missing cf_maker_fit_results metadata; rerun stages.cf3d=true");
    }
    auto metadata = nlohmann::json::parse(text->GetString().Data());
    const auto& identity = metadata.at("identity");
    if (metadata.at("schema_version") != 1 || identity.at("model") != "lcms_gaussian_v1") {
        throw std::runtime_error("unsupported saved-fit schema or model");
    }
    const auto& config = identity.at("fit");
    (void)checked_options(config);
    const double q_max = config.at("q_max").get<double>();
    if (!std::isfinite(q_max) || q_max <= 0) {
        throw std::runtime_error("invalid saved fit range");
    }
    const auto& freeze = config.at("freeze");
    if (!freeze.is_array() || freeze.size() != 7) {
        throw std::runtime_error("invalid saved parameter freeze array");
    }
    const auto& names = identity.at("binning").at("names");
    if (!names.is_array() || names.empty()) {
        throw std::runtime_error("invalid saved bin names");
    }
    for (const auto& name : names) {
        (void)name.get<std::string>();
    }
    for (const auto& selection :
         {std::make_pair("charges", 2), std::make_pair("centralities", 4)}) {
        const auto& selected = identity.at("selection").at(selection.first);
        if (!selected.is_array() || selected.empty()) {
            throw std::runtime_error("invalid saved selection");
        }
        std::vector<int> seen;
        for (const auto& value : selected) {
            const int index = value.get<int>();
            if (index < 0 || index >= selection.second ||
                std::find(seen.begin(), seen.end(), index) != seen.end()) {
                throw std::runtime_error("invalid or duplicate saved selection index");
            }
            seen.push_back(index);
        }
    }
    return metadata;
}

FitAudit audit_fit(TFile& input, const nlohmann::json& metadata, int charge, int centrality,
                   int bin)
{
    FitAudit audit;
    audit.charge = charge;
    audit.centrality = centrality;
    audit.bin = bin;
    try {
        const auto& identity = metadata.at("identity");
        const auto& config = identity.at("fit");
        for (const auto& parameter : config.at("freeze")) {
            audit.free_parameters += parameter.is_null();
        }
        constexpr const char* kCharges[] = {"pos", "neg"};
        constexpr const char* kCentralities[] = {"0-10", "10-30", "30-50", "50-80"};
        audit.histogram_name = std::string("CF (charge=") + kCharges[charge] +
                               ", centrality=" + kCentralities[centrality] + ", " +
                               identity.at("input_type").get<std::string>() + "=" +
                               identity.at("binning").at("names").at(bin).get<std::string>() + ")";
        const auto& saved = metadata.at("fit_grid").at(charge).at(centrality).at(bin);
        audit.attempts = saved.at("attempts").get<int>();
        audit.status = saved.at("status").get<int>();
        audit.cov_status = saved.at("cov_status").get<int>();
        audit.ok = saved.at("ok").get<bool>();
        audit.at_limit = saved.at("at_limit").get<bool>();
        if (audit.attempts <= 0) {
            throw std::runtime_error("requested fit has no saved fit attempt");
        }
        const auto radii = saved.at("r").get<std::array<double, 6>>();
        audit.lambda = saved.at("lambda").get<double>();
        audit.saved_chi2 = saved.at("chi2").get<double>();
        audit.saved_ndf = saved.at("ndf").get<int>();
        if (saved.contains("chi2_comparison")) {
            const auto& comparison = saved.at("chi2_comparison");
            if (comparison.at("root_chi2").is_number()) {
                audit.root_chi2 = comparison.at("root_chi2").get<double>();
                audit.root_ndf = comparison.at("root_ndf").get<int>();
                audit.root_available = std::isfinite(audit.root_chi2);
            }
        }
        if (!std::isfinite(audit.lambda) || !std::isfinite(audit.saved_chi2) ||
            !std::all_of(radii.begin(), radii.end(), [](double v) { return std::isfinite(v); })) {
            throw std::runtime_error("saved model parameters or chi2 are non-finite");
        }
        audit.histogram = dynamic_cast<TH3D*>(input.Get(audit.histogram_name.c_str()));
        if (!audit.histogram) {
            throw std::runtime_error("missing saved TH3D CF");
        }
        const double q_max = config.at("q_max").get<double>();
        const bool use_range = checked_options(config).find('R') != std::string::npos;
        const auto x_range = selected_bins(*audit.histogram->GetXaxis(), use_range, q_max);
        const auto y_range = selected_bins(*audit.histogram->GetYaxis(), use_range, q_max);
        const auto z_range = selected_bins(*audit.histogram->GetZaxis(), use_range, q_max);
        for (int ix = x_range.first; ix <= x_range.second; ++ix) {
            const double q_out = audit.histogram->GetXaxis()->GetBinCenter(ix);
            for (int iy = y_range.first; iy <= y_range.second; ++iy) {
                const double q_side = audit.histogram->GetYaxis()->GetBinCenter(iy);
                for (int iz = z_range.first; iz <= z_range.second; ++iz) {
                    const double q_long = audit.histogram->GetZaxis()->GetBinCenter(iz);
                    const double value = audit.histogram->GetBinContent(ix, iy, iz);
                    const double sigma = audit.histogram->GetBinError(ix, iy, iz);
                    if (!std::isfinite(value) || !std::isfinite(sigma)) {
                        throw std::runtime_error("non-finite saved CF content or error");
                    }
                    if (sigma <= 0) {
                        continue;
                    }
                    const double q_r_q =
                        radii[0] * radii[0] * q_out * q_out +
                        radii[1] * radii[1] * q_side * q_side +
                        radii[2] * radii[2] * q_long * q_long + 2 * radii[3] * q_out * q_side +
                        2 * radii[4] * q_out * q_long + 2 * radii[5] * q_side * q_long;
                    constexpr double kHc2 = 0.197 * 0.197;
                    const double model = 1 + audit.lambda * std::exp(-q_r_q / kHc2);
                    const long double pull = (static_cast<long double>(value) - model) / sigma;
                    const long double contribution = pull * pull;
                    if (!std::isfinite(pull) || !std::isfinite(contribution)) {
                        throw std::runtime_error("non-finite model or chi2 contribution");
                    }
                    audit.rows.push_back({ix, iy, iz, q_out, q_side, q_long, value, model, sigma,
                                          pull, contribution});
                    audit.chi2 += contribution;
                }
            }
        }
        if (!std::isfinite(audit.chi2)) {
            throw std::runtime_error("non-finite total chi2");
        }
        const double tolerance = 1e-8 * std::max(1.0, std::abs(audit.saved_chi2));
        audit.comparison_ok = std::abs(audit.chi2 - audit.saved_chi2) <= tolerance &&
                              std::max(0LL, audit.ndf()) == audit.saved_ndf;
        if (audit.root_available) {
            const double root_tolerance = 1e-8 * std::max(1.0, std::abs(audit.root_chi2));
            audit.root_comparison_ok = std::abs(audit.chi2 - audit.root_chi2) <= root_tolerance &&
                                       std::max(0LL, audit.ndf()) == audit.root_ndf;
        }
    } catch (const std::exception& error) {
        audit.error = error.what();
        audit.rows.clear();
        audit.chi2 = std::numeric_limits<long double>::quiet_NaN();
    }
    return audit;
}

std::string csv_text(const std::string& value)
{
    std::string result = "\"";
    for (char character : value) {
        if (character == '"') {
            result += '"';
        }
        result += character;
    }
    return result + '"';
}

using CsvFile = std::unique_ptr<std::FILE, decltype(&std::fclose)>;

CsvFile open_csv(const std::filesystem::path& path)
{
    const int descriptor = ::open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0666);
    if (descriptor < 0) {
        throw std::runtime_error("cannot create output (existing files are preserved): " +
                                 path.string() + ": " + std::strerror(errno));
    }
    std::FILE* file = ::fdopen(descriptor, "w");
    if (!file) {
        ::close(descriptor);
        throw std::runtime_error("cannot open CSV stream: " + path.string());
    }
    return CsvFile(file, &std::fclose);
}

void finish_csv(CsvFile& file)
{
    const bool failed = std::ferror(file.get()) != 0;
    const int close_status = std::fclose(file.release());
    if (failed || close_status != 0) {
        throw std::runtime_error("cannot finish CSV output");
    }
}

class PdfOutput
{
  public:
    explicit PdfOutput(std::filesystem::path path)
        : path_(std::move(path)), canvas_("manual_chi2_audit", "Chi-square audit", 1800, 1100)
    {
        auto reservation = open_csv(path_);
        finish_csv(reservation);
        canvas_.SetCanvasSize(1800, 1100);
        canvas_.Print((path_.string() + "[").c_str());
    }

    ~PdfOutput()
    {
        if (open_) {
            canvas_.Print((path_.string() + "]").c_str());
        }
    }

    TCanvas& canvas()
    {
        return canvas_;
    }

    void finish()
    {
        canvas_.Print((path_.string() + "]").c_str());
        open_ = false;
        std::ifstream file(path_, std::ios::binary);
        std::array<char, 4> signature{};
        file.read(signature.data(), signature.size());
        if (!file || std::string(signature.data(), signature.size()) != "%PDF") {
            throw std::runtime_error("PDF output was not created: " + path_.string());
        }
    }

  private:
    std::filesystem::path path_;
    TCanvas canvas_;
    bool open_ = true;
};

long double largest_share(const FitAudit& audit, std::size_t count)
{
    if (!audit.error.empty() || audit.chi2 <= 0) {
        return std::numeric_limits<long double>::quiet_NaN();
    }
    std::vector<long double> largest;
    for (const auto& row : audit.rows) {
        largest.push_back(row.chi2);
    }
    count = std::min(count, largest.size());
    std::partial_sort(largest.begin(), largest.begin() + count, largest.end(),
                      std::greater<long double>());
    long double total = 0;
    for (std::size_t i = 0; i < count; ++i) {
        total += largest[i];
    }
    return 100 * total / audit.chi2;
}

void write_summary_row(std::FILE* output, const FitAudit& audit)
{
    const bool computed = audit.error.empty();
    const std::string error =
        audit.error.empty() && !audit.comparison_ok
            ? "manual chi2/ndf does not match saved fit statistics"
        : audit.error.empty() && audit.root_available && !audit.root_comparison_ok
            ? "manual chi2/ndf does not match ROOT fit statistics"
            : audit.error;
    std::fprintf(output, "%d,%d,%d,%s,%zu,%d,%.21Lg,", audit.charge, audit.centrality, audit.bin,
                 csv_text(audit.histogram_name).c_str(), audit.rows.size(), audit.free_parameters,
                 audit.chi2);
    if (computed) {
        std::fprintf(output, "%lld,", audit.ndf());
    } else {
        std::fprintf(output, "nan,");
    }
    const auto ratio = computed && audit.ndf() > 0 ? audit.chi2 / audit.ndf()
                                                   : std::numeric_limits<long double>::quiet_NaN();
    std::fprintf(output, "%.21Lg,%.17g,%d,%d,%d,%d,%d,%d,%d,%.17g,%.21Lg,%.21Lg,", ratio,
                 audit.saved_chi2, audit.saved_ndf, audit.comparison_ok, audit.status,
                 audit.cov_status, audit.ok, audit.at_limit, audit.attempts, audit.lambda,
                 largest_share(audit, 1), largest_share(audit, 5));
    if (audit.root_available) {
        std::fprintf(output, "%.17g,%d,%.21Lg,%d,", audit.root_chi2, audit.root_ndf,
                     audit.root_over_manual(), audit.root_comparison_ok);
    } else {
        std::fprintf(output, "nan,nan,nan,nan,");
    }
    std::fprintf(output, "%s\n", csv_text(error).c_str());
}

void write_detailed_rows(std::FILE* output, const FitAudit& audit)
{
    for (const auto& row : audit.rows) {
        std::fprintf(output,
                     "%d,%d,%d,%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                     "%.21Lg,%.21Lg\n",
                     audit.charge, audit.centrality, audit.bin,
                     audit.histogram->GetBin(row.ix, row.iy, row.iz), row.ix, row.iy, row.iz,
                     row.q_out, row.q_side, row.q_long, row.value, row.model, row.sigma, row.pull,
                     row.chi2);
    }
}

std::string audit_number(long double value, int significant_digits = 6)
{
    if (!std::isfinite(value)) {
        return "undefined";
    }
    std::array<char, 96> buffer{};
    std::snprintf(buffer.data(), buffer.size(), "%.*Lg", significant_digits, value);
    return buffer.data();
}

void draw_audit_text(double left, double bottom, double right, double top,
                     const std::vector<std::string>& lines, double text_size = 0.034, int color = 1)
{
    TPaveText text(left, bottom, right, top, "NDC");
    text.SetFillStyle(0);
    text.SetBorderSize(0);
    text.SetTextAlign(12);
    text.SetTextFont(42);
    text.SetTextSize(text_size);
    text.SetTextColor(color);
    for (const auto& line : lines) {
        text.AddText(line.c_str());
    }
    text.DrawClone();
}

std::vector<std::string> audit_wrapped_error(const std::string& error)
{
    constexpr std::size_t kWidth = 38;
    std::vector<std::string> lines;
    std::size_t start = 0;
    while (start < error.size()) {
        std::size_t end = std::min(start + kWidth, error.size());
        if (end < error.size()) {
            const auto space = error.rfind(' ', end);
            if (space != std::string::npos && space > start) {
                end = space;
            }
        }
        lines.push_back(error.substr(start, end - start));
        start = end;
        while (start < error.size() && error[start] == ' ') {
            ++start;
        }
    }
    return lines;
}

void draw_audit_unavailable(const FitAudit* audit, int centrality, bool cumulative)
{
    constexpr const char* kCentralities[] = {"0-10%", "10-30%", "30-50%", "50-80%"};
    std::vector<std::string> lines{std::string("Centrality ") + kCentralities[centrality],
                                   cumulative ? "Cumulative #chi^{2}: unavailable"
                                              : "Cell contributions: unavailable"};
    if (!audit) {
        lines.push_back("Not selected for this audit");
    } else {
        const auto explanation = audit_wrapped_error(audit->error);
        lines.insert(lines.end(), explanation.begin(), explanation.end());
    }
    draw_audit_text(0.08, 0.25, 0.96, 0.90, lines, 0.040, audit ? 632 : 1);
}

void draw_audit_contributions(const FitAudit& audit, int color)
{
    const auto count = audit.rows.size();
    int minimum_cell = audit.histogram->GetNcells();
    int maximum_cell = -1;
    TGraph graph;
    graph.SetName("manual_chi2_cell_contributions");
    graph.SetTitle(";ROOT global cell i;#chi^{2}_{i}");
    graph.SetLineColor(color);
    graph.SetLineWidth(1);
    graph.SetMarkerColor(color);
    graph.SetMarkerStyle(20);
    graph.SetMarkerSize(count < 1000 ? 0.7 : 0.25);
    double smallest_positive = std::numeric_limits<double>::infinity();
    double largest = 0;
    std::size_t zero_count = 0;
    for (const auto& row : audit.rows) {
        const int global_cell = audit.histogram->GetBin(row.ix, row.iy, row.iz);
        minimum_cell = std::min(minimum_cell, global_cell);
        maximum_cell = std::max(maximum_cell, global_cell);
        const double contribution = static_cast<double>(row.chi2);
        if (contribution > 0) {
            smallest_positive = std::min(smallest_positive, contribution);
            largest = std::max(largest, contribution);
            const int cell = audit.histogram->GetBin(row.ix, row.iy, row.iz);
            graph.SetPoint(graph.GetN(), cell, contribution);
        } else {
            ++zero_count;
        }
    }
    if (graph.GetN() == 0) {
        // ROOT cannot draw logarithmic axes for an identically zero sum.
        // Keep all zero contributions in the linear plot and in N_used.
        for (const auto& row : audit.rows) {
            graph.SetPoint(graph.GetN(), audit.histogram->GetBin(row.ix, row.iy, row.iz), 0);
        }
        graph.SetMinimum(0);
        graph.SetMaximum(1);
    } else {
        gPad->SetLogy();
        graph.SetMinimum(std::max(std::numeric_limits<double>::min(), smallest_positive * 0.5));
        graph.SetMaximum(largest * 2);
    }
    if (count > 0) {
        graph.GetXaxis()->SetLimits(std::max(0.0, static_cast<double>(minimum_cell) - 1),
                                    static_cast<double>(maximum_cell) + 1);
        graph.GetXaxis()->SetTitleSize(0.043);
        graph.GetYaxis()->SetTitleSize(0.043);
        graph.GetXaxis()->SetLabelSize(0.037);
        graph.GetYaxis()->SetLabelSize(0.037);
        graph.GetYaxis()->SetTitleOffset(1.55);
        graph.GetXaxis()->SetTitleOffset(1.25);
        graph.DrawClone("AP");
    } else {
        draw_audit_text(0.18, 0.27, 0.96, 0.57, {"No cells with usable errors"}, 0.040, 632);
    }

    const long long ndf = static_cast<long long>(count) - audit.free_parameters;
    const auto ratio = ndf > 0 ? audit_number(audit.chi2 / ndf) : "undefined";
    std::string flags = std::string("#lambda = ") + audit_number(audit.lambda);
    if (audit.at_limit) {
        flags += "  [parameter at limit]";
    }
    draw_audit_text(
        0.07, 0.685, 0.98, 0.985,
        {"#chi^{2} = " + audit_number(audit.chi2) + "; ndf = " + std::to_string(ndf),
         "#chi^{2}/ndf = " + ratio,
         "#chi^{2}_{ROOT}/#chi^{2}_{manual} = " + (std::isfinite(audit.root_over_manual())
                                                       ? audit_number(audit.root_over_manual())
                                                       : "undefined / unavailable"),
         "N_{used} = " + std::to_string(count) +
             "; N_{free} = " + std::to_string(audit.free_parameters),
         "status = " + std::to_string(audit.status) + "; cov = " +
             std::to_string(audit.cov_status) + "; attempts = " + std::to_string(audit.attempts),
         flags,
         std::string("Source fit: ") +
             ((audit.ok && !audit.at_limit && audit.saved_ndf > 0) ? "usable" : "UNUSABLE") +
             (audit.comparison_ok ? "; saved sum: MATCH" : "; saved sum: MISMATCH")},
        0.032);
    if (zero_count > 0 && largest > 0) {
        draw_audit_text(0.18, 0.645, 0.98, 0.683,
                        {std::to_string(zero_count) + " zero contributions omitted on log y"},
                        0.026);
    }
}

void draw_audit_cumulative(const FitAudit& audit, int color)
{
    if (!(audit.chi2 > 0)) {
        draw_audit_text(0.10, 0.30, 0.96, 0.85,
                        {"Cumulative fraction is undefined", "because total #chi^{2} = 0.",
                         "All cells remain included in N_{used}."},
                        0.040);
        return;
    }
    std::vector<const BinContribution*> ordered;
    ordered.reserve(audit.rows.size());
    for (const auto& row : audit.rows) {
        ordered.push_back(&row);
    }
    std::sort(ordered.begin(), ordered.end(), [&audit](const auto* left, const auto* right) {
        if (left->chi2 != right->chi2) {
            return left->chi2 > right->chi2;
        }
        return audit.histogram->GetBin(left->ix, left->iy, left->iz) <
               audit.histogram->GetBin(right->ix, right->iy, right->iz);
    });
    TGraph graph;
    graph.SetName("manual_chi2_cumulative_contributions");
    graph.SetTitle(";Number of cells k, largest #chi^{2}_{i} first;Accumulated #chi^{2} [%]");
    graph.SetLineColor(color);
    graph.SetLineWidth(2);
    graph.SetMarkerColor(color);
    graph.SetMarkerStyle(1);
    long double cumulative = 0;
    long double top_five = 0;
    for (std::size_t index = 0; index < ordered.size(); ++index) {
        cumulative += ordered[index]->chi2;
        if (index < 5) {
            top_five += ordered[index]->chi2;
        }
        graph.SetPoint(graph.GetN(), index + 1, static_cast<double>(100 * cumulative / audit.chi2));
    }
    graph.SetMinimum(0);
    graph.SetMaximum(105);
    graph.GetXaxis()->SetLimits(0.9, std::max(1.1, ordered.size() * 1.03));
    graph.GetXaxis()->SetTitleSize(0.043);
    graph.GetYaxis()->SetTitleSize(0.043);
    graph.GetXaxis()->SetLabelSize(0.037);
    graph.GetYaxis()->SetLabelSize(0.037);
    graph.GetYaxis()->SetTitleOffset(1.55);
    graph.GetXaxis()->SetTitleOffset(1.25);
    gPad->SetLogx(ordered.size() > 1);
    graph.DrawClone("AL");
    const double maximum_rank = std::max(1.1, ordered.size() * 1.03);
    for (const double level : {50.0, 90.0}) {
        TLine reference(0.9, level, maximum_rank, level);
        reference.SetLineStyle(3);
        reference.SetLineColor(920);
        reference.DrawClone();
    }
    draw_audit_text(0.08, 0.82, 0.98, 0.985,
                    {"Top " + std::to_string(std::min<std::size_t>(5, ordered.size())) +
                         " cells: " + audit_number(100 * top_five / audit.chi2, 5) +
                         "% of #chi^{2}",
                     "100 #times #Sigma_{j=1}^{k} #chi^{2}_{(j)} / #chi^{2}",
                     "All used cells included; descending contribution"},
                    0.032);
}

void draw_audit_page(TCanvas& canvas, const std::array<const FitAudit*, 4>& audits, int charge,
                     int bin, const nlohmann::json& identity, const char* pdf_path)
{
    canvas.Clear();
    canvas.SetCanvasSize(1800, 1100);
    canvas.SetFillColor(0);
    constexpr const char* kCentralities[] = {"0-10%", "10-30%", "30-50%", "50-80%"};
    constexpr int kColors[] = {601, 807, 418, 617};
    std::string bin_label = "bin " + std::to_string(bin);
    const auto& names = identity.at("binning").at("names");
    if (bin >= 0 && static_cast<std::size_t>(bin) < names.size()) {
        bin_label += " (" + identity.at("input_type").get<std::string>() + " = " +
                     names.at(bin).get<std::string>() + ")";
    }
    draw_audit_text(0.02, 0.966, 0.98, 0.998,
                    {"Manual #chi^{2} audit: charge " + std::to_string(charge) +
                     (charge == 0 ? " (positive), " : " (negative), ") + bin_label},
                    0.022);
    draw_audit_text(0.02, 0.935, 0.98, 0.965,
                    {"#chi^{2}_{i} = [(C_{2,i} - C_{2,i}^{fit}) / #sigma_{C_{2,i}}]^{2}; "
                     "#chi^{2} = #Sigma_{i} #chi^{2}_{i}; ndf = N_{used} - N_{free}. "
                     "Saved fit, no minimization."},
                    0.018);
    for (int centrality = 0; centrality < 4; ++centrality) {
        for (int row = 0; row < 2; ++row) {
            canvas.cd();
            const std::string name =
                "chi2_audit_pad_" + std::to_string(centrality) + "_" + std::to_string(row);
            const double left = 0.002 + 0.25 * centrality;
            const double right = 0.248 + 0.25 * centrality;
            const double bottom = row == 0 ? 0.488 : 0.010;
            const double top = row == 0 ? 0.932 : 0.480;
            auto* pad = new TPad(name.c_str(), "", left, bottom, right, top);
            pad->SetFillColor(0);
            pad->SetLeftMargin(0.19);
            pad->SetRightMargin(0.04);
            pad->SetBottomMargin(0.16);
            pad->SetTopMargin(row == 0 ? 0.36 : 0.21);
            pad->SetGrid();
            pad->Draw();
            pad->cd();
            const auto* audit = audits[centrality];
            if (!audit || !audit->error.empty()) {
                draw_audit_unavailable(audit, centrality, row == 1);
            } else if (row == 0) {
                draw_audit_contributions(*audit, kColors[centrality]);
            } else {
                draw_audit_cumulative(*audit, kColors[centrality]);
            }
            if (row == 0) {
                canvas.cd();
                TLatex heading;
                heading.SetNDC();
                heading.SetTextAlign(22);
                heading.SetTextFont(62);
                heading.SetTextSize(0.018);
                heading.DrawLatex(0.125 + 0.25 * centrality, 0.934,
                                  (std::string("Centrality ") + kCentralities[centrality]).c_str());
            }
        }
    }
    canvas.cd();
    canvas.Modified();
    canvas.Update();
    canvas.Print(pdf_path);
    canvas.Clear();
}

void run_batch(TFile& input, const nlohmann::json& metadata, const std::vector<int>& charges,
               const std::vector<int>& centralities, const std::vector<int>& bins,
               const char* detailed_path, const char* requested_pdf)
{
    std::filesystem::path pdf_path =
        requested_pdf && *requested_pdf
            ? requested_pdf
            : std::filesystem::path(input.GetName()).parent_path() /
                  (std::filesystem::path(input.GetName()).stem().string() + "_chi2.pdf");
    if (pdf_path.extension() != ".pdf") {
        throw std::invalid_argument("PDF output path must end in .pdf");
    }
    auto summary_path = pdf_path;
    summary_path.replace_extension(".csv");
    std::vector<std::filesystem::path> outputs{pdf_path, summary_path};
    if (detailed_path && *detailed_path) {
        outputs.emplace_back(detailed_path);
    }
    for (std::size_t i = 0; i < outputs.size(); ++i) {
        outputs[i] = std::filesystem::absolute(outputs[i]).lexically_normal();
        if (std::filesystem::exists(outputs[i])) {
            throw std::runtime_error("output exists; choose a new PDF/CSV path: " +
                                     outputs[i].string());
        }
        for (std::size_t j = 0; j < i; ++j) {
            if (outputs[i] == outputs[j]) {
                throw std::runtime_error("output paths must be distinct");
            }
        }
        std::filesystem::create_directories(outputs[i].parent_path());
    }
    auto summary = open_csv(summary_path);
    std::fprintf(
        summary.get(),
        "charge,centrality,bin,histogram,used_cells,free_parameters,chi2,ndf,"
        "chi2_ndf,saved_chi2,saved_ndf,comparison_ok,status,cov_status,ok,"
        "at_limit,attempts,lambda,top1_share_percent,top5_share_percent,root_chi2,root_ndf,"
        "root_over_manual,root_comparison_ok,error\n");
    CsvFile detailed(nullptr, &std::fclose);
    if (detailed_path && *detailed_path) {
        detailed = open_csv(detailed_path);
        std::fprintf(detailed.get(),
                     "charge,centrality,bin,cell,ix,iy,iz,q_out,q_side,q_long,C,model,"
                     "sigma,pull,chi2_contribution\n");
    }
    PdfOutput pdf(pdf_path);
    std::size_t processed = 0, errors = 0, source_unusable = 0;
    for (int charge : charges) {
        for (int bin : bins) {
            std::array<FitAudit, 4> page;
            std::array<const FitAudit*, 4> panels{};
            for (int centrality : centralities) {
                auto& audit = page[centrality];
                audit = audit_fit(input, metadata, charge, centrality, bin);
                panels[centrality] = &audit;
                ++processed;
                errors += !audit.error.empty() || !audit.comparison_ok ||
                          (audit.root_available && !audit.root_comparison_ok);
                source_unusable +=
                    !audit.ok || audit.at_limit || audit.attempts <= 0 || audit.saved_ndf <= 0;
                write_summary_row(summary.get(), audit);
                if (detailed) {
                    write_detailed_rows(detailed.get(), audit);
                }
                std::cout << "audit (ch=" << charge << ", centr=" << centrality << ", b=" << bin
                          << "): ";
                if (audit.error.empty()) {
                    print_summary(audit.chi2, audit.rows.size(), audit.free_parameters);
                    std::cout << "  status=" << audit.status << " covariance=" << audit.cov_status
                              << " at_limit=" << audit.at_limit
                              << " comparison=" << (audit.comparison_ok ? "OK" : "MISMATCH")
                              << " ROOT/manual=" << audit.root_over_manual()
                              << " ROOT_ndf=" << audit.root_ndf << '\n';
                } else {
                    std::cout << "ERROR: " << audit.error << '\n';
                }
            }
            draw_audit_page(pdf.canvas(), panels, charge, bin, metadata.at("identity"),
                            pdf_path.string().c_str());
        }
    }
    pdf.finish();
    finish_csv(summary);
    if (detailed) {
        finish_csv(detailed);
    }
    std::cout << "PDF=" << pdf_path << "\nsummary CSV=" << summary_path << '\n';
    std::cout << "processed=" << processed << " audit_errors=" << errors
              << " source_unusable=" << source_unusable << '\n';
    if (errors > 0) {
        throw std::runtime_error("some histograms could not be checked or did not match; "
                                 "see PDF and summary CSV (remaining histograms were processed)");
    }
}
} // namespace manual_chi2_detail

// All selected CFs, four centralities per PDF page, plus a summary CSV:
// root -l -b -q 'tools/manual_chi2.C("output/cf3d.root")'
// One CF and its per-cell CSV (the original explicit-index interface):
// root -l -b -q 'tools/manual_chi2.C("output/cf3d.root",0,0,0,"/tmp/chi2.csv")'
// Filter by charge/bin; -1 means all saved selections. Last argument sets the PDF path.
void manual_chi2(const char* cf_path = "", int charge = -1, int centrality = -1, int bin = -1,
                 const char* csv_path = "", const char* pdf_path = "")
{
    using namespace manual_chi2_detail;
    std::cout << std::setprecision(std::numeric_limits<long double>::max_digits10);
    if (!cf_path || !*cf_path) {
        std::cout << std::setprecision(12);
        std::vector<BinContribution> rows;
        long double chi2 = 0;
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
    gROOT->SetBatch(kTRUE);
    gStyle->SetGridColor(kGray + 1);
    gStyle->SetGridStyle(3);
    if (charge < -1 || charge >= 2 || centrality < -1 || centrality >= 4 || bin < -1) {
        throw std::invalid_argument("charge must be -1..1, centrality -1..3, and bin >= -1");
    }
    std::unique_ptr<TFile> input(TFile::Open(cf_path, "READ"));
    if (!input || input->IsZombie()) {
        throw std::runtime_error(std::string("cannot open ROOT file: ") + cf_path);
    }
    const auto metadata = read_metadata(*input);
    const auto& identity = metadata.at("identity");
    auto charges = identity.at("selection").at("charges").get<std::vector<int>>();
    auto centralities = identity.at("selection").at("centralities").get<std::vector<int>>();
    const auto filter = [](std::vector<int>& selected, int index) {
        if (index == -1) {
            return;
        }
        if (std::find(selected.begin(), selected.end(), index) == selected.end()) {
            throw std::invalid_argument("requested charge/centrality was not selected in this run");
        }
        selected = {index};
    };
    filter(charges, charge);
    filter(centralities, centrality);
    const auto bin_count = identity.at("binning").at("names").size();
    if (bin >= 0 && static_cast<std::size_t>(bin) >= bin_count) {
        throw std::out_of_range("bin is outside saved binning");
    }
    std::vector<int> bins;
    for (std::size_t index = 0; index < bin_count; ++index) {
        if (bin == -1 || static_cast<std::size_t>(bin) == index) {
            bins.push_back(static_cast<int>(index));
        }
    }
    if (charge == -1 || centrality == -1 || bin == -1 || (pdf_path && *pdf_path)) {
        run_batch(*input, metadata, charges, centralities, bins, csv_path, pdf_path);
        return;
    }
    const auto audit = audit_fit(*input, metadata, charge, centrality, bin);
    if (!audit.error.empty()) {
        throw std::runtime_error(audit.histogram_name + ": " + audit.error);
    }
    std::cout << "CF=" << audit.histogram_name << '\n'
              << "fit status=" << audit.status << " covariance=" << audit.cov_status
              << " ok=" << audit.ok << " at_limit=" << audit.at_limit << '\n';
    print_summary(audit.chi2, audit.rows.size(), audit.free_parameters);
    std::cout << "saved chi2=" << audit.saved_chi2 << " ndf=" << audit.saved_ndf
              << " delta_chi2=" << audit.chi2 - audit.saved_chi2 << '\n';
    if (audit.root_available) {
        std::cout << "ROOT chi2=" << audit.root_chi2 << " ndf=" << audit.root_ndf
                  << " ROOT/manual=" << audit.root_over_manual()
                  << " comparison=" << (audit.root_comparison_ok ? "OK" : "MISMATCH") << '\n';
    } else {
        std::cout << "ROOT/manual=unavailable (no original ROOT statistics were saved)\n";
    }
    print_largest_contributions(*audit.histogram, audit.rows, audit.chi2);
    write_csv(csv_path, audit.rows);
    if (!audit.comparison_ok || (audit.root_available && !audit.root_comparison_ok)) {
        throw std::runtime_error("manual chi2/ndf does not match saved fit statistics");
    }
    std::cout << "comparison=OK\n";
}
