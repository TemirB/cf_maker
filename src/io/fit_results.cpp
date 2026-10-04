#include "io/fit_results.h"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include <TDirectory.h>
#include <TFile.h>
#include <TH3D.h>
#include <TKey.h>
#include <TObjString.h>
#include <nlohmann/json.hpp>
#include <openssl/evp.h>

#include "core/correlation.h"
#include "io/input.h"
#include "io/output.h"

namespace
{
constexpr int kSchemaVersion = 1;
constexpr const char* kMetadataName = "cf_maker_fit_results";

std::string input_sha256(const std::string& path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("cannot fingerprint input file: " + path);
    }
    const std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> context(EVP_MD_CTX_new(),
                                                                          &EVP_MD_CTX_free);
    if (!context || EVP_DigestInit_ex(context.get(), EVP_sha256(), nullptr) != 1) {
        throw std::runtime_error("cannot initialize input SHA256");
    }
    std::array<char, 65536> buffer{};
    while (input.read(buffer.data(), static_cast<std::streamsize>(buffer.size())) ||
           input.gcount() > 0) {
        if (EVP_DigestUpdate(context.get(), buffer.data(),
                             static_cast<std::size_t>(input.gcount())) != 1) {
            throw std::runtime_error("cannot fingerprint input file: " + path);
        }
    }
    if (!input.eof()) {
        throw std::runtime_error("cannot read input file for SHA256: " + path);
    }
    std::array<unsigned char, EVP_MAX_MD_SIZE> digest{};
    unsigned int length = 0;
    if (EVP_DigestFinal_ex(context.get(), digest.data(), &length) != 1) {
        throw std::runtime_error("cannot finish input SHA256");
    }
    std::ostringstream result;
    result << std::hex << std::setfill('0');
    for (unsigned int i = 0; i < length; ++i) {
        result << std::setw(2) << static_cast<unsigned int>(digest[i]);
    }
    return result.str();
}

nlohmann::json fit_identity(const Config& cfg)
{
    TDirectory::TContext context;
    TFile input(cfg.input.file.c_str(), "READ");
    if (input.IsZombie()) {
        throw std::runtime_error("cannot read input statistics mode: " + cfg.input.file);
    }
    const auto statistics = correlation_statistics(input);
    const char* statistics_version = statistics == CorrelationStatistics::PairWeights
                                         ? "weighted-mean-source-roundoff-v2"
                                         : "fixed_reference_v1";
    nlohmann::json freeze = nlohmann::json::array();
    for (const auto& value : cfg.fit.freeze) {
        freeze.push_back(value ? nlohmann::json(*value) : nlohmann::json(nullptr));
    }
    // Increment the model version when the formula, units or parameter conventions change.
    return {{"model", "lcms_gaussian_v1"},
            {"statistics", statistics_version},
            {"retry_selection", "usable-first-v1"},
            {"input_type", cfg.input.type},
            {"binning", {{"values", cfg.binning.values}, {"names", cfg.binning.names}}},
            {"selection",
             {{"charges", cfg.selection.charges}, {"centralities", cfg.selection.centralities}}},
            {"fit",
             {{"q_max", cfg.fit.q_max},
              {"use_default_ip", cfg.fit.use_default_ip},
              {"radius_sq_min", cfg.fit.radius_sq_min},
              {"radius_sq_max", cfg.fit.radius_sq_max},
              {"cross_min", cfg.fit.cross_min},
              {"cross_max", cfg.fit.cross_max},
              {"lambda_min", cfg.fit.lambda_min},
              {"lambda_max", cfg.fit.lambda_max},
              {"options", cfg.fit.options},
              {"minimizer", cfg.fit.minimizer},
              {"retry_with_defaults", cfg.fit.retry_with_defaults},
              {"use_integral", cfg.fit.use_integral},
              {"minos_errors", cfg.fit.minos_errors},
              {"freeze", freeze}}}};
}

nlohmann::json encode_double(double value)
{
    if (std::isnan(value)) {
        return "nan";
    }
    if (std::isinf(value)) {
        return value > 0 ? "+inf" : "-inf";
    }
    return value;
}

double decode_double(const nlohmann::json& value)
{
    if (value.is_number()) {
        return value.get<double>();
    }
    if (value == "nan") {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (value == "+inf") {
        return std::numeric_limits<double>::infinity();
    }
    if (value == "-inf") {
        return -std::numeric_limits<double>::infinity();
    }
    throw std::runtime_error("invalid fit-result floating-point value");
}

template <std::size_t N> nlohmann::json encode_array(const std::array<double, N>& values)
{
    nlohmann::json result = nlohmann::json::array();
    for (const double value : values) {
        result.push_back(encode_double(value));
    }
    return result;
}

template <std::size_t N>
void decode_array(const nlohmann::json& values, std::array<double, N>& result)
{
    if (!values.is_array() || values.size() != N) {
        throw std::runtime_error("invalid fit-result array size");
    }
    for (std::size_t i = 0; i < N; ++i) {
        result[i] = decode_double(values[i]);
    }
}

nlohmann::json encode_fit(const FitResult& fit)
{
    return {{"r", encode_array(fit.r)},
            {"e_r", encode_array(fit.e_r)},
            {"lambda", encode_double(fit.lambda)},
            {"e_lambda", encode_double(fit.e_lambda)},
            {"chi2", encode_double(fit.chi2)},
            {"ndf", fit.ndf},
            {"p_value", encode_double(fit.p_value)},
            {"ok", fit.ok},
            {"status", fit.status},
            {"cov_status", fit.cov_status},
            {"at_limit", fit.at_limit},
            {"attempts", fit.attempts},
            {"corr", encode_array(fit.corr)}};
}

FitResult decode_fit(const nlohmann::json& value)
{
    FitResult fit;
    decode_array(value.at("r"), fit.r);
    decode_array(value.at("e_r"), fit.e_r);
    fit.lambda = decode_double(value.at("lambda"));
    fit.e_lambda = decode_double(value.at("e_lambda"));
    fit.chi2 = decode_double(value.at("chi2"));
    fit.ndf = value.at("ndf").get<int>();
    fit.p_value = decode_double(value.at("p_value"));
    fit.ok = value.at("ok").get<bool>();
    fit.status = value.at("status").get<int>();
    fit.cov_status = value.at("cov_status").get<int>();
    fit.at_limit = value.at("at_limit").get<bool>();
    fit.attempts = value.at("attempts").get<int>();
    decode_array(value.at("corr"), fit.corr);
    if (fit.attempts < 0 || (fit.ok && fit.attempts == 0)) {
        throw std::runtime_error("invalid fit-result attempt count");
    }
    return fit;
}

void check_grid_size(const Config& cfg, const FitGrid& grid)
{
    if (grid.size() != charge::kCount) {
        throw std::runtime_error("invalid saved fit-grid charge count");
    }
    for (const auto& centralities : grid) {
        if (centralities.size() != centrality::kCount) {
            throw std::runtime_error("invalid saved fit-grid centrality count");
        }
        for (const auto& bins : centralities) {
            if (bins.size() != static_cast<std::size_t>(cfg.binning.count)) {
                throw std::runtime_error("invalid saved fit-grid bin count");
            }
        }
    }
}
} // namespace

void write_fit_results(TFile& file, const Config& cfg)
{
    check_grid_size(cfg, cfg.fit_results);
    nlohmann::json grid = nlohmann::json::array();
    for (const auto& centralities : cfg.fit_results) {
        nlohmann::json charges = nlohmann::json::array();
        for (const auto& bins : centralities) {
            nlohmann::json results = nlohmann::json::array();
            for (const auto& fit : bins) {
                results.push_back(encode_fit(fit));
            }
            charges.push_back(std::move(results));
        }
        grid.push_back(std::move(charges));
    }
    const nlohmann::json metadata = {{"schema_version", kSchemaVersion},
                                     {"input_sha256", input_sha256(cfg.input.file)},
                                     {"identity", fit_identity(cfg)},
                                     {"fit_grid", std::move(grid)}};
    TObjString object(metadata.dump().c_str());
    write_output_object(file, object, kMetadataName, TObject::kOverwrite);
}

void read_fit_results(TFile& file, Config& cfg)
{
    try {
        TKey* key = file.GetKey(kMetadataName);
        const std::unique_ptr<TObject> object(key ? key->ReadObj() : nullptr);
        const auto* text = dynamic_cast<TObjString*>(object.get());
        if (!text) {
            throw std::runtime_error("missing complete fit metadata; rerun with stages.cf3d=true");
        }
        const auto metadata = nlohmann::json::parse(text->GetString().Data());
        if (metadata.at("schema_version") != kSchemaVersion) {
            throw std::runtime_error("unsupported saved-fit schema; rerun with stages.cf3d=true");
        }
        if (metadata.at("identity") != fit_identity(cfg)) {
            throw std::runtime_error(
                "saved fits use a different binning, selection or fit configuration");
        }
        if (metadata.at("input_sha256") != input_sha256(cfg.input.file)) {
            throw std::runtime_error("saved fits were produced from a different input file");
        }
        const auto& saved = metadata.at("fit_grid");
        if (!saved.is_array() || saved.size() != charge::kCount) {
            throw std::runtime_error("invalid saved fit-grid charge count");
        }
        FitGrid grid;
        for (const auto& charge_results : saved) {
            if (!charge_results.is_array() || charge_results.size() != centrality::kCount) {
                throw std::runtime_error("invalid saved fit-grid centrality count");
            }
            std::vector<std::vector<FitResult>> centralities;
            for (const auto& bin_results : charge_results) {
                if (!bin_results.is_array() ||
                    bin_results.size() != static_cast<std::size_t>(cfg.binning.count)) {
                    throw std::runtime_error("invalid saved fit-grid bin count");
                }
                std::vector<FitResult> bins;
                for (const auto& result : bin_results) {
                    bins.push_back(decode_fit(result));
                }
                centralities.push_back(std::move(bins));
            }
            grid.push_back(std::move(centralities));
        }
        for (const int ch : cfg.selection.charges) {
            for (const int centr : cfg.selection.centralities) {
                for (int b = 0; b < cfg.binning.count; ++b) {
                    // An unavailable fit has attempts=0 and deliberately remains unavailable.
                    if (grid[ch][centr][b].attempts == 0) {
                        continue;
                    }
                    const std::string name =
                        get_cf_name(ch, centr, cfg.input.type, cfg.binning.names[b]);
                    TKey* histogram_key = file.GetKey(name.c_str());
                    const std::unique_ptr<TObject> histogram(
                        histogram_key ? histogram_key->ReadObj() : nullptr);
                    auto* cf = dynamic_cast<TH3D*>(histogram.get());
                    if (cf) {
                        cf->SetDirectory(nullptr);
                    }
                    if (!cf) {
                        throw std::runtime_error("saved fit has no corresponding 3D CF: " + name);
                    }
                }
            }
        }
        cfg.fit_results = std::move(grid);
    } catch (const std::exception& error) {
        throw std::runtime_error("cannot resume from " + std::string(file.GetName()) + ": " +
                                 error.what());
    }
}
