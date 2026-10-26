#include "io/input.h"

#include <memory>
#include <string>

#include <TDirectory.h>
#include <TFile.h>
#include <TH3F.h>
#include <TKey.h>
#include <TString.h>

#include "core/binning.h"
#include "core/correlation.h"
#include "core/log.h"

namespace
{
std::unique_ptr<TH3D> read_histogram(TFile& file, const TString& name, const char* label)
{
    TKey* key = file.GetKey(name);
    const std::string prefix = std::string("[") + label + "] ";
    if (!key) {
        logging::warn(prefix + "NOT FOUND: " + name.Data());
        return nullptr;
    }
    // Read a fresh object so ownership is unambiguous with either AddDirectory
    // policy and an object already cached by TFile::Get remains untouched.
    std::unique_ptr<TObject> object(key->ReadObj());
    if (!object) {
        logging::warn(prefix + "CANNOT READ: " + name.Data());
        return nullptr;
    }
    if (auto* histogram = dynamic_cast<TH1*>(object.get())) {
        histogram->SetDirectory(nullptr);
    }
    // Profiles have different statistics and must not be treated as count histograms.
    if (object->IsA() != TH3D::Class() && object->IsA() != TH3F::Class()) {
        logging::warn(prefix + "UNSUPPORTED CLASS: " + name.Data() + " (" + object->ClassName() +
                      "; expected TH3F or TH3D)");
        return nullptr;
    }

    // TH3::Copy converts bin storage through the destination's virtual methods,
    // preserving axes, flow bins, Sumw2 and cached statistics without float casts.
    // Keep the copy out of gDirectory even when automatic registration is enabled.
    TDirectory::TContext context(nullptr);
    auto histogram = std::make_unique<TH3D>();
    const auto& source = *static_cast<const TH3*>(object.get());
    source.Copy(*histogram);
    histogram->SetStatOverflows(source.GetStatOverflows());
    histogram->SetDirectory(nullptr);
    set_moment_storage(*histogram, moment_storage(source));
    return histogram;
}
} // namespace

std::pair<TH3D*, TH3D*> get_hists(TFile* f, int ch, int centr, int bin)
{
    if (!f) {
        logging::warn("Cannot read input histograms: null file");
        return {nullptr, nullptr};
    }
    auto num = read_histogram(*f, Form("bp_%d_%d_num_%d", ch, centr, bin), "Num");
    if (!num) {
        return {nullptr, nullptr};
    }
    auto wei = read_histogram(*f, Form("bp_%d_%d_num_wei_%d", ch, centr, bin), "NumWei");
    if (!wei) {
        return {nullptr, nullptr};
    }
    return {num.release(), wei.release()};
}

std::string get_cf_name(int ch, int centr, const std::string& binType, const std::string& binName)
{
    return Form("CF (charge=%s, centrality=%s, %s=%s)", charge::kNames[ch],
                centrality::kNames[centr], binType.c_str(), binName.c_str());
}
