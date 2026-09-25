#include "io/input.h"

#include <string>

#include <TFile.h>
#include <TString.h>

#include "core/binning.h"
#include "core/log.h"

std::pair<TH3D*, TH3D*> get_hists(TFile* f, int ch, int centr, int bin)
{
    TString numName = Form("bp_%d_%d_num_%d", ch, centr, bin);
    TH3D* num = dynamic_cast<TH3D*>(f->Get(numName));

    if (!num) {
        logging::warn(std::string("[Num] NOT FOUND: ") + numName.Data());
        return {nullptr, nullptr};
    }

    TString weiName = Form("bp_%d_%d_num_wei_%d", ch, centr, bin);
    TH3D* wei = dynamic_cast<TH3D*>(f->Get(weiName));

    if (!wei) {
        logging::warn(std::string("[NumWei] NOT FOUND: ") + weiName.Data());
        return {nullptr, nullptr};
    }
    return {num, wei};
}

std::string get_cf_name(int ch, int centr, const std::string& binType, const std::string& binName)
{
    return Form("CF (charge=%s, centrality=%s, %s=%s)", charge::kNames[ch],
                centrality::kNames[centr], binType.c_str(), binName.c_str());
}
