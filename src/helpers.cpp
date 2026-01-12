#include "helpers.h"

#include <sstream>
#include <iomanip>

#include <TH3.h>
#include <TDirectory.h>

#include <sys/stat.h>
#include <sys/types.h>

int getIdx(int chIdx, int centIdx, int ktIdx, int yIdx) {
    return chIdx * (centralitySize * ktSize * rapiditySize)
         + centIdx * (ktSize * rapiditySize)
         + ktIdx * rapiditySize
         + yIdx;
}

std::string getPrefix(int chIdx, int centIdx, int yIdx) {
    std::ostringstream ss;
    ss << chIdx << "_" << 0 << "_" << centIdx << "_y"
       << std::setw(2) << std::setfill('0') << yIdx
       << "_num_";
    return ss.str();
}

bool sameBinning(const TH3D* a, const TH3D* b) {
    return a->GetNbinsX() == b->GetNbinsX()
        && a->GetNbinsY() == b->GetNbinsY()
        && a->GetNbinsZ() == b->GetNbinsZ()
        && a->GetXaxis()->GetXmin() == b->GetXaxis()->GetXmin()
        && a->GetXaxis()->GetXmax() == b->GetXaxis()->GetXmax()
        && a->GetYaxis()->GetXmin() == b->GetYaxis()->GetXmin()
        && a->GetYaxis()->GetXmax() == b->GetYaxis()->GetXmax()
        && a->GetZaxis()->GetXmin() == b->GetZaxis()->GetXmin()
        && a->GetZaxis()->GetXmax() == b->GetZaxis()->GetXmax();
}

TDirectory* getOrMakeDir(TDirectory* base, const std::string& name) {
    if (!base) return nullptr;
    if (auto d = base->GetDirectory(name.c_str())) return d;
    return base->mkdir(name.c_str());
}

void EnsureDir(const std::string& dir) {
    struct stat st;
    if (stat(dir.c_str(), &st) != 0) {
        mkdir(dir.c_str(), 0755);
    }
}