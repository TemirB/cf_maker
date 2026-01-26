#include "helpers.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <libgen.h>
#include <limits.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <openssl/sha.h>

#include <TH3.h>
#include <TTree.h>
#include <TFile.h>

#include "fit/types.h"

TH3D* getNum(TFile* f, int charge, int cent, int ktIdx) {
    TString name = Form("bp_%d_%d_num_%d", charge, cent, ktIdx);
    TH3D* h = (TH3D*) f->Get(name);

    if (!h) {
        std::cerr << "[getNum] NOT FOUND: " << name << std::endl;
        return nullptr;
    }
    return h;
};

TH3D* getNumWei(TFile* f, int charge, int cent, int ktIdx) {
    TString name = Form("bp_%d_%d_num_wei_%d", charge, cent, ktIdx);
    TH3D* h = (TH3D*) f->Get(name);

    if (!h) {
        std::cerr << "[getNumWei] NOT FOUND: " << name << std::endl;
        return nullptr;
    }
    return h;
};

std::string getCFName(int chIdx, int centIdx, int ktIdx) {
    return Form("CF_%d_%d_%d", chIdx, centIdx, ktIdx);
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

std::string GetExeDir() {
    char buf[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", buf, sizeof(buf)-1);
    buf[len] = '\0';
    return std::string(dirname(buf));
}

void EnsureDir(const std::string& dir) {
    struct stat st;
    if (stat(dir.c_str(), &st) != 0) {
        mkdir(dir.c_str(), 0755);
    }
}

std::string Sha256OfFile(const std::string& fname) {
    std::ifstream f(fname, std::ios::binary);
    SHA256_CTX ctx;
    SHA256_Init(&ctx);

    char buf[8192];
    while (f.read(buf, sizeof(buf)))
        SHA256_Update(&ctx, buf, f.gcount());
    if (f.gcount())
        SHA256_Update(&ctx, buf, f.gcount());

    unsigned char hash[SHA256_DIGEST_LENGTH];
    SHA256_Final(hash, &ctx);

    char out[65];
    for (int i=0;i<32;i++) sprintf(out+2*i,"%02x",hash[i]);
    out[64]=0;
    return out;
}

int FindCentralityIndex(const std::string& name) {
    for (int i=0;i<centralitySize;i++)
        if (centralityNames[i] == name)
            return i;
    return -1;
}

int FindKtIndex(double low, double high) {
    for (int i=0;i<ktSize;i++)
        if (fabs(ktValues[i]-low)<1e-6 &&
            fabs(ktValues[i+1]-high)<1e-6)
            return i;
    return -1;
}

bool IsBadFit(const FitResult& r) {
    if (!r.ok) return true;
    if (!r.IsFinite()) return true;
    if (!r.IsValid()) return true;
    // if (r.ndf <= 0) return true;
    // if (r.Chi2Ndf() > 3.0) return true;
    return false;
}

void CollectBadFits(
    FitGrid& fitRes,
    std::vector<BadFitPoint>& badPoints
) {
    for (int ch = 0; ch < chargeSize; ch++)
    for (int cent = 0; cent < centralitySize; cent++)
    for (int kt = 0; kt < ktSize; kt++){
        const FitResult& res = fitRes[ch][cent][kt];
        if (!IsBadFit(res)) continue;

        BadFitPoint p{};
        p.charge = ch;
        p.cent   = cent;
        p.kt     = kt;

        p.chi2ndf = (res.ndf > 0 ? res.Chi2Ndf() : -1);
        p.pvalue  = res.pvalue;
        p.lambda  = res.lambda;
        p.elambda = res.elambda;

        for (int i = 0; i < 3; i++) {
            p.R[i]  = res.R[i];
            p.eR[i] = res.eR[i];
        }

        p.ktVal = 0.5 * (ktValues[kt] + ktValues[kt + 1]);

        badPoints.push_back(p);
    }
}

TTree* WriteBadFitTree(TFile* f, const std::vector<BadFitPoint>& badPoints) {
    f->cd();

    TTree* t = new TTree("badFits", "Bad fit points");
    BadFitPoint p;

    t->Branch("charge", &p.charge);
    t->Branch("cent", &p.cent);
    t->Branch("kt", &p.kt);

    t->Branch("chi2ndf", &p.chi2ndf);
    t->Branch("pvalue", &p.pvalue);
    t->Branch("lambda", &p.lambda);
    t->Branch("elambda", &p.elambda);

    t->Branch("R", p.R, "R[3]/D");
    t->Branch("eR", p.eR, "eR[3]/D");
    t->Branch("ktVal", &p.ktVal);

    for (const auto& b : badPoints) {
        p = b;
        t->Fill();
    }

    t->Write();
    return t;
}
