#include "helpers/utils.h"

#include <string>
#include <libgen.h>
#include <sys/stat.h>
#include <openssl/sha.h>
#include <linux/limits.h>
#include <fstream>

#include <TString.h>

#include "helpers/constants.h"

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