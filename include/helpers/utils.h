#pragma once

#include <string>

int FindCentralityIndex(const std::string& name);
int FindKtIndex(double low, double high);

std::string GetExeDir();
void EnsureDir(const std::string& dir);

std::string Sha256OfFile(const std::string& fname);
