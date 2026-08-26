#pragma once

#include <map>
#include <string>
#include <vector>

enum class Binning : char
{
    Kt,
    Rapidity
};

class InitialParameters
{
  public:
    explicit InitialParameters(Binning binning, bool useDefaults = false);

    [[nodiscard]] double Get(int charge, const std::string& parameter, int centrality, int y) const;

  private:
    void fill_rapidity();
    void fill_kt();
    void fill_default_kt();
    void fill_default_rapidity();

    std::map<std::string, std::map<int, std::vector<double>>> pData_;
    std::map<std::string, std::map<int, std::vector<double>>> nData_;
};
