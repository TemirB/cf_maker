#include <map>
#include <vector>
#include <string>

class InitialParameters {
private:
    std::map<std::string, std::map<int, std::vector<double>>> pData_;
    std::map<std::string, std::map<int, std::vector<double>>> nData_;
public:
    InitialParameters();
    
    double get(int charge, std::string parameter, int centrality, int y);
};