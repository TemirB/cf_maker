#include <map>
#include <vector>
#include <string>

class InitialParameters {
private:
    std::map<std::string, std::map<int, std::vector<double>>> data_;
public:
    InitialParameters();
    
    double get(std::string parameter, int centrality, int kt);
};