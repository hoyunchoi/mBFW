#pragma once

#include <set>
#include <string>
#include <vector>

#include "../library/stringFormat.hpp"

namespace mBFW {
const std::string dataDirectory = "../data/mBFW/";
const std::vector<std::string> states = {"0A1", "A1A2", "A2B", "BC", "C1"};
const std::vector<std::string> pointTypes = {"m_a1", "m_a2", "m_b", "m_c", "m_inflection", "t_a1", "t_a2", "t_b", "t_c", "t_inflection"};
const std::set<int> networkSizeList = {10000, 20000, 40000, 80000, 160000, 320000, 640000, 1280000, 2560000, 5120000, 10240000};
const std::set<double> acceptanceThresholdList = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

namespace fileName {
inline const std::string base(const int&, const double&);

const std::string NG(const int&, const double&, const int&);

const std::string NGE(const int&, const double&, const unsigned&, const int&);

const std::string NGES(const int&, const double&, const unsigned&, const std::string&, const double&, const int&);
}  // namespace fileName

inline const std::string fileName::base(const int& t_networkSize, const double& t_acceptanceThreshold) {
    return "N" + to_stringWithExponent((double)t_networkSize, 1) + ",G" + to_stringWithPrecision(t_acceptanceThreshold, 1);
}

const std::string fileName::NG(const int& t_networkSize, const double& t_acceptanceThreshold, const int& t_coreNum = -1) {
    const std::string fileName = base(t_networkSize, t_acceptanceThreshold);
    return t_coreNum == -1 ? fileName + ".txt" : fileName + "-" + std::to_string(t_coreNum) + ".txt";
}

const std::string fileName::NGE(const int& t_networkSize, const double& t_acceptanceThreshold, const unsigned& t_ensembleSize, const int& t_coreNum = -1) {
    const std::string fileName = base(t_networkSize, t_acceptanceThreshold) + ",E" + std::to_string(t_ensembleSize);
    return t_coreNum == -1 ? fileName + ".txt" : fileName + "-" + std::to_string(t_coreNum) + ".txt";
}

const std::string fileName::NGES(const int& t_networkSize, const double& t_acceptanceThreshold, const unsigned& t_ensembleSize, const std::string& t_standard, const double& t_repeater, const int& t_coreNum = -1) {
    const std::string fileName = base(t_networkSize, t_acceptanceThreshold) + ",E" + std::to_string(t_ensembleSize) + "," + t_standard + to_stringWithPrecision(t_repeater, 6);
    return t_coreNum == -1 ? fileName + ".txt" : fileName + "-" + std::to_string(t_coreNum) + ".txt";
}
}  // namespace mBFW