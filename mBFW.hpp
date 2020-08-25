#pragma once

#include <iostream>
#include <vector>

#include "../library-Git/stringFormat.hpp"

namespace mBFW{
    //! Declaration of global variables defined at mBFW namespace
    static int networkSize;
    static double acceptanceThreshold;
    static int ensembleSize;
    static int coreNum;
    static const double m_a = 0.05;
    static double m_c;
    static double t_c;
    static std::vector<double> time_orderParameterDistribution;
    static std::vector<double> orderParameter_clusterSizeDistribution;
    static const std::vector<std::string> states = {"before", "during"};
    const std::string directory = "../data/mBFW/";

    //! File name conventions
    const std::string defaultFileName(const int& t_networkSize, const double& t_acceptanceThreshold, const int&t_ensembleSize, const int& t_coreNum=-1){
        const std::string filename = "N" + to_stringWithExponent((double)t_networkSize, 1) + ",G"+to_stringWithPrecision(t_acceptanceThreshold,1) + ",E"+std::to_string(t_ensembleSize);

        return t_coreNum==-1 ? filename + ".txt" : filename + "-"+std::to_string(t_coreNum)+".txt";
    }

    const std::string filename_time(const int& t_networkSize, const double& t_acceptanceThreshold, const int&t_ensembleSize, const double& t_time, const int& t_coreNum=-1){
        const std::string filename = "N" + to_stringWithExponent((double)t_networkSize, 1) + ",G"+to_stringWithPrecision(t_acceptanceThreshold,1) + ",E"+std::to_string(t_ensembleSize) + ",T"+to_stringWithPrecision(t_time,4);

        return t_coreNum==-1 ? filename + ".txt" : filename + "-"+std::to_string(t_coreNum)+".txt";
    }

    const std::string filename_orderParameter(const int& t_networkSize, const double& t_acceptanceThreshold, const int&t_ensembleSize, const double& t_orderParameter, const int& t_coreNum=-1){
        const std::string filename = "N" + to_stringWithExponent((double)t_networkSize, 1) + ",G"+to_stringWithPrecision(t_acceptanceThreshold,1) + ",E"+std::to_string(t_ensembleSize) + ",OP"+to_stringWithPrecision(t_orderParameter,4);

        return t_coreNum==-1 ? filename + ".txt" : filename + "-"+std::to_string(t_coreNum)+".txt";
    }
} //* End of namespace mBFW