#include <chrono>
#include <iostream>

#include "common.hpp"
#include "data.hpp"

int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::cout.tie(NULL);

    const int networkSize = std::stoi(argv[1]);
    const double acceptanceThreshold = std::stod(argv[2]);
    const bool deletion = true;

    //* Check input network size and acceptance threshold
    if (mBFW::networkSizeList.find(networkSize) == mBFW::networkSizeList.end()) {
        std::ofstream ERROR(mBFW::logDirectory + "ERROR.log", std::ios_base::app);
        ERROR << mBFW::fileName::NG(networkSize, acceptanceThreshold) << ": Not valid network size " << networkSize << "\n";
        ERROR.close();
        return -1;
    }
    if (mBFW::acceptanceThresholdList.find(acceptanceThreshold) == mBFW::acceptanceThresholdList.end()) {
        std::ofstream ERROR(mBFW::logDirectory + "ERROR.log", std::ios_base::app);
        ERROR << mBFW::fileName::NG(networkSize, acceptanceThreshold) << ": Not valid acceptance threshold " << acceptanceThreshold << "\n";
        ERROR.close();
        return -1;
    }

    //* Check list of each observables
    std::map<std::string, bool> checkList;
    checkList["ageDist"] = false;
    checkList["clusterSizeDist"] = false;
    checkList["deltaUpperBoundDist"] = false;
    checkList["interEventTime"] = false;
    checkList["interEventTimeDist"] = false;
    checkList["interEventTime_orderParameter"] = false;
    checkList["meanClusterSize"] = false;
    checkList["orderParameter"] = false;
    checkList["orderParameterDist"] = false;
    checkList["orderParameterVariance"] = false;

    //* Generate and Run mBFW::data model
    auto start = std::chrono::system_clock::now();
    mBFW::Data model(networkSize, acceptanceThreshold);
    // model.run(checkList, deletion);
    model.temp("clusterSizeDist");
    model.temp("orderParmeterDist");
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    std::ofstream logFile(mBFW::logDirectory + "/time.log", std::ios_base::app);
    std::cout << mBFW::fileName::NG(networkSize, acceptanceThreshold) << ": " << std::setprecision(6) << sec.count() << "seconds\n";
    logFile << mBFW::fileName::NG(networkSize, acceptanceThreshold) << ": " << std::setprecision(6) << sec.count() << "seconds\n";

    return 0;
}
