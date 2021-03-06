#include <chrono>
#include <iostream>

#include "common.hpp"
#include "data.hpp"

int main(int argc, char* argv[]) {
    const int networkSize = std::stoi(argv[1]);
    const double acceptanceThreshold = std::stod(argv[2]);
    const bool deletion = true;
    const bool writeLog = true;

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
    checkList["ageDist_time"] = false;
    checkList["clusterSizeDist"] = true;
    checkList["clusterSizeDist_time"] = false;
    checkList["deltaUpperBoundDist"] = false;
    checkList["deltaUpperBoundDist_time"] = false;
    checkList["interEventTime"] = false;
    checkList["interEventTime_orderParameter"] = false;
    checkList["interEventTimeDist"] = false;
    checkList["interEventTimeDist_time"] = false;
    checkList["meanClusterSize"] = false;
    checkList["netOrderParameter"] = false;
    checkList["netSecondMoment"] = false;
    checkList["orderParameter"] = false;
    checkList["orderParameter_interEventTime"] = false;
    checkList["orderParameterDist"] = false;
    checkList["secondMaximum"] = false;
    checkList["secondMoment"] = false;

    //* Generate and Run mBFW::data model
    auto start = std::chrono::system_clock::now();
    mBFW::Data model(networkSize, acceptanceThreshold, deletion, writeLog);
    model.run(checkList);
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    std::ofstream logFile(mBFW::logDirectory + "/time.log", std::ios_base::app);
    logFile << mBFW::fileName::NG(networkSize, acceptanceThreshold) << ": " << std::setprecision(6) << sec.count() << "seconds\n";
    std::cout << mBFW::fileName::NG(networkSize, acceptanceThreshold) << ": " << std::setprecision(6) << sec.count() << "seconds\n";

    return 0;
}
