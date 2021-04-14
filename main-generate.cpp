#include <chrono>
#include <fstream>
#include <iostream>

#include "common.hpp"
#include "generate.hpp"

int main(int argc, char* argv[]) {
    //* Get input parameters
    const int networkSize = std::stoi(argv[1]);
    const double acceptanceThreshold = std::stod(argv[2]);
    const unsigned ensembleSize = std::stoul(argv[3]);
    const int coreNum = std::stoi(argv[4]);
    const int randomEngineSeed = -1;    //* seed chosen by std::random_device()


    //* Check input network size and acceptance threshold
    if (mBFW::networkSizeList.find(networkSize) == mBFW::networkSizeList.end()) {
        std::ofstream ERROR("ERROR.log", std::ios_base::app);
        ERROR << mBFW::fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum) << ": Not valid network size " << networkSize << "\n";
        ERROR.close();
        return -1;
    }
    if (mBFW::acceptanceThresholdList.find(acceptanceThreshold) == mBFW::acceptanceThresholdList.end()) {
        std::ofstream ERROR("ERROR.log", std::ios_base::app);
        ERROR << mBFW::fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum) << ": Not valid acceptance threshold " << acceptanceThreshold << "\n";
        ERROR.close();
        return -1;
    }
    if (coreNum <= 0) {
        std::ofstream ERROR("ERROR.log", std::ios_base::app);
        ERROR << mBFW::fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum) << ": Not valid core number " << coreNum << "\n";
        ERROR.close();
        return -1;
    }

    //* Generate and Run mBFW::Generate model
    const auto start = std::chrono::system_clock::now();
    mBFW::Generate model(networkSize, acceptanceThreshold, coreNum, randomEngineSeed);
    std::cout << "temp1\n";
    model.run(ensembleSize);
    std::cout << "temp2\n";
    model.save();
    std::cout << "temp3\n";
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    std::ofstream logFile("time.log", std::ios_base::app);
    logFile << mBFW::fileName::NGE(networkSize, acceptanceThreshold, ensembleSize, coreNum) << ": " << std::setprecision(6) << sec.count() << "seconds\n";

    return 0;
}
