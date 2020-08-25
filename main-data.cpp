#include <chrono>

#include "data.hpp"

int main(){
    const int networkSize=10000;
    const double acceptanceThreshold = 0.5;
    const double logBinDelta = 0.1;
    std::vector<int> ensembleList = {10,10};

    auto start = std::chrono::system_clock::now();
    mBFW::data::setParameters(networkSize, acceptanceThreshold, ensembleList, logBinDelta);

    mBFW::data::average_orderParameter();
    mBFW::data::average_meanClusterSize();
    mBFW::data::average_secondGiant();
    mBFW::data::average_interEventTime();
    mBFW::data::average_deltaAcceptance();
    mBFW::data::average_orderParameterDistribution();
    mBFW::data::average_clusterSizeDistribution();
    mBFW::data::average_ageDistribution();
    mBFW::data::average_interEventTimeDistribution();
    mBFW::data::average_deltaUpperBoundDistribution();
    mBFW::data::average_deltaAcceptanceDistribution();
    mBFW::data::average_interEventTime_DeltaAcceptance();
    mBFW::data::average_upperBound_DeltaAcceptance();
    mBFW::data::average_deltaUpperbound_DeltaAcceptance();

    std::chrono::duration<double> sec = std::chrono::system_clock::now()-start;
    printf("%.6f second to process data\n", sec.count());

    return 0;

}