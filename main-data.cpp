#include "dataProcess.hpp"

int main(){
    const int networkSize=10000;
    const double acceptanceThreshold = 0.5;
    const double logBinDelta = 0.1;
    std::vector<int> ensembleList(1,2);

    mBFW::process::setParameters(networkSize, acceptanceThreshold, ensembleList, logBinDelta);
    // mBFW::process::f();


    // meanCore("meanClusterSize", ensembleList);
    // meanCore("orderParameter", ensembleList);
    // meanOrderParameterDistribution(ensembleList);
    // logBinClusterSizeDistribution(ensembleList,1);
    // logBinInterEventTimeDistribution(ensembleList);

    // logBinDeltaMDistribution(ensembleList);
    // logBinK_DeltaAcceptance(ensembleList);
    // logBinDeltaK_DeltaAcceptance(ensembleList);
    // logBinTime_DeltaAcceptance(ensembleList);
    // logBinAgeDistribution(ensembleList);
    // meanPeriodAcceptance_UpperBoundRatio(ensembleList);
    // meanPeriodAcceptance_DeltaK(ensembleList);
    // meanPeriodAcceptanceAreaDistribution(ensembleList);
    // meanAcceptanceDistribution(ensembleList);
    // logBinInterEventTime(ensembleList);
    // logBinInterEventTime_Acceptace(ensembleList);


    // meanCore("secondGiant", ensembleList);
    // logBinClusterSizeDistribution(ensembleList,2);
    // meanOrderParameterAfter(ensembleList);

    return 0;

}