#include "dataProcess.hpp"

extern const int networkSize=1280000;
extern const double g=0.5;
extern const double binSize=0.1;
extern const int precision=1e4;

int main(){
    std::vector<int> ensembleList(20,20000);
    // std::vector<int> ensembleList={10};
    // ensembleList[0]=400000;

    // meanCore("meanClusterSize", ensembleList);
    // meanCore("orderParameter", ensembleList);
    // meanOrderParameterDistribution(ensembleList);
    // logBinClusterSizeDistribution(ensembleList,1);
    logBinInterEventTimeDistribution(ensembleList);

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