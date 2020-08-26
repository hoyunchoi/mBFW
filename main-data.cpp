#include <chrono>

#include "data.hpp"

int main(){
    const int networkSize=10000;
    const double acceptanceThreshold = 0.5;
    const double logBinDelta = 0.1;
    std::vector<int> ensembleList(2,10);


    //* Determine which observables to calculate
    std::vector<bool> observables(15);
    observables[0] = false;      //! Order Parameter
    observables[1] = false;      //! Mean Cluster Size
    observables[2] = false;      //! Second Giant
    observables[3] = false;      //! Inter Event Time
    observables[4] = false;      //! Delta Acceptance
    observables[5] = false;      //! Order Parameter Distribution
    observables[6] = false;      //! Cluster Size Distribution
    observables[7] = false;      //! Age Distribution
    observables[8] = true;      //! Inter Event Time Distribution
    observables[9] = true;      //! Delta Upper Bound Distribution
    observables[10] = true;     //! Delta Acceptance Distribution
    observables[11] = true;     //! Inter Event Time vs Delta Acceptance
    observables[12] = true;     //! Upper Bound vs Delta Acceptance
    observables[13] = true;     //! Delta Upper Bound vs Delta Acceptance
    observables[14] = false;    //! Dynamics

    auto start = std::chrono::system_clock::now();
    mBFW::data::setParameters(networkSize, acceptanceThreshold, ensembleList, logBinDelta, observables);
    mBFW::data::average();

    // mBFW::data::average_orderParameter();
    // mBFW::data::average_meanClusterSize();
    // mBFW::data::average_secondGiant();
    // mBFW::data::average_interEventTime();
    // mBFW::data::average_deltaAcceptance();
    // mBFW::data::average_orderParameterDistribution();
    // mBFW::data::average_clusterSizeDistribution();
    // mBFW::data::average_ageDistribution();
    // mBFW::data::average_interEventTimeDistribution();
    // mBFW::data::average_deltaUpperBoundDistribution();
    // mBFW::data::average_deltaAcceptanceDistribution();
    // mBFW::data::average_interEventTime_DeltaAcceptance();
    // mBFW::data::average_upperBound_DeltaAcceptance();
    // mBFW::data::average_deltaUpperbound_DeltaAcceptance();

    std::chrono::duration<double> sec = std::chrono::system_clock::now()-start;
    printf("%.6f second to process data\n", sec.count());

    return 0;

}