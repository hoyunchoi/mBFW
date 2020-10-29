#include <iostream>
#include <chrono>

#include "data.hpp"

int main(){
    const int networkSize=160000;
    const double acceptanceThreshold = 0.5;
    const bool deleteFile = true;

    const double logBinDelta = 0.1;
    std::vector<int> ensembleList(20,80000);
    // ensembleList[0] = 200000;

    //* Determine which observables to calculate
    std::vector<bool> observables(15);
    observables[0] = false;      //! Order Parameter
    observables[1] = false;      //! Mean Cluster Size
    observables[2] = false;      //! Second Giant
    observables[3] = false;      //! Inter Event Time
    observables[4] = false;      //! Delta Acceptance
    observables[5] = false;      //! Order Parameter Distribution
    observables[6] = true;      //! Cluster Size Distribution
    observables[7] = false;      //! Age Distribution
    observables[8] = false;      //! Inter Event Time Distribution
    observables[9] = false;      //! Delta Upper Bound Distribution
    observables[10] = false;     //! Delta Acceptance Distribution
    observables[11] = false;     //! Inter Event Time vs Delta Acceptance
    observables[12] = false;     //! Upper Bound vs Delta Acceptance
    observables[13] = false;     //! Delta Upper Bound vs Delta Acceqptance
    observables[14] = false;    //! Dynamics

    auto start = std::chrono::system_clock::now();
    mBFW::data::setParameters(networkSize, acceptanceThreshold, ensembleList, logBinDelta, observables, deleteFile);
    mBFW::data::run();
    std::chrono::duration<double> sec = std::chrono::system_clock::now()-start;
    printf("%.6f second to process data\n", sec.count());

    return 0;

}