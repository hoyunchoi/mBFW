#include <chrono>

#include "generate.hpp"

int main(int argc, char *argv[]){
    //* Input parameters
    const int networkSize=std::stoul(argv[1]);
    const double acceptanceThreshold=std::stod(argv[2]);
    const int ensembleSize=std::stoul(argv[3]);
    const std::string machine=argv[4];
    const int coreNum=std::stoul(argv[5]);
    constexpr int randomEngineSeed = 0;

    //* Set precision
    double precision;
    acceptanceThreshold==0.2 ? precision=1e3 : precision=1e4;

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
    observables[8] = false;      //! Inter Event Time Distribution
    observables[9] = false;      //! Delta Upper Bound Distribution
    observables[10] = false;     //! Delta Acceptance Distribution
    observables[11] = false;     //! Inter Event Time vs Delta Acceptance
    observables[12] = false;     //! Upper Bound vs Delta Acceptance
    observables[13] = false;     //! Delta Upper Bound vs Delta Acceptance
    observables[14] = true;    //! Dynamics

    //* run mBFW
    auto start=std::chrono::system_clock::now();
    mBFW::generate::setParameters(networkSize, ensembleSize, acceptanceThreshold, precision, coreNum, randomEngineSeed, observables);
    mBFW::generate::run();
    std::chrono::duration<double> sec=std::chrono::system_clock::now()-start;
    printf(" %.6fs for N=%.1e, g=%.1f, ensemble=%d-%d at %s\n", sec.count(),(double)networkSize, acceptanceThreshold, ensembleSize, coreNum, machine.c_str());

    //* save parameters
    start = std::chrono::system_clock::now();
    mBFW::generate::save();
    sec=std::chrono::system_clock::now()-start;
    printf(" %0.6fs for saving\n", sec.count());

    return 0;
}
