#include <chrono>

#include "generate.hpp"

int main(int argc, char *argv[]){
    //* Input parameters
    const int networkSize=std::stoul(argv[1]);
    const double g=std::stod(argv[2]);
    const int ensembleSize=std::stoul(argv[3]);
    const std::string machine=argv[4];
    const int coreNum=std::stoul(argv[5]);
    const int randomEngineSeed = -1;

    //* Set precision and parameters
    double precision;
    if(g==0.2){
        precision=1e3;
    }
    else{
        precision=1e4;
    }

    //* Determine which observables to calculate


    //* run mBFW
    auto start=std::chrono::system_clock::now();
    mBFW::generate::setParameters(networkSize, ensembleSize, g, precision, coreNum, randomEngineSeed);
    mBFW::generate::run();
    std::chrono::duration<double> sec=std::chrono::system_clock::now()-start;
    printf(" %.6fs for N=%.1e, g=%.1f, ensemble=%d-%d at %s\n", sec.count(),(double)networkSize, g, ensembleSize, coreNum, machine.c_str());


    //* save parameters
    start = std::chrono::system_clock::now();

    mBFW::generate::save_orderParameter();
    mBFW::generate::save_meanClusterSize();
    mBFW::generate::save_secondGiant();
    mBFW::generate::save_interEventTime();
    mBFW::generate::save_deltaAcceptance();
    mBFW::generate::save_orderParameterDistribution();
    mBFW::generate::save_clusterSizeDistribution();
    mBFW::generate::save_ageDistribution();
    mBFW::generate::save_interEventTimeDistribution();
    mBFW::generate::save_deltaUpperBoundDistribution();
    mBFW::generate::save_deltaAcceptanceDistribution();
    mBFW::generate::save_interEventTime_DeltaAcceptance();
    mBFW::generate::save_upperBound_DeltaAcceptance();
    mBFW::generate::save_deltaUpperBound_DeltaAcceptance();
    mBFW::generate::save_dynamics();

    sec=std::chrono::system_clock::now()-start;
    printf(" %0.6fs for saving\n", sec.count());

    return 0;
}
