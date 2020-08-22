#include <chrono>

#include "mBFW.hpp"

int main(int argc, char *argv[]){
    //* Input parameters
    const int networkSize=std::stoul(argv[1]);
    const double g=std::stod(argv[2]);
    const int ensembleSize=std::stoul(argv[3]);
    const std::string machine=argv[4];
    const int coreNum=std::stoul(argv[5]);
    const int randomEngineSeed = 0;

    //* Set precision and parameters
    double precision;
    if(g==0.2){
        precision=1e3;
    }
    else{
        precision=1e4;
    }
    mBFW::setParameters(networkSize, ensembleSize, g, precision, coreNum, randomEngineSeed);

    //* run mBFW
    auto start=std::chrono::system_clock::now();
    mBFW::run();
    std::chrono::duration<double> sec=std::chrono::system_clock::now()-start;
    printf(" %0.6fs for N=%.1e, g=%.1f, ensemble=%d-%d at %s\n", sec.count(),(double)networkSize, g, ensembleSize, coreNum, machine.c_str());


    //* save parameters
    start = std::chrono::system_clock::now();
    
    mBFW::process::save_orderParameter();
    mBFW::process::save_meanClusterSize();
    mBFW::process::save_secondGiant();
    mBFW::process::save_interEventTime();
    mBFW::process::save_deltaAcceptance();
    mBFW::process::save_orderParameterDistribution();
    mBFW::process::save_clusterSizeDistribution();
    mBFW::process::save_ageDistribution();
    mBFW::process::save_interEventTimeDistribution();
    mBFW::process::save_deltaUpperBoundDistribution();
    mBFW::process::save_deltaAcceptanceDistribution();
    mBFW::process::save_interEventTime_DeltaAcceptance();
    mBFW::process::save_upperBound_DeltaAcceptance();
    mBFW::process::save_deltaUpperBound_DeltaAcceptance();
    mBFW::process::save_dynamics();
    
    sec=std::chrono::system_clock::now()-start;
    printf(" %0.6fs for saving\n", sec.count());

    return 0;
}
