#include <chrono>

#include "generate.hpp"

int main(int argc, char *argv[]){
    //* Input parameters
    const int networkSize=std::stoul(argv[1]);
    const double g=std::stod(argv[2]);
    const int ensembleSize=std::stoul(argv[3]);
    const std::string machine=argv[4];
    const int coreNum=std::stoul(argv[5]);
    constexpr int randomEngineSeed = -1;

    //* Set precision
    double precision;
    g==0.2 ? precision=1e3 : precision=1e4;

    //* Determine which observables to calculate
    const bool orderParameter = true;
    const bool meanClusterSize = true;
    const bool secondGiant = true;
    const bool interEventTime = true;
    const bool deltaAcceptance = true;
    const bool orderParameterDistribution = true;
    const bool clusterSizeDistribution = true;
    const bool ageDistribution = true;
    const bool interEventTimeDistribution = true;
    const bool deltaUpperBoundDistribution = true;
    const bool deltaAcceptanceDistribution = true;
    const bool interEventTime_DeltaAcceptance = true;
    const bool upperBound_DeltaAcceptance = true;
    const bool deltaUpperBound_DeltaAcceptance = true;
    const bool dynamics = false;

    //* run mBFW
    auto start=std::chrono::system_clock::now();
    mBFW::generate::setParameters(networkSize, ensembleSize, g, precision, coreNum, randomEngineSeed);
    mBFW::generate::run(orderParameter, meanClusterSize, secondGiant, interEventTime, deltaAcceptance, orderParameterDistribution, clusterSizeDistribution, ageDistribution, interEventTimeDistribution, deltaUpperBoundDistribution, deltaAcceptanceDistribution, interEventTime_DeltaAcceptance, upperBound_DeltaAcceptance, deltaUpperBound_DeltaAcceptance, dynamics);
    std::chrono::duration<double> sec=std::chrono::system_clock::now()-start;
    printf(" %.6fs for N=%.1e, g=%.1f, ensemble=%d-%d at %s\n", sec.count(),(double)networkSize, g, ensembleSize, coreNum, machine.c_str());

    //* save parameters
    start = std::chrono::system_clock::now();
    if (orderParameter){mBFW::generate::save_orderParameter();}
    if (meanClusterSize){mBFW::generate::save_meanClusterSize();}
    if (secondGiant){mBFW::generate::save_secondGiant();}
    if (interEventTime){mBFW::generate::save_interEventTime();}
    if (deltaAcceptance){mBFW::generate::save_deltaAcceptance();}
    if (orderParameterDistribution){mBFW::generate::save_orderParameterDistribution();}
    if (clusterSizeDistribution){mBFW::generate::save_clusterSizeDistribution();}
    if (ageDistribution){mBFW::generate::save_ageDistribution();}
    if (interEventTimeDistribution){mBFW::generate::save_interEventTimeDistribution();}
    if (deltaUpperBoundDistribution){mBFW::generate::save_deltaUpperBoundDistribution();}
    if (deltaAcceptanceDistribution){mBFW::generate::save_deltaAcceptanceDistribution();}
    if (interEventTime_DeltaAcceptance){mBFW::generate::save_interEventTime_DeltaAcceptance();}
    if (upperBound_DeltaAcceptance){mBFW::generate::save_upperBound_DeltaAcceptance();}
    if (deltaUpperBound_DeltaAcceptance){mBFW::generate::save_deltaUpperBound_DeltaAcceptance();}
    if (dynamics){mBFW::generate::save_dynamics(randomEngineSeed);}
    sec=std::chrono::system_clock::now()-start;
    printf(" %0.6fs for saving\n", sec.count());

    return 0;
}
