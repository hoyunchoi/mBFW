#include <chrono>

#include "mBFW.hpp"

int main(int argc, char *argv[]){
    //* Input parameters
    const int networkSize=std::stoul(argv[1]);
    const double g=std::stod(argv[2]);
    const int ensembleSize=std::stoul(argv[3]);
    const std::string machine=argv[4];
    const int coreNum=std::stoul(argv[5]);

    double precision;
    if(g==0.2){
        precision=1e3;
    }
    else{
        precision=1e4;
    }

    //* Do mBFW
    auto start=std::chrono::system_clock::now();
    mBFW(networkSize, g, ensembleSize, coreNum, precision);
    std::chrono::duration<double> sec=std::chrono::system_clock::now()-start;
    printf(" %0.1fs for N=%.1e, g=%.1f, ensemble=%d-%d at %s\n", sec.count(),(double)networkSize, g, ensembleSize, coreNum, machine.c_str());
    return 0;
}
