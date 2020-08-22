# pragma once
#include "/pds/pds172/hoyun1009/library-Git/stringFormat.hpp"

//* File names of observables
std::string defaultFileName(){
    return "N"+to_stringWithExponent((double)networkSize, 1)+"G"+to_stringWithPrecision(acceptanceThreshold,1)+"E"+to_stringWithExponent((double)ensembleSize)+"-"+std::to_string(coreNum)+".txt";
}

std::string filename_time(const double& t_time){
    return "N"+to_stringWithExponent((double)networkSize, 1)+"G"+to_stringWithPrecision(acceptanceThreshold,1)+"E"+to_stringWithExponent((double)ensembleSize)+"T"+to_stringWithPrecision(t_time,4)+"-"+std::to_string(coreNum)+".txt";
}

std::string filename_orderParameter(const double& t_orderParameter){
    return "N"+to_stringWithExponent((double)networkSize, 1)+"G"+to_stringWithPrecision(acceptanceThreshold,1)+"E"+to_stringWithExponent((double)ensembleSize)+"OP"+to_stringWithPrecision(t_orderParameter,4)+"-"+std::to_string(coreNum)+".txt";
}