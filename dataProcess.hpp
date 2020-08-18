#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "/pds/pds172/hoyun1009/library/CSV.hpp"
#include "/pds/pds172/hoyun1009/library/linearAlgebra.hpp"
#include "fileName.hpp"
#include "parameters.hpp"

extern const int networkSize;
extern const double g;
extern const double binSize;
extern const int precision;
extern const int degenerated;
const std::string folder="/pds/pds172/hoyun1009/mBFW/data/";


// Get Mean Value of every cores
std::vector<double> coreMean(const std::string &t_fileWOcoreNum, const int& t_coreSize){
    // Get default value
    std::vector<double> result;
    readCSV(t_fileWOcoreNum+"-0.txt",result);

    // Read next data and add it to temp
    for (int core=1; core<t_coreSize; ++core){
        std::vector<double> temp;
        readCSV(t_fileWOcoreNum+"-"+std::to_string(core)+".txt",temp);
        {
            using namespace linearAlgebra;
            result+=temp;
        }
    }
    {
        using namespace linearAlgebra;
        result/=t_coreSize;
    }
    return result;
}

// Add x value to vector and make it to matrix
std::vector<std::vector<double>> addXValueAndCompress(const std::vector<double> &t_input, const int &t_initialX=0){
    const int size=t_input.size();
    std::vector<std::vector<double>> result;
    result.reserve(size);
    // Only save nonzero values
    for (int i=0; i<size; ++i){
        if (t_input[i]>0){
            result.push_back(std::vector<double>{(double)i+t_initialX,t_input[i]});
        }
    }
    return result;
}

// Log Bin the Rawdata. First column of the raw data should be int x_value, increasing
std::vector<std::vector<double>> logBin(const int &t_size, const std::vector<std::vector<double>> &t_rawData, const double &t_binSize){
    // Find total number of bins
    int totalBinNum=0;
    while(pow(10,t_binSize*totalBinNum)<t_size){
        ++totalBinNum;
    }
    // Log Bin
    std::vector<std::vector<double>> temp;
    temp.resize(totalBinNum,std::vector<double>{0,0});
    for (int bin=0; bin<totalBinNum; ++bin){
        const int lowBin=(int)pow(10,t_binSize*bin);
        const int highBin=(int)pow(10,t_binSize*(bin+1));
        temp[bin][0]=sqrt((double)lowBin*highBin)/t_size;
        for (int i=0; i<t_rawData.size(); ++i){
            if (lowBin<=t_rawData[i][0] && t_rawData[i][0]<highBin){
                temp[bin][1]+=t_rawData[i][1];
            }
        }
        if (highBin-lowBin>0){
            temp[bin][1]/=(double)(highBin-lowBin);
        }
    }
    // Save binned data
    std::vector<std::vector<double>> result;
    for (auto &e : temp){
        if (e[1]>0){
            result.push_back(e);
        }
    }
    return result;
}

template<typename T, typename TT>
std::vector<std::vector<double>> logBin(const int &t_size, const std::map<T,TT> &t_rawData, const double &t_binSize){
    // Find total number of bins
    int totalBinNum=0;
    while(pow(10,t_binSize*totalBinNum)<t_size){
        ++totalBinNum;
    }
    // Log Bin
    std::vector<std::vector<double>> temp;
    temp.resize(totalBinNum,std::vector<double>{0,0});
    for (int bin=0; bin<totalBinNum; ++bin){
        const int lowBin=(int)pow(10,t_binSize*bin);
        const int highBin=(int)pow(10,t_binSize*(bin+1));
        // temp[bin][0]=sqrt((double)lowBin*highBin)/t_size;
        temp[bin][0]=sqrt((double)lowBin*highBin);

        for (auto it=t_rawData.begin(); it!=t_rawData.end(); ++it){
            if (lowBin<=it->first && it->first<highBin){
                temp[bin][1]+=it->second;
            }
        }
        if (highBin-lowBin>0){
            temp[bin][1]/=(double)(highBin-lowBin);
        }
    }
    // Save binned data
    std::vector<std::vector<double>> result;
    for (auto &e : temp){
        if (e[1]>0){
            result.push_back(e);
        }
    }
    return result;
}

void meanCore(const std::string &t_type, const std::vector<int>& t_ensembleList){
    int totalEnsemble=0;
    for (auto ensemble : t_ensembleList){
        totalEnsemble+=ensemble;
    }
    std::string readFile=folder+t_type+"/"+fileName(networkSize, g, t_ensembleList[0])+"-0.txt";
    std::vector<double> mean, temp;
    readCSV(readFile, mean);
    for (int core=1; core<t_ensembleList.size(); ++core){
        readFile=folder+t_type+"/"+fileName(networkSize, g, t_ensembleList[core])+"-"+std::to_string(core)+".txt";
        readCSV(readFile, temp);
        {
            using namespace linearAlgebra;
            mean+=temp;
        }
    }
    {
        using namespace linearAlgebra;
        mean/=t_ensembleList.size();
    }
    const std::string writeFile=folder+t_type+"/mean/"+fileName(networkSize, g, totalEnsemble)+".txt";
    const std::string writeFile2=folder+t_type+"/"+fileName(networkSize, g, totalEnsemble)+"-0.txt";
    writeCSV(writeFile, mean);
    writeCSV(writeFile2, mean);
}

void meanAcceptanceDistribution(const std::vector<int>& t_ensembleList){
    int totalEnsemble=0;
    for (auto ensemble : t_ensembleList){
        totalEnsemble+=ensemble;
    }
    std::string readFile=folder+"acceptanceDistribution/"+fileName(networkSize, g, t_ensembleList[0])+"-0.txt";
    std::map<double, double> acceptance;
    std::vector<std::vector<double>>temp;
    const int coreNum=t_ensembleList.size();
    for (int core=0; core<coreNum; ++core){
        readFile=folder+"acceptanceDistribution/"+fileName(networkSize, g, t_ensembleList[core])+"-"+std::to_string(core)+".txt";
        readCSV(readFile, temp);
        for (auto &e : temp){
            acceptance[e[0]]+=e[1];
        }
    }
    for (auto it=acceptance.begin(); it!=acceptance.end(); ++it){
        it->second/=coreNum;
    }

    const std::string writeFile=folder+"acceptanceDistribution/mean/"+fileName(networkSize, g, totalEnsemble)+".txt";
    const std::string writeFile2=folder+"acceptanceDistribution/"+fileName(networkSize, g, totalEnsemble)+"-0.txt";
    const int precision=10;
    writeCSV(writeFile, acceptance, precision);
    writeCSV(writeFile2, acceptance, precision);
}

void meanOrderParameterAfter(const std::vector<int>& t_ensembleList){
    int totalEnsemble=0;
    for (auto ensemble : t_ensembleList){
        totalEnsemble+=ensemble;
    }
    std::string readFile=folder+"orderParameterAfter/"+fileName(networkSize, g, t_ensembleList[0])+"-0.txt";
    std::string readSampledFile=folder+"orderParameterAfter/sampled/"+fileName(networkSize, g, t_ensembleList[0])+"-0.txt";
    std::vector<double> mean, temp;
    std::vector<int> sampledMean, sampledTemp;
    readCSV(readFile, mean);
    readCSV(readSampledFile, sampledMean);
    for (int core=1; core<t_ensembleList.size(); ++core){
        readFile=folder+"orderParameterAfter/"+fileName(networkSize, g, t_ensembleList[core])+"-"+std::to_string(core)+".txt";
        readSampledFile=folder+"orderParameterAfter/sampled/"+fileName(networkSize, g, t_ensembleList[core])+"-"+std::to_string(core)+".txt";
        readCSV(readFile, temp);
        readCSV(readSampledFile, sampledTemp);
        {
            using namespace linearAlgebra;
            mean+=temp;
            sampledMean+=sampledTemp;
        }
    }
    temp=mean;
    for (int i=0; i<temp.size(); ++i){
        if (sampledMean[i]!=0){
            temp[i]/=sampledMean[i];
        }
    }

    const std::string writeFile=folder+"orderParameterAfter/mean/"+fileName(networkSize, g, totalEnsemble)+".txt";
    const std::string writeFile2=folder+"orderParameterAfter/"+fileName(networkSize, g, totalEnsemble)+"-0.txt";
    const std::string writeSampledFile=folder+"orderParameterAfter/sampled/"+fileName(networkSize, g, totalEnsemble)+"-0.txt";
    writeCSV(writeFile, temp);
    writeCSV(writeFile2, mean);
    writeCSV(writeSampledFile, sampledMean);

}

void meanOrderParameterDistribution(const std::vector<int> t_ensembleList){
    int totalEnsemble=0;
    for (auto ensemble : t_ensembleList){
        totalEnsemble+=ensemble;
    }
    std::vector<double> orderParameterDistributionTime,remove;
    double m_c;
    std::tie(orderParameterDistributionTime, remove, m_c)=getParameters(networkSize, g);

    for (auto &time : orderParameterDistributionTime){
        std::string readFile=folder+"orderParameterDistribution/"+fileName(networkSize, g, t_ensembleList[0])+",t="+to_stringWithPrecision(time, 4)+"-0.txt";
        std::map<double, double> meanOrderParameterDistribution;
        std::map<double, double> rawOrderParameterDistribution;
        std::vector<std::vector<double>> temp;
        readCSV(readFile, temp);
        for (auto &distribution : temp){
            double roundedSize=round(distribution[0]*precision)/precision;
            meanOrderParameterDistribution[roundedSize]+=distribution[1];
            rawOrderParameterDistribution[distribution[0]]+=distribution[1];
        }

        for (int core=1; core<t_ensembleList.size(); ++core){
            readFile=folder+"orderParameterDistribution/"+fileName(networkSize, g, t_ensembleList[core])+",t="+to_stringWithPrecision(time, 4)+"-"+std::to_string(core)+".txt";
            readCSV(readFile, temp);
            for (auto &distribution : temp){
                double roundedSize=round(distribution[0]*precision)/precision;
                meanOrderParameterDistribution[roundedSize]+=distribution[1];
                rawOrderParameterDistribution[distribution[0]]+=distribution[1];
            }
        }

        for (auto it=meanOrderParameterDistribution.begin(); it!=meanOrderParameterDistribution.end(); ++it){
            it->second/=totalEnsemble;
        }
        const std::string writeFile=folder+"orderParameterDistribution/mean/"+fileName(networkSize, g, totalEnsemble)+",t="+to_stringWithPrecision(time,4)+".txt";
        const std::string writeFile2=folder+"orderParameterDistribution/"+fileName(networkSize, g, totalEnsemble)+",t="+to_stringWithPrecision(time,4)+"-0.txt";
        writeCSV(writeFile, meanOrderParameterDistribution);
        writeCSV(writeFile2, rawOrderParameterDistribution);
    }
}

void meanPeriodAcceptance_UpperBoundRatio(const std::vector<int> &t_ensembleList){
    int totalEnsemble=0;
    for (auto ensemble : t_ensembleList){
        totalEnsemble+=ensemble;
    }
    std::string readFile=folder+"periodAcceptance_UpperBoundRatio/"+fileName(networkSize, g, t_ensembleList[0])+"-0.txt";
    std::vector<std::vector<double>> mean, temp, result;
    readCSV(readFile, mean);
    for (int core=1; core<t_ensembleList.size(); ++core){
        readFile=folder+"periodAcceptance_UpperBoundRatio/"+fileName(networkSize, g, t_ensembleList[core])+"-"+std::to_string(core)+".txt";
        readCSV(readFile, temp);
        for (int i=0; i<mean.size(); ++i){
            mean[i][1] += temp[i][1];
        }
    }
    for (int i=0; i<mean.size(); ++i){
        mean[i][1]/=t_ensembleList.size();
        if (mean[i][1]!=0){
            result.push_back(mean[i]);
        }
    }
    const std::string writeFile=folder+"periodAcceptance_UpperBoundRatio/mean/"+fileName(networkSize, g, totalEnsemble)+".txt";
    const std::string writeFile2=folder+"periodAcceptance_UpperBoundRatio/"+fileName(networkSize, g, totalEnsemble)+"-0.txt";
    const int precision=10;
    writeCSV(writeFile, result, precision);
    writeCSV(writeFile2, mean, precision);
}

void meanPeriodAcceptance_DeltaK(const std::vector<int> &t_ensembleList){
    int totalEnsemble=0;
    for (auto ensemble : t_ensembleList){
        totalEnsemble+=ensemble;
    }
    std::string readFile=folder+"periodAcceptance_DeltaK/"+fileName(networkSize, g, t_ensembleList[0])+"-0.txt";
    std::vector<std::vector<double>> mean, temp, result;
    readCSV(readFile, mean);
    for (int core=1; core<t_ensembleList.size(); ++core){
        readFile=folder+"periodAcceptance_DeltaK/"+fileName(networkSize, g, t_ensembleList[core])+"-"+std::to_string(core)+".txt";
        readCSV(readFile, temp);
        for (int i=0; i<mean.size(); ++i){
            mean[i][1] += temp[i][1];
        }
    }
    for (int i=0; i<mean.size(); ++i){
        mean[i][1]/=t_ensembleList.size();
        if (mean[i][1]!=0){
            result.push_back(mean[i]);
        }
    }
    const std::string writeFile=folder+"periodAcceptance_DeltaK/mean/"+fileName(networkSize, g, totalEnsemble)+".txt";
    const std::string writeFile2=folder+"periodAcceptance_DeltaK/"+fileName(networkSize, g, totalEnsemble)+"-0.txt";
    const int precision=10;
    writeCSV(writeFile, result, precision);
    writeCSV(writeFile2, mean, precision);
}


void meanPeriodAcceptanceAreaDistribution(const std::vector<int> &t_ensembleList){
    int totalEnsemble=0;
    for (auto ensemble : t_ensembleList){
        totalEnsemble+=ensemble;
    }
    std::map<double, double> mean;
    std::vector<std::vector<double>> temp;
    // readCSV(readFile, mean);
    for (int core=0; core<t_ensembleList.size(); ++core){
        const std::string readFile=folder+"periodAcceptanceAreaDistribution/"+fileName(networkSize, g, t_ensembleList[0])+"-"+std::to_string(core)+".txt";
        readCSV(readFile, temp);
        for (int i=0; i<temp.size(); ++i){
            if (temp[i][1]!=0){
                mean[temp[i][0]]+=temp[i][1];
            }
        }
    }
    for (auto it=mean.begin(); it!=mean.end(); ++it){
        it->second/=t_ensembleList.size();
    }

    const std::string writeFile=folder+"periodAcceptanceAreaDistribution/mean/"+fileName(networkSize, g, totalEnsemble)+".txt";
    const std::string writeFile2=folder+"periodAcceptanceAreaDistribution/"+fileName(networkSize, g, totalEnsemble)+"-0.txt";
    const int precision=10;
    writeCSV(writeFile, mean, precision);
    writeCSV(writeFile2, mean, precision);
}

void logBinClusterSizeDistribution(const std::vector<int> t_ensembleList, const int &t_exclude){
    int totalEnsemble=0;
    for (auto ensemble : t_ensembleList){
        totalEnsemble+=ensemble;
    }
    int coreNum=t_ensembleList.size();
    std::vector<double> remove, observedOrderParameter;
    double m_c;
    std::tie(remove, observedOrderParameter, m_c)=getParameters(networkSize, g);


    std::map<double,int> sampled;
    std::string readFile=folder+"clusterSizeDistribution"+std::to_string(t_exclude)+"/sampled/"+fileName(networkSize, g, t_ensembleList[0])+"-0.txt";
    readCSV(readFile, sampled);
    for (int core=1; core<coreNum; ++core){
        std::map<double,int> temp;
        readFile=folder+"clusterSizeDistribution"+std::to_string(t_exclude)+"/sampled/"+fileName(networkSize, g, t_ensembleList[core])+"-"+std::to_string(core)+".txt";
        readCSV(readFile, temp);
        for (auto it=temp.begin(); it!=temp.end(); it++){
            sampled[it->first]+=it->second;
        }
    }
    std::string writeFile=folder+"clusterSizeDistribution"+std::to_string(t_exclude)+"/sampled/"+fileName(networkSize, g, totalEnsemble)+"-0.txt";
    writeCSV(writeFile,sampled);

    for (auto &currentOrderParameter : observedOrderParameter){
        std::string readFile=folder+"clusterSizeDistribution"+std::to_string(t_exclude)+"/"+fileName(networkSize, g, t_ensembleList[0])+",m="+to_stringWithPrecision(currentOrderParameter,4)+"-0.txt";
        std::map<int, double> clusterSizeDistribution;
        std::vector<std::vector<double>> temp;
        for (int core=0; core<coreNum; ++core){
            readFile=folder+"clusterSizeDistribution"+std::to_string(t_exclude)+"/"+fileName(networkSize, g, t_ensembleList[core])+",m="+to_stringWithPrecision(currentOrderParameter,4)+"-"+std::to_string(core)+".txt";
            readCSV(readFile, temp);
            for (auto &distribution : temp){
                clusterSizeDistribution[distribution[0]]+=distribution[1];
            }
        }
        for (auto it=clusterSizeDistribution.begin(); it!=clusterSizeDistribution.end(); ++it){
            it->second/=coreNum;
        }

        const auto binnedData=logBin(networkSize, clusterSizeDistribution, binSize);
        const std::string writeFile=folder+"clusterSizeDistribution"+std::to_string(t_exclude)+"/logBin/"+fileName(networkSize, g, totalEnsemble)+",m="+to_stringWithPrecision(currentOrderParameter,4)+".txt";
        const std::string writeFile2=folder+"clusterSizeDistribution"+std::to_string(t_exclude)+"/"+fileName(networkSize, g, totalEnsemble)+",m="+to_stringWithPrecision(currentOrderParameter,4)+"-0.txt";
        writeCSV(writeFile, binnedData);
        writeCSV(writeFile2, clusterSizeDistribution);
    }
}

void logBinInterEventTimeDistribution(const std::vector<int> t_ensembleList){
    int totalEnsemble=0;
    for (auto ensemble : t_ensembleList){
        totalEnsemble+=ensemble;
    }
    const int coreNum=t_ensembleList.size();
    const std::vector<std::string> type={"before", "after"};
    for (auto typeName : type){
        std::map<int, double> interEventTime;
        std::vector<std::vector<double>> temp;
        for (int core=0; core<t_ensembleList.size(); ++core){
            std::string readFile=folder+"interEventTimeDistribution/"+typeName+"/"+fileName(networkSize, g, t_ensembleList[core])+"-"+std::to_string(core)+".txt";
            readCSV(readFile, temp);
            for (auto & iet : temp){
                interEventTime[iet[0]]+=iet[1];
            }
        }

        for (auto it=interEventTime.begin(); it!=interEventTime.end(); ++it){
            it->second/=coreNum;
        }
        const auto binnedData=logBin(networkSize, interEventTime, binSize);
        const std::string writeFile=folder+"interEventTimeDistribution/"+typeName+"/logBin/"+fileName(networkSize, g, totalEnsemble)+".txt";
        const std::string writeFile2=folder+"interEventTimeDistribution/"+typeName+"/"+fileName(networkSize, g, totalEnsemble)+"-0.txt";
        writeCSV(writeFile, binnedData);
        writeCSV(writeFile2, interEventTime);
    }
}


void logBinInterEventTime_Acceptace(const std::vector<int> &t_ensembleList){
    int totalEnsemble=0;
    for (auto ensemble : t_ensembleList){
        totalEnsemble+=ensemble;
    }
    const int coreNum=t_ensembleList.size();
    std::map<int, double> interEventTime_Acceptance;
    std::vector<std::vector<double>> temp;
    for(int core=0; core<coreNum; ++core){
        std::string readFile=folder+"interEventTime_Acceptance/"+fileName(networkSize, g, t_ensembleList[core])+"-"+std::to_string(core)+".txt";
        readCSV(readFile, temp);
        for (auto &pa : temp){
            interEventTime_Acceptance[pa[0]]+=pa[1];
        }
    }
    for (auto it=interEventTime_Acceptance.begin(); it!=interEventTime_Acceptance.end(); ++it){
        it->second/=coreNum;
    }
    const auto binnedData=logBin(networkSize, interEventTime_Acceptance, binSize);
    const std::string writeFile=folder+"interEventTime_Acceptance/logBin/"+fileName(networkSize, g, totalEnsemble)+".txt";
    const std::string writeFile2=folder+"interEventTime_Acceptance/"+fileName(networkSize, g, totalEnsemble)+"-0.txt";
    writeCSV(writeFile, binnedData);
    writeCSV(writeFile2, interEventTime_Acceptance);
}


void logBinDeltaMDistribution(const std::vector<int> t_ensembleList){
    int totalEnsemble=0;
    for (auto ensemble : t_ensembleList){
        totalEnsemble+=ensemble;
    }
    const int coreNum=t_ensembleList.size();
    const std::vector<std::string> type={"before", "after"};
    for (auto typeName : type){
        std::map<int, double> deltaM;
        std::vector<std::vector<double>> temp;
        for (int core=0; core<t_ensembleList.size(); ++core){
            std::string readFile=folder+"deltaMDistribution/"+typeName+"/"+fileName(networkSize, g, t_ensembleList[core])+"-"+std::to_string(core)+".txt";
            readCSV(readFile, temp);
            for (auto & delta : temp){
                deltaM[delta[0]]+=delta[1];
            }
        }

        for (auto it=deltaM.begin(); it!=deltaM.end(); ++it){
            it->second/=coreNum;
        }
        const auto binnedData=logBin(networkSize, deltaM, binSize);
        const std::string writeFile=folder+"deltaMDistribution/"+typeName+"/logBin/"+fileName(networkSize, g, totalEnsemble)+".txt";
        const std::string writeFile2=folder+"deltaMDistribution/"+typeName+"/"+fileName(networkSize, g, totalEnsemble)+"-0.txt";
        writeCSV(writeFile, binnedData);
        writeCSV(writeFile2, deltaM);
    }
}

void logBinK_DeltaAcceptance(const std::vector<int> t_ensembleList){
    int totalEnsemble=0;
    for (auto ensemble : t_ensembleList){
        totalEnsemble+=ensemble;
    }
    const int coreNum=t_ensembleList.size();
    std::map<int, double> deltaAcceptance;
    std::vector<std::vector<double>> temp;
    for (int core=0; core<t_ensembleList.size(); ++core){
        std::string readFile=folder+"k_DeltaAcceptance/"+fileName(networkSize, g, t_ensembleList[core])+"-"+std::to_string(core)+".txt";
        readCSV(readFile, temp);
        for (auto &acceptance : temp){
            deltaAcceptance[acceptance[0]]+=acceptance[1];
        }
    }

    for (auto it=deltaAcceptance.begin(); it!=deltaAcceptance.end(); ++it){
        it->second/=coreNum;
    }
    const auto binnedData=logBin(networkSize, deltaAcceptance, binSize);
    const std::string writeFile=folder+"k_DeltaAcceptance/logBin/"+fileName(networkSize, g, totalEnsemble)+".txt";
    const std::string writeFile2=folder+"k_DeltaAcceptance/"+fileName(networkSize, g, totalEnsemble)+"-0.txt";
    writeCSV(writeFile, binnedData);
    writeCSV(writeFile2, deltaAcceptance);
}

void logBinDeltaK_DeltaAcceptance(const std::vector<int> t_ensembleList){
    int totalEnsemble=0;
    for (auto ensemble : t_ensembleList){
        totalEnsemble+=ensemble;
    }
    const int coreNum=t_ensembleList.size();
    std::map<int, double> deltaAcceptance;
    std::vector<std::vector<double>> temp;
    for (int core=0; core<t_ensembleList.size(); ++core){
        std::string readFile=folder+"deltaK_DeltaAcceptance/"+fileName(networkSize, g, t_ensembleList[core])+"-"+std::to_string(core)+".txt";
        readCSV(readFile, temp);
        for (auto &acceptance : temp){
            deltaAcceptance[acceptance[0]]+=acceptance[1];
        }
    }

    for (auto it=deltaAcceptance.begin(); it!=deltaAcceptance.end(); ++it){
        it->second/=coreNum;
    }
    const auto binnedData=logBin(networkSize, deltaAcceptance, binSize);
    const std::string writeFile=folder+"deltaK_DeltaAcceptance/logBin/"+fileName(networkSize, g, totalEnsemble)+".txt";
    const std::string writeFile2=folder+"deltaK_DeltaAcceptance/"+fileName(networkSize, g, totalEnsemble)+"-0.txt";
    writeCSV(writeFile, binnedData);
    writeCSV(writeFile2, deltaAcceptance);
}

void logBinTime_DeltaAcceptance(const std::vector<int> t_ensembleList){
    int totalEnsemble=0;
    for (auto ensemble : t_ensembleList){
        totalEnsemble+=ensemble;
    }
    const int coreNum=t_ensembleList.size();
    std::map<int, double> deltaAcceptance;
    std::vector<std::vector<double>> temp;
    for (int core=0; core<coreNum ; ++core){
        std::string readFile=folder+"time_DeltaAcceptance/"+fileName(networkSize, g, t_ensembleList[core])+"-"+std::to_string(core)+".txt";
        readCSV(readFile, temp);
        for (auto &acceptance : temp){
            deltaAcceptance[acceptance[0]]+=acceptance[1];
        }
    }

    for (auto it=deltaAcceptance.begin(); it!=deltaAcceptance.end(); ++it){
        it->second/=coreNum;
    }

    const double t_c=0.937;
    const int binNum=80;
    std::vector<double> exponent;
    {
        using namespace linearAlgebra;
        exponent=linspace(-8,0,binNum+1);
    }
    std::vector<double> minX(binNum);
    std::vector<double> maxX(binNum);
    std::vector<int> sampledBin(binNum);
    std::vector<std::vector<double>> tempBin(binNum);
    std::vector<std::vector<double>> binnedData;
    for(int i=0; i<binNum; ++i){
        minX[i]=pow(10,exponent[i]);
        maxX[i]=pow(10,exponent[i+1]);
        tempBin[i]=std::vector<double>{sqrt(minX[i]*maxX[i]),0};
    }
    for (auto it=deltaAcceptance.begin(); it!=deltaAcceptance.end(); ++it){
        const double x=t_c-(double)it->first/networkSize;
        for (int i=0; i<binNum; ++i){
            if (minX[i]<=x && x<maxX[i]){
                tempBin[i][1]+=it->second;
                ++sampledBin[i];
            }
        }
    }
    for (int i=0; i<binNum; ++i){
        if (sampledBin[i]!=0){
            tempBin[i][1]/=sampledBin[i];
            binnedData.push_back(tempBin[i]);
        }
    }
    const std::string writeFile=folder+"time_DeltaAcceptance/logBin/"+fileName(networkSize, g, totalEnsemble)+".txt";
    const std::string writeFile2=folder+"time_DeltaAcceptance/"+fileName(networkSize, g, totalEnsemble)+"-0.txt";
    writeCSV(writeFile, binnedData);
    writeCSV(writeFile2, deltaAcceptance);
}




void logBinAgeDistribution(const std::vector<int> t_ensembleList){
    int totalEnsemble=0;
    for (auto ensemble : t_ensembleList){
        totalEnsemble+=ensemble;
    }
    const int coreNum=t_ensembleList.size();
    const std::vector<std::string> type={"before", "after"};
    for (auto typeName : type){
        std::map<int, double> age;
        std::vector<std::vector<double>> temp;
        for (int core=0; core<t_ensembleList.size(); ++core){
            std::string readFile=folder+"ageDistribution/"+typeName+"/"+fileName(networkSize, g, t_ensembleList[core])+"-"+std::to_string(core)+".txt";
            readCSV(readFile, temp);
            for (auto &a : temp){
                age[a[0]]+=a[1];
            }
        }
        double sum=0;
        for (auto it=age.begin(); it!=age.end(); ++it){
            it->second/=coreNum;
            sum+=it->second;
        }
        for (auto it=age.begin(); it!=age.end(); ++it){
            it->second/=sum;
        }
        const auto binnedData=logBin(networkSize, age, binSize);
        const std::string writeFile=folder+"ageDistribution/"+typeName+"/logBin/"+fileName(networkSize, g, totalEnsemble)+".txt";
        const std::string writeFile2=folder+"ageDistribution/"+typeName+"/"+fileName(networkSize, g, totalEnsemble)+"-0.txt";
        writeCSV(writeFile, binnedData);
        writeCSV(writeFile2, age);
    }
}

void logBinInterEventTime(const std::vector<int> t_ensembleList){
    int totalEnsemble=0;
    for (auto ensemble : t_ensembleList){
        totalEnsemble+=ensemble;
    }
    const int coreNum=t_ensembleList.size();
    std::map<int, double> interEventTime;
    std::vector<std::vector<double>> temp;
    for (int core=0; core<coreNum; ++core){
        std::string readFile=folder+"interEventTime/"+fileName(networkSize, g, t_ensembleList[core])+"-"+std::to_string(core)+".txt";
        readCSV(readFile, temp);
        for (auto &iet : temp){
            interEventTime[iet[0]]+=iet[1];
        }
    }

    for (auto it=interEventTime.begin(); it!=interEventTime.end(); ++it){
        it->second/=coreNum;
    }

    const double t_c=0.937;
    std::vector<std::vector<double>> rawData;
    for (auto it=interEventTime.begin(); it!=interEventTime.end(); ++it){
        if (t_c>(double)it->first/networkSize){
            rawData.push_back(std::vector<double>{t_c-(double)it->first/networkSize,it->second});
        }
    }
    const int binNum=80;
    std::vector<double> exponent;
    {
        using namespace linearAlgebra;
        exponent=linspace(-8,0,binNum+1);
    }
    std::vector<double> minX(binNum);
    std::vector<double> maxX(binNum);
    std::vector<int> sampledBin(binNum);
    std::vector<std::vector<double>> tempBin(binNum);
    std::vector<std::vector<double>> binnedData;
    for(int i=0; i<binNum; ++i){
        minX[i]=pow(10,exponent[i]);
        maxX[i]=pow(10,exponent[i+1]);
        tempBin[i]=std::vector<double>{sqrt(minX[i]*maxX[i]),0};
    }
    for (auto &raw : rawData){
        for (int i=0; i<binNum; ++i){
            if (minX[i]<=raw[0] && raw[0]<maxX[i]){
                tempBin[i][1]+=raw[1];
                ++sampledBin[i];
            }
        }
    }
    for (int i=0; i<binNum; ++i){
        if (sampledBin[i]!=0){
            tempBin[i][1]/=sampledBin[i];
            binnedData.push_back(tempBin[i]);
        }
    }

    const std::string writeFile=folder+"interEventTime/logBin/"+fileName(networkSize, g, totalEnsemble)+".txt";
    const std::string writeFile2=folder+"interEventTime/"+fileName(networkSize, g, totalEnsemble)+"-0.txt";
    writeCSV(writeFile, binnedData);
    writeCSV(writeFile2, interEventTime);
}