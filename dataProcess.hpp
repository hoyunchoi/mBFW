#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "../library-Git/CSV.hpp"
#include "../library-Git/linearAlgebra.hpp"

#include "parameters.hpp"
#include "mBFW.hpp"

namespace mBFW::data{
    //* Declare input variables which are not declared at mBFW.hpp
    std::vector<int> ensembleList;
    int fileNum, totalEnsemble;
    const bool withoutCoreNum = true;
    double logBinDelta;
    int logBinNum;
    std::vector<double> exponent;
    std::vector<double> logBinMin, logBinMax, logBinned;
    std::vector<double> LogBinMin, LogBinMax, LogBinned;

    //* Set parameters for core average
    void setParameters(const int& t_networkSize, const double& t_acceptanceThreshold, const std::vector<int>& t_ensembleList, const double& t_logBinDelta){
        networkSize = t_networkSize;
        acceptanceThreshold = t_acceptanceThreshold;
        ensembleList = t_ensembleList;
        fileNum = t_ensembleList.size();
        totalEnsemble = std::accumulate(t_ensembleList.begin(), t_ensembleList.end(), 0);
        logBinDelta = t_logBinDelta;
        exponent = arange(-8, 0, t_logBinDelta);
        logBinNum = exponent.size()-1;
        logBinMin.resize(logBinNum);
        logBinMax.resize(logBinNum);
        logBinned.resize(logBinNum);
        for (int i=0; i<logBinNum; ++i){
            logBinMin[i] = pow(10, exponent[i]);
            logBinMax[i] = pow(10, exponent[i+1]);
            logBinned[i] = sqrt(logBinMin[i] * logBinMax[i]);
        }
        LogBinMin = logBinMin * 1e8;
        LogBinMax = logBinMax * 1e8;
        LogBinned = logBinned * 1e8;
        std::tie(time_orderParameterDistribution, orderParameter_clusterSizeDistribution, m_c, t_c) = getParameters(networkSize, acceptanceThreshold);
    }

    //! Order Parameter
    void average_orderParameter(){
        std::vector<double> average(networkSize);
        std::vector<double> temp(networkSize);
        for (int core=0; core<fileNum; ++core){
            coreNum = core;
            ensembleSize = ensembleList[core];
            const std::string readFile = directory + "orderParameter/" + defaultFileName();
            readCSV(readFile, temp);
            average += temp;
        }
        average /= fileNum;
        ensembleSize = totalEnsemble;
        coreNum = 0;
        const std::string writeFile1 = directory + "orderParameter/" + defaultFileName();
        const std::string writeFile2 = directory + "orderParameter/average/" + defaultFileName(withoutCoreNum);
        writeCSV(writeFile1, average);
        writeCSV(writeFile2, average);
    }

    //! Mean Cluster Size
    void average_meanClusterSize(){
        std::vector<double> average(networkSize);
        std::vector<double> temp(networkSize);
        for (int core=0; core<fileNum; ++core){
            coreNum = core;
            ensembleSize = ensembleList[core];
            const std::string readFile = directory + "meanClusterSize/" + defaultFileName();
            readCSV(readFile, temp);
            average += temp;
        }
        average /= fileNum;
        ensembleSize = totalEnsemble;
        coreNum = 0;
        const std::string writeFile1 = directory + "meanClusterSize/" + defaultFileName();
        const std::string writeFile2 = directory + "meanClusterSize/average/" + defaultFileName(withoutCoreNum);
        writeCSV(writeFile1, average);
        writeCSV(writeFile2, average);
    }

    //! Second Giant
    void average_meanClusterSize(){
        std::vector<double> average(networkSize);
        std::vector<double> temp(networkSize);
        for (int core=0; core<fileNum; ++core){
            coreNum = core;
            ensembleSize = ensembleList[core];
            const std::string readFile = directory + "secondGiant/" + defaultFileName();
            readCSV(readFile, temp);
            average += temp;
        }
        average /= fileNum;
        ensembleSize = totalEnsemble;
        coreNum = 0;
        const std::string writeFile1 = directory + "secondGiant/" + defaultFileName();
        const std::string writeFile2 = directory + "secondGiant/average/" + defaultFileName(withoutCoreNum);
        writeCSV(writeFile1, average);
        writeCSV(writeFile2, average);
    }

    //! Inter Event Time
    void average_interEventTime(){
        //* average
        std::map<double, double> total;
        std::map<double, int> sampled;
        for (auto state : states){
            std::map<double, double> average;
            std::map<double, int> sampledAverage;
            std::vector<std::vector<double>> temp;
            for (int core=0; core<fileNum; ++core){
                coreNum = core;
                ensembleSize = ensembleList[core];
                const std::string readFile = directory + "interEventTime/" + state + "/" + defaultFileName();
                readCSV(readFile, temp);
                for (auto e : temp){
                    average[e[0]] += e[1];
                    ++sampledAverage[e[0]];
                }
            }
            for (auto it=sampledAverage.begin(); it!=sampledAverage.end(); ++it){
                average[it->first] /= it->second;
            }
            ensembleSize = totalEnsemble;
            coreNum = 0;
            const std::string writeFile1 = directory + "interEventTime/" + state + "/" + defaultFileName();
            writeCSV(writeFile1, average);
            for (auto it=average.begin(); it!=average.end(); ++it){
                total[it->first] += it->second;
                ++sampled[it->first];
            }
        }
        for (auto it=sampled.begin(); it!=sampled.end(); ++it){
            total[it->first] /= it->second;
        }
        ensembleSize = totalEnsemble;
        const std::string writeFile2 = directory + "interEventTime/average/" + defaultFileName(withoutCoreNum);
        writeCSV(writeFile2, total);

        //* log bin w.r.t t_c-t
        std::map<double, double> result;
        std::map<double, int> sampledLogBin;
        for (auto it=total.begin(); it!=total.end(); ++it){
            bool stop = true;
            const double deltaT = t_c-it->first;
            for (int i=0; i<logBinNum; ++i){
                if (logBinMin[i] <= deltaT && deltaT < logBinMax[i]){
                    result[logBinned[i]] += it->second;
                    ++sampledLogBin[logBinned[i]];
                    stop = false;
                    break;
                }
            }
            if (stop){
                break;
            }
        }
        for (auto it=sampledLogBin.begin(); it!=sampledLogBin.end(); ++it){
            result[it->first] /= it->second;
        }
        ensembleSize = totalEnsemble;
        const std::string writeFile3 = directory + "interEventTime/logBin/" + defaultFileName(withoutCoreNum);
        writeCSV(writeFile3, result);
    }

    //! Delta Acceptance
    void average_deltaAcceptance(){
        //* average
        std::map<double, double> total;
        std::map<double, int> sampled;
        for (auto state : states){
            std::map<double, double> average;
            std::map<double, int> sampledAverage;
            std::vector<std::vector<double>> temp;
            for (int core=0; core<fileNum; ++core){
                coreNum = core;
                ensembleSize = ensembleList[core];
                const std::string readFile = directory + "deltaAcceptance/" + state + "/" + defaultFileName();
                readCSV(readFile, temp);
                for (auto e : temp){
                    average[e[0]] += e[1];
                    ++sampledAverage[e[0]];
                }
            }
            for (auto it=sampledAverage.begin(); it!=sampledAverage.end(); ++it){
                average[it->first] /= it->second;
            }
            ensembleSize = totalEnsemble;
            coreNum = 0;
            const std::string writeFile1 = directory + "deltaAcceptance/" + state + "/" + defaultFileName();
            writeCSV(writeFile1, average);
            for (auto it=average.begin(); it!=average.end(); ++it){
                total[it->first] += it->second;
                ++sampled[it->first];
            }
        }
        for (auto it=sampled.begin(); it!=sampled.end(); ++it){
            total[it->first] /= it->second;
        }
        ensembleSize = totalEnsemble;
        const std::string writeFile2 = directory + "deltaAcceptance/average/" + defaultFileName(withoutCoreNum);
        writeCSV(writeFile2, total);


        //* log bin w.r.t t_c-t
        std::map<double, double> result;
        std::map<double, int> sampledLogBin;
        for (auto it=total.begin(); it!=total.end(); ++it){
            bool stop = true;
            const double deltaT = t_c-it->first;
            for (int i=0; i<logBinNum; ++i){
                if (logBinMin[i] <= deltaT && logBinMax[i] < deltaT){
                    result[logBinned[i]] += it->second;
                    ++sampledLogBin[logBinned[i]];
                    stop = false;
                    break;
                }
            }
            if (stop){
                break;
            }
        }
        for (auto it=sampledLogBin.begin(); it!=sampledLogBin.end(); ++it){
            result[it->first] /= it->second;
        }
        ensembleSize = totalEnsemble;
        const std::string writeFile3 = directory + "deltaAcceptance/logBin/" + defaultFileName(withoutCoreNum);
        writeCSV(writeFile3, result);
    }

    //! Order Parameter Distribution
    void average_orderParameterDistribution(){
        for (const double& t : time_orderParameterDistribution){
            std::map<double, double> average, temp;
            for (int core=0; core<fileNum; ++core){
                coreNum = core;
                ensembleSize = ensembleList[core];
                const std::string readFile = directory + "orderParameterDistribution/" + filename_time(t);
                readCSV(readFile, temp);
                for (auto it=temp.begin(); it!=temp.end(); ++it){
                    average[it->first] += it->second;
                }
            }
            double tot = 0.0;
            for (auto it=average.begin(); it!=average.end(); ++it){
                tot += it->second;
            }
            for (auto it=average.begin(); it!=average.end(); ++it){
                it->second /= tot;
            }
            ensembleSize = totalEnsemble;
            coreNum = 0;
            const std::string writeFile1 = directory + "orderParameterDistribution/" + filename_time(t);
            const std::string writeFile2 = directory + "orderParameterDistribution/average/" + filename_time(t,withoutCoreNum);
            writeCSV(writeFile1, average);
            writeCSV(writeFile2, average);
        }
    }

    //! Cluster Size Distribution
    void average_clusterSizeDistribution(){
        for (const double& op : orderParameter_clusterSizeDistribution){
            //* average
            std::map<double, double> average, temp;
            for (int core=0; core<fileNum; ++core){
                coreNum = core;
                ensembleSize = ensembleList[core];
                const std::string readFile = directory + "clusterSizeDistribution/" + filename_orderParameter(op);
                readCSV(readFile, temp);
                for (auto it=temp.begin(); it!=temp.end(); ++it){
                    average[it->first] += it->second;
                }
            }
            double tot = 0.0;
            for (auto it=average.begin(); it!=average.end(); ++it){
                tot += it->second;
            }
            for (auto it=average.begin(); it!=average.end(); ++it){
                it->second /= tot;
            }
            ensembleSize = totalEnsemble;
            coreNum = 0;
            const std::string writeFile1 = directory + "clusterSizeDistribution/" + filename_orderParameter(op);
            const std::string writeFile2 = directory + "clusterSizeDistribution/average/" + filename_orderParameter(op,withoutCoreNum);
            writeCSV(writeFile1, average);
            writeCSV(writeFile2, average);

            //* Log Bin
            std::map<double, double> result;
            std::map<double, int> sampled;
            for (auto it=average.begin(); it!=average.end(); ++it){
                for (int i=0; i<logBinNum; ++i){
                    if (LogBinMin[i] <= it->first && it->first < LogBinMax[i]){
                        result[LogBinned[i]] += it->second;
                        ++sampled[LogBinned[i]];
                        break;
                    }
                }
            }
            for (auto it=sampled.begin(); it!=sampled.end(); ++it){
                result[it->first] /= it->second;
            }
            ensembleSize = totalEnsemble;
            const std::string writeFile3 = directory + "clusterSizeDistribution/logBin/" + filename_orderParameter(op, withoutCoreNum);
            writeCSV(writeFile3, result);
        }
    }

    //! Age Distribution
    void average_ageDistribution(){
        for (auto state: states){
            //* average
            std::map<int, double> average, temp;
            for (int core=0; core<fileNum; ++core){
                coreNum = core;
                ensembleSize = ensembleList[core];
                const std::string readFile = directory + "ageDistribution/" + state + "/" + defaultFileName();
                readCSV(readFile, temp);
                for (auto it=temp.begin(); it!=temp.end(); ++it){
                    average[it->first] += it->second;
                }
            }
            double tot = 0.0;
            for (auto it=average.begin(); it!=average.end(); ++it){
                tot += it->second;
            }
            for (auto it=average.begin(); it!=average.end(); ++it){
                it->second /= tot;
            }
            coreNum = 0;
            ensembleSize = totalEnsemble;
            const std::string writeFile1 = directory + "ageDistribution/" + state + "/" + defaultFileName();
            writeCSV(writeFile1, average);

            //* Log Bin
            std::map<double, double> result;
            std::map<double, int> sampled;
            for (auto it=average.begin(); it!=average.end(); ++it){
                for (int i=0; i<logBinNum; ++i){
                    if (LogBinMin[i] <= it->first && it->first < LogBinMax[i]){
                        result[LogBinned[i]] += it->second;
                        ++sampled[LogBinned[i]];
                        break;
                    }
                }
            }
            for (auto it=sampled.begin(); it!=sampled.end(); ++it){
                result[it->first] /= it->second;
            }
            ensembleSize = totalEnsemble;
            const std::string writeFile3 = directory + "ageDistribution/" + state +"/logBin/" + defaultFileName(withoutCoreNum);
            writeCSV(writeFile3, result);
        }
    }

    //! Inter Event Time Distribution
    void average_interEventTimeDistribution(){
        for (auto state : states){
            //* average
            std::map<int, double> average, temp;
            for (int core=0; core<fileNum; ++core){
                coreNum = core;
                ensembleSize = ensembleList[core];
                const std::string readFile = directory + "interEventTimeDistribution/" + state + "/" + defaultFileName();
                readCSV(readFile, temp);
                for (auto it=temp.begin(); it!=temp.end(); ++it){
                    average[it->first] += it->second;
                }
            }
            double tot = 0.0;
            for (auto it=average.begin(); it!=average.end(); ++it){
                tot += it->second;
            }
            for (auto it=average.begin(); it!=average.end(); ++it){
                it->second /= tot;
            }
            coreNum = 0;
            ensembleSize = totalEnsemble;
            const std::string writeFile1 = directory + "interEventTimeDistribution/" + state + "/" + defaultFileName();
            writeCSV(writeFile1, average);

            //* Log Bin
            std::map<double, double> result;
            std::map<double, int> sampled;
            for (auto it=average.begin(); it!=average.end(); ++it){
                for (int i=0; i<logBinNum; ++i){
                    if (LogBinMin[i] <= it->first && it->first < LogBinMax[i]){
                        result[LogBinned[i]] += it->second;
                        ++sampled[LogBinned[i]];
                        break;
                    }
                }
            }
            for (auto it=sampled.begin(); it!=sampled.end(); ++it){
                result[it->first] /= it->second;
            }
            ensembleSize = totalEnsemble;
            const std::string writeFile3 = directory + "interEventTimeDistribution/" + state +"/logBin/" + defaultFileName(withoutCoreNum);
            writeCSV(writeFile3, result);
        }
    }

    //! Delta Upper Bound Distribution
    void average_deltaUpperBoundDistribution(){
        for (auto state : states){
            //* average
            std::map<int, double> average, temp;
            for (int core=0; core<fileNum; ++core){
                coreNum = core;
                ensembleSize = ensembleList[core];
                const std::string readFile = directory + "deltaUpperBoundDistribution/" + state + "/" + defaultFileName();
                readCSV(readFile, temp);
                for (auto it=temp.begin(); it!=temp.end(); ++it){
                    average[it->first] += it->second;
                }
            }
            double tot = 0.0;
            for (auto it=average.begin(); it!=average.end(); ++it){
                tot += it->second;
            }
            for (auto it=average.begin(); it!=average.end(); ++it){
                it->second /= tot;
            }
            coreNum = 0;
            ensembleSize = totalEnsemble;
            const std::string writeFile1 = directory + "deltaUpperBoundDistribution/" + state + "/" + defaultFileName();
            writeCSV(writeFile1, average);

            //* Log Bin
            std::map<double, double> result;
            std::map<double, int> sampled;
            for (auto it=average.begin(); it!=average.end(); ++it){
                for (int i=0; i<logBinNum; ++i){
                    if (LogBinMin[i] <= it->first && it->first < LogBinMax[i]){
                        result[LogBinned[i]] += it->second;
                        ++sampled[LogBinned[i]];
                        break;
                    }
                }
            }
            for (auto it=sampled.begin(); it!=sampled.end(); ++it){
                result[it->first] /= it->second;
            }
            ensembleSize = totalEnsemble;
            const std::string writeFile3 = directory + "deltaUpperBoundDistribution/" + state +"/logBin/" + defaultFileName(withoutCoreNum);
            writeCSV(writeFile3, result);
        }
    }

    //! Delta Acceptance Distribution
    void average_deltaAcceptanceDistribution(){
        for (auto state : states){
            //* average
            std::map<double, double> average, temp;
            for (int core=0; core<fileNum; ++core){
                coreNum = core;
                ensembleSize = ensembleList[core];
                const std::string readFile = directory + "deltaAcceptanceDistribution/" + state + "/" + defaultFileName();
                readCSV(readFile, temp);
                for (auto it=temp.begin(); it!=temp.end(); ++it){
                    average[it->first] += it->second;
                }
            }
            double tot = 0.0;
            for (auto it=average.begin(); it!=average.end(); ++it){
                tot += it->second;
            }
            for (auto it=average.begin(); it!=average.end(); ++it){
                it->second /= tot;
            }
            const std::string writeFile1 = directory + "deltaAcceptanceDistribution/" + state + "/" + defaultFileName();
            writeCSV(writeFile1, average);

            //* log bin
            std::map<double, double> result;
            std::map<double, int> sampled;
            for (auto it=average.begin(); it!=average.end(); ++it){
                for (int i=0; i<logBinNum; ++i){
                    if (logBinMin[i] <= it->first && it->first < logBinMax[i]){
                        result[logBinned[i]] += it->second;
                        ++sampled[logBinned[i]];
                        break;
                    }
                }
            }
            for (auto it=sampled.begin(); it!=sampled.end(); ++it){
                result[it->first] /= it->second;
            }
            const std::string writeFile3 = directory + "deltaAcceptanceDistributuion/" + state + "/logBin/" + defaultFileName(withoutCoreNum);
            writeCSV(writeFile3, result);
        }
    }


    //! Inter Event Time vs Delta Acceptance
    void average_interEventTime_DeltaAcceptance(){
        //* average
        std::map<int, double> average, temp;
        std::map<int, int> average_sampled;
        for (int core=0; core<networkSize; ++core){
            coreNum = core;
            ensembleSize = ensembleList[core];
            const std::string readFile = directory + "interEVentTime_DeltaAcceptance/" + defaultFileName();
            readCSV(readFile, temp);
            for (auto it=temp.begin(); it!=temp.end(); ++it){
                average[it->first] += it->second;
                ++average_sampled[it->first];
            }
            for (auto it=average.begin(); it!=average.end(); ++it){
                it->second /= average_sampled[it->first];
            }
        }
        const std::string writeFile1 = directory + "interEventTime_DeltaAcceptance/" + defaultFileName();
        const std::string writeFile2 = directory + "interEventTime_DeltaAcceptance/average/" + defaultFileName(withoutCoreNum);
        writeCSV(writeFile1, average);
        writeCSV(writeFile2, average);

        //* Log Bin iet and log bin delta acceptance
        std::map<double, double> result;
        std::map<double, int> sampled;
        for (auto it=average.bein(); it!=average.end(); ++it){
            for (int i=0; i<logBinNum; ++i){
                if (LogBinMin[i] <= it->first && it->first < LogBinMax[i]){
                    result[LogBinned[i]] += it->second;
                    ++sampled[LogBinned[i]];
                    break;
                }
            }
        }
        for (auto it=result.begin(); it!=result.end(); ++it){
            it->second /= sampled[it->first];
        }
        const std::string writeFile3 = directory + "interEventTime_DeltaAcceptance/logBin/" + defaultFileName(withoutCoreNum);
        writeCSV(writeFile3, result);
    }

    //! Upper Bound vs Delta Acceptance


    //! Delta Upper Bound vs Delta Acceptance





}//* End of namespace mBFW::process


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
;    int totalBinNum=0;
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