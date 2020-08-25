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
    using namespace linearAlgebra;

    //! Declaration of variables only use at mBFW::data namespace
    std::vector<int> ensembleList;
    int fileNum;
    int totalEnsemble;
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
        exponent = linearAlgebra::arange(-8, 0, t_logBinDelta);
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
            const std::string readFile = directory + "orderParameter/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
            readCSV(readFile, temp);
            average += temp;
        }
        average /= fileNum;
        ensembleSize = totalEnsemble;
        coreNum = 0;
        const std::string writeFile1 = directory + "orderParameter/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
        const std::string writeFile2 = directory + "orderParameter/average/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize);
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
            const std::string readFile = directory + "meanClusterSize/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
            readCSV(readFile, temp);
            average += temp;
        }
        average /= fileNum;
        ensembleSize = totalEnsemble;
        coreNum = 0;
        const std::string writeFile1 = directory + "meanClusterSize/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
        const std::string writeFile2 = directory + "meanClusterSize/average/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize);
        writeCSV(writeFile1, average);
        writeCSV(writeFile2, average);
    }

    //! Second Giant
    void average_secondGiant(){
        std::vector<double> average(networkSize);
        std::vector<double> temp(networkSize);
        for (int core=0; core<fileNum; ++core){
            coreNum = core;
            ensembleSize = ensembleList[core];
            const std::string readFile = directory + "secondGiant/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
            readCSV(readFile, temp);
            average += temp;
        }
        average /= fileNum;
        ensembleSize = totalEnsemble;
        coreNum = 0;
        const std::string writeFile1 = directory + "secondGiant/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
        const std::string writeFile2 = directory + "secondGiant/average/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize);
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
                const std::string readFile = directory + "interEventTime/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
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
            const std::string writeFile1 = directory + "interEventTime/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
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
        const std::string writeFile2 = directory + "interEventTime/average/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize);
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
        const std::string writeFile3 = directory + "interEventTime/logBin/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize);
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
                const std::string readFile = directory + "deltaAcceptance/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
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
            const std::string writeFile1 = directory + "deltaAcceptance/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
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
        const std::string writeFile2 = directory + "deltaAcceptance/average/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize);
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
        const std::string writeFile3 = directory + "deltaAcceptance/logBin/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize);
        writeCSV(writeFile3, result);
    }

    //! Order Parameter Distribution
    void average_orderParameterDistribution(){
        for (const double& t : time_orderParameterDistribution){
            std::map<double, double> average, temp;
            for (int core=0; core<fileNum; ++core){
                coreNum = core;
                ensembleSize = ensembleList[core];
                const std::string readFile = directory + "orderParameterDistribution/" + filename_time(networkSize, acceptanceThreshold, ensembleSize, t, coreNum);
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
            const std::string writeFile1 = directory + "orderParameterDistribution/" + filename_time(networkSize, acceptanceThreshold, ensembleSize, t, coreNum);
            const std::string writeFile2 = directory + "orderParameterDistribution/average/" + filename_time(networkSize, acceptanceThreshold, ensembleSize, t);
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
                const std::string readFile = directory + "clusterSizeDistribution/" + filename_orderParameter(networkSize, acceptanceThreshold, ensembleSize, op, coreNum);
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
            const std::string writeFile1 = directory + "clusterSizeDistribution/" + filename_orderParameter(networkSize, acceptanceThreshold, ensembleSize, op, coreNum);
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
            const std::string writeFile3 = directory + "clusterSizeDistribution/logBin/" + filename_orderParameter(networkSize, acceptanceThreshold, ensembleSize, op, coreNum);
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
                const std::string readFile = directory + "ageDistribution/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
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
            const std::string writeFile1 = directory + "ageDistribution/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
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
            const std::string writeFile3 = directory + "ageDistribution/" + state +"/logBin/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize);
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
                const std::string readFile = directory + "interEventTimeDistribution/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
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
            const std::string writeFile1 = directory + "interEventTimeDistribution/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
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
            const std::string writeFile3 = directory + "interEventTimeDistribution/" + state +"/logBin/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize);
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
                const std::string readFile = directory + "deltaUpperBoundDistribution/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
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
            const std::string writeFile1 = directory + "deltaUpperBoundDistribution/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
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
            const std::string writeFile3 = directory + "deltaUpperBoundDistribution/" + state +"/logBin/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize);
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
                const std::string readFile = directory + "deltaAcceptanceDistribution/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
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
            const std::string writeFile1 = directory + "deltaAcceptanceDistribution/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
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
            ensembleSize = totalEnsemble;
            const std::string writeFile3 = directory + "deltaAcceptanceDistribution/" + state + "/logBin/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize);
            writeCSV(writeFile3, result);
        }
    }


    //! Inter Event Time vs Delta Acceptance
    void average_interEventTime_DeltaAcceptance(){
        //* average
        std::map<int, double> average, temp;
        std::map<int, int> average_sampled;
        for (int core=0; core<fileNum; ++core){
            coreNum = core;
            ensembleSize = ensembleList[core];
            const std::string readFile = directory + "interEventTime_DeltaAcceptance/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
            readCSV(readFile, temp);
            for (auto it=temp.begin(); it!=temp.end(); ++it){
                average[it->first] += it->second;
                ++average_sampled[it->first];
            }
        }
        for (auto it=average.begin(); it!=average.end(); ++it){
            it->second /= average_sampled[it->first];
        }
        ensembleSize = totalEnsemble;
        coreNum = 0;
        const std::string writeFile1 = directory + "interEventTime_DeltaAcceptance/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
        writeCSV(writeFile1, average);

        //* Log Bin iet and log bin delta acceptance
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
        for (auto it=result.begin(); it!=result.end(); ++it){
            it->second /= sampled[it->first];
        }
        ensembleSize = totalEnsemble;
        const std::string writeFile3 = directory + "interEventTime_DeltaAcceptance/logBin/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize);
        writeCSV(writeFile3, result);
    }

    //! Upper Bound vs Delta Acceptance
    void average_upperBound_DeltaAcceptance(){
        //* average
        std::map<int, double> average, temp;
        std::map<int, int> average_sampled;
        for (int core=0; core<fileNum; ++core){
            coreNum = core;
            ensembleSize = ensembleList[core];
            const std::string readFile = directory + "upperBound_DeltaAcceptance/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
            readCSV(readFile, temp);
            for (auto it=temp.begin(); it!=temp.end(); ++it){
                average[it->first] += it->second;
                ++average_sampled[it->first];
            }
        }
        for (auto it=average.begin(); it!=average.end(); ++it){
            it->second /=average_sampled[it->first];
        }
        ensembleSize = totalEnsemble;
        coreNum = 0;
        const std::string writeFile1 = directory + "upperBound_DeltaAcceptance/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
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
        for (auto it=result.begin(); it!=result.end(); ++it){
            it->second /= sampled[it->first];
        }
        ensembleSize = totalEnsemble;
        const std::string writeFile3 = directory + "upperBound_DeltaAcceptance/logBin/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize);
        writeCSV(writeFile3, result);
    }


    //! Delta Upper Bound vs Delta Acceptance
    void average_deltaUpperbound_DeltaAcceptance(){
        //* average
        std::map<int, double> average, temp;
        std::map<int, int> average_sampled;
        for (int core=0; core<fileNum; ++core){
            coreNum = core;
            ensembleSize = ensembleList[core];
            const std::string readFile = directory + "deltaUpperBound_DeltaAcceptance/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
            readCSV(readFile, temp);
            for (auto it=temp.begin(); it!=temp.end(); ++it){
                average[it->first] += it->second;
                ++average_sampled[it->first];
            }
        }
        for (auto it=average.begin(); it!=average.end(); ++it){
            it->second /=average_sampled[it->first];
        }
        ensembleSize = totalEnsemble;
        coreNum = 0;
        const std::string writeFile1 = directory + "deltaUpperBound_DeltaAcceptance/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
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
        for (auto it=result.begin(); it!=result.end(); ++it){
            it->second /= sampled[it->first];
        }
        ensembleSize = totalEnsemble;
        const std::string writeFile3 = directory + "deltaUpperBound_DeltaAcceptance/logBin/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize);
        writeCSV(writeFile3, result);
    }
}//* End of namespace mBFW::data