#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <cstdio>

#include "../library-Git/CSV.hpp"
#include "../library-Git/linearAlgebra.hpp"

#include "parameters.hpp"
#include "mBFW.hpp"

namespace mBFW::data{
    using namespace linearAlgebra;

    //*------------------------------------------- Declaration of varaibles used at mBFW::data ------------------------------------------------------
    std::vector<int> ensembleList;
    int fileNum;
    int totalEnsemble;
    bool deleteFile;
    int logBinNum;
    std::vector<double> double_min, double_value;
    std::vector<double> int_min, int_value;

    //*------------------------------------------- Set parameters for average, merge, log bin ------------------------------------------------------
    //* Set parameters for core average
    void setParameters(const int& t_networkSize, const double& t_acceptanceThreshold, const std::vector<int>& t_ensembleList, const double& t_logBinDelta, const std::vector<bool> t_observables, const bool t_deleteFile){
        //! Observables to be processed
        process_orderParameter = t_observables[0];
        process_meanClusterSize = t_observables[1];
        process_secondGiant = t_observables[2];
        process_interEventTime = t_observables[3];
        process_deltaAcceptance = t_observables[4];
        process_orderParameterDistribution = t_observables[5];
        process_clusterSizeDistribution = t_observables[6];
        process_ageDistribution = t_observables[7];
        process_interEventTimeDistribution = t_observables[8];
        process_deltaUpperBoundDistribution = t_observables[9];
        process_deltaAcceptanceDistribution = t_observables[10];
        process_interEventTime_DeltaAcceptance = t_observables[11];
        process_upperBound_DeltaAcceptance = t_observables[12];
        process_deltaUpperBound_DeltaAcceptance = t_observables[13];
        process_dynamics = t_observables[14];

        //! Input variables
        networkSize = t_networkSize;
        acceptanceThreshold = t_acceptanceThreshold;
        ensembleList = t_ensembleList;
        deleteFile = t_deleteFile;
        fileNum = t_ensembleList.size();
        totalEnsemble = std::accumulate(t_ensembleList.begin(), t_ensembleList.end(), 0);

        //! Log Binning
        const std::vector<double> exponent = arange(-8.0, 0.0, t_logBinDelta);
        logBinNum = exponent.size()-1;
        double_value.resize(logBinNum);
        int_value.resize(logBinNum);

        double_min = elementPow(10.0, exponent);
        int_min = double_min * 1e8;
        for (int i=0; i<logBinNum; ++i){
            double_value[i] = sqrt(double_min[i] * double_min[i+1]);
            int_value[i] = sqrt(int_min[i] * int_min[i+1]);
        }

        //! pre-defined parameters from "parameter.hpp"
        std::tie(m_c, t_c, time_orderParameterDistribution, orderParameter_clusterSizeDistribution) = mBFW::parameters::pre_defined(networkSize, acceptanceThreshold);
    }

    //*------------------------------------------- functions for data process ------------------------------------------------------
    //! average process
    //* t_c : only for type check
    template<typename T>
    const std::map<T, double> average(const std::string t_directory, const T& t_check){
        std::map<T, double> average, temp;
        std::map<T, int> sampledAverage;
        for (int core=0; core<fileNum; ++core){
            const std::string readFile = t_directory + defaultFileName(networkSize, acceptanceThreshold, ensembleList[core], core);
            readCSV(readFile, temp);
            average += temp;
            sampleNum(sampledAverage, temp);
            if (deleteFile){std::remove(readFile.c_str());}
        }
        average /= sampledAverage;
        return average;
    }

    //! average process of distribution
    //* for check point distribution, t_checkpoint
    //* t_check : only for type
    template <typename T>
    const std::map<T, double> averageDistribution(const std::string& t_observable, const std::string& t_directory, const T& t_check, const double& t_checkpoint){
        std::map<T, double> average, temp;
        for (int core=0; core<fileNum; ++core){
            const std::string readFile = t_directory + generalFileName(t_observable, networkSize, acceptanceThreshold, ensembleList[core], t_checkpoint, core);
            readCSV(readFile, temp);
            average += temp;
            if (deleteFile){std::remove(readFile.c_str());}
        }
        average /= fileNum;
        return average;
    }

    //! Log Binning
    template <typename T>
    const std::map<double, double> logBin(const std::map<T,double>& t_raw){
        //* Test whether T is double or int
        T test = 2;
        std::vector<double> min, value;
        if (3/test == 1.5){
            min = double_min;
            value = double_value;
        }
        else{
            min = int_min;
            value = int_value;
        }

        //* Log Binning
        std::map<double, double> binned;
        std::map<double, int> sampledLogBin;
        for (auto it=t_raw.begin(); it!=t_raw.end(); ++it){
            for (int i=0; i<logBinNum; ++i){
                if (min[i+1] > it->first && it->first){
                    binned[value[i]] += it->second;
                    ++sampledLogBin[value[i]];
                    break;
                }
            }
        }
        binned /= sampledLogBin;
        return binned;

    }

    //! time vs X
    void time_X(const std::string& t_observable){
        const std::string directory = rootDirectory + t_observable + "/";
        std::vector<double> average(networkSize);
        std::vector<double> temp(networkSize);
        for (int core=0; core<fileNum; ++core){
            const std::string readFile = directory + defaultFileName(networkSize, acceptanceThreshold, ensembleList[core], core);
            readCSV(readFile, temp);
            average += temp;
            if (deleteFile){std::remove(readFile.c_str());}
        }
        average /= fileNum;
        const std::string writeFile1 = directory + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble, 0);
        const std::string writeFile2 = directory + "average/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble);
        writeCSV(writeFile1, average);
        writeCSV(writeFile2, average);
    }

    //! time vs X with log binning w.r.t t_c-t
    void logBin_time_X(const std::string& t_observable){
        const std::string directory = rootDirectory + t_observable + "/";

        std::map<double, double> merge;
        std::map<double, int> sampledMerge;
        for (auto state : states){
            //* average
            const std::map<double, double> avg = average(directory + state + "/", 0.0);
            const std::string writeFile1 = directory + state + "/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble, 0);
            writeCSV(writeFile1, avg);

            //* merge
            merge += avg;
            sampleNum(sampledMerge, avg);
        }
        merge /= sampledMerge;
        const std::string writeFile3 = directory + "merge/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble);
        writeCSV(writeFile3, merge);

        //* log bin w.r.t t_c-t
        merge = minus_first(t_c, merge);
        const std::map<double, double> binned = logBin(merge);
        const std::string writeFile2 = directory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble);
        writeCSV(writeFile2, binned);
    }


    //! Check point distribution
    //* t_c : only for type check
    template <typename T>
    void checkPointDistribution(const std::string& t_observable, const T& t_check){
        const std::string directory = rootDirectory + t_observable + "/";
        std::set<double> checkPointList;
        if (t_observable == "orderParameterDistribution"){
            checkPointList = time_orderParameterDistribution;
        }
        else if (t_observable == "clusterSizeDistribution"){
            checkPointList = orderParameter_clusterSizeDistribution;
        }

        for (const double& checkPoint : checkPointList){
            //* average
            const std::map<T, double> avg = averageDistribution(t_observable, directory, t_check, checkPoint);
            const std::string writeFile1 = directory + generalFileName(t_observable, networkSize, acceptanceThreshold, totalEnsemble, checkPoint, 0);
            writeCSV(writeFile1, avg);

            //* Log Binning
            if (t_observable == "orderParameterDistribution"){
                const std::string writeFile2 = directory + "average/" + filename_time(networkSize, acceptanceThreshold, totalEnsemble, checkPoint);
                writeCSV(writeFile2, avg);
            }
            else if (t_observable == "clusterSizeDistribution"){
                const std::map<double, double> binned = logBin(avg);
                const std::string writeFile2 = directory + "logBin/" + filename_orderParameter(networkSize, acceptanceThreshold, totalEnsemble, checkPoint);
                writeCSV(writeFile2, binned);
            }
        }
    }

    //! Distribution of 'keys' distinguished by before and during jump
    //* t_c : only for type check
    template <typename T>
    void distribution(const std::string& t_observable, const T& t_check){
        for (auto state : states){
            const std::string directory = rootDirectory + t_observable + "/" + state + "/";

            //* average
            const std::map<T, double> avg = averageDistribution(t_observable, directory, t_check, 0.0);
            const std::string writeFile1 = directory + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble, 0);
            writeCSV(writeFile1, avg);

            //* Log binning
            const std::map<double, double> binned = logBin(avg);
            const std::string writeFile2 = directory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble);
            writeCSV(writeFile2, binned);
        }
    }

    //! X vs Delta Acceptance
    //* t_c : only for type check
    template <typename T>
    void X_deltaAcceptance(const std::string&  t_observable, const T& t_check){
        const std::string directory = rootDirectory + t_observable + "/";

        //* average
        const std::map<T, double> avg = average(directory, t_check);
        const std::string writeFile1 = directory + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble, 0);
        writeCSV(writeFile1, avg);

        //* Log Binning
        const std::map<double, double> binned = logBin(avg);
        const std::string writeFile2 = directory + "logBin/" + defaultFileName(networkSize, acceptanceThreshold, totalEnsemble);
        writeCSV(writeFile2, binned);

    }

    //*------------------------------------------- process the data ------------------------------------------------------
    void run(){
        //! Order Parameter
        if (process_orderParameter){time_X("orderParameter");}

        //! Mean Cluster Size
        if (process_meanClusterSize){time_X("meanClusterSize");}

        //! Second giant
        if (process_secondGiant){time_X("secondGiant");}

        //! Inter Event Time
        if (process_interEventTime){logBin_time_X("interEventTime");}

        //! Delta Acceptance
        if (process_deltaAcceptance){logBin_time_X("deltaAcceptance");}

        //! Order Parameter Distribution
        if (process_orderParameterDistribution){checkPointDistribution("orderParameterDistribution", 0.0);}

        //! Cluster Size Distribution
        if (process_clusterSizeDistribution){checkPointDistribution("clusterSizeDistribution", 0);}

        //! Age Distribution
        if (process_ageDistribution){distribution("ageDistribution", 0);}

        //! Inter Event Time Distribution
        if (process_interEventTimeDistribution){distribution("interEventTimeDistribution", 0);}

        //! Delta Upper bound Distribution
        if (process_deltaUpperBoundDistribution){distribution("deltaUpperBoundDistribution", 0);}

        //! Delta Acceptance Distribution
        if (process_deltaAcceptanceDistribution){distribution("deltaAcceptanceDistribution", 0.0);}

        //! Inter Event Time vs Delta Acceptance
        if (process_interEventTime_DeltaAcceptance){X_deltaAcceptance("interEventTime_DeltaAcceptance", 0);}

        //! Upper Bound vs Delta Acceptance
        if (process_upperBound_DeltaAcceptance){X_deltaAcceptance("upperBound_DeltaAcceptance", 0);}

        //! Delta Upper Bound vs Delta Acceptance
        if (process_deltaUpperBound_DeltaAcceptance){X_deltaAcceptance("deltaUpperBound_DeltaAcceptance", 0);}
    }
}//* End of namespace mBFW::data