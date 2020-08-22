#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <map>

#include "/pds/pds172/hoyun1009/library-Git/linearAlgebra.hpp"
#include "/pds/pds172/hoyun1009/library-Git/Networks.hpp"
#include "/pds/pds172/hoyun1009/library-Git/CSV.hpp"
#include "/pds/pds172/hoyun1009/library-Git/stringFormat.hpp"

#include "parameters.hpp"

namespace mBFW{
    //! Declaration of input variables
    int networkSize;
    int ensembleSize;
    int coreNum;
    int randomEngineSeed;
    double acceptanceThreshold;
    double precision;

    //! Declaration of pre-defined variables
    double degenerated;
    double m_c, t_c;
    std::vector<double> time_orderParameterDistribution, orderParameter_clusterSizeDistribution;
    std::vector<double> deltaAcceptanceExponent;
    std::vector<double> minDeltaAcceptance, maxDeltaAcceptance, logBinnedDeltaAcceptance;
    int deltaAcceptanceBinNum;
    const std::vector<std::string> states = {"before", "during"};
    const double m_a = 0.05;
    const double deltaAcceptanceBinDelta = 0.01;

    //! Declaration of output variables
    //* X[time] = value of observable X at time
    std::vector<double> orderParameter;
    std::vector<double> secondGiant;
    std::vector<double> meanClusterSize;
    std::map<std::string, std::vector<double>> interEventTime;
    std::map<std::string, std::vector<int>> sampledInterEventTime;
    std::map<std::string, std::vector<double>> deltaAcceptance;
    std::map<std::string, std::vector<int>> sampledDeltaAcceptance;

    //* orderParameterDistribution[time] : distribution of rounded order parameter at (rounded time)=(time)
    //* orderParameterDistribution[time][mcs] : number of cluster size of "maximum cluster size"
    std::map<double, std::vector<int>> orderParameterDistribution;

    //* clusterSizeDistribution[orderParameter] : cluster size distribution at (rounded order parameter)=(orderParameter)
    //* clusterSizeDistribution[orderParameter][cs] : number of cluster of size "cs"
    std::map<double, std::vector<int>> clusterSizeDistribution;

    //* interEventTimeDistribution["before"] : inter event time distribution before jump (order parameter < m_a)
    //* interEventTimeDistribution["during"] : inter event time distribution during jump (m_a < order parameter < m_c)
    std::map<std::string, std::vector<int>> interEventTimeDistribution;

    //* detaMDistribution["before"] : jump size of maximum cluster size distribution before jump
    //* detaMDistribution["during"] : jump size of maximum cluster size distribution during jump
    std::map<std::string, std::vector<int>> deltaUpperBoundDistribution;

    //* ageDistribution["before"] : age distribution before jump (order parameter < m_a)
    //* ageDistribution["during"] : age distribution during jump (m_a < order parameter < m_c)
    std::map<std::string, std::vector<int>> ageDistribution;

    //* deltaAcceptanceDistribution["before"] : max acceptance distribution befeore jump (order paramter < m_a)
    //* deltaAcceptanceDistribution["during"] : max acceptance distribution during jump (m_a < order paramter < m_c)
    std::map<std::string, std::vector<int>> deltaAcceptanceDistribution;

    //* X_maxAcceptance[x] : average max acceptance at X=x before jump
    //* sampledX_maxAcceptance[x] : number of samples for each X=x
    std::vector<double> interEventTime_DeltaAcceptance;
    std::vector<int> sampledInterEventTime_DeltaAcceptance;
    std::vector<double> upperBound_DeltaAcceptance;
    std::vector<int> sampledUpperBound_DeltaAcceptance;
    std::vector<double> deltaUpperBound_DeltaAcceptance;
    std::vector<int> sampledDeltaUpperBound_DeltaAcceptance;

    //* Dynamics
    std::vector<std::vector<int>> dynamics;

    //! Random Engine
    pcg32 randomEngine;
    std::uniform_int_distribution<int> nodeDistribution;

    void setParameters(const int& t_networkSize, const int& t_ensembleSize, const double& t_acceptanceThreshold, const double t_precision, const int& t_coreNum, const int& t_randomEngineSeed){
        //! Input variables
        networkSize = t_networkSize;
        ensembleSize = t_ensembleSize;
        coreNum = t_coreNum;
        acceptanceThreshold = t_acceptanceThreshold;
        randomEngineSeed = t_randomEngineSeed;
        t_networkSize < t_precision ? precision = t_networkSize : precision = t_precision;

        //! Pre-defined variables
        std::tie(time_orderParameterDistribution, orderParameter_clusterSizeDistribution, m_c, t_c) = getParameters(networkSize, acceptanceThreshold);
        degenerated=t_networkSize/t_precision;

        deltaAcceptanceExponent = linearAlgebra::arange(-8,0,deltaAcceptanceBinDelta);
        deltaAcceptanceBinNum = deltaAcceptanceExponent.size()-1;
        minDeltaAcceptance.resize(deltaAcceptanceBinNum);
        maxDeltaAcceptance.resize(deltaAcceptanceBinNum);
        logBinnedDeltaAcceptance.resize(deltaAcceptanceBinNum);
        for (int i=0; i<deltaAcceptanceBinNum; ++i){
            minDeltaAcceptance[i] = pow(10, deltaAcceptanceExponent[i]);
            maxDeltaAcceptance[i] = pow(10, deltaAcceptanceExponent[i+1]);
            logBinnedDeltaAcceptance[i] = sqrt(minDeltaAcceptance[i] * maxDeltaAcceptance[i]);
        }

        //! Output variables (Observables)
        //* time-X
        orderParameter.resize(t_networkSize,1);
        secondGiant.resize(t_networkSize,1);
        meanClusterSize.resize(t_networkSize,1);
        for (auto state : states){
            interEventTime[state].resize(t_networkSize);
            sampledInterEventTime[state].resize(t_networkSize);
            deltaAcceptance[state].resize(t_networkSize);
            sampledDeltaAcceptance[state].resize(t_networkSize);
        }

        //* Distributions
        for (const double& t : time_orderParameterDistribution){
            orderParameterDistribution[t].resize(t_networkSize);
        }
        for (const double& m : orderParameter_clusterSizeDistribution){
            clusterSizeDistribution[m].resize(t_networkSize);
        }
        for (auto state : states){
            interEventTimeDistribution[state].resize(t_networkSize);
            deltaUpperBoundDistribution[state].resize(t_networkSize);
            ageDistribution[state].resize(t_networkSize);
            deltaAcceptanceDistribution[state].resize(deltaAcceptanceBinNum);
        }

        //* X-DeltaAcceptance
        interEventTime_DeltaAcceptance.resize(t_networkSize);
        sampledInterEventTime_DeltaAcceptance.resize(t_networkSize);
        upperBound_DeltaAcceptance.resize(t_networkSize);
        sampledUpperBound_DeltaAcceptance.resize(t_networkSize);
        deltaUpperBound_DeltaAcceptance.resize(t_networkSize);
        sampledDeltaUpperBound_DeltaAcceptance.resize(t_networkSize);

        //* Dynamics
        dynamics.reserve(3*networkSize*ensembleSize);

        //! Random Engine
        randomEngineSeed == -1 ? randomEngine.seed((std::random_device())()) : randomEngine.seed(randomEngineSeed);
        nodeDistribution.param(std::uniform_int_distribution<int>::param_type(0, networkSize-1));

    } //* End of setParameters

    void run(){
        for (int ensemble=0; ensemble<ensembleSize; ++ensemble){
            //* Default values for one ensemble
            NZ_Network model(networkSize);
            Node node1, node2, root1, root2;
            int size1, size2;
            int time = 0;
            int trialTime = 0;
            int upperBound = 2;
            int periodTime = 0;
            int periodTrialTime = 0;
            int peakTime = 0;
            int peakTrialTime = 0;
            int updatedTime = 0; 
            double maxAcceptance = 0;
            double maxPeriodAcceptance = 0;
            bool findNewNodes = true;

            //! Dynamics only for small number of ensembles
            if (ensembleSize < 5){
                dynamics.emplace_back(std::vector<int>{time, trialTime, periodTime, periodTrialTime, model.getMaximumClusterSize(), upperBound});
            }

            //* Do rBFW algorithm until all clusters merge to one
            while (model.getMaximumClusterSize() < networkSize){
                //* Find new nodes
                if (findNewNodes){
                    //* Randomly choose new nodes
                    do {
                        node1 = nodeDistribution(randomEngine);
                        node2 = nodeDistribution(randomEngine);
                        root1 = model.getRoot(node1);
                        root2 = model.getRoot(node2);
                    } while(root1 == root2);

                    //* choose two clusters of each node
                    size1 = model.getClusterSize(root1);
                    size2 = model.getClusterSize(root2);
                }

                //* Merge two clusters, update time
                if (size1+size2 <= upperBound){
                    model.merge(root1, root2);
                    ++time;
                    ++trialTime;
                    ++periodTime;
                    ++periodTrialTime;
                    findNewNodes=true;                    
                    const double exactOrderParameter=model.getMaximumClusterSize()/(double)networkSize;

                    //* max acceptance
                    if ((double)time/trialTime > maxAcceptance){
                        peakTime = time;
                        peakTrialTime = trialTime;
                        maxAcceptance = (double)time/trialTime;
                    }

                    //! order parameter
                    orderParameter[time] += exactOrderParameter;

                    //! mean cluster size
                    meanClusterSize[time] += model.getMeanClusterSize();

                    //! second giant
                    secondGiant[time] += model.getSecondMaximumClusterSize()/(double)networkSize;

                    //! Age Distribution
                    if (exactOrderParameter < m_a){
                        for (auto& changedAge : model.getChangedAge()){
                            ageDistribution["before"][changedAge.first] += changedAge.second;
                        }
                    }
                    else if (exactOrderParameter < m_c){
                        for (auto& changedAge : model.getChangedAge()){
                            ageDistribution["during"][changedAge.first] += changedAge.second;
                        }
                    }

                    //! order parameter distribution
                    const double roundedTime=round(time/degenerated)/precision;
                    auto it = std::find(time_orderParameterDistribution.begin(), time_orderParameterDistribution.end(), roundedTime);
                    if (it != time_orderParameterDistribution.end()){
                        ++orderParameterDistribution[*it][model.getMaximumClusterSize()];
                    }

                    //* End of k-period
                    if (model.getDeltaMaximumClusterSize()){
                        const int currentInterEventTime = time-updatedTime;
                        const double currentDeltaAcceptance = maxAcceptance-acceptanceThreshold;
                        std::string currentState;

                        //! Cluster Size Distribution 
                        const double roundedOrderParameter=round(exactOrderParameter*precision)/precision;
                        auto it = std::find(orderParameter_clusterSizeDistribution.begin(), orderParameter_clusterSizeDistribution.end(), roundedOrderParameter);
                        if (it != orderParameter_clusterSizeDistribution.end()){
                            auto sortedCluster = model.getSortedCluster();
                            for (auto it2 = sortedCluster.begin(); it2!= sortedCluster.end(); ++it2){
                                clusterSizeDistribution[*it][it2->first] += it2->second;
                            }    
                        }

                        //* Before and During Jump
                        if (exactOrderParameter < m_c){
                            exactOrderParameter < m_a ? currentState = "before" : currentState = "during";

                            //! Inter Event Time Distribution
                            ++interEventTimeDistribution[currentState][currentInterEventTime];

                            //! Inter Event Time
                            interEventTime[currentState][time] += currentInterEventTime;
                            ++sampledInterEventTime[currentState][time];

                            //! Delta Acceptance
                            deltaAcceptance[currentState][time] += currentDeltaAcceptance;
                            ++sampledDeltaAcceptance[currentState][time];

                            //! Delta Upper Bound Distribution
                            ++deltaUpperBoundDistribution[currentState][model.getDeltaMaximumClusterSize()];

                            //! Delta Acceptance Distribution
                            for (int i=0; i<deltaAcceptanceBinNum; ++i){
                                if (minDeltaAcceptance[i] <= currentDeltaAcceptance && currentDeltaAcceptance < maxDeltaAcceptance[i]){
                                    ++deltaAcceptanceDistribution[currentState][i];
                                    break;
                                }
                            }

                            //! X vs Delta Acceptance
                            if (currentState == "before"){
                                interEventTime_DeltaAcceptance[currentInterEventTime] += currentDeltaAcceptance;
                                ++sampledInterEventTime_DeltaAcceptance[currentInterEventTime];
                                upperBound_DeltaAcceptance[upperBound] += currentDeltaAcceptance;
                                ++sampledUpperBound_DeltaAcceptance[upperBound];
                                deltaUpperBound_DeltaAcceptance[model.getDeltaMaximumClusterSize()] += currentDeltaAcceptance;
                                ++sampledDeltaUpperBound_DeltaAcceptance[model.getDeltaMaximumClusterSize()];
                            }
                        }

                        //* Initialize variable for new k-period
                        updatedTime = time;
                        periodTime = 0;
                        periodTrialTime = 0;
                        maxAcceptance=0;
                        maxPeriodAcceptance=0;
                    }//* End of updating k-period
                }
                else if ((double)time/trialTime < acceptanceThreshold){
                    ++upperBound;
                    findNewNodes=false;
                }
                else{
                    ++trialTime;
                    ++periodTrialTime;
                    findNewNodes=true;
                } 
             
                //! Dynamics only for small number of ensembles
                if (ensembleSize < 5){
                    dynamics.emplace_back(std::vector<int>{time, trialTime, periodTime, periodTrialTime, model.getMaximumClusterSize(), upperBound});
                }//* End of one step
            
            }//* End of network growing (one ensemble)
        } //* End of every ensembles
    } //* End of function run



    namespace process{
        //* Default directory for data of mBFW model
        const std::string directory = "/pds/pds172/hoyun1009/data/mBFW/";
        using namespace linearAlgebra;

        //* File names of observables
        std::string defaultFileName(const bool withoutCoreNum = false){
            std::string filename;
            if (withoutCoreNum){
                return "N"+to_stringWithExponent((double)networkSize, 1)+",G"+to_stringWithPrecision(acceptanceThreshold,1)+",E"+std::to_string(ensembleSize)+".txt";
            }
            else{
                return "N"+to_stringWithExponent((double)networkSize, 1)+",G"+to_stringWithPrecision(acceptanceThreshold,1)+",E"+std::to_string(ensembleSize)+"-"+std::to_string(coreNum)+".txt";
            }
        }

        std::string filename_time(const double& t_time, const bool withoutCoreNum = false){
            if (withoutCoreNum){
                return "N"+to_stringWithExponent((double)networkSize, 1)+",G"+to_stringWithPrecision(acceptanceThreshold,1)+",E"+std::to_string(ensembleSize)+",T"+to_stringWithPrecision(t_time,4)+".txt";
            }
            else{
                return "N"+to_stringWithExponent((double)networkSize, 1)+",G"+to_stringWithPrecision(acceptanceThreshold,1)+",E"+std::to_string(ensembleSize)+",T"+to_stringWithPrecision(t_time,4)+"-"+std::to_string(coreNum)+".txt";

            }
        }

        std::string filename_orderParameter(const double& t_orderParameter, const bool withoutCoreNum = false){
            if (withoutCoreNum){
                return "N"+to_stringWithExponent((double)networkSize, 1)+",G"+to_stringWithPrecision(acceptanceThreshold,1)+",E"+std::to_string((double)ensembleSize)+",OP"+to_stringWithPrecision(t_orderParameter,4)+".txt";

            }
            else{
                return "N"+to_stringWithExponent((double)networkSize, 1)+",G"+to_stringWithPrecision(acceptanceThreshold,1)+",E"+std::to_string((double)ensembleSize)+",OP"+to_stringWithPrecision(t_orderParameter,4)+"-"+std::to_string(coreNum)+".txt";
            }
        }

        //! Order Parameter
        void save_orderParameter(){
            const std::string fullFileName = directory + "orderParameter/" + defaultFileName();
            orderParameter /= ensembleSize;
            writeCSV(fullFileName, orderParameter);
        }

        //! Mean Cluster Size
        void save_meanClusterSize(){
            const std::string fullFileName = directory + "meanClusterSize/" + defaultFileName();
            meanClusterSize /= ensembleSize;
            writeCSV(fullFileName, meanClusterSize);
        }

        //! Second Giant
        void save_secondGiant(){
            const std::string fullFileName = directory + "secondGiant/" + defaultFileName();
            secondGiant /= ensembleSize;
            writeCSV(fullFileName, secondGiant);
        }

        //! Inter Event Time
        void save_interEventTime(){
            for (const auto& state : states){
                const std::string fullFileName = directory + "interEventTime/" + state + "/" + defaultFileName();

                //* save only useful data
                std::map<double, double> trimmed;
                for (int t=0; t<networkSize; ++t){
                    if (sampledInterEventTime[state][t]){
                        trimmed[(double)t/networkSize] = interEventTime[state][t]/sampledInterEventTime[state][t];
                    }
                }
                writeCSV(fullFileName, trimmed);
            }
        }

        //! Delta Acceptance
        void save_deltaAcceptance(){
            for (const auto& state : states){
                const std::string fullFileName = directory + "deltaAcceptance/" + state + "/" + defaultFileName();

                //* save only useful data
                std::map<double, double> trimmed;
                for (int t=0; t<networkSize; ++t){
                    if (sampledDeltaAcceptance[state][t]){
                        trimmed[(double)t/networkSize] = deltaAcceptance[state][t]/sampledDeltaAcceptance[state][t];
                    }
                }
                writeCSV(fullFileName, trimmed);
            }
        }

        //! Order Parameter Distribution
        void save_orderParameterDistribution(){
            for (const auto& t : time_orderParameterDistribution){
                const std::string fullFileName = directory + "orderParameterDistribution/" + filename_time(t);

                //* save only useful data
                std::map<double, double> trimmed;
                const double tot = std::accumulate(orderParameterDistribution[t].begin(), orderParameterDistribution[t].end(), 0);
                for (int mcs=0; mcs<networkSize; ++mcs){
                    if (orderParameterDistribution[t][mcs]){
                        trimmed[(double)mcs/networkSize] = orderParameterDistribution[t][mcs]/tot;
                    }
                }
                writeCSV(fullFileName, trimmed);
            }
        }

        //! Cluster Size Distribution
        void save_clusterSizeDistribution(){
            for (const auto& op : orderParameter_clusterSizeDistribution){
                const std::string fullFileName = directory + "clusterSizeDistribution/" + filename_orderParameter(op);
            
                //* save only useful data
                std::map<double, double> trimmed;
                const double tot = std::accumulate(clusterSizeDistribution[op].begin(), clusterSizeDistribution[op].end(), 0.0);
                for (int cs=0; cs<networkSize; ++cs){
                    if (clusterSizeDistribution[op][cs]){
                        trimmed[(double)cs/networkSize] = clusterSizeDistribution[op][cs]/tot;
                    }
                }
                writeCSV(fullFileName, trimmed);
            }
        }

        //! Age Distribution
        void save_ageDistribution(){
            for (const auto& state : states){
                const std::string fullFileName = directory + "ageDistribution/" + state + "/" + defaultFileName();

                //* save only useful data
                std::map<int, double> trimmed;
                const double tot = std::accumulate(ageDistribution[state].begin(), ageDistribution[state].end(), 0.0);
                for (int age=0; age<networkSize; ++age){
                    if (ageDistribution[state][age]){
                        trimmed[age] = ageDistribution[state][age]/tot;
                    }
                }
                writeCSV(fullFileName, trimmed);
            }
        }

        //! Inter Event Time Distribution
        void save_interEventTimeDistribution(){
            for (const auto& state : states){
                const std::string fullFileName = directory + "interEventTimeDistribution/" + state + "/" + defaultFileName();

                //* save only useful data
                std::map<int, double> trimmed;
                const double tot = std::accumulate(interEventTimeDistribution[state].begin(), interEventTimeDistribution[state].end(), 0.0);
                for (int iet=0; iet<networkSize; ++iet){
                    if (interEventTimeDistribution[state][iet]){
                        trimmed[iet] = interEventTimeDistribution[state][iet]/tot;
                    }
                }
                writeCSV(fullFileName, trimmed);
            }
        }

        //! Delta Upper Bound Distribution
        void save_deltaUpperBoundDistribution(){
            for (const auto& state : states){
                const std::string fullFileName = directory + "deltaUpperBoundDistribution/" + state + "/" + defaultFileName();

                //* save only useful data
                std::map<int, double> trimmed;
                const double tot = std::accumulate(deltaUpperBoundDistribution[state].begin(), deltaUpperBoundDistribution[state].end(), 0.0);
                for (int deltaK=0; deltaK<networkSize; ++deltaK){
                    if (deltaUpperBoundDistribution[state][deltaK]){
                        trimmed[deltaK] = deltaUpperBoundDistribution[state][deltaK]/tot;
                    }
                }
                writeCSV(fullFileName, trimmed);
            }
        }

        //! Delta Acceptance Distribution
        void save_deltaAcceptanceDistribution(){
            for (const auto& state : states){
                const std::string fullFileName = directory + "deltaAcceptanceDistribution/" + state + "/" + defaultFileName();

                //* save only useful data
                std::map<double, double> trimmed;
                const double tot = std::accumulate(deltaAcceptanceDistribution[state].begin(), deltaAcceptanceDistribution[state].end(), 0.0);
                for (int i=0; i<deltaAcceptanceBinNum; ++i){
                    if (deltaAcceptanceDistribution[state][i]){
                        trimmed[logBinnedDeltaAcceptance[i]] = deltaAcceptanceDistribution[state][i]/tot;
                    }
                }
                writeCSV(fullFileName, trimmed);
            }
        }

        //! Inter Event Time vs Delta Acceptance
        void save_interEventTime_DeltaAcceptance(){
            const std::string fullFileName = directory + "interEventTime_DeltaAcceptance/" + defaultFileName();

            //* save only useful data
            std::map<int, double> trimmed;
            for (int iet=0; iet<networkSize; ++iet){
                if (sampledInterEventTime_DeltaAcceptance[iet]){
                    trimmed[iet] = interEventTime_DeltaAcceptance[iet]/sampledInterEventTime_DeltaAcceptance[iet];
                }
            }
            writeCSV(fullFileName, trimmed);
        }

        //! Upper Bound vs Delta Acceptance
        void save_upperBound_DeltaAcceptance(){
            const std::string fullFileName = directory + "upperBound_DeltaAcceptance/" + defaultFileName();

            //* save only useful data
            std::map<int, double> trimmed;
            for (int k=0; k<networkSize; ++k){
                if (sampledUpperBound_DeltaAcceptance[k]){
                    trimmed[k] = upperBound_DeltaAcceptance[k]/sampledUpperBound_DeltaAcceptance[k];
                }
            }
            writeCSV(fullFileName, trimmed);
        }


        //! Delta Upper Bound vs Delta Acceptance
        void save_deltaUpperBound_DeltaAcceptance(){
            const std::string fullFileName = directory + "deltaUpperBound_DeltaAcceptance/" + defaultFileName();

            //* save only useful data
            std::map<int, double> trimmed;
            for (int deltaK=0; deltaK<networkSize; ++deltaK){
                if (sampledDeltaUpperBound_DeltaAcceptance[deltaK]){
                    trimmed[deltaK] = deltaUpperBound_DeltaAcceptance[deltaK]/sampledDeltaUpperBound_DeltaAcceptance[deltaK];
                }
            }
            writeCSV(fullFileName, trimmed);
        }

        //! Dynamics
        void save_dynamics(){
            const std::string fullFileName = directory + "dynamics/" + defaultFileName();
            writeCSV(fullFileName, dynamics);
        }
    } //* End of namespace mBFW::process
} //* End of namespace mBFW