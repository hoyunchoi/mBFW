#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>

#include "../library-Git/linearAlgebra.hpp"
#include "../library-Git/Networks.hpp"
#include "../library-Git/CSV.hpp"
#include "../library-Git/pcg_random.hpp"

#include "parameters.hpp"
#include "mBFW.hpp"

namespace mBFW::generate{
    //! Declaration of variables only used at mBFW::generate namespace
    int randomEngineSeed;
    double precision;
    double degenerated;

    //* variables for log binning
    std::vector<double> minDeltaAcceptance;
    std::vector<double> maxDeltaAcceptance;
    std::vector<double> logBinnedDeltaAcceptance;
    int deltaAcceptanceBinNum;

    //! Declaration of observables
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

        const double deltaAcceptanceBinDelta = 0.01;
        const std::vector<double> deltaAcceptanceExponent = linearAlgebra::arange(-8,0,deltaAcceptanceBinDelta);
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
        orderParameter.resize(t_networkSize);
        secondGiant.resize(t_networkSize);
        meanClusterSize.resize(t_networkSize);
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
    } //* End of function mBFW::generate::setParameters

    void run( const bool t_orderParameter, const bool t_meanClusterSize, const bool t_secondGiant, const bool t_interEventTime, const bool t_deltaAcceptance, const bool t_orderParameterDistribution, const bool t_clusterSizeDistribution, const bool t_ageDistribution, const bool t_interEventTimeDistribution, const bool t_deltaUpperBoundDistribution, const bool t_deltaAcceptanceDistribution, const bool t_interEventTime_DeltaAcceptance, const bool t_upperBound_DeltaAcceptance, const bool t_deltaUpperBound_DeltaAcceptance, const bool t_dynamics){
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

            //* initial condition
            orderParameter[0] += 1.0/networkSize;
            secondGiant[0] += 1.0/networkSize;
            meanClusterSize[0] += 1.0;

            //! Dynamics only for small number of ensembles
            if (t_dynamics && ensembleSize < 5){
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
                    if (t_orderParameter){
                    orderParameter[time] += exactOrderParameter;

                    }

                    //! mean cluster size
                    if (t_meanClusterSize){
                        meanClusterSize[time] += model.getMeanClusterSize();
                    }

                    //! second giant
                    if (t_secondGiant){
                        secondGiant[time] += model.getSecondMaximumClusterSize()/(double)networkSize;
                    }

                    //! Age Distribution
                    if (t_ageDistribution){
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
                    }

                    //! order parameter distribution
                    if (t_orderParameterDistribution){
                        const double roundedTime=round(time/degenerated)/precision;
                        auto it = std::find(time_orderParameterDistribution.begin(),    time_orderParameterDistribution.end(), roundedTime);
                        if (it != time_orderParameterDistribution.end()){
                            ++orderParameterDistribution[*it][model.getMaximumClusterSize()];
                        }
                    }

                    //* End of k-period
                    if (model.getDeltaMaximumClusterSize() && model.getMaximumClusterSize()>2){
                        const int currentInterEventTime = time-updatedTime;
                        const double currentDeltaAcceptance = maxAcceptance-acceptanceThreshold;
                        std::string currentState;

                        //! Cluster Size Distribution
                        if (t_clusterSizeDistribution){
                            const double roundedOrderParameter=round(exactOrderParameter*precision)/precision;
                            auto it = std::find(orderParameter_clusterSizeDistribution.begin(),     orderParameter_clusterSizeDistribution.end(), roundedOrderParameter);
                            if (it != orderParameter_clusterSizeDistribution.end()){
                                auto sortedCluster = model.getSortedCluster();
                                for (auto it2 = sortedCluster.begin(); it2!= sortedCluster.end(); ++it2){
                                    clusterSizeDistribution[*it][it2->first] += it2->second;
                                }
                            }
                        }

                        //* Before and During Jump
                        if (exactOrderParameter < m_c){
                            exactOrderParameter < m_a ? currentState = "before" : currentState = "during";

                            //! Inter Event Time Distribution
                            if (t_interEventTimeDistribution){
                                ++interEventTimeDistribution[currentState][currentInterEventTime];
                            }

                            //! Inter Event Time
                            if (t_interEventTime){
                                interEventTime[currentState][time] += currentInterEventTime;
                                ++sampledInterEventTime[currentState][time];
                            }

                            //! Delta Acceptance
                            if (t_deltaAcceptance){
                                deltaAcceptance[currentState][time] += currentDeltaAcceptance;
                                ++sampledDeltaAcceptance[currentState][time];
                            }

                            //! Delta Upper Bound Distribution
                            if (t_deltaUpperBoundDistribution){
                                ++deltaUpperBoundDistribution[currentState][model.getDeltaMaximumClusterSize()];
                            }

                            //! Delta Acceptance Distribution
                            if (t_deltaAcceptanceDistribution){
                                for (int i=0; i<deltaAcceptanceBinNum; ++i){
                                    if (minDeltaAcceptance[i] <= currentDeltaAcceptance && currentDeltaAcceptance < maxDeltaAcceptance[i]){
                                        ++deltaAcceptanceDistribution[currentState][i];
                                        break;
                                    }
                                }
                            }

                            //! X vs Delta Acceptance
                            if (currentState == "before"){
                                if (t_interEventTime_DeltaAcceptance){
                                    interEventTime_DeltaAcceptance[currentInterEventTime] += currentDeltaAcceptance;
                                    ++sampledInterEventTime_DeltaAcceptance[currentInterEventTime];
                                }
                                if (t_upperBound_DeltaAcceptance){
                                    upperBound_DeltaAcceptance[upperBound] += currentDeltaAcceptance;
                                    ++sampledUpperBound_DeltaAcceptance[upperBound];
                                }
                                if (t_deltaUpperBound_DeltaAcceptance){
                                    deltaUpperBound_DeltaAcceptance[model.getDeltaMaximumClusterSize()] += currentDeltaAcceptance;
                                    ++sampledDeltaUpperBound_DeltaAcceptance[model.getDeltaMaximumClusterSize()];
                                }
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
                if (ensembleSize < 5 && t_dynamics){
                    dynamics.emplace_back(std::vector<int>{time, trialTime, periodTime, periodTrialTime, model.getMaximumClusterSize(), upperBound});
                }//* End of one step
            }//* End of network growing (one ensemble)
        } //* End of every ensembles
    } //* End of function mBFW::generate::run

    //! Order Parameter
    void save_orderParameter(){
        using namespace linearAlgebra;
        const std::string fullFileName = directory + "orderParameter/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
        orderParameter /= ensembleSize;
        writeCSV(fullFileName, orderParameter);
    }

    //! Mean Cluster Size
    void save_meanClusterSize(){
        using namespace linearAlgebra;
        const std::string fullFileName = directory + "meanClusterSize/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
        meanClusterSize /= ensembleSize;
        writeCSV(fullFileName, meanClusterSize);
    }

    //! Second Giant
    void save_secondGiant(){
        using namespace linearAlgebra;
        const std::string fullFileName = directory + "secondGiant/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);
        secondGiant /= ensembleSize;
        writeCSV(fullFileName, secondGiant);
    }

    //! Inter Event Time
    void save_interEventTime(){
        for (const auto& state : states){
            const std::string fullFileName = directory + "interEventTime/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

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
            const std::string fullFileName = directory + "deltaAcceptance/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

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
            const std::string fullFileName = directory + "orderParameterDistribution/" + filename_time(networkSize, acceptanceThreshold, ensembleSize, t, coreNum);

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
            const std::string fullFileName = directory + "clusterSizeDistribution/" + filename_orderParameter(networkSize, acceptanceThreshold, ensembleSize, op, coreNum);

            //* save only useful data
            std::map<int, double> trimmed;
            const double tot = std::accumulate(clusterSizeDistribution[op].begin(), clusterSizeDistribution[op].end(), 0.0);
            for (int cs=0; cs<networkSize; ++cs){
                if (clusterSizeDistribution[op][cs]){
                    trimmed[cs] = clusterSizeDistribution[op][cs]/tot;
                }
            }
            writeCSV(fullFileName, trimmed);
        }
    }

    //! Age Distribution
    void save_ageDistribution(){
        for (const auto& state : states){
            const std::string fullFileName = directory + "ageDistribution/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

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
            const std::string fullFileName = directory + "interEventTimeDistribution/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

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
            const std::string fullFileName = directory + "deltaUpperBoundDistribution/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

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
            const std::string fullFileName = directory + "deltaAcceptanceDistribution/" + state + "/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

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
        const std::string fullFileName = directory + "interEventTime_DeltaAcceptance/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

        //* save only useful data
        std::map<int, double> trimmed;
        for (int iet=1; iet<networkSize; ++iet){
            if (sampledInterEventTime_DeltaAcceptance[iet] && interEventTime_DeltaAcceptance[iet]){
                trimmed[iet] = interEventTime_DeltaAcceptance[iet]/sampledInterEventTime_DeltaAcceptance[iet];
            }
        }
        writeCSV(fullFileName, trimmed);
    }

    //! Upper Bound vs Delta Acceptance
    void save_upperBound_DeltaAcceptance(){
        const std::string fullFileName = directory + "upperBound_DeltaAcceptance/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

        //* save only useful data
        std::map<int, double> trimmed;
        for (int k=2; k<networkSize; ++k){
            if (sampledUpperBound_DeltaAcceptance[k] && upperBound_DeltaAcceptance[k]){
                trimmed[k] = upperBound_DeltaAcceptance[k]/sampledUpperBound_DeltaAcceptance[k];
            }
        }
        writeCSV(fullFileName, trimmed);
    }

    //! Delta Upper Bound vs Delta Acceptance
    void save_deltaUpperBound_DeltaAcceptance(){
        const std::string fullFileName = directory + "deltaUpperBound_DeltaAcceptance/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, coreNum);

        //* save only useful data
        std::map<int, double> trimmed;
        for (int deltaK=1; deltaK<networkSize; ++deltaK){
            if (sampledDeltaUpperBound_DeltaAcceptance[deltaK] && deltaUpperBound_DeltaAcceptance[deltaK]){
                trimmed[deltaK] = deltaUpperBound_DeltaAcceptance[deltaK]/sampledDeltaUpperBound_DeltaAcceptance[deltaK];
            }
        }
        writeCSV(fullFileName, trimmed);
    }

    //! Dynamics
    void save_dynamics(const int& t_randomEngineSeed){
        const std::string fullFileName = directory + "dynamics/" + defaultFileName(networkSize, acceptanceThreshold, ensembleSize, -1, t_randomEngineSeed);
        writeCSV(fullFileName, dynamics);
    }
}//* End of namespace mBFW::generate