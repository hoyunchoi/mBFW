#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <map>

#include "/pds/pds172/hoyun1009/library/linearAlgebra.hpp"
#include "/pds/pds172/hoyun1009/library/randomNumber.hpp"
#include "/pds/pds172/hoyun1009/library/Networks.hpp"

// #include "Network.hpp"
#include "save.hpp"
#include "parameters.hpp"
#include "fileName.hpp"

namespace mBFW{
    //! Declaration of input variables
    int networkSize;
    int ensembleSize;
    int coreNum;
    double acceptanceThreshold;
    double precision;
    double degenerated;

    //! Declaration of pre-defined variables at "parameters.hpp"
    double m_c, m_a;
    std::vector<double> time_orderParameterDistribution, orderParameter_clusterSizeDistribution;
    double deltaAcceptanceBinDelta;
    std::vector<double> deltaAcceptanceExponent;
    std::vector<double> minExponent, maxExponent, logBinnedDeltaAcceptance;
    int deltaAcceptanceBinNum;
    const std::vector<std::string> states = {"before", "during"};

    //! Declaration of output variables
    //* X[time] = value of observable X at time
    std::vector<double> orderParameter;
    std::vector<double> secondGiant;
    std::vector<double> meanClusterSize;
    std::map<std::string, std::vector<double>> interEventTime;
    std::map<std::string, std::vector<int>> sampledInterEventTime;
    std::vector<double> deltaAcceptance;
    std::vector<int> sampledDeltaAcceptance;

    //* orderParameterDistribution[time] : distribution of rounded order parameter at (rounded time)=(time)
    std::map<double, std::vector<double>> orderParameterDistribution;

    //* clusterSizeDistribution[orderParameter] : cluster size distribution at (rounded order parameter)=(orderParameter)
    //* sampledClusterSizeDistribution[orderParameter] : number of samples for each rounded order parameter
    std::map<double, std::vector<double>> clusterSizeDistribution;
    std::map<double, int> sampledClusterSizeDistribution;

    //* interEventTimeDistribution["before"] : inter event time distribution before jump (order parameter < m_a)
    //* interEventTimeDistribution["during"] : inter event time distribution during jump (m_a < order parameter < m_c)
    std::map<std::string, std::vector<double>> interEventTimeDistribution;

    //* detaMDistribution["before"] : jump size of maximum cluster size distribution before jump
    //* detaMDistribution["during"] : jump size of maximum cluster size distribution during jump
    std::map<std::string, std::vector<double>> deltaMDistribution;

    //* ageDistribution["before"] : age distribution before jump (order parameter < m_a)
    //* ageDistribution["during"] : age distribution during jump (m_a < order parameter < m_c)
    std::map<std::string, std::vector<double>> ageDistribution;

    //* deltaAcceptanceDistribution["before"] : max acceptance distribution befeore jump (order paramter < m_a)
    //* deltaAcceptanceDistribution["during"] : max acceptance distribution during jump (m_a < order paramter < m_c)
    std::map<std::string, std::vector<double>> deltaAcceptanceDistribution;

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

    void setParameters(const int& t_networkSize, const int& t_ensembleSize, const double& t_acceptanceThreshold, const double t_precision, const int& t_coreNum){
        //! Input variables
        networkSize = t_networkSize;
        ensembleSize = t_ensembleSize;
        coreNum = t_coreNum;
        acceptanceThreshold = t_acceptanceThreshold;
        t_networkSize < t_precision ? precision = t_networkSize : precision = t_precision;

        //! Pre-defined variables
        std::tie(time_orderParameterDistribution, orderParameter_clusterSizeDistribution, m_c) = getParameters(networkSize, acceptanceThreshold);
        m_a = 0.05;
        degenerated=t_networkSize/t_precision;

        deltaAcceptanceBinDelta = 0.01;
        deltaAcceptanceExponent = linearAlgebra::arange(-8,0,deltaAcceptanceBinDelta);
        deltaAcceptanceBinNum = deltaAcceptanceExponent.size()-1;
        minExponent.resize(deltaAcceptanceBinNum);
        maxExponent.resize(deltaAcceptanceBinNum);
        logBinnedDeltaAcceptance.resize(deltaAcceptanceBinNum);
        for (int i=0; i<deltaAcceptanceBinNum; ++i){
            minExponent[i] = pow(10, deltaAcceptanceExponent[i]);
            maxExponent[i] = pow(10, deltaAcceptanceExponent[i+1]);
            logBinnedDeltaAcceptance[i] = sqrt(minExponent[i] * maxExponent[i]);
        }

        //! Output variables (Observables)
        //* time-X
        orderParameter.resize(t_networkSize);
        secondGiant.resize(t_networkSize);
        meanClusterSize.resize(t_networkSize);       
        deltaAcceptance.resize(t_networkSize);
        sampledDeltaAcceptance.resize(t_networkSize);
        for (auto state : states){
            interEventTime[state].resize(t_networkSize);
            sampledInterEventTime[state].resize(t_networkSize);
        }

        //* Distributions
        for (const double& t : time_orderParameterDistribution){
            orderParameterDistribution[t].resize(t_networkSize);
        }
        for (const double& m : orderParameter_clusterSizeDistribution){
            clusterSizeDistribution[m].resize(t_networkSize);
            sampledClusterSizeDistribution[m] = 0;
        }
        for (auto state : states){
            interEventTimeDistribution[state].resize(t_networkSize);
            deltaMDistribution[state].resize(t_networkSize);
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
    } //* End of setParameters

    void run(){
        for (int ensemble=0; ensemble<ensembleSize; ++ensemble){
            //* Default values for one ensemble
            NZ_Network model(networkSize);
            Node node1, node2, root1, root2;
            int size1, size2;
            int time, trialTime, upperBound=2;
            int periodTime, periodTrialTime, peakTime, peakTrialTime, updatedTime;
            double maxAcceptance, maxPeriodAcceptance;
            bool findNewNodes = true;

            //* Random integer generator for the algorithm
            SNU::CNRC::RandomIntGenerator randomNode(0, networkSize-1);

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
                        node1 = randomNode();
                        root1 = model.getRoot(node1);
                        node2 = randomNode();
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

                        //! Delta Acceptance
                        deltaAcceptance[time] += currentDeltaAcceptance;
                        ++sampledDeltaAcceptance[time];

                        //! Cluster Size Distribution w.r.t. order parameter
                        const double roundedOrderParameter=round(exactOrderParameter*precision)/precision;
                        auto it = std::find(orderParameter_clusterSizeDistribution.begin(), orderParameter_clusterSizeDistribution.end(), roundedOrderParameter);
                        if (it != orderParameter_clusterSizeDistribution.end()){
                            ++sampledClusterSizeDistribution[*it];
                            auto sortedCluster = model.getSortedCluster();
                            for (auto it2 = sortedCluster.begin(); it2!= sortedCluster.end(); ++it2){
                                clusterSizeDistribution[*it][it2->first] += it2->second/(double)networkSize;
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

                            //! Delta M Distribution
                            ++deltaMDistribution[currentState][model.getDeltaMaximumClusterSize()];

                            //! Delta Acceptance Distribution
                            for (int i=0; i<deltaAcceptanceBinNum; ++i){
                                if (minExponent[i] <= currentDeltaAcceptance && currentDeltaAcceptance < maxExponent[i]){
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

        //* normalization and save data
        {
            using namespace linearAlgebra;
            const std::string filename=fileName(networkSize, acceptanceThreshold, ensembleSize);

            //! Order Parameter
            orderParameter /= ensembleSize;
            // saveOrderParameter(filename, coreNum, orderParameter);

            //! Mean Cluster Size
            meanClusterSize /= ensembleSize;
            // saveMeanClusterSize(filename, coreNum, meanClusterSize);

            //! Second Giant
            secondGiant /= ensembleSize;
            // saveSecondGiant(filename, coreNum, secondGiant);

            //! Inter Event Time
            for (const auto& state : states){
                for (int t=0; t<networkSize; ++t){
                    if (sampledInterEventTime[state][t]){
                        interEventTime[state][t] /= sampledInterEventTime[state][t];
                    }
                }
            }
            // saveInterEventTime(filename, coreNum, interEventTime);

            //! Delta Acceptance
            for (int t=0; t<networkSize; ++t){
                if (sampledDeltaAcceptance[t]){
                    deltaAcceptance[t] /= sampledDeltaAcceptance[t];
                }
            }
            // saveDeltaAcceptance(filename, coreNum, deltaAcceptance);

            //! Order Parameter Distribution
            for (const auto& t : time_orderParameterDistribution){
                const double tot = std::accumulate(orderParameterDistribution[t].begin(), orderParameterDistribution[t].end(), 0.0);
                orderParameterDistribution[t] /= tot;
            }
            // saveOrderParameterDistribution(filename, coreNum, timeOfOrderParameterDistribution, orderParameterDistribution);

            //! Cluster Size Distribution
            for (const auto& op : orderParameter_clusterSizeDistribution){
                if (sampledClusterSizeDistribution[op]){
                    clusterSizeDistribution[op] /= sampledClusterSizeDistribution[op];
                }
            }
            // saveClusterSizeDistribution(filename, coreNum, orderParameterOfClusterSizeDistribution, clusterSizeDistribution1, sampledClusterSizeDistribution, 1);

            //! Age Distribution
            for (const auto& state : states){
                const double tot = std::accumulate(ageDistribution[state].begin(), ageDistribution[state].end(), 0.0);
                ageDistribution[state] /= tot;
            }
            // saveAgeDistribution(filename, coreNum, ageDistribution);

            //! Inter Event Time Distribution
            for (const auto& state : states){
                const double tot = std::accumulate(interEventTimeDistribution[state].begin(), interEventTimeDistribution[state].end(), 0.0);
                interEventTimeDistribution[state] /= tot;
            }
            // saveInterEventTimeDistribution(filename, coreNum, interEventTimeDistribution);
            
            //! Delta M Distribution
            for (const auto& state : states){
                const double tot = std::accumulate(deltaMDistribution[state].begin(), deltaMDistribution[state].end(), 0.0);
                deltaMDistribution[state] /= tot;
            }
            // saveDeltaMDistribution(filename, coreNum, deltaMDistribution);

            //! Delta Acceptance Distribution
            for (const auto& state : states){
                const double tot = std::accumulate(deltaAcceptanceDistribution[state].begin(), deltaAcceptanceDistribution[state].end(), 0.0);
                deltaAcceptanceDistribution[state] /= tot;
            }
            // saveAcceptanceDistribution(filename, coreNum, logBinnedDeltaAcceptance, acceptanceDistribution);

            //! Inter Event Time vs Delta Acceptance
            for (int i=0; i<networkSize; ++i){
                if (sampledInterEventTime_DeltaAcceptance[i]){
                    interEventTime_DeltaAcceptance[i] /= sampledInterEventTime_DeltaAcceptance[i];
                }
            }
            // saveInterEventTime_DeltaAcceptance(filename, core, interEventTime_DeltaAcceptance);

            //! Upper Bound vs Delta Acceptance
            for (int i=0; i<networkSize; ++i){
                if (sampledUpperBound_DeltaAcceptance[i]){
                    upperBound_DeltaAcceptance[i] /= sampledUpperBound_DeltaAcceptance[i];
                }
            }
            // saveUpperBound_DeltaAcceptance(filename, core, upperBound_DeltaAcceptance);

            //! Delta Upper Bound vs Delta Acceptance
            for (int i=0; i<networkSize; ++i){
                if (sampledDeltaUpperBound_DeltaAcceptance[i]){
                    deltaUpperBound_DeltaAcceptance[i] /= sampledDeltaUpperBound_DeltaAcceptance[i];
                }
            }
            // saveDeltaUpperBound_DeltaAcceptance(filename, core, deltaUpperBound_DeltaAcceptance);

            //! Dynamics
            // saveDynamics(filename, core, dynamics);
        }//* End of saving
    } //* End of function run
} //* End of namespace mBFW