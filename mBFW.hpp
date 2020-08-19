#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <map>

#include "/pds/pds172/hoyun1009/library/linearAlgebra.hpp"
#include "/pds/pds172/hoyun1009/library/randomNumber.hpp"

#include "Network.hpp"
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
    double maxAcceptanceBinDelta;
    std::vector<double> maxAcceptanceExponent;
    std::vector<double> minExponent, maxExponent, logBinnedMaxAcceptance;
    int maxAcceptanceBinNum

    //! Declaration of output variables
    //* X[time] = value of observable X at time
    std::vector<double> orderParameter;
    std::vector<double> secondGiant;
    std::vector<double> meanClusterSize;
    std::map<string, std::vector<double>> interEventTime;
    std::map<string, std::vector<int>> sampledInterEventTime;
    std::vector<double> maxAcceptance;
    std::vector<int> sampledMaxAcceptance;

    //* orderParameterDistribution[time] : distribution of rounded order parameter at (rounded time)=(time)
    std::map<double, std::vector<double>> orderParameterDistributon;

    //* clusterSizeDistribution[orderParameter] : cluster size distribution at (rounded order parameter)=(orderParameter)
    //* sampledClusterSizeDistribution[orderParameter] : number of samples for each rounded order parameter
    std::map<double, std::vector<double>> clusterSizeDistribution;
    std::map<double, int> sampledClusterSizeDistribution;

    //* interEventTimeDistribution['before'] : inter event time distribution before jump (order parameter < m_a)
    //* interEventTimeDistribution['during'] : inter event time distribution during jump (m_a < order parameter < m_c)
    std::map<std::string, std::vector<double>> interEventTimeDistribution;

    //* detaMDistribution['before'] : jump size of maximum cluster size distribution before jump
    //* detaMDistribution['during'] : jump size of maximum cluster size distribution during jump
    std::map<std::string, std::vector<double>> deltaMDistribution;

    //* ageDistribution['before'] : age distribution before jump (order parameter < m_a)
    //* ageDistribution['during'] : age distribution during jump (m_a < order parameter < m_c)
    std::map<std::string, std::vector<double>> ageDistribution;

    //* maxAcceptanceDistribution['before'] : max acceptance distribution befeore jump (order paramter < m_a)
    //* maxAcceptanceDistribution['during'] : max acceptance distribution during jump (m_a < order paramter < m_c)
    std::map<std::string, std::vector<double>> maxAcceptanceDistribution;

    //* X_maxAcceptance[x] : average max acceptance at X=x before jump
    //* sampledX_maxAcceptance[x] : number of samples for each X=x
    std::vector<double> interEventTime_maxAcceptance;
    std::vector<int> sampledInterEventTime_maxAcceptance;
    std::vector<double> upperBound_maxAcceptance;
    std::vector<int> sampledUpperBound_maxAcceptance;
    std::vector<double> deltaUpperBound_maxAcceptance;
    std::vector<int> sampledDeltaUpperbound_maxAcceptance;

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
        
        maxAcceptanceBinDelta = 0.01;
        maxAcceptanceExponent = linearAlgebra::arange(-8,0,maxAcceptanceBinDelta);
        maxAcceptanceBinNum = maxAcceptanceExponent.size()-1;
        minExponent.resize(maxAcceptanceBinNum);
        maxExponent.resize(maxAcceptanceBinNum);
        logBinnedMaxAcceptance.resize(maxAcceptanceBinNum);
        for (int i=0; i<maxAcceptanceBinNum; ++i){
            minExponent[i] = pow(10, maxAcceptanceExponent[i]);
            maxExponent[i] = pow(10, maxAcceptanceExponent[i+1]);
            logBinnedMaxAcceptance[i] = sqrt(minExponent[i] * maxExponent[i]);
        }

        //! Output variables (Observables)
        //* time-X
        orderParameter.resize(t_networkSize);
        secondGiant.resize(t_networkSize);
        meanClusterSize.resize(t_networkSize);
        interEventTime['before'].resize(t_networkSize);
        interEventTime['during'].resize(t_networkSize);
        sampledInterEventTime['before'].resize(t_networkSize);
        sampledInterEventTime['during'].resize(t_networkSize);
        maxAcceptance.resize(t_networkSize);
        sampledMaxAcceptance.resize(t_networkSize);

        //* Distributions
        for (const double& t : time_orderParameterDistribution){
            orderParameterDistribution[t].resize(t_netorkSize);
        }
        for (const double& m : orderParameter_clusterSizeDistribution){
            clusterSizeDistribution[m].resize(t_networkSize);
            sampledClusterSizeDistribution[m] = 0;
        }
        interEventTimeDistribution['before'].resize(t_networkSize);
        interEventTimeDistribution['during'].resize(t_networkSize);
        deltaMDistribution['before'].resize(t_networkSize);
        deltaMDistribution['during'].resize(t_networkSize);
        ageDistribution['before'].resize(t_networkSize);
        ageDistribution['during'].resize(t_networkSize);
        maxAcceptanceDistribution['before'].resize(maxAcceptanceBinNum);
        maxAcceptanceDistribution['during'].resize(maxAcceptanceBinNum);

        //* X_maxAcceptance
        interEventTime_maxAcceptance.resize(t_networkSize);
        sampledInterEventTime_maxAcceptance.resize(t_networkSize);
        upperBound_MaxAcceptance.resize(t_networkSize);
        sampledUpperBound_MaxAcceptance.resize(t_networkSize);
        deltaUpperBound_MaxAcceptance.resize(t_networkSize);
        sampledDeltaUpperbound_MaxAcceptance.resize(t_networkSize);
    }

    void run(){
        for (int ensemble=0; ensemble<ensembleSize; ++ensemble){
            //* Default values for one ensemble
            Network model(networkSize);
            int time, trialTime, upperBound=2;
            int periodTime, periodTrialTime, peakTime, peakTrialTime, updatedTime;
            double maxAcceptance, maxPeriodAcceptance;
            bool fineNewNode = true;

            //* Random integer generator for the algorithm
            SNU::CNRC::RandomIntGenerator randomNode(0,t_networkSize-1);

            //* Do rBFW algorithm until all clusters merge to one
            while (model.getMaximumCluster()<t_networkSize){
                //* Find new nodes
                if (findNewNodes){
                    //* Randomly choose new nodes
                    do {
                        node1 = randomNode();
                        root1 = model.getRoot(node1);
                        node2 = randomNode();
                        root2 = model.getRoot(node2);
                    }while (root1 == root2);

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
                    const double exactOrderParameter=model.getMaximumCluster()/(double)t_networkSize;

                    //* max acceptance
                    if ((double)time/trialTime > maxAcceptance){
                        peakTime = time;
                        peakTrialTime = trialTime;
                        maxAcceptance = (double)time/trialTime;
                    }

                    //! order parameter
                    orderParameter[time] += exactOrderParameter;

                    //! mean cluster size
                    meanClusterSize[time] += model.meanCluster();

                    //! second giant
                    secondGiant[time] += model.getSecondGiant()/(double)t_networkSize;

                    //! Age Distribution
                    if (exactOrderParameter < m_a){
                        for (auto& changedAge : model.getChangedAge()){
                            ageDistribution['before'][changedAge[0]] += changedAge[1];
                        }
                    }
                    else if (exactOrderParameter < m_c){
                        for (auto& changedAge : model.getChangedAge()){
                            ageDistribution['during'][changedAge[0]] += changedAge[1]
                        }
                    }

                    //! order parameter distribution
                    const double roundedTime=round(time/degenerated)/t_precision;
                    auto it = std::find(time_orderParameterDistribution.begin(), time_orderParameterDistribution.end(), roundedTime);
                    if (it != time_orderParameterDistribution.end()){
                        ++orderParameterDistributon[*it][model.getMaximumCluster()];
                    }


                    //* End of k-period
                    if (model.getMaximumClusterUpdated()){
                        //! Cluster Size Distribution w.r.t. order parameter
                        const double roundedOrderParameter=round(exactOrderParameter*t_precision)/t_precision;
                        auto it = std::find(orderParameter_clusterSizeDistribution.begin(), orderParameter_clusterSizeDistribution.end(), roundedOrderParameter);
                        if (it != orderParameter_clusterSizeDistribution.end()){
                            ++sampledClusterSizeDistribution[*it];
                            auto sortedCluster = model.getSortedCluster();
                            for (auto it2 = sortedCluster.begin(); it2!= sortedCluster.end(); ++it2){
                                clusterSizeDistribution[*it][it2->first] += it2->second/(double)networkSize;
                            }    
                        }

                        //* Before jump
                        if (exactOrderParameter < m_a){
                            int currentInterEventTime = time-updatedTime;

                            //! Inter Event Time Distribution
                            ++interEventTimeDistribution['before'][currentInterEventTime];

                            //! inter event time
                            interEventTime['before'][time] += currentInterEventTime;
                            ++sampledInterEventTime['before'][time];


                        }//* End of before jump

                        //* During jump
                        else if(exactOrderParameter<m_c){
                            int currentInterEventTime = time-updatedTime;

                            //! Inter Event Time Distribution
                            ++interEventTimeDistribution['during'][currentInterEventTime];

                            //! inter event time
                            interEventTime['during'][time] += currentInterEventTime;
                            ++sampledInterEventTime['during'][time];
                         
                        }//* End of during jump

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
                } //* End of one step
            
            
            }//* End of network growing




        }//* End of one ensemble

    }

}


void mBFW(const int &t_networkSize, const double &t_g, const int &t_ensembleSize, const int &t_coreNum, double t_precision){
    if (t_networkSize<t_precision){
        t_precision=t_networkSize;
    }
    const double degenerated=t_networkSize/t_precision;
    double m_c;
    std::vector<double> timeOfOrderParameterDistribution, orderParameterOfClusterSizeDistribution;
    std::tie(timeOfOrderParameterDistribution, orderParameterOfClusterSizeDistribution, m_c) = getParameters(t_networkSize, t_g);

    //* Cluster Size distribution at Order Parameter
    //* clusterSizeDistribution[i] : cluster size distribution at (rounded order parameter)=orderParameterOfClusterSizeDistribution[i]
    std::map<std::vector<double>> clusterSizeDistribution;
    for (const double orderParameter : orderParameterOfClusterSizeDistirbution){
        clusterSizeDistribution[orderParameter].resize(t_networkSize);
    }

    std::vector<int> sampledClusterSizeDistribution(orderParameterOfClusterSizeDistribution.size());

    //* Order Parameter Distribution
    //* orderParameterDistribution[i] : distribution of rounded order parameter at (rounded time)=timeOforderparameterDistribution[i]
    std::vector<std::vector<double>> orderParameterDistribution(timeOfOrderParameterDistribution.size());
    for (auto &e : orderParameterDistribution){
        e.resize(t_networkSize);
    }

    //* order Parameter and second Giant, mean cluster size at each time
    std::vector<double> orderParameter(t_networkSize);
    std::vector<double> secondGiant(t_networkSize);
    std::vector<double> meanClusterSize(t_networkSize);

    //* order parmeter after the jump
    std::vector<double> orderParameterAfter(t_networkSize);
    std::vector<int> sampledOrderParameterAfter(t_networkSize);

    //* Inter Event Time Distribution
    //* interEventTimeDistribution[0] : when m<m_a (before jump), interEventTimeDistribution[1] : during jump
    //* interEventTimeDistribution[0][i] : number of i inter event time
    std::vector<std::vector<double>> interEventTimeDistribution(2);
    interEventTimeDistribution[0].resize(t_networkSize);
    interEventTimeDistribution[1].resize(t_networkSize);
    int sampledPeriod=0;

    //* Inter Event Time vs acceptance
    std::vector<double> interEventTime_Acceptance(t_networkSize);
    std::vector<int> sampledInterEventTime_Acceptance(t_networkSize);

    //* Inter Event Time
    //* interEventTime[i] : average inter event time at time i
    std::vector<double> interEventTime, sampledInterEventTime;
    interEventTime.resize(t_networkSize);
    sampledInterEventTime.resize(t_networkSize);

    //* Delta m distribution
    //* deltaMDistribution[0] : when m<m_a (before jump), deltaMDistribution[1] : during jump
    std::vector<std::vector<double>>deltaMDistribution(2);
    deltaMDistribution[0].resize(t_networkSize);
    deltaMDistribution[1].resize(t_networkSize);

    //* Acceptance Distribution
    // int binNum=80;
    // std::vector<double> exponent;
    // {
    //     using namespace linearAlgebra;
    //     exponent=linspace(-8,0,binNum+1);
    // }
    // std::vector<double> minX(binNum);
    // std::vector<double> maxX(binNum);
    // std::vector<std::vector<double>> acceptanceDistribution(binNum);
    // for(int i=0; i<binNum; ++i){
    //     minX[i]=pow(10,exponent[i]);
    //     maxX[i]=pow(10,exponent[i+1]);
    //     acceptanceDistribution[i]=std::vector<double>{sqrt(minX[i]*maxX[i]),0};
    // }

    //* period acceptance - upper bound ratio
    // std::vector<std::vector<double>> periodAcceptance_UpperBoundRatio=acceptanceDistribution;
    // std::vector<int> sampledPeriodAcceptance_UpperBoundRatio(binNum);

    //* period acceptance - delta k
    // std::vector<std::vector<double>> periodAcceptance_DeltaK=acceptanceDistribution;
    // std::vector<int> sampledPeriodAcceptance_DeltaK(binNum);

    //* period acceptance area distribution
    // const int binNum=180;
    // std::vector<double> exponent;
    // std::vector<std::vector<double>> periodAcceptanceAreaDistribution(binNum);
    // {
    //     using namespace linearAlgebra;
    //     exponent=linspace(-8,10,binNum+1);
    // }
    // std::vector<double> minX(binNum);
    // std::vector<double> maxX(binNum);
    // for(int i=0; i<binNum; ++i){
    //     minX[i]=pow(10,exponent[i]);
    //     maxX[i]=pow(10,exponent[i+1]);
    //     periodAcceptanceAreaDistribution[i]=std::vector<double>{sqrt(minX[i]*maxX[i]),0};
    // }

    //* Delta Acceptance
    //* x_DeltaAcceptance[i] : delta acceptance at x
    //* sampledX_DeltaAcceptance[i] : sampled number at x
    std::vector<double> k_DeltaAcceptance(t_networkSize);
    std::vector<int> sampledK_DeltaAcceptance(t_networkSize);
    std::vector<double> deltaK_DeltaAcceptance(t_networkSize);
    std::vector<int> sampledDeltaK_DeltaAcceptance(t_networkSize);
    std::vector<double> time_DeltaAcceptance(t_networkSize);
    std::vector<int> sampledTime_DeltaAcceptance(t_networkSize);

    //* Age Distribution
    //* ageDistribution[0][i] : number of i-age at [0,t_a]
    //* ageDistribution[1][i] : number of i-age at [t_a,t_c]
    std::vector<std::vector<double>> ageDistribution(2);
    ageDistribution[0].resize(t_networkSize);
    ageDistribution[1].resize(t_networkSize);

    //* upper bound ratio
    std::vector<std::vector<int>> successiveK;

    //* time, trialTime, maximumClusterSize, upperBound
    std::vector<std::vector<int>> LUMK;

    //* time, trialTime, periodTime, periodTrialTime
    std::vector<std::vector<int>> LUPLPU;

    //* temp
    int op;
    //* Iterate the algorithm for t_ensembleSize times
    for (int ensemble=0; ensemble<t_ensembleSize; ++ensemble){
        //* Default values for one loop
        Network model(t_networkSize);
        int upperBound=2, time=0, trialTime=0, updatedTime=0, periodTime=0, periodTrialTime=0, maxTime=0, maxTrialTime=0; 
        int node1=0, node2=0, root1=0, root2=0, size1=0, size2=0;
        bool findNewNodes=true;
        double maxAcceptance=0,periodAcceptanceArea=0, maxPeriodAcceptance=0, area=0;

        //* Random integer generator for the algorithm
        SNU::CNRC::RandomIntGenerator randomNode(0,t_networkSize-1);

        //* Do rBFW algorithm until all clusters mere to one
        while (model.getMaximumCluster()<t_networkSize){
            //* Find new nodes
            if (findNewNodes){
                //* Randomly choose new nodes
                node1=randomNode();
                node2=randomNode();

                //* If sampled nodes are in same cluster, find new nodes
                while(model.getRoot(node1)==model.getRoot(node2)){
                    node1=randomNode();
                    node2=randomNode();
                }
                //* choose two clusters of each node
                root1=model.getRoot(node1);
                root2=model.getRoot(node2);
                size1=model.getClusterSize(root1);
                size2=model.getClusterSize(root2);
            }

            //* Merge two clusters. Time step
            if (size1+size2<=upperBound){
                model.merge(root1,root2);
                ++time;
                ++trialTime;
                ++periodTime;
                ++periodTrialTime;
                periodAcceptanceArea+=(double)periodTime/periodTrialTime;
                findNewNodes=true;
                const double exactOrderParameter=model.getMaximumCluster()/(double)t_networkSize;

                //* Max Acceptance
                if ((double)time/trialTime > maxAcceptance){
                    maxTime=time;
                    maxTrialTime=trialTime;
                    maxAcceptance=(double)time/trialTime;
                }

                //* Max period Acceptance
                if ((double)periodTime/periodTrialTime > maxPeriodAcceptance && (double)periodTime/periodTrialTime){
                    maxPeriodAcceptance=(double)periodTime/periodTrialTime;
                }

                //! single sample value
                LUPLPU.emplace_back(std::vector<int>{time, trialTime, periodTime, periodTrialTime, op});
                op = model.getMaximumCluster();

                //! Age Distribution
                // for (auto& changedAge : model.getChangedAge()){
                //     if (exactOrderParameter < m_a){
                //         ageDistribution[0][changedAge[0]]+=changedAge[1];
                //     }
                //     else if (exactOrderParameter < m_c){
                //         ageDistribution[1][changedAge[0]]+=changedAge[1];
                //     }
                // }

                //! Order Parameter
                // orderParameter[time]+=exactOrderParameter;

                //! Mean Cluster Size
                // meanClusterSize[time]+=model.meanCluster(time);

                //! Second Giant
                // secondGiant[time]+=model.getSecondGiant()/(double)t_networkSize;

                //! Order Parameter after jump
                // if (exactOrderParameter>0.86){
                //     orderParameterAfter[time]+=exactOrderParameter;
                //     sampledOrderParameterAfter[time]++;
                // }

                //* End of k-period
                if (model.getMaximumClusterUpdated()){
                    //! Cluster Size Distribution
                    // const double roundedOrderParameter=round(exactOrderParameter*t_precision)/t_precision;
                    // std::vector<double>::iterator it=std::find(orderParameterOfClusterSizeDistribution.begin(), orderParameterOfClusterSizeDistribution.end(), roundedOrderParameter);
                    // if (it!=orderParameterOfClusterSizeDistribution.end()){
                    //     auto temp=it-orderParameterOfClusterSizeDistribution.begin();
                    //     ++sampledClusterSizeDistribution[temp];
                    //     const std::map<int,int> sortedCluster=model.getSortedCluster();
                    //     for (auto it2=sortedCluster.begin(); it2!=sortedCluster.end(); ++it2){
                    //         clusterSizeDistribution1[temp][it2->first]+=it2->second/(double)t_networkSize;
                    //         clusterSizeDistribution2[temp][it2->first]+=it2->second/(double)t_networkSize; 
                    //     }
                    //     clusterSizeDistribution2[temp][model.getSecondGiant()]-=1/(double)t_networkSize;
                    // }

                    //! Order Parameter distribution
                    // const double roundedTime=round(time/degenerated)/t_precision;
                    // std::vector<double>::iterator it=std::find(timeOfOrderParameterDistribution.begin(), timeOfOrderParameterDistribution.end(), roundedTime);
                    // if (it!=timeOfOrderParameterDistribution.end()){
                    //     auto temp=it-timeOfOrderParameterDistribution.begin();
                    //     orderParameterDistribution[temp][model.getMaximumCluster()]++;
                    // }

                    //! Inter Event Time
                    // interEventTime[time]+=time-updatedTime;
                    // ++sampledInterEventTime[time];

                    //* Before jump
                    if (exactOrderParameter < m_a){
                        //! Inter Event Time Distribution
                        // ++interEventTimeDistribution[0][time-updatedTime];
                        // ++sampledPeriod;

                        //! Delta M Distribution
                        // ++deltaMDistribution[0][model.getDeltaMaximumCluster()];

                        //! Acceptance Distribution
                        // for (int i=0; i<binNum; ++i){
                        //     if (minX[i]<=maxAcceptance-t_g && maxAcceptance-t_g<maxX[i]){
                        //         acceptanceDistribution[i][1]++;
                        //         break;
                        //     }
                        // }

                        //! x_DeltaAcceptance 
                        // k_DeltaAcceptance[upperBound]+=maxAcceptance-t_g;
                        // sampledK_DeltaAcceptance[upperBound]++;
                        // deltaK_DeltaAcceptance[model.getDeltaMaximumCluster()]+=maxAcceptance-t_g;
                        // sampledDeltaK_DeltaAcceptance[model.getDeltaMaximumCluster()]++;
                        // time_DeltaAcceptance[maxTime]+=maxAcceptance-t_g;
                        // sampledTime_DeltaAcceptance[maxTime]++;

                        //! period acceptance-Upper Bound Ratio 
                        // for (int i=0; i<binNum; ++i){
                        //     if (minX[i]<=maxPeriodAcceptance-t_g && maxPeriodAcceptance-t_g<maxX[i]){
                        //         periodAcceptance_UpperBoundRatio[i][1]+=(double)upperBound/(upperBound-model.getDeltaMaximumCluster())-1;
                        //         ++sampledPeriodAcceptance_UpperBoundRatio[i];
                        //         break;
                        //     }
                        // }

                        //! period acceptance-delta K
                        // for (int i=0; i<binNum; ++i){
                        //     if (minX[i]<=maxPeriodAcceptance-t_g && maxPeriodAcceptance-t_g<maxX[i]){
                        //         periodAcceptance_DeltaK[i][1]+=model.getDeltaMaximumCluster();
                        //         ++sampledPeriodAcceptance_DeltaK[i];
                        //         break;
                        //     }
                        // }

                        //! Period acceptance area distribution
                        // double normalizedArea=periodAcceptanceArea/(maxPeriodAcceptance-t_g);
                        // for (int i=0; i<binNum; ++i){
                        //     // if (minX[i]<=periodAcceptanceArea && periodAcceptanceArea<maxX[i]){
                        //     if (minX[i] <= normalizedArea && normalizedArea < maxX[i]){
                        //         periodAcceptanceAreaDistribution[i][1]++;
                        //     }
                        // }

                        //! inter event time vs max acceptance
                        // interEventTime_Acceptance[time-updatedTime]+=(maxAcceptance-t_g);
                        // sampledInterEventTime_Acceptance[time-updatedTime]++;

                        //! Successive upper bound
                        // successiveK.emplace_back(std::vector<int>{upperBound-model.getDeltaMaximumCluster(),upperBound});
                    }//* End of before jump

                    //* During jump
                    else if(exactOrderParameter<m_c){
                        //! Inter Event Time distribution
                        // ++interEventTimeDistribution[1][time-updatedTime];

                        //! Delta M Distribution
                        // ++deltaMDistribution[1][model.getDeltaMaximumCluster()];

                        //! x_DeltaAcceptance
                        // k_DeltaAcceptance[upperBound]+=maxAcceptance-t_g;
                        // sampledK_DeltaAcceptance[upperBound]++;
                        // deltaK_DeltaAcceptance[model.getDeltaMaximumCluster()]+=maxAcceptance-t_g;
                        // sampledDeltaK_DeltaAcceptance[model.getDeltaMaximumCluster()]++;
                        // time_DeltaAcceptance[time]+=maxAcceptance-t_g;
                        // sampledTime_DeltaAcceptance[time]++;

                        //! period acceptance-Upper Bound Ratio
                        // for (int i=0; i<binNum; ++i){
                            // if (minX[i]<=periodAcceptanceArea-t_g && periodAcceptanceArea-t_g<maxX[i]){
                        //     if (minX[i]<=maxPeriodAcceptance-t_g && maxPeriodAcceptance-t_g<maxX[i]){
                        //         periodAcceptance_UpperBoundRatio[i][1]+=(double)upperBound/(upperBound-model.getDeltaMaximumCluster())-1;
                        //         ++sampledPeriodAcceptance_UpperBoundRatio[i];
                        //         break;
                        //     }
                        // }

                        //! period acceptance-delta K
                        // for (int i=0; i<binNum; ++i){
                        //     if (minX[i]<=maxPeriodAcceptance-t_g && maxPeriodAcceptance-t_g<maxX[i]){
                        //         periodAcceptance_DeltaK[i][1]+=model.getDeltaMaximumCluster();
                        //         ++sampledPeriodAcceptance_DeltaK[i];
                        //         break;
                        //     }
                        // }
                    }//* End of during jump

                    //* Initialize variable for new k-period
                    // LUPLPU.emplace_back(std::vector<int>{time, trialTime, periodTime, periodTrialTime, model.getMaximumCluster()});
                    updatedTime=time;
                    maxAcceptance=0;
                    periodTime=0;
                    periodTrialTime=0;
                    periodAcceptanceArea=0;
                    maxPeriodAcceptance=0;
                    area=0;
                }//* End of k-period
            }
            else if ((double)time/trialTime<t_g){
                ++upperBound;
                findNewNodes=false;
            }
            else{
                ++trialTime;
                ++periodTrialTime;
                findNewNodes=true;
            }

            //! single sample value
            // LUMK.emplace_back(std::vector<int>{time, trialTime, model.getMaximumCluster(), upperBound});

        //* End of one step
        }//* End of network growing
    }

    //* Write CSV files
    {
        using namespace linearAlgebra;
        const std::string filename=fileName(t_networkSize, t_g, t_ensembleSize);

        //! Order Parameter
        // orderParameter/=t_ensembleSize;
        // saveOrderParameter(filename, t_coreNum, orderParameter);

        //! Mean Cluster Size
        // meanClusterSize/=t_ensembleSize;
        // saveMeanClusterSize(filename, t_coreNum, meanClusterSize);

        //! Second Giant
        // secondGiant/=t_ensembleSize;
        // saveSecondGiant(filename, t_coreNum, secondGiant);

        //! Order Parameter after jump (NOT NORMALIZED)
        // saveOrderParameterAfter(filename,t_networkSize, t_ensembleSize, t_coreNum, orderParameterAfter, sampledOrderParameterAfter);

        //! Cluster Size Distribution
        // for (int i=0; i<orderParameterOfClusterSizeDistribution.size(); ++i){
        //     if (sampledClusterSizeDistribution[i]!=0){
        //         clusterSizeDistribution1[i]/=sampledClusterSizeDistribution[i];
                // clusterSizeDistribution2[i]/=sampledClusterSizeDistribution[i];
        //     }
        // }
        // saveClusterSizeDistribution(filename, t_coreNum, orderParameterOfClusterSizeDistribution, clusterSizeDistribution1, sampledClusterSizeDistribution, 1);
        // saveClusterSizeDistribution(filename, t_coreNum, orderParameterOfClusterSizeDistribution, clusterSizeDistribution2, sampledClusterSizeDistribution, 2);

        //! Inter Event Time Distribution
        // for (auto &iet : interEventTimeDistribution){
        //     iet/=sampledPeriod;
        // }
        // saveInterEventTimeDistribution(filename, t_coreNum, interEventTimeDistribution);

        //! Inter Event Time vs Acceptance
        // for (int i=0; i<sampledInterEventTime_Acceptance.size(); ++i){
        //     if(sampledInterEventTime_Acceptance[i]){
        //         interEventTime_Acceptance[i]/=sampledInterEventTime_Acceptance[i];
        //     }
        // }
        // saveInterEventTimeAcceptance(filename, t_coreNum, interEventTime_Acceptance);

        //! Order Parameter Distribution
        // for (int i=0; i<timeOfOrderParameterDistribution.size(); ++i){
        //     orderParameterDistribution[i]/=degenerated;
        // }
        // saveOrderParameterDistribution(filename, t_coreNum, timeOfOrderParameterDistribution, orderParameterDistribution);

        //! single sample value
        // saveLUMK(LUMK);
        saveLUPLPU(filename, LUPLPU);

        //! Delta M distribution
        // for (auto &deltaM : deltaMDistribution){
        //     deltaM/=t_ensembleSize;
        // }
        // saveDeltaMDistribution(filename, t_coreNum, deltaMDistribution);

        //! Acceptance Distribution
        // double sum=0;
        // for (auto &acceptance : acceptanceDistribution){
        //     sum+=acceptance[1];
        // }
        // for (auto &acceptance : acceptanceDistribution){
        //     acceptance[1]/=sum;
        // }
        // saveAcceptanceDistribution(filename, t_coreNum, acceptanceDistribution);

        //! Age Distribution
        // for (std::vector<double> &age : ageDistribution){
        //     age/=t_ensembleSize;
        // }
        // saveAgeDistribution(filename, t_coreNum, ageDistribution);

        //! x_Delta Acceptance
        // for (int upperBound=0; upperBound<t_networkSize; ++upperBound){
        //     if (sampledK_DeltaAcceptance[upperBound]){
        //         k_DeltaAcceptance[upperBound]/=sampledK_DeltaAcceptance[upperBound];
        //     }
        // }
        // saveK_DeltaAcceptance(filename, t_coreNum, k_DeltaAcceptance);
        // for (int deltaK=0; deltaK<t_networkSize; ++deltaK){
        //     if (sampledDeltaK_DeltaAcceptance[deltaK]){
        //         deltaK_DeltaAcceptance[deltaK]/=sampledDeltaK_DeltaAcceptance[deltaK];
        //     }
        // }
        // saveDeltaK_DeltaAcceptance(filename, t_coreNum, deltaK_DeltaAcceptance);
        // for (int time=0; time<t_networkSize; ++time){
        //     if (sampledTime_DeltaAcceptance[time]){
        //         time_DeltaAcceptance[time]/=sampledTime_DeltaAcceptance[time];
        //     }
        // }
        // saveTime_DeltaAcceptance(filename, t_coreNum, time_DeltaAcceptance);

        //! Successive k
        // saveSuccessiveK(successiveK);

        //! period acceptance vs upper bound ratio
        // for (int i=0; i<binNum; ++i){
        //     if (sampledPeriodAcceptance_UpperBoundRatio[i]){
        //         periodAcceptance_UpperBoundRatio[i][1]/=sampledPeriodAcceptance_UpperBoundRatio[i];
        //     }
        // }
        // savePeriodAcceptance_UpperBoundRatio(filename, t_coreNum, periodAcceptance_UpperBoundRatio);

        //! period acceptance - delta k
        // for (int i=0; i<binNum; ++i){
        //     if (sampledPeriodAcceptance_DeltaK[i]){
        //         periodAcceptance_DeltaK[i][1]/=sampledPeriodAcceptance_DeltaK[i];
        //     }
        // }
        // savePeriodAcceptance_DeltaK(filename, t_coreNum, periodAcceptance_DeltaK);

        //! Period acceptance area distribution
        // double sum=0;
        // for (int i=0; i<binNum; ++i){
        //     sum+=periodAcceptanceAreaDistribution[i][1];
        // }
        // for (auto &e : periodAcceptanceAreaDistribution){
        //     e[1]/=sum;
        // }
        // savePeriodAcceptanceAreaDistribution(filename, t_coreNum, periodAcceptanceAreaDistribution);

        //! inter event time
        // for (int i=0; i<t_networkSize; ++i){
        //     if (sampledInterEventTime[i]!=0){
        //         interEventTime[i]/=sampledInterEventTime[i];
        //     }
        // }
        // saveInterEventTime(filename, t_coreNum, interEventTime);
    }
}