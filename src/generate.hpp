#pragma once

#include <map>
#include <random>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "pcg_random.hpp"

//*-------------------------
#include <algorithm>
#include <cmath>
#include <unistd.h>

#include "CSV.hpp"
#include "NZ_Network.hpp"
#include "common.hpp"
#include "linearAlgebra.hpp"
#include "parameters.hpp"
//*-------------------------

/*
    Observables
    obs_orderParameter[t] : Order parameter at time t
    obs_netOrderParameter[t] : Order Parameter greater than m_c
    obs_singleOrderParameter[ensemble][op] : ordre parameter for single ensemble
    obs_secondMaximum[t] : Second maximum cluster at time t
    obs_secondMoment[t] : Second moment (op^2) at time t
    obs_meanClusterSize[t] : Mean cluster size (of finite clusters) at time t
    obs_interEventTime[t] : Inter event time finished at time t, sample number

    obs_ageDist["state"][age] : Number of age accumulated at interval "state"
    obs_ageDist_time["state_time"][age] : Number of age accumulated at interval "state_time"
    obs_interEventTimeDist["state"][iet] : Number of inter event time accumulated through interval "state"
    obs_deltaUpperBoundDist["state"][dk] : Number of delta upper bound accumulated through interval "state"
    obs_clusterSizeDist[op][cs] : Number of clusters with size cs when order parameter exceeds op
    obs_orderParameterDist[t][op] : Number of ensembles with order parameter op at time t

    obs_interEventTime_orderParameter[iet] : Order parameter when inter event time is iet, sample number
    obs_orderParameter_interEventTime[m] : Inter Event time when maximum cluster size is m, sample number

    obs_dynamics[u] : trialTime, time, maximumClusterSize, upperBound
    obs_preiodDynamics[u] : preiodTrialTime, periodTime
*/

namespace mBFW {
struct Generate {
  protected:
    //* Member variables
    int m_networkSize;
    double m_acceptanceThreshold;
    unsigned m_maxTime;
    unsigned m_ensembleSize;
    int m_coreNum;
    int m_randomEngineSeed;
    pcg32 m_randomEngine;
    std::uniform_int_distribution<int> m_nodeDistribution;
    std::map<std::string, int> m_points;
    std::vector<std::string> m_state;
    std::vector<std::string> m_state_time;
    std::vector<std::string> m_subSuperState;
    std::set<int> m_clusterSizeDist_orderParameter;
    std::set<int> m_clusterSizeDist_time;
    std::set<unsigned> m_orderParameterDist_time;

    //* Observables
    std::vector<double> obs_orderParameter;
    std::vector<double> obs_secondMaximum;
    std::vector<double> obs_secondMoment;
    std::vector<double> obs_meanClusterSize;

    std::map<std::string, std::vector<unsigned>> obs_subSuperEnsemble;
    std::map<std::string, std::vector<double>> obs_netOrderParameter;
    std::map<std::string, std::vector<double>> obs_netSecondMoment;
    std::vector<std::pair<unsigned, unsigned>> obs_interEventTime;

    std::map<std::string, std::vector<unsigned long long>> obs_ageDist;
    std::map<std::string, std::vector<unsigned long long>> obs_ageDist_time;
    std::map<std::string, std::vector<unsigned>> obs_interEventTimeDist;
    std::map<std::string, std::vector<unsigned>> obs_interEventTimeDist_time;
    std::map<std::string, std::vector<unsigned>> obs_deltaUpperBoundDist;
    std::map<std::string, std::vector<unsigned>> obs_deltaUpperBoundDist_time;

    std::map<int, std::map<int, unsigned long long>> obs_clusterSizeDist;
    std::map<int, std::map<int, unsigned long long>> obs_clusterSizeDist_time;
    std::map<int, std::map<int, unsigned>> obs_orderParameterDist;

    std::vector<std::pair<double, unsigned>> obs_interEventTime_orderParameter;
    std::vector<std::pair<double, unsigned>> obs_orderParameter_interEventTime;

    std::vector<std::vector<double>> obs_singleOrderParameter;
    std::vector<std::vector<int>> obs_dynamics;
    std::vector<std::vector<int>> obs_periodDynamics;

  public:
    Generate() {}
    Generate(const int&, const double&, const unsigned&, const int&, const int& t_randomEngineSeed = -1);

    //* Member functions
    void run(const unsigned&);
    void save() const;

  protected:
    void m_singleRun();
    const double m_getAcceptanceThreshold(const int&) const;
};

Generate::Generate(const int& t_networkSize,
                   const double& t_acceptanceThreshold,
                   const unsigned& t_maxTime,
                   const int& t_coreNum,
                   const int& t_randomEngineSeed)
    : m_networkSize(t_networkSize),
      m_acceptanceThreshold(t_acceptanceThreshold),
      m_maxTime(t_maxTime),
      m_coreNum(t_coreNum),
      m_randomEngineSeed(t_randomEngineSeed) {
    //* Get default values from parameter
    mBFW::Parameter parameter(t_networkSize, t_acceptanceThreshold);
    m_points = parameter.get_points();
    m_clusterSizeDist_orderParameter = parameter.get_clusterSizeDist_orderParameter();
    m_clusterSizeDist_time = parameter.get_clusterSizeDist_time();
    m_orderParameterDist_time = parameter.get_orderParameterDist_time();

    //* Initialize random variables
    t_randomEngineSeed == -1 ? m_randomEngine.seed((std::random_device())()) : m_randomEngine.seed(t_randomEngineSeed);
    m_nodeDistribution.param(std::uniform_int_distribution<int>::param_type(0, t_networkSize - 1));

    //* Initialize state vectors using points
    m_state.resize(m_maxTime);
    m_state_time.resize(m_maxTime);
    m_subSuperState.resize(m_maxTime);
    {
        for (int m = 0; m < m_points.at("m_a1"); ++m) {
            m_state[m] = "0_A1";
        }
        for (int m = m_points.at("m_a1"); m < m_points.at("m_a2"); ++m) {
            m_state[m] = "A1_A2";
        }
        for (int m = m_points.at("m_a2"); m < m_points.at("m_b"); ++m) {
            m_state[m] = "A2_B";
        }
        for (int m = m_points.at("m_b"); m < m_points.at("m_c"); ++m) {
            m_state[m] = "B_C";
        }
        for (int m = m_points.at("m_c"); m < t_networkSize; ++m) {
            m_state[m] = "C_1";
        }
        for (int t = 0; t < m_points.at("t_a1"); ++t) {
            m_state_time[t] = "0_A1";
        }
        for (int t = m_points.at("t_a1"); t < m_points.at("t_a2"); ++t) {
            m_state_time[t] = "A1_A2";
        }
        for (int t = m_points.at("t_a2"); t < m_points.at("t_b"); ++t) {
            m_state_time[t] = "A2_B";
        }
        for (int t = m_points.at("t_b"); t < m_points.at("t_c"); ++t) {
            m_state_time[t] = "B_C";
        }
        for (int t = m_points.at("t_c"); t < t_networkSize; ++t) {
            m_state_time[t] = "C_1";
        }
        const int subSuper = (int)(m_networkSize * parameter.get_subSuper() * 7.0 / 8.0);
        for (int m = 0; m < subSuper; ++m) {
            m_subSuperState[m] = "sub";
        }
        for (int m = subSuper; m < t_networkSize; ++m) {
            m_subSuperState[m] = "super";
        }
    }

    //* Initialize observables
    obs_orderParameter.assign(m_maxTime, 0.0);
    obs_secondMaximum.assign(m_maxTime, 0.0);
    obs_meanClusterSize.assign(m_maxTime, 0.0);
    obs_secondMoment.assign(m_maxTime, 0.0);
    obs_interEventTime.assign(m_maxTime, std::pair<unsigned, unsigned>{0, 0});

    // for (const std::string& state : std::set<std::string>{"sub", "super"}){
    //     obs_subSuperEnsemble[state].assign(t_networkSize, 0);
    //     obs_netOrderParameter[state].assign(t_networkSize, 0.0);
    //     obs_netSecondMoment[state].assign(t_networkSize, 0.0);
    // }

    for (const std::string& state : mBFW::states) {
        obs_ageDist[state].assign(m_maxTime, 0);
        obs_interEventTimeDist[state].assign(m_maxTime, 0);
        obs_deltaUpperBoundDist[state].assign(m_maxTime, 0);
        // obs_ageDist_time[state].assign(t_networkSize, 0);
        // obs_interEventTimeDist_time[state].assign(t_networkSize, 0);
        // obs_deltaUpperBoundDist_time[state].assign(t_networkSize, 0);
    }

    for (const int& op : m_clusterSizeDist_orderParameter) {
        obs_clusterSizeDist[op] = std::map<int, unsigned long long>{};
    }
    // for (const int& t : m_clusterSizeDist_time) {
    //     obs_clusterSizeDist_time[t] = std::map<int, unsigned long long>{};
    // }

    for (const int& t : m_orderParameterDist_time) {
        obs_orderParameterDist[t] = std::map<int, unsigned>{};
    }

    // obs_interEventTime_orderParameter.assign(t_networkSize, std::pair<double, unsigned>{0.0, 0});
    // obs_orderParameter_interEventTime.assign(t_networkSize, std::pair<double, unsigned>{0.0, 0});
}

const double Generate::m_getAcceptanceThreshold(const int& t_upperBound) const {
    return m_acceptanceThreshold + std::pow(2.0 * t_upperBound, -0.5);
}

void Generate::m_singleRun() {
    //* Initialize for single ensemble
    NZ_Network network(m_networkSize);
    unsigned time = 0;
    int periodTime = 0;
    unsigned trialTime = 0;
    int periodTrialTime = 0;
    int eventTime = 0;
    int upperBound = 2;
    bool findNewLink = true;
    int root1, root2, node1, node2, newSize;
    std::set<int> findingClusterSizeDist = m_clusterSizeDist_orderParameter;
    std::set<int> newFindingClusterSizeDist = m_clusterSizeDist_orderParameter;
    std::set<int> findingClusterSizeDist_time = m_clusterSizeDist_time;
    std::set<unsigned> findingOrderParameterDist = m_orderParameterDist_time;
    std::vector<double> singleOrderParameter(m_networkSize, 0.0);

    //! Observables at time=0 state
    {
        obs_orderParameter[time] += 1.0 / m_networkSize;
        obs_secondMaximum[time] += 1.0 / m_networkSize;
        obs_secondMoment[time] += std::pow(1.0 / m_networkSize, 2.0);
        obs_meanClusterSize[time] += 1.0;
        // ++obs_subSuperEnsemble["sub"][time];
        // obs_netOrderParameter["sub"][time] += 1.0 / m_networkSize;
        // obs_netSecondMoment["sub"][time] += std::pow(1.0 / m_networkSize, 2.0);
        // singleOrderParameter[time] = 1.0 / m_networkSize;
    }

    //* Do m-BFW algorithm until all clusters are merged to one
    while (time < m_maxTime - 1) {
        //* Find new link to add
        if (findNewLink) {
            do {
                node1 = m_nodeDistribution(m_randomEngine);
                node2 = m_nodeDistribution(m_randomEngine);
            } while (node1 == node2);
            root1 = network.getRoot(m_nodeDistribution(m_randomEngine));
            root2 = network.getRoot(m_nodeDistribution(m_randomEngine));
            newSize = network.getSize(root1) + network.getSize(root2);
        }

        //* Chosen link is accepted: time, trial time increased
        if (newSize <= upperBound) {
            if (root1 == root2) {
                ++network.linkSize;
            } else {
                network.merge(root1, root2);
            }
            ++time;
            ++trialTime;
            ++periodTime;
            ++periodTrialTime;
            findNewLink = true;

            //* Get basic data from network
            const int maximumClusterSize = network.maximumClusterSize;
            const double orderParameter = maximumClusterSize / (double)m_networkSize;
            const int deltaMaximumClusterSize = network.deltaMaximumClusterSize;

            //! Sub_Super Ensemble
            {
                // ++obs_subSuperEnsemble[m_subSuperState[maximumClusterSize-1]][time];
            }

            //! (Net) Order Parameter
            {
                obs_orderParameter[time] += orderParameter;
                // obs_netOrderParameter[m_subSuperState[maximumClusterSize-1]][time] += orderParameter;
            }

            //! Single Order Parameter
            {
                // singleOrderParameter[time] = orderParameter;
            }

            //! (Period) Dynamics
            {
                // obs_dynamics.emplace_back(std::vector<int>{trialTime, time, maximumClusterSize, upperBound});
                // obs_periodDynamics.emplace_back(std::vector<int>{periodTrialTime, periodTime});
            }

            //! Second Maximum
            {
                obs_secondMaximum[time] += network.getSecondMaximumClusterSize() / (double)m_networkSize;
            }

            //! (Net) Second Moment
            {
                obs_secondMoment[time] += std::pow(orderParameter, 2.0);
                // obs_netSecondMoment[m_subSuperState[maximumClusterSize-1]][time] += std::pow(orderParameter, 2.0);
            }

            //! Mean Cluster Size
            {
                obs_meanClusterSize[time] += network.getMeanClusterSize();
            }

            //! Age Distribution (time)
            {
                for (const std::pair<unsigned long long, int>& changedAge : network.changedAge) {
                    obs_ageDist.at(m_state[maximumClusterSize - 1])[changedAge.first] += changedAge.second;
                }
                // for (const std::pair<unsigned long long, int>& changedAge : network.changedAge) {
                //     obs_ageDist_time.at(m_state_time[time])[changedAge.first] += changedAge.second;
                // }
            }

            //! Order Parameter Distribution
            {
                if (time == *findingOrderParameterDist.begin()) {
                    ++obs_orderParameterDist.at(time)[maximumClusterSize];
                    findingOrderParameterDist.erase(time);
                }
            }

            //! Cluster Size Distribution time
            {
                // if (time == *findingClusterSizeDist_time.begin()) {
                //     const std::map<int, int> sortedCluster = network.getSortedCluster();
                //     for (const std::pair<int, int>& cluster : sortedCluster) {
                //         obs_clusterSizeDist_time.at(time)[cluster.first] += cluster.second;
                //     }
                //     findingClusterSizeDist_time.erase(time);
                // }
            }

            //* Maximum cluster size of network is updated <=> Upper bound of m-BFW model is updated right before
            if (deltaMaximumClusterSize && time >= 2) {
                //* Get basic data
                const int interEventTime = time - eventTime;

                //! Inter Event Time
                {
                    obs_interEventTime[time].first += interEventTime;
                    ++obs_interEventTime[time].second;
                }

                //! Inter Event Time Distribution (time)
                {
                    ++obs_interEventTimeDist.at(m_state[maximumClusterSize - 1])[interEventTime];
                    // ++obs_interEventTimeDist_time.at(m_state_time[time])[interEventTime];
                }

                //! Delta Upper Bound Distribution (time)
                {
                    ++obs_deltaUpperBoundDist.at(m_state[maximumClusterSize - 1])[deltaMaximumClusterSize];
                    // ++obs_deltaUpperBoundDist_time.at(m_state_time[time])[deltaMaximumClusterSize];
                }

                //! Cluster Size Distribution
                {
                    findingClusterSizeDist = newFindingClusterSizeDist;
                    for (const int& op : findingClusterSizeDist) {
                        if (op < maximumClusterSize) {
                            const std::map<int, int> sortedCluster = network.getSortedCluster();
                            for (const std::pair<int, int>& cluster : sortedCluster) {
                                obs_clusterSizeDist.at(op)[cluster.first] += cluster.second;
                            }
                            newFindingClusterSizeDist.erase(op);
                        }
                    }
                }

                //! Inter Event Time - Order Parameter
                {
                    // obs_interEventTime_orderParameter[interEventTime].first += orderParameter;
                    // ++obs_interEventTime_orderParameter[interEventTime].second;
                    // obs_orderParameter_interEventTime[maximumClusterSize - 1].first += (double)interEventTime;
                    // ++obs_orderParameter_interEventTime[maximumClusterSize - 1].second;
                }

                //* Update event time
                eventTime = time;

            } //* End of maximum cluster size updated

        } //* End of accepting links

        //* Chosen link is rejected: only trial time increased
        else if ((double)time / trialTime > m_getAcceptanceThreshold(upperBound)) {
            ++trialTime;
            findNewLink = true;

            //! (Period) Dynamics
            {
                // ++periodTrialTime;
                // obs_dynamics.emplace_back(std::vector<int>{trialTime,
                //                                            time,
                //                                            network.maximumClusterSize,
                //                                            upperBound});
            }
        } //* End of rejecting links

        //* Upper bound is changing. Going to accept chosen link right after
        else {
            upperBound = newSize;
            findNewLink = false;

            //! (Period) Dynamics
            {
                // periodTime = 0;
                // periodTrialTime = 0;
                // obs_dynamics.emplace_back(std::vector<int>{trialTime,
                //                                            time,
                //                                            network.maximumClusterSize,
                //                                            upperBound});
            }
        } //* End of updating upper bound
    }

    //! Single Order Parameter
    {
        // obs_singleOrderParameter.emplace_back(singleOrderParameter);
    }
}

void Generate::run(const unsigned& t_ensembleSize) {
    m_ensembleSize = t_ensembleSize;
    //! Single Order Parameter
    {
        // obs_singleOrderParameter.reserve(t_ensembleSize);
    }
    for (unsigned ensemble = 0; ensemble < t_ensembleSize; ++ensemble) {
        m_singleRun();
    }
}

void Generate::save() const {
    using namespace linearAlgebra;
    const int precision = -1; //* Maximum precision
    const std::string NG = fileName::NG(m_networkSize, m_acceptanceThreshold, m_coreNum);
    const std::string NGE = fileName::NGE(m_networkSize, m_acceptanceThreshold, m_ensembleSize, m_coreNum);

    //! Order Parameter
    {
        const std::string directory = dataDirectory + "orderParameter/";
        CSV::generateDirectory(directory);
        const std::vector<double> orderParameter = obs_orderParameter / (double)m_ensembleSize;
        CSV::write(directory + NGE, orderParameter, precision);
    }

    //! Net Order Parameter
    {
        // for (const std::string& state : std::set<std::string>{"sub", "super"}){
        //     const std::string directory = dataDirectory + "netOrderParameter_mc78/" + state + "/";
        //     CSV::generateDirectory(directory);
        //     std::map<std::pair<int, unsigned>, double> trimmed;
        //     for (int t = 0; t < m_networkSize; ++t) {
        //         if (obs_subSuperEnsemble.at(state)[t]){
        //             const unsigned ensemble = obs_subSuperEnsemble.at(state)[t];
        //             trimmed[std::pair<int, unsigned>{t, ensemble}] = obs_netOrderParameter.at(state)[t] / (double)ensemble;
        //         }
        //     }
        //     CSV::write(directory + NG, trimmed, precision);
        // }
    }

    //! Single Order Parameter
    {
        // const std::string directory = dataDirectory + "singleOrderParameter/";
        // CSV::generateDirectory(directory);
        // for (unsigned ensemble=0; ensemble<m_ensembleSize; ++ensemble){
        //     CSV::write(directory + fileName::base(m_networkSize, m_acceptanceThreshold) + ",E" + std::to_string(ensemble) + ".txt", obs_singleOrderParameter[ensemble]);
        // }
    }

    //! (Period) Dynamics
    {
        // const std::string directory = dataDirectory + "dynamics/";
        // CSV::generateDirectory(directory);
        // CSV::write(directory + fileName::NG(m_networkSize, m_acceptanceThreshold, m_randomEngineSeed), obs_dynamics, precision);

        // const std::string periodDirectory = dataDirectory + "periodDynamics/";
        // CSV::generateDirectory(periodDirectory);
        // CSV::write(periodDirectory + fileName::NG(m_networkSize, m_acceptanceThreshold, m_randomEngineSeed), obs_periodDynamics, precision);
    }

    //! Second Maximum
    {
        const std::string directory = dataDirectory + "secondMaximum/";
        CSV::generateDirectory(directory);
        const std::vector<double> secondMaximum = obs_secondMaximum / (double)m_ensembleSize;
        CSV::write(directory + NGE, secondMaximum, precision);
    }

    //! Second Moment
    {
        const std::string directory = dataDirectory + "secondMoment/";
        CSV::generateDirectory(directory);
        const std::vector<double> secondMoment = obs_secondMoment / (double)m_ensembleSize;
        CSV::write(directory + NGE, secondMoment, precision);
    }

    //! Net Second Moment
    {
        // for (const std::string& state : std::set<std::string>{"sub", "super"}){
        //     const std::string directory = dataDirectory + "netSecondMoment_mc78/" + state + "/";
        //     CSV::generateDirectory(directory);
        //     std::map<std::pair<int, unsigned>, double> trimmed;
        //     for (int t = 0; t < m_networkSize; ++t) {
        //         if (obs_subSuperEnsemble.at(state)[t]){
        //             const unsigned ensemble = obs_subSuperEnsemble.at(state)[t];
        //             trimmed[std::pair<int, unsigned>{t, ensemble}] = obs_netSecondMoment.at(state)[t] / (double)ensemble;
        //         }
        //     }
        //     CSV::write(directory + NG, trimmed, precision);
        // }
    }

    //! Mean Cluster Size
    {
        const std::string directory = dataDirectory + "meanClusterSize/";
        CSV::generateDirectory(directory);
        std::vector<double> meanClusterSize = obs_meanClusterSize / (double)m_ensembleSize;
        meanClusterSize.back() = 0.0;
        CSV::write(directory + NGE, meanClusterSize, precision);
    }

    //! Inter Event Time
    {
        const std::string directory = dataDirectory + "interEventTime/";
        CSV::generateDirectory(directory);
        std::map<int, double> trimmed;
        for (int t = 0; t < m_networkSize; ++t) {
            if (obs_interEventTime[t].second) {
                trimmed[t] = obs_interEventTime[t].first / (double)obs_interEventTime[t].second;
            }
        }
        CSV::write(directory + NGE, trimmed, precision);
    }

    //! Age Distribution (time)
    {
        for (const std::string& state : mBFW::states) {
            const std::string directory = dataDirectory + "ageDist/" + state + "/";
            CSV::generateDirectory(directory);
            const std::vector<double> ageDist(obs_ageDist.at(state).begin(), obs_ageDist.at(state).end());
            CSV::write(directory + NGE, ageDist / (double)m_ensembleSize, precision);
        }

        // const std::vector<std::string> states = m_acceptanceThreshold == 1.0 ? std::vector<std::string>{"0_A1", "C_1"} : mBFW::states;
        // for (const std::string& state : states){
        //     const std::string directory = dataDirectory + "ageDist_time/" + state + "/";
        //     CSV::generateDirectory(directory);
        //     const std::vector<double> ageDist_time(obs_ageDist_time.at(state).begin(), obs_ageDist_time.at(state).end());
        //     CSV::write(directory + NGE, ageDist_time / (double)m_ensembleSize, precision);
        // }
    }

    //! Inter Event Time Distribution (time)
    {
        for (const std::string& state : mBFW::states) {
            const std::string directory = dataDirectory + "interEventTimeDist/" + state + "/";
            CSV::generateDirectory(directory);
            std::map<int, double> trimmed;
            for (int iet = 0; iet < m_networkSize; ++iet) {
                if (obs_interEventTimeDist.at(state)[iet]) {
                    trimmed[iet] = obs_interEventTimeDist.at(state)[iet] / (double)m_ensembleSize;
                }
            }
            CSV::write(directory + NGE, trimmed, precision);
        }

        // for (const std::string& state : mBFW::states) {
        //     const std::string directory = dataDirectory + "interEventTimeDist_time/" + state + "/";
        //     CSV::generateDirectory(directory);
        //     std::map<int, double> trimmed_time;
        //     for (int iet = 0; iet < m_networkSize; ++iet) {
        //         if (obs_interEventTimeDist_time.at(state)[iet]) {
        //             trimmed_time[iet] = obs_interEventTimeDist_time.at(state)[iet] / (double)m_ensembleSize;
        //         }
        //     }
        //     CSV::write(directory + NGE, trimmed_time, precision);
        // }
    }

    //! Delta Upper Bound Distribution (time)
    {
        for (const std::string& state : mBFW::states) {
            const std::string directory = dataDirectory + "deltaUpperBoundDist/" + state + "/";
            CSV::generateDirectory(directory);
            std::map<int, double> trimmed;
            for (int dk = 0; dk < m_networkSize; ++dk) {
                if (obs_deltaUpperBoundDist.at(state)[dk]) {
                    trimmed[dk] = obs_deltaUpperBoundDist.at(state)[dk] / (double)m_ensembleSize;
                }
            }
            CSV::write(directory + NGE, trimmed, precision);
        }

        // for (const std::string& state : mBFW::states) {
        //     const std::string directory = dataDirectory + "deltaUpperBoundDist_time/" + state + "/";
        //     CSV::generateDirectory(directory);
        //     std::map<int, double> trimmed_time;
        //     for (int dk = 0; dk < m_networkSize; ++dk) {
        //         if (obs_deltaUpperBoundDist_time.at(state)[dk]) {
        //             trimmed_time[dk] = obs_deltaUpperBoundDist_time.at(state)[dk] / (double)m_ensembleSize;
        //         }
        //     }
        //     CSV::write(directory + NGE, trimmed_time, precision);
        // }
    }

    //! Cluster Size Distribution
    {
        const std::string directory = dataDirectory + "clusterSizeDist/";
        CSV::generateDirectory(directory);
        std::vector<std::vector<double>> totalData;
        totalData.reserve(m_clusterSizeDist_orderParameter.size());
        for (const int& op : m_clusterSizeDist_orderParameter) {
            std::vector<double> clusterSizeDist;
            clusterSizeDist.reserve(2 + 2 * obs_clusterSizeDist.at(op).size());
            clusterSizeDist.emplace_back((double)op);
            clusterSizeDist.emplace_back((double)m_ensembleSize);
            for (const std::pair<int, unsigned long long>& csd : obs_clusterSizeDist.at(op)) {
                clusterSizeDist.emplace_back((double)csd.first);
                clusterSizeDist.emplace_back(csd.second / (double)m_ensembleSize);
            }
            totalData.emplace_back(clusterSizeDist);
        }
        CSV::write(directory + NG, totalData, precision);
    }

    //! Cluster Size Distribution time
    {
        // const std::string directory = dataDirectory + "clusterSizeDist_time/";
        // CSV::generateDirectory(directory);
        // std::vector<std::vector<double>> totalData_time;
        // totalData_time.reserve(m_clusterSizeDist_time.size());
        // for (const int& t : m_clusterSizeDist_time) {
        //     std::vector<double> clusterSizeDist_time;
        //     clusterSizeDist_time.reserve(2 + 2 * obs_clusterSizeDist_time.at(t).size());
        //     clusterSizeDist_time.emplace_back((double)t);
        //     clusterSizeDist_time.emplace_back((double)m_ensembleSize);
        //     for (const std::pair<int, unsigned long long>& csd : obs_clusterSizeDist_time.at(t)) {
        //         clusterSizeDist_time.emplace_back((double)csd.first);
        //         clusterSizeDist_time.emplace_back(csd.second / (double)m_ensembleSize);
        //     }
        //     totalData_time.emplace_back(clusterSizeDist_time);
        // }
        // CSV::write(directory + NG, totalData_time, precision);
    }

    //! Order Parameter Distribution
    {
        const std::string directory = dataDirectory + "orderParameterDist/";
        CSV::generateDirectory(directory);
        std::vector<std::vector<double>> totalData;
        totalData.reserve(m_orderParameterDist_time.size());
        for (const int& t : m_orderParameterDist_time) {
            std::vector<double> orderParameterDist;
            orderParameterDist.reserve(2 + 2 * obs_orderParameterDist.at(t).size());
            orderParameterDist.emplace_back((double)t);
            orderParameterDist.emplace_back((double)m_ensembleSize);
            for (const std::pair<int, unsigned>& opd : obs_orderParameterDist.at(t)) {
                orderParameterDist.emplace_back((double)opd.first);
                orderParameterDist.emplace_back(opd.second / (double)m_ensembleSize);
            }
            totalData.emplace_back(orderParameterDist);
        }
        CSV::write(directory + NG, totalData, precision);
    }

    //! Inter Event Time - Order Parameter
    {
        // const std::string directory = dataDirectory + "interEventTime_orderParameter/";
        // CSV::generateDirectory(directory);
        // std::map<int, double> trimmed;
        // for (int iet = 0; iet < m_networkSize; ++iet) {
        //     if (obs_interEventTime_orderParameter[iet].second) {
        //         trimmed[iet] = obs_interEventTime_orderParameter[iet].first / (double)obs_interEventTime_orderParameter[iet].second;
        //     }
        // }
        // CSV::write(directory + NGE, trimmed, precision);

        // const std::string directory2 = dataDirectory + "orderParameter_interEventTime/";
        // CSV::generateDirectory(directory2);
        // std::map<int, double> trimmed2;
        // for (int m = 0; m < m_networkSize; ++m) {
        //     if (obs_orderParameter_interEventTime[m].second) {
        //         trimmed2[m] = obs_orderParameter_interEventTime[m].first / (double)obs_orderParameter_interEventTime[m].second;
        //     }
        // }
        // CSV::write(directory2 + NGE, trimmed, precision);
    }
} //* End of function mBFW::Generate::save

} // namespace mBFW
