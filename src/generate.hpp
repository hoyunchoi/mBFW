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

#include "CSV.hpp"
#include "NZ_Network.hpp"
#include "common.hpp"
#include "linearAlgebra.hpp"
#include "parameters.hpp"
//*-------------------------

/*
    Observables
    obs_orderParameter[t] : Order parameter at time t
    obs_secondMoment[t] : Second moment (op^2) at time t
    obs_meanClusterSize[t] : Mean cluster size (of finite clusters) at time t
    obs_interEventTime[t] : Inter event time finished at time t, sample number
    obs_netOrderParameter[t] : Order Parameter greater than m_c
    obs_singleOrderParameter[ensemble][op] : ordre parameter for single ensemble

    obs_ageDist["state"][age] : Number of age accumulated at interval "state"
    obs_interEventTimeDist["state"][iet] : Number of inter event time accumulated through interval "state"
    obs_deltaUpperBoundDist["state"][dk] : Number of delta upper bound accumulated through interval "state"
    obs_clusterSizeDist[op][cs] : Number of clusters with size cs when order parameter exceeds op
    obs_orderParameterDist[t][op] : Number of ensembles with order parameter op at time t

    obs_interEventTime_orderParameter[iet] : Order parameter when inter event time is iet, sample number
    obs_orderParameter_interEventTime[m] : Inter Event time when maximum cluster size is m, sample number
*/

namespace mBFW {
struct Generate {
   protected:
    //* Member variables
    int m_networkSize;
    double m_acceptanceThreshold;
    unsigned m_ensembleSize;
    int m_coreNum;
    pcg32 m_randomEngine;
    std::uniform_int_distribution<int> m_nodeDistribution;
    std::map<std::string, int> m_points;
    std::vector<std::string> m_state;
    std::vector<std::string> m_state_time;
    std::vector<std::string> m_sub_super;
    std::set<int> m_clusterSizeDist_orderParameter;
    std::set<int> m_orderParameterDist_time;

    //* Observables
    std::vector<double> obs_orderParameter;
    std::vector<double> obs_secondMoment;
    std::vector<double> obs_meanClusterSize;
    std::map<std::string, std::vector<std::pair<double, unsigned>>> obs_netOrderParameter;
    std::vector<std::vector<double>> obs_singleOrderParameter;
    std::vector<std::pair<unsigned, unsigned>> obs_interEventTime;
    std::map<std::string, std::vector<unsigned long long>> obs_ageDist;
    std::map<std::string, std::vector<unsigned>> obs_interEventTimeDist;
    std::map<std::string, std::vector<unsigned>> obs_deltaUpperBoundDist;
    std::map<int, std::map<int, unsigned long long>> obs_clusterSizeDist;
    std::map<int, std::map<int, unsigned>> obs_orderParameterDist;
    std::vector<std::pair<double, unsigned>> obs_interEventTime_orderParameter;
    std::vector<std::pair<double, unsigned>> obs_orderParameter_interEventTime;

   public:
    Generate() {}
    Generate(const int&, const double&, const int&, const int& t_randomEngineSeed = -1);

    //* Member functions
    void run(const unsigned&);
    void save() const;

   protected:
    void m_singleRun();
};

Generate::Generate(const int& t_networkSize, const double& t_acceptanceThreshold, const int& t_coreNum, const int& t_randomEngineSeed) : m_networkSize(t_networkSize), m_acceptanceThreshold(t_acceptanceThreshold), m_coreNum(t_coreNum) {
    //* Get default values from parameter
    mBFW::Parameter parameter(t_networkSize, t_acceptanceThreshold);
    m_points = parameter.get_points();
    m_clusterSizeDist_orderParameter = parameter.get_clusterSizeDist_orderParameter();
    m_orderParameterDist_time = parameter.get_orderParameterDist_time();

    //* Initialize random variables
    t_randomEngineSeed == -1 ? m_randomEngine.seed((std::random_device())()) : m_randomEngine.seed(t_randomEngineSeed);
    m_nodeDistribution.param(std::uniform_int_distribution<int>::param_type(0, t_networkSize - 1));

    //* Initialize state vectors using points
    m_state.resize(t_networkSize);
    m_state_time.resize(t_networkSize);
    m_sub_super.resize(t_networkSize);
    {
        for (int m = 0; m < m_points.at("m_a1"); ++m) {
            m_state[m] = "0_A1";
            m_sub_super[m] = "sub";
        }
        for (int m = m_points.at("m_a1"); m < m_points.at("m_a2"); ++m) {
            m_state[m] = "A1_A2";
            m_sub_super[m] = "sub";
        }
        for (int m = m_points.at("m_a2"); m < m_points.at("m_b"); ++m) {
            m_state[m] = "A2_B";
            m_sub_super[m] = "sub";
        }
        for (int m = m_points.at("m_b"); m < m_points.at("m_c"); ++m) {
            m_state[m] = "B_C";
            m_sub_super[m] = "sub";
        }
        for (int m = m_points.at("m_c"); m < t_networkSize; ++m) {
            m_state[m] = "C_1";
            m_sub_super[m] = "super";
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
    }

    //* Initialize observables
    // obs_orderParameter.assign(t_networkSize, 0.0);
    // obs_meanClusterSize.assign(t_networkSize, 0.0);
    // obs_secondMoment.assign(t_networkSize, 0.0);
    // obs_interEventTime.assign(t_networkSize, std::pair<unsigned, unsigned>{0, 0});

    for (const std::string& state : std::set<std::string>{"sub", "super"}){
        obs_netOrderParameter[state].assign(t_networkSize, std::pair<double, unsigned>{0.0, 0});
    }

    // for (const std::string& state : mBFW::states) {
    //     obs_ageDist[state].assign(t_networkSize, 0);
    //     obs_interEventTimeDist[state].assign(t_networkSize, 0);
    //     obs_deltaUpperBoundDist[state].assign(t_networkSize, 0);
    // }
    // for (const double& op : m_clusterSizeDist_orderParameter) {
    //     obs_clusterSizeDist[op] = std::map<int, unsigned long long>{};
    // }
    // for (const int& t : m_orderParameterDist_time) {
    //     obs_orderParameterDist[t] = std::map<int, unsigned>{};
    // }

    // obs_interEventTime_orderParameter.assign(t_networkSize, std::pair<double, unsigned>{0.0, 0});
    // obs_orderParameter_interEventTime.assign(t_networkSize, std::pair<double, unsigned>{0.0, 0});
}

void Generate::m_singleRun() {
    //* Initialize for single ensemble
    NZ_Network network(m_networkSize);
    int time = 0;
    int trialTime = 0;
    int eventTime = 0;
    int upperBound = 2;
    bool findNewLink = true;
    int root1, root2, newSize;
    std::set<int> findingClusterSizeDist = m_clusterSizeDist_orderParameter;
    std::set<int> newFindingClusterSizeDist = m_clusterSizeDist_orderParameter;
    std::set<int> findingOrderParameterDist = m_orderParameterDist_time;
    std::vector<double> orderParameter_vec(m_networkSize, 0.0);

    //! Observables at time=0 state
    // obs_orderParameter[time] += 1.0 / m_networkSize;
    // obs_secondMoment[time] += std::pow(1.0 / m_networkSize, 2.0);
    // obs_meanClusterSize[time] += 1.0;
    orderParameter_vec[time] = 1.0 / m_networkSize;

    //* Do m-BFW algorithm until all clusters are merged to one
    while (time < m_networkSize - 1) {
        //* Find new link to add
        if (findNewLink) {
            do {
                root1 = network.getRoot(m_nodeDistribution(m_randomEngine));
                root2 = network.getRoot(m_nodeDistribution(m_randomEngine));
            } while (root1 == root2);
            newSize = network.getSize(root1) + network.getSize(root2);
        }

        //* Chosen link is accepted: time, trial time increased
        if (newSize <= upperBound) {
            network.merge(root1, root2);
            ++time;
            ++trialTime;
            findNewLink = true;

            //* Get basic data from network
            const int maximumClusterSize = network.maximumClusterSize;
            const double orderParameter = (double)maximumClusterSize / m_networkSize;
            const int deltaMaximumClusterSize = network.deltaMaximumClusterSize;

            //! Order Parameter
            {
                // obs_orderParameter[time] += orderParameter;
            }

            //! Net Order Parameter
            {
                // obs_netOrderParameter[m_sub_super[maximumClusterSize-1]][time].first += orderParameter;
                // ++obs_netOrderParameter[m_sub_super[maximumClusterSize-1]][time].second;
            }

            //! Single Order Parameter
            {
                orderParameter_vec[time] = orderParameter;
            }

            //! Second Moment
            {
                // obs_secondMoment[time] += std::pow(orderParameter, 2.0);
            }

            //! Mean Cluster Size
            {
                // obs_meanClusterSize[time] += network.getMeanClusterSize();
            }

            //! Age Distribution
            {
                // for (const std::pair<unsigned long long, int>& changedAge : network.changedAge) {
                //     obs_ageDist.at(m_state[maximumClusterSize - 1])[changedAge.first] += changedAge.second;
                // }
            }

            //! Order Parameter Distribution
            {
                // if (time == *findingOrderParameterDist.begin()) {
                //     ++obs_orderParameterDist.at(time)[maximumClusterSize];
                //     findingOrderParameterDist.erase(time);
                // }
            }

            //* Maximum cluster size of network is updated <=> Upper bound of m-BFW model is updated right before
            if (deltaMaximumClusterSize && time >= 2) {
                //* Get basic data
                const int interEventTime = time - eventTime;

                //! Inter Event Time
                {
                    // obs_interEventTime[time].first += interEventTime;
                    // ++obs_interEventTime[time].second;
                }

                //! Inter Event Time Distribution
                {
                    // ++obs_interEventTimeDist.at(m_state[maximumClusterSize - 1])[interEventTime];
                }

                //! Delta Upper Bound Distribution
                {
                    // ++obs_deltaUpperBoundDist.at(m_state[maximumClusterSize - 1])[deltaMaximumClusterSize];
                }

                //! Cluster Size Distribution
                {
                    // findingClusterSizeDist = newFindingClusterSizeDist;
                    // for (const int& op : findingClusterSizeDist) {
                    //     if (op < maximumClusterSize) {
                    //         const std::map<int, int> sortedCluster = network.getSortedCluster();
                    //         for (const std::pair<int, int>& cluster : sortedCluster) {
                    //             obs_clusterSizeDist.at(op)[cluster.first] += cluster.second;
                    //         }
                    //         newFindingClusterSizeDist.erase(op);
                    //     }
                    // }
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

            }  //* End of maximum cluster size updated

        }  //* End of accepting links

        //* Chosen link is rejected: only trial time increased
        else if ((double)time / trialTime > m_acceptanceThreshold) {
            ++trialTime;
            findNewLink = true;
        }  //* End of rejecting links

        //* Upper bound is changing. Going to accept chosen link right after
        else {
            upperBound = newSize;
            findNewLink = false;
        }  //* End of updating upper bound
    }

    //! Order Parameter single
    {
        obs_singleOrderParameter.emplace_back(orderParameter_vec);
    }
}

void Generate::run(const unsigned& t_ensembleSize) {
    m_ensembleSize = t_ensembleSize;
    obs_singleOrderParameter.reserve(t_ensembleSize);
    for (unsigned ensemble = 0; ensemble < t_ensembleSize; ++ensemble) {
        m_singleRun();
    }
}

void Generate::save() const {
    using namespace linearAlgebra;
    const int precision = -1;  //* Maximum precision
    const std::string NG = fileName::NG(m_networkSize, m_acceptanceThreshold, m_coreNum);
    const std::string NGE = fileName::NGE(m_networkSize, m_acceptanceThreshold, m_ensembleSize, m_coreNum);

    //! Order Parameter
    {
        // const std::string directory = dataDirectory + "orderParameter/";
        // CSV::generateDirectory(directory);
        // const std::vector<double> orderParameter = obs_orderParameter / (double)m_ensembleSize;
        // CSV::write(directory + NGE, orderParameter, precision);
    }

    //! Net Order Parameter
    {
        // for (const std::string& state : std::set<std::string>{"sub", "super"}){
        //     const std::string directory = dataDirectory + "netOrderParameter/" + state + "/";
        //     CSV::generateDirectory(directory);
        //     std::map<std::pair<int, unsigned>, double> trimmed;
        //     for (int t = 0; t < m_networkSize; ++t) {
        //         if (obs_netOrderParameter.at(state)[t].second) {
        //             trimmed[std::pair<int, unsigned>{t, obs_netOrderParameter.at(state)[t].second}] = obs_netOrderParameter.at(state)[t].first / (double)obs_netOrderParameter.at(state)[t].second;
        //         }
        //     }
        //     CSV::write(directory + NG, trimmed, precision);
        // }
    }

    //! Single Order Parameter
    {
        const std::string directory = dataDirectory + "singleOrderParameter/";
        CSV::generateDirectory(directory);
        for (unsigned ensemble=0; ensemble<m_ensembleSize; ++ensemble){
            CSV::write(directory + fileName::base(m_networkSize, m_acceptanceThreshold) + ",E" + std::to_string(ensemble) + ".txt", obs_singleOrderParameter[ensemble]);
        }
    }

    //! Second Moment
    {
        // const std::string directory = dataDirectory + "orderParameterVariance/";
        // CSV::generateDirectory(directory);
        // const std::vector<double> orderParameterVariance = elementPow(obs_secondMoment / (double)m_ensembleSize - elementPow(obs_orderParameter / (double)m_ensembleSize, 2.0), 0.5);
        // CSV::write(directory + NGE, orderParameterVariance, precision);
    }

    //! Mean Cluster Size
    {
        // const std::string directory = dataDirectory + "meanClusterSize/";
        // CSV::generateDirectory(directory);
        // std::vector<double> meanClusterSize = obs_meanClusterSize / (double)m_ensembleSize;
        // meanClusterSize.back() = 0.0;
        // CSV::write(directory + NGE, meanClusterSize, precision);
    }

    //! Inter Event Time
    {
        // const std::string directory = dataDirectory + "interEventTime/";
        // CSV::generateDirectory(directory);
        // std::map<int, double> trimmed;
        // for (int t = 0; t < m_networkSize; ++t) {
        //     if (obs_interEventTime[t].second) {
        //         trimmed[t] = obs_interEventTime[t].first / (double)obs_interEventTime[t].second;
        //     }
        // }
        // CSV::write(directory + NGE, trimmed, precision);
    }

    //! Age Distribution
    {
        // for (const std::string& state : mBFW::states) {
        //     const std::string directory = dataDirectory + "ageDist/" + state + "/";
        //     CSV::generateDirectory(directory);
        //     const std::vector<double> ageDist(obs_ageDist.at(state).begin(), obs_ageDist.at(state).end());
        //     CSV::write(directory + NGE, ageDist / (double)m_ensembleSize, precision);
        // }
    }

    //! Inter Event Time Distribution
    {
        // for (const std::string& state : mBFW::states) {
        //     const std::string directory = dataDirectory + "interEventTimeDist/" + state + "/";
        //     CSV::generateDirectory(directory);
        //     std::map<int, double> trimmed;
        //     for (int iet = 0; iet < m_networkSize; ++iet) {
        //         if (obs_interEventTimeDist.at(state)[iet]) {
        //             trimmed[iet] = obs_interEventTimeDist.at(state)[iet] / (double)m_ensembleSize;
        //         }
        //     }
        //     CSV::write(directory + NGE, trimmed, precision);
        // }
    }

    //! Delta Upper Bound Distribution
    {
        // for (const std::string& state : mBFW::states) {
        //     const std::string directory = dataDirectory + "deltaUpperBoundDist/" + state + "/";
        //     CSV::generateDirectory(directory);
        //     std::map<int, double> trimmed;
        //     for (int dk = 0; dk < m_networkSize; ++dk) {
        //         if (obs_deltaUpperBoundDist.at(state)[dk]) {
        //             trimmed[dk] = obs_deltaUpperBoundDist.at(state)[dk] / (double)m_ensembleSize;
        //         }
        //     }
        //     CSV::write(directory + NGE, trimmed, precision);
        // }
    }

    //! Cluster Size Distribution
    {
        // const std::string directory = dataDirectory + "clusterSizeDist/";
        // CSV::generateDirectory(directory);
        // std::vector<std::vector<double>> totalData;
        // totalData.reserve(m_clusterSizeDist_orderParameter.size());
        // for (const int& op : m_clusterSizeDist_orderParameter) {
        //     std::vector<double> clusterSizeDist;
        //     clusterSizeDist.reserve(2 + 2 * obs_clusterSizeDist.at(op).size());
        //     clusterSizeDist.emplace_back((double)op);
        //     clusterSizeDist.emplace_back((double)m_ensembleSize);
        //     for (const std::pair<int, unsigned long long>& csd : obs_clusterSizeDist.at(op)) {
        //         clusterSizeDist.emplace_back((double)csd.first);
        //         clusterSizeDist.emplace_back(csd.second / (double)m_ensembleSize);
        //     }
        //     totalData.emplace_back(clusterSizeDist);
        // }
        // CSV::write(directory + NG, totalData, precision);
    }

    //! Order Parameter Distribution
    {
        // const std::string directory = dataDirectory + "orderParameterDist/";
        // CSV::generateDirectory(directory);
        // std::vector<std::vector<double>> totalData;
        // totalData.reserve(m_orderParameterDist_time.size());
        // for (const int& t : m_orderParameterDist_time) {
        //     std::vector<double> orderParameterDist;
        //     orderParameterDist.reserve(2 + 2 * obs_orderParameterDist.at(t).size());
        //     orderParameterDist.emplace_back((double)t);
        //     orderParameterDist.emplace_back((double)m_ensembleSize);
        //     for (const std::pair<int, unsigned>& opd : obs_orderParameterDist.at(t)) {
        //         orderParameterDist.emplace_back((double)opd.first);
        //         orderParameterDist.emplace_back(opd.second / (double)m_ensembleSize);
        //     }
        //     totalData.emplace_back(orderParameterDist);
        // }
        // CSV::write(directory + NG, totalData, precision);
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
}  //* End of function mBFW::Generate::save

}  // namespace mBFW
