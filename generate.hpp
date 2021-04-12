#pragma once

#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "../library/CSV.hpp"
#include "../library/linearAlgebra.hpp"
#include "../library/pcg_random.hpp"
#include "NZ_Network.hpp"
#include "common.hpp"
#include "parameters.hpp"

/*
    Observables
    obs_orderParameter[t] : Order parameter at time t
    obs_secondMoment[t] : Second moment (op^2) at time t
    obs_meanClusterSize[t] : Mean cluster size (of finite clusters) at time t
    obs_interEventTime[t] : Inter event time finished at time t, sample number

    obs_ageDist["state"][age] : Number of age accumulated at interval "state"
    obs_interEventTimeDist["state"][iet] : Number of inter event time accumulated through interval "state"
    obs_deltaUpperBoundDist["state"][dk] : Number of delta upper bound accumulated through interval "state"
    obs_clusterSizeDist[op][cs] : Number of clusters with size cs when order parameter exceeds op
    obs_orderParameterDist[t][op] : Number of ensembles with order parameter op at time t

    obs_interEventTime_orderParameter[iet] : Order parameter when inter event time is iet, sample number
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
    std::set<int> m_clusterSizeDist_orderParameter;
    std::set<int> m_orderParameterDist_time;

    //* Observables
    std::vector<double> obs_orderParameter;
    std::vector<double> obs_secondMoment;
    std::vector<double> obs_meanClusterSize;
    std::vector<std::pair<unsigned, unsigned>> obs_interEventTime;
    std::map<std::string, std::vector<unsigned long long>> obs_ageDist;
    std::map<std::string, std::vector<unsigned>> obs_interEventTimeDist;
    std::map<std::string, std::vector<unsigned>> obs_deltaUpperBoundDist;
    std::map<int, std::map<int, unsigned long long>> obs_clusterSizeDist;
    std::map<int, std::map<int, unsigned>> obs_orderParameterDist;
    std::vector<std::pair<double, unsigned>> obs_interEventTime_orderParameter;

   public:
    Generate() {}
    Generate(const int&, const double&, const int&, const int&);

    //* Member functions
    void run(const unsigned&);
    void save() const;

   protected:
    const std::string m_getState(const int&) const;
    const std::string m_getState_time(const int&) const;
    void m_singleRun();
};

Generate::Generate(const int& t_networkSize, const double& t_acceptanceThreshold, const int& t_coreNum, const int& t_randomEngineSeed = -1) : m_networkSize(t_networkSize), m_acceptanceThreshold(t_acceptanceThreshold), m_coreNum(t_coreNum) {
    //* Get default values from parameter
    mBFW::Parameter parameter(t_networkSize, t_acceptanceThreshold);
    m_points = parameter.get_points();
    m_clusterSizeDist_orderParameter = parameter.get_clusterSizeDist_orderParameter();
    m_orderParameterDist_time = parameter.get_orderParameterDist_time();

    //* Initialize random variables
    t_randomEngineSeed == -1 ? m_randomEngine.seed((std::random_device())()) : m_randomEngine.seed(t_randomEngineSeed);
    m_nodeDistribution.param(std::uniform_int_distribution<int>::param_type(0, t_networkSize - 1));

    //* Initialize observables
    obs_orderParameter.assign(t_networkSize, 0.0);
    obs_meanClusterSize.assign(t_networkSize, 0.0);
    obs_secondMoment.assign(t_networkSize, 0.0);
    obs_interEventTime.assign(t_networkSize, std::pair<unsigned, unsigned>{0, 0});

    for (const std::string& state : mBFW::states) {
        obs_ageDist[state].assign(t_networkSize, 0);
        obs_interEventTimeDist[state].assign(t_networkSize, 0);
        obs_deltaUpperBoundDist[state].assign(t_networkSize, 0);
    }
    for (const double& op : m_clusterSizeDist_orderParameter) {
        obs_clusterSizeDist[op] = std::map<int, unsigned long long>{};
    }
    for (const int& t : m_orderParameterDist_time) {
        obs_orderParameterDist[t] = std::map<int, unsigned>{};
    }

    obs_interEventTime_orderParameter.assign(t_networkSize, std::pair<double, unsigned>{0.0, 0});
}

const std::string Generate::m_getState(const int& t_maximumClusterSize) const {
    if (t_maximumClusterSize < m_points.at("m_a1")) {
        return "0A1";
    } else if (t_maximumClusterSize < m_points.at("m_a2")){
        return "A1A2";
    } else if (t_maximumClusterSize < m_points.at("m_b")) {
        return "A2B";
    } else if (t_maximumClusterSize < m_points.at("m_c")) {
        return "BC";
    } else {
        return "C1";
    }
}

const std::string Generate::m_getState_time(const int& t_time) const {
    if (t_time < m_points.at("t_a1")) {
        return "0A1";
    } else if (t_time < m_points.at("t_a2")){
        return "A1A2";
    } else if (t_time < m_points.at("t_b")) {
        return "A2B";
    } else if (t_time < m_points.at("t_c")) {
        return "BC";
    } else {
        return "C1";
    }
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
    std::string state = "0A1";
    std::string state_time = "0A1";
    std::set<int> findingClusterSizeDist = m_clusterSizeDist_orderParameter;
    std::set<int> newFindingClusterSizeDist = m_clusterSizeDist_orderParameter;
    std::set<int> findingOrderParameterDist = m_orderParameterDist_time;

    //* Observables at time=0 state
    obs_orderParameter[time] += 1.0 / m_networkSize;
    obs_secondMoment[time] += std::pow(1.0 / m_networkSize, 2.0);
    obs_meanClusterSize[time] += 1.0;

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
            state_time = m_getState_time(time);

            //* Get basic data from network
            const int maximumClusterSize = network.maximumClusterSize;
            const double orderParameter = (double)maximumClusterSize / m_networkSize;
            const int deltaMaximumClusterSize = network.deltaMaximumClusterSize;

            //! Order Parameter
            {
                obs_orderParameter[time] += orderParameter;
            }
            //! Second Moment
            {
                obs_secondMoment[time] += std::pow(orderParameter, 2.0);
            }

            //! Mean Cluster Size
            {
                obs_meanClusterSize[time] += network.getMeanClusterSize();
            }

            //! Age Distribution
            {
                for (const std::pair<unsigned long long, int>& changedAge : network.changedAge) {
                    obs_ageDist.at(state)[changedAge.first] += changedAge.second;
                }
            }

            //! Order Parameter Distribution
            {
                if (time == *findingOrderParameterDist.begin()) {
                    ++obs_orderParameterDist.at(time)[maximumClusterSize];
                    findingOrderParameterDist.erase(time);
                }
            }

            //* Maximum cluster size of network is updated <=> Upper bound of m-BFW model is updated right before
            if (deltaMaximumClusterSize && time >= 2) {
                state = m_getState(maximumClusterSize);
                //* Get basic data
                const int interEventTime = time - eventTime;

                //! Inter Event Time
                {
                    obs_interEventTime[time].first += interEventTime;
                    ++obs_interEventTime[time].second;
                }

                //! Inter Event Time Distribution
                {
                    ++obs_interEventTimeDist.at(state)[interEventTime];
                }

                //! Delta Upper Bound Distribution
                {
                    ++obs_deltaUpperBoundDist.at(state)[deltaMaximumClusterSize];
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
                    obs_interEventTime_orderParameter[interEventTime].first += orderParameter;
                    ++obs_interEventTime_orderParameter[interEventTime].second;
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
}

void Generate::run(const unsigned& t_ensembleSize) {
    m_ensembleSize = t_ensembleSize;
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
        const std::string directory = dataDirectory + "orderParameter/";
        CSV::generateDirectory(directory);
        const std::vector<double> orderParameter = obs_orderParameter / (double)m_ensembleSize;
        CSV::write(directory + NGE, orderParameter, precision);
    }

    //! Second Moment
    {
        const std::string directory = dataDirectory + "orderParameterVariance/";
        CSV::generateDirectory(directory);
        const std::vector<double> orderParameterVariance = elementPow(obs_secondMoment / (double)m_ensembleSize - elementPow(obs_orderParameter / (double)m_ensembleSize, 2.0), 0.5);
        CSV::write(directory + NGE, orderParameterVariance, precision);
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

    //! Age Distribution
    {
        for (const std::string& state : mBFW::states) {
            const std::string directory = dataDirectory + "ageDist/" + state + "/";
            CSV::generateDirectory(directory);
            const std::vector<double> ageDist(obs_ageDist.at(state).begin(), obs_ageDist.at(state).end());
            CSV::write(directory + NGE, ageDist / (double)m_ensembleSize, precision);
        }
    }

    //! Inter Event Time Distribution
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
    }

    //! Delta Upper Bound Distribution
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
        const std::string directory = dataDirectory + "interEventTime_orderParameter/";
        CSV::generateDirectory(directory);
        std::map<int, double> trimmed;
        for (int iet = 0; iet < m_networkSize; ++iet) {
            if (obs_interEventTime_orderParameter[iet].second) {
                trimmed[iet] = obs_interEventTime_orderParameter[iet].first / (double)obs_interEventTime_orderParameter[iet].second;
            }
        }
        CSV::write(directory + NGE, trimmed, precision);
    }
}  //* End of function mBFW::Generate::save

}  // namespace mBFW
