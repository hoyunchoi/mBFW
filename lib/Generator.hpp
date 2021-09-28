#pragma once

#include <map>
#include <random>
#include <vector>

#include "CSV.hpp"
#include "CheckList.hpp"
#include "NZ_Network.hpp"
#include "Name.hpp"
#include "linearAlgebra.hpp"
#include "pcg_random.hpp"

namespace mBFW {
constexpr CheckList<std::string_view, bool, observables.size()> check_list{{observables}};

struct Generator {
   protected:
    //* Basic member variables
    Name m_name;
    int m_network_size;
    double m_acceptance_threshold;
    unsigned m_ensemble_size;
    int m_core_num;

    //* Random variables
    int m_random_engine_seed;
    pcg32 m_random_engine;
    std::uniform_int_distribution<int> m_node_distribution;

    //* Observable
    std::vector<double> obs_order_parameter;
    std::vector<double> obs_second_maximum;

   public:
    Generator() {}
    Generator(const Name&, const int& = -1);

    void initialize();
    void run();
    void save() const;

   protected:
    void _run();
    double _gk(const double&) const;
};

Generator::Generator(const Name& t_name,
                     const int& t_random_engine_seed) : m_name(t_name),
                                                        m_network_size(t_name.network_size),
                                                        m_acceptance_threshold(t_name.acceptance_threshold),
                                                        m_ensemble_size(t_name.ensemble_size),
                                                        m_core_num(t_name.core_num),
                                                        m_random_engine_seed(t_random_engine_seed) {
    //* Initialize random variables
    m_random_engine_seed == -1
        ? m_random_engine.seed((std::random_device())())
        : m_random_engine.seed(t_random_engine_seed);
    m_node_distribution.param(std::uniform_int_distribution<int>::param_type(0, m_network_size - 1));

    initialize();
}

void Generator::initialize() {
    /*
        Initialize obesrvables according to check list
    */
    if (check_list.at("order_parameter")) {
        obs_order_parameter.assign(m_network_size, 0.0);
    }
    if (check_list.at("second_maximum")) {
        obs_second_maximum.assign(m_network_size, 0.0);
    }
}

void Generator::run() {
    /*
        do _run for ensemble times
    */
    for (unsigned ensemble = 0; ensemble < m_ensemble_size; ++ensemble) {
        _run();
    }
}

void Generator::save() const {
    using namespace linearAlgebra;
    const int precision = -1;  //* Maximum precision at writing to csv file

    //! Order parameter
    if (check_list.at("order_parameter")) {
        const std::string directory = data_dir + "order_parameter/";
        CSV::generateDirectory(directory);
        const std::vector<double> order_parameter = obs_order_parameter / (double)m_ensemble_size;
        CSV::write(directory + m_name.get_NGE(), order_parameter, precision);
    }
    //! Second Maximum
    if (check_list.at("second_maximum")) {
        const std::string directory = data_dir + "second_maximum/";
        CSV::generateDirectory(directory);
        const std::vector<double> second_maximum = obs_second_maximum / (double)m_ensemble_size;
        CSV::write(directory + m_name.get_NGE(), second_maximum, precision);
    }
}

double Generator::_gk(const double& t_upperBound) const {
    return 0.5 + sqrt(0.5 / t_upperBound);
}

void Generator::_run() {
    //* Declaration of parameters for single ensemble
    int time = 0;
    int trial_time = 0;
    int upper_bound = 2;
    bool find_new_link = true;
    NZ_Network network(m_network_size);
    int root1, root2, newSize;

    //! Observables at time=0
    {
        if (check_list.at("order_parameter")) {
            obs_order_parameter[time] += 1.0 / m_network_size;
        }
        if (check_list.at("second_maximum")) {
            obs_second_maximum[time] += 1.0 / m_network_size;
        }
    }

    //* Do mBFW algorithm until all clusters are merged into one
    while (network.maximum_cluster_size < m_network_size) {
        //* Choose new candidate link
        if (find_new_link) {
            do {
                root1 = network.get_root(m_node_distribution(m_random_engine));
                root2 = network.get_root(m_node_distribution(m_random_engine));
            } while (root1 == root2);
            newSize = network.get_size(root1) + network.get_size(root2);
        }

        //* --------------------- Chosen link is accepted ---------------------
        if (newSize <= upper_bound) {
            network.merge(root1, root2);
            ++time;
            ++trial_time;
            find_new_link = true;

            //* Get basic data from network
            const int maximum_cluster_size = network.maximum_cluster_size;
            const double order_parameter = maximum_cluster_size / (double)m_network_size;
            const double second_maximum = network.get_second_maximum_cluster_size() / (double)m_network_size;
            const int delta_maximum_cluster_size = network.delta_maximum_cluster_size;

            //! Order Parameter
            if (check_list.at("order_parameter")) {
                obs_order_parameter[time] += order_parameter;
            }
            //! Second Maximum
            if (check_list.at("second_maximum")) {
                obs_second_maximum[time] += second_maximum;
            }

        }  //* --------------------- End of accepting candidate ---------------------

        //* --------------------- Chosen link is rejected ---------------------
        else if (time / (double)trial_time > _gk(upper_bound)) {
            ++trial_time;
            find_new_link = true;
        }  //* --------------------- End of rejecting candidate ---------------------

        //* --------------------- Upper bound is changing ---------------------
        else {
            ++upper_bound;
            find_new_link = true;
        }  //* --------------------- End of updating upper bound ---------------------
    }
}

}  // namespace mBFW