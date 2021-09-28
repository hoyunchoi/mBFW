#pragma once

#include <map>
#include <random>
#include <utility>
#include <vector>

#include "CSV.hpp"
#include "CheckList.hpp"
#include "NZ_Network.hpp"
#include "Name.hpp"
#include "linearAlgebra.hpp"
#include "pcg_random.hpp"

/*
    Observables
    obs_maximum_cluster_size[t] : Order parameter at time t
    obs_second_maximum_cluster_size[t] : Second maximum cluster at time t
    obs_second_moment[t] : Second moment (op^2) at time t
    obs_mean_cluster_size[t] : Mean cluster size (of finite clusters) at time t
    obs_inter_event_time[t] : Inter event time finished at time t, sample number

    obs_max_delta_order_parameter[e]: pair of <t, delta OP> for each enesmble e
*/

namespace mBFW {
constexpr CheckList<std::string_view, bool, observables.size()> check_list{{observables}};

struct Generator {
   public:
    //* Basic member variables
    Name name;
    int network_size;
    double acceptance_threshold;
    unsigned ensemble_size;
    int core_num;

   protected:
    //* Random variables
    int _random_engine_seed;
    pcg32 _random_engine;
    std::uniform_int_distribution<int> _node_distribution;

    //* Observable
    std::vector<unsigned long long> obs_maximum_cluster_size;
    std::vector<unsigned long long> obs_second_maximum_cluster_size;
    std::vector<double> obs_second_moment;
    std::vector<double> obs_mean_cluster_size;
    std::vector<std::pair<unsigned long long, unsigned>> obs_inter_event_time;

    std::vector<std::vector<double>> obs_max_delta_order_parameter;

   public:
    Generator() {}
    Generator(const Name&, const int& = -1);

    void initialize();
    void run();
    void save() const;

   protected:
    void _run();
};

Generator::Generator(const Name& t_name,
                     const int& t_random_engine_seed) : name(t_name),
                                                        network_size(t_name.network_size),
                                                        acceptance_threshold(t_name.acceptance_threshold),
                                                        ensemble_size(t_name.ensemble_size),
                                                        core_num(t_name.core_num),
                                                        _random_engine_seed(t_random_engine_seed) {
    //* Initialize random variables
    _random_engine_seed == -1
        ? _random_engine.seed((std::random_device())())
        : _random_engine.seed(t_random_engine_seed);
    _node_distribution.param(std::uniform_int_distribution<int>::param_type(0, network_size - 1));

    initialize();
}

void Generator::initialize() {
    /*
        Initialize obesrvables according to check list
    */
    //* Observable with plain average
    if (check_list.at("maximum_cluster_size")) {
        obs_maximum_cluster_size.assign(network_size, 0);
    }
    if (check_list.at("second_maximum_cluster_size")) {
        obs_second_maximum_cluster_size.assign(network_size, 0);
    }
    if (check_list.at("mean_cluster_size")) {
        obs_mean_cluster_size.assign(network_size, 0.0);
    }
    if (check_list.at("second_moment")) {
        obs_second_moment.assign(network_size, 0.0);
    }

    //* Observable with selective average
    if (check_list.at("inter_event_time")) {
        obs_inter_event_time.assign(network_size, std::pair<unsigned long long, unsigned>{0, 0});
    }

    //* Accumulated Observable
    if (check_list.at("max_delta_order_parameter")) {
        obs_max_delta_order_parameter.reserve(ensemble_size);
    }
}

void Generator::run() {
    /*
        do _run for ensemble times
    */
    for (unsigned ensemble = 0; ensemble < ensemble_size; ++ensemble) {
        _run();
    }
}

void Generator::_run() {
    //* Declaration of parameters for single ensemble
    int time = 0;
    int trial_time = 0;
    int upper_bound = 2;
    bool find_new_link = true;
    NZ_Network network(network_size);
    int root1, root2, newSize;

    //! Inter Event Time
    int event_time = 0;

    //! Max Delta Order Parameter
    double max_delta_order_parameter_t = 0.0;
    double max_delta_order_parameter_m = 0.0;

    //* Do mBFW algorithm until all clusters are merged into one
    while (network.maximum_cluster_size < network_size) {
        //* Choose new candidate link
        if (find_new_link) {
            do {
                root1 = network.get_root(_node_distribution(_random_engine));
                root2 = network.get_root(_node_distribution(_random_engine));
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
            const double order_parameter = maximum_cluster_size / (double)network_size;
            const int second_maximum_cluster_size = network.get_second_maximum_cluster_size();
            const int delta_maximum_cluster_size = network.delta_maximum_cluster_size;

            //! Order Parameter
            if (check_list.at("maximum_cluster_size")) {
                obs_maximum_cluster_size[time] += maximum_cluster_size;
            }
            //! Second Maximum
            if (check_list.at("second_maximum_cluster_size")) {
                obs_second_maximum_cluster_size[time] += second_maximum_cluster_size;
            }
            //! Mean Cluster Size
            if (check_list.at("mean_cluster_size")) {
                obs_mean_cluster_size[time] += network.get_mean_cluster_size();
            }
            //! Second Moment
            if (check_list.at("second_moment")) {
                obs_second_moment[time] += pow(order_parameter, 2.0);
            }

            //* ------------ Maximum cluster size updated ------------
            if (delta_maximum_cluster_size && time >= 2) {
                const int inter_event_time = time - event_time;
                event_time = time;

                //! Inter Event Time
                if (check_list.at("inter_event_time")) {
                    obs_inter_event_time[time].first += inter_event_time;
                    ++obs_inter_event_time[time].second;
                }

                //! Max Delta Order Parameter
                if (check_list.at("max_delta_order_parameter")) {
                    if (max_delta_order_parameter_m < delta_maximum_cluster_size / (double)network_size) {
                        max_delta_order_parameter_t = time / (double)network_size;
                        max_delta_order_parameter_m = delta_maximum_cluster_size / (double)network_size;
                    }
                }
            }  //* ---------- End of maximum cluster size update ------------

        }  //* --------------------- End of accepting candidate ---------------------

        //* --------------------- Chosen link is rejected ---------------------
        else if (time / (double)trial_time > acceptance_threshold) {
            ++trial_time;
            find_new_link = true;
        }  //* --------------------- End of rejecting candidate ---------------------

        //* --------------------- Upper bound is changing ---------------------
        else {
            upper_bound = newSize;
            find_new_link = true;
        }  //* --------------------- End of updating upper bound ---------------------
    }      //* End of iteration

    //! Max Delta Order Parameter
    if (check_list.at("max_delta_order_parameter")) {
        obs_max_delta_order_parameter.emplace_back(std::vector<double>{max_delta_order_parameter_t, max_delta_order_parameter_m});
    }
}

void Generator::save() const {
    using namespace linearAlgebra;
    const int precision = -1;  //* Maximum precision at writing to csv file

    //! maximum cluster size
    if (check_list.at("maximum_cluster_size")) {
        const std::string directory = data_dir + "maximum_cluster_size/";
        CSV::generateDirectory(directory);
        std::vector<double> order_parameter(network_size);
        for (int i = 1; i < network_size; ++i) {
            order_parameter[i] = obs_maximum_cluster_size[i] / (double)ensemble_size;
        }
        order_parameter[0] = 1.0;

        CSV::write(directory + name.get_NGE(), order_parameter, precision);
    }
    //! Second Maximum
    if (check_list.at("second_maximum_cluster_size")) {
        const std::string directory = data_dir + "second_maximum_cluster_size/";
        CSV::generateDirectory(directory);
        std::vector<double> second_maximum_cluster_size(network_size);
        for (int i = 1; i < network_size; ++i) {
            second_maximum_cluster_size[i] = obs_second_maximum_cluster_size[i] / (double)ensemble_size;
        }
        second_maximum_cluster_size[0] = 1.0;

        CSV::write(directory + name.get_NGE(), second_maximum_cluster_size, precision);
    }
    //! Mean Cluster Size
    if (check_list.at("mean_cluster_size")) {
        const std::string directory = data_dir + "mean_cluster_size/";
        CSV::generateDirectory(directory);
        std::vector<double> mean_cluster_size = obs_mean_cluster_size / (double)ensemble_size;
        mean_cluster_size[0] = 1.0;

        CSV::write(directory + name.get_NGE(), mean_cluster_size, precision);
    }
    //! Second Moment
    if (check_list.at("second_moment")) {
        const std::string directory = data_dir + "second_moment/";
        CSV::generateDirectory(directory);
        std::vector<double> second_moment = obs_second_moment / (double)ensemble_size;
        second_moment[0] = 1.0 / pow(network_size, 2.0);

        CSV::write(directory + name.get_NGE(), second_moment, precision);
    }

    //! Inter Event Time
    if (check_list.at("inter_event_time")) {
        const std::string directory = data_dir + "inter_event_time/";
        CSV::generateDirectory(directory);
        std::map<int, double> inter_event_time;
        for (int t = 0; t < network_size; ++t) {
            if (obs_inter_event_time[t].second) {
                inter_event_time[t] = obs_inter_event_time[t].first / (double)obs_inter_event_time[t].second;
            }
        }

        CSV::write(directory + name.get_NGE(), inter_event_time, precision);
    }

    //! Max Delta Order Parameter
    if (check_list.at("max_delta_order_parameter")) {
        const std::string directory = data_dir + "max_delta_order_parameter/";
        CSV::generateDirectory(directory);

        CSV::write(directory + name.get_NGE(), obs_max_delta_order_parameter, precision);
    }
}

}  // namespace mBFW