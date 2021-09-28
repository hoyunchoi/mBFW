#pragma once

#include <array>
#include <set>
#include <string>
#include <utility>

namespace mBFW {
//* Path
const std::string data_dir = "/pds/pds11/hoyun/mBFW/data/";
const std::string log_dir = "/pds/pds11/hoyun/mBFW/log/";
const std::string data_prefix = ".dat";

//* Generator parameter
const std::set<unsigned> network_size_set = {10000, 20000, 40000, 80000, 160000, 320000, 640000, 1280000, 2560000, 5120000, 10240000};
const std::set<double> acceptance_threshold_set = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

//* Checklist
constexpr std::array<std::pair<std::string_view, bool>, 6> observables{{{"order_parameter", true},
                                                                        {"second_maximum", true},
                                                                        {"second_moment", true},
                                                                        {"mean_cluster_size", true},
                                                                        {"inter_event_time", true},
                                                                        {"max_delta_order_parameter", true}}};

const bool check_network_size(const unsigned& t_network_size) {
    /*
        Return true only if input network size is valid
    */
    if (network_size_set.find(t_network_size) == network_size_set.end()) {
        return false;
    }
    return true;
}

const bool check_acceptance_threshold(const double& t_acceptance_threshold) {
    /*
        Return true only if input acceptanceThreshold is valid
    */
    if (acceptance_threshold_set.find(t_acceptance_threshold) == acceptance_threshold_set.end()) {
        return false;
    }
    return true;
}

}  // namespace mBFW