#pragma once

#include <string>

#include "Common.hpp"
#include "StringFormat.hpp"

namespace mBFW {
struct Name {
    //* Member variables
   public:
    int network_size;
    double acceptance_threshold;
    unsigned ensemble_size;
    int core_num;

    std::string NG;
    std::string NGE;
    std::string core_prefix;

    //* Member methods
   public:
    Name() {}
    Name(const int&, const double&, const unsigned&, const int&);

    const std::string get_NG() const;
    const std::string get_NGE() const;
    const std::string get_NGES(const std::string&, const double&) const;
};

Name::Name(const int& t_network_size,
           const double& t_acceptance_threshold,
           const unsigned& t_ensemble_size,
           const int& t_core_num) : network_size(t_network_size),
                                   acceptance_threshold(t_acceptance_threshold),
                                   ensemble_size(t_ensemble_size),
                                   core_num(t_core_num) {
    NG = "N" + to_stringWithExponent(t_network_size, 1) + ",G" + to_stringWithPrecision(t_acceptance_threshold, 1);
    NGE = NG + ",E" + std::to_string(ensemble_size);
    core_prefix = (t_core_num == -1) ? "" : "-" + std::to_string(t_core_num);
}

const std::string Name::get_NG() const {
    return NG + core_prefix + data_prefix;
}

const std::string Name::get_NGE() const {
    return NGE + core_prefix + data_prefix;
}

const std::string Name::get_NGES(const std::string& t_standard, const double& t_repeater) const {
    return NGE + "," + t_standard + to_stringWithPrecision(t_repeater) + core_prefix + data_prefix;
}

}  // namespace mBFW