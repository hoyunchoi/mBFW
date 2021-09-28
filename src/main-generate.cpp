#include <chrono>
#include <iostream>
#include <string>

#include "CheckList.hpp"
#include "Common.hpp"
#include "Generator.hpp"
#include "Logger.hpp"
#include "Name.hpp"

int main(int argc, char* argv[]) {
    //* Get input parameters
    const int network_size = std::stoi(argv[1]);
    const double acceptance_threshold = std::stod(argv[2]);
    const unsigned ensemble_size = std::stoul(argv[3]);
    const int core_num = std::stoi(argv[4]);
    const int random_engine_seed = -1;    //* seed chosen by std::random_device()

    //* Generate Name instance
    mBFW::Name name(network_size,
                    acceptance_threshold,
                    ensemble_size,
                    core_num);

    //* Generate log files
    mBFW::Logger error_logger("error");
    mBFW::Logger time_logger("time");

    //* Check input values
    if (not mBFW::check_network_size(network_size)) {
        error_logger.log(name.NGE + ": Not valid network size");
        return -1;
    }
    if (not mBFW::check_acceptance_threshold(acceptance_threshold)) {
        error_logger.log(name.NGE + ": Not valid acceptance threshold");
        return -1;
    }
    if (core_num <= 0) {
        error_logger.log(name.NGE + ": Not valid core number " + std::to_string(core_num));
        return -1;
    }

    //* Generate generator instance
    const auto start = std::chrono::system_clock::now();
    mBFW::Generator generator(name, random_engine_seed);
    generator.run();
    generator.save();
    std::chrono::duration<double> sec = std::chrono::system_clock::now() - start;
    time_logger.log(name.NGE + ": " + std::to_string(sec.count()));

    return 0;
}
