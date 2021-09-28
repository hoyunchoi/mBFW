#pragma once

#include <fstream>
#include <string>

#include "Common.hpp"

namespace mBFW {

struct Logger {
   protected:
    //* Member variables
    std::string _log_path;

   public:
    //* Member methods
    Logger() {}
    Logger(const std::string&);

    void log(const std::string&);
};

Logger::Logger(const std::string& t_name) {
    /*
        Define logger with name
    */
    _log_path = log_dir + t_name + ".log";
}

void Logger::log(const std::string& t_message) {
    /*
        Print input message to log file
        All log should be 'append'
    */
    std::ofstream logFile(_log_path, std::ios_base::app);
    logFile << t_message << "\n";
    logFile.close();
}

}  // namespace mBFW