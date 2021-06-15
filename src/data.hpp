#pragma once

#include <filesystem>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "CSV.hpp"
#include "common.hpp"
#include "linearAlgebra.hpp"

namespace mBFW {
struct Data {
  protected:
    //* Member variables
    int m_networkSize;
    double m_acceptanceThreshold;
    bool m_deletion;
    bool m_writeLog;
    std::string m_target;

  public:
    //* Member functions
    Data() {}
    Data(const int&, const double&, const bool, const bool);
    void run(const std::map<std::string, bool>&);

    template <typename T>
    void continuousAverage(const std::string&, const T&) const;
    void continuousAverage_repeater(const std::string&) const;

    template <typename T, typename TT>
    void discreteAverage(const std::string&, const std::map<T, TT>&) const;

    void netObservable(const std::string&) const;

  protected:
    //* Directory and fild handling
    const std::string m_defineAdditionalDirectory(const std::string&, const std::string&) const;
    const std::set<std::string> m_findTargetFileNameList(const std::string&, const std::string&) const;
    const unsigned m_extractEnsemble(const std::string&) const;
    void m_conditionallyDeleteFile(const std::string&) const;

    //* Data handling
    template <typename T>
    std::tuple<T, unsigned> m_continuousAvgFile(const std::string&, const std::set<std::string>&, const T&) const;

    template <typename T, typename TT>
    std::tuple<std::map<T, TT>, unsigned> m_discreteAvgFile(const std::string&, const std::set<std::string>&, const std::map<T, TT>&) const;
};

Data::Data(const int& t_networkSize,
           const double& t_acceptanceThreshold,
           const bool t_deletion,
           const bool t_writeLog)
    : m_networkSize(t_networkSize),
      m_acceptanceThreshold(t_acceptanceThreshold),
      m_deletion(t_deletion),
      m_writeLog(t_writeLog) {
    m_target = fileName::base(t_networkSize, t_acceptanceThreshold);
}

void Data::run(const std::map<std::string, bool>& t_checkList) {
    if (t_checkList.at("ageDist")) {
        for (const std::string& state : mBFW::states) {
            continuousAverage("ageDist/" + state, std::vector<double>{});
        }
    }
    if (t_checkList.at("ageDist_time")) {
        const std::vector<std::string> states = m_acceptanceThreshold == 1.0 ? std::vector<std::string>{"0_A1", "C_1"} : mBFW::states;
        for (const std::string& state : states) {
            continuousAverage("ageDist_time/" + state, std::vector<double>{});
        }
    }
    if (t_checkList.at("clusterSizeDist")) {
        continuousAverage_repeater("clusterSizeDist");
    }
    if (t_checkList.at("clusterSizeDist_time")) {
        continuousAverage_repeater("clusterSizeDist_time");
    }
    if (t_checkList.at("deltaUpperBoundDist")) {
        for (const std::string& state : mBFW::states) {
            continuousAverage("deltaUpperBoundDist/" + state, std::map<int, double>{});
        }
    }
    if (t_checkList.at("deltaUpperBoundDist_time")) {
        for (const std::string& state : mBFW::states) {
            continuousAverage("deltaUpperBoundDist_time/" + state, std::map<int, double>{});
        }
    }
    if (t_checkList.at("interEventTime")) {
        discreteAverage("interEventTime", std::map<int, double>{});
    }
    if (t_checkList.at("interEventTime_orderParameter")) {
        discreteAverage("interEventTime_orderParameter", std::map<int, double>{});
    }
    if (t_checkList.at("interEventTimeDist")) {
        for (const std::string& state : mBFW::states) {
            continuousAverage("interEventTimeDist/" + state, std::map<int, double>{});
        }
    }
    if (t_checkList.at("interEventTimeDist_time")) {
        for (const std::string& state : mBFW::states) {
            continuousAverage("interEventTimeDist_time/" + state, std::map<int, double>{});
        }
    }
    if (t_checkList.at("meanClusterSize")) {
        continuousAverage("meanClusterSize", std::vector<double>{});
    }
    if (t_checkList.at("netOrderParameter")) {
        for (const std::string& state : std::set<std::string>{"sub", "super"}) {
            netObservable("netOrderParameter_mc78/" + state);
        }
    }
    if (t_checkList.at("netSecondMoment")) {
        for (const std::string& state : std::set<std::string>{"sub", "super"}) {
            netObservable("netSecondMoment_mc78/" + state);
        }
    }
    if (t_checkList.at("orderParameter")) {
        continuousAverage("orderParameter", std::vector<double>{});
    }
    if (t_checkList.at("orderParameter_interEventTime")) {
        discreteAverage("orderParameter_interEventTime", std::map<int, double>{});
    }
    if (t_checkList.at("orderParameterDist")) {
        continuousAverage_repeater("orderParameterDist");
    }
    if (t_checkList.at("secondMaximum")) {
        continuousAverage("secondMaximum", std::vector<double>{});
    }
    if (t_checkList.at("secondMoment")) {
        continuousAverage("secondMoment", std::vector<double>{});
    }
}

const std::string Data::m_defineAdditionalDirectory(const std::string& t_baseDirectory,
                                                    const std::string& t_additionalDirectoryName) const {
    const std::string averageDirectory = t_baseDirectory + t_additionalDirectoryName + "/";
    CSV::generateDirectory(averageDirectory);
    return averageDirectory;
}

const std::set<std::string> Data::m_findTargetFileNameList(const std::string& t_directory,
                                                           const std::string& t_target) const {
    namespace fs = std::filesystem;
    std::set<std::string> targetFileNameList;
    for (const auto& file : fs::directory_iterator(t_directory)) {
        const std::string fileName = file.path().filename();
        if (fileName.find(t_target) != fileName.npos) {
            targetFileNameList.emplace(fileName);
        }
    }
    return targetFileNameList;
}

const unsigned Data::m_extractEnsemble(const std::string& t_fileName) const {
    std::string temp = t_fileName.substr(t_fileName.find("E") + 1);
    temp = temp.substr(0, temp.find_first_of(",-"));
    return std::stoul(temp);
}

void Data::m_conditionallyDeleteFile(const std::string& t_deletionFile) const {
    if (m_deletion) {
        if (m_writeLog) {
            std::cout << "Deleting file " << t_deletionFile << "\n";
        }
        CSV::deleteFile(t_deletionFile);
    }
}

template <typename T>
std::tuple<T, unsigned> Data::m_continuousAvgFile(const std::string& t_directory,
                                                  const std::set<std::string>& t_fileNameList, const T& t_format) const {
    using namespace linearAlgebra;
    unsigned totalEnsemble = 0;
    for (const std::string& fileName : t_fileNameList) {
        totalEnsemble += m_extractEnsemble(fileName);
    }
    T average;
    for (const std::string& fileName : t_fileNameList) {
        T temp;
        CSV::read(t_directory + fileName, temp);
        const double ratio = m_extractEnsemble(fileName) / (double)totalEnsemble;
        average += temp * ratio;
    }
    return std::make_tuple(average, totalEnsemble);
}

//? T=std::pair<int,unsigned>, TT=double
template <typename T, typename TT>
std::tuple<std::map<T, TT>, unsigned> Data::m_discreteAvgFile(const std::string& t_directory,
                                                              const std::set<std::string>& t_fileNameList,
                                                              const std::map<T, TT>& t_format) const {
    unsigned totalEnsemble = 0;
    std::map<T, unsigned> totalEnsembleMap;
    std::map<std::string, std::map<T, TT>> totalData;

    for (const std::string& fileName : t_fileNameList) {
        const unsigned ensemble = m_extractEnsemble(fileName);
        totalEnsemble += ensemble;
        std::map<T, TT> temp;
        CSV::read(t_directory + fileName, temp);
        totalData[fileName] = temp;
        for (const std::pair<T, TT>& e : temp) {
            totalEnsembleMap[e.first] += ensemble;
        }
    }

    std::map<T, TT> average;
    for (const std::string& fileName : t_fileNameList) {
        const unsigned ensemble = m_extractEnsemble(fileName);
        for (const std::pair<T, TT>& e : totalData.at(fileName)) {
            average[e.first] += e.second * (double)ensemble / totalEnsembleMap.at(e.first);
        }
    }
    return std::make_tuple(average, totalEnsemble);
}

template <typename T>
void Data::continuousAverage(const std::string& t_type,
                             const T& t_format) const {
    //* Define Directories and get target files
    const std::string directory = m_defineAdditionalDirectory(mBFW::dataDirectory, t_type);
    const std::set<std::string> fileNameList = m_findTargetFileNameList(directory, m_target);

    //* Check the number of files
    if (fileNameList.empty()) {
        std::ofstream ERROR(mBFW::logDirectory + "/ERROR.log", std::ios_base::app);
        ERROR << m_target << ": No file at " << directory << "\n";
        ERROR.close();
        exit(1);
    } else if (fileNameList.size() == 1) {
        if (m_writeLog) {
            std::cout << "Passing file " << directory + *fileNameList.begin() << "\n";
        }
        return;
    }

    //* Average raw data
    const auto [average, totalEnsemble] = m_continuousAvgFile(directory, fileNameList, t_format);

    //* Write the file
    const std::string newFileName = fileName::NGE(m_networkSize, m_acceptanceThreshold, totalEnsemble, 0);
    if (m_writeLog) {
        std::cout << "Writing file " << directory + newFileName << "\n";
    }
    CSV::write(directory + newFileName, average, -1);

    //* Delete previous averaged and trimmed data after successfully writing
    for (const std::string& fileName : fileNameList) {
        if (fileName != newFileName) {
            m_conditionallyDeleteFile(directory + fileName);
        }
    }

    //* void return
    return;
}

void Data::continuousAverage_repeater(const std::string& t_type) const {
    //* Define Directories and get target files
    const std::string directory = m_defineAdditionalDirectory(mBFW::dataDirectory, t_type);
    std::set<std::string> fileNameList = m_findTargetFileNameList(directory, m_target);

    //* Check the number of files
    if (fileNameList.empty()) {
        std::ofstream ERROR(mBFW::logDirectory + "ERROR.log", std::ios_base::app);
        ERROR << m_target << ": No file at " << directory << "\n";
        ERROR.close();
        exit(1);
    } else if (fileNameList.size() == 1) {
        if (m_writeLog) {
            std::cout << "Passing file " << directory + *fileNameList.begin() << "\n";
        }
        return;
    }

    //* Check the standard of input type
    const std::string standard = (t_type == "orderParameterDist" || t_type == "clusterSizeDist_time") ? "T" : "OP";

    //* Get total ensemble size for each repeater value
    std::map<std::string, std::vector<std::vector<double>>> totalData;
    std::map<int, unsigned> totalEnsembleMap;
    for (const std::string& fileName : fileNameList) {
        std::vector<std::vector<double>> temp;
        CSV::read(directory + fileName, temp);
        totalData[fileName] = temp;
        for (const std::vector<double>& data : temp) {
            totalEnsembleMap[(int)data[0]] += (unsigned)data[1];
        }
    }

    //* Average for each repeater value
    std::map<int, std::map<int, double>> averageMap;
    for (const std::string& fileName : fileNameList) {
        for (const std::vector<double>& data : totalData.at(fileName)) {
            const int repeater = (int)data[0];
            const double ratio = data[1] / (double)totalEnsembleMap.at(repeater);
            for (unsigned i = 2; i < data.size(); i += 2) {
                averageMap[repeater][(int)data[i]] += data[i + 1] * ratio;
            }
        }
    }

    //* Change averageMap to averageVector
    std::vector<std::vector<double>> average;
    average.reserve(averageMap.size());
    for (const std::pair<int, std::map<int, double>>& avg : averageMap) {
        std::vector<double> currentData;
        currentData.reserve(2 + 2 * avg.second.size());
        currentData.emplace_back((double)avg.first);
        currentData.emplace_back((double)totalEnsembleMap.at(avg.first));
        for (const std::pair<int, double>& e : avg.second) {
            currentData.emplace_back((double)e.first);
            currentData.emplace_back(e.second);
        }
        average.emplace_back(currentData);
    }

    //* Write the file
    const std::string newFileName = fileName::NG(m_networkSize, m_acceptanceThreshold, 0);
    if (m_writeLog) {
        std::cout << "Writing file " << directory + newFileName << "\n";
    }
    CSV::write(directory + newFileName, average, -1);

    //* Delete previous averaged and trimmed data after successfully writing
    for (const std::string& fileName : fileNameList) {
        if (fileName != newFileName) {
            m_conditionallyDeleteFile(directory + fileName);
        }
    }

    //* Get single file directory and file lists
    const std::string singleDirectory = m_defineAdditionalDirectory(directory, "single");
    std::set<std::string> singleFileNameList = m_findTargetFileNameList(singleDirectory, m_target);

    //* Seperate to multiple files containing single repeater
    std::set<std::string> newSingleFileNameList;
    for (const std::pair<int, std::map<int, double>>& single : averageMap) {
        const int repeater = single.first;
        const unsigned totalEnsemble = totalEnsembleMap.at(repeater);
        const std::string newSingleFileName = fileName::NGES(m_networkSize,
                                                             m_acceptanceThreshold,
                                                             totalEnsemble,
                                                             standard,
                                                             repeater / (double)m_networkSize, 0);
        newSingleFileNameList.emplace(newSingleFileName);

        if (singleFileNameList.find(newSingleFileName) != singleFileNameList.end()) {
            if (m_writeLog) {
                std::cout << "Passing file " << singleDirectory + newSingleFileName << "\n";
            }
            continue;
        } else {
            if (m_writeLog) {
                std::cout << "Writing file " << singleDirectory + newSingleFileName << "\n";
            }
            CSV::write(singleDirectory + newSingleFileName, single.second, -1);
        }
    }

    //* Delete old single files if they are updated
    for (const std::string& singleFileName : singleFileNameList) {
        if (newSingleFileNameList.find(singleFileName) == newSingleFileNameList.end()) {
            m_conditionallyDeleteFile(singleDirectory + singleFileName);
        }
    }

    //* void return
    return;
}

template <typename T, typename TT>
void Data::discreteAverage(const std::string& t_type,
                           const std::map<T, TT>& t_format) const {
    //* Define directories and target files
    const std::string directory = m_defineAdditionalDirectory(mBFW::dataDirectory, t_type);
    const std::set<std::string> fileNameList = m_findTargetFileNameList(directory, m_target);

    //* Check the number of files
    if (fileNameList.empty()) {
        std::ofstream ERROR(mBFW::logDirectory + "ERROR.log", std::ios_base::app);
        ERROR << m_target << ": No file at " << directory << "\n";
        ERROR.close();
        exit(1);
    } else if (fileNameList.size() == 1) {
        if (m_writeLog) {
            std::cout << "Passing file " << directory + *fileNameList.begin() << "\n";
        }
        return;
    }

    //* Average raw data
    const auto [average, totalEnsembleSize] = m_discreteAvgFile(directory, fileNameList, t_format);

    //* Write the file
    const std::string newFileName = fileName::NGE(m_networkSize, m_acceptanceThreshold, totalEnsembleSize, 0);
    if (m_writeLog) {
        std::cout << "Writing file " << directory + newFileName << "\n";
    }
    CSV::write(directory + newFileName, average, -1);

    //* Delete previous averaged and trimmed data after successfully writing
    for (const std::string& fileName : fileNameList) {
        if (fileName != newFileName) {
            m_conditionallyDeleteFile(directory + fileName);
        }
    }

    //* void return
    return;
}

void Data::netObservable(const std::string& t_type) const {
    //* Define directories and target files
    const std::string directory = m_defineAdditionalDirectory(mBFW::dataDirectory, t_type);
    const std::set<std::string> fileNameList = m_findTargetFileNameList(directory, m_target);

    //* Check the number of files
    if (fileNameList.empty()) {
        std::ofstream ERROR(mBFW::logDirectory + "ERROR.log", std::ios_base::app);
        ERROR << m_target << ": No file at " << directory << "\n";
        ERROR.close();
        exit(1);
    } else if (fileNameList.size() == 1) {
        if (m_writeLog) {
            std::cout << "Passing file " << directory + *fileNameList.begin() << "\n";
        }
        return;
    }

    //* Get total ensemble size for each time
    std::map<std::string, std::vector<std::vector<double>>> totalData;
    std::map<int, unsigned> totalEnsembleMap;
    for (const std::string& fileName : fileNameList) {
        std::vector<std::vector<double>> temp;
        CSV::read(directory + fileName, temp);
        totalData[fileName] = temp;
        for (const std::vector<double>& data : temp) {
            totalEnsembleMap[(int)data[0]] += (unsigned)data[1];
        }
    }

    //* Average for each time
    std::map<std::pair<int, unsigned>, double> average;
    for (const std::string& fileName : fileNameList) {
        for (const std::vector<double>& data : totalData.at(fileName)) {
            const int time = (int)data[0];
            const double ratio = data[1] / (double)totalEnsembleMap.at(time);
            average[std::pair<int, unsigned>{time, totalEnsembleMap.at(time)}] += data[2] * ratio;
        }
    }

    //* Write the file
    const std::string newFileName = fileName::NG(m_networkSize, m_acceptanceThreshold, 0);
    if (m_writeLog) {
        std::cout << "Writing file " << directory + newFileName << "\n";
    }
    CSV::write(directory + newFileName, average, -1);

    //* Delete previous averaged and trimmed data after successfully writing
    for (const std::string& fileName : fileNameList) {
        if (fileName != newFileName) {
            m_conditionallyDeleteFile(directory + fileName);
        }
    }
    return;
}
} // namespace mBFW
