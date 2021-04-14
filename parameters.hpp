#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "../library/linearAlgebra.hpp"
#include "common.hpp"

namespace mBFW {
struct Parameter {
   protected:
    //* Member varialbes
    int m_networkSize;
    double m_acceptanceThreshold;
    std::vector<double> m_clusterSizeDist_orderParameter;
    std::vector<double> m_clusterSizeDist_time;
    std::vector<double> m_orderParameterDist_time;
    std::map<std::string, double> m_points;

   public:
    //* Member functions
    Parameter() {}
    Parameter(const int&, const double&);

    const std::map<std::string, int> get_points() const;
    const std::set<int> get_clusterSizeDist_orderParameter() const;
    const std::set<int> get_clusterSizeDist_time() const;
    const std::set<int> get_orderParameterDist_time() const;

   protected:
    void m_set_points();
    void m_set_clusterSizeDist_orderParameter();
    void m_set_clusterSizeDist_time();
    void m_set_orderParameterDist_time();
    const std::set<int> m_doubleVec2intSet(const std::vector<double>&) const;
};

Parameter::Parameter(const int& t_networkSize, const double& t_acceptanceThreshold) : m_networkSize(t_networkSize), m_acceptanceThreshold(t_acceptanceThreshold) {
    m_set_points();
    m_set_clusterSizeDist_orderParameter();
    m_set_clusterSizeDist_time();
    m_set_orderParameterDist_time();
}

void Parameter::m_set_points() {
    //* Read file
    const std::string pointFileName = mBFW::dataDirectory + "points/" + mBFW::fileName::NG(m_networkSize, m_acceptanceThreshold);
    std::ifstream readFile(pointFileName);
    std::string line;
    while (getline(readFile, line)) {
        for (const std::string& type : mBFW::pointTypes){
            if (line.find(type) != line.npos){
                m_points[type] = std::stod(line.substr(line.find(": ") + 2));
            }
        }
    }

    //* Check if every points are set
    if (m_points.size() != mBFW::pointTypes.size()) {
        std::ofstream ERROR("ERROR.log", std::ios_base::app);
        ERROR << "Problem at reading " << pointFileName << "\n";
        ERROR.close();
        exit(1);
    }
    return;
}

void Parameter::m_set_clusterSizeDist_orderParameter() {
    std::vector<double> extra;
    //* Default values
    m_clusterSizeDist_orderParameter = m_networkSize >= 1280000 ? linearAlgebra::elementPow(10.0, linearAlgebra::linspace(-6.0, -4.0, 100)) : linearAlgebra::elementPow(10.0, linearAlgebra::linspace(-4.0, -2.0, 100));
    extra = linearAlgebra::arange(0.01, 0.99, 0.01);
    m_clusterSizeDist_orderParameter.insert(m_clusterSizeDist_orderParameter.end(), extra.begin(), extra.end());

    //* Additional values
    if (m_acceptanceThreshold == 0.2) {
        extra = {0.971, 0.972, 0.973, 0.974, 0.975, 0.976, 0.977, 0.978, 0.979, 0.981, 0.982, 0.983, 0.984, 0.985, 0.986, 0.987, 0.988, 0.989};
        if (m_networkSize >= 2560000) {
            extra = {0.971, 0.972, 0.973, 0.974, 0.975, 0.976, 0.977};
        }
    } else if (m_acceptanceThreshold == 0.3) {
        extra = {0.928, 0.929, 0.931, 0.932, 0.933, 0.934, 0.935, 0.936, 0.937, 0.938, 0.939, 0.941, 0.942, 0.943, 0.944, 0.945, 0.946, 0.947, 0.948, 0.949};
        if (m_networkSize >= 2560000) {
            extra = {0.928, 0.929, 0.931, 0.932, 0.933, 0.934, 0.935};
        }
    } else if (m_acceptanceThreshold == 0.4) {
        extra = {0.877, 0.878, 0.879, 0.881, 0.882, 0.883, 0.884, 0.885, 0.886, 0.887, 0.888, 0.889, 0.891, 0.892, 0.893, 0.894, 0.995, 0.896, 0.897, 0.898, 0.899};
        if (m_networkSize >= 2560000) {
            extra = {0.877, 0.878, 0.879, 0.881, 0.882, 0.883, 0.884};
        }
    } else if (m_acceptanceThreshold == 0.5) {
        extra = {0.811, 0.812, 0.813, 0.814, 0.815, 0.816, 0.817, 0.818, 0.819, 0.821, 0.822, 0.823, 0.824, 0.825, 0.826, 0.827, 0.828, 0.829, 0.831, 0.832, 0.833, 0.834, 0.835, 0.836, 0.837, 0.838, 0.839};
        if (m_networkSize >= 2560000) {
            extra = {0.811, 0.812, 0.813, 0.814, 0.815, 0.816, 0.817, 0.818, 0.819};
        }
    } else if (m_acceptanceThreshold == 0.6) {
        extra = {0.739, 0.741, 0.742, 0.743, 0.744, 0.745, 0.746, 0.747, 0.748, 0.749, 0.751, 0.752, 0.753, 0.754, 0.755, 0.756, 0.757, 0.758, 0.759, 0.761, 0.762, 0.763, 0.764, 0.765};
        if (m_networkSize >= 2560000) {
            extra = {0.739, 0.741, 0.742, 0.743, 0.744, 0.745, 0.746, 0.747, 0.748};
        }
    } else if (m_acceptanceThreshold == 0.7) {
        extra = {0.649, 0.651, 0.652, 0.653, 0.654, 0.655, 0.656, 0.657, 0.658, 0.659, 0.661, 0.662, 0.663, 0.664, 0.665, 0.667, 0.668, 0.669, 0.671, 0.672, 0.673, 0.674, 0.675};
        if (m_networkSize >= 2560000) {
            extra = {0.649, 0.651, 0.652, 0.653, 0.654, 0.655, 0.656, 0.657, 0.658, 0.659};
        }
    } else if (m_acceptanceThreshold == 0.8) {
        extra = {0.545, 0.546, 0.547, 0.548, 0.549, 0.551, 0.552, 0.553, 0.554, 0.555, 0.556, 0.557, 0.558, 0.559, 0.561, 0.562, 0.563, 0.564, 0.565, 0.566, 0.567, 0.568, 0.569};
        if (m_networkSize >= 2560000) {
            extra = {0.545, 0.546, 0.547, 0.548, 0.549, 0.551, 0.552, 0.553, 0.554, 0.555};
        }
    } else if (m_acceptanceThreshold == 0.9) {
        extra = {0.401, 0.402, 0.403, 0.404, 0.405, 0.406, 0.407, 0.408, 0.409, 0.411, 0.412, 0.413, 0.414, 0.415, 0.416, 0.417, 0.418, 0.419, 0.421, 0.422, 0.423, 0.424, 0.425, 0.426, 0.427, 0.428, 0.429};
        if (m_networkSize >= 2560000) {
            extra = {0.401, 0.402, 0.403, 0.404, 0.405, 0.406, 0.407, 0.408, 0.409, 0.411, 0.412, 0.413};
        }
    }
    m_clusterSizeDist_orderParameter.insert(m_clusterSizeDist_orderParameter.end(), extra.begin(), extra.end());

    return;
}

void Parameter::m_set_clusterSizeDist_time() {
    if (m_acceptanceThreshold == 0.2) {
        m_clusterSizeDist_time = {0.98, 0.985, 0.99, 0.991, 0.992, 0.993, 0.994, 0.995, 0.996, 0.997, 0.998, 0.999};
        m_clusterSizeDist_time = {0.9971, 0.9972, 0.9973, 0.9974, 0.9975, 0.9976, 0.9977, 0.9978, 0.9979};
    } else if (m_acceptanceThreshold == 0.3) {
        m_clusterSizeDist_time = {0.95, 0.96, 0.97, 0.98, 0.981, 0.982, 0.983, 0.984, 0.985, 0.986, 0.987, 0.989, 0.99, 0.995};
    } else if (m_acceptanceThreshold == 0.4) {
        m_clusterSizeDist_time = {0.92, 0.94, 0.96, 0.961, 0.962, 0.963, 0.964, 0.965, 0.966, 0.967, 0.968, 0.969, 0.97, 0.98, 0.99};
    } else if (m_acceptanceThreshold == 0.5) {
        m_clusterSizeDist_time = {0.9, 0.91, 0.92, 0.93, 0.931, 0.932, 0.933, 0.934, 0.935, 0.936, 0.937, 0.938, 0.939, 0.94, 0.941, 0.942, 0.943, 0.944, 0.945, 0.946, 0.947, 0.948, 0.949, 0.95, 0.96};
    } else if (m_acceptanceThreshold == 0.6) {
        m_clusterSizeDist_time = {0.84, 0.86, 0.88, 0.89, 0.896, 0.898, 0.90, 0.902, 0.904, 0.905, 0.906, 0.908, 0.91, 0.92, 0.94};
    } else if (m_acceptanceThreshold == 0.7) {
        m_clusterSizeDist_time = {0.8, 0.82, 0.83, 0.84, 0.846, 0.848, 0.85, 0.852, 0.854, 0.856, 0.858, 0.86, 0.87};
    } else if (m_acceptanceThreshold == 0.8) {
        m_clusterSizeDist_time = {0.74, 0.76, 0.78, 0.79, 0.792, 0.794, 0.796, 0.798, 0.80, 0.802, 0.804, 0.806, 0.808, 0.81, 0.83};
    } else if (m_acceptanceThreshold == 0.9) {
        m_clusterSizeDist_time = {0.66, 0.68, 0.69, 0.7, 0.71, 0.712, 0.714, 0.716, 0.718, 0.72, 0.722, 0.724, 0.726, 0.728, 0.73, 0.74};
    } else if (m_acceptanceThreshold == 1.0) {
        m_clusterSizeDist_time = {0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55};
    }

    return;
}

void Parameter::m_set_orderParameterDist_time() {
    std::vector<double> extra;

    //* Default values: near t_a
    m_orderParameterDist_time = linearAlgebra::arange(m_points.at("t_a1") - 0.01, m_points.at("t_a1") + 0.01, 0.0001);
    m_orderParameterDist_time.emplace_back(m_points.at("t_a2"));
    m_orderParameterDist_time.emplace_back(m_points.at("t_b"));
    m_orderParameterDist_time.emplace_back(m_points.at("t_c"));

    //* Extra values
    if (m_acceptanceThreshold == 0.2) {
        extra = {0.990, 0.991, 0.992, 0.993, 0.994, 0.995, 0.996, 0.997, 0.998, 0.999};
    } else if (m_acceptanceThreshold == 0.3) {
        extra = {0.975, 0.976, 0.977, 0.978, 0.979, 0.980, 0.981, 0.982, 0.983, 0.984, 0.985, 0.986, 0.987, 0.988};
    } else if (m_acceptanceThreshold == 0.4) {
        extra = {0.955, 0.956, 0.957, 0.958, 0.959, 0.960, 0.961, 0.962, 0.963, 0.964, 0.965, 0.966, 0.967, 0.968, 0.969, 0.970};
    } else if (m_acceptanceThreshold == 0.5) {
        extra = {0.925, 0.926, 0.927, 0.928, 0.929, 0.930, 0.931, 0.932, 0.933, 0.934, 0.935, 0.936, 0.937, 0.938, 0.939, 0.940, 0.941, 0.942};
    } else if (m_acceptanceThreshold == 0.6) {
        extra = {0.888, 0.889, 0.890, 0.891, 0.892, 0.893, 0.894, 0.895, 0.896, 0.897, 0.898, 0.899, 0.900, 0.901, 0.902, 0.903, 0.904, 0.905, 0.906, 0.907};
    } else if (m_acceptanceThreshold == 0.7) {
        extra = {0.840, 0.841, 0.842, 0.843, 0.844, 0.845, 0.846, 0.847, 0.848, 0.849, 0.850, 0.851, 0.852, 0.853, 0.854, 0.855, 0.856, 0.857, 0.858, 0.859, 0.860};
    } else if (m_acceptanceThreshold == 0.8) {
        extra = {0.778, 0.779, 0.780, 0.781, 0.782, 0.783, 0.784, 0.785, 0.786, 0.787, 0.788, 0.789, 0.790, 0.791, 0.792, 0.793, 0.794, 0.795, 0.796, 0.797, 0.798, 0.799, 0.800, 0.801, 0.802, 0.803};
    } else if (m_acceptanceThreshold == 0.9) {
        extra = {0.689, 0.690, 0.691, 0.692, 0.693, 0.694, 0.695, 0.696, 0.697, 0.698, 0.699, 0.700, 0.701, 0.702, 0.703, 0.704, 0.705, 0.706, 0.707, 0.708, 0.709, 0.710, 0.711, 0.712, 0.713, 0.714, 0.715, 0.716, 0.717, 0.718, 0.719, 0.720, 0.721, 0.722};
    }
    m_orderParameterDist_time.insert(m_orderParameterDist_time.end(), extra.begin(), extra.end());

    return;
}

const std::set<int> Parameter::m_doubleVec2intSet(const std::vector<double>& t_doubleVec) const {
    std::set<int> intSet;
    for (const double& e : t_doubleVec) {
        intSet.emplace_hint(intSet.end(), (int)(e * m_networkSize));
    }
    return intSet;
}

const std::map<std::string, int> Parameter::get_points() const {
    std::map<std::string, int> result;
    for (const std::pair<std::string, double>& e : m_points) {
        result[e.first] = (int)(e.second * m_networkSize);
    }
    return result;
}

const std::set<int> Parameter::get_clusterSizeDist_orderParameter() const {
    return m_doubleVec2intSet(m_clusterSizeDist_orderParameter);
}

const std::set<int> Parameter::get_clusterSizeDist_time() const {
    return m_doubleVec2intSet(m_clusterSizeDist_time);
}

const std::set<int> Parameter::get_orderParameterDist_time() const {
    return m_doubleVec2intSet(m_orderParameterDist_time);
}

}  // namespace mBFW
