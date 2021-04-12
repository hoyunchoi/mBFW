from dataProcess import *
import pandas as pd
import glob
import numpy as np
import sys
sys.path.append("../library/")


dataDirectory = "../data/mBFW_v3/"
states = ["0A", "AG", "GC", "C1"]
observables = set()

#* Continuous observables
observables.add("orderParameter")
observables.add("meanClusterSize")
observables.add("orderParameterVariance")

observables.add("orderParameter_trial")
observables.add("meanClusterSize_trial")
observables.add("orderParameterVariance_trial")

#* Discrete observables
observables.add("interEventTime")
observables.add("dotOrderParameter")
observables.add("noRestriction")

#* Basic distributions
for standard in ["", "_exact", "_time"]:
    observables.add("clusterSizeDist" + standard)
observables.add("orderParameterDist")


#* Observables distinguished by intervals
for state in states:
    observables.add("ageDist/" + state)
    observables.add("ageDist_time/" + state)
    observables.add("interEventTimeDist/" + state)
    observables.add("interEventTimeDist_time/" + state)
    observables.add("deltaUpperBoundDist/" + state)
    observables.add("deltaUpperBoundDist_time/" + state)
    # observables.add("interEventTime_deltaUpperBound/" + state)

#* 3D sampled and its derivatives
observables.add("sampled_deltaUpperBound_interEventTime")
observables.add("sampled_upperBound_interEventTime")
observables.add("sampled_time_interEventTime")
observables.add("deltaUpperBound_interEventTime_tot")
observables.add("interEventTime_deltaUpperBound_tot")
observables.add("interEventTime_upperBound_tot")
observables.add("upperBound_interEventTime_tot")
observables.add("time_interEventTime_tot")
observables.add("interEventTime_time_tot")

#* Other observables
observables.add("points")
observables.add("dynamics")
observables.add("periodDynamics")


absolutePathList = {}
for observable in observables:
    if ("clusterSizeDist" in observable) or (observable == "orderParameterDist"):
        absolutePathList[observable] = dataDirectory + observable + "/single/"
    else:
        absolutePathList[observable] = dataDirectory + observable + "/"


# * CSV Reader


def readCSV(t_fileName):
    data = pd.read_csv(t_fileName, sep=',', header=None)
    data = data.values.transpose()
    if (len(data) == 0):
        return None
    elif (len(data) == 1):
        return data[0]
    else:
        return tuple([row for row in data])


# * File Name Convections
def __NG__(t_networkSize, t_acceptanceThreshold):
    return "N{:.1e},G{:.1f}*-0.txt".format(t_networkSize, t_acceptanceThreshold)


def __NGT__(t_networkSize, t_acceptanceThreshold, t_time):
    return "N{:.1e},G{:.1f}*,T{:.4f}*-0.txt".format(t_networkSize, t_acceptanceThreshold, t_time)


def __NGOP__(t_networkSize, t_acceptanceThreshold, t_orderParameter):
    return "N{:.1e},G{:.1f}*,OP{:.4f}*-0.txt".format(t_networkSize, t_acceptanceThreshold, t_orderParameter)


# * Get the order parameter/time value in directory
def extractRepeater(t_type, t_networkSize, t_acceptanceThreshold):
    repeaterList = set()
    if (t_type == "clusterSizeDist_time" or t_type == "orderParameterDist"):
        target = "T"
    elif (t_type == "clusterSizeDist" or t_type == "clusterSizeDist_exact"):
        target = "OP"
    for file in glob.glob(absolutePathList[t_type] + __NG__(t_networkSize, t_acceptanceThreshold)):
        repeaterList.add(float(file[file.find(target) + len(target): file.find("-0.txt")]))
    return np.array(list(sorted(repeaterList)))

# * Read various points


def readPoint(t_type, t_networkSize, t_acceptanceThreshold):
    with open(absolutePathList["points"] + "N{:.1e},G{:.1f}.txt".format(t_networkSize, t_acceptanceThreshold)) as file:
        point = 0
        for line in file.readlines():
            if t_type in line:
                point = float(line[line.find(": ") + 2:])
    return point


def readPoints(t_networkSize, t_acceptanceThreshold):
    points = {}
    with open(absolutePathList["points"] + "N{:.1e},G{:.1f}.txt".format(t_networkSize, t_acceptanceThreshold)) as file:
        content = file.readlines()
        for line in content:
            if "t_a" in line:
                points["t_a"] = float(line[line.find(": ") + 2:])
            elif "m_a" in line:
                points["m_a"] = float(line[line.find(": ") + 2:])
            elif "inf_ta" in line:
                points["t_a2"] = float(line[line.find(": ") + 2:])
            elif "inf_ma" in line:
                points["m_a2"] = float(line[line.find(": ") + 2:])
            elif "inflection_t" in line:
                points["inflection_t"] = float(line[line.find(": ") + 2:])
            elif "inflection_m" in line:
                points["inflection_m"] = float(line[line.find(": ") + 2:])
            elif "t_c_csd" in line:
                points["t_c"] = float(line[line.find(": ") + 2:])
            elif "m_c_csd" in line:
                points["m_c"] = float(line[line.find(": ") + 2:])
            elif "t_g" in line:
                points["t_g"] = float(line[line.find(": ") + 2:])
            elif "m_g" in line:
                points["m_g"] = float(line[line.find(": ") + 2:])
            elif "t_c_var" in line:
                points["t_c_var"] = float(line[line.find(": ") + 2:])
            elif "m_c_var" in line:
                points["m_c_var"] = float(line[line.find(": ") + 2:])
            elif "t_c_mcs" in line:
                points["t_c_mcs"] = float(line[line.find(": ") + 2:])
            elif "m_c_mcs" in line:
                points["m_c_mcs"] = float(line[line.find(": ") + 2:])
    return points

# * Read Observables


def read(t_type, t_networkSize, t_acceptanceThreshold, t_reapeater=None):
    # * Read time-accumulated distributions
    if (t_type == "clusterSizeDist_time" or t_type == "orderParameterDist"):
        file = glob.glob(absolutePathList[t_type] + __NGT__(t_networkSize, t_acceptanceThreshold, t_reapeater))

    # * Read orderparameter-accumulated distributions
    elif (t_type == "clusterSizeDist" or t_type == "clusterSizeDist_exact"):
        file = glob.glob(absolutePathList[t_type] + __NGOP__(t_networkSize, t_acceptanceThreshold, t_reapeater))

    # * Read general data
    else:
        file = glob.glob(absolutePathList[t_type] + __NG__(t_networkSize, t_acceptanceThreshold))

    # * Check found files
    if (len(file) != 1):
        print("There is problem at reading " + t_type + " at N={:.1e}".format(t_networkSize) + ", G={:.1f}".format(t_acceptanceThreshold))
        return
    return readCSV(file[0])


def opList2tList(orderParameter, opList):
    tList = np.zeros_like(opList)
    index = 0
    for t, value in enumerate(orderParameter):
        while value > opList[index]:
            tList[index] = t
            index += 1
            if (index >= len(opList)):
                return tList


def op2t(orderParameter, op):
    for t, value in enumerate(orderParameter):
        if (value > op):
            return t


def getSubState(state):
    index = np.zeros(2, dtype=np.int8)
    for i, s in enumerate(states):
        if s[0] == state[0]:
            index[0] = i
        if s[1] == state[1]:
            index[1] = i
    return states[index[0]: index[1] + 1]


def find_inflection_ta(networkSize, orderParameter, delta=1e-4):
    time = np.arange(0.0, 1.0, 1 / networkSize)
    bin_t, bin_op = avgLinBin(time, orderParameter, delta=delta)
    inclination = (bin_op[1:] - bin_op[:-1]) / delta
    inflection_index = np.argmax(inclination)

    max_inclination, inflection_t, inflection_op = inclination[inflection_index], bin_t[inflection_index], bin_op[inflection_index]
    inflection_ta = inflection_t - (inflection_op / max_inclination)
    inflection_ma = orderParameter[int(inflection_ta * networkSize)]
    return inflection_ta, inflection_ma, inflection_t, inflection_op


if __name__ == "__main__":
    print("This is a module readData.py")
    points = {}
    for g in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        for N in [20000, 160000, 1280000, 10240000]:
            points[N, g] = readPoints(N, g)

    for g in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
        print("g={:.1f}".format(g))
        print("t_a: ", end='')
        for N in [20000, 160000, 1280000, 10240000]:
            print("{:.4f}".format(points[N, g]["t_a"]), end='\t')
        print()
        print("t_b: ", end='')
        for N in [20000, 160000, 1280000, 10240000]:
            print("{:.4f}".format(points[N, g]["t_g"]), end='\t')
        print()
        print("t_c: ", end='')
        for N in [20000, 160000, 1280000, 10240000]:
            print("{:.4f}".format(points[N, g]["t_c"]), end='\t')
        print("\n-----------------------")
