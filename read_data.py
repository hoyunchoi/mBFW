import sys
sys.path.append("../library/")
import pandas as pd
import glob
import numpy as np
from dataProcess import *


dataDirectory = "../data/mBFW/"
states = ["0A1", "A1A2", "A2G", "GC", "C1"]
point_type = ["t_a1", "m_a1", "t_a2", "m_a2", "t_b", "m_b", "t_c", "m_c", "t_inflection", "m_inflection"]

observables = set()
observables.add("points")

observables.add("orderParameter")
observables.add("meanClusterSize")
observables.add("orderParameterVariance")
observables.add("interEventTime")

#* Observables distinguished by intervals
for state in states:
    observables.add("ageDist/" + state)
    observables.add("interEventTimeDist/" + state)
    observables.add("deltaUpperBoundDist/" + state)

observables.add("clusterSizeDist")
observables.add("orderParameterDist")

observables.add("interEventTime_orderParameter")

absolutePathList = {}
for observable in observables:
    if (observable == "clusterSizeDist") or (observable == "orderParameterDist"):
        absolutePathList[observable] = dataDirectory + observable + "/single/"
    else:
        absolutePathList[observable] = dataDirectory + observable + "/"


# * CSV Reader


def readCSV(fileName):
    data = pd.read_csv(fileName, sep=',', header=None, dtype=np.double)
    data = data.values.transpose()
    if (len(data) == 0):
        return None
    elif (len(data) == 1):
        return data[0]
    else:
        return tuple([row for row in data])


# * File Name Convections
def __NG__(networkSize, acceptanceThreshold):
    return "N{:.1e},G{:.1f}*-0.txt".format(networkSize, acceptanceThreshold)


def __NGT__(networkSize, acceptanceThreshold, t_time):
    return "N{:.1e},G{:.1f}*,T{:.4f}*-0.txt".format(networkSize, acceptanceThreshold, t_time)


def __NGOP__(networkSize, acceptanceThreshold, t_orderParameter):
    return "N{:.1e},G{:.1f}*,OP{:.4f}*-0.txt".format(networkSize, acceptanceThreshold, t_orderParameter)


# * Get the order parameter/time value in directory
def extractRepeater(type, networkSize, acceptanceThreshold):
    repeaterList = set()
    if (type == "clusterSizeDist_time" or type == "orderParameterDist"):
        target = "T"
    elif (type == "clusterSizeDist" or type == "clusterSizeDist_exact"):
        target = "OP"
    for file in glob.glob(absolutePathList[type] + __NG__(networkSize, acceptanceThreshold)):
        repeaterList.add(float(file[file.find(target) + len(target): file.find("-0.txt")]))
    return np.array(list(sorted(repeaterList)))

# * Read various points

def readPoints(networkSize, acceptanceThreshold):
    points = {}
    with open(absolutePathList["points"] + "N{:.1e},G{:.1f}.txt".format(networkSize, acceptanceThreshold), 'r') as file:
        content = file.readlines()
        for line in content:
            for type in point_type:
                if type in line:
                    points[type] = np.double(line[line.find(": ") + 2:])
                    break
    return points


# * Read Observables


def read(type, networkSize, acceptanceThreshold, t_reapeater=None):
    # * Read time-accumulated distributions
    if (type == "clusterSizeDist_time" or type == "orderParameterDist"):
        file = glob.glob(absolutePathList[type] + __NGT__(networkSize, acceptanceThreshold, t_reapeater))

    # * Read orderparameter-accumulated distributions
    elif (type == "clusterSizeDist" or type == "clusterSizeDist_exact"):
        file = glob.glob(absolutePathList[type] + __NGOP__(networkSize, acceptanceThreshold, t_reapeater))

    # * Read general data
    else:
        file = glob.glob(absolutePathList[type] + __NG__(networkSize, acceptanceThreshold))

    # * Check found files
    if (len(file) != 1):
        print("There is problem at reading " + type + " at N={:.1e}".format(networkSize) + ", G={:.1f}".format(acceptanceThreshold))
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


def get_ta_inflection(networkSize, orderParameter, delta=1e-4):
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