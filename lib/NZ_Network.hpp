#pragma once

#include <cmath>
#include <map>
#include <utility>
#include <vector>

struct NZ_Node {
   public:
    int index;
    int parent{-1};
    int birth{0};

   public:
    NZ_Node() {}
    NZ_Node(const int& t_index) : index(t_index){};
};

struct NZ_Network {
   public:
    int networkSize{0};
    unsigned long long linkSize{0};
    int maximumClusterSize{1};
    int deltaMaximumClusterSize{0};
    std::vector<std::pair<unsigned long long, int>> changedAge;  //* changedAge[0,1] : {age, size} of two merged roots

   protected:
    std::map<int, int> m_sortedCluster;  //* sortedCluster[size] : number of cluster of 'size'
    std::vector<NZ_Node> m_nodes;

   public:
    //* Constructor
    NZ_Network() {}
    NZ_Network(const int&);

    //* Get the root of input node
    const int getRoot(const int&);

    //* Merge two clusters
    //* Two input index should be different root
    void merge(const int&, const int&);
    const int getSize(const int&) const;
    const int getSecondMaximumClusterSize() const;
    const std::map<int, int> getSortedCluster(const int& t_excludeNum = 1) const;
    const double getMeanClusterSize() const;

   protected:
    void m_updateAge(const int&, const int&, const int&, const int&);
    void m_updateSortedCluster(const int&, const int&);
    void m_updateMaximumClusterSize(const int&);
};


NZ_Network::NZ_Network(const int& t_networkSize) : networkSize(t_networkSize) {
    //* Generate default nodes
    m_nodes.reserve(t_networkSize);
    for (int index = 0; index < t_networkSize; ++index) {
        NZ_Node node(index);
        m_nodes.emplace_back(node);
    }

    //* Initialize sorted cluster and age
    m_sortedCluster[1] = t_networkSize;
    changedAge.assign(2, std::pair<int, int>{0, 0});
}

const int NZ_Network::getRoot(const int& t_index) {
    if (m_nodes[t_index].parent < 0) {
        return t_index;
    }
    return m_nodes[t_index].parent = getRoot(m_nodes[t_index].parent);
}

const int NZ_Network::getSize(const int& t_root) const {
    return -1 * m_nodes[t_root].parent;
}

void NZ_Network::m_updateAge(const int& t_root1, const int& t_size1, const int& t_root2, const int& t_size2) {
    changedAge[0] = std::pair<unsigned long long, int>{linkSize - m_nodes[t_root1].birth, t_size1};
    changedAge[1] = std::pair<unsigned long long, int>{linkSize - m_nodes[t_root2].birth, t_size2};
    m_nodes[t_root1].birth = linkSize;
}

void NZ_Network::m_updateSortedCluster(const int& t_size1, const int& t_size2) {
    --m_sortedCluster[t_size1];
    // if (!m_sortedCluster[t_size1]) {
    //     m_sortedCluster.erase(t_size1);
    // }
    --m_sortedCluster[t_size2];
    // if (!m_sortedCluster[t_size2]) {
    //     m_sortedCluster.erase(t_size2);
    // }
    ++m_sortedCluster[t_size1 + t_size2];
}

void NZ_Network::m_updateMaximumClusterSize(const int& t_newSize) {
    if (maximumClusterSize < t_newSize) {
        deltaMaximumClusterSize = t_newSize - maximumClusterSize;
        maximumClusterSize = t_newSize;
    } else {
        deltaMaximumClusterSize = 0;
    }
}

void NZ_Network::merge(const int& t_root1, const int& t_root2) {
    //* Update link size
    ++linkSize;

    //* Get size of two clusters
    const int size1 = getSize(t_root1);
    const int size2 = getSize(t_root2);
    const int newSize = size1 + size2;

    //* Optional observables
    m_updateMaximumClusterSize(newSize);
    m_updateAge(t_root1, size1, t_root2, size2);
    m_updateSortedCluster(size1, size2);

    //* Update parents
    m_nodes[t_root1].parent -= size2;
    m_nodes[t_root2].parent = t_root1;
}

const int NZ_Network::getSecondMaximumClusterSize() const {
    for (auto it = m_sortedCluster.rbegin(); it != m_sortedCluster.rend(); ++it) {
        if (it->first != maximumClusterSize || it->second > 1) {
            return it->first;
        }
    }
    return 1;
}

const std::map<int, int> NZ_Network::getSortedCluster(const int& t_excludeNum) const {
    std::map<int, int> result = m_sortedCluster;

    //* Exclue maximum cluster
    --result[maximumClusterSize];
    // if (!result[maximumClusterSize]) {
    //     result.erase(maximumClusterSize);
    // }
    //* If exclude number = 2, exclude second maximum cluster
    if (t_excludeNum == 2) {
        --result[getSecondMaximumClusterSize()];
    }

    for (auto it=result.begin(); it!= result.end(); ){
        it->second ? ++it : result.erase(it++);
    }

    return result;
}

const double NZ_Network::getMeanClusterSize() const {
    //* Get first moment of cluster size excluding maximum cluster
    const int firstMoment = networkSize - maximumClusterSize;

    //* Get second moment of cluster size excluding maximum cluster
    double secondMoment = std::pow(maximumClusterSize, 2.0) * (m_sortedCluster.at(maximumClusterSize) - 1);
    for (auto it = ++m_sortedCluster.rbegin(); it != m_sortedCluster.rend(); ++it) {
        secondMoment += std::pow(it->first, 2.0) * it->second;
    }
    return secondMoment / firstMoment;
}