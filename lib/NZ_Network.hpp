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
    int network_size{0};
    unsigned long long link_size{0};
    int maximum_cluster_size{1};
    int delta_maximum_cluster_size{0};
    std::vector<std::pair<unsigned long long, int>> changed_age; //* changed_age[0,1] : {age, size} of two merged roots

  protected:
    std::map<int, int> m_sorted_cluster; //* sorted_cluster[size] : number of cluster of 'size'
    std::vector<NZ_Node> m_nodes;

  public:
    //* Constructor
    NZ_Network() {}
    NZ_Network(const int&);

    //* Get the root of input node
    const int get_root(const int&);

    //* Merge two clusters
    //* Two input index should be different root
    void merge(const int&, const int&);
    const int get_size(const int&) const;
    const int get_second_maximum_cluster_size() const;
    const std::map<int, int> get_sorted_cluster(const int& t_exclude_num = 1) const;
    const double get_mean_cluster_size() const;

  protected:
    void _update_age(const int&, const int&, const int&, const int&);
    void _update_sorted_cluster(const int&, const int&);
    void _update_maximum_cluster_size(const int&);
};

NZ_Network::NZ_Network(const int& t_network_size)
    : network_size(t_network_size) {
    //* Generate default nodes
    m_nodes.reserve(t_network_size);
    for (int index = 0; index < t_network_size; ++index) {
        NZ_Node node(index);
        m_nodes.emplace_back(node);
    }

    //* Initialize sorted cluster and age
    m_sorted_cluster[1] = t_network_size;
    changed_age.assign(2, std::pair<int, int>{0, 0});
}

const int NZ_Network::get_root(const int& t_index) {
    if (m_nodes[t_index].parent < 0) {
        return t_index;
    }
    return m_nodes[t_index].parent = get_root(m_nodes[t_index].parent);
}

const int NZ_Network::get_size(const int& t_root) const {
    return -1 * m_nodes[t_root].parent;
}

void NZ_Network::_update_age(const int& t_root1,
                             const int& t_size1,
                             const int& t_root2,
                             const int& t_size2) {
    changed_age[0] = std::pair<unsigned long long, int>{link_size - m_nodes[t_root1].birth, t_size1};
    changed_age[1] = std::pair<unsigned long long, int>{link_size - m_nodes[t_root2].birth, t_size2};
    m_nodes[t_root1].birth = link_size;
}

void NZ_Network::_update_sorted_cluster(const int& t_size1,
                                       const int& t_size2) {
    --m_sorted_cluster[t_size1];
    if (m_sorted_cluster[t_size1] == 0) {
        m_sorted_cluster.erase(t_size1);
    }
    --m_sorted_cluster[t_size2];
    if (m_sorted_cluster[t_size2] == 0) {
        m_sorted_cluster.erase(t_size2);
    }
    ++m_sorted_cluster[t_size1 + t_size2];
}

void NZ_Network::_update_maximum_cluster_size(const int& t_newSize) {
    if (maximum_cluster_size < t_newSize) {
        delta_maximum_cluster_size = t_newSize - maximum_cluster_size;
        maximum_cluster_size = t_newSize;
    } else {
        delta_maximum_cluster_size = 0;
    }
}

void NZ_Network::merge(const int& t_root1,
                       const int& t_root2) {
    //* Update link size
    ++link_size;

    //* Get size of two clusters
    const int size1 = get_size(t_root1);
    const int size2 = get_size(t_root2);
    const int newSize = size1 + size2;

    //* Optional observables
    _update_maximum_cluster_size(newSize);
    _update_age(t_root1, size1, t_root2, size2);
    _update_sorted_cluster(size1, size2);

    //* Update parents
    m_nodes[t_root1].parent -= size2;
    m_nodes[t_root2].parent = t_root1;
}

const int NZ_Network::get_second_maximum_cluster_size() const {
    for (auto it = m_sorted_cluster.rbegin(); it != m_sorted_cluster.rend(); ++it) {
        if (it->first != maximum_cluster_size || it->second > 1) {
            return it->first;
        }
    }
    return 1;
}

const std::map<int, int> NZ_Network::get_sorted_cluster(const int& t_exclude_num) const {
    std::map<int, int> result = m_sorted_cluster;

    //* Exclue maximum cluster
    --result[maximum_cluster_size];
    if (result[maximum_cluster_size] == 0) {
        result.erase(maximum_cluster_size);
    }

    //* If exclude number = 2, exclude second maximum cluster
    if (t_exclude_num == 2) {
        const int second_maximum_cluster_size = get_second_maximum_cluster_size();
        --result[second_maximum_cluster_size];
        if (result[second_maximum_cluster_size]) {
            result.erase(second_maximum_cluster_size);
        }
    }
    return result;
}

const double NZ_Network::get_mean_cluster_size() const {
    //* Get first moment of cluster size excluding maximum cluster
    const int first_moment = network_size - maximum_cluster_size;

    //* Get second moment of cluster size excluding maximum cluster
    double second_moment = std::pow(maximum_cluster_size, 2.0) * (m_sorted_cluster.at(maximum_cluster_size) - 1);
    for (auto it = ++m_sorted_cluster.rbegin(); it != m_sorted_cluster.rend(); ++it) {
        second_moment += std::pow(it->first, 2.0) * it->second;
    }
    return second_moment / first_moment;
}