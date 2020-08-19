#pragma once
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

//* Network class with merging cluster by Newman Ziff algorithm
class Network{
private:
    //* Size of Network
    int m_size;

    //* Number of link
    int m_linkNum{0};

    //* Maximum cluster size of the network
    int m_maximumCluster{1};
    int m_secondGiant{1};
    int m_deltaMaximumCluster{0};

    //* Store the parent of each nodes.
    std::vector<int> m_parent;

    //* Flag that Maximum Cluster Size is increased
    bool m_maximumClusterUpdated{false};

    //* sortedCluster[i] : number of i size cluster
    std::map<int,int> m_sortedCluster;

    //* age[i]: age of each node
    //* changedAge[0] : {age,size} of root1
    std::vector<int> m_birth;
    std::vector<std::vector<int>> m_changedAge;

public:
    Network(const int &t_size)
    : m_size(t_size)
    {
        //! Make every node to root node with size 1
        m_parent.resize(t_size,-1);

        //! Initialize Sorted Cluster
        m_sortedCluster[1]=t_size;

        //! Initialize age and accumulated age
        m_birth.resize(t_size);
        m_changedAge.resize(2,std::vector<int>{0,0});
    }

    //* Simple get functions
    int getMaximumCluster() const {return m_maximumCluster;}
    int getSecondGiant() const {return m_secondGiant;}
    bool getMaximumClusterUpdated() const {return m_maximumClusterUpdated;}
    int getClusterSize(const int &t_root) const {return -m_parent[t_root];}
    int getDeltaMaximumCluster() const {return m_deltaMaximumCluster;}
    std::vector<std::vector<int>> getChangedAge() const {return m_changedAge;}
    std::map<int,int> getSortedCluster(const int &excludeNum=1) const {
        std::map<int,int> result=m_sortedCluster;

        //! exclude maximum cluster
        --result[m_maximumCluster];

        //! exclude second giant
        if (excludeNum-1){
            --result[m_secondGiant];
        }
        return result;
    }

    //* get the root of input node
    int getRoot(const int &t_node){
        //! Input node is root
        if (m_parent[t_node]<0){return t_node;}

        //! Find the parent node until the parent is root
        else {return m_parent[t_node]=getRoot(m_parent[t_node]);}
    }

    //* Merge two clusters
    void merge(const int &t_root1, const int &t_root2){
        //! upate line number
        m_linkNum++;

        //! Save age and reset birth time
        m_changedAge[0]=std::vector<int>{m_linkNum-m_birth[t_root1], -m_parent[t_root1]};
        m_changedAge[1]=std::vector<int>{m_linkNum-m_birth[t_root2], -m_parent[t_root2]};
        m_birth[t_root1]=m_linkNum;
        m_birth[t_root2]=m_linkNum;

        //! Get the size of each clusters
        const int size1=-m_parent[t_root1];
        const int size2=-m_parent[t_root2];
        const int newSize=size1+size2;

        //! Update Parent
        m_parent[t_root1]-=size2;
        m_parent[t_root2]=t_root1;

        //! Update Sorted Cluster
        --m_sortedCluster[size1];
        if (m_sortedCluster[size1]==0){
            m_sortedCluster.erase(size1);
        }
        --m_sortedCluster[size2];
        if (m_sortedCluster[size2]==0){
            m_sortedCluster.erase(size2);
        }
        ++m_sortedCluster[newSize];

        //! find maximum cluster
        if (m_maximumCluster<newSize){
            m_deltaMaximumCluster=newSize-m_maximumCluster;
            m_maximumCluster=newSize;
            m_maximumClusterUpdated=true;
        }
        else{
            m_maximumClusterUpdated=false;
        }

        //! find second giant cluster
        for (auto it=m_sortedCluster.rbegin(); it!=m_sortedCluster.rend(); ++it){
            if (it->first!=m_maximumCluster || it->second>1){
                m_secondGiant=it->first;
                break;
            }
        }
    }

    //* Mean cluster size
    double meanCluster() const{
        const int first=m_size-m_maximumCluster;
        double second=0;
        for (auto it=m_sortedCluster.begin(); it!=m_sortedCluster.end(); ++it){
            second+=pow(it->first,2)*it->second;
        }
        //! exclude infinite size cluster
        second-=pow(m_maximumCluster,2);

        return second/first;
    }
};