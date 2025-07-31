#pragma once
#include "Cluster.h"
#include <map>
#include <vector>

class ClusterGraph {
public:
    std::map<std::string, Cluster> clusters;
    std::map<std::string, std::map<std::string, double>> edgeWeights; // adjacency: [src][dst] = weight

    void addCluster(const Cluster& c) {
        clusters[c.id] = c;
    }

    void addEdge(const std::string& from, const std::string& to, double weight) {
        edgeWeights[from][to] = weight;
        clusters[from].edges.insert(to);
    }

    std::vector<std::string> getNeighbors(const std::string& cid) const {
    std::vector<std::string> neighbors;
    if (edgeWeights.count(cid)) {
        for (const auto& kv : edgeWeights.at(cid)) {
            neighbors.push_back(kv.first);
        }
    }
    return neighbors;
}


    bool hasCluster(const std::string& cid) const {
        return clusters.find(cid) != clusters.end();
    }
};
