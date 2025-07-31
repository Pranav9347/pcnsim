#pragma once
#include <string>
#include <set>
#include <map>

struct Cluster {
    std::string id;
    std::string CIP;
    std::set<std::string> members;
    std::set<std::string> edges; // clusterIDs of neighboring clusters

    Cluster() {}

    Cluster(std::string cid, std::string cip) : id(cid), CIP(cip) {
        members.insert(cip);
    }
};
