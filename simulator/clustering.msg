message ClusterInfoMsg {
    string sender;
    string clusterID;
    string CIP;
    int numNodes;
    int numEdges;
    int mergedNumNodes;
    int mergedNumEdges;
}

message ClusterUpdateMsg {
    string newClusterID;
    string newCIP;
}

message MergeRequestMsg {
    string sender;
    string clusterID;
}

message MergeAckMsg {
    string sender;
    string acceptedCluster;
}

message MergeNackMsg {
    string sender;
}
