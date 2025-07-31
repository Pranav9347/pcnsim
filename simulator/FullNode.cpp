#include "globals.h"
#include "crypto.h"
#include "baseMessage_m.h"
#include "payment_m.h"
#include "PaymentChannel.h"
#include "invoice_m.h"
#include "commitmentSigned_m.h"
#include "revokeAndAck_m.h"
#include "paymentRefused_m.h"
#include "HTLC.h"
#include "FullNode.h"
#include "clustering_m.h"
#include "clusterGraph.h"
#include "Cluster.h"
std::map<std::string, std::string> FullNode::nodeToCluster;
ClusterGraph FullNode::clusterGraph;
// Define module and initialize random number generator
Define_Module(FullNode);

/***********************************************************************************************************************/
/* OMNETPP FUNCTIONS                                                                                                   */
/***********************************************************************************************************************/

void FullNode::sendClusterInfoTo(const std::string& neighbor) {
    ClusterInfoMsg* msg = new ClusterInfoMsg();

    msg->setSender(getName());
    msg->setClusterID(clusterID.c_str());
    msg->setCIP(CIP.c_str());

    int clusterNodeCount = clusterMembers.size();  // Real number of nodes in this cluster
    int intraClusterEdgeCount = countClusterEdges(clusterMembers);  // Actual EC(vi)

    msg->setNumNodes(clusterNodeCount);
    msg->setNumEdges(intraClusterEdgeCount);

    // These fields are not used anymore, can be set to zero or removed
    msg->setMergedNumNodes(0);
    msg->setMergedNumEdges(0);

    send(msg, "out", rtable[neighbor]);
}



ClusterInfo FullNode::waitForClusterInfoFrom(const std::string& neighbor) {
    ClusterInfoMsg msg = receivedClusterInfos[neighbor];  // Filled in handleMessage()

    ClusterInfo ci(
        msg.getSender(),
        msg.getClusterID(),
        msg.getCIP(),
        msg.getNumNodes(),
        msg.getNumEdges(),
        msg.getMergedNumNodes(),
        msg.getMergedNumEdges()
    );
    return ci;
}
bool FullNode::areNeighbors(const std::string& nodeA, const std::string& nodeB) {
    if (nameToPCs.find(nodeA) == nameToPCs.end()) return false;
    return nameToPCs[nodeA].count(nodeB) > 0;
}


int FullNode::countInterClusterEdges(const std::set<std::string>& clusterA, const std::set<std::string>& clusterB) {
    int count = 0;
    for (const std::string& nodeA : clusterA) {
        if (nameToPCs.count(nodeA) == 0) continue;

        for (const auto& [neighbor, pc] : nameToPCs[nodeA]) {
            if (clusterB.count(neighbor)) {
                count++;
                EV << "[DEBUG] Inter-cluster edge: " << nodeA << " -- " << neighbor << "\n";
            }
        }
    }
    return count;
}


int FullNode::countAllInterClusterEdges() {
    std::set<std::pair<std::string, std::string>> countedEdges;
    int count = 0;

    for (const auto& [nodeA, neighbors] : nameToPCs) {
        for (const auto& [nodeB, pc] : neighbors) {
            // Check only one direction to avoid double-counting
            if (nodeA < nodeB) {
                std::string clusterA = nodeToCluster[nodeA];
                std::string clusterB = nodeToCluster[nodeB];

                if (clusterA != clusterB) {
                    countedEdges.insert({nodeA, nodeB});
                    count++;
                }
            }
        }
    }

    return count;
}


double FullNode::computeDeltaS(
    const std::set<std::string>& A_members,
    const std::set<std::string>& B_members
) {
    int V = _localTopology->getNumNodes();
    if (V == 0) return 0.0; // Avoid division by zero

    // 1. Calculate the change in inter-cluster edges: |L| - |L'|
    // This is simply the number of edges between cluster A and B that will become intra-cluster.
    int L_change = countInterClusterEdges(A_members, B_members);

    // 2. Calculate the change in the weighted intra-cluster edge term.
    int V_A = A_members.size();
    int E_A = countClusterEdges(A_members);

    int V_B = B_members.size();
    int E_B = countClusterEdges(B_members);

    std::set<std::string> merged_members = A_members;
    merged_members.insert(B_members.begin(), B_members.end());
    int V_merged = merged_members.size();
    int E_merged = countClusterEdges(merged_members);

    double weighted_edges_before = static_cast<double>(V_A * E_A + V_B * E_B);
    double weighted_edges_after = static_cast<double>(V_merged * E_merged);

    // 3. Combine the terms according to the paper's formula
    double deltaS = L_change + (weighted_edges_before - weighted_edges_after) / V;

    EV << "[DEBUG] ΔS with " << (B_members.empty() ? "N/A" : *B_members.begin()) << ":\n"
       << "  L_change: " << L_change << "\n"
       << "  V_A=" << V_A << ", E_A=" << E_A << "\n"
       << "  V_B=" << V_B << ", E_B=" << E_B << "\n"
       << "  V_merged=" << V_merged << ", E_merged=" << E_merged << "\n"
       << "  Weighted Term Change: " << (weighted_edges_before - weighted_edges_after) / V << "\n"
       << "  Final ΔS = " << deltaS << "\n";

    return deltaS;
}





void FullNode::sendMergeRequest(const std::string& targetCIP) {
    MergeRequestMsg* msg = new MergeRequestMsg();
    msg->setSender(getName());
    msg->setClusterID(clusterID.c_str());
    send(msg, "out", rtable[targetCIP]);
    EV << "[Cluster] Sent MergeRequest to " << targetCIP << "\n";
}

// In FullNode.cpp
// In FullNode.cpp
void FullNode::handleMergeRequest(MergeRequestMsg* msg) {
    std::string proposerName = msg->getSender();
    EV << "[Cluster] Received MergeRequest from " << proposerName << "\n";

    // Condition 1: Check if I am still an active CIP and haven't already accepted a merge this round.
    if (CIP != getName() || hasMergedThisRound) {
        EV << "[Cluster] REJECTING " << proposerName << " because I am no longer available to merge.\n";
        MergeNackMsg* nackMsg = new MergeNackMsg();
        nackMsg->setSender(getName());
        send(nackMsg, "out", rtable[proposerName]);
        delete msg;
        return;
    }

    // Condition 2: Check if merging with the proposer is beneficial (i.e., has a positive deltaS).
    if (deltaSByNeighbor.count(proposerName) && deltaSByNeighbor.at(proposerName) > 0) {
        EV << "[Cluster] ACCEPTING proposal from " << proposerName << ". Locking myself and sending accept.\n";
        hasMergedThisRound = true; // Lock myself from other offers this round.

        // Send the acceptance message back to the proposer.
        MergeAcceptMsg* acceptMsg = new MergeAcceptMsg();
        acceptMsg->setSender(getName());
        send(acceptMsg, "out", rtable[proposerName]);
    } else {
        // Reject if they are not a beneficial partner.
        EV << "[Cluster] REJECTING " << proposerName << " because the merge is not beneficial.\n";
        MergeNackMsg* nackMsg = new MergeNackMsg();
        nackMsg->setSender(getName());
        send(nackMsg, "out", rtable[proposerName]);
    }

    delete msg;
}

void FullNode::handleMergeAck(MergeAckMsg* msg) {
    std::string sender = msg->getSender();
    EV << "[Cluster] Received MergeAck from " << sender << "\n";
    // Store merge intent (you can refine this with a set/map)
    mergeResponses.insert(sender);
    delete msg;
}
void FullNode::handleMergeNack(MergeNackMsg* msg) {
    std::string rejectorName = msg->getSender();

    if (sentMergeRequestTo == rejectorName) {
        EV << "[Cluster] My merge proposal to " << rejectorName << " was REJECTED.\n";
        sentMergeRequestTo = ""; // Reset state.
    }

    delete msg;
}
// In FullNode.cpp
void FullNode::handleMergeAccept(MergeAcceptMsg* msg) {
    std::string acceptorName = msg->getSender();

    // NEW CHECK: If I have already merged this round, I must ignore this acceptance.
    if (hasMergedThisRound) {
        EV << "[Cluster] Received MergeAccept from " << acceptorName << ", but I have already merged. Ignoring.\n";
        delete msg;
        return;
    }

    if (sentMergeRequestTo == acceptorName) {
        EV << "[Cluster] My merge proposal to " << acceptorName << " was ACCEPTED!\n";
        EV << "[Cluster] I will now finalize the merge.\n";

        // NEW: Lock myself from any further merge activity this round.
        hasMergedThisRound = true;

        std::string newCIP = std::max(std::string(getName()), acceptorName);
        mergeClusters(getName(), acceptorName, newCIP);

        sentMergeRequestTo = "";
    } else {
        EV << "[Cluster] Received an unexpected MergeAccept from " << acceptorName << ". Ignoring.\n";
    }

    delete msg;
}

// In FullNode.cpp
void FullNode::mergeClusters(const std::string& nodeA_name, const std::string& nodeB_name, const std::string& mergedID) {
    EV << "[Cluster] Finalizing merge between " << nodeA_name << " and " << nodeB_name << " into new cluster: " << mergedID << "\n";

    // --- Step 1: Get the current cluster IDs for the two nodes ---
    // It's possible one of the nodes is already part of a larger cluster from a previous round.
    std::string clusterA_id = nodeToCluster[nodeA_name];
    std::string clusterB_id = nodeToCluster[nodeB_name];

    // If they are already in the same cluster, something is wrong. Abort.
    if (clusterA_id == clusterB_id) {
        EV << "[WARNING] Attempted to merge nodes that are already in the same cluster. Aborting merge.\n";
        return;
    }

    // --- Step 2: Get the full member lists from the static map ---
    std::set<std::string> membersOfA = clusterIDToMembers[clusterA_id];
    std::set<std::string> membersOfB = clusterIDToMembers[clusterB_id];

    // --- Step 3: Create the new merged set of members ---
    std::set<std::string> allNewMembers = membersOfA;
    allNewMembers.insert(membersOfB.begin(), membersOfB.end());

    // --- Step 4: Update the static maps (the single source of truth) ---
    // Point every member of the new cluster to the new cluster ID.
    for (const std::string& memberNode : allNewMembers) {
        nodeToCluster[memberNode] = mergedID;
    }

    // Update the cluster-to-members map.
    clusterIDToMembers[mergedID] = allNewMembers;
    // Remove the old, now-empty clusters.
    clusterIDToMembers.erase(clusterA_id);
    clusterIDToMembers.erase(clusterB_id);

    // --- Step 5: The node that initiated the merge becomes the new CIP and schedules the next round ---
    // This node's local state needs to be updated.
    this->clusterID = mergedID;
    this->CIP = mergedID;
    this->clusterMembers = allNewMembers;

    EV << "[Cluster] Merge complete. New cluster '" << mergedID << "' has " << allNewMembers.size() << " members.\n";

    // The node that performed the merge schedules the next round of clustering.
    scheduleAt(simTime() + SimTime(0.1, SIMTIME_S), new cMessage("nextClusteringRound"));

    // --- Step 6: Broadcast the update to all other members of the new cluster ---
    // The acceptor node and other members need to be told about their new cluster ID and CIP.
    for (const std::string& memberNode : allNewMembers) {
        if (memberNode != getName()) { // Don't send to self
            ClusterUpdateMsg* msg = new ClusterUpdateMsg();
            msg->setNewClusterID(mergedID.c_str());
            msg->setNewCIP(mergedID.c_str());
            send(msg, "out", rtable[memberNode]);
        }
    }
}


// Implementing the cluster formation algorithm (Algorithm 1 from the CRP paper)
int FullNode::countClusterEdges(const std::set<std::string>& members) {
    int count = 0;
    for (const std::string& n : members) {
        for (const auto& [neighbor, pc] : nameToPCs[n]) {
            if (members.count(neighbor)) {
                count++;
            }
        }
    }
    return count / 2;  // Each edge counted twice
}

void FullNode::runClusteringRound() {
    hasMergedThisRound = false; // Reset the flag at the start of each round
    deltaSByNeighbor.clear(); // Add this line to clear old results
    std::string myName = getName();

    if (CIP != myName) {
        EV << "[Cluster] Skipping clustering round — not a CIP: " << myName << "\n";
        return;
    }

    // Clear old state
    receivedClusterInfos.clear();
    mergeResponses.clear();

    // Step 1: Send EC(vi) to all neighbor CIPs
    for (const std::string& neighbor : interClusterNeighbors) {
        sendClusterInfoTo(neighbor);
    }

    // Now we wait for all ClusterInfoMsg to arrive
    // and proceed in onAllClusterInfoReceived()
}

// In FullNode.cpp

// In FullNode.cpp
void FullNode::onAllClusterInfoReceived() {
    std::string myName = getName();

    // --- 1. Calculate deltaS for each neighbor who responded ---
    // This loop correctly iterates over the received info, not the results map.
    for (const auto& [neighborName, ciMsg] : receivedClusterInfos) {
        // Get the members of the neighbor's cluster
        std::string neighborClusterID = ciMsg.getClusterID();
        std::set<std::string> neighborClusterMembers = clusterIDToMembers[neighborClusterID];

        // Calculate the benefit of merging with this neighbor
        double deltaS = computeDeltaS(clusterMembers, neighborClusterMembers);

        EV << "[Cluster] ΔS with " << neighborName << " (" << neighborClusterID << ") = " << deltaS << "\n";

        // Store the result if it's a positive benefit
        if (deltaS > 0) {
            deltaSByNeighbor[neighborName] = deltaS;
        }
    }

    // --- 2. Find the single best merge partner from the calculated results ---
    std::string bestMergeTarget = "";
    double maxDeltaS = 1e-9; // Use a small epsilon to avoid floating point issues

    for (const auto& pair : deltaSByNeighbor) {
        const std::string& currentNeighbor = pair.first;
        double currentDeltaS = pair.second;

        if (currentDeltaS > maxDeltaS) {
            maxDeltaS = currentDeltaS;
            bestMergeTarget = currentNeighbor;
        } else if (std::abs(currentDeltaS - maxDeltaS) < 1e-9) {
            if (currentNeighbor < bestMergeTarget) {
                bestMergeTarget = currentNeighbor;
            }
        }
    }

    // --- 3. If a good partner is found, send a proposal ---
    if (bestMergeTarget.empty()) {
        EV << "[Cluster] No suitable merge partners found this round for " << myName << ".\n";
        return;
    }

    EV << "[Cluster] Final Decision: Best merge target for " << myName << " is " << bestMergeTarget
       << " with ΔS = " << maxDeltaS << "\n";

    // Send the proposal and remember who we sent it to
    sendMergeRequest(bestMergeTarget);
    sentMergeRequestTo = bestMergeTarget;
}


void FullNode::initialize() {   
    std::string myName = getName();
    clusterID = myName;   // Initially, every node is its own cluster
    CIP = myName;         // Each node starts as its own Cluster Information Proxy (CIP)

    EV << "[Init] Node: " << myName << " initialized with clusterID: " << clusterID << " and acts as its own CIP.\n";

    // Step 1: Setup local topology view
    _localTopology = globalTopology;

    // Step 2: Copy workload from global structure
    std::map<std::string, std::vector<std::tuple<std::string, double, simtime_t>>> localPendingPayments = pendingPayments;

    // Step 3: Initialize Payment Channels and statistics
    for (auto& neighborToPCs : nameToPCs[myName]) {
        std::string neighborName = neighborToPCs.first;
        auto pcTuple = neighborToPCs.second;

        double capacity               = std::get<0>(pcTuple);
        double fee                    = std::get<1>(pcTuple);
        double quality                = std::get<2>(pcTuple);
        int maxAcceptedHTLCs         = std::get<3>(pcTuple);
        int numHTLCs                 = 0;
        double HTLCMinimumMsat       = std::get<4>(pcTuple);
        double channelReserveSatoshis= std::get<5>(pcTuple);
        cGate* localGate             = std::get<6>(pcTuple);
        cGate* neighborGate          = std::get<7>(pcTuple);

        PaymentChannel pc = PaymentChannel(capacity, fee, quality, maxAcceptedHTLCs, numHTLCs,
                                           HTLCMinimumMsat, channelReserveSatoshis, localGate, neighborGate);

        _paymentChannels[neighborName] = pc;

        // Register signal/statistic for capacity
        std::string signalName = myName + "-to-" + neighborName + ":capacity";
        simsignal_t signal = registerSignal(signalName.c_str());
        _signals[signalName] = signal;
        emit(_signals[signalName], capacity);

        cProperty* statTemplate = getProperties()->get("statisticTemplate", "pcCapacities");
        getEnvir()->addResultRecorders(this, signal, signalName.c_str(), statTemplate);

        // Build initial routing table
        rtable[neighborName] = localGate->getIndex();

        // TEMPORARY (before clustering round): Treat all neighbors as inter-cluster
        if (neighborName != myName) {
            interClusterNeighbors.insert(neighborName);
        }
    }

    // Step 4: Register self as initial cluster member (only self at start)
    clusterMembers.insert(myName);
    knownClusterMembers[clusterID] = clusterMembers;
    clusterIDToMembers[clusterID].insert(getName());
    nodeToCluster[myName] = clusterID;

    // Step 5: Run clustering round
    runClusteringRound();

    // Step 6: After clustering, populate intra-cluster neighbor set
    intraClusterNeighbors.clear();
    for (const std::string& member : knownClusterMembers[clusterID]) {
        if (member != myName) {
            intraClusterNeighbors.insert(member);
        }
    }
    intraClusterNeighbors.insert(myName); // include self for Dijkstra

    // Debug: Print intra-cluster neighbors
    EV << "[CRP] Intra-cluster neighbors of " << myName << ": ";
    for (const auto& n : intraClusterNeighbors) {
        EV << n << " ";
    }
    EV << endl;

    // Step 7: Register module-wide stats
    initPerModuleStatistics();

    // Step 8: Schedule payments (workload)
    if (pendingPayments.find(myName) != pendingPayments.end()) {
        std::vector<std::tuple<std::string, double, simtime_t>> myWorkload = pendingPayments[myName];

        for (const auto& paymentTuple : myWorkload) {
            std::string srcName = std::get<0>(paymentTuple);
            double value        = std::get<1>(paymentTuple);
            simtime_t time      = std::get<2>(paymentTuple);

            char msgname[100];
            sprintf(msgname, "%s-to-%s;value:%0.1f", srcName.c_str(), myName.c_str(), value);

            Payment* trMsg = new Payment(msgname);
            trMsg->setSource(srcName.c_str());
            trMsg->setDestination(myName.c_str());
            trMsg->setValue(value);
            trMsg->setHopCount(0);

            BaseMessage* baseMsg = new BaseMessage();
            baseMsg->setMessageType(TRANSACTION_INIT);
            baseMsg->setHopCount(0);
            baseMsg->encapsulate(trMsg);

            scheduleAt(simTime() + time + 10.0, baseMsg);
            _isFirstSelfMessage = true;
        }
    } else {
        EV << "No workload found for " << myName << ".\n";
    }
    EV << "[DEBUG] nameToPCs[" << getName() << "] contains:\n";
for (const auto& [neighbor, pc] : nameToPCs[getName()]) {
    EV << "  " << getName() << " -- " << neighbor << "\n";
}

}


// FullNode.cpp (continued)



void FullNode::handleMessage(cMessage *msg) {
    // Check for clustering-related messages first
    if (strcmp(msg->getName(), "nextClusteringRound") == 0) {
        delete msg;
        runClusteringRound();  // ← Your method that starts clustering again
        return;
    }
    if (ClusterInfoMsg* ci = dynamic_cast<ClusterInfoMsg*>(msg)) {
    std::string sender = ci->getSender();
    std::string clusterID = ci->getClusterID();

    receivedClusterInfos[sender] = *ci;

    clusterIDToMembers[clusterID].insert(sender);
    knownClusterMembers[clusterID].insert(sender);

    EV << "[DEBUG] Learned that " << sender << " is in cluster " << clusterID << "\n";

    if (receivedClusterInfos.size() == interClusterNeighbors.size()) {
        onAllClusterInfoReceived();
    }

    delete ci;
    return;
    }


    if (ClusterUpdateMsg* cu = dynamic_cast<ClusterUpdateMsg*>(msg)) {
        clusterID = cu->getNewClusterID();
        CIP = cu->getNewCIP();
        EV << "[Cluster] Updated to new clusterID: " << clusterID << ", CIP: " << CIP << "\n";
        delete cu;
        return;
    }


    if (MergeRequestMsg* mr = dynamic_cast<MergeRequestMsg*>(msg)) {
        handleMergeRequest(mr);
        return;
    }
    if (MergeAcceptMsg* ma = dynamic_cast<MergeAcceptMsg*>(msg)) {
        handleMergeAccept(ma);
        return;
    }
    if (MergeAckMsg* ma = dynamic_cast<MergeAckMsg*>(msg)) {
        handleMergeAck(ma);
        return;
    }

    if (MergeNackMsg* mn = dynamic_cast<MergeNackMsg*>(msg)) {
        handleMergeNack(mn);
        return;
    }

    // Proceed with PCN BaseMessage logic
    BaseMessage *baseMsg = check_and_cast<BaseMessage *>(msg);

    switch (baseMsg->getMessageType()) {
        case TRANSACTION_INIT:       initHandler(baseMsg); break;
        case INVOICE:                invoiceHandler(baseMsg); break;
        case UPDATE_ADD_HTLC:        updateAddHTLCHandler(baseMsg); break;
        case UPDATE_FULFILL_HTLC:    updateFulfillHTLCHandler(baseMsg); break;
        case UPDATE_FAIL_HTLC:       updateFailHTLCHandler(baseMsg); break;
        case PAYMENT_REFUSED:        paymentRefusedHandler(baseMsg); break;
        case COMMITMENT_SIGNED:      commitSignedHandler(baseMsg); break;
        case REVOKE_AND_ACK:         revokeAndAckHandler(baseMsg); break;
    }
}

std::vector<std::string> FullNode::intraClusterDijkstra(const std::string& src, const std::string& dst) {
    std::map<std::string, double> dist;
    std::map<std::string, std::string> prev;
    std::set<std::string> visited;

    // Initialize distances to infinity and prev to empty
    for (const std::string& node : intraClusterNeighbors) {
        dist[node] = std::numeric_limits<double>::infinity();
        prev[node] = "";
    }
    dist[src] = 0;

    auto minDistanceNode = [&dist, &visited]() -> std::string {
        double minDist = std::numeric_limits<double>::infinity();
        std::string minNode = "";
        for (const auto& [node, d] : dist) {
            if (visited.find(node) == visited.end() && d < minDist) {
                minDist = d;
                minNode = node;
            }
        }
        return minNode;
    };

    while (visited.size() < intraClusterNeighbors.size()) {
        std::string u = minDistanceNode();
        if (u == "") break;
        visited.insert(u);

        // For each neighbor of u within the same cluster
        for (const auto& [neighbor, pc] : nameToPCs[u]) {
            if (intraClusterNeighbors.find(neighbor) == intraClusterNeighbors.end()) continue;

            double capacity = std::get<0>(pc);  // 0th element is _capacity

            double weight = 1.0 / capacity;

            if (dist[u] + weight < dist[neighbor]) {
                dist[neighbor] = dist[u] + weight;
                prev[neighbor] = u;
            }
        }
    }

    // Reconstruct path from src to dst
    std::vector<std::string> path;
    std::string current = dst;
    while (current != "" && current != src) {
        path.insert(path.begin(), current);
        current = prev[current];
    }

    if (current == src) {
        path.insert(path.begin(), src);
    } else {
        path.clear();  // No path found
    }

    return path;
}
std::vector<std::vector<std::string>> FullNode::yensKShortestPaths(
    const ClusterGraph& graph,
    const std::string& srcCluster,
    const std::string& dstCluster,
    int K)
{
    using Path = std::vector<std::string>;
    using WeightedPath = std::pair<double, Path>; // cost, path

    auto dijkstraCluster = [&](const std::string& src, const std::string& dst) -> WeightedPath {
        std::map<std::string, double> dist;
        std::map<std::string, std::string> prev;
        std::set<std::string> visited;

        for (const auto& kv : graph.clusters) {
            dist[kv.first] = std::numeric_limits<double>::infinity();
            prev[kv.first] = "";
        }
        dist[src] = 0;

        while (!visited.count(dst)) {
            std::string u;
            double minDist = std::numeric_limits<double>::infinity();
            for (const auto& kv : dist) {
                if (!visited.count(kv.first) && kv.second < minDist) {
                    minDist = kv.second;
                    u = kv.first;
                }
            }
            if (u.empty()) break;

            visited.insert(u);
            for (const auto& neighbor : graph.getNeighbors(u)) {
                double weight = graph.edgeWeights.at(u).at(neighbor);
                if (dist[u] + weight < dist[neighbor]) {
                    dist[neighbor] = dist[u] + weight;
                    prev[neighbor] = u;
                }
            }
        }

        Path path;
        std::string current = dst;
        while (!current.empty()) {
            path.insert(path.begin(), current);
            current = prev[current];
        }

        return { dist[dst], path };
    };

    std::vector<WeightedPath> A;  // shortest paths found
    std::set<Path> pathSet;
    std::vector<WeightedPath> B;  // candidate paths

    A.push_back(dijkstraCluster(srcCluster, dstCluster));
    pathSet.insert(A[0].second);

    for (int k = 1; k < K; ++k) {
        const Path& lastPath = A[k - 1].second;
        for (size_t i = 0; i < lastPath.size() - 1; ++i) {
            Path rootPath(lastPath.begin(), lastPath.begin() + i + 1);

            // Temporarily remove edges
            ClusterGraph modifiedGraph = graph;
            for (const auto& p : A) {
                const Path& other = p.second;
                if (other.size() > i && std::equal(rootPath.begin(), rootPath.end(), other.begin())) {
                    std::string u = other[i];
                    std::string v = other[i + 1];
                    modifiedGraph.edgeWeights[u].erase(v);
                }
            }

            auto spurPath = dijkstraCluster(rootPath.back(), dstCluster);
            if (!spurPath.second.empty()) {
                Path totalPath = rootPath;
                totalPath.insert(totalPath.end(), spurPath.second.begin() + 1, spurPath.second.end());
                if (pathSet.find(totalPath) == pathSet.end()) {
                    double totalCost = 0.0;
                    for (size_t j = 0; j < totalPath.size() - 1; ++j) {
                        totalCost += graph.edgeWeights.at(totalPath[j]).at(totalPath[j + 1]);
                    }
                    B.push_back({ totalCost, totalPath });
                    pathSet.insert(totalPath);
                }
            }
        }

        if (B.empty()) break;

        std::sort(B.begin(), B.end(), [](const WeightedPath& a, const WeightedPath& b) {
            return a.first < b.first;
        });

        A.push_back(B.front());
        B.erase(B.begin());
    }

    std::vector<Path> result;
    for (const auto& wp : A) result.push_back(wp.second);
    return result;
}
std::vector<PathMetrics> FullNode::probeCandidatePaths(const std::vector<std::vector<std::string>>& clusterPaths, double paymentAmount) {
    std::vector<PathMetrics> results;

    for (const auto& path : clusterPaths) {
        double totalFee = 0.0;
        double minCapacity = std::numeric_limits<double>::infinity();
        double totalQuality = 0.0;

        for (size_t i = 0; i < path.size() - 1; ++i) {
            std::string fromCluster = path[i];
            std::string toCluster = path[i + 1];

            // Simulated: you can later replace this with a real message to CIP of fromCluster
            double fee = clusterGraph.edgeWeights[fromCluster][toCluster];   // simulate fee
            double capacity = 100.0;   // dummy capacity, query actual if known
            double quality = 1.0;      // dummy latency or reliability

            totalFee += fee;
            totalQuality += quality;
            minCapacity = std::min(minCapacity, capacity);
        }

        results.emplace_back(path, totalFee, minCapacity, totalQuality);
    }

    return results;
}
PathMetrics FullNode::selectOptimalRouteLP(const std::vector<PathMetrics>& candidatePaths, double paymentAmount) {
    double bestScore = std::numeric_limits<double>::infinity();
    PathMetrics bestPath({}, 0, 0, 0);

    for (const auto& p : candidatePaths) {
        if (p.minCapacity < paymentAmount) continue;

        // Score can be weighted: fee + α × quality
        double score = p.totalFee + 0.1 * p.totalQuality;

        if (score < bestScore) {
            bestScore = score;
            bestPath = p;
        }
    }

    return bestPath;
}


void FullNode::refreshDisplay() const {

    for(auto& it : _paymentChannels) {
        char buf[30];
        std::string neighborName = it.first;
        std::string neighborPath = "PCN." + neighborName;
        float capacity = it.second.getCapacity();
        cGate *gate = it.second.getLocalGate();
        cChannel *channel = gate->getChannel();
        sprintf(buf, "%0.1f\n", capacity);
        channel->getDisplayString().setTagArg("t", 0, buf);
        channel->getDisplayString().setTagArg("t", 1, "l");
    }
}

void FullNode::finish() {

    std::string myName = getName();
    int countCompleted = 0;
    int countFailed = 0;
    int countCanceled = 0;

    for (auto const & payment :_myPayments) {
        if (payment.second == "COMPLETED")
            countCompleted++;
        else if (payment.second == "FAILED")
            countFailed++;
        else if (payment.second == "CANCELED")
            countCanceled++;
        else
            continue;
    }

    int countTotal = countCompleted + countFailed + countCanceled;

    if (countTotal != 0 ) {
        EV << "------------------ Statistics for node " + myName + "------------------\n";
        double goodput = (double(countCompleted)/double(countFailed+countCompleted));
        EV << "COMPLETED/FAILED/CANCELED: " + std::to_string(countCompleted) + "/" + std::to_string(countFailed) + "/" + std::to_string(countCanceled) + "\n";
        EV << "Goodput: " +  std::to_string(goodput) + "\n";
    }
    // In void FullNode::finish()

// Only have one node print the final cluster state to avoid duplicate output
if (strcmp(getName(), "node0") == 0) {
        EV << "------------------ Final Cluster State (printed by " << getName() << ") ------------------\n";
        // ADD THIS DEBUG LOOP
    // END OF DEBUG LOOP
        // Reverse the map to group nodes by cluster ID
        std::map<std::string, std::vector<std::string>> finalClusters;
        for (const auto& pair : nodeToCluster) {
            std::string nodeName = pair.first;
            std::string clusterId = pair.second;
            finalClusters[clusterId].push_back(nodeName);
        }

        // Print each cluster and its members
        EV << "\n--- Final Cluster Information ---\n";
        EV << "[DEBUG] Final number of clusters: " << finalClusters.size() << "\n";
for (const auto& pair : finalClusters) {
    std::string clusterId = pair.first;
    const auto& members = pair.second; // 'members' is the vector of strings

    // Build a string of all members in the list
    std::string memberList;
    for (const std::string& member : members) {
        memberList += member + " ";
    }

    // Now print the formatted output
    EV << "Cluster ID: " << clusterId << " (" << members.size() << " members): " << memberList << "\n";
}
        EV << "--------------------------------------------------------------------------------\n";
    }

}


/***********************************************************************************************************************/
/* ROUTING FUNCTIONS                                                                                                   */
/***********************************************************************************************************************/

std::string FullNode::minDistanceNode (std::map<std::string, double> distances, std::map<std::string, bool> visited) {
    // Selects the node with minimum distance out of a distances list while disregarding visited nodes

    // Initialize min value
    double minDist = INT_MAX;
    std::string minName = "";

    for (auto& node : distances) {
        std::string nodeName = node.first;
        double dist = node.second;
        bool isNodeVisited = visited.find(nodeName)->second;
        if (!isNodeVisited && dist <= minDist) {
            minDist = dist;
            minName = nodeName;
        }
    }
    if (minDist == INT_MAX){
        return NULL;
    }
    else{
        return minName;
    }

}

std::vector<std::string> FullNode::getPath (std::map<std::string, std::string> parents, std::string target) {
    // Return the path to source given a target and its parent nodes

    std::vector<std::string> path;
    std::string node = target;

    // Recursively traverse the parents path
    while (parents[node] != node) {
        path.insert(path.begin(), node);
        node = parents[node];
    }

    // Set ourselves as the source node
    path.insert(path.begin(), getName());

    return path;
}

std::vector<std::string> FullNode::runCRPRouting(
    const std::string& srcNode,
    const std::string& dstNode,
    double paymentAmount)
{
    std::string srcCluster = nodeToCluster[srcNode];
    std::string dstCluster = nodeToCluster[dstNode];

    EV << "[CRP] Routing payment of " << paymentAmount << " from " << srcNode << " (" << srcCluster
       << ") to " << dstNode << " (" << dstCluster << ")\n";

    // Step 1: Get K shortest cluster paths using Yen's algorithm
    int K = 3;
    auto clusterPaths = yensKShortestPaths(clusterGraph, srcCluster, dstCluster, K);
    if (clusterPaths.empty()) {
        EV << "[CRP] No cluster-level paths found\n";
        return {};
    }

    // Step 2: Probe each path for fee, quality, capacity
    auto metrics = probeCandidatePaths(clusterPaths, paymentAmount);
    if (metrics.empty()) {
        EV << "[CRP] No feasible candidate paths found\n";
        return {};
    }

    // Step 3: Select optimal path using LP/scoring
    auto bestPathMetrics = selectOptimalRouteLP(metrics, paymentAmount);
    const auto& clusterPath = bestPathMetrics.clusterPath;

    // Step 4: For each cluster segment, compute intra-cluster Dijkstra
    std::vector<std::string> fullNodePath;

    for (size_t i = 0; i < clusterPath.size(); ++i) {
        std::string currentCluster = clusterPath[i];
        std::string fromNode, toNode;

        if (i == 0) {
            fromNode = srcNode;
        } else {
            fromNode = getCIP(clusterPath[i]); // entry to cluster
        }

        if (i == clusterPath.size() - 1) {
            toNode = dstNode;
        } else {
            toNode = getCIP(clusterPath[i + 1]); // exit from current cluster
        }

        std::vector<std::string> subPath = intraClusterDijkstra(fromNode, toNode);
        if (subPath.empty()) {
            EV << "[CRP] No intra-cluster path from " << fromNode << " to " << toNode << "\n";
            return {};
        }

        if (!fullNodePath.empty()) subPath.erase(subPath.begin()); // avoid duplicate hops
        fullNodePath.insert(fullNodePath.end(), subPath.begin(), subPath.end());
    }

    EV << "[CRP] Selected node path: ";
    for (auto& node : fullNodePath) EV << node << " ";
    EV << "\n";

    return fullNodePath;
}

std::string FullNode::getCIP(const std::string& clusterID) {
    if (clusterGraph.hasCluster(clusterID)) {
        return clusterGraph.clusters[clusterID].CIP;
    }
    return "";
}

std::vector<std::string> FullNode::dijkstraWeightedShortestPath (std::string src, std::string target, std::map<std::string, std::vector<std::pair<std::string, std::vector<double> > > > graph) {
    // This function returns the Dijkstra's shortest path from a source to some target given an adjacency matrix

    int numNodes = graph.size();
    std::map<std::string, double> distances;
    std::map<std::string, bool> visited;
    std::map<std::string, std::string> parents; // a.k.a nodeToParent

    // Initialize distances as infinite and visited array as false
    for (auto & node : graph) {
        distances[node.first] = INT_MAX;
        visited[node.first] = false;
        parents[node.first] = "";
    }

    // Set source node distance to zero and parent to itself
    distances[src] = 0;
    parents[src] = src;

    for (int i = 0; i < numNodes-1; i++) {
        std::string node = minDistanceNode(distances, visited);
        visited[node] = true;

        std::vector<std::pair<std::string,std::vector<double>>>::iterator it;

        // Update distance value of neighbor nodes of the current node
        for (it = graph[node].begin(); it != graph[node].end(); it++){
            std::string neighbor = it->first;
            std::vector<double> weightVector = it->second;
            double capacity = weightVector[0];
            double fee = weightVector[1];
            double linkQuality = weightVector[2];

            // Define weight as a combination of parameters in the weightVector
            double linkWeight = 1/capacity;

            if (!visited[neighbor]) {
                if(distances[node] + linkWeight < distances[neighbor]) {
                    parents[neighbor] = node;
                    distances[neighbor] = distances[node] + linkWeight;
                }
            }
        }
    }
    return getPath(parents, target);
}


/***********************************************************************************************************************/
/* MESSAGE HANDLERS                                                                                                    */
/***********************************************************************************************************************/

void FullNode::initHandler (BaseMessage *baseMsg) {

    Payment *initMsg = check_and_cast<Payment *> (baseMsg->decapsulate());
    EV << "TRANSACTION_INIT received. Starting payment "<< initMsg->getName() << "\n";

    // Create ephemeral communication channel with the payment source
    std::string srcName = initMsg->getSource();
    std::string srcPath = "PCN." + srcName;
    double value = initMsg->getValue();
    cModule* srcMod = findModuleByPath(srcPath.c_str());
    if (srcMod == nullptr) {
        EV_ERROR << "TRANSACTION_INIT failed: Could not find module " << srcPath << "\n";
        delete baseMsg;
        return;
    }

    cGate* myGate = this->getOrCreateFirstUnconnectedGate("out", 0, false, true);
    cGate* srcGate = srcMod->getOrCreateFirstUnconnectedGate("in", 0, false, true);
    cDelayChannel *tmpChannel = cDelayChannel::create("tmpChannel");
    tmpChannel->setDelay(100);
    myGate->connectTo(srcGate, tmpChannel);

    // Create invoice and send it to the payment source
    Invoice *invMsg = generateInvoice(srcName, value);
    baseMsg->setMessageType(INVOICE);
    baseMsg->encapsulate(invMsg);
    baseMsg->setName("INVOICE");
    send(baseMsg, myGate);

    // Close ephemeral connection
    myGate->disconnect();
}

void FullNode::invoiceHandler (BaseMessage *baseMsg) {

    Invoice *invMsg = check_and_cast<Invoice *> (baseMsg->decapsulate());
    EV << "INVOICE received. Payment hash: " << invMsg->getPaymentHash() << "\n";

    std::string dstName = invMsg->getDestination();
    std::string myName = getName();
    std::string paymentHash = invMsg->getPaymentHash();
    int htlcType = UPDATE_ADD_HTLC;
    std::string htlcId = createHTLCId(paymentHash, htlcType);
    double value = invMsg->getValue();

    // Find route to destination
    std::vector<std::string> path = runCRPRouting(myName, dstName, value);

if (path.empty()) {
    EV << "[CRP] No path found for invoice to " << dstName << ". Canceling.\n";
    _myPayments[paymentHash] = "CANCELED";
    _countCanceled++;
    _paymentGoodputAll = double(_countCompleted)/double(_countCompleted + _countFailed + _countCanceled);
    emit(_signals["canceledPayments"], _countCanceled);
    emit(_signals["paymentGoodputAll"], _paymentGoodputAll);
    return;
}

std::string firstHop = path[1];


    // If payment is larger than our capacity in the outbound payment channel, mark is as canceled and return
   if (!hasCapacityToForward(firstHop, value)) {
       _myPayments[paymentHash] = "CANCELED";
       EV << "WARNING: Canceling payment " + paymentHash + " on node " + myName + " due to insufficient funds in the first hop.\n";

       _countCanceled++;
       _paymentGoodputAll = double(_countCompleted)/double(_countCompleted + _countFailed + _countCanceled);

       emit(_signals["canceledPayments"], _countCanceled);
       emit(_signals["paymentGoodputAll"], _paymentGoodputAll);

       return;
   }

   // Add payment into payment list and set status = pending
   _myPayments[paymentHash] = "PENDING";

   // Print route
    std::string printPath = "Full route to destination: ";
    for (auto hop: path)
        printPath = printPath + hop + ", ";
    printPath += "\n";
    EV << printPath;

    //Create HTLC
    EV << "Creating HTLC to kick off the payment process \n";
    BaseMessage *newMessage = new BaseMessage();
    newMessage->setDestination(dstName.c_str());
    newMessage->setMessageType(UPDATE_ADD_HTLC);
    newMessage->setHopCount(1);
    newMessage->setHops(path);
    newMessage->setName("UPDATE_ADD_HTLC");
    newMessage->setDisplayString("i=block/encrypt;is=s");

    UpdateAddHTLC *firstUpdateAddHTLC = new UpdateAddHTLC();
    firstUpdateAddHTLC->setHtlcId(htlcId.c_str());
    firstUpdateAddHTLC->setSource(myName.c_str());
    firstUpdateAddHTLC->setPaymentHash(paymentHash.c_str());
    firstUpdateAddHTLC->setValue(value);

    HTLC *firstHTLC = new HTLC(firstUpdateAddHTLC);
    _paymentChannels[firstHop].setPendingHTLC(htlcId, firstHTLC);
    _paymentChannels[firstHop].setLastPendingHTLCFIFO(firstHTLC);
    _paymentChannels[firstHop].setPreviousHopUp(htlcId, myName);

    newMessage->encapsulate(firstUpdateAddHTLC);
    cGate *gate = _paymentChannels[firstHop].getLocalGate();

    //Sending HTLC out
    EV << "Sending HTLC to " + firstHop + " with payment hash " + paymentHash + "\n";
    send(newMessage, gate);

}

void FullNode::updateAddHTLCHandler (BaseMessage *baseMsg) {


    EV << "UPDATE_ADD_HTLC received at " + std::string(getName()) + " from " + std::string(baseMsg->getSenderModule()->getName()) + ".\n";
    std::string myName = getName();

    // If the message is a self message, it means we already attempted to commit changes but failed because the batch size was insufficient. So we wait for the timeout.
    // Otherwise, we attempt to commit normally.
    if (!baseMsg->isSelfMessage()){

        // Decapsulate message and get path
        UpdateAddHTLC *updateAddHTLCMsg = check_and_cast<UpdateAddHTLC *> (baseMsg->decapsulate());
        std::string dstName = baseMsg->getDestination();
        std::vector<std::string> path = baseMsg->getHops();
        std::string sender = baseMsg->getSenderModule()->getName();
        std::string paymentHash = updateAddHTLCMsg->getPaymentHash();
        int htlcType = UPDATE_ADD_HTLC;
        std::string htlcId = updateAddHTLCMsg->getHtlcId();
        double value = updateAddHTLCMsg->getValue();

        if ((myName == "node2") && (sender == "node0"))
        {
            // Set breakpoint here
            int x = 1;
        }


        // Create new HTLC in the backward direction and set it as pending
        HTLC *htlcBackward = new HTLC(updateAddHTLCMsg);
        EV << "Storing UPDATE_ADD_HTLC from node " + sender + " as pending.\n";
        EV << "Payment hash:" + paymentHash + ".\n";
        _paymentChannels[sender].setPendingHTLC(htlcId, htlcBackward);
        _paymentChannels[sender].setLastPendingHTLCFIFO(htlcBackward);
        _paymentChannels[sender].setPreviousHopUp(htlcId, sender);

        // If I'm the destination, trigger commit immediately and return
        if (dstName == this->getName()){
            EV << "Payment reached its destination. Not forwarding.\n";

            // Store base message for retrieving path on the first fulfill message later
            _myStoredMessages[paymentHash] = baseMsg;

            if (!tryCommitTxOrFail(sender, false)){
                EV << "Setting timeout for node " + myName + "\n";
                scheduleAt((simTime() + SimTime(500,SIMTIME_MS)),baseMsg);
            }
            return;
        }

        // If I'm not the destination, forward the message to the next hop in the UPSTREAM path
        std::string nextHop = path[baseMsg->getHopCount() + 1];
        std::string previousHop = path[baseMsg->getHopCount()-1];


        // Check if we have sufficient funds before forwarding
        if (!hasCapacityToForward(nextHop, value)) {
            // Not enough capacity to forward payment. Remove pending HTLCs and send a PAYMENT_REFUSED
            // message to the previous hop.
            _paymentChannels[sender].removePendingHTLC(htlcId);
            _paymentChannels[sender].removeLastPendingHTLCFIFO();
            _paymentChannels[sender].removePreviousHopUp(htlcId);

            BaseMessage *newMessage = new BaseMessage();
            newMessage->setDestination(previousHop.c_str());
            newMessage->setMessageType(PAYMENT_REFUSED);
            newMessage->setHopCount(baseMsg->getHopCount()-1);
            newMessage->setHops(path);
            newMessage->setName("PAYMENT_REFUSED");
            newMessage->setDisplayString("i=status/stop");

            PaymentRefused *paymentRefusedMsg = new PaymentRefused();
            paymentRefusedMsg->setPaymentHash(paymentHash.c_str());
            paymentRefusedMsg->setErrorReason("INSUFFICIENT CAPACITY");
            paymentRefusedMsg->setValue(value);

            newMessage->encapsulate(paymentRefusedMsg);

            cGate *gate = _paymentChannels[previousHop].getLocalGate();
            EV << "Sending PAYMENT_REFUSED to " + path[(newMessage->getHopCount())] + " with payment hash " + paymentRefusedMsg->getPaymentHash() + "\n";
            send(newMessage, gate);

        } else {
            // Enough funds. Forward HTLC.
            EV << "Creating HTLC to kick off the payment process \n";
            BaseMessage *newMessage = new BaseMessage();
            newMessage->setDestination(dstName.c_str());
            newMessage->setMessageType(UPDATE_ADD_HTLC);
            newMessage->setHopCount(baseMsg->getHopCount() + 1);
            newMessage->setHops(path);
            newMessage->setDisplayString("i=block/encrypt;is=s");
            newMessage->setName("UPDATE_ADD_HTLC");

            UpdateAddHTLC *newUpdateAddHTLC = new UpdateAddHTLC();
            newUpdateAddHTLC->setHtlcId(htlcId.c_str());
            newUpdateAddHTLC->setSource(myName.c_str());
            newUpdateAddHTLC->setPaymentHash(paymentHash.c_str());
            newUpdateAddHTLC->setValue(value);

            HTLC *htlcForward = new HTLC(newUpdateAddHTLC);

            // Add HTLC as pending in the forward direction and set previous hop as ourselves
            _paymentChannels[nextHop].setPendingHTLC(htlcId, htlcForward);
            _paymentChannels[nextHop].setLastPendingHTLCFIFO(htlcForward);
            _paymentChannels[nextHop].setPreviousHopUp(htlcId, myName);

            newMessage->encapsulate(newUpdateAddHTLC);

            cGate *gate = _paymentChannels[nextHop].getLocalGate();

            //Sending HTLC out
            EV << "Sending HTLC to " + path[(newMessage->getHopCount())] + " with payment hash " + paymentHash + "\n";
            send(newMessage, gate);

            if (!tryCommitTxOrFail(sender, false)){
                EV << "Setting timeout for node " + myName + ".\n";
                scheduleAt((simTime() + SimTime(500,SIMTIME_MS)),baseMsg);
            }
        }

    } else {
        // The message is the result of a timeout.
        EV << myName + " timeout expired. Creating commit.\n";
        std::vector<std::string> path = baseMsg->getHops();
        std::string previousHop = path[baseMsg->getHopCount()-1];
        std::string dstName = baseMsg->getDestination();
        tryCommitTxOrFail(previousHop, true);
    }
}

void FullNode::updateFulfillHTLCHandler (BaseMessage *baseMsg) {

    EV << "UPDATE_FULFILL_HTLC received at " + std::string(getName()) + " from " + std::string(baseMsg->getSenderModule()->getName()) + ".\n";
    std::string myName = getName();

    // If the message is a self message, it means we already attempted to commit changes but failed because the batch size was insufficient. So we wait for the timeout.
    // Otherwise, we attempt to commit normally.
    if (!baseMsg->isSelfMessage()){

        // Decapsulate message, get path, and preimage
        UpdateFulfillHTLC *fulfillHTLCMsg = check_and_cast<UpdateFulfillHTLC *> (baseMsg->decapsulate());
        std::string dstName = baseMsg->getDestination();
        std::vector<std::string> path = baseMsg->getHops();
        std::string sender = baseMsg->getSenderModule()->getName();
        std::string paymentHash = fulfillHTLCMsg->getPaymentHash();
        std::string preImage = fulfillHTLCMsg->getPreImage();
        double value = fulfillHTLCMsg->getValue();
        int htlcType = UPDATE_FULFILL_HTLC;
        std::string htlcId = fulfillHTLCMsg->getHtlcId();

         // Verify preimage
         if (sha256(preImage) != paymentHash){
             throw std::invalid_argument("ERROR: Failed to fulfill HTLC. Different hash value.");
         }

         // Create new HTLC in the backward direction and set it as pending
         HTLC *htlcBackward = new HTLC(fulfillHTLCMsg);
         EV << "Storing UPDATE_FULFILL_HTLC from node " + sender + " as pending.\n";
         EV << "Payment hash:" + paymentHash + ".\n";
         _paymentChannels[sender].setPendingHTLC(htlcId, htlcBackward);
         _paymentChannels[sender].setLastPendingHTLCFIFO(htlcBackward);
         _paymentChannels[sender].setPreviousHopDown(htlcId, sender);

         // If we are the destination, just try to commit the payment and return
         if (getName() == dstName) {
             EV << "Payment fulfillment has reached the payment's origin. Trying to commit...\n";
             if (!tryCommitTxOrFail(sender, false)) {
                 EV << "Setting timeout for node " + myName + "\n";
                 scheduleAt((simTime() + SimTime(500,SIMTIME_MS)),baseMsg);
             }
             return;
         }

        // If we're not the destination, forward the message to the next hop in the DOWNSTREAM path
        std::string nextHop = path[baseMsg->getHopCount()-1];
        std::string previousHop = path[baseMsg->getHopCount()+1];

        EV << "Forwarding UPDATE_FULFILL_HTLC in the downstream direction...\n";
        BaseMessage *newMessage = new BaseMessage();
        newMessage->setDestination(dstName.c_str());
        newMessage->setMessageType(UPDATE_FULFILL_HTLC);
        newMessage->setHopCount(baseMsg->getHopCount()-1);
        newMessage->setHops(path);
        newMessage->setDisplayString("i=block/decrypt;is=s");
        newMessage->setName("UPDATE_FULFILL_HTLC");

        UpdateFulfillHTLC *forwardFulfillHTLC = new UpdateFulfillHTLC();
        forwardFulfillHTLC->setHtlcId(htlcId.c_str());
        forwardFulfillHTLC->setPaymentHash(paymentHash.c_str());
        forwardFulfillHTLC->setPreImage(preImage.c_str());
        forwardFulfillHTLC->setValue(value);

        // Set UPDATE_FULFILL_HTLC as pending and invert the previous hop (now we're going downstream)
        HTLC *forwardBaseHTLC  = new HTLC(forwardFulfillHTLC);
        _paymentChannels[nextHop].setPendingHTLC(htlcId, forwardBaseHTLC);
        _paymentChannels[nextHop].setLastPendingHTLCFIFO(forwardBaseHTLC);
        _paymentChannels[nextHop].setPreviousHopDown(htlcId, myName);

        newMessage->encapsulate(forwardFulfillHTLC);

        cGate *gate = _paymentChannels[nextHop].getLocalGate();

        //Sending HTLC out
        EV << "Sending preimage " + preImage + " to " + path[(newMessage->getHopCount())] + " for payment hash " + paymentHash + "\n";
        send(newMessage, gate);

        // Try to commit
        if (!tryCommitTxOrFail(sender, false)){
            EV << "Setting timeout for node " + myName + ".\n";
            scheduleAt((simTime() + SimTime(500,SIMTIME_MS)),baseMsg);
        }

    } else {
        // The message is the result of a timeout.
        EV << myName + " timeout expired. Creating commit.\n";
        std::vector<std::string> path = baseMsg->getHops();
        std::string previousHop = path[baseMsg->getHopCount()+1];
        tryCommitTxOrFail(previousHop, true);
    }
}

void FullNode::updateFailHTLCHandler (BaseMessage *baseMsg) {

    EV << "UPDATE_FAIL_HTLC received at " + std::string(getName()) + " from " + std::string(baseMsg->getSenderModule()->getName()) + ".\n";
    std::string myName = getName();
    std::vector<std::string> path = baseMsg->getHops();
    std::string dstName = baseMsg->getDestination();

    // If the message is a self message, it means we already attempted to commit changes but failed because the batch size was insufficient. So we wait for the timeout.
    // Otherwise, we attempt to commit normally.
    if (!baseMsg->isSelfMessage()) {

        // Decapsulate message, get path, and preimage
        UpdateFailHTLC *failHTLCMsg = check_and_cast<UpdateFailHTLC *> (baseMsg->decapsulate());
        std::string sender = baseMsg->getSenderModule()->getName();
        std::string paymentHash = failHTLCMsg->getPaymentHash();
        std::string errorReason = failHTLCMsg->getErrorReason();
        double value = failHTLCMsg->getValue();
        int htlcType = UPDATE_FAIL_HTLC;
        std::string htlcId = failHTLCMsg->getHtlcId();

        // Create new HTLC in the backward direction and set it as pending
        HTLC *htlcBackward = new HTLC(failHTLCMsg);
        EV << "Storing UPDATE_FAIL_HTLC from node " + sender + " as pending.\n";
        EV << "Payment hash:" + paymentHash + ".\n";
        _paymentChannels[sender].setPendingHTLC(htlcId, htlcBackward);
        _paymentChannels[sender].setLastPendingHTLCFIFO(htlcBackward);
        _paymentChannels[sender].setPreviousHopDown(htlcId, sender);

        // If we are the destination, just try to commit and return
        if (getName() == dstName) {
            EV << "Payment fail has reached the payment's origin. Trying to commit...\n";
            if (!tryCommitTxOrFail(sender, false)) {
                EV << "Setting timeout for node " + myName + "\n";
                scheduleAt((simTime() + SimTime(500,SIMTIME_MS)),baseMsg);
            }
            return;
        }

        // If we're not the destination, forward the message to the next hop in the DOWNSTREAM path
        std::string nextHop = path[baseMsg->getHopCount()-1];
        std::string previousHop = path[baseMsg->getHopCount()+1];

        EV << "Forwarding UPDATE_FAIL_HTLC in the downstream direction...\n";
        BaseMessage *newMessage = new BaseMessage();
        newMessage->setDestination(dstName.c_str());
        newMessage->setMessageType(UPDATE_FAIL_HTLC);
        newMessage->setHopCount(baseMsg->getHopCount()-1);
        newMessage->setHops(path);
        newMessage->setDisplayString("i=status/stop");
        newMessage->setName("UPDATE_FAIL_HTLC");

        UpdateFailHTLC *forwardFailHTLC = new UpdateFailHTLC();
        forwardFailHTLC->setHtlcId(htlcId.c_str());
        forwardFailHTLC->setPaymentHash(paymentHash.c_str());
        forwardFailHTLC->setErrorReason(errorReason.c_str());
        forwardFailHTLC->setValue(value);

        // Set UPDATE_FAIL_HTLC as pending and invert the previous hop (now we're going downstream)
        HTLC *forwardBaseHTLC  = new HTLC(forwardFailHTLC);
        _paymentChannels[nextHop].setPendingHTLC(htlcId, forwardBaseHTLC);
        _paymentChannels[nextHop].setLastPendingHTLCFIFO(forwardBaseHTLC);
        //_paymentChannels[nextHop].removePreviousHopUp(htlcId);
        _paymentChannels[nextHop].setPreviousHopDown(htlcId, myName);

        newMessage->encapsulate(forwardFailHTLC);

        cGate *gate = _paymentChannels[nextHop].getLocalGate();

        //Sending HTLC out
        EV << "Sending UPDATE_FAIL_HTLC to " + path[(newMessage->getHopCount())] + "for payment hash " + paymentHash + "\n";
        send(newMessage, gate);

        // Try to commit
        if (!tryCommitTxOrFail(sender, false)){
            EV << "Setting timeout for node " + myName + ".\n";
            scheduleAt((simTime() + SimTime(500,SIMTIME_MS)),baseMsg);
        }

    } else {
        // The message is the result of a timeout.
        EV << myName + " timeout expired. Creating commit.\n";
        std::string previousHop = path[baseMsg->getHopCount()+1];
        tryCommitTxOrFail(previousHop, true);
    }
}

void FullNode::paymentRefusedHandler (BaseMessage *baseMsg) {

    std::string myName = getName();

    EV << "PAYMENT_REFUSED received at " + std::string(getName()) + " from " + std::string(baseMsg->getSenderModule()->getName()) + ".\n";

    PaymentRefused *paymentRefusedMsg = check_and_cast<PaymentRefused *> (baseMsg->getEncapsulatedPacket());
    std::string sender = baseMsg->getSenderModule()->getName();
    std::string paymentHash = paymentRefusedMsg->getPaymentHash();
    std::string errorReason = paymentRefusedMsg->getErrorReason();
    double value = paymentRefusedMsg->getValue();

    if ((myName == "node0"))
    {
        // Set breakpoint here
        int x = 1;
    }

    // If it's a self message, we know we are waiting for an UPDATE_ADD_HTLC to be committed before sending
    if (baseMsg->isSelfMessage()) {

        std::vector<std::string> path = baseMsg->getHops();
        std::string nextHop = path[baseMsg->getHopCount()-1];

        // Create dummy UPDATE_ADD_HTLC to use in lookup
        HTLC *tempHTLC = new HTLC();
        tempHTLC->setPaymentHash(paymentHash);
        tempHTLC->setType(UPDATE_ADD_HTLC);

        // Check if the UPDATE_ADD_HTLC has been committed
        if (!_paymentChannels[nextHop].isCommittedHTLC(tempHTLC)) {
            EV << "Waiting to send first UPDATE_FAIL_HTLC of payment " + paymentHash + ".\n";
            scheduleAt((simTime() + SimTime(500,SIMTIME_MS)),baseMsg);
        } else {
            EV << "Sending delayed first UPDATE_FAIL_HTLC for payment " + paymentHash + ".\n";

            int htlcType = UPDATE_FAIL_HTLC;
            std::string htlcId = createHTLCId(paymentHash, htlcType);

            UpdateFailHTLC *failHTLC = new UpdateFailHTLC();
            failHTLC->setHtlcId(htlcId.c_str());
            failHTLC->setPaymentHash(paymentHash.c_str());
            failHTLC->setValue(value);
            failHTLC->setErrorReason(errorReason.c_str());

            HTLC *baseHTLC = new HTLC(failHTLC);
            sendFirstFailHTLC(baseHTLC, nextHop);
        }
        return;
    }

    // Create temporary HTLC for status check
    HTLC *tempHTLC = new HTLC();
    tempHTLC->setPaymentHash(paymentHash);
    tempHTLC->setType(UPDATE_ADD_HTLC);

    EV << "Payment " + paymentHash +  "has been refused at node " + sender + ". Error reason: " + errorReason + ". Undoing updates...\n";

    if(!_paymentChannels[sender].isPendingHTLC(tempHTLC)) {
        // If we don't find the HTLC in our pending list, we look into our committed HTLCs list.
        if (!_paymentChannels[sender].isInFlight(tempHTLC)) {
            // The received HTLC is neither pending nor in flight. Something unexpected happened...
            throw std::invalid_argument( "ERROR: Unknown PAYMENT_REFUSED received!" );
        } else {
            // The HTLC is in flight. This cannot happen.
            throw std::invalid_argument( "ERROR: Refused HTLC is already in flight!" );
        }
    } else {

        int htlcType = UPDATE_ADD_HTLC;
        std::string htlcId = createHTLCId(paymentHash, htlcType);
        HTLC* addHTLC = _paymentChannels[sender].getPendingHTLC(htlcId);

        // The HTLC is still pending so we need to remove it upstream before the next commitment and trigger UPDATE_FAIL_HTLC downstream.
        _paymentChannels[sender].removePendingHTLC(htlcId);
        _paymentChannels[sender].removePendingHTLCFIFOByValue(addHTLC);
        _paymentChannels[sender].removePreviousHopUp(htlcId);

        // We should also look in our HTLCs waiting for ack (in case our payment has been refused after we sent
        // a commitment signed message with it)
        std::map<int, std::vector <HTLC *>> allHTLCsWaitingForAck = _paymentChannels[sender].getAllHTLCsWaitingForAck();

        for (const auto & pair : allHTLCsWaitingForAck) {
            int ackId = pair.first;
            std::vector<HTLC*> htlcVector = pair.second;
            for (const auto & htlc : htlcVector) {
                if (createHTLCId(htlc->getPaymentHash(), htlc->getType()) == htlcId)
                    _paymentChannels[sender].removeHTLCFromWaitingForAck(ackId, htlc);
            }
        }

        std::vector<std::string> path = baseMsg->getHops();

        // If I'm not the origin of the payment, trigger update fail downstream
        if (myName != path[0]) {
            // Store base message for retrieving path on the first fail message
            _myStoredMessages[paymentHash] = baseMsg;

            std::string nextHop = path[baseMsg->getHopCount()-1];

            // Create dummy UPDATE_ADD_HTLC to use in lookup
            HTLC *addHTLC = new HTLC();
            addHTLC->setPaymentHash(paymentHash);
            addHTLC->setType(UPDATE_ADD_HTLC);

            // If the corresponding UPDATE_ADD_HTLC has not been committed in the next downstream hop, wait to send
            // the UPDATE_FULFILL_HTLC. This way we avoid sending a fail message for a pending payment.
            if (_paymentChannels[nextHop].isCommittedHTLC(addHTLC)) {
                EV << "Waiting to send first UPDATE_FAIL_HTLC of payment " + paymentHash + ".\n";
                scheduleAt((simTime() + SimTime(500,SIMTIME_MS)),baseMsg);

            } else {
                UpdateFailHTLC *failHTLC = new UpdateFailHTLC();
                failHTLC->setHtlcId(htlcId.c_str());
                failHTLC->setPaymentHash(paymentHash.c_str());
                failHTLC->setValue(value);
                failHTLC->setErrorReason(errorReason.c_str());

                HTLC *baseHTLC = new HTLC(failHTLC);
                sendFirstFailHTLC(baseHTLC, nextHop);
            }
        }
    }
}

void FullNode::commitSignedHandler (BaseMessage *baseMsg) {
    EV << "COMMITMENT_SIGNED received at " + std::string(getName()) + " from " + std::string(baseMsg->getSenderModule()->getName()) + ".\n";
    commitmentSigned *commitMsg = check_and_cast<commitmentSigned *>(baseMsg->decapsulate());

    std::string sender = baseMsg->getSenderModule()->getName();
    HTLC *htlc = NULL;
    std::string myName = getName();
    std::string paymentHash;
    unsigned short index = 0;
    std::vector<HTLC *> HTLCs = commitMsg->getHTLCs();
    size_t numberHTLCs = HTLCs.size();
    std::vector<HTLC *> sortedHTLCs = this->getSortedPendingHTLCs(HTLCs, sender);

    // Iterate through the sorted HTLC list and attempt to commit them
    for (const auto & htlc : sortedHTLCs) {

        paymentHash = htlc->getPaymentHash();
        double value = htlc->getValue();
        int htlcType = htlc->getType();

        // Skip HTLC if it has already been committed
        if (_paymentChannels[sender].isCommittedHTLC(htlc)) {
            EV << "WARNING: Skipped " + std::to_string(htlcType) + " with paymentHash " + paymentHash + " on node " + myName + ".\n";
            continue;
        }

        switch(htlc->getType()) {
            case UPDATE_ADD_HTLC: {
                commitUpdateAddHTLC(htlc, sender);
                break;
            }
            case UPDATE_FULFILL_HTLC: {
                commitUpdateFulfillHTLC(htlc, sender);
                break;
            }
            case UPDATE_FAIL_HTLC: {
                commitUpdateFailHTLC(htlc, sender);
                break;
            }
        }
     }

    std::string signalName = myName +"-to-" + sender + ":capacity";
    emit(_signals[signalName], _paymentChannels[sender]._capacity);

    revokeAndAck *ack = new revokeAndAck();
    ack->setAckId(commitMsg->getId());
    ack->setHTLCs(sortedHTLCs);

    BaseMessage *newMessage = new BaseMessage();
    newMessage->setDestination(sender.c_str());
    newMessage->setMessageType(REVOKE_AND_ACK);
    newMessage->setHopCount(0);
    newMessage->setName("REVOKE_AND_ACK");
    //newMessage->setUpstreamDirection(!baseMsg->getUpstreamDirection());

    newMessage->encapsulate(ack);

    cGate *gate = _paymentChannels[sender].getLocalGate();

    //Sending pre image out
    EV << "Sending ack to " + sender + "with id " + std::to_string(commitMsg->getId()) + "\n";
    send(newMessage, gate);
}

void FullNode::revokeAndAckHandler (BaseMessage *baseMsg) {

    EV << "REVOKE_AND_ACK received at " + std::string(getName()) + " from " + std::string(baseMsg->getSenderModule()->getName()) + ".\n";
    revokeAndAck *ackMsg = check_and_cast<revokeAndAck *> (baseMsg->decapsulate());

    std::string sender = baseMsg->getSenderModule()->getName();
    int ackId = ackMsg->getAckId();

    std::vector<HTLC *> HTLCs = _paymentChannels[sender].getHTLCsWaitingForAck(ackId);
    HTLC *htlc;
    std::string myName = getName();
    std::string paymentHash;
    size_t index = 0;
    std::vector<HTLC *> sortedHTLCs = this->getSortedPendingHTLCs(HTLCs, sender);

    // Iterate through the sorted HTLC list and attempt to commit them
    for (const auto & htlc : sortedHTLCs) {

        paymentHash = htlc->getPaymentHash();
        double value = htlc->getValue();
        int htlcType = htlc->getType();

        if (_paymentChannels[sender].isCommittedHTLC(htlc)) {
            EV << "WARNING: Skipped " + std::to_string(htlcType) + " with paymentHash " + paymentHash + " on node " + myName + ".\n";
            continue;
        }

        switch(htlc->getType()) {
            case UPDATE_ADD_HTLC: {
                commitUpdateAddHTLC(htlc, sender);
                break;
            }
            case UPDATE_FULFILL_HTLC: {
                commitUpdateFulfillHTLC(htlc, sender);
                break;
            }
            case UPDATE_FAIL_HTLC: {
                commitUpdateFailHTLC(htlc, sender);
                break;
            }
        }
    }
    _paymentChannels[sender].removeHTLCsWaitingForAck(ackId);
    _paymentChannels[sender].setWaitingForAck(false);
}


/***********************************************************************************************************************/
/* HTLC SENDERS                                                                                                        */
/***********************************************************************************************************************/

void FullNode::sendFirstFulfillHTLC (HTLC *htlc, std::string firstHop) {
    // This function creates and sends an UPDATE_FULFILL_HTLC to the first hop in the downstream direction, triggering the beginning of payment completion

    EV << "Payment reached its destination. Releasing preimage... \n";

    //Get the stored pre image
    std::string htlcId = htlc->getHtlcId();
    std::string paymentHash = htlc->getPaymentHash();
    std::string preImage = _myPreImages[paymentHash];
    BaseMessage *storedBaseMsg = _myStoredMessages[paymentHash];
    std::vector<std::string> path = storedBaseMsg->getHops();
    std::string myName = getName();
    int htlcType = UPDATE_FULFILL_HTLC;

    //std::string htlcId = createHTLCId(paymentHash, htlcType);

    //Generate an UPDATE_FULFILL_HTLC message
    BaseMessage *newMessage = new BaseMessage();
    newMessage->setDestination(path[0].c_str());
    newMessage->setMessageType(UPDATE_FULFILL_HTLC);
    newMessage->setHopCount(storedBaseMsg->getHopCount() - 1);
    newMessage->setHops(storedBaseMsg->getHops());
    newMessage->setName("UPDATE_FULFILL_HTLC");
    newMessage->setDisplayString("i=block/decrypt;is=s");

    UpdateFulfillHTLC *firstFulfillHTLC = new UpdateFulfillHTLC();
    firstFulfillHTLC->setHtlcId(htlcId.c_str());
    firstFulfillHTLC->setPaymentHash(paymentHash.c_str());
    firstFulfillHTLC->setPreImage(preImage.c_str());
    firstFulfillHTLC->setValue(htlc->getValue());

    // Set UPDATE_FULFILL_HTLC as pending and invert the previous hop (now we're going downstream)
    HTLC *baseHTLC  = new HTLC(firstFulfillHTLC);
    _paymentChannels[firstHop].setPendingHTLC(htlcId, baseHTLC);
    _paymentChannels[firstHop].setLastPendingHTLCFIFO(baseHTLC);
    //_paymentChannels[firstHop].removePreviousHopUp(htlcId);
    _paymentChannels[firstHop].setPreviousHopDown(htlcId, myName);

    newMessage->encapsulate(firstFulfillHTLC);

    cGate *gate = _paymentChannels[firstHop].getLocalGate();

    _myPreImages.erase(paymentHash);
    _myStoredMessages.erase(paymentHash);

    //Sending HTLC out
    EV << "Sending pre image " + preImage + " to " + path[(newMessage->getHopCount()-1)] + "for payment hash " + paymentHash + "\n";
    send(newMessage, gate);
}

void FullNode::sendFirstFailHTLC (HTLC *htlc, std::string firstHop) {
    // This function creates and sends an UPDATE_FAIL_HTLC to the first hop in the downstream direction, triggering the beginning of payment failure

    //Get the stored base message
    std::string myName = getName();
    std::string htlcId = htlc->getHtlcId();
    std::string paymentHash = htlc->getPaymentHash();
    BaseMessage *storedBaseMsg = _myStoredMessages[paymentHash];
    std::vector<std::string> failPath = storedBaseMsg->getHops();
    double value = htlc->getValue();
    int htlcType = UPDATE_FAIL_HTLC;

    EV << "Initializing downstream unlocking of HTLCs for payment " + paymentHash + "... \n";

    //Generate an UPDATE_FAIL_HTLC message
    BaseMessage *newMessage = new BaseMessage();
    newMessage->setDestination(failPath[0].c_str());
    newMessage->setMessageType(UPDATE_FAIL_HTLC);
    newMessage->setHopCount(storedBaseMsg->getHopCount()-1);
    newMessage->setHops(failPath);
    newMessage->setName("UPDATE_FAIL_HTLC");
    newMessage->setDisplayString("i=status/stop");

    UpdateFailHTLC *firstFailHTLC = new UpdateFailHTLC();
    firstFailHTLC->setHtlcId(htlcId.c_str());
    firstFailHTLC->setPaymentHash(paymentHash.c_str());
    firstFailHTLC->setValue(htlc->getValue());
    firstFailHTLC->setErrorReason(htlc->getErrorReason().c_str());

    // Set UPDATE_FAIL_HTLC as pending and invert the previous hop (now we're going downstream)
    HTLC *baseHTLC  = new HTLC(firstFailHTLC);
    _paymentChannels[firstHop].setPendingHTLC(htlcId, baseHTLC);
    _paymentChannels[firstHop].setLastPendingHTLCFIFO(baseHTLC);
    _paymentChannels[firstHop].setPreviousHopDown(htlcId, myName);

    newMessage->encapsulate(firstFailHTLC);

    cGate *gate = _paymentChannels[firstHop].getLocalGate();

    _myPreImages.erase(paymentHash);
    _myStoredMessages.erase(paymentHash);

    //Sending HTLC out
    EV << "Sending first UPDATE_FAIL_HTLC to " + failPath[(newMessage->getHopCount())] + "for payment hash " + paymentHash + "\n";
    send(newMessage, gate);
}


/***********************************************************************************************************************/
/* HTLC COMMITTERS                                                                                                     */
/***********************************************************************************************************************/

void FullNode::commitUpdateAddHTLC (HTLC *htlc, std::string neighbor) {

    std::string myName = getName();
    std::string paymentHash = htlc->getPaymentHash();
    int htlcType = htlc->getType();
    std::string htlcId = htlc->getHtlcId();
    //std::string htlcId = createHTLCId(paymentHash, htlcType);
    std::string previousHop = _paymentChannels[neighbor].getPreviousHopUp(htlcId);

    EV << "Committing UPDATE_ADD_HTLC on channel " + myName + "->" + neighbor + " with payment hash " + paymentHash + "...\n";


    // If our neighbor is the HTLC's previous hop, we should commit but not set inFlight (that's the neighbors's responsibility)
    if (_paymentChannels[neighbor].getPreviousHopUp(htlcId) == neighbor) {

        // If we are the destination, commit and trigger first UPDATE_FULFILL_HTLC function
        if (!_myPreImages[paymentHash].empty()) {
            commitHTLC(htlc, neighbor);
            sendFirstFulfillHTLC(htlc, neighbor);

        // Otherwise just commit
        } else {
            commitHTLC(htlc, neighbor);
        }
    // If our neighbor is the HTLC's next hop, we must set it as in flight and decrement the channel balance
    } else if (_paymentChannels[neighbor].getPreviousHopUp(htlcId) == myName) {
        setInFlight(htlc, neighbor);
        commitHTLC(htlc, neighbor);

    // If either case is satisfied, this is unexpected behavior
    } else {
        throw std::invalid_argument("ERROR: Could not commit UPDATE_ADD_HTLC. Reason: previousHop unknown.");
    }
}

void FullNode::commitUpdateFulfillHTLC (HTLC *htlc, std::string neighbor) {

    std::string myName = getName();
    std::string htlcId = htlc->getHtlcId();
    std::string paymentHash = htlc->getPaymentHash();
    std::string previousHop = _paymentChannels[neighbor].getPreviousHopDown(htlcId);
    double value = htlc->getValue();
    int htlcType = htlc->getType();
    //std::string htlcId = createHTLCId(paymentHash, htlcType);

    EV << "Committing UPDATE_FULFILL_HTLC on channel " + myName + "->" + neighbor + " with payment hash " + paymentHash + "...\n";

    // If our neighbor is the fulfill's previous hop, we must remove the in flight HTLCs
    if (_paymentChannels[neighbor].getPreviousHopDown(htlcId) == neighbor) {
        _paymentChannels[neighbor].removeInFlight(htlcId);
        commitHTLC(htlc, neighbor);

        // If we are the destination, the payment has completed successfully
        if(_myPayments[paymentHash] == "PENDING") {
            bubble("Payment completed!");
            EV << "Payment " + paymentHash + " completed!\n";

            _myPayments[paymentHash] = "COMPLETED";
            _countCompleted++;
            _paymentGoodputSent = double(_countCompleted)/double(_countCompleted + _countFailed);
            _paymentGoodputAll = double(_countCompleted)/double(_countCompleted + _countFailed + _countCanceled);

            emit(_signals["completedPayments"], _countCompleted);
            emit(_signals["paymentGoodputSent"], _paymentGoodputSent);
            emit(_signals["paymentGoodputAll"], _paymentGoodputAll);
        }

        // If we are the fulfill's previous hop, we should claim our money
    } else if (_paymentChannels[neighbor].getPreviousHopDown(htlcId) == myName) {
        tryUpdatePaymentChannel(neighbor, value, true);
        commitHTLC(htlc, neighbor);

    // If either case is satisfied, this is unexpected behavior
    } else {
        throw std::invalid_argument("ERROR: Could not commit UPDATE_FULFILL_HTLC. Reason: previousHop unknown.");
    }
}

void FullNode::commitUpdateFailHTLC (HTLC *htlc, std::string neighbor) {

    std::string myName = getName();
    std::string htlcId = htlc->getHtlcId();
    std::string paymentHash = htlc->getPaymentHash();
    std::string previousHop = _paymentChannels[neighbor].getPreviousHopDown(htlcId);
    double value = htlc->getValue();
    int htlcType = htlc->getType();

    EV << "Committing UPDATE_FAIL_HTLC on channel " + myName + "->" + neighbor + " with payment hash " + paymentHash + "...\n";

    // If our neighbor is the fail's previous hop, we should we must remove the in flight HTLCs and claim our money back
    if (_paymentChannels[neighbor].getPreviousHopDown(htlcId) == neighbor) {
        _paymentChannels[neighbor].removeInFlight(htlcId);
        tryUpdatePaymentChannel(neighbor, value, true);
        commitHTLC(htlc, neighbor);

        EV << "Claimed HTLC back. Value: " + std::to_string(value) + "\n.";

        // If we are the destination, the payment has failed
        if(_myPayments[paymentHash] == "PENDING") {
            bubble("Payment failed!");
            EV << "Payment " + paymentHash + " failed!\n";
            _myPayments[paymentHash] = "FAILED";

            _countFailed++;
            _paymentGoodputSent = double(_countCompleted)/double(_countCompleted + _countFailed);
            _paymentGoodputAll = double(_countCompleted)/double(_countCompleted + _countFailed + _countCanceled);

            emit(_signals["failedPayments"], _countFailed);
            emit(_signals["paymentGoodputSent"], _paymentGoodputSent);
            emit(_signals["paymentGoodputAll"], _paymentGoodputAll);
        }

    // If our neighbor is the fails's next hop, just remove from pending (the updates have been applied in the sender node)
    } else if (_paymentChannels[neighbor].getPreviousHopDown(htlcId) == myName) {
        commitHTLC(htlc, neighbor);

    // If either case is satisfied, this is unexpected behavior
    } else {
        throw std::invalid_argument("ERROR: Could not commit UPDATE_FAIL_HTLC. Reason: previousHop unknown.");
    }
}

void FullNode::commitHTLC (HTLC *htlc, std::string neighbor) {
    // Removes HTLC from pending list and adds it to the commited HTLCs

    std::string htlcId = htlc->getHtlcId();
    _paymentChannels[neighbor].removePendingHTLC(htlcId);
    _paymentChannels[neighbor].removePendingHTLCFIFOByValue(htlc);
    _paymentChannels[neighbor].setCommittedHTLC(htlcId, htlc);
    _paymentChannels[neighbor].setLastCommittedHTLCFIFO(htlc);
}


/***********************************************************************************************************************/
/* STATISTICS                                                                                                          */
/***********************************************************************************************************************/

//void FullNode::initStatistics() {
//    // Register signals and statistics for producing simulaton results
//
//    // Per channel statistics
//
//    std::string signalName = myName +"-to-" + neighborName + ":capacity";
//    simsignal_t signal = registerSignal(signalName.c_str());
//    _signals[signalName] = signal;
//    emit(_signals[signalName], _paymentChannels[neighborName]._capacity);
//
//    std::string statisticName = myName +"-to-" + neighborName + ":capacity";
//    cProperty *statisticTemplate = getProperties()->get("statisticTemplate", "pcCapacities");
//    getEnvir()->addResultRecorders(this, signal, statisticName.c_str(), statisticTemplate);
//}

void FullNode::initPerModuleStatistics () {
    // Registers all statistics that do not depend on payment channels

    simsignal_t completedPayments = registerSignal("completedPayments");
    _signals["completedPayments"] = completedPayments;

    simsignal_t failedPayments = registerSignal("failedPayments");
    _signals["failedPayments"] = failedPayments;

    simsignal_t canceledPayments = registerSignal("canceledPayments");
    _signals["canceledPayments"] = canceledPayments;

    simsignal_t paymentGoodputSent = registerSignal("paymentGoodputSent");
    _signals["paymentGoodputSent"] = paymentGoodputSent;

    simsignal_t paymentGoodputAll = registerSignal("paymentGoodputAll");
    _signals["paymentGoodputAll"] = paymentGoodputAll;
}


/***********************************************************************************************************************/
/* UTIL FUNCTIONS                                                                                                      */
/***********************************************************************************************************************/

bool FullNode::tryUpdatePaymentChannel (std::string nodeName, double value, bool increase) {
    // Helper function to update payment channels. If increase = true, the function attemps to increase
    // the capacity of the channel. Else, check whether the channel has enough capacity to process the payment.

    if(increase) {
        _paymentChannels[nodeName].increaseCapacity(value);
        return true;
    } else {
        double capacity = _paymentChannels[nodeName]._capacity;
        if(capacity - value < 0)
            return false;
        else {
            _paymentChannels[nodeName].decreaseCapacity(value);
            return true;
        }
    }
}

bool FullNode::hasCapacityToForward  (std::string nodeName, double value) {
    // Helper function that calculates the payment channel capacity after applying the pending HTLCs and checks if the node has sufficient funds to forward a payment.

    std::deque< HTLC*> pendingHTLCsFIFO = _paymentChannels[nodeName].getPendingHTLCsFIFO();
    std::string myName = getName();

    // Calculate the capacity after applying pending HTLCs
    double capacity = _paymentChannels[nodeName].getCapacity();
    for (const auto & htlc : pendingHTLCsFIFO) {
        std::string htlcId  = htlc->getHtlcId();
        int htlcType = htlc->getType();

        // If it's an add update, subtract value from capacity if we are the previous hop uptstream
        // (because we'll have less money when we commmit it)
        if (htlcType == UPDATE_ADD_HTLC) {
            if (_paymentChannels[nodeName].getPreviousHopUp(htlcId) == myName) {
                capacity -= htlc->getValue();
                if (capacity <= 0)
                    return false;
            }
        } else if (htlcType == UPDATE_FAIL_HTLC) {
            // If it's a fail update, add value to capacity if we are not the previous hop downstream
            // (because we'll recover money when we commmit it)
            if (_paymentChannels[nodeName].getPreviousHopDown(htlcId) == nodeName) {
                capacity += htlc->getValue();
            }
        } else {}; // If it's a fulfill update, do nothing (fulfills don't change the capacity in the upstream direction)
    }

    // Check if the capacity would become negative after forwarding the next payment
    if ((capacity - value) <= 0)
        return false;
    else
        return true;
}

bool FullNode::tryCommitTxOrFail(std::string sender, bool timeoutFlag) {
    /***********************************************************************************************************************/
    /* tryCommitOrFail verifies whether the pending transactions queue has reached the defined commitment batch size       */
    /* or not. If the batch is full, tryCommitOrFail creates a commitment_signed message and sends it to the node that     */
    /* shares the payment channel where the batch is full. If the batch has not yet reached its limit, tryCommitOrFail     */
    /* does nothing.                                                                                                       */
    /***********************************************************************************************************************/

    unsigned short int index = 0;
    HTLC *htlc = NULL;
    std::vector<HTLC *> HTLCVector;
    std::string paymentHash;
    bool through = false;
    std::string myName = getName();


    EV << "Entered tryCommitTxOrFail. Current batch size: " + std::to_string(_paymentChannels[sender].getPendingBatchSize()) + "\n";

    if (_paymentChannels[sender].getPendingBatchSize() >= COMMITMENT_BATCH_SIZE || timeoutFlag == true) {
        for (const auto & htlc : _paymentChannels[sender].getPendingHTLCsFIFO()) {
            HTLCVector.push_back(htlc);
        }

        EV << "Setting through to true\n";
        through = true;
        _paymentChannels[sender].setWaitingForAck(true);

        commitmentSigned *commitTx = new commitmentSigned();
        commitTx->setHTLCs(HTLCVector);
        commitTx->setId(localCommitCounter);

        _paymentChannels[sender].setHTLCsWaitingForAck(localCommitCounter, HTLCVector);

        localCommitCounter += 1;
        //int gateIndex = rtable[sender];
        cGate *gate = _paymentChannels[sender].getLocalGate();

        BaseMessage *baseMsg = new BaseMessage();
        baseMsg->setDestination(sender.c_str());
        baseMsg->setMessageType(COMMITMENT_SIGNED);
        baseMsg->setHopCount(0);
        baseMsg->setName("COMMITMENT_SIGNED");

        baseMsg->encapsulate(commitTx);

        EV << "Sending Commitment Signed from node " + std::string(getName()) + "to " + sender + ". localCommitCounter: " + std::to_string(localCommitCounter) + "\n";

        send (baseMsg, gate);
    }
    return through;
}

Invoice* FullNode::generateInvoice(std::string srcName, double value) {

    std::string preImage;
    std::string preImageHash;
    preImage = generatePreImage();
    preImageHash = sha256(preImage);

    EV<< "Generated pre image " + preImage + " with hash " + preImageHash + "\n";

    _myPreImages[preImageHash] = preImage;

    Invoice *invoice = new Invoice();
    invoice->setSource(srcName.c_str());
    invoice->setDestination(getName());
    invoice->setValue(value);
    invoice->setPaymentHash(preImageHash.c_str());

    return invoice;
}

void FullNode::setInFlight(HTLC *htlc, std::string nextHop) {
    // Sets payment in flight and removes from pending

    std::string htlcId = htlc->getHtlcId();
    std::string paymentHash = htlc->getPaymentHash();
    int htlcType = htlc->getType();

    // If payment is already in flight, do nothing.
    if (_paymentChannels[nextHop].getInFlight(htlcId)) {
        EV << "Payment " + paymentHash + "already in flight! Ignoring...\n";

    // Else, try to put the HTLC in flight and decrease capacity on the forward direction
    } else {
        if (!tryUpdatePaymentChannel(nextHop, htlc->getValue(), false)) {
            // Insufficient funds. Trigger UPDATE_FAIL_HTLC
            throw std::invalid_argument("ERROR: Could not commit UPDATE_ADD_HTLC. Reason: Insufficient funds.");
        }
        _paymentChannels[nextHop].setInFlight(htlcId, htlc);
        EV << "Payment hash " + paymentHash + " set in flight.\n";
    }
}

bool FullNode::isInFlight(HTLC *htlc, std::string nextHop) {
    // Checks if payment is already in flight

    if(!_paymentChannels[nextHop].getInFlight(htlc->getHtlcId()))
        return false;
    else
        return true;
}

std::vector <HTLC *> FullNode::getSortedPendingHTLCs (std::vector<HTLC *> HTLCs, std::string neighbor) {
    // Util function that receies a vector of HTLCs and sorts them according to the local order
    // (also discards HTLCs that are not in the pending list)
    std::deque<HTLC *> pendingHTLCsFIFO = _paymentChannels[neighbor].getPendingHTLCsFIFO();
    std::vector<HTLC *> sortedHTLCs;

    for (const auto & pendingHTLC : pendingHTLCsFIFO) {
        for (const auto & htlc : HTLCs) {
            std::string htlcId = htlc->getHtlcId();
           if (htlcId == pendingHTLC->getHtlcId()) {
               sortedHTLCs.push_back(htlc);
           }
        }
    }
    return sortedHTLCs;
}

std::string FullNode::createHTLCId (std::string paymentHash, int htlcType) {
    return paymentHash + ":" + std::to_string(htlcType);
}

