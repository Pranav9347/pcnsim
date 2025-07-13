#ifndef __PCNSIM_FULLNODE_H
#define __PCNSIM_FULLNODE_H

#include <omnetpp.h>
#include <map>
#include <set>
#include <vector>
#include <tuple>
#include <string>

#include "PaymentChannel.h"
#include "clustering_m.h"
#include "clusterGraph.h"
#pragma once

#include <string>
#include <vector>

// ClusterInfo struct for cluster metrics exchange
struct ClusterInfo {
    std::string sender;
    std::string clusterID;
    std::string CIP;
    int numNodes;
    int numEdges;
    int mergedNumNodes;
    int mergedNumEdges;

    ClusterInfo() {}

    ClusterInfo(std::string _sender, std::string _clusterID, std::string _cip,
                int _n, int _e, int _mN, int _mE)
        : sender(_sender), clusterID(_clusterID), CIP(_cip),
          numNodes(_n), numEdges(_e), mergedNumNodes(_mN), mergedNumEdges(_mE) {}
};

// PathMetrics for evaluating candidate cluster paths
struct PathMetrics {
    std::vector<std::string> clusterPath;
    double totalFee;
    double minCapacity;
    double totalQuality;

    PathMetrics(std::vector<std::string> path, double fee, double cap, double quality)
        : clusterPath(std::move(path)), totalFee(fee), minCapacity(cap), totalQuality(quality) {}
};


using namespace omnetpp;
class FullNode : public cSimpleModule {

    protected:
        // Protected data structures
        std::map<std::string, std::string> _myPreImages; // paymentHash to preImage
        std::map<std::string, std::string> _myInFlights; // paymentHash to nodeName (who owes me)
        std::map<std::string, std::string> _myPayments; //paymentHash to status (PENDING, COMPLETED, FAILED, or CANCELED)
        std::map<std::string, BaseMessage *> _myStoredMessages; // paymentHash to baseMsg (for finding reverse path)
        std::map<std::string, cModule*> _senderModules; // paymentHash to Module

        // Omnetpp functions
        virtual void initialize() override;
        virtual void handleMessage(cMessage *msg) override;
        virtual void refreshDisplay() const override;
        virtual void finish() override;

        // Routing functions
        virtual std::vector<std::string> getPath(std::map<std::string, std::string> parents, std::string target);
        virtual std::string minDistanceNode (std::map<std::string, double> distances, std::map<std::string, bool> visited);
        virtual std::vector<std::string> dijkstraWeightedShortestPath (std::string src, std::string target, std::map<std::string, std::vector<std::pair<std::string, std::vector<double> > > > graph);

        // Message handlers
        virtual void initHandler (BaseMessage *baseMsg);
        virtual void invoiceHandler (BaseMessage *baseMsg);
        virtual void updateAddHTLCHandler (BaseMessage *baseMsg);
        virtual void updateFulfillHTLCHandler (BaseMessage *baseMsg);
        virtual void updateFailHTLCHandler (BaseMessage *baseMsg);
        virtual void paymentRefusedHandler (BaseMessage *baseMsg);
        virtual void commitSignedHandler (BaseMessage *baseMsg);
        virtual void revokeAndAckHandler (BaseMessage *baseMsg);

        // HTLC senders
        virtual void sendFirstFulfillHTLC (HTLC *htlc, std::string firstHop);
        virtual void sendFirstFailHTLC (HTLC *htlc, std::string firstHop);

        // HTLC committers
        virtual void commitUpdateAddHTLC (HTLC *htlc, std::string neighbor);
        virtual void commitUpdateFulfillHTLC (HTLC *htlc, std::string neighbor);
        virtual void commitUpdateFailHTLC (HTLC *htlc, std::string neighbor);
        virtual void commitHTLC(HTLC *htlc, std::string neighbor);

        // Statistics
        // virtual void initStatistics();
        virtual void initPerModuleStatistics();

        // Util functions
        virtual bool tryUpdatePaymentChannel (std::string nodeName, double value, bool increase);
        virtual bool hasCapacityToForward (std::string nodeName, double value);
        virtual bool tryCommitTxOrFail (std::string, bool);
        virtual Invoice* generateInvoice (std::string srcName, double value);
        virtual void setInFlight (HTLC *htlc, std::string nextHop);
        virtual bool isInFlight (HTLC *htlc, std::string nextHop);
        virtual std::vector <HTLC *> getSortedPendingHTLCs (std::vector<HTLC *> HTLCs, std::string neighbor);
        virtual std::string createHTLCId (std::string paymentHash, int htlcType);

    public:
        // Public data structures
        bool _isFirstSelfMessage;
        cTopology *_localTopology;
        int localCommitCounter;
        typedef std::map<std::string, int> RoutingTable;  // neighborName to gateIndex
        RoutingTable rtable;
        std::map<std::string, PaymentChannel> _paymentChannels; // neighborName to PaymentChannel
        std::map<std::string, int> _signals; // myName to signal

        // Statistic-related variables
        int _countCompleted = 0;
        int _countFailed = 0;
        int _countCanceled = 0;
        double _paymentGoodputSent = 0;
        double _paymentGoodputAll = 0;
        // --- Cluster-Based Routing Protocol ---
int countClusterEdges(const std::set<std::string>& members);
void runClusteringRound();
void sendClusterInfoTo(const std::string& neighbor);
ClusterInfo waitForClusterInfoFrom(const std::string& neighbor);
double computeDeltaS(const ClusterInfo& local, const ClusterInfo& remote);
void sendMergeRequest(const std::string& targetCIP);
void handleMergeRequest(MergeRequestMsg* msg);
void handleMergeAck(MergeAckMsg* msg);
void handleMergeNack(MergeNackMsg* msg);
bool receivedMutualMergeRequest(const std::string& targetCIP);
void mergeClusters(const std::string& a, const std::string& b, const std::string& mergedID);
std::string getCIP(const std::string& clusterID);

std::vector<std::string> runCRPRouting(const std::string& srcNode, const std::string& dstNode, double paymentAmount);
std::vector<std::string> intraClusterDijkstra(const std::string& src, const std::string& dst);
std::vector<std::vector<std::string>> yensKShortestPaths(const ClusterGraph& graph, const std::string& srcCluster, const std::string& dstCluster, int K = 3);
std::vector<PathMetrics> probeCandidatePaths(const std::vector<std::vector<std::string>>& clusterPaths, double paymentAmount);
PathMetrics selectOptimalRouteLP(const std::vector<PathMetrics>& candidatePaths, double paymentAmount);

// --- Clustering State ---
std::string clusterID;
std::string CIP;
std::set<std::string> clusterMembers;
std::set<std::string> intraClusterNeighbors;
std::set<std::string> interClusterNeighbors;
std::map<std::string, std::set<std::string>> knownClusterMembers;
std::map<std::string, ClusterInfoMsg> receivedClusterInfos;
std::set<std::string> mergeResponses;

// --- Static CRP Globals ---
static std::map<std::string, std::string> nodeToCluster;
static ClusterGraph clusterGraph;


};
#endif