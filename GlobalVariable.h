#pragma once
#include <vector>
using namespace std;

extern const double NOPATH;
//extern const int N_AS; // the number of AS's in the topology

// here we have three const parameters related with paths.
// 1) kK, which is used in the YenTopKShortestPathsAlg to initialize the topo table
// 2) kPath, which is the maximum number of paths for one source and one sink in the topo table
// 3) kMaxComputePath, which is the maximum number of paths used to compute the aggregated path from
//    one source to one sink when updating the AdvertisedTable
// extern const int kPath;
extern const int kK;
//extern const int kMaxComputePath;

extern double prob_generate_traffic; // the prob of a node generating traffic
extern double prob_within_AS; // the prob of the destination of a connection is within the AS
extern double prob_Google; // the prob of an AS using google's TE algorithm
extern const double delta;
extern const double multiplier; // multiply each link's bw
//extern const double loadC; // each commodity's demand multiply with this coefficient

extern int my_commodity_demand[];
extern vector<int> commodity_demand; // the possible demands of commodities
extern double my_prob_demand[];
extern vector<double> prob_demand; //the prob of each demand
