#pragma once
#include <vector>
using namespace std;

extern const double NOPATH;
extern const int N_AS; // the number of AS's in the topology
extern double prob_generate_traffic; // the prob of a node generating traffic
extern double prob_within_AS; // the prob of the destination of a connection is within the AS
extern double prob_Google; // the prob of an AS using google's TE algorithm
extern const double delta;

extern int my_commodity_demand[];
extern vector<int> commodity_demand; // the possible demands of commodities
extern double my_prob_demand[];
extern vector<double> prob_demand; //the prob of each demand
