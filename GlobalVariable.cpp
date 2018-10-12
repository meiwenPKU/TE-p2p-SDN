/*
 * GlobalVariable.cpp
 *
 *  Created on: Apr 21, 2016
 *      Author: hemin
 */
#include "GlobalVariable.h"


const double NOPATH = -10;
const int N_AS = 3 ; // the number of AS's in the topology
const int kPath = 3; // the maximum paths between a source and a sink in a topology table
const int kK = 3; //k shortest paths used to initialize topo table
const int kMaxComputePath = 3;
double prob_generate_traffic = 0.5; // the prob of a node generating traffic
double prob_within_AS = 0.1; // the prob of the destination of a connection is within the AS
double prob_Google = 0.5;
const double delta = 0.0001;
const double multiplier = 1;
const double loadC = 3.99;

int my_commodity_demand[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
//int my_commodity_demand[] = {10, 20};
vector<int> commodity_demand (my_commodity_demand,my_commodity_demand + sizeof(my_commodity_demand)/sizeof(int)); // the possible demands of commodities
double my_prob_demand[] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
//double my_prob_demand[] = {0.5,0.5};
vector<double> prob_demand (my_prob_demand, my_prob_demand + sizeof(my_prob_demand)/sizeof(double)); //the prob of each demand
