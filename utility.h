#pragma once
#include <set>
#include <map>
#include <queue>
#include <string>
#include <cmath>
#include <sstream>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <memory>
#include <unistd.h>
#include "GraphElements.h"
#include "Graph.h"
#include "Commodity.h"
#include "DijkstraShortestPathAlg.h"
#include "YenTopKShortestPathsAlg.h"
#include "FordFulkersonAlg.h"
//hemin
#include "TopoTable.h"
#include "GlobalVariable.h"
using namespace std;

struct setcomp
{
  bool operator()(const pair<BasePath,double>& lhs, const pair<BasePath,double>& rhs)const
  {
    if (lhs.first == rhs.first)
    {
      return false;
    }
    else
    {
      return true;
    }
  }
};

void process_mem_usage(double& vm_usage, double& resident_set);

void GenerateCommodity(vector<Graph*> ASes);

/*
 * map the vertex in the network to the vertex in one AS
 * if the graph is constructed based on several ASes,
 * this map is used to map the node in the AS to the node id of the constructed graph.
 * mpNodeID[graph_id*max_Nodes+vertex_id] = new vertex_id in the constructed graph
 */
BaseVertex* NetworkToAS(Graph& network, vector<Graph*> ASes, BaseVertex* ori_vertex);

double BWFunction(Commodity* commodity, double fair_share);

/*
 * select the value according to the prob distribution
 *
 */
int SelectValue(vector<int>& value, vector<double>& prob_dist);

void testDijkstraGraph(Graph* my_graph_pt);

void testYenAlg(Graph& my_graph);

/*
 * 1. find the shortest path for each commodity
 * 2. find the shortest path among all the paths find in step 1
 * 3. serve this path with the capacity of this path
 */
void Max_Throughput_TE(Graph* AS, map<Commodity*, set<pair<BasePath,double>, setcomp> >* Allocation);


/*
 * Google's TE optimization algorithm, please refer to S. Jain etc [B4]
 * 1. find the shortest path for each commodity
 * 2. find the bottleneck edge (with minimum fair share at its capacity)
 * 3. remove the bottleneck edge, repeat 1~2 until no more bottleneck edge is found or each commodity is satisfied
 */
void Google_TE_Optimization(Graph* AS, map<Commodity*, set<pair<BasePath,double>, setcomp> >* Allocation);

/*
 * Google's TE optimization algorithm, please refer to S. Jain etc [B4]
 * 1. find the shortest path for each commodity
 * 2. find the bottleneck edge (with minimum fair share at its capacity)
 * 3. remove the bottleneck edge, repeat 1~2 until no more bottleneck edge is found or each commodity is satisfied
 */
void Google_TE_Optimization_Benchmark(Graph* AS, map<Commodity*, set<pair<BasePath,double>, setcomp> >* Allocation);
