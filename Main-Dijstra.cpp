/************************************************************************/
/* This script is used to measure the scenario of applying Dijkstra
/************************************************************************/

#include <limits>
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

void testDijkstraGraph(Graph* my_graph_pt)
{
  DijkstraShortestPathAlg shortest_path_alg(my_graph_pt);
  BasePath* result =
    shortest_path_alg.get_shortest_path(
      my_graph_pt->get_vertex(6,my_graph_pt->get_graphID()), my_graph_pt->get_vertex(3,my_graph_pt->get_graphID()));
  result->PrintOut(cout);
  result =
    shortest_path_alg.get_shortest_path(
      my_graph_pt->get_vertex(10,my_graph_pt->get_graphID()), my_graph_pt->get_vertex(5,my_graph_pt->get_graphID()));
  result->PrintOut(cout);
  result =
    shortest_path_alg.get_shortest_path(
      my_graph_pt->get_vertex(2,my_graph_pt->get_graphID()), my_graph_pt->get_vertex(0,my_graph_pt->get_graphID()));
  result->PrintOut(cout);
}

void testYenAlg(Graph& my_graph)
{
  YenTopKShortestPathsAlg yenAlg(my_graph, my_graph.get_vertex(10,my_graph.get_graphID()),
                                 my_graph.get_vertex(5,my_graph.get_graphID()));

  int i=0;
  while(yenAlg.has_next())
  {
    ++i;
    yenAlg.next()->PrintOut(cout);
  }

// 	System.out.println("Result # :"+i);
// 	System.out.println("Candidate # :"+yenAlg.get_cadidate_size());
// 	System.out.println("All generated : "+yenAlg.get_generated_path_size());

}

int Graph::Graph_ID = 0;

/*
 * select the value according to the prob distribution
 *
 */

int SelectValue(vector<int>& value, vector<double>& prob_dist)
{
  double prob = (rand()%101)/100.0;
  double sum = 0;
  int index = 0;
  for (vector<double>::iterator it = prob_dist.begin(); it != prob_dist.end(); ++it)
  {
    if (prob >= sum && prob <= sum+*it)
    {
      break;
    }
    sum += *it;
    index++;
  }
  return index;
}

/*
 * \brief find the shortest path for each commodity, and transmit these commodities
 * \param network
 * \param v_commodity: the vector of commodities
 */
void Dijkstra(Graph* network, double& Total_Thr, double& max_cost)
{
  DijkstraShortestPathAlg shortest_path_alg(network);
  for (vector<Commodity*>::iterator it = network->m_vCommodity.begin(); it != network->m_vCommodity.end(); ++it)
  {
    BasePath* result = shortest_path_alg.get_shortest_path((*it)->source_,(*it)->sink_);
    if (result->length() > 0)
    {
      //---find the bottleneck link in the path---
      double BW_path = network->get_edge_BW(result->GetVertex(0),result->GetVertex(1))-
                       network->get_edge_UsedBW(result->GetVertex(0), result->GetVertex(1));
      if (result->Weight() > max_cost )
      {
        max_cost = result->Weight();
      }
      for (int i = 1; i < result->length()-1; ++i)
      {
        int edge_BW = network->get_edge_BW(result->GetVertex(i),result->GetVertex(i+1))
                      - network->get_edge_UsedBW(result->GetVertex(i), result->GetVertex(i+1));
        if (BW_path < edge_BW)
        {
          BW_path = edge_BW;
        }
      }
      //---transmit the traffic along this path
      double Thr = min(BW_path,double((*it)->demand_));
      Total_Thr += Thr;
      for (int i = 0; i < result->length()-1; ++i)
      {
        double usedBW = network->get_edge_UsedBW(result->GetVertex(i), result->GetVertex(i+1))+Thr;
        network->set_edge_UsedBW(result->GetVertex(i),result->GetVertex(i+1),usedBW);
        if (usedBW + delta > network->get_edge_BW(result->GetVertex(i), result->GetVertex(i+1)))
        {
          pair<int, int> removed_edge(network->get_vertex_code(result->GetVertex(i)),
                                      network->get_vertex_code(result->GetVertex(i+1)));
          network->remove_edge(removed_edge);
        }
      }
      cout << "Throughput: " << Thr << endl;
      result->PrintOut(cout);
    }
  }
  cout << "Total throughput: " << Total_Thr << endl;
  cout << "Maximal cost: " << max_cost << endl;
}

/*
 * generate commodities in the network
 * Different from v1, this version will generate a commodity between any two nodes in the network with a probability
 */
void GenerateCommodityV2(Graph* network)
{
  cout << "commodities: source node --> sink node: demand" << endl;
  srand(0);
  int num_Nodes = network->get_vertex_num() / N_AS; // number of switches in one AS
  for (int i = 0; i < network->get_vertex_num(); i++)
  {
    for (int j = 0; j< network->get_vertex_num(); j++)
    {
      if (i==j)
      {
        continue;
      }
      int demand = commodity_demand[SelectValue(commodity_demand,prob_demand)];
      int source_vertex_id = i % num_Nodes;
      int des_vertex_id = j % num_Nodes;
      int source_AS = i / num_Nodes;
      int des_AS = j / num_Nodes;
      Commodity* pt_commodity = new Commodity(network->get_vertex(source_vertex_id,source_AS),
                                              network->get_vertex(des_vertex_id,des_AS),demand);
      network->m_vCommodity.push_back(pt_commodity);
    }
  }
}

/*
 * generate commodities in the network
 * For every node, it generates a commodity with prob_generate_traffic
 * If it generates a commodity, its destination can be either in the same AS or different AS
 */
void GenerateCommodityV1(Graph* network, int seed)
{
  cout << "commodities: source node --> sink node: demand" << endl;
  srand(0);
  int num_Nodes = network->get_vertex_num()/N_AS; // number of switches in one AS
  for (int i = 0; i < network->get_vertex_num(); i++)
  {
    if (rand()%101/100.0 < prob_generate_traffic)
    {
      int source_vertex_id = i % num_Nodes;
      int source_AS = i / num_Nodes;
      int des_vertex_id, des_AS;
      int demand = commodity_demand[SelectValue(commodity_demand,prob_demand)];
      if (rand()%101/100.0 < prob_within_AS)
      {
        // the destination is in the same AS
        des_vertex_id = rand()% num_Nodes;
        while (des_vertex_id == source_vertex_id)
        {
          des_vertex_id = rand() % num_Nodes;
        }
        des_AS = source_AS;
      }
      else
      {
        // the destination is in another AS
        // select a random AS
        des_AS = rand()%N_AS;
        while (des_AS == source_AS)
        {
          des_AS = rand() % N_AS;
        }
        des_vertex_id = rand()%num_Nodes;
      }
      int des_vertex_global_id = des_vertex_id + num_Nodes * des_AS;

      Commodity* pt_commodity = new Commodity(network->get_vertex(i,source_AS),
                                              network->get_vertex(des_vertex_global_id,des_AS),demand);
      network->m_vCommodity.push_back(pt_commodity);
      cout << "(" << source_AS << ", " << source_vertex_id << ") --> (" << des_AS << ", " << des_vertex_id << "): " << demand << endl;
    }
  }
}


int main(...)
{
  //--------create the topology-----------------
  string file_name = "data/test_10AS";
  // build the entire network
  Graph* network = new Graph(file_name, true);

  int time = 0;
  int runTime = 1;
  double Throughput = 0;
  double Aver_cost = 0;
  for (time = 0; time < runTime; time++)
  {
    /*
     * initialize the network
     */
    for(vector<Commodity*>::iterator it_com = network->m_vCommodity.begin();
        it_com != network->m_vCommodity.end(); ++it_com)
    {
      delete (*it_com);
    }
    //generate commodities
    cout << "time = " << time << endl;
    network->m_vCommodity.clear();
    GenerateCommodityV1(network, time);
    for (map<int,double>::iterator it_used = network->m_mpEdgeCodeUsedBW.begin();
         it_used != network->m_mpEdgeCodeUsedBW.end(); ++it_used)
    {
      network->set_edge_UsedBW(it_used->first,0);
    }
    Dijkstra(network, Throughput, Aver_cost);
  }
}
