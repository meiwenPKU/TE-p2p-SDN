/************************************************************************/
/* This script is used to measure the performance of naive exchange protocol
 * In this naive protocol, all k-best routes info are advertised to the neighbors
 * And we will evaulate the computation overhead, msg overhead, and memory overhead
 * compared with our sdni-TE protocol
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

double BWFunction(Commodity* commodity, double fair_share)
{
  if (fair_share > commodity->demand_)
    return commodity->demand_;
  else
    return fair_share;
}

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
 * 1. find the shortest path for each commodity
 * 2. find the shortest path among all the paths find in step 1
 * 3. serve this path with the capacity of this path
 */
void Max_Throughput_TE(Graph* AS, map<Commodity*, set<pair<BasePath,double>, setcomp> >* Allocation)
{
  DijkstraShortestPathAlg shortest_path_alg(AS);
  double current_fair_share = 0;
  int num_iteration = 0;
  BasePath* shortest_path;
  Commodity* target_com;
  double min_cost = 1000000000;
  double Thr;


  while (true)
  {
    /*
     * find the shortest path (with the minimal weight) for each commodity
     */
    min_cost = 1000000000;
    for (vector<Commodity*>::iterator it = AS->m_vCommodity.begin(); it != AS->m_vCommodity.end(); ++it)
    {
      if ((*it)->Allocated_ == NOPATH ||
          (*it)->Allocated_ == (*it)->demand_ || (*it)->isSaturated_)
      {
        continue;
      }
      BasePath* result = shortest_path_alg.get_shortest_path((*it)->source_, (*it)->sink_);
      if (result->length() > 0)
      {
        // if the path only contains one vertex
        if (result->GetVertex(0) == result->GetLastVertex())
        {
          cout << "wrong in Max_Throughput_TE(): the path only contains one node" << endl;
        }
        else
        {
          if (min_cost > result->Weight())
          {
            min_cost = result->Weight();
            shortest_path = result;
            target_com = *it;
          }
        }
      }
      else
      {
        if ((*it)->Allocated_==0)
        {
          (*it)->Allocated_ = NOPATH;
        }
        else
        {
          (*it)->isSaturated_ = true;
        }
      }
    }
    /*
    		if (target_com->source_->getID() == 0 && target_com->sink_->getID()==271)
    		{
    			cout << "debug here" << endl;
    		}
    */

    if (min_cost == 1000000000)
    {
      break;
    }

    //determine the throughput of this path
    Thr = min(AS->get_path_BW(shortest_path),target_com->demand_ - target_com->Allocated_);
    if (Thr < 0.1)
    {
      target_com->isSaturated_ = true;
      continue;
    }

    //transmit the shortest path
    for (int i = 0; i < shortest_path->length() -1; ++i)
    {
      double newUsedBW = AS->get_edge_UsedBW(shortest_path->GetVertex(i), shortest_path->GetVertex(i+1)) + Thr;
      AS->set_edge_UsedBW(shortest_path->GetVertex(i), shortest_path->GetVertex(i+1), newUsedBW);
      if (newUsedBW + delta > AS->get_edge_BW(shortest_path->GetVertex(i), shortest_path->GetVertex(i+1)))
      {
        pair<int, int> removed_edge(AS->get_vertex_code(shortest_path->GetVertex(i)),
                                    AS->get_vertex_code(shortest_path->GetVertex(i+1)));
        AS->remove_edge(removed_edge);
      }
    }

    /*
     * update the allocated bandwidth for each commodity
     */
    target_com->Allocated_ += Thr;

    /*
     * update the return value
     */
    // check whether this path is already in the Allocation
    bool isFind = false;
    for (set<pair<BasePath,double> >::iterator it_set = Allocation->at(target_com).begin();
         it_set != Allocation->at(target_com).end(); ++it_set)
    {
      if (it_set->first == *shortest_path)
      {
        isFind = true;
        BasePath temp_path = it_set->first;
        double oldAllocated = it_set->second;
        double extraAllocated = Thr;
        Allocation->at(target_com).erase(it_set);
        Allocation->at(target_com).insert(pair<BasePath,double>(temp_path,oldAllocated + extraAllocated));
        break;
      }
    }
    if (!isFind)
    {
      double extraAllocated = Thr;
      Allocation->at(target_com).insert(pair<BasePath,double>(*shortest_path,extraAllocated));
    }
    num_iteration++;
  }
  cout << "iteration number = " << num_iteration << endl;
}


/*
 * Google's TE optimization algorithm, please refer to S. Jain etc [B4]
 * 1. find the shortest path for each commodity
 * 2. find the bottleneck edge (with minimum fair share at its capacity)
 * 3. remove the bottleneck edge, repeat 1~2 until no more bottleneck edge is found or each commodity is satisfied
 */
void Google_TE_Optimization(Graph* AS, map<Commodity*, set<pair<BasePath,double>, setcomp> >* Allocation)
{
  DijkstraShortestPathAlg shortest_path_alg(AS);
  double current_fair_share = 0;


  while (true)
  {
    vector<BasePath*> v_prefer_path;
    map<BasePath*, Commodity*> mp_PathToCommodity;

    /*
     * find the shortest path (with the minimal weight) for each commodity
     */
    for (vector<Commodity*>::iterator it = AS->m_vCommodity.begin(); it != AS->m_vCommodity.end(); ++it)
    {
      if ((*it)->Allocated_ == NOPATH ||
          (*it)->Allocated_ == (*it)->demand_ || (*it)->isSaturated_)
      {
        continue;
      }


      /*
       * for the sink d, ....->a1->a2->d, if a1 and a2 are out of the AS, remove the edge a1->a2
       */
      int sinkNodeID = (*it)->sink_->getID();
      int oriSinkNodeID,oriSinkASID;
      AS->get_original_id(sinkNodeID,oriSinkNodeID,oriSinkASID);
      set<pair<int,int> > RemovedEdge;
      if (oriSinkASID != AS->get_graphID())
      {
        set<BaseVertex*> vertex_set;
        AS->get_precedent_vertices((*it)->sink_,vertex_set);
        for (set<BaseVertex*>::iterator it_pre = vertex_set.begin(); it_pre != vertex_set.end(); ++it_pre)
        {
          int preNodeID = (*it_pre)->getID();
          int preOriASID, preOriNodeID;
          AS->get_original_id(preNodeID,preOriNodeID,preOriASID);
          if (preOriASID != AS->get_graphID())
          {
            set<BaseVertex*> pre_vertex_set;
            AS->get_precedent_vertices(*it_pre,pre_vertex_set);
            for (set<BaseVertex*>::iterator it_pre2 = pre_vertex_set.begin(); it_pre2 != pre_vertex_set.end(); ++it_pre2)
            {
              int pre2NodeID = (*it_pre2)->getID();
              int pre2OriASID, pre2OriNodeID;
              AS->get_original_id(pre2NodeID,pre2OriNodeID,pre2OriASID);
              if (pre2OriASID != AS->get_graphID())
              {
                RemovedEdge.insert(pair<int,int>(pre2NodeID,preNodeID));
                AS->remove_edge(pair<int,int>(pre2NodeID,preNodeID));
              }
            }
          }
        }
      }

      BasePath* result = shortest_path_alg.get_shortest_path((*it)->source_, (*it)->sink_);

      //recover the removed edges
      for (set<pair<int,int> >::iterator  it_remove = RemovedEdge.begin(); it_remove != RemovedEdge.end(); ++it_remove)
      {
        AS->recover_removed_edge(*it_remove);
      }

      if (result->length() > 0)
      {
        // if the path only contains one vertex
        if (result->GetVertex(0) == result->GetLastVertex())
        {
          cout << "wrong in Google_TE_Optimization(): the path only contains one node" << endl;
        }
        else
        {
          v_prefer_path.push_back(result);
          mp_PathToCommodity.insert(pair<BasePath*, Commodity*>(result,*it));
          for (int i = 0; i < result->length()-1; ++i)
          {
            int edgeCode = AS->get_edge_code(result->GetVertex(i), result->GetVertex(i+1));
            (AS->m_mpEdgeCodeCommodity[edgeCode])->insert((*it));
          }
        }
      }
      else
      {
        if ((*it)->Allocated_==0)
        {
          (*it)->Allocated_ = NOPATH;
        }
        else
        {
          (*it)->isSaturated_ = true;
        }
      }
    }

    //cout << "prefer paths size = " << v_prefer_path.size() << endl;
    if (v_prefer_path.size()==0)
    {
      break;
    }

    /*
     * for each edge i, we find the minimal fair share which makes the edge saturated.
     *      i.e.,          sum{j=1 to m}{fj(xi)-fj(s)} = ci
     * where xi is the fair share of edge i, m is the number of commodities carried by edge i,
     * fj(x) is the bandwidth function for commodity j, s is the current fair share level of edge i
     * ci is the left capacity of edge i
     *
     * Then the bottleneck of the graph is the edge with the minimal xi.
     *
     * In this implementation, we set fj(x)=min{x,dj}, where dj is the demand of the commodity of j
     * To solve this problem, we use the following algorithms:
     * for each edge: c= ci+ sum{j=1 to m}{fj(s)}
     * 1) compute N[i] = the number of commodities have demand of di
     * 2) Initilize F(0) = (0,m), where F(j)[0]+F(j)[1]*x = sum{j=1 to m}{fj(x)}, d(j-1) <= x < dj
     * 3) F(j)[0] = F(j-1)[0] + N[j-1]*d(j-1), F(j)[1] = F(j-1)[1]-N[j-1]
     * 4) Find the maximal j such that F(j)<=c and F(j+1)>c
     * 5) x = (c-F(j)[0])/F(j)[1]
     * x is the minimal fair share making this edge saturated.
     */
    double min_fair_share = max_BW;
    bool findBottle = false;
    for (map<int, set<Commodity*>*>::iterator it = AS->m_mpEdgeCodeCommodity.begin();
         it != AS->m_mpEdgeCodeCommodity.end(); ++it)
    {
      if (it->second->size() == 0)
      {
        continue;
      }
      map<double, int> mp_demand;
      double c = AS->get_edge_BW(it->first) - AS->get_edge_UsedBW(it->first); // left capacity of this edge
      for (set<Commodity*>::iterator it_commodity = (*(it->second)).begin();
           it_commodity != (*(it->second)).end(); ++it_commodity)
      {
        c += BWFunction((*it_commodity),current_fair_share);
        if (mp_demand.find((*it_commodity)->demand_) != mp_demand.end())
        {
          mp_demand.at((*it_commodity)->demand_) += 1;
        }
        else
        {
          mp_demand.insert(pair<double,int>((*it_commodity)->demand_, 1));
        }
      }

      double d_less = 0; //F(j)[0], the sum of the demands of commodities, where the demand<=fair share level
      int N_large = it->second->size();// F(j)[1], the num of the commodities, whose demand>fair share level
      //int N_large = AS->m_vCommodity.size(); // F(j)[1]

      for (map<double, int>::iterator it_map = mp_demand.begin(); it_map != mp_demand.end(); ++it_map)
      {
        double share;
        if (d_less + N_large * it_map->first >= c)
        {
          share = (c - d_less)/N_large;
        }
        else
        {
          d_less += it_map->first * it_map->second;
          N_large -= it_map->second;
          share = it_map->first;
        }
        if (share < min_fair_share)
        {
          findBottle = true;
          min_fair_share = share;
        }
      }
    }

    /*
     * According to the minimal fair share, update the used bandwidth of each edge
     */
    for (vector<BasePath*>::iterator it = v_prefer_path.begin(); it != v_prefer_path.end(); ++it)
    {
      for (int i = 0; i < (*it)->length()-1; ++i)
      {
        double newUsedBW = AS->get_edge_UsedBW((*it)->GetVertex(i), (*it)->GetVertex(i+1))
                           + BWFunction(mp_PathToCommodity.at(*it),min_fair_share)
                           - BWFunction(mp_PathToCommodity.at(*it),current_fair_share);
        AS->set_edge_UsedBW((*it)->GetVertex(i), (*it)->GetVertex(i+1), newUsedBW);
        if (newUsedBW + delta > AS->get_edge_BW((*it)->GetVertex(i), (*it)->GetVertex(i+1)))
        {
          pair<int, int> removed_edge(AS->get_vertex_code((*it)->GetVertex(i)), AS->get_vertex_code((*it)->GetVertex(i+1)));
          AS->remove_edge(removed_edge);
        }
      }
    }

    /*
     * update the allocated bandwidth for each commodity
     */
    int numSatisfiedCommodity = 0;
    for (vector<Commodity*>::iterator it = AS->m_vCommodity.begin(); it != AS->m_vCommodity.end(); ++it)
    {
      if ((*it)->Allocated_ == NOPATH ||
          (*it)->Allocated_ >= (*it)->demand_ || (*it)->isSaturated_)
      {
        numSatisfiedCommodity ++;
        continue;
      }
      (*it)->Allocated_ += BWFunction((*it),min_fair_share) - BWFunction((*it),current_fair_share);
      if ((*it)->Allocated_ >= (*it)->demand_)
      {
        numSatisfiedCommodity ++;
      }
    }

    /*
     * update the return value
     */
    for (vector<BasePath*>::iterator it = v_prefer_path.begin(); it != v_prefer_path.end(); ++it)
    {
      // check whether this path is already in the Allocation
      bool isFind = false;
      for (set<pair<BasePath,double> >::iterator it_set = Allocation->at(mp_PathToCommodity.at(*it)).begin();
           it_set != Allocation->at(mp_PathToCommodity.at(*it)).end(); ++it_set)
      {
        if (it_set->first == *(*it))
        {
          isFind = true;
          BasePath temp_path = it_set->first;
          double oldAllocated = it_set->second;
          double extraAllocated = BWFunction(mp_PathToCommodity.at(*it),min_fair_share)
                                  - BWFunction(mp_PathToCommodity.at(*it),current_fair_share);
          Allocation->at(mp_PathToCommodity.at(*it)).erase(it_set);
          Allocation->at(mp_PathToCommodity.at(*it)).insert(pair<BasePath,double>(temp_path,oldAllocated + extraAllocated));
          break;
        }
      }
      if (!isFind)
      {
        double extraAllocated = BWFunction(mp_PathToCommodity.at(*it),min_fair_share)
                                - BWFunction(mp_PathToCommodity.at(*it),current_fair_share);
        Allocation->at(mp_PathToCommodity.at(*it)).insert(pair<BasePath,double>(*(*it),extraAllocated));
      }
    }

    current_fair_share = min_fair_share;
    for (map<int, set<Commodity*>*>::iterator it = AS->m_mpEdgeCodeCommodity.begin();
         it != AS->m_mpEdgeCodeCommodity.end(); ++it)
    {
      delete it->second;
      it->second = new set<Commodity*>();
    }

    /*
     * if all commodities are satisfied or no more bottleneck is found, exit the algorithm
     */
    //cout << "number of satisfied commodities = " << numSatisfiedCommodity << endl;
    if (numSatisfiedCommodity == AS->m_vCommodity.size() || !findBottle)
    {
      break;
    }
  }
}


void naive_TE(vector<Graph*>& ASes, vector<InterGraph*>& InterAS)
{
  /*
   * for each AS, generate its own topology table
   */
  for (vector<Graph*>::iterator it = ASes.begin(); it != ASes.end(); ++it)
  {
    (*it)->ComputeTopoTable();
  }

  /*
   * advertise topology table until all ASes' topology tables are stable
   */
  int timeStamp;
  for (timeStamp = 0; timeStamp < 1000000; timeStamp++)
  {
    bool notStable = false;
    //----------------exchange the topology table between neighboring AS--------------
    for (vector<Graph*>::iterator it = ASes.begin(); it != ASes.end(); ++it)
    {
      for (vector<BaseVertex*>::iterator it_border = (*it)->m_vtBorderVertices.begin();
           it_border != (*it)->m_vtBorderVertices.end(); ++it_border) // border switch in the AS
      {

        vector<InterGraph*>::iterator Inter_AS;
        for ( Inter_AS = InterAS.begin(); Inter_AS != InterAS.end(); ++Inter_AS)
        {
          set<BaseVertex*> vertex_set;
          if (((*Inter_AS)->m_left->get_graphID()) == (*it)->get_graphID() ||
              ((*Inter_AS)->m_right->get_graphID()) == (*it)->get_graphID())
          {
            int neighborASid; // this is the neighboring AS
            if (((*Inter_AS)->m_left->get_graphID()) == (*it)->get_graphID())
            {
              neighborASid = (*Inter_AS)->m_right->get_graphID();
            }
            else
            {
              neighborASid = (*Inter_AS)->m_left->get_graphID();
            }

            (*Inter_AS)->get_out_vertices(*it_border,vertex_set);
            if (vertex_set.size() == 0)
            {
              continue;
            }
            for (vector<Graph*>::iterator it_AS = ASes.begin(); it_AS != ASes.end(); ++it_AS)
            {
              //find the AS the peer switch belonging to
              if ((*it_AS)->get_graphID() == neighborASid)
              {
                for (set<BaseVertex*>::iterator it_peer = vertex_set.begin(); it_peer != vertex_set.end(); ++it_peer) // the peer switch of the border switch
                {
                  //debug
                  cout << "(" << (*it_border)->getGraphID() << "," << (*it_border)->getID() << ")-->"
                       << "(" << (*it_peer)->getGraphID() << "," << (*it_peer)->getID() << ")" << endl;

                  for (vector<TopoTableEntry>::iterator it_topoEntry = (*it_AS)->m_TopoTable.m_vEntry.begin();
                       it_topoEntry != (*it_AS)->m_TopoTable.m_vEntry.end(); ++it_topoEntry)
                  {
                    if (it_topoEntry->m_isExchanged.find((*it)->get_graphID()) != it_topoEntry->m_isExchanged.end())
                      continue;
                    if (it_topoEntry->m_source == (*it_peer))
                    {
                      bool isCycle = false;
                      for (vector<int>::iterator it_as = it_topoEntry->m_vASPath.begin();
                           it_as != it_topoEntry->m_vASPath.end(); ++it_as)
                      {
                        if (*it_as == (*it)->get_graphID())
                        {
                          isCycle = true;
                          break;
                        }
                      }
                      if (isCycle)
                      {
                        continue;
                      }
                      bool localStable;
                      it_topoEntry->m_isExchanged.insert((*it)->get_graphID());
		      // insert the new entry into the topo table


                      localStable = (*it)->UpdateTopoTable(TopoTableEntry(*it_border, it_topoEntry->m_sink, it_topoEntry->m_source,
                                                           it_topoEntry->m_weight + (*Inter_AS)->get_edge_weight(*it_border,*it_peer),
                                                           min(it_topoEntry->m_BW,(*Inter_AS)->get_edge_BW(*it_border,*it_peer)),it_topoEntry->m_vASPath));
                      notStable = localStable || notStable;
                    }
                  }
                }
                break;
              }
            }
          }
        }
      }
    }
    if (!notStable)
    {
      break;
    }
  }
  cout << "time stamp = " << timeStamp << endl;
}

void GenerateCommodity(vector<Graph*> ASes)
{
  int GenerateMethod = 1; // Generate specified commodity if 0, otherwise generate random commodities
  cout << "commodities: source node --> sink node: demand" << endl;
  if (GenerateMethod == 0)
  {
    int s_AS_id, s_node_id, d_AS_id, d_node_id;
    double demand;
    //----commodity 1-----
    s_AS_id = 0;
    s_node_id = 1;
    d_AS_id = 1;
    d_node_id = 0;
    demand = 4;
    Commodity* pt_commodity1 = new Commodity(ASes[s_AS_id]->get_vertex(s_node_id,ASes[s_AS_id]->get_graphID()),
        ASes[d_AS_id]->get_vertex(d_node_id,ASes[d_AS_id]->get_graphID()),demand);
    ASes[s_AS_id]->m_vCommodity.push_back(pt_commodity1);
    cout << "(" << s_AS_id << ", " << s_node_id << ") --> (" << d_AS_id << ", " << d_node_id << "): " << demand << endl;


    //----commodity 2-----
    s_AS_id = 0;
    s_node_id = 3;
    d_AS_id = 1;
    d_node_id = 1;
    demand = 3;
    Commodity* pt_commodity2 = new Commodity(ASes[s_AS_id]->get_vertex(s_node_id,ASes[s_AS_id]->get_graphID()),
        ASes[d_AS_id]->get_vertex(d_node_id,ASes[d_AS_id]->get_graphID()),demand);
    ASes[s_AS_id]->m_vCommodity.push_back(pt_commodity2);
    cout << "(" << s_AS_id << ", " << s_node_id << ") --> (" << d_AS_id << ", " << d_node_id << "): " << demand << endl;


    //----commodity 3-----
    s_AS_id = 2;
    s_node_id = 2;
    d_AS_id = 1;
    d_node_id = 3;
    demand = 6;
    Commodity* pt_commodity3 = new Commodity(ASes[s_AS_id]->get_vertex(s_node_id,ASes[s_AS_id]->get_graphID()),
        ASes[d_AS_id]->get_vertex(d_node_id,ASes[d_AS_id]->get_graphID()),demand);
    ASes[s_AS_id]->m_vCommodity.push_back(pt_commodity3);
    cout << "(" << s_AS_id << ", " << s_node_id << ") --> (" << d_AS_id << ", " << d_node_id << "): " << demand << endl;


    //----commodity 4-----
    s_AS_id = 1;
    s_node_id = 3;
    d_AS_id = 1;
    d_node_id = 0;
    demand = 10;
    Commodity* pt_commodity4 = new Commodity(ASes[s_AS_id]->get_vertex(s_node_id,ASes[s_AS_id]->get_graphID()),
        ASes[d_AS_id]->get_vertex(d_node_id,ASes[d_AS_id]->get_graphID()),demand);
    ASes[s_AS_id]->m_vCommodity.push_back(pt_commodity4);
    cout << "(" << s_AS_id << ", " << s_node_id << ") --> (" << d_AS_id << ", " << d_node_id << "): " << demand << endl;

  }
  else
  {
    srand(0);
    vector<Commodity*> v_commodity;
    for (int i=0; i < N_AS; i++)
    {
      for (int j = 0; j < ASes[i]->get_vertex_num(); ++j)
      {
        if (rand()%101/100.0 < prob_generate_traffic)
        {
          if (rand()%101/100.0 < prob_within_AS)
          {
            //the destination is in the same AS
            int desti = rand()%(ASes[i]->get_vertex_num());
            while (desti == j)
            {
              desti = rand()%(ASes[i]->get_vertex_num());
            }
            double demand = commodity_demand[SelectValue(commodity_demand,prob_demand)];
            demand *= loadC;
	    Commodity* pt_commodity = new Commodity(ASes[i]->get_vertex(j,ASes[i]->get_graphID()),ASes[i]->get_vertex(desti,ASes[i]->get_graphID()),demand);
            ASes[i]->m_vCommodity.push_back(pt_commodity);
            v_commodity.push_back(pt_commodity);
            cout << "(" << i << ", " << j << ") --> (" << i << ", " << desti << "): " << demand << endl;
          }
          else
          {
            // the destination is in another AS
            int desti_AS = rand()%N_AS;
            while (desti_AS == i)
            {
              desti_AS = rand()%N_AS;
            }
            int desti = rand()%(ASes[desti_AS]->get_vertex_num());
            double demand = commodity_demand[SelectValue(commodity_demand,prob_demand)];
            demand *= loadC;
	    Commodity* pt_commodity = new Commodity(ASes[i]->get_vertex(j,ASes[i]->get_graphID()),
                                                    ASes[desti_AS]->get_vertex(desti,ASes[desti_AS]->get_graphID()),demand);
            ASes[i]->m_vCommodity.push_back(pt_commodity);
            v_commodity.push_back(pt_commodity);
            cout << "(" << i << ", " << j << ") --> (" << desti_AS << ", " << desti << "): " << demand << endl;
          }
        }
      }
    }
  }
}
/*
 * map the vertex in the network to the vertex in one AS
 * if the graph is constructed based on several ASes,
 * this map is used to map the node in the AS to the node id of the constructed graph.
 * mpNodeID[graph_id*max_Nodes+vertex_id] = new vertex_id in the constructed graph
 */
BaseVertex* NetworkToAS(Graph& network, vector<Graph*> ASes, BaseVertex* ori_vertex)
{
  for (map<int, int>::iterator it = network.mpNodeID.begin(); it != network.mpNodeID.end(); ++it)
  {
    if (it->second == ori_vertex->getID())
    {
      int vertexID = it->first % max_Nodes;
      int graphID = (it->first - vertexID)/max_Nodes;
      return ASes[graphID]->get_vertex(vertexID,graphID);
    }
  }
  cout << "wrong: can not map the vertex in the network to the AS" << endl;
  return NULL;
}

int main(...)
{
  //--------create the topology-----------------
  vector<Graph*> ASes;
  string file_name = "data/test_6Degree10AS";

  for (int i = 0; i < N_AS; i++)
  {
    Graph* temp = new Graph(file_name);
    ASes.push_back(temp);
  }

  vector<InterGraph*> InterASes;

  for (vector<Graph*>::iterator itl = ASes.begin(); itl != ASes.end(); ++itl)
  {
    for (vector<Graph*>::iterator itr = itl + 1; itr != ASes.end(); ++itr)
    {
      InterGraph* temp = new InterGraph((*itl),(*itr),file_name);
      if (temp->m_nVertexNum != 0)
      {
        InterASes.push_back(temp);
      }
    }
  }
  cout << InterASes[0]->m_nVertexNum << endl;

  SDNi_TE(ASes,InterASes);

  /*
   * construct a new graph for each AS based on the topology table
   * 1. for every vertex in the topology table but not in the AS, we add a new edge from its peer to it
   * 2. the weight and bandwidth of the new edge equal to the corresponding item in the topology table
   */
  vector<Graph*> VirtualAS;
  for (vector<Graph*>::iterator it = ASes.begin(); it != ASes.end(); ++it)
  {
    Graph* temp = new Graph();
    (*it)->ConstructVirtualGraph(temp,file_name);
    VirtualAS.push_back(temp);
  }


  ofstream result;
  result.open("commodityTrace.txt");
  int time = 0;
  int runTime = 30;
  double Throughput = 0;
  double Aver_cost = 0;
  for (time = 0; time < runTime; time++)
  {
    //generate commodities
    GenerateCommodity(ASes);
    cout << "time = " << time << endl;
    // each AS do the traffic engineering for its own commodities based on the topology
    vector<map<Commodity*, set<pair<BasePath,double>,setcomp> >* > v_Allocation;
    double totalSendThr = 0;
    for (vector<Graph*>::iterator it = ASes.begin(); it != ASes.end(); ++it)
    {
      result << "AS " << (*it)->get_graphID() << "\n";
      int index = it - ASes.begin();
      /*
       * initialize the network
       */
      VirtualAS[index]->m_vCommodity.clear();
      for (map<int,double>::iterator it_used = VirtualAS[index]->m_mpEdgeCodeUsedBW.begin();
           it_used != VirtualAS[index]->m_mpEdgeCodeUsedBW.end(); ++it_used)
      {
        VirtualAS[index]->set_edge_UsedBW(it_used->first,0);
      }

      map<Commodity*, set<pair<BasePath,double>,setcomp> >* Allocation = new map<Commodity*, set<pair<BasePath,double>,setcomp> >();
      for (vector<Commodity*>::iterator it_commodity = (*it)->m_vCommodity.begin();
           it_commodity != (*it)->m_vCommodity.end(); ++it_commodity)
      {
        set<pair<BasePath,double>,setcomp> empty_set;

        int source_id = (*it_commodity)->source_->getGraphID()*max_Nodes
                        +(*it_commodity)->source_->getID();
        int sink_id = (*it_commodity)->sink_->getGraphID()*max_Nodes
                      +(*it_commodity)->sink_->getID();

        if (VirtualAS[index]->mpNodeID.find(source_id) == VirtualAS[index]->mpNodeID.end())
        {
          cout << "wrong in MainP.cpp: the source node is not find" << endl;
          exit(EXIT_FAILURE);
        }
        if (VirtualAS[index]->mpNodeID.find(sink_id) == VirtualAS[index]->mpNodeID.end())
        {
          cout << "Commodity (" << (*it_commodity)->source_->getGraphID() << "," << (*it_commodity)->source_->getID() << ") --> (" << (*it_commodity)->sink_->getGraphID() << "," << (*it_commodity)->sink_->getID() << ") is not reachable" << endl;
          continue;
        }
        int s_vertexID = VirtualAS[index]->mpNodeID.at(source_id);
        int d_vertexID = VirtualAS[index]->mpNodeID.at(sink_id);
        totalSendThr += (*it_commodity)->demand_;
        BaseVertex* p_source = VirtualAS[index]->get_vertex(s_vertexID,VirtualAS[index]->get_graphID());
        BaseVertex* p_sink = VirtualAS[index]->get_vertex(d_vertexID,VirtualAS[index]->get_graphID());
        Commodity* map_com = new Commodity(p_source,p_sink, (*it_commodity)->demand_);
        Allocation->insert(pair<Commodity*, set<pair<BasePath,double>,setcomp> >(map_com,empty_set));
        VirtualAS[index]->m_vCommodity.push_back(map_com);
      }
      cout << "The total send feasible throughput = " << totalSendThr << endl;
      //Google_TE_Optimization(VirtualAS[index], Allocation);
      //Max_Throughput_TE(VirtualAS[index], Allocation);

      if (rand()%101/100.0 < prob_Google)
      {
        Google_TE_Optimization(VirtualAS[index], Allocation);
      }
      else
      {
        Max_Throughput_TE(VirtualAS[index], Allocation);
      }




      /*
       * print out the allocation results
       */
      for (map<Commodity*, set<pair<BasePath,double>,setcomp> >::iterator it_com = Allocation->begin();
           it_com != Allocation->end(); ++it_com)
      {
        it_com->first->Print(result);
        //debug
        double totalAllocated = 0;
        for (set<pair<BasePath,double>,setcomp>::iterator it_set = it_com->second.begin();
             it_set != it_com->second.end(); ++ it_set)
        {
          it_set->first.PrintOut(result);
          totalAllocated += it_set->second;
        }
        if (it_com->first->Allocated_ == NOPATH)
        {
          it_com->first->Print(cout);
          cout << "No path for this commodity" << endl;
          cout << "***************************" << endl;

        }
        else if (abs(it_com->first->Allocated_ - totalAllocated) > delta)
        {
          it_com->first->Print(cout);
          cout << "wrong in main(): allocation does not match" << endl;
          cout << "***************************" << endl;
          exit (EXIT_FAILURE);
        }
      }
      v_Allocation.push_back(Allocation);
      VirtualAS[index]->recover_removed_edges();
    }
    result << "****************************************\n";

    vector<vector<Commodity*> > NewCommodity(ASes.size(),vector<Commodity*>());

    /*
     * update the commodities and packets sent
     */
    int index;
    for (vector<map<Commodity*, set<pair<BasePath,double>,setcomp> >* >::iterator it = v_Allocation.begin();
         it != v_Allocation.end(); ++it)
    {
      index = it - v_Allocation.begin();
      result << "AS " << index << "\n";
      cout << "*************************************************************" << endl;
      cout << "AS " << index << "\n";
      map<Commodity*, set<pair<BasePath,double>,setcomp> >::iterator it_map;
      /*
       * if the commodity's destination and source are in the same AS, then the commodity arrives the destination.
       * Otherwise, the commodity enter into transit AS
       */
      for (it_map = (*(*it)).begin(); it_map != (*(*it)).end(); ++it_map)
      {
        if (it_map->first->Allocated_ == NOPATH)
        {
          continue;
        }
        it_map->first->Print(result);

        int oriNodeID, oriASID;
        VirtualAS[index]->get_original_id(it_map->first->source_->getID(),oriNodeID,oriASID);
        cout << "(" << oriASID << "," << oriNodeID << ")-->";
        VirtualAS[index]->get_original_id(it_map->first->sink_->getID(),oriNodeID,oriASID);
        cout << "(" << oriASID << "," << oriNodeID << ");";
        cout << "demand = " << it_map->first->demand_ << "; Allocated =" << it_map->first->Allocated_ << endl;

        BaseVertex* MappedNode = NetworkToAS(*(VirtualAS[index]),ASes,it_map->first->sink_);
        if (MappedNode->getGraphID() == ASes[index]->get_graphID())
        {
          // the commodity arrives at the destination
          //Throughput += it_map->first->Allocated_;
          for (set<pair<BasePath,double>,setcomp>::iterator it_set = it_map->second.begin();
               it_set != it_map->second.end(); ++ it_set)
          {
            Throughput += it_set->second;
            Aver_cost += it_set->first.Weight();
            it_set->first.PrintOut(result);
            cout << "allocation to this path = " << it_set->second << endl;
            VirtualAS[index]->printPath(&(it_set->first));
          }
        }
        else
        {
          // the commodity enters into a transit AS
          for (set<pair<BasePath,double>,setcomp>::iterator it_set = it_map->second.begin();
               it_set != it_map->second.end(); ++it_set)
          {
            // find the first switch along the path which does not belong to the AS
            BaseVertex* p_border;
            cout << "allocation to this path = " << it_set->second << endl;
            VirtualAS[index]->get_original_id(it_set->first.GetVertex(0)->getID(),oriNodeID,oriASID);
            cout << "(" << oriASID << "," << oriNodeID << ")";
            int len = it_set->first.length();
            for (int i = 1; i < len; ++i)
            {
              BaseVertex* Pref_MappedNode = NetworkToAS(*(VirtualAS[index]),ASes,it_set->first.GetVertex(i-1));
              BaseVertex* MappedNode = NetworkToAS(*(VirtualAS[index]),ASes,it_set->first.GetVertex(i));

              VirtualAS[index]->get_original_id(it_set->first.GetVertex(i)->getID(),oriNodeID,oriASID);
              cout << "-->(" << oriASID << "," << oriNodeID << ")";

              Aver_cost += VirtualAS[index]->get_edge_weight(it_set->first.GetVertex(i-1),
                           it_set->first.GetVertex(i));
              if (MappedNode->getGraphID() != ASes[index]->get_graphID())
              {
                p_border = MappedNode;
                break;
              }
            }
            cout << endl;
            if (p_border == NetworkToAS(*(VirtualAS[index]),ASes,it_set->first.GetLastVertex()))
            {
              // the commodity arrives at the destination
              Throughput += it_set->second;
            }
            else
            {
              // generate the commodity in the neighboring AS
              NewCommodity[p_border->getGraphID()].push_back(new Commodity(p_border,
                  NetworkToAS(*(VirtualAS[index]),ASes,it_set->first.GetLastVertex()),
                  it_set->second));
            }
          }
        }
        cout << "Throughput=" << Throughput << "; Aver_cost=" << Aver_cost << endl;
        cout << "---------------------------------------" << endl;
      }
    }
    for (vector<Graph*>::iterator it = ASes.begin(); it != ASes.end(); ++it)
    {
      // clear the commodities in each AS
      (*it)->m_vCommodity.clear();
      int index = it - ASes.begin();
      for (vector<Commodity*>::iterator it_c = NewCommodity[index].begin();
           it_c != NewCommodity[index].end(); ++it_c)
      {
        (*it)->m_vCommodity.push_back(*it_c);
      }
    }
  }
  result.close();
}
