/************************************************************************/
/* This version is from Main-Benchmark.cpp, but use dijistra algo to determine
 * the allocation
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

void Dijkstra(Graph* AS, map<Commodity*, set<pair<BasePath,double>, setcomp> >* Allocation)
{
  DijkstraShortestPathAlg shortest_path_alg(AS);
  BasePath* shortest_path;
  Commodity* target_com;
  double Thr;


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
    BasePath* shortest_path = shortest_path_alg.get_shortest_path((*it)->source_, (*it)->sink_);
    if (shortest_path->length() > 0)
    {
      // if the path only contains one vertex
      if (shortest_path->GetVertex(0) == shortest_path->GetLastVertex())
      {
        cout << "wrong in Dijkstra(): the path only contains one node" << endl;
      }
      else
      {
        target_com = *it;
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
      continue;
    }

    //determine the throughput of this path
    Thr = min(AS->get_path_BW(shortest_path), target_com->demand_ - target_com->Allocated_);
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
  }
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


  /*
   * construct a global view of the entire network for each AS
   */
  vector<Graph*> VirtualAS;
  for (vector<Graph*>::iterator it = ASes.begin(); it != ASes.end(); ++it)
  {
    Graph* temp = new Graph();
    (*it)->ConstructGlobalView(temp,file_name);
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
    for (vector<Graph*>::iterator it = ASes.begin(); it != ASes.end(); ++it)
    {
      result << "AS " << (*it)->get_graphID() << "\n";
      int index = it - ASes.begin();
      /*
       * initialize the network
       */
      for(vector<Commodity*>::iterator it_com = VirtualAS[index]->m_vCommodity.begin();
          it_com != VirtualAS[index]->m_vCommodity.end(); ++it_com)
      {
        delete (*it_com);
      }
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
        int source_AS = (*it_commodity)->source_->getGraphID();
        int source_id = source_AS*max_Nodes + (*it_commodity)->source_->getID();
        int sink_AS = (*it_commodity)->sink_->getGraphID();
        int sink_id = sink_AS * max_Nodes + (*it_commodity)->sink_->getID();

        if (VirtualAS[index]->mpNodeID.find(source_id) == VirtualAS[index]->mpNodeID.end())
        {
          cout << "wrong in Main-Benchmark.cpp: the source node is not find" << endl;
          exit(EXIT_FAILURE);
        }
        if (VirtualAS[index]->mpNodeID.find(sink_id) == VirtualAS[index]->mpNodeID.end())
        {
          cout << "the destination is not reachable" << endl;
          continue;
        }
        int s_vertexID = VirtualAS[index]->mpNodeID.at(source_id);
        int d_vertexID = VirtualAS[index]->mpNodeID.at(sink_id);

        BaseVertex* p_source = VirtualAS[index]->get_vertex(s_vertexID,source_AS);
        BaseVertex* p_sink = VirtualAS[index]->get_vertex(d_vertexID,sink_AS);
        Commodity* map_com = new Commodity(p_source,p_sink, (*it_commodity)->demand_);
        Allocation->insert(pair<Commodity*, set<pair<BasePath,double>,setcomp> >(map_com,empty_set));
        VirtualAS[index]->m_vCommodity.push_back(map_com);
      }

      Dijkstra(VirtualAS[index], Allocation);
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
    // release the memory of v_allocation
    for (vector<map<Commodity*, set<pair<BasePath,double>,setcomp> >* >::iterator it = v_Allocation.begin();
         it != v_Allocation.end(); ++it)
    {
      delete (*it);
    }
    for (vector<Graph*>::iterator it = ASes.begin(); it != ASes.end(); ++it)
    {
      // clear the commodities in each AS
      for (vector<Commodity*>::iterator it_commodity = (*it)->m_vCommodity.begin();
           it_commodity != (*it)->m_vCommodity.end(); ++it_commodity)
      {
        delete (*it_commodity);
      }
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
