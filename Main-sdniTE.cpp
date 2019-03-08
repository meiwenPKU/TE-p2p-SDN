/************************************************************************/
/* This script is used to measure the performance of SDNi_TE
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
#include "utility.h"
using namespace std;

int Graph::Graph_ID = 0;

void SDNi_TE(vector<Graph*>& ASes, vector<InterGraph*>& InterAS)
{
  /*
   * for each AS, generate its own topology table
   */
  for (vector<Graph*>::iterator it = ASes.begin(); it != ASes.end(); ++it)
  {
    (*it)->ComputeTopoTable();
  }

  /*
   * for each AS, compute its advertised topology table
   */
  for (vector<Graph*>::iterator it = ASes.begin(); it != ASes.end(); ++it)
  {
    (*it)->ComputeAdvertisedTable();
  }

  /*
   * exchange advertised topology table until all ASes' topology tables are stable
   */
  int timeStamp;
  int totalMsgExchange = 0;
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
          if (((*Inter_AS)->m_left->get_graphID()) == (*it)->get_graphID())
          {
            int neighborASid = (*Inter_AS)->m_right->get_graphID();; // this is the neighboring AS
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
                  // cout << "(" << (*it_border)->getGraphID() << "," << (*it_border)->getID() << ")-->"
                  //      << "(" << (*it_peer)->getGraphID() << "," << (*it_peer)->getID() << ")" << endl;
                  for (vector<TopoTableEntry>::iterator it_topoEntry = (*it_AS)->m_AdvertisedTopoTable.m_vEntry.begin();
                       it_topoEntry != (*it_AS)->m_AdvertisedTopoTable.m_vEntry.end(); ++it_topoEntry)
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
                      totalMsgExchange += 1;
                      it_topoEntry->m_isExchanged.insert((*it)->get_graphID());
                      auto new_entry = TopoTableEntry(
                          *it_border, it_topoEntry->m_sink,
                          it_topoEntry->m_source,
                          it_topoEntry->m_weight + (*Inter_AS)->get_edge_weight(
                                                       *it_border, *it_peer),
                          min(it_topoEntry->m_BW,
                              (*Inter_AS)->get_edge_BW(*it_border, *it_peer)),
                          it_topoEntry->m_vASPath);
                      // std::stringstream buffer;
                      // new_entry.printEntry(buffer);
                      // cout << buffer.str();
                      localStable = (*it)->UpdateTopoTable(new_entry);
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
    cout << "time stamp = " << timeStamp << "; total msg overhead = " << totalMsgExchange << endl;
    if (!notStable)
    {
      break;
    }
  }

}

int main(int argc, char** argv)
{
  // deal with the argument
  if (argc != 5){
    cout << "wrong number of arguments" << endl;
    return 0;
  }

  string file_name = argv[1];
  int N_AS = atoi(argv[2]);
  int kPath = atoi(argv[3]);
  double loadC = atof(argv[4]);

  //--------create the topology-----------------
  vector<Graph*> ASes;

  for (int i = 0; i < N_AS; i++)
  {
    Graph* temp = new Graph(file_name, kPath);
    ASes.push_back(temp);
  }

  vector<InterGraph*> InterASes;

  for (vector<Graph*>::iterator itl = ASes.begin(); itl != ASes.end(); ++itl)
  {
    for (vector<Graph*>::iterator itr = ASes.begin(); itr != ASes.end(); ++itr)
    {
      InterGraph* temp = new InterGraph((*itl),(*itr),file_name);
      if (temp->m_nVertexNum != 0)
      {
        InterASes.push_back(temp);
      }
    }
  }
  SDNi_TE(ASes,InterASes);

  /*
   * construct a new graph for each AS based on the topology table
   * 1. for every vertex in the topology table but not in the AS, we add a new edge from its peer to it
   * 2. the weight and bandwidth of the new edge equal to the corresponding item in the topology table
   */
  vector<Graph*> VirtualAS;
  int numTopoEntries = 0;
  int numAdvEntries = 0;
  for (vector<Graph*>::iterator it = ASes.begin(); it != ASes.end(); ++it)
  {
    Graph* temp = new Graph(kPath);
    //cout << "The topo table of AS " << (*it)->get_graphID() << endl;
    //(*it)->printTopoTable();
    numTopoEntries += (*it)->m_TopoTable.m_nEntry;
    numAdvEntries += (*it)->m_AdvertisedTopoTable.m_nEntry;
    //(*it)->printAdvertisedTable();
    (*it)->ConstructVirtualGraph(temp,file_name);
    //cout << "The built AS " << (*it)->get_graphID() << endl;
    //temp->printGraph();
    VirtualAS.push_back(temp);
  }
  cout << "number of topotable entries = " << numTopoEntries << ", number of AdvertisedTable = " << numAdvEntries << endl;


  //ofstream result;
  //result.open("commodityTrace.txt");
  int time = 0;
  int runTime = 30;
  double Throughput = 0;
  double Aver_cost = 0;
  double vm, rss;
  vector<map<Commodity*, set<pair<BasePath,double>,setcomp> > > v_Allocation;
  for (time = 0; time < runTime; time++)
  {
    //cout << "size of v_Allocation = " << sizeof(map<Commodity*, set<pair<BasePath,double>,setcomp> >) * v_Allocation.size() << endl;
    v_Allocation.clear();
    // get the memory footprint
    //process_mem_usage(vm, rss);
    //cout << "Before one iteratin, VM: " << vm << "; RSS: " << rss << endl;
    //generate commodities
    GenerateCommodity(ASes, loadC, N_AS);
    // get the memory footprint
    //process_mem_usage(vm, rss);
    //cout << "After generating commodities, VM: " << vm << "; RSS: " << rss << endl;

    cout << "time = " << time << endl;
    // each AS do the traffic engineering for its own commodities based on the topology
    double totalSendThr = 0;
    for (vector<Graph*>::iterator it = ASes.begin(); it != ASes.end(); ++it)
    {
      //result << "AS " << (*it)->get_graphID() << "\n";
      int index = it - ASes.begin();
      /*
       * initialize the network
       */
      VirtualAS[index]->clearCommodities();
      for (map<int,double>::iterator it_used = VirtualAS[index]->m_mpEdgeCodeUsedBW.begin();
           it_used != VirtualAS[index]->m_mpEdgeCodeUsedBW.end(); ++it_used)
      {
        VirtualAS[index]->set_edge_UsedBW(it_used->first,0);
      }

      map<Commodity*, set<pair<BasePath,double>,setcomp> > Allocation;
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
          //cout << "Commodity (" << (*it_commodity)->source_->getGraphID() << "," << (*it_commodity)->source_->getID() << ") --> (" << (*it_commodity)->sink_->getGraphID() << "," << (*it_commodity)->sink_->getID() << ") is not reachable" << endl;
          continue;
        }
        int s_vertexID = VirtualAS[index]->mpNodeID.at(source_id);
        int d_vertexID = VirtualAS[index]->mpNodeID.at(sink_id);
        totalSendThr += (*it_commodity)->demand_;
        BaseVertex* p_source = VirtualAS[index]->get_vertex(s_vertexID,VirtualAS[index]->get_graphID());
        BaseVertex* p_sink = VirtualAS[index]->get_vertex(d_vertexID,VirtualAS[index]->get_graphID());
        Commodity* map_com = new Commodity(p_source,p_sink, (*it_commodity)->demand_);
        Allocation[map_com] = empty_set;
        VirtualAS[index]->m_vCommodity.push_back(map_com);
      }
      //cout << "The total send feasible throughput = " << totalSendThr << endl;
      //Google_TE_Optimization(VirtualAS[index], &Allocation);
      // get the memory footprint
      //process_mem_usage(vm, rss);
      //cout << "Before TE, VM: " << vm << "; RSS: " << rss << endl;
      //Max_Throughput_TE(VirtualAS[index], &Allocation);
      // get the memory footprint
      //process_mem_usage(vm, rss);
      //cout << "After TE, VM: " << vm << "; RSS: " << rss << endl;

      if (rand()%101/100.0 < prob_Google)
      {
        Google_TE_Optimization(VirtualAS[index], &Allocation);
      }
      else
      {
        Max_Throughput_TE(VirtualAS[index], &Allocation);
      }


      /*
       * print out the allocation results
       */
      for (auto it_com = Allocation.begin();it_com != Allocation.end(); ++it_com)
      {
        //it_com->first->Print(result);
        double totalAllocated = 0;
        for (set<pair<BasePath,double>,setcomp>::iterator it_set = it_com->second.begin();
             it_set != it_com->second.end(); ++ it_set)
        {
          //it_set->first.PrintOut(result);
          totalAllocated += it_set->second;
        }
        if (it_com->first->Allocated_ == NOPATH)
        {
          //it_com->first->Print(cout);
          //cout << "No path for this commodity" << endl;
          //cout << "***************************" << endl;
        }
        else if (abs(it_com->first->Allocated_ - totalAllocated) > delta)
        {
          //it_com->first->Print(cout);
          //cout << "wrong in main(): allocation does not match" << endl;
          //cout << "***************************" << endl;
          exit (EXIT_FAILURE);
        }
      }
      v_Allocation.push_back(Allocation);
      VirtualAS[index]->recover_removed_edges();
    }
    //result << "****************************************\n";

    vector<vector<Commodity*> > NewCommodity(ASes.size(),vector<Commodity*>());

    /*
     * update the commodities and packets sent
     */
    int index;
    for (auto it = v_Allocation.begin(); it != v_Allocation.end(); ++it)
    {
      index = it - v_Allocation.begin();
      //result << "AS " << index << "\n";
      //cout << "*************************************************************" << endl;
      //cout << "AS " << index << "\n";
      map<Commodity*, set<pair<BasePath,double>,setcomp> >::iterator it_map;
      /*
       * if the commodity's destination and source are in the same AS, then the commodity arrives the destination.
       * Otherwise, the commodity enter into transit AS
       */
      for (it_map = (*it).begin(); it_map != (*it).end(); ++it_map)
      {
        if (it_map->first->Allocated_ == NOPATH)
        {
          continue;
        }
        //it_map->first->Print(result);

        int oriNodeID, oriASID;
        VirtualAS[index]->get_original_id(it_map->first->source_->getID(),oriNodeID,oriASID);
        //cout << "(" << oriASID << "," << oriNodeID << ")-->";
        VirtualAS[index]->get_original_id(it_map->first->sink_->getID(),oriNodeID,oriASID);
        //cout << "(" << oriASID << "," << oriNodeID << ");";
        //cout << "demand = " << it_map->first->demand_ << "; Allocated =" << it_map->first->Allocated_ << endl;

        BaseVertex* MappedNode = NetworkToAS(*(VirtualAS[index]),ASes,it_map->first->sink_);
        if (MappedNode->getGraphID() == ASes[index]->get_graphID())
        {
          // the commodity arrives at the destination
          //Throughput += it_map->first->Allocated_;
          for (set<pair<BasePath,double>,setcomp>::iterator it_set = it_map->second.begin();
               it_set != it_map->second.end(); ++ it_set)
          {
            Throughput += it_set->second;
            Aver_cost += it_set->second * it_set->first.Weight();
            //it_set->first.PrintOut(result);
            //cout << "allocation to this path = " << it_set->second << endl;
            //VirtualAS[index]->printPath(&(it_set->first));
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
            //cout << "allocation to this path = " << it_set->second << endl;
            VirtualAS[index]->get_original_id(it_set->first.GetVertex(0)->getID(),oriNodeID,oriASID);
            //cout << "(" << oriASID << "," << oriNodeID << ")";
            int len = it_set->first.length();
            for (int i = 1; i < len; ++i)
            {
              BaseVertex* Pref_MappedNode = NetworkToAS(*(VirtualAS[index]),ASes,it_set->first.GetVertex(i-1));
              BaseVertex* MappedNode = NetworkToAS(*(VirtualAS[index]),ASes,it_set->first.GetVertex(i));

              VirtualAS[index]->get_original_id(it_set->first.GetVertex(i)->getID(),oriNodeID,oriASID);
              //cout << "-->(" << oriASID << "," << oriNodeID << ")";

              Aver_cost += VirtualAS[index]->get_edge_weight(it_set->first.GetVertex(i-1),
                           it_set->first.GetVertex(i)) * it_set->second;
              if (MappedNode->getGraphID() != ASes[index]->get_graphID())
              {
                p_border = MappedNode;
                break;
              }
            }
            //cout << endl;
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
        //cout << "Throughput=" << Throughput << "; Aver_cost=" << Aver_cost << endl;
        //cout << "---------------------------------------" << endl;
      }
    }
    cout << "Throughput=" << Throughput << "; Aver_cost=" << Aver_cost << endl;
    // get the memory footprint
    //process_mem_usage(vm, rss);
    //cout << "After updating commodities, VM: " << vm << "; RSS: " << rss << endl;

    for (vector<Graph*>::iterator it = ASes.begin(); it != ASes.end(); ++it)
    {
      // clear the commodities in each AS
      (*it)->clearCommodities();
      int index = it - ASes.begin();
      for (vector<Commodity*>::iterator it_c = NewCommodity[index].begin();
           it_c != NewCommodity[index].end(); ++it_c)
      {
        (*it)->m_vCommodity.push_back(*it_c);
      }
    }
    // get the memory footprint
    //process_mem_usage(vm, rss);
    //cout << "After memory release, VM: " << vm << "; RSS: " << rss << endl;
  }
  // clear the memory for graphs
  for (auto g = ASes.begin(); g != ASes.end(); ++g){
    (*g)->clear();
  }
  for (auto g = InterASes.begin(); g != InterASes.end(); ++g){
    (*g)->clear();
  }
  for (auto g = VirtualAS.begin(); g != VirtualAS.end(); ++g){
    (*g)->clear();
  }
  //result.close();
}
