/************************************************************************/
/* difference between v6 and v5:
 * In this version, we build the advertised table and then exchange the table between
 * neighboring ASes, based on the exchanged information, the AS get the view of the
 * whole network, which can be used to do traffic engineering
 *
 * difference between v6 and v4:
 * In this version,we include AS path in the topoTable entry, which are used to
 * prevent from routing loops
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

int main(int argc, char** argv)
{
  // deal with the argument
  if (argc != 5){
    cout << "wrong number of arguments" << endl;
    return;
  }

  string file_name = argv[1];
  int N_AS = int(argv[2]);
  int kPath = int(argv[3]);
  double loadC = double(argv[4]);
  //--------create the topology-----------------
  vector<Graph*> ASes;

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
    GenerateCommodity(ASes, loadC, N_AS);
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
      //Google_TE_Optimization_Benchmark(VirtualAS[index], Allocation);
      //Max_Throughput_TE(VirtualAS[index], Allocation);


      if (rand()%101/100.0 < prob_Google)
      {
      	Google_TE_Optimization_Benchmark(VirtualAS[index], Allocation);
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
            // cost is defined as the sum of the costs of the links in one path
            // but not cost*flows
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
