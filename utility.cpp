#include "utility.h"

void process_mem_usage(double& vm_usage, double& resident_set)
{
    vm_usage     = 0.0;
    resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> vsize >> rss;
    }

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}

bool IsASCycle(vector<int>& ASpath, int graphID){
  for (auto it_as = ASpath.begin(); it_as != ASpath.end(); ++it_as)
  {
    if (*it_as == graphID)
    {
      return true;
    }
  }
  return false;
}

void GenerateCommodity(vector<Graph*> ASes, double loadC, int N_AS)
{
  int GenerateMethod = 1; // Generate specified commodity if 0, otherwise generate random commodities
  //cout << "commodities: source node --> sink node: demand" << endl;
  if (GenerateMethod == 0)
  {
    // this method is to generate explicit commodities such that debugging can be simplified
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
    //cout << "(" << s_AS_id << ", " << s_node_id << ") --> (" << d_AS_id << ", " << d_node_id << "): " << demand << endl;


    //----commodity 2-----
    s_AS_id = 0;
    s_node_id = 3;
    d_AS_id = 1;
    d_node_id = 1;
    demand = 3;
    Commodity* pt_commodity2 = new Commodity(ASes[s_AS_id]->get_vertex(s_node_id,ASes[s_AS_id]->get_graphID()),
        ASes[d_AS_id]->get_vertex(d_node_id,ASes[d_AS_id]->get_graphID()),demand);
    ASes[s_AS_id]->m_vCommodity.push_back(pt_commodity2);
    //cout << "(" << s_AS_id << ", " << s_node_id << ") --> (" << d_AS_id << ", " << d_node_id << "): " << demand << endl;


    //----commodity 3-----
    s_AS_id = 2;
    s_node_id = 2;
    d_AS_id = 1;
    d_node_id = 3;
    demand = 6;
    Commodity* pt_commodity3 = new Commodity(ASes[s_AS_id]->get_vertex(s_node_id,ASes[s_AS_id]->get_graphID()),
        ASes[d_AS_id]->get_vertex(d_node_id,ASes[d_AS_id]->get_graphID()),demand);
    ASes[s_AS_id]->m_vCommodity.push_back(pt_commodity3);
    //cout << "(" << s_AS_id << ", " << s_node_id << ") --> (" << d_AS_id << ", " << d_node_id << "): " << demand << endl;


    //----commodity 4-----
    s_AS_id = 1;
    s_node_id = 3;
    d_AS_id = 1;
    d_node_id = 0;
    demand = 10;
    Commodity* pt_commodity4 = new Commodity(ASes[s_AS_id]->get_vertex(s_node_id,ASes[s_AS_id]->get_graphID()),
        ASes[d_AS_id]->get_vertex(d_node_id,ASes[d_AS_id]->get_graphID()),demand);
    ASes[s_AS_id]->m_vCommodity.push_back(pt_commodity4);
    //cout << "(" << s_AS_id << ", " << s_node_id << ") --> (" << d_AS_id << ", " << d_node_id << "): " << demand << endl;

  }
  else
  {
    // generate random commodities
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
            //cout << "(" << i << ", " << j << ") --> (" << i << ", " << desti << "): " << demand << endl;
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
            //cout << "(" << i << ", " << j << ") --> (" << desti_AS << ", " << desti << "): " << demand << endl;
          }
        }
      }
    }
  }
}

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
  //cout << "wrong: can not map the vertex in the network to the AS" << endl;
  return NULL;
}

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
}


/*
 * 1. find the shortest path for each commodity
 * 2. find the shortest path among all the paths find in step 1
 * 3. serve this path with the capacity of this path
 */
void Max_Throughput_TE(Graph* AS, map<Commodity*, set<pair<BasePath,double>, setcomp> >* Allocation)
{
  double vm, rss;
  double oldvm, oldrss, newvm, newrss;
  //process_mem_usage(vm, rss);
  //cout << "Before call func, VM: " << vm << "; RSS: " << rss << endl;
  DijkstraShortestPathAlg shortest_path_alg(AS);
  double current_fair_share = 0;
  int num_iteration = 0;
  BasePath shortest_path;
  Commodity* target_com;
  double min_cost = 1000000000;
  double Thr;

  //process_mem_usage(vm, rss);
  //cout << "Before iterating, VM: " << vm << "; RSS: " << rss << endl;

  while (true)
  {
    //process_mem_usage(vm, rss);
    //cout << "Before finding the shortest path, VM: " << vm << "; RSS: " << rss << endl;
    /*
     * find the shortest path (with the minimal weight) for each commodity
     */
    min_cost = 1000000000;
    //double pathvm = 0;
    //double pathrss = 0;
    for (vector<Commodity*>::iterator it = AS->m_vCommodity.begin(); it != AS->m_vCommodity.end(); ++it)
    {
      if ((*it)->Allocated_ == NOPATH ||
          (*it)->Allocated_ == (*it)->demand_ || (*it)->isSaturated_)
      {
        continue;
      }
      //process_mem_usage(oldvm, oldrss);
      BasePath* result = shortest_path_alg.get_shortest_path((*it)->source_, (*it)->sink_);
      //process_mem_usage(newvm, newrss);
      //pathvm += newvm - oldvm;
      //pathrss += newrss - oldrss;
      if (result->length() > 0)
      {
        // if the path only contains one vertex
        if (result->GetVertex(0) == result->GetLastVertex())
        {
          //cout << "wrong in Max_Throughput_TE(): the path only contains one node" << endl;
        }
        else
        {
          if (min_cost > result->Weight())
          {
            min_cost = result->Weight();
            shortest_path = *result;
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
      //process_mem_usage(oldvm, oldrss);
      //std::cout << "before deleting, # BasePath objects = " << result->counter << '\n';
      delete result;
      result = NULL;
      //std::cout << "after deleting, # BasePath objects = " << result->counter << '\n';
      //process_mem_usage(newvm, newrss);
      //cout << "release vm = " << newvm - oldvm << "; release rss = " << newrss - oldrss << endl;
    }
    //cout << "Path vm = " << pathvm << "; path rss = " << pathrss << endl;

    if (min_cost == 1000000000)
    {
      break;
    }

    //determine the throughput of this path
    Thr = min(AS->get_path_BW(&shortest_path),target_com->demand_ - target_com->Allocated_);
    if (Thr < 0.1)
    {
      target_com->isSaturated_ = true;
      continue;
    }
    //process_mem_usage(vm, rss);
    //cout << "Before transmiting the shortest path, VM: " << vm << "; RSS: " << rss << endl;

    //transmit the shortest path
    for (int i = 0; i < shortest_path.length() -1; ++i)
    {
      double newUsedBW = AS->get_edge_UsedBW(shortest_path.GetVertex(i), shortest_path.GetVertex(i+1)) + Thr;
      AS->set_edge_UsedBW(shortest_path.GetVertex(i), shortest_path.GetVertex(i+1), newUsedBW);
      if (newUsedBW + delta > AS->get_edge_BW(shortest_path.GetVertex(i), shortest_path.GetVertex(i+1)))
      {
        pair<int, int> removed_edge(AS->get_vertex_code(shortest_path.GetVertex(i)),
                                    AS->get_vertex_code(shortest_path.GetVertex(i+1)));
        AS->remove_edge(removed_edge);
      }
    }

    /*
     * update the allocated bandwidth for each commodity
     */
    target_com->Allocated_ += Thr;

    //process_mem_usage(vm, rss);
    //cout << "Before updating the allocation, VM: " << vm << "; RSS: " << rss << endl;

    /*
     * update the return value
     */
    // check whether this path is already in the Allocation
    bool isFind = false;
    for (set<pair<BasePath,double> >::iterator it_set = Allocation->at(target_com).begin();
         it_set != Allocation->at(target_com).end(); ++it_set)
    {
      if (it_set->first == shortest_path)
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
      Allocation->at(target_com).insert(pair<BasePath,double>(shortest_path,extraAllocated));
    }
    num_iteration++;
    //process_mem_usage(vm, rss);
    //cout << "After one iteration, VM: " << vm << "; RSS: " << rss << endl;
  }
  //cout << "iteration number = " << num_iteration << endl;
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
          //cout << "wrong in Google_TE_Optimization(): the path only contains one node" << endl;
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

    ////cout << "prefer paths size = " << v_prefer_path.size() << endl;
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

    //TODO optimize the condition to exit such that less iterations are needed
    if (min_fair_share == current_fair_share){
      for (auto it = v_prefer_path.begin(); it != v_prefer_path.end(); ++it){
        delete (*it);
      }
      v_prefer_path.clear();
      break;
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

    // release memory allocated to BasePath object
    for (auto it = v_prefer_path.begin(); it != v_prefer_path.end(); ++it){
      delete (*it);
    }
    v_prefer_path.clear();

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
    if (numSatisfiedCommodity == AS->m_vCommodity.size() || !findBottle)
    {
      break;
    }
  }
}

/*
 * Google's TE optimization algorithm, please refer to S. Jain etc [B4]
 * 1. find the shortest path for each commodity
 * 2. find the bottleneck edge (with minimum fair share at its capacity)
 * 3. remove the bottleneck edge, repeat 1~2 until no more bottleneck edge is found or each commodity is satisfied
 */
void Google_TE_Optimization_Benchmark(Graph* AS, map<Commodity*, set<pair<BasePath,double>, setcomp> >* Allocation)
{
  DijkstraShortestPathAlg shortest_path_alg(AS);
  double current_fair_share = 0;
  //int numIterations = 0;
  while (true)
  {
    //numIterations ++;
    //cout << "start iteration = " << numIterations << endl;
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

      BasePath* result = shortest_path_alg.get_shortest_path((*it)->source_, (*it)->sink_);
      if (result->length() > 0)
      {
        // if the path only contains one vertex
        if (result->GetVertex(0) == result->GetLastVertex())
        {
          //cout << "wrong in Google_TE_Optimization(): the path only contains one node" << endl;
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

    //TODO optimize the condition to exit such that less iterations are needed
    if (min_fair_share == current_fair_share){
      for (auto it = v_prefer_path.begin(); it != v_prefer_path.end(); ++it){
        delete (*it);
      }
      v_prefer_path.clear();
      break;
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
      double preAllocate = (*it)->Allocated_;
      (*it)->Allocated_ += BWFunction((*it),min_fair_share) - BWFunction((*it),current_fair_share);
      // debug
      //(*it)->Print(cout);
      //cout << "allocation is increased by " << (*it)->Allocated_ - preAllocate << endl;
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

    // release memory allocated to BasePath object
    for (auto it = v_prefer_path.begin(); it != v_prefer_path.end(); ++it){
      delete (*it);
    }
    v_prefer_path.clear();

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
    if (numSatisfiedCommodity == AS->m_vCommodity.size() || !findBottle)
    {
      break;
    }
  }
  //cout << "number of iterations = " << numIterations << endl;
}
