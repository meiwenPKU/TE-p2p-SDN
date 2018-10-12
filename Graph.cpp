///////////////////////////////////////////////////////////////////////////////
///  Graph.cpp
///  <TODO: insert file description here>
///
///  @remarks <TODO: insert remarks here>
///
///  @author Yan Qi @date 8/18/2010
///
///  $Id: Graph.cpp 65 2010-09-08 06:48:36Z yan.qi.asu $
///////////////////////////////////////////////////////////////////////////////

#include <limits>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "GraphElements.h"
#include "Graph.h"
#include "YenTopKShortestPathsAlg.h"
#include "FordFulkersonAlg.h"
#include "GlobalVariable.h"

const double Graph::DISCONNECT = (numeric_limits<double>::max)();
const double InterGraph::INTER_DISCONNECT = (numeric_limits<double>::max)();

void Graph::ConstructGlobalView(Graph* graph, const string& AS_file)
{
  const char* file_name = AS_file.c_str();
  ifstream ifs(file_name);
  if (!ifs)
  {
    cerr << "The file " << file_name << " can not be opened!" << endl;
    exit(1);
  }

  int num_AS, num_Nodes;
  ifs >> num_AS >> num_Nodes;

  graph->set_graphID(this->get_graphID());
  graph->set_vertex_num(num_AS*num_Nodes);
  /*
   * read the files for constructing the ASes
   */
  int start_vertex, start_inAS_vertex, end_vertex, end_inAS_vertex, start_AS, end_AS;
  double edge_weight, edge_BW;
  while (ifs >> start_vertex)
  {
    ifs >> end_vertex >> edge_weight >> edge_BW >> start_AS >> end_AS;
    edge_BW *= multiplier;
    start_inAS_vertex = start_vertex - num_Nodes * start_AS; // the vertex number in one AS
    end_inAS_vertex = end_vertex - num_Nodes * end_AS;

    int start_ori_id = start_inAS_vertex + start_AS * max_Nodes;
    if (graph->mpNodeID.find(start_ori_id) == graph->mpNodeID.end())
    {
      graph->mpNodeID.insert(pair<int,int>(start_ori_id,start_vertex));
    }

    int end_ori_id = end_inAS_vertex + end_AS * max_Nodes;
    if (graph->mpNodeID.find(end_ori_id) == graph->mpNodeID.end())
    {
      graph->mpNodeID.insert(pair<int,int>(end_ori_id,end_vertex));
    }

    ///3.2.1 construct the vertices
    BaseVertex* start_vertex_pt = graph->get_vertex(graph->mpNodeID.at(start_ori_id),start_AS);
    BaseVertex* end_vertex_pt = graph->get_vertex(graph->mpNodeID.at(end_ori_id),end_AS);

    int edgeCode = graph->get_edge_code(start_vertex_pt, end_vertex_pt);

    graph->m_mpEdgeCodeWeight[edgeCode] = edge_weight;
    graph->m_mpEdgeCodeBW[edgeCode] = edge_BW;
    graph->m_mpEdgeCodeUsedBW[edgeCode] = 0;
    set<Commodity*>* p_commodity = new set<Commodity*>();
    graph->m_mpEdgeCodeCommodity[edgeCode] = p_commodity;

    ///3.2.3 update the fan-in or fan-out variables
    //// Fan-in
    graph->get_vertex_set_pt(end_vertex_pt, graph->m_mpFaninVertices)->insert(start_vertex_pt);

    //// Fan-out
    graph->get_vertex_set_pt(start_vertex_pt, graph->m_mpFanoutVertices)->insert(end_vertex_pt);

  }
  ifs.close();
}

void Graph::ConstructVirtualGraphNaive(Graph* graph,const string& AS_file)
{
  graph->set_graphID(this->get_graphID());
  // count the number of the vertices
  int numVertices = this->get_vertex_num();
  set<BaseVertex*> s_outVertices;
  for (vector<TopoTableEntry>::iterator it = m_TopoTable.m_vEntry.begin();
       it != m_TopoTable.m_vEntry.end(); ++it)
  {
    if (it->m_sink->getGraphID() != this->get_graphID())
    {
      s_outVertices.insert(it->m_next);
      s_outVertices.insert(it->m_sink);
    }
  }
  numVertices += s_outVertices.size();
  graph->set_vertex_num(numVertices);

  const char* file_name = AS_file.c_str();
  ifstream ifs(file_name);
  if (!ifs)
  {
    cerr << "The file " << file_name << " can not be opened!" << endl;
    exit(1);
  }

  int num_AS, num_Nodes;
  ifs >> num_AS >> num_Nodes;

  if (num_Nodes != this->get_vertex_num())
  {
    cout << "wrong in Graph::ConstructVirtualGraph(): unmatched vertex number" << endl;
  }

  /*
   * read the files for constructing the ASes
   */
  int start_vertex, end_vertex, start_AS, end_AS;
  double edge_weight, edge_BW;
  int graphID = this->get_graphID();
  int nodeID = 0;
  while (ifs >> start_vertex)
  {
    ifs >> end_vertex >> edge_weight >> edge_BW >> start_AS >> end_AS;
    edge_BW *= multiplier;
    start_vertex = start_vertex - num_Nodes * start_AS;
    end_vertex = end_vertex - num_Nodes * end_AS;

    if (start_AS != graphID || end_AS != graphID)
      continue;

    int start_ori_id = start_vertex + graphID * max_Nodes;
    if (graph->mpNodeID.find(start_ori_id) == graph->mpNodeID.end())
    {
      graph->mpNodeID.insert(pair<int,int>(start_ori_id,nodeID));
      //cout << start_vertex << " " << graphID << " " << nodeID << endl;
      nodeID++;
    }

    int end_ori_id = end_vertex + graphID * max_Nodes;
    if (graph->mpNodeID.find(end_ori_id) == graph->mpNodeID.end())
    {
      graph->mpNodeID.insert(pair<int,int>(end_ori_id,nodeID));
      nodeID++;
    }

    ///3.2.1 construct the vertices
    BaseVertex* start_vertex_pt = graph->get_vertex(graph->mpNodeID.at(start_ori_id),graphID);
    BaseVertex* end_vertex_pt = graph->get_vertex(graph->mpNodeID.at(end_ori_id),graphID);

    int edgeCode = graph->get_edge_code(start_vertex_pt, end_vertex_pt);

    graph->m_mpEdgeCodeWeight[edgeCode] = edge_weight;
    graph->m_mpEdgeCodeBW[edgeCode] = edge_BW;
    graph->m_mpEdgeCodeUsedBW[edgeCode] = 0;
    set<Commodity*>* p_commodity = new set<Commodity*>();
    graph->m_mpEdgeCodeCommodity[edgeCode] = p_commodity;

    ///3.2.3 update the fan-in or fan-out variables
    //// Fan-in
    graph->get_vertex_set_pt(end_vertex_pt, graph->m_mpFaninVertices)->insert(start_vertex_pt);

    //// Fan-out
    graph->get_vertex_set_pt(start_vertex_pt, graph->m_mpFanoutVertices)->insert(end_vertex_pt);

  }
  ifs.close();
  // build the links to the nodes in other ASes
  for (vector<TopoTableEntry>::iterator it_topo = m_TopoTable.m_vEntry.begin();
       it_topo != m_TopoTable.m_vEntry.end(); ++it_topo)
  {
    if (it_topo->m_sink->getGraphID() != this->get_graphID())
    {
      int start_vertex, end_vertex, start_AS, end_AS;
      double edge_weight, edge_BW;
      if (it_topo->m_source != it_topo->m_next)
      {
        start_vertex = it_topo->m_source->getID();
        end_vertex = it_topo->m_next->getID();
        start_AS = it_topo->m_source->getGraphID();
        end_AS = it_topo->m_next->getGraphID();

        for (vector<TopoTableEntry>::iterator it_topo1 = m_TopoTable.m_vEntry.begin();
             it_topo1 != m_TopoTable.m_vEntry.end(); ++it_topo1)
        {
          if (it_topo1->m_source == it_topo->m_source &&
              it_topo1->m_next == it_topo1->m_sink &&
              it_topo1->m_sink == it_topo->m_next)
          {
            edge_weight = it_topo1->m_weight;
            edge_BW = it_topo1->m_BW;
            break;
          }
        }
        int end_ori_id = end_vertex + end_AS * max_Nodes;
        if (graph->mpNodeID.find(end_ori_id) == graph->mpNodeID.end())
        {
          graph->mpNodeID.insert(pair<int,int>(end_ori_id,nodeID));
          nodeID++;
        }

        ///3.1 construct the vertices
        BaseVertex* start_vertex_pt = graph->get_vertex(graph->mpNodeID.at(start_vertex+start_AS * max_Nodes),start_AS);
        BaseVertex* end_vertex_pt = graph->get_vertex(graph->mpNodeID.at(end_vertex+end_AS * max_Nodes),end_AS);


        ///3.2 add the edge weight
        //// note that the duplicate edge would overwrite the one occurring before.
        int edgeCode = graph->get_edge_code(start_vertex_pt, end_vertex_pt);
        // debug
        if (graph->m_mpEdgeCodeBW.find(edgeCode) != graph->m_mpEdgeCodeBW.end() && graph->m_mpEdgeCodeBW[edgeCode] != edge_BW){
          cout << "wrong in Graph::ConstructVirtualGraphNaive(): multiple links between two nodes" << endl;
        }
        graph->m_mpEdgeCodeWeight[edgeCode] = edge_weight;
        graph->m_mpEdgeCodeBW[edgeCode] = edge_BW;
        graph->m_mpEdgeCodeUsedBW[edgeCode] = 0;
        set<Commodity*>* p_commodity = new set<Commodity*>();
        graph->m_mpEdgeCodeCommodity[edgeCode] = p_commodity;
        ///3.3 update the fan-in or fan-out variables
        //// Fan-in
        graph->get_vertex_set_pt(end_vertex_pt, graph->m_mpFaninVertices)->insert(start_vertex_pt);
        //// Fan-out
        graph->get_vertex_set_pt(start_vertex_pt, graph->m_mpFanoutVertices)->insert(end_vertex_pt);
      }
      if (it_topo->m_next != it_topo->m_sink)
      {
        start_vertex = it_topo->m_next->getID();
        end_vertex = it_topo->m_sink->getID();
        start_AS = it_topo->m_next->getGraphID();
        end_AS = it_topo->m_sink->getGraphID();

        edge_weight = it_topo->m_weight - edge_weight;
        edge_BW = min (it_topo->m_BW,edge_BW);

        int end_ori_id = end_vertex + end_AS * max_Nodes;
        if (graph->mpNodeID.find(end_ori_id) == graph->mpNodeID.end())
        {
          graph->mpNodeID.insert(pair<int,int>(end_ori_id,nodeID));
          nodeID++;
        }

        ///3.1 construct the vertices
        BaseVertex* start_vertex_pt = graph->get_vertex(graph->mpNodeID.at(start_vertex+start_AS * max_Nodes),start_AS);
        BaseVertex* end_vertex_pt = graph->get_vertex(graph->mpNodeID.at(end_vertex+end_AS * max_Nodes),end_AS);


        ///3.2 add the edge weight
        //// note that when there are multiple edges between two nodes, say (w1, b1), (w2, b2), ..., (wn, bn)
        //// then the aggregated edge will be (max(w1, w2, ..., wn), b1+b2+...+bn)
        int edgeCode = graph->get_edge_code(start_vertex_pt, end_vertex_pt);
        if (graph->m_mpEdgeCodeBW.find(edgeCode) == graph->m_mpEdgeCodeBW.end()){
          graph->m_mpEdgeCodeWeight[edgeCode] = edge_weight;
          graph->m_mpEdgeCodeBW[edgeCode] = edge_BW;
          graph->m_mpEdgeCodeUsedBW[edgeCode] = 0;
          set<Commodity*>* p_commodity = new set<Commodity*>();
          graph->m_mpEdgeCodeCommodity[edgeCode] = p_commodity;
        } else {
          graph->m_mpEdgeCodeWeight[edgeCode] = graph->m_mpEdgeCodeWeight[edgeCode] > edge_weight ? graph->m_mpEdgeCodeWeight[edgeCode] : edge_weight;
          graph->m_mpEdgeCodeBW[edgeCode] += edge_BW;
        }

        ///3.3 update the fan-in or fan-out variables
        //// Fan-in
        graph->get_vertex_set_pt(end_vertex_pt, graph->m_mpFaninVertices)->insert(start_vertex_pt);
        //// Fan-out
        graph->get_vertex_set_pt(start_vertex_pt, graph->m_mpFanoutVertices)->insert(end_vertex_pt);
      }
    }
  }
  graph->m_nEdgeNum = graph->m_mpEdgeCodeWeight.size();
}

void Graph::ConstructVirtualGraph(Graph* graph,const string& AS_file)
{
  graph->set_graphID(this->get_graphID());
  // count the number of the vertices
  int numVertices = this->get_vertex_num();
  set<BaseVertex*> s_outVertices;
  for (vector<TopoTableEntry>::iterator it = m_TopoTable.m_vEntry.begin();
       it != m_TopoTable.m_vEntry.end(); ++it)
  {
    if (it->m_sink->getGraphID() != this->get_graphID())
    {
      s_outVertices.insert(it->m_next);
      s_outVertices.insert(it->m_sink);
    }
  }
  numVertices += s_outVertices.size();
  graph->set_vertex_num(numVertices);

  const char* file_name = AS_file.c_str();
  ifstream ifs(file_name);
  if (!ifs)
  {
    cerr << "The file " << file_name << " can not be opened!" << endl;
    exit(1);
  }

  int num_AS, num_Nodes;
  ifs >> num_AS >> num_Nodes;

  if (num_Nodes != this->get_vertex_num())
  {
    cout << "wrong in Graph::ConstructVirtualGraph(): unmatched vertex number" << endl;
  }

  /*
   * read the files for constructing the ASes
   */
  int start_vertex, end_vertex, start_AS, end_AS;
  double edge_weight, edge_BW;
  int graphID = this->get_graphID();
  int nodeID = 0;
  while (ifs >> start_vertex)
  {
    ifs >> end_vertex >> edge_weight >> edge_BW >> start_AS >> end_AS;
    edge_BW *= multiplier;
    start_vertex = start_vertex - num_Nodes * start_AS;
    end_vertex = end_vertex - num_Nodes * end_AS;

    if (start_AS != graphID || end_AS != graphID)
      continue;

    int start_ori_id = start_vertex + graphID * max_Nodes;
    if (graph->mpNodeID.find(start_ori_id) == graph->mpNodeID.end())
    {
      graph->mpNodeID.insert(pair<int,int>(start_ori_id,nodeID));
      //cout << start_vertex << " " << graphID << " " << nodeID << endl;
      nodeID++;
    }

    int end_ori_id = end_vertex + graphID * max_Nodes;
    if (graph->mpNodeID.find(end_ori_id) == graph->mpNodeID.end())
    {
      graph->mpNodeID.insert(pair<int,int>(end_ori_id,nodeID));
      //cout << end_vertex << " " << graphID << " " << nodeID << endl;
      nodeID++;
    }

    ///3.2.1 construct the vertices
    BaseVertex* start_vertex_pt = graph->get_vertex(graph->mpNodeID.at(start_ori_id),graphID);
    BaseVertex* end_vertex_pt = graph->get_vertex(graph->mpNodeID.at(end_ori_id),graphID);

    int edgeCode = graph->get_edge_code(start_vertex_pt, end_vertex_pt);

    graph->m_mpEdgeCodeWeight[edgeCode] = edge_weight;
    graph->m_mpEdgeCodeBW[edgeCode] = edge_BW;
    graph->m_mpEdgeCodeUsedBW[edgeCode] = 0;
    set<Commodity*>* p_commodity = new set<Commodity*>();
    graph->m_mpEdgeCodeCommodity[edgeCode] = p_commodity;

    ///3.2.3 update the fan-in or fan-out variables
    //// Fan-in
    graph->get_vertex_set_pt(end_vertex_pt, graph->m_mpFaninVertices)->insert(start_vertex_pt);

    //// Fan-out
    graph->get_vertex_set_pt(start_vertex_pt, graph->m_mpFanoutVertices)->insert(end_vertex_pt);

  }
  ifs.close();

  for (vector<TopoTableEntry>::iterator it_topo = m_TopoTable.m_vEntry.begin();
       it_topo != m_TopoTable.m_vEntry.end(); ++it_topo)
  {
    if (it_topo->m_sink->getGraphID() != this->get_graphID())
    {
      int start_vertex, end_vertex, start_AS, end_AS;
      double edge_weight, edge_BW;
      if (it_topo->m_source != it_topo->m_next)
      {
        start_vertex = it_topo->m_source->getID();
        end_vertex = it_topo->m_next->getID();
        start_AS = it_topo->m_source->getGraphID();
        end_AS = it_topo->m_next->getGraphID();

        for (vector<TopoTableEntry>::iterator it_topo1 = m_TopoTable.m_vEntry.begin();
             it_topo1 != m_TopoTable.m_vEntry.end(); ++it_topo1)
        {
          if (it_topo1->m_source == it_topo->m_source &&
              it_topo1->m_next == it_topo1->m_sink &&
              it_topo1->m_sink == it_topo->m_next)
          {
            edge_weight = it_topo1->m_weight;
            edge_BW = it_topo1->m_BW;
            break;
          }
        }
        int end_ori_id = end_vertex + end_AS * max_Nodes;
        if (graph->mpNodeID.find(end_ori_id) == graph->mpNodeID.end())
        {
          graph->mpNodeID.insert(pair<int,int>(end_ori_id,nodeID));
          //cout << end_vertex << " " << graphID << " " << nodeID << endl;
          nodeID++;
        }

        ///3.1 construct the vertices
        BaseVertex* start_vertex_pt = graph->get_vertex(graph->mpNodeID.at(start_vertex+start_AS * max_Nodes),start_AS);
        BaseVertex* end_vertex_pt = graph->get_vertex(graph->mpNodeID.at(end_vertex+end_AS * max_Nodes),end_AS);


        ///3.2 add the edge weight
        //// note that the duplicate edge would overwrite the one occurring before.
        int edgeCode = graph->get_edge_code(start_vertex_pt, end_vertex_pt);
        graph->m_mpEdgeCodeWeight[edgeCode] = edge_weight;
        graph->m_mpEdgeCodeBW[edgeCode] = edge_BW;
        graph->m_mpEdgeCodeUsedBW[edgeCode] = 0;
        set<Commodity*>* p_commodity = new set<Commodity*>();
        graph->m_mpEdgeCodeCommodity[edgeCode] = p_commodity;
        ///3.3 update the fan-in or fan-out variables
        //// Fan-in
        graph->get_vertex_set_pt(end_vertex_pt, graph->m_mpFaninVertices)->insert(start_vertex_pt);
        //// Fan-out
        graph->get_vertex_set_pt(start_vertex_pt, graph->m_mpFanoutVertices)->insert(end_vertex_pt);
      }
      if (it_topo->m_next != it_topo->m_sink)
      {
        start_vertex = it_topo->m_next->getID();
        end_vertex = it_topo->m_sink->getID();
        start_AS = it_topo->m_next->getGraphID();
        end_AS = it_topo->m_sink->getGraphID();

        edge_weight = it_topo->m_weight - edge_weight;
        edge_BW = min (it_topo->m_BW,edge_BW);

        int end_ori_id = end_vertex + end_AS * max_Nodes;
        if (graph->mpNodeID.find(end_ori_id) == graph->mpNodeID.end())
        {
          graph->mpNodeID.insert(pair<int,int>(end_ori_id,nodeID));
          //cout << end_vertex << " " << graphID << " " << nodeID << endl;
          nodeID++;
        }

        ///3.1 construct the vertices
        BaseVertex* start_vertex_pt = graph->get_vertex(graph->mpNodeID.at(start_vertex+start_AS * max_Nodes),start_AS);
        BaseVertex* end_vertex_pt = graph->get_vertex(graph->mpNodeID.at(end_vertex+end_AS * max_Nodes),end_AS);


        ///3.2 add the edge weight
        //// note that the duplicate edge would overwrite the one occurring before.
        int edgeCode = graph->get_edge_code(start_vertex_pt, end_vertex_pt);
        graph->m_mpEdgeCodeWeight[edgeCode] = edge_weight;
        graph->m_mpEdgeCodeBW[edgeCode] = edge_BW;
        graph->m_mpEdgeCodeUsedBW[edgeCode] = 0;
        set<Commodity*>* p_commodity = new set<Commodity*>();
        graph->m_mpEdgeCodeCommodity[edgeCode] = p_commodity;
        ///3.3 update the fan-in or fan-out variables
        //// Fan-in
        graph->get_vertex_set_pt(end_vertex_pt, graph->m_mpFaninVertices)->insert(start_vertex_pt);
        //// Fan-out
        graph->get_vertex_set_pt(start_vertex_pt, graph->m_mpFanoutVertices)->insert(end_vertex_pt);
      }
    }
  }
  graph->m_nEdgeNum = graph->m_mpEdgeCodeWeight.size();
}


InterGraph::InterGraph(Graph* left_graph, Graph* right_graph, const string& input_file_name)
{
  //-----------import data for edges connecting neighboring AS from the file------
  const char* file_name = input_file_name.c_str();
  m_left = left_graph;
  m_right = right_graph;

  //1. Check the validity of the file
  ifstream ifs(file_name);
  if (!ifs)
  {
    cerr << "The file " << file_name << " can not be opened!" << endl;
    exit(1);
  }

  //2. Reset the members of the class
  clear();


  int num_AS, num_Nodes;
  ifs >> num_AS >> num_Nodes;


  /*
   * read the files for constructing the ASes
   */
  int start_vertex, end_vertex, start_AS, end_AS;
  double edge_weight, edge_BW;

  //int nodeID = 0;
  while (ifs >> start_vertex)
  {
    ifs >> end_vertex >> edge_weight >> edge_BW >> start_AS >> end_AS;
    edge_BW *= multiplier;
    start_vertex = start_vertex - num_Nodes * start_AS;
    end_vertex = end_vertex - num_Nodes * end_AS;


    if (!((m_left->get_graphID() == start_AS && m_right->get_graphID() == end_AS)
          || (m_left->get_graphID() == end_AS && m_right->get_graphID() == start_AS)))
    {
      continue;
    }

    BaseVertex* start_vertex_pt = NULL;
    BaseVertex* end_vertex_pt = NULL;
    if (m_left->get_graphID() == start_AS)
    {
      start_vertex_pt = m_left->get_vertex(start_vertex,start_AS);
      end_vertex_pt = m_right->get_vertex(end_vertex,end_AS);

      if (!(m_left->FindBorderVertex(start_vertex_pt)))
      {
        m_left->m_vtBorderVertices.push_back(start_vertex_pt);
      }
      if (!(m_right->FindBorderVertex(end_vertex_pt)))
      {
        m_right->m_vtBorderVertices.push_back(end_vertex_pt);
      }
    }
    else
    {
      start_vertex_pt = m_right->get_vertex(start_vertex,start_AS);
      end_vertex_pt = m_left->get_vertex(end_vertex,end_AS);

      if (!(m_left->FindBorderVertex(end_vertex_pt)))
      {
        m_left->m_vtBorderVertices.push_back(end_vertex_pt);
      }
      if (!(m_right->FindBorderVertex(start_vertex_pt)))
      {
        m_right->m_vtBorderVertices.push_back(start_vertex_pt);
      }
    }

    m_mpVertexIndex[start_vertex_pt->getGraphID()*max_Nodes + start_vertex_pt->getID()] = start_vertex_pt;
    m_mpVertexIndex[end_vertex_pt->getGraphID()*max_Nodes + end_vertex_pt->getID()] = end_vertex_pt;

    if (!find_vertex(start_vertex_pt))
    {
      m_vtVertices.push_back(start_vertex_pt);
    }

    if (!find_vertex(end_vertex_pt))
    {
      m_vtVertices.push_back(end_vertex_pt);
    }
    ///3.2 add the edge weight
    //// note that the duplicate edge would overwrite the one occurring before.
    m_mpEdgeCodeWeight[get_edge_code(start_vertex_pt, end_vertex_pt)] = edge_weight;
    m_mpEdgeCodeBW[get_edge_code(start_vertex_pt, end_vertex_pt)] = edge_BW;
    m_mpEdgeCodeUsedBW[get_edge_code(start_vertex_pt, end_vertex_pt)] = 0;
    ///3.3 update the fan-in or fan-out variables
    //// Fan-in
    get_vertex_set_pt(end_vertex_pt, m_mpFaninVertices)->insert(start_vertex_pt);
    //// Fan-out
    get_vertex_set_pt(start_vertex_pt, m_mpFanoutVertices)->insert(end_vertex_pt);
  }

  m_nEdgeNum = m_mpEdgeCodeWeight.size();
  m_nVertexNum = m_vtVertices.size();
  ifs.close();
}

Graph::Graph( const string& file_name )
{
  _import_from_file(file_name);
}

Graph::Graph(const string& file_name, bool buildEntireNetwork)
{
  if (!buildEntireNetwork)
  {
    cerr << "buildEntireNetwork should be set as True" << endl;
    exit(1);
  }
  else
  {
    const char* input_file_name = file_name.c_str();

    //1. Check the validity of the file
    ifstream ifs(input_file_name);
    if (!ifs)
    {
      cerr << "The file " << input_file_name << " can not be opened!" << endl;
      exit(1);
    }
    clear();

    m_graphID = -1; // set as -1 which means this is the entire network not single AS network
    int num_AS, num_Nodes;
    ifs >> num_AS >> num_Nodes;
    m_nVertexNum = num_AS*num_Nodes;

    /*
     * read the file for constructing the entire network
     */
    int start_vertex, start_inAS_vertex, end_vertex, end_inAS_vertex, start_AS, end_AS;
    double edge_weight, edge_BW;

    int nodeID = 0;
    while (ifs >> start_vertex)
    {
      ifs >> end_vertex >> edge_weight >> edge_BW >> start_AS >> end_AS;
      edge_BW *= multiplier;
      start_inAS_vertex = start_vertex - num_Nodes * start_AS; // the vertex number in one AS
      end_inAS_vertex = end_vertex - num_Nodes * end_AS;
      ///3.2.1 construct the vertices
      BaseVertex* start_vertex_pt = get_vertex(start_vertex,start_AS);
      BaseVertex* end_vertex_pt = get_vertex(end_vertex,end_AS);

      int edgeCode = get_edge_code(start_vertex_pt, end_vertex_pt);

      m_mpEdgeCodeWeight[edgeCode] = edge_weight;
      m_mpEdgeCodeBW[edgeCode] = edge_BW;
      m_mpEdgeCodeUsedBW[edgeCode] = 0;
      set<Commodity*>* p_commodity = new set<Commodity*>();
      m_mpEdgeCodeCommodity[edgeCode] = p_commodity;

      ///3.2.3 update the fan-in or fan-out variables
      //// Fan-in
      get_vertex_set_pt(end_vertex_pt, m_mpFaninVertices)->insert(start_vertex_pt);

      //// Fan-out
      get_vertex_set_pt(start_vertex_pt, m_mpFanoutVertices)->insert(end_vertex_pt);
    }
    m_nEdgeNum = m_mpEdgeCodeWeight.size();
    if (m_vtVertices.size() != m_nVertexNum)
    {
      cout << "wrong happened: vertex number does not match. The right number should be "
           << m_nVertexNum << ", but gets "<< m_vtVertices.size() << endl;
    }
    ifs.close();
  }
}

Graph::Graph( const Graph* graph )
{
  //
  m_TopoTable = graph->m_TopoTable;
  m_AdvertisedTopoTable = graph->m_AdvertisedTopoTable;

  //vector
  m_vPathTable = graph->m_vPathTable;
  m_vCommodity = graph->m_vCommodity;
  m_vtVertices.assign(graph->m_vtVertices.begin(),graph->m_vtVertices.end());
  m_vtBorderVertices.assign(graph->m_vtBorderVertices.begin(),graph->m_vtBorderVertices.end());

  //variable
  m_nVertexNum = graph->m_nVertexNum;
  m_nEdgeNum = graph->m_nEdgeNum;
  m_maxBW = max_BW;
  m_graphID = graph->m_graphID;

  //map
  m_mpFaninVertices.insert(graph->m_mpFaninVertices.begin(),graph->m_mpFaninVertices.end());
  m_mpFanoutVertices.insert(graph->m_mpFanoutVertices.begin(),graph->m_mpFanoutVertices.end());
  m_mpEdgeCodeWeight.insert(graph->m_mpEdgeCodeWeight.begin(),graph->m_mpEdgeCodeWeight.end());
  m_mpVertexIndex.insert(graph->m_mpVertexIndex.begin(),graph->m_mpVertexIndex.end());
  m_mpEdgeCodeBW.insert(graph->m_mpEdgeCodeBW.begin(),graph->m_mpEdgeCodeBW.end());
  m_mpEdgeCodeUsedBW.insert(graph->m_mpEdgeCodeUsedBW.begin(),graph->m_mpEdgeCodeUsedBW.end());
  m_mpEdgeCodeCommodity.insert(graph->m_mpEdgeCodeCommodity.begin(),graph->m_mpEdgeCodeCommodity.end());
  mpNodeID.insert(graph->mpNodeID.begin(),graph->mpNodeID.end());

  //set
  m_stRemovedVertexIds = graph->m_stRemovedVertexIds;
  m_stRemovedEdge = graph->m_stRemovedEdge;
}

Graph::~Graph(void)
{
  clear();
}

///////////////////////////////////////////////////////////////////////////////
///  public  _import_from_file
///  Construct the graph by importing the edges from the input file.
///
///  @param [in]       file_name const std::string &    The input graph file
///
///  This function doesn't return a value
///
///  @remarks The format of the file is as follows:
///   1. The first line has an integer as the number of vertices of the graph
///   2. Each line afterwards contains a directed edge in the graph:
///		     starting point, ending point and the weight of the edge.
///		 These values are separated by 'white space'.
///
///  @see <TODO: insert text here>
///
///  @author Yan Qi @date 5/29/2010
///////////////////////////////////////////////////////////////////////////////
void Graph::_import_from_file( const string& input_file_name )
{
  const char* file_name = input_file_name.c_str();

  //1. Check the validity of the file
  ifstream ifs(file_name);
  if (!ifs)
  {
    cerr << "The file " << file_name << " can not be opened!" << endl;
    exit(1);
  }
  clear();

  m_graphID = Graph_ID;
  Graph_ID++;

  int num_AS, num_Nodes;
  ifs >> num_AS >> num_Nodes;
  m_nVertexNum = num_Nodes;

  /*
   * read the files for constructing the ASes
   */
  int start_vertex, end_vertex, start_AS, end_AS;
  double edge_weight, edge_BW;

  int nodeID = 0;
  while (ifs >> start_vertex)
  {
    ifs >> end_vertex >> edge_weight >> edge_BW >> start_AS >> end_AS;
    edge_BW *= multiplier;
    start_vertex = start_vertex - num_Nodes * start_AS;
    end_vertex = end_vertex - num_Nodes * end_AS;

    if (start_AS != m_graphID || end_AS != m_graphID)
      continue;

    ///3.2.1 construct the vertices
    BaseVertex* start_vertex_pt = get_vertex(start_vertex,m_graphID);
    BaseVertex* end_vertex_pt = get_vertex(end_vertex,m_graphID);

    ///3.2.2 add the edge weight
    //// note that the duplicate edge would overwrite the one occurring before.
    m_mpEdgeCodeWeight[get_edge_code(start_vertex_pt, end_vertex_pt)] = edge_weight;

    //hemin
    m_mpEdgeCodeBW[get_edge_code(start_vertex_pt, end_vertex_pt)] = edge_BW;
    m_mpEdgeCodeUsedBW[get_edge_code(start_vertex_pt, end_vertex_pt)] = 0;

    ///3.2.3 update the fan-in or fan-out variables
    //// Fan-in
    get_vertex_set_pt(end_vertex_pt, m_mpFaninVertices)->insert(start_vertex_pt);

    //// Fan-out
    get_vertex_set_pt(start_vertex_pt, m_mpFanoutVertices)->insert(end_vertex_pt);

  }

  m_nEdgeNum = m_mpEdgeCodeWeight.size();
  if (m_vtVertices.size() != m_nVertexNum)
  {
    cout << "wrong happened: vertex number does not match" << endl;
  }
  ifs.close();
}

void Graph::set_graphID(int graphID)
{
  m_graphID = graphID;
}

void Graph::set_vertex_num(int numVertex)
{
  m_nVertexNum = numVertex;
}

bool Graph::FindBorderVertex(BaseVertex* vertex)
{
  for (vector<BaseVertex*>::iterator it = m_vtBorderVertices.begin();
       it != m_vtBorderVertices.end(); ++it)
  {
    if (*it == vertex)
      return true;
  }
  return false;
}

bool InterGraph::find_vertex(BaseVertex* vertex)
{
  for (vector<BaseVertex*>::iterator it = m_vtVertices.begin(); it != m_vtVertices.end(); ++it)
  {
    if (*it == vertex)
    {
      return true;
    }
  }
  return false;
}

BaseVertex* Graph::get_vertex( int node_id, int graphID )
{
  if (m_stRemovedVertexIds.find(node_id) != m_stRemovedVertexIds.end())
  {
    return NULL;
  }
  else
  {
    BaseVertex* vertex_pt = NULL;
    const map<int, BaseVertex*>::iterator pos = m_mpVertexIndex.find(node_id);
    if (pos == m_mpVertexIndex.end())
    {
      int vertex_id = m_vtVertices.size();
      vertex_pt = new BaseVertex();
      vertex_pt->setID(node_id);
      vertex_pt->setGraphID(graphID);
      m_mpVertexIndex[node_id] = vertex_pt;
      m_vtVertices.push_back(vertex_pt);
    }
    else
    {
      vertex_pt = pos->second;
    }
    return vertex_pt;
  }
}

BaseVertex* InterGraph::get_vertex( int node_id, int graphID )
{
  BaseVertex* vertex_pt = NULL;
  const map<int, BaseVertex*>::iterator pos = m_mpVertexIndex.find(node_id);
  if (pos == m_mpVertexIndex.end())
  {
    int vertex_id = m_vtVertices.size();
    vertex_pt = new BaseVertex();
    vertex_pt->setID(node_id);
    vertex_pt->setGraphID(graphID);
    m_mpVertexIndex[node_id] = vertex_pt;

    m_vtVertices.push_back(vertex_pt);
  }
  else
  {
    vertex_pt = pos->second;
  }

  return vertex_pt;
}

void Graph::printGraph()
{
  cout << "Source    Sink    Weight    Bandwidth" << endl;
  for (auto it = m_mpFanoutVertices.begin(); it != m_mpFanoutVertices.end(); ++it){
    for (auto it_set = it->second->begin(); it_set != it->second->end(); ++it_set){
      auto edgeCode = get_edge_code(it->first, *it_set);
      cout << it->first->getID() << "    " << (*it_set)->getID() << "    " << m_mpEdgeCodeWeight.at(edgeCode) << "    " << m_mpEdgeCodeBW.at(edgeCode) << endl;
    }
  }
}

void Graph::clear()
{
  m_nEdgeNum = 0;
  m_nVertexNum = 0;
  m_graphID = 0;

  for(map<BaseVertex*, set<BaseVertex*>*>::const_iterator pos=m_mpFaninVertices.begin();
      pos!=m_mpFaninVertices.end(); ++pos)
  {
    delete pos->second;
  }
  m_mpFaninVertices.clear();

  for(map<BaseVertex*, set<BaseVertex*>*>::const_iterator pos=m_mpFanoutVertices.begin();
      pos!=m_mpFanoutVertices.end(); ++pos)
  {
    delete pos->second;
  }
  m_mpFanoutVertices.clear();

  for (map<int, set<Commodity*>*>::const_iterator pos = m_mpEdgeCodeCommodity.begin();
       pos != m_mpEdgeCodeCommodity.end(); ++pos)
  {
    delete pos->second;
  }
  m_mpEdgeCodeCommodity.clear();




  m_mpEdgeCodeWeight.clear();
  //--hemin---
  m_mpEdgeCodeBW.clear();
  m_mpEdgeCodeUsedBW.clear();
  mpNodeID.clear();
  //----------

  //clear the list of vertices objects
  for_each(m_vtVertices.begin(), m_vtVertices.end(), DeleteFunc<BaseVertex>());
  m_vtVertices.clear();

  for_each(m_vPathTable.begin(), m_vPathTable.end(), DeleteFunc<BasePath>());
  m_vPathTable.clear();

  for_each(m_vCommodity.begin(), m_vCommodity.end(), DeleteFunc<Commodity>());
  m_vCommodity.clear();

  m_vtBorderVertices.clear();
  m_mpVertexIndex.clear();

  m_stRemovedVertexIds.clear();
  m_stRemovedEdge.clear();
}

void InterGraph::clear()
{
  m_nEdgeNum = 0;
  m_nVertexNum = 0;

  for(map<BaseVertex*, set<BaseVertex*>*>::const_iterator pos=m_mpFaninVertices.begin();
      pos!=m_mpFaninVertices.end(); ++pos)
  {
    delete pos->second;
  }
  m_mpFaninVertices.clear();

  for(map<BaseVertex*, set<BaseVertex*>*>::const_iterator pos=m_mpFanoutVertices.begin();
      pos!=m_mpFanoutVertices.end(); ++pos)
  {
    delete pos->second;
  }
  m_mpFanoutVertices.clear();

  m_mpEdgeCodeWeight.clear();
  m_mpEdgeCodeBW.clear();
  m_mpEdgeCodeUsedBW.clear();

  //clear the list of vertices objects
  for_each(m_vtVertices.begin(), m_vtVertices.end(), DeleteFunc<BaseVertex>());
  m_vtVertices.clear();
  m_mpVertexIndex.clear();

}

int Graph::get_edge_code( const BaseVertex* start_vertex_pt, const BaseVertex* end_vertex_pt ) const
{
  /// Note that the computation below works only if
  /// the result is smaller than the maximum of an integer!
  return start_vertex_pt->getID()*m_nVertexNum+end_vertex_pt->getID();
}

int InterGraph::get_edge_code( const BaseVertex* start_vertex_pt, const BaseVertex* end_vertex_pt ) const
{
  /// Note that the computation below works only if
  /// the result is smaller than the maximum of an integer!
  int vertex_num;
  if (start_vertex_pt->getGraphID() == m_left->get_graphID())
  {
    vertex_num = m_right->get_vertex_num();
  }
  else if (start_vertex_pt->getGraphID() == m_right->get_graphID())
  {
    vertex_num = m_left->get_vertex_num();
  }
  else
  {
    cout << "wrong happened in InterGraph::get_edge_code, graph id does not match" << endl;
    return -1;
  }
  return start_vertex_pt->getID()*vertex_num+end_vertex_pt->getID();
}

set<BaseVertex*>* Graph::get_vertex_set_pt( BaseVertex* vertex_, map<BaseVertex*, set<BaseVertex*>*>& vertex_container_index )
{
  BaseVertexPt2SetMapIterator pos = vertex_container_index.find(vertex_);

  if(pos == vertex_container_index.end())
  {
    set<BaseVertex*>* vertex_set = new set<BaseVertex*>();
    pair<BaseVertexPt2SetMapIterator,bool> ins_pos =
      vertex_container_index.insert(make_pair(vertex_, vertex_set));

    pos = ins_pos.first;
  }

  return pos->second;
}

set<BaseVertex*>* InterGraph::get_vertex_set_pt( BaseVertex* vertex_, map<BaseVertex*, set<BaseVertex*>*>& vertex_container_index )
{
  BaseVertexPt2SetMapIterator pos = vertex_container_index.find(vertex_);

  if(pos == vertex_container_index.end())
  {
    set<BaseVertex*>* vertex_set = new set<BaseVertex*>();
    pair<BaseVertexPt2SetMapIterator,bool> ins_pos =
      vertex_container_index.insert(make_pair(vertex_, vertex_set));

    pos = ins_pos.first;
  }

  return pos->second;
}

double Graph::get_edge_weight( const BaseVertex* source, const BaseVertex* sink )
{
  int source_id = source->getID();
  int sink_id = sink->getID();

  if (m_stRemovedVertexIds.find(source_id) != m_stRemovedVertexIds.end()
      || m_stRemovedVertexIds.find(sink_id) != m_stRemovedVertexIds.end()
      || m_stRemovedEdge.find(make_pair(source_id, sink_id)) != m_stRemovedEdge.end())
  {
    return DISCONNECT;
  }
  else
  {
    return get_original_edge_weight(source, sink);
  }
}

double InterGraph::get_edge_weight( const BaseVertex* source, const BaseVertex* sink )
{
  map<int, double>::const_iterator pos =
    m_mpEdgeCodeWeight.find(get_edge_code(source, sink));

  if (pos != m_mpEdgeCodeWeight.end())
  {
    return pos->second;
  }
  else
  {
    return INTER_DISCONNECT;
  }
}

//---hemin-----------
double Graph::get_edge_BW( const BaseVertex* source, const BaseVertex* sink )
{
  int source_id = source->getID();
  int sink_id = sink->getID();

  if (m_stRemovedVertexIds.find(source_id) != m_stRemovedVertexIds.end()
      || m_stRemovedVertexIds.find(sink_id) != m_stRemovedVertexIds.end()
      || m_stRemovedEdge.find(make_pair(source_id, sink_id)) != m_stRemovedEdge.end())
  {
    return DISCONNECT;
  }
  else
  {
    return get_original_edge_BW(source, sink);
  }
}

double InterGraph::get_edge_BW( const BaseVertex* source, const BaseVertex* sink )
{
  map<int, double>::const_iterator pos =
    m_mpEdgeCodeBW.find(get_edge_code(source, sink));

  if (pos != m_mpEdgeCodeBW.end())
  {
    return pos->second;
  }
  else
  {
    return INTER_DISCONNECT;
  }
}

double Graph::get_edge_BW(int code)
{
  return m_mpEdgeCodeBW.at(code);
}

void Graph::set_edge_BW( const BaseVertex* source, const BaseVertex* sink, double BW )
{
  int source_id = source->getID();
  int sink_id = sink->getID();

  if (m_stRemovedVertexIds.find(source_id) != m_stRemovedVertexIds.end()
      || m_stRemovedVertexIds.find(sink_id) != m_stRemovedVertexIds.end()
      || m_stRemovedEdge.find(make_pair(source_id, sink_id)) != m_stRemovedEdge.end())
  {

  }
  else
  {
    set_original_edge_BW(source, sink, BW);
  }
}

void Graph::set_edge_BW(int code, double BW)
{
  m_mpEdgeCodeBW.at(code) = BW;
}

double Graph::get_edge_UsedBW( const BaseVertex* source, const BaseVertex* sink )
{
  int source_id = source->getID();
  int sink_id = sink->getID();

  if (m_stRemovedVertexIds.find(source_id) != m_stRemovedVertexIds.end()
      || m_stRemovedVertexIds.find(sink_id) != m_stRemovedVertexIds.end()
      || m_stRemovedEdge.find(make_pair(source_id, sink_id)) != m_stRemovedEdge.end())
  {
    return DISCONNECT;
  }
  else
  {
    return get_original_edge_UsedBW(source, sink);
  }
}

double Graph::get_edge_UsedBW(int code)
{
  return m_mpEdgeCodeUsedBW.at(code);
}

void Graph::set_edge_UsedBW( const BaseVertex* source, const BaseVertex* sink, double UsedBW )
{
  int source_id = source->getID();
  int sink_id = sink->getID();

  if (m_stRemovedVertexIds.find(source_id) != m_stRemovedVertexIds.end()
      || m_stRemovedVertexIds.find(sink_id) != m_stRemovedVertexIds.end()
      || m_stRemovedEdge.find(make_pair(source_id, sink_id)) != m_stRemovedEdge.end())
  {

  }
  else
  {
    set_original_edge_UsedBW(source, sink, UsedBW);
  }
}

void Graph::set_edge_UsedBW(int code, double UsedBW)
{
  m_mpEdgeCodeUsedBW.at(code) = UsedBW;
}


//-------------------

void Graph::get_adjacent_vertices( BaseVertex* vertex, set<BaseVertex*>& vertex_set )
{
  int starting_vt_id = vertex->getID();

  if (m_stRemovedVertexIds.find(starting_vt_id) == m_stRemovedVertexIds.end())
  {
    set<BaseVertex*>* vertex_pt_set = get_vertex_set_pt(vertex, m_mpFanoutVertices);
    for(set<BaseVertex*>::const_iterator pos=(*vertex_pt_set).begin();
        pos != (*vertex_pt_set).end(); ++pos)
    {
      int ending_vt_id = (*pos)->getID();
      if (m_stRemovedVertexIds.find(ending_vt_id) != m_stRemovedVertexIds.end()
          || m_stRemovedEdge.find(make_pair(starting_vt_id, ending_vt_id)) != m_stRemovedEdge.end())
      {
        continue;
      }
      //
      vertex_set.insert(*pos);
    }
  }
}

void InterGraph::get_out_vertices( BaseVertex* vertex, set<BaseVertex*>& vertex_set )
{
  int starting_vt_id = vertex->getID();

  BaseVertexPt2SetMapIterator pos = m_mpFanoutVertices.find(vertex);
  if (pos == m_mpFanoutVertices.end())
  {
    return;
  }
  else
  {
    for (set<BaseVertex*>::const_iterator pos_set = (*pos->second).begin();
         pos_set != (*pos->second).end(); ++pos_set)
    {
      vertex_set.insert(*pos_set);
    }
  }
}

void InterGraph::get_in_vertices( BaseVertex* vertex, set<BaseVertex*>& vertex_set )
{
  int starting_vt_id = vertex->getID();

  BaseVertexPt2SetMapIterator pos = m_mpFaninVertices.find(vertex);
  if (pos == m_mpFaninVertices.end())
  {
    return;
  }
  else
  {
    for (set<BaseVertex*>::const_iterator pos_set = (*pos->second).begin();
         pos_set != (*pos->second).end(); ++pos_set)
    {
      vertex_set.insert(*pos_set);
    }
  }
}

void Graph::get_precedent_vertices( BaseVertex* vertex, set<BaseVertex*>& vertex_set )
{
  if (m_stRemovedVertexIds.find(vertex->getID()) == m_stRemovedVertexIds.end())
  {
    int ending_vt_id = vertex->getID();
    set<BaseVertex*>* pre_vertex_set = get_vertex_set_pt(vertex, m_mpFaninVertices);
    for(set<BaseVertex*>::const_iterator pos=(*pre_vertex_set).begin();
        pos != (*pre_vertex_set).end(); ++pos)
    {
      int starting_vt_id = (*pos)->getID();
      if (m_stRemovedVertexIds.find(starting_vt_id) != m_stRemovedVertexIds.end()
          || m_stRemovedEdge.find(make_pair(starting_vt_id, ending_vt_id)) != m_stRemovedEdge.end())
      {
        continue;
      }
      //
      vertex_set.insert(*pos);
    }
  }
}

void Graph::get_original_id(int newNodeID, int& oriNodeID, int& oriGraphID)
{
  for (map<int,int>::iterator it_m = mpNodeID.begin();
       it_m != mpNodeID.end(); ++it_m)
  {
    if (it_m->second == newNodeID)
    {
      oriGraphID = it_m->first / max_Nodes;
      oriNodeID = it_m->first % max_Nodes;
      break;
    }
  }
}

double Graph::get_original_edge_weight( const BaseVertex* source, const BaseVertex* sink )
{
  map<int, double>::const_iterator pos =
    m_mpEdgeCodeWeight.find(get_edge_code(source, sink));

  if (pos != m_mpEdgeCodeWeight.end())
  {
    return pos->second;
  }
  else
  {
    return DISCONNECT;
  }
}


//-------------------hemin-----------------------------
void Graph::printTopoTable()
{
  m_TopoTable.printTable();
}

void Graph::printAdvertisedTable()
{
  m_AdvertisedTopoTable.printTable();
}

void Graph::printTableEntry(TopoTableEntry& entry)
{
  entry.printEntry();
}

void Graph::printPath(const BasePath* path)
{
  std::cout << "Cost: " << path->Weight() << " Length: " << path->length() << std::endl;
  for (int i = 0; i < path->length(); ++i)
  {
    BaseVertex* temp = path->GetVertex(i);
    int oriNodeID, oriASID;
    get_original_id(temp->getID(),oriNodeID, oriASID);
    std::cout << "(" << oriASID << "," << oriNodeID << ")";
    if (i != path->length() - 1)
    {
      std::cout << "->";
    }
  }
  std::cout << std::endl;
}
double Graph::get_original_edge_BW( const BaseVertex* source, const BaseVertex* sink )
{
  map<int, double>::const_iterator pos =
    m_mpEdgeCodeBW.find(get_edge_code(source, sink));

  if (pos != m_mpEdgeCodeBW.end())
  {
    return pos->second;
  }
  else
  {
    return DISCONNECT;
  }
}

void Graph::set_original_edge_BW( const BaseVertex* source, const BaseVertex* sink, double BW )
{
  map<int, double>::const_iterator pos =
    m_mpEdgeCodeBW.find(get_edge_code(source, sink));

  if (pos != m_mpEdgeCodeBW.end())
  {
    m_mpEdgeCodeBW.at(get_edge_code(source, sink)) = BW;
  }
  else
  {

  }
}

double Graph::get_original_edge_UsedBW( const BaseVertex* source, const BaseVertex* sink )
{
  map<int, double>::const_iterator pos =
    m_mpEdgeCodeUsedBW.find(get_edge_code(source, sink));

  if (pos != m_mpEdgeCodeUsedBW.end())
  {
    return pos->second;
  }
  else
  {
    return DISCONNECT;
  }
}

void Graph::set_original_edge_UsedBW( const BaseVertex* source, const BaseVertex* sink, double UsedBW )
{
  map<int, double>::const_iterator pos =
    m_mpEdgeCodeUsedBW.find(get_edge_code(source, sink));

  if (pos != m_mpEdgeCodeUsedBW.end())
  {
    m_mpEdgeCodeUsedBW.at(get_edge_code(source, sink)) = UsedBW;
  }
  else
  {

  }
}

int Graph::get_edge_num()
{
  return this->m_nEdgeNum;
}

int Graph::get_vertex_num()
{
  return this->m_nVertexNum;
}

int Graph::get_graphID()
{
  return this->m_graphID;
}

double Graph::get_path_BW(BasePath* path)
{
  double BW_path = get_edge_BW(path->GetVertex(0),path->GetVertex(1));
  for (int i = 1; i < path->length()-1; ++i)
  {
    double edge_BW = get_edge_BW(path->GetVertex(i),path->GetVertex(i+1));
    if (BW_path > edge_BW)
    {
      BW_path = edge_BW;
    }
  }
  return BW_path;
}

void Graph::UpdateCommodity(Commodity* commodity)
{
  for (vector<Commodity*>::iterator it_commodity = m_vCommodity.begin(); it_commodity != m_vCommodity.end(); ++it_commodity)
  {
    if ((*it_commodity)->sink_ == commodity->sink_ && (*it_commodity)->source_ == commodity->source_)
    {
      (*it_commodity)->demand_ += commodity->demand_;
      return;
    }
  }
  m_vCommodity.push_back(commodity);
}


void Graph::ComputeTopoTable()
{
  for (int j = 0; j < m_vtBorderVertices.size(); ++j)
  {
    for (int i = 0; i < get_vertex_num(); ++i)
    {
      BaseVertex* p_source = m_vtBorderVertices[j];
      BaseVertex* p_sink = get_vertex(i,get_graphID());
      if(p_source != p_sink)
      {
        // the destination is different from the source
        vector<BasePath*> kPathList;
        YenTopKShortestPathsAlg yenAlg(*this);
        yenAlg.get_shortest_paths(p_source,p_sink, kK, kPathList);
        if (kPathList.size()==0)
        {
          continue;
        }
        for (vector<BasePath*>::iterator it_path = kPathList.begin();
             it_path != kPathList.end(); ++it_path)
        {
          (*it_path)->set_BW(get_path_BW((*it_path)));
          /*
           * if the source vertex is a border vertex, then insert to the topoTable
           */
          vector<int> ASPath;
          ASPath.push_back(m_graphID);
          m_TopoTable.Insert(TopoTableEntry(p_source,p_sink,p_source,
                                            (*it_path)->Weight(), (*it_path)->get_BW(),ASPath));
        }
        m_vPathTable.insert(m_vPathTable.end(),kPathList.begin(),kPathList.end());
      }
      else
      {
        /*
         * if the source vertex is a border vertex, then insert to the topoTable
         */
        vector<BaseVertex*> vertex_list;
        vertex_list.push_back(p_source);
        vertex_list.push_back(p_sink);
        m_vPathTable.push_back(new BasePath(vertex_list, 0));
        vector<int> ASPath;
        ASPath.push_back(m_graphID);
        m_TopoTable.Insert(TopoTableEntry(p_source,p_sink,p_source,0, max_BW, ASPath));
      }
    }
  }
}
/*
 * Update the topology table in the naive protocol
 * In this protocol, there is no AdvertisedTable, so no need to update it
 */
bool Graph::UpdateTopoTableNaive(TopoTableEntry entry)
{
  /*
   * the updated entry should be from the neighboring AS
   */
  if (entry.m_source == entry.m_next)
  {
    cout << "wrong happen in Graph::UpdateTopoTable" << endl;
    return false;
  }

  vector<TopoTableEntry>::iterator it_max_Entry = m_TopoTable.m_vEntry.begin();
  int num_path = 0;
  for (vector<TopoTableEntry>::iterator it = m_TopoTable.m_vEntry.begin();
       it != m_TopoTable.m_vEntry.end(); ++it)
  {
    if (it->m_source == entry.m_source && it->m_sink == entry.m_sink
        && it->m_next == entry.m_next)
    {
      return false;
    }
    if (it->m_source == entry.m_source && it->m_sink == entry.m_sink)
    {
      num_path ++;
      if (it_max_Entry->m_weight < it->m_weight)
      {
        it_max_Entry = it;
      }
    }
  }
  if (num_path < m_TopoTable.m_nK)
  {
    entry.m_vASPath.push_back(m_graphID);
    m_TopoTable.Insert(entry);
    return true;
  }
  else if (it_max_Entry->m_weight > entry.m_weight)
  {
    m_TopoTable.Delete(it_max_Entry);
    entry.m_vASPath.push_back(m_graphID);
    m_TopoTable.Insert(entry);
    return true;
  }
  return false;
}
/*
 * the input entry has already included the edge between two border switches
 */
bool Graph::UpdateTopoTable(TopoTableEntry entry)
{
  /*
   * the updated entry should be from the neighboring AS
   */
  if (entry.m_source == entry.m_next)
  {
    cout << "wrong happen in Graph::UpdateTopoTable" << endl;
  }

  vector<TopoTableEntry>::iterator it_max_Entry = m_TopoTable.m_vEntry.begin();
  int num_path = 0;
  for (vector<TopoTableEntry>::iterator it = m_TopoTable.m_vEntry.begin();
       it != m_TopoTable.m_vEntry.end(); ++it)
  {
    if (it->m_source == entry.m_source && it->m_sink == entry.m_sink
        && it->m_next == entry.m_next)
    {
      return false;
    }
    if (it->m_source == entry.m_source && it->m_sink == entry.m_sink)
    {
      num_path ++;
      if (it_max_Entry->m_weight < it->m_weight)
      {
        it_max_Entry = it;
      }
    }
  }
  if (num_path < m_TopoTable.m_nK)
  {
    entry.m_vASPath.push_back(m_graphID);
    m_TopoTable.Insert(entry);
    UpdateAdvertisedTable(entry);
    return true;
  }
  else if (it_max_Entry->m_weight > entry.m_weight)
  {
    m_TopoTable.Delete(it_max_Entry);
    entry.m_vASPath.push_back(m_graphID);
    m_TopoTable.Insert(entry);
    UpdateAdvertisedTable(entry);
    return true;
  }
  return false;
}

void Graph::UpdateAdvertisedTable(TopoTableEntry entry)
{
  // find whether the sink of the entry is in the advertise table
  bool isFind = false;
  for (vector<TopoTableEntry>::iterator it = m_AdvertisedTopoTable.m_vEntry.begin();
       it != m_AdvertisedTopoTable.m_vEntry.end(); ++it)
  {
    if (it->m_sink == entry.m_sink)
    {
      isFind = true;
      break;
    }
  }
  if (isFind)
  {
    vector<TopoTableEntry> v_insertTopoEntry;
    for (vector<BaseVertex*>::iterator it_border = m_vtBorderVertices.begin();
         it_border != m_vtBorderVertices.end(); ++it_border)
    {
      vector<vector<TopoTableEntry>::iterator> it_BorderToBorder;
      vector<vector<TopoTableEntry>::iterator> it_BorderToSink;
      double k_max = 0;
      vector<TopoTableEntry>::iterator it_k_maxB2B; // the iterator of entry from the border node to another border node
      vector<TopoTableEntry>::iterator it_k_maxB2S; // the iterator of entry from the border to the sink
      int index_max;
      // find the k-shortest path from the border switch to the destination
      for (vector<TopoTableEntry>::iterator it = m_TopoTable.m_vEntry.begin();
           it != m_TopoTable.m_vEntry.end(); ++it)
      {
        if (it->m_sink == entry.m_sink)
        {
          for (vector<TopoTableEntry>::iterator it1 = m_TopoTable.m_vEntry.begin();
               it1 != m_TopoTable.m_vEntry.end(); ++it1)
          {
            if (it1->m_sink == it->m_source && it1->m_source == *it_border)
            {
              double totalWeight = it1->m_weight + it->m_weight;
              if (it_BorderToBorder.size() < kMaxComputePath)
              {
                it_BorderToBorder.push_back(it1);
                it_BorderToSink.push_back(it);
                if (totalWeight > k_max)
                {
                  k_max = totalWeight;
                  it_k_maxB2B = it1;
                  it_k_maxB2S = it;
                  index_max = it_BorderToBorder.size() - 1;
                }
              }
              else
              {
                if (totalWeight < k_max)
                {
                  //delete the k-max path, and insert the new path
                  it_BorderToBorder[index_max] = it1;
                  it_BorderToSink[index_max] = it;
                  // find the k-max path in the new table
                  k_max = 0;
                  for (int i = 0; i < kMaxComputePath; ++i)
                  {
                    if (it_BorderToBorder[i]->m_weight + it_BorderToSink[i]->m_weight > k_max)
                    {
                      k_max = it_BorderToBorder[i]->m_weight + it_BorderToSink[i]->m_weight;
                      it_k_maxB2B = it_BorderToBorder[i];
                      it_k_maxB2S = it_BorderToSink[i];
                      index_max = i;
                    }
                  }
                }
              }
            }
          }
        }
      }

      if (it_BorderToBorder.size() == 0)
        continue;

      map<int,BaseVertex*> indexToVertex;
      map<BaseVertex*, int> vertexToIndex;
      int numVertices = 0;
      for (vector<vector<TopoTableEntry>::iterator>::iterator it_de = it_BorderToBorder.begin();
           it_de != it_BorderToBorder.end(); ++it_de)
      {
        for (vector<BasePath*>::iterator it_path = m_vPathTable.begin();
             it_path != m_vPathTable.end(); ++it_path)
        {
          if ((*it_path)->GetVertex(0) == (*(*it_de)).m_source &&
              (*it_path)->GetLastVertex() == (*(*it_de)).m_sink &&
              (*it_path)->Weight() == (*(*it_de)).m_weight)
          {
            for (int i = 0; i < (*it_path)->length(); ++i)
            {
              bool isAppear = false;
              for (map<BaseVertex*, int>::iterator it_vertex = vertexToIndex.begin();
                   it_vertex != vertexToIndex.end(); ++it_vertex)
              {
                if (it_vertex->first == (*it_path)->GetVertex(i))
                {
                  isAppear = true;
                  break;
                }
              }
              if (!isAppear)
              {
                if (vertexToIndex.find((*it_path)->GetVertex(i)) == vertexToIndex.end())
                {
                  vertexToIndex[(*it_path)->GetVertex(i)] = numVertices;
                  indexToVertex[numVertices] = (*it_path)->GetVertex(i);
                  numVertices++;
                }
              }
            }
            break;
          }
        }
      }

      for (vector<vector<TopoTableEntry>::iterator>::iterator it_de = it_BorderToSink.begin();
           it_de != it_BorderToSink.end(); ++it_de)
      {
        bool isAppear = false;
        for (map<BaseVertex*, int>::iterator it_vertex = vertexToIndex.begin();
             it_vertex != vertexToIndex.end(); ++it_vertex)
        {
          if (it_vertex->first == (*(*it_de)).m_next)
          {
            isAppear = true;
            break;
          }
        }
        if (!isAppear)
        {
          if (vertexToIndex.find((*(*it_de)).m_next) == vertexToIndex.end())
          {
            vertexToIndex[(*(*it_de)).m_next] = numVertices;
            indexToVertex[numVertices] = (*(*it_de)).m_next;
            numVertices++;
          }
        }
      }


      if (vertexToIndex.find(entry.m_sink) == vertexToIndex.end())
      {
        vertexToIndex[entry.m_sink] = numVertices;
        indexToVertex[numVertices] = entry.m_sink;
        numVertices++;
      }

      edge** graph = new edge*[numVertices];
      for (int i = 0; i < numVertices; ++i)
      {
        graph[i] = new edge[numVertices];
      }
      for (int i = 0; i < numVertices; ++i)
      {
        for (int j = 0; j < numVertices; ++j)
        {
          graph[i][j].flow = 0;
          graph[i][j].capacity = 0;
          graph[i][j].cost = 0;
        }
      }

      for (vector<vector<TopoTableEntry>::iterator>::iterator it_de = it_BorderToBorder.begin();
           it_de != it_BorderToBorder.end(); ++it_de)
      {
        for (vector<BasePath*>::iterator it_path = m_vPathTable.begin();
             it_path != m_vPathTable.end(); ++it_path)
        {
          if ((*it_path)->GetVertex(0) == (*(*it_de)).m_source &&
              (*it_path)->GetLastVertex() == (*(*it_de)).m_sink &&
              (*it_path)->Weight() == (*(*it_de)).m_weight)
          {
            for (int i = 0; i < (*it_path)->length()-1; ++i)
            {
              BaseVertex* p_source = (*it_path)->GetVertex(i);
              BaseVertex* p_sink = (*it_path)->GetVertex(i+1);
              if (p_source == p_sink)
                continue;
              graph[vertexToIndex.at(p_source)][vertexToIndex.at(p_sink)].capacity =
                get_edge_BW(p_source,p_sink);
              graph[vertexToIndex.at(p_source)][vertexToIndex.at(p_sink)].cost =
                get_edge_weight(p_source,p_sink);
            }
            break;
          }
        }
      }

      for (vector<vector<TopoTableEntry>::iterator>::iterator it_de = it_BorderToSink.begin();
           it_de != it_BorderToSink.end(); ++it_de)
      {
        for (vector<TopoTableEntry>::iterator it_first = m_TopoTable.m_vEntry.begin();
             it_first != m_TopoTable.m_vEntry.end(); ++it_first)
        {
          if (it_first->m_source == (*(*it_de)).m_source &&
              it_first->m_sink == (*(*it_de)).m_next)
          {
            BaseVertex* p_source = it_first->m_source;
            BaseVertex* p_sink = it_first->m_sink;
            if (p_source == p_sink)
              continue;
            graph[vertexToIndex.at(p_source)][vertexToIndex.at(p_sink)].capacity = it_first->m_BW;
            graph[vertexToIndex.at(p_source)][vertexToIndex.at(p_sink)].cost = it_first->m_weight;

            if (p_sink == entry.m_sink)
              continue;
            graph[vertexToIndex.at(p_sink)][vertexToIndex.at(entry.m_sink)].capacity = (*(*it_de)).m_BW;
            graph[vertexToIndex.at(p_sink)][vertexToIndex.at(entry.m_sink)].cost = (*(*it_de)).m_weight - it_first->m_weight;

          }
        }
      }

      FordFulkersonAlg myAlg(graph, numVertices,
                             vertexToIndex.at(*it_border), vertexToIndex.at(entry.m_sink));
      double max_flow = 0;
      double max_cost = 0;
      myAlg.MaxFlow(max_flow,max_cost);

      vector<int> v_ASpath((*(*it_BorderToSink.begin())).m_vASPath.begin(),
                           (*(*it_BorderToSink.begin())).m_vASPath.end());
      vector<int>::iterator it_path;
      for (vector<vector<TopoTableEntry>::iterator>::iterator it_de = it_BorderToSink.begin() + 1;
           it_de != it_BorderToSink.end(); ++it_de)
      {
        vector<int> v_temp_path(N_AS);
        it_path = set_intersection(v_ASpath.begin(), v_ASpath.end(),
                                   (*(*it_de)).m_vASPath.begin(), (*(*it_de)).m_vASPath.end(),
                                   v_temp_path.begin());
        v_temp_path.resize(it_path-v_temp_path.begin());
        v_ASpath = v_temp_path;
      }


      v_insertTopoEntry.push_back(TopoTableEntry(*it_border,entry.m_sink,NULL,max_cost,max_flow,v_ASpath));
      //m_AdvertisedTopoTable.Insert(TopoTableEntry(*it_border,entry.m_sink,NULL,max_cost,max_flow,v_ASpath));
    }
    for (vector<TopoTableEntry>::iterator it = v_insertTopoEntry.begin();
         it != v_insertTopoEntry.end(); ++it)
    {
      for (vector<TopoTableEntry>::iterator it_ad = m_AdvertisedTopoTable.m_vEntry.begin();
           it_ad != m_AdvertisedTopoTable.m_vEntry.end(); ++it_ad)
      {
        if (it_ad->m_source == it->m_source && it_ad->m_sink == it->m_sink)
        {
          m_AdvertisedTopoTable.Delete(it_ad);
          m_AdvertisedTopoTable.Insert(*it);
        }
      }
    }
  }
  else
  {
    vector<TopoTableEntry> v_insertTopoEntry;
    for (vector<TopoTableEntry>::iterator it = m_AdvertisedTopoTable.m_vEntry.begin();
         it != m_AdvertisedTopoTable.m_vEntry.end(); ++it)
    {
      if (it->m_sink == entry.m_source)
      {
        TopoTableEntry insertEntry = entry;
        insertEntry.m_source = it->m_source;
        insertEntry.m_next = NULL;
        insertEntry.m_weight += it->m_weight;
        insertEntry.m_BW = min(it->m_BW,entry.m_BW);
        v_insertTopoEntry.push_back(insertEntry);
      }
    }
    for (vector<TopoTableEntry>::iterator it = v_insertTopoEntry.begin();
         it != v_insertTopoEntry.end(); ++it)
    {
      m_AdvertisedTopoTable.Insert(*it);
    }
  }
}

/*
 * Compute the advertised table according to the topology table during initializing the AS based on the following way:
 * 1) for each entry e(source, next, sink) in topoTable
 * 2) if source=next, this is interior path,
 *    2.1) find all the paths in routing table where source and sink matches with e
 *    2.2) extract the nodes and edges in these paths to construct a new subgraph G'
 *    2.3) in G', compute the max flow from source to sink by FordFulkerson Alg
 *    2.4) based on the maximal flow, compute the cost of the flows c= sum{e in G'}{f(e)c(e)}/maxFlow,
 *         where f(e) is the flow on the edge e, c(e) is the cost of passing e, maxFlow is the size of max flow in G'
 *    2.5) insert the entry(source, sink, c, maxFlow) into the advertised table
 * 3) if source!=next, this is cross-AS path,
 *    3.1) sum the BW of all these entries whose source and sink matches with e, say sum_BW
 *    3.2) compute the expected weight, say aver_weight, where the probability of each entry is its BW/sum_BW
 *    3.3) Insert the entry (source, sink, aver_weight, sum_BW) into the advertised table
 */

void Graph::ComputeAdvertisedTable()
{

  for (vector<TopoTableEntry>::iterator it_topoEntry = m_TopoTable.m_vEntry.begin();
       it_topoEntry != m_TopoTable.m_vEntry.end(); ++it_topoEntry)
  {
    bool isComputed = false;
    for (vector<TopoTableEntry>::iterator it_advertise = m_AdvertisedTopoTable.m_vEntry.begin();
         it_advertise != m_AdvertisedTopoTable.m_vEntry.end(); ++it_advertise)
    {
      if (it_advertise->m_source == it_topoEntry->m_source &&
          it_advertise->m_sink == it_topoEntry->m_sink)
      {
        isComputed = true;
        break;
      }
    }
    if (!isComputed)
    {
      if (it_topoEntry->m_source == it_topoEntry->m_next)
      {
        map<int,BaseVertex*> indexToVertex;
        map<BaseVertex*, int> vertexToIndex;
        int numVertices = 0;

        for (vector<BasePath*>::iterator it_path = m_vPathTable.begin();
             it_path != m_vPathTable.end(); ++it_path)
        {
          if ((*it_path)->GetVertex(0) == it_topoEntry->m_source &&
              (*it_path)->GetLastVertex() == it_topoEntry->m_sink)
          {
            for (int i = 0; i < (*it_path)->length(); ++i)
            {
              bool isAppear = false;
              for (map<BaseVertex*, int>::iterator it_vertex = vertexToIndex.begin();
                   it_vertex != vertexToIndex.end(); ++it_vertex)
              {
                if (it_vertex->first == (*it_path)->GetVertex(i))
                {
                  isAppear = true;
                  break;
                }
              }
              if (!isAppear)
              {
                vertexToIndex[(*it_path)->GetVertex(i)] = numVertices;
                indexToVertex[numVertices] = (*it_path)->GetVertex(i);
                numVertices++;
              }
            }
          }
        }
        if (numVertices == 1)
        {
          m_AdvertisedTopoTable.Insert(TopoTableEntry(it_topoEntry->m_source, it_topoEntry->m_sink,
                                       NULL, 0, max_BW, it_topoEntry->m_vASPath));
        }
        else if (numVertices == 2)
        {
          m_AdvertisedTopoTable.Insert(TopoTableEntry(it_topoEntry->m_source, it_topoEntry->m_sink,
                                       NULL, it_topoEntry->m_weight, it_topoEntry->m_BW,
                                       it_topoEntry->m_vASPath));
        }
        else
        {
          edge** graph = new edge*[numVertices];
          for (int i = 0; i < numVertices; ++i)
          {
            graph[i] = new edge[numVertices];
          }
          for (int i = 0; i < numVertices; ++i)
          {
            for (int j = 0; j < numVertices; ++j)
            {
              graph[i][j].flow = 0;
              graph[i][j].capacity = 0;
              graph[i][j].cost = 0;
            }
          }

          for (vector<BasePath*>::iterator it_path = m_vPathTable.begin();
               it_path != m_vPathTable.end(); ++it_path)
          {
            if ((*it_path)->GetVertex(0) == it_topoEntry->m_source &&
                (*it_path)->GetLastVertex() == it_topoEntry->m_sink)
            {
              for (int i = 0; i < (*it_path)->length()-1; ++i)
              {
                BaseVertex* p_source = (*it_path)->GetVertex(i);
                BaseVertex* p_sink = (*it_path)->GetVertex(i+1);
                graph[vertexToIndex.at(p_source)][vertexToIndex.at(p_sink)].capacity =
                  get_edge_BW(p_source,p_sink);
                graph[vertexToIndex.at(p_source)][vertexToIndex.at(p_sink)].cost =
                  get_edge_weight(p_source,p_sink);
              }
            }
          }
          FordFulkersonAlg myAlg(graph, numVertices,
                                 vertexToIndex.at(it_topoEntry->m_source), vertexToIndex.at(it_topoEntry->m_sink));
          double max_flow = 0;
          double max_cost = 0;
          myAlg.MaxFlow(max_flow,max_cost);
          m_AdvertisedTopoTable.Insert(TopoTableEntry(it_topoEntry->m_source, it_topoEntry->m_sink,
                                       NULL, max_cost, max_flow, it_topoEntry->m_vASPath));
        }
      }
      else
      {
        cout << "wrong in Graph::ComputeAdvertisedTable(): the next and the source are in different ASes" << endl;
      }
    }
  }
}


int Graph::get_vertex_code(BaseVertex* vertex)
{
  for (map<int, BaseVertex*>::iterator it = m_mpVertexIndex.begin(); it != m_mpVertexIndex.end(); ++it)
  {
    if (it->second == vertex)
    {
      return it->first;
    }
  }
  return -1;
}


void Graph::add_vertex(BaseVertex* vertex)
{
  BaseVertex* vertex_pt = NULL;

  if (vertex->getGraphID() != m_graphID)
  {
    // add a vertex from another graph into the graph
    int ori_id = vertex->getID() + vertex->getGraphID() * max_Nodes;
    if (mpNodeID.find(ori_id) != mpNodeID.end())
    {
      /*
       * if the vertex has already been added into the graph, just ignore
       */
      return;
    }
    else
    {
      /*
       * if the vertex is not added yet, we create a new vertex with a new node id
       */
      int vertex_id = m_vtVertices.size();
      vertex_pt = new BaseVertex();
      vertex_pt->setID(vertex_id);
      vertex_pt->setGraphID(this->get_graphID());
      m_vtVertices.push_back(vertex_pt);
      m_mpVertexIndex.insert(pair<int,BaseVertex*>(vertex_id, vertex_pt));

      int ori_id = vertex->getID() + vertex->getGraphID() * max_Nodes;
      mpNodeID.insert(pair<int,int>(ori_id,vertex_id));
    }
  }
  else
  {
    // add a new vertex to the graph
    const map<int, BaseVertex*>::iterator pos = m_mpVertexIndex.find(vertex->getID());
    if (pos != m_mpVertexIndex.end())
    {
      /*
       * if there already exits a vertex with the same node id,
       * we create a new vertex with a new node id
       */
      int vertex_id = m_vtVertices.size();
      vertex_pt = new BaseVertex();
      vertex_pt->setID(vertex_id);
      vertex_pt->setGraphID(vertex->getGraphID());
      m_vtVertices.push_back(vertex_pt);
      m_mpVertexIndex.insert(pair<int,BaseVertex*>(vertex_id, vertex_pt));

    }
    else
    {
      m_vtVertices.push_back(vertex);
      m_mpVertexIndex.insert(pair<int,BaseVertex*>(vertex->getID(), vertex));
    }
  }
}

void Graph::add_edge(BaseVertex* start_vertex_pt, BaseVertex* end_vertex_pt, double edge_weight, double edge_BW)
{
  m_mpEdgeCodeWeight[get_edge_code(start_vertex_pt, end_vertex_pt)] = edge_weight;

  m_mpEdgeCodeBW[get_edge_code(start_vertex_pt, end_vertex_pt)] = edge_BW;
  m_mpEdgeCodeUsedBW[get_edge_code(start_vertex_pt, end_vertex_pt)] = 0;

  get_vertex_set_pt(end_vertex_pt, m_mpFaninVertices)->insert(start_vertex_pt);

  get_vertex_set_pt(start_vertex_pt, m_mpFanoutVertices)->insert(end_vertex_pt);
}
//-------------------------------------------------------------------
