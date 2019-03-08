///////////////////////////////////////////////////////////////////////////////
///  Graph.h
///  <TODO: insert file description here>
///
///  @remarks <TODO: insert remarks here>
///
///  @author Yan Qi @date 8/18/2010
///
///  $Id: Graph.h 65 2010-09-08 06:48:36Z yan.qi.asu $
///////////////////////////////////////////////////////////////////////////////


#pragma once
#include "TopoTable.h"
#include "Commodity.h"
#include <map>
#include <set>

using namespace std;

const double max_BW = 1000000;
const int max_Nodes = 100; // the maximal nodes per AS
const int max_AS = 20; //  the maximal number of ASes

class Path : public BasePath
{
public:

  Path(const std::vector<BaseVertex*>& vertex_list, double weight):BasePath(vertex_list,weight) {}

  // display the content
  void PrintOut(std::ostream& out_stream) const
  {
    out_stream << "Cost: " << m_dWeight << " Length: " << m_vtVertexList.size() << std::endl;
    for(std::vector<BaseVertex*>::const_iterator pos=m_vtVertexList.begin(); pos!=m_vtVertexList.end(); ++pos)
    {
      out_stream << (*pos)->getID() << " ";
    }
    out_stream << std::endl <<  "*********************************************" << std::endl;
  }
};



class Graph
{
public: // members

  const static double DISCONNECT;
  static int Graph_ID;

  typedef set<BaseVertex*>::iterator VertexPtSetIterator;
  typedef map<BaseVertex*, set<BaseVertex*>*>::iterator BaseVertexPt2SetMapIterator;

  //hemin yang
  vector<BaseVertex*> m_vtBorderVertices;
  map<pair<BaseVertex*, BaseVertex*>, pair<double, double>> m_BorderEdges; // m_BorderEdges[<src, sink>]=<bw, cost>
  TopoTable m_TopoTable;  // store the information from the border switch to all the other switches in the whole network
  TopoTable m_AdvertisedTopoTable;
  vector<BasePath*> m_vPathTable; // the k shortest pathes from each vertex in the AS to all the other vertex
  double m_maxBW;
  map<int, set<Commodity*>*> m_mpEdgeCodeCommodity;
  map<int, double> m_mpEdgeCodeUsedBW;
  vector<Commodity*> m_vCommodity;
  map<int,int> mpNodeID; // if the graph is constructed based on several ASes,
  // this map is used to map the node in the AS to the node id of the constructed graph.
  // mpNodeID[graph_id*max_Nodes+vertex_id] = new vertex_id in the constructed graph

protected: // members

  // Basic information
  map<BaseVertex*, set<BaseVertex*>*> m_mpFanoutVertices;
  map<BaseVertex*, set<BaseVertex*>*> m_mpFaninVertices;
  map<int, double> m_mpEdgeCodeWeight;
  vector<BaseVertex*> m_vtVertices;

  map<int, double> m_mpEdgeCodeBW;
  //map<int, double> m_mpEdgeCodeUsedBW;


  int m_nEdgeNum;
  int m_nVertexNum;
  int m_graphID; // the id of the graph

  map<int, BaseVertex*> m_mpVertexIndex;

  // Members for graph modification
  set<int> m_stRemovedVertexIds;
  set<pair<int,int> > m_stRemovedEdge;

public:

  // Constructors and Destructor
  Graph(const string& file_name, int kPath = 0);
  Graph(const Graph* rGraph);
  Graph(int kPath = 0):m_nEdgeNum(0),m_nVertexNum(0),m_graphID(0), m_maxBW(0){
    m_TopoTable.m_nK = kPath;
    m_AdvertisedTopoTable.m_nK = kPath;
  }
  /*---hemin-----
   * \brief this constructor is to construct a network containing several ASs
   * \param file_name, the txt file containing all the information needed for the internet
   * \param graphID, the id of this graph
   */
  Graph(const string& file_name, bool buildEntireNetwork);

  /*
   * construct a virtual graph based on the topoTable
   */
  void ConstructVirtualGraph(Graph* graph, const string& AS_file);

  /*
   * construct a global view of the entire network
   */
  void ConstructGlobalView(Graph* graph, const string& AS_file);


  ~Graph(void);

  void clear();
  void clearCommodities();

  //----hemin------
  int get_vertex_num();
  int get_edge_num();
  int get_graphID();
  void set_graphID(int graphID);
  void set_vertex_num(int numVertex);
  /*
   * if the graph is constructed based on several graphs, this function is used to
   * get the original node id and graph id of the node.
   */
  void get_original_id(int newNodeID, int& oriNodeID, int& oriGraphID);

  void printGraph();
  void printTopoTable();
  void printAdvertisedTable();
  void printTableEntry(TopoTableEntry& entry);
  void printPath(const BasePath* path);

  /*
   * get the bandwidth of the bottleneck of one path
   */
  double get_path_BW(BasePath* path);

  /*
   * find whether the vertex is a border vertex
   */
  bool FindBorderVertex(BaseVertex* vertex);

  /*
   * Insert a new commodity into the vector of commodity
   * 1) check whether there is a commodity share same source and sink as the inserted one
   * 2) if no, insert the new commodity to the vector
   * 3) else, plus the demand of the new commodity to the old one
   */
  void UpdateCommodity(Commodity* commodity);
  /*
   * compute the topology table
   * 1) compute the k shortest paths from the border switch to all the other switches in the AS
   */
  void ComputeTopoTable();
  /*
   * given an topology table entry from neighboring AS, the AS update its own topology table in the following way
   * 1) find the entries having the same destination as the received entry
   * 2) if the number of such entries < kPath, insert the elements into the topotable
   * 3) if not, find the entry with maximal weight
   * 4) if the maximal weight > the received entry's weight,
   *    delete the entry with the maximal weight and insert the received entry
   * \param element the input topotable entry
   * \return return true if the topotable is updated, otherwise return false
   */
  bool UpdateTopoTable(TopoTableEntry element);
  bool UpdateTopoTableNaive(TopoTableEntry element);

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
   *    3.1) sum the BW of all these entries in the topoTable whose source and sink matches with e, say sum_BW
   *    3.2) compute the expected weight, say aver_weight, where the probability of each entry is its BW/sum_BW
   *    3.3) Insert the entry (source, sink, aver_weight, sum_BW) into the advertised table
   */
  void ComputeAdvertisedTable();

  /*
   * Compute the advertisement table for the navie protocol (BGP-Addpath) during initialization. The advertisement table is exactly same as the topology table
   */
  void ComputeAdvertisedTableNaive();

  /*
   * Update the advertised table if the entries of topology table for source-sink pair are updated
   * 1) find the entries whose source and sink match the input
   * 2) recompute the BW and weight for this pair according to the way in ComputeAdvertisedTable()
   * 3) if there exists an entry with same source and sink as the inputs, update the entry's BW and weight
   * 4) otherwise, insert a new entry (source, sink, aver_weight, BW) to the AdvertisedTable
   */
  void UpdateAdvertisedTable(TopoTableEntry element);

  //---------------------------

  BaseVertex* get_vertex(int node_id, int graphID);

  int get_edge_code(const BaseVertex* start_vertex_pt, const BaseVertex* end_vertex_pt) const;

  set<BaseVertex*>* get_vertex_set_pt(BaseVertex* vertex_, map<BaseVertex*, set<BaseVertex*>*>& vertex_container_index);

  double get_original_edge_weight(const BaseVertex* source, const BaseVertex* sink);
  //---hemin----
  /*
   * get the vertex index of a specific vertex
   * If no such vertex exist, return -1
   */
  int get_vertex_code(BaseVertex* vertex);

  double get_original_edge_BW(const BaseVertex* source, const BaseVertex* sink);
  void set_original_edge_BW(const BaseVertex* source, const BaseVertex* sink, double BW);
  double get_original_edge_UsedBW(const BaseVertex* source, const BaseVertex* sink);
  void set_original_edge_UsedBW(const BaseVertex* source, const BaseVertex* sink, double UsedBW);
  //--------------

  double get_edge_weight(const BaseVertex* source, const BaseVertex* sink);

  //----hemin-----------
  double get_edge_BW(const BaseVertex* source, const BaseVertex* sink);
  double get_edge_BW(int code);
  void set_edge_BW(const BaseVertex* source, const BaseVertex* sink, double BW);
  void set_edge_BW(int code, double BW);

  double get_edge_UsedBW(const BaseVertex* source, const BaseVertex* sink);
  double get_edge_UsedBW(int code);
  void set_edge_UsedBW(const BaseVertex* source, const BaseVertex* sink, double UsedBW);
  void set_edge_UsedBW(int code, double UsedBW);

  //--------------------

  void get_adjacent_vertices(BaseVertex* vertex, set<BaseVertex*>& vertex_set);
  void get_precedent_vertices(BaseVertex* vertex, set<BaseVertex*>& vertex_set);

  /// Methods for changing graph
  void remove_edge(const pair<int,int> edge)
  {
    m_stRemovedEdge.insert(edge);
  }

  void remove_vertex(const int vertex_id)
  {
    m_stRemovedVertexIds.insert(vertex_id);
  }

  void recover_removed_edges()
  {
    m_stRemovedEdge.clear();
  }

  void recover_removed_vertices()
  {
    m_stRemovedVertexIds.clear();
  }

  void recover_removed_edge(const pair<int,int> edge)
  {
    m_stRemovedEdge.erase(m_stRemovedEdge.find(edge));
  }

  void recover_removed_vertex(int vertex_id)
  {
    m_stRemovedVertexIds.erase(m_stRemovedVertexIds.find(vertex_id));
  }

  //--hemin---
  void add_vertex(BaseVertex* vertex);
  void add_edge(BaseVertex* start_vertex_pt, BaseVertex* end_vertex_pt, double edge_weight, double edge_BW);

private:
  void _import_from_file(const std::string& file_name, int kPath);

};

class InterGraph
{

public:
  const static double INTER_DISCONNECT;
  Graph* m_left;
  Graph* m_right;
  map<BaseVertex*, set<BaseVertex*>*> m_mpFanoutVertices;
  map<BaseVertex*, set<BaseVertex*>*> m_mpFaninVertices;
  map<int, double> m_mpEdgeCodeWeight;
  map<int, double> m_mpEdgeCodeBW;
  map<int, double> m_mpEdgeCodeUsedBW;

  map<int, BaseVertex*> m_mpVertexIndex;
  vector<BaseVertex*> m_vtVertices;

  int m_nEdgeNum;
  int m_nVertexNum;


  typedef set<BaseVertex*>::iterator VertexPtSetIterator;
  typedef map<BaseVertex*, set<BaseVertex*>*>::iterator BaseVertexPt2SetMapIterator;

  InterGraph(Graph* left_graph, Graph* right_graph, const string& file_name);
  bool find_vertex(BaseVertex*);
  /*
   * get the vertices directly connected from a vertex
   */
  void get_out_vertices(BaseVertex* vertex, set<BaseVertex*>& vertex_set);
  /*
   * get the vertices directly connected to a vertex
   */
  void get_in_vertices(BaseVertex* vertex, set<BaseVertex*>& vertex_set);

  double get_edge_weight(const BaseVertex* source, const BaseVertex* sink);

  double get_edge_BW(const BaseVertex* source, const BaseVertex* sink);
  double get_edge_BW(int code);

  void set_edge_BW(const BaseVertex* source, const BaseVertex* sink, double BW);
  void set_edge_BW(int code, double BW);

  double get_edge_UsedBW(const BaseVertex* source, const BaseVertex* sink);
  double get_edge_UsedBW(int code);

  void set_edge_UsedBW(const BaseVertex* source, const BaseVertex* sink, double UsedBW);
  void set_edge_UsedBW(int code, double UsedBW);

  set<BaseVertex*>* get_vertex_set_pt(BaseVertex* vertex_, map<BaseVertex*, set<BaseVertex*>*>& vertex_container_index);
  int get_edge_code(const BaseVertex* start_vertex_pt, const BaseVertex* end_vertex_pt) const;
  BaseVertex* get_vertex(int node_id, int graphID);
  void clear();

};
