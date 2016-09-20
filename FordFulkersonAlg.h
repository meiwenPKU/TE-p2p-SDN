/*
 * FordFulkersonAlg.h
 *
 *  Created on: Apr 5, 2016
 *      Author: Hemin
 */


#pragma once
#include <cstdlib>
#include <stdio.h>
using namespace std;

struct edge{
	double capacity;
	double cost;
	double flow;
};

class FordFulkersonAlg
{
public:
	int m_Source;
	int m_Sink;
	edge **m_graph; // the adjacent matrix of the graph
	double **m_Rgraph; // the residual graph
	int m_numVertices;
	FordFulkersonAlg(edge** graph, int numVertices, int s, int t):
		m_Source(s),m_Sink(t),m_graph(graph),m_numVertices(numVertices){
		 m_Rgraph = new double*[numVertices];
		 for (int i = 0; i < numVertices; ++i)
		 {
			 m_Rgraph[i] = new double[numVertices];
		 }
	}

	bool bfs(int parent[]);
	/*
	 * maxflow is the max flow from the source to the sink,
	 * and the cost is the cost to transmit these flows in the network
	 */
	void MaxFlow(double& maxflow, double& cost);
};

