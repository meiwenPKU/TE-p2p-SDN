/*
 * FordFulkersonAlg.cpp
 *
 *  Created on: Apr 5, 2016
 *      Author: Hemin
 */
#include "FordFulkersonAlg.h"
#include <memory.h>
#include <queue>
#include <iostream>
const double max_BW = 1000000;
using namespace std;

bool FordFulkersonAlg::bfs(int parent[])
{
    // Create a visited array and mark all vertices as not visited
    bool visited[m_numVertices];
    memset(visited, 0, sizeof(visited));

    // Create a queue, enqueue source vertex and mark source vertex
    // as visited
    queue <int> q;
    q.push(m_Source);
    visited[m_Source] = true;
    parent[m_Source] = -1;

    // Standard BFS Loop
    while (!q.empty())
    {
        int u = q.front();
        q.pop();

        for (int v=0; v<m_numVertices; v++)
        {
            if (visited[v]==false && m_Rgraph[u][v] > 0)
            {
                q.push(v);
                parent[v] = u;
                visited[v] = true;
            }
        }
    }
    // If we reached sink in BFS starting from source, then return
    // true, else false
    return (visited[m_Sink] == true);
}

void FordFulkersonAlg::MaxFlow(double& max_flow, double& cost)
{
    int u, v;

    // Create a residual graph and fill the residual graph with
    // given capacities in the original graph as residual capacities
    // in residual graph
    for (u = 0; u < m_numVertices; u++)
        for (v = 0; v < m_numVertices; v++)
             m_Rgraph[u][v] = m_graph[u][v].capacity;

    int parent[m_numVertices];  // This array is filled by BFS and to store path
    max_flow = 0;
    // Augment the flow while tere is path from source to sink
    while (bfs(parent))
    {
        // Find minimum residual capacity of the edges along the
        // path filled by BFS. Or we can say find the maximum flow
        // through the path found.
        double path_flow = max_BW;
        for (v=m_Sink; v!=m_Source; v=parent[v])
        {
            u = parent[v];
            path_flow = min(path_flow, m_Rgraph[u][v]);
        }

        // update residual capacities of the edges and reverse edges
        // along the path
        for (v=m_Sink; v != m_Source; v=parent[v])
        {
            u = parent[v];
            m_Rgraph[u][v] -= path_flow;
            m_Rgraph[v][u] += path_flow;
        }
        // Add path flow to overall flow
        max_flow += path_flow;
    }

    // compute the cost of the maximal flow
    for (u = 0; u < m_numVertices; ++u)
    {
    	for (v = 0; v < m_numVertices; ++v)
    	{
    		if (m_Rgraph[u][v] > 0)
    		{
    			if (m_graph[u][v].capacity > 0 & m_graph[v][u].capacity <= 0)
    			{
    				m_graph[u][v].flow = m_graph[u][v].capacity - m_Rgraph[u][v];
    			}
    			else if (m_graph[u][v].capacity <= 0 & m_graph[v][u].capacity > 0)
    			{
    				m_graph[v][u].flow = m_Rgraph[u][v];

    			}
    			else if (m_graph[u][v].capacity > 0 & m_graph[v][u].capacity > 0)
    			{
            cout << "wrong: there are two reverse edges between two nodes" << endl;
    			}
    		}
    	}
    }
    cost = 0;
    for (u = 0; u < m_numVertices; ++u)
    {
    	for (v = 0; v < m_numVertices; ++v)
    	{
    		if (m_graph[u][v].flow > 0)
    		{
    			cost += m_graph[u][v].flow * m_graph[u][v].cost;
    		}
    	}
    }
    cost  = cost / max_flow;
}
