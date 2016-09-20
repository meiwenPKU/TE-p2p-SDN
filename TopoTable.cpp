/*
 * TopoTable.cpp
 *
 *  Created on: Apr 15, 2016
 *      Author: Hemin
 */
#include "TopoTable.h"
#include "GraphElements.h"
#include "TopoTableEntry.h"
#include<vector>
using namespace std;
void TopoTable::printTable() const{
	cout << "**********************************" << endl;
	cout << "The number of entries = " << m_nEntry <<
			"; Maximal number of paths for each source-sink pair = " << m_nK <<  endl;
	cout << "source    next    sink    weight    bandwidth    ASPath" << endl;
	for (vector<TopoTableEntry>::const_iterator it = m_vEntry.begin(); it != m_vEntry.end(); ++it)
	{
		cout << "(" << it->m_source->getGraphID() << "," << it->m_source->getID() << ")    "
				<< "(" << it->m_next->getGraphID() << "," << it->m_next->getID() << ")    "
				<< "(" << it->m_sink->getGraphID() << "," << it->m_sink->getID() << ")    "
				<< it->m_weight << "   " << it->m_BW << "   ";
		for (vector<int>::const_iterator it_as = it->m_vASPath.begin(); it_as != it->m_vASPath.end(); ++it_as)
		{
			cout << *it_as << ",";
		}
		cout << endl;
		cout << "**************************************" << endl;
	}
}
