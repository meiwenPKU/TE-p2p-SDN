/*
 * TopoTable.h
 *
 *  Created on: Mar 8, 2016
 *      Author: hemin
 */
#pragma once
#include "GraphElements.h"
#include "TopoTableEntry.h"
#include<vector>
using namespace std;
const int kPath = 3; // the maximum number of paths for each commodity, k may be different for different traffic
                 // (i.e., self->other, self->self, other->self)
class TopoTable{
public:
	int m_nEntry; // the number of entries in the topology table
	int m_nK; // the maximal number of entries for the same (source-sink) pair
    vector<TopoTableEntry> m_vEntry;

	TopoTable(){
		m_nEntry = 0;
		m_nK = kPath;
	}
	TopoTable(const TopoTable& topoTable ){
		m_nEntry = topoTable.m_nEntry;
		m_nK = topoTable.m_nK;
		m_vEntry = topoTable.m_vEntry;
	}
	TopoTable(vector<TopoTableEntry> vEntry){
		m_vEntry.insert(m_vEntry.begin(), vEntry.begin(), vEntry.end());
		m_nEntry = vEntry.size();
		m_nK = kPath;
	}

    void operator=(const TopoTable &topoTable)
	{
		m_nEntry = topoTable.m_nEntry;
		m_nK = topoTable.m_nK;
		m_vEntry = topoTable.m_vEntry;
	}


	void Insert(TopoTableEntry entry){
		m_vEntry.push_back(entry);
		m_nEntry++;
	}
	void Delete(vector<TopoTableEntry>::iterator it){
		m_vEntry.erase(it);
		m_nEntry--;
	}
	void printTable() const{
		cout << "The number of entries = " << m_nEntry <<
				"; Maximal number of paths for each source-sink pair = " << m_nK <<  endl;
		cout << "***************************************" << endl;
		cout << "source    next    sink    weight  BW    ASPath    ExchangeSet" << endl;
		for (vector<TopoTableEntry>::const_iterator it = m_vEntry.begin(); it != m_vEntry.end(); ++it)
		{
			if (it->m_next == NULL)
			{
				cout << "(" << it->m_source->getGraphID() << "," << it->m_source->getID() << ")    "
					 << "(N,N)    "
					 << "(" << it->m_sink->getGraphID() << "," << it->m_sink->getID() << ")    "
					 << it->m_weight << "   " << it->m_BW << "   ";
			}
			else
			{
				cout << "(" << it->m_source->getGraphID() << "," << it->m_source->getID() << ")    "
					 << "(" << it->m_next->getGraphID() << "," << it->m_next->getID() << ")    "
					 << "(" << it->m_sink->getGraphID() << "," << it->m_sink->getID() << ")    "
					 << it->m_weight << "   " << it->m_BW << "   ";
			}

			for (vector<int>::const_iterator it_as = it->m_vASPath.begin(); it_as != it->m_vASPath.end(); ++it_as)
			{
				cout << *it_as << ",";
			}

            cout << "    ";
			for (set<int>::const_iterator it_set = it->m_isExchanged.begin(); it_set != it->m_isExchanged.end(); ++it_set)
			{
				cout << *it_set << ",";
			}
			cout << endl;
	   }
	   cout << "**************************************" << endl;
	}
};



