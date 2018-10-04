/*
 * TopoTableEntry.h
 *
 *  Created on: Mar 15, 2016
 *      Author: Hemin
 */

#ifndef TOPO_TABLE_ENTRY_H
#define TOPO_TABLE_ENTRY_H
#include<stdint.h>
#include <set>
using namespace std;

class TopoTableEntry{
public:
	BaseVertex* m_source;
	BaseVertex* m_sink;
	BaseVertex* m_next; // if m_next != m_source, then this entry is for inter-AS path, and m_source and m_next
	                    // are in different ASes
	vector<int> m_vASPath;
	double m_weight;
	double m_BW;
	set<int> m_isExchanged; // the set of ASes which the entry has been exchanged with
	TopoTableEntry(BaseVertex* source, BaseVertex* sink, BaseVertex* next, double weight, double BW):
		m_source(source), m_sink(sink), m_next(next),m_weight(weight), m_BW(BW){

	}
	TopoTableEntry(BaseVertex* source, BaseVertex* sink, BaseVertex* next, double weight, double BW, vector<int> ASPath):
			m_source(source), m_sink(sink), m_next(next),m_weight(weight), m_BW(BW){
		m_vASPath = ASPath;
    }
	void printEntry() const {
		if (m_next == NULL)
		{
			cout << "(" << m_source->getGraphID() << "," << m_source->getID() << ")    "
				 << "(N,N)    "
				 << "(" << m_sink->getGraphID() << "," << m_sink->getID() << ")    "
				 << m_weight << "   " << m_BW << "   ";
		}
		else
		{
			cout << "(" << m_source->getGraphID() << "," << m_source->getID() << ")    "
				 << "(" << m_next->getGraphID() << "," << m_next->getID() << ")    "
				 << "(" << m_sink->getGraphID() << "," << m_sink->getID() << ")    "
				 << m_weight << "   " << m_BW << "   ";
		}

		for (vector<int>::const_iterator it_as = m_vASPath.begin(); it_as != m_vASPath.end(); ++it_as)
		{
			cout << *it_as << ",";
		}

        cout << "    ";
		for (set<int>::const_iterator it_set = m_isExchanged.begin(); it_set != m_isExchanged.end(); ++it_set)
		{
			cout << *it_set << ",";
		}
		cout << endl;
	}

};


#endif
