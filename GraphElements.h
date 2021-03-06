///////////////////////////////////////////////////////////////////////////////
///  GraphElements.h
///  <TODO: insert file description here>
///
///  @remarks <TODO: insert remarks here>
///
///  @author Yan Qi @date 5/28/2010
///
///  $Id: GraphElements.h 65 2010-09-08 06:48:36Z yan.qi.asu $
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <iostream>
#include <vector>


template<class T>
class WeightGreater
{
public:
  // Determine priority.
  bool operator()(const T& a, const T& b) const
  {
    return a.Weight() > b.Weight();
  }

  bool operator()(const T* a, const T* b) const
  {
    return a->Weight() > b->Weight();
  }
};

template<class T>
class WeightLess
{
public:
  // Determine priority.
  bool operator()(const T& a, const T& b) const
  {
    return a.Weight() < b.Weight();
  }

  bool operator()(const T* a, const T* b) const
  {
    return a->Weight() < b->Weight();
  }
};

//////////////////////////////////////////////////////////////////////////
// A class for the object deletion
//////////////////////////////////////////////////////////////////////////
template<class T>
class DeleteFunc
{
public:
  void operator()(const T* it) const
  {
    delete it;
  }
};



/**************************************************************************
*  BaseVertex
*  <TODO: insert class description here>
*
*
*  @remarks <TODO: insert remarks here>
*
*  @author Yan Qi @date 6/6/2010
**************************************************************************/
class BaseVertex
{
  int m_nID;
  double m_dWeight;

  //hemin
  int m_graphID; // the id of the AS/graph containing this vertex

public:
  //---------hemin---------------
  BaseVertex(int nID, int graphID):m_nID(nID),m_graphID(graphID),m_dWeight(0)
  {

  }
  BaseVertex():m_nID(0),m_graphID(0),m_dWeight(0) {}
  bool operator==(const BaseVertex &other) const
  {
    if (this->m_nID != other.m_nID)
    {
      return false;
    }
    if (this->m_graphID != other.m_graphID)
    {
      return false;
    }
    return true;
  }

  //------------------------------

  int getID() const
  {
    return m_nID;
  }
  void setID(int ID_)
  {
    m_nID = ID_;
  }

  //--hemin---
  int getGraphID() const
  {
    return m_graphID;
  }
  void setGraphID(int graphID_)
  {
    m_graphID = graphID_;
  }

  double Weight() const
  {
    return m_dWeight;
  }
  void Weight(double val)
  {
    m_dWeight = val;
  }

  void PrintOut(std::ostream& out_stream)
  {
    out_stream << "(" << m_graphID << ", " << m_nID % 30 << ")";
  }
};


/**************************************************************************
*  BasePath
*  <TODO: insert class description here>
*
*
*  @remarks <TODO: insert remarks here>
*
*  @author Yan Qi @date 6/6/2010
**************************************************************************/
#pragma once

class BasePath
{
protected:

  int m_nLength;
  double m_dWeight;
  double m_BW;
  std::vector<BaseVertex*> m_vtVertexList;

public:
  BasePath(const std::vector<BaseVertex*>& vertex_list, double weight)
    :m_dWeight(weight),m_BW(0)
  {
    m_vtVertexList.assign(vertex_list.begin(), vertex_list.end());
    m_nLength = m_vtVertexList.size();
  }
  ~BasePath(void) {}

  bool operator==(const BasePath &other) const
  {

    if (this->m_nLength != other.m_nLength)
    {
      return false;
    }
    for (int i = 0; i < m_nLength; ++i)
    {
      if (*m_vtVertexList[i] == *(other.m_vtVertexList[i]))
      {
        continue;
      }
      else
      {
        return false;
      }
    }
    return true;
  }


  double Weight() const
  {
    return m_dWeight;
  }
  void Weight(double val)
  {
    m_dWeight = val;
  }
  double get_BW() const
  {
    return m_BW;
  }
  void set_BW(double bw)
  {
    m_BW = bw;
  }

  int length() const
  {
    return m_nLength;
  }

  BaseVertex* GetVertex(int i) const
  {
    return m_vtVertexList.at(i);
  }
  BaseVertex* GetLastVertex() const
  {
    return m_vtVertexList.at(m_nLength-1);
  }

  void SetVertex(int i, BaseVertex* newVertex)
  {
    m_vtVertexList[i] = newVertex;
  }

  bool SubPath(std::vector<BaseVertex*>& sub_path, BaseVertex* ending_vertex_pt)
  {

    for (std::vector<BaseVertex*>::const_iterator pos = m_vtVertexList.begin();
         pos != m_vtVertexList.end(); ++pos)
    {
      if (*pos != ending_vertex_pt)
      {
        sub_path.push_back(*pos);
      }
      else
      {
        //break;
        return true;
      }
    }

    return false;
  }

  // display the content
  void PrintOut(std::ostream& out_stream) const
  {
    out_stream << "Cost: " << m_dWeight << " Length: " << m_vtVertexList.size() << std::endl;
    for(std::vector<BaseVertex*>::const_iterator pos=m_vtVertexList.begin(); pos!=m_vtVertexList.end(); ++pos)
    {
      (*pos)->PrintOut(out_stream);
      if (pos != m_vtVertexList.end() - 1)
      {
        out_stream << "->";
      }
    }
    out_stream << std::endl <<  "*********************************************" << std::endl;
  }
};
