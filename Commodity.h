/*
 * Commodity.h
 *
 *  Created on: Mar 27, 2016
 *      Author: Hemin
 */
#pragma once

#ifndef COMMODITY_H
#define COMMODITY_H
#include "GraphElements.h"
const double priority = 1 ; // the priority of the commodity
class Commodity
{
public:
  BaseVertex* source_; // the source of this commodity
  BaseVertex* sink_;   // the sink of this commodity
  double demand_;         // the demand of this commodity
  double Allocated_;      // the allocated bandwidth for this commodity
  double priority_;    // the priority of this commodity, i.e., the slope of the bandwidth function, see Google B4
  bool isSaturated_;   // indicate whether the commodity is saturated

  Commodity(BaseVertex* source,BaseVertex* sink, double demand):
    source_(source),
    sink_(sink),
    priority_(priority),
    demand_(demand),
    Allocated_(0),
    isSaturated_(false)
  {
  }
  void Print(std::ostream& out_stream)
  {
    out_stream << "(" << source_->getGraphID() << "," << source_->getID() << ") -->"
               << "(" << sink_->getGraphID() << "," << sink_->getID() << "):"
               << "demand = " << demand_ << "; Allocated = " << Allocated_ << "; priority = " << priority_ << std::endl;
  }
};

#endif



