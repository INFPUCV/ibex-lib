//============================================================================
//                                  I B E X
// File        : ibex_CellSet.h
// Author      : Gilles Chabert, Jordan Ninin
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Sep 12, 2014
//============================================================================

#ifndef __IBEX_CELL_SET_H__
#define __IBEX_CELL_SET_H__

#include "ibex_Random.h"
#include "ibex_CellBuffer.h"
#include <set>
// #include "../strategy/ibex_Cell.h"

namespace ibex {


template<class T>
class CellSet : public CellBuffer {
public:

	CellSet();

	void flush();

  /** Return the size of the buffer. */
  unsigned int size() const;

  /** Return true if the buffer is empty. */
  bool empty() const;

  /** push a new cell on the stack. */
  void push(Cell* cell);

  /** Pop a cell from the stack and return it.*/
  Cell* pop();

  /** Return the next box (but does not pop it).*/
  Cell* top() const;

private:
	/* Set of Cells */
	std::set<Cell*, T> cset;

};

struct minLB {
  bool operator() (const Cell* c1, const Cell* c2) const
  {
	  if(c1->get<CellBS>().lb != c2->get<CellBS>().lb) return (c1->get<CellBS>().lb < c2->get<CellBS>().lb);
	  if(c1->get<CellBS>().depth != c2->get<CellBS>().depth) return (c1->get<CellBS>().depth < c2->get<CellBS>().depth);
	  return (c1->get<CellBS>().id > c2->get<CellBS>().id);
  }
};

template class CellSet<minLB>;

#endif // __IBEX_CELL_SET_H__
