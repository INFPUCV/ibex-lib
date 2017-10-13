//============================================================================
//                                  I B E X
// File        : ibex_CellNSSet.h
// Author      : Matias Campusano, Damir Aliquintui, Ignacio Araya
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Oct 12, 2017
// Last Update : Oct 12, 2017
//============================================================================

#ifndef __IBEX_CELL_NS_SET_H__
#define __IBEX_CELL_NS_SET_H__

#include "ibex_CellBufferOptim.h"
#include <set>
#include <list>

namespace ibex {

	struct maxsize {
	  bool operator() (const Cell* c1, const Cell* c2) const
	  {
		  int n = c1->box.size();

		  float areaC1 = c1->box[n-1].diam()*c1->box[n-2].diam();
		  float areaC2 = c2->box[n-1].diam()*c2->box[n-2].diam();
		  return (areaC1 > areaC2);
	  }
	};

	class CellNSSet: public CellBufferOptim {
		public:
			CellNSSet();

			virtual ~CellNSSet();

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


			  virtual double minimum() const;

			  virtual void contract(double loup);

		private:
			/* Set of Cells */
			  typename std::set<Cell*, maxsize> nondset;
			/* Stack of cells */
			  typename std::list<Cell*> dset;
	};


} // namespace ibex

#endif /* __IBEX_CELL_NS_SET_H__ */
