/*
 * ibex_DistanceSortedCellBufferMOP.h
 *
 *  Created on: 20 oct. 2017
 *      Author: iaraya
 */

#ifndef OPTIM_SRC_STRATEGY_IBEX_DISTANCESORTEDCELLBUFFERMOP_H_
#define OPTIM_SRC_STRATEGY_IBEX_DISTANCESORTEDCELLBUFFERMOP_H_

#include "ibex_CellBuffer.h"
#include "ibex_CellBufferOptim.h"
#include "ibex_CellFeasibleDiving.h"
#include <queue>
#include <map>

using namespace std;

namespace ibex {


/** \ingroup strategy
 *
 * \brief Buffer which selects next the box maximizing the distance to the non dominated set.
 */
class DistanceSortedCellBufferMOP : public CellBufferOptim {
 public:

   virtual void add_backtrackable(Cell& root){
     root.add<CellBS>();
   }
   
  /** Flush the buffer.
   * All the remaining cells will be *deleted* */
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

  /**
	* \brief Return the minimum value of the heap
	*
	*/
  virtual double minimum() const {
	  cout << "DistanceSortedCellBufferMOP::minimum is not implemented" << endl;
	  exit(0);
	  return 0.0;
  }

	/**
	 * \brief Contract the buffer using the UB
	 */
	virtual void contract(double loup){
		  cout << "DistanceSortedCellBufferMOP::contract is not implemented" << endl;
		  exit(0);
	}


	mutable std::priority_queue<Cell*, std::vector<Cell*>, max_distance > cells;


};





} // end namespace ibex
#endif //  /* OPTIM_SRC_STRATEGY_IBEX_DISTANCESORTEDCELLBUFFERMOP_H_ */
