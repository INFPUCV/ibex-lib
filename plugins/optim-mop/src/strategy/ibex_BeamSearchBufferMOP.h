/*
 * ibex_BeamSearchBufferMOP.h
 *
 *  Created on: 20 oct. 2017
 *      Author: matias y pablo
 */

#ifndef OPTIM_SRC_STRATEGY_IBEX_BEAMSEARCHBUFFERMOP_H_
#define OPTIM_SRC_STRATEGY_IBEX_BEAMSEARCHBUFFERMOP_H_



#include "ibex_CellMOP.h"
#include "ibex_CellSet.h"
#include "ibex_NDS.h"
#include <queue>
#include <map>

using namespace std;

namespace ibex {

/**
 * Criteria for bi-objective problems
 */
struct max_distanceBeam {


	bool operator() (const Cell* c1, const Cell* c2){
	   int n = c1->box.size();
	   if(c1->get<CellMOP>().ub_distance != c2->get<CellMOP>().ub_distance)
		   return (c1->get<CellMOP>().ub_distance < c2->get<CellMOP>().ub_distance);
	   else if(c1->box[n-2].lb() >= c2->box[n-2].lb() && c1->box[n-1].lb() >= c2->box[n-1].lb()) return true;
	   else return false;
	}

};

struct min_distanceBeam {


	bool operator() (const Cell* c1, const Cell* c2){
	   int n = c1->box.size();
	   if(c1->get<CellMOP>().ub_distance != c2->get<CellMOP>().ub_distance)
		   return (c1->get<CellMOP>().ub_distance > c2->get<CellMOP>().ub_distance);
	   else if(c1->box[n-2].lb() <= c2->box[n-2].lb() && c1->box[n-1].lb() <= c2->box[n-1].lb()) return true;
	   else return false;
	}

};

/** \ingroup strategy
 *
 * \brief Buffer which selects next the box maximizing the distance to the non dominated set.
 */
class BeamSearchBufferMOP : public CellBufferOptim {
 public:

   void set(NDS_seg& nds) {
		 this->nds=&nds;
	 }

   virtual void add_backtrackable(Cell& root){
     root.add<CellMOP>();
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
	  cout << "BeamSearchBufferMOP::minimum is not implemented" << endl;
	  exit(0);
	  return 0.0;
  }

	/**
	 * \brief Contract the buffer using the UB
	 */
	virtual void contract(double loup){
		  cout << "BeamSearchBufferMOP::contract is not implemented" << endl;
		  exit(0);
	}


	/**
	 * A heap data structure for keeping the cells sorted by distance
	 */

	mutable std::priority_queue<Cell*, std::vector<Cell*>, max_distanceBeam > globalBuffer;
    mutable std::priority_queue<Cell*, std::vector<Cell*>, max_distanceBeam > currentBuffer;
    mutable std::multiset <Cell*, min_distanceBeam> nextBuffer;

  NDS_seg* nds;
  private:
	int cont = 0;
};





} // end namespace ibex
#endif //  /* OPTIM_SRC_STRATEGY_IBEX_BeamSearchBufferMOP_H_ */
