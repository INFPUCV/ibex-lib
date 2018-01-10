/*
 * ibex_CellMOP.h
 *
 *  Created on: 10 ene. 2018
 *      Author: iaraya
 */

#ifndef OPTIM_MOP_SRC_STRATEGY_IBEX_CELLMOP_H_
#define OPTIM_MOP_SRC_STRATEGY_IBEX_CELLMOP_H_

#include "ibex_CellBuffer.h"
#include "ibex_CellBufferOptim.h"

namespace ibex {

class CellMOP : public Backtrackable {
public:
	/**
	 * \brief Constructor for the root node (followed by a call to init_root).
	 */
	CellMOP() : depth(0), id(0), a(0.0), w_lb(POS_INFINITY), ub_distance(POS_INFINITY) {}

	/**
	 * \brief Copy constructor
	 */

	CellMOP(const CellMOP& c) : depth(c.depth+1), id(++nb_cells),
			a(c.a), w_lb(c.w_lb), ub_distance(c.ub_distance) { }

	/**
	 * \brief Duplicate the structure into the left/right nodes
	 */
	std::pair<Backtrackable*,Backtrackable*> down(){
		CellMOP* c1=new CellMOP(*this);
		CellMOP* c2=new CellMOP(*this);

		return std::pair<Backtrackable*,Backtrackable*>(c1,c2);
	}

	static int nb_cells;
	static Interval y1_init;
	static Interval y2_init;

    /**unique identifier for comparisons*/
    int id;

	/** depth of the node **/
	int depth;

	/** MOP: after filtering we know that z1+a*z2 > w_lb and we can
	 * use this information for filtering**/
	double a;
	double w_lb;


	/** MOP: distance of the box to the current non dominated set (UB) */
	double ub_distance;

};


} /* namespace ibex */

#endif /* OPTIM_MOP_SRC_STRATEGY_IBEX_CELLMOP_H_ */
