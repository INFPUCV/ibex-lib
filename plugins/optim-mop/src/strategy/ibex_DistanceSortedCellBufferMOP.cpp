/*
 * ibex_DistanceSortedCellBufferMOP.cpp
 *
 *  Created on: 20 oct. 2017
 *      Author: iaraya
 */

#include "ibex_DistanceSortedCellBufferMOP.h"
#include "ibex_OptimizerMOP.h"
#include <algorithm>    // std::min_element, std::max_element

#ifndef cdata
#define cdata ((BxpMOPData*) c->prop[BxpMOPData::id])
#endif


namespace ibex {



	map< pair <double, double>, IntervalVector >* max_distance::UB=NULL;


	void DistanceSortedCellBufferMOP::flush() {
		while (!cells.empty()) {
			delete pop();
		}
	}

	unsigned int DistanceSortedCellBufferMOP::size() const {
		return cells.size();
	}

	bool DistanceSortedCellBufferMOP::empty() const {
		return cells.empty();
	}

	void DistanceSortedCellBufferMOP::push(Cell* c) {
		double dist = nds->distance(c->box);
		if(dist < cdata->ub_distance){
			cdata->ub_distance=dist;

		}
		cells.push(c);
	}

	Cell* DistanceSortedCellBufferMOP::pop() {
        Cell* c = top();
        cells.pop();
		return c;
	}


	Cell* DistanceSortedCellBufferMOP::top() const {
		Cell* c = cells.top();
		if(!c) return NULL;
		//Obtain dist from cell top
		double dist=nds->distance(c->box);
		//we update the distance and reinsert the element
		while(dist < cdata->ub_distance){
			cells.pop();
			cdata->ub_distance=dist;
			cells.push(c);
			c = cells.top();
			dist=nds->distance(c->box);
		}
		return c;
	}

} // end namespace ibex
