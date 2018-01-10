/*
 * ibex_DistanceSortedCellBufferMOP.cpp
 *
 *  Created on: 20 oct. 2017
 *      Author: iaraya
 */

#include "ibex_DistanceSortedCellBufferMOP.h"
#include "ibex_OptimizerMOP.h"
#include <algorithm>    // std::min_element, std::max_element

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

	void DistanceSortedCellBufferMOP::push(Cell* cell) {
		double dist=OptimizerMOP::distance2(cell);
		if(dist < cell->get<CellMOP>().ub_distance )
			cell->get<CellMOP>().ub_distance=dist;
		cells.push(cell);
	}

	Cell* DistanceSortedCellBufferMOP::pop() {

        Cell* c = top();
        cells.pop();

		return c;
	}

	Cell* DistanceSortedCellBufferMOP::top() const {
		Cell* c = cells.top();
		if(!c) return NULL;


		double dist=OptimizerMOP::distance2(c);

		//we update the distance and reinsert the element
		while(dist < c->get<CellMOP>().ub_distance){
			cells.pop();
			c->get<CellMOP>().ub_distance=dist;
			cells.push(c);
			c = cells.top();
			dist=OptimizerMOP::distance2(c);
		}

		//cout << "dist:" << dist << endl;
		return c;
	}

} // end namespace ibex
