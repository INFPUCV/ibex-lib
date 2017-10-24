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



	void DistanceSortedCellBufferMOP::flush() {
		while (!cells.empty()) {
			delete cells.top();
			cells.pop();
		}
	}

	unsigned int DistanceSortedCellBufferMOP::size() const {
		return cells.size();
	}

	bool DistanceSortedCellBufferMOP::empty() const {
		return cells.empty();
	}

	void DistanceSortedCellBufferMOP::push(Cell* cell) {
		cell->get<CellBS>().ub_distance=max_distance::distance(cell->box);
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
		while(dist!=c->get<CellBS>().ub_distance){
			cells.pop();
			c->get<CellBS>().ub_distance=dist;
			cells.push(c);
			c = cells.top();
			dist=OptimizerMOP::distance2(c);
		}

		cout << "dist:" << dist << endl;
		return c;
	}

} // end namespace ibex

