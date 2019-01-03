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



	//map< pair <double, double>, IntervalVector >* max_distance::UB=NULL;


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
		double dist=nds->distance(cell);
		if(dist < cell->get<CellMOP>().ub_distance && !OptimizerMOP::_hv)
			cell->get<CellMOP>().ub_distance=dist;
		cells.push(cell);
	}

	Cell* DistanceSortedCellBufferMOP::pop() {

        Cell* c = top();
        cells.pop();

		return c;
	}

  int counter=0;
	Cell* DistanceSortedCellBufferMOP::top() const {

		Cell* c = cells.top();
		if(!c) return NULL;

		if (OptimizerMOP::_hv) return c;

		double dist=nds->distance(c);

		//we update the distance and reinsert the element
		while(dist < c->get<CellMOP>().ub_distance){
			cells.pop();
			c->get<CellMOP>().ub_distance=dist;
			cells.push(c);
			c = cells.top();
			dist=nds->distance(c);
		}

    counter ++;
		//cout << counter  <<":" <<  c->get<CellMOP>().ub_distance << endl;

		return c;
	}

} // end namespace ibex
