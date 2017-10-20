/*
 * ibex_DistanceSortedCellBufferMOP.cpp
 *
 *  Created on: 20 oct. 2017
 *      Author: iaraya
 */

#include "ibex_DistanceSortedCellBufferMOP.h"
#include <algorithm>    // std::min_element, std::max_element

namespace ibex {



	void DistanceSortedCellBufferMOP::flush() {
		while (!cells.empty()) {
			delete cells.front();
			cells.pop_front();
		}
	}

	unsigned int DistanceSortedCellBufferMOP::size() const {
		return cells.size();
	}

	bool DistanceSortedCellBufferMOP::empty() const {
		return cells.empty();
	}

	void DistanceSortedCellBufferMOP::push(Cell* cell) {
		cells.push_back(cell);
	}

	Cell* DistanceSortedCellBufferMOP::pop() {
		std::list<Cell*>::iterator it_max = std::min_element(cells.begin(),cells.end(),max_distance());
		Cell* c = *it_max;
		cells.erase(it_max);
		return c;
	}

	Cell* DistanceSortedCellBufferMOP::top() const {
		std::list<Cell*>::const_iterator it_max = std::min_element(cells.begin(),cells.end(),max_distance());
		return *it_max;
	}

	} // end namespace ibex

