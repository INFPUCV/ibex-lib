//============================================================================
//                                  I B E X                                   
// File        : CellDatas
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : May 10, 2012
// Last Update : May 10, 2012
//============================================================================

#ifndef __IBEX_CELL_DATA_H__
#define __IBEX_CELL_DATA_H__

#include "ibex_Cell.h"
#include <set>

namespace ibex {

/**
 * \defgroup strategy Strategies
 */

/** \ingroup strategy
 *
 * \brief Representation of the search space.
 *
 * This representation includes default data (current box) and data related to
 * user-defined contractors or bisectors. A different cell is associated to each
 * node and cell construction/inheritance can be controlled (see #ibex::Backtrackable).
 *
 * The cell on its own contains the minimum of information associated to the actual search space.
 * Besides the current box (the search space), this minimum information includes, e.g., the number
 * of the last bisected variable (other fields might be added with future releases).
 *
 * The amount of information contained in a cell can be arbitrarily augmented thanks to the
 * "data registration" technique (see #ibex::Contractor::require()).
 */
class CellData: public Cell {
public:

	/**
	 * \brief Create the root cell.
	 *
	 * \param box - Box (passed by copy).
	 */
	CellData(const IntervalVector& box);




	std::set<int> HC4;
	std::set<int> Acid;
	std::set<int> Compo;

	/**
	 * \set<int>::iterator it;
	myset.insert(0);
	myset.insert(1);

	for (it=myset.begin(); it!=myset.end(); ++it)
	    cout << ' ' << *it;
	 */

private:
	/* A constant to be used when no variable has been split yet (root cell). */
	//static const int ROOT_CELL;
};

} // end namespace ibex

#endif // __IBEX_CELL_DATA_H__
