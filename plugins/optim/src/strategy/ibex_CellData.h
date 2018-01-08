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

#include "ibex_Random.h"
#include "ibex_Backtrackable.h"
#include "ibex_CellBuffer.h"
#include <map>
#include <queue>


using namespace std;

namespace ibex {

	class CellData : public Backtrackable {

public:
		/**
		 * \brief Constructor for the root node (followed by a call to init_root).
		 */
	CellData();

		/**
		 * \brief Copy constructor
		 */

	CellData(const CellData& c);

	std::pair<Backtrackable*,Backtrackable*> down();


		std::set<int> HC4;
		std::set<int> Acid;
		std::set<int> Compo;

};

} // end namespace ibex

#endif // __IBEX_CELL_DATA_H__
