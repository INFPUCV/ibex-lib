//============================================================================
//                                  I B E X                                   
// File        : ibex_Cell.cpp
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : May 10, 2012
// Last Update : May 10, 2012
//============================================================================

#include "ibex_CellData.h"
#include "ibex_Backtrackable.h"

namespace ibex {

 CellData::CellData() {}

 		/**
 		 * \brief Copy constructor
 		 */

 CellData::CellData(const CellData& c): HC4(c.HC4), ACID(c.ACID), COMPO(c.COMPO) { }


	/**
	 * \brief Duplicate the structure into the left/right nodes
	 */
	std::pair<Backtrackable*,Backtrackable*> CellData::down(){
		CellData* c1= new CellData(*this);
		CellData* c2=new CellData(*this);

		return std::pair<Backtrackable*,Backtrackable*>(c1,c2);
	}

} // end namespace ibex
