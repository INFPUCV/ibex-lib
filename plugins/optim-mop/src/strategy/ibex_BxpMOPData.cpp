//============================================================================
//                                  I B E X
// File        : ibex_BxpMOPData.cpp
// Author      : Ignacio Araya
// License     : See the LICENSE file
// Created     : Oct 18, 2014
// Last Update : Jul 05, 2018
//============================================================================

#include "ibex_BxpMOPData.h"
#include "ibex_OptimizerMOP.h"
#include "ibex_Id.h"


namespace ibex {

int BxpMOPData::nb_cells = 0;
const long BxpMOPData::id = next_id();

Map<long,false>& BxpMOPData::ids() {
	static Map<long,false> _ids;
	return _ids;
}

BxpMOPData::BxpMOPData() : Bxp(id), idd(0), a(0.0), w_lb(POS_INFINITY), ub_distance(POS_INFINITY){
	yn_init.resize(OptimizerMOP::nb_ObjFunc);
	for(int i=0; i<OptimizerMOP::nb_ObjFunc; i++) yn_init.set_empty();
}

BxpMOPData::BxpMOPData(const BxpMOPData& e) : Bxp(e.id), idd(e.idd), a(e.a), w_lb(e.w_lb), ub_distance(e.ub_distance) {

}

BxpMOPData::~BxpMOPData() {

}
} // end namespace ibex
