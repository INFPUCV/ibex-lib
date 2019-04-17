/*
 * ibex_pyPlotter.cpp
 *
 *  Created on: 10 ene. 2018
 *      Author: iaraya
 */

#include "ibex_pyPlotter.h"
#include "ibex_OptimizerMOP.h"
#include <iostream>

#ifndef cdata
#define cdata ((BxpMOPData*) c->prop[BxpMOPData::id])
#endif


#include "ibex_OptimizerMOP.h"

namespace ibex {


void py_Plotter::offline_plot(map< pair <double, double>, IntervalVector, struct sorty2 >& NDS,
 map< pair <double, double>, IntervalVector, struct sorty2 >* NDS2, const char* output_file){
	ofstream output;
	output.open(output_file);

	output << "[";

	map< pair <double, double>, IntervalVector > :: iterator ub=NDS.begin();
	for(;ub!=NDS.end();ub++)
		output << "(" << ub->first.first << "," << ub->first.second << "),";

  output << "]" << endl;

  if(NDS2){
		output << "[";
		ub=NDS2->begin();
		for(;ub!=NDS2->end();ub++)
			output << "(" << ub->first.first << "," << ub->first.second << "),";

	  output << "]" << endl;
  }else
		output << "[]" << endl;

	output.close();

}

} /* namespace ibex */
