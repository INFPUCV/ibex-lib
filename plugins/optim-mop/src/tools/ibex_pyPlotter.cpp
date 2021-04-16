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


void py_Plotter::offline_plot(set< Vector*, sort_rp >& NDS,
 set< Vector*, sort_rp >* NDS2, const char* output_file, IntervalVector* focus){
	ofstream output;
	output.open(output_file);

	output << "[";

	set< Vector* > :: iterator ub=NDS.begin();
	for(;ub!=NDS.end();ub++){
		if(!focus || (*focus).contains(**ub)){
			output << "(" << (**ub)[0] << "," << (**ub)[1] << "),";

      
      //output  << "(" << ub->first[0] << " ; " << ub->first[1] << ")_" <<
        //        ((ub->second.n==1)? ub->second.x1:0.0) << ",";
    }
	}

  output << "]" << endl;

  if(NDS2){
		output << "[";
		ub=NDS2->begin();
		for(;ub!=NDS2->end();ub++){
			if(!focus || (*focus).contains(**ub))
				output << "(" << (**ub)[0] << "," << (**ub)[1] << "),";
		}

	  output << "]" << endl;
  }else
		output << "[]" << endl;

	output.close();

}


void py_Plotter::offline_plot(list<vector<double> > &upperList, 
list<vector<double> > &lowerList, const char* output_file){
  ofstream output;
  output.open(output_file);

  output << "[";
  for(vector<double> v : upperList)
	output << "(" << v[0] << "," << v[1] << "),";
  output << "]" << endl;

   output << "[";
  for(vector<double> v : lowerList)
	output << "(" << v[0] << "," << v[1] << "),";
  output << "]" << endl;

  output.close();

}


} /* namespace ibex */
