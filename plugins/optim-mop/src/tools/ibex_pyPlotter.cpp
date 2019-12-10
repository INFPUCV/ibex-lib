/*
 * ibex_pyPlotter.cpp
 *
 *  Created on: 10 ene. 2018
 *      Author: iaraya
 */

#include "ibex_pyPlotter.h"
#include "ibex_OptimizerMOP.h"
#include <unordered_map>
#include <iostream>

#ifndef cdata
#define cdata ((BxpMOPData*) c->prop[BxpMOPData::id])
#endif


#include "ibex_OptimizerMOP.h"

namespace ibex {


void py_Plotter::offline_plot(map< Vector, NDS_data, struct sorty2 >& NDS,
    map< Vector, NDS_data, struct sorty2 >* NDS2, const char* output_file, IntervalVector* focus){
	ofstream output;
	output.open(output_file);

	output << "[";

	map< Vector, NDS_data > :: iterator ub=NDS.begin();
	for(;ub!=NDS.end();ub++){
		if(!focus || (*focus).contains(ub->first))
			output << "(" << ub->first[0] << "," << ub->first[1] << "),";
	}

  output << "]" << endl;

  if(NDS2){
		output << "[";
		ub=NDS2->begin();
		for(;ub!=NDS2->end();ub++){
			if(!focus || (*focus).contains(ub->first))
				output << "(" << ub->first[0] << "," << ub->first[1] << "),";
		}

	  output << "]" << endl;
  }else
		output << "[]" << endl;

	output.close();

}

void py_Plotter::offline_plot(unordered_map< Vector, NDS_X, MyHashFunction > & NDS,
	unordered_map< Vector, NDS_X, MyHashFunction >* NDS2, const char* output_file, const char* output_file_solutions, IntervalVector* focus){
	ofstream output, output2;
	output.open(output_file);


	//Print pareto front
	unordered_map< Vector, NDS_X, MyHashFunction >:: iterator ub=NDS.begin();
	for(;ub!=NDS.end();ub++){
		if(!focus || (*focus).contains(ub->first)){
			for(auto i=0; i<ub->first.size(); i++){
				output <<std::setprecision(10)<< ub->first[i];
				if(i+1!=ub->first.size()) output<<" ";
			}
//			output << ub->first[0] << ";" << ub->first[1] <<  ";" <<ub->first[2];
			output << "\n";
		}else{
			output <<"entraaaaa";
		}
	}
	output.close();


	//Print Pareto optimal set
		output2.open(output_file_solutions);
		ub=NDS.begin();
		for(;ub!=NDS.end();ub++){
				if(!focus || (*focus).contains(ub->first)){
//					auto i = ub->second.x1->begin();
					Vector aux(*ub->second.x1);
					for(int i=0; i<aux.size(); i++){
						output2 << aux[i];
						if(i+1!=aux.size()) output2<<";";
					}
					output2 << "\n";
				}else{
					output2 <<"entraaaaa";
				}
		}
		output2.close();


}
} /* namespace ibex */
