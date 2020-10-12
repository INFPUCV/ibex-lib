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
						if(i+1!=aux.size()) output2<<" ";
					}
					output2 << "\n";
				}else{
					output2 <<"entraaaaa";
				}
		}
		output2.close();


}

void py_Plotter::timer_save(list<pair<double, pair<double, double>> > t_checkDominance, list<pair<double, pair<double, double>> > t_distance, list<pair<double, pair<double, double>> > t_upperbounding){
	ofstream output_dominance, output_distance, o_upperbounding;

	//timer check dominance inside
	std::string nombre = "tCheckDominance.txt";
	output_dominance.open(nombre.c_str());
	output_dominance<<"t_exec t_oper n_sol\n";
	for(auto tiempo = t_checkDominance.begin(); tiempo != t_checkDominance.end(); ++tiempo){
		output_dominance<<tiempo->first<<" "<<tiempo->second.first<<" "<< tiempo ->second.second;
		output_dominance << "\n";
	}

	output_dominance.close();

	nombre = "tDistance.txt";
	output_distance.open(nombre.c_str());
	output_distance <<"t_exec t_oper n_nodes\n";
	for(auto data = t_distance.begin(); data != t_distance.end(); ++data){
		output_distance<<data->first<<" "<< data->second.first<<" "<< data ->second.second;
		output_distance << "\n";
	}

	output_distance.close();


	nombre = "tUpperBounding.txt";
	o_upperbounding.open(nombre.c_str());
	o_upperbounding<<"t_exec t_oper n_sol\n";
	for(auto tiempo = t_upperbounding.begin(); tiempo != t_upperbounding.end(); ++tiempo){
		o_upperbounding<<tiempo->first<<" "<<tiempo->second.first<<" "<< tiempo ->second.second;
		o_upperbounding << "\n";
	}

	o_upperbounding.close();
}


void py_Plotter::upper_envelope_save(list<pair<IntervalVector,Vector>> upper_envelope){
	ofstream output;
	std::string nombre = "upper_envelope.txt";
	output.open(nombre.c_str());
	output<<"box"<<endl;
	for(auto it=upper_envelope.begin(); it != upper_envelope.end(); ++it){
		//output<<it->first<<endl;
		for(auto it2 = it->first.begin(); it2 != it->first.end(); ++it2){
			Interval itv(*it2);
			output << itv.lb()<<" "<<itv.ub();
			if(it->first.size()+1 != it->first.size()) output<<" ";
		}
		output << "\n";
	}

	output<<"hyperplane"<<endl;
	for(auto it=upper_envelope.begin(); it != upper_envelope.end(); ++it){
		//output<<it->second<<endl;
		for(auto it2 = it->second.begin(); it2 != it->second.end(); ++it2){
			output << *it2;
			if(it->first.size()+1 != it->first.size()) output<<" ";
		}
		output << "\n";
	}

	output.close();

}


void py_Plotter::timer_save_component(const char* nombre, list<pair<double, pair<double, double>> > data){
	ofstream output_dominance;

	//timer check dominance inside
	//std::string nombre = "tCheckDominance.txt";
	output_dominance.open(nombre);
	output_dominance<<"t_exec t_oper n_sol/nb_nodes\n";
	for(auto tiempo = data.begin(); tiempo != data.end(); ++tiempo){
		output_dominance<<tiempo->first<<" "<<tiempo->second.first<<" "<< tiempo ->second.second;
		output_dominance << "\n";
	}

	output_dominance.close();
}





} /* namespace ibex */
