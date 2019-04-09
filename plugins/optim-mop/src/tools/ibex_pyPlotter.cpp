/*
 * ibex_pyPlotter.cpp
 *
 *  Created on: 10 ene. 2018
 *      Author: iaraya
 */

#include "ibex_pyPlotter.h"
#include "ibex_OptimizerMOP.h"
#include <iostream>


#include "ibex_OptimizerMOP.h"

namespace ibex {

int py_Plotter::n=0;

void py_Plotter::plot_add_ub(pair<double, double> eval){
	std::cout << "add ub: {\"pts\": (" << eval.first << "," << eval.second << ")}" << endl;
}

void py_Plotter::plot_del_ub(pair<double, double> eval){
	std::cout << "del ub: {\"pts\": (" << eval.first << "," << eval.second << ")}" << endl;
}

void py_Plotter::plot_add_lb(Cell* c){
	std::cout << "add lb: {\"id\":" << c->get<CellMOP>().id;
	std::cout << ", 'pts':(" << c->box[n].lb() << "," <<  c->box[n+1].lb() << ")";
	std::cout << "}" << endl;

}

void py_Plotter::plot_add_box(Cell* c){
	std::cout << "add: {\"id\":" << c->get<CellMOP>().id;
	std::cout << ", 'pts':(" << c->box[n].lb() << "," <<  c->box[n+1].lb() << "),";
	std::cout << "'diam_x': " <<  c->box[n].diam() << ",'diam_y': " << c->box[n+1].diam();
	std::cout << ", 'pA':(" << c->box[n].lb() <<"," <<  (((c)->get<CellMOP>().w_lb-c->box[n].lb())/(c)->get<CellMOP>().a)   << "),";
	std::cout << "'pB':(" << (c->get<CellMOP>().w_lb-c->get<CellMOP>().a*c->box[n+1].lb()) <<"," <<  c->box[n+1].lb()  << ")";
	std::cout << "}" << endl;
}

void py_Plotter::plot_del_box(Cell* c){
	std::cout << "del: {\"id\":" << c->get<CellMOP>().id;
	std::cout << "}" << endl;
}
void py_Plotter::offline_plot(Cell* c, map< pair <double, double>, IntervalVector >& NDS){
	cout << "print plot:" << NDS.size() << endl;
	ofstream output;
	output.open("output.txt");

	output << "[";

	map< pair <double, double>, IntervalVector > :: iterator ub=NDS.begin();
	for(;ub!=NDS.end();ub++){
		output << "(" << ub->first.first << "," << ub->first.second << "),";
	}
output << "]" << endl;

	output.close();

}

void py_Plotter::offline_plot(Cell* c, map< pair <double, double>, IntervalVector, struct sorty2 >& NDS,
 map< pair <double, double>, IntervalVector, struct sorty2 >* NDS2){
	//cout << "print plotX2:" << NDS.size() << endl;
	ofstream output;
	output.open("output2.txt");

	output << "[";

	map< pair <double, double>, IntervalVector > :: iterator ub=NDS.begin();
	for(;ub!=NDS.end();ub++){
		//cout << "(" << ub->first.first << "," << ub->first.second << ")" << endl;
		output << "(" << ub->first.first << "," << ub->first.second << "),";
	}
    output << "]" << endl;

  if(NDS2){
		output << "[";
		ub=NDS2->begin();
		for(;ub!=NDS2->end();ub++){
			//cout << "(" << ub->first.first << "," << ub->first.second << ")" << endl;
			output << "(" << ub->first.first << "," << ub->first.second << "),";
		}
	  output << "]" << endl;
  }else{
		output << "[]" << endl;
	}
	output.close();

}

void py_Plotter::offline_plot(Cell* c, map< pair <double, double>, IntervalVector, struct sorty2 >& NDS,
		std::vector< pair <double, double> > rectaUB,
		std::vector< pair <double, double> > functionPoly
		){
	cout << "print plotX1:" << NDS.size() << endl;
	ofstream output;
	output.open("output2.txt");
	output << "[";

	map< pair <double, double>, IntervalVector > :: iterator ub=NDS.begin();
	for(;ub!=NDS.end();ub++){
		output << "(" << ub->first.first << "," << ub->first.second << "),";
	}
output << "]" << endl;

		output << "[";
		for (int i=0;i<rectaUB.size();i++) {
			output << "(" << rectaUB[i].first << "," << rectaUB[i].second << "),";
		}
		output << "]" << endl;

		output << "[";
		for (int i=0;i<functionPoly.size();i++) {
			output << "(" << functionPoly[i].first << "," << functionPoly[i].second << "),";
		}
		output << "]" << endl;

	output.close();
	// system("python3 plot.py");
}

/**
 * Hamburger plot
 *
 */
void py_Plotter::offline_plot(Cell* c,
		map< pair <double, double>, IntervalVector, struct sorty2 > NDS,
		std::vector< pair <double, double> > rectaUB,
		std::vector< pair <double, double> > functionPoly_origin,
		std::vector< pair <double, double> > functionPoly
		){

	cout << "print plot Hamburger:" << NDS.size() << endl;
	ofstream output;
	output.open("outputH.txt");
	output << "[";
	map< pair <double, double>, IntervalVector, sorty2 >::iterator ub=NDS.begin();
	// NDS_seg::iterator ub=NDS.begin();
	for(;ub!=NDS.end();ub++){
		output << "(" << ub->first.first << "," << ub->first.second << "),";
	}
output << "]" << endl;


		output << "[";
		for (int i=0;i<rectaUB.size();i++) {
			output << "(" << rectaUB[i].first << "," << rectaUB[i].second << "),";
		}
		output << "]" << endl;

		output << "[";
		for (int i=0;i<functionPoly_origin.size();i++) {
			output << "(" << functionPoly_origin[i].first << "," << functionPoly_origin[i].second << "),";
		}
		output << "]" << endl;

		output << "[";
		for (int i=0;i<functionPoly.size();i++) {
			output << "(" << functionPoly[i].first << "," << functionPoly[i].second << "),";
		}
		output << "]" << endl;

	output.close();
	// system("python3 plot.py");
}

} /* namespace ibex */
