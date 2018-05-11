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
	//set<  Cell* > :: iterator cell=buffer_cells.begin();
/*
	output << "(";
	if(c){
		output << "{'pts':(" << c->box[n].lb() << "," <<  c->box[n+1].lb() << "),";
		output << "'diam_x': " <<  c->box[n].diam() << ",'diam_y': " << c->box[n+1].diam()<< ",";
		output << "'pA':(" << c->box[n].lb() <<"," <<  (((c)->get<CellMOP>().w_lb-c->box[n].lb())/(c)->get<CellMOP>().a)   << "),";
		output << "'pB':(" << (c->get<CellMOP>().w_lb-c->get<CellMOP>().a*c->box[n+1].lb()) <<"," <<  c->box[n+1].lb()  << ")";
		output << "},";
  }

	for(;cell!=buffer_cells.end();cell++){
		//if(distance2(*cell) < 0){continue;}

		output << "{'pts':(" << (*cell)->box[n].lb() << "," <<  (*cell)->box[n+1].lb() << "),";
		output << "'diam_x': " <<  (*cell)->box[n].diam() << ",'diam_y': " <<  (*cell)->box[n+1].diam() << ",";
		output << "'pA':(" << (*cell)->box[n].lb() <<"," <<  (((*cell)->get<CellMOP>().w_lb-(*cell)->box[n].lb())/(*cell)->get<CellMOP>().a)   << "),";
		output << "'pB':(" << ((*cell)->get<CellMOP>().w_lb-(*cell)->get<CellMOP>().a*(*cell)->box[n+1].lb()) <<"," <<  (*cell)->box[n+1].lb()  << ")";
		output << "},";
	}
	output << ")" << endl;
*/
	output << "[";

	map< pair <double, double>, IntervalVector > :: iterator ub=NDS.begin();
	for(;ub!=NDS.end();ub++){
		output << "(" << ub->first.first << "," << ub->first.second << "),";
	}
output << "]" << endl;

	output.close();
	// system("python3 plot.py");
}

void py_Plotter::offline_plot(Cell* c, map< pair <double, double>, IntervalVector, struct sorty2 >& NDS){
	cout << "print plot:" << NDS.size() << endl;
	ofstream output;
	output.open("output.txt");
	//set<  Cell* > :: iterator cell=buffer_cells.begin();
/*
	output << "(";
	if(c){
		output << "{'pts':(" << c->box[n].lb() << "," <<  c->box[n+1].lb() << "),";
		output << "'diam_x': " <<  c->box[n].diam() << ",'diam_y': " << c->box[n+1].diam()<< ",";
		output << "'pA':(" << c->box[n].lb() <<"," <<  (((c)->get<CellMOP>().w_lb-c->box[n].lb())/(c)->get<CellMOP>().a)   << "),";
		output << "'pB':(" << (c->get<CellMOP>().w_lb-c->get<CellMOP>().a*c->box[n+1].lb()) <<"," <<  c->box[n+1].lb()  << ")";
		output << "},";
  }

	for(;cell!=buffer_cells.end();cell++){
		//if(distance2(*cell) < 0){continue;}

		output << "{'pts':(" << (*cell)->box[n].lb() << "," <<  (*cell)->box[n+1].lb() << "),";
		output << "'diam_x': " <<  (*cell)->box[n].diam() << ",'diam_y': " <<  (*cell)->box[n+1].diam() << ",";
		output << "'pA':(" << (*cell)->box[n].lb() <<"," <<  (((*cell)->get<CellMOP>().w_lb-(*cell)->box[n].lb())/(*cell)->get<CellMOP>().a)   << "),";
		output << "'pB':(" << ((*cell)->get<CellMOP>().w_lb-(*cell)->get<CellMOP>().a*(*cell)->box[n+1].lb()) <<"," <<  (*cell)->box[n+1].lb()  << ")";
		output << "},";
	}
	output << ")" << endl;
*/
	output << "[";

	map< pair <double, double>, IntervalVector > :: iterator ub=NDS.begin();
	for(;ub!=NDS.end();ub++){
		output << "(" << ub->first.first << "," << ub->first.second << "),";
	}
output << "]" << endl;

	output.close();
	// system("python3 plot.py");
}

void py_Plotter::offline_plot(Cell* c, map< pair <double, double>, IntervalVector, struct sorty2 >& NDS,
		std::vector< pair <double, double> > rectaUB,
		std::vector< pair <double, double> > functionPoly
		){
	cout << "print plot:" << NDS.size() << endl;
	ofstream output;
	output.open("output.txt");
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

} /* namespace ibex */
