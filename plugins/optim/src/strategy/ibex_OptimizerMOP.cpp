//                                  I B E X                                   
// File        : ibex_Optimizer.cpp
// Author      : Gilles Chabert, Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : May 14, 2012
// Last Update : December 24, 2012
//============================================================================

#include "ibex_OptimizerMOP.h"
#include "ibex_Timer.h"
#include "ibex_Function.h"
#include "ibex_NoBisectableVariableException.h"
#include "ibex_Backtrackable.h"
#include "ibex_OptimData.h"

#include <float.h>
#include <stdlib.h>
#include <iomanip>

using namespace std;

namespace ibex {

const double OptimizerMOP::default_eps_x = 1e-5;

void OptimizerMOP::write_ext_box(const IntervalVector& box, IntervalVector& ext_box) {
	int i2=0;
	for (int i=0; i<n; i++,i2++) {
		if (i2>=n) i2++; // skip goal variables
		ext_box[i2]=box[i];
	}
}

void OptimizerMOP::read_ext_box(const IntervalVector& ext_box, IntervalVector& box) {
	int i2=0;
	for (int i=0; i<n; i++,i2++) {
		if (i2>=n) i2++; // skip goal variable
		box[i]=ext_box[i2];
	}
}

OptimizerMOP::OptimizerMOP(int n, Ctc& ctc, Bsc& bsc, LoupFinder& finder, CellBufferOptim& buffer,
		int goal_var, int goal_var2, double eps_x) : n(n),
                				ctc(ctc), bsc(bsc), UB_finder(finder), buffer(buffer),
                				eps_x(eps_x), trace(false), timeout(-1), status(SUCCESS),
                				time(0), nb_cells(0) {

	if (trace) cout.precision(12);
}

OptimizerMOP::~OptimizerMOP() {

}

//TODO: Implementar algoritmo sencillo para encontrar soluciones no dominadas y actualizar UB
bool OptimizerMOP::update_UB(const IntervalVector& box) {

	/*
	 1. Tomar punto medio de la caja mid(box)
	 2. Verificar si es factible usando las restricciones (ctrs)
	 3. Si es factible, evaluar el punto usando funciones objetivo (goal[0] y goal[1])
	 4. Insertar en mapa UB (si es no dominada) y actualizar eliminar soluciones dominadas de UB
	 5. Si el mapa UB fue modificado retornar true, si no false
	 */

}

//TODO: Insertar box en mapa UB (si es no dominada) y actualizar LB eliminando soluciones dominadas por box
void OptimizerMOP::update_LB_with_epsboxes(const IntervalVector& box) {
	double f1 = box[n].lb();
	double f2 = box[n+1].lb();


}



void OptimizerMOP::handle_cell(Cell& c, const IntervalVector& init_box ){

	contract_and_bound(c, init_box);

	if (c.box.is_empty()) {
		delete &c;
	} else {
		buffer.push(&c);

		nb_cells++;
	}
}

void OptimizerMOP::contract_and_bound(Cell& c, const IntervalVector& init_box) {

	//TODO: descartar c.box si es dominada por UB set

	//TODO: contract (y1,y2) with UB??



	/*================ contract x with f(x)=y1, f(x)=y2 and g(x)<=0 ================*/

	ctc.contract(c.box);

	if (c.box.is_empty()) return;

	/*========================= update loup =============================*/

	//IntervalVector tmp_box(n);
	//read_ext_box(c.box,tmp_box);

	bool loup_ch=update_UB(c.box);


	/*====================================================================*/
	if ( c.box.max_diam()<=eps_x ) {
		update_LB_with_epsboxes(c.box);
		c.box.set_empty();
		return;
	}

	// the current extended box in the cell is updated
	//write_ext_box(tmp_box,c.box);
}

OptimizerMOP::Status OptimizerMOP::optimize(const IntervalVector& init_box) {

	nb_cells=0;

	buffer.flush();
	LB.clear();
	UB.clear();

	//the box in cells have the n original variables plus the two objective variables (y1 and y2)
	Cell* root=new Cell(IntervalVector(n+2));

	write_ext_box(init_box,root->box);

	// add data required by the bisector
	bsc.add_backtrackable(*root);

	// add data required by the buffer
	buffer.add_backtrackable(*root);

	// add data required by optimizer + KKT contractor
//	root->add<EntailedCtr>();
//	//root->add<Multipliers>();
//	entailed=&root->get<EntailedCtr>();
//	entailed->init_root(user_sys,sys);

	time=0;
	Timer::start();
	handle_cell(*root,init_box);

	try {
		//TODO: Stopping criteria considering a precision UB-LB
		while (!buffer.empty()) {

		  if (trace >= 2) cout << buffer;

			Cell *c = buffer.top();

			try {
				pair<IntervalVector,IntervalVector> boxes=bsc.bisect(*c);

				pair<Cell*,Cell*> new_cells=c->bisect(boxes.first,boxes.second);

				buffer.pop();
				delete c; // deletes the cell.

				handle_cell(*new_cells.first, init_box);
				handle_cell(*new_cells.second, init_box);

				time_limit_check(); // TODO: not reentrant

			}
			catch (NoBisectableVariableException& ) {
				//boxes con size<eps son guardadas en LB
				cout << "Error: NoBisectableVariableException en OptimizerMOP" << endl;
			}
		}
	}
	catch (TimeOutException& ) {
		status = TIME_OUT;
		return status;
	}

	Timer::stop();
	time+= Timer::VIRTUAL_TIMELAPSE();

	if (LB.empty())
		status=INFEASIBLE;
	else
		status=SUCCESS;

	return status;
}

void OptimizerMOP::report(bool verbose) {

	if (!verbose) {
		cout << get_status() << endl;
		cout << "LB:" << endl;
		map< pair <double, double>, IntervalVector > :: iterator lb=LB.begin();

		for(;lb!=LB.end();lb++)
			cout << "(" << lb->first.first << "," << lb->first.second << ")" << endl;

		cout << "UB:" << endl;
		map< pair <double, double>, IntervalVector > :: iterator ub=UB.begin();

		for(;ub!=UB.end();ub++)
			cout << "(" << ub->first.first << "," << ub->first.second << ")" << endl;

		cout << endl << get_time() << " " << get_nb_cells() << endl;
		return;
	}

	switch(status) {
	case SUCCESS: cout << "\033[32m" << " optimization successful!" << endl;
	break;
	case INFEASIBLE: cout << "\033[31m" << " infeasible problem" << endl;
	break;
	case NO_FEASIBLE_FOUND: cout << "\033[31m" << " no feasible point found (the problem may be infesible)" << endl;
	break;
	case UNBOUNDED_OBJ: cout << "\033[31m" << " possibly unbounded objective (f*=-oo)" << endl;
	break;
	case TIME_OUT: cout << "\033[31m" << " time limit " << timeout << "s. reached " << endl;
	break;
	case UNREACHED_PREC: cout << "\033[31m" << " unreached precision" << endl;
	}

	cout << "\033[0m" << endl;

	cout << " cpu time used: " << time << "s." << endl;
	cout << " number of cells: " << nb_cells << endl;
}

void OptimizerMOP::time_limit_check () {
	if (timeout<=0) return;
	Timer::stop();
	time += Timer::VIRTUAL_TIMELAPSE();
	if (time >=timeout ) throw TimeOutException();
	Timer::start();
}

} // end namespace ibex
