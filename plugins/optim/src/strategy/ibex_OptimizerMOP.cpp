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

const double OptimizerMOP::default_eps_z = 1e-5;

OptimizerMOP::OptimizerMOP(int n, const Array<NumConstraint>& ctrs, const Function &f1,  const Function &f2,
		Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder, double eps_x, double eps_z) : n(n),
                				ctc(ctc), bsc(bsc), buffer(buffer), ctrs(ctrs), goal1(f1), goal2(f2),
								finder(finder),
                				eps_x(eps_x), eps_z(eps_z), trace(false), timeout(-1), status(SUCCESS),
                				time(0), nb_cells(0) {

	if (trace) cout.precision(12);
}

OptimizerMOP::~OptimizerMOP() {

}

Interval OptimizerMOP::eval_goal(const Function& goal, Vector& x){
	//the objectives are set to 0.0
	x[n]=0.0;
	x[n+1]=0.0;
	return goal.eval(x);
}


bool OptimizerMOP::update_UB(const IntervalVector& box) {

	list<Vector> feasible_points;

	//We attempt to find two feasible points which minimize both objectives
	//and the middle point between them
	IntervalVector box2(box); box2.resize(n);
	finder.find(box2,feasible_points);

	bool new_ub=false;
	list<Vector>::iterator it=feasible_points.begin();

	for(;it!=feasible_points.end();it++){
		Vector vec=*it;
		vec.resize(n+2);

		//3. Se evalua el punto usando funciones objetivo (goal1 y goal2)
		pair< double, double> eval = make_pair(eval_goal(goal1,vec).ub(), eval_goal(goal2,vec).ub());

		//4. Insertar en mapa UB (si es no dominada) y actualizar eliminar soluciones dominadas de UB
		bool insert_ub=true;
		for(std::map<pair<double, double>, Vector>::iterator it=UB.begin(); it!=UB.end(); ++it ){
			if (eval.first >= it->first.first && eval.second >= it->first.second){
				insert_ub=false;
				break;
			}
		}

		if(insert_ub){
			for(std::map<pair<double, double>, Vector>::iterator it=UB.begin(); it!=UB.end(); ){
				if (eval.first <= it->first.first && eval.second <= it->first.second){
					std::map<pair<double, double>, Vector>::iterator aux = it;
					++aux;
					UB.erase(it);
					it = aux;
				}
				else {
					++it;
				}
			}

			UB.insert(make_pair(eval, vec));
			new_ub = true;
		}
		//cout << UB.size() << endl;
		//5. Si el mapa UB fue modificado retornar true, si no false

	}
	return new_ub;

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
	double z1, z2;
	map< pair <double, double>, Vector >:: iterator ent1;
	for(ent1 = UB.begin(); ent1!= UB.end() ; ent1++) {
		z1 = ent1->first.first; // pair 1
		z2 = ent1->first.second; // pair 2
		c.box[n].lb(); // valor minimo
		c.box[n].ub(); // valor maximo
		// se elimina c si un UB es dominante de c
		if(z1 < c.box[n].lb() && z2 < c.box[n+1].lb()) {
			c.box.set_empty();
			return;
		}
	}


	//TODO (DA, MC): contract c.box[n] && c.box[n+1] with UB points



	/*================ contract x with f(x)=y1, f(x)=y2 and g(x)<=0 ================*/

	ctc.contract(c.box);

	if (c.box.is_empty()) return;

	/*========================= update loup =============================*/



	bool loup_ch=update_UB(c.box);


	/*====================================================================*/
	double diamX = 0.0, valueBox;
	int i;
	for (i=0; i < n; i++) {
		valueBox = c.box[i].diam();
		if(diamX < valueBox) {
			diamX = valueBox;
		}
	}
	// Metodo de termino para las restricciones
	if ( diamX<=eps_x ) {
		//se guarda c.box en lista de soluciones (Sout)
		Sout.push_back(c.box);
		c.box.set_empty();
		return;
	}

	double diamZ = 0;
	for (i=n; i < n+2; i++) {
		valueBox = c.box[i].diam();
		if(diamZ < valueBox) {
			diamZ = valueBox;
		}
	}
	// Metodo de termino para las funciones objetivo
	if ( diamZ<=eps_z ) {
		//se guarda c.box en lista de soluciones (Sout)
		Sout.push_back(c.box);
		c.box.set_empty();
		return;
	}

}

OptimizerMOP::Status OptimizerMOP::optimize(const IntervalVector& init_box) {

	nb_cells=0;

	buffer.flush();
	//LB.clear();
	UB.clear();

	//the box in cells have the n original variables plus the two objective variables (y1 and y2)
	Cell* root=new Cell(IntervalVector(n+2));

	root->box=init_box;

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
		/*
		 * TODO (IAZ): Stopping criteria considering a precision UB-LB
		 * Agregar restriccion a*z1+b*z2 = w para mejorar precision LB
		 */
		while (!buffer.empty()) {

		  if (trace >= 2) cout << buffer;

			/**
			 * TODO  (DA, MC): Extender CellBuffer (CellFeasibleDivingMOP) para que seleccione
			 * cajas con (z1=box[n].lb() ; z2=box[n+1].lb()) no dominados. Puede ser la caja con mayor dimension
			 * en alguno de los objetivos.
			 * Implementar otros metodos para comparar.
			 */

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
				exit(0);
			}
		}
	}
	catch (TimeOutException& ) {
		status = TIME_OUT;
		return status;
	}

	Timer::stop();
	time+= Timer::VIRTUAL_TIMELAPSE();

	if (Sout.empty())
		status=INFEASIBLE;
	else
		status=SUCCESS;

	return status;
}

void OptimizerMOP::report(bool verbose) {

	if (!verbose) {
		cout << get_status() << endl;
		cout << "LB:" << endl;
		list<  IntervalVector > :: iterator sol=Sout.begin();

		for(;sol!=Sout.end();sol++)
			cout << "(" << (*sol)[n].lb() << "," << (*sol)[n+1].lb() << ")" << endl;

		cout << "UB:" << endl;
		map< pair <double, double>, Vector > :: iterator ub=UB.begin();

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
