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
#include "ibex_CellSet.h"

#include <float.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>

using namespace std;

namespace ibex {

const double OptimizerMOP::default_eps = 1e-5;


OptimizerMOP::OptimizerMOP(int n, const Array<NumConstraint>& ctrs, const Function &f1,  const Function &f2,
		Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder, double eps) : n(n),
                				ctc(ctc), bsc(bsc), buffer(buffer), ctrs(ctrs), goal1(f1), goal2(f2),
								finder(finder), eps(eps), trace(false), timeout(-1), status(SUCCESS),
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


bool OptimizerMOP::update_UB(const IntervalVector& box, int np) {

	list<Vector> feasible_points;

	//We attempt to find two feasible points which minimize both objectives
	//and the middle point between them
	IntervalVector box2(box); box2.resize(n);
	finder.find(box2,feasible_points,np);

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


	}
	//5. Si el mapa UB fue modificado retornar true, si no false
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

	/*====== filtering using the line z1+a*z2>w_lb and the ub_set ======*/

	//Conservativo: ok!
	if (c.get<CellBS>().depth>0 && discard_test(c.get<CellBS>().w_lb, c.get<CellBS>().a, c.box[n], c.box[n+1])){
		c.box.set_empty(); return;
	}

	/*====================================================================*/


	double z1, z2;
	pair <double, double> valueZ1;
	pair <double, double> valueZ2;

	valueZ1.first = NEG_INFINITY;
	valueZ2.second = NEG_INFINITY;

	IntervalVector box2=c.box;
	//cout << c.box[n] << "," << c.box[n+1] << endl;


	map< pair <double, double>, Vector >:: iterator ent1;
	for(ent1 = UB.begin(); ent1 != UB.end() ; ent1++) {
		z1 = (Interval(ent1->first.first)-Interval(eps)).ub(); // pair 1
		z2 = (Interval(ent1->first.second)-Interval(eps)).ub(); // pair 2
		if(z1 < c.box[n].lb() && z2 < c.box[n+1].lb()) {
			c.box.set_empty();
			return;
		}
		// contract c.box[n] && c.box[n+1] with PNS points
		if(z1 < c.box[n].lb()) {
			if(valueZ1.first == NEG_INFINITY || valueZ1.first < z1) {

				valueZ1 = ent1->first;
			}
		}
		if(z2 < c.box[n+1].lb()) {
			if(valueZ2.second == NEG_INFINITY || valueZ2.second < z2) {
				valueZ2 = ent1->first;
			}
		}
	}

	// contract c.box[n] && c.box[n+1] with PNS points
	if(valueZ1.first != NEG_INFINITY && valueZ1.second < c.box[n+1].ub() && valueZ1.second > c.box[n+1].lb()) {
		c.box[n+1] = Interval(c.box[n+1].lb(),valueZ1.second);
	}

	if(valueZ2.second != NEG_INFINITY && valueZ2.first < c.box[n].ub() && valueZ2.first > c.box[n].lb()) {
		c.box[n] = Interval(c.box[n].lb(), valueZ2.first);
	}



	/*================ contract x with f(x)=y1, f(x)=y2 and g(x)<=0 ================*/

	//adding the constraint z1+a*z2=w for filtering w
	IntervalVector box3(c.box);
	box3.resize(n+4);

	box3[n+2] = Interval(NEG_INFINITY, POS_INFINITY); /* w */
	if(c.box[n+1].diam() < POS_INFINITY)
		box3[n+3] = c.box[n].diam()/c.box[n+1].diam(); /* a */
	else
		box3[n+3] = 1.0;

	//the contraction is performed
	ctc.contract(box3);

	c.box=box3;
	c.box.resize(n+2);
	if (c.box.is_empty()) return;

	c.get<CellBS>().a = box3[n+3].mid();
	c.get<CellBS>().w_lb = box3[n+2].lb();



	/*========================= update loup =============================*/



	bool loup_ch=update_UB(c.box, 3);




	/*====================================================================*/
	/*double diamX = 0.0, valueBox;
	int i;
	for (i=0; i < n; i++) {
		valueBox = c.box[i].diam();
		if(diamX < valueBox) {
			diamX = valueBox;
		}
	}*/
	// Metodo de termino para las restricciones
	/*if ( diamX<=eps_x ) {
		//se guarda c.box en lista de soluciones (Sout)
		Sout.push_back(c.box);
		c.box.set_empty();
		return;
	}*/

	//we compute the min diameter between z1, z2 and w
	/*double diamZ = dist_w;
	for (i=n; i < n+2; i++) {
		valueBox = c.box[i].diam();
		if(diamZ > valueBox) {
			diamZ = valueBox;
		}
	}
	// Metodo de termino para las funciones objetivo
	if ( diamZ<=eps_z ) {
		//se guarda c.box en lista de soluciones (Sout)
		Sout.push_back(c.box);
		c.box.set_empty();
		return;
	}*/

}

OptimizerMOP::Status OptimizerMOP::optimize(const IntervalVector& init_box) {

	nb_cells=0;

	buffer.flush();
	//LB.clear();
	UB.clear();
	//the first point
	UB.insert(make_pair(make_pair(NEG_INFINITY,POS_INFINITY), Vector(1)));
	//the last point
	UB.insert(make_pair(make_pair(POS_INFINITY,NEG_INFINITY), Vector(1)));

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
	CellBS::z1_init=root->box[n];
	CellBS::z2_init=root->box[n+1];

	try {
		/** Criterio de termino: todas los nodos filtrados*/
		while (!buffer.empty()) {
		  plot();
		  getchar();
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

void OptimizerMOP::plot(){
	ofstream output;
	output.open("output.txt");
	list<  IntervalVector > :: iterator sol=Sout.begin();

	output << "(";
	for(;sol!=Sout.end();sol++){
		output << "{'pts':(" << (*sol)[n].lb() << "," << (*sol)[n+1].lb() << "),";
		output << "'diam_x': " << (*sol)[n].diam() << ",'diam_y': " << (*sol)[n+1].diam();
		output << "},";
	}
	output << ")" << endl;
	output << "[";
	map< pair <double, double>, Vector > :: iterator ub=UB.begin();
	for(;ub!=UB.end();ub++){
		output << "(" << ub->first.first << "," << ub->first.second << "),";
	}
	output << "]" << endl;
	output.close();
	// system("python3 plot.py");
}

void OptimizerMOP::report(bool verbose) {

	if (!verbose) {
		cout << get_status() << endl;
		cout << "LB:" << endl;
		list<  IntervalVector > :: iterator sol=Sout.begin();


		for(;sol!=Sout.end();sol++){
			cout << "(" << (*sol)[n].lb() << "," << (*sol)[n+1].lb() << ")" << endl;
		}

		cout << "UB:" << endl;
		map< pair <double, double>, Vector > :: iterator ub=UB.begin();

		for(;ub!=UB.end();ub++){
			cout << "(" << ub->first.first << "," << ub->first.second << ")" << endl;
		}

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
