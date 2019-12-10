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
#include <float.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <set>
//#include "ibex_CellSet.h"

#ifndef cdata
#define cdata ((BxpMOPData*) c->prop[BxpMOPData::id])
#endif

using namespace std;

namespace ibex {

const double OptimizerMOP::default_eps=0.01;

bool OptimizerMOP::_plot = false;
double OptimizerMOP::_min_ub_dist = 1e-7;
bool OptimizerMOP::_cy_upper =false;
bool OptimizerMOP::cy_contract_var = false;
bool OptimizerMOP::_eps_contract = false;
double OptimizerMOP::_rh = 0.1;

int OptimizerMOP::nb_ObjFunc = 5; //MAKE NUMBER OF OBJECTIVE FUNCTION Todo BxpMOPData use it
bool OptimizerMOP::imprimir = false;

/*NEW CONSTRUCTOR,
 * WIP to adapt goals, eliminate goal1 and goal2
 *
 */
OptimizerMOP::OptimizerMOP(int n,Array<const Function> fObj,
		Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder,
		Mode nds_mode, Mode split_mode, double eps, double rel_eps) : n(n),
                				ctc(ctc), bsc(bsc), buffer(buffer),
								goals(fObj), finder(finder), trace(false), timeout(-1), status(SUCCESS),
                				time(0), nb_cells(0), eps(eps), nds_mode(nds_mode), split_mode(split_mode),
												rel_eps(rel_eps){

	if (trace) cout.precision(12);
}

OptimizerMOP::~OptimizerMOP() {

}


Interval OptimizerMOP::eval_goal(const Function& goal, const IntervalVector& x, int n, int nFuncObj){
	//the objectives are set to 0.0
	IntervalVector xz(x);
	xz.resize(n+nFuncObj);
	for(int i=0; i<nFuncObj; i++) xz[n+i]=0.0;
//	std::cout<<"x = "<<x<<endl;
//	std::cout <<"goal = "<<goal<<endl;
//	//std::cout <<"xz(n+nFuncObj) ="<<xz<<endl;
//	std::cout<<"goal.eval(xz) = "<<goal.eval(xz)<<endl;
//	std::cout<<"\n";
	return goal.eval(xz);
}

IntervalVector OptimizerMOP::deriv_goal(const Function& goal, const IntervalVector& x, int n){
	//the objectives are set to 0.0
	IntervalVector xz(x);
	xz.resize(n+2);

	xz[n]=0.0;
	xz[n+1]=0.0;
	IntervalVector g(goal.gradient(xz));
	g.resize(n);
	return g;
}

bool OptimizerMOP::upper_bounding(const IntervalVector& box) {
	if(imprimir) std::cout<<"\t \n Upper Bounding ++++++++++++++ \n";
	//We attempt to find two feasible points which minimize both objectives
	//and the middle point between them
	IntervalVector box2(box); box2.resize(n);
	IntervalVector xa(n), xb(n);
	finder.clear();

	list< pair <double, double> > points;
	list< pair< pair< double, double> , pair< double, double> > > segments;


	Vector mid=box2.mid();
	if(imprimir) std::cout<<"box mid = "<<mid<<endl;
	if (finder.norm_sys.is_inner(mid)){
		Vector v(goals.size());
		for(int i=0; i<goals.size(); i++) {
			//evaluation of function n in mid point of box of variables (x1, x2, x3, ...)
			v[i] = eval_goal(goals[i], mid, n, goals.size()).ub();
		}
		if(imprimir) std::cout <<"eval fn(mid): v= "<<v<<endl;
		//addPoint ndsh2
		ndsh2.addPoint(v, NDS_X(mid));

	}

	return true;

}

//Funtion is not executed (already commented in branch and bound method)
void OptimizerMOP::discard_generalized_monotonicty_test(IntervalVector& box, const IntervalVector& initbox){
//	IntervalVector grad_f1= goal1.gradient(box);
//	IntervalVector grad_f2= goal2.gradient(box);
	IntervalVector grad_f1= goals[0].gradient(box);
	IntervalVector grad_f2= goals[1].gradient(box);

	IntervalVector new_box(box);

	//bool discard=false;

	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			if(grad_f1[i].lb() > 0.0 && grad_f2[j].lb()>0.0 &&
				    (i==j || (Interval(grad_f2[i].lb()) / grad_f2[j].lb() -  Interval(grad_f1[i].lb()) / grad_f1[j].lb()).lb()  > 0.0) ){

				    if(i==j) new_box[i] = box[i].lb(); //simple monotonicity test

					if( box[i].lb() != initbox[i].lb() && box[j].lb() != initbox[j].lb() ){
						if(is_inner_facet(box,i,box[i].lb()) && is_inner_facet(box,j,box[j].lb())){
							box.set_empty();
							return;
						}
					}
			}else if(grad_f1[i].ub() < 0.0 && grad_f2[j].lb()>0.0 &&
					(Interval(grad_f1[i].ub()) / grad_f1[j].lb() -  Interval(grad_f2[i].ub()) / grad_f2[j].lb()).lb()  > 0.0) {

				    if( box[i].ub() == initbox[i].ub() && box[j].lb() != initbox[j].lb() ) {
						if(is_inner_facet(box,i,box[i].ub()) && is_inner_facet(box,j,box[j].lb())){
							box.set_empty();
							return;
						}
					}

			}else if(grad_f1[i].lb() > 0.0 && grad_f2[j].ub()<0.0 &&
					(Interval(grad_f1[i].lb()) / grad_f1[j].ub() -  Interval(grad_f2[i].lb()) / grad_f2[j].ub()).lb()  > 0.0) {

				    if( box[i].lb() != initbox[i].lb() && box[j].ub() != initbox[j].ub() )
				    {
						if(is_inner_facet(box,i,box[i].lb()) && is_inner_facet(box,j,box[j].ub())){
							box.set_empty();
							return;
						}
					}

			 }else if(grad_f1[i].ub() < 0.0 && grad_f2[j].ub() < 0.0 &&
					(i==j || (Interval(grad_f2[i].ub()) / grad_f2[j].ub() -  Interval(grad_f1[i].ub()) / grad_f1[j].ub()).lb()  > 0.0) ){

					if(i==j) new_box[i] = box[i].ub(); //simple monotonicity test

					if( box[i].ub() != initbox[i].ub() && box[j].ub() != initbox[j].ub() )
					{
						if(is_inner_facet(box,i,box[i].ub()) && is_inner_facet(box,j,box[j].ub())){
							box.set_empty();
							return;
						}
					}
			}
		}
	}

	IntervalVector bb=new_box;
	bb.resize(n);

	if(finder.norm_sys.is_inner(bb))
		box=new_box;

}


void OptimizerMOP::dominance_peeler2(IntervalVector& box, list < Vector >& inpoints){
	Vector firstp=inpoints.front();
	Vector lastp=inpoints.back();
	// contract c.box[n]
	if(firstp[1] < box[1].ub())
		box[1] = Interval(box[1].lb(),firstp[1]);
	// contract c.box[n+1]
	if(lastp[0] < box[0].ub() )
		box[0] = Interval(box[0].lb(), lastp[0]);

}

//Todo Dominance Peeler for n objective functions
void OptimizerMOP::dominance_peeler_n(IntervalVector& boxy, list < pair< Vector, int > >& cuttingPoints){

	for(auto itr = cuttingPoints.begin(); itr != cuttingPoints.end(); ++itr){
		if(imprimir) std::cout<<"itr->first = "<<itr->first<<endl;
		if(imprimir) std::cout<<"itr->second = "<<itr->second<<"\n";
		boxy[itr->second] = Interval(boxy[itr->second].lb(), itr->first[itr->second]);
	}
}

/*
 * Retorna un vector de intervalos de la caja de variables (cell) asociada con las funciones objetivos
 *
 * Eg:
 * 	Variables x1, x2, x3, z1, z2, z3, z4 		(4 objectives function)
 *
 * 	cell box = ([0, 5] ; [0, 3] ; [0, 3] ; [0, 62.55000000000001] ; [1, 50] ; [0, 8] ; [-3, 3])
 * 	box_y = ([0, 62.55000000000001] ; [1, 50] ; [0, 8] ; [-3, 3])
 */
IntervalVector OptimizerMOP::get_box_y(IntervalVector &box){
	IntervalVector box_y(goals.size());
	int n=box.size()-goals.size();
	for(int i=0; i< goals.size(); i++) box_y[i]=box[n+i];
	return box_y;
}
/*
 * Reemplaza los valores de box asociado a las funciones objetivos con los valores de box_y
 */
void OptimizerMOP::set_box_y(IntervalVector &box, IntervalVector &box_y){
	for(int i=0; i<box_y.size(); i++) box[n+i] = box_y[i];
	//cout<<"set box_y = "<<get_box_y(box)<<endl;
}

void OptimizerMOP::contract_and_bound(Cell& c, const IntervalVector& init_box) {
	if(imprimir) std::cout<<"\t \nContract and Bound ++++++++++++++ \n\n";

	//Obtain box from cell
	IntervalVector boxy=get_box_y(c.box);
	if(imprimir) std::cout<<"initial boxy = "<<boxy<<endl;

//TODO- Return a list of point non-dominated by box_y.lb()
	list< pair< Vector, int >> cuttingPoints = ndsh2.cutting_points(boxy.lb(), boxy.ub());

//TODO Dominance Peeler contract box_y
	if(imprimir) cout<<"Old DomPeel boxy = "<<boxy<<endl;
	dominance_peeler_n(boxy, cuttingPoints);
	set_box_y(c.box, boxy);
	if(imprimir) cout<<"New DomPeel boxy = "<<boxy<<endl;

	//discard_generalized_monotonicty_test(c.box, init_box);

	if (c.box.is_empty()) return;

	if(cy_contract_var){
		//cy_contract2(c,inner_segments);
	}else{
		//Apply contraction with acidhc4 on the box (Default)
		ctc.contract(c.box);
	}
	if (c.box.is_empty()) return;
}

//WIP add array of goals
void OptimizerMOP::pre_optimize(const IntervalVector& init_box, Cell* root){
	//the box in cells have the n original variables plus the two objective variables (y1 and y2)
//	root->box=init_box;
//	root->prop.add(new BxpMOPData());

	nb_cells=0;
	buffer.flush();
	ndsH.clear();

	ndsh2.nObjFunc = goals.size();
	ndsh2.clear();
	MOPData = new BxpMOPData();

	root->box=init_box;
	root->prop.add(MOPData);

	// add data required by the cell buffer
	buffer.add_property(init_box, root->prop);

	// add data required by the bisector
	bsc.add_property(init_box, root->prop);

	// add data required by the contractor
	ctc.add_property(init_box, root->prop);


	//Evaluate every objective function into initial domains value of original variables
	for(int i=0; i<goals.size();i++){
		Interval goal_evaluation = eval_goal(goals[i], root->box, n, goals.size());
		MOPData->yn_init.put(i, goal_evaluation);
	}

	//Init min feasible value for each objective
	double inf = POS_INFINITY;
	yn_ub.resize(goals.size());
	for(int i=0; i<goals.size();i++) yn_ub.set_ref(i,inf);

	time=0;
	buffer.push(root);
    max_dist_eps = NEG_INFINITY;
}

//WIP to add goals
OptimizerMOP::Status OptimizerMOP::optimize(const IntervalVector& init_box) {

	status=SUCCESS;

	//the box in cells have the n original variables plus the objective variables (y1, y2, ..., yn)
	Cell* root=new Cell(IntervalVector(n+goals.size()));
	pre_optimize(init_box, root);
	Timer timer;
	timer.start();
	set<Cell*> cells;
	cells.insert(root);

	IntervalVector focus(goals.size());
	//Set focus, utilizado en OptimizerMOPserver, pyPlotter y en el calculo de eps cuando rel_eps>0
	for(int i=0; i<goals.size();i++) {
		//focus[i] = BxpMOPData::yn_init.operator [](i);
		focus[i] = MOPData->yn_init.operator [](i);
	}

	try {
		bool server_pause=false;
		while (!buffer.empty()) {

			if(buffer.empty()) break;

			Cell *c = buffer.top();
			//std::cout<<"ub_distance = "<<cdata->ub_distance<<endl;
			if(imprimir) std::cout <<"-------------------------------------------------------------------------------------------\nCurrent box cell"<<c->box<<endl;
			//if(rel_eps>0.0)	eps=std::max(focus[0].diam(),focus[1].diam())*rel_eps;

			buffer.pop();
			cells.erase(c);

//TODO si la cell con la mayor distancia encontrada (ub_distance) es menor a eps, elimina esa cell y toma la siguiente del buffer
			if(cdata->ub_distance <= eps){
				std::cout<<"ub_distance final = "<<cdata->ub_distance<<endl;
				delete c;

				if(dynamic_cast<DistanceSortedCellBufferMOP*>(&buffer)) break;
				else continue;
			}

			nb_cells++;

//TODO dominance Peeler and ndsh.nondominatedPoints
			contract_and_bound(*c, init_box);

			if (c->box.is_empty()) {
				delete c;
				continue;
			}
			if(imprimir) std::cout <<"Current box cell contracted"<<c->box<<endl;


			upper_bounding(c->box);

			pair<Cell*,Cell*> new_cells;
			bool atomic_box=false;
			try {
				new_cells=pair<Cell*,Cell*>(bsc.bisect(*c));
			}
			catch (NoBisectableVariableException& ) {
				atomic_box=true;
			}

			double dist2 = 0.0 ;
			if(!atomic_box) {
//TODO Distance ndsh2
				dist2=ndsh2.distance(c->box);
			}

//Todo Se elimina la caja
			if(atomic_box || dist2 < eps){
				if(imprimir) std::cout<<"caja eliminada  dist< eps or atomic_box"<<endl;
				if(new_cells.first){
					delete new_cells.first;
					delete new_cells.second;
				}

			delete c;

			continue;

			}

			delete c; // deletes the cell.

			buffer.push(new_cells.first);
			cells.insert(new_cells.first);

			buffer.push(new_cells.second);
			cells.insert(new_cells.second);

			if(imprimir) std::cout<<"new cell bisected, box 1 = "<<new_cells.first->box<<endl;
			if(imprimir) std::cout<<"new cell bisected, box 2 = "<<new_cells.second->box<<endl;


			if (timeout>0) timer.check(timeout); // TODO: not reentrant, JN: done
			time = timer.get_time();


		}
	}
	catch (TimeOutException& ) {
		status = TIME_OUT;
		cout << "timeout" << endl;
	}

	if(!buffer.empty()){
		Cell *c = buffer.top();
		std::cout<<"ub_distance final = "<<cdata->ub_distance<<endl;
	}


	timer.stop();
	time = timer.get_time();

//	py_Plotter::offline_plot(ndsH.NDS2, NULL, "output2.txt");
	py_Plotter::offline_plot(ndsh2.NDS, NULL, "output2.txt", "solutions.txt");
	return status;
}

void OptimizerMOP::report(bool verbose) {

	if (!verbose) {
     cout << endl 	<< "time 	#nodes 		|Y|" << endl;
		 //cout << get_time() << " " << get_nb_cells() << " " << ndsH.size() << endl;
		 cout << get_time() << " " << get_nb_cells() << " " << ndsh2.size() <<endl;
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

	cout << " cpu time used: " << get_time() << "s." << endl;
	cout << " number of cells: " << get_nb_cells() << endl;
	cout<< "number of cells not executed: "<< buffer.size()<<endl;

	cout << " number of solutions: "  << ndsh2.size() << endl;
	for(auto ub : ndsh2.NDS){
		//cout << "(" << ub.first[0] << "," << ub.first[1] << ")"  << endl;
		cout<<"(";
		for(auto i=0; i<ub.first.size(); i++){
			cout<<ub.first[i];
			if(i+1!=ub.first.size()) cout<<", ";
		}
		cout<<")"<<endl;
	}



}


} // end namespace ibex
