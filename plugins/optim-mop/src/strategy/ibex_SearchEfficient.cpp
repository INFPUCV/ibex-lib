//                                  I B E X
// File        : ibex_Optimizer.cpp
// Author      : Gilles Chabert, Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : May 14, 2012
// Last Update : December 24, 2012
//============================================================================

#include "ibex_SearchEfficient.h"
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

const double SearchEfficient::default_eps=0.01;

bool SearchEfficient::_plot = false;
double SearchEfficient::_min_ub_dist = 1e-7;
bool SearchEfficient::_cy_upper =false;
bool SearchEfficient::cy_contract_var = false;
bool SearchEfficient::_eps_contract = false;
double SearchEfficient::_rh = 0.1;
double* manhattan::efficient = NULL;
int manhattan::mode = 0;

template class CellSet<manhattan>;

SearchEfficient::SearchEfficient(int n, const Function &f1,  const Function &f2,
		Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder,
		 double eps, double rel_eps) : n(n),
                				ctc(ctc), bsc(bsc), buffer(buffer), goal1(f1), goal2(f2),
								finder(finder), trace(false), timeout(-1), status(SUCCESS),
                				time(0), nb_cells(0), eps(eps), efficient_solution(n) {

	if (trace) cout.precision(12);
}


SearchEfficient::~SearchEfficient() {

}


Interval SearchEfficient::eval_goal(const Function& goal, const IntervalVector& x, int n){
	//the objectives are set to 0.0
	IntervalVector xz(x);
	xz.resize(n+2);

	xz[n]=0.0;
	xz[n+1]=0.0;
	return goal.eval(xz);
}

IntervalVector SearchEfficient::deriv_goal(const Function& goal, const IntervalVector& x, int n){
	//the objectives are set to 0.0
	IntervalVector xz(x);
	xz.resize(n+2);

	xz[n]=0.0;
	xz[n+1]=0.0;
	IntervalVector g(goal.gradient(xz));
	g.resize(n);
	return g;
}

bool SearchEfficient::upper_bounding(const IntervalVector& box, const IntervalVector& init_box, opt_mode mode) {

	//We attempt to find two feasible points which minimize both objectives
	//and the middle point between them
	IntervalVector box2(box); box2.resize(n);
	IntervalVector xa(n), xb(n);
	finder.clear();

	list< pair <double, double> > points;
	list< pair< pair< double, double> , pair< double, double> > > segments;

	Vector mid=box2.mid();
	if (finder.norm_sys.is_inner(mid)){
		Vector v(2);
		v[0]=eval_goal(goal1,mid,n).ub();

		v[1]=eval_goal(goal2,mid,n).ub();



		if((mode==MINF1 && efficient[0] > v[0]) ||
			 (mode==MINF2 && efficient[1] > v[1]) ||
			 (mode==EFFICIENT && efficient[0] >= v[0] && efficient [1] >= v[1]))
		{
			//bajo el ub de la caja
			if(v[0] <= init_box[n].ub() && v[1] <= init_box[n+1].ub()){
			   efficient[0]=std::max(v[0], init_box[n].lb());
			   efficient[1]=std::max(v[1], init_box[n+1].lb());
			   efficient_solution = mid;
			}


		}
	}


	try{
		while(true){
			xa = finder.find(box2,box2,POS_INFINITY).first;

      //from phase 2, the points are in the line segment, thus more centered
			if(mode==EFFICIENT && finder.get_phase()<2)  continue;

			Vector v(2);
			v[0]=eval_goal(goal1,xa,n).ub();
			//cout << "dentro del while eterno (?):" << v[0] << endl;
			v[1]=eval_goal(goal2,xa,n).ub();
			//cout << "dentro del while eterno (?):"<< v[1] << endl;
			//verificar que el vector este dentro de la caja actual
			if (finder.norm_sys.is_inner(xa.mid())){
	         	if((mode==MINF1 && efficient[0] > v[0]) ||
				   (mode==MINF2 && efficient[1] > v[1]) ||
				   (mode==EFFICIENT && efficient[0] >= v[0] && efficient [1] >= v[1]))
				{
					//bajo el ub de la caja
					if(v[0] <= init_box[n].ub() && v[1] <= init_box[n+1].ub()){
					   efficient[0]= std::max(v[0], init_box[n].lb());
					   efficient[1]= std::max(v[1], init_box[n+1].lb());
					   efficient_solution = xa.mid();
				  }
				}
		    }
		}
	}catch (LoupFinder::NotFound& ) {
		return true;
	}


	return true;

}

void SearchEfficient::contract_and_bound(Cell& c, const IntervalVector& init_box, opt_mode mode) {

  //cout << efficient[0] << ", " << efficient[1] << endl;
	//cout << c.box << endl;
	switch (mode) {
		case EFFICIENT:
					if(efficient[0] < c.box[n].lb() || efficient[1] < c.box[n+1].lb()){
						c.box.set_empty();
						return;
					}
					if(efficient[0] < c.box[n].ub())
						c.box[n] = Interval(c.box[n].lb(),efficient[0]);

					if(efficient[1] < c.box[n+1].ub())
						c.box[n+1] = Interval(c.box[n+1].lb(),efficient[1]);
					break;

		case MINF1://disminuye y
					if(efficient[0] < c.box[n].lb()){
						c.box.set_empty();
						return;
					}
					if(efficient[0] < c.box[n].ub())
						c.box[n]=Interval(c.box[n].lb(),efficient[0]);
					break;
		case MINF2://disminuye x
				if(efficient[1] < c.box[n+1].lb()){
					c.box.set_empty();
					return;
				}
				if(efficient[1] < c.box[n+1].ub())
					c.box[n+1]=Interval(c.box[n+1].lb(),efficient[1]);
				break;
	}

	ctc.contract(c.box);


	if (c.box.is_empty()){
		return;
	}



}

void SearchEfficient::pre_optimize(const IntervalVector& init_box, Cell* root){
	//the box in cells have the n original variables plus the two objective variables (y1 and y2)
	root->box=init_box;
	root->prop.add(new BxpMOPData());

	nb_cells=0;
	buffer.flush();

	root->box=init_box;

	root->prop.add(new BxpMOPData());

	// add data required by the cell buffer
	buffer.add_property(init_box, root->prop);

	// add data required by the bisector
	bsc.add_property(init_box, root->prop);

	// add data required by the contractor
	ctc.add_property(init_box, root->prop);

	BxpMOPData::y1_init=eval_goal(goal1, root->box, n);
	BxpMOPData::y2_init=eval_goal(goal2, root->box, n);

	//cout << "esto es algo que ya olvide (?)" << endl;
	//cout << BxpMOPData::y1_init.ub() << endl;
	//cout << BxpMOPData::y2_init.lb() << endl;

	y1_ub.first=POS_INFINITY;
	y2_ub.second=POS_INFINITY;

	time=0;

	buffer.push(root);

    max_dist_eps = NEG_INFINITY;
}

SearchEfficient::Status SearchEfficient::optimize(const IntervalVector& init_box, opt_mode mode) {

	status=SUCCESS;
	//cout << "STATUSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS:  " << status << endl;

	Cell* root=new Cell(IntervalVector(n+2));
	pre_optimize(init_box, root);

	Timer timer;
	timer.start();

	set<Cell*> cells;
	cells.insert(root);



	efficient[0]= POS_INFINITY;
	efficient[1]= POS_INFINITY;

	manhattan::efficient=efficient;
	manhattan::mode = mode;

	try {
		bool server_pause=false;
		while (!buffer.empty()) {

			if(buffer.empty()) break;

			Cell *c = buffer.top();

			//cout << "distancia " << endl;
			//cout << max_area::distance(c) << endl;

			buffer.pop();
			cells.erase(c);

			if(manhattan::distance(c) < eps){
				delete c;
				continue;
			}

			nb_cells++;

			//EL DOMINIO ORIGINAL
			//	cout << "El dominio orifinal" << c->box[n].lb() << " , " << c->box[n].ub() << endl;
			contract_and_bound(*c, init_box, mode);
			//LUEGO DE HABER CORTADO LA CAJA
			//cout << "El dominio luego de haber cortado la caja" << c->box[n].lb() << " , " << c->box[n].ub() << endl;

			if (c->box.is_empty()) {
				delete c;
				continue;
			}

			upper_bounding(c->box, init_box, mode);

			if(manhattan::distance(c) < eps){
				delete c;
				continue;
			}

			pair<Cell*,Cell*> new_cells;
			bool atomic_box=false;
			try {
				new_cells=pair<Cell*,Cell*>(bsc.bisect(*c));
			}catch (NoBisectableVariableException& ) {
				if(new_cells.first){
					delete new_cells.first;
					delete new_cells.second;
				}
				continue;
			}

			delete c;

			buffer.push(new_cells.first);
			cells.insert(new_cells.first);

			buffer.push(new_cells.second);
			cells.insert(new_cells.second);


			if (timeout>0) timer.check(timeout); // TODO: not reentrant, JN: done
			time = timer.get_time();

		}
	}
	catch (TimeOutException& ) {
		status = TIME_OUT;
		//cout << "timeout" << endl;
	}


	timer.stop();
	time = timer.get_time();

//	py_Plotter::offline_plot(efficient.NDS2, NULL, "output2.txt");
	return status;
}




} // end namespace ibex
