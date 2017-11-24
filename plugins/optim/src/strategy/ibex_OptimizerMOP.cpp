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

const double OptimizerMOP::default_abs_eps = 1e-7;
const double OptimizerMOP::default_rel_eps = 0.1;
bool OptimizerMOP::_plot = false;
int OptimizerMOP::_nb_ub_sols = 3;
double OptimizerMOP::_min_ub_dist = 1e-7;

map< pair <double, double>, IntervalVector > OptimizerMOP::UB;


OptimizerMOP::OptimizerMOP(int n, const Array<NumConstraint>& ctrs, const Function &f1,  const Function &f2,
		Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder, double rel_eps, double abs_eps) : n(n),
                				ctc(ctc), bsc(bsc), buffer(buffer), ctrs(ctrs), goal1(f1), goal2(f2),
								finder(finder), trace(false), timeout(-1), status(SUCCESS),
                				time(0), nb_cells(0), rel_eps(rel_eps), abs_eps(abs_eps), nb_sols(0),
								y1_max(NEG_INFINITY), y2_max(NEG_INFINITY) {

	if (trace) cout.precision(12);
}

OptimizerMOP::~OptimizerMOP() {

}


Interval OptimizerMOP::eval_goal(const Function& goal, IntervalVector& x){
	//the objectives are set to 0.0
	IntervalVector xz(x);
	xz.resize(n+2);

	xz[n]=0.0;
	xz[n+1]=0.0;
	return goal.eval(xz);
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

		IntervalVector vec=*it;

		//3. Se evalua el punto usando funciones objetivo (goal1 y goal2)
		pair< double, double> eval = make_pair(eval_goal(goal1,vec).ub(), eval_goal(goal2,vec).ub());

		//4. Insertar en mapa UB (si es no dominada) y actualizar eliminar soluciones dominadas de UB
		std::map<pair<double, double>, IntervalVector>::iterator it2 = UB.lower_bound(eval);

		//there is an equivalent point
		if(it2->first == eval) continue;
		it2--;
		//it is dominated by the previous ub point
		if(eval.second >= it2->first.second)	continue;



		/**** UB correction ****/
		if(finder.ub_correction(vec.mid(), vec)){
			//cout << "pre:" << eval.first << "," << eval.second << endl;
			eval = make_pair(eval_goal(goal1,vec).ub(), eval_goal(goal2,vec).ub());
			//cout << "corrected:" << eval.first << "," << eval.second << endl;
		}
		else continue;


		it2= UB.lower_bound(eval);

		//there is an equivalent point
		if(it2->first == eval) continue;
		it2--;
		//it is dominated by the previous ub point
		if(eval.second >= it2->first.second)	continue;

		/**** end UB correction ****/

		bool domine=false;
		for(it2++; it2!=UB.end(); ){
			if(eval.second > it2->first.second) break;
			std::map<pair<double, double>, IntervalVector>::iterator aux = it2;
			++aux;
			UB.erase(it2);
			it2 = aux;
			domine=true;
		}


		//the point is inserted in UB only if its distance to the neighbor points is greater than (abs_eps/2.0)
		if(domine || std::min(it2->first.first - eval.first,  eval.second - it2->first.second) >= _min_ub_dist){
			//it is not dominated and we remove the new dominated points


			if(eval.first < y1_ub.first) y1_ub=eval;
			if(eval.second < y2_ub.second) y2_ub=eval;

			UB.insert(make_pair(eval, vec));
			new_ub = true;
		}else{
			it2--;
			if( std::min(eval.first - it2->first.first,  it2->first.second - eval.second) >= _min_ub_dist ){
				//it is not dominated and we remove the new dominated points

				if(eval.first < y1_ub.first) y1_ub=eval;
				if(eval.second < y2_ub.second) y2_ub=eval;

				UB.insert(make_pair(eval, vec));
				new_ub = true;
			}

		}
		if(trace) cout << eval.first  <<"," << eval.second << endl;


		//UB.insert(make_pair(eval, vec));
		//new_ub = true;


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
		buffer_cells.insert(&c);

		nb_cells++;
	}
}

void OptimizerMOP::contract_and_bound(Cell& c, const IntervalVector& init_box) {

	/*====================================================================*/

	double z1, z2;
	pair <double, double> valueZ1;
	pair <double, double> valueZ2;

	valueZ1.first = NEG_INFINITY;
	valueZ2.second = NEG_INFINITY;

	IntervalVector box2=c.box;
	//cout << c.box[n] << "," << c.box[n+1] << endl;


	map< pair <double, double>, IntervalVector >:: iterator ent1=UB.lower_bound(make_pair(c.box[n].lb(),c.box[n+1].lb()));
    ent1--;
	for(; ent1 != UB.end() ; ent1++) {
		z1 = ent1->first.first; // pair 1
		z2 = ent1->first.second; // pair 2
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

		if(c.box[n].lb() > z2) break;
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

	map< pair <double, double>, IntervalVector >::iterator it = UB.upper_bound(make_pair(c.box[n].lb(), NEG_INFINITY));
	it--;
	map< pair <double, double>, IntervalVector >::iterator it2 = UB.lower_bound(make_pair(c.box[n].ub(), NEG_INFINITY));

    if(it->first.first != NEG_INFINITY && it2->first.second != NEG_INFINITY && it->first!=it2->first){
    	box3[n+3] = (it2->first.first-it->first.first)/(it->first.second-it2->first.second);
    }else if(c.box[n+1].diam() < POS_INFINITY)
		box3[n+3] = c.box[n].diam()/c.box[n+1].diam(); /* a */
	else
		box3[n+3] = 1.0;

  //setting w_ub with the UB points
	  double w_ub=NEG_INFINITY;

    if(UB.size()==2)  w_ub = POS_INFINITY;
		else{
			it = UB.lower_bound(make_pair(c.box[n].lb(), NEG_INFINITY));
			it--;
			while(it!=UB.end()){
				pair <double, double> p = it->first; it++;
				if(it==UB.end() || p.second < c.box[n+1].lb()) break;
				pair <double, double> p2 = it->first;
				pair <double, double> pmax= make_pair(p2.first, p.second);
        //cout << "pmax:" << pmax.first << "," << pmax.second << endl;
				if(pmax.first==POS_INFINITY || pmax.second==POS_INFINITY)
				   w_ub = POS_INFINITY;
				else{
		  		double ww = ( Interval(pmax.first) + box3[n+3]*Interval(pmax.second) ).ub();
		  		if(w_ub < ww )  w_ub = ww;
				}
			}
		}

		//cout << w_ub << endl;
		//box3[n+2] = Interval(NEG_INFINITY, POS_INFINITY); /* w */
		box3[n+2] = Interval(NEG_INFINITY, w_ub); /* w */


	//the contraction is performed
	ctc.contract(box3);


	c.box=box3;
	c.box.resize(n+2);
	if (c.box.is_empty()) return;

	c.get<CellBS>().a = box3[n+3].mid();
	c.get<CellBS>().w_lb = box3[n+2].lb();



	/*========================= update loup =============================*/

	bool loup_ch=update_UB(c.box, _nb_ub_sols);
	if(trace && loup_ch && _plot) { plot(&c);  getchar(); }

}

OptimizerMOP::Status OptimizerMOP::optimize(const IntervalVector& init_box) {

	nb_cells=0;

	nb_sols=0;

	buffer.flush();
	buffer_cells.clear();
	//LB.clear();
	UB.clear();
	//the first point
	UB.insert(make_pair(make_pair(NEG_INFINITY,POS_INFINITY), Vector(1)));
	//the last point
	UB.insert(make_pair(make_pair(POS_INFINITY,NEG_INFINITY), Vector(1)));

	//the box in cells have the n original variables plus the two objective variables (y1 and y2)
	Cell* root=new Cell(IntervalVector(n+2));

	root->box=init_box;

	y1=eval_goal(goal1, root->box);
	y2=eval_goal(goal2, root->box);
  //Initial points
	CellBS::y1_init=y1;
	CellBS::y2_init=y2;
	abs_eps=rel_eps*CellBS::y1_init.diam();

	y1_ub.first=POS_INFINITY;
	y2_ub.second=POS_INFINITY;

	y1_max=NEG_INFINITY;
	y2_max=NEG_INFINITY;


	// add data required by the bisector
	bsc.add_backtrackable(*root);

	// add data required by the buffer
	buffer.add_backtrackable(*root);

	time=0;
	Timer timer;
	timer.start();

	handle_cell(*root,init_box);


	try {
		/** Criterio de termino: todas los nodos filtrados*/
		while (!buffer.empty()) {

		  if (trace >= 2) cout << buffer;


			Cell *c = buffer.pop();
			buffer_cells.erase(c);

			//cout << "abs_eps:" << abs_eps << endl;
			double dist=distance2(c);

			if(dist < abs_eps){
				if( dist>0 ) nb_sols++;
				if(trace) cout << "dist: " << dist << ":" << abs_eps << endl;

				if(LB.empty()) {
					if(_plot) {plot(c);  getchar();}
					y2_max=y1_ub.second;
					y1_max=y2_ub.first;
				}

                //for obtaining the nadir point (y1max, y2max) used to compute the hypervolume
				if(y1_max < c->box[n].ub() &&  c->box[n+1].lb()<y2_ub.second)
					y1_max=c->box[n].ub();
				if(y2_max < c->box[n+1].ub() &&  c->box[n].lb()<y1_ub.first)
					y2_max=c->box[n+1].ub();

				double ya2 = ((c)->get<CellBS>().w_lb-c->box[n].lb())/(c)->get<CellBS>().a;
				if(trace) cout << "ya2:" << ya2 << endl;
				if(ya2 > c->box[n+1].lb() && ya2 < c->box[n+1].ub()){
					double yb1 = (c)->get<CellBS>().w_lb-((c)->get<CellBS>().a*c->box[n+1].lb());
					if(trace) cout << "yb1:" << yb1 << endl;
					if(yb1 > c->box[n].lb() && yb1 < c->box[n].ub()){
						    insert_lb_segment( point2(c->box[n].lb(),ya2),
						    point2(yb1 ,  c->box[n+1].lb()) );
					}else{
						double yb2=((c)->get<CellBS>().w_lb-c->box[n].ub())/(c)->get<CellBS>().a;
						insert_lb_segment( point2(c->box[n].lb(),ya2),
								point2(c->box[n].ub() ,  yb2 ) );
					}

				}else if(ya2 >= c->box[n+1].ub()){
					double ya1=(c)->get<CellBS>().w_lb-((c)->get<CellBS>().a*c->box[n+1].ub());
					double yb1 = (c)->get<CellBS>().w_lb-((c)->get<CellBS>().a*c->box[n+1].lb());
					if(trace) cout << "ya1:" << ya1 << endl;
					if(trace) cout << "yb1:" << yb1 << endl;
                    if(yb1 > c->box[n].lb() && yb1 < c->box[n].ub()){

						insert_lb_segment( point2(ya1, c->box[n+1].ub()),
									point2(yb1 ,  c->box[n+1].lb() ) );
					}else{

						double yb2=((c)->get<CellBS>().w_lb-c->box[n].ub())/(c)->get<CellBS>().a;
						insert_lb_segment( point2(ya1, c->box[n+1].ub()),
									point2(c->box[n].ub() ,  yb2) );
					}


				}else
				  insert_lb_segment(point2(c->box[n].lb(),c->box[n+1].lb()),point2(c->box[n].lb(),c->box[n+1].lb()));

				//trace=0;
                if(_plot && buffer.empty())
        	       plot(c);



				delete c; continue;
			}

			//cout << c->get<CellBS>().ub_distance << ":" << abs_eps << endl;

		 //	plot(c); //getchar();

			try {
				pair<IntervalVector,IntervalVector> boxes=bsc.bisect(*c);

				pair<Cell*,Cell*> new_cells=c->bisect(boxes.first,boxes.second);

				delete c; // deletes the cell.

				handle_cell(*new_cells.first, init_box);
				handle_cell(*new_cells.second, init_box);

				if (timeout>0) timer.check(timeout); // TODO: not reentrant, JN: done
				time = timer.get_time();

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
		cout << "timeout" << endl;
		cout << "lb-hypervolume:" << compute_lb_hypervolume() << endl;
		cout << "ub-hypervolume:" << compute_ub_hypervolume() << endl;

		cout << "diff-hypervolume:" << (compute_lb_hypervolume()-compute_ub_hypervolume())/y1.diam() << endl;

		return status;
	}


	timer.stop();
	time = timer.get_time();


	cout << "lb-hypervolume:" << compute_lb_hypervolume() << endl;
	cout << "ub-hypervolume:" << compute_ub_hypervolume() << endl;

	cout << "diff-hypervolume:" << (compute_lb_hypervolume()-compute_ub_hypervolume())/y1.diam() << endl;

	//if (Sout.empty())
		status=INFEASIBLE;
	//else
		status=SUCCESS;

	return status;
}

void OptimizerMOP::plot(Cell* c){
	ofstream output;
	output.open("output.txt");
	set<  Cell* > :: iterator cell=buffer_cells.begin();


	output << "(";
	if(c){
		output << "{'pts':(" << c->box[n].lb() << "," <<  c->box[n+1].lb() << "),";
		output << "'diam_x': " <<  c->box[n].diam() << ",'diam_y': " << c->box[n+1].diam()<< ",";
		output << "'pA':(" << c->box[n].lb() <<"," <<  (((c)->get<CellBS>().w_lb-c->box[n].lb())/(c)->get<CellBS>().a)   << "),";
		output << "'pB':(" << (c->get<CellBS>().w_lb-c->get<CellBS>().a*c->box[n+1].lb()) <<"," <<  c->box[n+1].lb()  << ")";
		output << "},";
  }

	for(;cell!=buffer_cells.end();cell++){
		if(distance2(*cell) < 0){continue;}

		output << "{'pts':(" << (*cell)->box[n].lb() << "," <<  (*cell)->box[n+1].lb() << "),";
		output << "'diam_x': " <<  (*cell)->box[n].diam() << ",'diam_y': " <<  (*cell)->box[n+1].diam() << ",";
		output << "'pA':(" << (*cell)->box[n].lb() <<"," <<  (((*cell)->get<CellBS>().w_lb-(*cell)->box[n].lb())/(*cell)->get<CellBS>().a)   << "),";
		output << "'pB':(" << ((*cell)->get<CellBS>().w_lb-(*cell)->get<CellBS>().a*(*cell)->box[n+1].lb()) <<"," <<  (*cell)->box[n+1].lb()  << ")";
		output << "},";
	}
	output << ")" << endl;


	output << "[";

	map< pair <double, double>, IntervalVector > :: iterator ub=UB.begin();
	for(;ub!=UB.end();ub++){
		output << "(" << ub->first.first << "," << ub->first.second << "),";
	}
output << "]" << endl;

	output << "[";
	set< point2 > :: iterator lb=LB.begin();
	for(;lb!=LB.end();lb++){
		if(lb->x.mid() >= 1e9)  output << "(inf,";
		else if(lb->x.mid() <= -1e9)  output << "(-inf,";
		else output << "(" << lb->x.mid() << ",";

		if(lb->y.mid() >= 1e9)  output << "inf),";
		else if(lb->y.mid() <= -1e9)  output << "-inf),";
		else output << lb->y.mid() << "),";
	}

	output << "]" << endl;
	output.close();
	// system("python3 plot.py");
}

void OptimizerMOP::report(bool verbose) {

	if (!verbose) {
		cout << get_status() << endl;
		cout << "LB:" << endl;
		//list<  IntervalVector > :: iterator sol=Sout.begin();


		//for(;sol!=Sout.end();sol++){
			//cout << "(" << (*sol)[n].lb() << "," << (*sol)[n+1].lb() << ")" << endl;
		//}

		cout << "UB:" << UB.size() << endl;
		map< pair <double, double>, IntervalVector > :: iterator ub=UB.begin();
/*
		for(;ub!=UB.end();ub++){
			cout << ub->second << endl;
			cout << "(" << ub->first.first << "," << ub->first.second << ")" << endl;
		}
*/
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


} // end namespace ibex
