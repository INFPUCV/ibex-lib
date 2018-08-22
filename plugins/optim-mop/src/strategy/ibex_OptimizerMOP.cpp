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
#include <iostream>
//#include "ibex_CellSet.h"

using namespace std;

namespace ibex {

const double OptimizerMOP::default_eps=0.01;

bool OptimizerMOP::_plot = false;
double OptimizerMOP::_min_ub_dist = 1e-7;
bool OptimizerMOP::_cy_upper =false;
//bool OptimizerMOP::_hv =false;
bool OptimizerMOP::cy_contract_var = false;
bool OptimizerMOP::_eps_contract = false;



OptimizerMOP::OptimizerMOP(int n, const Function &f1,  const Function &f2,
		Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder, Mode nds_mode, double eps) : n(n),
                				ctc(ctc), bsc(bsc), buffer(buffer), goal1(f1), goal2(f2),
								finder(finder), trace(false), timeout(-1), status(SUCCESS),
                				time(0), nb_cells(0), eps(eps), nds_mode(nds_mode)
								 {

	py_Plotter::n=n;

	if (trace) cout.precision(12);
}


OptimizerMOP::~OptimizerMOP() {

}


Interval OptimizerMOP::eval_goal(const Function& goal, const IntervalVector& x, int n){
	//the objectives are set to 0.0
	IntervalVector xz(x);
	xz.resize(n+2);

	xz[n]=0.0;
	xz[n+1]=0.0;
	return goal.eval(xz);
}

IntervalVector OptimizerMOP::deriv_goal(const Function& goal, const IntervalVector& x, int n){
	//the objectives are set to 0.0
	IntervalVector xz(x);
	xz.resize(n+2);

	xz[n]=0.0;
	xz[n+1]=0.0;
	return goal.gradient(xz);
}



bool OptimizerMOP::update_NDS2(const IntervalVector& box) {

	//We attempt to find two feasible points which minimize both objectives
	//and the middle point between them
	IntervalVector box2(box); box2.resize(n);
	IntervalVector xa(n), xb(n);
	finder.clear();

	list< pair <double, double> > points;
	list< pair< pair< double, double> , pair< double, double> > > segments;

	try{
		xa = finder.find(box2,box2,POS_INFINITY).first;

	}catch (LoupFinder::NotFound& ) {
		return false;
	}

	try{
		xb = finder.find(box2,box2,POS_INFINITY).first;
	}catch (LoupFinder::NotFound& ) {
		ndsH.addPoint(make_pair(eval_goal(goal1,xa,n).ub(), eval_goal(goal2,xa,n).ub()));
		return true;
	}

	hamburger(xa, xb);

	return true;

}

void OptimizerMOP::discard_generalized_monotonicty_test(IntervalVector& box, const IntervalVector& initbox){
	IntervalVector grad_f1= goal1.gradient(box);
	IntervalVector grad_f2= goal2.gradient(box);

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

void OptimizerMOP::cy_contract2(Cell& c, list <pair <double,double> >& inpoints){
	IntervalVector& box = c.box;
	IntervalVector box3(box);
	box3.resize(n+4);
	box3[n+3] = 1.0;

	pair <double, double> firstp=inpoints.front();
	pair <double, double> lastp=inpoints.back();

	//Pent of cy

	if(firstp!=lastp)
		box3[n+3] = (lastp.first-firstp.first)/(firstp.second-lastp.second);
	if(box3[n+3]==0.0 || box3[n+3]==1.0 || box3[n+3].is_empty())
		box3[n+3] = box[n].diam()/box[n+1].diam(); // a

   //cout <<box3[n+3]  << endl;

   //setting w_ub with the NDS points
	 double w_ub=POS_INFINITY;

   if(_cy_upper){
			w_ub=NEG_INFINITY;
			for(auto pmax:inpoints){
				  if(pmax.first == POS_INFINITY || pmax.second == POS_INFINITY) {
						w_ub=POS_INFINITY;
						break;
					}

					double ww;
					if(_eps_contract)
						ww = ( Interval(pmax.first)-eps + box3[n+3]*(Interval(pmax.second)-eps) ).ub();
					else
						ww = ( Interval(pmax.first) + box3[n+3]*(Interval(pmax.second)) ).ub();
				  if(w_ub < ww )  w_ub = ww;
			}
	 }

	box3[n+2] = Interval(NEG_INFINITY, w_ub); // w
   // cout << 	box3[n+2] <<endl;
	//the contraction is performed
	ctc.contract(box3);
	c.get<CellMOP>().a = box3[n+3].mid();
	c.get<CellMOP>().w_lb = box3[n+2].lb();

	box=box3;
	box.resize(n+2);
}

void OptimizerMOP::dominance_peeler2(IntervalVector& box, list <pair <double,double> >& inpoints){
	/*=================Dominance peeler for NDS2 and cy =========*/

	pair <double, double> firstp=inpoints.front();
	pair <double, double> lastp=inpoints.back();

	// contract c.box[n]
	if(firstp.second < box[n+1].ub())
		box[n+1] = Interval(box[n+1].lb(),firstp.second);

	// contract c.box[n+1]
	if(lastp.first < box[n].ub() )
		box[n] = Interval(box[n].lb(), lastp.first);

}


void OptimizerMOP::contract_and_bound(Cell& c, const IntervalVector& init_box) {


	//if(nds_mode == HAMBURGER){
		list<pair <double,double> > inner_segments = ndsH.non_dominated_points(c.box[n].lb(), c.box[n+1].lb());
		dominance_peeler2(c.box,inner_segments);

		if(cy_contract_var)
			cy_contract2(c,inner_segments);
		else
			ctc.contract(c.box);


	/*}else if(nds_mode==POINTS || nds_mode==SEGMENTS){
		dominance_peeler(c.box);
		//discard_generalized_monotonicty_test(c.box, init_box);

		if(cy_contract_var){
			cy_contract(c);
		}else
			ctc.contract(c.box);
	}*/

	if (c.box.is_empty()) return;



}

OptimizerMOP::Status OptimizerMOP::optimize(const IntervalVector& init_box) {

	status=SUCCESS;

	nb_cells=0;

	buffer.flush();

	/*if(nds_mode == POINTS || nds_mode == SEGMENTS){
		NDS.clear();
		//the first point
		NDS.insert(make_pair(make_pair(NEG_INFINITY,POS_INFINITY), Vector(1)));
		//the last point
		NDS.insert(make_pair(make_pair(POS_INFINITY,NEG_INFINITY), Vector(1)));
	}

	if( nds_mode == SEGMENTS)
		nds.clear();

	if(nds_mode == HAMBURGER)*/

	ndsH.clear();

	//the box in cells have the n original variables plus the two objective variables (y1 and y2)
	Cell* root=new Cell(IntervalVector(n+2));

	root->box=init_box;

	CellMOP::y1_init=eval_goal(goal1, root->box, n);
	CellMOP::y2_init=eval_goal(goal2, root->box, n);

	y1_ub.first=POS_INFINITY;
	y2_ub.second=POS_INFINITY;


	// add data required by the bisector
	bsc.add_backtrackable(*root);

	// add data required by the buffer
	buffer.add_backtrackable(*root);

	time=0;
	Timer timer;
	timer.start();



	//handle_cell(*root,init_box);
	buffer.push(root);
	if(_plot) py_Plotter::plot_add_box(root);

	int iter=0;

	try {
		/** Criterio de termino: todas los nodos filtrados*/
		while (!buffer.empty()) {
		  if(_plot) {
			  cout << "iter:" << iter << endl;
			  cout << "buffer_size:" << buffer.size() << endl;
		  }
		  iter++;


		  //if((nds_mode == POINTS || nds_mode==SEGMENTS) && _plot) {py_Plotter::offline_plot(NULL, NDS);}

		  if (trace >= 2) cout << buffer;

			Cell *c = buffer.pop();
			if(c->get<CellMOP>().ub_distance < eps) break;

			if(_plot) py_Plotter::plot_del_box(c);

			nb_cells++;
			contract_and_bound(*c, init_box);

			if (c->box.is_empty()) {
				delete c;
				continue;
			}


			//if(nds_mode==POINTS)
			//	update_NDS(c->box);

			//if(nds_mode==SEGMENTS || nds_mode == HAMBURGER)
				update_NDS2(c->box);



			pair<IntervalVector,IntervalVector>* boxes=NULL;
			bool atomic_box=false;
			try {
				boxes=new pair<IntervalVector,IntervalVector>(bsc.bisect(*c));
			}
			catch (NoBisectableVariableException& ) {
				atomic_box=true;
			}


        	double dist=0.0;
        	if(!atomic_box) dist= NDS_seg::distance(c);

        	if(dist < eps || atomic_box){
        		if(dist <0.0){
        			delete c;
        			continue;
        		}

        		if(_plot) py_Plotter::plot_add_lb(c);

        		//TODO: Puede que sea lo mismo que la linea 659 (if(dist <0.0))

						/*
        		if(nds_mode==POINTS || nds_mode==SEGMENTS){
        			map< pair <double, double>, IntervalVector >:: iterator ent1=NDS.upper_bound(make_pair(c->box[n].lb(),c->box[n+1].lb()));
        			ent1--;
        			if(ent1->first.second <= c->box[n+1].lb()){
        				delete c;
        				continue;
        			}
        		}*/

        		if(boxes) delete boxes;
				delete c; continue;

			}

      /** Improvement for avoiding big boxes when lb1 < y1_ub or lb2< y2_ub*/
		 IntervalVector left(c->box);
		  if(c->box[n].lb() < y1_ub.first && c->box[n].ub() > y1_ub.first &&
				  (c->box[n].ub()-y1_ub.first)*(c->box[n+1].ub()-y1_ub.second) <  (c->box[n].diam())*(c->box[n+1].diam()) ){

			 left[n]=Interval(c->box[n].lb(),y1_ub.first);
			 c->box[n]=Interval(y1_ub.first,c->box[n].ub());
		  }else left.set_empty();

         pair<Cell*,Cell*> new_cells;
			if(!left.is_empty())
				  new_cells=c->bisect(left,c->box);
			else{
				IntervalVector bottom(c->box);
		  	if(c->box[n+1].lb() < y2_ub.second && c->box[n+1].ub() > y2_ub.second  &&
		  			(c->box[n].ub()-y2_ub.first)*(c->box[n+1].ub()-y2_ub.second) <  (c->box[n].diam())*(c->box[n+1].diam()) ) {
					bottom[n+1]=Interval(c->box[n+1].lb(),y2_ub.second);
					c->box[n+1]=Interval(y2_ub.second,c->box[n+1].ub());
				}else bottom.set_empty();

				if(!bottom.is_empty())
				  new_cells=c->bisect(bottom,c->box);
				else
				 new_cells=c->bisect(boxes->first,boxes->second); //originally we should do only do this
			}
      /****/

    	delete boxes;
			delete c; // deletes the cell.

			buffer.push(new_cells.first);
			if(_plot) py_Plotter::plot_add_box(new_cells.first);

			buffer.push(new_cells.second);

			if(_plot) py_Plotter::plot_add_box(new_cells.second);

			if (timeout>0) timer.check(timeout); // TODO: not reentrant, JN: done
			time = timer.get_time();


		}
	}
	catch (TimeOutException& ) {
		status = TIME_OUT;
		cout << "timeout" << endl;

		//return status;
	}


	timer.stop();
	time = timer.get_time();


	// if(_plot) (nds_mode==POINTS)?
	py_Plotter::offline_plot(NULL, ndsH.NDS2);
	return status;
}


  // https://docs.google.com/document/d/1oXQhagd1dgZvkbPs34B4Nvye_GqA8lFxGZNQxqYwEgo/edit
void OptimizerMOP::hamburger(const IntervalVector& aIV, const IntervalVector& bIV) {

	IntervalVector xa=aIV;
	IntervalVector xb=bIV;

	PFunction pf(goal1, goal2, xa, xb);

	Node_t n_init (Interval(0,1), 0.0, POS_INFINITY);
	std::priority_queue<Node_t, vector<Node_t> > n;
	if(process_node(pf, n_init)) n.push(n_init);

	//cout << "INIT HAMBURGER" << endl;
	while(n.size() > 0) {
		Node_t nt = n.top();
		n.pop();
		if(nt.dist < eps) continue;

		//cout << "dist:" << nt.dist << endl;
		Node_t n1( Interval(nt.t.lb(), nt.b), 0.0, POS_INFINITY);
		if(process_node(pf, n1)){ n.push(n1);}

		Node_t n2( Interval(nt.b, nt.t.ub()), 0.0, POS_INFINITY);
		if(process_node(pf, n2)) n.push(n2);

	}
	//cout << "END HAMBURGER" << endl;

	//getchar();
}

bool OptimizerMOP::process_node(PFunction& pf, Node_t& n_t) {
    // https://docs.google.com/document/d/1oXQhagd1dgZvkbPs34B4Nvye_GqA8lFxGZNQxqYwEgo/edit

	Interval t=n_t.t;

  if(!t.is_bisectable())
		 return false;

	// get extreme points
	IntervalVector ft_lb = pf.get_point(t.lb());
	IntervalVector ft_ub = pf.get_point(t.ub());
	Interval ya1=ft_lb[0];
	Interval ya2=ft_lb[1];
	Interval yb1=ft_ub[0];
	Interval yb2=ft_ub[1];

	// add extreme points
	ndsH.addPoint(ft_lb);
	ndsH.addPoint(ft_ub);

	if(nds_mode==POINTS) return false;

	//dominated extremal points
	//if(ndsH.is_dominated(make_pair(ya1.ub(),ya2.ub())) &&  ndsH.is_dominated(make_pair(yb1.ub(),yb2.ub())))
		//return false;

  //too-close points
	if( fabs(ya1.ub() - yb1.ub()) < eps && fabs(ya2.ub() - yb2.ub()) < eps ) {
		if(_plot){
			std::vector< pair <double, double> > curve_y;
			pf.get_curve_y( curve_y );
			std::vector< pair <double, double> > curve_y_origin;
			pf.get_curve_y( curve_y_origin );
			std::vector< pair <double, double> > rectaUB;
			py_Plotter::offline_plot(NULL, ndsH.NDS2, rectaUB, curve_y_origin, curve_y);
			cout << "y1 and y2 <<< eps" << endl;
			getchar();
		}
		return false;
	}

	if(ya1.ub() > yb1.ub() || ya2.ub() < yb2.ub()) {
		Interval aux = ya1;
		ya1 = yb1;
		yb1 = aux;
		aux = ya2;
		ya2 = yb2;
		yb2 = aux;
	}

	// m ← getSlope(n.t)
	Interval m = (yb2-ya2)/(yb1-ya1);
	Interval m_horizontal = Interval(0);
	Interval m_vertical = Interval(POS_INFINITY);


	// get minf1 and minf2
	pair<double, double> c1_t1 = make_pair(POS_INFINITY,0);
	pair<double, double> c2_t2 = make_pair(POS_INFINITY,0);
	if(nds_mode==HAMBURGER){
		c1_t1 = pf.optimize(m_vertical, PFunction::MIN, POS_INFINITY, t);
		c2_t2 = pf.optimize(m_horizontal, PFunction::MIN, POS_INFINITY, t);
	}

	// check if the point dominating the curve is dominated by ndsH
	if(ndsH.is_dominated(pair<double,double>(c1_t1.first, c2_t2.first)))
			return false;

	// set bisection of inter
	list<double> v;

	//we save candidate points for bisection (t1 and t2)
	v.push_back(c1_t1.second);
	v.push_back(c2_t2.second);


	if(m.ub() < 0) {
		pair<double, double> c3_t3 = pf.optimize(m, PFunction::MAX, POS_INFINITY, t);
		pair<double, double> c4_t4;
  	if(nds_mode==HAMBURGER)
			 pair<double, double> c4_t4 = pf.optimize(m, PFunction::MIN, POS_INFINITY, t);

    ndsH.addSegment(make_pair(((ya2-c3_t3.first)/m).ub(),ya2.ub()),
										make_pair(yb1.ub(),(yb1*m+c3_t3.first).ub()));

		if(nds_mode==HAMBURGER)
			n_t.dist=ndsH.distance(c1_t1.first,c2_t2.first,m.mid(),c4_t4.first);

    v.push_back(c3_t3.second);
		v.push_back(c4_t4.second);

	}else if(nds_mode==HAMBURGER)
		n_t.dist=ndsH.distance(c1_t1.first,c2_t2.first);

	double min_dist=2.0;
	double vv;
	for(auto point:v){
		 double dist=std::abs(point-t.mid());
		 if(dist<min_dist){
			 min_dist=dist;
			 vv=point;
		 }
	}

	//if the bisection point is too close of a bound we bisect in the middle
	if( std::min(vv-t.lb(),t.ub()-vv)/t.diam() < 0.1 ) vv=t.mid();
	n_t.b=vv;


	if(_plot){
		// plot curve
		std::vector< pair <double, double> > curve_y;
		pf.get_curve_y( curve_y );
		std::vector< pair <double, double> > curve_y_origin;
		pf.get_curve_y( curve_y_origin );
		std::vector< pair <double, double> > rectaUB;
		py_Plotter::offline_plot(NULL, ndsH.NDS2, rectaUB, curve_y_origin, curve_y);
		cout << "v" << endl;
		getchar();
	}

	if(nds_mode!=HAMBURGER) return false;
	return true;
}


void OptimizerMOP::report(bool verbose) {

	if (!verbose) {
        cout << endl 	<< "time 	#nodes 		|Y|" << endl;
		cout << get_time() << " " << get_nb_cells() << " " << ndsH.size() <<  endl;
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
	/*if(nds_mode==SEGMENTS){
		cout << " number of solutions: "  << nds.size() << endl;
		for(auto ub : nds.NDS2)
			cout << "(" << ub.first.first << "," << ub.first.second << "): " << ub.second.mid() << endl;
	}
	if(nds_mode==HAMBURGER){*/
		cout << " number of solutions: "  << ndsH.size() << endl;
		for(auto ub : ndsH.NDS2)
			 cout << "(" << ub.first.first << "," << ub.first.second << "): " << ub.second.mid() << endl;
	/*}else{
		cout << " number of solutions: "  << NDS.size() << endl;
		for(auto ub : NDS)
			cout << "(" << ub.first.first << "," << ub.first.second << "): " << ub.second.mid() << endl;
	}*/

}


} // end namespace ibex
