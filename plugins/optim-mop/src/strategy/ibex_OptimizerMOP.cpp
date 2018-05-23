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

map< pair <double, double>, IntervalVector > OptimizerMOP::NDS;
map< pair <double, double>, IntervalVector, struct sorty2 > OptimizerMOP::NDS2;


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


bool OptimizerMOP::is_dominated(pair< double, double>& eval){
	std::map<pair<double, double>, IntervalVector>::iterator it2 = NDS.lower_bound(eval);

	//there is an equivalent point
	if(it2->first == eval) return true;
	it2--;
	//it is dominated by the previous ub point
	if(eval.second >= it2->first.second) return true;

	return false;
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
		update_NDS_pt(xa);

	}catch (LoupFinder::NotFound& ) {
		return false;
	}

	try{
		xb = finder.find(box2,box2,POS_INFINITY).first;
		update_NDS_pt(xb);
	}catch (LoupFinder::NotFound& ) {
		addPointtoNDS(make_pair(eval_goal(goal1,xa,n).ub(), eval_goal(goal2,xa,n).ub()));
		return true;
	}
	cout << "nb2" << endl;
	py_Plotter::offline_plot(NULL, NDS);
	add_upper_segment(xa, xb);


	return true;

}

bool OptimizerMOP::update_NDS_pt(IntervalVector& vec) {
		bool new_ub=false;
	//3. Se evalua el punto usando funciones objetivo (goal1 y goal2)
	pair< double, double> eval = make_pair(eval_goal(goal1,vec,n).ub(), eval_goal(goal2,vec,n).ub());

	//cout << eval.first << "," << eval.second << endl;
	//4. Insertar en mapa NDS (si es no dominada) y actualizar eliminar soluciones dominadas de NDS

	if (is_dominated(eval)) return false;

	/**** NDS correction ****/
	/*
	if(finder.ub_correction(vec.mid(), vec)){
		eval = make_pair(eval_goal(goal1,vec,n).ub(), eval_goal(goal2,vec,n).ub());
	}
	else  return false;;

	if (is_dominated(eval))  return false;;
*/

	/**** end NDS correction ****/

	bool domine=false;
	std::map<pair<double, double>, IntervalVector>::iterator it2 = NDS.lower_bound(eval);


	for(; it2!=NDS.end(); ){

		if(eval.second > it2->first.second) break;
		std::map<pair<double, double>, IntervalVector>::iterator aux = it2;
		++aux;
		if(_plot)	py_Plotter::plot_del_ub(it2->first);

		NDSy.erase(it2->first);
		NDS.erase(it2);
		it2 = aux;
		domine=true;
	}


	//the point is inserted in NDS only if its distance to the neighbor points is greater than (abs_eps/2.0)
	if(domine || std::min(it2->first.first - eval.first,  eval.second - it2->first.second) >= _min_ub_dist*eps){
		//it is not dominated and we remove the new dominated points

		if(eval.first < y1_ub.first) y1_ub=eval;
		if(eval.second < y2_ub.second) y2_ub=eval;

		NDS.insert(make_pair(eval, vec));
		NDSy.insert(make_pair(eval, vec));
		//cout << "passed" << endl;
		new_ub = true;
	}else{
		it2--;
		if( std::min(eval.first - it2->first.first,  it2->first.second - eval.second) >= _min_ub_dist*eps ){
			//it is not dominated and we remove the new dominated points

			if(eval.first < y1_ub.first) y1_ub=eval;
			if(eval.second < y2_ub.second) y2_ub=eval;

			NDS.insert(make_pair(eval, vec));
			NDSy.insert(make_pair(eval, vec));
			//cout << "passed" << endl;
			new_ub = true;
		}
	}

	if(_plot && new_ub) py_Plotter::plot_add_ub(eval);
	if(trace) {cout << eval.first  <<"," << eval.second << "(" << NDS.size() << ")" << endl;}

}

bool OptimizerMOP::update_NDS(const IntervalVector& box) {

	list<Vector> feasible_points;

	//We attempt to find two feasible points which minimize both objectives
	//and the middle point between them
	IntervalVector box2(box); box2.resize(n);
  finder.clear();

	bool new_ub=false;
	bool flag=true;


	int i=0;

	while(flag){
		i++;

		IntervalVector vec(n);

		try{
			vec = finder.find(box2,box2,POS_INFINITY).first;
		}catch (LoupFinder::NotFound& ) {
			break;
			vec = box2.mid();
			if(!finder.norm_sys.is_inner(vec)) break;
			flag=false;
		}

		cout << "nb:"  << i << endl;
		new_ub= update_NDS_pt(vec);

	}

    //cout << i << endl;


	return new_ub;

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

//TODO: Juntar con non_dominated_segments
list<pair <double,double> > OptimizerMOP::extremal_non_dominated(const IntervalVector& box){
	int n=box.size()-2;
	list <pair <double,double> > inpoints;

	pair <double, double> firstp;
	pair <double, double> lastp;

	//first point in the box (v11)
	map< pair <double, double>, IntervalVector >:: iterator ent1=NDS2.upper_bound(make_pair(box[n].lb(),NEG_INFINITY));
	while(ent1->first.second > box[n].ub()) ent1++;
	if(ent1->first.first > box[n].lb()) return inpoints; //no point in the box

	pair <double, double> v11 = ent1->first; ent1--;
    pair <double, double> v10 = ent1->first;

	//TODO: interseccion conservativa: upperPointIntersection
    firstp = pointIntersection( v10, v11, make_pair(box[n].lb(),v11.second),  make_pair(box[n].lb(),v10.second));
   	if(firstp.second <= box[n+1].lb() ){
   		cout << "error: the box is dominated" << endl;
   		exit(0);
   	}

	inpoints.push_back(firstp);

	//last point in the box (v21)
	while(ent1->first.second > box[n+1].lb() || ent1->first.first > box[n].lb())
		ent1++;

	pair <double, double> v21 = ent1->first;
	ent1++;
	pair <double, double> v20 = ent1->first;


	lastp = pointIntersection( v20, v21, make_pair(v20.first,box[n+1].lb()),  make_pair(v21.first,box[n+1].lb()));

	inpoints.push_back(lastp);

	return inpoints;
}


list<pair <double,double> > OptimizerMOP::non_dominated_segments(IntervalVector& box){
	list <pair <double,double> > inpoints;

	pair <double, double> firstp;
	pair <double, double> lastp;

	//first point in the box (v11)
	map< pair <double, double>, IntervalVector >:: iterator ent1=NDS2.upper_bound(make_pair(box[n].lb(),NEG_INFINITY));
	while(ent1->first.second > box[n].ub()) ent1++;
	if(ent1->first.first > box[n].lb()) return inpoints; //no point in the box

	pair <double, double> v11 = ent1->first; ent1--;
    pair <double, double> v10 = ent1->first;

	//TODO: interseccion conservativa: upperPointIntersection
    firstp = pointIntersection( v10, v11, make_pair(box[n].lb(),v11.second),  make_pair(box[n].lb(),v10.second));
   	if(firstp.second <= box[n+1].lb() ){
   		box.set_empty();
   		return inpoints;
   	}

	if(firstp.second > box[n+1].ub())
		firstp = pointIntersection( v10, v11, make_pair(v10.first,box[n+1].ub()),  make_pair(v11.first,box[n+1].ub()));

	//valueZ1 is a point in the bounds lb(y1) or ub(y2) of the box delimiting the NDS in the box
	inpoints.push_back(firstp);

	//last point in the box (v21)
	while(ent1->first.second > box[n+1].lb() || ent1->first.first > box[n].lb()){
		inpoints.push_back(ent1->first);
		ent1++;
	}
	pair <double, double> v21 = ent1->first;
	ent1++;
	pair <double, double> v20 = ent1->first;


	lastp = pointIntersection( v20, v21, make_pair(v20.first,box[n+1].lb()),  make_pair(v21.first,box[n+1].lb()));

	if(lastp.first > box[n].ub() )
		lastp = pointIntersection( v20, v21, make_pair(box[n].ub(),v21.second),  make_pair(box[n].ub(),v20.second));

	//valueZ2 is a point in the bounds ub(y1) or lb(y2) of the box delimiting the NDS in the box
	inpoints.push_back(lastp);

	return inpoints;
}

void OptimizerMOP::cy_contract2(Cell& c, list <pair <double,double> >& inpoints){
	IntervalVector& box = c.box;
	IntervalVector box3(box);
	box3.resize(n+4);

	pair <double, double> firstp=inpoints.front();
	pair <double, double> lastp=inpoints.back();

	//Pent of cy
    if(firstp!=lastp)
    	box3[n+3] = (lastp.first-firstp.first)/(firstp.second-lastp.second);
    else
		box3[n+3] = box[n].diam()/box[n+1].diam(); // a

     //setting w_ub with the NDS points
  	 double w_ub;

      if(_cy_upper){
		w_ub=NEG_INFINITY;
		for(auto pmax:inpoints){
			double ww;
			if(_eps_contract)
				ww = ( Interval(pmax.first)-eps + box3[n+3]*(Interval(pmax.second)-eps) ).ub();
			else
				ww = ( Interval(pmax.first) + box3[n+3]*(Interval(pmax.second)) ).ub();
			if(w_ub < ww )  w_ub = ww;
		}
  	  }

	box3[n+2] = Interval(NEG_INFINITY, w_ub); // w
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


double OptimizerMOP::distance22(const Cell* c){
	cout << 1 << endl;
	double max_dist=NEG_INFINITY;

	int n=c->box.size();

	Interval z1 = c->box[n-2];
	Interval z2 = c->box[n-1];

	double a = c->get<CellMOP>().a;
	double w_lb = c->get<CellMOP>().w_lb;

	map< pair <double, double>, IntervalVector >::iterator it = NDS2.lower_bound(make_pair(z1.lb(),POS_INFINITY)); //NDS.begin();
	//it--;


	list<pair <double,double> > inner_segments= extremal_non_dominated(c->box);

	if(inner_segments.size()>=2){
		pair <double,double> p1 = inner_segments.front();
		double dist = std::min (p1.first - z1.lb(), p1.second - z2.lb());
		//Damir's distance
		if(cy_contract_var && w_lb!=POS_INFINITY)
		  dist = std::min(dist, (Interval(p1.second)-(Interval(w_lb) - (Interval(p1.first) - Interval(p1.second)))/(Interval(a)+1.0)).ub());

		if(dist > max_dist) max_dist=dist;

		p1 = inner_segments.back();
		dist = std::min (p1.first - z1.lb(), p1.second - z2.lb());
		//Damir's distance
		if(cy_contract_var && w_lb!=POS_INFINITY)
		  dist = std::min(dist, (Interval(p1.second)-(Interval(w_lb) - (Interval(p1.first) - Interval(p1.second)))/(Interval(a)+1.0)).ub());

		if(dist > max_dist) max_dist=dist;
	}

	for(;it!=NDS2.end(); it++){
		pair <double, double> pmax= it->first;


		if(pmax.first==POS_INFINITY) pmax.first=CellMOP::y1_init.ub();
		if(pmax.second==POS_INFINITY) pmax.second=CellMOP::y2_init.ub();

		//cout << "z:" << z1.lb() << "," << z2.lb() << "  -->  ";
		//cout << "pmax:" << pmax.first << "," << pmax.second << endl;

		//el punto esta dentro de la zona de interes
		if(pmax.first >= z1.lb() && pmax.second >= z2.lb()){
			double dist = std::min (pmax.first - z1.lb(), pmax.second - z2.lb());
			//Damir's distance
			if(cy_contract_var && w_lb!=POS_INFINITY)
			  dist = std::min(dist, (Interval(pmax.second)-(Interval(w_lb) - (Interval(pmax.first) - Interval(pmax.second)))/(Interval(a)+1.0)).ub());

			if(dist > max_dist) max_dist=dist;
		}else break;
	}

	cout << "z:" << z1.lb() << "," << z2.lb() << "  -->  " << max_dist << endl;
	return max_dist;
}

void OptimizerMOP::dominance_peeler(IntervalVector& box){
	/*=================Dominance peeler ==================*/

	double z1, z2;
	pair <double, double> valueZ1;
	pair <double, double> valueZ2;

	valueZ1.first = NEG_INFINITY;
	valueZ2.second = NEG_INFINITY;


	map< pair <double, double>, IntervalVector >:: iterator ent1=NDS.upper_bound(make_pair(box[n].lb(),POS_INFINITY /*box[n+1].lb()*/));
    ent1--;

    //z1 < box[n].lb()
	z2 = ent1->first.second; // pair 2
	//the box is dominated
	if(z2 <= box[n+1].lb()){
		box.set_empty();
		return;
	}
	valueZ1 = ent1->first;

	// contract c.box[n] && c.box[n+1] with PNS points
	if(valueZ1.second > box[n+1].lb() && valueZ1.second < box[n+1].ub() ) {
		box[n+1] = Interval(box[n+1].lb(),valueZ1.second);
	}


	map< pair <double, double>, IntervalVector, sorty>:: iterator ent2=NDSy.lower_bound(make_pair(NEG_INFINITY,box[n+1].lb()));
	if(ent2==NDSy.end()) return;

	valueZ2 = ent2->first;
	//cout << valueZ1.first <<","<<valueZ1.second << " ; " << valueZ2.first <<","<<valueZ2.second << endl;

	if(valueZ2.first > box[n].lb() && valueZ2.first <= box[n].ub() ) {
		box[n] = Interval(box[n].lb(), valueZ2.first);
	}

}


void OptimizerMOP::cy_contract(Cell& c){
	IntervalVector& box = c.box;
	IntervalVector box3(box);
	box3.resize(n+4);
	map< pair <double, double>, IntervalVector >::iterator it = NDS.upper_bound(make_pair(box[n].lb(), NEG_INFINITY));
	it--;
	map< pair <double, double>, IntervalVector >::iterator it2 = NDS.lower_bound(make_pair(box[n].ub(), NEG_INFINITY));

    if(it->first.first != NEG_INFINITY && it2->first.second != NEG_INFINITY && it->first!=it2->first){
    	box3[n+3] = (it2->first.first-it->first.first)/(it->first.second-it2->first.second);
    }else if(box[n+1].diam() < POS_INFINITY)
		box3[n+3] = box[n].diam()/box[n+1].diam(); // a
	else
		box3[n+3] = 1.0;

  //setting w_ub with the NDS points
	  double w_ub=POS_INFINITY;
   //TODO::revisar
    if(_cy_upper){
	    if(NDS.size()==2)  w_ub = POS_INFINITY;
			else{
				it = NDS.lower_bound(make_pair(box[n].lb(), POS_INFINITY));
				it--;
				 w_ub=NEG_INFINITY;
				while(it!=NDS.end()){
					pair <double, double> p = it->first; it++;
					if(it==NDS.end() || p.second < box[n+1].lb()) break;
					pair <double, double> p2 = it->first;
					pair <double, double> pmax= make_pair(p2.first, p.second);

			 		if(pmax.first==POS_INFINITY || pmax.second==POS_INFINITY)
					   w_ub = POS_INFINITY;
					else{
						double ww;
						if(_eps_contract)
							ww = ( Interval(pmax.first)-eps + box3[n+3]*(Interval(pmax.second)-eps) ).ub();
						else
							ww = ( Interval(pmax.first) + box3[n+3]*(Interval(pmax.second)) ).ub();
			  		if(w_ub < ww )  w_ub = ww;
					}
				}
			}
	  }

		box3[n+2] = Interval(NEG_INFINITY, w_ub); // w
		//the contraction is performed
		ctc.contract(box3);
		c.get<CellMOP>().a = box3[n+3].mid();
		c.get<CellMOP>().w_lb = box3[n+2].lb();

		box=box3;
		box.resize(n+2);

}

void OptimizerMOP::contract_and_bound(Cell& c, const IntervalVector& init_box) {


	if(false){
		cout << 1 << endl;
		list<pair <double,double> > inner_segments; //= non_dominated_segments(c.box);

		cout << 2 << endl;
		//if(inner_segments.size()>=2) dominance_peeler2(c.box,inner_segments);


		if(cy_contract_var && inner_segments.size()>=2){
			cy_contract2(c,inner_segments);
		}else
			ctc.contract(c.box);
		cout << 4 << endl;

	}else if(nds_mode==POINTS || nds_mode==SEGMENTS){
		dominance_peeler(c.box);
		discard_generalized_monotonicty_test(c.box, init_box);

		if(cy_contract_var){
			cy_contract(c);
		}else
			ctc.contract(c.box);
	}

	if (c.box.is_empty()) return;



}

OptimizerMOP::Status OptimizerMOP::optimize(const IntervalVector& init_box) {

	status=SUCCESS;

	nb_cells=0;


	buffer.flush();

	if(nds_mode == POINTS || nds_mode == SEGMENTS){
		NDS.clear();
		//the first point
		NDS.insert(make_pair(make_pair(NEG_INFINITY,POS_INFINITY), Vector(1)));
		//the last point
		NDS.insert(make_pair(make_pair(POS_INFINITY,NEG_INFINITY), Vector(1)));
	}

	if( nds_mode == SEGMENTS){
		NDS2.clear();
		//the first point
		NDS2.insert(make_pair(make_pair(NEG_INFINITY,POS_INFINITY), Vector(1)));
		//the middle point
		NDS2.insert(make_pair(make_pair(POS_INFINITY,POS_INFINITY), Vector(1)));
		//the last point
		NDS2.insert(make_pair(make_pair(POS_INFINITY,NEG_INFINITY), Vector(1)));
	}


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
			if(_plot) py_Plotter::plot_del_box(c);



			nb_cells++;


			contract_and_bound(*c, init_box);

			if (c->box.is_empty()) {
				delete c;
				continue;
			}

			cout << "updateNDS" << endl;

			//if(nds_mode==POINTS)
			//	update_NDS(c->box);

			if(nds_mode==SEGMENTS)
				update_NDS2(c->box);
			cout << "end" << endl;



			pair<IntervalVector,IntervalVector>* boxes=NULL;
			bool atomic_box=false;
			try {
				boxes=new pair<IntervalVector,IntervalVector>(bsc.bisect(*c));
			}
			catch (NoBisectableVariableException& ) {
				atomic_box=true;
			}


        	double dist=0.0;
        	if(!atomic_box) dist= (nds_mode==POINTS || nds_mode==SEGMENTS)? distance2(c) : distance22(c);
        	cout << "dist:"  << dist << endl;

        	if(dist < eps || atomic_box){
        		if(dist <0.0){
        			delete c;
        			continue;
        		}

        		if(_plot) py_Plotter::plot_add_lb(c);

        		//TODO: Puede que sea lo mismo que la linea 659 (if(dist <0.0))

        		if(nds_mode==POINTS || nds_mode==SEGMENTS){
        			map< pair <double, double>, IntervalVector >:: iterator ent1=NDS.upper_bound(make_pair(c->box[n].lb(),c->box[n+1].lb()));
        			ent1--;
        			if(ent1->first.second <= c->box[n+1].lb()){
        				delete c;
        				continue;
        			}
        		}




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


	// if(_plot) (nds_mode==POINTS)? py_Plotter::offline_plot(NULL, NDS) : py_Plotter::offline_plot(NULL, NDS2);
	return status;
}

void OptimizerMOP::report(bool verbose) {

	if (!verbose) {
        cout << endl 	<< "time 	#nodes 		|Y|" << endl;
		cout << get_time() << " " << get_nb_cells() << " " << ((nds_mode==POINTS)? NDS.size():NDS2.size()) <<  endl;
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
	if(nds_mode==SEGMENTS){
		cout << " number of solutions: "  << NDS2.size() << endl;
		for(auto ub : NDS2)
			cout << "(" << ub.first.first << "," << ub.first.second << "): " << ub.second.mid() << endl;
	}else{
		cout << " number of solutions: "  << NDS.size() << endl;
		for(auto ub : NDS)
			cout << "(" << ub.first.first << "," << ub.first.second << "): " << ub.second.mid() << endl;
	}

}

double OptimizerMOP::distance2(const Cell* c){
	double max_dist=NEG_INFINITY;

	int n=c->box.size();

	Interval z1 = c->box[n-2];
	Interval z2 = c->box[n-1];

	double a = c->get<CellMOP>().a;
	double w_lb = c->get<CellMOP>().w_lb;

	map< pair <double, double>, IntervalVector >::iterator it = NDS.lower_bound(make_pair(z1.lb(),-NEG_INFINITY)); //NDS.begin();
	it--;

	for(;it!=NDS.end(); ){
		pair <double, double> p = it->first; it++;
		if(it==NDS.end()) break;
		pair <double, double> p2 = it->first;

		pair <double, double> pmax= make_pair(p2.first, p.second);
		if(pmax.first==POS_INFINITY) pmax.first=CellMOP::y1_init.ub();
		if(pmax.second==POS_INFINITY) pmax.second=CellMOP::y2_init.ub();

		//el punto esta dentro de la zona de interes
		if(pmax.first >= z1.lb() && pmax.second >= z2.lb()){
			double dist = std::min (pmax.first - z1.lb(), pmax.second - z2.lb());
			//Damir's distance
			if(cy_contract_var && w_lb!=POS_INFINITY)
			  dist = std::min(dist, (Interval(pmax.second)-(Interval(w_lb) - (Interval(pmax.first) - Interval(pmax.second)))/(Interval(a)+1.0)).ub());

			if(dist > max_dist) max_dist=dist;
		}else break;
	}

	return max_dist;
}

} // end namespace ibex
