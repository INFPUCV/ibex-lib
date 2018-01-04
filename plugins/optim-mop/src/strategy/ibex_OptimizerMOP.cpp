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
#include "ibex_CellSet.h"

using namespace std;

namespace ibex {

const double OptimizerMOP::default_eps=0.01;

bool OptimizerMOP::_plot = false;
int OptimizerMOP::_nb_ub_sols = 3;
double OptimizerMOP::_min_ub_dist = 1e-7;
bool OptimizerMOP::_cy_upper =false;
bool OptimizerMOP::_hv =false;
bool OptimizerMOP::cy_contract_var = false;
bool OptimizerMOP::_eps_contract = false;

map< pair <double, double>, IntervalVector > OptimizerMOP::UB;

OptimizerMOP::OptimizerMOP(int n, const Array<NumConstraint>& ctrs, const Function &f1,  const Function &f2,
		Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder,double eps) : n(n),
                				ctc(ctc), bsc(bsc), buffer(buffer), ctrs(ctrs), goal1(f1), goal2(f2),
								finder(finder), trace(false), timeout(-1), status(SUCCESS),
                				time(0), nb_cells(0), nb_sols(0), eps(eps),
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
		if(eval.second >= it2->first.second) continue;


		/**** UB correction ****/
		if(finder.ub_correction(vec.mid(), vec)){
			eval = make_pair(eval_goal(goal1,vec).ub(), eval_goal(goal2,vec).ub());
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
			if(_plot){
				OptimizerMOP::plot_del_ub(it2->first);
			}
			UBy.erase(it2->first);
			UB.erase(it2);
			it2 = aux;
			domine=true;
		}


		//the point is inserted in UB only if its distance to the neighbor points is greater than (abs_eps/2.0)
		if(domine || std::min(it2->first.first - eval.first,  eval.second - it2->first.second) >= _min_ub_dist*eps){
			//it is not dominated and we remove the new dominated points

			if(eval.first < y1_ub.first) y1_ub=eval;
			if(eval.second < y2_ub.second) y2_ub=eval;

			UB.insert(make_pair(eval, vec));
			if(_plot){
				OptimizerMOP::plot_add_ub(eval);
			}
			UBy.insert(make_pair(eval, vec));
			new_ub = true;
		}else{
			it2--;
			if( std::min(eval.first - it2->first.first,  it2->first.second - eval.second) >= _min_ub_dist*eps ){
				//it is not dominated and we remove the new dominated points

				if(eval.first < y1_ub.first) y1_ub=eval;
				if(eval.second < y2_ub.second) y2_ub=eval;
				if(_plot){
					OptimizerMOP::plot_add_ub(eval);
				}
				UB.insert(make_pair(eval, vec));
				UBy.insert(make_pair(eval, vec));
				new_ub = true;
			}

		}
		if(trace) {cout << eval.first  <<"," << eval.second << "(" << UB.size() << ")" << endl;}


		//UB.insert(make_pair(eval, vec));
		//new_ub = true;


	}

	//5. Si el mapa UB fue modificado retornar true, si no false
	return new_ub;

}


void OptimizerMOP::dominance_peeler(IntervalVector& box){
	/*=================Dominance peeler ==================*/

	double z1, z2;
	pair <double, double> valueZ1;
	pair <double, double> valueZ2;

	valueZ1.first = NEG_INFINITY;
	valueZ2.second = NEG_INFINITY;


	map< pair <double, double>, IntervalVector >:: iterator ent1=UB.upper_bound(make_pair(box[n].lb(),POS_INFINITY /*box[n+1].lb()*/));
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


	map< pair <double, double>, IntervalVector, sorty>:: iterator ent2=UBy.lower_bound(make_pair(NEG_INFINITY,box[n+1].lb()));
	if(ent2==UBy.end()) return;

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
	map< pair <double, double>, IntervalVector >::iterator it = UB.upper_bound(make_pair(box[n].lb(), NEG_INFINITY));
	it--;
	map< pair <double, double>, IntervalVector >::iterator it2 = UB.lower_bound(make_pair(box[n].ub(), NEG_INFINITY));

    if(it->first.first != NEG_INFINITY && it2->first.second != NEG_INFINITY && it->first!=it2->first){
    	box3[n+3] = (it2->first.first-it->first.first)/(it->first.second-it2->first.second);
    }else if(box[n+1].diam() < POS_INFINITY)
		box3[n+3] = box[n].diam()/box[n+1].diam(); // a
	else
		box3[n+3] = 1.0;

  //setting w_ub with the UB points
	  double w_ub=POS_INFINITY;
   //TODO::revisar
    if(_cy_upper){
	    if(UB.size()==2)  w_ub = POS_INFINITY;
			else{
				it = UB.lower_bound(make_pair(box[n].lb(), POS_INFINITY));
				it--;
				 w_ub=NEG_INFINITY;
				while(it!=UB.end()){
					pair <double, double> p = it->first; it++;
					if(it==UB.end() || p.second < box[n+1].lb()) break;
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
		c.get<CellBS>().a = box3[n+3].mid();
		c.get<CellBS>().w_lb = box3[n+2].lb();

		box=box3;
		box.resize(n+2);

}

void OptimizerMOP::contract_and_bound(Cell& c, const IntervalVector& init_box) {


	/*map< pair <double, double>, IntervalVector >:: iterator ent1=UB.upper_bound(make_pair(c.box[n].lb(),c.box[n+1].lb()));
	ent1--;
	if(ent1->first.second <= c.box[n+1].lb()){
		c.box.set_empty();
		return;
	}*/

  if(trace && _plot ) { plot(&c);  cout << "init" << endl; getchar(); }

	dominance_peeler(c.box);
	if(!c.box.is_empty() && trace && _plot ) { plot(&c);  cout << "dominance_peeler" << endl;  getchar();}
	discard_generalized_monotonicty_test(c.box, init_box);

	if(cy_contract_var){
		cy_contract(c);
	}else
		ctc.contract(c.box);

	if(!c.box.is_empty() && trace && _plot ) { plot(&c);  cout << "contreacted" << endl; getchar(); }

	if (c.box.is_empty()) return;



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

	CellBS::y1_init=eval_goal(goal1, root->box);
	CellBS::y2_init=eval_goal(goal2, root->box);

	y1_ub.first=POS_INFINITY;
	y2_ub.second=POS_INFINITY;

	y1_max=CellBS::y1_init.ub();
	y2_max=CellBS::y2_init.ub();


	// add data required by the bisector
	bsc.add_backtrackable(*root);

	// add data required by the buffer
	buffer.add_backtrackable(*root);

	time=0;
	Timer timer;
	timer.start();

	//handle_cell(*root,init_box);
	buffer.push(root);
	if(_plot){
		buffer_cells.insert(root);
		OptimizerMOP::plot_add_box(root);
	}

	top_dist=POS_INFINITY;

	try {
		/** Criterio de termino: todas los nodos filtrados*/
		while (!buffer.empty()) {

		  if (trace >= 2) cout << buffer;


			Cell *c = buffer.pop();
			if(_plot){
				OptimizerMOP::plot_del_box(c);
			}
			top_dist=c->get<CellBS>().ub_distance;

			if(_plot) buffer_cells.erase(c);
			nb_cells++;


			contract_and_bound(*c, init_box);

			if (c->box.is_empty()) {
				if(_plot && buffer.empty())
					plot(c);
				delete c;
				continue;
			}

			bool loup_ch=update_UB(c->box, _nb_ub_sols);




			pair<IntervalVector,IntervalVector>* boxes=NULL;
			bool atomic_box=false;
			try {
				boxes=new pair<IntervalVector,IntervalVector>(bsc.bisect(*c));
			}
			catch (NoBisectableVariableException& ) {
				atomic_box=true;
			}


        	double dist=0.0;
        	if(!atomic_box /*&& eps>0.0*/) dist=distance2(c);

        	if(trace) cout << "distance:" << dist << endl;
        	if(trace && loup_ch && _plot  && dist>=0) { cout << c->box[n] << ";" << c->box[n+1 ] << endl;   ; plot(c);  getchar(); }

        	if(dist < eps || atomic_box){
        		if(dist <0.0){
        			delete c;
        			continue;
        		}

        		if(_plot){
        			OptimizerMOP::plot_add_lb(c);
        		}

        		map< pair <double, double>, IntervalVector >:: iterator ent1=UB.upper_bound(make_pair(c->box[n].lb(),c->box[n+1].lb()));
        		ent1--;
        		if(ent1->first.second <= c->box[n+1].lb()){
        			if(_plot && buffer.empty()) plot(c);
        			delete c;
        			continue;
        		}


				if(_plot && buffer.empty() )
					plot(c);

				nb_sols++;

				if(_hv && dist>=0) updateLB(c, dist);

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
			if(_plot){
				OptimizerMOP::plot_add_box(new_cells.first);
				buffer_cells.insert(new_cells.first);
			}
			buffer.push(new_cells.second);

			if(_plot){
				OptimizerMOP::plot_add_box(new_cells.second);
				buffer_cells.insert(new_cells.second);
			}


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
	if(_hv){
		cout << "lb-hypervolume:" << compute_lb_hypervolume() << endl;
		cout << "ub-hypervolume:" << compute_ub_hypervolume() << endl;

		cout << "diff-hypervolume:" << (compute_lb_hypervolume()-compute_ub_hypervolume())/(CellBS::y1_init.diam()*CellBS::y2_init.diam())<< endl;
	}

	status=SUCCESS;

	return status;
}

void OptimizerMOP::plot_add_ub(pair<double, double> eval){
	std::cout << "add ub: {\"pts\": (" << eval.first << "," << eval.second << ")}" << endl;
}

void OptimizerMOP::plot_del_ub(pair<double, double> eval){
	std::cout << "del ub: {\"pts\": (" << eval.first << "," << eval.second << ")}" << endl;
}

void OptimizerMOP::plot_add_lb(Cell* c){
	std::cout << "add lb: {\"id\":" << c->get<CellBS>().id;
	std::cout << ", 'pts':(" << c->box[n].lb() << "," <<  c->box[n+1].lb() << ")";
	std::cout << "}" << endl;

}

void OptimizerMOP::plot_add_box(Cell* c){
	std::cout << "add: {\"id\":" << c->get<CellBS>().id;
	std::cout << ", 'pts':(" << c->box[n].lb() << "," <<  c->box[n+1].lb() << "),";
	std::cout << "'diam_x': " <<  c->box[n].diam() << ",'diam_y': " << c->box[n+1].diam();
	std::cout << ", 'pA':(" << c->box[n].lb() <<"," <<  (((c)->get<CellBS>().w_lb-c->box[n].lb())/(c)->get<CellBS>().a)   << "),";
	std::cout << "'pB':(" << (c->get<CellBS>().w_lb-c->get<CellBS>().a*c->box[n+1].lb()) <<"," <<  c->box[n+1].lb()  << ")";
	std::cout << "}" << endl;
}

void OptimizerMOP::plot_del_box(Cell* c){
	std::cout << "del: {\"id\":" << c->get<CellBS>().id;
	std::cout << "}" << endl;
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
		//cout << get_status() << endl;
		//cout << "#nb_sols:" << nb_sols << endl;
		//list<  IntervalVector > :: iterator sol=Sout.begin();


		//for(;sol!=Sout.end();sol++){
			//cout << "(" << (*sol)[n].lb() << "," << (*sol)[n+1].lb() << ")" << endl;
		//}

		//cout << "|Y|:" << UB.size() << endl;
		map< pair <double, double>, IntervalVector > :: iterator ub=UB.begin();
/*
		for(;ub!=UB.end();ub++){
			cout << ub->second << endl;
			cout << "(" << ub->first.first << "," << ub->first.second << ")" << endl;
		}
*/
        cout << endl 	<< "time #nodes |Y| #sols" << endl;
		cout << endl 	<< get_time() << " " << top_dist << " " << get_nb_cells() << " " << UB.size() <<  endl;
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

double OptimizerMOP::distance2(const Cell* c){
	double max_dist=NEG_INFINITY;

	int n=c->box.size();

	Interval z1 = c->box[n-2];
	Interval z2 = c->box[n-1];

	double a = c->get<CellBS>().a;
	double w_lb = c->get<CellBS>().w_lb;

	map< pair <double, double>, IntervalVector >::iterator it = UB.lower_bound(make_pair(z1.lb(),-NEG_INFINITY)); //UB.begin();
	it--;

	for(;it!=UB.end(); ){
		pair <double, double> p = it->first; it++;
		if(it==UB.end()) break;
		pair <double, double> p2 = it->first;

		pair <double, double> pmax= make_pair(p2.first, p.second);
		if(pmax.first==POS_INFINITY) pmax.first=CellBS::y1_init.ub();
		if(pmax.second==POS_INFINITY) pmax.second=CellBS::y2_init.ub();

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
