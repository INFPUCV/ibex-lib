/*
 * ibex_OptimizerMOPserver.cpp
 *
 *  Created on: Sep 25, 2019
 *      Author: iaraya
 */

#include "ibex_Timer.h"
#include "ibex_Function.h"
#include "ibex_NoBisectableVariableException.h"
#include <float.h>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <set>

#include "ibex_OptimizerMOPinteractive.h"

#ifndef cdata
#define cdata ((BxpMOPData*) c->prop[BxpMOPData::id])
#endif

namespace ibex {


OptimizerMOP_I::OptimizerMOP_I(int n, const Function &f1,  const Function &f2,
		Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder,
		Mode nds_mode, Mode split_mode) :
		OptimizerMOP(n, f1, f2, ctc, bsc, buffer, finder, nds_mode, split_mode),
		refpoint(2), istatus(NONE), initbox(n) {

}

void OptimizerMOP_I::plot(){
	NDS_seg LBseg;
	for(auto cc:cells)	LBseg.add_lb(*cc);
	for(auto cc:paused_cells) LBseg.add_lb(*cc);

	//py_Plotter::offline_plot(ndsH.NDS2, &LBseg.NDS2, "output2.txt");
}

void OptimizerMOP_I::plot(list<vector<double> > &upperList, list<vector<double> > &lowerList){
	//py_Plotter::offline_plot(upperList, lowerList, "output2.txt");
}

list  < pair < bool, Vector> > OptimizerMOP_I::changes_lower_envelope(int nb_changes){
	if(changes_lower.empty()){
		NDS_seg LBseg_new;
		for(auto cc:cells)	LBseg_new.add_lb(*cc);
		for(auto cc:paused_cells) LBseg_new.add_lb(*cc);
		list  < pair< Vector, NDS_data> > add_changes;
		list  <  pair< Vector, NDS_data> > rem_changes;
		std::set_difference( LBseg_new.NDS2.begin(), LBseg_new.NDS2.end(),
	    LBseg.begin(), LBseg.end(),
	    std::back_inserter(add_changes), sorty2() );

		std::set_difference( LBseg.begin(), LBseg.end(),
		    LBseg_new.NDS2.begin(), LBseg_new.NDS2.end(),
		    std::back_inserter(rem_changes), sorty2() );

		for (auto ch : rem_changes)
					changes_lower.push_back(make_pair(false, ch.first));

		for (auto ch : add_changes)
			changes_lower.push_back(make_pair(true, ch.first));

	  LBseg = LBseg_new.NDS2;
	}

  list < pair < bool, Vector> > ch ;
	if(nb_changes==-1){
		ch=changes_lower;
		changes_lower.clear();
	}else{
		for(int i=0; i<nb_changes && !changes_lower.empty() ;i++){
			ch.push_back(changes_lower.front());
			changes_lower.pop_front();
    }
  }

	return ch;
}

void OptimizerMOP_I::write_envelope(string output_file){
	NDS_seg LBaux;
	NDS_seg UBaux=ndsH;
	for(auto cc:cells)	LBaux.add_lb(*cc);
	for(auto cc:paused_cells) LBaux.add_lb(*cc);
	//py_Plotter::offline_plot(UBaux.NDS2,  &LBaux.NDS2, output_file.c_str());
}

IntervalVector OptimizerMOP_I::load(const IntervalVector& init_box) {
	this->initbox=initbox;

	Cell* root=new Cell(IntervalVector(n+2)); //crea nodo raiz
	pre_optimize(init_box, root);
	cells.clear();
	cells.insert(root);
  istatus=PAUSED;

	IntervalVector focus(2);
	focus[0]=BxpMOPData::y1_init;
	focus[1]=BxpMOPData::y2_init;
	refpoint[0] = POS_INFINITY;
	refpoint[1] = POS_INFINITY;

	paused_cells.clear();
	current_precision = POS_INFINITY;
	return focus;
	//return _optimize(init_box, focus);
}


void OptimizerMOP_I::update_refpoint(Vector refpoint, double eps){
	cout << "cells:" << cells.size() << endl;
	cout << "paused_cells:" << paused_cells.size() << endl;

	this->refpoint = refpoint;

	set<Cell*> paused_cells_to_erase;	

	for(auto c:paused_cells){

		IntervalVector boxy=get_boxy(c->box,n);
		if(refpoint[0] > boxy[0].lb() && refpoint[1] > boxy[1].lb() &&  cdata->ub_distance > eps){

			buffer.push(c);
			cells.insert(c);
			paused_cells_to_erase.insert(c);
			//paused_cells.erase(c);
		}
	}

	for(auto c:paused_cells_to_erase){
		paused_cells.erase(c);
	}

	if(!cells.empty()) istatus=PAUSED;
	else istatus=REACHED_PRECISION;
}

OptimizerMOP_I::IStatus OptimizerMOP_I::run(int maxiter, double eps) {
	

	if(current_precision < eps) return REACHED_PRECISION;
	else update_refpoint(refpoint, eps);


  current_precision = eps;


  for(int iter=0; iter<maxiter;){
	  
		time = timer.get_time();

    	if(buffer.empty() && paused_cells.empty()) return FINISHED;
		if(buffer.empty()) return REACHED_PRECISION;

		Cell* c=NULL;

		c = buffer.top();
		buffer.pop();
		cells.erase(c);

    	//we verify that the box_lb dominates the ref_point, otherwise it is paused
		IntervalVector boxy = get_boxy(c->box,n);
		
		

		if(refpoint[0] < boxy[0].lb() || refpoint[1] < boxy[1].lb() ){
			paused_cells.insert(c);
			continue;
		}

    	//we verify that the boxes in the search tree are not eps-dominated
    	current_precision=cdata->ub_distance;
		
		if(current_precision <= eps){
			 if(current_precision < 0.0) delete c;
			 else paused_cells.insert(c);

			 while(!buffer.empty()) paused_cells.insert(buffer.pop());
			 cells.clear();
			 return REACHED_PRECISION;
		}

    iter++;
		nb_cells++;

    //Contraction
		contract_and_bound(*c, initbox);
		if (c->box.is_empty()) { delete c; continue; }

		//Aquí se debería imprimir (para visualización):
		//c->id, c->parent->id, c->dist

    //upper_bounding
		upper_bounding(c->box);

    //Discarding by using distance and epsilon
		double dist=ndsH.distance(c);
    if(dist < 0.0) {delete c; continue; }
		if(dist <= eps){ paused_cells.insert(c); continue; }

    /********************* bisection *************/
		pair<Cell*,Cell*> new_cells;
		try {
			new_cells=pair<Cell*,Cell*>(bsc.bisect(*c));
			delete c; // deletes the cell. Se mantiene par visualización
		}
		catch (NoBisectableVariableException& ) {
			paused_cells.insert(c);
			cout << "blocked:" << c->box <<endl;
			continue;
			throw NoBisectableVariableException();
		}

		buffer.push(new_cells.first);
		cells.insert(new_cells.first);

		buffer.push(new_cells.second);
		cells.insert(new_cells.second);



    /************* end of bisection *************/

	}

  

  time = timer.get_time();
  return PAUSED;
}


} /* namespace ibex */
