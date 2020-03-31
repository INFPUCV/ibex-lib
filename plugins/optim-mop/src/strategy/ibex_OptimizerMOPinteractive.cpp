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

void OptimizerMOP_I::write_envelope(string output_file){
	NDS_seg LBaux;
	NDS_seg UBaux=ndsH;
	for(auto cc:cells)	LBaux.add_lb(*cc);
	for(auto cc:paused_cells) LBaux.add_lb(*cc);
	py_Plotter::offline_plot(UBaux.NDS2,  &LBaux.NDS2, output_file.c_str());
}


IntervalVector OptimizerMOP_I::load(const IntervalVector& init_box, string filename) {

	nb_cells=0;
	buffer.flush();
	ndsH.clear();

  timer.start();
	this->initbox=initbox;

	y1_ub.first=POS_INFINITY;
	y2_ub.second=POS_INFINITY;
	time=0;

	cells.clear();

	load_state_from_file(filename, init_box);

	for(auto c:cells){
		buffer.push(c);
	}

	IntervalVector focus(2);
	focus[0]=BxpMOPData::y1_init;
	focus[1]=BxpMOPData::y2_init;

	paused_cells.clear();
	current_precision = POS_INFINITY;
	return focus;

}

IntervalVector OptimizerMOP_I::load(const IntervalVector& init_box) {
	this->initbox=initbox;

	Cell* root=new Cell(IntervalVector(n+2)); //crea nodo raiz
	pre_optimize(init_box, root);
	cells.clear();
	cells.insert(root);
  istatus=READY;

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

void OptimizerMOP_I::update_refpoint(Vector& refpoint){
	this->refpoint = refpoint;
  //we move from paused_cells to cells, all the boxes with lb dominated by refpoint
	for(auto cc:paused_cells){
		IntervalVector boxy=get_boxy(cc->box,n);

		if(refpoint[0] > boxy[0].lb() && refpoint[1] > boxy[1].lb() ){
			buffer.push(cc);
			cells.insert(cc);
			paused_cells.erase(cc);
		}
	}

	if(!cells.empty()) istatus=READY;
	else istatus=STOPPED;
}

OptimizerMOP_I::IStatus OptimizerMOP_I::run(int maxiter, double eps) {
  current_precision = eps;


  for(int iter=0; iter<maxiter;){
		time = timer.get_time();

    if(buffer.empty() && paused_cells.empty()) return FINISHED;
		if(buffer.empty()) return STOPPED;

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
			 if(current_precision <= 0.0) delete c;
			 else paused_cells.insert(c);

			 while(!buffer.empty()) paused_cells.insert(buffer.pop());
			 cells.clear();
			 return STOPPED;
		}

    iter++;
		nb_cells++;

    //Contraction
		contract_and_bound(*c, initbox);
		if (c->box.is_empty()) { delete c; continue; }

    //upper_bounding
		upper_bounding(c->box);

    //Discarding by using distance and epsilon
		double dist=ndsH.distance(c);
    if(dist <= 0.0) {delete c; continue; }
		if(dist <= eps){ paused_cells.insert(c); continue; }

    /********************* bisection *************/
		pair<Cell*,Cell*> new_cells;
		try {
			new_cells=pair<Cell*,Cell*>(bsc.bisect(*c));
			delete c; // deletes the cell.
		}
		catch (NoBisectableVariableException& ) {
			throw NoBisectableVariableException();
		}

		buffer.push(new_cells.first);
		cells.insert(new_cells.first);

		buffer.push(new_cells.second);
		cells.insert(new_cells.second);
    /************* end of bisection *************/

	}


  time = timer.get_time();
	return READY;
}


} /* namespace ibex */
