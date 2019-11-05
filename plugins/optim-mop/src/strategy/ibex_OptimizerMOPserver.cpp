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

#include "ibex_OptimizerMOPserver.h"

#ifndef cdata
#define cdata ((BxpMOPData*) c->prop[BxpMOPData::id])
#endif

namespace ibex {

string OptimizerMOP_S::instructions_file="";
string OptimizerMOP_S::output_file="";

OptimizerMOP_S::OptimizerMOP_S(int n, const Function &f1,  const Function &f2,
		Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder,
		Mode nds_mode, Mode split_mode, double eps, double rel_eps) :
		OptimizerMOP(n, f1, f2, ctc, bsc, buffer, finder, nds_mode, split_mode, eps, rel_eps), sstatus(STAND_BY_SEARCH){

}


void OptimizerMOP_S::zoom(bool out, set<Cell*>& cells, set<Cell*>& paused_cells, IntervalVector& focus, ifstream& myfile){
	sstatus=FOCUS_SEARCH;
	double y1_lb,y1_ub,y2_lb,y2_ub;
	myfile >> y1_lb >> y1_ub;
	myfile >> y2_lb >> y2_ub;
	focus[0] = Interval(y1_lb,y1_ub);
	focus[1] = Interval(y2_lb,y2_ub);

	if(out){
		focus[0]=BxpMOPData::y1_init;
		focus[1]=BxpMOPData::y2_init;
		update_focus(cells, paused_cells, focus);
	}

	for(auto cc:paused_cells){
		buffer.push(cc);
		cells.insert(cc);
	}
	paused_cells.clear();
	max_dist_eps=NEG_INFINITY;
}

void OptimizerMOP_S::get_solution(ifstream& myfile){

	string output_file;
	double y1,y2;
	string y1_str, y2_str;
	myfile >> output_file;
	myfile >> y1_str >> y2_str;
	y1 = stod(y1_str);
	y2 = stod(y2_str);


	Vector y(2); y[0]=y1; y[1]=y2;

	pair<Vector, NDS_data> data = ndsH.get(y);
	cout << "y:" << data.first << endl;
	if(data.second.x1) cout << *data.second.x1 << endl;
	if(data.second.x2) cout << *data.second.x2 << endl;

	Vector* v=NULL;
	Vector realy(2);
	if(!data.second.x2 || data.second.x1 == data.second.x2){
		realy[0]=eval_goal(goal1, *data.second.x1, data.second.x1->size()).ub();
		realy[1]=eval_goal(goal2, *data.second.x1, data.second.x1->size()).ub();
		if( realy[0] < y[0] + eps && realy[1] < y[1] + eps)
		  v=new Vector(*data.second.x1);
	}else{
		PFunction pf(goal1, goal2, *data.second.x1, *data.second.x2);
		v=pf.find_feasible(y, 1e-8);
	}

	ofstream output, output_tmp;
	output.open(output_file,ios_base::app);
	output_tmp.open("output.tmp");
	if(v){
		realy[0]=eval_goal(goal1, *v, v->size()).ub();
		realy[1]=eval_goal(goal2, *v, v->size()).ub();
    output << y1_str << " " << y2_str << endl;
		output << *v << endl;
		output << realy << endl;
		output_tmp << *v << endl;
		output_tmp << realy << endl;
		delete v;
	}else {
		output << "not found" << endl;
		output_tmp << "not found" << endl;
	}
	output.close();
	output_tmp.close();


}



void OptimizerMOP_S::read_instructions(set<Cell*>& cells, set<Cell*>& paused_cells, IntervalVector& focus){


	//se lee el archivo de instrucciones y se elimina
	string line; ifstream myfile;
	myfile.open(instructions_file);
	if (myfile.is_open()){
		string instruction;
		while(!myfile.eof()){
			myfile >> instruction ;
			if(instruction=="zoom_in" || instruction=="zoom_out"){
				zoom(instruction == "zoom_out", cells, paused_cells, focus, myfile);

			}else if(instruction=="upper_envelope"){
				get_solution(myfile);

			}else if(instruction=="pause"){
				 if(sstatus==SEARCH) sstatus=STAND_BY_SEARCH;
				 else if(sstatus==FOCUS_SEARCH) sstatus=STAND_BY_FOCUS;
			}else if(instruction=="continue"){
				 if(sstatus==STAND_BY_SEARCH) sstatus=SEARCH;
				 else if(sstatus==STAND_BY_FOCUS) sstatus=FOCUS_SEARCH;
			}else if(instruction=="finish"){
				 sstatus=FINISHED;
			}
		}

		myfile.close();
		rename(instructions_file.c_str(), (instructions_file+".old").c_str());
	}

}

void OptimizerMOP_S::write_envelope(set<Cell*>& cells, set<Cell*>& paused_cells, IntervalVector& focus){
	//escritura de archivos
	//dormir 1 segundo y lectura de instrucciones
	cout << "escritura de archivo" << endl;

	NDS_seg LBaux;
	NDS_seg UBaux=ndsH;

	update_focus(cells, paused_cells, focus);

	for(auto cc:cells)	LBaux.add_lb(*cc);
	for(auto cc:paused_cells) LBaux.add_lb(*cc);

	//se escribe el archivo de salida
	IntervalVector focus2(2);
	focus2[0]=BxpMOPData::y1_init;
	focus2[1]=BxpMOPData::y2_init;
	update_focus(cells, paused_cells, focus2);

	py_Plotter::offline_plot(UBaux.NDS2,  &LBaux.NDS2, output_file.c_str(), &focus2);
}

void OptimizerMOP_S::write_status(double rel_prec){
	ofstream output;
	output.open( (output_file+".state").c_str());
	switch(sstatus){
		case STAND_BY_SEARCH:
		case STAND_BY_FOCUS: output << "STAND_BY" ; break;
		case REACHED_PRECISION: output << "REACHED_PRECISION" ; break;
		case SEARCH: output << "SEARCH" ; break;
		case FOCUS_SEARCH: output << "FOCUS_SEARCH" ; break;
		case FINISHED: output << "FINISHED" ; break;
	}
	output << "," << rel_prec*100 << endl;
	output.close();
}

void OptimizerMOP_S::update_focus(set<Cell*>& cells, set<Cell*>& paused_cells, IntervalVector& focus){

	IntervalVector new_focus(2);
	new_focus.set_empty();

	for(auto cc:cells){
		IntervalVector boxy=get_boxy(cc->box,n);
		if(new_focus.is_empty())
			new_focus=boxy;
		else new_focus|=boxy;
	}

	for(auto cc:paused_cells){
		IntervalVector boxy=get_boxy(cc->box,n);
		if(new_focus.is_empty())
			new_focus=boxy;
		else new_focus|=boxy;
	}

	focus&=new_focus;

}


OptimizerMOP_S::Status OptimizerMOP_S::optimize(const IntervalVector& init_box) {

	sstatus=SEARCH;

	Cell* root=new Cell(IntervalVector(n+2));
	pre_optimize(init_box, root);

	Timer timer;
	timer.start();

	Timer timer_stand_by;

	set<Cell*> cells;
	set<Cell*> paused_cells;
	cells.insert(root);

	IntervalVector focus(2);
	focus[0]=BxpMOPData::y1_init;
	focus[1]=BxpMOPData::y2_init;

	int iter = 0;
	try {
		bool server_pause=false;
		while (!buffer.empty() || !paused_cells.empty()) {
			iter++;
			Cell *c = buffer.top();

			if(iter%10==0) server_pause=true;
			while(buffer.empty() || sstatus==REACHED_PRECISION || sstatus==STAND_BY_FOCUS || sstatus==STAND_BY_SEARCH || server_pause){
        if (sstatus == SEARCH || sstatus== FOCUS_SEARCH) timer_stand_by.restart();
				if(server_pause) {
			    	cout << "buffer size:" << buffer.size() << endl;
			    	cout << "eps:" << eps << endl;
						write_envelope(cells, paused_cells, focus);
				}
				sleep(4);
				read_instructions(cells, paused_cells, focus);
				write_status(std::max(max_dist_eps,cdata->ub_distance));

				if(sstatus == FINISHED || (buffer.empty() && paused_cells.empty()) || timer_stand_by.get_time()>1000 ){
					 sstatus = FINISHED;
					 write_status(std::max(max_dist_eps,cdata->ub_distance));
					 exit(0);
				}
				if(buffer.empty()) sstatus = REACHED_PRECISION;
				server_pause=false;
			}


			if(rel_eps>0.0)	eps=std::max(focus[0].diam(),focus[1].diam())*rel_eps;

			buffer.pop();
			cells.erase(c);


			//server stuff
			IntervalVector boxy(2); boxy[0]=c->box[n]; boxy[1]=c->box[n+1];
			list< Vector > inner_segments = ndsH.non_dominated_points(focus.lb());
			dominance_peeler2(focus,inner_segments);

			if(focus[0].ub()<boxy[0].lb() || focus[1].ub()<boxy[1].lb() ){
				paused_cells.insert(c);
				continue;
			}




			if(cdata->ub_distance <= eps){
				 //server stuff: the cell is paused
				 while(!buffer.empty()){
				   paused_cells.insert(buffer.pop());
				   if(max_dist_eps<cdata->ub_distance) max_dist_eps=cdata->ub_distance;
				 }
				 cells.clear();
				 continue;
				 //break;
			}



			nb_cells++;
			contract_and_bound(*c, init_box);

			if (c->box.is_empty()) {
				delete c;
				continue;
			}

			upper_bounding(c->box);

			pair<Cell*,Cell*> new_cells;
			bool atomic_box=false;
			try {
				new_cells=pair<Cell*,Cell*>(bsc.bisect(*c));
			}
			catch (NoBisectableVariableException& ) {
				atomic_box=true;
			}

			double dist=0.0;
			if(!atomic_box) dist=ndsH.distance(c);

			//if the box is dominated it is removed, otherwise it is paused
			if(dist < eps || atomic_box){

				if(new_cells.first){
					delete new_cells.first;
					delete new_cells.second;
				}

				if(dist>=0.0){
					paused_cells.insert(c);
				}else delete c;

				continue;
			}

			delete c; // deletes the cell.

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

		cout << "timeout" << endl;
	}

	timer.stop();
	time = timer.get_time();

	return status;
}


} /* namespace ibex */
