//                                  I B E X                                   
// File        : ibex_OptimizerANN.cpp
// Author      : Gilles Chabert, Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : May 14, 2012
// Last Update : December 24, 2012
//============================================================================

#include "ibex_OptimizerANN.h"
#include "ibex_Timer.h"
#include "ibex_Function.h"
#include "ibex_NoBisectableVariableException.h"
#include "ibex_Backtrackable.h"
#include "ibex_OptimData.h"

#include "ibex_CtcCompo.h"

#include "ibex_IntervalVector.h"

#include "ibex_CellData.h"

#include "ibex_ANN.h"

#include <map>

#include <float.h>
#include <stdlib.h>
#include <iomanip>

#include "ibex_BitSet.h"
#include "ibex_CtcFixPoint.h"
#include "ibex_CtcCompo.h"
#include "ibex_CtcPolytopeHull.h"

using namespace std;

namespace ibex {

const double OptimizerANN::default_eps_x = 0;
const double OptimizerANN::default_rel_eps_f = 1e-03;
const double OptimizerANN::default_abs_eps_f = 1e-07;
const double OptimizerANN::default_threshold = 0.5;
const int OptimizerANN::default_trainingdata = 2000;

void OptimizerANN::write_ext_box(const IntervalVector& box, IntervalVector& ext_box) {
	int i2=0;
	for (int i=0; i<n; i++,i2++) {
		if (i2==goal_var) i2++; // skip goal variable
		ext_box[i2]=box[i];
	}
}

void OptimizerANN::read_ext_box(const IntervalVector& ext_box, IntervalVector& box) {
	int i2=0;
	for (int i=0; i<n; i++,i2++) {
		if (i2==goal_var) i2++; // skip goal variable
		box[i]=ext_box[i2];
	}
}



OptimizerANN::OptimizerANN(int n, CtcCompo& ctc, Bsc& bsc, LoupFinder& finder,
		CellBufferOptim& buffer,
		int goal_var, double eps_x, double rel_eps_f, double abs_eps_f, double threshold, double trainingdata) :
                				n(n), goal_var(goal_var),
                				ctc(ctc), bsc(bsc), loup_finder(finder), buffer(buffer),
                				eps_x(eps_x), rel_eps_f(rel_eps_f), abs_eps_f(abs_eps_f),
                				trace(0), timeout(-1),
                				status(SUCCESS),
                				//kkt(normalized_user_sys),
						uplo(NEG_INFINITY), uplo_of_epsboxes(POS_INFINITY), loup(POS_INFINITY),
                				loup_point(n), initial_loup(POS_INFINITY), loup_changed(false),
                                                time(0), nb_cells(0), ann("trainingData.txt", n+1), threshold(threshold), trainingdata(trainingdata) {
	if (trace) cout.precision(12);

	/*
	vector<double> inputVals;
	cout << "training " << endl;
	// training ANN
	for(int i=0; i < 7000; i++) {
		ann.trainingNeuron(inputVals, inputVals);

	}
	cout << "testing " << endl;
	// testing ANN
	for(int i=0; i < 10; i++) {
		ann.testingNeuron(inputVals);
	}
	*/

}


OptimizerANN::~OptimizerANN() {

}

// compute the value ymax (decreasing the loup with the precision)
// the heap and the current box are contracted with y <= ymax
double OptimizerANN::compute_ymax() {
	double ymax = loup - rel_eps_f*fabs(loup);
	if (loup - abs_eps_f < ymax)
		ymax = loup - abs_eps_f;
	return ymax;
}

bool OptimizerANN::update_loup(const IntervalVector& box) {

	try {
		pair<IntervalVector,double> p=loup_finder.find(box,loup_point,loup);
		loup_point = p.first;
		loup = p.second;

		if (trace) {
			cout << "                    ";
			cout << "\033[32m loup= " << loup << "\033[0m" << endl;
//			cout << " loup point=";
//			if (loup_finder.rigorous())
//				cout << loup_point << endl;
//			else
//				cout << loup_point.lb() << endl;
		}
		return true;

	} catch(LoupFinder::NotFound&) {
		return false;
	}
}

//bool OptimizerANN::update_entailed_ctr(const IntervalVector& box) {
//	for (int j=0; j<m; j++) {
//		if (entailed->normalized(j)) {
//			continue;
//		}
//		Interval y=sys.ctrs[j].f.eval(box);
//		if (y.lb()>0) return false;
//		else if (y.ub()<=0) {
//			entailed->set_normalized_entailed(j);
//		}
//	}
//	return true;
//}

void OptimizerANN::update_uplo() {
	double new_uplo=POS_INFINITY;

	if (! buffer.empty()) {
		new_uplo= buffer.minimum();
		if (new_uplo > loup) {
			cout << " loup = " << loup << " new_uplo=" << new_uplo << endl;
			ibex_error("OptimizerANN: new_uplo>loup (please report bug)");
		}
		if (new_uplo < uplo) {
			cout << "uplo= " << uplo << " new_uplo=" << new_uplo << endl;
			ibex_error("OptimizerANN: new_uplo<uplo (please report bug)");
		}

		// uplo <- max(uplo, min(new_uplo, uplo_of_epsboxes))
		if (new_uplo < uplo_of_epsboxes) {
			if (new_uplo > uplo) {
				uplo = new_uplo;

				if (trace)
					cout << "\033[33m uplo= " << uplo << "\033[0m" << endl;
			}
		}
		else uplo = uplo_of_epsboxes;
	}
	else if (buffer.empty() && loup != POS_INFINITY) {
		// empty buffer : new uplo is set to ymax (loup - precision) if a loup has been found
		new_uplo=compute_ymax(); // not new_uplo=loup, because constraint y <= ymax was enforced
		//    cout << " new uplo buffer empty " << new_uplo << " uplo " << uplo << endl;

		double m = (new_uplo < uplo_of_epsboxes) ? new_uplo :  uplo_of_epsboxes;
		if (uplo < m) uplo = m; // warning: hides the field "m" of the class
		// note: we always have uplo <= uplo_of_epsboxes but we may have uplo > new_uplo, because
		// ymax is strictly lower than the loup.
	}

}

void OptimizerANN::update_uplo_of_epsboxes(double ymin) {

	// the current box cannot be bisected.  ymin is a lower bound of the objective on this box
	// uplo of epsboxes can only go down, but not under uplo : it is an upperbound for uplo,
	//that indicates a lowerbound for the objective in all the small boxes
	// found by the precision criterion
	assert (uplo_of_epsboxes >= uplo);
	assert(ymin >= uplo);
	if (uplo_of_epsboxes > ymin) {
		uplo_of_epsboxes = ymin;
		if (trace) {
			cout << " unprocessable tiny box: now uplo<=" << setprecision(12) <<  uplo_of_epsboxes << " uplo " << uplo << endl;
		}
	}
}

void OptimizerANN::handle_cell(Cell& c, const IntervalVector& init_box ){

	contract_and_bound(c, init_box);

	if(!c.box.is_empty() && father > 0) {
		cout << "\"" << father << "\":f0 -> \"" << c.get<CellData>().id << "\":f0 [" << endl;
		cout << "id = " << id << endl;
		cout << "];" << endl;
		id++;
	}

	if (c.box.is_empty()) {
		delete &c;
	} else {
		buffer.push(&c);
	}
}

void OptimizerANN::contract_and_bound(Cell& c, const IntervalVector& init_box) {
	/*======================== contract y with y<=loup ========================*/
	Interval& y=c.box[goal_var];

	double ymax;
	if (loup==POS_INFINITY) ymax = POS_INFINITY;
	// ymax is slightly increased to favour subboxes of the loup
	// TODO: useful with double heap??
	else ymax = compute_ymax()+1.e-15;

	y &= Interval(NEG_INFINITY,ymax);

	if (y.is_empty()) {
		c.box.set_empty();
		return;
	}

	// ctc.contract(c.box);
	IntervalVector boxOld = c.box;
	cout << "\"" << c.get<CellData>().id << "\" [" << endl;
	cout << "label = \" <f0> " << c.get<CellData>().id;

	/*
	c.get<CellData>().HC4.clear();
	c.get<CellData>().ACID.clear();
	c.get<CellData>().COMPO.clear();
	*/

	// contractor HC4 y ACID
	// HC4 -> 0
	if(!c.box.is_empty()) {
		boxOld = c.box;
		try {
			ctc.list[0].contract(c.box);
		}
		catch(Exception& e) { // ibex exceptions
			cout << "Error " << 0 << endl;
			throw e;
		}
		catch (std::exception& e) { // other exceptions
			cout << "Error " << 0 << endl;
			throw e;
		}
		catch (...) {
			ibex_error("contract: cannot handle exception");
		}
		for(int i=0; i < c.box.size(); i++) {
			if(c.box.is_empty()) c.get<CellData>().HC4.push_back(0);
			else if(boxOld[i].diam() != c.box[i].diam()) c.get<CellData>().HC4.push_back(1);
			else  c.get<CellData>().HC4.push_back(0);
		}
		if(c.box.is_empty())  c.get<CellData>().HC4.push_back(1);
		else  c.get<CellData>().HC4.push_back(0);
		/*
		for (it=c.get<CellData>().HC4.begin(); it!=c.get<CellData>().HC4.end(); ++it)
			cout << *it << ".0 ";
		*/

		cout << "|";
		for(int i=0; i < c.box.size(); i++) {
			if(c.box.is_empty()) cout << 0 << " ";
			else if(boxOld[i].diam() != c.box[i].diam()) cout << 1 << " ";
			else  cout << 0 << " ";
		}
		if(c.box.is_empty())  cout << 1 << " ";
		else  cout << 0 << " ";
	}

	cout << "|";
	// ACID -> 1
	if(!c.box.is_empty()) {
		boxOld = c.box;
		if (c.box.is_empty()) return;
		try {
			ctc.list[1].contract(c.box);
		}
		catch(Exception& e) { // ibex exceptions
			cout << "Error " << 1 << endl;
			throw e;
		}
		catch (std::exception& e) { // other exceptions
			cout << "Error " << 1 << endl;
			throw e;
		}
		catch (...) {
			ibex_error("contract: cannot handle exception");
		}
		for(int i=0; i < c.box.size(); i++) {
			if(c.box.is_empty()) c.get<CellData>().ACID.push_back(0);
			else if(boxOld[i].diam() != c.box[i].diam()) c.get<CellData>().ACID.push_back(1);
			else  c.get<CellData>().ACID.push_back(0);
		}
		if(c.box.is_empty())  c.get<CellData>().ACID.push_back(1);
		else  c.get<CellData>().ACID.push_back(0);
		/*
		for (it=c.get<CellData>().ACID.begin(); it!=c.get<CellData>().ACID.end(); ++it)
			cout << *it << ".0 ";
		*/
		for(int i=0; i < c.box.size(); i++) {
			if(c.box.is_empty()) cout << 0 << " ";
			else if(boxOld[i].diam() != c.box[i].diam()) cout << 1 << " ";
			else  cout << 0 << " ";
		}
		if(c.box.is_empty())  cout << 1 << " ";
		else  cout << 0 << " ";
	}

	cout << "|";
	// COMPO -> 2
	if(!c.box.is_empty()) {
		// training ANN with COMPO
		boxOld = c.box;
		if(c.get<CellData>().id < trainingdata) {
			if (c.box.is_empty()) return;
			try {
				ctc.list[2].contract(c.box);
			}
			catch(Exception& e) { // ibex exceptions
				cout << "Error " << 2 << endl;
				throw e;
			}
			catch (std::exception& e) { // other exceptions
				cout << "Error " << 2 << endl;
				throw e;
			}
			catch (...) {
				ibex_error("contract: cannot handle exception");
			}
			for(int i=0; i < c.box.size(); i++) {
				if(c.box.is_empty()) c.get<CellData>().COMPO.push_back(0);
				else if(boxOld[i].diam() != c.box[i].diam()) c.get<CellData>().COMPO.push_back(1);
				else  c.get<CellData>().COMPO.push_back(0);
			}
			if(c.box.is_empty())  c.get<CellData>().COMPO.push_back(1);
			else  c.get<CellData>().COMPO.push_back(0);


			vector<double> inputVals, targetVals;
			//cout << endl << "in: ";
			vector<int>::iterator it;
			if(c.get<CellData>().ACID.size() > 0) {
				for (it=c.get<CellData>().ACID.begin(); it!=c.get<CellData>().ACID.end(); ++it) {
					//cout << *it << ".0 ";
					inputVals.push_back(*it);
				}
			}
			//cout << endl;

			//cout << "out: ";
			int aux = 0;
			for (it=c.get<CellData>().COMPO.begin(); it!=c.get<CellData>().COMPO.end(); ++it) {
				aux += (int) *it;
				//cout << *it << ".0 ";
				targetVals.push_back(*it);
			}
			//cout << endl;

			ann.trainingNeuron(inputVals, targetVals);

			//cout << "input size " << inputVals.size() <<  endl;
			//cout << "output size " << targetVals.size() <<  endl;

		// testing ANN with COMPO
		} else {
			vector<double> inputVals, targetVals, resultsVals;

			//cout << endl << "in: ";
			vector<int>::iterator it;
			if(c.get<CellData>().ACID.size() > 0) {
				for (it=c.get<CellData>().ACID.begin(); it!=c.get<CellData>().ACID.end(); ++it) {
					//cout << *it << ".0 ";
					inputVals.push_back(*it);
				}
			}
			//cout << endl;

			resultsVals = ann.testingNeuron(inputVals, inputVals);
			//cout << "out: ";
			int aux = 0;
			for (int i=0; i< resultsVals.size(); i++) {
				//cout << resultsVals[i] << " ";
				targetVals.push_back(*it);
			}
			//cout << endl;

			BitSet contractors(resultsVals.size()-1);


			// contrae los valores sobre threshold
			int contract = 0;
			for(int i=0; i<resultsVals.size()-1;i++) {
				if(resultsVals[i] > threshold) {
					contractors.add(i);
					contract++;
				}
			}
			// Si no contrae nada y el ultimo valor de resultsVals
			// esta sobre threshold se contrae uno de forma aleatoria
			if(contract==0 && resultsVals[resultsVals.size()-1] > threshold) contractors.add(resultsVals.size()-1);


			/*
			if(resultsVals[resultsVals.size()-1] > 0.3) {
				for(int i=0; i<resultsVals.size()-1;i++) {
					contractors.add(i);
				}
			} else {
				for(int i=0; i<resultsVals.size()-1;i++) {
					if(resultsVals[i] > 0.2) contractors.add(i);
				}
			}
			*/

			/*
			if(resultsVals[resultsVals.size()-1] > 0.7) {
				contractors.add(0);
			} else {
				for(int i=0; i<resultsVals.size()-1;i++) {
					if(resultsVals[i] > 0.5) contractors.add(i);
				}
			}
			*/



			// contract all
			/*
			for(int i=0; i<resultsVals.size()-1;i++) {
				contractors.add(i);
			}
			*/


			// COMPO -> 2
			boxOld = c.box;
			try {
				CtcFixPoint* fixpoint = dynamic_cast<CtcFixPoint*> (&ctc.list[2]);
				CtcCompo* compo = dynamic_cast<CtcCompo*> (&fixpoint->ctc);
				CtcPolytopeHull* poly = dynamic_cast<CtcPolytopeHull*> (&compo->list[0]);
				poly->set_contracted_vars(contractors);
				poly->contract(c.box);

			}
			catch(Exception& e) { // ibex exceptions
				cout << "Error " << 2 << endl;
				throw e;
			}
			catch (std::exception& e) { // other exceptions
				cout << "Error " << 2 << endl;
				throw e;
			}
			catch (...) {
				ibex_error("contract: cannot handle exception");
			}
			for(int i=0; i < c.box.size(); i++) {
				if(c.box.is_empty()) c.get<CellData>().COMPO.push_back(0);
				else if(boxOld[i].diam() != c.box[i].diam()) c.get<CellData>().COMPO.push_back(1);
				else  c.get<CellData>().COMPO.push_back(0);
			}
			if(c.box.is_empty())  c.get<CellData>().COMPO.push_back(1);
			else  c.get<CellData>().COMPO.push_back(0);


		}

		for(int i=0; i < c.box.size(); i++) {
			if(c.box.is_empty()) cout << 0 << " ";
			else if(boxOld[i].diam() != c.box[i].diam()) cout << 1 << " ";
			else  cout << 0 << " ";
		}
		if(c.box.is_empty())  cout << 1 << " ";
		else  cout << 0 << " ";
	}




	cout << "\"" << endl;
	cout << "shape = \"record\"" << endl;
	cout << "];" << endl;


	if (c.box.is_empty()) return;

	//cout << " [contract]  x after=" << c.box << endl;
	//cout << " [contract]  y after=" << y << endl;
	/*====================================================================*/

	/*========================= update loup =============================*/

	IntervalVector tmp_box(n);
	read_ext_box(c.box,tmp_box);

//	entailed = &c.get<EntailedCtr>();
//	if (!update_entailed_ctr(tmp_box)) {
//		c.box.set_empty();
	//		return;
//	}

	bool loup_ch=update_loup(tmp_box);

	// update of the upper bound of y in case of a new loup found
	if (loup_ch) y &= Interval(NEG_INFINITY,compute_ymax());

	//TODO: should we propagate constraints again?

	loup_changed |= loup_ch;

	if (y.is_empty()) { // fix issue #44
		c.box.set_empty();
		return;
	}

	/*====================================================================*/
	// Note: there are three different cases of "epsilon" box,
	// - NoBisectableVariableException raised by the bisector (---> see optimize(...)) which
	//   is independent from the OptimizerANN
	// - the width of the box is less than the precision given to the OptimizerANN ("prec" for the original variables
	//   and "goal_abs_prec" for the goal variable)
	// - the extended box has no bisectable domains (if prec=0 or <1 ulp)
	if ((tmp_box.max_diam()<=eps_x && y.diam() <=abs_eps_f) || !c.box.is_bisectable()) {
		update_uplo_of_epsboxes(y.lb());
		c.box.set_empty();
		return;
	}

	// ** important: ** must be done after upper-bounding
	//kkt.contract(tmp_box);

	if (tmp_box.is_empty()) {
		c.box.set_empty();
	} else {
		// the current extended box in the cell is updated
		write_ext_box(tmp_box,c.box);
	}
}

OptimizerANN::Status OptimizerANN::optimize(const IntervalVector& init_box, double obj_init_bound) {

	loup=obj_init_bound;

	// Just to initialize the "loup" for the buffer
	// TODO: replace with a set_loup function
	buffer.contract(loup);

	uplo=NEG_INFINITY;
	uplo_of_epsboxes=POS_INFINITY;

	nb_cells=0;

	buffer.flush();

	Cell* root=new Cell(IntervalVector(n+1));

	//cout << root->box.size() << " " << n << endl;
	//getchar();

	root->add<CellData>();
	root->get<CellData>().id = id;

	write_ext_box(init_box,root->box);

	// add data required by the bisector
	bsc.add_backtrackable(*root);
	// bsc.add_backtrackable(*root);

	// add data required by the buffer
	buffer.add_backtrackable(*root);
	// buffer.add_backtrackable(*root);

	loup_changed=false;
	initial_loup=obj_init_bound;

	// TODO: no loup-point if handle_cell contracts everything
	loup_point=init_box;
	time=0;
	Timer timer;
	timer.start();

	handle_cell(*root,init_box);
	id++;

	update_uplo();


	try {
	     while (!buffer.empty()) {

			loup_changed=false;
			// for double heap , choose randomly the buffer : top  has to be called before pop
			// Cell *c = buffer.top();
			Cell *c = buffer.top();

			if (trace >= 2) cout << " current box " << c->box << endl;

			try {

				father = c->get<CellData>().id;

				pair<IntervalVector,IntervalVector> boxes=bsc.bisect(*c);

				pair<Cell*,Cell*> new_cells=c->bisect(boxes.first,boxes.second);


				buffer.pop();
				delete c; // deletes the cell.

				nb_cells+=2;  // counting the cells handled ( in previous versions nb_cells was the number of cells put into the buffer after being handled)

				new_cells.first->get<CellData>().id = id;
				handle_cell(*new_cells.first, init_box);

				new_cells.second->get<CellData>().id = id;
				handle_cell(*new_cells.second, init_box);

				if (uplo_of_epsboxes == NEG_INFINITY) {
					cout << " possible infinite minimum " << endl;
					break;
				}
				if (loup_changed) {
					// In case of a new upper bound (loup_changed == true), all the boxes
					// with a lower bound greater than (loup - goal_prec) are removed and deleted.
					// Note: if contraction was before bisection, we could have the problem
					// that the current cell is removed by contractHeap. See comments in
					// older version of the code (before revision 284).

					double ymax=compute_ymax();

					buffer.contract(ymax);
				
					//cout << " now buffer is contracted and min=" << buffer.minimum() << endl;

					// TODO: check if happens. What is the return code in this case?
					if (ymax <= NEG_INFINITY) {
						if (trace) cout << " infinite value for the minimum " << endl;
						break;
					}
				}
				update_uplo();
				if (timeout>0) timer.check(timeout); // TODO: not reentrant, JN: done
				time = timer.get_time();

			}
			catch (NoBisectableVariableException& ) {
				update_uplo_of_epsboxes((c->box)[goal_var].lb());
				buffer.pop();
				delete c; // deletes the cell.
				update_uplo(); // the heap has changed -> recalculate the uplo (eg: if not in best-first search)

			}
		}
	}
	catch (TimeOutException& ) {
		status = TIME_OUT;
		return status;
	}

	timer.stop();
	time = timer.get_time();

	if (uplo_of_epsboxes == POS_INFINITY && (loup==POS_INFINITY || (loup==initial_loup && abs_eps_f==0 && rel_eps_f==0)))
		status=INFEASIBLE;
	else if (loup==initial_loup)
		status=NO_FEASIBLE_FOUND;
	else if (uplo_of_epsboxes == NEG_INFINITY)
		status=UNBOUNDED_OBJ;
	else if (get_obj_rel_prec()>rel_eps_f && get_obj_abs_prec()>abs_eps_f)
		status=UNREACHED_PREC;
	else
		status=SUCCESS;

	return status;
}

void OptimizerANN::report(bool verbose) {

	if (!verbose) {
		cout << get_status() << endl;
		cout << get_uplo() << ' ' << get_loup() << endl;
		for (int i=0; i<n; i++) {
			if (i>0) cout << ' ';
			cout << get_loup_point()[i].lb();
			if (loup_finder.rigorous())
				cout << ' ' << get_loup_point()[i].ub();
		}
		cout << endl << get_time() << " " << get_nb_cells() << endl;
		return;
	}

	switch(status) {
	case SUCCESS: cout << "\033[32m" << " optimization successful!" << endl;
	break;
	case INFEASIBLE: cout << "\033[31m" << " infeasible problem" << endl;
	break;
	case NO_FEASIBLE_FOUND: cout << "\033[31m" << " no feasible point found (the problem may be infeasible)" << endl;
	break;
	case UNBOUNDED_OBJ: cout << "\033[31m" << " possibly unbounded objective (f*=-oo)" << endl;
	break;
	case TIME_OUT: cout << "\033[31m" << " time limit " << timeout << "s. reached " << endl;
	break;
	case UNREACHED_PREC: cout << "\033[31m" << " unreached precision" << endl;
	}

	cout << "\033[0m" << endl;

	// No solution found and optimization stopped with empty buffer  before the required precision is reached => means infeasible problem
	if (buffer.empty() && uplo_of_epsboxes == POS_INFINITY && (loup==POS_INFINITY || (loup==initial_loup && abs_eps_f==0 && rel_eps_f==0))) {
		cout << " infeasible problem " << endl;
	} else {
		cout << " best bound in: [" << uplo << "," << loup << "]" << endl;

		double rel_prec=get_obj_rel_prec();
		double abs_prec=get_obj_abs_prec();

		cout << " relative precision obtained on objective function: " << rel_prec << " " <<
				(rel_prec <= rel_eps_f? " [passed]" : " [failed]") << endl;

		cout << " absolute precision obtained on objective function: " << abs_prec << " " <<
				(abs_prec <= abs_eps_f? " [passed]" : " [failed]") << endl;

		if (loup==initial_loup)
			cout << " no feasible point found " << endl;
		else {
			cout << " best feasible point: ";

			if (loup_finder.rigorous())
				cout << loup_point << endl;
			else
				cout << loup_point.lb() << endl;
		}
	}
	cout << " cpu time used: " << time << "s." << endl;
	cout << " number of cells: " << nb_cells << endl;
}



} // end namespace ibex
