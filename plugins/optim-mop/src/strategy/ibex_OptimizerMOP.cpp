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
bool OptimizerMOP::imprimir = true;



void print_sol(Vector& v){
	cout << "feasible_solution:";
	for(int i=0; i<v.size(); i++)
		cout << v[i] << " ";
	cout << endl;
}






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

void OptimizerMOP::polytope_search(IntervalVector box2){
	vector<IntervalVector> polytope_vertex;
	IntervalVector point(n);
	try{
		finder.clear();
		//encuentra m puntos factibles, donde m es el
		//numero de funciones objetivo
		for(int i=0; i<goals.size(); i++){
			point = finder.find2(box2,box2,POS_INFINITY).first;
			polytope_vertex.push_back(point);
		}

		//Se busca el hull de los puntos.
		IntervalVector box_hull(n);
		for(int j = 0 ; j < n ; j++)
			box_hull[j] = Interval(polytope_vertex[0][j].lb(),polytope_vertex[0][j].ub());
		for (int i = 1 ; i < polytope_vertex.size() ; i++)
			for(int j = 0 ; j < n ; j++){
				if (box_hull[j].lb() > polytope_vertex[i][j].lb())
					box_hull[j] = Interval(polytope_vertex[i][j].lb(),box_hull[j].ub());
				if (box_hull[j].ub() < polytope_vertex[i][j].ub())
					box_hull[j] = Interval(box_hull[j].lb(),polytope_vertex[i][j].ub());
			}

		IntervalVector input_box(box_hull);
		input_box.resize(n+goals.size());

		for (int i = 0 ; i < goals.size() ; i++)
			input_box[n+i]=0.0;

		IntervalVector corner = box_hull.lb();
		IntervalVector g_corners(goals.size());

		for(int i = 0 ; i < goals.size() ; i++)
			g_corners[i] = eval_goal(goals[i],corner,n,goals.size());

		vector<IntervalVector> J_eval(goals.size());

		for(int i = 0 ; i < goals.size() ; i++)
			J_eval[i] = goals[i].gradient(input_box);

		vector<pair <Vector, int> > linear_obj_f;

		/*linealizacion de los objetivos*/
		for(int i = 0 ; i < goals.size() ; i++)
			linear_obj_f.push_back(linearize_goal(box_hull,corner,J_eval[i],g_corners[i], n));

		IntervalVector p_point(goals.size());
		vector <IntervalVector> p_points;

		std::vector<std::vector<float>> vector_points;
		vector<float> aux;

		/*evaluacion de los objetivos*/
		for (int k = 0 ; k < goals.size() ; k++){
			for (int i = 0 ; i < goals.size() ; i++)
				p_point[i] = eval_fl(linear_obj_f[i], polytope_vertex[k]);
			p_points.push_back(p_point);
		}

		for (int i = 0 ; i < p_points.size() ; i++){
			for (int j = 0 ; j < p_points[i].size() ; j++)
				aux.push_back(p_points[i][j].ub());
			vector_points.push_back(aux);
			aux.clear();
		}
		/*ordenando puntos para poder borrar los repetidos*/
		std::sort(vector_points.begin(), vector_points.end());
		vector_points.erase(std::unique(vector_points.begin(), vector_points.end()), vector_points.end());

		char filename[ ] = "polytope_points.txt";
		fstream thefile;
		/*los puntos se guardan en el archivo*/
		thefile.open(filename, fstream::app);
		thefile <<"(";
		for (int i = 0 ; i < vector_points.size() ; i++){
			thefile <<"(";
			for (int j = 0 ; j < vector_points[i].size() ; j++){
				thefile << vector_points[i][j];
				if (j < vector_points[i].size()-1) thefile << ", ";
				else thefile << ")";
			}
			if (i < vector_points.size()-1) thefile << ";";
			else thefile << ")\n";
		}
		thefile.close();

	}

	catch (LoupFinder::NotFound& ) {

	}
}


bool OptimizerMOP::upper_bounding(const IntervalVector& box, double timer) {
	//We attempt to find two feasible points which minimize both objectives
	//and the middle point between them
	IntervalVector box2(box);
	IntervalVector box2_y = OptimizerMOP::get_box_y(box2);
	box2.resize(n);
	IntervalVector xa(n), xb(n);
	finder.clear();

	Vector mid=box2.mid();


	if (finder.norm_sys.is_inner(mid) && nds_mode==POINTS){
		Vector v(goals.size());
		for(int i=0; i<goals.size(); i++) {
			//evaluation of function n in mid point of box of variables (x1, x2, x3, ...)
			v[i] = eval_goal(goals[i], mid, n, goals.size()).ub();
		}


		//addPoint ndsh2
		ndsh2.addPoint(v, NDS_X(mid));
		print_sol(v);
		

	}

	/*testing Victor*/
//	polytope_search(box2);

	if(nds_mode==SEGMENTS) {
		list<Vector> simplex_points_eval;
		list <pair <Vector, Vector> > points;
		list <pair <Vector, Vector> > simplex_vertex; //Save vertex simplex inside of polytope

		simplex_points_eval.clear();
		points. clear();





			try{
				simplex_vertex.clear();
				int phase = 0;

				//Add midpoint
				Vector v(goals.size());
				for(int i=0; i<goals.size(); i++)	v[i] = eval_goal(goals[i], mid, n, goals.size()).ub();
				if(finder.norm_sys.is_inner(mid)){
					points.push_back(make_pair(mid, v));
					ndsh2.addPoint(v, NDS_X(mid));
				}

				//Search all the points inside of simplex
				while(true){
					xa = finder.find2(box2,box2,POS_INFINITY).first;
					if(simplex_vertex.back().first != xa.mid()){
						//Evaluation point in simplex
						Vector v(goals.size());
						for(int i=0; i<goals.size(); i++) v[i] = eval_goal(goals[i], xa.mid(), n, goals.size()).ub();

						//Save external vertex of simplex for testing
						if(phase < goals.size())
							simplex_vertex.push_back(make_pair(xa.mid(), v));

						//Check dominance of points
						finder.add_nondominated_point(make_pair(xa.mid(),v), points);

					}
					++phase;
				}

			}catch (LoupFinder::NotFound& ) {

				
				if(points.size() >= 3){

					/*
					 * Add points to ndsh2
					 */

					//N POINTS
					if(finder.get_nds_simplex() == finder.get_nb_sol()){
						for(auto itr = points.begin(); itr != points.end(); ++itr){
							ndsh2.addPoint(itr->second, NDS_X(itr->first));
							print_sol(itr->second);
						}
					}

					//D POINTS
					if(finder.get_nds_simplex() != finder.get_nb_sol()){
						for(auto itr = simplex_vertex.begin(); itr != simplex_vertex.end(); ++itr ){
							ndsh2.addPoint(itr->second, NDS_X(itr->first));
							print_sol(itr->second);
						}
					}


				}else{
					//actualizar upper envelope cuando solo se agrega el punto medio

					Vector v(goals.size());
					for(int i=0; i<goals.size(); i++)	v[i] = eval_goal(goals[i], mid, n, goals.size()).ub();
					if(finder.norm_sys.is_inner(mid)){
						update_upper_envelope_mid(v);
						print_sol(v);
					}


				}
				return true;
			}

		}

	return true;

}

void OptimizerMOP::update_upper_envelope(IntervalVector box_y, Vector hyperplane_coefficients){

	if(ndsh2.is_dominated(box_y.lb())) return;

	//por cada punto se verifica que la caja no este dominada, si lo esta se elimina
	for(auto itr_point = ndsh2.NDS.begin(); itr_point != ndsh2.NDS.end(); ++itr_point){

		auto it_box = upper_envelope.begin();
		while ( it_box != upper_envelope.end()){
			if( ndsh2.is_dominated(it_box->first.lb(), itr_point->first) ){
				it_box = upper_envelope.erase(it_box);
			}else
				++it_box;
		}

	}


	upper_envelope.push_back(make_pair(box_y, hyperplane_coefficients));
}


void OptimizerMOP::update_upper_envelope2(IntervalVector box_y, Vector hyperplane_coefficients, list <pair <Vector, Vector> > points_box){

	//check lower bound dominance wrt ndsh
	if(ndsh2.is_dominated(box_y.lb())) return;

	//por cada punto se verifica que la caja no este dominada, si lo esta se elimina
	for(auto itr_point = points_box.begin(); itr_point != points_box.end(); ++itr_point){

		auto it_box = upper_envelope.begin();
		while ( it_box != upper_envelope.end()){
			if( ndsh2.is_dominated(it_box->first.lb(), itr_point->second) ){
				it_box = upper_envelope.erase(it_box);
			}else
				++it_box;
		}
	}


	upper_envelope.push_back(make_pair(box_y, hyperplane_coefficients));
}


void OptimizerMOP::update_upper_envelope_mid(Vector mid){
	auto it = upper_envelope.begin();
	while ( it != upper_envelope.end()){
		if( ndsh2.is_dominated(it->first.lb(), mid) ){
			it = upper_envelope.erase(it);
		}else
			++it;
	}
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

/*
 * Dominance Peeler for n objective functions
 */
void OptimizerMOP::dominance_peeler_n(IntervalVector& boxy, list < pair< Vector, int > >& cuttingPoints){
	for(auto itr = cuttingPoints.begin(); itr != cuttingPoints.end(); ++itr)
		boxy[itr->second] = Interval(boxy[itr->second].lb(), itr->first[itr->second]);
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
}

void OptimizerMOP::contract_and_bound(Cell& c, const IntervalVector& init_box) {
	//Obtain box from cell
	IntervalVector boxy=get_box_y(c.box);

	//Return a list of point non-dominated by box_y.lb()
	list< pair< Vector, int >> cuttingPoints = ndsh2.cutting_points(boxy.lb(), boxy.ub());

	// Dominance Peeler contract box_y
	dominance_peeler_n(boxy, cuttingPoints);
	set_box_y(c.box, boxy);

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

void OptimizerMOP::pre_optimize(const IntervalVector& init_box, Cell* root){
	//the box in cells have the n original variables plus the two objective variables (y1 and y2)

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

OptimizerMOP::Status OptimizerMOP::optimize(const IntervalVector& init_box) {

	status=SUCCESS;

	//the box in cells have the n original variables plus the objective variables (y1, y2, ..., yn)
	Cell* root=new Cell(IntervalVector(n+goals.size()));
	pre_optimize(init_box, root);
	Timer timer;
	timer.start();
	cells_map cells;
	int id=0;
	cells.insert(make_pair(root,id++));
	//map<Cell*, int> idbox;
	

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

			
			//if(rel_eps>0.0)	eps=std::max(focus[0].diam(),focus[1].diam())*rel_eps;
			buffer.pop();
			cells.erase(c);

			if(cdata->ub_distance <= eps){
				delete c;
				if(dynamic_cast<DistanceSortedCellBufferMOP*>(&buffer)) break;
				else continue;
			}
			nb_cells++;

			contract_and_bound(*c, init_box);

			if (c->box.is_empty()) {
				delete c;
				continue;
			}

			upper_bounding(c->box, timer.get_time());

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
				// Distance ndsh2
				dist2=ndsh2.distance(c->box);
			}

			if(atomic_box || dist2 < eps){
				if(new_cells.first){
					delete new_cells.first;
					delete new_cells.second;
				}

			delete c;

			continue;

			}

			delete c; // deletes the cell.

			buffer.push(new_cells.first);
			cells.insert(make_pair(new_cells.first,id++));

			buffer.push(new_cells.second);
			cells.insert(make_pair(new_cells.second,id++));

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
		cout << "left_boxes:" << endl;
		for(auto cc : cells){
			for(int i=0; i<goals.size(); i++)
				cout << cc.first->box[n+i].lb() << " " << cc.first->box[n+i].ub() << " ";
			cout << endl;
		}
		cout << "end" << endl;
	}



	timer.stop();
	time = timer.get_time();

//	py_Plotter::offline_plot(ndsH.NDS2, NULL, "output2.txt");
	py_Plotter::offline_plot(ndsh2.NDS, NULL, "eval_solutions.txt", "solutions.txt");

	//save upper envelope hyperplanes
	py_Plotter::upper_envelope_save(upper_envelope);

	return status;
}

void OptimizerMOP::report(bool verbose) {

	if (!verbose) {
		cout << endl 	<< "time 	#nodes 		|Y|		|H|" << endl;
		cout << get_time() << " " << get_nb_cells() << " " << ndsh2.size() <<" "<< upper_envelope.size()<<endl;

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

}


} // end namespace ibex
