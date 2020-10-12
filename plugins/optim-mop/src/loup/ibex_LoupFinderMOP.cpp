/*
 * ibex_LoupFinderMOP.cpp
 *
 *  Created on: 5 oct. 2017
 *      Author: iaraya
 */

#include "ibex_LoupFinderMOP.h"
#include "ibex_PdcHansenFeasibility.h"
#include "ibex_FncActivation.h"
#include "ibex.h"
#include <random>
#include <chrono>
#include <bits/stdc++.h>
#include "ibex_NDS2.h"
#include "ibex_Interval.h"
#include <sstream>



namespace ibex {

 double LoupFinderMOP::_weight2=0.0;


 LoupFinderMOP::LoupFinderMOP(const System& sys, const Array<const Function> goals, double eqeps, int nb_sol, int nds_simplex) :
 		sys(sys), norm_sys(sys,eqeps), lr(norm_sys,LinearizerXTaylor::RESTRICT, LinearizerXTaylor::INF),
 		lp_solver(sys.nb_var, std::max((sys.nb_var)*3,100/*LPSolver::default_max_iter*/)),
 		goals(goals), has_equality(false), nb_sol(nb_sol), nds_simplex(nds_simplex), phase(0), vec1(norm_sys.nb_var), vec2(norm_sys.nb_var) {

 	if (sys.nb_ctr>0)
 		// ==== check if the system contains equalities ====
 		for (int i=0; i<sys.f_ctrs.image_dim(); i++) {
 			if (sys.ops[i]==EQ) {
 				(bool&) has_equality = true;
 				break;
 			}
 		}
 }

bool LoupFinderMOP::ub_correction(Vector p, IntervalVector& res){
    //if(!norm_sys.is_inner(p)) return false;

	if (!has_equality && norm_sys.is_inner(p)){
		res=IntervalVector(p);
		return true;
	}

	//sys is the original system (with equations)

	FncActivation af(sys,p,NormalizedSystem::default_eps_h);

	if (af.image_dim()==0) {
		res=IntervalVector(p);
		return true;
	}
	//cout << "gere" << endl;

	IntervalVector epsbox(p);

	//epsbox.inflate(NormalizedSystem::default_eps_h);
	//PdcHansenFeasibility pdc(af, false);
	// ====================================================
	// solution #2: we call Hansen test in inflating mode.
	PdcHansenFeasibility pdc(af, true);
	// ====================================================
	//cout << epsbox << endl;
	//cout << af.eval(0,epsbox) << endl;
	if (pdc.test(epsbox)==YES) {
		//note: don't call is_inner because it would check again all equalities (which is useless
		// and perhaps wrong as the solution box may violate the relaxed inequality (still, very unlikely))
		bool satisfy_inequalities=true;
		for (int j=0; j<sys.f_ctrs.image_dim(); j++) {
			if (!af.activated()[j]){

				if (((sys.ops[j]==LEQ || sys.ops[j]==LT)
					  && sys.f_ctrs.eval(j,pdc.solution()).ub()>0)
						||
					((sys.ops[j]==GEQ || sys.ops[j]==GT)
					  && sys.f_ctrs.eval(j,pdc.solution()).lb()<0)) {

					/* TODO: && !entailed->original(j)*/
						satisfy_inequalities=false;
						cout << "not corrected" << endl;
						break;
					}
			}
		}
		if (satisfy_inequalities) {
			res = pdc.solution();
			//cout << "corrected" << endl;
			return true;
		}
	}
	//===========================================================
	return false;
}

std::pair<IntervalVector, double> LoupFinderMOP::find(const IntervalVector& box, const IntervalVector& lp, double l){
	int n=norm_sys.nb_var;
	//phase 0 or 1: call to simplex
    if(phase < nb_sol && phase<=1 && lp_solver.default_max_bound > box.max_diam() ){
		lp_solver.clean_ctrs();
		lp_solver.set_bounds(box);

		IntervalVector box2(box);
		//box2.resize(n+2);
		//box2[n]=0.0; box2[n+1]=0.0;
		box2.resize(n+goals.size());
		for(int i=n; i<n+goals.size(); i++) box2[n] = 0.0;


		IntervalVector ig= (phase==0 && (nb_sol>1 || rand()%2==0))?
						(goals[0].gradient(box2.mid())+ _weight2*goals[1].gradient(box2.mid())) :
						(goals[1].gradient(box2.mid())+ _weight2*goals[0].gradient(box2.mid()));




		if (ig.is_empty()){ // unfortunately, at the midpoint the function is not differentiable
			std::cout<<"at the midpoint the function is not differentiable"<<endl;
			phase = 0;
			throw NotFound(); // not a big deal: wait for another box...
		}
		ig.resize(n);
		Vector g=ig.mid();
		// set the objective coefficient
		// TODO: replace with lp_solver.set_obj(g) when implemented
		for (int j=0; j<n; j++)
			//lp_solver.set_obj(g);
			lp_solver.set_obj_var(j,g[j]);

		int count = lr.linearize(box,lp_solver);

		if (count==-1) {
			lp_solver.clean_ctrs();
			phase = 0;
			throw NotFound();
		}

		LPSolver::Status_Sol stat = lp_solver.solve();

		if (stat == LPSolver::OPTIMAL) {
			//the linear solution is mapped to intervals
			Vector loup_point(lp_solver.get_primal_sol());

			//correct the point
			for(int i=0;i<box.size();i++){
				if(loup_point[i] < box[i].lb())  loup_point[i] = box[i].lb(); //se encuentra dentro caja?
				if(loup_point[i] > box[i].ub())  loup_point[i] = box[i].ub();
			}

			if(phase==0) vec1=loup_point;
			if(phase==1) {
				if(vec1==vec2){
				    phase=0;
				    throw NotFound();
				 }
				vec2=loup_point;
			}
			phase ++;

			return make_pair(loup_point, 0.0);

		}else{
			phase = 0;
			throw NotFound();
		}
    }else if(phase>=2 && phase < nb_sol){
    	Vector vec3 = vec1 + (((double)phase-1.0)/((double)nb_sol))*(vec2-vec1);
		for(int i=0;i<box.size();i++){
		    if(vec3[i] < box[i].lb())  vec3[i] = box[i].lb();
			if(vec3[i] > box[i].ub())  vec3[i] = box[i].ub();
		}
		phase ++;
		return make_pair(vec3, 0.0);
    }

    phase=0;
    throw NotFound();

}

//Find modified to multiples objectives
std::pair<IntervalVector, double> LoupFinderMOP::find2(const IntervalVector& box, const IntervalVector& lp, double l){
	int n=norm_sys.nb_var;
	//std::cout<<"PHASE = "<<phase<<endl;
//Find coeficientes of every
	if(phase < nb_sol && phase<goals.size() && lp_solver.default_max_bound > box.max_diam() ){
		lp_solver.clean_ctrs();
		lp_solver.set_bounds(box);

		IntervalVector box2(box);
		box2.resize(n+goals.size());
		for(int i=n; i<n+goals.size(); i++) box2[n] = 0.0;

		//Sumatory of _weight2
		IntervalVector sum_coefficients(n+goals.size(),0);

		//Sum of the gradient of the objective functions evaluated at the midpoint
		int itr=phase;
		bool check = true;
		while(check){
			if(itr+1 > goals.size()-1 ){
				itr = 0;
			}else{
				itr++;
			}
			if(itr == phase) check=false;
			if(check) {
				sum_coefficients = goals[itr].gradient(box2.mid()) + sum_coefficients;
			}
		}

		//Coeficcient calculator
		IntervalVector ig = goals[phase].gradient(box2.mid())+ _weight2*sum_coefficients;

		if (ig.is_empty()){ // unfortunately, at the midpoint the function is not differentiable
			phase = 0;
			throw NotFound(); // not a big deal: wait for another box...
		}
		ig.resize(n);
		Vector g=ig.mid();

		// Set the objective coefficient
		for (int j=0; j<n; j++) lp_solver.set_obj_var(j,g[j]);

		int count = lr.linearize(box,lp_solver);

		if (count==-1) {
			//Linearizer problem
			lp_solver.clean_ctrs();
			phase = 0;
			throw NotFound();
		}

		LPSolver::Status_Sol stat = lp_solver.solve();

		if (stat == LPSolver::OPTIMAL) {
			//the linear solution is mapped to intervals
			Vector loup_point(lp_solver.get_primal_sol());

			//correct the point
			for(int i=0;i<box.size();i++){
				if(loup_point[i] < box[i].lb())  loup_point[i] = box[i].lb();
				if(loup_point[i] > box[i].ub())  loup_point[i] = box[i].ub();
			}

			if(phase == 0){
				solutions.clear(); //Clean all previous solutions
				solutions.push_back(make_pair(loup_point, phase));
			}
			else{
				if(!isSaved(loup_point)){
					solutions.push_back(make_pair(loup_point, phase));
				}
				if(phase == goals.size()-1 && solutions.size() == 1){
					//Existe una solucion factible en todas las funciones objetivos
					phase = 0;
					throw NotFound();
				}
			}
			phase ++;



			return make_pair(loup_point, 0.0);

			}else{
				phase = 0;
				throw NotFound();
			}

	 }else if(phase>=goals.size() && phase < nb_sol){
		 //Punto nuevo
		 Vector point(n, 0.0);

		 /*
		  * Generate random number between [0,1]
		  */
		 std::mt19937_64 rng;
		 // initialize the random number generator with time-dependent seed
		 uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
		 std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed>>32)};
		 rng.seed(ss);
		 // initialize a uniform distribution between 0 and 1
		 std::uniform_real_distribution<double> unif(0, 1);


		 /*
		  * Generate random points in n-Simplex
		  */

		 //Generacion de vertices base en n-simplex
		 double *randomNumber = new double [solutions.size()];
		 randomNumber[0] = unif(rng); //primer vertice
		 randomNumber[1] = 1 - randomNumber[0]; //segundo vertice

		 //Calculo de punto con respecto a los otros vertices
		 auto itr = solutions.begin();
		 point = randomNumber[0] * itr->first;
		 itr++;
		 point = point + randomNumber[1] * itr->first;
		 itr++;
		 double random;
		 for(int i=2; itr != solutions.end();i++){
			 if(i==2){
				 random = sqrtf(unif(rng));
			 }
			 else{
				 random = std::pow( unif(rng), (double)(1.0/i) );
			 }
			 point = random*point + (1-random)* itr->first;
			 itr++;
		 }

		 delete [] randomNumber;
		 phase ++;
		 return make_pair(point, 0.0);
	     }

	     phase=0;
	     throw NotFound();

}



bool LoupFinderMOP::isSaved(Vector loup_point){
	double min;
	int isequal=0;

	for(auto itr = solutions.begin(); itr != solutions.end(); ++itr) {
		isequal=0;
		for(int i=0; i<loup_point.size();i++){
			if( fabs(loup_point[i]) < fabs(itr->first[i]) ) min=fabs(loup_point[i]);
			else min = fabs(itr->first[i]);

			if(min < 1){
				if(fabs(itr->first[i] - loup_point[i]) <= 1e-8) isequal++;//test same value eps_x
			}else{
				if(fabs(itr->first[i] - loup_point[i])/ min <= 1e-8) isequal++;
			}
		}
		if(isequal == loup_point.size()){
			return true;
		}
	}
	return false;
}

void LoupFinderMOP::add_external_vertex(Vector base_point, list<Vector> list_points, list<Vector>& external_vertex, IntervalVector box_y){

	//check dominance base_point
	for(auto itr = list_points.begin(); itr != list_points.end(); ++itr){
		Vector comparison_point(*itr);
		NDS_XY a;
		if(comparison_point != base_point && a.is_dominated(base_point, comparison_point)) return;
	}

	//Find vertex points
	for(int dim=0; dim<base_point.size(); dim++){
		Vector vertex_point(base_point);
		vertex_point[dim] = box_y.ub()[dim]; //se inicia el punto de vertice en la dimesion del UB de la caja

		//Busqueda de cada vertice externo con respecto a los otros puntos
		for(auto itr = list_points.begin(); itr != list_points.end(); ++itr){
			Vector comparison_point(*itr);
			if( comparison_point != base_point){

				if(vertex_point[dim] > comparison_point[dim] && base_point[dim] < comparison_point[dim]){
					vertex_point[dim] = comparison_point[dim];
				}
			}
		}
	external_vertex.push_back(vertex_point);
	}
}

void LoupFinderMOP::add_nondominated_point(pair<Vector,Vector> point, list<pair<Vector,Vector>>& list_points){
	//Check if point is dominated
	bool point_is_dominated = false;
	for(auto itr=list_points.begin(); itr != list_points.end(); ++itr){
		if(is_dominated(point.second, itr->second)){
			point_is_dominated = true;
			itr = --list_points.end();
		}
	}
	if(!point_is_dominated){
		//Delete all points dominated by this new point
		auto it = list_points.begin();
		while(it != list_points.end()){
			if(is_dominated(it->second, point.second) ){
				it = list_points.erase(it);
			}else{
				++it;
			}
		}
		list_points.push_back(point);
	}

}


//Return True si y_prime es dominado por y
	bool LoupFinderMOP::is_dominated(const Vector& y_prime, const Vector& y){
		for(int i=0; i < y_prime.size(); i++){
			if(y[i] <= y_prime[i]){
			}else return false;
		}
		return true;
	}

Vector LoupFinderMOP::find_minimum(IntervalVector box, Vector hyperplane_coefficients, int dim_goal){
	Vector point(0);
	try{
		double prec1=1e-8; // precision
		// Build the system to find minimum point on the hyperplane inside of box:
		SystemFactory opt_factory;

		Variable xn[goals.size()]; // not to be confused with xn()
		Variable xn_ctr[goals.size()];
		Array<const ExprSymbol> vars(goals.size());
		Array<const ExprSymbol> vars_ctr(goals.size());
		for (int i=0; i<goals.size(); i++){
			vars.set_ref(i,xn[i]);
			vars_ctr.set_ref(i,xn_ctr[i]);
		}
		opt_factory.add_var(vars);


		/**
		 * Assign goal
		 */
		string goal_expr = to_string(hyperplane_coefficients[hyperplane_coefficients.size()-1]);
		int itr=dim_goal;
		bool check = true;
		const ExprNode* func_node=&( vars[itr]*0 - hyperplane_coefficients[hyperplane_coefficients.size()-1]);

		while(check){
			if(itr+1 > hyperplane_coefficients.size()-2 ){
				itr = 0;
			}else{
				itr++;
			}
			if(itr == dim_goal) {
				func_node = & (*func_node / hyperplane_coefficients[itr]);
				check=false;
			}
			if(check) {
				func_node = & (*func_node - vars[itr]*hyperplane_coefficients[itr]);
			}
		}

		Function f(vars, *func_node, "fu");
		opt_factory.add_goal(f);
		//std::cout<<"goal_expr = "<<str<<endl;


        /*
         * Contraint assign
         */
		const ExprNode* exp_node_ctr= &(vars[0]*hyperplane_coefficients[0]);

		for(int i = 1; i< goals.size(); i++){
			exp_node_ctr = &(*exp_node_ctr + vars[i]*hyperplane_coefficients[i]);
		}
		exp_node_ctr = &(*exp_node_ctr + hyperplane_coefficients[hyperplane_coefficients.size()-1]);

		//std::cout<<"ctr_expr = "<<*exp_node_ctr<<endl;

		ExprCtr num_constraint(*exp_node_ctr, CmpOp::EQ);
		opt_factory.add_ctr(num_constraint);

		System sys(opt_factory);

		for(int dim = 0; dim< box.size(); dim++){
			Interval aux(box[dim].lb() , box[dim].ub());
			sys.box[dim]=aux;
		}

		DefaultOptimizer optim(sys,prec1);
		optim.optimize(sys.box);
		//std::cout<<"status = "<<optim.get_status()<<endl;

		if(optim.get_status() == Optimizer::Status::SUCCESS){
			return optim.get_loup_point().mid();
		}

		return point;

	}catch(ibex::SyntaxError& e) {
	  cout << e << endl;
	  std::cout<<"error----"<<endl;
	  return point;
	}


}

void LoupFinderMOP::add_external_vertex2(pair<Vector,Vector> point, list <pair<Vector,Vector>> list_points, list<Vector>& external_vertex, IntervalVector box_y){

	//check dominance base_point
	for(auto itr = list_points.begin(); itr != list_points.end(); ++itr){
		Vector comparison_point(itr->second);
		NDS_XY a;
		if(comparison_point != point.second && a.is_dominated(point.second, comparison_point)) return;
	}
	//Find vertex points
	for(int dim=0; dim<point.second.size(); dim++){
		Vector vertex_point(point.second);
		vertex_point[dim] = box_y.ub()[dim]; //se inicia el punto de vertice en la dimesion del UB de la caja
		//busqueda de cada vertice externo con respecto a los otros puntos
		for(auto itr = list_points.begin(); itr != list_points.end(); ++itr){
			Vector comparison_point(itr->second);
			if( comparison_point != point.second){
				if(vertex_point[dim] > comparison_point[dim] && point.second[dim] < comparison_point[dim]){
					vertex_point[dim] = comparison_point[dim];
				}
			}
		}
	external_vertex.push_back(vertex_point);
	}

}

void LoupFinderMOP::print_hyperplane_intersection(IntervalVector caja_comp, Vector hyp_comp, IntervalVector caja_add, Vector hyp_add, IntervalVector caja_int, list<Vector> puntos_min_add, list<Vector> puntos_min_comp){
	//Imprimir
	std::cout<<"\n----imprimir"<<endl;

	std::cout<<"caja_comp"<<endl;
	for(int i = 0 ;i<caja_comp.size();i++)
		std::cout<<caja_comp[i].lb()<<" "<<caja_comp[i].ub()<<endl;

	std::cout<<"hyperplano_comp"<<endl;
	for(int i=0; i < hyp_comp.size(); i++){
		std::cout<<hyp_comp[i];
		if(hyp_comp.size() + 1 != hyp_comp.size()) std::cout<<" ";
	}
	std::cout<<"\n";


	std::cout<<"caja_add"<<endl;
	for(int i = 0 ;i<caja_add.size();i++){
		std::cout<<std::setprecision(15)<<fixed<<caja_add[i].lb();
	    std::cout<<" ";
	    std::cout<<setprecision(15)<<caja_add[i].ub()<<endl;
	}

	std::cout<<"hyperplano_add"<<endl;
	for(int i=0; i < hyp_add.size(); i++){
		std::cout<<hyp_add[i];
		if(hyp_add.size() + 1 != hyp_add.size()) std::cout<<" ";
	}
	std::cout<<"\n";


	std::cout<<"caja_interseccion"<<endl;
	for(int i = 0 ;i<caja_int.size();i++)
		std::cout<<caja_int[i].lb()<<" "<<caja_int[i].ub()<<endl;

	std::cout<<"puntos_minimos_add"<<endl;
	for(auto itr = puntos_min_add.begin(); itr != puntos_min_add.end(); ++itr ){
		for(auto itr2 = itr->begin(); itr2 != itr->end(); ++itr2) {
			std::cout<<*itr2;
			if( itr->size()+1 != itr->size()) std::cout<<" ";
		}
		std::cout<<"\n";
	}

	std::cout<<"puntos_minimos_comp"<<endl;
	for(auto itr = puntos_min_comp.begin(); itr != puntos_min_comp.end(); ++itr ){
		for(auto itr2 = itr->begin(); itr2 != itr->end(); ++itr2) {
			std::cout<<*itr2;
			if( itr->size()+1 != itr->size()) std::cout<<" ";
		}
		std::cout<<"\n";
	}

}

//Version 1.0
bool LoupFinderMOP::is_dominated_hyperplane_add(list< pair<IntervalVector, Vector> > upper_envelope, IntervalVector box_y, Vector hyperplane_coefficients){
	//Check if hyperplane is dominated, compare every hyperplane with this new hyperplane
	for(auto itr = upper_envelope.begin(); itr != upper_envelope.end(); ++itr){
		std::cout<<"-----"<<endl;
		std::cout<<"\n box_comp = "<<itr->first<<endl;
		std::cout<<"\n box = "<<box_y<<endl;

		//Check box intersect to verify hyperplane dominated
		bool intersect = true;
		for(int dim = 0 ; dim < box_y.size(); dim++){

			if( itr->first.lb()[dim] <= box_y.lb()[dim] &&
				box_y.lb()[dim] <= itr->first.ub()[dim]){
			}else{
				intersect = false;
				dim = box_y.size();
				std::cout<<"-----no existe una interseccion valida en la caja check hyperplano"<<endl;
			}
		}

		IntervalVector box_its(goals.size());
		list<Vector> min_Points_Hyp_add;
		list<Vector> min_Points_Hyp_comp;
		if(intersect){
			//Obtain box intersect
			for(int dim = 0 ;dim < box_y.size(); dim++){
				double box_lb, box_ub;
				//std::cout<<"box_y[dim].lb() = "<<box_y[dim].lb()<<endl;
				//box_its[dim].lb() = box_y[dim].lb();
				box_lb = box_y[dim].lb();

				if(box_y[dim].ub() <= itr->first.ub()[dim]){
					box_ub = box_y[dim].ub();
				}else{
					box_ub = itr->first.ub()[dim];
				}
				Interval aux(box_lb, box_ub);
				box_its[dim] = aux;
				//std::cout<<"box_its = "<<box_its<<endl;
			}

			std::cout<<"\n intersection box = "<<box_its<<endl;


			//Puntos que minimizan a las funciones en hyperplano dentro de la caja intersectada
			min_Points_Hyp_add.clear();
			min_Points_Hyp_comp.clear();
			for(int dim_goal=0; dim_goal<box_its.size(); dim_goal++){
				//Search minimum points of every function in the hyperplane to add inside the box intersection
				Vector point_min = find_minimum(box_its, hyperplane_coefficients, dim_goal);
				if(point_min.size() != 0){
					min_Points_Hyp_add.push_back(point_min);
				}
				//Search minimum points of every function in the hyperplane to comp inside the box intersection
				Vector point_min2 = find_minimum(box_its, itr->second, dim_goal);
				if(point_min2.size() != 0){
					min_Points_Hyp_comp.push_back(point_min2);
					min_Points_Hyp_add.push_back(point_min2);
				}

			}

			if(min_Points_Hyp_add.size() != 0 || min_Points_Hyp_comp.size() != 0){
			//Verificar que el hiperplano se  encuentra dominado
				double nb_hypadd_dominated = 0;
				double nb_hypadd_not_dominated = 0;
				for(auto itr_minpoints=min_Points_Hyp_add.begin(); itr_minpoints!=min_Points_Hyp_add.end(); ++itr_minpoints){
					Vector min_point(*itr_minpoints);
					//verificar si c es mayor o menor en cada minimo de ambos hyperplanos
					double pos_hyp_comp = 0;
					double pos_hyp_add = 0;
					for(int dim_hyp = 0; dim_hyp<itr_minpoints->size(); dim_hyp++){
						pos_hyp_comp += itr->second[dim_hyp]*min_point[dim_hyp];
						pos_hyp_add += hyperplane_coefficients[dim_hyp]*min_point[dim_hyp];
					}

					pos_hyp_comp += itr->second[goals.size()];
					pos_hyp_add  += hyperplane_coefficients[goals.size()];

					if(pos_hyp_comp <= pos_hyp_add)
						nb_hypadd_dominated++;
					else
						nb_hypadd_not_dominated++;


				}
				print_hyperplane_intersection(itr->first, itr->second, box_y, hyperplane_coefficients, box_its, min_Points_Hyp_add, min_Points_Hyp_comp);
				//Verificar si se encuentra dominado completamente en aquellos puntos
				if(nb_hypadd_dominated != 0 && nb_hypadd_not_dominated == 0 ){
					std::cout<<"add domina a comp parece, verificar"<<endl;

				}

				if(nb_hypadd_dominated == 0 && nb_hypadd_not_dominated != 0){
					std::cout<<"add domina3"<<endl;
					print_hyperplane_intersection(itr->first, itr->second, box_y, hyperplane_coefficients, box_its, min_Points_Hyp_add, min_Points_Hyp_comp);
					return true;
				}

				if(nb_hypadd_dominated != 0 && nb_hypadd_not_dominated != 0 ){
					std::cout<<"parece interseccion"<<endl;
				}

			}

		}

		//Imprimir
		//finder.print_hyperplane_intersection(itr->first, itr->second, box_y, hyperplane_coefficients, box_its, min_Points_Hyp_add, min_Points_Hyp_comp);
	}

	return false;

}

//hyperplano_comp domina a hiperplano_prime V 2.0
bool LoupFinderMOP::is_dominated_hyperplane(IntervalVector box_comp, Vector hyperplane_coefficients_comp, IntervalVector box_prime, Vector hyperplane_coefficients_prime){

	//Check box intersect to verify hyperplane dominated
	bool intersect = true;
	for(int dim = 0 ; dim < box_comp.size(); dim++){

		if( box_comp.lb()[dim] <= box_prime.lb()[dim] &&
				box_prime.lb()[dim] <= box_comp.ub()[dim]){
		}else{
			return false;
		}
	}

	IntervalVector box_its(goals.size());
	list<Vector> min_Points_Hyp_add;
	list<Vector> min_Points_Hyp_comp;

	//Obtain box intersect
	for(int dim = 0 ;dim < box_prime.size(); dim++){
		double box_lb, box_ub;
		box_lb = box_prime[dim].lb();

		if(box_prime[dim].ub() <= box_comp.ub()[dim]){
			box_ub = box_prime[dim].ub();
		}else{
			box_ub = box_comp.ub()[dim];
		}
		Interval aux(box_lb, box_ub);
		box_its[dim] = aux;
	}



	//Puntos que minimizan a las funciones en hyperplano dentro de la caja intersectada
	min_Points_Hyp_add.clear();
	min_Points_Hyp_comp.clear();
	for(int dim_goal=0; dim_goal<box_its.size(); dim_goal++){
		//Search minimum points of every function in the hyperplane to add inside the box intersection
		Vector point_min = find_minimum(box_its, hyperplane_coefficients_prime, dim_goal);
		if(point_min.size() != 0){
			min_Points_Hyp_add.push_back(point_min);
		}
		//Search minimum points of every function in the hyperplane to comp inside the box intersection
		Vector point_min2 = find_minimum(box_its, hyperplane_coefficients_comp, dim_goal);
		if(point_min2.size() != 0){
			min_Points_Hyp_comp.push_back(point_min2);
			min_Points_Hyp_add.push_back(point_min2);
		}

	}

	if(min_Points_Hyp_add.size() != 0){
	//Verificar que el hiperplano se  encuentra dominado
		double nb_hypadd_dominated     = 0;
		double nb_hypadd_not_dominated = 0;

		for(auto itr_minpoints=min_Points_Hyp_add.begin(); itr_minpoints!=min_Points_Hyp_add.end(); ++itr_minpoints){
			Vector min_point(*itr_minpoints);
			//verificar si c es mayor o menor en cada minimo de ambos hyperplanos
			double pos_hyp_comp = 0;
			double pos_hyp_add = 0;
			for(int dim_hyp = 0; dim_hyp<itr_minpoints->size(); dim_hyp++){
				pos_hyp_comp += hyperplane_coefficients_comp[dim_hyp]*min_point[dim_hyp];
				pos_hyp_add += hyperplane_coefficients_prime[dim_hyp]*min_point[dim_hyp];
			}
			pos_hyp_comp += hyperplane_coefficients_comp[goals.size()];
			pos_hyp_add  += hyperplane_coefficients_prime[goals.size()];

			if(pos_hyp_comp <= pos_hyp_add)
				nb_hypadd_dominated++;
			else
				nb_hypadd_not_dominated++;
		}

		//Check simple dominance between boxes with hyperplanes
		if(nb_hypadd_dominated == 0 && nb_hypadd_not_dominated != 0) return true;
	}

	return false;
}

int LoupFinderMOP::get_nds_simplex(){
	return nds_simplex;
}

int LoupFinderMOP::get_nb_sol(){
	return nb_sol;
}

} /* namespace ibex */
