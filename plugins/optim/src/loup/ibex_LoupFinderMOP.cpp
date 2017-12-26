/*
 * ibex_LoupFinderMOP.cpp
 *
 *  Created on: 5 oct. 2017
 *      Author: iaraya
 */

#include "ibex_LoupFinderMOP.h"
#include "ibex_PdcHansenFeasibility.h"
#include "ibex_FncActivation.h"


namespace ibex {

 double LoupFinderMOP::_weight2=0.0;

//TODO: remove this recipe for the argument of the max number of iterations of the LP solver
LoupFinderMOP::LoupFinderMOP(const System& sys, const Function& goal1, const Function& goal2, double eqeps) :
		sys(sys), norm_sys(sys,eqeps), lr(norm_sys,LinearizerXTaylor::RESTRICT),
		lp_solver(sys.nb_var, std::max((sys.nb_var)*3,LPSolver::default_max_iter)), goal1(goal1), goal2(goal2), has_equality(false) {

	if (sys.nb_ctr>0)
		// ==== check if the system contains equalities ====
		for (int i=0; i<sys.f_ctrs.image_dim(); i++) {
			if (sys.ops[i]==EQ) {
				(bool&) has_equality = true;
				break;
			}
		}

//	nb_simplex=0;
//	diam_simplex=0;
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

void LoupFinderMOP::find(const IntervalVector& box, list<Vector>& candidate_points, int nn) {

	int n=norm_sys.nb_var;

  if(nn>1 && (lp_solver.default_limit_diam_box.contains(box.max_diam()))){
	  lp_solver.clean_ctrs();
	  lp_solver.set_bounds(box);


    for(int i=0; i<2; i++){

    	IntervalVector box2(box);
    	box2.resize(n+2);
    	box2[n]=0.0; box2[n+1]=0.0;
		IntervalVector ig= (i==0)? (goal1.gradient(box2.mid())+ _weight2*goal2.gradient(box2.mid())) : (goal2.gradient(box2.mid())+ _weight2*goal1.gradient(box2.mid()));

		if (ig.is_empty()) // unfortunately, at the midpoint the function is not differentiable
			continue; // not a big deal: wait for another box...

		ig.resize(n);

		Vector g=ig.mid();

		// set the objective coefficient
		// TODO: replace with lp_solver.set_obj(g) when implemented
		for (int j=0; j<n; j++)
			lp_solver.set_obj_var(j,g[j]);

		int count = lr.linearize(box,lp_solver);

		if (count==-1) {
			lp_solver.clean_ctrs();
			return;
		}

		//lp_solver.write_file("holo.txt");


		LPSolver::Status_Sol stat = lp_solver.solve();

		if (stat == LPSolver::OPTIMAL) {
			//the linear solution is mapped to intervals
			Vector loup_point(n);
			lp_solver.get_primal_sol(loup_point);


			//std::cout << " simplex result " << loup_point << std::endl;

			//std::cout << box << endl;

			//correct the point
			for(int i=0;i<box.size();i++){
				if(loup_point[i] < box[i].lb())  loup_point[i] = box[i].lb();
				if(loup_point[i] > box[i].ub())  loup_point[i] = box[i].ub();
			}

            if(candidate_points.size()==0 || candidate_points.front()!=loup_point )
			   candidate_points.push_back(loup_point);
		}
    }
	}
    //we add the feasible inner points
    if(candidate_points.size()==2 && nn>2){
    	double step = 1.0/((double)nn-1.0);

    	double r=0.0;
    	for(int i=0; i<nn; i++, r+=step){
    		if(i==nn-1) r=1.0;
    		list<Vector>::iterator it = candidate_points.begin();
    		Vector vec1 = *it; it++;
    		Vector vec2 = *it;
    		Vector vec3 = (1.0-r)*vec1 + r*vec2;
			  for(int i=0;i<box.size();i++){
				  if(vec3[i] < box[i].lb())  vec3[i] = box[i].lb();
				  if(vec3[i] > box[i].ub())  vec3[i] = box[i].ub();
			  }

    		candidate_points.push_back(vec3);


    	}
    	candidate_points.pop_front();
    	candidate_points.pop_front();

    }else{

		Vector vec=box.mid();
		if(norm_sys.is_inner(vec)) {
			candidate_points.push_back(vec);
		}
	}


}

} /* namespace ibex */
