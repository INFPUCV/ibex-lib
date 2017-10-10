/*
 * ibex_LoupFinderMOP.cpp
 *
 *  Created on: 5 oct. 2017
 *      Author: iaraya
 */

#include "ibex_LoupFinderMOP.h"

namespace ibex {

//TODO: remove this recipe for the argument of the max number of iterations of the LP solver
LoupFinderMOP::LoupFinderMOP(const System& sys, const Function& goal1, const Function& goal2) : sys(sys), lr(sys,LinearizerXTaylor::RESTRICT),
		lp_solver(sys.nb_var, std::max((sys.nb_var)*3,LinearSolver::default_max_iter)), goal1(goal1), goal2(goal2) {
//	nb_simplex=0;
//	diam_simplex=0;
}

void LoupFinderMOP::find(const IntervalVector& box, list<Vector>& feasible_points) {

	if (!(lp_solver.default_limit_diam_box.contains(box.max_diam())))
		return;

	int n=sys.nb_var;

	lp_solver.clean_ctrs();

	lp_solver.set_bounds(box);



    for(int i=0; i<2; i++){

    	IntervalVector box2(box);
    	box2.resize(n+2);
    	box2[n]=0.0; box2[n+1]=0.0;
		IntervalVector ig= (i==0)? goal1.gradient(box2.mid()) : goal2.gradient(box2.mid());

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


		LinearSolver::Status_Sol stat = lp_solver.solve();

		if (stat == LinearSolver::OPTIMAL) {
			//the linear solution is mapped to intervals and evaluated
			Vector loup_point(n);
			lp_solver.get_primal_sol(loup_point);


			//std::cout << " simplex result " << loup_point << std::endl;

			//std::cout << box << endl;

			if (!box.contains(loup_point)) continue;

			feasible_points.push_back(loup_point);
		}
    }

    //we add the feasible middle point

    if(feasible_points.size()==2){
    	list<Vector>::iterator it = feasible_points.begin();
    	Vector vec1 = *it; it++;
    	Vector vec2 = *it;
    	Vector vec3 = (vec1+vec2);
    	vec3 *= 0.5;
    	if (box.contains(vec3)) feasible_points.push_back(vec3);
    }

    return;
}

} /* namespace ibex */
