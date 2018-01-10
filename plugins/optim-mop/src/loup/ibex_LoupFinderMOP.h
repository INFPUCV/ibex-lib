/*
 * ibex_LoupFinderMOP.h
 *
 *  Created on: 5 oct. 2017
 *      Author: iaraya
 */

#ifndef OPTIM_SRC_LOUP_IBEX_LOUPFINDERMOP_H_
#define OPTIM_SRC_LOUP_IBEX_LOUPFINDERMOP_H_

#include "ibex_Vector.h"
#include "ibex_System.h"
#include "ibex_LinearizerXTaylor.h"
#include "ibex_NormalizedSystem.h"

#include <list>

using namespace std;

namespace ibex {

/**
 * \ingroup optim
 *
 * \brief Upper-bounding algorithm based on XTaylor restriction.
 *
 * The algorithm builds an inner (feasible) polytope inside the
 * current box (see #LinearizerXTaylor) and then minimizes the TWO
 * linear approximation of the goal functions on this polytope via
 * a LP solver. The resulting points are verified a posteriori to
 * be feasible (wrt nonlinear constraint).
 *
 * \note Only works with inequality constraints.
 */
class LoupFinderMOP {

public:

	/**
	 * \brief Create the algorithm for a given system.
	 *
	 * The system is an inequality system of constraints.
	 * Goal functions have the form: f1 - z1  and f2 - z2.
	 *
	 * \param sys         - The NLP problem.
	 * \param goal1
	 * \param goal2
	 */
	LoupFinderMOP(const System& sys, const Function& goal1, const Function& goal2, double eqeps=NormalizedSystem::default_eps_h);


	/**
	 * \brief Correct the solution p by using a Hansen feasibility test with eps-inflation
	 */
	bool ub_correction(Vector p, IntervalVector& res);

	/*
	 * find up to n feasible points inside a inner-polytope
	 */
	void find(const IntervalVector& box, list<Vector>& feasible_points, int n=2);

	/**
	 * \brief The real NLP problem.
	 */
	const System& sys;

	/**
	 * \brief The relaxed NLP problem for finding feasible points
	 */
	const NormalizedSystem norm_sys;

	/**
	 * \brief Objective function f1
	 * Functions have the form: f1 - z1  and f2 - z2. Thus, in order to
	 * evaluate them we have to set z1 and z2 to [0,0].
	 */
	const Function& goal1;

	/**
	 * \brief Objective function f2
	 */
	const Function& goal2;

	/**
	 * Weight of the secondary objective function for the linear program
	 */
	double static _weight2;
protected:

	/**
	 * \brief True iff there is an equality.
	 */
	const bool has_equality;

	/** Linearization technique. */
	LinearizerXTaylor lr;

	/** linear solver */
	LPSolver lp_solver;


};

} /* namespace ibex */

#endif /* OPTIM_SRC_LOUP_IBEX_LOUPFINDERMOP_H_ */
