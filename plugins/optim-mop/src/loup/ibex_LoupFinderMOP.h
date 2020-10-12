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
#include "ibex_LoupFinder.h"
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
class LoupFinderMOP : public LoupFinder{

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

	//LoupFinderMOP(const System& sys, const Function& goal1, const Function& goal2, const Array<const Function> goals, double eqeps=NormalizedSystem::default_eps_h, int nb_sol=50);

	LoupFinderMOP(const System& sys, const Array<const Function> goals, double eqeps=NormalizedSystem::default_eps_h, int nb_sol=50, int nds_simplex=1);


	/**
	 * \brief Correct the solution p by using a Hansen feasibility test with eps-inflation
	 */
	bool ub_correction(Vector p, IntervalVector& res);

	/**
	 * \brief Find a candidate solution for the non-dominated set.
	 *
	 * \param box        - the box where the solution is searched
	 * \param loup_point - not used
	 * \param loup       - not used
	 * \return             a candidate solution <x,f(x)> (may be not feasible)
	 * \throws             NotFound in case of failure.
	 */
	virtual std::pair<IntervalVector, double> find(const IntervalVector& box, const IntervalVector& loup_point, double loup=POS_INFINITY);



	virtual std::pair<IntervalVector, double> find2(const IntervalVector& box, const IntervalVector& loup_point, double loup=POS_INFINITY);

	/**
	 * \brief True if equalities are accepted.
	 *
	 * This function is abstract and may be overriden in the subclass.
	 *
	 * By default: return true.
	 */
	virtual bool rigorous() const {
		return false;
	}

	virtual void clear() {
		lp_solver.clean_ctrs();
		phase=0;
	}

	/**
	 * \brief Delete this.
	 */
	virtual ~LoupFinderMOP() {}


	/**
	 * \brief The relaxed NLP problem for finding feasible points
	 */
	const NormalizedSystem norm_sys;

	/**
	 * Weight of the secondary objective function for the linear program
	 */
	double static _weight2;


	bool isSaved(Vector loup_point);

	void add_external_vertex(Vector base_point, list<Vector> simplex_points, list<Vector>& external_vertex, IntervalVector box_y);

	void add_nondominated_point(pair<Vector,Vector> point, list<pair<Vector,Vector>>& list);

	bool is_dominated(const Vector& y_prime, const Vector& y);

	int get_nds_simplex();

	int get_nb_sol();

	//testing another method
	void add_external_vertex2(pair<Vector,Vector> base_point, list<pair<Vector,Vector>> list_points, list<Vector> & external_vertex, IntervalVector box_y);

	Vector find_minimum(IntervalVector box, Vector hyperplane_coefficients, int dim_goal);

	bool is_dominated_hyperplane_add(list< pair<IntervalVector, Vector> > upper_envelope, IntervalVector box_y, Vector hyperplane_coefficients);

    bool is_dominated_hyperplane(IntervalVector box, Vector hyperplane_coefficients, IntervalVector box_prime, Vector hyperplane_coefficients_prime);

    void print_hyperplane_intersection(IntervalVector caja_comp, Vector hyp_comp, IntervalVector caja_add, Vector hyp_add, IntervalVector caja_int, list<Vector> puntos_min_add, list<Vector> puntos_min_comp);

protected:

	/**
	 * \brief The real NLP problem.
	 */
	const System& sys;



	/**
	 * \brief Objective function fi
	 * Functions have the form: fi - zi. Thus, in order to
	 * evaluate them we have to set zi [0,0].
	 */

	/*
	 * Array of objective funcion fn
	 */
	const Array<const Function> goals;


	/**
	 * \brief True iff there is an equality.
	 */
	const bool has_equality;

	/** Linearization technique. */
	LinearizerXTaylor lr;

	/** linear solver */
	LPSolver lp_solver;

private:

	//number of solutions to find
	int nb_sol;

	//number of solution to be added
	int nds_simplex;

	//0: means the first solution of the polytope (or the midpoint)
	//1: means the second solution of the polytope
	//>=2: means the solutions in the line
	int phase;
	Vector vec1; //the first solution of the poltytope
	Vector vec2; //the second solution of the polytope


	//The index of the objective function and its solution is saved
	list <pair<Vector, int>> solutions; //Solution of each objective function in the polytope



};

} /* namespace ibex */

#endif /* OPTIM_SRC_LOUP_IBEX_LOUPFINDERMOP_H_ */
