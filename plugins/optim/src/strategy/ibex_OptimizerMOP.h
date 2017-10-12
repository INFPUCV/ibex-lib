//============================================================================
//                                  I B E X
// File        : ibex_Optimizer.h
// Author      : Matias Campusano, Damir Aliquintui, Ignacio Araya
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Sep 24, 2017
// Last Update : Sep 24, 2017
//============================================================================

#ifndef __IBEX_OPTIMIZERMOP_H__
#define __IBEX_OPTIMIZERMOP_H__

#include "ibex_Ctc.h"
#include "ibex_Bsc.h"
#include "ibex_LoupFinderMOP.h"
#include "ibex_CellBufferOptim.h"
//#include "ibex_EntailedCtr.h"
#include "ibex_CtcKhunTucker.h"

#include <set>
#include <map>
#include <list>

using namespace std;
namespace ibex {

/**
 * \defgroup optim IbexOpt
 */

/**
 * \ingroup optim
 *
 * \brief Global MOP Optimizer.
 *
 */
class OptimizerMOP {

public:

	/**
	 * \brief Return status of the optimizer
	 *
	 * See comments for optimize(...) below.
	 */
	typedef enum {SUCCESS, INFEASIBLE, NO_FEASIBLE_FOUND, UNBOUNDED_OBJ, TIME_OUT, UNREACHED_PREC} Status;

	/**
	 *  \brief Create an optimizer.
	 *
	 * Inputs:
	 *   \param n        - number of variables or the <b>original system</b>
	 *   \param ctc      - contractor for <b>extended<b> boxes (of size n+2)
	 *   \param bsc      - bisector for <b>extended<b> boxes (of size n+2)
	 *   \param finder   - upper-bounding procedure for the original system (n-sized boxes)
	 *   \param buffer   - buffer for <b>extended<b> boxes (of size n+2)
	 *   \param goal_var - index of the first goal variable in an extended box.
   *   	 \param goal_var2- index of the second goal variable in an extended box.
	 *
	 * And optionally:
	 *   \param eps_x         - absolute precision for the boxes (bisection control)
	 *
	 * \warning The optimizer relies on the contractor \a ctc to contract the domain of the goal variable
	 *          and increase the uplo. If this contractor never contracts this goal variable,
	 *          the optimizer will only rely on the evaluation of f and will be very slow.
	 *
	 * We are assuming that the objective variables are n and n+1
	 *
	 */
	OptimizerMOP(int n, const Array<NumConstraint>& ctcs, const Function &f1,  const Function &f2,
			Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder, double eps_x=default_eps_x, double eps_z=default_eps_z);

	/**
	 * \brief Delete *this.
	 */
	virtual ~OptimizerMOP();

	/**
	 * \brief Run the optimization.
	 *
	 * \param init_box             The initial box
	 * \param obj_init_bound       (optional) can be set when an initial upper bound of the objective minimum is known a priori.
	 *                             This bound can be obtained, e.g., by a local solver. This is equivalent to (but more practical
	 *                             than) adding a constraint f(x)<=obj_init_bound.
	 *
	 * \return SUCCESS             if the global minimum (with respect to the precision required) has been found.
	 *                             In particular, at least one feasible point has been found, less than obj_init_bound, and in the
	 *                             time limit.
	 *
	 *         INFEASIBLE          if no feasible point exist less than obj_init_bound. In particular, the function returns INFEASIBLE
	 *                             if the initial bound "obj_init_bound" is LESS than the true minimum (this case is only possible if
	 *                             goal_abs_prec and goal_rel_prec are 0). In the latter case, there may exist feasible points.
	 *
	 *         NO_FEASIBLE_FOUND   if no feasible point could be found less than obj_init_bound. Contrary to INFEASIBLE,
	 *                             infeasibility is not proven here. Warning: this return value is sensitive to the abs_eps_f and
	 *                             rel_eps_f parameters. The upperbounding makes the optimizer only looking for points less than
	 *                             min { (1-rel_eps_f)*obj_init_bound, obj_init_bound - abs_eps_f }.
	 *
	 *         UNBOUNDED_OBJ       if the objective function seems unbounded (tends to -oo).
	 *
	 *         TIMEOUT             if time is out.
	 *
	 *         UNREACHED_PREC      if the search is over but the resulting interval [uplo,loup] does not satisfy the precision
	 *                             requirements. There are several possible reasons: the goal function may be too pessimistic
	 *                             or the constraints function may be too pessimistic with respect to the precision requirement
	 *                             (which can be too stringent). This results in tiny boxes that can neither be contracted nor
	 *                             used as new loup candidates. Finally, the eps_x parameter may be too large.
	 *
	 */
	Status optimize(const IntervalVector& init_box);

	/* =========================== Output ============================= */

	/**
	 * \brief Displays on standard output a report of the last call to optimize(...).
	 *
	 * Information provided:
	 * <ul><li> interval of the cost  [uplo,loup]
	 *     <li> the best feasible point found
	 *     <li> total running time
	 *     <li> total number of cells (~boxes) created during the exploration
	 * </ul>
	 */
	void report(bool verbose=true);

	/**
	 * \brief Get the status.
	 *
	 * \return the status of last call to optimize(...).
	 */
	Status get_status() const;

	/**
	 * \brief Get the "UB" set of the pareto front.
	 *
	 * \return the UB of the last call to optimize(...).
	 */
	map< pair <double, double>, Vector > get_UB() const;


	/**
	 * \brief Get the time spent.
	 *
	 * \return the total CPU time of last call to optimize(...)
	 */
	double get_time() const;

	/**
	 * \brief Get the number of cells.
	 *
	 * \return the number of cells generated by the last call to optimize(...).
	 */
	double get_nb_cells() const;

	/**
	 * \brief Get the number of solutions.
	 */
	int get_nb_sol(){
		return Sout.size();
	}

	/* =========================== Settings ============================= */

	/**
	 * \brief Number of variables.
	 */
	const int n;

	/**
	 * \brief Objective functions
	 * Functions have the form: f1 - z1  and f2 - z2. Thus, in order to
	 * evaluate them we have to set z1 and z2 to [0,0].
	 */
	const Function& goal1;
	const Function& goal2;

	/**
	 * \brief Constraints
	 */
	const Array<NumConstraint>& ctrs;

	/**
	 * \brief Contractor for the extended system.
	 *
	 * The extended system:
	 * (y=f(x), g_1(x)<=0,...,g_m(x)<=0).
	 */
	Ctc& ctc;

	/**
	 * \brief Bisector.
	 *
	 * Must work on extended boxes.
	 */
	Bsc& bsc;

	/**
	 * Cell buffer.
	 */
	CellBuffer& buffer;

	/**
	 * \brief LoupFinder
	 */
	LoupFinderMOP& finder;

	/** Precision (bisection control constraints) */
	const double eps_x;

	/** Default bisection precision: 1e-5 */
	static const double default_eps_x;

	/** Precision (bisection control objective functions) */
	const double eps_z;

	/** Default bisection precision: 1e-5 */
	static const double default_eps_z;

	/**
	 * \brief Trace activation flag.
	 *
	 * The value can be fixed by the user.
	 * - 0 (by default): nothing is printed
	 * - 1:              prints every loup/uplo update.
	 * - 2:              prints also each handled node (warning: can generate very
	 *                   long trace).
	 */
	int trace;

	/**
	 * \brief Time limit.
	 *
	 * Maximum CPU time used by the strategy.
	 * This parameter allows to bound time consumption.
	 * The value can be fixed by the user.
	 */
	double timeout;


protected:

	/**
	 * \brief Main procedure for processing a box.
	 *
	 * <ul>
	 * <li> contract the cell box and try to find a new loup (see contract_and_bound)
	 * <li> push the cell into the buffer or delete the cell in case of empty box detected.
	 * </ul>
	 *
	 */
	void handle_cell(Cell& c, const IntervalVector& init_box);

	/**
	 * \brief Contract and bound procedure for processing a box.
	 *
	 * <ul>
	 * <li> contract the cell's box w.r.t the "loup",
	 * <li> contract with the contractor ctc,
	 * <li> search for a new loup,
	 * <li> (optional) call the first order contractor
	 * </ul>
	 *
	 */
	void contract_and_bound(Cell& c, const IntervalVector& init_box);


	/**
	 * \brief compute the max distance between the line z1 + a*z2 = w_lb and the UB step function
	 */
	double max_distance(double w_lb, double a, Interval z1, Interval z2){
		double max_dist=-1.0;

		//cout << z1 << ";" << z2 << endl;
		map< pair <double, double>, Vector >::iterator it = UB.begin();
		for(;it!=UB.end(); ){
			pair <double, double> p = it->first; it++;
			if(it==UB.end()) break;
			pair <double, double> p2 = it->first;

			//la solucion esta dentro de la caja
			if( (z1.contains(p.first) && z2.contains(p.second)) || (z1.contains(p2.first) && z2.contains(p2.second)) ){
				//cout << "(" << p.first << "," << p.second << ") ; (" << p2.first << "," << p2.second << ")" << endl;

				pair <double, double> pmax= make_pair( std::min(p2.first,z1.ub()), std::min(p.second,z2.ub()));

				double dist= std::min( (pmax.first + a*pmax.second) - w_lb, std::min(pmax.first-z1.lb(),  pmax.second-z2.lb()));
				//cout << dist << endl;
				if(dist > max_dist) max_dist=dist;
			}
		}

		if(max_dist == -1.0) return std::min( z1.ub() + a*z2.ub() - w_lb, std::min(z1.diam(), z2.diam()));
		else return max_dist;
	}

	/**
	 * TODO: Funcion que genera segmentos LB una vez terminada la busqueda
	 */
	void generate_LB(){
		std::list< IntervalVector >::iterator it = Sout.begin();

		for(;it!=Sout.end(); it++){
		//	update_LB()
		}
	}



	/**
	 * TODO: (DA, MC) Funcion que agrega segmentos al set LB
	 * \brief Update the set LB by adding a segment
	 */
	void update_LB(const pair<double, double> p1, const pair<double, double> p2){
        /** Idea:
		* 1. Encontrar segmentos de LB que corten el segmento p1-p2
		* 2. Encontrar puntos de interseccion
		* 3. Actualizar LB
		*/
	}

	/**
	 * \brief Update the set LB by adding two segments (p.first, +inf)--(p.first, p.second) and (p.first, p.second)--(+inf, p.second)
	 */
	void update_LB(const pair<double, double> p){
		update_LB(make_pair(p.first, POS_INFINITY), p);
		update_LB(p, make_pair(POS_INFINITY, p.second));
	}

	/**
	 * \brief Main procedure for updating the loup.
	 */
	bool update_UB(const IntervalVector& box, int n);

	/**
	 * \brief Check time is not out.
	 */
	void time_limit_check();


private:

	/**
	 * \brief Evaluate the goal in the point x
	 */
	Interval eval_goal(const Function& goal, Vector& x);

	/** Currently entailed constraints */
	//EntailedCtr* entailed;

	//!! warning: sys.box should be properly set before call to constructor !!
	//CtcKhunTucker kkt;

	/* Remember return status of the last optimization. */
	Status status;

	/** The set of final solutions */
	std::list< IntervalVector > Sout;

	/** The current upper bounds (f1(x), f2(x)) of the pareto front associated
	 * to its corresponding  point x
	 * If the loup-finder is rigorous, x may be a (non-degenerated) box. */
	map< pair <double, double>, Vector > UB;

	/**
	 * A set of points denoting the segments related to the lowerbound of the
	 * pareto front.
	 */
	std::set< pair <double, double> > LB;


	/** True if loup has changed in the last call to handle_cell(..) */
	//bool loup_changed;

	/* CPU running time of the current optimization. */
	double time;

	/** Number of cells pushed into the heap (which passed through the contractors) */
	int nb_cells;
};

inline OptimizerMOP::Status OptimizerMOP::get_status() const { return status; }

inline map< pair <double, double>, Vector > OptimizerMOP::get_UB() const { return UB; }

inline double OptimizerMOP::get_time() const { return time; }

inline double OptimizerMOP::get_nb_cells() const { return nb_cells; }


} // end namespace ibex

#endif // __IBEX_OPTIMIZER_H__
