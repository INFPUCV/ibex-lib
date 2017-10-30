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
#include "ibex_CellSet.h"
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
	 *   \param f1	     - the objective function f1
     *	 \param f2       - the objective function f2
	 *   \param ctc      - contractor for <b>extended<b> boxes (of size n+2)
	 *   \param bsc      - bisector for <b>extended<b> boxes (of size n+2)
	 *   \param buffer   - buffer for <b>extended<b> boxes (of size n+2)
	 *   \param finder   - the finder of ub solutions
	 *   \param eps	     - the required precision
	 *

	 *
	 * \warning The optimizer relies on the contractor \a ctc to contract the domain of the goal variable
	 *          and increase the uplo. If this contractor never contracts this goal variable,
	 *          the optimizer will only rely on the evaluation of f and will be very slow.
	 *
	 * We are assuming that the objective variables are n and n+1
	 *
	 */
	OptimizerMOP(int n, const Array<NumConstraint>& ctcs, const Function &f1,  const Function &f2,
			Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder, double rel_eps=default_rel_eps, double abs_eps=default_abs_eps);

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

	void plot(Cell* current);

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
	map< pair <double, double>, IntervalVector >& get_UB()  { return UB; }


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

	/** Precision (bisection control objective functions) */
	double abs_eps;

	const double rel_eps;

	/** Default absolute precision: 1e-7 */
	static const double default_abs_eps;

	/** Default relative precision: 0.1 */
	static const double default_rel_eps;

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


	/**
	 * \brief returns true if the box+z1 + a*z2 > w_lb is dominated by the ub_set
	 */
	static double distance2(const Cell* c){
		//if(c->get<CellBS>().ub_distance != POS_INFINITY) return c->get<CellBS>().ub_distance;
		double max_dist=NEG_INFINITY;
		if(UB.size()==2) return POS_INFINITY;

		int n=c->box.size();

		Interval z1 = c->box[n-2];
		Interval z2 = c->box[n-1];
		double a = c->get<CellBS>().a;
		double w_lb = c->get<CellBS>().w_lb;

		//TODO: optimize this
		map< pair <double, double>, IntervalVector >::iterator it = UB.begin();


		for(;it!=UB.end(); ){
			pair <double, double> p = it->first; it++;
			if(it==UB.end()) break;
			pair <double, double> p2 = it->first;

			pair <double, double> pmax= make_pair(p2.first, p.second);
			//cout << "pmax: (" << pmax.first <<"," << pmax.second << ")" << endl;



			//el punto esta dentro de la zona de interes
			if(pmax.first > z1.lb() && pmax.second > z2.lb()){
				double dist = std::min (pmax.first - z1.lb(), pmax.second - z2.lb());
				//here we add the distance to the line
			  dist = std::min (dist, (Interval(pmax.first) + Interval(a)*Interval(pmax.second) - Interval(w_lb)).ub() );

				if(dist > max_dist) max_dist=dist;

			}
		}

		return max_dist;
	}

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



private:

	/**
	 * \brief Evaluate the goal in the point x
	 */
	Interval eval_goal(const Function& goal, IntervalVector& x);


	/** Currently entailed constraints */
	//EntailedCtr* entailed;

	//!! warning: sys.box should be properly set before call to constructor !!
	//CtcKhunTucker kkt;

	/* Remember return status of the last optimization. */
	Status status;

	/** The set of final solutions */
	std::list< IntervalVector > Sout;

  /** The cells in the buffer for plotting
	 * the set should be updated each time the real buffer is popped
	 * and pushed.
	 */
	set< Cell* > buffer_cells;

	/** The current upper bounds (f1(x), f2(x)) of the pareto front associated
	 * to its corresponding  point x
	 */
	static map< pair <double, double>, IntervalVector > UB;

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

inline double OptimizerMOP::get_time() const { return time; }

inline double OptimizerMOP::get_nb_cells() const { return nb_cells; }




} // end namespace ibex

#endif // __IBEX_OPTIMIZER_H__
