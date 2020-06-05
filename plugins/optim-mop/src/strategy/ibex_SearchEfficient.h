//============================================================================
//                                  I B E X
// File        : ibex_Optimizer.h
// Author      : Matias Campusano, Damir Aliquintui, Ignacio Araya
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Sep 24, 2017
// Last Update : Sep 24, 2017
//============================================================================

#ifndef __IBEX_SEARCHEFFICIENT_H__
#define __IBEX_SEARCHEFFICIENT_H__

#include "ibex_Ctc.h"
#include "ibex_Bsc.h"
#include "ibex_LoupFinderMOP.h"
#include "ibex_CtcKhunTucker.h"
#include "ibex_DistanceSortedCellBufferMOP.h"
#include "ibex_pyPlotter.h"
#include "ibex_PFunction.h"
#include "ibex_NDS.h"
#include "ibex_OptimizerMOP.h"
#include "vector.h"

#include <set>
#include <map>
#include <list>
#include <stack>
#include <math.h>
#include <iostream>
#include <unistd.h>

#include "ibex_BxpMOPData.h"
//#include "ibex_DistanceSorted.h"

using namespace std;
namespace ibex {

  struct manhattan{
  	static double distance(const Cell* c){
  		int n = c->box.size()-2;
      double dist = 0.0;
      if(mode == 0 || mode == 1 )
         dist += std::max(efficient[0]-c->box[n].lb(),0.0);
      if(mode == 0 || mode == 2 )
         dist += std::max(efficient[1]-c->box[n+1].lb(),0.0);

  		return dist;
  	}


  	bool operator() (const Cell* c1, const Cell* c2){
  		int n = c1->box.size();
  		double dist1=distance(c1);
  		double dist2=distance(c2);

  		return dist1>dist2;

  	}

  	static double* efficient;
    static int mode;
  };



/**
 * \brief Global biObjetive Optimizer (ibexMOP).
 *
 * This class is an implementation of a global optimization algorithm for biObjective problems
 * described in https://github.com/INFPUCV/ibex-lib/tree/master/plugins/optim-mop by Araya et al.
 *
 * \remark In all the comments of this class, "NDS" means "Non Dominated Set"
 * related to the objectives f1 and f2
 */
class SearchEfficient{

public:

	/**
	 * \brief Return status of the optimizer
	 *
	 * See comments for optimize(...) below.
	 */
	typedef enum {SUCCESS, INFEASIBLE, NO_FEASIBLE_FOUND, UNBOUNDED_OBJ, TIME_OUT, UNREACHED_PREC} Status;

	//typedef enum {POINTS, SEGMENTS, HAMBURGER, /* splitting strategies */ MIDPOINT, MAXDIST, ALL} Mode;


	/**
	 *  \brief Create an optimizer.
	 *
	 * Inputs:
	 *   \param n        		   number of variables of the <b>original system</b>
	 *   \param f1	       		   the objective function f1
     *	 \param f2       		   the objective function f2
	 *   \param ctc       		   contractor for <b>extended<b> boxes (of size n+2)
	 *   \param bsc       		   bisector for <b>extended<b> boxes (of size n+2)
	 *   \param buffer    		   buffer for <b>extended<b> boxes (of size n+2)
	 *   \param finder    		   the finder of ub solutions
	 *   \param eps	      		   the required precision
	 *

	 *
	 * \warning The optimizer relies on the contractor \a ctc to contract the domain of the goal variable
	 *          and increase the uplo. If this contractor never contracts this goal variable,
	 *          the optimizer will only rely on the evaluation of f and will be very slow.
	 *
	 * We are assuming that the objective variables are n and n+1
	 *
	 */
	SearchEfficient(int n, const Function &f1,  const Function &f2,
			Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder,
			 double eps=default_eps, double rel_eps=0.0);

	/**
	 * \brief Delete *this.
	 */
	virtual ~SearchEfficient();



	void pre_optimize(const IntervalVector& init_box, Cell* root);

	/**
	 * \brief Run the optimization.
	 *
	 * \param init_box             The initial box
	 *
	 * \return SUCCESS             if the global minimum (with respect to the precision required) has been found.
	 *                             In particular, at least one feasible point has been found, less than obj_init_bound, and in the
	 *                             time limit.
	 **
	 *         TIMEOUT             if time is out.
	 *
	 */
  typedef enum{EFFICIENT, MINF1, MINF2}opt_mode;
	Status optimize(const IntervalVector& init_box, opt_mode mode=EFFICIENT);



	/* =========================== Output ============================= */

	/**
	 * \brief Displays on standard output a report of the last call to optimize(...).
	 *
	 * Information provided:
	 * <ul><li> total running time
	 *     <li> totl number of cells (boxes)
	 *     <li> total number of best non-dominated solutions found
	 *     <li> the set of non-dominated solutions found
	 * </ul>
	 */
	void report(bool verbose=true){
		if (!verbose) {
     cout << endl 	<< "time 	#nodes 		|Y|" << endl;
		 //cout << get_time() << " " << get_nb_cells() << " " << ndsH.size() << endl;
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

	/*cout << "\033[0m" << endl;

	cout << " cpu time used: " << get_time() << "s." << endl;
	cout << " number of cells: " << get_nb_cells() << endl;

	cout << " number of solutions: "  << ndsH.size() << endl;
	for(auto ub : ndsH.NDS2)
		 cout << "(" << ub.first[0] << "," << ub.first[1] << ")"  << endl;*/


}



	/**
	 * \brief Get the status.
	 *
	 * \return the status of last call to optimize(...).
	 */
	Status get_status() const;


	//std::set< point2 >& get_LB()  { return LB; }

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




	/*
    * \brief Get the box of potential solutions
    *
    * Inputs:
    *    \param v 				       a
    *    \param n		               a
    */

	static IntervalVector get_boxy(IntervalVector& v, int n){
		IntervalVector boxy(2);
		boxy[0]=v[n];
		boxy[1]=v[n+1];
		return boxy;
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

	/** Required precision for the envelope */
	double eps;




	/** Default precision: 0.01 */
	static const double default_eps;



	/**
	 * \brief Trace activation flag.
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

	//pair <double, double> y1ref;
	//pair <double, double> y2ref;

	//Interval y1refi;
	//Interval y2refi;

	/* ======== Some other parameters of the solver =========== */

	//if true: Save a file to be plotted by plot.py (default value= false).
	bool static _plot;

	//Min distance between two non dominated points to be considered, expressed as a fraction of eps (default value= 0.1)
	double static _min_ub_dist;

	//True: the solver uses the upper envelope of the cy contract for contraction
	static bool _cy_upper;

	//True: the solver uses the lower envelope of the cy contract in contraction
	static bool cy_contract_var;

	//True: the solver reduces the search spaces by reducing the NDS vectors in (eps, eps)
	static bool _eps_contract;

	//Termination criteria for the hamburger algorithm (dist < rh*ini_dist)
	static double _rh;

  //max distance of cells rejected by eps-close
   double max_dist_eps;

	//NDS mode: POINTS or SEGMENTS
	//Mode nds_mode;
	//Mode split_mode;


	/**
	 * \brief Evaluate the goal in the point x
	 */
	static Interval eval_goal(const Function& goal, const IntervalVector& x, int n);

	/**
	 * \brief Gradient of the goal in the point x
	 */
	static IntervalVector deriv_goal(const Function& goal, const IntervalVector& x, int n);

	/*double distance(const Cell* c){
		return NDS_seg::distance(c);
	}*/

	double efficient[2];

protected:


	/**
	 * \brief return the first and last points dominating the lb of the box
	 */
	static list<pair <double,double> > extremal_non_dominated(const IntervalVector& box);



	/**
	 * \brief Contract and bound procedure for processing a box.
	 *
	 * <ul>
	 * <li> contract the cell using #dominance_peeler() and #discard_generalized_monotonicty_test(),
	 * <li> contract with the contractor ctc,
	 * </ul>
	 *
	 */
	void contract_and_bound(Cell& c, const IntervalVector& init_box, opt_mode mode=EFFICIENT);

	//SearchEfficient* search_left_box( const int n,
	//								double inf_x, double sup_y, double inf_y, double sup_x);

	//SearchEfficient* search_right_box(const IntervalVector& init_box, const int n,
	//								double sup_y, double inf_x, double inf_y, double sup_x);

    /**
     *  \brief returns true if the facet orthogonal to the i direction of the box is feasible.
     *
     *  See #discard_generalized_monotonicty_test()
     */
    bool is_inner_facet(IntervalVector box, int i, Interval bound){
    	box.resize(n);
    	box[i]=bound;
    	return finder.norm_sys.is_inner(box);
    }


    /**
     * \brief Implements the method for discarding boxes proposed in [J. Fernandez and B. Toth,
     * "Obtaining the efficient set of nonlinear biobjective optimiza-tion  problems  via  interval
     * branch-and-bound  methods" (2009)]
     */
	void discard_generalized_monotonicty_test(IntervalVector& box, const IntervalVector& initbox);



	/**
	 * \brief Main procedure for updating the NDS.
	 * <ul>
	 * <li> finds two points xa and xb in a polytope using the #LoupFinderMOP,
	 * <li> finds a segment passing over the points (f1(x),f2(x)) in the segment xa-xb,
	 * <li> add the segment to NDS
	 * </ul>
	 */
	bool upper_bounding(const IntervalVector& box, const IntervalVector& init_box, opt_mode mode=EFFICIENT);


	/* ======== Some other parameters of the solver =========== */

	/**
	 * min feasible value found for each objective
	 */
    pair <double, double> y1_ub, y2_ub;

	/* Remember return status of the last optimization. */
	Status status;



	/** The current non-dominated set sorted by increasing x */
	//static map< pair <double, double>, IntervalVector, sorty2 > NDS2;

	/* CPU running time of the current optimization. */
	double time;

	/** Number of cells pushed into the heap (which passed through the contractors) */
	int nb_cells;

};

inline SearchEfficient::Status SearchEfficient::get_status() const { return status; }

inline double SearchEfficient::get_time() const { return time; }

inline double SearchEfficient::get_nb_cells() const { return nb_cells; }




} // end namespace ibex

#endif // __IBEX_OPTIMIZER_H__
