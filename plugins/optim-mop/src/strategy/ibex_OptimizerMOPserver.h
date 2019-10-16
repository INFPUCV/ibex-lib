/*
 * ibex_OptimizerMOPserver.h
 *
 *  Created on: Sep 25, 2019
 *      Author: iaraya
 */

#ifndef __IBEX_OPTIMIZERMOP_S_H__
#define __IBEX_OPTIMIZERMOP_S_H__

#include "ibex_OptimizerMOP.h"

namespace ibex {

class OptimizerMOP_S : public OptimizerMOP {
public:


    typedef enum {STAND_BY_FOCUS, STAND_BY_SEARCH, REACHED_PRECISION, SEARCH, FOCUS_SEARCH, FINISHED} ServerStatus;

//	OptimizerMOP_S(int n, const Function &f1,  const Function &f2,
//			Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder,
//			Mode nds_mode=POINTS, Mode split_mode=MIDPOINT, double eps=default_eps, double rel_eps=0.0);

	virtual ~OptimizerMOP_S() { }


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
	virtual Status optimize(const IntervalVector& init_box);

	void write_status(double rel_prec);

	/*
    * \brief Print the region of solutions
    *
    *
    *
    * Inputs:
    *    \param cells 				   a
    *    \param paused_cells		   a
    *    \param focus 				   a
    */

 	void write_envelope(set<Cell*>& cells, set<Cell*>& paused_cells, IntervalVector& focus);


	/*
    * \brief Read the instruccions of work
    *
    * after reading a instruction, this delete it from the file
    *
    * Inputs:
    *    \param cells 				   a
    *    \param paused_cells		   a
    *    \param focus 				   a
    */

	void read_instructions(set<Cell*>& cells, set<Cell*>& paused_cells, IntervalVector& focus);


	/*
    * \brief Update the focus of solution
    *
    * This take in count the hull of the region and the found solutions
    *
    * Inputs:
    *    \param cells 				   a
    *    \param paused_cells		   a
    *    \param focus 				   a
    */

	void update_focus(set<Cell*>& cells, set<Cell*>& paused_cells, IntervalVector& focus);


	void zoom(bool out, set<Cell*>& cells, set<Cell*>& paused_cells, IntervalVector& focus, ifstream& myfile);

	void get_solution(ifstream& myfile);

  //server variables
	static string instructions_file;
	static string output_file;

	ServerStatus sstatus;

};

} /* namespace ibex */

#endif /* PLUGINS_OPTIM_MOP_SRC_STRATEGY_IBEX_OPTIMIZERMOPSERVER_H_ */
