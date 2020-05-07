/*
 * ibex_OptimizerMOPserver.h
 *
 *  Created on: Sep 25, 2019
 *      Author: iaraya
 */

#ifndef __IBEX_OPTIMIZERMOPINTERACTIVE_S_H__
#define __IBEX_OPTIMIZERMOPINTERACTIVE_S_H__

#include "ibex_OptimizerMOP.h"
#include "ibex_Timer.h"

#ifndef cdata
#define cdata ((BxpMOPData*) c->prop[BxpMOPData::id])
#endif

namespace ibex {

class OptimizerMOP_I : public OptimizerMOP {
public:

	OptimizerMOP_I(int n, const Function &f1,  const Function &f2,
			Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder,
			Mode nds_mode=POINTS, Mode split_mode=MIDPOINT);

	virtual ~OptimizerMOP_I() { }


	/**
	 * \brief Load the initial_box (x) in the optimizer and returns the initial y
	 */
	virtual IntervalVector load(const IntervalVector& init_box);

  /**
	 * \brief Load the nodes in the search tree and returns the image hull y
	 */
  virtual IntervalVector load(const IntervalVector& init_box, string filename);

  /**
  * \brief perform maxiter iterations and returns SUCCESS (if the search is )
  */

  typedef enum {NONE, READY, STOPPED, FINISHED} IStatus;

  /*Funciones para interactuar con la api */

  void save_state_in_file(string filename);

  void load_state_from_file(string filename, const IntervalVector& init_box);

	IStatus run(int maxiter=5, double eps=1e-1);

  void update_refpoint(Vector& refpoint);

 	void write_envelope(string output_file);

  /****************************************/

  IntervalVector initbox;

  IStatus istatus;
  set<Cell*> cells;
  set<Cell*> paused_cells;
  Vector refpoint;

  double current_precision;

  Timer timer;

};

} /* namespace ibex */

#endif /* PLUGINS_OPTIM_MOP_SRC_STRATEGY_IBEX_OPTIMIZERMOPSERVER_H_ */
