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

  void update_refpoint(Vector& refpoint, double eps);

	void plot();

  /*Return a list of pairs. The first element indicates if the second has
	to be inserted or removed from the upper_envelope */
	list  < pair < bool, Vector> > changes_upper_envelope(int nb_changes=-1){
		return ndsH.get_and_clear_changes(nb_changes);
	}

	/* Return a list of pairs. The first element indicates if the second has
	to be inserted or removed from the lower_envelope */
	list  < pair < bool, Vector> > changes_lower_envelope(int nb_changes=-1);

	void lower_envelope_tostring(list<vector<double> > &lowerList, double eps){
		//the envelope is generated
		NDS_seg LBseg;
		for(auto cc:cells)	LBseg.add_lb(*cc);
		for(auto cc:paused_cells) LBseg.add_lb(*cc);
		LBseg.to_string(lowerList, eps, false);
	}

	void upper_envelope_tostring(list<vector<double> > &upperList, double eps){
		ndsH.to_string(upperList, eps, true);
	}

 	void write_envelope(string output_file);

  /****************************************/

  IntervalVector initbox;

  IStatus istatus;
  set<Cell*> cells;
  set<Cell*> paused_cells;
  Vector refpoint;

  double current_precision;

	map< Vector, NDS_data, sorty2 > LBseg;

  Timer timer;
	list  < pair < bool, Vector> > changes_lower;

};

} /* namespace ibex */

#endif /* PLUGINS_OPTIM_MOP_SRC_STRATEGY_IBEX_OPTIMIZERMOPSERVER_H_ */
