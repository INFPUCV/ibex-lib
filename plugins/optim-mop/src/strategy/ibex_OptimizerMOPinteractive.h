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
  IStatus run(int maxiter=5, double eps=1e-1);

  void save_state_in_file(string filename){
     //write object into the file
     // Object to write in file
     ofstream file_obj;

     // Opening file in append mode
     file_obj.open(filename, ios::out);
     file_obj.write((char*) &BxpMOPData::y1_init, sizeof(BxpMOPData::y1_init));
     file_obj.write((char*) &BxpMOPData::y2_init, sizeof(BxpMOPData::y2_init));

     int len = cells.size();
     file_obj.write((char*) &len, sizeof(int));
     for (auto c:cells){
       file_obj.write((char*) &c->box[0], sizeof(c->box[0])*c->box.size());
       file_obj.write((char*) &cdata->a, sizeof(cdata->a));
       file_obj.write((char*) &cdata->w_lb, sizeof(cdata->w_lb));
       file_obj.write((char*) &cdata->ub_distance, sizeof(cdata->ub_distance));
     }

     len = paused_cells.size();
     file_obj.write((char*) &len, sizeof(int));
     for (auto c:paused_cells){
       file_obj.write((char*) &c->box[0], sizeof(c->box[0])*c->box.size());
       file_obj.write((char*) &cdata->a, sizeof(cdata->a));
       file_obj.write((char*) &cdata->w_lb, sizeof(cdata->w_lb));
       file_obj.write((char*) &cdata->ub_distance, sizeof(cdata->ub_distance));
     }

     len = ndsH.size();

     file_obj.write((char*) &len, sizeof(int));
     for (auto elem:ndsH.NDS2){
       file_obj.write((char*) &elem.first[0], sizeof(elem.first[0])*elem.first.size());
       file_obj.write((char*) &elem.second.n, sizeof(elem.second.n));
       if(elem.second.n >= 1) file_obj.write((char*) &elem.second.x1[0], sizeof(elem.second.x1[0])*elem.second.x1.size());
       if(elem.second.n == 2) file_obj.write((char*) &elem.second.x2[0], sizeof(elem.second.x2[0])*elem.second.x2.size());
     }

     for (auto p : ndsH.NDS2) cout << p.first << endl;

     file_obj.close();
  }

  void load_state_from_file(string filename, const IntervalVector& init_box){
    // Object to read from file
    ifstream file_obj;

    // Opening file in input mode
    file_obj.open(filename, ios::in);
    file_obj.read((char*)&BxpMOPData::y1_init, sizeof(BxpMOPData::y1_init));
    file_obj.read((char*)&BxpMOPData::y2_init, sizeof(BxpMOPData::y2_init));

    int len;
    file_obj.read((char*) &len, sizeof(int));
    for (int i=0; i<len; i++){
      Cell* c=new Cell(init_box);
      file_obj.read((char*)&c->box[0], sizeof(c->box[0])*c->box.size());
      c->prop.add(new BxpMOPData());
    	buffer.add_property(c->box, c->prop);
    	bsc.add_property(c->box, c->prop);
    	ctc.add_property(c->box, c->prop);

      file_obj.read((char*) &cdata->a, sizeof(cdata->a));
      file_obj.read((char*) &cdata->w_lb, sizeof(cdata->w_lb));
      file_obj.read((char*) &cdata->ub_distance, sizeof(cdata->ub_distance));
      cells.insert(c);
    }

    file_obj.read((char*) &len, sizeof(int));
    for (int i=0; i<len; i++){
      Cell* c=new Cell(init_box);
      file_obj.read((char*)&c->box[0], sizeof(c->box[0])*c->box.size());
      c->prop.add(new BxpMOPData());
      buffer.add_property(c->box, c->prop);
    	bsc.add_property(c->box, c->prop);
    	ctc.add_property(c->box, c->prop);

      file_obj.read((char*) c->prop[BxpMOPData::id], sizeof(BxpMOPData));
      cells.insert(c);
    }

    file_obj.read((char*) &len, sizeof(int));
    ndsH.NDS2.clear ();

    for (int i=0; i<len; i++){
      Vector y(2);
      int nn;

      file_obj.read((char*) &y[0], sizeof(y[0])*2);
      file_obj.read((char*) &nn, sizeof(int));
      Vector x1(1);
      Vector x2(1); //or n+2?

      if(nn>=1) {x1.resize(n); file_obj.read((char*) &x1[0], sizeof(x1[0])*(n));}
      if(nn==2) {x2.resize(n); file_obj.read((char*) &x2[0], sizeof(x2[0])*(n));}

      ndsH.NDS2[y]=NDS_data(x1,x2);
      ndsH.NDS2[y].n=nn;
    }
    cout.precision(17);
    for (auto p : ndsH.NDS2) cout << p.first << endl;


    file_obj.close();
  }

  void update_refpoint(Vector& refpoint);

 	void write_envelope(string output_file);


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
