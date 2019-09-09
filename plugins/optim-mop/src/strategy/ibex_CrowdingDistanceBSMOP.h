/*
 * ibex_CrowdingDistanceBSMOP.h
 *
 *  Created on: 20 oct. 2017
 *      Author: matias y pablo
 */

#ifndef OPTIM_SRC_STRATEGY_IBEX_CROWDINGDISTANCEBSMOP_H_
#define OPTIM_SRC_STRATEGY_IBEX_CROWDINGDISTANCEBSMOP_H_



#include "ibex_CellMOP.h"
#include "ibex_CellSet.h"
#include "ibex_NDS.h"
#include "ibex_OptimizerMOP.h"
#include <queue>
#include <map>


using namespace std;

namespace ibex {

/**
 * Criteria for bi-objective problems
 */


//Struct that contains Cell and crowding distance
	typedef struct CDBox{
		Cell* C;
		double crowding_distance;
        int size=0;
	} CDBox;



    //Ordena por objetivo y1, y2
    struct crowding_distanceBeam {
        bool operator() (const Cell* c1, const Cell* c2){
            int n = c1->box.size();
            if(c1->box[n-2].lb() <= c2->box[n-2].lb() && c1->box[n-1].lb() <= c2->box[n-1].lb()) return true;
            else return false;
        }
    };

    struct sortByCrowdingDistance{
        bool operator()(const CDBox* c1, const CDBox* c2){
            if(c1->crowding_distance <= c2->crowding_distance) return true;
            else return false;
        }

    };


struct max_distanceCrowding {


	bool operator() (const Cell* c1, const Cell* c2){
	   int n = c1->box.size();
	   if(c1->get<CellMOP>().ub_distance != c2->get<CellMOP>().ub_distance)
		   return (c1->get<CellMOP>().ub_distance < c2->get<CellMOP>().ub_distance);
	   else if(c1->box[n-2].lb() >= c2->box[n-2].lb() && c1->box[n-1].lb() >= c2->box[n-1].lb()) return true;
	   else return false;
	}

};

struct objective_sort {


	bool operator() (const Cell* c1, const Cell* c2){

	    int n = c1->box.size();
	    if (c1->box[n-2].lb() < c2->box[n-2].lb()) return true;
	    return false;

	}

};

// struct min_distanceBeam {


// 	bool operator() (const Cell* c1, const Cell* c2){
// 	   int n = c1->box.size();
// 	   if(c1->get<CellMOP>().ub_distance != c2->get<CellMOP>().ub_distance)
// 		   return (c1->get<CellMOP>().ub_distance > c2->get<CellMOP>().ub_distance);
// 	   else if(c1->box[n-2].lb() < c2->box[n-2].lb() && c1->box[n-1].lb() < c2->box[n-1].lb()) return true;
// 	   else return false;
// 	}

// };

/** \ingroup strategy
 *
 * \brief Buffer which selects next the box maximizing the distance to the non dominated set.
 */
class CrowdingDistanceBSMOP : public CellBufferOptim {
 public:

	//RECORDAR CAMBIAR NEXTBUFFERSIZE POR CURRENTBUFFERMAXSIZE (EN EL MAIN), TAMBIEN CON EL AUXILIAR
	CrowdingDistanceBSMOP(int currentBufferMaxSize, double bs_tolerance, bool crowding) : CellBufferOptim(),
	depth(0), currentBufferMaxSize(currentBufferMaxSize), bs_tolerance(bs_tolerance), crowding(crowding){

	}

   int currentBufferMaxSize=8;
   int currentBufferSizeAux=0;
   double bs_tolerance=0.5;
   int T=1;
   bool crowding=false;
   int bs_level=0;
   int max_level=0;

   //static int nn;
   int iterBS=0;
    
        //static std::multiset<Cell*,max_distanceBeam> getCrowdingDistance(
   void crowdingDistance(
   	        std::multiset<Cell*, max_distanceCrowding>& nextBuffer,
   	        std::priority_queue<Cell*, std::vector<Cell*>, max_distanceCrowding >& currentBuffer,
   	        std::priority_queue<Cell*, std::vector<Cell*>, max_distanceCrowding >& globalBuffer,
   	        int currentBufferMaxSize
   	        );

    static bool dominates(Cell* a, Cell* c);



    void nonDominatedSort(
            std::multiset<Cell*, max_distanceCrowding>& nextBuffer,
            std::priority_queue<Cell*, std::vector<Cell*>, max_distanceCrowding >& currentBuffer,
            std::priority_queue<Cell*, std::vector<Cell*>, max_distanceCrowding >& globalBuffer,
            int currentBufferMaxSize
            );

    void extractNonDominated(
        std::multiset<Cell*, max_distanceCrowding>& nextBuffer, std::multiset<Cell*, max_distanceCrowding>& nonDominated);





   void set(NDS_seg& nds, double& pruebaprom, int&depthMayor) {
		 this->nds=&nds;
		 this->pruebaprom=&pruebaprom;
		 this->depthMayor=&depthMayor;
	 }

   virtual void add_backtrackable(Cell& root){
     root.add<CellMOP>();
   }

  /** Flush the buffer.
   * All the remaining cells will be *deleted* */
  void flush();

  /** Return the size of the buffer. */
  unsigned int size() const;

  /** Return true if the buffer is empty. */
  bool empty() const;

  void bsPerformance(int iter);

  /** push a new cell on the stack. */
  void push(Cell* cell);

  /** Pop a cell from the stack and return it.*/
  Cell* pop();

  /** Return the next box (but does not pop it).*/
  Cell* top() const;

  /**
	* \brief Return the minimum value of the heap
	*
	*/
  virtual double minimum() const {
	  cout << "BeamSearchBufferMOP::minimum is not implemented" << endl;
	  exit(0);
	  return 0.0;
  }

	/**
	 * \brief Contract the buffer using the UB
	 */
	virtual void contract(double loup){
		  cout << "BeamSearchBufferMOP::contract is not implemented" << endl;
		  exit(0);
	}


	/**
	 * A heap data structure for keeping the cells sorted by distance
	 */

	mutable std::priority_queue<Cell*, std::vector<Cell*>, max_distanceCrowding > globalBuffer;
    mutable std::priority_queue<Cell*, std::vector<Cell*>, max_distanceCrowding > currentBuffer;
	//mutable std::priority_queue<Cell*, std::vector<Cell*>, max_distanceBeam > nextBuffer;
    mutable std::multiset <Cell*, max_distanceCrowding> nextBuffer;

  NDS_seg* nds=NULL;
  double* pruebaprom=NULL;

  bool bs_state=false;
  bool bs_performance=false;
  int* depthMayor=0;
  int depth=0, depthTotal=0;
  double depthPromedio=0;

  private:
	int cont = 0, cantBeam = 0;
	double hv=0,hv0=0,hvE=0,initial_reduction=0.0,mejora=0;
	double mejor_mejora=0.0;
};






} // end namespace ibex
#endif //  /* OPTIM_SRC_STRATEGY_IBEX_BeamSearchBufferMOP_H_ */
