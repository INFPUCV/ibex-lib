#ifndef __IBEX_CELL_FEASIBLE_DIVING_H__
#define __IBEX_CELL_FEASIBLE_DIVING_H__

#include "ibex_CellBufferOptim.h"
#include "ibex_CellSet.h"
#include "ibex_ExtendedSystem.h"
#include "ibex_CellCostFunc.h"
#include "ibex_CellNSSet.h"

#include <map>

namespace ibex {


/**
 * Criteria for bi-objective problems
 */
struct max_distance {
	/**
	 * \brief distance from the box to the non dominated set
	 */
	double distance(const IntervalVector& b){
	   int n=b.size();
     double min_dist = NEG_INFINITY;
     map< pair <double, double>, Vector >::const_iterator it_lb=UB->lower_bound(make_pair(b[n-2].lb(),POS_INFINITY));

     for (;it_lb!=UB->end(); it_lb++){
  	   pair <double, double> z = it_lb->first;
  	   double dist = std::max(z.first -  b[n-2].lb(), z.second - b[n-1].lb());
  	   if(dist < min_dist) min_dist=dist;
  	   if(z.second <= b[n-1].lb()) break;
     }

     return min_dist;
	}

	bool operator() (const Cell* c1, const Cell* c2){
       return (distance(c1->box) > distance(c2->box));
    }

	static map< pair <double, double>, Vector >* UB;
};

  //  TODO: verificar que sea CellSet en vez de CellBuffer (por que?, ya no me acuerdo)
  // TODO: Agregar descripcion (comentarios) a las funciones
  template<class T>
  class CellFeasibleDiving : public CellBufferOptim {
  public:
    CellFeasibleDiving(CellBufferOptim& cset /*const ExtendedSystem& sys*/);

    ~CellFeasibleDiving();

    virtual void add_backtrackable(Cell& root);

    virtual void flush();

    virtual unsigned int size() const;

  	virtual bool empty() const;

  	virtual void push(Cell* cell);

  	virtual Cell* pop();

  	virtual Cell* top() const;

  	virtual std::ostream& print(std::ostream& os) const;

  	virtual double minimum() const;

  	virtual void contract(double loup);

    CellBufferOptim& bufferset;

    protected:

      	/**
      	 * The system
      	 */
      	//const ExtendedSystem& sys;

        Cell *cl;
		Cell *cr;
  };

} // namespace ibex

#endif // __IBEX_CELL_FEASIBLE_DIVING_H__
