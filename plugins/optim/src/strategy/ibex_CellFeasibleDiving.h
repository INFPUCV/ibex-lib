#ifndef __IBEX_CELL_FEASIBLE_DIVING_H__
#define __IBEX_CELL_FEASIBLE_DIVING_H__

#include "ibex_CellBufferOptim.h"
#include "ibex_CellSet.h"
#include "ibex_ExtendedSystem.h"
#include "ibex_CellCostFunc.h"

namespace ibex {


  //  TODO: verificar que sea CellSet en vez de CellBuffer
  class CellFeasibleDiving: public CellBufferOptim {
  public:
    CellFeasibleDiving(const ExtendedSystem& sys);

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

    CellBuffer& bufferset;

    protected:

      	/**
      	 * The system
      	 */
      	const ExtendedSystem& sys;

        Cell *cl=NULL,*cr=NULL;

        double loup_lb;
  };


} // namespace ibex

#endif // __IBEX_CELL_FEASIBLE_DIVING_H__
