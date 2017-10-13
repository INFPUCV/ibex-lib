#ifndef __IBEX_CELL_FEASIBLE_DIVING_H__
#define __IBEX_CELL_FEASIBLE_DIVING_H__

#include "ibex_CellBufferOptim.h"
#include "ibex_CellSet.h"
#include "ibex_ExtendedSystem.h"
#include "ibex_CellCostFunc.h"
#include "ibex_CellNSSet.h"

namespace ibex {


  //  TODO: verificar que sea CellSet en vez de CellBuffer
  // TODO: Agregar descripcion (comentarios) a las funciones
  template<class T>
  class CellFeasibleDiving: public CellBufferOptim {
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
