#ifndef __IBEX_CELL_FEASIBLE_DIVING_H__
#define __IBEX_CELL_FEASIBLE_DIVING_H__

#include "ibex_CellBufferOptim.h"
#include "ibex_CellSet.h"
#include "ibex_ExtendedSystem.h"

namespace ibex {

  class CellFeasibleDiving: CellSet<Cell>, CellBufferOptim {
  public:
    CellFeasibleDiving();

    ~CellFeasibleDiving();

    virtual void add_backtrackable(Cell& root);

    void flush();

  	unsigned int size() const;

  	bool empty() const;

  	void push(Cell* cell);

  	Cell* pop();

  	Cell* top() const;

  	std::ostream& print(std::ostream& os) const;

  	virtual double minimum() const;

  	virtual void contract(double loup);
  }

} // namespace ibex

#endif // __IBEX_CELL_FEASIBLE_DIVING_H__
