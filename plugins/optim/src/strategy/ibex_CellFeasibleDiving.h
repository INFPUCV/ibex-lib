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
    CellFeasibleDiving(const ExtendedSystem& sys, int crit2_pr=50,
  			CellCostFunc::criterion crit2=CellCostFunc::UB);

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

    CellBuffer& bufferset;

    protected:

      	/**
      	 * The system
      	 */
      	const ExtendedSystem& sys;

        Cell *cl=NULL,*cr=NULL;


  };

  /*================================== inline implementations ========================================*/

    // functions about CellFeasibleDiving
  inline CellFeasibleDiving::CellFeasibleDiving(const ExtendedSystem& sys, int crit2_pr, CellCostFunc::criterion crit2) :
  		bufferset(*new CellSet<minLB>),
  		sys(sys) {
  }

  inline CellFeasibleDiving::~CellFeasibleDiving() { }

// TODO: verificar si cellset tiene add_backtrackable
  inline void CellFeasibleDiving::add_backtrackable(Cell& root) {
      // root.add<CellBS>();
  }

   inline std::ostream& CellFeasibleDiving::print(std::ostream& os) const
   {	os << "==============================================================================\n";
     os << " first cell " << " size " <<"test" << " top " << "test" << std::endl;
       os << " second cell " << " size " << "test" << " top " << "test" ;
       return  os << std::endl;
   }

    // functions about CellSet
  inline void CellFeasibleDiving::push(Cell* cell) {
      // TODO: imprimir erro cuando ningun nodo sea nulo
      // cell->get<CellBS>().lb=cell->box[cell->box.size()-1].lb();
      if(cl == NULL) {
        cl = cell;
      } else {
        cr = cell;
      }
  }

  inline Cell* CellFeasibleDiving::top() const {
      if(cr == NULL && cl == NULL) {
        return bufferset.top();
      } else if(cr != NULL && cl == NULL) {
        return cr;
      } else if(cr == NULL && cl != NULL) {
        return cl;
      } else if(cr->box[cr->box.size()-1].lb() < cl->box[cl->box.size()-1].lb()) {
        return cr;
      } else {
        return cl;
      }
  }

  inline Cell* CellFeasibleDiving::pop() {
      Cell *c;
      if(cr == NULL && cl == NULL) {
          c = bufferset.pop();
      } else if(cr != NULL && cl == NULL) {
          c = cr;
          cr = NULL;
      } else if(cr == NULL && cl != NULL) {
          c = cl;
          cl = NULL;
      } else if(cr->box[cr->box.size()-1].lb() < cl->box[cl->box.size()-1].lb()) {
          c = cr;
          cr = NULL;
      } else {
          c = cl;
          cl = NULL;
      }
      return c;
  }

  inline bool CellFeasibleDiving::empty() const {
      if(cr == NULL && cl == NULL) {
          return bufferset.empty();
      } else {
          return false;
      }
  }

  inline unsigned int CellFeasibleDiving::size() const { bufferset.size(); }

  inline void CellFeasibleDiving::flush() { bufferset.flush(); }


  // functions about CellBufferOptim
  inline void CellFeasibleDiving::contract(double new_loup) {
    // TODO: contraer cr y cl
  	// bufferset.contract(new_loup);
  }

  inline double CellFeasibleDiving::minimum() const     { return 0; }


} // namespace ibex

#endif // __IBEX_CELL_FEASIBLE_DIVING_H__
