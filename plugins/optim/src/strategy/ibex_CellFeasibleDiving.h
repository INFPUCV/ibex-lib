#ifndef __IBEX_CELL_FEASIBLE_DIVING_H__
#define __IBEX_CELL_FEASIBLE_DIVING_H__

#include "ibex_CellBufferOptim.h"
#include "ibex_CellSet.h"
#include "ibex_ExtendedSystem.h"
#include "ibex_CellCostFunc.h"

namespace ibex {

    int CellBS::nb_cells=0;
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

        double loup_lb;
  };

  /*================================== inline implementations ========================================*/


    // functions about CellFeasibleDiving
  inline CellFeasibleDiving::CellFeasibleDiving(const ExtendedSystem& sys, int crit2_pr, CellCostFunc::criterion crit2) :
  		bufferset(*new CellSet<minLB>),
      loup_lb(POS_INFINITY),
  		sys(sys) {
  }

  inline CellFeasibleDiving::~CellFeasibleDiving() { }

// TODO: verificar si cellset tiene add_backtrackable
  inline void CellFeasibleDiving::add_backtrackable(Cell& root) {
      root.add<CellBS>();
  }

   inline std::ostream& CellFeasibleDiving::print(std::ostream& os) const
   {	os << "==============================================================================\n";
     os << " first cell " << " size " << size() << " top " << minimum() << std::endl;
       return  os << std::endl;
   }

    // functions about CellSet
  inline void CellFeasibleDiving::push(Cell* cell) {
      // TODO: imprimir erro cuando ningun nodo sea nulo
      if(cell->box[cell->box.size()-1].lb() < loup_lb) {
          cell->get<CellBS>().lb=cell->box[cell->box.size()-1].lb();
          std::cout << cell->get<CellBS>().lb << std::endl;
          if(cl == NULL) {
            cl = cell;
          } else {
            cr = cell;
          }
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
          bufferset.push(cl);
          cl = NULL;
      } else {
          c = cl;
          cl = NULL;
          bufferset.push(cr);
          cr = NULL;
      }
      CellFeasibleDiving::contract(loup_lb);
      return c;
  }

  inline bool CellFeasibleDiving::empty() const {
      if(cr == NULL && cl == NULL) {
          return bufferset.empty();
      } else {
          return false;
      }
  }

  inline unsigned int CellFeasibleDiving::size() const {
    unsigned int size = bufferset.size();
    if(cr != NULL) {
      size++;
    }
    if(cl != NULL) {
      size++;
    }
    return size;
  }

  inline void CellFeasibleDiving::flush() {
      if(cr != NULL)
          delete cr;
      if(cl != NULL)
          delete cl;
      cr = NULL;
      cl = NULL;
      bufferset.flush();
  }


  // functions about CellBufferOptim
  inline void CellFeasibleDiving::contract(double new_loup) {
      loup_lb = new_loup;
      // TODO: Ver cuando se hace el contract guardar el new_loup para no agregar otros
      if(new_loup < minimum()) {
        CellFeasibleDiving::flush();
      }
  }

  inline double CellFeasibleDiving::minimum() const  {
      Cell *cellBuffer = bufferset.top(), *cellTop = CellFeasibleDiving::top();
      if(cellTop == NULL) {
        return POS_INFINITY;
      }
      if(cellBuffer != NULL && cellBuffer->box[cellBuffer->box.size()-1].lb() < cellTop->box[cellTop->box.size()-1].lb()) {
          return cellBuffer->box[cellBuffer->box.size()-1].lb();
      } else {
          return cellTop->box[cellTop->box.size()-1].lb();
      }
  }


} // namespace ibex

#endif // __IBEX_CELL_FEASIBLE_DIVING_H__
