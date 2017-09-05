//============================================================================
//                                  I B E X
// File        : ibex_DefaultSolver.cpp
// Author      : Matias Campusano, Damir Aliquintui, Ignacio Araya
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Sep 05, 2017
// Last Update : Sep 05, 2017
//============================================================================

#include "ibex_CellFeasibleDiving.h"

using namespace std;

namespace ibex {

	int CellBS::nb_cells=0;

	/*================================== inline implementations ========================================*/


	    // functions about CellFeasibleDiving
	CellFeasibleDiving::CellFeasibleDiving(const ExtendedSystem& sys) :
	  		bufferset(*new CellSet<minLB>),
	      loup_lb(POS_INFINITY),
	  		sys(sys) {
	}

	CellFeasibleDiving::~CellFeasibleDiving() { }

	// TODO: verificar si cellset tiene add_backtrackable
	void CellFeasibleDiving::add_backtrackable(Cell& root) {
	      root.add<CellBS>();
	}

	std::ostream& CellFeasibleDiving::print(std::ostream& os) const
	   {	os << "==============================================================================\n";
	     os << " first cell " << " size " << size() << " top " << minimum() << std::endl;
	       return  os << std::endl;
	   }

	    // functions about CellSet
	  void CellFeasibleDiving::push(Cell* cell) {
	      // TODO: imprimir erro cuando ningun nodo sea nulo
	      if(cell->box[cell->box.size()-1].lb() < loup_lb) {
	          cell->get<CellBS>().lb=cell->box[cell->box.size()-1].lb();
	         // std::cout << cell->get<CellBS>().lb << std::endl;
	          if(cl == NULL) {
	            cl = cell;
	          } else {
	            cr = cell;
	          }
	      }
	  }

	  Cell* CellFeasibleDiving::top() const {
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

	  Cell* CellFeasibleDiving::pop() {
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

	  bool CellFeasibleDiving::empty() const {
	      if(cr == NULL && cl == NULL) {
	          return bufferset.empty();
	      } else {
	          return false;
	      }
	  }

	  unsigned int CellFeasibleDiving::size() const {
	    unsigned int size = bufferset.size();
	    if(cr != NULL) {
	      size++;
	    }
	    if(cl != NULL) {
	      size++;
	    }
	    return size;
	  }

	  void CellFeasibleDiving::flush() {
	      if(cr != NULL)
	          delete cr;
	      if(cl != NULL)
	          delete cl;
	      cr = NULL;
	      cl = NULL;
	      bufferset.flush();
	  }


	  // functions about CellBufferOptim
	  void CellFeasibleDiving::contract(double new_loup) {
	      loup_lb = new_loup;
	      // TODO: Ver cuando se hace el contract guardar el new_loup para no agregar otros
	      if(new_loup < minimum()) {
	        CellFeasibleDiving::flush();
	      }
	  }

	  double CellFeasibleDiving::minimum() const  {
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



} // end namespace ibex
