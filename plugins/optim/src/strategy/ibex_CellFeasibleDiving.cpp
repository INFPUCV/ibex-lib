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
	template<class T>
	CellFeasibleDiving<T>::CellFeasibleDiving(CellBufferOptim& cset) :
			bufferset(cset)/*, sys(sys)*/, cl(NULL), cr(NULL) {
	}

	template<class T>
	CellFeasibleDiving<T>::~CellFeasibleDiving() { }

	// TODO: verificar si cellset tiene add_backtrackable
	template<class T>
	void CellFeasibleDiving<T>::add_backtrackable(Cell& root) {
	      root.add<CellBS>();
	}

	template<class T>
	std::ostream& CellFeasibleDiving<T>::print(std::ostream& os) const
	{	os << "==============================================================================\n";
		 os << " first cell " << " size " << size() << " top " << minimum() << std::endl;
		   return  os << std::endl;
	}

		// functions about CellSet
	template<class T>
	void CellFeasibleDiving<T>::push(Cell* cell) {

		if(cl == NULL)
		 cl = cell;
		else if(cr == NULL)
		 cr = cell;
		else
		 ibex_error("CellFeasibleDiving: triple push error");
	}
	  template<class T>
	  Cell* CellFeasibleDiving<T>::top() const {

	      if(cr == NULL && cl == NULL) {
	        return bufferset.top();
	      } else if(cl == NULL) {
	        return cr;
	      } else if(cr == NULL) {
	        return cl;
	      } else if(T()(cr, cl)) {
	        return cr;
	      } else {
	        return cl;
	      }
	  }

	  template<class T>
	  Cell* CellFeasibleDiving<T>::pop() {
	      Cell *c;
	      if(cr == NULL && cl == NULL) {
	          c = bufferset.pop();
	      } else if(cl == NULL) {
	          c = cr;
	          cr = NULL;
	      } else if(cr == NULL) {
	          c = cl;
	          cl = NULL;
	      } else if(T()(cr, cl)) {
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

	      return c;
	  }

	  template<class T>
	  bool CellFeasibleDiving<T>::empty() const {
	      if(cr == NULL && cl == NULL) {
	          return bufferset.empty();
	      } else {
	          return false;
	      }
	  }

	  template<class T>
	  unsigned int CellFeasibleDiving<T>::size() const {
	    unsigned int size = bufferset.size();
	    if(cr != NULL) {
	      size++;
	    }
	    if(cl != NULL) {
	      size++;
	    }
	    return size;
	  }

	  template<class T>
	  void CellFeasibleDiving<T>::flush() {
	      if(cr != NULL)
	          delete cr;
	      if(cl != NULL)
	          delete cl;
	      cr = NULL;
	      cl = NULL;
	      bufferset.flush();
	  }


	  // functions about CellBufferOptim
	  template<class T>
	  void CellFeasibleDiving<T>::contract(double new_loup) {
		  bufferset.contract(new_loup);
		  if(cr && cr->box[cr->box.size()-1].lb() > new_loup){
			  delete cr; cr=NULL;
		  }

		  if(cl && cl->box[cl->box.size()-1].lb() > new_loup){
			  delete cl; cl=NULL;
		  }
	  }

	  template<class T>
	  double CellFeasibleDiving<T>::minimum() const  {
		  double min = bufferset.minimum();
		  if(cr && cr->box[cr->box.size()-1].lb() < min)
			  min=cr->box[cr->box.size()-1].lb();

		  if(cl && cl->box[cl->box.size()-1].lb() < min)
			  min=cl->box[cl->box.size()-1].lb();

		  return min;
	  }

	  template class CellFeasibleDiving<minLB>;
	  template class CellFeasibleDiving<maxsize>;

} // end namespace ibex
