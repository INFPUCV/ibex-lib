//============================================================================
//                                  I B E X
// File        : ibex_CellSet.h
// Author      : Gilles Chabert, Jordan Ninin
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Sep 12, 2014
//============================================================================

#ifndef __IBEX_CELL_SET_H__
#define __IBEX_CELL_SET_H__

#include "ibex_Random.h"
#include "ibex_CellBuffer.h"
#include <set>
#include <map>
#include <queue>
// #include "ibex_Set.h"
// #include "../strategy/ibex_Cell.h"

using namespace std;

namespace ibex {

	class CellBS : public Backtrackable {
	public:
		/**
		 * \brief Constructor for the root node (followed by a call to init_root).
		 */
		CellBS() : depth(0), id(0), a(0.0), w_lb(0.0), ub_distance(POS_INFINITY) {}

		/**
		 * \brief Copy constructor
		 */

		CellBS(const CellBS& c) : depth(c.depth+1), id(nb_cells++),
				a(c.a), w_lb(c.w_lb), ub_distance(c.ub_distance) { }

		/**
		 * \brief Duplicate the structure into the left/right nodes
		 */
		std::pair<Backtrackable*,Backtrackable*> down(){
			CellBS* c1=new CellBS(*this);
			CellBS* c2=new CellBS(*this);

			return std::pair<Backtrackable*,Backtrackable*>(c1,c2);
		}

		static int nb_cells;
		static Interval y1_init;
		static Interval y2_init;

	    /**unique identifier for comparisons*/
	    int id;

		/** depth of the node **/
		int depth;

		/** MOP: after filtering we know that z1+a*z2 > w_lb and we can
		 * use this information for filtering**/
		double a;
		double w_lb;


		/** MOP: distance of the box to the current non dominated set (UB) */
		double ub_distance;

	};


	/**
	 * Criteria for bi-objective problems
	 */
	struct max_distance {


		bool operator() (const Cell* c1, const Cell* c2){
		   int n = c1->box.size();
		   if(c1->get<CellBS>().ub_distance != c2->get<CellBS>().ub_distance)
			   return (c1->get<CellBS>().ub_distance < c2->get<CellBS>().ub_distance);
		   else if(c1->box[n-2].lb() >= c2->box[n-2].lb() && c1->box[n-1].lb() >= c2->box[n-1].lb()) return true;
		   else return false;
	    }

		static map< pair <double, double>, IntervalVector >* UB;
	};

	template<class T>
	class CellSet : public CellBufferOptim {
	public:

	  CellSet();

	  void flush();

	    virtual void add_backtrackable(Cell& root){
	    	root.add<CellBS>();
	    }

	  /** Return the size of the buffer. */
	  unsigned int size() const;

	  /** Return true if the buffer is empty. */
	  bool empty() const;

	  /** push a new cell on the stack. */
	  void push(Cell* cell);

	  /** Pop a cell from the stack and return it.*/
	  Cell* pop();

	  /** Return the next box (but does not pop it).*/
	  Cell* top() const;


	  virtual double minimum() const;

	  virtual void contract(double loup);

	private:
		/* Set of Cells */
		typename std::multiset<Cell*, T> cset;

	};

	/**
	 * This criterion corresponds to the SR1 criterion of the paper of Fernandez&Toth (2007) for bi-objective problems
	 * Authors says that in this way the curve is generated from left-top to right-bottom
	 */
	struct minLB {
	  bool operator() (const Cell* c1, const Cell* c2) const
	  {
		  int n = c1->box.size();

		  if(c1->box[n-1].lb() != c2->box[n-1].lb()) return (c1->box[n-1].lb() < c2->box[n-1].lb());
		  return (c1->get<CellBS>().depth < c2->get<CellBS>().depth);

	  }
	};

	/**
	 * OC1 in https://tel.archives-ouvertes.fr/tel-01146856/document
	 * min y1.lb
	 */
	struct OC1 {
	  bool operator() (const Cell* c1, const Cell* c2) const
	  {
		  int n = c1->box.size();

		  if(c1->box[n-1].lb() != c2->box[n-1].lb()) return (c1->box[n-1].lb() < c2->box[n-1].lb());
		  if(c1->box[n-2].lb() != c2->box[n-2].lb()) return (c1->box[n-2].lb() < c2->box[n-2].lb());
		  return (c1->get<CellBS>().depth < c2->get<CellBS>().depth);
		  //return (c1->get<CellBS>().id > c2->get<CellBS>().id);
	  }
	};

	/**
	 * OC2 in https://tel.archives-ouvertes.fr/tel-01146856/document
	 * min y2.lb
	 */
	struct OC2 {
	  bool operator() (const Cell* c1, const Cell* c2) const
	  {
		  int n = c1->box.size();

		  if(c1->box[n-2].lb() != c2->box[n-2].lb()) return (c1->box[n-2].lb() < c2->box[n-2].lb());
		  if(c1->box[n-1].lb() != c2->box[n-1].lb()) return (c1->box[n-1].lb() < c2->box[n-1].lb());
		  return (c1->get<CellBS>().depth < c2->get<CellBS>().depth);
		  //return (c1->get<CellBS>().id > c2->get<CellBS>().id);
	  }
	};

	/**
	 * OC3 in https://tel.archives-ouvertes.fr/tel-01146856/document
	 * increasing value of y1.lb + y2.lb
	 */
	struct OC3 {
	  bool operator() (const Cell* c1, const Cell* c2) const
	  {
		  int n = c1->box.size();

		  if(c1->box[n-1].lb() + c1->box[n-2].lb() != c2->box[n-1].lb() + c2->box[n-2].lb())
			  return (c1->box[n-1].lb() + c1->box[n-2].lb() < c2->box[n-1].lb() + c2->box[n-2].lb());
		  return (c1->get<CellBS>().depth < c2->get<CellBS>().depth);
	  }
	};

	/**
	 * OC4 in https://tel.archives-ouvertes.fr/tel-01146856/document
	 * decreasing value of hypervolume of the point y
	 */
	struct OC4 {
	  bool operator() (const Cell* c1, const Cell* c2) const
	  {
		  int n = c1->box.size();

		  double hyper1=(CellBS::y1_init.ub()-c1->box[n-1].lb())*(CellBS::y2_init.ub()-c1->box[n-2].lb());
		  double hyper2=(CellBS::y1_init.ub()-c1->box[n-1].lb())*(CellBS::y2_init.ub()-c1->box[n-2].lb());

		  if(hyper1 != hyper2) return (hyper1 < hyper2);
		  return (c1->get<CellBS>().depth < c2->get<CellBS>().depth);
	  }
	};

	/**
	 * Criteria for bi-objective problems used in the paper by Martin et al. (2016)
	 */
	struct weighted_sum {
	  bool operator() (const Cell* c1, const Cell* c2) const
	  {
		  int n = c1->box.size();
		  double c1_ev= (c1->box[n-2].lb()-CellBS::y1_init.lb())/CellBS::y1_init.diam() +
				  (c1->box[n-1].lb()-CellBS::y2_init.lb())/CellBS::y2_init.diam();

		  double c2_ev= (c2->box[n-2].lb()-CellBS::y1_init.lb())/CellBS::y1_init.diam() +
				  (c2->box[n-1].lb()-CellBS::y2_init.lb())/CellBS::y2_init.diam();

		  return c1_ev < c2_ev;
	  }
	};



	template<class T>
	CellSet<T>::CellSet() {
	}

	template<class T>
	void CellSet<T>::flush() {
		while (!cset.empty()) {
			delete *cset.begin();
			cset.erase(cset.begin());
		}
	}

	template<class T>
	unsigned int CellSet<T>::size() const {
		return cset.size();
	}

	template<class T>
	bool CellSet<T>::empty() const {
		return cset.empty();
	}

	template<class T>
	void CellSet<T>::push(Cell* cell) {
		if (capacity>0 && size() == capacity) throw CellBufferOverflow();
		cset.insert(cell);
	}

	template<class T>
	Cell* CellSet<T>::pop() {
		Cell* c = *cset.begin();
		cset.erase(cset.begin());
		return c;
	}

	template<class T>
	Cell* CellSet<T>::top() const{
		return *cset.begin();
	}

	template<class T>
	double CellSet<T>::minimum() const{
	      if(size()==0)
	        return POS_INFINITY;
	      else
	    	return top()->box[top()->box.size()-1].lb();
	}

	template<class T>
	void CellSet<T>::contract(double loup){
		  typename std::set<Cell*, T>::iterator it= cset.begin();

		  while(it!= cset.end()){
			  if( (*it)->box[(*it)->box.size()-1].lb() > loup ){
				  typename std::set<Cell*, T>::iterator it2=it; it2++;
				  cset.erase(it);
				  it=it2;
			  }else it++;
		  }
	}


	template class CellSet<OC1>;
	template class CellSet<OC2>;
	template class CellSet<OC3>;
	template class CellSet<OC4>;
	template class CellSet<minLB>;
	template class CellSet<weighted_sum>;
}
#endif // __IBEX_CELL_SET_H__
