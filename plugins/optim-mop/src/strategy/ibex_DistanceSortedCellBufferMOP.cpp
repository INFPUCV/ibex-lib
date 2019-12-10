/*
 * ibex_DistanceSortedCellBufferMOP.cpp
 *
 *  Created on: 20 oct. 2017
 *      Author: iaraya
 */

#include "ibex_DistanceSortedCellBufferMOP.h"
#include "ibex_OptimizerMOP.h"
#include <algorithm>    // std::min_element, std::max_element

#ifndef cdata
#define cdata ((BxpMOPData*) c->prop[BxpMOPData::id])
#endif


namespace ibex {



	map< pair <double, double>, IntervalVector >* max_distance::UB=NULL;


	void DistanceSortedCellBufferMOP::flush() {
		while (!cells.empty()) {
			delete pop();
		}
	}

	unsigned int DistanceSortedCellBufferMOP::size() const {
		return cells.size();
	}

	bool DistanceSortedCellBufferMOP::empty() const {
		return cells.empty();
	}

	void DistanceSortedCellBufferMOP::push(Cell* c) {
		if(OptimizerMOP::imprimir)std::cout<<"\t PUSH+++++++++++++++++"<<endl;

		double dist = nds->distance(c->box);

		//std::cout<<"dist = "<<dist<<endl;
		//std::cout <<"ub_distance = "<<cdata->ub_distance<<endl;

		if(dist < cdata->ub_distance){
			cdata->ub_distance=dist;
			//std::cout <<"new ub_distance = "<<cdata->ub_distance<<endl;
		}
		cells.push(c);

		//todo test aux priority queue
//		cells_aux.push(c);

	}

	Cell* DistanceSortedCellBufferMOP::pop() {
		if(OptimizerMOP::imprimir) std::cout<<"\t POP+++++++++++++++++"<<endl;
        Cell* c = top();
        cells.pop();

		return c;
	}

  int counter=0;
	Cell* DistanceSortedCellBufferMOP::top() const {
		if(OptimizerMOP::imprimir) std::cout<<"\t TOP+++++++++++++++++"<<endl;
		Cell* c = cells.top();
		if(!c) return NULL;

		//Obtain dist from cell top
		double dist=nds->distance(c->box);

		//we update the distance and reinsert the element
		while(dist < cdata->ub_distance){
			cells.pop();
			cdata->ub_distance=dist;
			if(OptimizerMOP::imprimir) std::cout<<"PUSH Normal"<<endl;
			cells.push(c);
			c = cells.top();
			dist=nds->distance(c->box);
		}

    counter ++;
//		cout << counter  <<":" <<  c->get<CellMOP>().ub_distance << endl;

		return c;
	}

} // end namespace ibex
