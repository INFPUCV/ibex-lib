/*
 * ibex_BeamSearchBufferMOP.cpp
 *
 *  Created on: 20 oct. 2017
 *      Author: matias y pablo
 */

#include "ibex_BeamSearchBufferMOP.h"
#include "ibex_OptimizerMOP.h"
#include <algorithm>    // std::min_element, std::max_element

namespace ibex {

	void BeamSearchBufferMOP::flush() {
		while (!globalBuffer.empty()) {
			delete pop();
		}
	}

	unsigned int BeamSearchBufferMOP::size() const {
		return globalBuffer.size();
	}

	bool BeamSearchBufferMOP::empty() const {
		return globalBuffer.empty();
	}

	void BeamSearchBufferMOP::push(Cell* cell) {
		cout << nds << endl;
        double dist=nds->distance(cell);
		cout << dist << endl;
		int delta=0,i=0;
		if(dist < cell->get<CellMOP>().ub_distance)
		{	

			cell->get<CellMOP>().ub_distance=dist;
		}
 
        
        if(globalBuffer.empty() && currentBuffer.empty() && nextBuffer.empty()){
		
			globalBuffer.push(cell);
		
		}else{

			nextBuffer.insert(cell);
			if(nextBuffer.size()>4){
	
				globalBuffer.push(*nextBuffer.end());
	
			}
		}  
  

		

	}

	Cell* BeamSearchBufferMOP::pop() {
		Cell* c = NULL;
		std::multiset <Cell*>::iterator it;
		//sacar de current

		if(currentBuffer.empty() && !nextBuffer.empty()){

				for(it = nextBuffer.begin(); it != nextBuffer.end(); ++it){
					
					currentBuffer.push(*it);
					
					it = nextBuffer.erase(it);
					
				}
				cout << "salida for" << endl;
			}else if(currentBuffer.empty() && nextBuffer.empty()){	
				
					c = globalBuffer.top();
					
					globalBuffer.pop();
					
			}else if(!currentBuffer.empty()){
					
					c = currentBuffer.top();
					currentBuffer.pop();

			} 
	
		return c;
	}

  int counter1=0;
	Cell* BeamSearchBufferMOP::top() const {

		Cell* c = globalBuffer.top();
		if(!c) return NULL;

		if (OptimizerMOP::_hv) return c;

		double dist=nds->distance(c);

		//we update the distance and reinsert the element
		while(dist < c->get<CellMOP>().ub_distance){
			globalBuffer.pop();
			c->get<CellMOP>().ub_distance=dist;
			globalBuffer.push(c);
			c = globalBuffer.top();
			dist=nds->distance(c);
		}

    	counter1 ++;
		//cout << counter1  <<":" <<  c->get<CellMOP>().ub_distance << endl;

		return c;
	}

} // end namespace ibex
