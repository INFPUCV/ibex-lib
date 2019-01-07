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
		return (globalBuffer.empty() && currentBuffer.empty() && nextBuffer.empty());
	}

	void BeamSearchBufferMOP::push(Cell* cell) {
		//cout << nds << endl;
        double dist=nds->distance(cell);
		std::multiset <Cell*>::iterator it;
		//cout << dist << endl;
		int delta=0,i=0;
		if(dist < cell->get<CellMOP>().ub_distance)
		{	

			cell->get<CellMOP>().ub_distance=dist;
		}
 
        //cout << cont << endl;
        if(globalBuffer.empty() && nextBuffer.empty() && cont==0){
			
			globalBuffer.push(cell);
			//cout << "tamaño global cuando los 3 estan vacios" << endl;
			//cout << globalBuffer.size() << endl;
			cont=1;
		
		}else{
			
			nextBuffer.insert(cell);

			//cout << "tamaño next cuando el global tiene elementos " << endl;
			//cout << globalBuffer.size() << endl;
			while(nextBuffer.size()>4){
			/*	it = nextBuffer.end();
				globalBuffer.push(*it);

				nextBuffer.erase(*it);
				JELP
*/
				globalBuffer.push(*nextBuffer.begin())	;
				nextBuffer.erase(nextBuffer.begin());
			}
		}  		
		iter++;
		cout << iter << endl;
	}

	Cell* BeamSearchBufferMOP::pop() {
		Cell* c = NULL;
		std::multiset <Cell*>::iterator it;
		//sacar de current
		if(currentBuffer.empty() && !nextBuffer.empty()){
				
			while(!nextBuffer.empty()){

				//cout << "entrada while" << endl;
				it = nextBuffer.begin();
				currentBuffer.push(*it);	
				nextBuffer.erase(it);
			}
			//cout << "salida while" << endl;
		}
		
		if(currentBuffer.empty() && nextBuffer.empty()){

			c = globalBuffer.top();
			globalBuffer.pop();
					
		}else if(!currentBuffer.empty()){
			
			c = currentBuffer.top();
			currentBuffer.pop();

		}else{
			cout << "error" << endl;
		 	exit;
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
