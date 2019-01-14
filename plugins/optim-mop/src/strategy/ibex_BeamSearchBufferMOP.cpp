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

	int BeamSearchBufferMOP::nextBufferSize = 4;

	void BeamSearchBufferMOP::flush() {
		while (!globalBuffer.empty()) {
			delete pop();
		}
		while (!currentBuffer.empty()) {
			delete pop();
		}
		while (!nextBuffer.empty()) {
			nextBuffer.erase(nextBuffer.begin());
		}
	}

	unsigned int BeamSearchBufferMOP::size() const {
		return (globalBuffer.size()+currentBuffer.size()+nextBuffer.size());
	}

	bool BeamSearchBufferMOP::empty() const {
		return (globalBuffer.empty() && currentBuffer.empty() && nextBuffer.empty());
	}

	void BeamSearchBufferMOP::push(Cell* cell) {
		
        double dist=nds->distance(cell);
		std::multiset <Cell*>::iterator it;

		int delta=0,i=0;
		if(dist < cell->get<CellMOP>().ub_distance)
		{	

			cell->get<CellMOP>().ub_distance=dist;
		}
 
        //primera iteracion
        if(globalBuffer.empty() && nextBuffer.empty() && cont==0){
			
			globalBuffer.push(cell);
			//se cambia el valor del flag cont para no entrar nuevamente
			cont=1;
		
		}else{
			
			nextBuffer.insert(cell);

			//Si el NextBuffer sobrepasa la capacidad maxima, se borran y se mueven al global
			//cout << "size next antes while: " << nextBuffer.size() << endl;
			while(nextBuffer.size()> nextBufferSize){
				
				globalBuffer.push(*nextBuffer.begin());
				nextBuffer.erase(nextBuffer.begin());
				// cout << "size next en while: " << nextBuffer.size() << endl;
				// cout << "size global en while: " << globalBuffer.size() << endl;
			}
		}  	

		// cout << "tamaño global: " << globalBuffer.size() << endl;
		// cout << "tamaño current: " << currentBuffer.size() << endl;
		// cout << "tamaño next: " << nextBuffer.size() << endl;
		//getchar();	
	}

	Cell* BeamSearchBufferMOP::pop() {
		Cell* c = NULL;
		std::multiset <Cell*>::iterator it;
		
		//SI el current esta vacio y el next tiene elementos, se pasan del next al current
		if(currentBuffer.empty() && !nextBuffer.empty()){
			//cantBeam++;
			//cout << "BeamSearch: " << cantBeam << endl;

			while(!nextBuffer.empty()){

				//it = nextBuffer.begin();
				//double distNextBegin = nds->distance(*nextBuffer.begin());
				//cout << "distancia primero: " << distNextBegin << endl;
				/*if(nextBuffer.size()>1){
					it++;
					double distNextEnd = nds->distance(*it);
					//cout << "distancia siguiente: " << distNextEnd << endl;
				}*/

				currentBuffer.push(*nextBuffer.begin());	
				nextBuffer.erase(nextBuffer.begin());
				
			}
		}
		
		//Si current y next estan vacios, se popea del global
		if(currentBuffer.empty() && nextBuffer.empty() && !globalBuffer.empty()){

			c = globalBuffer.top();
			globalBuffer.pop();
			//cantBeam++;
			//cout << "BeamSearch: " << cantBeam << endl;
			//int p = c->get<CellMOP>().depth;
			//cout << "Profundidad: " << p << endl;
			//getchar();
					
		}else if(!currentBuffer.empty()){
			//si current tiene elementos, siempre se sacan de current
			c = currentBuffer.top();
			currentBuffer.pop();

		}else{
			cout << "error" << endl;
		 	exit;
		} 
		// cout << "tamaño next 2: " << nextBuffer.size() << endl;
		// cout << "tamaño current 2: " << currentBuffer.size() << endl;
		//getchar();
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
