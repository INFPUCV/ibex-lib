/*
 * ibex_BeamSearchBufferMOP.cpp
 *
 *  Created on: 20 oct. 2017
 *     Authors: matias y pablo
 */

#include "ibex_BeamSearchBufferMOP.h"
#include "ibex_OptimizerMOP.h"
#include "ibex_CrowdingDistance.h"
#include <algorithm>    // std::min_element, std::max_element


namespace ibex { 

	int BeamSearchBufferMOP::currentBufferMaxSize = 8;
	int BeamSearchBufferMOP::nn = 0;

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
			Cell* c=NULL;

		}  	

	}

	Cell* BeamSearchBufferMOP::pop() {
		Cell *c = NULL, *c2 = NULL;
		std::multiset <Cell*>::iterator it;
		
		//SI el current esta vacio y el next tiene elementos, se pasan del next al current
		//TODO: aqui se debe aplicar non dominated sorting + crowding distance
		if(currentBuffer.empty() && !nextBuffer.empty()){

			c=*nextBuffer.begin();
			++c;
			c2=c;

			while(!nextBuffer.empty()){
				std::cout << "BEGIN: " << endl;
				for(auto c:nextBuffer) cout << c->box[c->box.size()-1].lb() << "," << c->box[c->box.size()-2].lb() << endl;

				NyuCrowdingDistance::getCrowdingDistance(nextBuffer, currentBuffer, globalBuffer, currentBufferMaxSize); //al current



				std::cout << "getCrowdingDistance done" << endl;
			}


		}
		
		//Si current y next estan vacios, se popea del global
		if(currentBuffer.empty() && nextBuffer.empty() && !globalBuffer.empty()){
			
			// Reset de archivo
			ofstream myfile4;
			myfile4.open("global.txt");

			c = globalBuffer.top();

			myfile4 << c->box[nn-1] << "\n" << c->box[nn-2] << "\n";
			myfile4.close();

			
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
		// cout << "tama������o next 2: " << nextBuffer.size() << endl;
		// cout << "tama������o current 2: " << currentBuffer.size() << endl;
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
