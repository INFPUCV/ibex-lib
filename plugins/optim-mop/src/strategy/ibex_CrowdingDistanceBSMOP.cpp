/*
 * ibex_CrowdingDistanceBSMOP.cpp
 *
 *  Created on: 20 oct. 2017
 *     Authors: matias y pablo
 */

#include "ibex_CrowdingDistanceBSMOP.h"
#include "ibex_OptimizerMOP.h"
#include "ibex_NDS.h"
#include <algorithm>    // std::min_element, std::max_element


namespace ibex { 

	int CrowdingDistanceBSMOP::nn = 0;


	void CrowdingDistanceBSMOP::flush() {
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

	unsigned int CrowdingDistanceBSMOP::size() const {
		return (globalBuffer.size()+currentBuffer.size()+nextBuffer.size());
	}

	void CrowdingDistanceBSMOP::bsPerformance(int iter){

	}

	bool CrowdingDistanceBSMOP::empty() const {
		return (globalBuffer.empty() && currentBuffer.empty() && nextBuffer.empty());
	}

	void CrowdingDistanceBSMOP::push(Cell* cell) {
		
		if(global_hv) {
			aux=nds->hypervolume(CellMOP::y1_init,CellMOP::y2_init).mid();

			initial_reduction=aux-aux2;
			//cout << initial_reduction << endl;

			//cout << "hv global:" << initial_reduction << endl; 
			global_hv=false;

		}
        double dist=nds->distance(cell);
		std::multiset <Cell*>::iterator it;

		int delta=0,i=0;
		if(dist < cell->get<CellMOP>().ub_distance)
		{	

			cell->get<CellMOP>().ub_distance=dist;
		}
 
        //primera iteracion, se usa un flag porque para la segunda iteracion vuelven a estar todos vacios, por lo que en vez
		//de hacer el push en el next, lo vuelve a hacer en el global
        if(globalBuffer.empty() && nextBuffer.empty() && /*currentBuffer.empty()*/ cont==0){
			
			globalBuffer.push(cell);
			//se cambia el valor del flag cont para no entrar nuevamente
			cont=1;
		
		}else{
			
			//next no tiene tama��o maximo
			nextBuffer.insert(cell);
			Cell* c=NULL;

			//cout << "iter en bs: " << iterBS << endl;
			//cout << T << endl;
			//cout << iterBS%T << endl;
			/*if(iterBS%T==0){
				currentBufferSizeAux=currentBufferMaxSize;
			}else{
				currentBufferSizeAux=0;
			}*/
		}  		
	}

	//NDS + crowding distance se realiza cuando no hay elementos en el currentBuffer
	Cell* CrowdingDistanceBSMOP::pop() {

		//cout << iterBS << endl;
		iterBS++;

		Cell *c = NULL, *c2 = NULL;
		std::multiset <Cell*>::iterator it;
		
		//SI el current esta vacio y el next tiene elementos, se pasan del next al current
		if(currentBuffer.empty() && !nextBuffer.empty()){
			depth++;
			if(depth>*depthMayor){
				*depthMayor=depth;
			}
			depthTotal++;//=depthTotal+depth;
			//	getchar();

			c=*nextBuffer.begin();
			++c;
			c2=c;
			while(!nextBuffer.empty()){
				if(crowding){
					

					//Esto ahora no se usa, se utiliza el nonDominatedSorting en la funcion del crowding distance
					//if(nextBuffer.size()>4)removeDominated(nextBuffer, globalBuffer);

					std::cout << "BEGIN: " << endl;
					for(auto c:nextBuffer) cout << c->box[c->box.size()-1].lb() << "," << c->box[c->box.size()-2].lb() << endl;
					//if(nextBuffer.size() >= nextBufferSize){

					int sizeDiffMax = currentBufferMaxSize - currentBuffer.size();
					//cout << currentBuffer.size() + nextBuffer.size()<< endl;
					//getchar();

					nonDominatedSort(nextBuffer,currentBuffer, globalBuffer, currentBufferMaxSize); //al current

				}else{
					while(currentBuffer.size()<=currentBufferMaxSize){
						c=*nextBuffer.begin();
						currentBuffer.push(*nextBuffer.begin());
						nextBuffer.erase(nextBuffer.begin());
					}
				}
			}

			//*****esto quizas puede quedar igual
			/*if(!currentBuffer.empty()){
				//intento de hv
				//aux2 es antes de trabajar la caja y aux despues de trabajar la caja


				aux=nds->hypervolume(CellMOP::y1_init,CellMOP::y2_init).mid();
				//cout <<  aux << endl;
				if(aux-aux2!=0){
					
					if(initial_reduction==0) initial_reduction=aux-aux2;

					double mejora=(aux-aux2)/initial_reduction;

					if(mejora>=bs_tolerance*mejor_mejora || aux-aux2==initial_reduction){

						if(mejora>mejor_mejora) {
							if(mejor_mejora!=0.0) {T=1; bs_performance=true;}
							else bs_performance=false;

							mejor_mejora=mejora;
						}



					//if(mejora>=errorBS || aux-aux2==initial_reduction){
						//errorBS=mejora;

					}else{
						if(!bs_performance && currentBufferSizeAux!=0){
							T=T*2;
						}

						if(!currentBuffer.empty()){
							while(!currentBuffer.empty()){
				
								globalBuffer.push(currentBuffer.top());
								currentBuffer.pop();

							}
						}
						//en teoria nunca entro aqui
						if(!nextBuffer.empty()){
							while(!nextBuffer.empty()){
				
								c=*nextBuffer.begin();
								globalBuffer.push(*nextBuffer.begin());
					
								nextBuffer.erase(nextBuffer.begin());

							}
						}
					}
				}

				aux2=aux;

			}*/

		}
		
		//Si current y next estan vacios, se popea del global
		if(currentBuffer.empty() && nextBuffer.empty() && !globalBuffer.empty()){

			depth=0;
			mejor_mejora=0;
			aux2=nds->hypervolume(CellMOP::y1_init,CellMOP::y2_init).mid();


			c=top();
			//c = globalBuffer.top();
			global_hv=true;

			cout << globalBuffer.size() << endl;
			//AQUI SE CAE
			globalBuffer.pop();
			cantBeam++;

			if(cantBeam!=0 && depthTotal!=0){
				depthPromedio=depthTotal/(double)cantBeam;
				*pruebaprom=depthPromedio;
				//cout << *pruebaprom <<endl;
			} 
					
		}else if(!currentBuffer.empty()){

			//si current tiene elementos, siempre se sacan de current
			c = currentBuffer.top();
			currentBuffer.pop();
			//cout << "llego aca? " << currentBuffer.size() << endl;

		}else{
			cout << "error" << endl;
		} 
		return c;
	}

	//no deberia ser mas complejo que solo hacerlo con el global?
	Cell* CrowdingDistanceBSMOP::top() const {

		Cell* c = globalBuffer.top();

		if (OptimizerMOP::_hv) return c;

		cout << "aqui es" << endl;
		cout << globalBuffer.size() << endl;

		//AQUI SE CAE
		double dist=nds->distance(c);
		//we update the distance and reinsert the element

		while(dist < c->get<CellMOP>().ub_distance){
			globalBuffer.pop();
			c->get<CellMOP>().ub_distance=dist;
			globalBuffer.push(c);
			c = globalBuffer.top();
			dist=nds->distance(c);
		}

		//cout << counter1  <<":" <<  c->get<CellMOP>().ub_distance << endl;

		return c;
	}

	void CrowdingDistanceBSMOP::crowdingDistance(
	        std::multiset<Cell*, max_distanceCrowding>& nonDominated,
	        std::priority_queue<Cell*, std::vector<Cell*>, max_distanceCrowding >& currentBuffer,
	        std::priority_queue<Cell*, std::vector<Cell*>, max_distanceCrowding >& globalBuffer,
	        int currentBufferMaxSize
	        ){

		Cell *c = NULL;	

        int i = 0; //%TEMP% Used just to print stuff
        //First, we take out the dominated ones.

        int returnSize=currentBufferMaxSize - currentBuffer.size();

        //removeDominated(nextBuffer, currentBuffer); //We do this now in the BeamSearchBufferMOP Class.
		std::multiset<CDBox*, sortByCrowdingDistance> cdSet;  //We now declare this before the while
        while(returnSize>0 && nonDominated.size()>0){
            
			cdSet.clear(); //We clear cdSet each iteration.

            for(std::multiset<Cell*, crowding_distanceBeam>::iterator it = nonDominated.begin();
            it != nonDominated.end();
            it++)
            {

                CDBox* cdBox= (CDBox*)malloc(sizeof(CDBox));
                cdBox->C = *it;

                //If they are the first/last one, set their Crowding Distances to Infinite
                if(*it == *nonDominated.begin() || *next(it) == *nonDominated.end()){
                    cdBox->crowding_distance = std::numeric_limits<double>::infinity();
                }
                else{ //If not, then calculate them.
                    int n = (*it)->box.size();
                    Cell* first = *nonDominated.begin();
                    Cell* last = *prev(nonDominated.end());


                    // Here we calculate the crowding distance. We won't use a for, since we already know that there's only 2 objectives (y1, y2)
                    //Primer Objetivo
                    cdBox->crowding_distance += ((*std::prev(it))->box[n-1].lb() - (*std::next(it))->box[n-1].lb())/((first->box[n-1].lb() - last->box[n-1].lb()));
                    //Segundo objetivo
                    cdBox->crowding_distance += ((*std::prev(it))->box[n-2].lb() - (*std::next(it))->box[n-2].lb())/((first->box[n-2].lb() - last->box[n-2].lb()));
                }

                //Inserts the CDBox node into the set
                cdSet.insert(cdBox);
            }
            
			//////////////////////////////
			//OLD METHOD:
            //Here we return the first element in the multiset, we do this 'returnSize' times, we calculate the crowding distance every time we do the loop.
            /*CDBox* cdbox = *(cdSet.begin());
            currentBuffer.push(cdbox->C);
            nextBuffer.erase(cdbox->C);*/
			//////////////////////////////

			/* New Method:
			*	We "remove" only the last (lowest Crowding Distance)
			*	and then it's sent into the globalBuffer.
			*/
			CDBox* cdbox = *(--cdSet.end());

			globalBuffer.push(cdbox->C); //Inserts the cell into the global Buffer.

            std::cout << i << ", " << cdbox->crowding_distance << ", "
                    << cdbox->C->box[cdbox->C->box.size()-1].lb() << "," << cdbox->C->box[cdbox->C->box.size()-2].lb()<< endl;
            i++;

            returnSize--; //If the returnSize == 0, then the cdSet isn't cleared.
        }
		//If we have stuff in the cdSet, it means that we need to put them into the nextBuffer!
		//Since the lowest crowding distance ones are already in the globalBuffer, we won't need to remove them from here. We will just have to add stuff into the nextBuffer.
		while(!cdSet.empty()){
			CDBox* cdbox = *(cdSet.begin());
			currentBuffer.push(cdbox->C);
			cdSet.erase(cdSet.begin());
		}
		//After putting everything into the currentBuffer, we clear cdSet.
		cdSet.clear(); //I'm not sure if this deletes the pointers from everything, TODO: CHECK that.



	}

	void CrowdingDistanceBSMOP::nonDominatedSort(
        std::multiset<Cell*, max_distanceCrowding>& nextBuffer,
        std::priority_queue<Cell*, std::vector<Cell*>, max_distanceCrowding >& currentBuffer,
        std::priority_queue<Cell*, std::vector<Cell*>, max_distanceCrowding >& globalBuffer,
        int currentBufferMaxSize
        ){

		int returnSize = currentBufferMaxSize - currentBuffer.size();

		//nonDominatedSort

		//1. en cada iteracion se extraen los nodos no dominados del nextBuffer y se guardan en nonDominated
		//2. se verifica que los no dominados quepan en currentBuffer, si caben se vuelve a 1
		//3. Si no caben se llama al crowdingDistance para filtrar el nonDominated
		while(returnSize>0 && nextBuffer.size()>0){

			std::multiset<Cell*, max_distanceCrowding> nonDominated;
			std::list< std::multiset<Cell*, max_distanceCrowding>::iterator >::iterator findIter;
			extractNonDominated(nextBuffer, nonDominated);

			//Si la cantidad de no dominados es menor o igual que el tamaño disponible del current, se pasan todos y se borran del next buffer
			if(nonDominated.size() < returnSize){
				for(auto el:nonDominated)
					currentBuffer.push(el);

				//sino, en el conjunto no dominado extraído, hay que realizar el crowding distance hasta tener la cantidad
				//que necesitamos, los sobrantes se van al global
			}else
				crowdingDistance(nonDominated,currentBuffer,globalBuffer,currentBufferMaxSize);


			returnSize = currentBufferMaxSize - currentBuffer.size();
		}

		for(auto el:nextBuffer)
			globalBuffer.push(el);

		nextBuffer.clear();
    }

    bool CrowdingDistanceBSMOP::isDominated(Cell* a, Cell* b){
        int n = a->box.size();
        return (a->box[n-1].lb() <= b->box[n-1].lb() && a->box[n-2].lb() <= b->box[n-2].lb());
    }

    //Mueve los elementos no dominados de nextBuffer a nonDominated
	//se podria mandar como parametro el tamaño que admite el current

	//Por lo que estoy viendo, se pasan siempre todas las cajas del next al nonDominated
    void CrowdingDistanceBSMOP::extractNonDominated(
    std::multiset<Cell*, max_distanceCrowding>& nextBuffer, std::multiset<Cell*, max_distanceCrowding>& nonDominated){

    	std::multiset<Cell*, max_distanceCrowding>::iterator it=nextBuffer.begin();

		cout << "nextbuffer: " << endl;
		for(auto b : nextBuffer){
            cout << b->box[nn-1] << "\n" << b->box[nn-2] << "\n";
        }

        for(; it!=nextBuffer.end(); ){
        	Cell* a=*it;
        	bool dominated=false;
            for(auto b : nextBuffer){
                if(a != b && isDominated(a, b)){
                	dominated=true;
                	break;
                }
            }

            if(!dominated){
            	nonDominated.insert(*it);
            	it = nextBuffer.erase(it);
            }else
            	it++;
        }

		cout << "nondominated: " << endl;
		for(auto b : nonDominated){
            cout << b->box[nn-1] << "\n" << b->box[nn-2] << "\n";
        }

		cout << "nextbuffer (de nuevo): " << endl;
		for(auto b : nextBuffer){
            cout << b->box[nn-1] << "\n" << b->box[nn-2] << "\n";
        }
		cout << nextBuffer.size() << endl;

    }

    //Removes the dominated cells from nextBuffers and adds them into globalBuffer.
    void CrowdingDistanceBSMOP::removeDominated(
    std::multiset<Cell*, max_distanceCrowding>& nextBuffer,
    std::priority_queue<Cell*, std::vector<Cell*>, max_distanceCrowding >& globalBuffer
    ){
        //std::multiset<Cell*, max_distanceCrowding> nextBufferCopy;
        //Here we take the dominateds out of nextBuffer and insert them into globalBuffer.

		for(auto itA = nextBuffer.begin(); itA != nextBuffer.end(); ++itA){
			for(auto itB = nextBuffer.begin(); itB != nextBuffer.end(); ++itB){
				if(itA != itB && isDominated(*itA, *itB)){
					globalBuffer.push(*itA);
					nextBuffer.erase(*itA++);
					break;
				}
			}
		}

    }

} // end namespace ibex
