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
#include <iostream>
#include <string>




namespace ibex { 

	int CrowdingDistanceBSMOP::bs_level=0;


	void CrowdingDistanceBSMOP::flush() {
		while (!globalBuffer.empty()) {
			delete pop();
		}
		while (!currentBuffer.empty()) {
			delete pop();
		}
		nextBuffer.clear();
		intersectedCells.clear();
		g_pausedCells.clear();
		c_pausedCells.clear();
		n_pausedCells.clear();
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
		if(bs_level==0) {
			
			//hv=nds->hypervolume(CellMOP::y1_init,CellMOP::y2_init).mid();
			initial_reduction=hv-hv0;
			//cout << initial_reduction << endl;
			ofstream myfile;
			int n = cell->box.size();
  			myfile.open ("global.txt");//, std::ios_base::app);
			myfile << cell->box[n-1] << "\n" << cell->box[n-2] << "\n";
			myfile.close();

			//getchar();

		}

		/*if(initial_reduction<0){
			cout << "error" << endl;
			getchar();
		}*/

        double dist=nds->distance(cell);

		if(dist < cell->get<CellMOP>().ub_distance)
		{	

			cell->get<CellMOP>().ub_distance=dist;
		}
 
		nextBuffer.insert(cell);

		/*if(iterBS%T==0){
			currentBufferSizeAux=currentBufferMaxSize;
		}else{
			currentBufferSizeAux=0;
		}*/
		 		
	}

	//NDS + crowding distance se realiza cuando no hay elementos en el currentBuffer
	Cell* CrowdingDistanceBSMOP::pop() {

		Cell *c,*c2 = NULL;
		
		//SI el current esta vacio y el next tiene elementos, se pasan del next al current
		if(currentBuffer.empty() && !nextBuffer.empty() /*&& bs_state*/){
			depth++;
			if(depth>*depthMayor){
				*depthMayor=depth;
			}
			depthTotal++;
			bs_level++;
			if(nextBuffer.size()<=currentBufferMaxSize){
				while(!nextBuffer.empty()){
					c=*nextBuffer.begin();
					currentBuffer.push(*nextBuffer.begin());
					nextBuffer.erase(nextBuffer.begin());
				}
			}else{
				nonDominatedSort(nextBuffer,currentBuffer, globalBuffer, currentBufferMaxSize); //al current
				if(nextBuffer.size()>0){
					for(auto el:nextBuffer)
						globalBuffer.push(el);

					nextBuffer.clear();
				}
			}

		}
		
		/*if(bs_level>maxBSHeight){
			while(!nextBuffer.empty()){
				c=*nextBuffer.begin();
				globalBuffer.push(*nextBuffer.begin());
				nextBuffer.erase(nextBuffer.begin());
			}
			while(!currentBuffer.empty()){

				globalBuffer.push(currentBuffer.top());
				currentBuffer.pop();

			}
		}*/

		double x1,x2,y1,y2;
		std::ifstream selectedSpace;
		selectedSpace.open("selectedSpace.txt");
 		char output[100];
		vector<string> points;
 		if (selectedSpace.is_open()) {
 			while (!selectedSpace.eof()) {


    			selectedSpace >> output;
   				//cout<<output;
				char * pch;
     			pch = strtok (output," ,");
      			while (pch != NULL){
        			//printf ("%s\n",pch);
					points.push_back(pch);
        			pch = strtok (NULL, " , ");
      			}


 			}
		}
		selectedSpace.close();

		//for (int i=0; i<points.size(); i++)     
        //	cout << points[i] << "\n";

		if(points.size()==4){
			//cout << "inter " << intersectedCells.size() << endl;
			//cout << "g_paused " << g_pausedCells.size() << endl;
			//cout << "c_paused " << c_pausedCells.size() << endl;
			//cout << "n_paused " << n_pausedCells.size() << endl;
			//cout << "global " << globalBuffer.size() << endl;
			//cout << "current " << currentBuffer.size() << endl;
			//cout << "next " << nextBuffer.size() << endl;
			//cout << "dentro if point size" << endl;
			x1 = ::atof(points[0].c_str());
			x2 = ::atof(points[1].c_str());
			y1 = ::atof(points[2].c_str());
			y2 = ::atof(points[3].c_str());
			//cout << "nextbuffer size: " << nextBuffer.size() << endl;
			//cout << "current size: " << currentBuffer.size() << endl;
			Point l2 = {x1, y2}, r2 = {x2, y1};
			for(auto nextrec: nextBuffer){
				cout << "recorriendo next" << endl;
				int n = nextrec->box.size();
				Point l1 = {nextrec->box[n-1].lb(), nextrec->box[n-2].ub()}, r1 = {nextrec->box[n-1].ub(), nextrec->box[n-2].lb()}; 				
				 
				//cout << l1.x << "," << l1.y << endl;
				//cout << r1.x << "," << r1.y << endl;
				//cout << l2.x << "," << l2.y << endl;
				//cout << r2.x << "," << r2.y << endl;
				if (doOverlap(l1, r1, l2, r2)){
					//cout << "Rectangles Overlap" << endl;
					intersectedCells.insert(nextrec);
					//getchar();
				}else{
					//cout << "Rectangles don't Overlap" << endl;
					n_pausedCells.insert(nextrec);
				}	
			}
			nextBuffer.clear();

			while(!currentBuffer.empty()){
				Cell *currentrec= NULL;
				cout << "recorriendo current" << endl;
				currentrec=currentBuffer.top();
				currentBuffer.pop();
				//cout << c->box << endl;
				int n = currentrec->box.size();
				Point l1 = {currentrec->box[n-1].lb(), currentrec->box[n-2].ub()}, r1 = {currentrec->box[n-1].ub(), currentrec->box[n-2].lb()}; 

				if (doOverlap(l1, r1, l2, r2)){
					//cout << "Rectangles Overlap" << endl;
					intersectedCells.insert(currentrec);
					//getchar();
				}else{
					//cout << "Rectangles don't Overlap" << endl;
					c_pausedCells.insert(currentrec);
				}

			}

			while(!globalBuffer.empty()){
				Cell *globalrec= NULL;
				cout << "recorriendo global" << endl;
				globalrec=globalBuffer.top();
				globalBuffer.pop();
				int n = globalrec->box.size();
				Point l1 = {globalrec->box[n-1].lb(), globalrec->box[n-2].ub()}, r1 = {globalrec->box[n-1].ub(), globalrec->box[n-2].lb()}; 

				if (doOverlap(l1, r1, l2, r2)){
					//cout << "Rectangles Overlap" << endl;
					intersectedCells.insert(globalrec);
					//getchar();
				}else{
					//cout << "Rectangles don't Overlap" << endl;
					g_pausedCells.insert(globalrec);

				}
			}

			if(!intersectedCells.empty()){
				cout << "intersected not empty" << endl;
				while (!globalBuffer.empty()) {
					globalBuffer.pop();
				}
				while (!currentBuffer.empty()) {
					currentBuffer.pop();
				}
				nextBuffer.clear();
			}else{
				cout << "intersected empty" << endl;
				for(auto el:g_pausedCells)
					globalBuffer.push(el);
				g_pausedCells.clear();
				
				for(auto el:c_pausedCells)
					currentBuffer.push(el);
				c_pausedCells.clear();

				for(auto el:n_pausedCells)
					nextBuffer.insert(el);
				n_pausedCells.clear();
			}
		
  		}

		for(auto caja: nextBuffer){

		}

		while(!intersectedCells.empty()){
			c=*intersectedCells.begin();
			intersectedCells.erase(intersectedCells.begin());
			c2=*intersectedCells.begin();
			*c.box&=*c2.box;
		}

		//cout << "inter " << intersectedCells.size() << endl;
		//cout << "g_paused " << g_pausedCells.size() << endl;
		//cout << "c_paused " << c_pausedCells.size() << endl;
		//cout << "n_paused " << n_pausedCells.size() << endl;
		//cout << "global " << globalBuffer.size() << endl;
		//cout << "current " << currentBuffer.size() << endl;
		//cout << "next " << nextBuffer.size() << endl;

		//cout << "Value of x1 : " << x1 << endl; 
		//cout << "Value of x2 : " << x2 << endl; 
		//cout << "Value of y1 : " << y1 << endl; 
		//cout << "Value of y2 : " << y2 << endl; 

		//cout << "pausa" << endl;
		//getchar();








		//Si el archivo contiene algo
		//if(selectedSpace.peek() != std::ifstream::traits_type::eof()){
			
		//}



		//Si current y next estan vacios, se popea del global
		if(currentBuffer.empty() && nextBuffer.empty() && !globalBuffer.empty()){

			depth=0;
			mejor_mejora=0;
			//hvE=nds->hypervolume(CellMOP::y1_init,CellMOP::y2_init).mid();

			mejora=hvE-hv0;
			/*if(mejora>0 && initial_reduction/mejora){
				T=1;
			}else{
				T=T*2;
			}*/

			if(iterBS%T==0){
				bs_state=true;
			}else{
				bs_state=false;
			}

			c=top();
			//c = globalBuffer.top();

			bs_level=0;
			iterBS++;
			//hv0=nds->hypervolume(CellMOP::y1_init,CellMOP::y2_init).mid();



			globalBuffer.pop();
			cantBeam++;
			if(cantBeam!=0 && depthTotal!=0){
				depthPromedio=depthTotal/(double)cantBeam;
				*pruebaprom=depthPromedio;
			} 
					
		}else if(!currentBuffer.empty()){

			//si current tiene elementos, siempre se sacan de current
			c = currentBuffer.top();
			currentBuffer.pop();

		}else if(!intersectedCells.empty()){

			c = *intersectedCells.begin();
			intersectedCells.erase(intersectedCells.begin());

			if(intersectedCells.empty()){

				for(auto el:g_pausedCells)
					globalBuffer.push(el);
				g_pausedCells.clear();
				
				for(auto el:c_pausedCells)
					currentBuffer.push(el);
				c_pausedCells.clear();

				for(auto el:n_pausedCells)
					nextBuffer.insert(el);
				n_pausedCells.clear();

			}

		}
		return c;
	}

	//no deberia ser mas complejo que solo hacerlo con el global?
	Cell* CrowdingDistanceBSMOP::top() const {

		Cell* c = globalBuffer.top();

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

		//cout << counter1  <<":" <<  c->get<CellMOP>().ub_distance << endl;

		return c;
	}

	void CrowdingDistanceBSMOP::nonDominatedSort(
        std::multiset<Cell*, max_distanceCrowding>& nextBuffer,
        std::priority_queue<Cell*, std::vector<Cell*>, max_distanceCrowding >& currentBuffer,
        std::priority_queue<Cell*, std::vector<Cell*>, max_distanceCrowding >& globalBuffer,
        int currentBufferMaxSize
        ){
		std::multiset<Cell*, max_distanceCrowding> nonDominated;
		//nonDominatedSort

		//1. en cada iteracion se extraen los nodos no dominados del nextBuffer y se guardan en nonDominated
		//2. se verifica que los no dominados quepan en currentBuffer, si caben se vuelve a 1
		//3. Si no caben se llama al crowdingDistance para filtrar el nonDominated
		while(currentBuffer.size()<currentBufferMaxSize && nextBuffer.size()>0){
			
			extractNonDominated(nextBuffer, nonDominated);

			//Si la cantidad de no dominados es menor o igual que el tama������o disponible del current, se pasan todos y se borran del next buffer
			if(nonDominated.size()>0 && (nonDominated.size()+currentBuffer.size())<=currentBufferMaxSize){
				for(auto el:nonDominated)
					currentBuffer.push(el);

				nonDominated.clear();
				//sino, en el conjunto no dominado extra������do, hay que realizar el crowding distance hasta tener la cantidad
				//que necesitamos, los sobrantes se van al global
			}else{
				crowdingDistance(nonDominated,currentBuffer,globalBuffer,currentBufferMaxSize);
			}
		}

		
    }

    void CrowdingDistanceBSMOP::extractNonDominated(
    std::multiset<Cell*, max_distanceCrowding>& nextBuffer, std::multiset<Cell*, max_distanceCrowding>& nonDominated){

    	std::multiset<Cell*, max_distanceCrowding>::iterator it=nextBuffer.begin();

        for(; it!=nextBuffer.end(); ){
        	Cell* a=*it;
        	bool dominated=false;
            for(auto b : nextBuffer){
				int x = a->box.size();
				int y = b->box.size();
				//cout << "A: " << a->box[x-1].lb() << " , " << a->box[x-2].lb() << "\n";
				//cout << "B: " << b->box[y-1].lb() << " , " << b->box[y-2].lb() << "\n";
                if(a != b && dominates(b, a)){
					//cout << "dominated true" << endl;
					//getchar();
                	dominated=true; 
                	break;
                }
				//getchar();
            }

            if(!dominated){
            	nonDominated.insert(*it);
            	it = nextBuffer.erase(it);
            }else{
            	it++;
			}
        }

		if(nonDominated.size()==0){
			Cell* c=NULL;
            while(!nextBuffer.empty()){
				c=*nextBuffer.begin();
				nonDominated.insert(*nextBuffer.begin());
				nextBuffer.erase(nextBuffer.begin());
			}
		}

		/*cout << "nondominated: " << endl;
		for(auto b : nonDominated){
			int n = b->box.size();
            cout << b->box[n-1].lb() << " , " << b->box[n-2].lb() << "\n";
        }

		cout << "nextbuffer (de nuevo): " << endl;
		for(auto b : nextBuffer){
			int n = b->box.size();
            cout << b->box[n-1].lb() << " , " << b->box[n-2].lb() << "\n";
        }
		cout << nextBuffer.size() << endl;*/
		//getchar();

    }

	void CrowdingDistanceBSMOP::crowdingDistance(
	        std::multiset<Cell*, max_distanceCrowding>& nonDominated,
	        std::priority_queue<Cell*, std::vector<Cell*>, max_distanceCrowding >& currentBuffer,
	        std::priority_queue<Cell*, std::vector<Cell*>, max_distanceCrowding >& globalBuffer,
	        int currentBufferMaxSize
	        ){
		
		Cell *c = NULL;	
		std::multiset <CDBox*>::iterator it;
		std::multiset<CDBox*, sortByCrowdingDistance> cdSet;

        int i = 0; //%TEMP% Used just to print stuff
        //First, we take out the dominated ones.

        int returnSize=nonDominated.size() - cdSet.size();

        while(nonDominated.size()>currentBufferMaxSize && nonDominated.size()>0){
            
			//cdSet.clear(); //We clear cdSet each iteration.
			std::multiset<Cell*, crowding_distanceBeam>::iterator erase_it = nonDominated.end();

            for(std::multiset<Cell*, crowding_distanceBeam>::iterator it = nonDominated.begin();
            it != nonDominated.end();
            it++)
            {
            	double min_cd=POS_INFINITY;

                //CDBox* cdBox= (CDBox*)malloc(sizeof(CDBox)); new CDBox();
                //cdBox->C = *it;
            	double crowding_distance=0.0;

                //If they are the first/last one, set their Crowding Distances to Infinite
                if(*it == *nonDominated.begin() || *next(it) == *nonDominated.end()){
                    crowding_distance = POS_INFINITY; //std::numeric_limits<double>::infinity();
                }
                else{ //If not, then calculate them.
                    int n = (*it)->box.size();
                    Cell* first = *nonDominated.begin();
                    Cell* last = *prev(nonDominated.end());

                    // Here we calculate the crowding distance. We won't use a for, since we already know that there's only 2 objectives (y1, y2)
                    //Primer Objetivo
                    crowding_distance += ((*std::prev(it))->box[n-1].lb() - (*std::next(it))->box[n-1].lb())/((first->box[n-1].lb() - last->box[n-1].lb()));
                    //Segundo objetivo
                    crowding_distance += ((*std::prev(it))->box[n-2].lb() - (*std::next(it))->box[n-2].lb())/((first->box[n-2].lb() - last->box[n-2].lb()));
                }

                if(min_cd >= crowding_distance){
                	min_cd=crowding_distance;
                	erase_it = it;
                }

                //Inserts the CDBox node into the set
                //cdSet.insert(cdBox);
            }

            //cout << "nondominated size antes " << nonDominated.size() << endl;
            globalBuffer.push(*erase_it);
            nonDominated.erase(erase_it);
			//cout << "nondominated size despues " << nonDominated.size() << endl;
			//getchar();
            //std::cout << i << ", " << cdBox->crowding_distance << ", "
              //      << cdBox->C->box[cdBox->C->box.size()-1].lb() << "," << cdBox->C->box[cdBox->C->box.size()-2].lb()<< endl;
            //i++;

            returnSize--; //If the returnSize == 0, then the cdSet isn't cleared.
        }
		//If we have stuff in the nonDominated, it means that we need to put them into the nextBuffer!
		//Since the lowest crowding distance ones are already in the globalBuffer, we won't need to remove them from here. We will just have to add stuff into the nextBuffer.
		while(!nonDominated.empty()){
			currentBuffer.push(*nonDominated.begin());
			nonDominated.erase(nonDominated.begin());
		}


	}

	bool CrowdingDistanceBSMOP::dominates(Cell* a, Cell* b){
        int n = a->box.size();
		int m = b->box.size();
        return (a->box[n-1].lb() <= b->box[m-1].lb() && a->box[n-2].lb() <= b->box[m-2].lb());
    }

	bool CrowdingDistanceBSMOP::doOverlap(Point l1, Point r1, Point l2, Point r2){ 
    	// If one rectangle is on left side of other 
		if (l1.x > r2.x || l2.x > r1.x) 
			return false; 
	
		// If one rectangle is above other 
		if (l1.y < r2.y || l2.y < r1.y) 
			return false; 
	
		return true; 
	} 

} // end namespace ibex
