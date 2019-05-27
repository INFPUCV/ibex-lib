/*      Crowding Distance Implementation
 *                  v0.0
 *  
 *              Author: Nyuku
 */

#include <iostream>
#include "ibex_CrowdingDistance.h"

using namespace std;

namespace ibex{


    std::multiset<Cell*, max_distanceBeam> NyuCrowdingDistance::getCrowdingDistance(std::multiset<Cell*, max_distanceBeam>& nextBuffer, int currentBuffer){

        std::multiset<Cell*, crowding_distanceBeam> cdBuffer;
        std::multiset<CDBox*, sortByCrowdingDistance> cdSet;

        //Leaving these commented out, we are now working as an object!!!
        //std::multiset<Cell*, crowding_distanceBeam> cdBuffer;
        //std::multiset<CDBox*, sortByCrowdingDistance> cdSet;
        std::multiset<Cell*, max_distanceBeam> ret; //Returning set

        //First, we insert nextBuffer items into cdBuffer, which we will use to insert into
        //cdSet, which will be the one we will use to insert into the current Buffer.
        //AKA: Sort(Solution,Objective)
        //TODO: Optimize this, use less structures
        for(std::multiset<Cell*, max_distanceBeam>::iterator it = nextBuffer.begin(); 
        it != nextBuffer.end(); it++)
        {
            cdBuffer.insert(*it);
        }

        for(std::multiset<Cell*, crowding_distanceBeam>::iterator it = cdBuffer.begin();
        it != cdBuffer.end(); it++)
        {

            CDBox* cdBox= (CDBox*)malloc(sizeof(CDBox));
            cdBox->C = *it;

            //If they are the first/last one, set their Crowding Distances to Infinite
            if(*it == *cdBuffer.begin() || *next(it) == *cdBuffer.end()){
                std::cout << "SAME" << endl;
                cdBox->crowding_distance = std::numeric_limits<double>::infinity();
            }
            else{ //If not, then calculate them.
                int n = (*it)->box.size();
                Cell* first = *cdBuffer.begin();
                Cell* last = *prev(cdBuffer.end());


                /* Here we calculate the crowding distance. We won't use a for, since we already know that there's only 2 objectives (y1, y2)*/
                //Primer Objetivo
                cdBox->crowding_distance += ((*std::prev(it))->box[n-1].lb() - (*std::next(it))->box[n-1].lb())/((first->box[n-1].lb() - last->box[n-1].lb()));
                //Segundo objetivo
                cdBox->crowding_distance += ((*std::prev(it))->box[n-2].lb() - (*std::next(it))->box[n-2].lb())/((first->box[n-2].lb() - last->box[n-2].lb()));
            }

            //Inserts the CDBox node into the set
            cdSet.insert(cdBox);
        }
        
        //Returns a set with a size of -currentBuffer- Cells. The first -currentBuffer- cells in cdBuffer get returned (First and last + currentBuffer-2 other ones)
        int i=0;
        for(auto cdbox: cdSet){
        	ret.insert(cdbox->C);
            nextBuffer.erase(cdbox->C); //Elimina la cell del nextBuffer, ya que solo interesa que este en el return!
            std::cout << i << ": " << cdbox->crowding_distance << " and " << cdbox->C->box[cdbox->C->box.size()] << endl;
        	i++;
        	if(i>=currentBuffer) break;
        }
        getchar();
        return ret;
    }
}
