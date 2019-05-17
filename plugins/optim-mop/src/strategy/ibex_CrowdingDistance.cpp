/*      Crowding Distance Implementation
 *                  v0.0
 *  
 *              Author: Nyuku
 */

#include <iostream>
#include "ibex_CrowdingDistance.h"


namespace ibex{

    //Struct that contains Cell and crowding distance
    typedef struct CDBox{
        Cell* C;
        double crowding_distance;
    } CDBox;

    int NyuCrowdingDistance::currentBuffer = 4;

    std::multiset<Cell*, max_distanceBeam> NyuCrowdingDistance::getCrowdingDistance(std::multiset nextBuffer){
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
            NyuCrowdingDistance::cdBuffer.insert(*it);
        }

        for(std::multiset<Cell*, crowding_distanceBeam>::iterator it = NyuCrowdingDistance::cdBuffer.begin();
        it != NyuCrowdingDistance::cdBuffer.end(); it++)
        {

            CDBox* cdBox= (CDBox*)malloc(sizeof(CDBox));
            cdBox->C = *it;
            //If they are the first/last one, set their Crowding Distances to Infinite
            if(it == NyuCrowdingDistance::cdBuffer.begin || it == NyuCrowdingDistance::cdBuffer.end()){
                cdBox.crowding_distance = INF;
            }
            else{ //If not, then calculate them.
                int n = it->box.size();
                Cell* first = (Cell*)NyuCrowdingDistance::cdBuffer.end();
                Cell* last = (Cell*)NyuCrowdingDistance::cdBuffer.begin();

                /* Here we calculate the crowding distance. We won't use a for, since we already know that there's only 2 objectives (y1, y2)*/
                //Primer Objetivo
                cdBox.crowding_distance += ((CDBox*)(std::prev(it))->box[n-1].lb() - ((CDBox*)(std::next(it))->box[n-1].lb()))/((first->box[first->box.size()].lb() - last->box[first->box.size()].lb()));
                //Segundo objetivo
                cdBox.crowding_distance += ((CDBox*)(std::prev(it))->box[n].lb() - ((CDBox*)(std::next(it))->box[n].lb()))/((first->box[first->box.size()].lb() - last->box[first->box.size()].lb()));
            }

            //Inserts the CDBox node into the set
            NyuCrowdingDistance::cdSet.insert(cdBox);
            std::cout << cdBox.crowding_distance << endl; //%temp%
        }
        getchar(); //%temp%
        
        //Returns a set with a size of -currentBuffer- Cells. The first -currentBuffer- cells in cdBuffer get returned (First and last + currentBuffer-2 other ones)
        std::multiset<Cell*, max_distanceBeam>::iterator it = NyuCrowdingDistance::cdBuffer.begin();
        for(i = 0; i<NyuCrowdingDistance::currentBuffer; i++){
            ret.insert(((CDBox*)(it)).C);
            it++;
        }
        return ret;
    }
}