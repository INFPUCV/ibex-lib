/*      Crowding Distance Implementation
 *                  v0.0
 *  
 *              Author: Nyuku
 */

#include <iostream>
#include <set>
#include <iterator>
#include <algorithm>
#include <limits>
#include "ibex_CellMOP.h"
#include "ibex_CellSet.h"
#include "ibex_BeamSearchBufferMOP.h"

#define INF std::numeric_limits<double>::infinity();

namespace ibex{

    typedef struct CDBox{
        Cell* C;
        int rank;
        double crowding_distance;
    } CDBox;

    struct crowding_distanceBeam {
        bool operator() (const Cell* c1, const Cell* c2){
            int n = c1->box.size();
            if(c1->box[n-2].lb() <= c2->box[n-2].lb() && c1->box[n-1].lb() <= c2->box[n-1].lb()) return true;
            else return false;
        }
    };

    struct sortByCrowdingDistance{
        bool operator()(const CDBox* c1, const CDBox* c2){
            if(c1.crowding_distance <= c2.crowding_distance) return true;
            else return false;
        }

    }


    void getCrowdingDistance(std::multiset nextBuffer){
        std::multiset<Cell*, crowding_distanceBeam> cdBuffer;
        std::multiset<CDBox*, sortByCrowdingDistance> cdSet;

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
            if(it == cdBuffer.begin || it == cdBuffer.end()){
                cdBox.crowding_distance = INF;
            }
            else{ //If not, then calculate them.
                int n = it->box.size();
                Cell* first = (Cell*)cdBuffer.end();
                Cell* last = (Cell*)cdBuffer.begin();

                //TODO: where the fuck do i get the objective numbers from!?
                cdBox.crowding_distance = it->box[n-2].lb() - it->box[n-1].lb()
                cdBox.crowding_distance /= (first->box[first->box.size()].lb() - last->box[0].lb());
            }

            //Inserts the CDBox node into the set
            cdSet.insert(cdBox);
            std::cout << cdBox.crowding_distance << endl; //%temp%
        }
        getchar(); //%temp%
        //TODO: Return a set with X amount of Cells.
    }
}