#ifndef CROWDING_DISTANCE_NYU_H_
#define CROWDING_DISTANCE_NYU_H_

#include "ibex_CellMOP.h"
#include "ibex_CellSet.h"
#include "ibex_NDS.h"
#include "ibex_BeamSearchBufferMOP.h"
#include <queue>
#include <map>
#include <set>
#include <iterator>
#include <algorithm>
#include <limits>



namespace ibex {

	//Struct that contains Cell and crowding distance
	typedef struct CDBox{
		Cell* C;
		double crowding_distance;
	} CDBox;



    //Ordena por objetivo y1, y2
    struct crowding_distanceBeam {
        bool operator() (const Cell* c1, const Cell* c2){
            int n = c1->box.size();
            if(c1->box[n-2].lb() <= c2->box[n-2].lb() && c1->box[n-1].lb() <= c2->box[n-1].lb()) return true;
            else return false;
        }
    };

    struct sortByCrowdingDistance{
        bool operator()(const CDBox* c1, const CDBox* c2){
            if(c1->crowding_distance <= c2->crowding_distance) return true;
            else return false;
        }

    };

    class NyuCrowdingDistance : public CellBufferOptim{
        public:
    	static std::multiset<Cell*,max_distanceBeam> getCrowdingDistance(std::multiset<Cell*, max_distanceBeam>& nextBuffer, int currentBuffer);

    };



}

#endif
