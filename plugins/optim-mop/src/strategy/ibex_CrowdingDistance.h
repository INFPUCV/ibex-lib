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

#define INF std::numeric_limits<double>::infinity();

namespace ibex {

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
            if(c1.crowding_distance <= c2.crowding_distance) return true;
            else return false;
        }

    };

    class NyuCrowdingDistance : public CellBufferOptim{
        public:
            static int currentBuffer;
            mutable std::multiset<Cell*, crowding_distanceBeam> cdBuffer;
            mutable std::multiset<CDBox*, sortByCrowdingDistance> cdSet;

            std::multiset getCrowdingDistance(std::multiset nextBuffer);

    }



}
