/*
 * ibex_NDS.h
 *
 *  Created on: 24 may. 2018
 *      Author: iaraya
 */

#include "ibex_IntervalVector.h"
#include <map>
#include <list>
#include <unordered_set>
#include "ibex_BxpMOPData.h"
#include "ibex_OptimizerMOP.h"

#ifndef OPTIM_MOP_SRC_STRATEGY_IBEX_NDSRP_H_
#define OPTIM_MOP_SRC_STRATEGY_IBEX_NDSRP_H_

using namespace std;

namespace ibex {

class Point : public Vector{
    public:
        double hv = 0.0;
        Point* next=NULL;
        Point* prev=NULL;


        void push_before(Point* next){
            if(next->prev) next->prev->next=this;
            this->prev=next->prev;
            next->prev=this;
            this->next=next;
        }

        Point(double a, double b, Point* next=NULL) : Vector(2) {
            (*this)[0]=a; (*this)[1]=b;
            if(next) push_before(next);
        }

        //create and add the point before next
        Point(const Vector &v, Point* next=NULL) : Vector(v){
            if(next) push_before(next);
            //update hv
        }

        //remove the point
        ~Point(){
            if(prev){
                prev->next->prev=prev;
                prev->next=next;
            }
        }

};

/**
 * comparation function for sorting NDS by increasing x and decreasing by y
 */
struct sort_rp{
	bool operator()(const Vector* y1, const Vector* y2){
		if((*y1)[0] != (*y2)[0])
			return (*y1)[0]<(*y2)[0];
		return (*y1)[1]>(*y2)[1];

	}
};

/**
 * \brief Segment based non-dominated set
 */
class NDSrp {
public:

	virtual void clear();

	virtual void NDS_clear(){
        for(auto p:NDS)
            delete p;
        
		NDS.clear();
	}

	virtual void NDS_insert(const Vector& p){
        set<Vector*>::iterator next = NDS.lower_bound((Vector*)&p);
        Point* pp = new Point(p, (Point*)(*next));

        //set<Point>::iterator curr = NDS.insert(Point(p)).first;
        //next->prev->next=&(*curr);
        //curr->prev=next->prev;
        //next->prev=&(*curr);
        //curr->next=&(*next);
	}

	virtual void NDS_erase(std::set<Vector*>::iterator it){
		Vector* p= *it;
		NDS.erase(it);
        delete p;
	}

    inline int size() const{
		return NDS.size();
	}

    /**
	 * \brief return true if new_p is dominated by the NDS
	 */
	bool is_dominated(const Vector& new_p);

	/**
	* Add a point in the NDS structure
	*/
	void addPoint(const Vector& y);

	/**
	* Add a segment in the NDS structure. Return 0 if the segment did not modify the NDS
	*/
	bool addSegment(const pair< Vector, Vector>& y1y2);

    //returns a list of points non-dominated by lb
	list< Vector > non_dominated_points(const Vector& lb, bool cutting_points=true);

    /**
        * \brief Returns the point intersecting two segments. Otherwise it throw
        * a NoIntersectionException
        * It is conservative, that is:
        * 1) if there are intersection it should return a point dominated by the real intersection
        * 2) if there are no intersection it may return an exception or a point dominated by one segment
        * TODO: revise with colinear generated examples
	*/
	static Vector pointIntersection(const Vector& p0, const Vector& p1, const Vector& p2, const Vector& p3);

    //Distance from a cell (cy_box) to the NDS
    double distance(const Cell* c);

	// Distance from segment (yA,yB) to the NDS
	double distance(Vector& yA, Vector& yB, double m=POS_INFINITY, double c=POS_INFINITY);

    //returns the segment yA--yB of the line y_2=m*y_2+c dominated by lb
    static pair <Vector, Vector> get_segment(const Vector& lb, double m=POS_INFINITY, double c=POS_INFINITY);

	// The current non-dominated set sorted by increasing y1
	set< Vector*, sort_rp > NDS;

    struct NoIntersectionException : public exception {
       const char * what () const throw () {
          return "NoIntersectionException";
        }
    };
};


} /* namespace ibex */

#endif /* OPTIM_MOP_SRC_STRATEGY_IBEX_NDS_H_ */
