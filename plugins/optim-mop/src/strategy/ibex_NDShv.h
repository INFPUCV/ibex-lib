/*
 * ibex_NDS.h
 *
 *  Created on: 24 may. 2018
 *      Author: iaraya
 */

#include "ibex_IntervalVector.h"
#include <map>
#include <list>
#include <set>
#include <unordered_set>
#include "ibex_BxpMOPData.h"


#ifndef OPTIM_POINT_H_
#define OPTIM_POINT_H_
//#include "ibex_NDSrp.h"

using namespace std;

namespace ibex {

class Point : public Vector{
    public:
        double hv_contribution = POS_INFINITY;
        //static Vector ref;
        bool up_point;
        

        Point* next=NULL;
        Point* prev=NULL;

        void push_before(Point* next){
            if(next->prev) next->prev->next=this;
            this->prev=next->prev;
            next->prev=this;
            this->next=next;
        }

        Point(double a, double b) : Vector(2) {
            (*this)[0]=a; (*this)[1]=b;
        }

        //create and add the point before next
        Point(const Vector &v) : Vector(v){
        }

        //remove the point
        ~Point(){}

        static double compute_m(const Vector &p1, const Vector &p2);

        /*
        set the variable up_point which indicates if the point is located
        up or down de line segment connecting the previous and next point
        */
        void compute_location();

        /*
        compute how much the hypervolume is reduced if the point is removed
        return a new position for the next vector in case of removal
        */
        double compute_hv_contribution() {
            Vector tmp_next(1);
            double hv = _compute_hv_contribution(tmp_next);
            if (isnan(hv))
                hv_contribution = POS_INFINITY;
            else hv_contribution = hv;
            return hv_contribution;
        }
        
        double _compute_hv_contribution(Vector& tmp_next);
        static double get_hv_contribution(const Vector& p, const Vector& prev, const Vector& next);

        //Area of a triangle given its vertices
        static double compute_area(const Vector& x, const Vector& y, const Vector& z);

        bool is_upper() const;
        static bool is_upper(const Vector& y, const Vector&x, const Vector& z);

};

struct sort_hv{

    bool operator()(const Point* y1, const Point* y2){
		if( y1->hv_contribution == y2->hv_contribution)
			return y1 < y2;

        if(y2->hv_contribution == POS_INFINITY) return true;
        if(y1->hv_contribution == POS_INFINITY) return false;


		return (y1->hv_contribution < y2->hv_contribution);
    }
};

struct sort_rp2{
	bool operator()(const Vector* y1, const Vector* y2){
		if((*y1)[0] != (*y2)[0])
			return (*y1)[0]<(*y2)[0];
		return (*y1)[1]>(*y2)[1];

	}
};

class NDShv : public set<Point*, sort_hv>, public set< Vector*, sort_rp2 > {
public:
    NDShv();
    virtual ~NDShv() {}

    virtual void insert(const Vector& p);

    virtual void erase(Vector* p);
    virtual void erase(Point* p);
    virtual void erase(std::set<Vector*>::iterator it);
    void force_erase(Point* p);

    Point* front();
    double pop_front();

    void update(Point* p);
    void update(Point* p, Vector& v);
    
    int size() const{
        return set<Vector*, sort_rp2>::size();
    }

    inline set<Vector*>::iterator lower_bound(Vector* v){
        return set< Vector*, sort_rp2 >::lower_bound(v);
    }

    inline set<Vector*>::iterator upper_bound(Vector* v){
        return set< Vector*, sort_rp2 >::upper_bound(v);
    }

    inline set<Vector*>::iterator begin(){
        return set< Vector*, sort_rp2 >::begin();
    }

    inline set<Vector*>::iterator end(){
        return set< Vector*, sort_rp2 >::end();
    }
    int insertions;
    static double min_hv;
};

} /* namespace ibex */

#endif /* OPTIM_POINT_H_ */