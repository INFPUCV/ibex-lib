/*
 * ibex_NDS.h
 *
 *  Created on: 24 may. 2018
 *      Author: iaraya
 */

#include "ibex_IntervalVector.h"
#include "ibex_CellMOP.h"
#include "ibex_pyPlotter.h"
#include <map>
#include <list>

#ifndef OPTIM_MOP_SRC_STRATEGY_IBEX_NDS_H_
#define OPTIM_MOP_SRC_STRATEGY_IBEX_NDS_H_

using namespace std;

namespace ibex {


/**
 * comparation function for sorting NDS2 by increasing x and decreasing by y
 */
struct sorty2{
	bool operator()(const pair<double,double> p1, const pair<double,double> p2){
		if(p1.first != p2.first)
			return p1.first<p2.first;
		return p1.second>p2.second;

	}
};

/**
 * \brief Segment based non-dominated set
 */
class NDS_seg {
public:


	NDS_seg() {
		NDS2.clear();
	};

	virtual ~NDS_seg() { };


	void clear(){
		NDS2.clear();
		//the first point
		NDS2.insert(make_pair(make_pair(NEG_INFINITY,POS_INFINITY), Vector(1)));
		//the middle point
		NDS2.insert(make_pair(make_pair(POS_INFINITY,POS_INFINITY), Vector(1)));
		//the last point
		NDS2.insert(make_pair(make_pair(POS_INFINITY,NEG_INFINITY), Vector(1)));
	}

	int size(){
		return NDS2.size();
	}

	/**
	 * \brief return true if new_p is dominated by the NDS
	 */
	bool is_dominated(vector<double>& new_p);
	 bool is_dominated(pair< double, double> new_p);

  /**
	* Add a point in the NDS structure
	*/
	void addPoint(pair< double, double> eval);

	/**
	* Add a segment in the NDS structure
	*/
	void addSegment(pair< double, double> p1, pair< double, double> p2);

    struct NoIntersectionException : public exception {
       const char * what () const throw () {
          return "NoIntersectionException";
       }
    };

	//TODO: Juntar con non_dominated_segments
	static list<pair <double,double> > extremal_non_dominated(const IntervalVector& box){
		int n=box.size()-2;
		list <pair <double,double> > inpoints;

		pair <double, double> firstp;
		pair <double, double> lastp;

		//first point in the box (v11)
		map< pair <double, double>, IntervalVector >:: iterator ent1=NDS2.upper_bound(make_pair(box[n].lb(),NEG_INFINITY));
		while(ent1->first.second > box[n].ub()) ent1++;
		if(ent1->first.first > box[n].lb()) return inpoints; //no point in the box

		pair <double, double> v11 = ent1->first; ent1--;
	    pair <double, double> v10 = ent1->first;

		//TODO: interseccion conservativa: upperPointIntersection
	    cout << 1 << endl;
	    firstp = pointIntersection( v10, v11, make_pair(box[n].lb(),v11.second),  make_pair(box[n].lb(),v10.second));
	    cout << 2 << endl;
	   	if(firstp.second <= box[n+1].lb() ){
	   		cout << "error: the box is dominated" << endl;
	   		exit(0);
	   	}

		inpoints.push_back(firstp);

		//last point in the box (v21)
		while(ent1->first.second > box[n+1].lb() || ent1->first.first > box[n].lb())
			ent1++;

		pair <double, double> v21 = ent1->first;
		ent1++;
		pair <double, double> v20 = ent1->first;

	    cout << 3 << endl;
		lastp = pointIntersection( v20, v21, make_pair(v20.first,box[n+1].lb()),  make_pair(v21.first,box[n+1].lb()));
	    cout << 4 << endl;

		inpoints.push_back(lastp);

		return inpoints;
	}

  static void remove_infinity(pair<double, double>& p);

	static void add_infinity(pair<double, double>& p);

  /**
	* \brief Returns the point intersecting two segments. Otherwise it throw
	* a NoIntersectionException
	* It is conservative, that is:
	* 1) if there are intersection it should return a point dominated by the real intersection
	* 2) if there are no intersection it may return an exception or a point dominated by one segment
	* TODO: revise with colinear generated examples
	*/
	static pair<double, double> pointIntersection(pair<double, double> p0, pair<double, double> p1,
			pair<double, double> p2, pair<double, double> p3);

	static list<pair <double,double> > non_dominated_segments(IntervalVector& box, int n){

		list <pair <double,double> > inpoints;

		pair <double, double> firstp;
		pair <double, double> lastp;

		//first point in the box (v11)
		map< pair <double, double>, IntervalVector >:: iterator ent1=NDS2.upper_bound(make_pair(box[n].lb(),NEG_INFINITY));
		while(ent1->first.second > box[n].ub()) ent1++;
		if(ent1->first.first > box[n].lb()) return inpoints; //no point in the box

		pair <double, double> v11 = ent1->first; ent1--;
	    pair <double, double> v10 = ent1->first;

		//TODO: interseccion conservativa: upperPointIntersection
	    try{
	    	firstp = pointIntersection( v10, v11, make_pair(box[n].lb(),v11.second),  make_pair(box[n].lb(),v10.second));
	    }catch(NoIntersectionException& e){
	    	return inpoints;
	    }

	   	if(firstp.second <= box[n+1].lb() ){
	   		box.set_empty();
	   		return inpoints;
	   	}

	    cout << 1 << endl;
		if(firstp.second > box[n+1].ub())
			firstp = pointIntersection( v10, v11, make_pair(v10.first,box[n+1].ub()),  make_pair(v11.first,box[n+1].ub()));
	    cout << 2 << endl;

		//valueZ1 is a point in the bounds lb(y1) or ub(y2) of the box delimiting the NDS in the box
		inpoints.push_back(firstp);

		//last point in the box (v21)
		while(ent1->first.second > box[n+1].lb() || ent1->first.first > box[n].lb()){
			inpoints.push_back(ent1->first);
			ent1++;
		}
		pair <double, double> v21 = ent1->first;
		ent1++;
		pair <double, double> v20 = ent1->first;


		try{
			lastp = pointIntersection( v20, v21, make_pair(v20.first,box[n+1].lb()),  make_pair(v21.first,box[n+1].lb()));
			if(lastp.first > box[n].ub() )
				lastp = pointIntersection( v20, v21, make_pair(box[n].ub(),v21.second),  make_pair(box[n].ub(),v20.second));
	    }catch(NoIntersectionException& e){

	    }



		//valueZ2 is a point in the bounds ub(y1) or lb(y2) of the box delimiting the NDS in the box
		inpoints.push_back(lastp);

		return inpoints;
	}


	static double distance(const Cell* c){
		//TODO: esto deberia ser parametro
		bool cy_contract_var=false;

		double max_dist=NEG_INFINITY;

		int n=c->box.size();

		Interval z1 = c->box[n-2];
		Interval z2 = c->box[n-1];

		double a = c->get<CellMOP>().a;
		double w_lb = c->get<CellMOP>().w_lb;

		map< pair <double, double>, IntervalVector >::iterator it = NDS2.lower_bound(make_pair(z1.lb(),POS_INFINITY)); //NDS.begin();
		//it--;


		list<pair <double,double> > inner_segments= extremal_non_dominated(c->box);

		if(inner_segments.size()>=2){
			pair <double,double> p1 = inner_segments.front();
			double dist = std::min (p1.first - z1.lb(), p1.second - z2.lb());
			//Damir's distance
			if(cy_contract_var && w_lb!=POS_INFINITY)
			  dist = std::min(dist, (Interval(p1.second)-(Interval(w_lb) - (Interval(p1.first) - Interval(p1.second)))/(Interval(a)+1.0)).ub());

			if(dist > max_dist) max_dist=dist;

			p1 = inner_segments.back();
			dist = std::min (p1.first - z1.lb(), p1.second - z2.lb());
			//Damir's distance
			if(cy_contract_var && w_lb!=POS_INFINITY)
			  dist = std::min(dist, (Interval(p1.second)-(Interval(w_lb) - (Interval(p1.first) - Interval(p1.second)))/(Interval(a)+1.0)).ub());

			if(dist > max_dist) max_dist=dist;
		}

		for(;it!=NDS2.end(); it++){
			pair <double, double> pmax= it->first;


			if(pmax.first==POS_INFINITY) pmax.first=CellMOP::y1_init.ub();
			if(pmax.second==POS_INFINITY) pmax.second=CellMOP::y2_init.ub();

			//cout << "z:" << z1.lb() << "," << z2.lb() << "  -->  ";
			//cout << "pmax:" << pmax.first << "," << pmax.second << endl;

			//el punto esta dentro de la zona de interes
			if(pmax.first >= z1.lb() && pmax.second >= z2.lb()){
				double dist = std::min (pmax.first - z1.lb(), pmax.second - z2.lb());
				//Damir's distance
				if(cy_contract_var && w_lb!=POS_INFINITY)
				  dist = std::min(dist, (Interval(pmax.second)-(Interval(w_lb) - (Interval(pmax.first) - Interval(pmax.second)))/(Interval(a)+1.0)).ub());

				if(dist > max_dist) max_dist=dist;
			}else break;
		}

		cout << "z:" << z1.lb() << "," << z2.lb() << "  -->  " << max_dist << endl;
		return max_dist;
	}


	/** The current non-dominated set sorted by increasing x */
	static map< pair <double, double>, IntervalVector, sorty2 > NDS2;
};

} /* namespace ibex */

#endif /* OPTIM_MOP_SRC_STRATEGY_IBEX_NDS_H_ */
