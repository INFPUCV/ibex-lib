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

	double hypervolume(const Interval& y1, const Interval& y2) const{
		double hv=0.0;
		double prev1=y1.lb();
		double prev2=y2.ub();
		for(auto ndp:NDS2){
			double next1=std::min(ndp.first.first,y1.ub());
			double next2=std::min(y2.ub(),ndp.first.second);
			if(next1 > prev1){
				hv+=(next1-prev1)*(next2 + (prev2-next2)/2.0);
				//cout << prev1 << "," << prev2 << endl;
				//cout << next1 << "," << next2 << endl;
			}

			if(next1>=prev1 && next2<=prev2){
				prev1=next1;
				prev2=next2;
			}
		}
		return hv;
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

	void addPoint(IntervalVector& y){
		addPoint(pair<double,double>(y[0].ub(), y[1].ub()));
	};

	/**
	* Add a segment in the NDS structure
	*/
	void addSegment(pair< double, double> p1, pair< double, double> p2);

    struct NoIntersectionException : public exception {
       const char * what () const throw () {
          return "NoIntersectionException";
       }
    };

	//returns a list of points non-dominated by lb
	static list<pair <double,double> > non_dominated_points(double lbx, double lby){
		list <pair <double,double> > inpoints;

		pair <double, double> firstp;
		pair <double, double> lastp;

		//first potential dominated point
		map< pair <double, double>, IntervalVector >:: iterator ent1=NDS2.upper_bound(make_pair(lbx,NEG_INFINITY));
		pair <double, double> v11 = ent1->first; ent1--;

		//last point before x=lbx
	    pair <double, double> v10 = ent1->first;

		//x-point cutting lbx
	    firstp = pointIntersection( v10, v11, make_pair(lbx,v11.second),  make_pair(lbx,v10.second));
	   	if(firstp.second <= lby ){
	   		return inpoints; //empty list
	   	}

		inpoints.push_back(firstp);

		ent1++;
		//points dominated by lb
		while(ent1->first.second > lby){
			inpoints.push_back(ent1->first);
			ent1++;
		}

		//last potential dominated point
		ent1--;
		pair <double, double> v21 = ent1->first;
		//first point after y=lby
		ent1++;
		pair <double, double> v20 = ent1->first;

		//x-point cutting lby
		lastp = pointIntersection( v20, v21, make_pair(v20.first,lby),  make_pair(v21.first,lby));

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

  static bool _trace;
	static double distance(const Cell* c){
		int n=c->box.size();

		double a = c->get<CellMOP>().a;
		double w_lb = c->get<CellMOP>().w_lb;
    //cout << "a:" << a << endl;
		//cout << "w_lb:" << w_lb << endl;

    double dist;
		_trace=true;
		if(a!=0)
		  dist= distance(c->box[n-2].lb(),c->box[n-1].lb(),-1/a, w_lb/a);
    else
		   dist= distance(c->box[n-2].lb(),c->box[n-1].lb());
		_trace=false;
		return dist;
	}



	// m in [-oo, 0]
	static double distance(double lbx, double lby, double m=POS_INFINITY, double c=POS_INFINITY){
		double max_dist=NEG_INFINITY;

		Interval Ay=lby;
		Interval Bx=lbx;

		if(m!=POS_INFINITY){
			Ay = Interval(m)*lbx+c;  // Ax=lbx
			Bx = (Interval(lby)-c)/m; // By=lby
			if(Ay.lb() < lby){ Ay=lby; }
			if(Bx.lb() < lbx){ Bx=lbx; }
		}


		//	cout << "A:" << lbx << "," << Ay.mid() << endl;
		//	cout << "B:" << Bx.mid() << "," << lby << endl;

		list<pair <double,double> > inner_segments= non_dominated_points(lbx, lby);
		pair <double,double>* p0=NULL;

		bool Adist=false;
		bool Bdist=false;

		//cout << "dist inner-size:" << inner_segments.size() << endl;

		for(auto p : inner_segments){
			if(p.first==POS_INFINITY && p.second==POS_INFINITY) return POS_INFINITY;

			Interval dist;
			//up-left point
			if(p.first-lbx < (p.second-Ay).ub() || p.second==POS_INFINITY){
				dist=p.first-Interval(lbx);
				//cout << "dist-UL:" << p.first  << "," << p.second << "-->" <<  dist.ub() << endl;
			}
			//bottom-right point
			else if(p.second-lby < (p.first-Bx).ub() || p.first==POS_INFINITY){
				dist=p.second-Interval(lby);
				//cout << "dist-BR:" << p.first  << "," << p.second << "-->" <<  dist.ub() << endl;
				if(!Bdist && p0){
					Interval mm=NEG_INFINITY;
					if(p.first-p0->first != 0)
						mm= (Interval(p.second)-p0->second)/(Interval(p.first)-p0->first);

						//cout << p.first-p0->first << endl;
					//cout << "dist-B (m):" << mm << endl;
					if(mm.lb() > -1 && mm.lb() < 0.0 ){
						Interval cc= p.second - mm*p.first;
						//cout << "dist-B:" << (mm*lbx - Ay + cc)/(1.0-mm) << endl;
						dist=std::max(dist.ub(), ((mm*lbx - Ay + cc)/(1.0-mm)).lb());
						Bdist=true;
					}
				}
			}
			//cy-45-degree zone
			else{
				dist= -(m*p.first - p.second+c)/(1.0-m);
			//	cout << "A:" << lbx << "," << Ay.mid() << endl;
			//	cout << "B:" << Bx.mid() << "," << lby << endl;
				//cout << "dist-IN:" << p.first  << "," << p.second << "-->" <<  dist.ub() << endl;
				if(!Adist && p0){
					Interval mm=NEG_INFINITY;
					if(p.first-p0->first != 0)
						mm= (Interval(p.second)-p0->second)/(Interval(p.first)-p0->first);

					if(mm.lb() <-1 && mm.lb()>NEG_INFINITY){
						Interval cc= p.second - mm*p.first;
						//cout << "dist-A:" << (mm*Bx - lby + cc)/(1.0-mm) << endl;
						dist=std::max(dist.ub(), ((mm*Bx - lby + cc)/(1.0-mm)).lb());
						Adist=true;
					}
				}
			}


			if(dist.ub()>max_dist){

				max_dist=dist.ub();
			}

			if(p0) delete p0;
			p0=new pair <double,double>(p);

		}
		if(p0) delete p0;
		return max_dist;
	}


	/** The current non-dominated set sorted by increasing x */
	static map< pair <double, double>, IntervalVector, sorty2 > NDS2;
};

} /* namespace ibex */

#endif /* OPTIM_MOP_SRC_STRATEGY_IBEX_NDS_H_ */
