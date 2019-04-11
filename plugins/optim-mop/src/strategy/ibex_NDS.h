/*
 * ibex_NDS.h
 *
 *  Created on: 24 may. 2018
 *      Author: iaraya
 */

#include "ibex_IntervalVector.h"
#include "ibex_pyPlotter.h"
#include <map>
#include <list>
#include "ibex_BxpMOPData.h"

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

	pair<double,double> lb(){
		if(NDS2.size()>3){
			auto it=NDS2.begin();
			it++;
			double min1=it->first.first;
			it=NDS2.end();
			it--;
			it--;
			double min2=it->first.second;
			return make_pair(min1,min2);
		}
		else return make_pair(POS_INFINITY, POS_INFINITY);
	}

	pair<double,double> nadir(){
		if(NDS2.size()>3){
			auto it=NDS2.begin();
			it++;it++;
			double min2=it->first.second;
			it=NDS2.end();
			it--;
			it--;it--;
			double min1=it->first.first;
			return make_pair(min1,min2);
		}
		else return make_pair(NEG_INFINITY, NEG_INFINITY);
	}

	Interval hypervolume(const Interval& y1, const Interval& y2) const{
		Interval hv=0.0;
		double prev1=y1.lb();
		double prev2=y2.ub();
		for(auto ndp:NDS2){
			if(ndp.first.first == NEG_INFINITY) continue;
			if(ndp.first.second == NEG_INFINITY) continue;

			double next1=std::min(ndp.first.first,y1.ub());
			double next2=std::min(y2.ub(),ndp.first.second);

			if(next1 > prev1){
				Interval hvv=(Interval(next1)-prev1)*( (y2.ub()	-Interval(prev2))  +  (Interval(prev2)-next2)/2.0);
				hv+=hvv;

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
	* Add a segment in the NDS structure. Return 0 if the segment did not modify the NDS
	*/
	bool addSegment(pair< double, double> p1, pair< double, double> p2);

    struct NoIntersectionException : public exception {
       const char * what () const throw () {
          return "NoIntersectionException";
       }
    };

	//returns a list of points non-dominated by lb
	list<pair <double,double> > non_dominated_points(double lbx, double lby){
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

  void remove_infinity(pair<double, double>& p);

	void add_infinity(pair<double, double>& p);

  /**
	* \brief Returns the point intersecting two segments. Otherwise it throw
	* a NoIntersectionException
	* It is conservative, that is:
	* 1) if there are intersection it should return a point dominated by the real intersection
	* 2) if there are no intersection it may return an exception or a point dominated by one segment
	* TODO: revise with colinear generated examples
	*/
	pair<double, double> pointIntersection(pair<double, double> p0, pair<double, double> p1,
			pair<double, double> p2, pair<double, double> p3);

  static bool _trace;
	double distance(const Cell* c){
		int n=c->box.size();


		double a = ((BxpMOPData*) c->prop[BxpMOPData::id])->a;
		double w_lb = ((BxpMOPData*) c->prop[BxpMOPData::id])->w_lb;


		double dist;
		IntervalVector points(4);
		if(a!=0){
			points = get_points(c->box[n-2].lb(),c->box[n-1].lb(),-1/a, w_lb/a);
			dist= distance(points,-1/a, w_lb/a);
		}else{
			points = get_points(c->box[n-2].lb(),c->box[n-1].lb());
			dist= distance(points);
		}

		return std::max(0.0,dist);
	}

  IntervalVector get_points(double lbx, double lby, double m=POS_INFINITY, double c=POS_INFINITY){
		double max_dist=NEG_INFINITY;

		Interval Ay=lby;
		Interval Bx=lbx;

		if(m!=POS_INFINITY){
			Ay = Interval(m)*lbx+c;  // Ax=lbx
			Bx = (Interval(lby)-c)/m; // By=lby
			if(Ay.lb() < lby){ Ay=lby; }
			if(Bx.lb() < lbx){ Bx=lbx; }
		}

		IntervalVector points(4);
		points[0]=lbx;
		points[1]=Ay;
		points[2]=Bx;
		points[3]=lby;
		return points;
	}


	// m in [-oo, 0]
	double distance(IntervalVector& points, double m=POS_INFINITY, double c=POS_INFINITY){
		Interval Ax=points[0];
		Interval Ay=points[1];
		Interval Bx=points[2];
		Interval By=points[3];

		double max_dist=NEG_INFINITY;

		//	cout << "A:" << lbx << "," << Ay.mid() << endl;
		//	cout << "B:" << Bx.mid() << "," << By << endl;

		list<pair <double,double> > inner_segments= non_dominated_points(Ax.mid(), By.mid());
		pair <double,double>* p0=NULL;

		bool Adist=false;
		bool Bdist=false;

		//cout << "dist inner-size:" << inner_segments.size() << endl;

		for(auto p : inner_segments){
			if(p.first==POS_INFINITY && p.second==POS_INFINITY) return POS_INFINITY;

			Interval dist;
			//up-left point
			if((p.first-Ax).lb() < (p.second-Ay).ub() || p.second==POS_INFINITY){
				dist=p.first-Interval(Ax);
				//cout << "dist-UL:" << p.first  << "," << p.second << "-->" <<  dist.ub() << endl;
			}
			//bottom-right point
			else if((p.second-By).lb() < (p.first-Bx).ub() || p.first==POS_INFINITY){
				dist=p.second-Interval(By);
				//cout << "dist-BR:" << p.first  << "," << p.second << "-->" <<  dist.ub() << endl;
				if(!Bdist && p0){
					Interval mm=NEG_INFINITY;
					if(p.first-p0->first != 0)
						mm= (Interval(p.second)-p0->second)/(Interval(p.first)-p0->first);

						//cout << p.first-p0->first << endl;
					//cout << "dist-B (m):" << mm << endl;
					if(mm.lb() > -1 && mm.lb() < 0.0 ){
						Interval cc= p.second - mm*p.first;
						//cout << "dist-B:" << (mm*Ax - Ay + cc)/(1.0-mm) << endl;
						dist=std::max(dist.ub(), ((mm*Ax - Ay + cc)/(1.0-mm)).lb());
						Bdist=true;
					}
				}
			}
			//cy-45-degree zone
			else{
				dist= -(m*p.first - p.second+c)/(1.0-m);
			//	cout << "A:" << Ax << "," << Ay.mid() << endl;
			//	cout << "B:" << Bx.mid() << "," << By << endl;
				//cout << "dist-IN:" << p.first  << "," << p.second << "-->" <<  dist.ub() << endl;
				if(!Adist && p0){
					Interval mm=NEG_INFINITY;
					if(p.first-p0->first != 0)
						mm= (Interval(p.second)-p0->second)/(Interval(p.first)-p0->first);

					if(mm.lb() <-1 && mm.lb()>NEG_INFINITY){
						Interval cc= p.second - mm*p.first;
						//cout << "dist-A:" << (mm*Bx - By + cc)/(1.0-mm) << endl;
						dist=std::max(dist.ub(), ((mm*Bx - By + cc)/(1.0-mm)).lb());
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
	map< pair <double, double>, IntervalVector, sorty2 > NDS2;
};

} /* namespace ibex */

#endif /* OPTIM_MOP_SRC_STRATEGY_IBEX_NDS_H_ */
