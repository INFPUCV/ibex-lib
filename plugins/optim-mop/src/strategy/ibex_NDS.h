/*
 * ibex_NDS.h
 *
 *  Created on: 24 may. 2018
 *      Author: iaraya
 */

#include "ibex_IntervalVector.h"
#include "ibex_CellMOP.h"
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


	void addPoint(pair< double, double> eval);

	void addSegment(pair< double, double> eval1, pair< double, double> eval2) {

		if(eval1.first == eval2.first  &&  eval1.second == eval2.second ){
			addPoint(eval1);
			return;
		}



		cout << "addVectortoNDS-------" << endl;
		std::map<pair<double, double>, IntervalVector>::iterator aux, it1 = --NDS2.lower_bound(eval1);
		pair< double, double> first, second, point, inter_last = make_pair(NEG_INFINITY, NEG_INFINITY);
		std::map< pair <double, double>, IntervalVector, sorty2 > DS2;
		std::map< pair <double, double>, IntervalVector, sorty2 > NDS_points;

		DS2.insert(*it1);
		it1++;

		bool flagDS2 = false;
		for(;it1 != NDS2.end();) {
			if(it1->first.second < eval1.second and it1->first.second < eval2.second) flagDS2= true;
			aux = it1;
			++aux;
			DS2.insert(*it1);
			NDS2.erase(it1);
			it1 = aux;
			if(flagDS2) break;
		}

		cout << "DS2" << endl;
		for(it1 = DS2.begin();it1!=DS2.end();++it1) {
			cout << "(" << it1->first.first << "," << it1->first.second << ") " << endl;
		}
		cout << "end DS2" << endl;

		it1 = DS2.begin();
		first = it1->first;
		it1++;

		// false cuando la recta pasa por fuera.
		/*
		bool flag = false;

		for(;it1 != DS2.end(); ++it1) {
			second = it1->first;
			cout << "(" << first.first << "," << first.second << ") ";
			cout << "(" << second.first << "," << second.second << ")" << endl;
			point = pointIntersection(first, second, eval1, eval2);

			cout << "intersection (" << point.first << "," << point.second << ")" << endl;

			// puntos muy cercanos son los mismos
			if(fabs(inter_last.first - point.first) > 1e-7 && fabs(inter_last.second - point.second) > 1e-7
					&& first.first <= point.first + 1e-4 && point.first - 1e-4 <= second.first  // por x
					&& second.second <= point.second + 1e-4 && point.second - 1e-4 <= first.second) { // por y
				cout << "intersection (" << point.first << "," << point.second << ")" << endl;
				inter_last = point;

				NDS_points.insert(make_pair(point,it1->second));

				if(flag) flag = false;
				else flag = true;
			}else{
				// si pasa por fuera se agregan los punto
				if(!flag) {
					NDS_points.insert(make_pair(second,it1->second));
				}
			}

			first = second;
		}
		*/
		bool flag = false;
		double epsilonRaro = 1e-7;
		for(;it1 != DS2.end(); ++it1) {
			second = it1->first;
			cout << "(" << first.first << "," << first.second << ") ";
			cout << "(" << second.first << "," << second.second << ")" << endl;
			point = pointIntersection(first, second, eval1, eval2);

			cout << "intersection (" << point.first << "," << point.second << ")" << endl;

			// puntos muy cercanos son los mismos
			if(fabs(inter_last.first - point.first) > epsilonRaro && fabs(inter_last.second - point.second) > epsilonRaro) {

				if( ((first.first <= second.first && first.first - epsilonRaro <= point.first && point.first <= second.first + epsilonRaro) ||
						(second.first <= first.first  && second.first - epsilonRaro <= point.first && point.first <= first.first + epsilonRaro)) &&
						((first.second <= second.second && first.second - epsilonRaro <= point.second && point.second <= second.second + epsilonRaro) ||
						(second.second <= first.second  && second.second - epsilonRaro <= point.second && point.second <= first.second + epsilonRaro)) &&
						((eval1.first <= eval2.first && eval1.first - epsilonRaro <= point.first && point.first <= eval2.first + epsilonRaro) ||
						(eval2.first <= eval1.first  && eval2.first - epsilonRaro <= point.first && point.first <= eval1.first + epsilonRaro)) &&
						((eval1.second <= eval2.second && eval1.second - epsilonRaro <= point.second && point.second <= eval2.second + epsilonRaro) ||
						(eval2.second <= eval1.second  && eval2.second - epsilonRaro <= point.second && point.second <= eval1.second + epsilonRaro))) {

					cout << (first.first <= second.first && first.first <= point.first && point.first <= second.first) <<
							" " << (second.first <= first.first  && second.first <= point.first && point.first <= first.first) <<
							" " << (first.second <= second.second && first.second <= point.second && point.second <= second.second) <<
							" " << (second.second <= first.second  && second.second <= point.second && point.second <= first.second) << endl;

					cout << "intersection (" << point.first << "," << point.second << ")" << endl;
					inter_last = point;

					NDS_points.insert(make_pair(point,it1->second));
				} else {
					cout << "intersection error" << endl;
				}

				//if(flag) flag = false;
				//else flag = true;

			//}else{
				// si pasa por fuera se agregan los punto
				//if(!flag) {
				//	NDS_points.insert(make_pair(second,it1->second));
				//}
			} else {
				cout << "point very close to the after point" << endl;
			}

			first = second;
		}

		cout << "------------" << endl;

		for(it1 = NDS_points.begin();it1 != NDS_points.end();++it1) {
			cout << "newpoint (" << it1->first.first << "," << it1->first.second << ")" << endl;

			if(it1->first.first != it1->first.first  || it1->first.second != it1->first.second ) getchar();
			NDS2.insert(*it1);
		}

		cout << "insert points" << endl;

		for(it1 = DS2.begin();it1!=DS2.end();++it1) {
			if( ((it1->first.first <= eval1.first && it1->first.first <= eval2.first) &&
					(it1->first.second >= eval1.second && it1->first.second >= eval2.second)) ||
				((it1->first.second <= eval1.second && it1->first.second <= eval2.second) &&
						(it1->first.first >= eval1.first && it1->first.first >= eval2.first)) ) {
				cout << "fuera" << endl;
				cout << "(" << it1->first.first << "," << it1->first.second << ") " << endl;
				NDS2.insert(*it1);
				// addPointtoNDS(make_pair(it1->first.first,it1->first.second));
			} else {
				double m = (eval1.second - eval2.second)/(eval1.first - eval2.first);
				if( (eval1.second - m*eval1.first) >= (it1->first.second - m*it1->first.first)) {
					cout << "dentro" << endl;
					cout << "(" << it1->first.first << "," << it1->first.second << ") " << endl;
					NDS2.insert(*it1);
				} else {
					cout << "eliminado" << endl;
					cout << "(" << it1->first.first << "," << it1->first.second << ") " << endl;
				}
			}
		}

	}

    /**
     * Obtiene la pendiente del segmento
     * @param first  - Vector inicial del segmento
     * @return              - Rectorna la pendiente del segmento
     */
    static double getSlopeSegment(pair<double, double> first, pair<double, double> last){
        double slope1, slope2, Ffirst, Fsecond, Lfirst, Lsecond;

        if((last.second == NEG_INFINITY and first.second == NEG_INFINITY) or
        		(last.second == POS_INFINITY and first.second == POS_INFINITY))
        	return 0;
        if((last.first == NEG_INFINITY and first.first == NEG_INFINITY) or
        		(last.first == POS_INFINITY and first.first == POS_INFINITY))
        	return POS_INFINITY;
        //cout << (last.first == first.first) << endl;
        if(last.second == first.second) return 0;
        if(last.first == first.first) return POS_INFINITY;
        slope1 = last.second - first.second;
        //cout << slope1 << endl;
        slope2 = last.first - first.first;
        //cout << slope2 << endl;
        return slope1/slope2;

    }


	//******************************************************************************************************************
	//Intersecciones y puntos
	/**
	 * @param v10 punto inicial del primer segmento
	 * @param v11 punto final del primer segmento
	 * @param v20 punto inicial del segnundo segmento
	 * @param v21 punto final del segundo segmento
	 * @return Punto o vector de interseccion entre los segmentos evaluados,
	 *         retorna null si no hay interseccion entre los segmentos
	 */
	static pair<double, double> pointIntersection(
			pair<double, double> v10, pair<double, double> v11,
			pair<double, double> v20, pair<double, double> v21){
		pair<double, double> interseccion = make_pair(0.0 ,0.0);

		double m = getSlopeSegment(v10 , v11);
		//cout << "m " << m << endl;
		double n = getSlopeSegment(v20, v21);
		double c = v10.second - v10.first * m;
		double d = v20.second - v20.first * n;

		/*
		cout << "v10 (" << v10.first << "," << v10.second << ")" << endl;
		cout << "v11 (" << v11.first << "," << v11.second << ")" << endl;
		cout << "v20 (" << v20.first << "," << v20.second << ")" << endl;
		cout << "v21 (" << v21.first << "," << v21.second << ")" << endl;
		cout << v10.second << " " << v10.first << " " << m << endl;
		cout << "c " << c << " d " << d << endl;
		*/

		//Cuando el segmento es vertical su pendiente es infinito
		if(m == POS_INFINITY and n == 0) {
			//cout << "1111" << endl;
			if((v20.first <= v10.first and v10.first <= v21.first) or
					(v21.first <= v10.first and v10.first <= v20.first))
				interseccion.first = v10.first;
			if((v10.second <= v20.second and v20.second <= v11.second) or
					(v11.second <= v20.second and v20.second <= v10.second))
				interseccion.second = v20.second;
		}
		else if (m == 0 and n == POS_INFINITY){
			//cout << "22222" << endl;
			if((v10.first <= v20.first and v20.first <= v11.first) or
					(v11.first <= v20.first and v20.first <= v10.first))
				interseccion.first = v20.first;
			if((v20.second <= v10.second and v10.second <= v21.second) or
					(v21.second <= v10.second and v10.second <= v20.second))
				interseccion.second = v10.second;
		}
		else if (v11.first - v10.first == 0){
			//cout << "33333" << endl;
			interseccion.first = v11.first;
			interseccion.second = n*v11.first + d;
		}
		//Cuando el segmento es vertical su pendiente es infinito
		else if (v21.first - v20.first == 0){
			//cout << "444" << endl;
			//cout << "m " << m << " v21.first " << v21.first << endl;
			interseccion.first = v21.first;
			interseccion.second = m*v21.first + c;
		}
		else{
			//cout << "5555" << endl;
			interseccion.first = (d - c) / (m - n);
			interseccion.second = m*interseccion.first + c;
		}
		return interseccion;
	}


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
	    firstp = pointIntersection( v10, v11, make_pair(box[n].lb(),v11.second),  make_pair(box[n].lb(),v10.second));
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


		lastp = pointIntersection( v20, v21, make_pair(v20.first,box[n+1].lb()),  make_pair(v21.first,box[n+1].lb()));

		inpoints.push_back(lastp);

		return inpoints;
	}


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
	    firstp = pointIntersection( v10, v11, make_pair(box[n].lb(),v11.second),  make_pair(box[n].lb(),v10.second));
	   	if(firstp.second <= box[n+1].lb() ){
	   		box.set_empty();
	   		return inpoints;
	   	}

		if(firstp.second > box[n+1].ub())
			firstp = pointIntersection( v10, v11, make_pair(v10.first,box[n+1].ub()),  make_pair(v11.first,box[n+1].ub()));

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


		lastp = pointIntersection( v20, v21, make_pair(v20.first,box[n+1].lb()),  make_pair(v21.first,box[n+1].lb()));

		if(lastp.first > box[n].ub() )
			lastp = pointIntersection( v20, v21, make_pair(box[n].ub(),v21.second),  make_pair(box[n].ub(),v20.second));

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
