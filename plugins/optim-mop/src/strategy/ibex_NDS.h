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

	/**
	 * \brief return true if new_p is dominated by the NDS
	 */
	bool is_dominated(vector<double>& new_p);


	void addPoint(pair< double, double> eval);

	void addSegment(pair< double, double> p1, pair< double, double> p2) {

		if(p1.first == p2.first  &&  p1.second == p2.second ){
			addPoint(p1);
			return;
		}


		//se insertan en DS2 todos los puntos que se ubican entre p1 y p2 incluyendo en anterior en y1 y el siguiente en y2
		cout << "addVectortoNDS-------" << endl;
		std::map<pair<double, double>, IntervalVector>::iterator aux, it1 = --NDS2.lower_bound(p1);
		pair< double, double> second, point, inter_last = make_pair(NEG_INFINITY, NEG_INFINITY);
		list< pair <double, double>> DS2;
		std::map< pair <double, double>, IntervalVector, sorty2 > NDS_points;

		DS2.push_back(it1->first);
		it1++;

		bool flagDS2 = false;
		for(;it1 != NDS2.end();) {
			if(it1->first.second < p1.second and it1->first.second < p2.second) flagDS2= true;
			aux = it1;
			++aux;
			DS2.push_back(it1->first);
			NDS2.erase(it1);
			it1 = aux;
			if(flagDS2) break;
		}



		bool flag = false;
		double epsilonRaro = 1e-7;

		pair< double, double> prev = DS2.front();
		for(auto second:DS2){
			if(prev==second) continue;

			cout << "(" << prev.first << "," << prev.second << ") ";
			cout << "(" << second.first << "," << second.second << ")" << endl;

			//deberia retornan algo para indicar si existe la interseccion
			try{
				point = pointIntersection(prev, second, p1, p2);
			}catch(NoIntersectionException& e) {
				continue;
			}

			//luego el punto se agrega si hay interseccion
			//NDS_points.insert(make_pair(point,IntervalVector(1)));

			cout << "intersection (" << point.first << "," << point.second << ")" << endl;


			if(fabs(inter_last.first - point.first) > epsilonRaro && fabs(inter_last.second - point.second) > epsilonRaro) {

				if( ((prev.first <= second.first && prev.first - epsilonRaro <= point.first && point.first <= second.first + epsilonRaro) ||
						(second.first <= prev.first  && second.first - epsilonRaro <= point.first && point.first <= prev.first + epsilonRaro)) &&
						((prev.second <= second.second && prev.second - epsilonRaro <= point.second && point.second <= second.second + epsilonRaro) ||
						(second.second <= prev.second  && second.second - epsilonRaro <= point.second && point.second <= prev.second + epsilonRaro)) &&
						((p1.first <= p2.first && p1.first - epsilonRaro <= point.first && point.first <= p2.first + epsilonRaro) ||
						(p2.first <= p1.first  && p2.first - epsilonRaro <= point.first && point.first <= p1.first + epsilonRaro)) &&
						((p1.second <= p2.second && p1.second - epsilonRaro <= point.second && point.second <= p2.second + epsilonRaro) ||
						(p2.second <= p1.second  && p2.second - epsilonRaro <= point.second && point.second <= p1.second + epsilonRaro))) {

					cout << (prev.first <= second.first && prev.first <= point.first && point.first <= second.first) <<
							" " << (second.first <= prev.first  && second.first <= point.first && point.first <= prev.first) <<
							" " << (prev.second <= second.second && prev.second <= point.second && point.second <= second.second) <<
							" " << (second.second <= prev.second  && second.second <= point.second && point.second <= prev.second) << endl;

					cout << "intersection (" << point.first << "," << point.second << ")" << endl;
					inter_last = point;

					NDS_points.insert(make_pair(point,IntervalVector(1)));
				} else {
					cout << "intersection error" << endl;
				}
			} else {
				cout << "point very close to the after point" << endl;
			}

			prev = second;
		}

		cout << "------------" << endl;

		for(it1 = NDS_points.begin();it1 != NDS_points.end();++it1) {
			cout << "newpoint (" << it1->first.first << "," << it1->first.second << ")" << endl;

			if(it1->first.first != it1->first.first  || it1->first.second != it1->first.second ) getchar();
			NDS2.insert(*it1);
		}

		cout << "insert points" << endl;


		for(auto aux_p:DS2){
			if( ((aux_p.first <= p1.first && aux_p.first <= p2.first) &&
					(aux_p.second >= p1.second && aux_p.second >= p2.second)) ||
				((aux_p.second <= p1.second && aux_p.second <= p2.second) &&
						(it1->first.first >= p1.first && aux_p.first >= p2.first)) ) {
				cout << "fuera" << endl;
				cout << "(" << aux_p.first << "," << aux_p.second << ") " << endl;
				NDS2.insert(make_pair(aux_p,IntervalVector(1)));
				// addPointtoNDS(make_pair(it1->first.first,it1->first.second));
			} else {
				double m = (p1.second - p2.second)/(p1.first - p2.first);
				if( (p1.second - m*p1.first) >= (aux_p.second - m*aux_p.first)) {
					cout << "dentro" << endl;
					cout << "(" << aux_p.first << "," << aux_p.second << ") " << endl;
					NDS2.insert(make_pair(aux_p,IntervalVector(1)));
				} else {
					cout << "eliminado" << endl;
					cout << "(" << aux_p.first << "," << aux_p.second << ") " << endl;
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

	static pair<double, double> pointIntersection2(pair<double, double> p0, pair<double, double> p1,
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
