/*
 * ibex_NDS.cpp
 *
 *  Created on: 24 may. 2018
 *      Author: iaraya
 */

#include "ibex_NDS.h"

namespace ibex {

	 map< pair <double, double>, IntervalVector, sorty2 > NDS_seg::NDS2;
	 bool NDS_seg::_trace;

	bool NDS_seg::is_dominated(pair< double, double> new_p){
		if(new_p.first == POS_INFINITY && new_p.second == POS_INFINITY) return false;

		std::map<pair<double, double>, IntervalVector>::iterator it1 = --NDS2.lower_bound(new_p);
		// std::map<pair<double, double>, IntervalVector>::iterator it1 = NDS2.begin();
		pair< double, double> point1, point2;
		point1 = it1->first;
		it1++;
		point2 = it1->first;

		// Se comprueba que no sea dominado por el anterior al lower_bound (puede ocurrir)
		if(point1.first <= new_p.first && point1.second <= new_p.second){
			//cout << "is_dom: prev_dom" << endl;
			return true;
		}

		// Se comprueba que no sea dominado por el lower_bound
		if(point2.first <= new_p.first && point2.second <= new_p.second ){
			//cout << "is_dom: lb_dom" << endl;
			return true;
		}


		// comprobar que no este dominado por la recta que forma los dos puntos anteriores
		if(!(new_p.first <= point1.first && new_p.second <= point1.second ) &&
				!(new_p.first <= point2.first && new_p.second <= point2.second)) {

			//pendiente de los dos puntos
			Interval m = (Interval(point2.second)-Interval(point1.second))/
					(Interval(point2.first)-Interval(point1.first));

			// se obtiene el c de la funcion de los puntos
			Interval c = Interval(point1.second) - m*Interval(point1.first);

			// se obtiene el c del nuevo punto
			Interval cEval = Interval(new_p.second) - m*Interval(new_p.first);


			if(cEval.lb() > c.ub()){
				return true;
			}


		}

		return false;

	}


	bool NDS_seg::is_dominated(vector<double>& p){
			pair< double, double> new_p=make_pair(p[0],p[1]);
			return NDS_seg::is_dominated(new_p);
	}



	void NDS_seg::addSegment(pair< double, double> p1, pair< double, double> p2) {
		//cout << "add_segment:" << p1.first << "," << p1.second << " -- "	<< p2.first << "," << p2.second  <<endl;

		if(p1.first == p2.first  &&  p1.second == p2.second ){
			addPoint(p1);
			return;
		}

		//se insertan en DS2 todos los puntos que se ubican entre p1 y p2 incluyendo en anterior en y1 y el siguiente en y2
		std::map<pair<double, double>, IntervalVector>::iterator aux, it1 = --NDS2.lower_bound(p1);
		pair< double, double> second, point, inter_last = make_pair(NEG_INFINITY, NEG_INFINITY);
		list< pair <double, double>> DS2;

		Interval m = (Interval(p1.second) - Interval(p2.second))/
								(Interval(p1.first) - Interval(p2.first));
		double c_ub = (Interval(p1.second) - m*Interval(p1.first)).ub();

		DS2.push_back(it1->first);
		it1++;

		bool flagDS2 = false;
		for(;it1 != NDS2.end();) {
			if(it1->first.second < p1.second and it1->first.second < p2.second) flagDS2= true;
			DS2.push_back(it1->first);
			if(flagDS2) break;

      //se elimina el punto si es dominado por el segmento
			if( c_ub < (Interval(it1->first.second) - m*Interval(it1->first.first)).lb()){
				aux = it1; ++aux;
				NDS2.erase(it1);
				it1 = aux;
			} else it1++;
		}

    //se intersecta el segmento con los segmentos de la NDS
		//se agregan las intersecciones en NDS
		pair< double, double> prev = DS2.front();
		for(auto second:DS2){
			if(prev==second) continue;

			try{
				point = pointIntersection(prev, second, p1, p2);
				NDS2.insert(make_pair(point,IntervalVector(1)));
			}catch(NoIntersectionException& e) {
			}

			prev = second;
		}

		// getchar();
		//std::vector< pair <double, double> > curve_y;
		//std::vector< pair <double, double> > rectaUB;
		//py_Plotter::offline_plot(NULL, NDS2, rectaUB, curve_y);
	}

	//TODO: think about how to associate the solutions to segments and points in the NDS
	void NDS_seg::addPoint(pair< double, double> new_p) {
		vector<double> p(2);
		p[0]=new_p.first;
		p[1]=new_p.second;

		//cout << "add_point:" << p[0] << "," << p[1] << endl;
        //cout << "add_point: lp:" << (--NDS2.end())->first.first << "," << (--NDS2.end())->first.second << endl;

		if (is_dominated(p)) return;
		//cout << "add_point: non_dominated" << endl;

		// Removes from NDS the points dominated by p
		// Then, adds the new point between the corresponding NDS points and adds new ones
		// intersecting the old segments

		std::map<pair<double, double>, IntervalVector>::iterator it1 = NDS2.lower_bound(new_p); // Se llega al nodo izquierdo del nodo eval, deberia estar mas arriba
		std::map<pair<double, double>, IntervalVector>::iterator aux;

		pair <double, double> first_dom=it1->first;
		it1--;
		pair <double, double> last_dom=it1->first;
    bool first=true;
		for(;it1 != NDS2.end();) {
			// termina cuando it1 no este dentro de los rangos dominados del punto a agregar
			if(it1->first.second < p[1]) break;

			// comprueba si esta dominado el punto para agregarlo a DS2
			if(p[0] <= it1->first.first and p[1] <= it1->first.second) {
				aux = it1;
				++aux;
				if(first) first_dom=it1->first;
				first=false;
				last_dom=it1->first;

				NDS2.erase(it1);
				it1 = aux;
			} else ++it1;
		}

		std::map<pair<double, double>, IntervalVector>::iterator it2 = --NDS2.lower_bound(new_p);
		it1 = it2;
		it2++;

		pair<double, double> aux_p;

		aux_p = make_pair(p[0], POS_INFINITY);
		pair<double, double> intersection1 = pointIntersection(it1->first, first_dom, new_p, aux_p);
		aux_p = make_pair(POS_INFINITY, p[1]);
		pair<double, double> intersection2 = pointIntersection(last_dom, it2->first, new_p, aux_p);
		//cout << "add_point: lp:" << (--NDS2.end())->first.first << "," << (--NDS2.end())->first.second << endl;
		// se agregan el punto y los dos obtenidos anteriormente
		NDS2.insert(make_pair(new_p, IntervalVector(1)));
		NDS2.insert(make_pair(intersection1, IntervalVector(1)));
		NDS2.insert(make_pair(intersection2, IntervalVector(1)));


		//std::vector< pair <double, double> > curve_y;
		//std::vector< pair <double, double> > rectaUB;
		//py_Plotter::offline_plot(NULL, NDS2, rectaUB, curve_y);
		// getchar();


	}


  void NDS_seg::remove_infinity(pair<double, double>& p){
		if(p.first==NEG_INFINITY) p.first=-1e20;
		if(p.second==NEG_INFINITY) p.second=-1e20;
		if(p.first==POS_INFINITY) p.first=1e20;
		if(p.second==POS_INFINITY) p.second=1e20;
	}

	void NDS_seg::add_infinity(pair<double, double>& p){
		if(p.first==-1e20) p.first=NEG_INFINITY;
		if(p.second==-1e20) p.second=NEG_INFINITY;
		if(p.first==1e20) p.first=POS_INFINITY;
		if(p.second==1e20) p.second=POS_INFINITY;
	}

	// Returns 1 if the lines intersect, otherwise 0. In addition, if the lines
	// intersect the intersection point may be stored in the floats i_x and i_y.
	pair<double, double> NDS_seg::pointIntersection(pair<double, double> p0, pair<double, double> p1,
			pair<double, double> p2, pair<double, double> p3)
	{
		remove_infinity(p0);
		remove_infinity(p1);
		remove_infinity(p2);
		remove_infinity(p3);

		Interval p0_x=p0.first;
		Interval p0_y=p0.second;
		Interval p1_x=p1.first;
		Interval p1_y=p1.second;
		Interval p2_x=p2.first;
		Interval p2_y=p2.second;
		Interval p3_x=p3.first;
		Interval p3_y=p3.second;

		Interval i_x, i_y;
		pair<double,double> i ;

		//Cuando el segmento es vertical su pendiente es infinito
		if(p0_x==p1_x && p2_y==p3_y) {
			if( ((p2_x.mid() <= p0_x.mid() and p0_x.mid() <= p3_x.mid()) or
					(p3_x.mid() <= p0_x.mid() and p0_x.mid() <= p2_x.mid()))
			   && ((p0_y.mid() <= p2_y.mid() and p2_y.mid() <= p1_y.mid()) or
					(p1_y.mid() <= p2_y.mid() and p2_y.mid() <= p0_y.mid())) )

				i = make_pair(p0_x.mid(),p2_y.mid());

			else throw NoIntersectionException();
		}
		else if (p0_y==p1_y && p2_x==p3_x){

			if( ((p0_x.mid() <= p2_x.mid() and p2_x.mid() <= p1_x.mid()) or
					(p1_x.mid() <= p2_x.mid() and p2_x.mid() <= p1_x.mid())) &&
				((p2_y.mid() <= p0_y.mid() and p0_y.mid() <= p3_y.mid()) or
					(p3_y.mid() <= p0_y.mid() and p0_y.mid() <= p2_y.mid())))

					i = make_pair(p2_x.mid(),p0_y.mid());

			else throw NoIntersectionException();
		}else{
		  Interval r_x, r_y, s_x, s_y, u, t, p_x=p0_x, p_y=p0_y, q_x=p2_x, q_y=p2_y;
	    r_x = p1_x - p0_x;     r_y = p1_y - p0_y;
	    s_x = p3_x - p2_x;     s_y = p3_y - p2_y;

			Interval rxs = -s_x * r_y + r_x * s_y;
			if(!rxs.contains(0)){
	    	//u = (-r_y * (p_x - q_x) + r_x * (p_y - q_y)) / rxs;  //(p-q) x r /rxs
 	    	t = ( s_x * (p_y - q_y) - s_y * (p_x - q_x)) / rxs;  //(p-q) x s /rxs
			}else if((-r_y * (p_x - q_x) + r_x * (p_y - q_y)).contains(0)) {
				//colinear
				Interval rxr = (r_x*r_x + r_y*r_y);
				Interval t0= (q_x-p_x)*r_x + (q_y-p_y)*r_y / (r_x*r_x + r_y*r_y); // (q − p) · r / (r · r)
				Interval t1= t0 + s_x*r_x + s_y*r_y / (r_x*r_x + r_y*r_y); // t0 + s · r / (r · r)

				t = Interval(0,1);
				t &= Interval(std::min(t0.lb(),t1.lb()),std::max(t0.ub(),t1.ub()));
			}


	    if (/*u.ub() >= 0 && u.lb() <= 1 &&*/ t.ub() >= 0 && t.lb() <= 1)
	    {
	        // Collision detected
	        i_x = p_x + (t * r_x);
	        i_y = p_y + (t * r_y);

					if (p0_x==p1_x) i_x=p0_x;
					if (p2_x==p3_x) i_x=p2_x;
					if (p0_y==p1_y) i_y=p0_y;
					if (p2_y==p3_y) i_y=p2_y;

					i = make_pair(i_x.ub(),i_y.ub());

	    }else{
				throw NoIntersectionException();
			}
		}

		add_infinity(i);

		return i;



	}

} /* namespace ibex */
