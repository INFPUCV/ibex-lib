/*
 * ibex_NDS.cpp
 *
 *  Created on: 24 may. 2018
 *      Author: iaraya
 */

#include "ibex_NDS.h"

namespace ibex {

	 map< pair <double, double>, IntervalVector, sorty2 > NDS_seg::NDS2;

	bool NDS_seg::is_dominated(vector<double>& p){
		pair< double, double> new_p=make_pair(p[0],p[1]);

		std::map<pair<double, double>, IntervalVector>::iterator it1 = --NDS2.lower_bound(new_p);
		// std::map<pair<double, double>, IntervalVector>::iterator it1 = NDS2.begin();
		pair< double, double> point1, point2;
		point1 = it1->first;
		it1++;
		point2 = it1->first;

		cout << "punto (" << new_p.first << "," << new_p.second << ")" << endl;

		// Se comprueba que no sea dominado por el anterior al lower_bound (puede ocurrir)
		if(point1.first <= new_p.first && point1.second <= new_p.second){
			return true;
		}

		// Se comprueba que no sea dominado por el lower_bound
		if(point2.first <= new_p.first && point2.second <= new_p.second )
			return true;


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

			if(cEval.ub() > c.lb())
				return true;


		} return false;

	}

	//TODO: think about how to associate the solutions to segments and points in the NDS
	void NDS_seg::addPoint(pair< double, double> new_p) {
		vector<double> p(2);
		p[0]=new_p.first;
		p[1]=new_p.second;

		if (is_dominated(p)) return;

		// Removes from NDS the points dominated by p
		// Then, adds the new point between the corresponding NDS points and adds new ones
		// intersecting the old segments

		std::map<pair<double, double>, IntervalVector>::iterator it1 = NDS2.lower_bound(new_p); // Se llega al nodo izquierdo del nodo eval, deberia estar mas arriba
		std::map<pair<double, double>, IntervalVector>::iterator aux;

		pair <double, double> first_dom=it1->first;
		it1--;
		pair <double, double> last_dom=it1->first;

		for(;it1 != NDS2.end();) {
			// termina cuando it1 no este dentro de los rangos dominados del punto a agregar
			if(it1->first.second < p[1]) break;

			// comprueba si esta dominado el punto para agregarlo a DS2
			if(p[0] <= it1->first.first and p[1] <= it1->first.second) {
				aux = it1;
				++aux;
				first_dom=it1->first;
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
		cout << 0 << endl;
		pair<double, double> intersection1 = pointIntersection(it1->first, first_dom, new_p, aux_p);
		cout << intersection1.first  << "," <<  intersection1.second << endl;
		aux_p = make_pair(POS_INFINITY, p[1]);
		cout << 10 << endl;
		cout << it2->first.first  << "," <<  it2->first.second << endl;
		cout << last_dom.first  << "," <<  last_dom.second << endl;
		cout << new_p.first  << "," <<  new_p.second << endl;
		cout << aux_p.first  << "," <<  aux_p.second << endl;
		pair<double, double> intersection2 = pointIntersection(last_dom, it2->first, new_p, aux_p);
		cout << intersection2.first  << "," <<  intersection2.second << endl;

		// se agregan el punto y los dos obtenidos anteriormente
		NDS2.insert(make_pair(new_p, IntervalVector(1)));
		NDS2.insert(make_pair(intersection1, IntervalVector(1)));
		NDS2.insert(make_pair(intersection2, IntervalVector(1)));

	}


	// Returns 1 if the lines intersect, otherwise 0. In addition, if the lines
	// intersect the intersection point may be stored in the floats i_x and i_y.
	pair<double, double> NDS_seg::pointIntersection2(pair<double, double> p0, pair<double, double> p1,
			pair<double, double> p2, pair<double, double> p3)
	{
		Interval p0_x=p0.first;
		Interval p0_y=p0.second;
		Interval p1_x=p1.first;
		Interval p1_y=p1.second;
		Interval p2_x=p2.first;
		Interval p2_y=p2.second;
		Interval p3_x=p3.first;
		Interval p3_y=p3.second;

		Interval i_x, i_y;

		cout << p0_x << "," << p0_y << endl;
		cout << p1_x << "," << p1_y << endl;
		cout << p2_x << "," << p2_y << endl;
		cout << p3_x << "," << p3_y << endl;

		//Cuando el segmento es vertical su pendiente es infinito
		if(p0_x==p1_x && p2_y==p3_y) {
			if( ((p2_x.mid() <= p0_x.mid() and p0_x.mid() <= p3_x.mid()) or
					(p3_x.mid() <= p0_x.mid() and p0_x.mid() <= p2_x.mid()))
			   && ((p0_y.mid() <= p2_y.mid() and p2_y.mid() <= p1_y.mid()) or
					(p1_y.mid() <= p2_y.mid() and p2_y.mid() <= p0_y.mid())) )

				return make_pair(p0_x.mid(),p2_y.mid());

			else throw NoIntersectionException();
		}
		else if (p0_y==p1_y && p2_x==p3_x){

			if( ((p0_x.mid() <= p2_x.mid() and p2_x.mid() <= p1_x.mid()) or
					(p1_x.mid() <= p2_x.mid() and p2_x.mid() <= p1_x.mid())) &&
				((p2_y.mid() <= p0_y.mid() and p0_y.mid() <= p3_y.mid()) or
					(p3_y.mid() <= p0_y.mid() and p0_y.mid() <= p2_y.mid())))
				return make_pair(p2_x.mid(),p0_y.mid());

			else throw NoIntersectionException();
		}



		Interval s1_x, s1_y, s2_x, s2_y, sn, tn, sd, td, t;
	    s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
	    s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

	    sn = -s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y);
	    sd = -s2_x * s1_y + s1_x * s2_y;
	    tn =  s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x);
	    td = -s2_x * s1_y + s1_x * s2_y;



	    if (sn.ub() >= 0 && sn.lb() <= sd.ub() && tn.ub() >= 0 && tn.lb() <= td.ub())
	    {
	        // Collision detected
	        t = tn / td;
            i_x = p0_x + (tn * s1_x);
            i_y = p0_y + (tn * s1_y);

            cout << i_x << endl;
            cout << i_y << endl;
	        return make_pair(i_x.ub(),i_y.ub());
	    }

	    throw NoIntersectionException();

	}

} /* namespace ibex */
