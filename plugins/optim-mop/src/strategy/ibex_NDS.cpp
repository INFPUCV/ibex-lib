/*
 * ibex_NDS.cpp
 *
 *  Created on: 24 may. 2018
 *      Author: iaraya
 */

#include "ibex_NDS.h"

namespace ibex {

	 map< pair <double, double>, IntervalVector, sorty2 > NDS_seg::NDS2;

	void NDS_seg::addPoint(pair< double, double> eval) {

		std::map<pair<double, double>, IntervalVector>::iterator it1 = --NDS2.lower_bound(eval);
		// std::map<pair<double, double>, IntervalVector>::iterator it1 = NDS2.begin();
		pair< double, double> point1, point2;
		point1 = it1->first;
		it1++;
		point2 = it1->first;

		cout << "punto (" << eval.first << "," << eval.second << ")" << endl;

		// Se comprueba que no sea dominado por el anterior al lower_bound (puede ocurrir)
		if(point1.first <= eval.first && point1.second <= eval.second){
			return;
		}

		// Se comprueba que no sea dominado por el lower_bound
		if(point2.first <= eval.first && point2.second <= eval.second )
			return;


		// comprobar que no este dominado por la recta que forma los dos puntos anteriores
		if(!(eval.first <= point1.first && eval.second <= point1.second ) &&
				!(eval.first <= point2.first && eval.second <= point2.second)) {

			//pendiente de los dos puntos
			Interval m = (Interval(point2.second)-Interval(point1.second))/
					(Interval(point2.first)-Interval(point1.first));

			// se obtiene el c de la funcion de los puntos
			Interval c = Interval(point1.second) - m*Interval(point1.first);

			// se obtiene el c del nuevo punto
			Interval cEval = Interval(eval.second) - m*Interval(eval.first);

			if(cEval.ub() > c.lb())
				return;


		} //else eval domina uno de los puntos


		//TODO: think about how to associate the solutions to segments and points in the NDS
		IntervalVector vec(1);

		// insertar en NDS2
		// Agregar a la lista del DS
		//- se revisa si it1 es el comienzo, si no lo es se retrocede uno
		//- se revisa si es dominado el it1, en el caso que lo sea se elimina y se guarda en el set
		//- lo anterior se realiza hasta que el eje y del it1 sea menor que el eval

		it1 = --NDS2.lower_bound(eval); // Se llega al nodo izquierdo del nodo eval, deberia estar mas arriba
		if(it1 != NDS2.begin()) it1--; // se retrocede 1 si no es el primero

		// agrega todos los puntos dominados por el punto a DS2 y los elimina de NDS2
		std::map<pair<double, double>, IntervalVector>::iterator aux;
		std::map<pair <double, double>, IntervalVector, sorty2 > DS2;
		std::map<pair<double, double>, IntervalVector>::iterator beginit, endit;
		std::map<pair<double, double>, IntervalVector>::iterator it2 = --NDS2.lower_bound(eval);
		for(;it1 != NDS2.end();) {
			// termina cuando it1 no este dentro de los rangos dominados del punto a agregar
			if(it1->first.second < eval.second) break;

			// comprueba si esta dominado el punto para agregarlo a DS2
			if(eval.first <= it1->first.first and eval.second <= it1->first.second) {
				aux = it1;
				++aux;
				DS2.insert(*it1);
				NDS2.erase(it1);
				it1 = aux;
			} else ++it1;
		}

		if(DS2.size() > 0) {
			beginit = DS2.begin();
			endit = --DS2.end();
		} else {
			endit = it2;
			it2++;
			beginit = it2;
		}


		// se obtienen los 2 nuevos puntos
		it2 = --NDS2.lower_bound(eval);
		it1 = it2++;

		pair<double, double> eval2;

		eval2 = make_pair(eval.first, POS_INFINITY);
		pair<double, double> intersection1 = pointIntersection(it1->first, beginit->first, eval, eval2);
		cout << intersection1.first  << "," <<  intersection1.second << endl;
		eval2 = make_pair(POS_INFINITY, eval.second);
		pair<double, double> intersection2 = pointIntersection(it2->first, endit->first, eval, eval2);
		cout << intersection2.first  << "," <<  intersection2.second << endl;

		// se agregan el punto y los dos obtenidos anteriormente
		NDS2.insert(make_pair(eval, vec));
		NDS2.insert(make_pair(intersection1, vec));
		NDS2.insert(make_pair(intersection2, vec));

	}
} /* namespace ibex */
