/*
 * ibex_NDS.cpp
 *
 *  Created on: 24 may. 2018
 *      Author: iaraya
 */

#include "ibex_NDS.h"
#include "ibex_NDS2.h"
#include "ibex_OptimizerMOP.h"

namespace ibex {

	 //map< pair <double, double>, IntervalVector, sorty2 > NDS_seg::NDS2;
	 bool NDS_seg::_trace;

	bool NDS_seg::is_dominated(const Vector& new_p){
		if(new_p[0] == POS_INFINITY && new_p[1] == POS_INFINITY) return false;

		std::map<Vector, NDS_data >::iterator it1 = --NDS2.lower_bound(new_p);

		Vector point1 = it1->first;
		it1++;
		Vector point2 = it1->first;

		// Se comprueba que no sea dominado por el anterior al lower_bound (puede ocurrir)
		if(point1[0] <= new_p[0] && point1[1] <= new_p[1])
			return true;


		// Se comprueba que no sea dominado por el lower_bound
		if(point2[0] <= new_p[0] && point2[1] <= new_p[1] )
			return true;

		// comprobar que no este dominado por la recta que forma los dos puntos anteriores
		if(!(new_p[0] <= point1[0] && new_p[1] <= point1[1] ) &&
				!(new_p[0] <= point2[0] && new_p[1] <= point2[1])) {

			//pendiente de los dos puntos
			Interval m = (Interval(point2[1])-Interval(point1[1]))/
					(Interval(point2[0])-Interval(point1[0]));

			// se obtiene el c de la funcion de los puntos
			Interval c = Interval(point1[1]) - m*Interval(point1[0]);

			// se obtiene el c del nuevo punto
			Interval cEval = Interval(new_p[1]) - m*Interval(new_p[0]);

			if(cEval.lb() > c.ub())
				return true;
		}

		return false;
	}

	bool NDS_seg::addSegment(const pair<Vector, Vector>& y1y2, const NDS_data& data) {

		const Vector& y1=y1y2.first;
		const Vector& y2=y1y2.second;

		if(y1y2.first == y1y2.second){
			addPoint(y1y2.first, data);
			return true;
		}

		//se insertan en DS2 todos los puntos que se ubican entre y1 y y2 incluyendo el anterior a y1 y el siguiente a y2
		std::map<Vector, NDS_data >::iterator aux, it1 = --NDS2.lower_bound(y1);


		Interval m = (Interval(y1[1]) - Interval(y2[1]))/
								(Interval(y1[0]) - Interval(y2[0]));
		double c_ub = (Interval(y1[1]) - m*Interval(y1[0])).ub();


		//segmentos en el rango
		list< pair<Vector, NDS_data> > DS2;
		DS2.push_back(*it1);
		it1++;

		for(;it1 != NDS2.end();) {

			DS2.push_back(*it1);
			if(it1->first[1] < y1[1] && it1->first[1] < y2[1]){
				 break;
			 }

			//se elimina el punto si es dominado por el segmento
			if( c_ub < (Interval(it1->first[1]) - m*Interval(it1->first[0])).lb()
			 && y1[0] < it1->first[0] && y2[1] < it1->first[1]){
				aux = it1; ++aux;
				NDS2.erase(it1);
				it1 = aux;

			} else it1++;
		}

		//se intersecta el segmento con los segmentos de la NDS
		//se agregan las intersecciones en NDS
		int intersections=0;
		pair<Vector, NDS_data> prev = DS2.front();
		for(auto next:DS2){
			if(prev.first==next.first) continue;

			try{
				double m2=-1e100;
				if(prev.first[0]-next.first[0]!=0)
				   m2=(prev.first[1]-next.first[1])/(prev.first[0]-next.first[0]);

				Vector point = pointIntersection(prev.first, next.first, y1, y2);
				if(m2<m.mid())
				   NDS2.insert(make_pair(point,prev.second));
				else
				   NDS2.insert(make_pair(point,data));
				intersections++;
			}catch(NoIntersectionException& e) {  }

			prev.first = next.first;
			prev.second = next.second;
		}

		return (intersections>0);

	}

	void NDS_seg::addPoint(const Vector& new_y, const NDS_data& data) {
		if (is_dominated(new_y)) return;

		// Removes from NDS the points dominated by p
		// Then, adds the new point between the corresponding NDS points and adds new ones
		// intersecting the old segments

		std::map<Vector, NDS_data >::iterator it1 = NDS2.lower_bound(new_y); // Se llega al nodo izquierdo del nodo eval, deberia estar mas arriba
		std::map<Vector, NDS_data >::iterator aux;

		Vector first_dom=it1->first;
		it1--;
		Vector last_dom=it1->first;
		bool first=true;
		NDS_data prev_data;
		NDS_data next_data;
		for(;it1 != NDS2.end();) {
			// termina cuando it1 no este dentro de los rangos dominados del punto a agregar
			if(it1->first[1] < new_y[1]) break;

			// comprueba si esta dominado el punto
			if(new_y[0] <= it1->first[0] && new_y[1] <= it1->first[1]) {
				aux = it1;
				++aux;
				if(first) {
					prev_data=it1->second;
					first_dom=it1->first;
				}
				first=false;
				last_dom=it1->first;
				next_data = it1->second;

				NDS2.erase(it1);
				it1 = aux;
			} else ++it1;
		}

		std::map<Vector, NDS_data >::iterator it2 = --NDS2.lower_bound(new_y);
		it1 = it2;
		it2++;

		Vector aux_y(2);

		aux_y[0]=new_y[0]; aux_y[1]=POS_INFINITY;
		Vector intersection1 = pointIntersection(it1->first, first_dom, new_y, aux_y);
		aux_y[0]=POS_INFINITY; aux_y[1]=new_y[1];
		Vector intersection2 = pointIntersection(last_dom, it2->first, new_y, aux_y);

		// se agregan el punto y los dos obtenidos anteriormente
		NDS2.insert(make_pair(new_y, data));
		NDS2.insert(make_pair(intersection1, NDS_data()));
		NDS2.insert(make_pair(intersection2, next_data));
	}


	// Returns 1 if the lines intersect, otherwise 0. In addition, if the lines
	// intersect the intersection point may be stored in the floats i_x and i_y.
	Vector NDS_seg::pointIntersection(const Vector& p0, const Vector& p1, const Vector& p2, const Vector& p3)
	{

		double p0_x=p0[0];
		double p0_y=p0[1];
		double p1_x=p1[0];
		double p1_y=p1[1];
		double p2_x=p2[0];
		double p2_y=p2[1];
		double p3_x=p3[0];
		double p3_y=p3[1];


		Interval i_x, i_y;
		pair<double,double> i ;


		//Cuando el segmento es vertical su pendiente es infinito
		if(p0_x==p1_x && p2_y==p3_y) {
			if( ((p2_x <= p0_x and p0_x <= p3_x) or
					(p3_x <= p0_x and p0_x <= p2_x))
			   && ((p0_y <= p2_y and p2_y <= p1_y) or
					(p1_y <= p2_y and p2_y <= p0_y)) )

				i = make_pair(p0_x,p2_y);

			else throw NoIntersectionException();
		}
		else if (p0_y==p1_y && p2_x==p3_x){

			if( ((p0_x <= p2_x and p2_x <= p1_x) or
					(p1_x <= p2_x and p2_x <= p1_x)) &&
				((p2_y <= p0_y and p0_y <= p3_y) or
					(p3_y <= p0_y and p0_y <= p2_y)))

					i = make_pair(p2_x,p0_y);

			else throw NoIntersectionException();
		}else{

      if(p0_y==POS_INFINITY) p0_y=std::max(p2_y,p3_y);
			if(p0_y==POS_INFINITY) p0_y=std::min(p2_y,p3_y);
			if(p1_y==POS_INFINITY) p1_y=std::max(p2_y,p3_y);
			if(p1_y==POS_INFINITY) p1_y=std::min(p2_y,p3_y);
			if(p0_x==POS_INFINITY) p0_x=std::max(p2_x,p3_x);
			if(p0_x==POS_INFINITY) p0_x=std::min(p2_x,p3_x);
			if(p1_x==POS_INFINITY) p1_x=std::max(p2_x,p3_x);
			if(p1_x==POS_INFINITY) p1_x=std::min(p2_x,p3_x);
			if(p2_y==POS_INFINITY) p2_y=std::max(p0_y,p1_y);
			if(p3_y==POS_INFINITY) p3_y=std::max(p0_y,p1_y);
			if(p2_x==POS_INFINITY) p2_x=std::max(p0_x,p1_x);
			if(p3_x==POS_INFINITY) p3_x=std::max(p0_x,p1_x);

		  Interval r_x, r_y, s_x, s_y, u, t, p_x=p0_x, p_y=p0_y, q_x=p2_x, q_y=p2_y;
	    r_x = p1_x - p_x;     r_y = p1_y - p_y;
	    s_x = p3_x - q_x;     s_y = p3_y - q_y;

			Interval rxs = -s_x * r_y + r_x * s_y;


			if(!rxs.contains(0)){
	    	//u = (-r_y * (p_x - q_x) + r_x * (p_y - q_y)) / rxs;  //(p-q) x r /rxs
 	    	t = ( s_x * (p_y - q_y) - s_y * (p_x - q_x)) / rxs;  //(p-q) x s /rxs
			}else if((-r_y * (p_x - q_x) + r_x * (p_y - q_y)).contains(0)) {
				//colinear
				Interval rxr = (r_x*r_x + r_y*r_y);
				Interval t0= (q_x-p_x)*r_x + (q_y-p_y)*r_y / (r_x*r_x + r_y*r_y); // (q ��� p) �� r / (r �� r)
				Interval t1= t0 + s_x*r_x + s_y*r_y / (r_x*r_x + r_y*r_y); // t0 + s �� r / (r �� r)

				t = Interval(0,1);
				t &= Interval(std::min(t0.lb(),t1.lb()),std::max(t0.ub(),t1.ub()));
			}else //parallel
				throw NoIntersectionException();


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

		Vector v(2); v[0]=i.first; v[1]=i.second;
		return v;

	}

//TODO: Handle of ndsh2

//Return True si y esta dominado por algun punto de NDS
	bool NDS_XY::is_dominated(const Vector& new_point_y){
		for(auto itr = NDS.begin(); itr != NDS.end(); ++itr) {
			if(is_dominated(new_point_y, itr->first)) return true;
		}
		return false;
	}

//Return True si y_prime es dominado por y
	bool NDS_XY::is_dominated(const Vector& y_prime, const Vector& y){
		for(int i=0; i < y_prime.size(); i++){
			if(y[i] <= y_prime[i]){
			}else return false;
		}
		return true;
	}

	void NDS_XY::addPoint(const Vector& new_y, const NDS_X& data){
		if (is_dominated(new_y)) return;

		//Removes from NDS the points dominated by new_y
		auto it = NDS.begin();
		while ( it != NDS.end()){
			if(is_dominated(it->first, new_y)){
				it = NDS.erase(it);
			}else
				++it;
		}
		NDS.insert(make_pair(new_y, data));
	}


	IntervalVector get_box_y(const IntervalVector &box, int nFuncObj){
		IntervalVector box_y(nFuncObj);
		int n=box.size()-nFuncObj;
		for(int i=0; i< nFuncObj; i++) box_y[i]=box[n+i];
		return box_y;
	}
	double NDS_XY::distance(const IntervalVector& box){
		Vector y( get_box_y(box, nObjFunc).lb() );
		//Descarta la caja si boxy.lb se encuentra dominado por un punto de ndsh
		if (is_dominated(y)) {
			return -1;
		}

		/*
		 * Para cada punto en ndsh
		 * 		Primero se calcula la distancia entre ndsh[i] con boxy_lb en cada dimension
		 * 		Luego, para cada distancia obtenida en cada dimension (conponent_distance) se verifica,
		 * 			Si la caja boxy_lb + component_distance se encuentra dominada por el punto ndsh[i] (la distancia es factible)
		 * 		Para cada distancia factible se selecciona la menor entre ellas (minDistance)
		 * Finalmente  entre todas las minDistance obtenidas, se selecciona el mayor valor entre ellas (maxDistance)
		 */
		double maxDistance=POS_INFINITY;
		for(auto itr = NDS.begin(); itr != NDS.end(); ++itr){
			double PartialMaxDistance = NEG_INFINITY;
			for(int i=0; i<itr->first.size(); i++){
				//calculo de distancia de cada dimension
				double component_distance = itr->first[i] - y[i];
				if(component_distance > 0){
					//Obtain minDistance
					if(component_distance > PartialMaxDistance) {
						PartialMaxDistance = component_distance;
					}
				}
			}
			//Save maxDistance
			if( maxDistance > PartialMaxDistance ) maxDistance = PartialMaxDistance;
		}


		return maxDistance;
	}

	list< pair< Vector, int> > NDS_XY::cutting_points(const Vector& boxy_lb, const Vector& boxy_ub){

		//Initialize best cutting points
		Array< Vector > cp(boxy_lb.size());
		Vector init(boxy_lb.size(), POS_INFINITY);
		for(int i=0; i<boxy_lb.size(); i++) cp.set_ref(i, init);

		for(auto it=NDS.begin(); it != NDS.end(); ++it){
			int sum = 0;
			int dim = 0;

			for(int i=0; i< boxy_lb.size(); i++){
				//We take the points outside of the box with these two conditions
				//first condition
				if( boxy_ub[i] > it->first[i]) {
					//second condition
					if(boxy_lb[i] >= it->first[i]) {
						sum++;
					}
					else dim = i;
				}else {
					i = boxy_lb.size();
					sum = 0;
				}
			}
			//We add the points only if for every component of the point ndsh hold
			//- Satisfy the first condition
			//- Satisfy the second condition nFuncObj-1 times
			if (sum == boxy_lb.size() - 1) {
				Vector aux = cp[dim];
				if(it->first[dim] < aux[dim]) cp.set_ref(dim, (Vector &)it->first);
			}
		}

		list< pair<Vector, int>> cutting_points;
		for(int i=0; i<boxy_lb.size(); i++){
			double aux = POS_INFINITY;
			Vector aux2 = cp[i];
			if(aux2[0] != aux){
				cutting_points.push_back(make_pair(cp[i], i));
			}
		}

		return cutting_points;
	}

} /* namespace ibex */
