/*
 * ibex_NDS.cpp
 *
 *  Created on: 13 apr. 2021
 *      Author: iaraya
 */

#include "ibex_NDSrp.h"

namespace ibex {

	void NDSrp::clear(){
		NDS_clear();
		//the first point
		Point* first=new Point(NEG_INFINITY,POS_INFINITY);
		Point* middle=new Point(POS_INFINITY,POS_INFINITY);
		Point* last=new Point(POS_INFINITY,NEG_INFINITY);
		middle->push_before(last);
		first->push_before(middle);

        NDS.insert(last);
		NDS.insert(middle);
		NDS.insert(first);
	}

	bool NDSrp::is_dominated(const Vector& new_p){
        if(new_p[0] == POS_INFINITY && new_p[1] == POS_INFINITY) return false;

		std::set<Vector*>::iterator it1 = --NDS.lower_bound((Vector*)&new_p);

		Vector point1 = (**it1);
		it1++;
		Vector point2 = (**it1);

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


	void NDSrp::addPoint(const Vector& new_y){
        if (is_dominated(new_y)) return;
		// Removes from NDS the points dominated by p
		// Then, adds the new point between the corresponding NDS points and adds new ones
		// intersecting the old segments
		auto it1 = NDS.lower_bound((Vector*)&new_y); // Se llega al nodo izquierdo del nodo eval, deberia estar mas arriba
		std::set<Vector *>::iterator aux;

		Vector first_dom=(**it1);
		it1--;
		Vector last_dom=(**it1);
		bool first=true;

		for(;it1 != NDS.end();) {
			// termina cuando it1 no este dentro de los rangos dominados del punto a agregar
			if((**it1)[1] < new_y[1]) break;

			// comprueba si esta dominado el punto
			if(new_y[0] <= (**it1)[0] && new_y[1] <= (**it1)[1]) {
				aux = it1;
				++aux;
				if(first) {
					first_dom=(**it1);
				}
				first=false;
				last_dom=(**it1);
				NDS_erase(it1);
				it1 = aux;
			} else ++it1;
		}

		auto it2 = --NDS.lower_bound((Vector*)&new_y);
		it1 = it2;
		it2++;

		Vector aux_y(2);

		aux_y[0]=new_y[0]; aux_y[1]=POS_INFINITY;

		try{
			Vector intersection1 = pointIntersection(**it1, first_dom, new_y, aux_y);
			aux_y[0]=POS_INFINITY; aux_y[1]=new_y[1];
			Vector intersection2 = pointIntersection(last_dom, **it2, new_y, aux_y);

			// se agregan el punto y los dos obtenidos anteriormente
			NDS_insert(new_y);
			NDS_insert(intersection1);
			NDS_insert(intersection2);
		}catch(NoIntersectionException& e) {  }

    }

	/**
	* Add a segment in the NDS structure. Return 0 if the segment did not modify the NDS
	*/
	bool NDSrp::addSegment(const pair< Vector, Vector>& y1y2){
        const Vector& y1=y1y2.first;
		const Vector& y2=y1y2.second;

		if(y1y2.first == y1y2.second){
			addPoint(y1y2.first);
			return true;
		}

		//se insertan en DS2 todos los puntos que se ubican entre y1 y y2 incluyendo el anterior a y1 y el siguiente a y2
		std::set<Vector* >::iterator aux, it1 = --NDS.lower_bound((Vector*)&y1);

		Interval m = (Interval(y1[1]) - Interval(y2[1]))/
								(Interval(y1[0]) - Interval(y2[0]));
		double c_ub = (Interval(y1[1]) - m*Interval(y1[0])).ub();

		//segmentos en el rango
		list< Vector > DS2;
		DS2.push_back(**it1);
		it1++;

		for(;it1 != NDS.end();) {

			DS2.push_back(**it1);
			if((**it1)[1] < y1[1] && (**it1)[1] < y2[1]){
				 break;
			 }

			//se elimina el punto si es dominado por el segmento
			if( c_ub < (Interval((**it1)[1]) - m*Interval((**it1)[0])).lb()
			 && y1[0] < (**it1)[0] && y2[1] < (**it1)[1]){
				aux = it1; ++aux;
				NDS_erase(it1);
				it1 = aux;

			} else it1++;
		}

		//se intersecta el segmento con los segmentos de la NDS
		//se agregan las intersecciones en NDS
		int intersections=0;
		Vector prev = DS2.front();
		for(auto next:DS2){
			if(prev==next) continue; //esto ocurre la primera vez

			try{
				double m2=-1e100;
				if(prev[0]-next[0]!=0)
				   m2=(prev[1]-next[1])/(prev[0]-next[0]);

				Vector point = pointIntersection(prev, next, y1, y2);
				if(m2<m.mid())
				   NDS_insert(point);
				else
				   NDS_insert(point);
				intersections++;
			}catch(NoIntersectionException& e) {  }

			prev = next;
		}

		return (intersections>0);
    }

    //returns a list of points non-dominated by lb
	list< Vector > NDSrp::non_dominated_points(const Vector& lb, bool cutting_points){
        list < Vector > inpoints;

		//first potential dominated point
		Vector v(2); v[0]=lb[0]; v[1]=NEG_INFINITY;
		auto ent1=NDS.upper_bound((Vector*)&v);
		Vector v11 = (**ent1); ent1--;
		//last point before x=lbx
		Vector v10 = (**ent1);

		Vector v1(2); v1[0]=lb[0]; v1[1]=v11[1];
		Vector v2(2); v2[0]=lb[0]; v2[1]=v10[1];

		//x-point cutting lbx
		if(cutting_points){
			Vector firstp(2);
		    firstp = pointIntersection( v10, v11, v1, v2);
		   	if(firstp[1] <= lb[1] ){
		   		return inpoints; //empty list
		   	}
		   inpoints.push_back(firstp);
	    }

		ent1++;
		//points dominated by lb
		while((**ent1)[1] > lb[1]){
			inpoints.push_back((**ent1));
			ent1++;
		}

		if(cutting_points){
	    Vector lastp(2);
			//last potential dominated point
			ent1--;
			Vector v21 = (**ent1);
			//first point after y=lby
			ent1++;
			Vector v20 = (**ent1);

			//x-point cutting lby
		    v1[0]=v20[0]; v1[1]=lb[1];
		    v2[0]=v21[0]; v2[1]=lb[1];
			lastp = pointIntersection( v20, v21, v1, v2);
			inpoints.push_back(lastp);
	    }

		return inpoints;
    }

    /**
        * \brief Returns the point intersecting two segments. Otherwise it throw
        * a NoIntersectionException
        * It is conservative, that is:
        * 1) if there are intersection it should return a point dominated by the real intersection
        * 2) if there are no intersection it may return an exception or a point dominated by one segment
        * TODO: revise with colinear generated examples
	*/
	Vector NDSrp::pointIntersection(const Vector& p0, const Vector& p1, const Vector& p2, const Vector& p3){
        double p0_x=p0[0]; double p1_x=p1[0]; double p2_x=p2[0]; double p3_x=p3[0];
		double p0_y=p0[1]; double p1_y=p1[1]; double p2_y=p2[1]; double p3_y=p3[1];


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
		}else if (p0_y==p1_y && p2_x==p3_x){

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
	        }else throw NoIntersectionException();
		}

		Vector v(2); v[0]=i.first; v[1]=i.second;
		return v;
    }

    double NDSrp::distance(const Cell* c){
		IntervalVector box_y=OptimizerMOP::get_box_y(c);
		double a = ((BxpMOPData*) c->prop[BxpMOPData::id])->a;
		double w_lb = ((BxpMOPData*) c->prop[BxpMOPData::id])->w_lb;


		double dist;
		IntervalVector points(4);
		if(a!=0){
			pair<Vector, Vector> points = get_segment(box_y.lb(),-1.0/a, w_lb/a);
			dist= distance(points.first, points.second, -1.0/a, w_lb/a);
		}else{
			pair<Vector, Vector> points = get_segment(box_y.lb());
			dist= distance(points.first, points.second);
		}

		return std::max(0.0,dist);
    }

	// m in [-oo, 0]
	double NDSrp::distance(Vector& yA, Vector& yB, double m, double c){
		Interval Ax=yA[0]; Interval Ay=yA[1];
		Interval Bx=yB[0]; Interval By=yB[1];

		double max_dist=NEG_INFINITY;

		Vector lb(2); lb[0]=yA[0], lb[1]=yB[1];
		list< Vector> inner_segments= non_dominated_points(lb);

		Vector* p0=NULL;

		bool Adist=false;
		bool Bdist=false;

		for(auto p : inner_segments){
			if(p[0]==POS_INFINITY && p[1]==POS_INFINITY) return POS_INFINITY;

			Interval dist;
			//up-left point
			if((p[0]-Ax).lb() < (p[1]-Ay).ub() || p[1]==POS_INFINITY){
				//cout << "up-left: " << p[0] << "," << p[1] << endl;
				dist=p[0]-Interval(Ax);
				//cout << dist << endl;
			}
			//bottom-right point
			else if((p[1]-By).lb() < (p[0]-Bx).ub() || p[0]==POS_INFINITY){
				//cout << "bottom-right: " << p[0] << "," << p[1] << endl;
				dist=p[1]-Interval(By);
				//cout << dist << endl;
				if(!Bdist && p0){
					Interval mm=NEG_INFINITY;
					if(p[0]-(*p0)[0] != 0)
						mm= (Interval(p[1])-(*p0)[1])/(Interval(p[0])-(*p0)[0]);

                     //cout << mm << endl;
					if(/*mm.lb() > -1 &&*/ mm.lb() < 0.0 ){
						Interval cc= p[1] - mm*p[0];
						//cout << ((mm*Ax - Ay + cc)/(1.0-mm)).lb() << endl;
						dist=std::max(dist.ub(), ((mm*Ax - Ay + cc)/(1.0-mm)).lb());
						Bdist=true;
					}
				}
			}
			//cy-45-degree zone
			else{
				//cout << "45 zone: " << p[0] << "," << p[1] << endl;
				dist= -(m*p[0] - p[1]+c)/(1.0-m);
				if(!Adist && p0){
					Interval mm=NEG_INFINITY;
					if(p[0]-(*p0)[0] != 0)
						mm= (Interval(p[1])-(*p0)[1])/(Interval(p[0])-(*p0)[0]);

					if(/*mm.lb() <-1 &&*/ mm.lb()>NEG_INFINITY){
						Interval cc= p[1] - mm*p[0];
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
			p0=new Vector(p);

		}

		if(p0) delete p0;
		return max_dist;
    }

    //returns the segment yA--yB of the line y_2=m*y_2+c dominated by lb
    pair <Vector, Vector> NDSrp::get_segment(const Vector& lb, double m, double c){
        double max_dist=NEG_INFINITY;

		Interval Ay=lb[1];
		Interval Bx=lb[0];

		if(m!=POS_INFINITY){
			Ay = Interval(m)*lb[0]+c;  // Ax=lbx
			Bx = (Interval(lb[1])-c)/m; // By=lby
			if(Ay.lb() < lb[1]){ Ay=lb[1]; }
			if(Bx.lb() < lb[0]){ Bx=lb[0]; }
		}

		Vector yA(2);
		Vector yB(2);

		yA[0]=lb[0];
		yA[1]=Ay.lb();
		yB[0]=Bx.lb();
		yB[1]=lb[1];
		return make_pair(yA,yB);
    }

}