/*
 * ibex_Point.cpp
 *
 *  Created on: 13 apr. 2021
 *      Author: iaraya
 */

#include "ibex_NDShv.h"

namespace ibex {

        double NDShv::min_hv = 0.0;

        double Point::compute_m(const Vector &p1, const Vector &p2){
            return (p2[1]-p1[1])/(p2[0]-p1[0]);
        }

        bool Point::is_upper() const{
            return is_upper(*this, *prev, *next);
        }

        bool Point::is_upper(const Vector& y, const Vector&x, const Vector& z){
            double m = compute_m(x,z);
            double b = x[1] - (m*x[0]);

            double valor = (m*y[0]) + b;

            //si valor es mayor, quiere decir que el punto que se evalua queda por debajo
            if (valor>y[1])
                return false;

            //por el contrario, si es menor quiere decir que el punto queda por arriba
            else 
                return true;
        }

        

        double Point::compute_area(const Vector& x, const Vector& y, const Vector& z){
            return std::abs((x[0]*(y[1]-z[1]) + y[0]*(z[1]-x[1]) + z[0]*(x[1]-y[1]))/2);

        }

        Point lineIntersection(double a1, double a2, double c1, double c2)
        {
            return Point(POS_INFINITY, POS_INFINITY); //remove

            // Line AB represented as a1x + b1y = c1
            double b1 = 1;
        
            // Line CD represented as a2x + b2y = c2
            double b2 = 1;
        
            double determinant = a1*b2 - a2*b1;
            

            if (determinant == 0)
            {
                // The lines are parallel. This is simplified
                // by returning a pair of FLT_MAX
                return Point(POS_INFINITY, POS_INFINITY);
            }
            else
            {
                double x = (b2*c1 - b1*c2)/determinant;
                double y = (a1*c2 - a2*c1)/determinant;
                return Point(x, y);
            }
        }

        double Point::get_hv_contribution(const Vector& p, const Vector& prev, const Vector& next){
            if( p[0]==NEG_INFINITY || p[0]==POS_INFINITY || 
                    p[1]==NEG_INFINITY || p[1]==POS_INFINITY  ||
                    prev[1]==POS_INFINITY || next[0]==POS_INFINITY ) return POS_INFINITY;
            
            //if (next.next->next == NULL) return POS_INFINITY;
            if(!is_upper(p,prev,next)) return compute_area(prev,p,next);
            else return POS_INFINITY;
        }


        /*
        compute how much the hypervolume is reduced if the point is removed
        */
        double Point::_compute_hv_contribution(Vector& tmp_next){
            if( (*this)[0]==NEG_INFINITY || (*this)[0]==POS_INFINITY || 
                    (*this)[1]==NEG_INFINITY || (*this)[1]==POS_INFINITY  ||
                    (*prev)[1]==POS_INFINITY || (*next)[0]==POS_INFINITY ) return POS_INFINITY;
            
            if (this->next->next->next == NULL) return POS_INFINITY;
            
            if(is_upper()) up_point=true;
            else up_point=false;

            //si el punto queda por debajo del segmento entra aca
            if (!up_point)
                return compute_area(*prev,*this,*next);
            else{
                return POS_INFINITY;
                Vector &c=*this, &p=*prev, &n=*next;
                //si el punto es el penultimo entra aca retorna POS_INFINITY (para mantener cosas simples)
                /*if (this->next->next->next == NULL){
                    return POS_INFINITY;
                    Vector& pp=*prev->prev;
                    
                    double m1 = compute_m(pp,p);
                    double m2 = compute_m(c,n);

                    //si la recta del final es vertical entra acá
                    if (m2 == NEG_INFINITY) {
                        double b1 = p[1] - (m1*p[0]);
                        Vector x={c[0],m1*c[0] + b1};
                        return compute_area(p,x,c);
                    }

                    //si la pendiente del antepenultimo segmento es 0 entra aca 
                    if (m1 == 0){
                       double b2 = c[1] - (m2*c[0]);
                       Vector x={(p[1]-b2)/m2,p[1]};
                       return compute_area(p,x,c);
                    }

                    //si no es ninguno de los 2 casos calcula normalmente la interseccion entre los 2 segmentos
                    //b1 = y2[pto-1] - (m1*y1[pto-1])
                    //b2 = y2[pto+1] - (m2*y1[pto+1])
                    double b1 = p[1] - (m1*p[0]);
                    double b2 = n[1] - (m2*n[0]);

                    Point x=lineIntersection(-m1,b1,-m2,b2);
                    
                    return compute_area(p,x,c);
                }*/
                
                if(!next->is_upper()) return POS_INFINITY;
                
                Vector &nn =*next->next;

                double m1 = compute_m(p,c);
                double m2 = compute_m(n,nn);

                //si la segunda pendiente es vertical entra aca
                if (m2 == POS_INFINITY || m2 == NEG_INFINITY){
                    Vector x(2);
                    x[0] = n[0];
                    double b1 = c[1] - (m1*c[0]);
                    x[1] = ((Interval(m1)*Interval(x[0])) + Interval(b1)).ub();
                    tmp_next=Vector(x);
                    return compute_area(c,x,n);
                }

                //si la primera pendiente es horizontal entra aca
                if (m1 == 0){
                    Vector x(2);
                    x[1] = c[1];
                    double b2 = n[1] - (m2*n[0]);
                    x[0] = ((Interval(x[1])-Interval(b2))/Interval(m2)).ub();
                    tmp_next=Vector(x);
                    return compute_area(c,x,n);
                }

                //si no se cumplen las 2 anteriores continua por aca
                double b1 = c[1] - (m1*c[0]);
                double b2 = n[1] - (m2*n[0]);

                Point x=lineIntersection(-m1,-m2,b1,b2);
                tmp_next=Vector(x);
                return compute_area(c,x,n);
            }

        }

    NDShv::NDShv() : insertions(0){
        Point* first=new Point(NEG_INFINITY,POS_INFINITY);
		Point* middle=new Point(POS_INFINITY,POS_INFINITY);
		Point* last=new Point(POS_INFINITY,NEG_INFINITY);
		middle->push_before(last);
		first->push_before(middle);

        set<Vector*, sort_rp2>::insert(last);
		set<Vector*, sort_rp2>::insert(middle);
		set<Vector*, sort_rp2>::insert(first);


		set<Point*, sort_hv>::insert(last);
		set<Point*, sort_hv>::insert(middle);
		set<Point*, sort_hv>::insert(first);
    }

    void NDShv::insert(const Vector& p){
        if(set<Vector*, sort_rp2>::find((Vector*)&p) != set<Vector*, sort_rp2>::end() ) return;
        

        set<Vector*>::iterator next = set< Vector*, sort_rp2 >::lower_bound((Vector*) &p);
        if(NDShv::min_hv>0.0 && Point::get_hv_contribution(*((Vector*) &p),*((Point*)*next)->prev,**next) < NDShv::min_hv) return;

        Point* pp = new Point(p);
        pp->push_before((Point*)*next); 
        update((Point*)*next);      
        update(pp->prev);
        pp->compute_hv_contribution();

        set<Vector*, sort_rp2>::insert(pp);
        set<Point*, sort_hv>::insert(pp);
        if( set< Vector*, sort_rp2 >::size() != set<Point*, sort_hv>::size()) exit(0);

        insertions++;

        if(insertions%10000==0){
            cout << size() << "-->";
            while(front()->hv_contribution < NDShv::min_hv) pop_front();
            cout << size() << endl;
        }


    }

    void NDShv::erase(std::set<Vector*>::iterator it){
        Point* p=(Point*) *it;
        erase(p);
    }

    void NDShv::erase(Point* p){       
        p->next->prev=p->prev;
        p->prev->next=p->next;
        update(p->prev);
        update(p->next);

        set< Vector*, sort_rp2 >::erase(p);
        set<Point*, sort_hv>::erase(p);

        if( set< Vector*, sort_rp2 >::size() != set<Point*, sort_hv>::size()) exit(0);

        delete p;
    }

    void NDShv::erase(Vector* p){
        Point* pp = (Point*) (*set< Vector*, sort_rp2 >::find(p));
        erase(pp);
    }

    void NDShv::update(Point* p){

        if((*p)[0] == POS_INFINITY || (*p)[0] == NEG_INFINITY) return;
        if((*p)[1] == POS_INFINITY || (*p)[1] == NEG_INFINITY) return;

        set<Point*, sort_hv>::erase(p);
        p->compute_hv_contribution();
        set<Point*, sort_hv>::insert(p);

        
        if( set< Vector*, sort_rp2 >::size() != set<Point*, sort_hv>::size()) exit(0);
 
    }

   void NDShv::update(Point* p, Vector& v){
       if(set< Vector*, sort_rp2 >::find(&v) != set< Vector*, sort_rp2 >::end()) return;
       

        set< Vector*, sort_rp2 >::erase(p);
        (*p)[0]=v[0];
        (*p)[1]=v[1];

        set< Vector*, sort_rp2 >::insert(p);
        if( set< Vector*, sort_rp2 >::size() != set<Point*, sort_hv>::size()){
            exit(0);
        }
 
    }

    

    Point* NDShv::front(){
        return *set<Point*, sort_hv>::begin();
    }

    double NDShv::pop_front(){

        Point* p=NDShv::front();
        double r=p->hv_contribution;
        force_erase(p);
        return r;
    }

    //removing only down_points
    void NDShv::force_erase(Point* p){
        //Vector next(1);
        //p->_compute_hv_contribution(next);

        
        //if(next.size()==2)
          //  update(p->next,next);
        

        erase(p);

    }

}