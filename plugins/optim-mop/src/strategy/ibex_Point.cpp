/*
 * ibex_Point.cpp
 *
 *  Created on: 13 apr. 2021
 *      Author: iaraya
 */

#include "ibex_NDSrp.h"

namespace ibex {

        double Point::compute_m(Vector &p1, Vector &p2){
            return (p2[1]-p1[1])/(p2[0]-p1[0]);
        }

        bool Point::is_upper(){
            Vector z=*next; Vector y=*this; Vector x=*prev;

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

        

        double Point::compute_area(Vector& x, Vector& y, Vector& z){
            return std::abs((x[0]*(y[1]-z[1]) + y[0]*(z[1]-x[1]) + z[0]*(x[1]-y[1]))/2);

        }

        Point lineIntersection(double a1, double a2, double c1, double c2)
        {
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

        /*
        compute how much the hypervolume is reduced if the point is removed
        */
        double Point::compute_hv_contribution(){
            if( (*this)[0]==NEG_INFINITY || (*this)[0]==POS_INFINITY || 
                    (*this)[1]==NEG_INFINITY || (*this)[1]==POS_INFINITY) return POS_INFINITY;
            
            if(is_upper()) up_point=true;
            else up_point=false;

            //si el punto queda por debajo del segmento entra aca
            if (!up_point)
                return compute_area(*prev,*this,*next);
            else{
                Vector &c=*this, &p=*prev, &n=*next;
                //si el punto es el penultimo entra aca retorna POS_INFINITY (para mantener cosas simples)
                if (this->next->next->next == NULL){
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

                    /*A.x = B?
                    A = np.array([[-m1,1],[-m2,1]])
                    B = np.array([b1,b2])
                
                    if np.linalg.det(A) == 0:
                        x,y,z,w = np.linalg.lstsq(A,B)
                    else:
                        x = np.linalg.solve(A,B)
                    */
                    Point x=lineIntersection(-m1,b1,-m2,b2);
                    //pointIntersection(const Vector& p0, const Vector& p1, const Vector& p2, const Vector& p3)

                    return compute_area(p,x,c);
                }

                //No debería ir antes?
                if(!next->is_upper()) return POS_INFINITY;
                
                Vector &nn =*next->next;

                double m1 = compute_m(p,c);
                double m2 = compute_m(n,nn);

                //si la segunda pendiente es vertical entra aca
                if (m2 == POS_INFINITY || m2 == NEG_INFINITY){
                    Vector x(2);
                    x[0] = n[0];
                    double b1 = c[1] - (m1*c[0]);
                    x[1] = (m1*x[0]) + b1;
                    return compute_area(c,x,n);
                }

                //si la primera pendiente es horizontal entra aca
                if (m1 == 0){
                    Vector x(2);
                    x[1] = c[1];
                    double b2 = n[1] - (m2*n[0]);
                    x[0] = (x[1]-b2)/m2;
                    return compute_area(c,x,n);
                }

                //si no se cumplen las 2 anteriores continua por aca
                double b1 = c[1] - (m1*c[0]);
                double b2 = n[1] - (m2*n[0]);

                Point x=lineIntersection(-m1,-m2,b1,b2);
                cout << c << x << n << endl;
                return compute_area(c,x,n);
            }

        }
}