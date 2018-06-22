/*
 * ibex_PFunction.cpp
 *
 *  Created on: 23 may. 2018
 *      Author: iaraya
 */

#include "ibex_PFunction.h"
#include "ibex_OptimizerMOP.h"

namespace ibex {

PFunction::PFunction(const Function& f1, const Function& f2, const IntervalVector& xa, const IntervalVector& xb):
		f1(f1),f2(f2), xa(xa), xb(xb) {}

Interval PFunction::eval(const Interval& t, const Interval& m, bool minimize) const{
	IntervalVector xt = xa+t*(xb-xa);
	Interval result;
	if(m.is_empty()) {
		result = -OptimizerMOP::eval_goal(f1,xt,  xt.size());
	} else result = OptimizerMOP::eval_goal(f2,xt, xt.size()) - m*OptimizerMOP::eval_goal(f1,xt,  xt.size());
	if( (!minimize && m.ub() > 0) || (minimize && m.ub() <= 0) || (!minimize && m.is_empty()) ) return Interval(-result.ub(), -result.lb());
	else return result;
}

Interval PFunction::deriv(const Interval& t, const Interval& m, bool minimize) const{
	IntervalVector xt = xa+t*(xb-xa);
	IntervalVector g1 = OptimizerMOP::deriv_goal(f1, xt, xt.size());
	IntervalVector g2 = OptimizerMOP::deriv_goal(f2, xt, xt.size());
	Interval result;
	if(m.is_empty()) {
		return (-g1)*(xb-xa);
	} else return (g2-m*g1)*(xb-xa);
}

IntervalVector PFunction::get_point(const Interval& t) const{
	IntervalVector y(2);
	IntervalVector xt = xa+t*(xb-xa);
	y[0]=OptimizerMOP::eval_goal(f1,xt,  xt.size());
	y[1]=OptimizerMOP::eval_goal(f2,xt,  xt.size());
	return y;
}



void PFunction::get_curve_y(std::vector< pair <double, double> >& curve_y ){
	Interval a;
	IntervalVector interVector = IntervalVector(2);
	double value;
	int max_iterations = 100;
	for(int i=0;i <= max_iterations; i++) {
		a = Interval((double) i/max_iterations);
		//cout << a << endl;
		interVector = get_point(a);
		//cout << interVector[0].ub() << "," <<  interVector[1].ub() << endl;
		curve_y.push_back(make_pair(interVector[0].ub(), interVector[1].ub()));
	}
}


/**
 * \brief minimize/maximize the function pf: f1(t)+w*f2(t)
 * returning the best solution found t and its the lb/ub of its evaluation
 * input m, minimize, max_c=max_value
 */
pair<double, double> PFunction::optimize(const Interval& m, bool minimize, double max_c){
	//TODO: we are assuming minimize=false;

	// TODO: ver cuando es conveniente realizar Newton
	// 1. Revisar que cuando la derivada sea vacia no se contracte
	// 2. Revisar que cunado la derivada sea muy alta no haga nada
	// 3. Revisar que newton si supera el max c no siga buscando y no se contracte
	// 4. Distancia minima entre puntos ya e yb

	if(minimize && max_c != POS_INFINITY) max_c = -max_c;

	Interval derivate = deriv(Interval(0,1), m, minimize);

	double t_final;
	double epsilon = 0.0003;
	double lb = NEG_INFINITY;
	stack<Interval> pila;
	pila.push(Interval(0,1));
	Interval inter, left, right;
	double point_t, point_c, t_before, error, min_interval, max_diam;
	Interval y_r, y_c, y_l;
	double lb_interval;
	// global values
	error = 1e-4;
	min_interval = 1e-4;
	max_diam = 1e-3;

	// pila
	int iter = 0;
	double t_temp;


	while(!pila.empty() and lb < max_c) {
		cout << pila.size() <<endl;

		inter = pila.top();
		pila.pop();

		iter++;

		derivate = deriv(inter, m, minimize);

		// lowerbounding
		y_r=eval(inter.lb(), m, minimize);
		y_c=eval(inter.mid(), m, minimize);
		y_l= eval(inter.ub(), m, minimize);

		lb_interval = y_r.ub();
		t_temp = inter.lb();
		if(lb_interval < y_c.ub()) {
			lb_interval = y_c.ub();
			t_temp = inter.mid();
		}
		if(lb_interval < y_l.ub()) {
			lb_interval = y_l.ub();
			t_temp = inter.ub();
		}

		if(fabs(lb_interval) < 1 && lb_interval + epsilon > lb) {
			t_final = t_temp;
			lb = lb_interval+epsilon;
		}
		else if(fabs(lb_interval) >= 1 && lb_interval + fabs(lb_interval)*epsilon > lb) {
			t_final = t_temp;
			lb = lb_interval + fabs(lb_interval)*epsilon;
		}

		// Contract if derivate is not empty
		if(!derivate.is_empty()) {
			// contract Newton from left
			point_t = inter.lb();
			point_c = eval(point_t, m, minimize).ub();
			t_before = NEG_INFINITY;
			while(derivate.ub() > 0 and point_t - t_before > error and point_t < inter.ub()) {
				t_before = point_t;

				if(0 == derivate.ub())
					point_t = POS_INFINITY;
				else{
					cout << (fabs(lb) < 1) << endl;
					point_t = (lb - point_c)/derivate.ub() + t_before;
					point_c = eval(point_t, m, minimize).ub();
				}
				// error en caso de que c sea mayor al lb+epsilon
				if(point_t < inter.ub() and point_c > lb+epsilon) {
					//cout << "ERRROR LEFT: point right is greater than lb " << endl;
					//getchar();
					break;
					//exit(-1);
				}

				//cout << "point_t " << point_t << endl;
			}

			// Se elimina el intervalo ya que no contiene una solucion mejor a lb+epsilon
			if(point_t >= inter.ub()) {
				continue;
			} else {
				//se contracta el intervalo si point_t > 0
				if(point_t > 0)
					inter = Interval(point_t, inter.ub());
			}

			//cout << "contract left inter diam " << inter.diam() << " " << inter << endl;

			// contract Newton from right
			point_t = inter.ub();
			point_c = eval(point_t, m, minimize).ub();
			t_before = NEG_INFINITY;


			while(derivate.lb() < 0 and t_before - point_t > error and point_t > inter.lb()) {
				t_before = point_t;

				if(0 == derivate.lb())
					point_t = NEG_INFINITY;
				else{
					point_t = t_before - (lb - point_c)/derivate.lb();
					point_c = eval(point_t, m, minimize).ub();
				}

				// error en caso de que c sea mayor al lb+epsilon
				if(point_t > inter.lb() and point_c > lb+epsilon) {
					//cout << "ERRROR RIGHT: point right is greater than lb" << endl;
					//getchar();
					break;
				}
			}

			// Se remueve
			if(point_t <= inter.lb()) {
				continue;
			} else {
				inter = Interval(inter.lb(), point_t);
			}

			//cout << "contract right inter diam " << inter.diam() << " " << inter << endl;
		}


		// bisect interval and push in stack
		if(inter.is_bisectable() and inter.diam() > max_diam) {
			pair<Interval,Interval> bsc = inter.bisect(0.5);
			//cout << "bsc1 diam " << bsc.first.diam() << " " << bsc.first << endl;
			//cout << "bsc2 diam " << bsc.second.diam() << " " << bsc.second << endl;
			pila.push(bsc.first);
			pila.push(bsc.second);
		}
		//cout << "iteracion " << iter << " pila " << pila.size() << endl;
	}
	//cout << "derivate check2 " << derivate << endl;



	//if (derivate.is_empty()) {
		//if(minimize) return min(eval(0, m, minimize), eval(1, m, minimize)).lb();
		//else
	//		return max(eval(0, m, minimize), eval(1, m, minimize)).ub();
	//}

	// if(minimize) lb = max(eval(0, m, minimize), eval(1, m, minimize)).lb() - (lb - max(eval(0, m, minimize), eval(1, m, minimize)).lb());

	if( (!minimize && m.ub() > 0) || (minimize && m.ub() <= 0) || (!minimize && m.is_empty()) ) return make_pair(-lb, t_final);
	else return make_pair(lb, t_final);
}

} /* namespace ibex */
