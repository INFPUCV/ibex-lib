/*
 * ibex_PFunction.cpp
 *
 *  Created on: 23 may. 2018
 *      Author: iaraya
 */

#include "ibex_PFunction.h"
#include "ibex_OptimizerMOP.h"

namespace ibex {

bool PFunction::MIN=true;
bool PFunction::MAX=false;

double PFunction::_min_newton_step=0.25; //25% of the current diameter
double PFunction::_min_diam=0.15; //15% of the initial diameter
double PFunction::_eps_opt=1e-7;

PFunction::PFunction(const Function& f1, const Function& f2, const IntervalVector& xa, const IntervalVector& xb):
		f1(f1),f2(f2), xa(xa), xb(xb) {}

/**
 * convert pf.t to t in inter
 * t = inter.lb() + pf.t*(inter.ub() - inter.lb());
 */
void PFunction::contract_curve(const Interval& t) {
	if(t.is_empty() || t.lb() < 0 || t. ub() > 1) return;

	IntervalVector a = xa+t.lb()*(xb-xa);
	IntervalVector b = xa+t.ub()*(xb-xa);

	xa = a;
	xb = b;
}

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

bool PFunction::newton_rcontract(const Interval& m, bool minimize, Interval& inter, const Interval& derivate, double lb){
	double point_t = inter.ub();
	double point_c = eval(point_t, m, minimize).ub();
	double t_before = NEG_INFINITY;
	while(derivate.lb() < 0 && t_before - point_t > _min_newton_step*inter.diam() && point_t > inter.lb()) {

		t_before = point_t;

		point_t = t_before - (lb - point_c)/derivate.lb();
		point_c = eval(point_t, m, minimize).ub();
	}

	// TODO: puede pasar esto? error en caso de que c sea mayor al lb+epsilon
	if(point_t > inter.lb() and point_c > lb + _eps_opt){
		cout << "error: point_c > lb+epsilon" << endl;
		exit(0);
	}

	// Se elimina el intervalo ya que no contiene una solucion mejor a lb+epsilon
	if(point_t <= inter.lb())
		return false;
	 else
		inter = Interval(inter.lb(), point_t);


	return true;
}

bool PFunction::newton_lcontract(const Interval& m, bool minimize, Interval& inter, const Interval& derivate, double lb){
	double point_t = inter.lb();
	double point_c = eval(point_t, m, minimize).ub();
	double t_before = NEG_INFINITY;

	while(derivate.ub() > 0 && point_t - t_before > _min_newton_step*inter.diam() && point_t < inter.ub()) {
		t_before = point_t;

		point_t = (lb - point_c)/derivate.ub() + t_before;
		point_c = eval(point_t, m, minimize).ub();
	}

	// TODO: puede pasar esto? error en caso de que c sea mayor al lb+epsilon
	if(point_t < inter.ub() and point_c > lb+_eps_opt){
		cout << "error: point_c > lb+eps" << endl;
		exit(0);
	}

	// Se elimina el intervalo ya que no contiene una solucion mejor a lb+epsilon
	if(point_t >= inter.ub())
		return false;
	 else {
		//se contracta el intervalo si point_t > 0
		inter = Interval(point_t, inter.ub());
	}

	return true;
}

/**
 * \brief minimize/maximize the function pf: f1(t)+w*f2(t)
 * returning the best solution found t and its the lb/ub of its evaluation
 * input m, minimize, max_c=max_value
 */
pair<double, double> PFunction::optimize(const Interval& m, bool minimize, double max_c, Interval init){
	double diam0=init.diam();

	if(init.is_empty()) init=Interval(0,1);
	//TODO: we are assuming minimize=false;

	// TODO: ver cuando es conveniente realizar Newton
	// 1. Revisar que cuando la derivada sea vacia no se contracte
	// 2. Revisar que cunado la derivada sea muy alta no haga nada
	// 3. Revisar que newton si supera el max c no siga buscando y no se contracte
	// 4. Distancia minima entre puntos ya e yb

	//if(minimize && max_c != POS_INFINITY) max_c = -max_c;

	Interval derivate = deriv(init, m, minimize);

	double t_final;
	double lb = NEG_INFINITY;
	stack<Interval> pila;
	pila.push(init);
	Interval inter, left, right;
	double point_t, point_c, t_before;
	Interval y_r, y_c, y_l;
	double lb_interval;


	// pila
	double t_temp;

	int iter=0;
	while(!pila.empty() and lb < max_c) {
		// cout << pila.size() <<endl;

		inter = pila.top();
		pila.pop();

		derivate = deriv(inter, m, minimize);

		/************ lowerbounding ***************/
		//TODO: no se podrá optimizar esto?
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

		if(fabs(lb_interval) < 1 && lb_interval + _eps_opt > lb) {
			t_final = t_temp;
			lb = lb_interval+_eps_opt;
		}
		else if(fabs(lb_interval) >= 1 && lb_interval + fabs(lb_interval)*_eps_opt > lb) {
			t_final = t_temp;
			lb = lb_interval + fabs(lb_interval)*_eps_opt;
		}
		/***************************************/

		// Contract if derivate is not empty
		if(!derivate.is_empty()) {

			// contract Newton from left
			if(!newton_lcontract(m,  minimize, inter, derivate, lb))
				continue;

			if(!newton_rcontract(m,  minimize, inter, derivate, lb))
				continue;

		}


		// bisect interval and push in stack
		if(inter.is_bisectable() and inter.diam() > _min_diam*diam0) {
			pair<Interval,Interval> bsc = inter.bisect(0.5);
			pila.push(bsc.first);
			pila.push(bsc.second);
		}
		iter++;
	}
	//cout << iter << endl;

	if( (!minimize && m.ub() > 0) || (minimize && m.ub() <= 0) || (minimize && m.is_empty()) ) return make_pair(-lb, t_final);
	else return make_pair(lb, t_final);
}

} /* namespace ibex */
