/*
 * ibex_PFunction.h
 *
 *  Created on: 23 may. 2018
 *      Author: iaraya
 */

#include "ibex_Function.h"


#ifndef OPTIM_MOP_SRC_STRATEGY_IBEX_PFUNCTION_H_
#define OPTIM_MOP_SRC_STRATEGY_IBEX_PFUNCTION_H_


using namespace std;

namespace ibex {

/**
 * Parameterized function f(t) ← f1(xt) - m*f2(xt)
 * xt = xa + t*(xb-xa)
 */
class PFunction{

public:
	PFunction(const Function& f1, const Function& f2, const Interval& m, const IntervalVector& xa, const IntervalVector& xb);

	Interval eval(const Interval& t) const;

	Interval deriv(const Interval& t) const;

	IntervalVector get_point(const Interval& t) const;

	/**
	 * \brief gets the image of the segment line xa-xb
	 */
	void get_curve_y(std::vector< pair <double, double> >& curve_y );

	/**
	 * \brief minimize/maximize the function pf: f1(t)+w*f2(t)
	 * returning the best solution found t and its the lb/ub of its evaluation
	 */
	double optimize(Interval max_c, bool minimize=false);

private:

	const Function& f1;
	const Function& f2;
	Interval m;
	IntervalVector xa;
	IntervalVector xb;
};

} /* namespace ibex */

#endif /* OPTIM_MOP_SRC_STRATEGY_IBEX_PFUNCTION_H_ */
