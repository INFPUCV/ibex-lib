//============================================================================
//                                  I B E X
// File        : ibex_Optimizer.h
// Author      : Matias Campusano, Damir Aliquintui, Ignacio Araya
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Sep 24, 2017
// Last Update : Sep 24, 2017
//============================================================================

#ifndef __IBEX_OPTIMIZERMOP_H__
#define __IBEX_OPTIMIZERMOP_H__

#include "ibex_Ctc.h"
#include "ibex_Bsc.h"
#include "ibex_LoupFinderMOP.h"
#include "ibex_CellMOP.h"
#include "ibex_CtcKhunTucker.h"
#include "ibex_DistanceSortedCellBufferMOP.h"
#include "ibex_pyPlotter.h"

#include <set>
#include <map>
#include <list>
#include <stack>
#include <math.h>
//#include "ibex_DistanceSorted.h"

using namespace std;
namespace ibex {

/**
 * comparation function for sorting NDS by decreasing y
 */
struct sorty{
	bool operator()(const pair<double,double> p1, const pair<double,double> p2){
		return p1.second>p2.second;
	}
};


/**
 * comparation function for sorting NDS2 by increasing x and decreasing by y
 */
struct sorty2{
	bool operator()(const pair<double,double> p1, const pair<double,double> p2){
		if(p1.first != p2.first)
			return p1.first<p2.first;
		return p1.second>p2.second;

	}
};

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

private:

	const Function& f1;
	const Function& f2;
	Interval m;
	IntervalVector xa;
	IntervalVector xb;
};

class Node_t{
public:
	Node_t(Interval t, Interval ft) : t(t), ft(ft) { }

	friend bool operator>(Node_t& n1, Node_t& n2){
		return n1.ft.ub() > n2.ft.ub();
	}

	Interval t;
	Interval ft;
};


/**
 * \brief Global biObjetive Optimizer (ibexMOP).
 *
 * This class is an implementation of a global optimization algorithm for biObjective problems
 * described in https://github.com/INFPUCV/ibex-lib/tree/master/plugins/optim-mop by Araya et al.
 *
 * \remark In all the comments of this class, "NDS" means "Non Dominated Set"
 * related to the objectives f1 and f2
 */
class OptimizerMOP {

public:

	/**
	 * \brief Return status of the optimizer
	 *
	 * See comments for optimize(...) below.
	 */
	typedef enum {SUCCESS, INFEASIBLE, NO_FEASIBLE_FOUND, UNBOUNDED_OBJ, TIME_OUT, UNREACHED_PREC} Status;

	typedef enum {POINTS, SEGMENTS} Mode;

	/**
	 *  \brief Create an optimizer.
	 *
	 * Inputs:
	 *   \param n        - number of variables of the <b>original system</b>
	 *   \param f1	     - the objective function f1
     *	 \param f2       - the objective function f2
	 *   \param ctc      - contractor for <b>extended<b> boxes (of size n+2)
	 *   \param bsc      - bisector for <b>extended<b> boxes (of size n+2)
	 *   \param buffer   - buffer for <b>extended<b> boxes (of size n+2)
	 *   \param finder   - the finder of ub solutions
	 *   \param eps	     - the required precision
	 *

	 *
	 * \warning The optimizer relies on the contractor \a ctc to contract the domain of the goal variable
	 *          and increase the uplo. If this contractor never contracts this goal variable,
	 *          the optimizer will only rely on the evaluation of f and will be very slow.
	 *
	 * We are assuming that the objective variables are n and n+1
	 *
	 */
	OptimizerMOP(int n, const Function &f1,  const Function &f2,
			Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder, Mode nds_mode=POINTS, double eps=default_eps);

	/**
	 * \brief Delete *this.
	 */
	virtual ~OptimizerMOP();

	/**
	 * \brief Run the optimization.
	 *
	 * \param init_box             The initial box
	 *
	 * \return SUCCESS             if the global minimum (with respect to the precision required) has been found.
	 *                             In particular, at least one feasible point has been found, less than obj_init_bound, and in the
	 *                             time limit.
	 **
	 *         TIMEOUT             if time is out.
	 *
	 */
	Status optimize(const IntervalVector& init_box);

	/* =========================== Output ============================= */

	/**
	 * \brief Displays on standard output a report of the last call to optimize(...).
	 *
	 * Information provided:
	 * <ul><li> total running time
	 *     <li> totl number of cells (boxes)
	 *     <li> total number of best non-dominated solutions found
	 *     <li> the set of non-dominated solutions found
	 * </ul>
	 */
	void report(bool verbose=true);


	/**
	 * \brief Get the status.
	 *
	 * \return the status of last call to optimize(...).
	 */
	Status get_status() const;

	/**
	 * \brief Get the "UB" set of the pareto front.
	 *
	 * \return the UB of the last call to optimize(...).
	 */
	map< pair <double, double>, IntervalVector >& get_UB()  { return NDS; }

	//std::set< point2 >& get_LB()  { return LB; }

	/**
	 * \brief Get the time spent.
	 *
	 * \return the total CPU time of last call to optimize(...)
	 */
	double get_time() const;

	/**
	 * \brief Get the number of cells.
	 *
	 * \return the number of cells generated by the last call to optimize(...).
	 */
	double get_nb_cells() const;

	/**
	 * \brief returns the distance from the box to the current NDS
	 */
	static double distance2(const Cell* c);

	/* =========================== Settings ============================= */

	/**
	 * \brief Number of variables.
	 */
	const int n;

	/**
	 * \brief Objective functions
	 * Functions have the form: f1 - z1  and f2 - z2. Thus, in order to
	 * evaluate them we have to set z1 and z2 to [0,0].
	 */
	const Function& goal1;
	const Function& goal2;

	/**
	 * \brief Contractor for the extended system.
	 *
	 * The extended system:
	 * (y=f(x), g_1(x)<=0,...,g_m(x)<=0).
	 */
	Ctc& ctc;

	/**
	 * \brief Bisector.
	 *
	 * Must work on extended boxes.
	 */
	Bsc& bsc;

	/**
	 * Cell buffer.
	 */
	CellBuffer& buffer;

	/**
	 * \brief LoupFinder
	 */
	LoupFinderMOP& finder;

	/** Required precision for the envelope */
	double eps;


	/** Default precision: 0.01 */
	static const double default_eps;

	/**
	 * \brief Trace activation flag.
	 */
	int trace;

	/**
	 * \brief Time limit.
	 *
	 * Maximum CPU time used by the strategy.
	 * This parameter allows to bound time consumption.
	 * The value can be fixed by the user.
	 */
	double timeout;


	/* ======== Some other parameters of the solver =========== */

	//if true: Save a file to be plotted by plot.py (default value= false).
	bool static _plot;

	//Min distance between two non dominated points to be considered, expressed as a fraction of eps (default value= 0.1)
	double static _min_ub_dist;

	//True: the solver uses the upper envelope of the cy contract for contraction
	static bool _cy_upper;

	//True: the solver uses the lower envelope of the cy contract in contraction
	static bool cy_contract_var;

	//True: the solver reduces the search spaces by reducing the NDS vectors in (eps, eps)
	static bool _eps_contract;

	//NDS mode: POINTS or SEGMENTS
	Mode nds_mode;

	/**
	 * \brief Evaluate the goal in the point x
	 */
	static Interval eval_goal(const Function& goal, const IntervalVector& x, int n);

	/**
	 * \brief Gradient of the goal in the point x
	 */
	static IntervalVector deriv_goal(const Function& goal, const IntervalVector& x, int n);

protected:
	/**
	 * The contraction using y+cy
	 */
	void cy_contract(Cell& c);

	void cy_contract2(Cell& c, list <pair <double,double> >& inpoints);

	/**
	 * \brief return a set of non-dominated segments of the box
	 */
	list<pair <double,double> > non_dominated_segments(IntervalVector& box);

	double distance22(const Cell* c);

	/**
	 * \brief Contract and bound procedure for processing a box.
	 *
	 * <ul>
	 * <li> contract the cell using #dominance_peeler() and #discard_generalized_monotonicty_test(),
	 * <li> contract with the contractor ctc,
	 * </ul>
	 *
	 */
	void contract_and_bound(Cell& c, const IntervalVector& init_box);

	/**
	 * \brief The box is reduced using the NDS
	 *
	 * Details are given in [Martin, B. et al.
	 * Constraint propagation using dominance in interval
	 * Branch & Bound for nonlinear biobjective optimization (2017)]
	 */
	void dominance_peeler(IntervalVector& box);

	/**
	 * \brief The box is reduced using NDS2
	 *
	 * Adaptation of dominance_peeler to NDS2
	 * Returns the set of non-dominated segments in box
	 */
	void dominance_peeler2(IntervalVector &box, list <pair <double,double> >& inpoints);

    /**
     *  \brief returns true if the facet orthogonal to the i direction of the box is feasible.
     *
     *  See #discard_generalized_monotonicty_test()
     */
    bool is_inner_facet(IntervalVector box, int i, Interval bound){
    	box.resize(n);
    	box[i]=bound;
    	return finder.norm_sys.is_inner(box);
    }

    /**
     * \brief return true is the pair is dominated by some NDS point, false otherwise
     */
    bool is_dominated(pair< double, double>& eval);

    /**
     * \brief Implements the method for discarding boxes proposed in [J. Fernandez and B. Toth,
     * "Obtaining the efficient set of nonlinear biobjective optimiza-tion  problems  via  interval
     * branch-and-bound  methods" (2009)]
     */
	void discard_generalized_monotonicty_test(IntervalVector& box, const IntervalVector& initbox);


	/**
	 * \brief Main procedure for updating the NDS.
	 * <ul>
	 * <li> finds two points in a polytope using the #LoupFinderMOP,
	 * <li> generate n equi-distant points between the two points,
	 * <li> correct the points using a Hansen feasibility test (See #PdcHansenFeasibility)
	 * <li> add the non-dominated vectors to NDS
	 * </ul>
	 */
	bool update_NDS(const IntervalVector& box);


	/**
	 * \brief Main procedure for updating the NDS.
	 * <ul>
	 * <li> finds two points xa and xb in a polytope using the #LoupFinderMOP,
	 * <li> finds a segment passing over the points (f1(x),f2(x)) in the segment xa-xb,
	 * <li> add the segment to NDS
	 * </ul>
	 */
	bool update_NDS2(const IntervalVector& box);


	/**
	 * \brief minimize/maximize the function pf: f1(t)+w*f2(t)
	 * returning the best solution found t and its the lb/ub of its evaluation
	 */
	pair<double, double> optimize_pf(PFunction& pf, bool minimize=false){

		if(minimize){
			cout << "minimize f1(t)+w*f2(t) is not implemented yet!" << endl;
			exit(0);
		}
		//deep-first search
		stack<Node_t> nodes;
		//std::priority_queue<Interval, std::vector<Interval> > nodes;

		double LB=NEG_INFINITY, eps=0.01;
		double UB=NEG_INFINITY;
		double best_t=0.0;

		Interval t=Interval(0.0,1.0);

		nodes.push(Node_t(t,pf.eval(t)));

		int i=0;
		while(!nodes.empty()){
			i++;
			Node_t n= nodes.top(); nodes.pop();

			//the nodes is removed if its ub is lower than LB+eps
			if(n.ft.ub() < LB+eps) {UB=std::max(UB, n.ft.ub()); continue;}

			//the size of the node is too small, then the UB is updated
			if(n.t.diam() < 0.01) {
				if(n.ft.ub() > UB) UB=n.ft.ub();
				continue;
			}

			//we search for a better solution in the midpoint
			double probing = pf.eval(n.t.mid()).ub();
			if(probing > LB) {
				best_t=n.t.mid();
				LB=probing;
			}


			//newton step
			Interval d(pf.deriv(n.t));

			while(true){
				Interval y0 = pf.eval(n.t.lb());
				if(y0.is_empty()) break;
				if(y0.ub() > LB) {LB=y0.ub(); best_t=n.t.lb(); break;}

				if(d.ub()==0.0) break;
				Interval x= (Interval(LB) - y0)/Interval(d.ub());

				//contract t
				if(x.lb()>eps)
					n.t=Interval((n.t.lb()+x).lb(),n.t.ub());
				else break;
			}

			while(true){
				cout << "hola mundo este e un comentario para realizar pruebas " << nodes.size() << endl;
				Interval y0 = pf.eval(n.t.ub());
				if(y0.is_empty()) break;
				if(y0.ub() > LB) {LB=y0.ub(); best_t=n.t.ub(); break;}

				if(d.lb()==0.0) break;
				Interval x= (Interval(LB) - y0)/Interval(-d.lb());

				//contract t
				if(x.ub()>eps)
					n.t=Interval(n.t.lb(),(n.t.ub()-x).ub());
				else break;
			}

			//bisection
			Interval tl = Interval(n.t.lb(), n.t.mid());
			Interval tr = Interval(n.t.mid(), n.t.ub());

			nodes.push(Node_t(Interval(tl),pf.eval(tl)));
			nodes.push(Node_t(Interval(tr),pf.eval(tr)));
		}

		cout << i << endl;
		return make_pair(best_t, UB);

	}

	void addVectortoNDS(pair< double, double> eval1, pair< double, double> eval2) {
		cout << "addVectortoNDS-------" << endl;
		std::map<pair<double, double>, IntervalVector>::iterator aux, it1 = --NDS2.lower_bound(eval1);
		pair< double, double> first, second, point, inter_last = make_pair(NEG_INFINITY, NEG_INFINITY);
		std::map< pair <double, double>, IntervalVector, sorty2 > DS2;
		std::map< pair <double, double>, IntervalVector, sorty2 > NDS_points;

		DS2.insert(*it1);
		it1++;

		bool flagDS2 = false;
		for(;it1 != NDS2.end();) {
			if(it1->first.second < eval1.second and it1->first.second < eval2.second) flagDS2= true;
			aux = it1;
			++aux;
			DS2.insert(*it1);
			NDS2.erase(it1);
			it1 = aux;
			if(flagDS2) break;
		}

		cout << "DS2" << endl;
		for(it1 = DS2.begin();it1!=DS2.end();++it1) {
			cout << "(" << it1->first.first << "," << it1->first.second << ") " << endl;
		}
		cout << "end DS2" << endl;

		it1 = DS2.begin();
		first = it1->first;
		it1++;

		// false cuando la recta pasa por fuera.
		/*
		bool flag = false;

		for(;it1 != DS2.end(); ++it1) {
			second = it1->first;
			cout << "(" << first.first << "," << first.second << ") ";
			cout << "(" << second.first << "," << second.second << ")" << endl;
			point = pointIntersection(first, second, eval1, eval2);

			cout << "intersection (" << point.first << "," << point.second << ")" << endl;

			// puntos muy cercanos son los mismos
			if(fabs(inter_last.first - point.first) > 1e-7 && fabs(inter_last.second - point.second) > 1e-7
					&& first.first <= point.first + 1e-4 && point.first - 1e-4 <= second.first  // por x
					&& second.second <= point.second + 1e-4 && point.second - 1e-4 <= first.second) { // por y
				cout << "intersection (" << point.first << "," << point.second << ")" << endl;
				inter_last = point;

				NDS_points.insert(make_pair(point,it1->second));

				if(flag) flag = false;
				else flag = true;
			}else{
				// si pasa por fuera se agregan los punto
				if(!flag) {
					NDS_points.insert(make_pair(second,it1->second));
				}
			}

			first = second;
		}
		*/
		bool flag = false;
		for(;it1 != DS2.end(); ++it1) {
			second = it1->first;
			cout << "(" << first.first << "," << first.second << ") ";
			cout << "(" << second.first << "," << second.second << ")" << endl;
			point = pointIntersection(first, second, eval1, eval2);

			cout << "intersection (" << point.first << "," << point.second << ")" << endl;

			// puntos muy cercanos son los mismos
			if(fabs(inter_last.first - point.first) > 1e-7 && fabs(inter_last.second - point.second) > 1e-7) {

				if( ((first.first <= second.first && first.first <= point.first && point.first <= second.first) ||
						(second.first <= first.first  && second.first <= point.first && point.first <= first.first)) &&
						((first.second <= second.second && first.second <= point.second && point.second <= second.second) ||
						(second.second <= first.second  && second.second <= point.second && point.second <= first.second))) {

					cout << (first.first <= second.first && first.first <= point.first && point.first <= second.first) <<
							" " << (second.first <= first.first  && second.first <= point.first && point.first <= first.first) <<
							" " << (first.second <= second.second && first.second <= point.second && point.second <= second.second) <<
							" " << (second.second <= first.second  && second.second <= point.second && point.second <= first.second) << endl;

					cout << "intersection (" << point.first << "," << point.second << ")" << endl;
					inter_last = point;

					NDS_points.insert(make_pair(point,it1->second));
				} else {
					cout << "intersection error" << endl;
				}

				//if(flag) flag = false;
				//else flag = true;

			//}else{
				// si pasa por fuera se agregan los punto
				//if(!flag) {
				//	NDS_points.insert(make_pair(second,it1->second));
				//}
			} else {
				cout << "point very close to the after point" << endl;
			}

			first = second;
		}

		cout << "------------" << endl;

		for(it1 = NDS_points.begin();it1 != NDS_points.end();++it1) {
			cout << "newpoint (" << it1->first.first << "," << it1->first.second << ")" << endl;

			if(it1->first.first != it1->first.first  || it1->first.second != it1->first.second ) getchar();
			NDS2.insert(*it1);
		}

		cout << "insert points" << endl;

		for(it1 = DS2.begin();it1!=DS2.end();++it1) {
			if( ((it1->first.first <= eval1.first && it1->first.first <= eval2.first) &&
					(it1->first.second >= eval1.second && it1->first.second >= eval2.second)) ||
				((it1->first.second <= eval1.second && it1->first.second <= eval2.second) &&
						(it1->first.first >= eval1.first && it1->first.first >= eval2.first)) ) {
				cout << "fuera" << endl;
				cout << "(" << it1->first.first << "," << it1->first.second << ") " << endl;
				NDS2.insert(*it1);
				// addPointtoNDS(make_pair(it1->first.first,it1->first.second));
			} else {
				double m = (eval1.second - eval2.second)/(eval1.first - eval2.first);
				if( (eval1.second - m*eval1.first) >= (it1->first.second - m*it1->first.first)) {
					cout << "dentro" << endl;
					cout << "(" << it1->first.first << "," << it1->first.second << ") " << endl;
					addPointtoNDS(make_pair(it1->first.first,it1->first.second));
				} else {
					cout << "eliminado" << endl;
					cout << "(" << it1->first.first << "," << it1->first.second << ") " << endl;
				}
			}
		}



	}

	void addPointtoNDS(pair< double, double> eval) {
		std::map<pair<double, double>, IntervalVector>::iterator it1 = --NDS2.lower_bound(eval);
		// std::map<pair<double, double>, IntervalVector>::iterator it1 = NDS2.begin();
		pair< double, double> point1, point2;
		point1 = it1->first;
		it1++;
		point2 = it1->first;

		cout << "punto (" << eval.first << "," << eval.second << ")" << endl;
		// Se comprueba que no sea dominado por el anterior al lower_bound
		if(point1 == eval or (point1.first <= eval.first and point1.second <= eval.second) ) {
			cout << "punto dominado ant_lower_bound por (" << point1.first
					<< "," << point1.second << ")" << endl;
			if(_plot) py_Plotter::offline_plot(NULL, NDS2);
			//getchar();
			return;
		}
		// Se comprueba que no sea dominado por el lower_bound
		//it1++;
		if(point2 == eval or (point2.first <= eval.first and point2.second <= eval.second) ) {
			cout << "punto dominado lower_bound por (" << point2.first
					<< "," << point2.second << ")" << endl;
			if(_plot) py_Plotter::offline_plot(NULL, NDS2);
			//getchar();
			return;
		}

		// comprobar que no este dominado por la recta que forma los dos puntos anteriores
		// solo se comprueba si Eval no domina a los puntos
		if(!(eval.first <= point1.first and eval.second <= point1.second ) and
				!(eval.first <= point2.first and eval.second <= point2.second)) {
			//cout << "point 1: (" << point1.first << "," << point1.second << ")" << endl;
			//cout << "point 2: (" << point2.first << "," << point2.second << ")" << endl;
			//pendiente de los dos puntos
			float m = (point2.second-point1.second)/(point2.first-point1.first);
			//cout << "pendiente = " << m << endl;
			// se obtiene el c de la funcion de los puntos
			float c = point1.second - m*point1.first;
			//cout << "c = " << c << endl;
			// se obtiene el c del nuevo punto
			float cEval = eval.second - m*eval.first;
			//cout << "cEval = " << cEval << endl;
			if(cEval > c) {
				cout << "eval dominado por recta (" << point1.first
						<< "," << point1.second << ") y (" << point2.first
						<< "," << point2.second << ")" << endl;
				if(_plot) py_Plotter::offline_plot(NULL, NDS2);
				//getchar();
				return;
			}
			cout << "eval no es dominado por recta (" << point1.first
						<< "," << point1.second << ") y (" << point2.first
						<< "," << point2.second << ")" << endl;
		} else {
			cout << "eval domina uno de los puntos (" << point1.first
						<< "," << point1.second << ") y (" << point2.first
						<< "," << point2.second << ")" << endl;
		}


		IntervalVector vec(n);

		// insertar en NDS2
		// Agregar a la lista del DS
		//- se revisa si it1 es el comienzo, si no lo es se retrocede uno
		//- se revisa si es dominado el it1, en el caso que lo sea se elimina y se guarda en el set
		//- lo anterior se realiza hasta que el eje y del it1 sea menor que el eval

		it1 = --NDS2.lower_bound(eval); // Se llega al nodo izquierdo del nodo eval
		if(it1 != NDS2.begin()) it1--; // se retrocede 1 si no es el primero
		// agrega todos los puntos dominados por el punto a DS2 y los elimina de NDS2
		std::map<pair<double, double>, IntervalVector>::iterator aux;
		std::map< pair <double, double>, IntervalVector, sorty2 > DS2;
		std::map<pair<double, double>, IntervalVector>::iterator beginit, endit;
		std::map<pair<double, double>, IntervalVector>::iterator it2 = --NDS2.lower_bound(eval);
		for(;it1 != NDS2.end();) {
			// termina cuando it1 no este dentro de los rangos dominados del punto a agregar
			if(it1->first.second < eval.second) break;
			// comprueba si esta dominado el punto para agregarlo a DS2
			if(eval.first <= it1->first.first and eval.second <= it1->first.second) {
				cout << "punto dominado (" << it1->first.first
							<< "," << it1->first.second << ")" << endl;
				aux = it1;
				++aux;
				DS2.insert(*it1);
				NDS2.erase(it1);
				it1 = aux;
			} else ++it1;
		}
		cout << "Cantidad de puntos dominados " << DS2.size() << endl;

		if(DS2.size() > 0) {

			for(it1 =NDS2.begin();it1 != NDS2.end();++it1) {
					cout << "puntos no dominados (" << it1->first.first
								<< "," << it1->first.second << ")" << endl;
				}

			for(it1 =DS2.begin();it1 != DS2.end();++it1) {
					cout << "puntos dominados (" << it1->first.first
								<< "," << it1->first.second << ")" << endl;
				}
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


		cout << "recta superior " << endl;
		cout << "it1 (" << it1->first.first << "," << it1->first.second << ")" << endl;
		cout << "beginit (" << beginit->first.first << "," << beginit->first.second << ")" << endl;
		cout << "recta inferior " << endl;
		cout << "it2 (" << it2->first.first << "," << it2->first.second << ")" << endl;
		cout << "endit (" << endit->first.first << "," << endit->first.second << ")" << endl;

		pair<double, double> eval2;

		// cout << "eval (" << eval.first << "," << eval.second << ")" << endl;
		eval2 = make_pair(eval.first, POS_INFINITY);
		pair<double, double> intersection1 = pointIntersection(it1->first, beginit->first, eval, eval2);
		cout << "intersection1 (" << intersection1.first << "," << intersection1.second << ")" << endl;
		eval2 = make_pair(POS_INFINITY, eval.second);
		pair<double, double> intersection2 = pointIntersection(it2->first, endit->first, eval, eval2);
		cout << "intersection2 (" << intersection2.first << "," << intersection2.second << ")" << endl;
		// se agregan el punto y los dos obtenidos anteriormente
		NDS2.insert(make_pair(eval, vec));
		NDS2.insert(make_pair(intersection1, vec));
		NDS2.insert(make_pair(intersection2, vec));

		//NDS2.insert(make_pair(eval, vec));
		if(_plot) py_Plotter::offline_plot(NULL, NDS2);
		cout << "addPointtoNDS puntos no dominados " << NDS2.size() << endl;
		// getchar();
	}
	//******************************************************************************************************************
	//Intersecciones y puntos
	/**
	 * @param v10 punto inicial del primer segmento
	 * @param v11 punto final del primer segmento
	 * @param v20 punto inicial del segnundo segmento
	 * @param v21 punto final del segundo segmento
	 * @return Punto o vector de interseccion entre los segmentos evaluados,
	 *         retorna null si no hay interseccion entre los segmentos
	 */
	pair<double, double> pointIntersection(
			pair<double, double> v10, pair<double, double> v11,
			pair<double, double> v20, pair<double, double> v21){
		pair<double, double> interseccion = make_pair(0.0 ,0.0);

		double m = getSlopeSegment(v10 , v11);
		//cout << "m " << m << endl;
		double n = getSlopeSegment(v20, v21);
		double c = v10.second - v10.first * m;
		double d = v20.second - v20.first * n;

		/*
		cout << "v10 (" << v10.first << "," << v10.second << ")" << endl;
		cout << "v11 (" << v11.first << "," << v11.second << ")" << endl;
		cout << "v20 (" << v20.first << "," << v20.second << ")" << endl;
		cout << "v21 (" << v21.first << "," << v21.second << ")" << endl;
		cout << v10.second << " " << v10.first << " " << m << endl;
		cout << "c " << c << " d " << d << endl;
		*/

		//Cuando el segmento es vertical su pendiente es infinito
		if(m == POS_INFINITY and n == 0) {
			//cout << "1111" << endl;
			if((v20.first <= v10.first and v10.first <= v21.first) or
					(v21.first <= v10.first and v10.first <= v20.first))
				interseccion.first = v10.first;
			if((v10.second <= v20.second and v20.second <= v11.second) or
					(v11.second <= v20.second and v20.second <= v10.second))
				interseccion.second = v20.second;
		}
		else if (m == 0 and n == POS_INFINITY){
			//cout << "22222" << endl;
			if((v10.first <= v20.first and v20.first <= v11.first) or
					(v11.first <= v20.first and v20.first <= v10.first))
				interseccion.first = v20.first;
			if((v20.second <= v10.second and v10.second <= v21.second) or
					(v21.second <= v10.second and v10.second <= v20.second))
				interseccion.second = v10.second;
		}
		else if (v11.first - v10.first == 0){
			//cout << "33333" << endl;
			interseccion.first = v11.first;
			interseccion.second = n*v11.first + d;
		}
		//Cuando el segmento es vertical su pendiente es infinito
		else if (v21.first - v20.first == 0){
			//cout << "444" << endl;
			//cout << "m " << m << " v21.first " << v21.first << endl;
			interseccion.first = v21.first;
			interseccion.second = m*v21.first + c;
		}
		else{
			//cout << "5555" << endl;
			interseccion.first = (d - c) / (m - n);
			interseccion.second = m*interseccion.first + c;
		}
		return interseccion;
	}

    /**
     * Obtiene la pendiente del segmento
     * @param first  - Vector inicial del segmento
     * @return              - Rectorna la pendiente del segmento
     */
    double getSlopeSegment(pair<double, double> first, pair<double, double> last){
        double slope1, slope2, Ffirst, Fsecond, Lfirst, Lsecond;
    	//cout << "first " << first.first << " " << first.second << endl;
    	//cout << "last " << last.first << " " << last.second << endl;
        if((last.second == NEG_INFINITY and first.second == NEG_INFINITY) or
        		(last.second == POS_INFINITY and first.second == POS_INFINITY))
        	return 0;
        if((last.first == NEG_INFINITY and first.first == NEG_INFINITY) or
        		(last.first == POS_INFINITY and first.first == POS_INFINITY))
        	return POS_INFINITY;
        //cout << (last.first == first.first) << endl;
        if(last.second == first.second) return 0;
        if(last.first == first.first) return POS_INFINITY;
        slope1 = last.second - first.second;
        //cout << slope1 << endl;
        slope2 = last.first - first.first;
        //cout << slope2 << endl;
        return slope1/slope2;

    }

	/**
	 * \brief Finds the lower segment dominated by (f1(x),f2(x)) for some point in the line xa-xb
	 */
	void dominated_segment(const IntervalVector& xa, const IntervalVector& xb){

		// TODO: ver cuando es conveniente realizar Newton
		// 1. Revisar que cuando la derivada sea vacia no se contracte
		// 2. Revisar que cunado la derivada sea muy alta no haga nada
		// 3. Revisar que newton si supera el max c no siga buscando y no se contracte
		// 4. Distancia minima entre puntos ya e yb

		// TODO: Graficar resultados del algoritmo (segmento y curva)

		Interval ya1=OptimizerMOP::eval_goal(goal1,xa,n);
		Interval ya2=OptimizerMOP::eval_goal(goal2,xa,n);
		Interval yb1=OptimizerMOP::eval_goal(goal1,xb,n);
		Interval yb2=OptimizerMOP::eval_goal(goal2,xb,n);

		cout << "found two points" << endl;
		cout << "xa: " << xa << endl;
		cout << "xb: " << xb << endl;
		cout << "ya1: " << ya1 << endl;
		cout << "ya2: " << ya2 << endl;
		cout << "yb1: " << yb1 << endl;
		cout << "yb2: " << yb2 << endl;

		Interval m = (yb2-ya2)/(yb1-ya1);
		//Interval m = (yb1-ya1)/(yb2-ya2);
		PFunction pf(goal1, goal2, m, xa, xb);

		// maximo valor de c con el punto (yb1, ya2)  de la funcion f2 = m*f1 + c
		Interval max_c, min_c, min_c2;
		max_c = ya2 - (m*yb1);
		min_c = ya2 - (m*ya1); // deberia ser lo mismo que pf.eval(1)
		min_c2 = yb2 - (m*yb1); // deberia ser lo mismo que pf.eval(1)

		cout << "m: " << m << endl;
		Interval derivate;
		//cout << "derivate (min, max): " << derivate << endl;

		double t1=0.0, t2=0.5, t3=1.0;
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
		while(!pila.empty() and lb < max_c.ub()) {
			inter = pila.top();
			pila.pop();

			iter++;
			cout << "iteracion " << iter << " pila " << pila.size() << endl;
			cout << "inter diam " << inter.diam() << " " << inter << endl;

			derivate = pf.deriv(inter);
			cout << "derivate " << derivate << endl;
			if (derivate.is_empty()) break;

			// lowerbounding
			y_r=pf.eval(inter.lb());
			y_c=pf.eval(inter.mid());
			y_l= pf.eval(inter.ub());

			lb_interval = max(max(y_r, y_c), y_l).ub();
			//if(lb_interval > lb) lb = lb_interval;
			// epsilon relativo: if |lb|<1: lb+eps, otherwise lb + |lb|*eps
			if(fabs(lb_interval) < 1 && lb_interval + epsilon > lb)
				lb = lb_interval+epsilon;
			else if(fabs(lb_interval) >= 1 && lb_interval + fabs(lb_interval)*epsilon > lb)
				lb = lb_interval + fabs(lb_interval)*epsilon;


			// if derivate is empty the segment should not be created
			if(derivate.is_empty()) {
				cout << "derivate empty -> remove interval" << endl;
				//getchar();
				break;
			}
			// contract Newton from left
			point_t = inter.lb();
			point_c = pf.eval(point_t).ub();
			t_before = NEG_INFINITY;

			while(derivate.ub() > 0 and point_t - t_before > error and point_t < inter.ub()) {
				t_before = point_t;

				if(0 == derivate.ub())
					point_t = POS_INFINITY;
				else{
					cout << (fabs(lb) < 1) << endl;
					point_t = (lb - point_c)/derivate.ub() + t_before;
					point_c = pf.eval(point_t).ub();
				}
				// error en caso de que c sea mayor al lb+epsilon
				if(point_t < inter.ub() and point_c > lb+epsilon) {
					cout << "ERRROR LEFT: point right is greater than lb " << endl;
					//getchar();
					break;
					//exit(-1);
				}

				cout << "point_t " << point_t << endl;
			}

			// Se elimina el intervalo ya que no contiene una solucion mejor a lb+epsilon
			if(point_t >= inter.ub()) {
				continue;
			} else {
				//se contracta el intervalo si point_t > 0
				if(point_t > 0)
					inter = Interval(point_t, inter.ub());
			}

			cout << "contract left inter diam " << inter.diam() << " " << inter << endl;

			// contract Newton from right
			point_t = inter.ub();
			point_c = pf.eval(point_t).ub();
			t_before = NEG_INFINITY;

			while(derivate.lb() < 0 and t_before - point_t > error and point_t > inter.lb()) {
				t_before = point_t;

				if(0 == derivate.lb())
					point_t = NEG_INFINITY;
				else{
					point_t = t_before - (lb - point_c)/derivate.lb();
					point_c = pf.eval(point_t).ub();
				}

				// error en caso de que c sea mayor al lb+epsilon
				if(point_t > inter.lb() and point_c > lb+epsilon) {
					cout << "ERRROR RIGHT: point right is greater than lb" << endl;
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

			cout << "contract right inter diam " << inter.diam() << " " << inter << endl;

			// bisect interval and push in stack
			if(inter.is_bisectable() and inter.diam() > max_diam) {
				pair<Interval,Interval> bsc = inter.bisect(0.5);
				cout << "bsc1 diam " << bsc.first.diam() << " " << bsc.first << endl;
				cout << "bsc2 diam " << bsc.second.diam() << " " << bsc.second << endl;
				pila.push(bsc.first);
				pila.push(bsc.second);
			}
			cout << "iteracion " << iter << " pila " << pila.size() << endl;
		}

		/*
		// step method
		double step=1e-2;
		double tinf=0.0;
		double min = POS_INFINITY;
		double max = NEG_INFINITY;
		Interval maxinterval(0);
		while(tinf<1.0){
			double tsup=tinf+step;
			if(tsup>1.0) tsup=1.0;
			double ub = pf.eval(Interval(tinf,tsup)).ub();
			double lb = pf.eval(Interval(tinf,tsup)).lb();
			if(ub > max){maxinterval=Interval(tinf,tsup); max=ub;}
			if(lb < min) min=lb;
			tinf=tsup;

		}
		*/


		// pair <double,double> d = optimize_pf(pf, false);

		// TODO: Que hacer cuando no existaderivada o sea muy grande
		if(derivate.is_empty()) {
			cout << "derivate empty -> remove interval" << endl;
			cout << "Que hacer cuando no existaderivada o sea muy grande" << endl;
			// epsilon relativo: if |lb|<1: lb+eps, otherwise lb + |lb|*eps
			if(fabs(pf.eval(1.0).ub()) < 1)
				lb = pf.eval(1.0).ub()+epsilon;
			else if(fabs(pf.eval(1.0).ub()) >= 1)
				lb = pf.eval(1.0).ub() + fabs(pf.eval(1.0).ub())*epsilon;
			// getchar();
		}

		cout << "optim Newton:" << lb <<  endl;


		// cout << "optim:" << d.second <<  endl;

		//TODO: como obtener los maximosy minimos de c?
		cout << "como obtener los maximosy minimos de c?" << endl;
		cout << "max c " << max_c << endl;
		cout << "min c " << min_c << endl;
		cout << "min c2 " << min_c2 << endl;

		cout << "min y max " << pf.eval(Interval(0,1)) << endl;
		cout << "min " << pf.eval(1.0) << endl;

		/*
		if(fabs(d.second-lb) > 1) {
			cout << "Diferencia  entre optim newton y optim" << endl;
			// getchar();
		}
		*/



		// obtiene los dos puntos para generar la recta obtenida con el metodo de Newton
		cout << "ya1, ya2: " << ya2.ub() << "," << ya1.ub() << endl;
		cout << "yb1, yb2: " << yb2.ub() << "," << yb1.ub() << endl;
		cout << "optim Newton: " << lb <<  endl;
		cout << "m: " << m.ub() << endl;
		cout << "max c: " << max_c.ub() << endl;

		Interval x1, y1, x2, y2;
		/*
		// si lb es 0 la recta pasa por ya y yb
		if(lb == 0) {
			y1 = ya2;
			x1 = ya1;
			y2 = yb2;
			x2 = yb1;
			// cout << "point1: " << x1.ub() << "," << y1.ub() << endl;
			// cout << "point2: " << x2.ub() << "," << y2.ub() << endl;
		} else if(lb < max_c.ub()) { // si c < max_c existe una recta
			// primer punto (x1, y1)
			x1 = ya1;
			y1 = (x1 - lb)/m;
			// segundo punto (x2, y2)
			y2 = yb2;
			x2 = m*y2 + lb;
		}
		*/

		cout << "punto1 (" << ya1.ub() << "," << (ya1*m) + lb << ")" << endl;
		cout << "punto2 (" << yb1.ub() << "," << (yb1*m) + lb << ")" << endl;
		// corte horizontal
		y1 = ya2;
		x1 = (y1-lb)/m;
		// corte vertical
		x2 = yb1;
		y2 = (x2*m) + lb;

		cout << "ya1, ya2: " << ya1.ub() << "," << ya2.ub() << endl;
		cout << "yb1, yb2: " << yb1.ub() << "," << yb2.ub() << endl;


		cout << "Se agregan el punto la recta de X ---" << endl;
		cout << "ya1, ya2: " << ya1.ub() << "," << ya2.ub() << endl;
		addPointtoNDS(make_pair(ya1.ub(), ya2.ub()));
		//getchar();
		cout << "Se agregan el punto la recta de X ---" << endl;
		cout << "yb1, yb2: " << yb1.ub() << "," << yb2.ub() << endl;
		addPointtoNDS(make_pair(yb1.ub(), yb2.ub()));
		cout << lb << " " << max_c.ub() << endl;
		// getchar();

		// Si lb no esta entre los rangos permitidos no se agrega nada
		// if(lb < 0 || lb >= max_c.ub()) return;

		// obtiene funcion
		std::vector< pair <double, double> > functionPoly;
		Interval a;
		IntervalVector interVector = IntervalVector(2);
		double value;
		int max_iterations = 10;
		for(int i=0;i <= max_iterations; i++) {
			a = Interval((double) i/max_iterations);
			//cout << a << endl;
			interVector = pf.get_point(a);
			//cout << interVector[0].ub() << "," <<  interVector[1].ub() << endl;
			functionPoly.push_back(make_pair(interVector[0].ub(), interVector[1].ub()));
		}


		cout << "Se agrega una recta o punta---" << endl;
		std::vector< pair <double, double> > rectaUB;
		if(lb == 0 or (lb != 0 and lb < max_c.ub()) ) {
			if(x1.ub() == x2.ub() and  y1.ub() == y2.ub()) {
				cout << "point: " << x1.ub() << "," << y1.ub() << endl;
				rectaUB.push_back(make_pair(x1.ub(), y1.ub()));
				rectaUB.push_back(make_pair(x1.ub(), y1.ub()));
				cout << "antes" << endl;
				if(_plot) py_Plotter::offline_plot(NULL, NDS2, rectaUB, functionPoly);
				getchar();
				addPointtoNDS(make_pair(x1.ub(), y1.ub()));
				//getchar();
			}else {
				rectaUB.push_back(make_pair(x1.ub(),y1.ub()));
				rectaUB.push_back(make_pair(x2.ub(), y2.ub()));
				cout << "antes" << endl;
				if(_plot) py_Plotter::offline_plot(NULL, NDS2, rectaUB, functionPoly);
				getchar();

				//if(_plot) py_Plotter::offline_plot(NULL, NDS2);
				cout << "se va a gregar el vector "<< endl;
				cout << "point1: " << x1.ub() << "," << y1.ub() << endl;
				cout << "point2: " << x2.ub() << "," << y2.ub() << endl;
				cout << NDS2.size() << endl;
				std::map<pair<double, double>, IntervalVector>::iterator it2 = NDS2.begin();
				for(it2=NDS2.begin(); it2!=NDS2.end(); ++it2){
					cout << "point (" << it2->first.first << "," << it2->first.second << ")" << endl;
				}
				if(x1.ub() != x1.ub() || y1.ub() != y1.ub()) {
					//cout << "optim:" << d.second <<  endl;
					cout << "max c " << max_c << endl;
					cout << "min c" << min_c << endl;

					cout << "min " << pf.eval(Interval(0,1)) << endl;
					cout << "max " << pf.eval(1.0) << endl;
					cout << "bad point 1" << endl;
					// getchar();
				}
				addPointtoNDS(make_pair(x1.ub(),y1.ub()));
				if(x2.ub() != x2.ub() || y2.ub() != y2.ub()) {
					//cout << "optim:" << d.second <<  endl;
					cout << "max c " << max_c << endl;
					cout << "min c" << min_c << endl;

					cout << "min " << pf.eval(Interval(0,1)) << endl;
					cout << "max " << pf.eval(1.0) << endl;
					cout << "bad point 2" << endl;
					// getchar();
				}
				addPointtoNDS(make_pair(x2.ub(), y2.ub())); //error
				addVectortoNDS(make_pair(x1.ub(),y1.ub()), make_pair(x2.ub(), y2.ub()));
				//if(_plot) py_Plotter::offline_plot(NULL, NDS2);
				cout << "point1: " << x1.ub() << "," << y1.ub() << endl;
				cout << "point2: " << x2.ub() << "," << y2.ub() << endl;
				cout << NDS2.size() << endl;

				// if(_plot) py_Plotter::offline_plot(NULL, NDS2);
				/*
				for(it2=NDS2.begin(); it2!=NDS2.end(); ++it2){
					cout << "point (" << it2->first.first << "," << it2->first.second << ")" << endl;
				}*/
				// getchar();
			}
		} else {
			rectaUB.push_back(make_pair(0,0));
			rectaUB.push_back(make_pair(0,0));
		}

		//if(_plot) py_Plotter::offline_plot(NULL, NDS);
		// cout << "Sin NDS2 plot NDS" << endl;
		//getchar();


		if(_plot) py_Plotter::offline_plot(NULL, NDS2);
		cout << "Sin NDS2 plot NDS2" << endl;
		// NDS2, recta, funcion
		cout << "optim Newton:" << lb <<  endl;
		cout << "max c " << max_c.ub() << endl;
		cout << "min c " << min_c.ub() << endl;

		cout << "pendiente " << m.ub() << endl;

		cout << "punto1  (" << ya1.ub() << "," << ya2.ub() << ")" << endl;
		cout << "punto2 (" << yb1.ub() << "," << yb2.ub() << ")" << endl;

		// NDS2 listo
		cout << "NDS2 " << endl;
		map< pair <double, double>, IntervalVector > :: iterator ub=NDS2.begin();
		for(;ub!=NDS2.end();ub++){
			cout << "(" << ub->first.first << "," << ub->first.second << ")" << endl;
		}
		// recta listo
		cout << "recta " << endl;
		for (int i=0;i<rectaUB.size();i++) {
			cout << rectaUB[i].first << " " << rectaUB[i].second << endl;
		}
		// funcion listo
		cout << "funcion " << endl;
		for (int i=0;i<functionPoly.size();i++) {
			cout << functionPoly[i].first << " " << functionPoly[i].second << endl;
		}
		if(_plot) py_Plotter::offline_plot(NULL, NDS2, rectaUB, functionPoly);
		getchar();

	}


private:



	/**
	 * min feasible value found for each objective
	 */
    pair <double, double> y1_ub, y2_ub;

	/* Remember return status of the last optimization. */
	Status status;



	//TODO: this should not be static

	/** The current non-dominated set sorted by increasing x */
	static map< pair <double, double>, IntervalVector > NDS;

	/** The current non-dominated set sorted by decreasing y */
	map< pair <double, double>, IntervalVector, sorty > NDSy;

	/** The current non-dominated set sorted by increasing x */
	static map< pair <double, double>, IntervalVector, sorty2 > NDS2;

	/* CPU running time of the current optimization. */
	double time;

	/** Number of cells pushed into the heap (which passed through the contractors) */
	int nb_cells;

};

inline OptimizerMOP::Status OptimizerMOP::get_status() const { return status; }

inline double OptimizerMOP::get_time() const { return time; }

inline double OptimizerMOP::get_nb_cells() const { return nb_cells; }




} // end namespace ibex

#endif // __IBEX_OPTIMIZER_H__
