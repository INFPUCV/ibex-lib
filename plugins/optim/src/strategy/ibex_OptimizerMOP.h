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
#include "ibex_CellBufferOptim.h"
#include "ibex_CellSet.h"
//#include "ibex_EntailedCtr.h"
#include "ibex_CtcKhunTucker.h"

#include <set>
#include <map>
#include <list>

using namespace std;
namespace ibex {

class point2{
public:
	Interval x;
	Interval y;

	point2()  : x(0), y(0)  {}
	point2(double x, double y) : x(Interval(x)), y(Interval(y)) { }
	point2(Interval x, Interval y) : x(x), y(y) { }

	point2 operator+(const point2& p2) const {
		return point2( x+p2.x, y+p2.y );
	}

	point2 operator-(const point2& p2) const {
		return point2( x-p2.x, y-p2.y );
	}

	Interval operator*(const point2& p2) const {
		return ( x*p2.y - y*p2.x );
	}

	bool operator<(const point2& p2) const {
		if(x.mid()!=p2.x.mid()) return x.mid() < p2.x.mid();
		if(y.mid()!=p2.y.mid()) return y.mid() > p2.y.mid();
		return false;
	}

};
/**
 * \defgroup optim IbexOpt
 */

/**
 * \ingroup optim
 *
 * \brief Global MOP Optimizer.
 *
 */
class OptimizerMOP {

public:

	/**
	 * \brief Return status of the optimizer
	 *
	 * See comments for optimize(...) below.
	 */
	typedef enum {SUCCESS, INFEASIBLE, NO_FEASIBLE_FOUND, UNBOUNDED_OBJ, TIME_OUT, UNREACHED_PREC} Status;

	/**
	 *  \brief Create an optimizer.
	 *
	 * Inputs:
	 *   \param n        - number of variables or the <b>original system</b>
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
	OptimizerMOP(int n, const Array<NumConstraint>& ctcs, const Function &f1,  const Function &f2,
			Ctc& ctc, Bsc& bsc, CellBufferOptim& buffer, LoupFinderMOP& finder, double rel_eps=default_rel_eps, double abs_eps=default_abs_eps);

	/**
	 * \brief Delete *this.
	 */
	virtual ~OptimizerMOP();

	/**
	 * \brief Run the optimization.
	 *
	 * \param init_box             The initial box
	 * \param obj_init_bound       (optional) can be set when an initial upper bound of the objective minimum is known a priori.
	 *                             This bound can be obtained, e.g., by a local solver. This is equivalent to (but more practical
	 *                             than) adding a constraint f(x)<=obj_init_bound.
	 *
	 * \return SUCCESS             if the global minimum (with respect to the precision required) has been found.
	 *                             In particular, at least one feasible point has been found, less than obj_init_bound, and in the
	 *                             time limit.
	 *
	 *         INFEASIBLE          if no feasible point exist less than obj_init_bound. In particular, the function returns INFEASIBLE
	 *                             if the initial bound "obj_init_bound" is LESS than the true minimum (this case is only possible if
	 *                             goal_abs_prec and goal_rel_prec are 0). In the latter case, there may exist feasible points.
	 *
	 *         NO_FEASIBLE_FOUND   if no feasible point could be found less than obj_init_bound. Contrary to INFEASIBLE,
	 *                             infeasibility is not proven here. Warning: this return value is sensitive to the abs_eps_f and
	 *                             rel_eps_f parameters. The upperbounding makes the optimizer only looking for points less than
	 *                             min { (1-rel_eps_f)*obj_init_bound, obj_init_bound - abs_eps_f }.
	 *
	 *         UNBOUNDED_OBJ       if the objective function seems unbounded (tends to -oo).
	 *
	 *         TIMEOUT             if time is out.
	 *
	 *         UNREACHED_PREC      if the search is over but the resulting interval [uplo,loup] does not satisfy the precision
	 *                             requirements. There are several possible reasons: the goal function may be too pessimistic
	 *                             or the constraints function may be too pessimistic with respect to the precision requirement
	 *                             (which can be too stringent). This results in tiny boxes that can neither be contracted nor
	 *                             used as new loup candidates. Finally, the eps_x parameter may be too large.
	 *
	 */
	Status optimize(const IntervalVector& init_box);

	/* =========================== Output ============================= */

	/**
	 * \brief Displays on standard output a report of the last call to optimize(...).
	 *
	 * Information provided:
	 * <ul><li> interval of the cost  [uplo,loup]
	 *     <li> the best feasible point found
	 *     <li> total running time
	 *     <li> total number of cells (~boxes) created during the exploration
	 * </ul>
	 */
	void report(bool verbose=true);

	void plot(Cell* current);

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
	map< pair <double, double>, IntervalVector >& get_UB()  { return UB; }

	std::set< point2 >& get_LB()  { return LB; }

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


	int get_nb_sols() const {return nb_sols;}
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
	 * \brief Constraints
	 */
	const Array<NumConstraint>& ctrs;

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

	/** Precision (bisection control objective functions) */
	double abs_eps;

	const double rel_eps;

	/** Default absolute precision: 1e-7 */
	static const double default_abs_eps;

	/** Default relative precision: 0.1 */
	static const double default_rel_eps;

	/**
	 * \brief Trace activation flag.
	 *
	 * The value can be fixed by the user.
	 * - 0 (by default): nothing is printed
	 * - 1:              prints every loup/uplo update.
	 * - 2:              prints also each handled node (warning: can generate very
	 *                   long trace).
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


	/**
	 * \brief returns true if the box+z1 + a*z2 > w_lb is dominated by the ub_set
	 */
	static double distance2(const Cell* c){
		//if(c->get<CellBS>().ub_distance != POS_INFINITY) return c->get<CellBS>().ub_distance;
		double max_dist=NEG_INFINITY;
		if(UB.size()==2) return POS_INFINITY;

		int n=c->box.size();

		Interval z1 = c->box[n-2];
		Interval z2 = c->box[n-1];
		double a = c->get<CellBS>().a;
		double w_lb = c->get<CellBS>().w_lb;
		
		//TODO: optimize this
		map< pair <double, double>, IntervalVector >::iterator it = UB.begin();


		for(;it!=UB.end(); ){
			pair <double, double> p = it->first; it++;
			if(it==UB.end()) break;
			pair <double, double> p2 = it->first;

			pair <double, double> pmax= make_pair(p2.first, p.second);
//			cout << "pmax: (" << pmax.first <<"," << pmax.second << ")" << endl;
//			cout << "z: (" << z1.lb() <<"," << z2.lb() << ")" << endl;
//			cout << "a: " << a << endl;
//			cout << "w_lb: " << w_lb << endl;



			//el punto esta dentro de la zona de interes
			if(pmax.first > z1.lb() && pmax.second > z2.lb()){
				double dist = std::min (pmax.first - z1.lb(), pmax.second - z2.lb());
				//here we add the distance to the line
			    //dist = std::min (dist, (Interval(pmax.first) + Interval(a)*Interval(pmax.second) - Interval(w_lb)).ub() );

				// distancia Damir
				dist = std::min(dist, z2.lb() - ( (w_lb - pmax.first + pmax.second)/(a+1.0) ));
//				cout << "damir " << pmax.second - (w_lb - pmax.first - pmax.second)/(a - 1) << endl;
//				cout << "damir2 " <<  z2.lb() - ( (w_lb - pmax.first + pmax.second)/(a+1.0) ) << endl;
				//Damir's distance
				// dist = std::min(dist, (Interval(pmax.second)-(Interval(w_lb) - (Interval(pmax.first) - Interval(pmax.second)))/(Interval(a)+1.0)).ub());

//				cout << "nacho " << (Interval(pmax.second)-(Interval(w_lb) - (Interval(pmax.first) - Interval(pmax.second)))/(Interval(a)-1.0)).ub() << endl;
				if(dist > max_dist) max_dist=dist;

			}
		}

		return max_dist;
	}


	Interval y1,y2;

  //TODO: make it conservative!
	void insert_lb_segment(point2 p1, point2 p2){
      //trace=1;
		  if(trace) cout << "p1-p2: (" << p1.x.mid() << "," << p1.y.mid() << ") --> (" << p2.x.mid() << "," << p2.y.mid() << ")" << endl;
	    if(LB.size()==0){
			LB.insert(point2(y1.lb(),y2.ub()));
	    	LB.insert(point2(y1.ub(),y2.ub()));
	    	LB.insert(point2(y1.ub(),y2.lb()));

	    }

		point2 p1_p = point2(p1.x,y2.ub());
		point2 p2_p = point2(y1.ub(),p2.y);

		//point2 p1_p = point2(p1.x,1e10);
		//point2 p2_p = point2(1e10,p2.y);

		std::set< point2 > new_points;
		std::set< point2 >::iterator it=LB.upper_bound(p1);
		it--;

		point2 v1(it->x, it->y);
		it++;
		point2 v2(it->x, it->y);

		bool in = false;


		point2 s;
		if (intersect(v1, v2, p1_p,  p1, s)) {
			 if(trace) cout << "s1: (" << s.x << "," << s.y << ")" << endl;
			 if(s.x.lb()==NEG_INFINITY) exit(0);

      if(( (s.x==p1.x && s.y==p1.y) && ((p2-p1)*(v2-v1)).lb() <= 0 )){
				new_points.insert(s);
				new_points.insert(p1);
				in = true;  if(trace) cout << "in"  << endl;
				//it--; LB.erase(v2); it++; v1 = v2; v2=*it;
			}


			else if(s.x!=p1.x || s.y!=p1.y){
				new_points.insert(s);
				new_points.insert(p1);
				in = true;  if(trace) cout << "in"  << endl;
				//it--; LB.erase(v2); //it++; v1 = v2; v2=*it;
			}
			//it++; v1 = v2; v2=*it;
		}



	    while(v1.y.mid() > p2.y.mid()){

	        if (intersect(v1,v2, p1, p2, s)){
	           if(trace) cout << "s2: (" << s.x << "," << s.y << ")" << endl;
						if(s.x.lb()==NEG_INFINITY) exit(0);

	          in=!in;
              if(trace) cout << ((in)? "in":"out")  << endl;
	          new_points.insert(point2(s.x.lb(),s.y.lb()));
	        }

	        if ( v2.y.mid() <= p2.y.mid() && intersect(v1,v2, p2, p2_p,s)){
	           if(trace) cout << "s3: (" << s.x << "," << s.y << ")" << endl;
						if(s.x.lb()==NEG_INFINITY) exit(0);
	          in = false;
						 if(trace) cout << "out" << endl;
	          new_points.insert(point2(s.x.lb(),s.y.lb()));
	          new_points.insert(p2);
	          break;
	        }

	        if(in){ it--; LB.erase(v2); }

	        v1 = v2;
	        it++;
	        v2 = *it;
	    }

	    LB.insert(new_points.begin(), new_points.end());

/*
		cout << "LB points:" << endl;
		for(it=LB.begin(); it!=LB.end(); it++){
			cout << "(" << it->x << "," << it->y << ")" << endl;
		}*/
	}

protected:

	/**
	 * \brief Main procedure for processing a box.
	 *
	 * <ul>
	 * <li> contract the cell box and try to find a new loup (see contract_and_bound)
	 * <li> push the cell into the buffer or delete the cell in case of empty box detected.
	 * </ul>
	 *
	 */
	void handle_cell(Cell& c, const IntervalVector& init_box);

	/**
	 * \brief Contract and bound procedure for processing a box.
	 *
	 * <ul>
	 * <li> contract the cell's box w.r.t the "loup",
	 * <li> contract with the contractor ctc,
	 * <li> search for a new loup,
	 * <li> (optional) call the first order contractor
	 * </ul>
	 *
	 */
	void contract_and_bound(Cell& c, const IntervalVector& init_box);


  /**
  * \brief intersect two segments and return the intersection res
  * https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
  */
  bool intersect(point2 p, point2 p2,
		point2 q,  point2 q2, point2& res){

	  	  if(trace) cout << "p-p2: (" << p.x.mid() << "," << p.y.mid() << ") --> (" << p2.x.mid() << "," << p2.y.mid() << ")" << endl;
	  	  if(trace) cout << "q-q2: (" << q.x.mid() << "," << q.y.mid() << ") --> (" << q2.x.mid() << "," << q2.y.mid() << ")" << endl;

        //orthogonal segments
	  	  if(p.x.mid()==p2.x.mid() && q.y.mid()==q2.y.mid()) {
	  		  if(q.x.mid()>p.x.mid() || q2.x.mid() < p.x.mid() ||
					q.y.mid() < p2.y.mid() || q.y.mid() > p.y.mid()) return false;
	  		  res.x=p.x;
	  		  res.y=q.y;
	  		  return true;
	  	  }

        //orthogonal segments
	  	  if(p.y.mid()==p2.y.mid() && q.x.mid()==q2.x.mid()){
	  		  if(p.x.mid() > q.x.mid() || p2.x.mid() < q.x.mid() ||
					p.y.mid() < q2.y.mid() || p.y.mid() > q.y.mid()) return false;
	  		  res.x=q.x;
	  		  res.y=p.y;
	  		  return true;
	  	  }

	  	  if( (p.x.mid()==p2.x.mid() && p.y.mid()==p2.y.mid()) ||
				 (q.x.mid()==q2.x.mid() && q.y.mid()==q2.y.mid())) return false;

      //  if(p2.x.mid()>=1e9) p2.x = std::max(std::max(q.x.ub(),q2.x.ub()),p.x.ub());
      //  if(p.y.mid()>=1e9) p.y = std::max(std::max(q.y.ub(),q2.y.ub()),p2.y.ub());
      //  if(q2.x.mid()>=1e9) q2.x = std::max(std::max(p.x.ub(),p2.x.ub()),q.x.ub());
      //  if(q.y.mid()>=1e9) q.y = std::max(std::max(p.y.ub(),p2.y.ub()),q2.y.ub());



			point2 r = p2-p;
			point2 s = q2-q;

			//now we find a solution for the equation p+tr = q+us,


      if( (r * s).lb() > -1e-8 && (r * s).ub() < 1e-8) {
				//segments are collinear
				if((q - p) * r == 0){
					//segments are intersecting
					if(p2.y.mid() > q2.y.mid() || (p2.y.mid() == q2.y.mid() && p2.x.mid()<q2.x.mid())) res=p2;
					else res=q2;
					return true;
				}else //segments have no intersection
				  return false;

			}

			Interval t = ((q - p) * s) / (r * s);
			Interval u = ((p - q) * r) / (s * r);

      //cout << (r * s) << endl;


			if (t.ub()>=0 && t.lb() <=1  && u.ub()>=0 && u.lb() <=1){
				res = p + point2(t*r.x,t*r.y);
				//cout << "res::(" << res.x << "," << res.y << ")"  << endl;
				return true;
			}

			return false;
		}

	/**
	 * \brief Main procedure for updating the loup.
	 */
	bool update_UB(const IntervalVector& box, int n);

	Interval compute_lb_hypervolume(){

    Interval volume(0.0);
    point2 lb1 = *LB.begin();

    std::set< point2 >::iterator lb2=LB.begin();

    cout << y1_max << "," << y2_max << endl;

		for(lb2++;lb2!=LB.end();lb2++){
          if( (lb1.y.mid() <= y2_max) && ( lb2->x.mid() <= y1_max ) )
        	  volume += (( y2_max - lb1.y ) + (lb1.y - lb2->y)/2.0) * ( lb2->x - lb1.x );

          else if( lb2->x.mid() > y1_max ){
        	  volume += ( y2_max - lb1.y )  * ( y1_max - lb1.x );
        	  break;
          }

          lb1=*lb2;
		}

		return volume;
	}

	Interval compute_ub_hypervolume(){

		//cout << y1_max << ";" <<y2_max << endl;
        Interval volume(0.0);
        pair <double, double> ub1 = UB.begin()->first;
        map< pair <double, double>, IntervalVector >::iterator _ub2= UB.begin();

		for(;_ub2!=UB.end();_ub2++){

			pair <double, double> ub2=_ub2->first;
			pair <double, double> ubx=make_pair(ub2.first, ub1.second);

          if( (ub1.second <= y2_max) && ( ub2.first <= y1_max ) )
        	  volume += (Interval(y2_max) - ub1.second ) * ( Interval(ubx.first) - ub1.first );

          else if( ub2.first > y1_max ){
        	  volume += ( Interval(y2_max) - ub1.second )  * ( Interval(y1_max) - ub1.first );
        	  break;
          }

          ub1=ub2;
		}

		return volume;
	}


private:

	/**
	 * \brief Evaluate the goal in the point x
	 */
	Interval eval_goal(const Function& goal, IntervalVector& x);

	/**
	 * min feasible value found for the objectives
	 */
    pair <double, double> y1_ub, y2_ub;

    /**
     * the max possible value for the objectives s.t. the other objective is minimized
     */
    double y1_max, y2_max;

	/** Currently entailed constraints */
	//EntailedCtr* entailed;

	//!! warning: sys.box should be properly set before call to constructor !!
	//CtcKhunTucker kkt;

	/* Remember return status of the last optimization. */
	Status status;

  /** The cells in the buffer for plotting
	 * the set should be updated each time the real buffer is popped
	 * and pushed.
	 */
	set< Cell* > buffer_cells;

	/** The current upper bounds (f1(x), f2(x)) of the pareto front associated
	 * to its corresponding  point x
	 */
	static map< pair <double, double>, IntervalVector > UB;

	/**
	 * A set of points denoting the segments related to the lowerbound of the
	 * pareto front.
	 */
	std::set< point2 > LB;

	int nb_sols;


	/** True if loup has changed in the last call to handle_cell(..) */
	//bool loup_changed;

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
