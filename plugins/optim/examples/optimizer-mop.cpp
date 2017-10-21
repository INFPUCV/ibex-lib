//============================================================================
//                                  I B E X                                   
// File        : optimizer04.cpp
// Author      : Gilles Chabert  Bertrand Neveu
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Jul 12, 2012
// Last Update : Jul 12, 2012
//============================================================================


#include "ibex.h"



#ifndef _IBEX_WITH_OPTIM_
#error "You need the plugin Optim to run this example."
#endif

const double default_relax_ratio = 0.2;

using namespace std;
using namespace ibex;
int main(int argc, char** argv){


	// ------------------------------------------------
	// Parameterized Optimizer (with a system loaded from a file, and choice of contractor, linearization  and bisector)
        // Load a problem to optimize
	// --------------------------
	try {

	if (argc<8) {
		cerr << "usage: optimizer04 filename filtering linear_relaxation bisection strategy prec timelimit "  << endl;
		exit(1);
	}

	//restricciones del sistema original + goal=NULL
	System ext_sys(argv[1]);

	SystemFactory fac2;

	Variable w;
	Variable a;

	fac2.add_var(w);
	fac2.add_var(a);

	fac2.add_var(ext_sys.args[ext_sys.nb_var-2]);
	fac2.add_var(ext_sys.args[ext_sys.nb_var-1]);
	fac2.add_ctr(ext_sys.args[ext_sys.nb_var-2] + a * ext_sys.args[ext_sys.nb_var-1] - w = 0);

	System _ext_sys(ext_sys, System(fac2));

	cout << _ext_sys << endl;

	string filtering = argv[2];
	string linearrelaxation= argv[3];
	string bisection= argv[4];
	string strategy= argv[5];
	int nbinput=5;

	double prec= atof(argv[nbinput+1]);
	double timelimit = atof(argv[nbinput+2]);
	double eqeps= 1.e-8;

	RNG::srand(atoi(argv[nbinput+3]));

	// the extended system 
	// restricciones del sistema original + variables objetivo y restricciones

	//ExtendedSystem ext_sys(sys, eqeps);

	SystemFactory fac;
	for(int i=0; i<ext_sys.nb_var-2; i++ )
		fac.add_var(ext_sys.args[i]);


	for(int j=2; j<ext_sys.nb_ctr; j++ )
		fac.add_ctr(ext_sys.ctrs[j]);


	System sys(fac);
	for(int i=0; i<sys.nb_var; i++ )
		sys.box[i] = ext_sys.box[i];


	cout << sys << endl;

	IntervalVector box = ext_sys.box.mid();
	box[sys.nb_var]=0;
	box[sys.nb_var+1]=0;

	cout << ext_sys.ctrs[0].f.eval(box) << endl;
	cout << ext_sys.ctrs[1].f.eval(box) << endl;

	LoupFinderMOP finder(sys, ext_sys.ctrs[0].f, ext_sys.ctrs[1].f);

	//NormalizedSystem norm_sys(sys,eqeps);
	//LoupFinderDefault loupfinder (norm_sys,true);
	//LoupFinderDefault loupfinder (norm_sys,false);

	CellBufferOptim* buffer;
	if(strategy=="minlb")
	  buffer = new CellSet<OC1>;
	else if(strategy=="weighted_sum")
	  buffer = new CellSet<weighted_sum>;
	else if(strategy=="NDSdist")
	  buffer = new DistanceSortedCellBufferMOP;
	else if(strategy=="NDSsize")
	  buffer = new CellFeasibleDiving<maxsize>(*new CellNSSet);
	else if(strategy=="diving-minlb")
      buffer = new CellFeasibleDiving<minLB>(*new CellSet<minLB>);
	else if(strategy=="diving-weighted_sum")
	  buffer = new CellFeasibleDiving<weighted_sum>(*new CellSet<weighted_sum>);
	else if(strategy=="diving-NDSdist")
	  buffer = new CellFeasibleDiving<max_distance>(*new DistanceSortedCellBufferMOP);
	else if(strategy=="diving-NDSsize")
	  buffer = new CellFeasibleDiving<maxsize>(*new CellNSSet);

	/*else
		buffer = new CellDoubleHeap  (ext_sys);*/

	//        cout << "file " << argv[1] << endl;

	// Build the bisection heuristic
	// --------------------------

	Bsc * bs;

	if (bisection=="roundrobin")
	  bs = new RoundRobin (prec);
	else if (bisection== "largestfirst")
          bs= new LargestFirst(prec);
	else if (bisection=="smearsum")
	  bs = new SmearSum(ext_sys,prec);
	else if (bisection=="smearmax")
	  bs = new SmearMax(ext_sys,prec);
	else if (bisection=="smearsumrel")
	  bs = new SmearSumRelative(ext_sys,0);
	else if (bisection=="smearmaxrel")
	  bs = new SmearMaxRelative(ext_sys,0);
	//else if (bisection=="lsmear")
	 // bs = new LSmear(ext_sys,prec);
	else {cout << bisection << " is not an implemented  bisection mode "  << endl; return -1;}

	// The contractors

	// the first contractor called
	CtcHC4 hc4(_ext_sys.ctrs,0.01,true);
	// hc4 inside acid and 3bcid : incremental propagation beginning with the shaved variable
	CtcHC4 hc44cid(_ext_sys.ctrs,0.1,true);
	// hc4 inside xnewton loop 
	CtcHC4 hc44xn (_ext_sys.ctrs,0.01,false);

	// The 3BCID contractor on all variables (component of the contractor when filtering == "3bcidhc4") 
	Ctc3BCid c3bcidhc4(hc44cid);
	// hc4 followed by 3bcidhc4 : the actual contractor used when filtering == "3bcidhc4" 
	CtcCompo hc43bcidhc4 (hc4, c3bcidhc4);

	// The ACID contractor (component of the contractor  when filtering == "acidhc4")
	CtcAcid acidhc4(_ext_sys,hc44cid,true);
	// hc4 followed by acidhc4 : the actual contractor used when filtering == "acidhc4" 
	CtcCompo hc4acidhc4 (hc4, acidhc4);

      

	Ctc* ctc;
	if (filtering == "hc4")
	  ctc= &hc4;
	else if
	  (filtering =="acidhc4")   
	  ctc= &hc4acidhc4;
	else if 
	  (filtering =="3bcidhc4")
	  ctc= &hc43bcidhc4;
	else {cout << filtering <<  " is not an implemented  contraction  mode "  << endl; return -1;}

	Linearizer* lr;
	if (linearrelaxation=="art")
	  lr= new LinearizerCombo(_ext_sys,LinearizerCombo::ART);
	else if  (linearrelaxation=="compo")
	  lr= new LinearizerCombo(_ext_sys,LinearizerCombo::COMPO);
	else if (linearrelaxation=="xn")
	  lr= new LinearizerXTaylor (_ext_sys, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM_OPP);
	//	else {cout << linearrelaxation  <<  " is not an implemented  linear relaxation mode "  << endl; return -1;}
	// fixpoint linear relaxation , hc4  with default fix point ratio 0.2
	CtcFixPoint* cxn;
	CtcPolytopeHull* cxn_poly;
	CtcCompo* cxn_compo;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn")
          {
		cxn_poly = new CtcPolytopeHull(*lr);
		cxn_compo =new CtcCompo(*cxn_poly, hc44xn);
		cxn = new CtcFixPoint (*cxn_compo, default_relax_ratio);
	  }
	//  the actual contractor  ctc + linear relaxation 
	Ctc* ctcxn;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn")
          ctcxn= new CtcCompo  (*ctc, *cxn); 
	else
	  ctcxn = ctc;

	// the optimizer : the same precision goalprec is used as relative and absolute precision
	OptimizerMOP o(sys.nb_var,sys.ctrs,ext_sys.ctrs[0].f,ext_sys.ctrs[1].f, *ctcxn,*bs,*buffer,finder,prec);
	max_distance::UB= &o.get_UB();


	//	cout << " sys.box " << sys.box << endl;

	// the trace 
	o.trace=0;

	// the allowed time for search
	o.timeout=timelimit;

	// the search itself 
	o.optimize(ext_sys.box);

	// printing the results     
	o.report(false);
       cout << o.get_time() << "  " << o.get_nb_cells() << " "<< o.get_nb_sol() << endl;

	//	if (filtering == "acidhc4"  )
	//cout    << " nbcidvar " <<  acidhc4.nbvar_stat() << endl;

	delete bs;
	delete buffer;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn") {
		delete lr;
	    delete ctcxn;
	    delete cxn;
	    delete cxn_poly;
	    delete cxn_compo;
	}



	return 0;
	
	}


	catch(ibex::SyntaxError& e) {
	  cout << e << endl;
	}
}
