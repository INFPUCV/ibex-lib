//============================================================================
//                                  I B E X
// File        : ibexmop.cpp
// Author      : Ignacio Araya, Matias Campusano, Damir Aliquintui
// Copyright   : PUCV (Chile)
// License     : See the LICENSE file
// Created     : Jan 01, 2017
// Last Update : Jan 01, 2017
//============================================================================


#include "ibex.h"


using namespace std;
using namespace ibex;
int main(int argc, char** argv){

	System ext_sys("plugins/optim/benchs/hard/ex6_2_10.bch");
	SystemFactory fac;

	for(int i=0; i<ext_sys.args.size(); i++ ){
		fac.add_var(ext_sys.args[i]);
	}



	for(int j=2; j<ext_sys.nb_ctr; j++ ){
		const ExprNode& e=ext_sys.f_ctrs.expr()[j];
		ExprCtr cc(e, ext_sys.ctrs[j].op);
		fac.add_ctr(cc);
		//ExprCtr cc(ext_sys.ctrs[j].f.expr(),  ext_sys.ctrs[j].op);
		//fac.add_ctr(cc);
		delete &e;
	}



	return 0;

}
