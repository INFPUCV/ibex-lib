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
#include "args.hxx"

// Server side C/C++ program to demonstrate Socket programming
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstring>
#include <set> 
#include <iterator>
#include <vector>


#ifndef _IBEX_WITH_OPTIM_MOP_
#error "You need the plugin Optim MOP to run this example."
#endif

const double default_relax_ratio = 0.2;

using namespace std;
using namespace ibex;


int main(int argc, char** argv){

    NDSrp nds;
    nds.clear();
    Vector v(2);
    auto it = nds.NDS.begin(); it++;
    nds.NDS_erase(it);

    v[0]=0; v[1]=POS_INFINITY; nds.NDS_insert(v);
    v[0]=0; v[1]=10; nds.NDS_insert(v);
    v[0]=1; v[1]=9.5; nds.NDS_insert(v);
    v[0]=2; v[1]=8; nds.NDS_insert(v);
    v[0]=3; v[1]=5; nds.NDS_insert(v);
    v[0]=4; v[1]=4.8; nds.NDS_insert(v);
    v[0]=5; v[1]=4.5; nds.NDS_insert(v);
    v[0]=6; v[1]=4; nds.NDS_insert(v);
    v[0]=7; v[1]=3; nds.NDS_insert(v);
    v[0]=8; v[1]=2.8; nds.NDS_insert(v);
    v[0]=9; v[1]=2.7; nds.NDS_insert(v);
    v[0]=10; v[1]=2.5; nds.NDS_insert(v);
    v[0]=POS_INFINITY; v[1]=2.5; nds.NDS_insert(v);

    py_Plotter::offline_plot(nds.NDS, NULL, "output2.txt");

    int i=0;
    for(auto ub : nds.NDS){
        if(i<3) {i++; continue;}
		cout << "(" << (*ub)[0] << ", " << (*ub)[1] << ")" << endl;
        cout << dynamic_cast<Point*>(ub)->compute_hv_contribution() << endl;
    }

	return 0;


}