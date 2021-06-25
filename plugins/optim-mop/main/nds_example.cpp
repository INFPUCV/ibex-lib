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

    NDShv nds;
    Vector v(2);
    v[0]=POS_INFINITY; v[1]=POS_INFINITY;
    nds.erase(&v);

    v[0]=0; v[1]=POS_INFINITY; nds.insert(v);
    v[0]=POS_INFINITY; v[1]=2.5; nds.insert(v);

    v[0]=0; v[1]=10; nds.insert(v);
    v[0]=2; v[1]=10; nds.insert(v);
    v[0]=2; v[1]=8; nds.insert(v);
    v[0]=3; v[1]=8; nds.insert(v);
    v[0]=3; v[1]=5; nds.insert(v);
    v[0]=4; v[1]=4.8; nds.insert(v);
    v[0]=5; v[1]=4.39; nds.insert(v);
    v[0]=6; v[1]=4; nds.insert(v);
    v[0]=7; v[1]=3; nds.insert(v);
    v[0]=8; v[1]=2.8; nds.insert(v);
    v[0]=9; v[1]=2.7; nds.insert(v);
    v[0]=10; v[1]=2.5; nds.insert(v);
    

    for(auto ub :  nds){
		cout << "(" << (*ub)[0] << ", " << (*ub)[1] << ")" << endl;
        cout << dynamic_cast<Point*>(ub)->hv_contribution << endl;
    }

    //nds.pop_front();
    //nds.pop_front();
    set<Vector*, sort_rp> aux;
    for(auto ub : nds){
		cout << "(" << (*ub)[0] << ", " << (*ub)[1] << ")" << endl;
        aux.insert(ub);
    }

    py_Plotter::offline_plot(aux, NULL, "output2.txt");

    while(  getchar() ){
        nds.pop_front();
        aux.clear();
        for(auto ub : nds)   aux.insert(ub);
        py_Plotter::offline_plot(aux, NULL, "output2.txt");
    }


    


	return 0;


}