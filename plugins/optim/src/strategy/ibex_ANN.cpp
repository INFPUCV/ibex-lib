//============================================================================
//                                  I B E X                                   
// File        : ibex_Cell.cpp
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : May 10, 2012
// Last Update : May 10, 2012
//============================================================================

#include "ibex_ANN.h"

#include <vector>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>

namespace ibex {

ANN::ANN(const string filename) : trainData(filename) {

	// obtiene la topologia [cant inputs, cant nodo escondido, cant output, cantidad de datos]
	trainData.getTopology(topology);
	cout << "topology: " << topology[3] << endl;
}


} // end namespace ibex
