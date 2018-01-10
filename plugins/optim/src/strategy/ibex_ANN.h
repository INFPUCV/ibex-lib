//============================================================================
//                                  I B E X                                   
// File        : CellDatas
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : May 10, 2012
// Last Update : May 10, 2012
//============================================================================

#ifndef __IBEX_ANN_H__
#define __IBEX_ANN_H__

#include <vector>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>

#include <ibex_TrainingData.h>


using namespace std;

namespace ibex {

	class ANN {

	public:
		ANN(const string filename);
		TrainingData trainData;

	private:
		vector<unsigned> topology;

		void showVectorVals(string label, vector<double> &v);

};

} // end namespace ibex

#endif // __IBEX_ANN_H__
