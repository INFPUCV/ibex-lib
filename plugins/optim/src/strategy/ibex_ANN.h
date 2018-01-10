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


using namespace std;

namespace ibex {

	class TrainingData {

	public:
		TrainingData(string filename);

		bool isEof(void);

		void getTopology(vector<unsigned> &topology);

		// Returns the number of input values read from the file:
		unsigned getNextInputs(vector<double> &inputVals);
		unsigned getTargetOutputs(vector<double> &targetOutputVals);

	private:
		ifstream m_trainingDataFile;

};

} // end namespace ibex

#endif // __IBEX_ANN_H__
