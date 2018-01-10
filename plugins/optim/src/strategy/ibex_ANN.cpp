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


#include "ibex_Network.h"

namespace ibex {

void ANN::showVectorVals(string label, vector<double> &v)
{
	cout << label << " ";
	for(unsigned i = 0; i < v.size(); ++i)
	{
		cout << v[i] << " ";
	}
	cout << endl;
}

ANN::ANN(const string filename) : trainData(filename) {

	// obtiene la topologia [cant inputs, cant nodo escondido, cant output, cantidad de datos]
	trainData.getTopology(topology);
	cout << "topology: " << topology[3] << endl;

	// crea la red neuronal con la topologia obtenida
	Network myNet(topology);

	vector<double> inputVals, targetVals, resultVals;
	int trainingPass = 0;
	int iter = 0;

	while(!trainData.isEof())
		{
			++trainingPass;
			cout << endl << "Pass" << trainingPass;

			// Get new input data and feed it forward:
			if(trainData.getNextInputs(inputVals) != topology[0])
				break;
			showVectorVals(": Inputs :", inputVals);
			myNet.feedForward(inputVals);

			// Collect the net's actual results:
			myNet.getResults(resultVals);
			showVectorVals("Outputs:", resultVals);

			// Train the net what the outputs should have been:
			trainData.getTargetOutputs(targetVals);
			showVectorVals("Targets:", targetVals);
			assert(targetVals.size() == topology.back());

			myNet.backProp(targetVals);

			// Report how well the training is working, average over recnet
			cout << "Net recent average error: "
			     << myNet.getRecentAverageError() << endl;

			iter++;
		}

		// get error ANN
		cout << "Net recent average error: "
	                     << myNet.getRecentAverageError() << endl;


}


} // end namespace ibex
