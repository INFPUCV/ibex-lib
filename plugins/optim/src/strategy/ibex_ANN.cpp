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

vector<double> ANN::trainingNeuron(vector<double> &inputVals, vector<double> &targetVals) {
	vector<double> resultVals;

	++trainingPass;
	cout << endl << "Pass " << trainingPass;

	// Get new input data and feed it forward:
	if(trainData.getNextInputs(inputVals) != topology[0])
		return resultVals;
	showVectorVals(": Inputs :", inputVals);
	myNet->feedForward(inputVals);

	// Collect the net's actual results:
	myNet->getResults(resultVals);
	showVectorVals("Outputs:", resultVals);

	// Train the net what the outputs should have been:
	trainData.getTargetOutputs(targetVals);
	showVectorVals("Targets:", targetVals);
	assert(targetVals.size() == topology.back());

	myNet->backProp(targetVals);

	// Report how well the training is working, average over recnet
	cout << "Net recent average error: "
		 << myNet->getRecentAverageError() << endl;

	return resultVals;
}


vector<double> ANN::testingNeuron(vector<double> &inputVals) {
	vector<double> resultVals, targetVals;
	++trainingPass;
	cout << endl << "Pass " << trainingPass;

	// Get new input data and feed it forward:
	if(trainData.getNextInputs(inputVals) != topology[0])
		return resultVals;
	showVectorVals(": Inputs :", inputVals);
	myNet->feedForward(inputVals);

	// Collect the net's actual results:
	myNet->getResults(resultVals);
	showVectorVals("Outputs:", resultVals);


	// Train the net what the outputs should have been:
	trainData.getTargetOutputs(targetVals);
	showVectorVals("Targets:", targetVals);
	assert(targetVals.size() == topology.back());

	// Report how well the training is working, average over recnet
	cout << "Net recent average error: "
		 << myNet->getRecentAverageError() << endl;

	return resultVals;
}


ANN::ANN(const string filename) : trainData(filename), trainingPass(0), totalTraining(10) {

	// obtiene la topologia [cant inputs, cant nodo escondido, cant output, cantidad de datos]
	unsigned arr[] = {8, 3, 1};
	topology = vector<unsigned>(arr, arr+3);
	trainData.getTopology(topology);

	// crea la red neuronal con la topologia obtenida
	myNet = new Network(topology);

	/*
	vector<double> inputVals, targetVals, resultVals;

	while(!trainData.isEof())
		{
			++trainingPass;
			cout << endl << "Pass " << trainingPass;

			// Get new input data and feed it forward:
			if(trainData.getNextInputs(inputVals) != topology[0])
				break;
			showVectorVals(": Inputs :", inputVals);
			myNet->feedForward(inputVals);

			// Collect the net's actual results:
			myNet->getResults(resultVals);
			showVectorVals("Outputs:", resultVals);

			// Train the net what the outputs should have been:
			trainData.getTargetOutputs(targetVals);
			showVectorVals("Targets:", targetVals);
			assert(targetVals.size() == topology.back());

			myNet->backProp(targetVals);

			// Report how well the training is working, average over recnet
			cout << "Net recent average error: "
			     << myNet->getRecentAverageError() << endl;

			if(trainingPass > totalTraining) break;
		}

		// get error ANN
		cout << "Net recent average error: "
	                     << myNet->getRecentAverageError() << endl;


		cout << "testing data" << endl;

		while(!trainData.isEof())
		{
			++trainingPass;
			cout << endl << "Pass " << trainingPass;

			// Get new input data and feed it forward:
			if(trainData.getNextInputs(inputVals) != topology[0])
				break;
			showVectorVals(": Inputs :", inputVals);
			myNet->feedForward(inputVals);

			// Collect the net's actual results:
			myNet->getResults(resultVals);
			showVectorVals("Outputs:", resultVals);

			// Train the net what the outputs should have been:
			trainData.getTargetOutputs(targetVals);
			showVectorVals("Targets:", targetVals);
			assert(targetVals.size() == topology.back());

			// Report how well the training is working, average over recnet
			cout << "Net recent average error: "
				 << myNet->getRecentAverageError() << endl;

			if(trainingPass > totalTraining + 5) break;
		}
		*/

}


} // end namespace ibex
