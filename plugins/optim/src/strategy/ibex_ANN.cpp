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

	// if(trainingPass > 2000) return testingNeuron(inputVals, targetVals);

	cout << "Training" << endl;

	++trainingPass;
	cout << "Pass " << trainingPass;

	// Get new input data and feed it forward:
	// if(trainData.getNextInputs(inputVals) != topology[0])
	// 	return resultVals;
	if(inputVals.size() != topology[0])
			return resultVals;
	showVectorVals(": Inputs :", inputVals);
	myNet->feedForward(inputVals);

	// Collect the net's actual results:
	myNet->getResults(resultVals);
	cout << "results size " << resultVals.size() << endl;
	showVectorVals("Outputs:", resultVals);

	// Train the net what the outputs should have been:
	// trainData.getTargetOutputs(targetVals);
	showVectorVals("Targets:", targetVals);
	assert(targetVals.size() == topology.back());

	myNet->backProp(targetVals);

	// Report how well the training is working, average over recnet
	cout << "Net recent average error: "
		 << myNet->getRecentAverageError() << endl;

	return resultVals;
}


vector<double> ANN::testingNeuron(vector<double> &inputVals, vector<double> &targetVals) {
	vector<double> resultVals;
	++trainingPass;

	cout << "Testing" << endl;

	cout << "Pass " << trainingPass;

	// Get new input data and feed it forward:
	// if(trainData.getNextInputs(inputVals) != topology[0])
	// 	return resultVals;
	if(inputVals.size() != topology[0])
			return resultVals;
	showVectorVals(": Inputs :", inputVals);
	myNet->feedForward(inputVals);

	// Collect the net's actual results:
	myNet->getResults(resultVals);
	showVectorVals("Outputs:", resultVals);


	// Train the net what the outputs should have been:
	// trainData.getTargetOutputs(targetVals);
	showVectorVals("Targets:", targetVals);
	assert(targetVals.size() == topology.back());

	// Report how well the training is working, average over recnet
	cout << "Net recent average error: "
		 << myNet->getRecentAverageError() << endl;

	int auxTarget=0, auxResults=0;
	// any contract result
	if(targetVals[targetVals.size()-1] == 1) {
		for(int i=0; i< targetVals.size();i++) {
			if(resultVals[i] > 0.5) auxResults++;
		}
		if(resultVals[targetVals.size()-1] > 0.5 || auxResults > 0) VAC++;
		else FAC++;
	} else {
		for(int i=0; i< targetVals.size();i++) {
			if(targetVals[i] > 0.5) auxTarget++;
		}
		// no contract
		if(auxTarget == 0) {
			for(int i=0; i< resultVals.size();i++) {
				if(resultVals[i] > 0.7) auxResults++;
			}
			if(auxResults == 0) VNC++;
			else FNC++;
		// contract
		} else {
			int i = 0;
			for(i=0; i< resultVals.size();i++) {
				if(resultVals[i] < 0.5 && targetVals[i] == 1) break;
			}
			if(i == resultVals.size()) VC++;
			else FC++;

		}
	}

	cout << "VNC " << VNC << " FNC " << FNC << endl;
	cout << "VC " << VC << " FC " << FC << endl;
	cout << "VAC " << VAC << " FAC " << FAC << endl;
	float verdadero, falso;
	verdadero = (float) (VNC+VC+VAC)/(VNC+VC+VAC+FNC+FC+FAC);
	falso = (float) (FNC+FC+FAC)/(VNC+VC+VAC+FNC+FC+FAC);
	cout << "verdadero " << verdadero << " falso " << falso << endl;
	// if(resultVals[targetVals.size()-1] > 0.5 || auxResults > 0) getchar();

	return resultVals;
}


ANN::ANN(const string filename) : trainData(filename), trainingPass(0), totalTraining(10) {

	// obtiene la topologia [cant inputs, cant nodo escondido, cant output, cantidad de datos]
	unsigned arr[] = {8, 5, 8};
	topology = vector<unsigned>(arr, arr+3);

	// trainData.getTopology(topology);

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
