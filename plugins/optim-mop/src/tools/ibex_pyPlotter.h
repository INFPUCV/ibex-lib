/*
 * ibex_pyPlotter.h
 *
 *  Created on: 10 ene. 2018
 *      Author: iaraya
 */

#ifndef OPTIM_MOP_SRC_TOOLS_IBEX_PYPLOTTER_H_
#define OPTIM_MOP_SRC_TOOLS_IBEX_PYPLOTTER_H_

#include <set>
#include <map>
#include "ibex_BxpMOPData.h"
#include "ibex_Cell.h"

using namespace std;
namespace ibex {

//TODO: more general!

/**
 * \brief Plotter for python3 .
 *
 * This class implements some methods for plotting the results on a
 * python3 program (plot.py and plot2.py)
 */

class py_Plotter {
public:

	/**
	 * \brief writes a file (output.txt) to be read by the python3 program plot.py
	 */
	static void offline_plot(map< pair <double, double>, IntervalVector, struct sorty2 >& NDS,
		map< pair <double, double>, IntervalVector, struct sorty2 >* NDS2, const char* output_file);

};

} /* namespace ibex */

#endif /* OPTIM_MOP_SRC_TOOLS_IBEX_PYPLOTTER_H_ */
