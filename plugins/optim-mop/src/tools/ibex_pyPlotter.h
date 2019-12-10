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
#include <unordered_map>
#include "ibex_BxpMOPData.h"
#include "ibex_Cell.h"
#include "ibex_NDS.h"
#include "ibex_NDS2.h"

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
	static void offline_plot(map< Vector, NDS_data, struct sorty2 >& NDS,
		map< Vector, NDS_data, struct sorty2 >* NDS2, const char* output_file, IntervalVector* focus=NULL);


	static void offline_plot(unordered_map< Vector, NDS_X, MyHashFunction >& NDS,
			unordered_map< Vector, NDS_X, MyHashFunction >* NDS2, const char* output_file_pareto, const char* output_file_solutions, IntervalVector* focus=NULL);
};

} /* namespace ibex */

#endif /* OPTIM_MOP_SRC_TOOLS_IBEX_PYPLOTTER_H_ */
