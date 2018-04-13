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
#include "ibex_CellMOP.h"

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
	static void offline_plot(Cell* current, map< pair <double, double>, IntervalVector >& NDS);
	static void offline_plot(Cell* current, map< pair <double, double>, IntervalVector, struct sorty2 >& NDS);

	/**
	* write line commands to be interpreted by the python3 program plot.py
	*/
	static void plot_add_ub(pair<double, double> eval);
	static void plot_del_ub(pair<double, double> eval);
	static void plot_add_lb(Cell* c);
	static void plot_add_box(Cell* c);
	static void plot_del_box(Cell* c);



	static int n;

};

} /* namespace ibex */

#endif /* OPTIM_MOP_SRC_TOOLS_IBEX_PYPLOTTER_H_ */
