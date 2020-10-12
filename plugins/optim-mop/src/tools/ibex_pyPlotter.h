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
#include <list>

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




	static void timer_save(list<pair<double, pair<double, double>> >  t_checkdominance, list<pair<double, pair<double, double>> >  t_distance, list<pair<double, pair<double, double>> >  t_upperbounding);


	static void timer_save_component(const char* nombre, list<pair<double, pair<double, double>> > data);

	static void upper_envelope_save(list< pair<IntervalVector, Vector> >upper_envelope);
};

} /* namespace ibex */

#endif /* OPTIM_MOP_SRC_TOOLS_IBEX_PYPLOTTER_H_ */
