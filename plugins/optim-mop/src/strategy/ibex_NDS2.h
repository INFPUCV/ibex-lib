/*
 * ibex_NDS2.h
 *
 *  Created on: 18 oct. 2019
 *      Author: javier
 */

#include "ibex_IntervalVector.h"
#include "ibex_pyPlotter.h"
#include <unordered_map>
#include <list>
#include "ibex_BxpMOPData.h"


#ifndef OPTIM_MOP_SRC_STRATEGY_IBEX_NDS2_H_
#define OPTIM_MOP_SRC_STRATEGY_IBEX_NDS2_H_

using namespace std;

namespace ibex {

class NDS_X{
public:
	NDS_X() : x1(NULL){ }

	NDS_X(const Vector& x1) : x1(new Vector(x1)) { }

	NDS_X(const NDS_X& other) {
		x1=NULL;
		if(other.x1) x1=new Vector(*other.x1);
	}

	NDS_X& operator=(NDS_X& other){
		if(this != &other){
			if(x1) delete x1;
			x1=NULL;
			if(other.x1) x1=new Vector(*other.x1);
		}
		return *this;
	}

	~NDS_X(){
		if(x1) delete x1;
	}
	Vector* x1;
};

class MyHashFunction{
public:
    std::size_t operator()(Vector const &vect) const{
        std::size_t seed;
        seed = *vect.raw();
        return seed;
    }
};

class NDS_XY{
public:
	NDS_XY() {
		clear();
	};

	virtual ~NDS_XY() { };

	/*
 	 * Clear Map NDS
 	 */
	void clear(){
		NDS.clear();
		if(nObjFunc != 0){
			Vector init_point(nObjFunc);
			for(int i=0; i<nObjFunc; i++) init_point[i] = POS_INFINITY;
			NDS.insert(make_pair(init_point, NDS_X()));
		}
	}

	/*
	 *Return Size of Map NDS
	 */
	int size() const{
		return NDS.size();
	}

	/*
	 * Check if new_point_y is dominated by NDS set
	 */
	bool is_dominated(const Vector& new_point_y);

	/*
	 * Check if point new_point_y is dominated by y
	 */
	bool is_dominated(const Vector& new_point_y, const Vector& y);

	/*
	 * Add Point y in set NDS and delete every point dominanted by y
	 */
	void addPoint(const Vector& y, const NDS_X& data=NDS_X());

	/*
	 * Return the minimum dimension of the maximum distance  of a point "y" to the "NDS set"
	 */
	double distance(const IntervalVector& y);

	/*
	 * Return list of the best cutting points for the Dominance Peeler
	 */
	list< pair <Vector, int> > cutting_points(const Vector& boxy_lb, const Vector& boxy_ub);

/*Testing functions
 * */
	void print_NDS(){
		cout << "NDS:" << endl;
		for(auto itr = NDS.begin(); itr != NDS.end(); ++itr) {
			for(int i= 0; i<itr->first.size(); i++ )
				cout <<itr->first[i] << " ";
			cout << endl;
			
			//cout<<" x= ";
			//cout<<*itr->second.x1;
			//cout<<"\n";
		}

	}

	void print_NDS_y(){
		for(auto itr = NDS.begin(); itr != NDS.end(); ++itr) {
			cout<<"Punto no dominado= "<<itr->first<<endl;
		}
	}


	unordered_map< Vector, NDS_X, MyHashFunction > NDS;
	int nObjFunc = 0;
};


}
#endif /* PLUGINS_OPTIM_MOP_SRC_STRATEGY_IBEX_NDS2_H_ */

