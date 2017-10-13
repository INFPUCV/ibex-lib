//============================================================================
//                                  I B E X
// File        : ibex_CellNSSet.cpp
// Author      : Matias Campusano, Damir Aliquintui, Ignacio Araya
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Oct 12, 2017
// Last Update : Oct 12, 2017
//============================================================================

#include "ibex_CellNSSet.h"

using namespace std;

namespace ibex {


	CellNSSet::CellNSSet() {

		}

	CellNSSet::~CellNSSet() {
		flush();
	}

	void CellNSSet::flush() {
		while (!dset.empty()) {
			delete *dset.begin();
			dset.erase(dset.begin());
		}
		while (!nondset.empty()) {
			delete *nondset.begin();
			nondset.erase(nondset.begin());
		}
	}

	unsigned int CellNSSet::size() const {
		return (dset.size()+nondset.size());
	}

	bool CellNSSet::empty() const {
		return nondset.empty();
	}

	void CellNSSet::push(Cell* cell) {
		// if (capacity>0 && size() == capacity) throw CellBufferOverflow();
		// TODO Verficar que n-1 y n-2 sean las cajas de las funciones objetivo
		bool dominated = false;
		int n;
		std::multiset<Cell*, maxsize>::iterator it;
		std::multiset<Cell*, maxsize>::iterator aux;
		for (it = nondset.begin(); it != nondset.end(); ) {
			n = (*it)->box.size();
			// si cell es dominado por un no dominado
			if((*it)->box[n-1].lb() < cell->box[n-1].lb() && (*it)->box[n-2].lb() < cell->box[n-2].lb()) {
				dominated = true;
				break;
			}
			// si cell domina a un no dominado
			if(cell->box[n-1].lb() < (*it)->box[n-1].lb() && cell->box[n-2].lb() < (*it)->box[n-2].lb()) {
				dset.push_front(*it);
				aux = it;
				++aux;
				nondset.erase(*it);
				it=aux;
			} else ++it;
		}
		if(dominated) dset.push_front(cell);
		else nondset.insert(cell);
	}

	Cell* CellNSSet::pop() {
		int n,n1;
		bool dominated = false;
		typename std::list<Cell*> domcells;
		std::list<Cell*>::iterator it;
		std::list<Cell*>::iterator it1;
		std::list<Cell*>::iterator aux;
		std::multiset<Cell*, maxsize>::iterator it2;
		std::multiset<Cell*, maxsize>::iterator aux2;
		Cell* c = *nondset.begin();
		std::cout << "NS pop " << c  << std::endl;
		nondset.erase(nondset.begin());
		// se agregan todos los cell dominados por c en dset a la lista domcells
		for (it = dset.begin(); it != dset.end();) {
			n = (*it)->box.size();
			if(c->box[n-1].lb() < (*it)->box[n-1].lb() && c->box[n-2].lb() < (*it)->box[n-2].lb()) {
				domcells.push_front(*it);
				aux = it;
				++aux;
				dset.erase(it);
				it=aux;
			} else ++it;
		}
		// se cambian todos los cell que esten dominados en domcells a la lista dset
		for (it = domcells.begin(); it != domcells.end(); ) {
			n = (*it)->box.size();
			for (it1 = domcells.begin(); it1 != domcells.end(); ++it1) {
				n1 = (*it1)->box.size();
				if((*it1)->box[n1-1].lb() < (*it)->box[n-1].lb() && (*it1)->box[n1-2].lb() < (*it)->box[n-2].lb()) {
					dset.push_front(*it);
					aux = it;
					++aux;
					domcells.erase(it);
					it=aux;
				} else ++it;
			}
		}
		// se agrega a nondom si un elemento de domcells no es dominado por ningun elemento de nonset
		for (it = domcells.begin(); it != domcells.end(); ++it) {
			for (it2 = nondset.begin(); it2 != nondset.end(); ) {
				// si es dominado el elemento de domcells
				if((*it2)->box[n-1].lb() < (*it)->box[n-1].lb() && (*it2)->box[n-2].lb() < (*it)->box[n-2].lb()) {
					dominated = true;
					break;
				}
				// si domina el elemento de domcells a uno de nondset
				if((*it)->box[n-1].lb() < (*it2)->box[n-1].lb() && (*it)->box[n-2].lb() < (*it2)->box[n-2].lb()) {
					dset.push_front(*it2);
					aux2 = it2;
					++aux2;
					nondset.erase(it2);
					it2=aux2;
				}else ++it2;
			}
			if(dominated) dset.push_back(*it);
			else nondset.insert(*it);


		}
		return c;
	}

	Cell* CellNSSet::top() const{
		return *nondset.begin();
	}

	double CellNSSet::minimum() const{
		// TODO: que realizar
		return POS_INFINITY;
	}

	void CellNSSet::contract(double loup){
		// TODO: que realizar
	}

} // end namespace ibex
