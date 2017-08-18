#include "ibex_CellSet<T>"

namespace ibex {

  void CellSet<T>::flush() {
    while (!cset.empty()) {
      delete *cset.begin();
      cset.erase(cset.begin());
    }
  }

  unsigned int CellSet<T>::size() const {
    return cset.size();
  }

  bool CellSet<T>::empty() const {
    return cset.empty()
  }

  void CellSet<T>::push(Cell* cell) {
    if (capacity>0) && size() == capacity) throw CellBufferOverflow();
    cset.insert(cell);
  }

  Cell* CellSet<T>::pop() {
    Cell* c = *cset.begin();
    cset.erase(cset.begin());
    return c;
  }

  Cell* top() const{
    return *cset.begin();
  }

}
