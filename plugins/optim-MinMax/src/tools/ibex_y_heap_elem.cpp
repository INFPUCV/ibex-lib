#include "ibex_y_heap_elem.h"


y_heap_elem::y_heap_elem(IntervalVector box,Interval pf,double  pu):box(box),pf(pf),pu(pu)
{};

y_heap_elem::y_heap_elem(const y_heap_elem& original):box(original.box),pf(original.pf),pu(original.pu)
{};


pair<y_heap_elem*,y_heap_elem*> y_heap_elem::bisect(IntervalVector box1,IntervalVector box2) {
    y_heap_elem *elem1 = new y_heap_elem(*this);
    y_heap_elem *elem2 = new y_heap_elem(*this);
    elem1->box = box1;
    elem2->box = box2;
    pair<y_heap_elem*,y_heap_elem*> p;
    p.first = elem1;
    p.second = elem2;
    return p;

}


double y_heap_costfub::cost(const y_heap_elem& elem) const {
    return  -elem.pf.ub();
}