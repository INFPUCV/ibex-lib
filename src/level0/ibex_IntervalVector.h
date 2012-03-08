/* ============================================================================
 * I B E X - Interval Vector definition
 * ============================================================================
 * Copyright   : Ecole des Mines de Nantes (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Gilles Chabert
 * Created     : Dec 05, 2011
 * ---------------------------------------------------------------------------- */

#ifndef __IBEX_INTERVAL_VECTOR_H__
#define __IBEX_INTERVAL_VECTOR_H__

#include <cassert>
#include <iostream>
#include "ibex_Interval.h"
#include "ibex_InvalidIntervalVectorOp.h"

namespace ibex {

class Domain;
class IntervalMatrix; // declared only for friendship

/**
 * \ingroup arith
 *
 * \brief Vector of Intervals
 *
 * By convention an empty vector has a dimension. A vector becomes empty
 * when one of its component becomes empty and all the components
 * are set to the empty Interval.
 */
class IntervalVector {

public:
	/**
	 * \brief Create [(-oo,+oo) ; ..; (-oo,+oo)]
	 *
	 * Create a n-sized vector. All the components are (-oo,+oo)
	 * \pre n>0
	 */
	IntervalVector(int n);

	/**
	 * \brief Create [x; ....; x]
	 *
	 * Create a IntervalVector of dimension \a n with
	 * all the components initialized to \a x.
	 * \pre n>0
	 */
	IntervalVector(int n, const Interval& x);

	/**
	 * \brief Create a copy of \a x.
	 */
	IntervalVector(const IntervalVector& x);

	/**
	 * \brief Create the IntervalVector [bounds[0][0],bounds[0][1]]x...x[bounds[n-1][0],bounds[n-1][1]]
	 *
	 * \param bounds an nx2 array of doubles
	 * \pre n>0
	 */
	IntervalVector(int n, double  bounds[][2]);

	/**
	 * \brief Create [empty; ...; empty]
	 *
	 * Create an empty IntervalVector of dimension \a n
	 * (all the components being empty Intervals)
	 *
	 * \pre n>0
	 */
	static IntervalVector empty(int n);

	/**
	 * \brief Delete this vector
	 */
	virtual ~IntervalVector();

	/**
	 * \brief Return the ith Interval
	 *
	 * A return a const reference to the
	 * i^th component (i starts from 0)
	 */
	const Interval& operator[](int i) const;

	/**
	 * \brief Return the ith Interval
	 *
	 * A return a non-const reference to the
	 * i^th component (i starts from 0)
	 */
	Interval& operator[](int i);

	/**
	 * \brief Set this IntervalVector to the empty IntervalVector
	 *
	 * The dimension remains the same.
	 */
	void set_empty();

	/**
	 * \brief Resize this IntervalVector.
	 *
	 * If the size is increased, the existing components are not
	 * modified and the new ones are set to (-inf,+inf), even if
	 * (*this) is the empty Interval (however, in this case, the status of
	 * (*this) remains "empty").
	 */
	void resize(int n2);

	/**
	 * \brief Return a subvector.
	 *
	 * \pre (*this) must not be empty
	 * \return [ (*this)[start_index]; ...; (*this)[end_index] ].
	 */
	IntervalVector subvector(int start_index, int end_index);

	/*
	 * \brief Set the lower bound of the ith component of this
	 *
	 * Set the bound unless the IntervalVector was empty in which case the ith component remains
	 * the empty Interval.
	 */
	//void set_lb(int i, double lb);

	/*
	 * \brief Set the upper bound of the ith component of this
	 *
	 * Set the bound unless the IntervalVector was empty in which case the ith component remains
	 * the empty Interval.
	 */
	//void set_ub(int i, double ub);

	/**
	 * \brief Assign the ith component to [lb,ub]
	 *
	 * Assign unless the IntervalVector was empty in which case the ith component remains
	 * the empty Interval.
	 */
	//void set(int i, double lb, double ub);

	/** \brief Assign the ith component to x:
	 *
	 * Assignment unless the IntervalVector was empty in which case the ith component remains
	 * the empty Interval. If x is the empty Interval, the IntervalVector becomes empty.
	 */
	//void set(int i, const Interval& x);

	/**
	 * \brief Assign this IntervalVector to x.
	 *
	 * \pre Dimensions of this and x must match. 
	 * \note Emptiness is overridden.
	 */
	IntervalVector& operator=(const IntervalVector& x);

	/**
	 * \brief "Load" the domain into this
	 */
	IntervalVector& operator=(const Domain& d);

	/**
	 * \brief Set *this to its intersection with x
	 *
	 * \return a reference to this.
	 * \throws InvalidIntervalVectorOp if the vectors do not have the same dimensions.
	 */
	IntervalVector& operator&=(const IntervalVector& x);

	/**
	 * \brief Set this IntervalVector to the hull of itself and another.
	 *
	 * \return a reference to this.
	 * \throws InvalidIntervalVectorOp if IntervalVectores do not have the same dimensions.
	 */
	IntervalVector& operator|=(const IntervalVector& x);

	/**
	 * \brief Return the intersection of this and x.
	 */
	IntervalVector operator&(const IntervalVector& x) const;

	/**
	 * \brief Return the hull of this & x.
	 */
	IntervalVector operator|(const IntervalVector& x) const;

	/**
	 * \brief Return true if the bounds of this IntervalVector match that of \a x.
	 */
	bool operator==(const IntervalVector& x) const;

	/**
	 * \brief Return true if one bounds of one component of *this differs from \a x.
	 */
	bool operator!=(const IntervalVector& x) const;

	/**
	 * \brief The dimension (number of components)
	 */
	int size() const;

	/**
	 * \brief Return the midpoint (a degenerated IntervalVector)
	 * \pre (*this) must be nonempty
	 */
	IntervalVector mid() const;

	/**
	 * \brief Return true iff this IntervalVector is empty
	 */
	bool is_empty() const;

	/**
	 * \brief Return true iff this IntervalVector is flat.
	 *
	 * An IntervalVector is "flat" if the radius is 0 on at least one dimension
	 * An empty interval vector is considered as flat.
	 */
	bool is_flat() const;

	/**
	 * \brief True iff this interval vector contains \a x.
	 *
	 * \pre Dimension of \a x must be equal to the dimension of (*this).
	 * \sa #ibex::Interval::contains(double) const.
	 */
	//bool contains(const Vector& x) const;

	/**
	 * \brief true iff this interval vector contains an infinite bound.
	 *
	 * \note An empty interval vector is always bounded.
	 */
	bool is_unbounded() const;

	/**
	 * \brief True iff this interval vector is a subset of \a x.
	 *
	 * \pre Dimension of \a x must be equal to the dimension of this vector.

	 * \note Always return true if this interval vector is empty.

	 * \sa #ibex::Interval::is_subset(const Interval&) const.
	 */
	bool is_subset(const IntervalVector& x) const;

	/**
	 * \brief True iff this interval vector is inside the interior of \a x.
	 *
	 * \pre Dimension of \a x must be equal to the dimension of this vector.
	 *
	 * \note return true if this interval vector is empty and \a x not.
	 *
	 * \sa #ibex::Interval::is_strict_subset(const Interval&) const.
	 */
	bool is_strict_subset(const IntervalVector& x) const;

	/**
	 * \brief True iff this interval vector is a superset of \a x.
	 *
	 * \pre Dimension of \a x must be equal to the dimension of this vector.

	 * \note Always return true if \a x is empty.

	 * \sa #ibex::Interval::is_superset(const Interval&) const.
	 */
	bool is_superset(const IntervalVector& x) const;

	/**
	 * \brief True iff \a x is inside the interior of (*this).
	 *
	 * \pre Dimension of \a x must be equal to the dimension of this vector.
	 *
	 * \note return true if x is empty and not (*this).
	 *
	 * \sa #ibex::Interval::is_strict_superset(const Interval&) const.
	 */
	bool is_strict_superset(const IntervalVector& x) const;

	/**
	 * \brief Return the index of a component with minimal/maximal diameter.
	 *
	 *  \param min true => minimal diameter
	 *  \throws InvalidIntervalVectorOp if the IntervalVector is empty.
	 */
	int extr_diam_index(bool min) const;

	/**
	 * \brief Return the maximal diameter among all the components.
	 *
	 *  \throws InvalidIntervalVectorOp if the IntervalVector is empty.
	 */
	double max_diam() const;

	/**
	 * \brief Return the minimal diameter among all the components.
	 *
	 * \throws InvalidIntervalVectorOp if the IntervalVector is empty.
	 */
	double min_diam() const;

	/**
	 * \brief Return the volume of this interval vector.
	 *
	 * \note Return \c POS_INFINITY if the vector is unbounded and not flat.
	 * \note Return 0 if the vector is flat and not unbounded.
	 * \warning If the interval vector is both flat and unbounded, the result is undefined.
	 * \sa #flat()
	 * \sa #unbounded()
	 */
	double volume() const;

	/**
	 * \brief Return the perimeter of this interval vector.
	 *
	 * \note Return \c POS_INFINITY if unbounded.
	 */
	double perimeter() const;

	/**
	 * \brief Return the relative distance with x.
	 *
	 * \return \f$\displaystyle \max_{i=1..n} rel_distance([this]_i, x_i)/diam([this]_i)\f$.
	 *
	 * \sa #ibex::distance(const IntervalVector& x1, const IntervalVector& x2).
	 * \sa #Interval::rel_distance(const Interval& x) const.
	 */
	double rel_distance(const IntervalVector& x) const;

	/**
	 * \brief The opposite of this.
	 */
	IntervalVector operator-() const;

	/**
	 * \brief Sum of this and x.
	 */
	IntervalVector operator+(const IntervalVector& x) const;

	/**
	 * \brief Add x to this.
	 */
	IntervalVector& operator+=(const IntervalVector& x);

	/**
	 * \brief Subtraction of this and x.
	 */
	IntervalVector operator-(const IntervalVector& x) const;

	/**
	 * \brief Subtract x to this.
	 */
	IntervalVector& operator-=(const IntervalVector& x);

	/**
	 * \brief Return *this \ y (set difference).
	 *
	 * Store the difference under the form of a union of non-overlapping IntervalVectors
	 * into \a result, and return the size of the union.
	 *
	 * <p>
	 * If the difference is empty, \a result is an array of one element set to the empty box.
	 * It is <b>not</b> a zero-sized array containing no element.
	 */
	int diff(const IntervalVector& y, IntervalVector*& result) const;

	/**
	 * \brief Return the complementary of *this.
	 *
	 * Store the complementary under the form of a union of non-overlapping IntervalVectors,
	 * into \a result, and return the size of the union.
	 *
	 * <p>
	 * If (*this) is the empty set with n components, the complementary of (*this) is
	 * the n-dimensional box (-oo,oo)x...(-oo,oo).
	 *
	 * <p>
	 * If the complementary is empty, \a result is an array of one element set to the empty box.
	 * It is <b>not</b> a zero-sized array containing no element.
	 */
	int complementary(IntervalVector*& result) const;

	/**
	 * \brief Return a random vector (degenerated box) inside *this.
	 *
	 * \pre (*this) must be nonempty.
	 */
	IntervalVector random() const;

private:
	friend class IntervalMatrix;

	IntervalVector() : n(0), vec(NULL) { } // for IntervalMatrix & complementary()

	int n;             // dimension (size of vec)
	Interval *vec;	   // vector of elements
};

/**
 * \brief Return the infinite Hausdorff distance (i.e. the maximum of the distance componentwise).
 *
 * \return \f$\displaystyle \max_{i=1..n} distance(x1_i, x2_i)\f$.
 *
 * \pre Dimension of \a x1 and \a x2 must be equal.
 *
 * \sa #ibex::distance(const Interval&, const Interval&) const.
 */
double distance(const IntervalVector& x1, const IntervalVector& x2);


/*================================== inline implementations ========================================*/

inline IntervalVector::IntervalVector(int n) : n(n), vec(new Interval[n]) {
	assert(n>=1);
	for (int i=0; i<n; i++) vec[i]=Interval::ALL_REALS;
}

inline IntervalVector::IntervalVector(int n, const Interval& x) : n(n), vec(new Interval[n]) {
	assert(n>=1);
	for (int i=0; i<n; i++) vec[i]=x;
}

inline IntervalVector::IntervalVector(const IntervalVector& x) : n(x.n), vec(new Interval[x.n]) {
	for (int i=0; i<n; i++) vec[i]=x[i];
}

inline IntervalVector::IntervalVector(int n, double bounds[][2]) : n(n), vec(new Interval[n]) {
	for (int i=0; i<n; i++)
		vec[i]=Interval(bounds[i][0],bounds[i][1]);
}

inline IntervalVector IntervalVector::empty(int n) {
	return IntervalVector(n, Interval::EMPTY_SET);
}

inline IntervalVector::~IntervalVector() {
	delete[] vec;
}

inline const Interval& IntervalVector::operator[](int i) const {
	assert(i>=0 && i<n);
	return vec[i];
}

inline Interval& IntervalVector::operator[](int i) {
	assert(i>=0 && i<n);
	return vec[i];
}

inline void IntervalVector::set_empty() {
	// warning: do not insert this test:
	//	 if (isEmpty()) return;
	// because we call set_empty() from set(...) and the first component
	// may be empty in an intermediate state

	for (int i=0; i<size(); i++)
		(*this)[i]=Interval::EMPTY_SET;
}

inline void IntervalVector::resize(int n2) {
	assert(n2>=1);
	if (n2==size()) return;

	Interval* newVec=new Interval[n2];
	int i=0;
	for (; i<size() && i<n2; i++)
		newVec[i]=vec[i];
	for (; i<n2; i++)
		newVec[i]=Interval::ALL_REALS;
	if (vec!=NULL) // vec==NULL happens when default constructor is used (n==0)
		delete[] vec;

	n   = n2;
	vec = newVec;
}

inline IntervalVector IntervalVector::subvector(int start_index, int end_index) {
	assert(!is_empty());
	assert(end_index>=0 && start_index>=0);

	if (end_index>=size() || start_index>end_index)
		throw NonRecoverableException("Invalid indices for IntervalVector::subvector");

	IntervalVector v(end_index-start_index+1);
	int j=0;
	for (int i=start_index; i<=end_index; i++) {
		v[j++]=(*this)[i];
	}
	return v;
}

inline IntervalVector& IntervalVector::operator=(const IntervalVector& x) {
	assert(size()==x.size()); // throw InvalidIntervalVectorOp("Cannot set a IntervalVector to a IntervalVector with different dimension");

	if (x.is_empty())
		set_empty();
	else
		// don't use "set(...)" because the test "is_empty()" called inside
		// may return prematurely in case "this" is empty.
		// use physical copy instead:
		for (int i=0; i<size(); i++)
			(*this)[i]=x[i];

	return *this;
}

inline IntervalVector& IntervalVector::operator&=(const IntervalVector& x)  {
	// dimensions are non zero henceforth
	if (size()!=x.size()) throw InvalidIntervalVectorOp("Cannot intersect IntervalVectores with different dimensions");

	if (is_empty()) return *this;
	if (x.is_empty()) { set_empty(); return *this; }

	for (int i=0; i<size(); i++) {
		(*this)[i] &= x[i];
		if ((*this)[i].is_empty()) {
			set_empty();
			return *this;
		}
	}
	return *this;
}

inline IntervalVector& IntervalVector::operator|=(const IntervalVector& x)  {
	// dimensions are non zero henceforth
	if (size()!=x.size()) throw InvalidIntervalVectorOp("Cannot make the hull of IntervalVectores with different dimensions");

	if (x.is_empty()) return *this;
	if (is_empty()) { *this=x; return *this; }

	for (int i=0; i<size(); i++) {
		(*this)[i] |= x[i];
	}
	return *this;
}

inline IntervalVector IntervalVector::operator&(const IntervalVector& x) const {
	return IntervalVector(*this) &= x;
}

inline IntervalVector IntervalVector::operator|(const IntervalVector& x) const {
	return IntervalVector(*this) |= x;
}

inline bool IntervalVector::operator==(const IntervalVector& x) const {
	if (n!=x.size()) return false;
	if (is_empty() || x.is_empty()) return is_empty() && x.is_empty();
	for (int i=0; i<n; i++)
		if ((*this)[i]!=(x[i])) return false;
	return true;
}

inline bool IntervalVector::operator!=(const IntervalVector& x) const {
	return !(*this==x);
}

inline int IntervalVector::size() const {
	return n;
}

inline IntervalVector IntervalVector::mid() const {
	assert(!is_empty());

	IntervalVector mV(size());
	for (int i=0; i<size(); i++) {
		Interval m =  (*this)[i].mid();
		mV[i]=m;
	}
	return mV;
}

inline bool IntervalVector::is_empty() const {
	return (*this)[0].is_empty();
}

inline bool IntervalVector::is_flat() const {
	if (is_empty()) return true;
	for (int i=0; i<size(); i++)
		if ((*this)[i].is_degenerated()) // don't use diam() because of roundoff
			return true;
	return false;
}

/*
inline bool IntervalVector::contains(const VECTOR& x) const {
  for (int i=1; i<=size(); i++)
    if (!(*this)(i).contains(x(i))) return false;
  return true;
}*/

inline bool IntervalVector::is_unbounded() const {
  if (is_empty()) return false;
  for (int i=0; i<size(); i++)
    if ((*this)[i].is_unbounded()) return true;
  return false;
}

inline bool IntervalVector::is_subset(const IntervalVector& x) const {
  if (is_empty()) return true;
  for (int i=0; i<size(); i++)
    if (!(*this)[i].is_subset(x[i])) return false;
  return true;
}

inline bool IntervalVector::is_strict_subset(const IntervalVector& x) const {
  if (is_empty() && !x.is_empty()) return true;
  if (x.is_empty()) return false;
  for (int i=0; i<size(); i++)
    if (!(*this)[i].is_strict_subset(x[i])) return false;
  return true;
}

inline bool IntervalVector::is_superset(const IntervalVector& x) const {
	return x.is_subset(*this);
}

inline bool IntervalVector::is_strict_superset(const IntervalVector& x) const {
	return x.is_strict_subset(*this);
}

inline int IntervalVector::extr_diam_index(bool min) const {
	double d=(*this)[0].diam();
	int selectedIndex=0;
	if (is_empty()) throw InvalidIntervalVectorOp("Diameter of an empty IntervalVector is undefined");
	for (int i=1; i<size(); i++) {
		double w=(*this)[i].diam();
		if (min? w<d : w>d) {
			selectedIndex=i;
			d=w;
		}
	}
	return selectedIndex;
}

inline double IntervalVector::max_diam() const {
	return (*this)[extr_diam_index(false)].diam();
}

inline double IntervalVector::min_diam() const {
	return (*this)[extr_diam_index(true)].diam();
}

inline IntervalVector IntervalVector::operator-() const {
	const int n=size();
	if (is_empty()) return empty(n);

	IntervalVector x(n);
	for (int i=0; i<n; i++) {
		x[i]=-(*this)[i];
	}
	return x;
}

inline IntervalVector IntervalVector::operator+(const IntervalVector& x) const {
	return IntervalVector(*this)+=x;
	/*const int n=size();
	assert(x.size()==n);
	if (is_empty() || x.is_empty()) return empty(n);
	else {
		IntervalVector y(n);
		for (int i=0; i<n; i++) {
			(y[i]=(*this)[i])+=x[i]; // faster than y[i]=(*this)[i]+x[i]
		}
		return y;
	}*/
}

inline IntervalVector& IntervalVector::operator+=(const IntervalVector& x) {
	const int n=size();
	assert(x.size()==n);
	if (!is_empty()) {
		if (x.is_empty()) { set_empty(); }
		else {
			for (int i=0; i<n; i++) {
				(*this)[i]+=x[i];
			}
		}
	}
	return *this;
}

inline IntervalVector IntervalVector::operator-(const IntervalVector& x) const {
	return IntervalVector(*this)-=x;
}

inline IntervalVector& IntervalVector::operator-=(const IntervalVector& x) {
	const int n=size();
	assert(x.size()==n);
	if (!is_empty()) {
		if (x.is_empty()) { set_empty(); }
		else {
			for (int i=0; i<n; i++) {
				(*this)[i]-=x[i];
			}
		}
	}
	return *this;
}

/** \brief Display the IntervalVector \a x
 */
inline std::ostream& operator<<(std::ostream& os, const IntervalVector& x) {
	os << "(";
	for (int i=0; i<x.size(); i++)
		os << x[i] << (i<x.size()-1? " ; " : "");
	os << ")";
	return os;
}

/*
inline void IntervalVector::set_lb(int i, double lb) {
	if (!is_empty()) {
		assert(lb<=vec[i].ub());
		vec[i]=Interval(lb, vec[i].ub());
	}
}

inline void IntervalVector::set_ub(int i, double ub) {
	if (!is_empty()) {
		assert(ub>=vec[i].lb());
		vec[i]=Interval(vec[i].lb(),ub);
	}
}

inline void IntervalVector::set(int i, double lb, double ub) {
	if (!is_empty()) {
		assert(ub>=lb);
		(*this)[i]=Interval(lb,ub);
	}
}


inline void IntervalVector::set(int i, const Interval& x) {
	if (x.is_empty()) set_empty();
	else set(i, x.lb(), x.ub());
}*/


} // end namespace

#endif /* _IBEX_INTERVAL_VECTOR_H_ */
