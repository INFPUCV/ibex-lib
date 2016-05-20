//============================================================================
//                                  I B E X                                   
// File        : ibex_CtrGenerator.h
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Jun 12, 2012
// Last Update : Jun 12, 2012
//============================================================================

#ifndef __IBEX_CTR_GENERATOR_H__
#define __IBEX_CTR_GENERATOR_H__

#include <vector>
#include <stack>
#include <utility>
#include "ibex_Scope.h"
#include "ibex_NumConstraint.h"
#include "ibex_P_Expr.h"
#include "ibex_Expr.h"
#include "ibex_Function.h"
#include "ibex_P_Source.h"

namespace ibex {

namespace parser {

class P_NumConstraint;
class P_OneConstraint;
class P_ConstraintList;
class P_ConstraintLoop;

class CtrGenerator {
public:
	void generate(std::stack<Scope>& scopes, const P_Source& source, std::vector<NumConstraint*>& dest);

	void visit(const P_NumConstraint& c);
	void visit(const P_OneConstraint& c);
	void visit(const P_ConstraintList& l);
	void visit(const P_ConstraintLoop& l);


protected:

	std::stack<Scope>& scopes;

	P_Source& source;

	std::vector<NumConstraint*>& ctrs;

};

} // end namespace parser
} // end namespace ibex

#endif // __IBEX_CTR_GENERATOR_H__
