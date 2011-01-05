// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/10/2009-- carved out of CbcBranchActual

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "CoinTypes.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcDummyBranchingObject.hpp"
#include "CbcBranchActual.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

//##############################################################################

// Default Constructor
CbcDummyBranchingObject::CbcDummyBranchingObject(CbcModel * model)
        : CbcBranchingObject(model, 0, 0, 0.5)
{
    setNumberBranchesLeft(1);
}


// Copy constructor
CbcDummyBranchingObject::CbcDummyBranchingObject ( const CbcDummyBranchingObject & rhs) : CbcBranchingObject(rhs)
{
}

// Assignment operator
CbcDummyBranchingObject &
CbcDummyBranchingObject::operator=( const CbcDummyBranchingObject & rhs)
{
    if (this != &rhs) {
        CbcBranchingObject::operator=(rhs);
    }
    return *this;
}
CbcBranchingObject *
CbcDummyBranchingObject::clone() const
{
    return (new CbcDummyBranchingObject(*this));
}


// Destructor
CbcDummyBranchingObject::~CbcDummyBranchingObject ()
{
}

/*
  Perform a dummy branch
*/
double
CbcDummyBranchingObject::branch()
{
    decrementNumberBranchesLeft();
    return 0.0;
}
// Print what would happen
void
CbcDummyBranchingObject::print()
{
    printf("Dummy branch\n");
}
/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type.
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int
CbcDummyBranchingObject::compareOriginalObject
(const CbcBranchingObject* /*brObj*/) const
{
    throw("must implement");
}

/** Compare the \c this with \c brObj. \c this and \c brObj must be os the
    same type and must have the same original object, but they may have
    different feasible regions.
    Return the appropriate CbcRangeCompare value (first argument being the
    sub/superset if that's the case). In case of overlap (and if \c
    replaceIfOverlap is true) replace the current branching object with one
    whose feasible region is the overlap.
*/
CbcRangeCompare
CbcDummyBranchingObject::compareBranchingObject
(const CbcBranchingObject* /*brObj*/, const bool /*replaceIfOverlap*/)
{
    throw("must implement");
}

