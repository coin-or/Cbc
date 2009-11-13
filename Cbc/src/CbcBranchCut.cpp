/* $Id$ */
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchCut.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"


/** Default Constructor

*/
CbcBranchCut::CbcBranchCut ()
        : CbcObject()
{
}

/* Constructor so model can be passed up
*/
CbcBranchCut::CbcBranchCut (CbcModel * model)
        : CbcObject(model)
{
}
// Copy constructor
CbcBranchCut::CbcBranchCut ( const CbcBranchCut & rhs)
        : CbcObject(rhs)

{
}

// Clone
CbcObject *
CbcBranchCut::clone() const
{
    return new CbcBranchCut(*this);
}

// Assignment operator
CbcBranchCut &
CbcBranchCut::operator=( const CbcBranchCut& /*rhs*/)
{
    return *this;
}

// Destructor
CbcBranchCut::~CbcBranchCut ()
{
}
double
CbcBranchCut::infeasibility(const OsiBranchingInformation * /*info*/,
                            int &preferredWay) const
{
    throw CoinError("Use of base class", "infeasibility", "CbcBranchCut");
    preferredWay = -1;
    return 0.0;
}

// This looks at solution and sets bounds to contain solution
/** More precisely: it first forces the variable within the existing
    bounds, and then tightens the bounds to fix the variable at the
    nearest integer value.
*/
void
CbcBranchCut::feasibleRegion()
{
}
/* Return true if branch created by object should fix variables
 */
bool
CbcBranchCut::boundBranch() const
{
    return false;
}
CbcBranchingObject *
CbcBranchCut::createCbcBranch(OsiSolverInterface * /*solver*/, const OsiBranchingInformation * /*info*/, int /*way*/)
{
    throw CoinError("Use of base class", "createCbcBranch", "CbcBranchCut");
    return new CbcCutBranchingObject();
}


/* Given valid solution (i.e. satisfied) and reduced costs etc
   returns a branching object which would give a new feasible
   point in direction reduced cost says would be cheaper.
   If no feasible point returns null
*/
CbcBranchingObject *
CbcBranchCut::preferredNewFeasible() const
{
    throw CoinError("Use of base class", "preferredNewFeasible", "CbcBranchCut");
    return new CbcCutBranchingObject();
}

/* Given valid solution (i.e. satisfied) and reduced costs etc
   returns a branching object which would give a new feasible
   point in direction opposite to one reduced cost says would be cheaper.
   If no feasible point returns null
*/
CbcBranchingObject *
CbcBranchCut::notPreferredNewFeasible() const
{
    throw CoinError("Use of base class", "notPreferredNewFeasible", "CbcBranchCut");
    return new CbcCutBranchingObject();
}

/*
  Bounds may be tightened, so it may be good to be able to refresh the local
  copy of the original bounds.
 */
void
CbcBranchCut::resetBounds()
{
}

