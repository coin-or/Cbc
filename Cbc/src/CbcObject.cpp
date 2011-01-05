// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/12/2009 carved from CbcBranchBase

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "OsiChooseVariable.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchBase.hpp"


// Default Constructor
CbcObject::CbcObject()
        : OsiObject(),
        model_(NULL),
        id_(-1),
        position_(-1),
        preferredWay_(0)
{
}

// Constructor from model
CbcObject::CbcObject(CbcModel * model)
        : OsiObject(),
        model_(model),
        id_(-1),
        position_(-1),
        preferredWay_(0)
{
}


// Destructor
CbcObject::~CbcObject ()
{
}

// Copy constructor
CbcObject::CbcObject ( const CbcObject & rhs)
        : OsiObject(rhs)
{
    model_ = rhs.model_;
    id_ = rhs.id_;
    position_ = rhs.position_;
    preferredWay_ = rhs.preferredWay_;
}

// Assignment operator
CbcObject &
CbcObject::operator=( const CbcObject & rhs)
{
    if (this != &rhs) {
        OsiObject::operator=(rhs);
        model_ = rhs.model_;
        id_ = rhs.id_;
        position_ = rhs.position_;
        preferredWay_ = rhs.preferredWay_;
    }
    return *this;
}

/* Returns floor and ceiling i.e. closest valid points
 */
void
CbcObject::floorCeiling(double & floorValue, double & ceilingValue, double value,
                        double tolerance) const
{
    if (fabs(floor(value + 0.5) - value) > tolerance) {
        floorValue = floor(value);
    } else {
        floorValue = floor(value + 0.5);
    }
    ceilingValue = floorValue + 1.0;
}
/* For the variable(s) referenced by the object,
      look at the current solution and set bounds to match the solution.
      Returns measure of how much it had to move solution to make feasible
*/
double
CbcObject::feasibleRegion(OsiSolverInterface * /*solver*/) const
{
    //assert (solver==model_->solver());
    CbcObject * fudge = const_cast<CbcObject *>(this);
    fudge->feasibleRegion();
    return 0.0;
}

/* For the variable(s) referenced by the object,
      look at the current solution and set bounds to match the solution.
      Returns measure of how much it had to move solution to make feasible
*/
double
CbcObject::feasibleRegion(OsiSolverInterface * /*solver*/,
                          const OsiBranchingInformation * /*info*/) const
{
    //assert (solver==model_->solver());
    CbcObject * fudge = const_cast<CbcObject *>(this);
    fudge->feasibleRegion();
    return 0.0;
}
/* Create a branching object and indicate which way to branch first.

      The branching object has to know how to create branches (fix
      variables, etc.)
*/
OsiBranchingObject *
CbcObject::createOsiBranch(OsiSolverInterface * solver,
                           const OsiBranchingInformation * info,
                           int way) const
{
    //assert (solver==model_->solver());
    CbcObject * fudge = const_cast<CbcObject *>(this);
    return fudge->createBranch(solver, info, way);
}
/* Create an OsiSolverBranch object

This returns NULL if branch not represented by bound changes
*/
OsiSolverBranch *
CbcObject::solverBranch() const
{
    return NULL;
}
/* Pass in information on branch just done and create CbcObjectUpdateData instance.
   If object does not need data then backward pointer will be NULL.
   Assumes can get information from solver */
CbcObjectUpdateData
CbcObject::createUpdateInformation(const OsiSolverInterface * /*solver*/,
                                   const CbcNode * /*node*/,
                                   const CbcBranchingObject * /*branchingObject*/)
{
    return CbcObjectUpdateData();
}

