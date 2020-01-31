/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
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
  : OsiObject()
  , model_(NULL)
  , id_(-1)
  , position_(-1)
  , preferredWay_(0)
{
}

// Constructor from model
CbcObject::CbcObject(CbcModel *model)
  : OsiObject()
  , model_(model)
  , id_(-1)
  , position_(-1)
  , preferredWay_(0)
{
}

// Destructor
CbcObject::~CbcObject()
{
}

// Copy constructor
CbcObject::CbcObject(const CbcObject &rhs)
  : OsiObject(rhs)
{
  model_ = rhs.model_;
  id_ = rhs.id_;
  position_ = rhs.position_;
  preferredWay_ = rhs.preferredWay_;
}

// Assignment operator
CbcObject &
CbcObject::operator=(const CbcObject &rhs)
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
void CbcObject::floorCeiling(double &floorValue, double &ceilingValue, double value,
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
  CbcObject *fudge = const_cast< CbcObject * >(this);
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
  CbcObject *fudge = const_cast< CbcObject * >(this);
  fudge->feasibleRegion();
  return 0.0;
}
/* Create a branching object and indicate which way to branch first.

      The branching object has to know how to create branches (fix
      variables, etc.)
*/
OsiBranchingObject *
CbcObject::createOsiBranch(OsiSolverInterface *solver,
  const OsiBranchingInformation *info,
  int way) const
{
  //assert (solver==model_->solver());
  CbcObject *fudge = const_cast< CbcObject * >(this);
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

// Default Constructor
CbcBranchingObject::CbcBranchingObject()
  : OsiBranchingObject()
{
  model_ = NULL;
  originalCbcObject_ = NULL;
  variable_ = -1;
  way_ = 0;
}

// Useful constructor
CbcBranchingObject::CbcBranchingObject(CbcModel *model, int variable, int way, double value)
  : OsiBranchingObject(model->solver(), value)
{
  model_ = model;
  originalCbcObject_ = NULL;
  variable_ = variable;
  way_ = way;
}

// Copy constructor
CbcBranchingObject::CbcBranchingObject(const CbcBranchingObject &rhs)
  : OsiBranchingObject(rhs)
{
  model_ = rhs.model_;
  originalCbcObject_ = rhs.originalCbcObject_;
  variable_ = rhs.variable_;
  way_ = rhs.way_;
  value_ = rhs.value_;
}

// Assignment operator
CbcBranchingObject &
CbcBranchingObject::operator=(const CbcBranchingObject &rhs)
{
  if (this != &rhs) {
    OsiBranchingObject::operator=(rhs);
    model_ = rhs.model_;
    originalCbcObject_ = rhs.originalCbcObject_;
    variable_ = rhs.variable_;
    way_ = rhs.way_;
  }
  return *this;
}

// Destructor
CbcBranchingObject::~CbcBranchingObject()
{
}
// Default Constructor
CbcBranchDecision::CbcBranchDecision()
  : object_(NULL)
  , model_(NULL)
  , chooseMethod_(NULL)
{
}

// Copy Constructor
CbcBranchDecision::CbcBranchDecision(const CbcBranchDecision &rhs)
  : object_(NULL)
  , model_(rhs.model_)
  , chooseMethod_(NULL)
{
  if (rhs.chooseMethod_)
    chooseMethod_ = rhs.chooseMethod_->clone();
}

CbcBranchDecision::~CbcBranchDecision()
{
  delete object_;
  delete chooseMethod_;
}
/* Compare N branching objects. Return index of best
   and sets way of branching in chosen object.

   This routine is used only after strong branching.
   This is reccommended version as it can be more sophisticated
*/

int CbcBranchDecision::bestBranch(CbcBranchingObject **objects, int numberObjects,
  int /*numberUnsatisfied*/,
  double *changeUp, int *numberInfeasibilitiesUp,
  double *changeDown, int *numberInfeasibilitiesDown,
  double /*objectiveValue*/)
{
  int bestWay = 0;
  int whichObject = -1;
  if (numberObjects) {
    initialize(objects[0]->model());
    CbcBranchingObject *bestObject = NULL;
    for (int i = 0; i < numberObjects; i++) {
      int betterWay = betterBranch(objects[i],
        bestObject,
        changeUp[i],
        numberInfeasibilitiesUp[i],
        changeDown[i],
        numberInfeasibilitiesDown[i]);
      if (betterWay) {
        bestObject = objects[i];
        bestWay = betterWay;
        whichObject = i;
      }
    }
    // set way in best
    if (whichObject >= 0)
      objects[whichObject]->way(bestWay);
  }
  return whichObject;
}
// Set (clone) chooseMethod
void CbcBranchDecision::setChooseMethod(const OsiChooseVariable &method)
{
  delete chooseMethod_;
  chooseMethod_ = method.clone();
}
// Default constructor
CbcConsequence::CbcConsequence()
{
}

// Destructor
CbcConsequence::~CbcConsequence()
{
}

// Copy constructor
CbcConsequence::CbcConsequence(const CbcConsequence & /*rhs*/)
{
}

// Assignment operator
CbcConsequence &
CbcConsequence::operator=(const CbcConsequence &rhs)
{
  if (this != &rhs) {
  }
  return *this;
}
// Default constructor
CbcObjectUpdateData::CbcObjectUpdateData()
  : object_(NULL)
  , way_(0)
  , objectNumber_(-1)
  , change_(0.0)
  , status_(0)
  , intDecrease_(0)
  , branchingValue_(0.0)
  , originalObjective_(COIN_DBL_MAX)
  , cutoff_(COIN_DBL_MAX)
{
}

// Useful constructor
CbcObjectUpdateData::CbcObjectUpdateData(CbcObject *object,
  int way,
  double change,
  int status,
  int intDecrease,
  double branchingValue)
  : object_(object)
  , way_(way)
  , objectNumber_(-1)
  , change_(change)
  , status_(status)
  , intDecrease_(intDecrease)
  , branchingValue_(branchingValue)
  , originalObjective_(COIN_DBL_MAX)
  , cutoff_(COIN_DBL_MAX)
{
}

// Destructor
CbcObjectUpdateData::~CbcObjectUpdateData()
{
}

// Copy constructor
CbcObjectUpdateData::CbcObjectUpdateData(const CbcObjectUpdateData &rhs)
  : object_(rhs.object_)
  , way_(rhs.way_)
  , objectNumber_(rhs.objectNumber_)
  , change_(rhs.change_)
  , status_(rhs.status_)
  , intDecrease_(rhs.intDecrease_)
  , branchingValue_(rhs.branchingValue_)
  , originalObjective_(rhs.originalObjective_)
  , cutoff_(rhs.cutoff_)
{
}

// Assignment operator
CbcObjectUpdateData &
CbcObjectUpdateData::operator=(const CbcObjectUpdateData &rhs)
{
  if (this != &rhs) {
    object_ = rhs.object_;
    way_ = rhs.way_;
    objectNumber_ = rhs.objectNumber_;
    change_ = rhs.change_;
    status_ = rhs.status_;
    intDecrease_ = rhs.intDecrease_;
    branchingValue_ = rhs.branchingValue_;
    originalObjective_ = rhs.originalObjective_;
    cutoff_ = rhs.cutoff_;
  }
  return *this;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
