// $Id$
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>

#include "CoinHelperFunctions.hpp"
//#include "CoinIndexedVector.hpp"
#include "ClpQuadraticObjective.hpp"
#include "ClpLinearObjective.hpp"
#include "ClpObjective.hpp"
//#include "ClpSimplex.hpp"
#include "ClpQuadInterface.hpp"
//#############################################################################
// Solve methods
//#############################################################################
void ClpQuadInterface::initialSolve()
{
  // save cutoff
  double cutoff = modelPtr_->dualObjectiveLimit();
  modelPtr_->setDualObjectiveLimit(1.0e50);
  modelPtr_->scaling(0);
  modelPtr_->setLogLevel(0);
  // solve with no objective to get feasible solution
  setBasis(basis_, modelPtr_);
  modelPtr_->dual();
  basis_ = getBasis(modelPtr_);
  modelPtr_->setDualObjectiveLimit(cutoff);
  if (modelPtr_->problemStatus()) {
    assert(modelPtr_->problemStatus() == 1);
    return;
  }
  ClpObjective *saveObjective = modelPtr_->objectiveAsObject();
  modelPtr_->setObjectivePointer(quadraticObjective_);
  //modelPtr_->setLogLevel(1);
  // Could load up any data into a solver
  modelPtr_->primal();
  modelPtr_->setDualObjectiveLimit(cutoff);
  if (modelPtr_->objectiveValue() > cutoff)
    modelPtr_->setProblemStatus(1);
  // zero reduced costs
  // Should not have to as convex
  //CoinZeroN(modelPtr_->dualRowSolution(),modelPtr_->numberRows());
  //CoinZeroN(modelPtr_->dualColumnSolution(),modelPtr_->numberColumns());
  modelPtr_->setObjectivePointer(saveObjective);
}

//-----------------------------------------------------------------------------
void ClpQuadInterface::resolve()
{
  initialSolve();
}

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpQuadInterface::ClpQuadInterface()
  : OsiClpSolverInterface()
{
  quadraticObjective_ = NULL;
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
OsiSolverInterface *
ClpQuadInterface::clone(bool CopyData) const
{
  if (CopyData) {
    return new ClpQuadInterface(*this);
  } else {
    printf("warning ClpQuadInterface clone with copyData false\n");
    return new ClpQuadInterface();
  }
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpQuadInterface::ClpQuadInterface(
  const ClpQuadInterface &rhs)
  : OsiClpSolverInterface(rhs)
{
  if (rhs.quadraticObjective_)
    quadraticObjective_ = rhs.quadraticObjective_->clone();
  else
    quadraticObjective_ = NULL;
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
ClpQuadInterface::~ClpQuadInterface()
{
  delete quadraticObjective_;
}

//-------------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpQuadInterface &
ClpQuadInterface::operator=(const ClpQuadInterface &rhs)
{
  if (this != &rhs) {
    OsiClpSolverInterface::operator=(rhs);
    if (rhs.quadraticObjective_)
      quadraticObjective_ = rhs.quadraticObjective_->clone();
    else
      quadraticObjective_ = NULL;
  }
  return *this;
}
//-------------------------------------------------------------------
// Real initializer
//-------------------------------------------------------------------
void ClpQuadInterface::initialize()
{
  // Save true objective and create a fake one
  delete quadraticObjective_;
  quadraticObjective_ = modelPtr_->objectiveAsObject();
  ClpLinearObjective *linearObjective = new ClpLinearObjective(NULL, modelPtr_->numberColumns());
  modelPtr_->setObjectivePointer(linearObjective);
}
// Get objective function value (can't use default)
double
ClpQuadInterface::getObjValue() const
{
  // first try easy way
  return modelPtr_->objectiveValue();
}
