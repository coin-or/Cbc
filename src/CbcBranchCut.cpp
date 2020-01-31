/* $Id$ */
// Copyright (C) 2004, International Business Machines
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
//#define CBC_DEBUG

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchCut.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

/** Default Constructor

*/
CbcBranchCut::CbcBranchCut()
  : CbcObject()
{
}

/* Constructor so model can be passed up
*/
CbcBranchCut::CbcBranchCut(CbcModel *model)
  : CbcObject(model)
{
}
// Copy constructor
CbcBranchCut::CbcBranchCut(const CbcBranchCut &rhs)
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
CbcBranchCut::operator=(const CbcBranchCut & /*rhs*/)
{
  return *this;
}

// Destructor
CbcBranchCut::~CbcBranchCut()
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
void CbcBranchCut::feasibleRegion()
{
}
/* Return true if branch created by object should fix variables
 */
bool CbcBranchCut::boundBranch() const
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
void CbcBranchCut::resetBounds()
{
}

// Default Constructor
CbcCutBranchingObject::CbcCutBranchingObject()
  : CbcBranchingObject()
{
  down_ = OsiRowCut();
  up_ = OsiRowCut();
  canFix_ = false;
}

// Useful constructor
CbcCutBranchingObject::CbcCutBranchingObject(CbcModel *model,
  OsiRowCut &down,
  OsiRowCut &up,
  bool canFix)
  : CbcBranchingObject(model, 0, -1, 0.0)
{
  down_ = down;
  up_ = up;
  canFix_ = canFix;
}

// Copy constructor
CbcCutBranchingObject::CbcCutBranchingObject(const CbcCutBranchingObject &rhs)
  : CbcBranchingObject(rhs)
{
  down_ = rhs.down_;
  up_ = rhs.up_;
  canFix_ = rhs.canFix_;
}

// Assignment operator
CbcCutBranchingObject &
CbcCutBranchingObject::operator=(const CbcCutBranchingObject &rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
    down_ = rhs.down_;
    up_ = rhs.up_;
    canFix_ = rhs.canFix_;
  }
  return *this;
}
CbcBranchingObject *
CbcCutBranchingObject::clone() const
{
  return (new CbcCutBranchingObject(*this));
}

// Destructor
CbcCutBranchingObject::~CbcCutBranchingObject()
{
}

/*
  Perform a branch by adjusting bounds and/or adding a cut. Note
  that each arm of the branch advances the object to the next arm by
  advancing the value of way_.

  Returns change in guessed objective on next branch
*/
double
CbcCutBranchingObject::branch()
{
  decrementNumberBranchesLeft();
  OsiRowCut *cut;
  if (way_ < 0) {
    cut = &down_;
    way_ = 1;
  } else {
    cut = &up_;
    way_ = -1; // Swap direction
  }
  // printf("CUT %s ", (way_ == -1) ? "up" : "down");
  // cut->print();
  // See if cut just fixes variables
  double lb = cut->lb();
  double ub = cut->ub();
  int n = cut->row().getNumElements();
  const int *column = cut->row().getIndices();
  const double *element = cut->row().getElements();
  OsiSolverInterface *solver = model_->solver();
  const double *upper = solver->getColUpper();
  const double *lower = solver->getColLower();
  double low = 0.0;
  double high = 0.0;
  for (int i = 0; i < n; i++) {
    int iColumn = column[i];
    double value = element[i];
    if (value > 0.0) {
      high += upper[iColumn] * value;
      low += lower[iColumn] * value;
    } else {
      high += lower[iColumn] * value;
      low += upper[iColumn] * value;
    }
  }
  // leave as cut
  //model_->setNextRowCut(*cut);
  //return 0.0;
  // assume cut was cunningly constructed so we need not worry too much about tolerances
  if (low + 1.0e-8 >= ub && canFix_) {
    // fix
    for (int i = 0; i < n; i++) {
      int iColumn = column[i];
      double value = element[i];
      if (value > 0.0) {
        solver->setColUpper(iColumn, lower[iColumn]);
      } else {
        solver->setColLower(iColumn, upper[iColumn]);
      }
    }
  } else if (high - 1.0e-8 <= lb && canFix_) {
    // fix
    for (int i = 0; i < n; i++) {
      int iColumn = column[i];
      double value = element[i];
      if (value > 0.0) {
        solver->setColLower(iColumn, upper[iColumn]);
      } else {
        solver->setColUpper(iColumn, lower[iColumn]);
      }
    }
  } else {
    // leave as cut
    model_->setNextRowCut(*cut);
  }
  return 0.0;
}
// Print what would happen
void CbcCutBranchingObject::print()
{
  OsiRowCut *cut;
  if (way_ < 0) {
    cut = &down_;
    printf("CbcCut would branch down");
  } else {
    cut = &up_;
    printf("CbcCut would branch up");
  }
  double lb = cut->lb();
  double ub = cut->ub();
  int n = cut->row().getNumElements();
  const int *column = cut->row().getIndices();
  const double *element = cut->row().getElements();
  if (n > 5) {
    printf(" - %d elements, lo=%g, up=%g\n", n, lb, ub);
  } else {
    printf(" - %g <=", lb);
    for (int i = 0; i < n; i++) {
      int iColumn = column[i];
      double value = element[i];
      printf(" (%d,%g)", iColumn, value);
    }
    printf(" <= %g\n", ub);
  }
}

// Return true if branch should fix variables
bool CbcCutBranchingObject::boundBranch() const
{
  return false;
}

/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type.
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int CbcCutBranchingObject::compareOriginalObject(const CbcBranchingObject *brObj) const
{
  const CbcCutBranchingObject *br = dynamic_cast< const CbcCutBranchingObject * >(brObj);
  assert(br);
  const OsiRowCut &r0 = way_ == -1 ? down_ : up_;
  const OsiRowCut &r1 = br->way_ == -1 ? br->down_ : br->up_;
  return r0.row().compare(r1.row());
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
CbcCutBranchingObject::compareBranchingObject(const CbcBranchingObject *brObj, const bool replaceIfOverlap)
{
  const CbcCutBranchingObject *br = dynamic_cast< const CbcCutBranchingObject * >(brObj);
  assert(br);
  OsiRowCut &r0 = way_ == -1 ? down_ : up_;
  const OsiRowCut &r1 = br->way_ == -1 ? br->down_ : br->up_;
  double thisBd[2];
  thisBd[0] = r0.lb();
  thisBd[1] = r0.ub();
  double otherBd[2];
  otherBd[0] = r1.lb();
  otherBd[1] = r1.ub();
  CbcRangeCompare comp = CbcCompareRanges(thisBd, otherBd, replaceIfOverlap);
  if (comp != CbcRangeOverlap || (comp == CbcRangeOverlap && !replaceIfOverlap)) {
    return comp;
  }
  r0.setLb(thisBd[0]);
  r0.setUb(thisBd[1]);
  return comp;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
