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
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchLotsize.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

/*
  CBC_PRINT 1 just does sanity checks - no printing
  Larger values of CBC_PRINT set various printing levels.  Larger
  values print more.
*/
//#define CBC_PRINT 1
// First/last variable to print info on
#if CBC_PRINT
// preset does all - change to x,x to just do x
static int firstPrint = 0;
static int lastPrint = 1000000;
static CbcModel *saveModel = NULL;
#endif
// Just for debug (CBC_PRINT defined in CbcBranchLotsize.cpp)
void
#if CBC_PRINT
CbcLotsize::printLotsize(double value, bool condition, int type) const
#else
CbcLotsize::printLotsize(double, bool, int) const
#endif
{
#if CBC_PRINT
  if (columnNumber_ >= firstPrint && columnNumber_ <= lastPrint) {
    int printIt = CBC_PRINT - 1;
    // Get details
    OsiSolverInterface *solver = saveModel->solver();
    double currentLower = solver->getColLower()[columnNumber_];
    double currentUpper = solver->getColUpper()[columnNumber_];
    int i;
    // See if in a valid range (with two tolerances)
    bool inRange = false;
    bool inRange2 = false;
    double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
    // increase if type 2
    if (type == 2) {
      integerTolerance *= 100.0;
      type = 0;
      printIt = 2; // always print
    }
    // bounds should match some bound
    int rangeL = -1;
    int rangeU = -1;
    if (rangeType_ == 1) {
      for (i = 0; i < numberRanges_; i++) {
        if (fabs(currentLower - bound_[i]) < 1.0e-12)
          rangeL = i;
        if (fabs(currentUpper - bound_[i]) < 1.0e-12)
          rangeU = i;
        if (fabs(value - bound_[i]) < integerTolerance)
          inRange = true;
        if (fabs(value - bound_[i]) < 1.0e8)
          inRange2 = true;
      }
    } else {
      for (i = 0; i < numberRanges_; i++) {
        if (fabs(currentLower - bound_[2 * i]) < 1.0e-12)
          rangeL = i;
        if (fabs(currentUpper - bound_[2 * i + 1]) < 1.0e-12)
          rangeU = i;
        if (value > bound_[2 * i] - integerTolerance && value < bound_[2 * i + 1] + integerTolerance)
          inRange = true;
        if (value > bound_[2 * i] - integerTolerance && value < bound_[2 * i + 1] + integerTolerance)
          inRange = true;
      }
    }
    assert(rangeL >= 0 && rangeU >= 0);
    bool abortIt = false;
    switch (type) {
      // returning from findRange (fall through to just check)
    case 0:
      if (printIt) {
        printf("findRange returns %s for column %d and value %g",
          condition ? "true" : "false", columnNumber_, value);
        if (printIt > 1)
          printf(" LP bounds %g, %g", currentLower, currentUpper);
        printf("\n");
      }
      // Should match
    case 1:
      if (inRange != condition) {
        printIt = 2;
        abortIt = true;
      }
      break;
      //
    case 2:
      break;
      //
    case 3:
      break;
      //
    case 4:
      break;
    }
  }
#endif
}
/** Default Constructor

*/
CbcLotsize::CbcLotsize()
  : CbcObject()
  , columnNumber_(-1)
  , rangeType_(0)
  , numberRanges_(0)
  , largestGap_(0)
  , bound_(NULL)
  , range_(0)
{
}

/** Useful constructor

  Loads actual upper & lower bounds for the specified variable.
*/
CbcLotsize::CbcLotsize(CbcModel *model,
  int iColumn, int numberPoints,
  const double *points, bool range)
  : CbcObject(model)
{
#if CBC_PRINT
  if (!saveModel)
    saveModel = model;
#endif
  assert(numberPoints > 0);
  columnNumber_ = iColumn;
  // and set id so can be used for branching
  id_ = iColumn;
  // sort ranges
  int *sort = new int[numberPoints];
  double *weight = new double[numberPoints];
  int i;
  if (range) {
    rangeType_ = 2;
  } else {
    rangeType_ = 1;
  }
  for (i = 0; i < numberPoints; i++) {
    sort[i] = i;
    weight[i] = points[i * rangeType_];
  }
  CoinSort_2(weight, weight + numberPoints, sort);
  numberRanges_ = 1;
  largestGap_ = 0;
  if (rangeType_ == 1) {
    bound_ = new double[numberPoints + 1];
    bound_[0] = weight[0];
    for (i = 1; i < numberPoints; i++) {
      if (weight[i] != weight[i - 1])
        bound_[numberRanges_++] = weight[i];
    }
    // and for safety
    bound_[numberRanges_] = bound_[numberRanges_ - 1];
    for (i = 1; i < numberRanges_; i++) {
      largestGap_ = CoinMax(largestGap_, bound_[i] - bound_[i - 1]);
    }
  } else {
    bound_ = new double[2 * numberPoints + 2];
    bound_[0] = points[sort[0] * 2];
    bound_[1] = points[sort[0] * 2 + 1];
    double hi = bound_[1];
    assert(hi >= bound_[0]);
    for (i = 1; i < numberPoints; i++) {
      double thisLo = points[sort[i] * 2];
      double thisHi = points[sort[i] * 2 + 1];
      assert(thisHi >= thisLo);
      if (thisLo > hi) {
        bound_[2 * numberRanges_] = thisLo;
        bound_[2 * numberRanges_ + 1] = thisHi;
        numberRanges_++;
        hi = thisHi;
      } else {
        //overlap
        hi = CoinMax(hi, thisHi);
        bound_[2 * numberRanges_ - 1] = hi;
      }
    }
    // and for safety
    bound_[2 * numberRanges_] = bound_[2 * numberRanges_ - 2];
    bound_[2 * numberRanges_ + 1] = bound_[2 * numberRanges_ - 1];
    for (i = 1; i < numberRanges_; i++) {
      largestGap_ = CoinMax(largestGap_, bound_[2 * i] - bound_[2 * i - 1]);
    }
  }
  delete[] sort;
  delete[] weight;
  range_ = 0;
}

// Copy constructor
CbcLotsize::CbcLotsize(const CbcLotsize &rhs)
  : CbcObject(rhs)

{
  columnNumber_ = rhs.columnNumber_;
  rangeType_ = rhs.rangeType_;
  numberRanges_ = rhs.numberRanges_;
  range_ = rhs.range_;
  largestGap_ = rhs.largestGap_;
  if (numberRanges_) {
    assert(rangeType_ > 0 && rangeType_ < 3);
    bound_ = new double[(numberRanges_ + 1) * rangeType_];
    memcpy(bound_, rhs.bound_, (numberRanges_ + 1) * rangeType_ * sizeof(double));
  } else {
    bound_ = NULL;
  }
}

// Clone
CbcObject *
CbcLotsize::clone() const
{
  return new CbcLotsize(*this);
}

// Assignment operator
CbcLotsize &
CbcLotsize::operator=(const CbcLotsize &rhs)
{
  if (this != &rhs) {
    CbcObject::operator=(rhs);
    columnNumber_ = rhs.columnNumber_;
    rangeType_ = rhs.rangeType_;
    numberRanges_ = rhs.numberRanges_;
    largestGap_ = rhs.largestGap_;
    delete[] bound_;
    range_ = rhs.range_;
    if (numberRanges_) {
      assert(rangeType_ > 0 && rangeType_ < 3);
      bound_ = new double[(numberRanges_ + 1) * rangeType_];
      memcpy(bound_, rhs.bound_, (numberRanges_ + 1) * rangeType_ * sizeof(double));
    } else {
      bound_ = NULL;
    }
  }
  return *this;
}

// Destructor
CbcLotsize::~CbcLotsize()
{
  delete[] bound_;
}
/* Finds range of interest so value is feasible in range range_ or infeasible
   between hi[range_] and lo[range_+1].  Returns true if feasible.
*/
bool CbcLotsize::findRange(double value) const
{
  assert(range_ >= 0 && range_ < numberRanges_ + 1);
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  int iLo;
  int iHi;
  double infeasibility = 0.0;
  if (rangeType_ == 1) {
    if (value < bound_[range_] - integerTolerance) {
      iLo = 0;
      iHi = range_ - 1;
    } else if (value < bound_[range_] + integerTolerance) {
#if CBC_PRINT
      printLotsize(value, true, 0);
#endif
      return true;
    } else if (value < bound_[range_ + 1] - integerTolerance) {
#ifdef CBC_PRINT
      printLotsize(value, false, 0);
#endif
      return false;
    } else {
      iLo = range_ + 1;
      iHi = numberRanges_ - 1;
    }
    // check lo and hi
    bool found = false;
    if (value > bound_[iLo] - integerTolerance && value < bound_[iLo + 1] + integerTolerance) {
      range_ = iLo;
      found = true;
    } else if (value > bound_[iHi] - integerTolerance && value < bound_[iHi + 1] + integerTolerance) {
      range_ = iHi;
      found = true;
    } else {
      range_ = (iLo + iHi) >> 1;
    }
    //points
    while (!found) {
      if (value < bound_[range_]) {
        if (value >= bound_[range_ - 1]) {
          // found
          range_--;
          break;
        } else {
          iHi = range_;
        }
      } else {
        if (value < bound_[range_ + 1]) {
          // found
          break;
        } else {
          iLo = range_;
        }
      }
      range_ = (iLo + iHi) >> 1;
    }
    if (value - bound_[range_] <= bound_[range_ + 1] - value) {
      infeasibility = value - bound_[range_];
    } else {
      infeasibility = bound_[range_ + 1] - value;
      if (infeasibility < integerTolerance)
        range_++;
    }
#ifdef CBC_PRINT
    printLotsize(value, (infeasibility < integerTolerance), 0);
#endif
    return (infeasibility < integerTolerance);
  } else {
    // ranges
    if (value < bound_[2 * range_] - integerTolerance) {
      iLo = 0;
      iHi = range_ - 1;
    } else if (value < bound_[2 * range_ + 1] + integerTolerance) {
#ifdef CBC_PRINT
      printLotsize(value, true, 0);
#endif
      return true;
    } else if (value < bound_[2 * range_ + 2] - integerTolerance) {
#ifdef CBC_PRINT
      printLotsize(value, false, 0);
#endif
      return false;
    } else {
      iLo = range_ + 1;
      iHi = numberRanges_ - 1;
    }
    // check lo and hi
    bool found = false;
    if (value > bound_[2 * iLo] - integerTolerance && value < bound_[2 * iLo + 2] - integerTolerance) {
      range_ = iLo;
      found = true;
    } else if (value >= bound_[2 * iHi] - integerTolerance) {
      range_ = iHi;
      found = true;
    } else {
      range_ = (iLo + iHi) >> 1;
    }
    //points
    while (!found) {
      if (value < bound_[2 * range_]) {
        if (value >= bound_[2 * range_ - 2]) {
          // found
          range_--;
          break;
        } else {
          iHi = range_;
        }
      } else {
        if (value < bound_[2 * range_ + 2]) {
          // found
          break;
        } else {
          iLo = range_;
        }
      }
      range_ = (iLo + iHi) >> 1;
    }
    if (value >= bound_[2 * range_] - integerTolerance && value <= bound_[2 * range_ + 1] + integerTolerance)
      infeasibility = 0.0;
    else if (value - bound_[2 * range_ + 1] < bound_[2 * range_ + 2] - value) {
      infeasibility = value - bound_[2 * range_ + 1];
    } else {
      infeasibility = bound_[2 * range_ + 2] - value;
    }
#ifdef CBC_PRINT
    printLotsize(value, (infeasibility < integerTolerance), 0);
#endif
    return (infeasibility < integerTolerance);
  }
}
/* Returns floor and ceiling
 */
void CbcLotsize::floorCeiling(double &floorLotsize, double &ceilingLotsize, double value,
  double /*tolerance*/) const
{
  bool feasible = findRange(value);
  if (rangeType_ == 1) {
    floorLotsize = bound_[range_];
    ceilingLotsize = bound_[range_ + 1];
    // may be able to adjust
    if (feasible && fabs(value - floorLotsize) > fabs(value - ceilingLotsize)) {
      floorLotsize = bound_[range_ + 1];
      ceilingLotsize = bound_[range_ + 2];
    }
  } else {
    // ranges
    assert(value >= bound_[2 * range_ + 1]);
    floorLotsize = bound_[2 * range_ + 1];
    ceilingLotsize = bound_[2 * range_ + 2];
  }
}
double
CbcLotsize::infeasibility(const OsiBranchingInformation * /*info*/,
  int &preferredWay) const
{
  OsiSolverInterface *solver = model_->solver();
  const double *solution = model_->testSolution();
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  double value = solution[columnNumber_];
  value = CoinMax(value, lower[columnNumber_]);
  value = CoinMin(value, upper[columnNumber_]);
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  /*printf("%d %g %g %g %g\n",columnNumber_,value,lower[columnNumber_],
      solution[columnNumber_],upper[columnNumber_]);*/
  assert(value >= bound_[0] - integerTolerance
    && value <= bound_[rangeType_ * numberRanges_ - 1] + integerTolerance);
  double infeasibility = 0.0;
  bool feasible = findRange(value);
  if (!feasible) {
    if (rangeType_ == 1) {
      if (value - bound_[range_] < bound_[range_ + 1] - value) {
        preferredWay = -1;
        infeasibility = value - bound_[range_];
      } else {
        preferredWay = 1;
        infeasibility = bound_[range_ + 1] - value;
      }
    } else {
      // ranges
      if (value - bound_[2 * range_ + 1] < bound_[2 * range_ + 2] - value) {
        preferredWay = -1;
        infeasibility = value - bound_[2 * range_ + 1];
      } else {
        preferredWay = 1;
        infeasibility = bound_[2 * range_ + 2] - value;
      }
    }
  } else {
    // always satisfied
    preferredWay = -1;
  }
  if (infeasibility < integerTolerance)
    infeasibility = 0.0;
  else
    infeasibility /= largestGap_;
#ifdef CBC_PRINT
  printLotsize(value, infeasibility, 1);
#endif
  return infeasibility;
}
/* Column number if single column object -1 otherwise,
   so returns >= 0
   Used by heuristics
*/
int CbcLotsize::columnNumber() const
{
  return columnNumber_;
}
// This looks at solution and sets bounds to contain solution
/** More precisely: it first forces the variable within the existing
    bounds, and then tightens the bounds to make sure the variable is feasible
*/
void CbcLotsize::feasibleRegion()
{
  OsiSolverInterface *solver = model_->solver();
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  const double *solution = model_->testSolution();
  double value = solution[columnNumber_];
  value = CoinMax(value, lower[columnNumber_]);
  value = CoinMin(value, upper[columnNumber_]);
  findRange(value);
  double nearest;
  if (rangeType_ == 1) {
    nearest = bound_[range_];
    solver->setColLower(columnNumber_, nearest);
    solver->setColUpper(columnNumber_, nearest);
  } else {
    // ranges
    solver->setColLower(columnNumber_, CoinMax(bound_[2 * range_], lower[columnNumber_]));
    solver->setColUpper(columnNumber_, CoinMin(bound_[2 * range_ + 1], upper[columnNumber_]));
    if (value > bound_[2 * range_ + 1])
      nearest = bound_[2 * range_ + 1];
    else if (value < bound_[2 * range_])
      nearest = bound_[2 * range_];
    else
      nearest = value;
  }
#ifdef CBC_PRINT
  // print details
  printLotsize(value, true, 2);
#endif
  // Scaling may have moved it a bit
  // Lotsizing variables could be a lot larger
#ifndef NDEBUG
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  assert(fabs(value - nearest) <= (100.0 + 10.0 * fabs(nearest)) * integerTolerance);
#endif
}
CbcBranchingObject *
CbcLotsize::createCbcBranch(OsiSolverInterface *solver, const OsiBranchingInformation * /*info*/, int way)
{
  //OsiSolverInterface * solver = model_->solver();
  const double *solution = model_->testSolution();
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  double value = solution[columnNumber_];
  value = CoinMax(value, lower[columnNumber_]);
  value = CoinMin(value, upper[columnNumber_]);
  assert(!findRange(value));
  return new CbcLotsizeBranchingObject(model_, columnNumber_, way,
    value, this);
}

/* Given valid solution (i.e. satisfied) and reduced costs etc
   returns a branching object which would give a new feasible
   point in direction reduced cost says would be cheaper.
   If no feasible point returns null
*/
CbcBranchingObject *
CbcLotsize::preferredNewFeasible() const
{
  OsiSolverInterface *solver = model_->solver();

  assert(findRange(model_->testSolution()[columnNumber_]));
  double dj = solver->getObjSense() * solver->getReducedCost()[columnNumber_];
  CbcLotsizeBranchingObject *object = NULL;
  double lo, up;
  if (dj >= 0.0) {
    // can we go down
    if (range_) {
      // yes
      if (rangeType_ == 1) {
        lo = bound_[range_ - 1];
        up = bound_[range_ - 1];
      } else {
        lo = bound_[2 * range_ - 2];
        up = bound_[2 * range_ - 1];
      }
      object = new CbcLotsizeBranchingObject(model_, columnNumber_, -1,
        lo, up);
    }
  } else {
    // can we go up
    if (range_ < numberRanges_ - 1) {
      // yes
      if (rangeType_ == 1) {
        lo = bound_[range_ + 1];
        up = bound_[range_ + 1];
      } else {
        lo = bound_[2 * range_ + 2];
        up = bound_[2 * range_ + 3];
      }
      object = new CbcLotsizeBranchingObject(model_, columnNumber_, -1,
        lo, up);
    }
  }
  return object;
}

/* Given valid solution (i.e. satisfied) and reduced costs etc
   returns a branching object which would give a new feasible
   point in direction opposite to one reduced cost says would be cheaper.
   If no feasible point returns null
*/
CbcBranchingObject *
CbcLotsize::notPreferredNewFeasible() const
{
  OsiSolverInterface *solver = model_->solver();

#ifndef NDEBUG
  double value = model_->testSolution()[columnNumber_];
  double nearest = floor(value + 0.5);
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  // Scaling may have moved it a bit
  // Lotsizing variables could be a lot larger
  assert(fabs(value - nearest) <= (10.0 + 10.0 * fabs(nearest)) * integerTolerance);
#endif
  double dj = solver->getObjSense() * solver->getReducedCost()[columnNumber_];
  CbcLotsizeBranchingObject *object = NULL;
  double lo, up;
  if (dj <= 0.0) {
    // can we go down
    if (range_) {
      // yes
      if (rangeType_ == 1) {
        lo = bound_[range_ - 1];
        up = bound_[range_ - 1];
      } else {
        lo = bound_[2 * range_ - 2];
        up = bound_[2 * range_ - 1];
      }
      object = new CbcLotsizeBranchingObject(model_, columnNumber_, -1,
        lo, up);
    }
  } else {
    // can we go up
    if (range_ < numberRanges_ - 1) {
      // yes
      if (rangeType_ == 1) {
        lo = bound_[range_ + 1];
        up = bound_[range_ + 1];
      } else {
        lo = bound_[2 * range_ + 2];
        up = bound_[2 * range_ + 3];
      }
      object = new CbcLotsizeBranchingObject(model_, columnNumber_, -1,
        lo, up);
    }
  }
  return object;
}

/*
  Bounds may be tightened, so it may be good to be able to refresh the local
  copy of the original bounds.
 */
void CbcLotsize::resetBounds(const OsiSolverInterface * /*solver*/)
{
}

// Default Constructor
CbcLotsizeBranchingObject::CbcLotsizeBranchingObject()
  : CbcBranchingObject()
{
  down_[0] = 0.0;
  down_[1] = 0.0;
  up_[0] = 0.0;
  up_[1] = 0.0;
}

// Useful constructor
CbcLotsizeBranchingObject::CbcLotsizeBranchingObject(CbcModel *model,
  int variable, int way, double value,
  const CbcLotsize *lotsize)
  : CbcBranchingObject(model, variable, way, value)
{
  int iColumn = lotsize->modelSequence();
  assert(variable == iColumn);
  down_[0] = model_->solver()->getColLower()[iColumn];
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  lotsize->floorCeiling(down_[1], up_[0], value, integerTolerance);
  up_[1] = model->getColUpper()[iColumn];
}
// Useful constructor for fixing
CbcLotsizeBranchingObject::CbcLotsizeBranchingObject(CbcModel *model,
  int variable, int way,
  double lowerValue,
  double upperValue)
  : CbcBranchingObject(model, variable, way, lowerValue)
{
  setNumberBranchesLeft(1);
  down_[0] = lowerValue;
  down_[1] = upperValue;
  up_[0] = lowerValue;
  up_[1] = upperValue;
}

// Copy constructor
CbcLotsizeBranchingObject::CbcLotsizeBranchingObject(const CbcLotsizeBranchingObject &rhs)
  : CbcBranchingObject(rhs)
{
  down_[0] = rhs.down_[0];
  down_[1] = rhs.down_[1];
  up_[0] = rhs.up_[0];
  up_[1] = rhs.up_[1];
}

// Assignment operator
CbcLotsizeBranchingObject &
CbcLotsizeBranchingObject::operator=(const CbcLotsizeBranchingObject &rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
    down_[0] = rhs.down_[0];
    down_[1] = rhs.down_[1];
    up_[0] = rhs.up_[0];
    up_[1] = rhs.up_[1];
  }
  return *this;
}
CbcBranchingObject *
CbcLotsizeBranchingObject::clone() const
{
  return (new CbcLotsizeBranchingObject(*this));
}

// Destructor
CbcLotsizeBranchingObject::~CbcLotsizeBranchingObject()
{
}

/*
  Perform a branch by adjusting the bounds of the specified variable. Note
  that each arm of the branch advances the object to the next arm by
  advancing the value of way_.

  Providing new values for the variable's lower and upper bounds for each
  branching direction gives a little bit of additional flexibility and will
  be easily extensible to multi-way branching.
*/
double
CbcLotsizeBranchingObject::branch()
{
  decrementNumberBranchesLeft();
  int iColumn = variable_;
  if (way_ < 0) {
#ifndef NDEBUG
    {
      double olb, oub;
      olb = model_->solver()->getColLower()[iColumn];
      oub = model_->solver()->getColUpper()[iColumn];
#ifdef CBC_DEBUG
      printf("branching down on var %d: [%g,%g] => [%g,%g]\n",
        iColumn, olb, oub, down_[0], down_[1]);
#endif
      assert(olb < down_[0] + 1.0e-7 && oub > down_[1] - 1.0e-7);
    }
#endif
    model_->solver()->setColLower(iColumn, down_[0]);
    model_->solver()->setColUpper(iColumn, down_[1]);
    way_ = 1;
  } else {
#ifndef NDEBUG
    {
      double olb, oub;
      olb = model_->solver()->getColLower()[iColumn];
      oub = model_->solver()->getColUpper()[iColumn];
#ifdef CBC_DEBUG
      printf("branching up on var %d: [%g,%g] => [%g,%g]\n",
        iColumn, olb, oub, up_[0], up_[1]);
#endif
      assert(olb < up_[0] + 1.0e-7 && oub > up_[1] - 1.0e-7);
    }
#endif
    model_->solver()->setColLower(iColumn, up_[0]);
    model_->solver()->setColUpper(iColumn, up_[1]);
    way_ = -1; // Swap direction
  }
  return 0.0;
}
// Print
void CbcLotsizeBranchingObject::print()
{
  int iColumn = variable_;
  if (way_ < 0) {
    {
      double olb, oub;
      olb = model_->solver()->getColLower()[iColumn];
      oub = model_->solver()->getColUpper()[iColumn];
      printf("branching down on var %d: [%g,%g] => [%g,%g]\n",
        iColumn, olb, oub, down_[0], down_[1]);
    }
  } else {
    {
      double olb, oub;
      olb = model_->solver()->getColLower()[iColumn];
      oub = model_->solver()->getColUpper()[iColumn];
      printf("branching up on var %d: [%g,%g] => [%g,%g]\n",
        iColumn, olb, oub, up_[0], up_[1]);
    }
  }
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
CbcLotsizeBranchingObject::compareBranchingObject(const CbcBranchingObject *brObj, const bool replaceIfOverlap)
{
  const CbcLotsizeBranchingObject *br = dynamic_cast< const CbcLotsizeBranchingObject * >(brObj);
  assert(br);
  double *thisBd = way_ == -1 ? down_ : up_;
  const double *otherBd = br->way_ == -1 ? br->down_ : br->up_;
  return CbcCompareRanges(thisBd, otherBd, replaceIfOverlap);
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
