// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/10/2009-- carved out of CbcBranchActual

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "CoinTypes.h"
#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcFollowOn.hpp"
#include "CbcBranchActual.hpp"
#include "CbcBranchCut.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

//##############################################################################

// Default Constructor
CbcFollowOn::CbcFollowOn()
  : CbcObject()
  , rhs_(NULL)
{
}

// Useful constructor
CbcFollowOn::CbcFollowOn(CbcModel *model)
  : CbcObject(model)
{
  assert(model);
  OsiSolverInterface *solver = model_->solver();
  matrix_ = *solver->getMatrixByCol();
  matrix_.removeGaps();
  matrix_.setExtraGap(0.0);
  matrixByRow_ = *solver->getMatrixByRow();
  int numberRows = matrix_.getNumRows();

  rhs_ = new int[numberRows];
  int i;
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  // Row copy
  const double *elementByRow = matrixByRow_.getElements();
  const int *column = matrixByRow_.getIndices();
  const CoinBigIndex *rowStart = matrixByRow_.getVectorStarts();
  const int *rowLength = matrixByRow_.getVectorLengths();
  for (i = 0; i < numberRows; i++) {
    rhs_[i] = 0;
    double value = rowLower[i];
    if (value == rowUpper[i]) {
      if (floor(value) == value && value >= 1.0 && value < 10.0) {
        // check elements
        bool good = true;
        for (CoinBigIndex j = rowStart[i]; j < rowStart[i] + rowLength[i]; j++) {
          int iColumn = column[j];
          if (!solver->isBinary(iColumn))
            good = false;
          double elValue = elementByRow[j];
          if (floor(elValue) != elValue || value < 1.0)
            good = false;
        }
        if (good)
          rhs_[i] = static_cast< int >(value);
      }
    }
  }
}

// Copy constructor
CbcFollowOn::CbcFollowOn(const CbcFollowOn &rhs)
  : CbcObject(rhs)
  , matrix_(rhs.matrix_)
  , matrixByRow_(rhs.matrixByRow_)
{
  int numberRows = matrix_.getNumRows();
  rhs_ = CoinCopyOfArray(rhs.rhs_, numberRows);
}

// Clone
CbcObject *
CbcFollowOn::clone() const
{
  return new CbcFollowOn(*this);
}

// Assignment operator
CbcFollowOn &
CbcFollowOn::operator=(const CbcFollowOn &rhs)
{
  if (this != &rhs) {
    CbcObject::operator=(rhs);
    delete[] rhs_;
    matrix_ = rhs.matrix_;
    matrixByRow_ = rhs.matrixByRow_;
    int numberRows = matrix_.getNumRows();
    rhs_ = CoinCopyOfArray(rhs.rhs_, numberRows);
  }
  return *this;
}

// Destructor
CbcFollowOn::~CbcFollowOn()
{
  delete[] rhs_;
}
// As some computation is needed in more than one place - returns row
int CbcFollowOn::gutsOfFollowOn(int &otherRow, int &preferredWay) const
{
  int whichRow = -1;
  otherRow = -1;
  int numberRows = matrix_.getNumRows();

  int i;
  // For sorting
  int *sort = new int[numberRows];
  int *isort = new int[numberRows];
  // Column copy
  //const double * element = matrix_.getElements();
  const int *row = matrix_.getIndices();
  const CoinBigIndex *columnStart = matrix_.getVectorStarts();
  const int *columnLength = matrix_.getVectorLengths();
  // Row copy
  const double *elementByRow = matrixByRow_.getElements();
  const int *column = matrixByRow_.getIndices();
  const CoinBigIndex *rowStart = matrixByRow_.getVectorStarts();
  const int *rowLength = matrixByRow_.getVectorLengths();
  OsiSolverInterface *solver = model_->solver();
  const double *columnLower = solver->getColLower();
  const double *columnUpper = solver->getColUpper();
  const double *solution = solver->getColSolution();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  int nSort = 0;
  for (i = 0; i < numberRows; i++) {
    if (rhs_[i]) {
      // check elements
      double smallest = 1.0e10;
      double largest = 0.0;
      int rhsValue = rhs_[i];
      int number1 = 0;
      int numberUnsatisfied = 0;
      for (CoinBigIndex j = rowStart[i]; j < rowStart[i] + rowLength[i]; j++) {
        int iColumn = column[j];
        double value = elementByRow[j];
        double solValue = solution[iColumn];
        if (columnLower[iColumn] != columnUpper[iColumn]) {
          smallest = std::min(smallest, value);
          largest = std::max(largest, value);
          if (value == 1.0)
            number1++;
          if (solValue < 1.0 - integerTolerance && solValue > integerTolerance)
            numberUnsatisfied++;
        } else {
          rhsValue -= static_cast< int >(value * floor(solValue + 0.5));
        }
      }
      if (numberUnsatisfied > 1) {
        if (smallest < largest) {
          // probably no good but check a few things
          assert(largest <= rhsValue);
          if (number1 == 1 && largest == rhsValue)
            printf("could fix\n");
        } else if (largest == rhsValue) {
          sort[nSort] = i;
          isort[nSort++] = -numberUnsatisfied;
        }
      }
    }
  }
  if (nSort > 1) {
    CoinSort_2(isort, isort + nSort, sort);
    CoinZeroN(isort, numberRows);
    double *other = new double[numberRows];
    CoinZeroN(other, numberRows);
    int *which = new int[numberRows];
    //#define COUNT
#ifndef COUNT
    bool beforeSolution = model_->getSolutionCount() == 0;
#endif
    for (int k = 0; k < nSort - 1; k++) {
      i = sort[k];
#if defined(COUNT) || !defined(NDEBUG)
      int numberUnsatisfied = 0;
#endif
      int n = 0;
      CoinBigIndex j;
      for (j = rowStart[i]; j < rowStart[i] + rowLength[i]; j++) {
        int iColumn = column[j];
        if (columnLower[iColumn] != columnUpper[iColumn]) {
          double solValue = solution[iColumn] - columnLower[iColumn];
          if (solValue < 1.0 - integerTolerance && solValue > integerTolerance) {
#if defined(COUNT) || !defined(NDEBUG)
            numberUnsatisfied++;
#endif
            for (CoinBigIndex jj = columnStart[iColumn]; jj < columnStart[iColumn] + columnLength[iColumn]; jj++) {
              int iRow = row[jj];
              if (rhs_[iRow]) {
                other[iRow] += solValue;
                if (isort[iRow]) {
                  isort[iRow]++;
                } else {
                  isort[iRow] = 1;
                  which[n++] = iRow;
                }
              }
            }
          }
        }
      }
      double total = 0.0;
      // Take out row
      double sumThis = other[i];
      other[i] = 0.0;
      assert(numberUnsatisfied == isort[i]);
      // find one nearest half if solution, one if before solution
      int iBest = -1;
      double dtarget = 0.5 * total;
#ifdef COUNT
      int target = (numberUnsatisfied + 1) >> 1;
      int best = numberUnsatisfied;
#else
      double best;
      if (beforeSolution)
        best = dtarget;
      else
        best = 1.0e30;
#endif
      for (j = 0; j < n; j++) {
        int iRow = which[j];
        double dvalue = other[iRow];
        other[iRow] = 0.0;
#ifdef COUNT
        int value = isort[iRow];
#endif
        isort[iRow] = 0;
        if (fabs(dvalue) < 1.0e-8 || fabs(sumThis - dvalue) < 1.0e-8)
          continue;
        if (dvalue < integerTolerance || dvalue > 1.0 - integerTolerance)
          continue;
#ifdef COUNT
        if (abs(value - target) < best && value != numberUnsatisfied) {
          best = abs(value - target);
          iBest = iRow;
          if (dvalue < dtarget)
            preferredWay = 1;
          else
            preferredWay = -1;
        }
#else
        if (beforeSolution) {
          if (fabs(dvalue - dtarget) > best) {
            best = fabs(dvalue - dtarget);
            iBest = iRow;
            if (dvalue < dtarget)
              preferredWay = 1;
            else
              preferredWay = -1;
          }
        } else {
          if (fabs(dvalue - dtarget) < best) {
            best = fabs(dvalue - dtarget);
            iBest = iRow;
            if (dvalue < dtarget)
              preferredWay = 1;
            else
              preferredWay = -1;
          }
        }
#endif
      }
      if (iBest >= 0) {
        whichRow = i;
        otherRow = iBest;
        break;
      }
    }
    delete[] which;
    delete[] other;
  }
  delete[] sort;
  delete[] isort;
  return whichRow;
}
double
CbcFollowOn::infeasibility(const OsiBranchingInformation * /*info*/,
  int &preferredWay) const
{
  int otherRow = 0;
  int whichRow = gutsOfFollowOn(otherRow, preferredWay);
  if (whichRow < 0)
    return 0.0;
  else
    return 2.0 * model_->getDblParam(CbcModel::CbcIntegerTolerance);
}

// This looks at solution and sets bounds to contain solution
void CbcFollowOn::feasibleRegion()
{
}

CbcBranchingObject *
CbcFollowOn::createCbcBranch(OsiSolverInterface *solver, const OsiBranchingInformation * /*info*/, int way)
{
  int otherRow = 0;
  int preferredWay;
  int whichRow = gutsOfFollowOn(otherRow, preferredWay);
  assert(way == preferredWay);
  assert(whichRow >= 0);
  int numberColumns = matrix_.getNumCols();

  // Column copy
  //const double * element = matrix_.getElements();
  const int *row = matrix_.getIndices();
  const CoinBigIndex *columnStart = matrix_.getVectorStarts();
  const int *columnLength = matrix_.getVectorLengths();
  // Row copy
  //const double * elementByRow = matrixByRow_.getElements();
  const int *column = matrixByRow_.getIndices();
  const CoinBigIndex *rowStart = matrixByRow_.getVectorStarts();
  const int *rowLength = matrixByRow_.getVectorLengths();
  //OsiSolverInterface * solver = model_->solver();
  const double *columnLower = solver->getColLower();
  const double *columnUpper = solver->getColUpper();
  //const double * solution = solver->getColSolution();
  int nUp = 0;
  int nDown = 0;
  int *upList = new int[numberColumns];
  int *downList = new int[numberColumns];
  CoinBigIndex j;
  for (j = rowStart[whichRow]; j < rowStart[whichRow] + rowLength[whichRow]; j++) {
    int iColumn = column[j];
    if (columnLower[iColumn] != columnUpper[iColumn]) {
      bool up = true;
      for (CoinBigIndex jj = columnStart[iColumn]; jj < columnStart[iColumn] + columnLength[iColumn]; jj++) {
        int iRow = row[jj];
        if (iRow == otherRow) {
          up = false;
          break;
        }
      }
      if (up)
        upList[nUp++] = iColumn;
      else
        downList[nDown++] = iColumn;
    }
  }
  //printf("way %d\n",way);
  // create object
  //printf("would fix %d down and %d up\n",nDown,nUp);
  CbcBranchingObject *branch
    = new CbcFixingBranchingObject(model_, way,
      nDown, downList, nUp, upList);
  delete[] upList;
  delete[] downList;
  return branch;
}

//##############################################################################

// Default Constructor
CbcFixingBranchingObject::CbcFixingBranchingObject()
  : CbcBranchingObject()
{
  numberDown_ = 0;
  numberUp_ = 0;
  downList_ = NULL;
  upList_ = NULL;
}

// Useful constructor
CbcFixingBranchingObject::CbcFixingBranchingObject(CbcModel *model,
  int way,
  int numberOnDownSide, const int *down,
  int numberOnUpSide, const int *up)
  : CbcBranchingObject(model, 0, way, 0.5)
{
  numberDown_ = numberOnDownSide;
  numberUp_ = numberOnUpSide;
  downList_ = CoinCopyOfArray(down, numberDown_);
  upList_ = CoinCopyOfArray(up, numberUp_);
}

// Copy constructor
CbcFixingBranchingObject::CbcFixingBranchingObject(const CbcFixingBranchingObject &rhs)
  : CbcBranchingObject(rhs)
{
  numberDown_ = rhs.numberDown_;
  numberUp_ = rhs.numberUp_;
  downList_ = CoinCopyOfArray(rhs.downList_, numberDown_);
  upList_ = CoinCopyOfArray(rhs.upList_, numberUp_);
}

// Assignment operator
CbcFixingBranchingObject &
CbcFixingBranchingObject::operator=(const CbcFixingBranchingObject &rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
    delete[] downList_;
    delete[] upList_;
    numberDown_ = rhs.numberDown_;
    numberUp_ = rhs.numberUp_;
    downList_ = CoinCopyOfArray(rhs.downList_, numberDown_);
    upList_ = CoinCopyOfArray(rhs.upList_, numberUp_);
  }
  return *this;
}
CbcBranchingObject *
CbcFixingBranchingObject::clone() const
{
  return (new CbcFixingBranchingObject(*this));
}

// Destructor
CbcFixingBranchingObject::~CbcFixingBranchingObject()
{
  delete[] downList_;
  delete[] upList_;
}
double
CbcFixingBranchingObject::branch()
{
  decrementNumberBranchesLeft();
  OsiSolverInterface *solver = model_->solver();
  const double *columnLower = solver->getColLower();
  int i;
  // *** for way - up means fix all those in up section
  if (way_ < 0) {
#ifdef FULL_PRINT
    printf("Down Fix ");
#endif
    //printf("Down Fix %d\n",numberDown_);
    for (i = 0; i < numberDown_; i++) {
      int iColumn = downList_[i];
      model_->solver()->setColUpper(iColumn, columnLower[iColumn]);
#ifdef FULL_PRINT
      printf("Setting bound on %d to lower bound\n", iColumn);
#endif
    }
    way_ = 1; // Swap direction
  } else {
#ifdef FULL_PRINT
    printf("Up Fix ");
#endif
    //printf("Up Fix %d\n",numberUp_);
    for (i = 0; i < numberUp_; i++) {
      int iColumn = upList_[i];
      model_->solver()->setColUpper(iColumn, columnLower[iColumn]);
#ifdef FULL_PRINT
      printf("Setting bound on %d to lower bound\n", iColumn);
#endif
    }
    way_ = -1; // Swap direction
  }
#ifdef FULL_PRINT
  printf("\n");
#endif
  return 0.0;
}
void CbcFixingBranchingObject::print()
{
  int i;
  // *** for way - up means fix all those in up section
  if (way_ < 0) {
    printf("Down Fix ");
    for (i = 0; i < numberDown_; i++) {
      int iColumn = downList_[i];
      printf("%d ", iColumn);
    }
  } else {
    printf("Up Fix ");
    for (i = 0; i < numberUp_; i++) {
      int iColumn = upList_[i];
      printf("%d ", iColumn);
    }
  }
  printf("\n");
}

/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type.
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int CbcFixingBranchingObject::compareOriginalObject(const CbcBranchingObject * /*brObj*/) const
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
CbcFixingBranchingObject::compareBranchingObject(const CbcBranchingObject * /*brObj*/, const bool /*replaceIfOverlap*/)
{
#ifdef JJF_ZERO //ndef NDEBUG
  const CbcFixingBranchingObject *br = dynamic_cast< const CbcFixingBranchingObject * >(brObj);
  assert(br);
#endif
  // If two FixingBranchingObject's have the same base object then it's pretty
  // much guaranteed
  throw("must implement");
}

//##############################################################################

//##############################################################################

// Default Constructor
CbcIdiotBranch::CbcIdiotBranch()
  : CbcObject()
{
  id_ = 1000000000 + CutBranchingObj;
}

// Useful constructor
CbcIdiotBranch::CbcIdiotBranch(CbcModel *model)
  : CbcObject(model)
{
  assert(model);
  id_ = 1000000000 + CutBranchingObj;
}

// Copy constructor
CbcIdiotBranch::CbcIdiotBranch(const CbcIdiotBranch &rhs)
  : CbcObject(rhs)
  , randomNumberGenerator_(rhs.randomNumberGenerator_)
  , savedRandomNumberGenerator_(rhs.savedRandomNumberGenerator_)
{
}

// Clone
CbcObject *
CbcIdiotBranch::clone() const
{
  return new CbcIdiotBranch(*this);
}

// Assignment operator
CbcIdiotBranch &
CbcIdiotBranch::operator=(const CbcIdiotBranch &rhs)
{
  if (this != &rhs) {
    CbcObject::operator=(rhs);
    randomNumberGenerator_ = rhs.randomNumberGenerator_;
    savedRandomNumberGenerator_ = rhs.savedRandomNumberGenerator_;
  }
  return *this;
}

// Destructor
CbcIdiotBranch::~CbcIdiotBranch()
{
}
double
CbcIdiotBranch::infeasibility(const OsiBranchingInformation *info,
  int &preferredWay) const
{
  randomNumberGenerator_ = savedRandomNumberGenerator_;
  double rhs = buildCut(info, 0, preferredWay).ub();
  double fraction = rhs - floor(rhs);
  if (fraction > 0.5)
    fraction = 1.0 - fraction;
  return fraction;
}

// This looks at solution and sets bounds to contain solution
void CbcIdiotBranch::feasibleRegion()
{
}

CbcBranchingObject *
CbcIdiotBranch::createCbcBranch(OsiSolverInterface *solver, const OsiBranchingInformation *info, int way)
{
  // down and up
  randomNumberGenerator_ = savedRandomNumberGenerator_;
  int preferredWay;
  OsiRowCut downCut = buildCut(info, 0, preferredWay);
  double rhs = downCut.ub();
  assert(rhs == downCut.lb());
  OsiRowCut upCut = downCut;
  downCut.setUb(floor(rhs));
  downCut.setLb(-COIN_DBL_MAX);
  upCut.setLb(ceil(rhs));
  upCut.setUb(COIN_DBL_MAX);
  CbcBranchingObject *branch
    = new CbcCutBranchingObject(model_, downCut, upCut, true);
  return branch;
}
// Initialize for branching
void CbcIdiotBranch::initializeForBranching(CbcModel *)
{
  savedRandomNumberGenerator_ = randomNumberGenerator_;
}
// Build "cut"
OsiRowCut
CbcIdiotBranch::buildCut(const OsiBranchingInformation *info, int /*type*/, int &preferredWay) const
{
  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();
  int *which = new int[numberIntegers];
  double *away = new double[numberIntegers];
  const double *solution = info->solution_;
  const double *lower = info->lower_;
  const double *upper = info->upper_;
  double integerTolerance = model_->getIntegerTolerance();
  //int nMax=std::min(4,numberIntegers/2);
  int n = 0;
  for (int i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    double value = solution[iColumn];
    value = std::max(value, lower[iColumn]);
    value = std::min(value, upper[iColumn]);
    //#define JUST_SMALL
#ifndef JUST_SMALL
    double nearest = floor(value + 0.5);
    if (fabs(value - nearest) > integerTolerance) {
      double random = 0.0; //0.25*randomNumberGenerator_.randomDouble();
      away[n] = -fabs(value - nearest) * (1.0 + random);
      which[n++] = iColumn;
    }
#else
    double distanceDown = value - floor(value);
    if (distanceDown > 10.0 * integerTolerance && distanceDown < 0.1) {
      double random = 0.0; //0.25*randomNumberGenerator_.randomDouble();
      away[n] = distanceDown * (1.0 + random);
      which[n++] = iColumn;
    }
#endif
  }
  CoinSort_2(away, away + n, which);
  OsiRowCut possibleCut;
  possibleCut.setUb(0.0);
  if (n > 1) {
    int nUse = 0;
    double useRhs = 0.0;
    double best = 0.0;
    double rhs = 0.0;
    double scaleFactor = 1.0;
    /* should be able to be clever to find best cut
       maybe don't do if away >0.45
       maybe don't be random
       maybe do complete enumeration on biggest k (where sum of k >= 1.0)
       maybe allow +-1
     */
    for (int i = 0; i < n; i++) {
      int iColumn = which[i];
      double value = solution[iColumn];
      value = std::max(value, lower[iColumn]);
      value = std::min(value, upper[iColumn]);
#define PLUS_MINUS
#ifndef PLUS_MINUS
      away[i] = 1.0;
      rhs += value;
#else
      if (value - floor(value) <= 0.5) {
        away[i] = 1.0;
        rhs += value;
      } else {
        away[i] = -1.0;
        rhs -= value;
      }
#endif
      double nearest = floor(rhs + 0.5);
      double distance = fabs(rhs - nearest);
      // scale back
      distance *= scaleFactor;
      scaleFactor *= 0.95;
      if (distance > best) {
        nUse = i + 1;
        best = distance;
        useRhs = rhs;
      }
    }
    // create cut
    if (nUse > 1) {
      possibleCut.setRow(nUse, which, away);
      possibleCut.setLb(useRhs);
      possibleCut.setUb(useRhs);
    }
  }
  delete[] which;
  delete[] away;
  return possibleCut;
}
#if 0
//##############################################################################

// Default Constructor
CbcBoundingBranchingObject::CbcBoundingBranchingObject()
        : CbcBranchingObject()
{
  branchCut_.setUb(0.0);
}

// Useful constructor
CbcBoundingBranchingObject::CbcBoundingBranchingObject (CbcModel * model,
							int way , const OsiRowCut * cut)
        : CbcBranchingObject(model, 0, way, 0.5)
{
  branchCut_ = *cut;
}

// Copy constructor
CbcBoundingBranchingObject::CbcBoundingBranchingObject ( const CbcBoundingBranchingObject & rhs) : CbcBranchingObject(rhs)
{
  branchCut_ = rhs.branchCut_;
}

// Assignment operator
CbcBoundingBranchingObject &
CbcBoundingBranchingObject::operator=( const CbcBoundingBranchingObject & rhs)
{
    if (this != &rhs) {
        CbcBranchingObject::operator=(rhs);
	branchCut_ = rhs.branchCut_;
    }
    return *this;
}
CbcBranchingObject *
CbcBoundingBranchingObject::clone() const
{
    return (new CbcBoundingBranchingObject(*this));
}


// Destructor
CbcBoundingBranchingObject::~CbcBoundingBranchingObject ()
{
}
double
CbcBoundingBranchingObject::branch()
{
    decrementNumberBranchesLeft();
    OsiSolverInterface * solver = model_->solver();
    const double * columnLower = solver->getColLower();
    int i;
    // *** for way - up means bound all those in up section
    if (way_ < 0) {
#ifdef FULL_PRINT
        printf("Down Bound ");
#endif
        //printf("Down Bound %d\n",numberDown_);
        for (i = 0; i < numberDown_; i++) {
            int iColumn = downList_[i];
            model_->solver()->setColUpper(iColumn, columnLower[iColumn]);
#ifdef FULL_PRINT
            printf("Setting bound on %d to lower bound\n", iColumn);
#endif
        }
        way_ = 1;	  // Swap direction
    } else {
#ifdef FULL_PRINT
        printf("Up Bound ");
#endif
        //printf("Up Bound %d\n",numberUp_);
        for (i = 0; i < numberUp_; i++) {
            int iColumn = upList_[i];
            model_->solver()->setColUpper(iColumn, columnLower[iColumn]);
#ifdef FULL_PRINT
            printf("Setting bound on %d to lower bound\n", iColumn);
#endif
        }
        way_ = -1;	  // Swap direction
    }
#ifdef FULL_PRINT
    printf("\n");
#endif
    return 0.0;
}
void
CbcBoundingBranchingObject::print()
{
  OsiRowCut cut = branchCut_;
  if (way_ < 0) {
    printf("Down Fix ");
    cut.setUb(floor(branchCut_.ub()));
    cut.setLb(-COIN_DBL_MAX);
  } else {
    printf("Up Fix ");
    cut.setLb(ceil(branchCut_.lb()));
    cut.setUb(COIN_DBL_MAX);
  }
  printf("\n");
  cut.print();
}

/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type.
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int
CbcBoundingBranchingObject::compareOriginalObject
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
CbcBoundingBranchingObject::compareBranchingObject
(const CbcBranchingObject* /*brObj*/, const bool /*replaceIfOverlap*/)
{
#ifdef JJF_ZERO //ndef NDEBUG
    const CbcBoundingBranchingObject* br =
        dynamic_cast<const CbcBoundingBranchingObject*>(brObj);
    assert(br);
#endif
    // If two BoundingBranchingObject's have the same base object then it's pretty
    // much guaranteed
    throw("must implement");
}

//##############################################################################

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
