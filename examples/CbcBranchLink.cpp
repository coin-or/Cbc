// $Id$
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "CoinPragma.hpp"
#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchLink.hpp"
#include "CoinError.hpp"
#include "CoinPackedMatrix.hpp"

// Default Constructor
CbcLink::CbcLink()
  : CbcObject()
  , weights_(NULL)
  , numberMembers_(0)
  , numberLinks_(0)
  , which_(NULL)
  , sosType_(1)
{
}

// Useful constructor (which are indices)
CbcLink::CbcLink(CbcModel *model, int numberMembers,
  int numberLinks, int first, const double *weights, int identifier)
  : CbcObject(model)
  , numberMembers_(numberMembers)
  , numberLinks_(numberLinks)
  , which_(NULL)
  , sosType_(1)
{
  id_ = identifier;
  if (numberMembers_) {
    weights_ = new double[numberMembers_];
    which_ = new int[numberMembers_ * numberLinks_];
    if (weights) {
      memcpy(weights_, weights, numberMembers_ * sizeof(double));
    } else {
      for (int i = 0; i < numberMembers_; i++)
        weights_[i] = i;
    }
    // weights must be increasing
    int i;
    for (i = 0; i < numberMembers_; i++)
      assert(i == 0 || weights_[i] > weights_[i - 1] + 1.0e-12);
    for (i = 0; i < numberMembers_ * numberLinks_; i++) {
      which_[i] = first + i;
    }
  } else {
    weights_ = NULL;
  }
}

// Useful constructor (which are indices)
CbcLink::CbcLink(CbcModel *model, int numberMembers,
  int numberLinks, int sosType, const int *which, const double *weights, int identifier)
  : CbcObject(model)
  , numberMembers_(numberMembers)
  , numberLinks_(numberLinks)
  , which_(NULL)
  , sosType_(sosType)
{
  id_ = identifier;
  if (numberMembers_) {
    weights_ = new double[numberMembers_];
    which_ = new int[numberMembers_ * numberLinks_];
    if (weights) {
      memcpy(weights_, weights, numberMembers_ * sizeof(double));
    } else {
      for (int i = 0; i < numberMembers_; i++)
        weights_[i] = i;
    }
    // weights must be increasing
    int i;
    for (i = 0; i < numberMembers_; i++)
      assert(i == 0 || weights_[i] > weights_[i - 1] + 1.0e-12);
    for (i = 0; i < numberMembers_ * numberLinks_; i++) {
      which_[i] = which[i];
    }
  } else {
    weights_ = NULL;
  }
}

// Copy constructor
CbcLink::CbcLink(const CbcLink &rhs)
  : CbcObject(rhs)
{
  numberMembers_ = rhs.numberMembers_;
  numberLinks_ = rhs.numberLinks_;
  sosType_ = rhs.sosType_;
  if (numberMembers_) {
    weights_ = CoinCopyOfArray(rhs.weights_, numberMembers_);
    which_ = CoinCopyOfArray(rhs.which_, numberMembers_ * numberLinks_);
  } else {
    weights_ = NULL;
    which_ = NULL;
  }
}

// Clone
CbcObject *
CbcLink::clone() const
{
  return new CbcLink(*this);
}

// Assignment operator
CbcLink &
CbcLink::operator=(const CbcLink &rhs)
{
  if (this != &rhs) {
    CbcObject::operator=(rhs);
    delete[] weights_;
    delete[] which_;
    numberMembers_ = rhs.numberMembers_;
    numberLinks_ = rhs.numberLinks_;
    sosType_ = rhs.sosType_;
    if (numberMembers_) {
      weights_ = CoinCopyOfArray(rhs.weights_, numberMembers_);
      which_ = CoinCopyOfArray(rhs.which_, numberMembers_ * numberLinks_);
    } else {
      weights_ = NULL;
      which_ = NULL;
    }
  }
  return *this;
}

// Destructor
CbcLink::~CbcLink()
{
  delete[] weights_;
  delete[] which_;
}

// Infeasibility - large is 0.5
double
CbcLink::infeasibility(int &preferredWay) const
{
  int j;
  int firstNonZero = -1;
  int lastNonZero = -1;
  OsiSolverInterface *solver = model_->solver();
  const double *solution = model_->testSolution();
  //const double * lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double weight = 0.0;
  double sum = 0.0;

  // check bounds etc
  double lastWeight = -1.0e100;
  int base = 0;
  for (j = 0; j < numberMembers_; j++) {
    for (int k = 0; k < numberLinks_; k++) {
      int iColumn = which_[base + k];
      //if (lower[iColumn])
      //throw CoinError("Non zero lower bound in CBCLink","infeasibility","CbcLink");
      if (lastWeight >= weights_[j] - 1.0e-7)
        throw CoinError("Weights too close together in CBCLink", "infeasibility", "CbcLink");
      double value = CoinMax(0.0, solution[iColumn]);
      sum += value;
      if (value > integerTolerance && upper[iColumn]) {
        // Possibly due to scaling a fixed variable might slip through
        if (value > upper[iColumn] + 1.0e-8) {
          // Could change to #ifdef CBC_DEBUG
#ifndef NDEBUG
          if (model_->messageHandler()->logLevel() > 1)
            printf("** Variable %d (%d) has value %g and upper bound of %g\n",
              iColumn, j, value, upper[iColumn]);
#endif
        }
        value = CoinMin(value, upper[iColumn]);
        weight += weights_[j] * value;
        if (firstNonZero < 0)
          firstNonZero = j;
        lastNonZero = j;
      }
    }
    base += numberLinks_;
  }
  double valueInfeasibility;
  preferredWay = 1;
  if (lastNonZero - firstNonZero >= sosType_) {
    // find where to branch
    assert(sum > 0.0);
    weight /= sum;
    valueInfeasibility = lastNonZero - firstNonZero + 1;
    valueInfeasibility *= 0.5 / ((double)numberMembers_);
    //#define DISTANCE
#ifdef DISTANCE
    assert(sosType_ == 1); // code up
    /* may still be satisfied.
       For LOS type 2 we might wish to move coding around
       and keep initial info in model_ for speed
    */
    int iWhere;
    bool possible = false;
    for (iWhere = firstNonZero; iWhere <= lastNonZero; iWhere++) {
      if (fabs(weight - weights_[iWhere]) < 1.0e-8) {
        possible = true;
        break;
      }
    }
    if (possible) {
      // One could move some of this (+ arrays) into model_
      const CoinPackedMatrix *matrix = solver->getMatrixByCol();
      const double *element = matrix->getMutableElements();
      const int *row = matrix->getIndices();
      const CoinBigIndex *columnStart = matrix->getVectorStarts();
      const int *columnLength = matrix->getVectorLengths();
      const double *rowSolution = solver->getRowActivity();
      const double *rowLower = solver->getRowLower();
      const double *rowUpper = solver->getRowUpper();
      int numberRows = matrix->getNumRows();
      double *array = new double[numberRows];
      CoinZeroN(array, numberRows);
      int *which = new int[numberRows];
      int n = 0;
      int base = numberLinks_ * firstNonZero;
      for (j = firstNonZero; j <= lastNonZero; j++) {
        for (int k = 0; k < numberLinks_; k++) {
          int iColumn = which_[base + k];
          double value = CoinMax(0.0, solution[iColumn]);
          if (value > integerTolerance && upper[iColumn]) {
            value = CoinMin(value, upper[iColumn]);
            for (int j = columnStart[iColumn]; j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              double a = array[iRow];
              if (a) {
                a += value * element[j];
                if (!a)
                  a = 1.0e-100;
              } else {
                which[n++] = iRow;
                a = value * element[j];
                assert(a);
              }
              array[iRow] = a;
            }
          }
        }
        base += numberLinks_;
      }
      base = numberLinks_ * iWhere;
      for (int k = 0; k < numberLinks_; k++) {
        int iColumn = which_[base + k];
        const double value = 1.0;
        for (int j = columnStart[iColumn]; j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          double a = array[iRow];
          if (a) {
            a -= value * element[j];
            if (!a)
              a = 1.0e-100;
          } else {
            which[n++] = iRow;
            a = -value * element[j];
            assert(a);
          }
          array[iRow] = a;
        }
      }
      for (j = 0; j < n; j++) {
        int iRow = which[j];
        // moving to point will increase row solution by this
        double distance = array[iRow];
        if (distance > 1.0e-8) {
          if (distance + rowSolution[iRow] > rowUpper[iRow] + 1.0e-8) {
            possible = false;
            break;
          }
        } else if (distance < -1.0e-8) {
          if (distance + rowSolution[iRow] < rowLower[iRow] - 1.0e-8) {
            possible = false;
            break;
          }
        }
      }
      for (j = 0; j < n; j++)
        array[which[j]] = 0.0;
      delete[] array;
      delete[] which;
      if (possible) {
        valueInfeasibility = 0.0;
        printf("possible %d %d %d\n", firstNonZero, lastNonZero, iWhere);
      }
    }
#endif
  } else {
    valueInfeasibility = 0.0; // satisfied
  }
  return valueInfeasibility;
}

// This looks at solution and sets bounds to contain solution
void CbcLink::feasibleRegion()
{
  int j;
  int firstNonZero = -1;
  int lastNonZero = -1;
  OsiSolverInterface *solver = model_->solver();
  const double *solution = model_->testSolution();
  const double *upper = solver->getColUpper();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double weight = 0.0;
  double sum = 0.0;

  int base = 0;
  for (j = 0; j < numberMembers_; j++) {
    for (int k = 0; k < numberLinks_; k++) {
      int iColumn = which_[base + k];
      double value = CoinMax(0.0, solution[iColumn]);
      sum += value;
      if (value > integerTolerance && upper[iColumn]) {
        weight += weights_[j] * value;
        if (firstNonZero < 0)
          firstNonZero = j;
        lastNonZero = j;
      }
    }
    base += numberLinks_;
  }
#ifdef DISTANCE
  if (lastNonZero - firstNonZero > sosType_ - 1) {
    /* may still be satisfied.
       For LOS type 2 we might wish to move coding around
       and keep initial info in model_ for speed
    */
    int iWhere;
    bool possible = false;
    for (iWhere = firstNonZero; iWhere <= lastNonZero; iWhere++) {
      if (fabs(weight - weights_[iWhere]) < 1.0e-8) {
        possible = true;
        break;
      }
    }
    if (possible) {
      // One could move some of this (+ arrays) into model_
      const CoinPackedMatrix *matrix = solver->getMatrixByCol();
      const double *element = matrix->getMutableElements();
      const int *row = matrix->getIndices();
      const CoinBigIndex *columnStart = matrix->getVectorStarts();
      const int *columnLength = matrix->getVectorLengths();
      const double *rowSolution = solver->getRowActivity();
      const double *rowLower = solver->getRowLower();
      const double *rowUpper = solver->getRowUpper();
      int numberRows = matrix->getNumRows();
      double *array = new double[numberRows];
      CoinZeroN(array, numberRows);
      int *which = new int[numberRows];
      int n = 0;
      int base = numberLinks_ * firstNonZero;
      for (j = firstNonZero; j <= lastNonZero; j++) {
        for (int k = 0; k < numberLinks_; k++) {
          int iColumn = which_[base + k];
          double value = CoinMax(0.0, solution[iColumn]);
          if (value > integerTolerance && upper[iColumn]) {
            value = CoinMin(value, upper[iColumn]);
            for (int j = columnStart[iColumn]; j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              double a = array[iRow];
              if (a) {
                a += value * element[j];
                if (!a)
                  a = 1.0e-100;
              } else {
                which[n++] = iRow;
                a = value * element[j];
                assert(a);
              }
              array[iRow] = a;
            }
          }
        }
        base += numberLinks_;
      }
      base = numberLinks_ * iWhere;
      for (int k = 0; k < numberLinks_; k++) {
        int iColumn = which_[base + k];
        const double value = 1.0;
        for (int j = columnStart[iColumn]; j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          double a = array[iRow];
          if (a) {
            a -= value * element[j];
            if (!a)
              a = 1.0e-100;
          } else {
            which[n++] = iRow;
            a = -value * element[j];
            assert(a);
          }
          array[iRow] = a;
        }
      }
      for (j = 0; j < n; j++) {
        int iRow = which[j];
        // moving to point will increase row solution by this
        double distance = array[iRow];
        if (distance > 1.0e-8) {
          if (distance + rowSolution[iRow] > rowUpper[iRow] + 1.0e-8) {
            possible = false;
            break;
          }
        } else if (distance < -1.0e-8) {
          if (distance + rowSolution[iRow] < rowLower[iRow] - 1.0e-8) {
            possible = false;
            break;
          }
        }
      }
      for (j = 0; j < n; j++)
        array[which[j]] = 0.0;
      delete[] array;
      delete[] which;
      if (possible) {
        printf("possible feas region %d %d %d\n", firstNonZero, lastNonZero, iWhere);
        firstNonZero = iWhere;
        lastNonZero = iWhere;
      }
    }
  }
#else
  assert(lastNonZero - firstNonZero < sosType_);
#endif
  base = 0;
  for (j = 0; j < firstNonZero; j++) {
    for (int k = 0; k < numberLinks_; k++) {
      int iColumn = which_[base + k];
      solver->setColUpper(iColumn, 0.0);
    }
    base += numberLinks_;
  }
  // skip
  base += numberLinks_;
  for (j = lastNonZero + 1; j < numberMembers_; j++) {
    for (int k = 0; k < numberLinks_; k++) {
      int iColumn = which_[base + k];
      solver->setColUpper(iColumn, 0.0);
    }
    base += numberLinks_;
  }
}

// Creates a branching object
CbcBranchingObject *
CbcLink::createCbcBranch(OsiSolverInterface * /*solver*/, const OsiBranchingInformation * /*info*/, int way)
{
  int j;
  const double *solution = model_->testSolution();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  OsiSolverInterface *solver = model_->solver();
  const double *upper = solver->getColUpper();
  int firstNonFixed = -1;
  int lastNonFixed = -1;
  int firstNonZero = -1;
  int lastNonZero = -1;
  double weight = 0.0;
  double sum = 0.0;
  int base = 0;
  for (j = 0; j < numberMembers_; j++) {
    for (int k = 0; k < numberLinks_; k++) {
      int iColumn = which_[base + k];
      if (upper[iColumn]) {
        double value = CoinMax(0.0, solution[iColumn]);
        sum += value;
        if (firstNonFixed < 0)
          firstNonFixed = j;
        lastNonFixed = j;
        if (value > integerTolerance) {
          weight += weights_[j] * value;
          if (firstNonZero < 0)
            firstNonZero = j;
          lastNonZero = j;
        }
      }
    }
    base += numberLinks_;
  }
  assert(lastNonZero - firstNonZero >= sosType_);
  // find where to branch
  assert(sum > 0.0);
  weight /= sum;
  int iWhere;
  double separator = 0.0;
  for (iWhere = firstNonZero; iWhere < lastNonZero; iWhere++)
    if (weight < weights_[iWhere + 1])
      break;
  if (sosType_ == 1) {
    // SOS 1
    separator = 0.5 * (weights_[iWhere] + weights_[iWhere + 1]);
  } else {
    // SOS 2
    if (iWhere == firstNonFixed)
      iWhere++;
    ;
    if (iWhere == lastNonFixed - 1)
      iWhere = lastNonFixed - 2;
    separator = weights_[iWhere + 1];
  }
  // create object
  CbcBranchingObject *branch;
  branch = new CbcLinkBranchingObject(model_, this, way, separator);
  return branch;
}
// Useful constructor
CbcLinkBranchingObject::CbcLinkBranchingObject(CbcModel *model,
  const CbcLink *set,
  int way,
  double separator)
  : CbcBranchingObject(model, set->id(), way, 0.5)
{
  set_ = set;
  separator_ = separator;
}

// Copy constructor
CbcLinkBranchingObject::CbcLinkBranchingObject(const CbcLinkBranchingObject &rhs)
  : CbcBranchingObject(rhs)
{
  set_ = rhs.set_;
  separator_ = rhs.separator_;
}

// Assignment operator
CbcLinkBranchingObject &
CbcLinkBranchingObject::operator=(const CbcLinkBranchingObject &rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
    set_ = rhs.set_;
    separator_ = rhs.separator_;
  }
  return *this;
}
CbcBranchingObject *
CbcLinkBranchingObject::clone() const
{
  return (new CbcLinkBranchingObject(*this));
}

// Destructor
CbcLinkBranchingObject::~CbcLinkBranchingObject()
{
}
double
CbcLinkBranchingObject::branch()
{
  decrementNumberBranchesLeft();
  int numberMembers = set_->numberMembers();
  int numberLinks = set_->numberLinks();
  const double *weights = set_->weights();
  const int *which = set_->which();
  OsiSolverInterface *solver = model_->solver();
  //const double * lower = solver->getColLower();
  //const double * upper = solver->getColUpper();
  // *** for way - up means fix all those in down section
  if (way_ < 0) {
    int i;
    for (i = 0; i < numberMembers; i++) {
      if (weights[i] > separator_)
        break;
    }
    assert(i < numberMembers);
    int base = i * numberLinks;
    ;
    for (; i < numberMembers; i++) {
      for (int k = 0; k < numberLinks; k++) {
        int iColumn = which[base + k];
        solver->setColUpper(iColumn, 0.0);
      }
      base += numberLinks;
    }
    way_ = 1; // Swap direction
  } else {
    int i;
    int base = 0;
    for (i = 0; i < numberMembers; i++) {
      if (weights[i] >= separator_) {
        break;
      } else {
        for (int k = 0; k < numberLinks; k++) {
          int iColumn = which[base + k];
          solver->setColUpper(iColumn, 0.0);
        }
        base += numberLinks;
      }
    }
    assert(i < numberMembers);
    way_ = -1; // Swap direction
  }
  return 0.0;
}
// Print what would happen
void CbcLinkBranchingObject::print()
{
  int numberMembers = set_->numberMembers();
  int numberLinks = set_->numberLinks();
  const double *weights = set_->weights();
  const int *which = set_->which();
  OsiSolverInterface *solver = model_->solver();
  const double *upper = solver->getColUpper();
  int first = numberMembers;
  int last = -1;
  int numberFixed = 0;
  int numberOther = 0;
  int i;
  int base = 0;
  for (i = 0; i < numberMembers; i++) {
    for (int k = 0; k < numberLinks; k++) {
      int iColumn = which[base + k];
      double bound = upper[iColumn];
      if (bound) {
        first = CoinMin(first, i);
        last = CoinMax(last, i);
      }
    }
    base += numberLinks;
  }
  // *** for way - up means fix all those in down section
  base = 0;
  if (way_ < 0) {
    printf("SOS Down");
    for (i = 0; i < numberMembers; i++) {
      if (weights[i] > separator_)
        break;
      for (int k = 0; k < numberLinks; k++) {
        int iColumn = which[base + k];
        double bound = upper[iColumn];
        if (bound)
          numberOther++;
      }
      base += numberLinks;
    }
    assert(i < numberMembers);
    for (; i < numberMembers; i++) {
      for (int k = 0; k < numberLinks; k++) {
        int iColumn = which[base + k];
        double bound = upper[iColumn];
        if (bound)
          numberFixed++;
      }
      base += numberLinks;
    }
  } else {
    printf("SOS Up");
    for (i = 0; i < numberMembers; i++) {
      if (weights[i] >= separator_)
        break;
      for (int k = 0; k < numberLinks; k++) {
        int iColumn = which[base + k];
        double bound = upper[iColumn];
        if (bound)
          numberFixed++;
      }
      base += numberLinks;
    }
    assert(i < numberMembers);
    for (; i < numberMembers; i++) {
      for (int k = 0; k < numberLinks; k++) {
        int iColumn = which[base + k];
        double bound = upper[iColumn];
        if (bound)
          numberOther++;
      }
      base += numberLinks;
    }
  }
  assert((numberFixed % numberLinks) == 0);
  assert((numberOther % numberLinks) == 0);
  printf(" - at %g, free range %d (%g) => %d (%g), %d would be fixed, %d other way\n",
    separator_, first, weights[first], last, weights[last], numberFixed / numberLinks,
    numberOther / numberLinks);
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
CbcLinkBranchingObject::compareBranchingObject(const CbcBranchingObject *brObj, const bool replaceIfOverlap)
{
  throw("must implement");
}
