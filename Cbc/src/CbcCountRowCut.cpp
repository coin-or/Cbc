/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif
#include <cassert>

#include "OsiRowCut.hpp"
#include "CbcModel.hpp"
#include "CbcCountRowCut.hpp"
#include "CbcNode.hpp"
//#define CHECK_CUT_COUNTS
// Default Constructor
CbcCountRowCut::CbcCountRowCut()
  : OsiRowCut()
  , owner_(NULL)
  , ownerCut_(-1)
  , numberPointingToThis_(0)
  , whichCutGenerator_(-1)
{
#ifdef CHECK_CUT_COUNTS
  printf("CbcCountRowCut default constructor %x\n", this);
#endif
}

// Copy Constructor
CbcCountRowCut::CbcCountRowCut(const OsiRowCut &rhs)
  : OsiRowCut(rhs)
  , owner_(NULL)
  , ownerCut_(-1)
  , numberPointingToThis_(0)
  , whichCutGenerator_(-1)
{
#ifdef CHECK_CUT_COUNTS
  printf("CbcCountRowCut constructor %x from RowCut\n", this);
#endif
}
// Copy Constructor
CbcCountRowCut::CbcCountRowCut(const OsiRowCut &rhs,
  CbcNodeInfo *info, int whichOne,
  int whichGenerator,
  int numberPointingToThis)
  : OsiRowCut(rhs)
  , owner_(info)
  , ownerCut_(whichOne)
  , numberPointingToThis_(numberPointingToThis)
  , whichCutGenerator_(whichGenerator)
{
#ifdef CHECK_CUT_COUNTS
  printf("CbcCountRowCut constructor %x from RowCut and info %d\n",
    this, numberPointingToThis_);
#endif
  //assert (!numberPointingToThis||numberPointingToThis==1000000000);
}
CbcCountRowCut::~CbcCountRowCut()
{
#ifdef CHECK_CUT_COUNTS
  printf("CbcCountRowCut destructor %x - references %d\n", this,
    numberPointingToThis_);
#endif
  // Look at owner and delete
  if (owner_)
    owner_->deleteCut(ownerCut_);
  ownerCut_ = -1234567;
}
// Increment number of references
void CbcCountRowCut::increment(int change)
{
  assert(ownerCut_ != -1234567);
  numberPointingToThis_ += change;
}

// Decrement number of references and return number left
int CbcCountRowCut::decrement(int change)
{
  assert(ownerCut_ != -1234567);
  // See if plausible number
  if (change < 900000000) {
    //assert(numberPointingToThis_>=change);
    assert(numberPointingToThis_ >= 0);
    if (numberPointingToThis_ < change) {
      assert(numberPointingToThis_ > 0);
      COIN_DETAIL_PRINT(printf("negative cut count %d - %d\n", numberPointingToThis_, change));
      change = numberPointingToThis_;
    }
    numberPointingToThis_ -= change;
  }
  return numberPointingToThis_;
}

// Set information
void CbcCountRowCut::setInfo(CbcNodeInfo *info, int whichOne)
{
  owner_ = info;
  ownerCut_ = whichOne;
}
// Returns true if can drop cut if slack basic
bool CbcCountRowCut::canDropCut(const OsiSolverInterface *solver, int iRow) const
{
  // keep if COIN_DBL_MAX otherwise keep if slack zero
  if (effectiveness() < 1.0e20) {
    return true;
  } else if (effectiveness() != COIN_DBL_MAX) {
    if (iRow >= solver->getNumRows())
      return true;
    const double *rowActivity = solver->getRowActivity();
    const double *rowLower = solver->getRowLower();
    const double *rowUpper = solver->getRowUpper();
    double tolerance;
    solver->getDblParam(OsiPrimalTolerance, tolerance);
    double value = rowActivity[iRow];
    if (value < rowLower[iRow] + tolerance || value > rowUpper[iRow] - tolerance)
      return false;
    else
      return true;
  } else {
    return false;
  }
}
static double multiplier[] = { 1.23456789e2, -9.87654321 };
static int hashCut(const OsiRowCut2 &x, int size)
{
  int xN = x.row().getNumElements();
  double xLb = x.lb();
  double xUb = x.ub();
  const int *xIndices = x.row().getIndices();
  const double *xElements = x.row().getElements();
  unsigned int hashValue;
  double value = 1.0;
  if (xLb > -1.0e10)
    value += xLb * multiplier[0];
  if (xUb < 1.0e10)
    value += xUb * multiplier[1];
  for (int j = 0; j < xN; j++) {
    int xColumn = xIndices[j];
    double xValue = xElements[j];
    int k = (j & 1);
    value += (j + 1) * multiplier[k] * (xColumn + 1) * xValue;
  }
  // should be compile time but too lazy for now
  union {
    double d;
    unsigned int i[2];
  } xx;
  if (sizeof(value) > sizeof(hashValue)) {
    assert(sizeof(value) == 2 * sizeof(hashValue));
    xx.d = value;
    hashValue = (xx.i[0] + xx.i[1]);
  } else {
    assert(sizeof(value) == sizeof(hashValue));
    xx.d = value;
    hashValue = xx.i[0];
  }
  return hashValue % (size);
}
static int hashCut2(const OsiRowCut2 &x, int size)
{
  int xN = x.row().getNumElements();
  double xLb = x.lb();
  double xUb = x.ub();
  const int *xIndices = x.row().getIndices();
  const double *xElements = x.row().getElements();
  unsigned int hashValue;
  double value = 1.0;
  if (xLb > -1.0e10)
    value += xLb * multiplier[0];
  if (xUb < 1.0e10)
    value += xUb * multiplier[1];
  for (int j = 0; j < xN; j++) {
    int xColumn = xIndices[j];
    double xValue = xElements[j];
    int k = (j & 1);
    value += (j + 1) * multiplier[k] * (xColumn + 1) * xValue;
  }
  // should be compile time but too lazy for now
  if (sizeof(value) > sizeof(hashValue)) {
    assert(sizeof(value) == 2 * sizeof(hashValue));
    union {
      double d;
      unsigned int i[2];
    } xx;
    xx.d = value;
    hashValue = (xx.i[0] + xx.i[1]);
  } else {
    assert(sizeof(value) == sizeof(hashValue));
    union {
      double d;
      unsigned int i[2];
    } xx;
    xx.d = value;
    hashValue = xx.i[0];
  }
  return hashValue % (size);
}
static bool same(const OsiRowCut2 &x, const OsiRowCut2 &y)
{
  int xN = x.row().getNumElements();
  int yN = y.row().getNumElements();
  bool identical = false;
  if (xN == yN) {
    double xLb = x.lb();
    double xUb = x.ub();
    double yLb = y.lb();
    double yUb = y.ub();
    if (fabs(xLb - yLb) < 1.0e-8 && fabs(xUb - yUb) < 1.0e-8) {
      const int *xIndices = x.row().getIndices();
      const double *xElements = x.row().getElements();
      const int *yIndices = y.row().getIndices();
      const double *yElements = y.row().getElements();
      int j;
      for (j = 0; j < xN; j++) {
        if (xIndices[j] != yIndices[j])
          break;
        if (fabs(xElements[j] - yElements[j]) > 1.0e-12)
          break;
      }
      identical = (j == xN);
    }
  }
  return identical;
}
static bool same2(const OsiRowCut2 &x, const OsiRowCut2 &y)
{
  int xN = x.row().getNumElements();
  int yN = y.row().getNumElements();
  bool identical = false;
  if (xN == yN) {
    double xLb = x.lb();
    double xUb = x.ub();
    double yLb = y.lb();
    double yUb = y.ub();
    if (fabs(xLb - yLb) < 1.0e-8 && fabs(xUb - yUb) < 1.0e-8) {
      const int *xIndices = x.row().getIndices();
      const double *xElements = x.row().getElements();
      const int *yIndices = y.row().getIndices();
      const double *yElements = y.row().getElements();
      int j;
      for (j = 0; j < xN; j++) {
        if (xIndices[j] != yIndices[j])
          break;
        if (fabs(xElements[j] - yElements[j]) > 1.0e-12)
          break;
      }
      identical = (j == xN);
    }
  }
  return identical;
}
CbcRowCuts::CbcRowCuts(int initialMaxSize, int hashMultiplier)
{
  numberCuts_ = 0;
  size_ = initialMaxSize;
  hashMultiplier_ = hashMultiplier;
  int hashSize = hashMultiplier_ * size_;
  if (size_) {
    rowCut_ = new OsiRowCut2 *[size_];
    hash_ = new CoinHashLink[hashSize];
  } else {
    rowCut_ = NULL;
    hash_ = NULL;
  }
  for (int i = 0; i < hashSize; i++) {
    hash_[i].index = -1;
    hash_[i].next = -1;
  }
  lastHash_ = -1;
}
CbcRowCuts::~CbcRowCuts()
{
  for (int i = 0; i < numberCuts_; i++)
    delete rowCut_[i];
  delete[] rowCut_;
  delete[] hash_;
}
CbcRowCuts::CbcRowCuts(const CbcRowCuts &rhs)
{
  numberCuts_ = rhs.numberCuts_;
  hashMultiplier_ = rhs.hashMultiplier_;
  size_ = rhs.size_;
  int hashSize = size_ * hashMultiplier_;
  lastHash_ = rhs.lastHash_;
  if (size_) {
    rowCut_ = new OsiRowCut2 *[size_];
    hash_ = new CoinHashLink[hashSize];
    for (int i = 0; i < hashSize; i++) {
      hash_[i] = rhs.hash_[i];
    }
    for (int i = 0; i < numberCuts_; i++) {
      if (rhs.rowCut_[i])
        rowCut_[i] = new OsiRowCut2(*rhs.rowCut_[i]);
      else
        rowCut_[i] = NULL;
    }
  } else {
    rowCut_ = NULL;
    hash_ = NULL;
  }
}
CbcRowCuts &
CbcRowCuts::operator=(const CbcRowCuts &rhs)
{
  if (this != &rhs) {
    for (int i = 0; i < numberCuts_; i++)
      delete rowCut_[i];
    delete[] rowCut_;
    delete[] hash_;
    numberCuts_ = rhs.numberCuts_;
    hashMultiplier_ = rhs.hashMultiplier_;
    size_ = rhs.size_;
    lastHash_ = rhs.lastHash_;
    if (size_) {
      rowCut_ = new OsiRowCut2 *[size_];
      int hashSize = size_ * hashMultiplier_;
      hash_ = new CoinHashLink[hashSize];
      for (int i = 0; i < hashSize; i++) {
        hash_[i] = rhs.hash_[i];
      }
      for (int i = 0; i < numberCuts_; i++) {
        if (rhs.rowCut_[i])
          rowCut_[i] = new OsiRowCut2(*rhs.rowCut_[i]);
        else
          rowCut_[i] = NULL;
      }
    } else {
      rowCut_ = NULL;
      hash_ = NULL;
    }
  }
  return *this;
}
void CbcRowCuts::eraseRowCut(int sequence)
{
  // find
  assert(sequence >= 0 && sequence < numberCuts_);
  OsiRowCut2 *cut = rowCut_[sequence];
  int hashSize = size_ * hashMultiplier_;
  int ipos = hashCut(*cut, hashSize);
  int found = -1;
  while (true) {
    int j1 = hash_[ipos].index;
    if (j1 >= 0) {
      if (j1 != sequence) {
        int k = hash_[ipos].next;
        if (k != -1)
          ipos = k;
        else
          break;
      } else {
        found = j1;
        break;
      }
    } else {
      break;
    }
  }
  assert(found >= 0);
  assert(hash_[ipos].index == sequence);
  // shuffle up
  while (hash_[ipos].next >= 0) {
    int k = hash_[ipos].next;
    hash_[ipos] = hash_[k];
    ipos = k;
  }
  hash_[ipos].index = -1;
  // move last to found
  numberCuts_--;
  assert(found == numberCuts_); // debug when fails
  if (numberCuts_ && found < numberCuts_) {
    ipos = hashCut(*rowCut_[numberCuts_], hashSize);
    while (true) {
      int j1 = hash_[ipos].index;
      if (j1 != numberCuts_) {
        int k = hash_[ipos].next;
        ipos = k;
        assert(ipos >= 0);
      } else {
        // change
        hash_[ipos].index = found;
        rowCut_[found] = rowCut_[numberCuts_];
        rowCut_[numberCuts_] = NULL;
        break;
      }
    }
  }
  delete cut;
  rowCut_[numberCuts_] = NULL;
  //assert (!rowCut_[numberCuts_-1]);
}
// Truncate
void CbcRowCuts::truncate(int numberAfter)
{
  if (numberAfter < 0 || numberAfter >= numberCuts_)
    return;
  for (int i = numberAfter; i < numberCuts_; i++) {
    delete rowCut_[i];
    rowCut_[i] = NULL;
  }
  numberCuts_ = numberAfter;
  int hashSize = size_ * hashMultiplier_;
  for (int i = 0; i < hashSize; i++) {
    hash_[i].index = -1;
    hash_[i].next = -1;
  }
  OsiRowCut2 **temp = new OsiRowCut2 *[size_];
  lastHash_ = -1;
  for (int i = 0; i < numberCuts_; i++) {
    temp[i] = rowCut_[i];
    int ipos = hashCut(*temp[i], hashSize);
    int found = -1;
    int jpos = ipos;
    while (true) {
      int j1 = hash_[ipos].index;
      if (j1 >= 0) {
        if (!same(*temp[i], *temp[j1])) {
          int k = hash_[ipos].next;
          if (k != -1)
            ipos = k;
          else
            break;
        } else {
          found = j1;
          break;
        }
      } else {
        break;
      }
    }
    if (found < 0) {
      assert(hash_[ipos].next == -1);
      if (ipos == jpos) {
        // first
        hash_[ipos].index = i;
      } else {
        // find next space
        while (true) {
          ++lastHash_;
          assert(lastHash_ < hashSize);
          if (hash_[lastHash_].index == -1)
            break;
        }
        hash_[ipos].next = lastHash_;
        hash_[lastHash_].index = i;
      }
    }
  }
  delete[] rowCut_;
  rowCut_ = temp;
}
// Return 0 if added, 1 if not, -1 if not added because of space
int CbcRowCuts::addCutIfNotDuplicate(const OsiRowCut &cut, int whichType)
{
  int hashSize = size_ * hashMultiplier_;
  bool globallyValid = cut.globallyValid();
  if (numberCuts_ == size_) {
    size_ = 2 * size_ + 100;
    hashSize = hashMultiplier_ * size_;
    OsiRowCut2 **temp = new OsiRowCut2 *[size_];
    delete[] hash_;
    hash_ = new CoinHashLink[hashSize];
    for (int i = 0; i < hashSize; i++) {
      hash_[i].index = -1;
      hash_[i].next = -1;
    }
    lastHash_ = -1;
    for (int i = 0; i < numberCuts_; i++) {
      temp[i] = rowCut_[i];
      int ipos = hashCut(*temp[i], hashSize);
      int found = -1;
      int jpos = ipos;
      while (true) {
        int j1 = hash_[ipos].index;

        if (j1 >= 0) {
          if (!same(*temp[i], *temp[j1])) {
            int k = hash_[ipos].next;
            if (k != -1)
              ipos = k;
            else
              break;
          } else {
            found = j1;
            break;
          }
        } else {
          break;
        }
      }
      if (found < 0) {
        assert(hash_[ipos].next == -1);
        if (ipos == jpos) {
          // first
          hash_[ipos].index = i;
        } else {
          // find next space
          while (true) {
            ++lastHash_;
            assert(lastHash_ < hashSize);
            if (hash_[lastHash_].index == -1)
              break;
          }
          hash_[ipos].next = lastHash_;
          hash_[lastHash_].index = i;
        }
      }
    }
    delete[] rowCut_;
    rowCut_ = temp;
  }
  if (numberCuts_ < size_) {
    double newLb = cut.lb();
    double newUb = cut.ub();
    CoinPackedVector vector = cut.row();
    int numberElements = vector.getNumElements();
    int *newIndices = vector.getIndices();
    double *newElements = vector.getElements();
    CoinSort_2(newIndices, newIndices + numberElements, newElements);
    int i;
    bool bad = false;
    for (i = 0; i < numberElements; i++) {
      double value = fabs(newElements[i]);
      if (value < 1.0e-12 || value > 1.0e12)
        bad = true;
    }
    if (bad)
      return 1;
    OsiRowCut2 newCut(whichType);
    newCut.setLb(newLb);
    newCut.setUb(newUb);
    newCut.setRow(vector);
    int ipos = hashCut(newCut, hashSize);
    int found = -1;
    int jpos = ipos;
    while (true) {
      int j1 = hash_[ipos].index;

      if (j1 >= 0) {
        if (!same(newCut, *rowCut_[j1])) {
          int k = hash_[ipos].next;
          if (k != -1)
            ipos = k;
          else
            break;
        } else {
          found = j1;
          break;
        }
      } else {
        break;
      }
    }
    if (found < 0) {
      assert(hash_[ipos].next == -1);
      if (ipos == jpos) {
        // first
        hash_[ipos].index = numberCuts_;
      } else {
        // find next space
        while (true) {
          ++lastHash_;
          assert(lastHash_ < hashSize);
          if (hash_[lastHash_].index == -1)
            break;
        }
        hash_[ipos].next = lastHash_;
        hash_[lastHash_].index = numberCuts_;
      }
      OsiRowCut2 *newCutPtr = new OsiRowCut2(whichType);
      newCutPtr->setLb(newLb);
      newCutPtr->setUb(newUb);
      newCutPtr->setRow(vector);
      newCutPtr->setGloballyValid(globallyValid);
      rowCut_[numberCuts_++] = newCutPtr;
      //printf("addedGlobalCut of size %d to %x - cuts size %d\n",
      //     cut.row().getNumElements(),this,numberCuts_);
      return 0;
    } else {
      return 1;
    }
  } else {
    return -1;
  }
}
// Return 0 if added, 1 if not, -1 if not added because of space
int CbcRowCuts::addCutIfNotDuplicateWhenGreedy(const OsiRowCut &cut, int whichType)
{
  int hashSize = size_ * hashMultiplier_;
  if (numberCuts_ == size_) {
    size_ = 2 * size_ + 100;
    hashSize = hashMultiplier_ * size_;
    OsiRowCut2 **temp = new OsiRowCut2 *[size_];
    delete[] hash_;
    hash_ = new CoinHashLink[hashSize];
    for (int i = 0; i < hashSize; i++) {
      hash_[i].index = -1;
      hash_[i].next = -1;
    }
    lastHash_ = -1;
    for (int i = 0; i < numberCuts_; i++) {
      temp[i] = rowCut_[i];
      int ipos = hashCut2(*temp[i], hashSize);
      int found = -1;
      int jpos = ipos;
      while (true) {
        int j1 = hash_[ipos].index;

        if (j1 >= 0) {
          if (!same2(*temp[i], *temp[j1])) {
            int k = hash_[ipos].next;
            if (k != -1)
              ipos = k;
            else
              break;
          } else {
            found = j1;
            break;
          }
        } else {
          break;
        }
      }
      if (found < 0) {
        assert(hash_[ipos].next == -1);
        if (ipos == jpos) {
          // first
          hash_[ipos].index = i;
        } else {
          // find next space
          while (true) {
            ++lastHash_;
            assert(lastHash_ < hashSize);
            if (hash_[lastHash_].index == -1)
              break;
          }
          hash_[ipos].next = lastHash_;
          hash_[lastHash_].index = i;
        }
      }
    }
    delete[] rowCut_;
    rowCut_ = temp;
  }
  if (numberCuts_ < size_) {
    double newLb = cut.lb();
    double newUb = cut.ub();
    CoinPackedVector vector = cut.row();
    int numberElements = vector.getNumElements();
    int *newIndices = vector.getIndices();
    double *newElements = vector.getElements();
    CoinSort_2(newIndices, newIndices + numberElements, newElements);
    int i;
    bool bad = false;
    for (i = 0; i < numberElements; i++) {
      double value = fabs(newElements[i]);
      if (value < 1.0e-12 || value > 1.0e12)
        bad = true;
    }
    if (bad)
      return 1;
    OsiRowCut2 newCut(whichType);
    newCut.setLb(newLb);
    newCut.setUb(newUb);
    newCut.setRow(vector);
    int ipos = hashCut2(newCut, hashSize);
    int found = -1;
    int jpos = ipos;
    while (true) {
      int j1 = hash_[ipos].index;

      if (j1 >= 0) {
        if (!same2(newCut, *rowCut_[j1])) {
          int k = hash_[ipos].next;
          if (k != -1)
            ipos = k;
          else
            break;
        } else {
          found = j1;
          break;
        }
      } else {
        break;
      }
    }
    if (found < 0) {
      assert(hash_[ipos].next == -1);
      if (ipos == jpos) {
        // first
        hash_[ipos].index = numberCuts_;
      } else {
        // find next space
        while (true) {
          ++lastHash_;
          assert(lastHash_ < hashSize);
          if (hash_[lastHash_].index == -1)
            break;
        }
        hash_[ipos].next = lastHash_;
        hash_[lastHash_].index = numberCuts_;
      }
      OsiRowCut2 *newCutPtr = new OsiRowCut2(whichType);
      newCutPtr->setLb(newLb);
      newCutPtr->setUb(newUb);
      newCutPtr->setRow(vector);
      rowCut_[numberCuts_++] = newCutPtr;
      //printf("addedGreedyGlobalCut of size %d to %p - cuts size %d\n",
      //     cut.row().getNumElements(),this,numberCuts_);
      return 0;
    } else {
      return 1;
    }
  } else {
    return -1;
  }
}
// Add in cuts as normal cuts and delete
void CbcRowCuts::addCuts(OsiCuts &cs)
{
  for (int i = 0; i < numberCuts_; i++) {
    cs.insert(*rowCut_[i]);
    delete rowCut_[i];
    rowCut_[i] = NULL;
  }
  numberCuts_ = 0;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
