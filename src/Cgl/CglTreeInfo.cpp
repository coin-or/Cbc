// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <climits>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>

#include "CoinPragma.hpp"
#include "CglTreeInfo.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinSort.hpp"
#include "CoinPackedMatrix.hpp"
#include "CglStored.hpp"
#include "OsiRowCut.hpp"

// Default constructor
CglTreeInfo::CglTreeInfo()
  : level(-1)
  , pass(-1)
  , formulation_rows(-1)
  , options(0)
  , inTree(false)
  , hasParent(0)
  , parentSolver(NULL)
  , originalColumns(NULL)
  , strengthenRow(NULL)
  , randomNumberGenerator(NULL)
{
}

// Copy constructor
CglTreeInfo::CglTreeInfo(const CglTreeInfo &rhs)
  : level(rhs.level)
  , pass(rhs.pass)
  , formulation_rows(rhs.formulation_rows)
  , options(rhs.options)
  , inTree(rhs.inTree)
  , hasParent(rhs.hasParent)
  , parentSolver(rhs.parentSolver)
  , originalColumns(rhs.originalColumns)
  , strengthenRow(rhs.strengthenRow)
  , randomNumberGenerator(rhs.randomNumberGenerator)
{
}
// Clone
CglTreeInfo *
CglTreeInfo::clone() const
{
  return new CglTreeInfo(*this);
}

// Assignment operator
CglTreeInfo &
CglTreeInfo::operator=(const CglTreeInfo &rhs)
{
  if (this != &rhs) {
    //CglCutGenerator::operator=(rhs);
    level = rhs.level;
    pass = rhs.pass;
    formulation_rows = rhs.formulation_rows;
    options = rhs.options;
    inTree = rhs.inTree;
    hasParent = rhs.hasParent;
    parentSolver = rhs.parentSolver;
    originalColumns = rhs.originalColumns;
    strengthenRow = rhs.strengthenRow;
    randomNumberGenerator = rhs.randomNumberGenerator;
  }
  return *this;
}

// Destructor
CglTreeInfo::~CglTreeInfo()
{
}

// Default constructor
CglTreeProbingInfo::CglTreeProbingInfo()
  : CglTreeInfo()
  , fixEntry_(NULL)
  , toZero_(NULL)
  , toOne_(NULL)
  , integerVariable_(NULL)
  , backward_(NULL)
  , fixingEntry_(NULL)
  , numberVariables_(0)
  , numberIntegers_(0)
  , maximumEntries_(0)
  , numberEntries_(-1)
{
}
// Constructor from model
CglTreeProbingInfo::CglTreeProbingInfo(const OsiSolverInterface *model)
  : CglTreeInfo()
  , fixEntry_(NULL)
  , toZero_(NULL)
  , toOne_(NULL)
  , integerVariable_(NULL)
  , backward_(NULL)
  , fixingEntry_(NULL)
  , numberVariables_(0)
  , numberIntegers_(0)
  , maximumEntries_(0)
  , numberEntries_(-1)
{
  numberVariables_ = model->getNumCols();
  // Too many ... but
  integerVariable_ = new int[numberVariables_];
  backward_ = new int[numberVariables_];
  int i;
  // Get integer types
  const char *columnType = model->getColType(true);
  for (i = 0; i < numberVariables_; i++) {
    backward_[i] = -1;
    if (columnType[i]) {
      if (columnType[i] == 1) {
        backward_[i] = numberIntegers_;
        integerVariable_[numberIntegers_++] = i;
      } else {
        backward_[i] = -2;
      }
    }
  }
  // Set up to arrays
  toOne_ = new int[numberIntegers_];
  toZero_ = new int[numberIntegers_ + 1];
  // zero out
  CoinZeroN(toOne_, numberIntegers_);
  CoinZeroN(toZero_, numberIntegers_ + 1);
}

// Copy constructor
CglTreeProbingInfo::CglTreeProbingInfo(const CglTreeProbingInfo &rhs)
  : CglTreeInfo(rhs)
  , fixEntry_(NULL)
  , toZero_(NULL)
  , toOne_(NULL)
  , integerVariable_(NULL)
  , backward_(NULL)
  , fixingEntry_(NULL)
  , numberVariables_(rhs.numberVariables_)
  , numberIntegers_(rhs.numberIntegers_)
  , maximumEntries_(rhs.maximumEntries_)
  , numberEntries_(rhs.numberEntries_)
{
  if (numberVariables_) {
    fixEntry_ = new CliqueEntry[maximumEntries_];
    memcpy(fixEntry_, rhs.fixEntry_, maximumEntries_ * sizeof(CliqueEntry));
    if (numberEntries_ < 0) {
      // in order
      toZero_ = CoinCopyOfArray(rhs.toZero_, numberIntegers_ + 1);
      toOne_ = CoinCopyOfArray(rhs.toOne_, numberIntegers_);
    } else {
      // not in order
      fixingEntry_ = CoinCopyOfArray(rhs.fixingEntry_, maximumEntries_);
    }
    integerVariable_ = CoinCopyOfArray(rhs.integerVariable_, numberIntegers_);
    backward_ = CoinCopyOfArray(rhs.backward_, numberVariables_);
  }
}
// Clone
CglTreeInfo *
CglTreeProbingInfo::clone() const
{
  return new CglTreeProbingInfo(*this);
}

// Assignment operator
CglTreeProbingInfo &
CglTreeProbingInfo::operator=(const CglTreeProbingInfo &rhs)
{
  if (this != &rhs) {
    CglTreeInfo::operator=(rhs);
    delete[] fixEntry_;
    delete[] toZero_;
    delete[] toOne_;
    delete[] integerVariable_;
    delete[] backward_;
    delete[] fixingEntry_;
    numberVariables_ = rhs.numberVariables_;
    numberIntegers_ = rhs.numberIntegers_;
    maximumEntries_ = rhs.maximumEntries_;
    numberEntries_ = rhs.numberEntries_;
    if (numberVariables_) {
      fixEntry_ = new CliqueEntry[maximumEntries_];
      memcpy(fixEntry_, rhs.fixEntry_, maximumEntries_ * sizeof(CliqueEntry));
      if (numberEntries_ < 0) {
        // in order
        toZero_ = CoinCopyOfArray(rhs.toZero_, numberIntegers_ + 1);
        toOne_ = CoinCopyOfArray(rhs.toOne_, numberIntegers_);
        fixingEntry_ = NULL;
      } else {
        // not in order
        fixingEntry_ = CoinCopyOfArray(rhs.fixingEntry_, maximumEntries_);
        toZero_ = NULL;
        toOne_ = NULL;
      }
      toZero_ = CoinCopyOfArray(rhs.toZero_, numberIntegers_ + 1);
      toOne_ = CoinCopyOfArray(rhs.toOne_, numberIntegers_);
      integerVariable_ = CoinCopyOfArray(rhs.integerVariable_, numberIntegers_);
      backward_ = CoinCopyOfArray(rhs.backward_, numberVariables_);
    } else {
      fixEntry_ = NULL;
      toZero_ = NULL;
      toOne_ = NULL;
      integerVariable_ = NULL;
      backward_ = NULL;
      fixingEntry_ = NULL;
    }
  }
  return *this;
}

// Destructor
CglTreeProbingInfo::~CglTreeProbingInfo()
{
  delete[] fixEntry_;
  delete[] toZero_;
  delete[] toOne_;
  delete[] integerVariable_;
  delete[] backward_;
  delete[] fixingEntry_;
}
static int outDupsEtc(int numberIntegers, int &numberCliques, int &numberMatrixCliques,
  CoinBigIndex *&cliqueStart, char *&cliqueType, CliqueEntry *&entry,
  int numberLastTime, int printit)
{
  bool allNew = false;
  int *whichP = new int[numberIntegers];
  int iClique;
  assert(sizeof(int) == 4);
  assert(sizeof(CliqueEntry) == 4);
  // If lots then get rid of short ones
#define KEEP_CLIQUES 10000
  if (numberCliques - numberMatrixCliques > KEEP_CLIQUES) {
    int *sort = new int[numberCliques];
    for (iClique = numberMatrixCliques; iClique < numberCliques; iClique++) {
      CoinBigIndex j = cliqueStart[iClique];
      int n = static_cast< int >(cliqueStart[iClique + 1] - j);
      sort[iClique] = n;
    }
    std::sort(sort + numberMatrixCliques, sort + numberCliques);
    int allow = sort[numberCliques - KEEP_CLIQUES];
    int nEqual = 0;
    for (iClique = numberCliques - KEEP_CLIQUES; iClique < numberCliques; iClique++) {
      if (sort[iClique] > allow)
        break;
      else
        nEqual++;
    }
    delete[] sort;
    CoinBigIndex j = cliqueStart[numberMatrixCliques];
    CoinBigIndex put = j;
    int nClique = numberMatrixCliques;
    for (iClique = numberMatrixCliques; iClique < numberCliques; iClique++) {
      CoinBigIndex end = cliqueStart[iClique + 1];
      int n = static_cast< int >(end - j);
      bool copy = false;
      if (n > allow) {
        copy = true;
      } else if (n == allow && nEqual) {
        copy = true;
        nEqual--;
      }
      if (copy) {
        cliqueType[nClique++] = cliqueType[iClique];
        for (; j < end; j++)
          entry[put++] = entry[j];
      }
      j = cliqueStart[iClique + 1];
      cliqueStart[nClique] = put;
    }
    numberCliques = nClique;
  }
  // sort
  for (iClique = 0; iClique < numberCliques; iClique++) {
    CoinBigIndex j = cliqueStart[iClique];
    int n = static_cast< int >(cliqueStart[iClique + 1] - j);
    for (int i = 0; i < n; i++)
      whichP[i] = sequenceInCliqueEntry(entry[i + j]);
    CoinSort_2(whichP, whichP + n, (reinterpret_cast< int * >(entry)) + j);
  }
  // lexicographic sort
  int *which = new int[numberCliques];
  CoinBigIndex *position = new CoinBigIndex[numberCliques];
  int *sort = new int[numberCliques];
  int *value = new int[numberCliques];
  for (iClique = 0; iClique < numberCliques; iClique++) {
    which[iClique] = iClique;
    sort[iClique] = sequenceInCliqueEntry(entry[cliqueStart[iClique]]);
    value[iClique] = sort[iClique];
    position[iClique] = 0;
  }
  CoinSort_2(sort, sort + numberCliques, which);
  int lastDone = -1;
  int nDup = 0;
  CoinBigIndex nSave = 0;
  while (lastDone < numberCliques - 1) {
    int jClique = lastDone + 1;
    int jFirst = jClique;
    int iFirst = which[jFirst];
    int iValue = value[iFirst];
    CoinBigIndex iPos = position[iFirst];
    jClique++;
    for (; jClique < numberCliques; jClique++) {
      int kClique = which[jClique];
      int jValue = value[kClique];
      if (jValue > iValue || position[kClique] < iPos)
        break;
    }
    if (jClique == jFirst + 1) {
      // done that bit
      lastDone++;
    } else {
      // use next bit to sort and then repeat
      int jLast = jClique;
      for (jClique = jFirst; jClique < jLast; jClique++) {
        int kClique = which[jClique];
        int iValue = value[kClique];
        // put at end if finished
        if (iValue < numberIntegers) {
          CoinBigIndex kPos = position[kClique] + 1;
          position[kClique] = kPos;
          kPos += cliqueStart[kClique];
          if (kPos == cliqueStart[kClique + 1]) {
            iValue = numberIntegers;
          } else {
            iValue = sequenceInCliqueEntry(entry[kPos]);
          }
          value[kClique] = iValue;
        }
        sort[jClique] = iValue;
      }
      CoinSort_2(sort + jFirst, sort + jLast, which + jFirst);
      // if duplicate mark and move on
      int iLowest = numberCliques;
      for (jClique = jFirst; jClique < jLast; jClique++) {
        int kClique = which[jClique];
        int iValue = value[kClique];
        if (iValue < numberIntegers)
          break;
        iLowest = std::min(iLowest, kClique);
      }
      if (jClique > jFirst) {
        // mark all apart from lowest number as duplicate and move on
        lastDone = jClique - 1;
        for (jClique = jFirst; jClique <= lastDone; jClique++) {
          int kClique = which[jClique];
          if (kClique != iLowest) {
            value[kClique] = -2;
            nDup++;
            nSave += cliqueStart[kClique + 1] - cliqueStart[kClique];
          }
        }
      }
    }
  }
  if (printit)
    printf("%d duplicates\n", nDup);
  // Now see if any subset
  int nOut = 0;
  for (int jClique = 0; jClique < numberCliques; jClique++) {
    if (value[jClique] != -2) {
      position[jClique] = cliqueStart[jClique];
      value[jClique] = sequenceInCliqueEntry(entry[cliqueStart[jClique]]);
    }
  }
  nSave = 0;
  int startLooking = 0;
  for (int jClique = 0; jClique < numberCliques; jClique++) {
    int kClique = which[jClique];
    if (value[kClique] == -2) {
      nOut++;
      nSave += cliqueStart[kClique + 1] - cliqueStart[kClique];
      if (jClique == startLooking)
        startLooking++;
      continue;
    }
    int kValue = value[kClique];
    for (int iiClique = startLooking; iiClique < jClique; iiClique++) {
      int iClique = which[iiClique];
      int iValue = value[iClique];
      if (iValue == -2 || iValue == numberIntegers) {
        if (iiClique == startLooking)
          startLooking++;
        continue;
      } else {
        if (kValue > static_cast< int >(sequenceInCliqueEntry(entry[cliqueStart[iClique + 1] - 1]))) {
          value[iClique] = numberIntegers;
          continue;
        }
      }
      if (iValue < kValue) {
        while (iValue < kValue) {
          CoinBigIndex iPos = position[iClique] + 1;
          position[iClique] = iPos;
          if (iPos == cliqueStart[iClique + 1]) {
            iValue = numberIntegers;
          } else {
            iValue = sequenceInCliqueEntry(entry[iPos]);
          }
          value[iClique] = iValue;
        }
      }
      if (iValue > kValue)
        continue; // not a candidate
      // See if subset (remember duplicates have gone)
      if (cliqueStart[iClique + 1] - position[iClique] > cliqueStart[kClique + 1] - cliqueStart[kClique]) {
        // could be subset ?
        CoinBigIndex offset = cliqueStart[iClique] - position[kClique];
        CoinBigIndex j;
        bool subset = true;
        // what about different fixes bool odd=false;
        for (j = cliqueStart[kClique] + 1; j < cliqueStart[kClique + 1]; j++) {
          int kColumn = sequenceInCliqueEntry(entry[j]);
          int iColumn = sequenceInCliqueEntry(entry[j + offset]);
          if (iColumn > kColumn) {
            subset = false;
          } else {
            while (iColumn < kColumn) {
              offset++;
              if (j + offset < cliqueStart[iClique + 1]) {
                iColumn = sequenceInCliqueEntry(entry[j + offset]);
              } else {
                subset = false;
                break;
              }
            }
          }
          if (!subset)
            break;
        }
        if (subset) {
          value[kClique] = -2;
          if (printit > 1)
            printf("clique %d is subset of %d\n", kClique, iClique);
          nOut++;
          break;
        }
      }
    }
  }
  if (nOut) {
    if (printit)
      printf("Can get rid of %d cliques\n", nOut);
    // make new copy
    int nNewC = numberCliques - nOut;
    int size = static_cast< int >(cliqueStart[numberCliques] - nSave);
    int n = 0;
    CoinBigIndex *start = new CoinBigIndex[nNewC + 1];
    char *type = new char[nNewC];
    start[0] = 0;
    CliqueEntry *entryC = new CliqueEntry[size];
    CoinBigIndex nel = 0;
    allNew = true;
    for (int jClique = 0; jClique < numberCliques; jClique++) {
      int kClique = which[jClique];
      if (value[kClique] != -2 && kClique < numberMatrixCliques) {
        if (kClique >= numberLastTime)
          allNew = false;
        int nn = static_cast< int >(cliqueStart[kClique + 1] - cliqueStart[kClique]);
        memcpy(entryC + nel, entry + cliqueStart[kClique], nn * sizeof(CliqueEntry));
        nel += nn;
        type[n++] = cliqueType[kClique];
        start[n] = nel;
      }
    }
    int nM = n;
    for (int jClique = 0; jClique < numberCliques; jClique++) {
      int kClique = which[jClique];
      if (value[kClique] != -2 && kClique >= numberMatrixCliques) {
        if (kClique >= numberLastTime)
          allNew = false;
        int nn = static_cast< int >(cliqueStart[kClique + 1] - cliqueStart[kClique]);
        memcpy(entryC + nel, entry + cliqueStart[kClique], nn * sizeof(CliqueEntry));
        nel += nn;
        type[n++] = cliqueType[kClique];
        start[n] = nel;
      }
    }
    // move
    numberCliques = n;
    numberMatrixCliques = nM;
    delete[] cliqueStart;
    cliqueStart = start;
    delete[] entry;
    entry = entryC;
    delete[] cliqueType;
    cliqueType = type;
    if (printit > 1) {
      for (int jClique = 0; jClique < numberCliques; jClique++) {
        printf("%d [ ", jClique);
        for (CoinBigIndex i = cliqueStart[jClique]; i < cliqueStart[jClique + 1]; i++)
          printf("%d(%d) ", sequenceInCliqueEntry(entry[i]), oneFixesInCliqueEntry(entry[i]));
        printf("]\n");
      }
    }
    if (printit)
      printf("%d matrix cliques and %d found by probing\n", numberMatrixCliques, numberCliques - numberMatrixCliques);
  }
  delete[] value;
  delete[] sort;
  delete[] which;
  delete[] position;
  delete[] whichP;
  if (!allNew)
    return nOut;
  else
    return -1;
}
OsiSolverInterface *
CglTreeProbingInfo::analyze(const OsiSolverInterface &si, int createSolver,
  int numberExtraCliques, const CoinBigIndex *starts,
  const CliqueEntry *entries, const char *type)
{
  if (!createSolver)
    return NULL;
  convert();
  if (!numberIntegers_)
    return NULL;
  bool alwaysDo = false;
  if (numberExtraCliques < 0) {
    alwaysDo = true;
    numberExtraCliques = 0;
  }
  bool printit = false;
  int numberCliques = 0;
  CoinBigIndex numberEntries = 0;
  CoinBigIndex *cliqueStart = NULL;
  CliqueEntry *entry = NULL;
  char *cliqueType = NULL;
  int *whichP = new int[numberIntegers_];
  int *whichM = new int[numberIntegers_];
  int *whichClique = NULL;
  int numberRows = si.getNumRows();
  int numberMatrixCliques = 0;
  const CoinPackedMatrix *rowCopy = si.getMatrixByRow();
  assert(numberRows && si.getNumCols());
  int iRow;
  const int *column = rowCopy->getIndices();
  const double *elementByRow = rowCopy->getElements();
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
  const int *rowLength = rowCopy->getVectorLengths();
  const double *lower = si.getColLower();
  const double *upper = si.getColUpper();
  const double *rowLower = si.getRowLower();
  const double *rowUpper = si.getRowUpper();
  for (int iPass = 0; iPass < 2; iPass++) {
    if (iPass) {
      CoinBigIndex numberExtraEntries = 0;
      if (numberExtraCliques)
        numberExtraEntries = starts[numberExtraCliques];
      cliqueStart = new CoinBigIndex[numberCliques + 1 + numberExtraCliques];
      cliqueStart[0] = 0;
      entry = new CliqueEntry[numberEntries + numberExtraEntries];
      cliqueType = new char[numberCliques + numberExtraCliques];
      whichClique = new int[numberEntries + numberExtraEntries];
      numberCliques = 0;
      numberEntries = 0;
    }
#if 1
    for (iRow = 0; iRow < numberRows; iRow++) {
      int numberP1 = 0, numberM1 = 0;
      //int numberTotal = 0;
      CoinBigIndex j;
      double upperValue = rowUpper[iRow];
      double lowerValue = rowLower[iRow];
      bool good = true;
      for (j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
        int iColumn = column[j];
        double value = elementByRow[j];
        if (upper[iColumn] - lower[iColumn] < 1.0e-8) {
          // fixed
          upperValue -= lower[iColumn] * value;
          lowerValue -= lower[iColumn] * value;
          continue;
        } else if (backward_[iColumn] < 0) {
          good = false;
          break;
        } else {
          iColumn = backward_[iColumn];
          //numberTotal++;
        }
        if (fabs(value) != 1.0) {
          good = false;
        } else if (value > 0.0) {
          assert(numberP1 < numberIntegers_);
          whichP[numberP1++] = iColumn;
          ;
        } else {
          assert(numberM1 < numberIntegers_);
          whichM[numberM1++] = iColumn;
        }
      }
      int iUpper = upperValue > INT_MAX ? INT_MAX : static_cast<int> (floor(upperValue+1.0e-5));
      int iLower = lowerValue < INT_MIN ? INT_MIN : static_cast<int> (ceil(lowerValue-1.0e-5));
      int state = 0;
      if (upperValue < 1.0e6) {
        if (iUpper == 1 - numberM1)
          state = 1;
        else if (iUpper == -numberM1)
          state = 2;
        else if (iUpper < -numberM1)
          state = 3;
        if (fabs(static_cast< double >(iUpper) - upperValue) > 1.0e-9)
          state = -1;
      }
      if (!state && lowerValue > -1.0e6) {
        if (-iLower == 1 - numberP1)
          state = -1;
        else if (-iLower == -numberP1)
          state = -2;
        else if (-iLower < -numberP1)
          state = -3;
        if (fabs(static_cast< double >(iLower) - lowerValue) > 1.0e-9)
          state = -1;
      }
      if (numberP1 + numberM1 < 2)
        state = -1;
      if (good && state > 0) {
        if (abs(state) == 3) {
          // infeasible
          printf("FFF Infeasible\n");
          ;
          //feasible=false;
          break;
        } else if (abs(state) == 2) {
          // we can fix all
          //numberFixed += numberP1+numberM1;
          printf("FFF can fix %d\n", numberP1 + numberM1);
        } else {
          for (j = 0; j < numberP1; j++) {
            int iColumn = whichP[j];
            if (iPass) {
              CliqueEntry temp;
              setOneFixesInCliqueEntry(temp, true);
              setSequenceInCliqueEntry(temp, iColumn);
              entry[numberEntries] = temp;
            }
            numberEntries++;
          }
          for (j = 0; j < numberM1; j++) {
            int iColumn = whichM[j];
            if (iPass) {
              CliqueEntry temp;
              setOneFixesInCliqueEntry(temp, false);
              setSequenceInCliqueEntry(temp, iColumn);
              entry[numberEntries] = temp;
            }
            numberEntries++;
          }
          if (iPass) {
            if (iLower != iUpper) {
              // slack
              cliqueType[numberCliques] = 'S';
            } else {
              cliqueType[numberCliques] = 'E';
            }
            cliqueStart[numberCliques + 1] = numberEntries;
          }
          numberCliques++;
        }
      }
    }
#endif
    if (numberExtraCliques) {
      CoinBigIndex numberExtraEntries = starts[numberExtraCliques];
      memcpy(entry + numberEntries, entries, numberExtraEntries * sizeof(CliqueEntry));
      for (int iClique = 0; iClique < numberExtraCliques; iClique++) {
        cliqueType[numberCliques] = type[iClique];
        numberCliques++;
        cliqueStart[numberCliques] = starts[iClique] + numberEntries;
      }
      numberEntries += numberExtraEntries;
    }
    numberMatrixCliques = numberCliques;
    //int size = toZero_[numberIntegers_];
    //char * used = new char [size];
    //memset(used,0,size);
    // find two cliques
    int nFix = 0;
    for (int iColumn = 0; iColumn < static_cast< int >(numberIntegers_); iColumn++) {
      int j;
      for (j = toZero_[iColumn]; j < toOne_[iColumn]; j++) {
        int jColumn = sequenceInCliqueEntry(fixEntry_[j]);
        // just look at ones beore (this also skips non 0-1)
        if (jColumn < iColumn) {
          int k;
          for (k = toZero_[jColumn]; k < toOne_[jColumn]; k++) {
            if (sequenceInCliqueEntry(fixEntry_[k]) == (iColumn)) {
              if (oneFixesInCliqueEntry(fixEntry_[j])) {
                if (oneFixesInCliqueEntry(fixEntry_[k])) {
                  if (printit && !iPass)
                    printf("%d to zero implies %d to one and %d to zero implies %d to one\n",
                      iColumn, jColumn, jColumn, iColumn);
                  //0-0 illegal
                  if (iPass) {
                    CliqueEntry temp;
                    setOneFixesInCliqueEntry(temp, false);
                    setSequenceInCliqueEntry(temp, iColumn);
                    entry[numberEntries] = temp;
                  }
                  numberEntries++;
                  if (iPass) {
                    CliqueEntry temp;
                    setOneFixesInCliqueEntry(temp, false);
                    setSequenceInCliqueEntry(temp, jColumn);
                    entry[numberEntries] = temp;
                  }
                  numberEntries++;
                  if (iPass) {
                    // slack
                    cliqueType[numberCliques] = 'S';
                    cliqueStart[numberCliques + 1] = numberEntries;
                  }
                  numberCliques++;
                } else {
                  if (printit && !iPass)
                    printf("%d to zero implies %d to one and %d to zero implies %d to zero\n",
                      iColumn, jColumn, jColumn, iColumn);
                  nFix++; // jColumn is 1
                }
              } else {
                if (oneFixesInCliqueEntry(fixEntry_[k])) {
                  if (printit && !iPass)
                    printf("%d to zero implies %d to zero and %d to zero implies %d to one\n",
                      iColumn, jColumn, jColumn, iColumn);
                  nFix++; // iColumn is 1
                } else {
                  if (printit && !iPass)
                    printf("%d to zero implies %d to zero and %d to zero implies %d to zero\n",
                      iColumn, jColumn, jColumn, iColumn);
                  nFix++; // jColumn=iColumn
                }
              }
            }
          }
          for (k = toOne_[jColumn]; k < toZero_[jColumn + 1]; k++) {
            if (sequenceInCliqueEntry(fixEntry_[k]) == (iColumn)) {
              if (oneFixesInCliqueEntry(fixEntry_[j])) {
                if (oneFixesInCliqueEntry(fixEntry_[k])) {
                  if (printit && !iPass)
                    printf("%d to zero implies %d to one and %d to one implies %d to one\n",
                      iColumn, jColumn, jColumn, iColumn);
                  nFix++; //iColumn is 1
                } else {
                  if (printit && !iPass)
                    printf("%d to zero implies %d to one and %d to one implies %d to zero\n",
                      iColumn, jColumn, jColumn, iColumn);
                  nFix++; // iColumn+jcolumn=1
                }
              } else {
                if (oneFixesInCliqueEntry(fixEntry_[k])) {
                  if (printit && !iPass)
                    printf("%d to zero implies %d to zero and %d to one implies %d to one\n",
                      iColumn, jColumn, jColumn, iColumn);
                  // 0-1 illegal
                  if (iPass) {
                    CliqueEntry temp;
                    setOneFixesInCliqueEntry(temp, false);
                    setSequenceInCliqueEntry(temp, iColumn);
                    entry[numberEntries] = temp;
                  }
                  numberEntries++;
                  if (iPass) {
                    CliqueEntry temp;
                    setOneFixesInCliqueEntry(temp, true);
                    setSequenceInCliqueEntry(temp, jColumn);
                    entry[numberEntries] = temp;
                  }
                  numberEntries++;
                  if (iPass) {
                    // slack
                    cliqueType[numberCliques] = 'S';
                    cliqueStart[numberCliques + 1] = numberEntries;
                  }
                  numberCliques++;
                } else {
                  if (printit && !iPass)
                    printf("%d to zero implies %d to zero and %d to one implies %d to zero\n",
                      iColumn, jColumn, jColumn, iColumn);
                  nFix++; // jColumn is 0
                }
              }
            }
          }
        }
      }
      for (j = toOne_[iColumn]; j < toZero_[iColumn + 1]; j++) {
        int jColumn = sequenceInCliqueEntry(fixEntry_[j]);
        if (jColumn < iColumn) {
          int k;
          for (k = toZero_[jColumn]; k < toOne_[jColumn]; k++) {
            if (sequenceInCliqueEntry(fixEntry_[k]) == (iColumn)) {
              if (oneFixesInCliqueEntry(fixEntry_[j])) {
                if (oneFixesInCliqueEntry(fixEntry_[k])) {
                  if (printit && !iPass)
                    printf("%d to one implies %d to one and %d to zero implies %d to one\n",
                      iColumn, jColumn, jColumn, iColumn);
                  nFix++; // jColumn is 1
                } else {
                  if (printit && !iPass)
                    printf("%d to one implies %d to one and %d to zero implies %d to zero\n",
                      iColumn, jColumn, jColumn, iColumn);
                  // 1-0 illegal
                  if (iPass) {
                    CliqueEntry temp;
                    setOneFixesInCliqueEntry(temp, true);
                    setSequenceInCliqueEntry(temp, iColumn);
                    entry[numberEntries] = temp;
                  }
                  numberEntries++;
                  if (iPass) {
                    CliqueEntry temp;
                    setOneFixesInCliqueEntry(temp, false);
                    setSequenceInCliqueEntry(temp, jColumn);
                    entry[numberEntries] = temp;
                  }
                  numberEntries++;
                  if (iPass) {
                    // slack
                    cliqueType[numberCliques] = 'S';
                    cliqueStart[numberCliques + 1] = numberEntries;
                  }
                  numberCliques++;
                }
              } else {
                if (oneFixesInCliqueEntry(fixEntry_[k])) {
                  if (printit && !iPass)
                    printf("%d to one implies %d to zero and %d to zero implies %d to one\n",
                      iColumn, jColumn, jColumn, iColumn);
                  nFix++; // iColumn+jColumn=1
                } else {
                  if (printit && !iPass)
                    printf("%d to one implies %d to zero and %d to zero implies %d to zero\n",
                      iColumn, jColumn, jColumn, iColumn);
                  nFix++; // iColumn is 0
                }
              }
            }
          }
          for (k = toOne_[jColumn]; k < toZero_[jColumn + 1]; k++) {
            if (sequenceInCliqueEntry(fixEntry_[k]) == (iColumn)) {
              if (oneFixesInCliqueEntry(fixEntry_[j])) {
                if (oneFixesInCliqueEntry(fixEntry_[k])) {
                  if (printit && !iPass)
                    printf("%d to one implies %d to one and %d to one implies %d to one\n",
                      iColumn, jColumn, jColumn, iColumn);
                  nFix++; // iColumn == jColumn
                } else {
                  if (printit && !iPass)
                    printf("%d to one implies %d to one and %d to one implies %d to zero\n",
                      iColumn, jColumn, jColumn, iColumn);
                  nFix++; // iColumn is 0
                }
              } else {
                if (oneFixesInCliqueEntry(fixEntry_[k])) {
                  if (printit && !iPass)
                    printf("%d to one implies %d to zero and %d to one implies %d to one\n",
                      iColumn, jColumn, jColumn, iColumn);
                  nFix++; // jColumn is 0
                } else {
                  if (printit && !iPass)
                    printf("%d to one implies %d to zero and %d to one implies %d to zero\n",
                      iColumn, jColumn, jColumn, iColumn);
                  // 1-1 illegal
                  if (iPass) {
                    CliqueEntry temp;
                    setOneFixesInCliqueEntry(temp, true);
                    setSequenceInCliqueEntry(temp, iColumn);
                    entry[numberEntries] = temp;
                  }
                  numberEntries++;
                  if (iPass) {
                    CliqueEntry temp;
                    setOneFixesInCliqueEntry(temp, true);
                    setSequenceInCliqueEntry(temp, jColumn);
                    entry[numberEntries] = temp;
                  }
                  numberEntries++;
                  if (iPass) {
                    // slack
                    cliqueType[numberCliques] = 'S';
                    cliqueStart[numberCliques + 1] = numberEntries;
                  }
                  numberCliques++;
                }
              }
            }
          }
        }
      }
    }
    if (!iPass)
      printf("%d cliques and %d fixed (%d already from matrix))\n",
        numberCliques - numberMatrixCliques, nFix, numberMatrixCliques);
  }
  int iClique;
  outDupsEtc(numberIntegers_, numberCliques, numberMatrixCliques,
    cliqueStart, cliqueType, entry, -1, printit ? 2 : 1);
  printf("%d matrix cliques and %d found by probing\n", numberMatrixCliques, numberCliques - numberMatrixCliques);
  int *zeroStart = new int[numberIntegers_ + 1];
  int *oneStart = new int[numberIntegers_];
  int *zeroCount = new int[numberIntegers_];
  int *oneCount = new int[numberIntegers_];
  char *mark = new char[numberIntegers_];
  memset(mark, 0, numberIntegers_);
  int nStrengthen = -1;
  int iPass = 0;
  while (nStrengthen && iPass < 20) {
    iPass++;
    int numberLastTime = numberCliques;
    int *count = new int[numberCliques];
    int i, iColumn;
    for (i = 0; i < numberCliques; i++) {
      count[i] = 0;
    }
    int *whichCount = new int[numberCliques];
    CoinZeroN(zeroCount, numberIntegers_);
    CoinZeroN(oneCount, numberIntegers_);
    for (iClique = 0; iClique < numberCliques; iClique++) {
      for (CoinBigIndex j = cliqueStart[iClique]; j < cliqueStart[iClique + 1]; j++) {
        int iColumn = static_cast< int >(sequenceInCliqueEntry(entry[j]));
        if (oneFixesInCliqueEntry(entry[j])) {
          oneCount[iColumn]++;
        } else {
          zeroCount[iColumn]++;
        }
      }
    }
    int j;
    zeroStart[0] = 0;
    cliqueStart[0] = 0;
    for (j = 0; j < numberIntegers_; j++) {
      int n;
      n = zeroCount[j];
      zeroCount[j] = 0;
      oneStart[j] = zeroStart[j] + n;
      n = oneCount[j];
      oneCount[j] = 0;
      zeroStart[j + 1] = oneStart[j] + n;
    }
    for (iClique = 0; iClique < numberCliques; iClique++) {
      for (CoinBigIndex j = cliqueStart[iClique]; j < cliqueStart[iClique + 1]; j++) {
        int iColumn = static_cast< int >(sequenceInCliqueEntry(entry[j]));
        if (oneFixesInCliqueEntry(entry[j])) {
          int k = oneCount[iColumn];
          oneCount[iColumn]++;
          int put = oneStart[iColumn] + k;
          whichClique[put] = iClique;
        } else {
          int k = zeroCount[iColumn];
          zeroCount[iColumn]++;
          int put = zeroStart[iColumn] + k;
          whichClique[put] = iClique;
        }
      }
    }
    nStrengthen = 0;
    CoinBigIndex numberEntries = cliqueStart[numberCliques];
    CoinBigIndex maximumEntries = numberEntries;
    int maximumCliques = numberCliques;
    for (iColumn = 0; iColumn < numberIntegers_; iColumn++) {
      int i;
      int n;
      int nCount = 0;
      n = 0;
      for (i = zeroStart[iColumn]; i < oneStart[iColumn]; i++) {
        int jClique = whichClique[i];
        //if (jClique<numberMatrixCliques)
        //continue;
        CoinBigIndex j = cliqueStart[jClique];
        //assert (cliqueStart[jClique+1]==j+2);
        for (; j < cliqueStart[jClique + 1]; j++) {
          CliqueEntry eJ = entry[j];
          int jColumn = sequenceInCliqueEntry(eJ);
          if (jColumn > iColumn && !mark[jColumn]) {
            mark[jColumn] = 1;
            whichP[n++] = jColumn;
            assert(n < numberIntegers_);
            if (oneFixesInCliqueEntry(eJ)) {
              for (int k = oneStart[jColumn]; k < zeroStart[jColumn + 1]; k++) {
                int jClique = whichClique[k];
                if (!count[jClique])
                  whichCount[nCount++] = jClique;
                count[jClique]++;
              }
            } else {
              for (int k = zeroStart[jColumn]; k < oneStart[jColumn]; k++) {
                int jClique = whichClique[k];
                if (!count[jClique])
                  whichCount[nCount++] = jClique;
                count[jClique]++;
              }
            }
          }
        }
      }
      std::sort(whichP, whichP + n);
      for (i = 0; i < nCount; i++) {
        int jClique = whichCount[i];
        int jCount = count[jClique];
        count[jClique] = 0;
        if (jCount == cliqueStart[jClique + 1] - cliqueStart[jClique]) {
          printf("Zero can extend %d [ ", jClique);
          for (CoinBigIndex i = cliqueStart[jClique]; i < cliqueStart[jClique + 1]; i++)
            printf("%d(%d) ", sequenceInCliqueEntry(entry[i]), oneFixesInCliqueEntry(entry[i]));
          printf("] by %d(0)\n", iColumn);
          nStrengthen++;
          if (numberEntries + jCount + 1 > maximumEntries) {
            maximumEntries = std::max(numberEntries + jCount + 1, (maximumEntries * 12) / 10 + 100);
            CliqueEntry *temp = new CliqueEntry[maximumEntries];
            memcpy(temp, entry, numberEntries * sizeof(CliqueEntry));
            delete[] entry;
            entry = temp;
            int *tempI = new int[maximumEntries];
            memcpy(tempI, whichClique, numberEntries * sizeof(int));
            delete[] whichClique;
            whichClique = tempI;
          }
          if (numberCliques == maximumCliques) {
            maximumCliques = (maximumCliques * 12) / 10 + 100;
            CoinBigIndex *temp = new CoinBigIndex[maximumCliques + 1];
            memcpy(temp, cliqueStart, (numberCliques + 1) * sizeof(CoinBigIndex));
            delete[] cliqueStart;
            cliqueStart = temp;
            char *tempT = new char[maximumCliques];
            memcpy(tempT, cliqueType, numberCliques);
            delete[] cliqueType;
            cliqueType = tempT;
          }
          CliqueEntry eI;
          eI.fixes = 0;
          setSequenceInCliqueEntry(eI, iColumn);
          setOneFixesInCliqueEntry(eI, false);
          entry[numberEntries++] = eI;
          whichM[0] = iColumn;
          int n = 1;
          for (CoinBigIndex i = cliqueStart[jClique]; i < cliqueStart[jClique + 1]; i++) {
            entry[numberEntries++] = entry[i];
            whichM[n++] = sequenceInCliqueEntry(entry[i]);
          }
          CoinSort_2(whichM, whichM + n, (reinterpret_cast< int * >(entry)) + numberEntries - n);
          cliqueType[numberCliques] = 'S';
          numberCliques++;
          cliqueStart[numberCliques] = numberEntries;
        }
      }
      for (i = 0; i < n; i++)
        mark[whichP[i]] = 0;
      nCount = 0;
      n = 0;
      for (i = oneStart[iColumn]; i < zeroStart[iColumn + 1]; i++) {
        int jClique = whichClique[i];
        //if (jClique<numberMatrixCliques)
        //continue;
        CoinBigIndex j = cliqueStart[jClique];
        //assert (cliqueStart[jClique+1]==j+2);
        for (; j < cliqueStart[jClique + 1]; j++) {
          CliqueEntry eJ = entry[j];
          int jColumn = sequenceInCliqueEntry(eJ);
          if (jColumn > iColumn && !mark[jColumn]) {
            mark[jColumn] = 1;
            whichP[n++] = jColumn;
            assert(n < numberIntegers_);
            if (oneFixesInCliqueEntry(eJ)) {
              for (int k = oneStart[jColumn]; k < zeroStart[jColumn + 1]; k++) {
                int jClique = whichClique[k];
                if (!count[jClique])
                  whichCount[nCount++] = jClique;
                count[jClique]++;
              }
            } else {
              for (int k = zeroStart[jColumn]; k < oneStart[jColumn]; k++) {
                int jClique = whichClique[k];
                if (!count[jClique])
                  whichCount[nCount++] = jClique;
                count[jClique]++;
              }
            }
          }
        }
      }
      std::sort(whichP, whichP + n);
      for (i = 0; i < nCount; i++) {
        int jClique = whichCount[i];
        int jCount = count[jClique];
        count[jClique] = 0;
        if (jCount == cliqueStart[jClique + 1] - cliqueStart[jClique]) {
#if 1
          if (printit) {
            printf("One can extend %d [ ", jClique);
            for (CoinBigIndex i = cliqueStart[jClique]; i < cliqueStart[jClique + 1]; i++)
              printf("%d(%d) ", sequenceInCliqueEntry(entry[i]), oneFixesInCliqueEntry(entry[i]));
            printf("] by %d(1)\n", iColumn);
          }
#endif
          nStrengthen++;
          if (numberEntries + jCount + 1 > maximumEntries) {
            maximumEntries = std::max(numberEntries + jCount + 1, (maximumEntries * 12) / 10 + 100);
            CliqueEntry *temp = new CliqueEntry[maximumEntries];
            memcpy(temp, entry, numberEntries * sizeof(CliqueEntry));
            delete[] entry;
            entry = temp;
            int *tempI = new int[maximumEntries];
            memcpy(tempI, whichClique, numberEntries * sizeof(int));
            delete[] whichClique;
            whichClique = tempI;
          }
          if (numberCliques == maximumCliques) {
            maximumCliques = (maximumCliques * 12) / 10 + 100;
            CoinBigIndex *temp = new CoinBigIndex[maximumCliques + 1];
            memcpy(temp, cliqueStart, (numberCliques + 1) * sizeof(CoinBigIndex));
            delete[] cliqueStart;
            cliqueStart = temp;
            char *tempT = new char[maximumCliques];
            memcpy(tempT, cliqueType, numberCliques);
            delete[] cliqueType;
            cliqueType = tempT;
          }
          CliqueEntry eI;
          eI.fixes = 0;
          setSequenceInCliqueEntry(eI, iColumn);
          setOneFixesInCliqueEntry(eI, true);
          entry[numberEntries++] = eI;
          whichM[0] = iColumn;
          int n = 1;
          for (CoinBigIndex i = cliqueStart[jClique]; i < cliqueStart[jClique + 1]; i++) {
            entry[numberEntries++] = entry[i];
            whichM[n++] = sequenceInCliqueEntry(entry[i]);
          }
          CoinSort_2(whichM, whichM + n, (reinterpret_cast< int * >(entry)) + numberEntries - n);
          cliqueType[numberCliques] = 'S';
          numberCliques++;
          cliqueStart[numberCliques] = numberEntries;
        }
      }
      for (i = 0; i < n; i++)
        mark[whichP[i]] = 0;
    }
    if (nStrengthen) {
      int numberDeleted = outDupsEtc(numberIntegers_, numberCliques, numberMatrixCliques,
        cliqueStart, cliqueType, entry, numberLastTime, printit ? 2 : 1);
      if (numberDeleted < 0 || (iPass > 1 && numberCliques - numberDeleted > 5000))
        nStrengthen = 0;
    }
    delete[] count;
    delete[] whichCount;
  }
#if 0
  if (numberCliques>numberMatrixCliques) {
    // should keep as cliques and also use in branching ??
    double * element = new double [numberIntegers_];
    for (iClique=numberMatrixCliques;iClique<numberCliques;iClique++) {
      int n=0;
      double rhs=1.0;
      for (int i=cliqueStart[iClique];i<cliqueStart[iClique+1];i++) {
	CliqueEntry eI=entry[i];
	int iColumn = integerVariable_[sequenceInCliqueEntry(eI)];
	whichP[n]=iColumn;
	if (oneFixesInCliqueEntry(eI)) {
	  element[n++]=1.0;
	} else {
	  element[n++]=-1.0;
	  rhs -= 1.0;
	}
      }
      stored->addCut(-COIN_DBL_MAX,rhs,n,whichP,element);
    } 
    delete [] element;
  }
#endif
  OsiSolverInterface *newSolver = NULL;
  if (numberCliques > numberMatrixCliques || alwaysDo) {
    newSolver = si.clone();
    // Delete all rows
    CoinBigIndex *start = new CoinBigIndex[std::max(numberRows, numberCliques + 1)];
    int i;
    int *start2 = reinterpret_cast< int * >(start);
    for (i = 0; i < numberRows; i++)
      start2[i] = i;
    newSolver->deleteRows(numberRows, start2);
    start[0] = 0;
    CoinBigIndex numberElements = cliqueStart[numberCliques];
    int *column = new int[numberElements];
    double *element = new double[numberElements];
    double *lower = new double[numberCliques];
    double *upper = new double[numberCliques];
    numberElements = 0;
    for (iClique = 0; iClique < numberCliques; iClique++) {
      double rhs = 1.0;
      for (CoinBigIndex i = cliqueStart[iClique]; i < cliqueStart[iClique + 1]; i++) {
        CliqueEntry eI = entry[i];
        int iColumn = integerVariable_[sequenceInCliqueEntry(eI)];
        column[numberElements] = iColumn;
        if (oneFixesInCliqueEntry(eI)) {
          element[numberElements++] = 1.0;
        } else {
          element[numberElements++] = -1.0;
          rhs -= 1.0;
        }
      }
      start[iClique + 1] = numberElements;
      assert(cliqueType[iClique] == 'S' || cliqueType[iClique] == 'E');
      if (cliqueType[iClique] == 'S')
        lower[iClique] = -COIN_DBL_MAX;
      else
        lower[iClique] = rhs;
      upper[iClique] = rhs;
    }
    newSolver->addRows(numberCliques, start, column, element, lower, upper);
    delete[] start;
    delete[] column;
    delete[] element;
    delete[] lower;
    delete[] upper;
  }
  delete[] mark;
  delete[] whichP;
  delete[] whichM;
  delete[] cliqueStart;
  delete[] entry;
  delete[] cliqueType;
  delete[] zeroStart;
  delete[] oneStart;
  delete[] zeroCount;
  delete[] oneCount;
  delete[] whichClique;
  return newSolver;
}
// Take action if cut generator can fix a variable (toValue -1 for down, +1 for up)
bool CglTreeProbingInfo::fixes(int variable, int toValue, int fixedVariable, bool fixedToLower)
{
  //printf("%d going to %d fixes %d at %g\n",variable,toValue,fixedVariable,fixedToValue);
  int intVariable = backward_[variable];
  if (intVariable < 0) // off as no longer in order FIX
    return true; // not 0-1 (well wasn't when constructor was called)
  int intFix = backward_[fixedVariable];
  if (intFix < 0)
    intFix = numberIntegers_ + fixedVariable; // not 0-1
  int fixedTo = fixedToLower ? 0 : 1;
  if (numberEntries_ == maximumEntries_) {
    // See if taking too much memory
    if (maximumEntries_ >= std::max(1000000, 10 * numberIntegers_))
      return false;
    maximumEntries_ += 100 + maximumEntries_ / 2;
    CliqueEntry *temp1 = new CliqueEntry[maximumEntries_];
    if (fixEntry_)
      memcpy(temp1, fixEntry_, numberEntries_ * sizeof(CliqueEntry));
    delete[] fixEntry_;
    fixEntry_ = temp1;
    int *temp2 = new int[maximumEntries_];
    if (fixingEntry_)
      memcpy(temp2, fixingEntry_, numberEntries_ * sizeof(int));
    delete[] fixingEntry_;
    fixingEntry_ = temp2;
  }
  CliqueEntry entry1;
  entry1.fixes = 0;
  setOneFixesInCliqueEntry(entry1, fixedTo != 0);
  setSequenceInCliqueEntry(entry1, intFix);
  fixEntry_[numberEntries_] = entry1;
  assert(toValue == -1 || toValue == 1);
  assert(fixedTo == 0 || fixedTo == 1);
  if (toValue < 0)
    fixingEntry_[numberEntries_++] = intVariable << 1;
  else
    fixingEntry_[numberEntries_++] = (intVariable << 1) | 1;
  return true;
}
// Initalizes fixing arrays etc - returns true if we want to save info
int CglTreeProbingInfo::initializeFixing(const OsiSolverInterface *model)
{
  if (numberEntries_ >= 0)
    return 2; // already got arrays
  else if (numberEntries_ == -2)
    return numberEntries_;
  delete[] fixEntry_;
  delete[] toZero_;
  delete[] toOne_;
  delete[] integerVariable_;
  delete[] backward_;
  delete[] fixingEntry_;
  numberVariables_ = model->getNumCols();
  // Too many ... but
  integerVariable_ = new int[numberVariables_];
  backward_ = new int[numberVariables_];
  numberIntegers_ = 0;
  int i;
  // Get integer types
  const char *columnType = model->getColType(true);
  for (i = 0; i < numberVariables_; i++) {
    backward_[i] = -1;
    if (columnType[i]) {
      if (columnType[i] == 1) {
        backward_[i] = numberIntegers_;
        integerVariable_[numberIntegers_++] = i;
      } else {
        backward_[i] = -2;
      }
    }
  }
  toZero_ = NULL;
  toOne_ = NULL;
  fixEntry_ = NULL;
  fixingEntry_ = NULL;
  maximumEntries_ = 0;
  numberEntries_ = 0;
  return 1;
}
// Converts to ordered and takes out duplicates
void CglTreeProbingInfo::convert()
{
  if (numberEntries_ >= 0) {
    CoinSort_2(fixingEntry_, fixingEntry_ + numberEntries_, fixEntry_);
    assert(!toZero_);
    toZero_ = new int[numberIntegers_ + 1];
    toOne_ = new int[numberIntegers_];
    toZero_[0] = 0;
    int n = 0;
    int put = 0;
    for (int intVariable = 0; intVariable < numberIntegers_; intVariable++) {
      int last = n;
      for (; n < numberEntries_; n++) {
        int value = fixingEntry_[n];
        int iVar = value >> 1;
        int way = value & 1;
        if (intVariable != iVar || way)
          break;
      }
      if (n > last) {
        // sort
        assert(sizeof(int) == 4);
        std::sort(reinterpret_cast< unsigned int * >(fixEntry_) + last,
          reinterpret_cast< unsigned int * >(fixEntry_) + n);
        CliqueEntry temp2;
        temp2.fixes = 0;
        setSequenceInCliqueEntry(temp2, numberVariables_ + 1);
        for (int i = last; i < n; i++) {
          if (sequenceInCliqueEntry(temp2) != sequenceInCliqueEntry(fixEntry_[i]) || oneFixesInCliqueEntry(temp2) || oneFixesInCliqueEntry(fixEntry_[i])) {
            temp2 = fixEntry_[i];
            fixEntry_[put++] = temp2;
          }
        }
      }
      toOne_[intVariable] = put;
      last = n;
      for (; n < numberEntries_; n++) {
        int value = fixingEntry_[n];
        int iVar = value >> 1;
        if (intVariable != iVar)
          break;
      }
      if (n > last) {
        // sort
        assert(sizeof(int) == 4);
        std::sort(reinterpret_cast< unsigned int * >(fixEntry_) + last,
          reinterpret_cast< unsigned int * >(fixEntry_) + n);
        CliqueEntry temp2;
        temp2.fixes = 0;
        setSequenceInCliqueEntry(temp2, numberVariables_ + 1);
        for (int i = last; i < n; i++) {
          if (sequenceInCliqueEntry(temp2) != sequenceInCliqueEntry(fixEntry_[i]) || oneFixesInCliqueEntry(temp2) || oneFixesInCliqueEntry(fixEntry_[i])) {
            temp2 = fixEntry_[i];
            fixEntry_[put++] = temp2;
          }
        }
        last = n;
      }
      toZero_[intVariable + 1] = put;
    }
    delete[] fixingEntry_;
    fixingEntry_ = NULL;
    numberEntries_ = -2;
  }
}
// Fix entries in a solver using implications
int CglTreeProbingInfo::fixColumns(OsiSolverInterface &si) const
{
  int nFix = 0;
  const double *lower = si.getColLower();
  const double *upper = si.getColUpper();
  bool feasible = true;
  for (int jColumn = 0; jColumn < static_cast< int >(numberIntegers_); jColumn++) {
    int iColumn = integerVariable_[jColumn];
    if (upper[iColumn] == 0.0) {
      int j;
      for (j = toZero_[jColumn]; j < toOne_[jColumn]; j++) {
        int kColumn = sequenceInCliqueEntry(fixEntry_[j]);
        kColumn = integerVariable_[kColumn];
        bool fixToOne = oneFixesInCliqueEntry(fixEntry_[j]);
        if (fixToOne) {
          if (lower[kColumn] == 0.0) {
            if (upper[kColumn] == 1.0) {
              si.setColLower(kColumn, 1.0);
              nFix++;
            } else {
              // infeasible!
              feasible = false;
            }
          }
        } else {
          if (upper[kColumn] == 1.0) {
            if (lower[kColumn] == 0.0) {
              si.setColUpper(kColumn, 0.0);
              nFix++;
            } else {
              // infeasible!
              feasible = false;
            }
          }
        }
      }
    } else if (lower[iColumn] == 1.0) {
      int j;
      for (j = toOne_[jColumn]; j < toZero_[jColumn + 1]; j++) {
        int kColumn = sequenceInCliqueEntry(fixEntry_[j]);
        kColumn = integerVariable_[kColumn];
        bool fixToOne = oneFixesInCliqueEntry(fixEntry_[j]);
        if (fixToOne) {
          if (lower[kColumn] == 0.0) {
            if (upper[kColumn] == 1.0) {
              si.setColLower(kColumn, 1.0);
              nFix++;
            } else {
              // infeasible!
              feasible = false;
            }
          }
        } else {
          if (upper[kColumn] == 1.0) {
            if (lower[kColumn] == 0.0) {
              si.setColUpper(kColumn, 0.0);
              nFix++;
            } else {
              // infeasible!
              feasible = false;
            }
          }
        }
      }
    }
  }
  if (!feasible) {
#ifdef COIN_DEVELOP
    printf("treeprobing says infeasible!\n");
#endif
    nFix = -1;
  }
  return nFix;
}
// Fix entries in a solver using implications for one variable
int CglTreeProbingInfo::fixColumns(int iColumn, int value, OsiSolverInterface &si) const
{
  assert(value == 0 || value == 1);
  int nFix = 0;
  const double *lower = si.getColLower();
  const double *upper = si.getColUpper();
  bool feasible = true;
  int jColumn = backward_[iColumn];
  if (jColumn < 0 || !toZero_)
    return 0;
  if (!value) {
    int j;
    for (j = toZero_[jColumn]; j < toOne_[jColumn]; j++) {
      int kColumn = sequenceInCliqueEntry(fixEntry_[j]);
      kColumn = integerVariable_[kColumn];
      bool fixToOne = oneFixesInCliqueEntry(fixEntry_[j]);
      if (fixToOne) {
        if (lower[kColumn] == 0.0) {
          if (upper[kColumn] == 1.0) {
            si.setColLower(kColumn, 1.0);
            nFix++;
          } else {
            // infeasible!
            feasible = false;
          }
        }
      } else {
        if (upper[kColumn] == 1.0) {
          if (lower[kColumn] == 0.0) {
            si.setColUpper(kColumn, 0.0);
            nFix++;
          } else {
            // infeasible!
            feasible = false;
          }
        }
      }
    }
  } else {
    int j;
    for (j = toOne_[jColumn]; j < toZero_[jColumn + 1]; j++) {
      int kColumn = sequenceInCliqueEntry(fixEntry_[j]);
      kColumn = integerVariable_[kColumn];
      bool fixToOne = oneFixesInCliqueEntry(fixEntry_[j]);
      if (fixToOne) {
        if (lower[kColumn] == 0.0) {
          if (upper[kColumn] == 1.0) {
            si.setColLower(kColumn, 1.0);
            nFix++;
          } else {
            // infeasible!
            feasible = false;
          }
        }
      } else {
        if (upper[kColumn] == 1.0) {
          if (lower[kColumn] == 0.0) {
            si.setColUpper(kColumn, 0.0);
            nFix++;
          } else {
            // infeasible!
            feasible = false;
          }
        }
      }
    }
  }
  if (!feasible) {
#ifdef COIN_DEVELOP
    printf("treeprobing says infeasible!\n");
#endif
    nFix = -1;
  }
  return nFix;
}
// Packs down entries
int CglTreeProbingInfo::packDown()
{
  convert();
  int iPut = 0;
  int iLast = 0;
  for (int jColumn = 0; jColumn < static_cast< int >(numberIntegers_); jColumn++) {
    int j;
    for (j = iLast; j < toOne_[jColumn]; j++) {
      int kColumn = sequenceInCliqueEntry(fixEntry_[j]);
      if (kColumn < numberIntegers_)
        fixEntry_[iPut++] = fixEntry_[j];
    }
    iLast = toOne_[jColumn];
    toOne_[jColumn] = iPut;
    for (j = iLast; j < toZero_[jColumn + 1]; j++) {
      int kColumn = sequenceInCliqueEntry(fixEntry_[j]);
      if (kColumn < numberIntegers_)
        fixEntry_[iPut++] = fixEntry_[j];
    }
    iLast = toZero_[jColumn + 1];
    toZero_[jColumn + 1] = iPut;
  }
  return iPut;
}
void CglTreeProbingInfo::generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
  const CglTreeInfo /*info*/) const
{
  const double *lower = si.getColLower();
  const double *upper = si.getColUpper();
  const double *colsol = si.getColSolution();
  CoinPackedVector lbs(false);
  CoinPackedVector ubs(false);
  int numberFixed = 0;
  char *fixed = NULL;
  for (int jColumn = 0; jColumn < static_cast< int >(numberIntegers_); jColumn++) {
    int iColumn = integerVariable_[jColumn];
    assert(iColumn >= 0 && iColumn < si.getNumCols());
    if (lower[iColumn] == 0.0 && upper[iColumn] == 1.0) {
      double value1 = colsol[iColumn];
      int j;
      for (j = toZero_[jColumn]; j < toOne_[jColumn]; j++) {
        int kColumn = sequenceInCliqueEntry(fixEntry_[j]);
        kColumn = integerVariable_[kColumn];
        assert(kColumn >= 0 && kColumn < si.getNumCols());
        assert(kColumn != iColumn);
        if (lower[kColumn] == 0.0 && upper[kColumn] == 1.0) {
          double value2 = colsol[kColumn];
          bool fixToOne = oneFixesInCliqueEntry(fixEntry_[j]);
          if (fixToOne) {
            if (value1 + value2 < 0.99999) {
              OsiRowCut rc;
              int index[2];
              double element[2];
              index[0] = iColumn;
              element[0] = 1.0;
              index[1] = kColumn;
              element[1] = 1.0;
              rc.setLb(1.0);
              rc.setUb(COIN_DBL_MAX);
              rc.setRow(2, index, element, false);
              cs.insertIfNotDuplicate(rc);
            }
          } else {
            if (value1 - value2 < -0.00001) {
              OsiRowCut rc;
              int index[2];
              double element[2];
              index[0] = iColumn;
              element[0] = 1.0;
              index[1] = kColumn;
              element[1] = -1.0;
              rc.setLb(0.0);
              rc.setUb(COIN_DBL_MAX);
              rc.setRow(2, index, element, false);
              cs.insertIfNotDuplicate(rc);
            }
          }
        }
      }
      for (j = toOne_[jColumn]; j < toZero_[jColumn + 1]; j++) {
        int kColumn = sequenceInCliqueEntry(fixEntry_[j]);
        kColumn = integerVariable_[kColumn];
        assert(kColumn >= 0 && kColumn < si.getNumCols());
        assert(kColumn != iColumn);
        if (lower[kColumn] == 0.0 && upper[kColumn] == 1.0) {
          double value2 = colsol[kColumn];
          bool fixToOne = oneFixesInCliqueEntry(fixEntry_[j]);
          if (fixToOne) {
            if (value1 - value2 > 0.00001) {
              OsiRowCut rc;
              int index[2];
              double element[2];
              index[0] = iColumn;
              element[0] = 1.0;
              index[1] = kColumn;
              element[1] = -1.0;
              rc.setLb(-COIN_DBL_MAX);
              rc.setUb(0.0);
              rc.setRow(2, index, element, false);
              cs.insertIfNotDuplicate(rc);
            }
          } else {
            if (value1 + value2 > 1.00001) {
              OsiRowCut rc;
              int index[2];
              double element[2];
              index[0] = iColumn;
              element[0] = 1.0;
              index[1] = kColumn;
              element[1] = 1.0;
              rc.setLb(-COIN_DBL_MAX);
              rc.setUb(1.0);
              rc.setRow(2, index, element, false);
              cs.insertIfNotDuplicate(rc);
            }
          }
        }
      }
    } else if (upper[iColumn] == 0.0) {
      for (int j = toZero_[jColumn]; j < toOne_[jColumn]; j++) {
        int kColumn01 = sequenceInCliqueEntry(fixEntry_[j]);
        int kColumn = integerVariable_[kColumn01];
        assert(kColumn >= 0 && kColumn < si.getNumCols());
        bool fixToOne = oneFixesInCliqueEntry(fixEntry_[j]);
        if (lower[kColumn] == 0.0 && upper[kColumn] == 1.0) {
          if (!fixed) {
            fixed = new char[numberIntegers_];
            memset(fixed, 0, numberIntegers_);
          }
          if (fixToOne) {
            if ((fixed[kColumn01] & 1) == 0) {
              fixed[kColumn01] |= 1;
              lbs.insert(kColumn, 1.0);
            }
          } else {
            if ((fixed[kColumn01] & 2) == 0) {
              fixed[kColumn01] |= 2;
              ubs.insert(kColumn, 0.0);
            }
          }
          numberFixed++;
        } else if ((fixToOne && upper[kColumn] == 0.0) || (!fixToOne && lower[kColumn] == 1.0)) {
          // infeasible!
          // generate infeasible cut and return
          OsiRowCut rc;
          rc.setLb(COIN_DBL_MAX);
          rc.setUb(0.0);
          cs.insertIfNotDuplicate(rc);
          //printf("IMPINFEAS!\n");
          return;
        }
      }
    } else {
      for (int j = toOne_[jColumn]; j < toZero_[jColumn + 1]; j++) {
        int kColumn01 = sequenceInCliqueEntry(fixEntry_[j]);
        int kColumn = integerVariable_[kColumn01];
        assert(kColumn >= 0 && kColumn < si.getNumCols());
        bool fixToOne = oneFixesInCliqueEntry(fixEntry_[j]);
        if (lower[kColumn] == 0.0 && upper[kColumn] == 1.0) {
          if (!fixed) {
            fixed = new char[numberIntegers_];
            memset(fixed, 0, numberIntegers_);
          }
          if (fixToOne) {
            if ((fixed[kColumn01] & 1) == 0) {
              fixed[kColumn01] |= 1;
              lbs.insert(kColumn, 1.0);
            }
          } else {
            if ((fixed[kColumn01] & 2) == 0) {
              fixed[kColumn01] |= 2;
              ubs.insert(kColumn, 0.0);
            }
          }
          numberFixed++;
        } else if ((fixToOne && upper[kColumn] == 0.0) || (!fixToOne && lower[kColumn] == 1.0)) {
          // infeasible!
          // generate infeasible cut and return
          OsiRowCut rc;
          rc.setLb(COIN_DBL_MAX);
          rc.setUb(0.0);
          cs.insertIfNotDuplicate(rc);
          //printf("IMPINFEAS!\n");
          return;
        }
      }
    }
  }
  if (numberFixed) {
    int feasible = true;
    for (int i = 0; i < numberIntegers_; i++) {
      if (fixed[i] == 3) {
        feasible = false;
        break;
      }
    }
    delete[] fixed;
    fixed = NULL;
    if (feasible) {
      //printf("IMP fixed %d\n",numberFixed);
      OsiColCut cc;
      cc.setUbs(ubs);
      cc.setLbs(lbs);
      cc.setEffectiveness(1.0);
      cs.insert(cc);
    } else {
      // infeasible!
      // generate infeasible cut
      OsiRowCut rc;
      rc.setLb(COIN_DBL_MAX);
      rc.setUb(0.0);
      cs.insertIfNotDuplicate(rc);
      //printf("IMPINFEAS!\n");
    }
  }

  if (fixed)
    delete[] fixed;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
