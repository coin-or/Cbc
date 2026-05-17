// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>
//#define CGL_DEBUG 2
#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "OsiRowCutDebugger.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CglStored.hpp"
#include "CglTreeInfo.hpp"
#include "CoinFinite.hpp"
//-------------------------------------------------------------------
// Generate Stored cuts
//-------------------------------------------------------------------
void CglStored::generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
  const CglTreeInfo /*info*/)
{
  // Get basic problem information
  const double *solution = si.getColSolution();
  int numberRowCuts = cuts_.sizeRowCuts();
  for (int i = 0; i < numberRowCuts; i++) {
    const OsiRowCut *rowCutPointer = cuts_.rowCutPtr(i);
    double violation = rowCutPointer->violated(solution);
    if (violation >= requiredViolation_)
      cs.insert(*rowCutPointer);
  }
  if (probingInfo_) {
    int number01 = probingInfo_->numberIntegers();
    const CliqueEntry *entry = probingInfo_->fixEntries();
    const int *toZero = probingInfo_->toZero();
    const int *toOne = probingInfo_->toOne();
    const int *integerVariable = probingInfo_->integerVariable();
    const double *lower = si.getColLower();
    const double *upper = si.getColUpper();
    OsiRowCut cut;
    int column[2];
    double element[2];
    for (int i = 0; i < number01; i++) {
      int iColumn = integerVariable[i];
      if (upper[iColumn] == lower[iColumn])
        continue;
      double value1 = solution[iColumn];
      for (int j = toZero[i]; j < toOne[i]; j++) {
        int jColumn = sequenceInCliqueEntry(entry[j]);
        if (jColumn < number01) {
          jColumn = integerVariable[jColumn];
          assert(jColumn >= 0);
          double value2 = solution[jColumn];
          if (oneFixesInCliqueEntry(entry[j])) {
            double violation = 1.0 - value1 - value2;
            if (violation > requiredViolation_) {
              //printf("XXX can do %d + %d >=1\n",iColumn,jColumn);
              cut.setLb(1.0);
              cut.setUb(COIN_DBL_MAX);
              column[0] = iColumn;
              element[0] = 1.0;
              column[1] = jColumn;
              element[1] = 1.0;
              cut.setEffectiveness(violation);
              cut.setRow(2, column, element, false);
              cs.insertIfNotDuplicate(cut);
            }
          } else {
            double violation = value2 - value1;
            if (violation > requiredViolation_) {
              //printf("XXX can do %d >= %d\n",iColumn,jColumn);
              cut.setLb(0.0);
              cut.setUb(COIN_DBL_MAX);
              column[0] = iColumn;
              element[0] = 1.0;
              column[1] = jColumn;
              element[1] = -1.0;
              cut.setEffectiveness(violation);
              cut.setRow(2, column, element, false);
              cs.insertIfNotDuplicate(cut);
            }
          }
        } else {
          jColumn -= number01; // not 0-1
          double value2 = solution[jColumn];
          double lowerValue = lower[jColumn];
          double upperValue = upper[jColumn];
          if (oneFixesInCliqueEntry(entry[j])) {
            double violation = upperValue - value1 * (upperValue - lowerValue) - value2;
            if (violation > requiredViolation_) {
              //printf("XXX can do %g*%d + %d >=%g\n",(upperValue-lowerValue),iColumn,jColumn,upperValue);
              cut.setLb(upperValue);
              cut.setUb(COIN_DBL_MAX);
              column[0] = iColumn;
              element[0] = upperValue - lowerValue;
              column[1] = jColumn;
              element[1] = 1.0;
              cut.setEffectiveness(violation);
              cut.setRow(2, column, element, false);
              cs.insertIfNotDuplicate(cut);
            }
          } else {
            double violation = value2 - value1 * (upperValue - lowerValue) - lowerValue;
            if (violation > requiredViolation_) {
              //printf("XXX can do %g*%d >= %d -%g\n",(upperValue-lowerValue),iColumn,jColumn,lowerValue);
              cut.setLb(-lowerValue);
              cut.setUb(COIN_DBL_MAX);
              column[0] = iColumn;
              element[0] = upperValue - lowerValue;
              column[1] = jColumn;
              element[1] = -1.0;
              cut.setEffectiveness(violation);
              cut.setRow(2, column, element, false);
              cs.insertIfNotDuplicate(cut);
            }
          }
        }
      }
      for (int j = toOne[i]; j < toZero[i + 1]; j++) {
        int jColumn = sequenceInCliqueEntry(entry[j]);
        if (jColumn < number01) {
          jColumn = integerVariable[jColumn];
          assert(jColumn >= 0);
          double value2 = solution[jColumn];
          if (oneFixesInCliqueEntry(entry[j])) {
            double violation = value1 - value2;
            if (violation > requiredViolation_) {
              //printf("XXX can do %d <= %d\n",iColumn,jColumn);
              cut.setLb(-COIN_DBL_MAX);
              cut.setUb(0.0);
              column[0] = iColumn;
              element[0] = 1.0;
              column[1] = jColumn;
              element[1] = -1.0;
              cut.setEffectiveness(violation);
              cut.setRow(2, column, element, false);
              cs.insertIfNotDuplicate(cut);
            }
          } else {
            double violation = value1 + value2 - 1.0;
            if (violation > requiredViolation_) {
              //printf("XXX can do %d + %d <=1\n",iColumn,jColumn);
              cut.setLb(-COIN_DBL_MAX);
              cut.setUb(1.0);
              column[0] = iColumn;
              element[0] = 1.0;
              column[1] = jColumn;
              element[1] = 1.0;
              cut.setEffectiveness(violation);
              cut.setRow(2, column, element, false);
              cs.insertIfNotDuplicate(cut);
            }
          }
        } else {
          jColumn -= number01; // not 0-1
          double value2 = solution[jColumn];
          double lowerValue = lower[jColumn];
          double upperValue = upper[jColumn];
          if (oneFixesInCliqueEntry(entry[j])) {
            double violation = lowerValue + (upperValue - lowerValue) * value1 - value2;
            if (violation > requiredViolation_) {
              //printf("XXX can do %g*%d <= %d -%g\n",(upperValue-lowerValue),iColumn,jColumn,lowerValue);
              cut.setLb(-COIN_DBL_MAX);
              cut.setUb(-lowerValue);
              column[0] = iColumn;
              element[0] = upperValue - lowerValue;
              column[1] = jColumn;
              element[1] = -1.0;
              cut.setEffectiveness(violation);
              cut.setRow(2, column, element, false);
              cs.insertIfNotDuplicate(cut);
            }
          } else {
            double violation = (upperValue - lowerValue) * value1 + value2 - upperValue;
            if (violation > requiredViolation_) {
              //printf("XXX can do %g*%d + %d <=%g\n",(upperValue-lowerValue),iColumn,jColumn,upperValue);
              cut.setLb(-COIN_DBL_MAX);
              cut.setUb(upperValue);
              column[0] = iColumn;
              element[0] = upperValue - lowerValue;
              column[1] = jColumn;
              element[1] = 1.0;
              cut.setEffectiveness(violation);
              cut.setRow(2, column, element, false);
              cs.insertIfNotDuplicate(cut);
            }
          }
        }
      }
    }
  }
}
// Add cuts
void CglStored::addCut(const OsiCuts &cs)
{
  int numberRowCuts = cs.sizeRowCuts();
  for (int i = 0; i < numberRowCuts; i++) {
    int numberCuts = cuts_.sizeRowCuts();
    cuts_.insert(*cs.rowCutPtr(i));
    // Keep when diving
    cuts_.rowCutPtr(numberCuts)->setEffectiveness(COIN_DBL_MAX);
  }
}
// Add a row cut
void CglStored::addCut(const OsiRowCut &cut)
{
  int numberCuts = cuts_.sizeRowCuts();
  cuts_.insert(cut);
  // Keep when diving
  cuts_.rowCutPtr(numberCuts)->setEffectiveness(COIN_DBL_MAX);
}
// Add a row cut from a packed vector
void CglStored::addCut(double lb, double ub, const CoinPackedVector &vector)
{
  OsiRowCut rc;
  rc.setRow(vector);
  rc.mutableRow().setTestForDuplicateIndex(false);
  rc.setLb(lb);
  rc.setUb(ub);
  int numberCuts = cuts_.sizeRowCuts();
  cuts_.insert(rc);
  // Keep when diving
  cuts_.rowCutPtr(numberCuts)->setEffectiveness(COIN_DBL_MAX);
}
// Add a row cut from elements
void CglStored::addCut(double lb, double ub, int size, const int *colIndices, const double *elements)
{
  OsiRowCut rc;
  rc.setRow(size, colIndices, elements, false);
  rc.setLb(lb);
  rc.setUb(ub);
  int numberCuts = cuts_.sizeRowCuts();
  cuts_.insert(rc);
  // Keep when diving
  cuts_.rowCutPtr(numberCuts)->setEffectiveness(COIN_DBL_MAX);
}

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
CglStored::CglStored(int numberColumns)
  : CglCutGenerator()
  , requiredViolation_(1.0e-5)
  , probingInfo_(NULL)
  , numberColumns_(numberColumns)
  , bestSolution_(NULL)
  , bounds_(NULL)
{
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CglStored::CglStored(const CglStored &source)
  : CglCutGenerator(source)
  , requiredViolation_(source.requiredViolation_)
  , probingInfo_(NULL)
  , cuts_(source.cuts_)
  , numberColumns_(source.numberColumns_)
  , bestSolution_(NULL)
  , bounds_(NULL)
{
  if (source.probingInfo_)
    probingInfo_ = new CglTreeProbingInfo(*source.probingInfo_);
  if (numberColumns_) {
    bestSolution_ = CoinCopyOfArray(source.bestSolution_, numberColumns_ + 1);
    bounds_ = CoinCopyOfArray(source.bounds_, 2 * numberColumns_);
  }
}
//-------------------------------------------------------------------
// Constructor from file
//-------------------------------------------------------------------
CglStored::CglStored(const char *fileName)
  : CglCutGenerator()
  , requiredViolation_(1.0e-5)
  , probingInfo_(NULL)
  , numberColumns_(0)
  , bestSolution_(NULL)
  , bounds_(NULL)
{
  FILE *fp = fopen(fileName, "rb");
  if (fp) {
#ifndef NDEBUG
    size_t numberRead;
#endif
    int maxInCut = 0;
    int *index = NULL;
    double *coefficient = NULL;
    double rhs[2];
    int n = 0;
    while (n >= 0) {
#ifndef NDEBUG
      numberRead = fread(&n, sizeof(int), 1, fp);
      assert(numberRead == 1);
#else
      fread(&n, sizeof(int), 1, fp);
#endif
      if (n < 0)
        break;
      if (n > maxInCut) {
        maxInCut = n;
        delete[] index;
        delete[] coefficient;
        index = new int[maxInCut];
        coefficient = new double[maxInCut];
      }
#ifndef NDEBUG
      numberRead = fread(rhs, sizeof(double), 2, fp);
      assert(numberRead == 2);
#else
      fread(rhs, sizeof(double), 2, fp);
#endif
      fread(index, sizeof(int), n, fp);
      fread(coefficient, sizeof(double), n, fp);
      OsiRowCut rc;
      rc.setRow(n, index, coefficient, false);
      rc.setLb(rhs[0]);
      rc.setUb(rhs[1]);
      cuts_.insert(rc);
    }
    delete[] index;
    delete[] coefficient;
    fclose(fp);
  }
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglStored::clone() const
{
  return new CglStored(*this);
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
CglStored::~CglStored()
{
  delete probingInfo_;
  delete[] bestSolution_;
  delete[] bounds_;
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CglStored &
CglStored::operator=(const CglStored &rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    requiredViolation_ = rhs.requiredViolation_;
    cuts_ = rhs.cuts_;
    delete probingInfo_;
    if (rhs.probingInfo_)
      probingInfo_ = new CglTreeProbingInfo(*rhs.probingInfo_);
    else
      probingInfo_ = NULL;
    delete[] bestSolution_;
    delete[] bounds_;
    bestSolution_ = NULL;
    bounds_ = NULL;
    numberColumns_ = rhs.numberColumns_;
    if (numberColumns_) {
      bestSolution_ = CoinCopyOfArray(rhs.bestSolution_, numberColumns_ + 1);
      bounds_ = CoinCopyOfArray(rhs.bounds_, 2 * numberColumns_);
    }
  }
  return *this;
}
// Save stuff
void CglStored::saveStuff(double bestObjective, const double *bestSolution,
  const double *lower, const double *upper)
{
  assert(numberColumns_);
  delete[] bestSolution_;
  delete[] bounds_;
  if (bestSolution) {
    bestSolution_ = new double[numberColumns_ + 1];
    memcpy(bestSolution_, bestSolution, numberColumns_ * sizeof(double));
    bestSolution_[numberColumns_] = bestObjective;
  } else {
    bestSolution_ = NULL;
  }
  bounds_ = new double[2 * numberColumns_];
  memcpy(bounds_, lower, numberColumns_ * sizeof(double));
  memcpy(bounds_ + numberColumns_, upper, numberColumns_ * sizeof(double));
}
// Best objective
double
CglStored::bestObjective() const
{
  if (bestSolution_)
    return bestSolution_[numberColumns_];
  else
    return COIN_DBL_MAX;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
