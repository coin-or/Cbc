// $Id$
// Copyright (C) 2004, International Business Machines
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
#include "CbcBranchFollow2.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"
// Default Constructor
CbcFollowOn2::CbcFollowOn2()
  : CbcObject()
  , rhs_(NULL)
  , maximumRhs_(1)
{
}

// Useful constructor
CbcFollowOn2::CbcFollowOn2(CbcModel *model)
  : CbcObject(model)
{
  assert(model);
  OsiSolverInterface *solver = model_->solver();
  matrix_ = *solver->getMatrixByCol();
  matrix_.removeGaps();
  matrixByRow_ = *solver->getMatrixByRow();
  int numberRows = matrix_.getNumRows();
  maximumRhs_ = 1;

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
      if (floor(value) == value && value >= 1.0 && value < 100.0) {
        // check elements
        bool good = true;
        for (int j = rowStart[i]; j < rowStart[i] + rowLength[i]; j++) {
          int iColumn = column[j];
          if (!solver->isInteger(iColumn))
            good = false;
          double elValue = elementByRow[j];
          if (floor(elValue) != elValue || elValue < 1.0)
            good = false;
        }
        if (good)
          rhs_[i] = (int)value;
      }
    }
  }
}

// Copy constructor
CbcFollowOn2::CbcFollowOn2(const CbcFollowOn2 &rhs)
  : CbcObject(rhs)
  , matrix_(rhs.matrix_)
  , matrixByRow_(rhs.matrixByRow_)
  , maximumRhs_(rhs.maximumRhs_)
{
  int numberRows = matrix_.getNumRows();
  rhs_ = CoinCopyOfArray(rhs.rhs_, numberRows);
}

// Clone
CbcObject *
CbcFollowOn2::clone() const
{
  return new CbcFollowOn2(*this);
}

// Assignment operator
CbcFollowOn2 &
CbcFollowOn2::operator=(const CbcFollowOn2 &rhs)
{
  if (this != &rhs) {
    CbcObject::operator=(rhs);
    delete[] rhs_;
    matrix_ = rhs.matrix_;
    matrixByRow_ = rhs.matrixByRow_;
    int numberRows = matrix_.getNumRows();
    rhs_ = CoinCopyOfArray(rhs.rhs_, numberRows);
    maximumRhs_ = rhs.maximumRhs_;
  }
  return *this;
}

// Destructor
CbcFollowOn2::~CbcFollowOn2()
{
  delete[] rhs_;
}
/* As some computation is needed in more than one place - returns row.
   Also returns other row and effective rhs (so we can know if cut)
*/
int CbcFollowOn2::gutsOfFollowOn2(int &otherRow, int &preferredWay,
  int &effectiveRhs) const
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
      int numberOther = 0;
      int numberUnsatisfied = 0;
      for (int j = rowStart[i]; j < rowStart[i] + rowLength[i]; j++) {
        int iColumn = column[j];
        double value = elementByRow[j];
        double solValue = solution[iColumn];
        if (columnLower[iColumn] != columnUpper[iColumn]) {
          smallest = CoinMin(smallest, value);
          largest = CoinMax(largest, value);
          if (value == 1.0)
            number1++;
          else
            numberOther++;
          if (fabs(floor(solValue + 0.5) - solValue) > integerTolerance)
            numberUnsatisfied++;
        } else {
          rhsValue -= (int)(value * floor(solValue + 0.5));
        }
      }
      if (numberUnsatisfied > 1) {
        if (smallest < largest) {
#if 0
	  if  (largest>rhsValue)
	    printf("could fix\n");
	  if (number1==1&&largest==rhsValue)
	    printf("could fix\n");
#endif
          if (rhsValue <= maximumRhs_ && 0) {
            // will mean a cut but worth trying
            sort[nSort] = i;
            isort[nSort++] = 100000 - numberUnsatisfied;
          }
        } else if (largest == rhsValue) {
          sort[nSort] = i;
          isort[nSort++] = -numberUnsatisfied;
        }
      }
    }
  }
  if (nSort > 1) {
    CoinSort_2(isort, isort + nSort, sort);
    assert(isort[1] < 0);
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
      int numberUnsatisfied = 0;
      int n = 0;
      int j;
      for (j = rowStart[i]; j < rowStart[i] + rowLength[i]; j++) {
        int iColumn = column[j];
        if (columnLower[iColumn] != columnUpper[iColumn]) {
          double solValue = solution[iColumn] - columnLower[iColumn];
          if (fabs(floor(solValue + 0.5) - solValue) > integerTolerance) {
            numberUnsatisfied++;
            for (int jj = columnStart[iColumn]; jj < columnStart[iColumn] + columnLength[iColumn]; jj++) {
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
        if (fabs(floor(dvalue + 0.5) - dvalue) < integerTolerance)
          continue;
        dvalue -= floor(dvalue);
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
        //printf("Rows %d (%d) and %d (%d)\n",whichRow,rhs_[whichRow],
        //     otherRow,rhs_[otherRow]);
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

// Infeasibility - large is 0.5
double
CbcFollowOn2::infeasibility(int &preferredWay) const
{
  int otherRow = 0;
  int effectiveRhs;
  int whichRow = gutsOfFollowOn2(otherRow, preferredWay, effectiveRhs);
  if (whichRow < 0) {
    return 0.0;
  } else {
    assert(whichRow != otherRow);
    return 2.0 * model_->getDblParam(CbcModel::CbcIntegerTolerance);
  }
}

// This looks at solution and sets bounds to contain solution
void CbcFollowOn2::feasibleRegion()
{
}

// Creates a branching object
CbcBranchingObject *
CbcFollowOn2::createBranch(int way)
{
  int otherRow = 0;
  int preferredWay;
  int effectiveRhs;
  int whichRow = gutsOfFollowOn2(otherRow, preferredWay, effectiveRhs);
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
  OsiSolverInterface *solver = model_->solver();
  const double *columnLower = solver->getColLower();
  const double *columnUpper = solver->getColUpper();
  //const double * solution = solver->getColSolution();
#if 0
  //printf("Rows %d (%d) and %d (%d)\n",whichRow,rhs_[whichRow],
  //     otherRow,rhs_[otherRow]);
  int nFree=0;
  int nOut=0;
  int nImplicit=0;
  int i;
  int rhsx[100];
  double * colUpper2 = new double [numberColumns];
  memcpy(rhsx,rhs_,matrix_.getNumRows()*sizeof(int));
  for ( i=0;i<numberColumns;i++) {
    colUpper2[i]=columnUpper[i];
    if (columnLower[i]==columnUpper[i]) {
      for (int jj=columnStart[i];jj<columnStart[i]+columnLength[i];jj++) {
        int iRow = row[jj];
        nOut += (int) floor(element[jj]*solution[i]+0.5);
        rhsx[iRow] -= (int) floor(element[jj]*solution[i]+0.5);
      }
    }
  }
  int nFixedBut=0;
  for ( i=0;i<numberColumns;i++) {
    if (columnLower[i]!=columnUpper[i]) {
      nFree++;
      bool nonzero=false;
      if (fabs(solution[i])>1.0e-5) {
        nonzero=true;
        //printf("column %d value %g ",i,solution[i]);
        for (int jj=columnStart[i];jj<columnStart[i]+columnLength[i];jj++) {
          //int iRow = row[jj];
          //printf("(%d,%g) ",iRow,element[jj]);
        }
        //printf("\n");
      } 
      bool fixed=false;
      for (int jj=columnStart[i];jj<columnStart[i]+columnLength[i];jj++) {
        int iRow = row[jj];
        if (element[jj]>rhsx[iRow]) 
          fixed=true;
      }
      if (fixed) {
        nImplicit++;
        colUpper2[i]=0.0;
        if (nonzero)
          nFixedBut++;
        assert (!columnLower[i]);
      }
    }
  }
  // See if anything odd
  char * check = new char[numberColumns];
  memset(check,0,numberColumns);
  int * which2 = new int[numberColumns];
  int numberRows=matrix_.getNumRows();
  for (i=0;i<numberRows;i++) {
    if (rhsx[i]==1) {
      int nn=0;
      int j,k;
      for (j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
	int iColumn = column[j];
	double value = elementByRow[j];
	if (columnLower[iColumn]!=colUpper2[iColumn]) {
          assert (value==1.0);
          check[iColumn]=1;
          which2[nn++]=iColumn;
        }
      }
      for ( k=i+1;k<numberRows;k++) {
        if (rhsx[k]==1) {
          int nn2=0;
          int nnsame=0;
          for (int j=rowStart[k];j<rowStart[k]+rowLength[k];j++) {
            int iColumn = column[j];
            double value = elementByRow[j];
            if (columnLower[iColumn]!=colUpper2[iColumn]) {
              assert (value==1.0);
              nn2++;
              if (check[iColumn])
                nnsame++;
            }
          }
          if (nnsame==nn2) {
            if (nn2<nn)
              printf("row %d strict subset of row %d, fix some in row %d\n",
                     k,i,i);
            else if (nn2==nn)
              printf("row %d identical to row %d\n",
                     k,i);
            else if (nn2>=nn)
              abort();
          } else if (nnsame==nn&&nn2>nn) {
              printf("row %d strict superset of row %d, fix some in row %d\n",
                     k,i,k);
          }
        }
      }
      for (k=0;k<nn;k++) 
        check[which2[k]]=0;
      
    }
  }
  delete [] check;
  delete [] which2;
  delete [] colUpper2;
  printf("%d free (but %d implicitly fixed of which %d nonzero), %d out of rhs\n",nFree,nImplicit,nFixedBut,nOut);
#endif
  int nUp = 0;
  int nDown = 0;
  int *upList = new int[numberColumns];
  int *downList = new int[numberColumns];
  int j;
  for (j = rowStart[whichRow]; j < rowStart[whichRow] + rowLength[whichRow]; j++) {
    int iColumn = column[j];
    if (columnLower[iColumn] != columnUpper[iColumn]) {
      bool up = true;
      for (int jj = columnStart[iColumn]; jj < columnStart[iColumn] + columnLength[iColumn]; jj++) {
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
  CbcBranchingObject *branch
    = new CbcFixingBranchingObject(model_, way,
      nDown, downList, nUp, upList);
  delete[] upList;
  delete[] downList;
  return branch;
}
