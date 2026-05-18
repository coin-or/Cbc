// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif

#include "CoinUtilsConfig.h"

#include <cassert>
#include <cstdio>

#include "CoinFactorization.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"
#include <stdio.h>
#include <iostream>
#if COIN_FACTORIZATION_DENSE_CODE == 1
// using simple lapack interface
extern "C"
{
  /** LAPACK Fortran subroutine DGETRS. */
  void COINUTILS_LAPACK_FUNC(dgetrs,DGETRS)(char *trans, cipfint *n,
                  cipfint *nrhs, const double *A, cipfint *ldA,
                  cipfint *ipiv, double *B, cipfint *ldB, ipfint *info,
                  int trans_len);
}
#elif COIN_FACTORIZATION_DENSE_CODE == 2
// C interface
enum CBLAS_ORDER { CblasRowMajor = 101,
  CblasColMajor = 102 };
enum CBLAS_TRANSPOSE { CblasNoTrans = 111,
  CblasTrans = 112 };
extern "C" {
int clapack_dgetrs(const enum CBLAS_ORDER Order,
  const enum CBLAS_TRANSPOSE Trans,
  const int N, const int NRHS,
  const double *A, const int lda, const int *ipiv,
  double *B, const int ldb);
}
#elif COIN_FACTORIZATION_DENSE_CODE == 3
// Intel compiler
#include "mkl_lapacke.h"
#endif
// For semi-sparse
#define BITS_PER_CHECK 8
#define CHECK_SHIFT 3
typedef unsigned char CoinCheckZero;

//:class CoinFactorization.  Deals with Factorization and Updates

//  getRowSpaceIterate.  Gets space for one Row with given length
//may have to do compression  (returns true)
//also moves existing vector
bool CoinFactorization::getRowSpaceIterate(int iRow,
  int extraNeeded)
{
  const int *numberInRow = numberInRowArray_;
  int number = numberInRow[iRow];
  int *COIN_RESTRICT startRow = startRowUArray_;
  int *COIN_RESTRICT indexColumn = indexColumnUArray_;
  int *COIN_RESTRICT convertRowToColumn = convertRowToColumnUArray_;
  int space = lengthAreaU_ - startRow[maximumRowsExtra_];
  int *COIN_RESTRICT nextRow = nextRowArray_;
  int *COIN_RESTRICT lastRow = lastRowArray_;
  if (space < extraNeeded + number + 2) {
    //compression
    int iRow = nextRow[maximumRowsExtra_];
    int put = 0;
    while (iRow != maximumRowsExtra_) {
      //move
      int get = startRow[iRow];
      int getEnd = startRow[iRow] + numberInRow[iRow];

      startRow[iRow] = put;
      int i;
      for (i = get; i < getEnd; i++) {
        indexColumn[put] = indexColumn[i];
        convertRowToColumn[put] = convertRowToColumn[i];
        put++;
      }
      iRow = nextRow[iRow];
    } /* endwhile */
    numberCompressions_++;
    startRow[maximumRowsExtra_] = put;
    space = lengthAreaU_ - put;
    if (space < extraNeeded + number + 2) {
      //need more space
      //if we can allocate bigger then do so and copy
      //if not then return so code can start again
      status_ = -99;
      return false;
    }
  }
  int put = startRow[maximumRowsExtra_];
  int next = nextRow[iRow];
  int last = lastRow[iRow];

  //out
  nextRow[last] = next;
  lastRow[next] = last;
  //in at end
  last = lastRow[maximumRowsExtra_];
  nextRow[last] = iRow;
  lastRow[maximumRowsExtra_] = iRow;
  lastRow[iRow] = last;
  nextRow[iRow] = maximumRowsExtra_;
  //move
  int get = startRow[iRow];

  int *COIN_RESTRICT indexColumnU = indexColumnUArray_;
  startRow[iRow] = put;
  while (number) {
    number--;
    indexColumnU[put] = indexColumnU[get];
    convertRowToColumn[put] = convertRowToColumn[get];
    put++;
    get++;
  } /* endwhile */
  //add four for luck
  startRow[maximumRowsExtra_] = put + extraNeeded + 4;
  return true;
}

//  getColumnSpaceIterateR.  Gets space for one extra R element in Column
//may have to do compression  (returns true)
//also moves existing vector
bool CoinFactorization::getColumnSpaceIterateR(int iColumn, double value,
  int iRow)
{
  CoinFactorizationDouble *COIN_RESTRICT elementR = elementR_ + lengthAreaR_;
  int *COIN_RESTRICT indexRowR = indexRowR_ + lengthAreaR_;
  int *COIN_RESTRICT startR = startColumnRArray_ + maximumPivots_ + 1;
  int *COIN_RESTRICT numberInColumnPlus = numberInColumnPlusArray_;
  int number = numberInColumnPlus[iColumn];
  //*** modify so sees if can go in
  //see if it can go in at end
  int *COIN_RESTRICT nextColumn = nextColumnArray_;
  int *COIN_RESTRICT lastColumn = lastColumnArray_;
  if (lengthAreaR_ - startR[maximumColumnsExtra_] < number + 1) {
    //compression
    int jColumn = nextColumn[maximumColumnsExtra_];
    int put = 0;
    while (jColumn != maximumColumnsExtra_) {
      //move
      int get;
      int getEnd;

      get = startR[jColumn];
      getEnd = get + numberInColumnPlus[jColumn];
      startR[jColumn] = put;
      int i;
      for (i = get; i < getEnd; i++) {
        indexRowR[put] = indexRowR[i];
        elementR[put] = elementR[i];
        put++;
      }
      jColumn = nextColumn[jColumn];
    }
    numberCompressions_++;
    startR[maximumColumnsExtra_] = put;
  }
  // Still may not be room (as iColumn was still in)
  if (lengthAreaR_ - startR[maximumColumnsExtra_] < number + 1)
    return false;

  int next = nextColumn[iColumn];
  int last = lastColumn[iColumn];

  //out
  nextColumn[last] = next;
  lastColumn[next] = last;

  int put = startR[maximumColumnsExtra_];
  //in at end
  last = lastColumn[maximumColumnsExtra_];
  nextColumn[last] = iColumn;
  lastColumn[maximumColumnsExtra_] = iColumn;
  lastColumn[iColumn] = last;
  nextColumn[iColumn] = maximumColumnsExtra_;

  //move
  int get = startR[iColumn];
  startR[iColumn] = put;
  int i = 0;
  for (i = 0; i < number; i++) {
    elementR[put] = elementR[get];
    indexRowR[put++] = indexRowR[get++];
  }
  //insert
  elementR[put] = value;
  indexRowR[put++] = iRow;
  numberInColumnPlus[iColumn]++;
  //add 4 for luck
  startR[maximumColumnsExtra_] = std::min(static_cast< int >(put + 4), lengthAreaR_);
  return true;
}
int CoinFactorization::checkPivot(double saveFromU,
  double oldPivot) const
{
  int status;
#define ALLOW_SMALL_PIVOTS
#ifdef ALLOW_SMALL_PIVOTS
#define SMALL_PIVOT 1.0e-9
#else
#define SMALL_PIVOT 1.0e-8
#endif
  if (fabs(saveFromU) > SMALL_PIVOT) {
    double checkTolerance;

    if (numberRowsExtra_ < numberRows_ + 2) {
      checkTolerance = 1.0e-5;
    } else if (numberRowsExtra_ < numberRows_ + 10) {
      checkTolerance = 1.0e-6;
    } else if (numberRowsExtra_ < numberRows_ + 50) {
      checkTolerance = 1.0e-8;
    } else {
      checkTolerance = 1.0e-10;
    }
    checkTolerance *= relaxCheck_;
    if (fabs(1.0 - fabs(saveFromU / oldPivot)) < checkTolerance) {
      status = 0;
    } else {
#if COIN_DEBUG
      std::cout << "inaccurate pivot " << oldPivot << " "
                << saveFromU << std::endl;
#endif
      if (fabs(fabs(oldPivot) - fabs(saveFromU)) < 1.0e-12 || fabs(1.0 - fabs(saveFromU / oldPivot)) < 1.0e-8) {
        status = 1;
      } else {
        status = 2;
      }
    }
  } else {
    //error
    if (fabs(1.0 - fabs(saveFromU / oldPivot)) < 1.0e-10) {
      status = 0;
    } else {
      status = 2;
#if COIN_DEBUG
      std::cout << "inaccurate pivot " << saveFromU / oldPivot
                << " " << saveFromU << std::endl;
#endif
    }
  }
  //if (status==2)
  //printf("status %d\n",status);
  return status;
}
#ifdef CLP_FACTORIZATION_INSTRUMENT
extern double externalTimeStart;
extern double timeInFactorize;
extern double timeInUpdate;
extern double timeInFactorizeFake;
extern double timeInUpdateFake1;
extern double timeInUpdateFake2;
extern double timeInUpdateTranspose;
extern double timeInUpdateFT;
extern double timeInUpdateTwoFT;
extern double timeInReplace;
extern int numberUpdate;
extern int numberUpdateTranspose;
extern int numberUpdateFT;
extern int numberUpdateTwoFT;
extern int numberReplace;
extern int currentLengthR;
extern int currentLengthU;
extern int currentTakeoutU;
extern double averageLengthR;
extern double averageLengthL;
extern double averageLengthU;
extern double scaledLengthDense;
extern double scaledLengthDenseSquared;
extern double scaledLengthL;
extern double scaledLengthR;
extern double scaledLengthU;
extern int startLengthU;
extern int endLengthU;
#endif
//  replaceColumn.  Replaces one Column to basis
//      returns 0=OK, 1=Probably OK, 2=singular, 3=no room
//partial update already in U
int CoinFactorization::replaceColumn(CoinIndexedVector *regionSparse,
  int pivotRow,
  double pivotCheck,
  bool checkBeforeModifying,
  double)
{
#ifdef CLP_FACTORIZATION_INSTRUMENT
  double startTimeX = CoinCpuTime();
#endif
  assert(numberU_ <= numberRowsExtra_);
  int *COIN_RESTRICT startColumnU = startColumnUArray_;
  int *COIN_RESTRICT startColumn;
  int *COIN_RESTRICT indexRow;
  CoinFactorizationDouble *COIN_RESTRICT element;

  //return at once if too many iterations
  if (numberColumnsExtra_ >= maximumColumnsExtra_) {
    return 5;
  }
  if (lengthAreaU_ < startColumnU[maximumColumnsExtra_]) {
    return 3;
  }

  int *COIN_RESTRICT numberInRow = numberInRowArray_;
  int *COIN_RESTRICT numberInColumn = numberInColumnArray_;
  int *COIN_RESTRICT numberInColumnPlus = numberInColumnPlusArray_;
  int realPivotRow;
  realPivotRow = pivotColumnArray_[pivotRow];
  //zeroed out region
  double *COIN_RESTRICT region = regionSparse->denseVector();

  element = elementUArray_;
  //take out old pivot column

  // If we have done no pivots then always check before modification
  if (!numberPivots_)
    checkBeforeModifying = true;

  totalElements_ -= numberInColumn[realPivotRow];
  CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;
  CoinFactorizationDouble oldPivot = pivotRegion[realPivotRow];
  // for accuracy check
  pivotCheck = pivotCheck / oldPivot;
#if COIN_DEBUG > 1
  int checkNumber = 1000000;
  //if (numberL_) checkNumber=-1;
  if (numberR_ >= checkNumber) {
    printf("pivot row %d, check %g - alpha region:\n",
      realPivotRow, pivotCheck);
    /*int i;
      for (i=0;i<numberRows_;i++) {
      if (pivotRegion[i])
      printf("%d %g\n",i,pivotRegion[i]);
  }*/
  }
#endif
  pivotRegion[realPivotRow] = 0.0;

  int saveEnd = startColumnU[realPivotRow]
    + numberInColumn[realPivotRow];
#ifdef CLP_FACTORIZATION_INSTRUMENT
  currentTakeoutU += numberInColumn[realPivotRow];
  currentTakeoutU += numberInRow[realPivotRow];
#endif
  // not necessary at present - but take no chances for future
  numberInColumn[realPivotRow] = 0;
  //get entries in row (pivot not stored)
  int numberNonZero = 0;
  int *COIN_RESTRICT indexColumn = indexColumnUArray_;
  int *COIN_RESTRICT convertRowToColumn = convertRowToColumnUArray_;
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  int *COIN_RESTRICT startRow = startRowUArray_;
  int start = 0;
  int end = 0;
#if COIN_DEBUG > 1
  if (numberR_ >= checkNumber)
    printf("Before btranu\n");
#endif

#if COIN_ONE_ETA_COPY
  if (convertRowToColumn) {
#endif
    start = startRow[realPivotRow];
    end = start + numberInRow[realPivotRow];

    int smallestIndex = numberRowsExtra_;
    if (!checkBeforeModifying) {
      for (int i = start; i < end; i++) {
        int iColumn = indexColumn[i];
        assert(iColumn < numberRowsExtra_);
        smallestIndex = std::min(smallestIndex, iColumn);
        int j = convertRowToColumn[i];

        region[iColumn] = element[j];
#if COIN_DEBUG > 1
        if (numberR_ >= checkNumber)
          printf("%d %g\n", iColumn, region[iColumn]);
#endif
        element[j] = 0.0;
        regionIndex[numberNonZero++] = iColumn;
      }
    } else {
      for (int i = start; i < end; i++) {
        int iColumn = indexColumn[i];
        smallestIndex = std::min(smallestIndex, iColumn);
        int j = convertRowToColumn[i];

        region[iColumn] = element[j];
#if COIN_DEBUG > 1
        if (numberR_ >= checkNumber)
          printf("%d %g\n", iColumn, region[iColumn]);
#endif
        regionIndex[numberNonZero++] = iColumn;
      }
    }
    //do BTRAN - finding first one to use
    regionSparse->setNumElements(numberNonZero);
    updateColumnTransposeU(regionSparse, smallestIndex);
#if COIN_ONE_ETA_COPY
  } else {
    // use R to save where elements are
    int *COIN_RESTRICT saveWhere = NULL;
    if (checkBeforeModifying) {
      if (lengthR_ + maximumRowsExtra_ + 1 >= lengthAreaR_) {
        //not enough room
        return 3;
      }
      saveWhere = indexRowR_ + lengthR_;
    }
    replaceColumnU(regionSparse, saveWhere,
      realPivotRow);
  }
#endif
  numberNonZero = regionSparse->getNumElements();
  CoinFactorizationDouble saveFromU = 0.0;

  int startU = startColumnU[numberColumnsExtra_];
  int *COIN_RESTRICT indexU = &indexRowUArray_[startU];
  CoinFactorizationDouble *COIN_RESTRICT elementU = &elementUArray_[startU];

  // Do accuracy test here if caller is paranoid
  if (checkBeforeModifying) {
    double tolerance = zeroTolerance_;
    int number = numberInColumn[numberColumnsExtra_];

    for (int i = 0; i < number; i++) {
      int iRow = indexU[i];
      //if (numberCompressions_==99&&lengthU_==278)
      //printf("row %d saveFromU %g element %g region %g\n",
      //       iRow,saveFromU,elementU[i],region[iRow]);
      if (fabs(elementU[i]) > tolerance) {
        if (iRow != realPivotRow) {
          saveFromU -= elementU[i] * region[iRow];
        } else {
          saveFromU += elementU[i];
        }
      }
    }
    //check accuracy
    int status = checkPivot(saveFromU, pivotCheck);
    if (status) {
      // restore some things
      pivotRegion[realPivotRow] = oldPivot;
      number = saveEnd - startColumnU[realPivotRow];
      totalElements_ += number;
      numberInColumn[realPivotRow] = number;
      regionSparse->clear();
      return status;
#if COIN_ONE_ETA_COPY
    } else if (convertRowToColumn) {
#else
    } else {
#endif
      // do what we would have done by now
      for (int i = start; i < end; i++) {
        int j = convertRowToColumn[i];
        element[j] = 0.0;
      }
#if COIN_ONE_ETA_COPY
    } else {
      // delete elements
      // used R to save where elements are
      int *COIN_RESTRICT saveWhere = indexRowR_ + lengthR_;
      CoinFactorizationDouble *COIN_RESTRICT element = elementUArray_;
      int n = saveWhere[0];
      for (int i = 0; i < n; i++) {
        int where = saveWhere[i + 1];
        element[where] = 0.0;
      }
      //printf("deleting els\n");
#endif
    }
  }
  // Now zero out column of U
  //take out old pivot column
  for (int i = startColumnU[realPivotRow]; i < saveEnd; i++) {
    element[i] = 0.0;
  }
  //zero out pivot Row (before or after?)
  //add to R
  startColumn = startColumnRArray_;
  indexRow = indexRowR_;
  element = elementR_;
  int l = lengthR_;
  int number = numberR_;

  startColumn[number] = l; //for luck and first time
  number++;
  startColumn[number] = l + numberNonZero;
  numberR_ = number;
  lengthR_ = l + numberNonZero;
  totalElements_ += numberNonZero;
  if (lengthR_ >= lengthAreaR_) {
    //not enough room
    regionSparse->clear();
    return 3;
  }
#if COIN_DEBUG > 1
  if (numberR_ >= checkNumber)
    printf("After btranu\n");
#endif
  for (int i = 0; i < numberNonZero; i++) {
    int iRow = regionIndex[i];
#if COIN_DEBUG > 1
    if (numberR_ >= checkNumber)
      printf("%d %g\n", iRow, region[iRow]);
#endif

    indexRow[l] = iRow;
    element[l] = region[iRow];
    l++;
  }
  int *COIN_RESTRICT nextRow;
  int *COIN_RESTRICT lastRow;
  int next;
  int last;
#if COIN_ONE_ETA_COPY
  if (convertRowToColumn) {
#endif
    //take out row
    nextRow = nextRowArray_;
    lastRow = lastRowArray_;
    next = nextRow[realPivotRow];
    last = lastRow[realPivotRow];

    nextRow[last] = next;
    lastRow[next] = last;
    numberInRow[realPivotRow] = 0;
#if COIN_DEBUG
    nextRow[realPivotRow] = 777777;
    lastRow[realPivotRow] = 777777;
#endif
#if COIN_ONE_ETA_COPY
  }
#endif
  //do permute
  permuteArray_[numberRowsExtra_] = realPivotRow;
  // and other way
  permuteBackArray_[realPivotRow] = numberRowsExtra_;
  permuteBackArray_[numberRowsExtra_] = -1;
  ;
  //and for safety
  permuteArray_[numberRowsExtra_ + 1] = 0;

  pivotColumnArray_[pivotRow] = numberRowsExtra_;
  pivotColumnBack()[numberRowsExtra_] = pivotRow;
  startColumn = startColumnU;
  indexRow = indexRowUArray_;
  element = elementUArray_;

  numberU_++;
  number = numberInColumn[numberColumnsExtra_];

  totalElements_ += number;
  lengthU_ += number;
  if (lengthU_ >= lengthAreaU_) {
    //not enough room
    regionSparse->clear();
    return 3;
  }

  saveFromU = 0.0;

  //put in pivot
  //add row counts

#if COIN_DEBUG > 1
  if (numberR_ >= checkNumber)
    printf("On U\n");
#endif
#if COIN_ONE_ETA_COPY
  if (convertRowToColumn) {
#endif
    for (int i = 0; i < number; i++) {
      int iRow = indexU[i];
#if COIN_DEBUG > 1
      if (numberR_ >= checkNumber)
        printf("%d %g\n", iRow, elementU[i]);
#endif

      //assert ( fabs ( elementU[i] ) > zeroTolerance_ );
      if (iRow != realPivotRow) {
        int next = nextRow[iRow];
        int iNumberInRow = numberInRow[iRow];
        int space;
        int put = startRow[iRow] + iNumberInRow;

        space = startRow[next] - put;
        if (space <= 0) {
          getRowSpaceIterate(iRow, iNumberInRow + 4);
          put = startRow[iRow] + iNumberInRow;
        }
        indexColumn[put] = numberColumnsExtra_;
        convertRowToColumn[put] = i + startU;
        numberInRow[iRow] = iNumberInRow + 1;
        saveFromU = saveFromU - elementU[i] * region[iRow];
      } else {
        //zero out and save
        saveFromU += elementU[i];
        elementU[i] = 0.0;
      }
    }
    //in at end
    last = lastRow[maximumRowsExtra_];
    nextRow[last] = numberRowsExtra_;
    lastRow[maximumRowsExtra_] = numberRowsExtra_;
    lastRow[numberRowsExtra_] = last;
    nextRow[numberRowsExtra_] = maximumRowsExtra_;
    startRow[numberRowsExtra_] = startRow[maximumRowsExtra_];
    numberInRow[numberRowsExtra_] = 0;
#if COIN_ONE_ETA_COPY
  } else {
    //abort();
    for (int i = 0; i < number; i++) {
      int iRow = indexU[i];
#if COIN_DEBUG > 1
      if (numberR_ >= checkNumber)
        printf("%d %g\n", iRow, elementU[i]);
#endif

      if (fabs(elementU[i]) > tolerance) {
        if (iRow != realPivotRow) {
          saveFromU = saveFromU - elementU[i] * region[iRow];
        } else {
          //zero out and save
          saveFromU += elementU[i];
          elementU[i] = 0.0;
        }
      } else {
        elementU[i] = 0.0;
      }
    }
  }
#endif
  //column in at beginning (as empty)
  int *COIN_RESTRICT nextColumn = nextColumnArray_;
  int *COIN_RESTRICT lastColumn = lastColumnArray_;
  next = nextColumn[maximumColumnsExtra_];
  lastColumn[next] = numberColumnsExtra_;
  nextColumn[maximumColumnsExtra_] = numberColumnsExtra_;
  nextColumn[numberColumnsExtra_] = next;
  lastColumn[numberColumnsExtra_] = maximumColumnsExtra_;
  //check accuracy - but not if already checked (optimization problem)
  int status = (checkBeforeModifying) ? 0 : checkPivot(saveFromU, pivotCheck);

  if (status != 2) {

    CoinFactorizationDouble pivotValue = 1.0 / saveFromU;

    pivotRegion[numberRowsExtra_] = pivotValue;
    //modify by pivot
    for (int i = 0; i < number; i++) {
      elementU[i] *= pivotValue;
    }
    maximumU_ = std::max(maximumU_, startU + number);
    numberRowsExtra_++;
    numberColumnsExtra_++;
    numberGoodU_++;
    numberPivots_++;
  }
  if (numberRowsExtra_ > numberRows_ + 50) {
    int extra = factorElements_ >> 1;

    if (numberRowsExtra_ > numberRows_ + 100 + numberRows_ / 500) {
      if (extra < 2 * numberRows_) {
        extra = 2 * numberRows_;
      }
    } else {
      if (extra < 5 * numberRows_) {
        extra = 5 * numberRows_;
      }
    }
    int added = totalElements_ - factorElements_;

    if (added > extra && added > (factorElements_) << 1 && !status
      && 3 * totalElements_ > 2 * (lengthAreaU_ + lengthAreaL_)) {
      status = 3;
      if (messageLevel_ & 4) {
        std::cout << "Factorization has " << totalElements_
                  << ", basis had " << factorElements_ << std::endl;
      }
    }
  }
  if (numberInColumnPlus && status < 2) {
    // we are going to put another copy of R in R
    CoinFactorizationDouble *COIN_RESTRICT elementR = elementR_ + lengthAreaR_;
    int *COIN_RESTRICT indexRowR = indexRowR_ + lengthAreaR_;
    int *COIN_RESTRICT startR = startColumnRArray_ + maximumPivots_ + 1;
    int pivotRow = numberRowsExtra_ - 1;
    for (int i = 0; i < numberNonZero; i++) {
      int iRow = regionIndex[i];
      assert(pivotRow > iRow);
      next = nextColumn[iRow];
      int space;
      if (next != maximumColumnsExtra_)
        space = startR[next] - startR[iRow];
      else
        space = lengthAreaR_ - startR[iRow];
      int numberInR = numberInColumnPlus[iRow];
      if (space > numberInR) {
        // there is space
        int put = startR[iRow] + numberInR;
        numberInColumnPlus[iRow] = numberInR + 1;
        indexRowR[put] = pivotRow;
        elementR[put] = region[iRow];
        //add 4 for luck
        if (next == maximumColumnsExtra_)
          startR[maximumColumnsExtra_] = std::min(static_cast< int >(put + 4), lengthAreaR_);
      } else {
        // no space - do we shuffle?
        if (!getColumnSpaceIterateR(iRow, region[iRow], pivotRow)) {
          // printf("Need more space for R\n");
          numberInColumnPlus_.conditionalDelete();
          numberInColumnPlusArray_ = NULL;
          regionSparse->clear();
          break;
        }
      }
      region[iRow] = 0.0;
    }
    regionSparse->setNumElements(0);
  } else {
    regionSparse->clear();
  }
#ifdef CLP_FACTORIZATION_INSTRUMENT
  numberReplace++;
  timeInReplace += CoinCpuTime() - startTimeX;
  currentLengthR = lengthR_;
  currentLengthU = lengthU_;
#endif
  return status;
}
#if ABOCA_LITE_FACTORIZATION
// Does btranU part of replaceColumn (skipping entries)
void CoinFactorization::replaceColumn1(CoinIndexedVector *regionSparse,
  int pivotRow)
{
  //return;
  //return at once if too many iterations
  if (numberColumnsExtra_ >= maximumColumnsExtra_) {
    return;
  }
  int *COIN_RESTRICT startColumnU = startColumnUArray_;
  if (lengthAreaU_ < startColumnU[maximumColumnsExtra_]) {
    return;
  }
  // If we have done no pivots then always check before modification
  if (!numberPivots_)
    return;

  assert(numberU_ <= numberRowsExtra_);
  CoinFactorizationDouble *COIN_RESTRICT element;
  int *COIN_RESTRICT numberInRow = numberInRowArray_;
  //int * COIN_RESTRICT numberInColumn = numberInColumnArray_;
  int realPivotRow;
  realPivotRow = pivotColumnArray_[pivotRow];
  //zeroed out region
  double *COIN_RESTRICT region = regionSparse->denseVector();

  element = elementUArray_;
  //get entries in row (pivot not stored)
  int numberNonZero = 0;
  int *COIN_RESTRICT indexColumn = indexColumnUArray_;
  int *COIN_RESTRICT convertRowToColumn = convertRowToColumnUArray_;
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  int *COIN_RESTRICT startRow = startRowUArray_;
  int start = 0;
  int end = 0;

#if COIN_ONE_ETA_COPY
  xxxxxx;
#endif
  start = startRow[realPivotRow];
  end = start + numberInRow[realPivotRow];

  int smallestIndex = numberRowsExtra_;
  for (int i = start; i < end; i++) {
    int iColumn = indexColumn[i];
    assert(iColumn < numberRowsExtra_);
    smallestIndex = std::min(smallestIndex, iColumn);
    int j = convertRowToColumn[i];

    region[iColumn] = element[j];
    //element[j] = 0.0;
    regionIndex[numberNonZero++] = iColumn;
  }
  //do BTRAN - finding first one to use
  regionSparse->setNumElements(numberNonZero);
  // split up later
  int number = regionSparse->getNumElements();
  int goSparse;
  // Guess at number at end
  if (sparseThreshold_ > 0) {
    if (btranAverageAfterU_) {
      double newNumber = number * btranAverageAfterU_;
      if (newNumber < sparseThreshold_)
        goSparse = 2;
      else if (newNumber < sparseThreshold2_)
        goSparse = 1;
      else
        goSparse = 0;
    } else {
      if (number < sparseThreshold_)
        goSparse = 2;
      else
        goSparse = 0;
    }
  } else {
    goSparse = 0;
  }
  // change and use second part of temp region
  goSparse = 0;
  switch (goSparse) {
  case 0: // densish
    //updateColumnTransposeUDensish(regionSparse,smallestIndex);
    {
      double tolerance = zeroTolerance_;
      int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
      const int *startRow = startRowUArray_;
      const int *convertRowToColumn = convertRowToColumnUArray_;
      const int *indexColumn = indexColumnUArray_;

      const CoinFactorizationDouble *COIN_RESTRICT element = elementUArray_;
      int last = numberU_;

      const int *numberInRow = numberInRowArray_;
      numberNonZero = 0;
      for (int i = smallestIndex; i < last; i++) {
        if (i == realPivotRow)
          continue; // skip - do we need to?
        CoinFactorizationDouble pivotValue = region[i];
        if (fabs(pivotValue) > tolerance) {
          int start = startRow[i];
          int numberIn = numberInRow[i];
          int end = start + (numberIn & (~1));
          for (int j = start; j < end; j += 2) {
            int iRow0 = indexColumn[j];
            int iRow1 = indexColumn[j + 1];
            int getElement0 = convertRowToColumn[j];
            int getElement1 = convertRowToColumn[j + 1];
            CoinFactorizationDouble value0 = element[getElement0];
            CoinFactorizationDouble value1 = element[getElement1];
            region[iRow0] -= value0 * pivotValue;
            region[iRow1] -= value1 * pivotValue;
          }
          if ((numberIn & 1) != 0) {
            int iRow = indexColumn[end];
            int getElement = convertRowToColumn[end];
            CoinFactorizationDouble value = element[getElement];
            region[iRow] -= value * pivotValue;
          }
          regionIndex[numberNonZero++] = i;
        } else {
          region[i] = 0.0;
        }
      }
      //set counts
      regionSparse->setNumElements(numberNonZero);
    }
    break;
  case 1: // middling
    updateColumnTransposeUSparsish(regionSparse, smallestIndex);
    break;
  case 2: // sparse
    updateColumnTransposeUSparse(regionSparse);
    break;
  }
}
// Does replaceColumn - having already done btranU
int CoinFactorization::replaceColumn2(CoinIndexedVector *regionSparse,
  int pivotRow,
  double pivotCheck)
{
#ifdef CLP_FACTORIZATION_INSTRUMENT
  double startTimeX = CoinCpuTime();
#endif
  assert(numberU_ <= numberRowsExtra_);
  int *COIN_RESTRICT startColumnU = startColumnUArray_;
  int *COIN_RESTRICT startColumn;
  int *COIN_RESTRICT indexRow;
  CoinFactorizationDouble *COIN_RESTRICT element;

  //return at once if too many iterations
  if (numberColumnsExtra_ >= maximumColumnsExtra_) {
    return 5;
  }
  if (lengthAreaU_ < startColumnU[maximumColumnsExtra_]) {
    return 3;
  }

  int *COIN_RESTRICT numberInRow = numberInRowArray_;
  int *COIN_RESTRICT numberInColumn = numberInColumnArray_;
  int *COIN_RESTRICT numberInColumnPlus = numberInColumnPlusArray_;
  int realPivotRow;
  realPivotRow = pivotColumnArray_[pivotRow];
  //zeroed out region
  double *COIN_RESTRICT region = regionSparse->denseVector();

  element = elementUArray_;
  //take out old pivot column

  totalElements_ -= numberInColumn[realPivotRow];
  CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;
  CoinFactorizationDouble oldPivot = pivotRegion[realPivotRow];
  // for accuracy check
  pivotCheck = pivotCheck / oldPivot;
#if COIN_DEBUG > 1
  int checkNumber = 1000000;
  //if (numberL_) checkNumber=-1;
  if (numberR_ >= checkNumber) {
    printf("pivot row %d, check %g - alpha region:\n",
      realPivotRow, pivotCheck);
    /*int i;
      for (i=0;i<numberRows_;i++) {
      if (pivotRegion[i])
      printf("%d %g\n",i,pivotRegion[i]);
  }*/
  }
#endif
  pivotRegion[realPivotRow] = 0.0;

  int saveEnd = startColumnU[realPivotRow]
    + numberInColumn[realPivotRow];
#ifdef CLP_FACTORIZATION_INSTRUMENT
  currentTakeoutU += numberInColumn[realPivotRow];
  currentTakeoutU += numberInRow[realPivotRow];
#endif
  // not necessary at present - but take no chances for future
  numberInColumn[realPivotRow] = 0;
  //get entries in row (pivot not stored)
  int numberNonZero = 0;
  int *COIN_RESTRICT indexColumn = indexColumnUArray_;
  int *COIN_RESTRICT convertRowToColumn = convertRowToColumnUArray_;
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  int *COIN_RESTRICT startRow = startRowUArray_;
  int start = 0;
  int end = 0;
#if COIN_DEBUG > 1
  if (numberR_ >= checkNumber)
    printf("Before btranu\n");
#endif

  start = startRow[realPivotRow];
  end = start + numberInRow[realPivotRow];

  for (int i = start; i < end; i++) {
    int iColumn = indexColumn[i];
    assert(iColumn < numberRowsExtra_);
    int j = convertRowToColumn[i];
    element[j] = 0.0;
  }
  numberNonZero = regionSparse->getNumElements();
  CoinFactorizationDouble saveFromU = 0.0;

  int startU = startColumnU[numberColumnsExtra_];
  int *COIN_RESTRICT indexU = &indexRowUArray_[startU];
  CoinFactorizationDouble *COIN_RESTRICT elementU = &elementUArray_[startU];

  // Now zero out column of U
  //take out old pivot column
  for (int i = startColumnU[realPivotRow]; i < saveEnd; i++) {
    element[i] = 0.0;
  }
  //zero out pivot Row (before or after?)
  //add to R
  startColumn = startColumnRArray_;
  indexRow = indexRowR_;
  element = elementR_;
  int l = lengthR_;
  int number = numberR_;

  startColumn[number] = l; //for luck and first time
  number++;
  startColumn[number] = l + numberNonZero;
  numberR_ = number;
  lengthR_ = l + numberNonZero;
  totalElements_ += numberNonZero;
  if (lengthR_ >= lengthAreaR_) {
    //not enough room
    regionSparse->clear();
    return 3;
  }
#if COIN_DEBUG > 1
  if (numberR_ >= checkNumber)
    printf("After btranu\n");
#endif
  for (int i = 0; i < numberNonZero; i++) {
    int iRow = regionIndex[i];
#if COIN_DEBUG > 1
    if (numberR_ >= checkNumber)
      printf("%d %g\n", iRow, region[iRow]);
#endif

    indexRow[l] = iRow;
    element[l] = region[iRow];
    l++;
  }
  int *COIN_RESTRICT nextRow;
  int *COIN_RESTRICT lastRow;
  int next;
  int last;
  //take out row
  nextRow = nextRowArray_;
  lastRow = lastRowArray_;
  next = nextRow[realPivotRow];
  last = lastRow[realPivotRow];

  nextRow[last] = next;
  lastRow[next] = last;
  numberInRow[realPivotRow] = 0;
#if COIN_DEBUG
  nextRow[realPivotRow] = 777777;
  lastRow[realPivotRow] = 777777;
#endif
  //do permute
  permuteArray_[numberRowsExtra_] = realPivotRow;
  // and other way
  permuteBackArray_[realPivotRow] = numberRowsExtra_;
  permuteBackArray_[numberRowsExtra_] = -1;
  ;
  //and for safety
  permuteArray_[numberRowsExtra_ + 1] = 0;

  pivotColumnArray_[pivotRow] = numberRowsExtra_;
  pivotColumnBack()[numberRowsExtra_] = pivotRow;
  startColumn = startColumnU;
  indexRow = indexRowUArray_;
  element = elementUArray_;

  numberU_++;
  number = numberInColumn[numberColumnsExtra_];

  totalElements_ += number;
  lengthU_ += number;
  if (lengthU_ >= lengthAreaU_) {
    //not enough room
    regionSparse->clear();
    return 3;
  }

  saveFromU = 0.0;

  //put in pivot
  //add row counts

#if COIN_DEBUG > 1
  if (numberR_ >= checkNumber)
    printf("On U\n");
#endif
  for (int i = 0; i < number; i++) {
    int iRow = indexU[i];
#if COIN_DEBUG > 1
    if (numberR_ >= checkNumber)
      printf("%d %g\n", iRow, elementU[i]);
#endif

    //assert ( fabs ( elementU[i] ) > zeroTolerance_ );
    if (iRow != realPivotRow) {
      int next = nextRow[iRow];
      int iNumberInRow = numberInRow[iRow];
      int space;
      int put = startRow[iRow] + iNumberInRow;

      space = startRow[next] - put;
      if (space <= 0) {
        getRowSpaceIterate(iRow, iNumberInRow + 4);
        put = startRow[iRow] + iNumberInRow;
      }
      indexColumn[put] = numberColumnsExtra_;
      convertRowToColumn[put] = i + startU;
      numberInRow[iRow] = iNumberInRow + 1;
      saveFromU = saveFromU - elementU[i] * region[iRow];
    } else {
      //zero out and save
      saveFromU += elementU[i];
      elementU[i] = 0.0;
    }
  }
  //in at end
  last = lastRow[maximumRowsExtra_];
  nextRow[last] = numberRowsExtra_;
  lastRow[maximumRowsExtra_] = numberRowsExtra_;
  lastRow[numberRowsExtra_] = last;
  nextRow[numberRowsExtra_] = maximumRowsExtra_;
  startRow[numberRowsExtra_] = startRow[maximumRowsExtra_];
  numberInRow[numberRowsExtra_] = 0;
  //column in at beginning (as empty)
  int *COIN_RESTRICT nextColumn = nextColumnArray_;
  int *COIN_RESTRICT lastColumn = lastColumnArray_;
  next = nextColumn[maximumColumnsExtra_];
  lastColumn[next] = numberColumnsExtra_;
  nextColumn[maximumColumnsExtra_] = numberColumnsExtra_;
  nextColumn[numberColumnsExtra_] = next;
  lastColumn[numberColumnsExtra_] = maximumColumnsExtra_;
  //check accuracy - but not if already checked (optimization problem)
  int status = checkPivot(saveFromU, pivotCheck);

  if (status != 2) {

    CoinFactorizationDouble pivotValue = 1.0 / saveFromU;

    pivotRegion[numberRowsExtra_] = pivotValue;
    //modify by pivot
    for (int i = 0; i < number; i++) {
      elementU[i] *= pivotValue;
    }
    maximumU_ = std::max(maximumU_, startU + number);
    numberRowsExtra_++;
    numberColumnsExtra_++;
    numberGoodU_++;
    numberPivots_++;
  }
  if (numberRowsExtra_ > numberRows_ + 50) {
    int extra = factorElements_ >> 1;

    if (numberRowsExtra_ > numberRows_ + 100 + numberRows_ / 500) {
      if (extra < 2 * numberRows_) {
        extra = 2 * numberRows_;
      }
    } else {
      if (extra < 5 * numberRows_) {
        extra = 5 * numberRows_;
      }
    }
    int added = totalElements_ - factorElements_;

    if (added > extra && added > (factorElements_) << 1 && !status
      && 3 * totalElements_ > 2 * (lengthAreaU_ + lengthAreaL_)) {
      status = 3;
      if (messageLevel_ & 4) {
        std::cout << "Factorization has " << totalElements_
                  << ", basis had " << factorElements_ << std::endl;
      }
    }
  }
  if (numberInColumnPlus && status < 2) {
    // we are going to put another copy of R in R
    CoinFactorizationDouble *COIN_RESTRICT elementR = elementR_ + lengthAreaR_;
    int *COIN_RESTRICT indexRowR = indexRowR_ + lengthAreaR_;
    int *COIN_RESTRICT startR = startColumnRArray_ + maximumPivots_ + 1;
    int pivotRow = numberRowsExtra_ - 1;
    for (int i = 0; i < numberNonZero; i++) {
      int iRow = regionIndex[i];
      assert(pivotRow > iRow);
      next = nextColumn[iRow];
      int space;
      if (next != maximumColumnsExtra_)
        space = startR[next] - startR[iRow];
      else
        space = lengthAreaR_ - startR[iRow];
      int numberInR = numberInColumnPlus[iRow];
      if (space > numberInR) {
        // there is space
        int put = startR[iRow] + numberInR;
        numberInColumnPlus[iRow] = numberInR + 1;
        indexRowR[put] = pivotRow;
        elementR[put] = region[iRow];
        //add 4 for luck
        if (next == maximumColumnsExtra_)
          startR[maximumColumnsExtra_] = std::min(static_cast< int >(put + 4), lengthAreaR_);
      } else {
        // no space - do we shuffle?
        if (!getColumnSpaceIterateR(iRow, region[iRow], pivotRow)) {
          // printf("Need more space for R\n");
          numberInColumnPlus_.conditionalDelete();
          numberInColumnPlusArray_ = NULL;
          regionSparse->clear();
          break;
        }
      }
      region[iRow] = 0.0;
    }
    regionSparse->setNumElements(0);
  } else {
    regionSparse->clear();
  }
#ifdef CLP_FACTORIZATION_INSTRUMENT
  numberReplace++;
  timeInReplace += CoinCpuTime() - startTimeX;
  currentLengthR = lengthR_;
  currentLengthU = lengthU_;
#endif
  return status;
}
#endif

//  updateColumnTranspose.  Updates one column transpose (BTRAN)
int CoinFactorization::updateColumnTranspose(CoinIndexedVector *regionSparse,
  CoinIndexedVector *regionSparse2)
  const
{
#ifdef CLP_FACTORIZATION_INSTRUMENT
  double startTimeX = CoinCpuTime();
#endif
  //zero region
  regionSparse->clear();
  double *COIN_RESTRICT region = regionSparse->denseVector();
  double *COIN_RESTRICT vector = regionSparse2->denseVector();
  int *COIN_RESTRICT index = regionSparse2->getIndices();
  int numberNonZero = regionSparse2->getNumElements();
  const int *pivotColumn = pivotColumnArray_;

  //move indices into index array
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  bool packed = regionSparse2->packedMode();
  if (packed) {
    for (int i = 0; i < numberNonZero; i++) {
      int iRow = index[i];
      double value = vector[i];
      iRow = pivotColumn[iRow];
      vector[i] = 0.0;
      region[iRow] = value;
      regionIndex[i] = iRow;
    }
  } else {
    for (int i = 0; i < numberNonZero; i++) {
      int iRow = index[i];
      double value = vector[iRow];
      vector[iRow] = 0.0;
      iRow = pivotColumn[iRow];
      region[iRow] = value;
      regionIndex[i] = iRow;
    }
  }
  regionSparse->setNumElements(numberNonZero);
  if (collectStatistics_) {
    numberBtranCounts_++;
    btranCountInput_ += static_cast< double >(numberNonZero);
  }
  if (!doForrestTomlin_) {
    // Do PFI before everything else
    updateColumnTransposePFI(regionSparse);
    numberNonZero = regionSparse->getNumElements();
  }
  //  ******* U
  // Apply pivot region - could be combined for speed
  CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;

  int smallestIndex = numberRowsExtra_;
  for (int j = 0; j < numberNonZero; j++) {
    int iRow = regionIndex[j];
    smallestIndex = std::min(smallestIndex, iRow);
    region[iRow] *= pivotRegion[iRow];
  }
  updateColumnTransposeU(regionSparse, smallestIndex);
  if (collectStatistics_)
    btranCountAfterU_ += static_cast< double >(regionSparse->getNumElements());
  //permute extra
  //row bits here
  updateColumnTransposeR(regionSparse);
  //  ******* L
  updateColumnTransposeL(regionSparse);
  numberNonZero = regionSparse->getNumElements();
  if (collectStatistics_) {
    btranCountAfterL_ += static_cast< double >(numberNonZero);
#ifdef CLP_FACTORIZATION_INSTRUMENT
    scaledLengthDense += numberDense_ * numberNonZero;
    scaledLengthDenseSquared += numberDense_ * numberDense_ * numberNonZero;
    scaledLengthL += lengthL_ * numberNonZero;
    scaledLengthR += lengthR_ * numberNonZero;
    scaledLengthU += lengthU_ * numberNonZero;
#endif
  }
  const int *permuteBack = pivotColumnBack();
  int number = 0;
  if (packed) {
    for (int i = 0; i < numberNonZero; i++) {
      int iRow = regionIndex[i];
      double value = region[iRow];
      region[iRow] = 0.0;
      //if (fabs(value)>zeroTolerance_) {
      iRow = permuteBack[iRow];
      vector[number] = value;
      index[number++] = iRow;
      //}
    }
  } else {
    for (int i = 0; i < numberNonZero; i++) {
      int iRow = regionIndex[i];
      double value = region[iRow];
      region[iRow] = 0.0;
      //if (fabs(value)>zeroTolerance_) {
      iRow = permuteBack[iRow];
      vector[iRow] = value;
      index[number++] = iRow;
      //}
    }
  }
  regionSparse->setNumElements(0);
  regionSparse2->setNumElements(number);
#ifdef COIN_DEBUG
  for (i = 0; i < numberRowsExtra_; i++) {
    assert(!region[i]);
  }
#endif
#ifdef CLP_FACTORIZATION_INSTRUMENT
  numberUpdateTranspose++;
  timeInUpdateTranspose += CoinCpuTime() - startTimeX;
  averageLengthR += lengthR_;
  averageLengthU += lengthU_;
  averageLengthL += lengthL_;
#endif
  return number;
}
#define TYPE_TWO_TRANSPOSE 1
#if TYPE_TWO_TRANSPOSE
#if ABOCA_LITE_FACTORIZATION
#include <cilk/cilk.h>
#endif
void CoinFactorization::updateOneColumnTranspose(CoinIndexedVector *regionWork,
  int &statistics) const
{
  int numberNonZero = regionWork->getNumElements();
  double *COIN_RESTRICT region = regionWork->denseVector();
  int *COIN_RESTRICT regionIndex = regionWork->getIndices();
  CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;
  if (!doForrestTomlin_) {
    // Do PFI before everything else
    updateColumnTransposePFI(regionWork);
    numberNonZero = regionWork->getNumElements();
  }
  //  ******* U
  // Apply pivot region - could be combined for speed

  int smallestIndex = numberRowsExtra_;
  for (int j = 0; j < numberNonZero; j++) {
    int iRow = regionIndex[j];
    smallestIndex = std::min(smallestIndex, iRow);
    region[iRow] *= pivotRegion[iRow];
  }
  updateColumnTransposeU(regionWork, smallestIndex);
  statistics = regionWork->getNumElements();
  //permute extra
  //row bits here
  updateColumnTransposeR(regionWork);
  //  ******* L
  updateColumnTransposeL(regionWork);
}
#endif
/* Updates two columns (BTRAN) from regionSparse2 and 3
   regionSparse starts as zero and is zero at end
   Note - if regionSparse2 packed on input - will be packed on output - same for 3
   NOTE NOTE - Now we assume 2 is packed and 3 is not
   If changes then need to checl setNumElements as that can change packed mode
*/
void CoinFactorization::updateTwoColumnsTranspose(CoinIndexedVector *regionSparse,
  CoinIndexedVector *regionSparse2,
  CoinIndexedVector *regionSparse3,
  int type) const
{
  /*
    0 - two separate
    1 - two separate but do permutes here
    2 - do in three chunks
    3 - do in three chunks changing coding
    11 - ABOCA_LITE_FACTORIZATION and as 1 (apart from permute)
    12 - ABOCA_LITE_FACTORIZATION and 2 ????
  */
#if TYPE_TWO_TRANSPOSE == 0
  updateColumnTranspose(regionSparse, regionSparse3);
  updateColumnTranspose(regionSparse, regionSparse2);
#else
#ifdef CLP_FACTORIZATION_INSTRUMENT
  double startTimeX = CoinCpuTime();
#endif
  const int *pivotColumn = pivotColumnArray_;
  //zero region
  CoinIndexedVector *regionWorkA = regionSparse;
  regionWorkA->clear();
  double *COIN_RESTRICT regionA = regionWorkA->denseVector();
  double *COIN_RESTRICT vectorA = regionSparse3->denseVector();
  int *COIN_RESTRICT indexA = regionSparse3->getIndices();
  int numberNonZeroA = regionSparse3->getNumElements();
  int *COIN_RESTRICT regionIndexA = regionWorkA->getIndices();

  //move indices into index array
  bool packedA = regionSparse3->packedMode();
  assert(!packedA);
  packedA = false;
  if (packedA) {
    for (int i = 0; i < numberNonZeroA; i++) {
      int iRow = indexA[i];
      double value = vectorA[i];
      iRow = pivotColumn[iRow];
      vectorA[i] = 0.0;
      regionA[iRow] = value;
      regionIndexA[i] = iRow;
    }
  } else {
    for (int i = 0; i < numberNonZeroA; i++) {
      int iRow = indexA[i];
      double value = vectorA[iRow];
      vectorA[iRow] = 0.0;
      iRow = pivotColumn[iRow];
      regionA[iRow] = value;
      regionIndexA[i] = iRow;
    }
  }
  regionWorkA->setNumElements(numberNonZeroA);
  CoinIndexedVector *regionWorkB = regionSparse3;
  double *COIN_RESTRICT regionB = regionWorkB->denseVector();
  double *COIN_RESTRICT vectorB = regionSparse2->denseVector();
  int *COIN_RESTRICT indexB = regionSparse2->getIndices();
  int numberNonZeroB = regionSparse2->getNumElements();
  int *COIN_RESTRICT regionIndexB = regionWorkB->getIndices();
  bool packedB = regionSparse2->packedMode();
  assert(packedB);
  packedB = true;
  if (packedB) {
    for (int i = 0; i < numberNonZeroB; i++) {
      int iRow = indexB[i];
      double value = vectorB[i];
      iRow = pivotColumn[iRow];
      vectorB[i] = 0.0;
      regionB[iRow] = value;
      regionIndexB[i] = iRow;
    }
  } else {
    for (int i = 0; i < numberNonZeroB; i++) {
      int iRow = indexB[i];
      double value = vectorB[iRow];
      vectorB[iRow] = 0.0;
      iRow = pivotColumn[iRow];
      regionB[iRow] = value;
      regionIndexB[i] = iRow;
    }
  }
  regionWorkB->setNumElements(numberNonZeroB);
  if (collectStatistics_) {
    numberBtranCounts_ += 2;
    btranCountInput_ += static_cast< double >(numberNonZeroA + numberNonZeroB);
  }
  int statistics[2];
#if TYPE_TWO_TRANSPOSE == 1
#ifdef ABOCA_LITE_FACTORIZATION
  if (!type) {
#endif
    CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;
    if (!doForrestTomlin_) {
      // Do PFI before everything else
      updateColumnTransposePFI(regionWorkA);
      numberNonZeroA = regionWorkA->getNumElements();
    }
    //  ******* U
    // Apply pivot region - could be combined for speed

    int smallestIndexA = numberRowsExtra_;
    for (int j = 0; j < numberNonZeroA; j++) {
      int iRow = regionIndexA[j];
      smallestIndexA = std::min(smallestIndexA, iRow);
      regionA[iRow] *= pivotRegion[iRow];
    }
    updateColumnTransposeU(regionWorkA, smallestIndexA);
    statistics[0] = regionWorkA->getNumElements();
    //permute extra
    //row bits here
    updateColumnTransposeR(regionWorkA);
    //  ******* L
    updateColumnTransposeL(regionWorkA);
    if (!doForrestTomlin_) {
      // Do PFI before everything else
      updateColumnTransposePFI(regionWorkB);
      numberNonZeroB = regionWorkB->getNumElements();
    }
    //  ******* U
    // Apply pivot region - could be combined for speed

    int smallestIndexB = numberRowsExtra_;
    for (int j = 0; j < numberNonZeroB; j++) {
      int iRow = regionIndexB[j];
      smallestIndexB = std::min(smallestIndexB, iRow);
      regionB[iRow] *= pivotRegion[iRow];
    }
    updateColumnTransposeU(regionWorkB, smallestIndexB);
    statistics[1] = regionWorkB->getNumElements();
    //permute extra
    //row bits here
    updateColumnTransposeR(regionWorkB);
    //  ******* L
    updateColumnTransposeL(regionWorkB);
#ifdef ABOCA_LITE_FACTORIZATION
  } else {
    regionWorkB->setCapacity(regionWorkB->capacity() | 0x80000000);
    cilk_spawn updateOneColumnTranspose(regionWorkA, statistics[0]);
    //cilk_sync;
    cilk_spawn updateOneColumnTranspose(regionWorkB, statistics[1]);
    cilk_sync;
    regionWorkB->setCapacity(regionWorkB->capacity() & 0x7fffffff);
    numberNonZeroA = regionWorkA->getNumElements();
    numberNonZeroA = regionWorkB->getNumElements();
  }
#endif
#elif TYPE_TWO_TRANSPOSE == 2
#ifdef ABOCA_LITE_FACTORIZATION
  if (!type) {
#endif
    CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;
    if (!doForrestTomlin_) {
      // Do PFI before everything else
      updateColumnTransposePFI(regionWorkA);
      numberNonZeroA = regionWorkA->getNumElements();
      updateColumnTransposePFI(regionWorkB);
      numberNonZeroB = regionWorkB->getNumElements();
    }
    //  ******* U
    // Apply pivot region - could be combined for speed
    int smallestIndexA = numberRowsExtra_;
    for (int j = 0; j < numberNonZeroA; j++) {
      int iRow = regionIndexA[j];
      smallestIndexA = std::min(smallestIndexA, iRow);
      regionA[iRow] *= pivotRegion[iRow];
    }
    int smallestIndexB = numberRowsExtra_;
    for (int j = 0; j < numberNonZeroB; j++) {
      int iRow = regionIndexB[j];
      smallestIndexB = std::min(smallestIndexB, iRow);
      regionB[iRow] *= pivotRegion[iRow];
    }
    updateColumnTransposeU(regionWorkA, smallestIndexA);
    updateColumnTransposeU(regionWorkB, smallestIndexB);
    statistics[0] = regionWorkA->getNumElements();
    statistics[1] = regionWorkB->getNumElements();
    //permute extra
    //row bits here
    updateColumnTransposeR(regionWorkA);
    updateColumnTransposeR(regionWorkB);
    //  ******* L
    updateColumnTransposeL(regionWorkA);
    updateColumnTransposeL(regionWorkB);
#ifdef ABOCA_LITE_FACTORIZATION
  } else {
    regionWorkB->setCapacity(regionWorkB->capacity() | 0x80000000);
    cilk_spawn updateOneColumnTranspose(regionWorkA, statistics[0]);
    //cilk_sync;
    cilk_spawn updateOneColumnTranspose(regionWorkB, statistics[1]);
    cilk_sync;
    regionWorkB->setCapacity(regionWorkB->capacity() & 0x7fffffff);
    numberNonZeroA = regionWorkA->getNumElements();
    numberNonZeroA = regionWorkB->getNumElements();
  }
#endif
#else
  CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;
  if (!doForrestTomlin_) {
    // Do PFI before everything else
    updateColumnTransposePFI(regionWorkA);
    numberNonZeroA = regionWorkA->getNumElements();
    updateColumnTransposePFI(regionWorkB);
    numberNonZeroB = regionWorkB->getNumElements();
  }
  //  ******* U
  // Apply pivot region - could be combined for speed
  int smallestIndexA = numberRowsExtra_;
  for (int j = 0; j < numberNonZeroA; j++) {
    int iRow = regionIndexA[j];
    smallestIndexA = std::min(smallestIndexA, iRow);
    regionA[iRow] *= pivotRegion[iRow];
  }
  int smallestIndexB = numberRowsExtra_;
  for (int j = 0; j < numberNonZeroB; j++) {
    int iRow = regionIndexB[j];
    smallestIndexB = std::min(smallestIndexB, iRow);
    regionB[iRow] *= pivotRegion[iRow];
  }
  updateColumnTransposeU(regionWorkA, smallestIndexA);
  updateColumnTransposeU(regionWorkB, smallestIndexB);
  statistics[0] = regionWorkA->getNumElements();
  statistics[1] = regionWorkB->getNumElements();
  //permute extra
  //row bits here
  updateColumnTransposeR(regionWorkA);
  updateColumnTransposeR(regionWorkB);
  //  ******* L
#ifdef COIN_FACTORIZATION_DENSE_CODE
  if (!numberDense_) {
#endif
    updateColumnTransposeL(regionWorkA);
    updateColumnTransposeL(regionWorkB);
#ifdef COIN_FACTORIZATION_DENSE_CODE
  } else {
    abort();
  }
#endif
#endif
  if (collectStatistics_) {
    btranCountAfterL_ += static_cast< double >(numberNonZeroA + numberNonZeroB);
    btranCountAfterU_ += static_cast< double >(statistics[0] + statistics[1]);
#ifdef CLP_FACTORIZATION_INSTRUMENT
    scaledLengthDense += numberDense_ * numberNonZeroA;
    scaledLengthDenseSquared += numberDense_ * numberDense_ * numberNonZeroA;
    scaledLengthL += lengthL_ * numberNonZeroA;
    scaledLengthR += lengthR_ * numberNonZeroA;
    scaledLengthU += lengthU_ * numberNonZeroA;
#endif
  }
  const int *permuteBack = pivotColumnBack();
  numberNonZeroA = regionWorkA->getNumElements();
  numberNonZeroB = regionWorkB->getNumElements();
  int number = 0;
  if (packedB) {
    for (int i = 0; i < numberNonZeroB; i++) {
      int iRow = regionIndexB[i];
      double value = regionB[iRow];
      regionB[iRow] = 0.0;
      //if (fabs(value)>zeroTolerance_) {
      iRow = permuteBack[iRow];
      vectorB[number] = value;
      indexB[number++] = iRow;
      //}
    }
  } else {
    for (int i = 0; i < numberNonZeroB; i++) {
      int iRow = regionIndexB[i];
      double value = regionB[iRow];
      regionB[iRow] = 0.0;
      //if (fabs(value)>zeroTolerance_) {
      iRow = permuteBack[iRow];
      vectorB[iRow] = value;
      indexB[number++] = iRow;
      //}
    }
  }
  //regionWorkB->setNumElements(0);
  regionSparse2->setNumElements(number);
  number = 0;
  if (packedA) {
    for (int i = 0; i < numberNonZeroA; i++) {
      int iRow = regionIndexA[i];
      double value = regionA[iRow];
      regionA[iRow] = 0.0;
      //if (fabs(value)>zeroTolerance_) {
      iRow = permuteBack[iRow];
      vectorA[number] = value;
      indexA[number++] = iRow;
      //}
    }
  } else {
    for (int i = 0; i < numberNonZeroA; i++) {
      int iRow = regionIndexA[i];
      double value = regionA[iRow];
      regionA[iRow] = 0.0;
      //if (fabs(value)>zeroTolerance_) {
      iRow = permuteBack[iRow];
      vectorA[iRow] = value;
      indexA[number++] = iRow;
      //}
    }
  }
  regionWorkA->setNumElements(0);
  regionSparse3->setNumElements(number);
#ifdef COIN_DEBUG
  for (i = 0; i < numberRowsExtra_; i++) {
    assert(!regionA[i]);
  }
#endif
#ifdef CLP_FACTORIZATION_INSTRUMENT
  numberUpdateTranspose += 2;
  timeInUpdateTranspose += CoinCpuTime() - startTimeX;
  averageLengthR += 2 * lengthR_;
  averageLengthU += 2 * lengthU_;
  averageLengthL += 2 * lengthL_;
#endif
#endif
}

/* Updates part of column transpose (BTRANU) when densish,
   assumes index is sorted i.e. region is correct */
void CoinFactorization::updateColumnTransposeUDensish(CoinIndexedVector *regionSparse,
  int smallestIndex) const
{
  double *COIN_RESTRICT region = regionSparse->denseVector();
  int numberNonZero = regionSparse->getNumElements();
  double tolerance = zeroTolerance_;

  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();

  const int *startRow = startRowUArray_;

  const int *convertRowToColumn = convertRowToColumnUArray_;
  const int *indexColumn = indexColumnUArray_;

  const CoinFactorizationDouble *COIN_RESTRICT element = elementUArray_;
  int last = numberU_;

  const int *numberInRow = numberInRowArray_;
  numberNonZero = 0;
  for (int i = smallestIndex; i < last; i++) {
    CoinFactorizationDouble pivotValue = region[i];
    if (fabs(pivotValue) > tolerance) {
      int start = startRow[i];
      int numberIn = numberInRow[i];
      int end = start + (numberIn & (~1));
      for (int j = start; j < end; j += 2) {
        int iRow0 = indexColumn[j];
        int iRow1 = indexColumn[j + 1];
        int getElement0 = convertRowToColumn[j];
        int getElement1 = convertRowToColumn[j + 1];
        CoinFactorizationDouble value0 = element[getElement0];
        CoinFactorizationDouble value1 = element[getElement1];
        region[iRow0] -= value0 * pivotValue;
        region[iRow1] -= value1 * pivotValue;
      }
      if ((numberIn & 1) != 0) {
        int iRow = indexColumn[end];
        int getElement = convertRowToColumn[end];
        CoinFactorizationDouble value = element[getElement];
        region[iRow] -= value * pivotValue;
      }
      regionIndex[numberNonZero++] = i;
    } else {
      region[i] = 0.0;
    }
  }
  //set counts
  regionSparse->setNumElements(numberNonZero);
}
/* Updates part of column transpose (BTRANU) when sparsish,
      assumes index is sorted i.e. region is correct */
void CoinFactorization::updateColumnTransposeUSparsish(CoinIndexedVector *regionSparse,
  int smallestIndex) const
{
  double *COIN_RESTRICT region = regionSparse->denseVector();
  int numberNonZero = regionSparse->getNumElements();
  double tolerance = zeroTolerance_;

  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();

  const int *startRow = startRowUArray_;

  const int *convertRowToColumn = convertRowToColumnUArray_;
  const int *indexColumn = indexColumnUArray_;

  const CoinFactorizationDouble *COIN_RESTRICT element = elementUArray_;
  int last = numberU_;

  const int *numberInRow = numberInRowArray_;

  // mark known to be zero
  int nInBig = sizeof(int) / sizeof(int);
#if ABOCA_LITE_FACTORIZATION == 0
#define sparseOffset 0
#else
  int sparseOffset = ((regionSparse->capacity() & 0x80000000) != 0) ? sparseOffset_ : 0;
  assert(!sparseOffset);
#endif
  CoinCheckZero *COIN_RESTRICT mark = reinterpret_cast< CoinCheckZero * >(sparseArray_ + (2 + nInBig) * maximumRowsExtra_ + sparseOffset);

  for (int i = 0; i < numberNonZero; i++) {
    int iPivot = regionIndex[i];
    int iWord = iPivot >> CHECK_SHIFT;
    int iBit = iPivot - (iWord << CHECK_SHIFT);
    if (mark[iWord]) {
      mark[iWord] = static_cast< CoinCheckZero >(mark[iWord] | (1 << iBit));
    } else {
      mark[iWord] = static_cast< CoinCheckZero >(1 << iBit);
    }
  }

  numberNonZero = 0;
  // Find convenient power of 2
  smallestIndex = smallestIndex >> CHECK_SHIFT;
  int kLast = last >> CHECK_SHIFT;
  // do in chunks

  for (int k = smallestIndex; k < kLast; k++) {
    unsigned int iMark = mark[k];
    if (iMark) {
      // something in chunk - do all (as imark may change)
      int i = k << CHECK_SHIFT;
      int iLast = i + BITS_PER_CHECK;
      for (; i < iLast; i++) {
        CoinFactorizationDouble pivotValue = region[i];
        if (fabs(pivotValue) > tolerance) {
          int start = startRow[i];
          int numberIn = numberInRow[i];
          int end = start + numberIn;
          for (int j = start; j < end; j++) {
            int iRow = indexColumn[j];
            int getElement = convertRowToColumn[j];
            CoinFactorizationDouble value = element[getElement];
            int iWord = iRow >> CHECK_SHIFT;
            int iBit = iRow - (iWord << CHECK_SHIFT);
            if (mark[iWord]) {
              mark[iWord] = static_cast< CoinCheckZero >(mark[iWord] | (1 << iBit));
            } else {
              mark[iWord] = static_cast< CoinCheckZero >(1 << iBit);
            }
            region[iRow] -= value * pivotValue;
          }
          regionIndex[numberNonZero++] = i;
        } else {
          region[i] = 0.0;
        }
      }
      mark[k] = 0;
    }
  }
  mark[kLast] = 0;
  for (int i = kLast << CHECK_SHIFT; i < last; i++) {
    CoinFactorizationDouble pivotValue = region[i];
    if (fabs(pivotValue) > tolerance) {
      int start = startRow[i];
      int numberIn = numberInRow[i];
      int end = start + numberIn;
      for (int j = start; j < end; j++) {
        int iRow = indexColumn[j];
        int getElement = convertRowToColumn[j];
        CoinFactorizationDouble value = element[getElement];

        region[iRow] -= value * pivotValue;
      }
      regionIndex[numberNonZero++] = i;
    } else {
      region[i] = 0.0;
    }
  }
#ifdef COIN_DEBUG
  for (int i = 0; i < maximumRowsExtra_; i++) {
    assert(!mark[i]);
  }
#endif
  //set counts
  regionSparse->setNumElements(numberNonZero);
}
/* Updates part of column transpose (BTRANU) when sparse,
   assumes index is sorted i.e. region is correct */
void CoinFactorization::updateColumnTransposeUSparse(
  CoinIndexedVector *regionSparse) const
{
  double *COIN_RESTRICT region = regionSparse->denseVector();
  int numberNonZero = regionSparse->getNumElements();
  double tolerance = zeroTolerance_;

  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  const int *startRow = startRowUArray_;

  const int *convertRowToColumn = convertRowToColumnUArray_;
  const int *indexColumn = indexColumnUArray_;

  const CoinFactorizationDouble *COIN_RESTRICT element = elementUArray_;

  const int *numberInRow = numberInRowArray_;

  // use sparse_ as temporary area
  // mark known to be zero
#if ABOCA_LITE_FACTORIZATION == 0
#define sparseOffset 0
#else
  int sparseOffset = ((regionSparse->capacity() & 0x80000000) != 0) ? sparseOffset_ : 0;
#endif
  int *COIN_RESTRICT stack = sparseArray_ + sparseOffset; /* pivot */
  int *COIN_RESTRICT list = stack + maximumRowsExtra_; /* final list */
  int *COIN_RESTRICT next = reinterpret_cast< int * >(list + maximumRowsExtra_); /* jnext */
  char *COIN_RESTRICT mark = reinterpret_cast< char * >(next + maximumRowsExtra_);
  int nList;
#ifdef COIN_DEBUG
  for (int i = 0; i < maximumRowsExtra_; i++) {
    assert(!mark[i]);
  }
#endif
  nList = 0;
  for (int k = 0; k < numberNonZero; k++) {
    int kPivot = regionIndex[k];
    stack[0] = kPivot;
    int j = startRow[kPivot] + numberInRow[kPivot] - 1;
    next[0] = j;
    int nStack = 1;
    while (nStack) {
      /* take off stack */
      kPivot = stack[--nStack];
      if (mark[kPivot] != 1) {
        j = next[nStack];
        if (j >= startRow[kPivot]) {
          kPivot = indexColumn[j--];
          /* put back on stack */
          next[nStack++] = j;
          if (!mark[kPivot]) {
            /* and new one */
            j = startRow[kPivot] + numberInRow[kPivot] - 1;
            stack[nStack] = kPivot;
            mark[kPivot] = 2;
            next[nStack++] = j;
          }
        } else {
          // finished
          list[nList++] = kPivot;
          mark[kPivot] = 1;
        }
      }
    }
  }
  numberNonZero = 0;
  for (int i = nList - 1; i >= 0; i--) {
    int iPivot = list[i];
    mark[iPivot] = 0;
    CoinFactorizationDouble pivotValue = region[iPivot];
    if (fabs(pivotValue) > tolerance) {
      int start = startRow[iPivot];
      int numberIn = numberInRow[iPivot];
      int end = start + numberIn;
      for (int j = start; j < end; j++) {
        int iRow = indexColumn[j];
        int getElement = convertRowToColumn[j];
        CoinFactorizationDouble value = element[getElement];
        region[iRow] -= value * pivotValue;
      }
      regionIndex[numberNonZero++] = iPivot;
    } else {
      region[iPivot] = 0.0;
    }
  }
  //set counts
  regionSparse->setNumElements(numberNonZero);
}
//  updateColumnTransposeU.  Updates part of column transpose (BTRANU)
//assumes index is sorted i.e. region is correct
//does not sort by sign
void CoinFactorization::updateColumnTransposeU(CoinIndexedVector *regionSparse,
  int smallestIndex) const
{
#if COIN_ONE_ETA_COPY
  int *COIN_RESTRICT convertRowToColumn = convertRowToColumnUArray_;
  if (!convertRowToColumn) {
    //abort();
    updateColumnTransposeUByColumn(regionSparse, smallestIndex);
    return;
  }
#endif
  int number = regionSparse->getNumElements();
  int goSparse;
  // Guess at number at end
  if (sparseThreshold_ > 0) {
    if (btranAverageAfterU_) {
      double newNumber = number * btranAverageAfterU_;
      if (newNumber < sparseThreshold_)
        goSparse = 2;
      else if (newNumber < sparseThreshold2_)
        goSparse = 1;
      else
        goSparse = 0;
    } else {
      if (number < sparseThreshold_)
        goSparse = 2;
      else
        goSparse = 0;
    }
  } else {
    goSparse = 0;
  }
  switch (goSparse) {
  case 0: // densish
    updateColumnTransposeUDensish(regionSparse, smallestIndex);
    break;
  case 1: // middling
    updateColumnTransposeUSparsish(regionSparse, smallestIndex);
    break;
  case 2: // sparse
    updateColumnTransposeUSparse(regionSparse);
    break;
  }
}

/*  updateColumnTransposeLDensish.
    Updates part of column transpose (BTRANL) dense by column */
void CoinFactorization::updateColumnTransposeLDensish(CoinIndexedVector *regionSparse) const
{
  double *COIN_RESTRICT region = regionSparse->denseVector();
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  int numberNonZero;
  double tolerance = zeroTolerance_;
  int base;
  int first = -1;

  numberNonZero = 0;
  //scan
  for (first = numberRows_ - 1; first >= 0; first--) {
    if (region[first])
      break;
  }
  if (first >= 0) {
    base = baseL_;
    const int *COIN_RESTRICT startColumn = startColumnLArray_;
    const int *COIN_RESTRICT indexRow = indexRowLArray_;
    const CoinFactorizationDouble *COIN_RESTRICT element = elementLArray_;
    int last = baseL_ + numberL_;

    if (first >= last) {
      first = last - 1;
    }
    for (int i = first; i >= base; i--) {
      int j;
      CoinFactorizationDouble pivotValue = region[i];
      for (j = startColumn[i]; j < startColumn[i + 1]; j++) {
        int iRow = indexRow[j];
        CoinFactorizationDouble value = element[j];
        pivotValue -= value * region[iRow];
      }
      if (fabs(pivotValue) > tolerance) {
        region[i] = pivotValue;
        regionIndex[numberNonZero++] = i;
      } else {
        region[i] = 0.0;
      }
    }
    //may have stopped early
    if (first < base) {
      base = first + 1;
    }
    if (base > 5) {
      int i = base - 1;
      CoinFactorizationDouble pivotValue = region[i];
      bool store = fabs(pivotValue) > tolerance;
      for (; i > 0; i--) {
        bool oldStore = store;
        CoinFactorizationDouble oldValue = pivotValue;
        pivotValue = region[i - 1];
        store = fabs(pivotValue) > tolerance;
        if (!oldStore) {
          region[i] = 0.0;
        } else {
          region[i] = oldValue;
          regionIndex[numberNonZero++] = i;
        }
      }
      if (store) {
        region[0] = pivotValue;
        regionIndex[numberNonZero++] = 0;
      } else {
        region[0] = 0.0;
      }
    } else {
      for (int i = base - 1; i >= 0; i--) {
        CoinFactorizationDouble pivotValue = region[i];
        if (fabs(pivotValue) > tolerance) {
          region[i] = pivotValue;
          regionIndex[numberNonZero++] = i;
        } else {
          region[i] = 0.0;
        }
      }
    }
  }
  //set counts
  regionSparse->setNumElements(numberNonZero);
}
/*  updateColumnTransposeLByRow.
    Updates part of column transpose (BTRANL) densish but by row */
void CoinFactorization::updateColumnTransposeLByRow(CoinIndexedVector *regionSparse) const
{
  double *COIN_RESTRICT region = regionSparse->denseVector();
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  int numberNonZero;
  double tolerance = zeroTolerance_;
  int first = -1;

  // use row copy of L
  const CoinFactorizationDouble *COIN_RESTRICT element = elementByRowLArray_;
  const int *startRow = startRowLArray_;
  const int *column = indexColumnLArray_;
  for (first = numberRows_ - 1; first >= 0; first--) {
    if (region[first])
      break;
  }
  numberNonZero = 0;
  for (int i = first; i >= 0; i--) {
    CoinFactorizationDouble pivotValue = region[i];
    if (fabs(pivotValue) > tolerance) {
      regionIndex[numberNonZero++] = i;
      int j;
      for (j = startRow[i + 1] - 1; j >= startRow[i]; j--) {
        int iRow = column[j];
        CoinFactorizationDouble value = element[j];
        region[iRow] -= pivotValue * value;
      }
    } else {
      region[i] = 0.0;
    }
  }
  //set counts
  regionSparse->setNumElements(numberNonZero);
}
// Updates part of column transpose (BTRANL) when sparsish by row
void CoinFactorization::updateColumnTransposeLSparsish(CoinIndexedVector *regionSparse) const
{
  double *COIN_RESTRICT region = regionSparse->denseVector();
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  int numberNonZero = regionSparse->getNumElements();
  double tolerance = zeroTolerance_;

  // use row copy of L
  const CoinFactorizationDouble *COIN_RESTRICT element = elementByRowLArray_;
  const int *startRow = startRowLArray_;
  const int *column = indexColumnLArray_;
  // mark known to be zero
  int nInBig = sizeof(int) / sizeof(int);
#if ABOCA_LITE_FACTORIZATION == 0
#define sparseOffset 0
#else
  int sparseOffset = ((regionSparse->capacity() & 0x80000000) != 0) ? sparseOffset_ : 0;
#endif
  CoinCheckZero *COIN_RESTRICT mark = reinterpret_cast< CoinCheckZero * >(sparseArray_ + (2 + nInBig) * maximumRowsExtra_ + sparseOffset);
  for (int i = 0; i < numberNonZero; i++) {
    int iPivot = regionIndex[i];
    int iWord = iPivot >> CHECK_SHIFT;
    int iBit = iPivot - (iWord << CHECK_SHIFT);
    if (mark[iWord]) {
      mark[iWord] = static_cast< CoinCheckZero >(mark[iWord] | (1 << iBit));
    } else {
      mark[iWord] = static_cast< CoinCheckZero >(1 << iBit);
    }
  }
  numberNonZero = 0;
  // First do down to convenient power of 2
  int jLast = (numberRows_ - 1) >> CHECK_SHIFT;
  jLast = (jLast << CHECK_SHIFT);
  for (int i = numberRows_ - 1; i >= jLast; i--) {
    CoinFactorizationDouble pivotValue = region[i];
    if (fabs(pivotValue) > tolerance) {
      regionIndex[numberNonZero++] = i;
      int j;
      for (j = startRow[i + 1] - 1; j >= startRow[i]; j--) {
        int iRow = column[j];
        CoinFactorizationDouble value = element[j];
        int iWord = iRow >> CHECK_SHIFT;
        int iBit = iRow - (iWord << CHECK_SHIFT);
        if (mark[iWord]) {
          mark[iWord] = static_cast< CoinCheckZero >(mark[iWord] | (1 << iBit));
        } else {
          mark[iWord] = static_cast< CoinCheckZero >(1 << iBit);
        }
        region[iRow] -= pivotValue * value;
      }
    } else {
      region[i] = 0.0;
    }
  }
  // and in chunks
  jLast = jLast >> CHECK_SHIFT;
  mark[jLast] = 0;
  for (int k = jLast - 1; k >= 0; k--) {
    unsigned int iMark = mark[k];
    if (iMark) {
      // something in chunk - do all (as imark may change)
      int iLast = k << CHECK_SHIFT;
      for (int i = iLast + BITS_PER_CHECK - 1; i >= iLast; i--) {
        CoinFactorizationDouble pivotValue = region[i];
        if (fabs(pivotValue) > tolerance) {
          regionIndex[numberNonZero++] = i;
          int j;
          for (j = startRow[i + 1] - 1; j >= startRow[i]; j--) {
            int iRow = column[j];
            CoinFactorizationDouble value = element[j];
            int iWord = iRow >> CHECK_SHIFT;
            int iBit = iRow - (iWord << CHECK_SHIFT);
            if (mark[iWord]) {
              mark[iWord] = static_cast< CoinCheckZero >(mark[iWord] | (1 << iBit));
            } else {
              mark[iWord] = static_cast< CoinCheckZero >(1 << iBit);
            }
            region[iRow] -= pivotValue * value;
          }
        } else {
          region[i] = 0.0;
        }
      }
      mark[k] = 0;
    }
  }
#ifdef COIN_DEBUG
  for (int i = 0; i < maximumRowsExtra_; i++) {
    assert(!mark[i]);
  }
#endif
  //set counts
  regionSparse->setNumElements(numberNonZero);
}
/*  updateColumnTransposeLSparse.
    Updates part of column transpose (BTRANL) sparse */
void CoinFactorization::updateColumnTransposeLSparse(CoinIndexedVector *regionSparse) const
{
  double *COIN_RESTRICT region = regionSparse->denseVector();
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  int numberNonZero = regionSparse->getNumElements();
  double tolerance = zeroTolerance_;

  // use row copy of L
  const CoinFactorizationDouble *COIN_RESTRICT element = elementByRowLArray_;
  const int *startRow = startRowLArray_;
  const int *column = indexColumnLArray_;
  // use sparse_ as temporary area
  // mark known to be zero
#if ABOCA_LITE_FACTORIZATION == 0
#define sparseOffset 0
#else
  int sparseOffset = ((regionSparse->capacity() & 0x80000000) != 0) ? sparseOffset_ : 0;
#endif
  int *COIN_RESTRICT stack = sparseArray_ + sparseOffset; /* pivot */
  int *COIN_RESTRICT list = stack + maximumRowsExtra_; /* final list */
  int *COIN_RESTRICT next = reinterpret_cast< int * >(list + maximumRowsExtra_); /* jnext */
  char *COIN_RESTRICT mark = reinterpret_cast< char * >(next + maximumRowsExtra_);
  int nList;
  int number = numberNonZero;
#ifdef COIN_DEBUG
  for (i = 0; i < maximumRowsExtra_; i++) {
    assert(!mark[i]);
  }
#endif
  nList = 0;
  for (int k = 0; k < number; k++) {
    int kPivot = regionIndex[k];
    if (!mark[kPivot] && region[kPivot]) {
      stack[0] = kPivot;
      int j = startRow[kPivot + 1] - 1;
      int nStack = 0;
      while (nStack >= 0) {
        /* take off stack */
        if (j >= startRow[kPivot]) {
          int jPivot = column[j--];
          /* put back on stack */
          next[nStack] = j;
          if (!mark[jPivot]) {
            /* and new one */
            kPivot = jPivot;
            j = startRow[kPivot + 1] - 1;
            stack[++nStack] = kPivot;
            mark[kPivot] = 1;
            next[nStack] = j;
          }
        } else {
          /* finished so mark */
          list[nList++] = kPivot;
          mark[kPivot] = 1;
          --nStack;
          if (nStack >= 0) {
            kPivot = stack[nStack];
            j = next[nStack];
          }
        }
      }
    }
  }
  numberNonZero = 0;
  for (int i = nList - 1; i >= 0; i--) {
    int iPivot = list[i];
    mark[iPivot] = 0;
    CoinFactorizationDouble pivotValue = region[iPivot];
    if (fabs(pivotValue) > tolerance) {
      regionIndex[numberNonZero++] = iPivot;
      int j;
      for (j = startRow[iPivot]; j < startRow[iPivot + 1]; j++) {
        int iRow = column[j];
        CoinFactorizationDouble value = element[j];
        region[iRow] -= value * pivotValue;
      }
    } else {
      region[iPivot] = 0.0;
    }
  }
  //set counts
  regionSparse->setNumElements(numberNonZero);
}
//  updateColumnTransposeL.  Updates part of column transpose (BTRANL)
void CoinFactorization::updateColumnTransposeL(CoinIndexedVector *regionSparse) const
{
  int number = regionSparse->getNumElements();
  if (!numberL_ && !numberDense_) {
    if (sparseArray_ || number < numberRows_)
      return;
  }
  int goSparse;
  // Guess at number at end
  // we may need to rethink on dense
  if (sparseThreshold_ > 0) {
    if (btranAverageAfterL_) {
      double newNumber = number * btranAverageAfterL_;
      if (newNumber < sparseThreshold_)
        goSparse = 2;
      else if (newNumber < sparseThreshold2_)
        goSparse = 1;
      else
        goSparse = 0;
    } else {
      if (number < sparseThreshold_)
        goSparse = 2;
      else
        goSparse = 0;
    }
  } else {
    goSparse = -1;
  }
#ifdef COIN_FACTORIZATION_DENSE_CODE
  if (numberDense_) {
    //take off list
    int lastSparse = numberRows_ - numberDense_;
    int number = regionSparse->getNumElements();
    double *COIN_RESTRICT region = regionSparse->denseVector();
    int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
    int i = 0;
    bool doDense = false;
    if (number <= numberRows_) {
      while (i < number) {
        int iRow = regionIndex[i];
        if (iRow >= lastSparse) {
          doDense = true;
          regionIndex[i] = regionIndex[--number];
        } else {
          i++;
        }
      }
    } else {
      for (i = numberRows_ - 1; i >= lastSparse; i--) {
        if (region[i]) {
          doDense = true;
          // numbers are all wrong - do a scan
          regionSparse->setNumElements(0);
          regionSparse->scan(0, lastSparse, zeroTolerance_);
          number = regionSparse->getNumElements();
          break;
        }
      }
      if (sparseThreshold_)
        goSparse = 0;
      else
        goSparse = -1;
    }
    if (doDense) {
      regionSparse->setNumElements(number);
#if COIN_FACTORIZATION_DENSE_CODE == 1
      char trans = 'T';
      int ione = 1;
      int info;
      COINUTILS_LAPACK_FUNC(dgetrs,DGETRS)(&trans,&numberDense_,&ione,
              denseAreaAddress_,&numberDense_,densePermute_,region+lastSparse,
              &numberDense_,&info,1);
#elif COIN_FACTORIZATION_DENSE_CODE==2
      clapack_dgetrs(CblasColMajor,CblasTrans,numberDense_,1,
                     denseAreaAddress_,numberDense_,densePermute_,
                     region+lastSparse,numberDense_);
#elif COIN_FACTORIZATION_DENSE_CODE==3
      LAPACKE_dgetrs(LAPACK_COL_MAJOR,'T',numberDense_,1,
                     denseAreaAddress_,numberDense_,densePermute_,
                     region+lastSparse,numberDense_);
#endif
      //and scan again
      if (goSparse > 0 || !numberL_)
        regionSparse->scan(lastSparse, numberRows_, zeroTolerance_);
    }
    if (!numberL_) {
      // could be odd combination of sparse/dense
      if (number > numberRows_) {
        regionSparse->setNumElements(0);
        regionSparse->scan(0, numberRows_, zeroTolerance_);
      }
      return;
    }
  }
#endif
  if (goSparse > 0 && regionSparse->getNumElements() > numberRows_)
    goSparse = 0;
  switch (goSparse) {
  case -1: // No row copy
    updateColumnTransposeLDensish(regionSparse);
    break;
  case 0: // densish but by row
    updateColumnTransposeLByRow(regionSparse);
    break;
  case 1: // middling(and by row)
    updateColumnTransposeLSparsish(regionSparse);
    break;
  case 2: // sparse
    updateColumnTransposeLSparse(regionSparse);
    break;
  }
}
#if COIN_ONE_ETA_COPY
/* Combines BtranU and delete elements
   If deleted is NULL then delete elements
   otherwise store where elements are
*/
void CoinFactorization::replaceColumnU(CoinIndexedVector *regionSparse,
  int *COIN_RESTRICT deleted,
  int internalPivotRow)
{
  double *COIN_RESTRICT region = regionSparse->denseVector();
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  double tolerance = zeroTolerance_;
  const int *startColumn = startColumnUArray_;
  const int *indexRow = indexRowUArray_;
  CoinFactorizationDouble *COIN_RESTRICT element = elementUArray_;
  int numberNonZero = 0;
  const int *numberInColumn = numberInColumnArray_;
  //const CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;
  bool deleteNow = true;
  if (deleted) {
    deleteNow = false;
    deleted++;
  }
  int nPut = 0;
  for (int i = std::max(numberSlacks_, internalPivotRow);
       i < numberU_; i++) {
    assert(!region[i]);
    CoinFactorizationDouble pivotValue = 0.0; //region[i]*pivotRegion[i];
    //printf("Epv %g reg %g pr %g\n",
    //   pivotValue,region[i],pivotRegion[i]);
    //pivotValue = region[i];
    for (int j = startColumn[i];
         j < startColumn[i] + numberInColumn[i]; j++) {
      int iRow = indexRow[j];
      CoinFactorizationDouble value = element[j];
      if (iRow != internalPivotRow) {
        pivotValue -= value * region[iRow];
      } else {
        assert(!region[iRow]);
        pivotValue += value;
        if (deleteNow)
          element[j] = 0.0;
        else
          deleted[nPut++] = j;
      }
    }
    if (fabs(pivotValue) > tolerance) {
      regionIndex[numberNonZero++] = i;
      region[i] = pivotValue;
    } else {
      region[i] = 0;
    }
  }
  if (!deleteNow) {
    deleted--;
    deleted[0] = nPut;
  }
  regionSparse->setNumElements(numberNonZero);
}
/* Updates part of column transpose (BTRANU) by column
   assumes index is sorted i.e. region is correct */
void CoinFactorization::updateColumnTransposeUByColumn(CoinIndexedVector *regionSparse,
  int smallestIndex) const
{
  //CoinIndexedVector temp = *regionSparse;
  //updateColumnTransposeUDensish(&temp,smallestIndex);
  double *COIN_RESTRICT region = regionSparse->denseVector();
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  double tolerance = zeroTolerance_;
  const int *startColumn = startColumnUArray_;
  const int *indexRow = indexRowUArray_;
  const CoinFactorizationDouble *COIN_RESTRICT element = elementUArray_;
  int numberNonZero = 0;
  const int *numberInColumn = numberInColumnArray_;
  const CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;

  for (int i = smallestIndex; i < numberSlacks_; i++) {
    double value = region[i];
    if (value) {
      //region[i]=-value;
      regionIndex[numberNonZero] = i;
      if (fabs(value) > tolerance)
        numberNonZero++;
      else
        region[i] = 0.0;
    }
  }
  for (int i = std::max(numberSlacks_, smallestIndex);
       i < numberU_; i++) {
    CoinFactorizationDouble pivotValue = region[i] * pivotRegion[i];
    //printf("pv %g reg %g pr %g\n",
    //   pivotValue,region[i],pivotRegion[i]);
    pivotValue = region[i];
    for (int j = startColumn[i];
         j < startColumn[i] + numberInColumn[i]; j++) {
      int iRow = indexRow[j];
      CoinFactorizationDouble value = element[j];
      pivotValue -= value * region[iRow];
    }
    if (fabs(pivotValue) > tolerance) {
      regionIndex[numberNonZero++] = i;
      region[i] = pivotValue;
    } else {
      region[i] = 0;
    }
  }
  regionSparse->setNumElements(numberNonZero);
  //double * region2 = temp.denseVector();
  //for (i=0;i<maximumRowsExtra_;i++) {
  //assert(fabs(region[i]-region2[i])<1.0e-4);
  //}
}
#endif
// Updates part of column transpose (BTRANR) when dense
void CoinFactorization::updateColumnTransposeRDensish(CoinIndexedVector *regionSparse) const
{
  double *COIN_RESTRICT region = regionSparse->denseVector();

#ifdef COIN_DEBUG
  regionSparse->checkClean();
#endif
  int last = numberRowsExtra_ - 1;

  const int *indexRow = indexRowR_;
  const CoinFactorizationDouble *COIN_RESTRICT element = elementR_;
  const int *startColumn = startColumnRArray_ - numberRows_;
  //move using permute_ (stored in inverse fashion)
  const int *permute = permuteArray_;

  for (int i = last; i >= numberRows_; i--) {
    int putRow = permute[i];
    CoinFactorizationDouble pivotValue = region[i];
    //zero out  old permuted
    region[i] = 0.0;
    if (pivotValue) {
      for (int j = startColumn[i]; j < startColumn[i + 1]; j++) {
        CoinFactorizationDouble value = element[j];
        int iRow = indexRow[j];
        region[iRow] -= value * pivotValue;
      }
      region[putRow] = pivotValue;
      //putRow must have been zero before so put on list ??
      //but can't catch up so will have to do L from end
      //unless we update lookBtran in replaceColumn - yes
    }
  }
}
// Updates part of column transpose (BTRANR) when sparse
void CoinFactorization::updateColumnTransposeRSparse(CoinIndexedVector *regionSparse) const
{
  double *COIN_RESTRICT region = regionSparse->denseVector();
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  int numberNonZero = regionSparse->getNumElements();
  double tolerance = zeroTolerance_;

#ifdef COIN_DEBUG
  regionSparse->checkClean();
#endif
  int last = numberRowsExtra_ - 1;

  const int *indexRow = indexRowR_;
  const CoinFactorizationDouble *COIN_RESTRICT element = elementR_;
  const int *startColumn = startColumnRArray_ - numberRows_;
  //move using permute_ (stored in inverse fashion)
  const int *permute = permuteArray_;

  // we can use sparse_ as temporary array
#if ABOCA_LITE_FACTORIZATION == 0
#define sparseOffset 0
#else
  int sparseOffset = ((regionSparse->capacity() & 0x80000000) != 0) ? sparseOffset_ : 0;
#endif
  int *COIN_RESTRICT spare = sparseArray_ + sparseOffset;
  for (int i = 0; i < numberNonZero; i++) {
    spare[regionIndex[i]] = i;
  }
  // still need to do in correct order (for now)
  for (int i = last; i >= numberRows_; i--) {
    int putRow = permute[i];
    assert(putRow <= i);
    CoinFactorizationDouble pivotValue = region[i];
    //zero out  old permuted
    region[i] = 0.0;
    if (pivotValue) {
      for (int j = startColumn[i]; j < startColumn[i + 1]; j++) {
        CoinFactorizationDouble value = element[j];
        int iRow = indexRow[j];
        CoinFactorizationDouble oldValue = region[iRow];
        CoinFactorizationDouble newValue = oldValue - value * pivotValue;
        if (oldValue) {
          if (newValue)
            region[iRow] = newValue;
          else
            region[iRow] = COIN_INDEXED_REALLY_TINY_ELEMENT;
        } else if (fabs(newValue) > tolerance) {
          region[iRow] = newValue;
          spare[iRow] = numberNonZero;
          regionIndex[numberNonZero++] = iRow;
        }
      }
      region[putRow] = pivotValue;
      // modify list
      int position = spare[i];
      regionIndex[position] = putRow;
      spare[putRow] = position;
    }
  }
  regionSparse->setNumElements(numberNonZero);
}

//  updateColumnTransposeR.  Updates part of column (FTRANR)
void CoinFactorization::updateColumnTransposeR(CoinIndexedVector *regionSparse) const
{
  if (numberRowsExtra_ == numberRows_)
    return;
  int numberNonZero = regionSparse->getNumElements();

  if (numberNonZero) {
    if (numberNonZero < (sparseThreshold_ << 2) || (!numberL_ && sparseArray_)) {
      updateColumnTransposeRSparse(regionSparse);
      if (collectStatistics_)
        btranCountAfterR_ += regionSparse->getNumElements();
    } else {
      updateColumnTransposeRDensish(regionSparse);
      // we have lost indices
      // make sure won't try and go sparse again
      if (collectStatistics_)
        btranCountAfterR_ += std::min((numberNonZero << 1), numberRows_);
      regionSparse->setNumElements(numberRows_ + 1);
    }
  }
}
//  makes a row copy of L
void CoinFactorization::goSparse()
{
  if (!sparseThreshold_) {
#define setSparseSmall 300
#define setSparseMinMedium 500
#define setSparseDivMedium 6
#define setSparseThresh2 2
    if (numberRows_ > setSparseSmall) {
      if (numberRows_ < 10000) {
        sparseThreshold_ =
	  std::min(numberRows_ / setSparseDivMedium, setSparseMinMedium);
        sparseThreshold2_ = numberRows_/setSparseThresh2;
      } else {
        sparseThreshold_ = 1000;
        sparseThreshold2_ = numberRows_ >> 2;
        sparseThreshold_ = 500;
        sparseThreshold2_ = std::max(sparseThreshold_, numberRows_ >> 3);
      }
      //printf("sparseThreshold %d threshold2 %d - numberDense %d\n",
      //     sparseThreshold_,sparseThreshold2_,numberDense_);
    } else {
      sparseThreshold_ = 0;
      sparseThreshold2_ = 0;
    }
  } else {
    if (!sparseThreshold_ && numberRows_ > 400) {
      sparseThreshold_ = std::min((numberRows_ - 300) / 9, 1000);
    }
    sparseThreshold2_ = sparseThreshold_;
  }
  if (!sparseThreshold_)
    return;
  // allow for stack, list, next and char map of mark
  int nRowIndex = (maximumRowsExtra_ + CoinSizeofAsInt(int) - 1) / CoinSizeofAsInt(char);
  int nInBig = static_cast< int >(sizeof(int) / sizeof(int));
  assert(nInBig >= 1);
#if ABOCA_LITE_FACTORIZATION == 0
  sparse_.conditionalNew((2 + nInBig) * maximumRowsExtra_ + nRowIndex);
  // zero out mark
  memset(sparse_.array() + (2 + nInBig) * maximumRowsExtra_,
    0, maximumRowsExtra_ * sizeof(char));
#else
  sparseOffset_ = (2 + nInBig) * maximumRowsExtra_ + nRowIndex;
  sparse_.conditionalNew(2 * sparseOffset_);
  // zero out mark
  memset(sparse_.array() + (2 + nInBig) * maximumRowsExtra_,
    0, maximumRowsExtra_ * sizeof(char));
  memset(sparse_.array() + sparseOffset_ + (2 + nInBig) * maximumRowsExtra_,
    0, maximumRowsExtra_ * sizeof(char));
#endif
  sparseArray_ = sparse_.array();
  elementByRowL_.conditionalDelete();
  indexColumnL_.conditionalDelete();
  startRowL_.conditionalNew(numberRows_ + 1);
  startRowLArray_ = startRowL_.array();
  if (lengthAreaL_) {
    elementByRowL_.conditionalNew(lengthAreaL_);
    indexColumnL_.conditionalNew(lengthAreaL_);
    elementByRowLArray_ = elementByRowL_.array();
    indexColumnLArray_ = indexColumnL_.array();
  }
  // counts
  int *COIN_RESTRICT startRowL = startRowLArray_;
  CoinZeroN(startRowL, numberRows_);
  const int *startColumnL = startColumnLArray_;
  CoinFactorizationDouble *COIN_RESTRICT elementL = elementLArray_;
  const int *indexRowL = indexRowLArray_;
  for (int i = baseL_; i < baseL_ + numberL_; i++) {
    for (int j = startColumnL[i]; j < startColumnL[i + 1]; j++) {
      int iRow = indexRowL[j];
      startRowL[iRow]++;
    }
  }
  // convert count to lasts
  int count = 0;
  for (int i = 0; i < numberRows_; i++) {
    int numberInRow = startRowL[i];
    count += numberInRow;
    startRowL[i] = count;
  }
  startRowL[numberRows_] = count;
  // now insert
  CoinFactorizationDouble *COIN_RESTRICT elementByRowL = elementByRowLArray_;
  int *COIN_RESTRICT indexColumnL = indexColumnLArray_;
  for (int i = baseL_ + numberL_ - 1; i >= baseL_; i--) {
    for (int j = startColumnL[i]; j < startColumnL[i + 1]; j++) {
      int iRow = indexRowL[j];
      int start = startRowL[iRow] - 1;
      startRowL[iRow] = start;
      elementByRowL[start] = elementL[j];
      indexColumnL[start] = i;
    }
  }
}

//  set sparse threshold
void CoinFactorization::sparseThreshold(int value)
{
  if (value > 0 && sparseThreshold_) {
    sparseThreshold_ = value;
    sparseThreshold2_ = sparseThreshold_;
  } else if (!value && sparseThreshold_) {
    // delete copy
    sparseThreshold_ = 0;
    sparseThreshold2_ = 0;
    elementByRowL_.conditionalDelete();
    startRowL_.conditionalDelete();
    indexColumnL_.conditionalDelete();
    sparse_.conditionalDelete();
  } else if (value > 0 && !sparseThreshold_) {
    if (value > 1)
      sparseThreshold_ = value;
    else
      sparseThreshold_ = 0;
    sparseThreshold2_ = sparseThreshold_;
    goSparse();
  }
}
void CoinFactorization::maximumPivots(int value)
{
  if (value > 0) {
    maximumPivots_ = value;
  }
}
void CoinFactorization::messageLevel(int value)
{
  if (value > 0 && value < 16) {
    messageLevel_ = value;
  }
}
void CoinFactorization::pivotTolerance(double value)
{
  if (value > 0.0 && value <= 1.0) {
    //if (value<pivotTolerance_) {
    //printf("reducing pivot tolerance from %g to %g\n",
    //	     pivotTolerance_,value);
    //}
    pivotTolerance_ = value;
  }
}
void CoinFactorization::zeroTolerance(double value)
{
  if (value > 0.0 && value < 1.0) {
    zeroTolerance_ = value;
  }
}
#ifndef COIN_FAST_CODE
void CoinFactorization::slackValue(double value)
{
  if (value >= 0.0) {
    slackValue_ = 1.0;
  } else {
    slackValue_ = -1.0;
  }
}
#endif
// Reset all sparsity etc statistics
void CoinFactorization::resetStatistics()
{

  /// Below are all to collect
  ftranCountInput_ = 0.0;
  ftranCountAfterL_ = 0.0;
  ftranCountAfterR_ = 0.0;
  ftranCountAfterU_ = 0.0;
  btranCountInput_ = 0.0;
  btranCountAfterU_ = 0.0;
  btranCountAfterR_ = 0.0;
  btranCountAfterL_ = 0.0;

  /// We can roll over factorizations
  numberFtranCounts_ = 0;
  numberBtranCounts_ = 0;

  /// While these are averages collected over last
  ftranAverageAfterL_ = 0.0;
  ftranAverageAfterR_ = 0.0;
  ftranAverageAfterU_ = 0.0;
  btranAverageAfterU_ = 0.0;
  btranAverageAfterR_ = 0.0;
  btranAverageAfterL_ = 0.0;
}
/*  getColumnSpaceIterate.  Gets space for one extra U element in Column
    may have to do compression  (returns true)
    also moves existing vector.
    Returns -1 if no memory or where element was put
    Used by replaceRow (turns off R version) */
int CoinFactorization::getColumnSpaceIterate(int iColumn, double value,
  int iRow)
{
  if (numberInColumnPlusArray_) {
    numberInColumnPlus_.conditionalDelete();
    numberInColumnPlusArray_ = NULL;
  }
  int *COIN_RESTRICT numberInRow = numberInRowArray_;
  int *COIN_RESTRICT numberInColumn = numberInColumnArray_;
  int *COIN_RESTRICT nextColumn = nextColumnArray_;
  int *COIN_RESTRICT lastColumn = lastColumnArray_;
  int number = numberInColumn[iColumn];
  int iNext = nextColumn[iColumn];
  int *COIN_RESTRICT startColumnU = startColumnUArray_;
  int *COIN_RESTRICT startRowU = startRowUArray_;
  int space = startColumnU[iNext] - startColumnU[iColumn];
  int put;
  int *COIN_RESTRICT convertRowToColumnU = convertRowToColumnUArray_;
  int *COIN_RESTRICT indexColumnU = indexColumnUArray_;
  CoinFactorizationDouble *COIN_RESTRICT elementU = elementUArray_;
  int *COIN_RESTRICT indexRowU = indexRowUArray_;
  if (space < number + 1) {
    //see if it can go in at end
    if (lengthAreaU_ - startColumnU[maximumColumnsExtra_] < number + 1) {
      //compression
      int jColumn = nextColumn[maximumColumnsExtra_];
      int put = 0;
      while (jColumn != maximumColumnsExtra_) {
        //move
        int get;
        int getEnd;

        get = startColumnU[jColumn];
        getEnd = get + numberInColumn[jColumn];
        startColumnU[jColumn] = put;
        int i;
        for (i = get; i < getEnd; i++) {
          CoinFactorizationDouble value = elementU[i];
          if (value) {
            indexRowU[put] = indexRowU[i];
            elementU[put] = value;
            put++;
          } else {
            numberInColumn[jColumn]--;
          }
        }
        jColumn = nextColumn[jColumn];
      }
      numberCompressions_++;
      startColumnU[maximumColumnsExtra_] = put;
      //space for cross reference
      int *COIN_RESTRICT convertRowToColumn = convertRowToColumnUArray_;
      int j = 0;
      int *COIN_RESTRICT startRow = startRowUArray_;

      int iRow;
      for (iRow = 0; iRow < numberRowsExtra_; iRow++) {
        startRow[iRow] = j;
        j += numberInRow[iRow];
      }
      factorElements_ = j;

      CoinZeroN(numberInRow, numberRowsExtra_);
      int i;
      for (i = 0; i < numberRowsExtra_; i++) {
        int start = startColumnU[i];
        int end = start + numberInColumn[i];

        int j;
        for (j = start; j < end; j++) {
          int iRow = indexRowU[j];
          int iLook = numberInRow[iRow];

          numberInRow[iRow] = iLook + 1;
          int k = startRow[iRow] + iLook;

          indexColumnU[k] = i;
          convertRowToColumn[k] = j;
        }
      }
    }
    // Still may not be room (as iColumn was still in)
    if (lengthAreaU_ - startColumnU[maximumColumnsExtra_] >= number + 1) {
      int next = nextColumn[iColumn];
      int last = lastColumn[iColumn];

      //out
      nextColumn[last] = next;
      lastColumn[next] = last;

      put = startColumnU[maximumColumnsExtra_];
      //in at end
      last = lastColumn[maximumColumnsExtra_];
      nextColumn[last] = iColumn;
      lastColumn[maximumColumnsExtra_] = iColumn;
      lastColumn[iColumn] = last;
      nextColumn[iColumn] = maximumColumnsExtra_;

      //move
      int get = startColumnU[iColumn];
      startColumnU[iColumn] = put;
      int i = 0;
      for (i = 0; i < number; i++) {
        CoinFactorizationDouble value = elementU[get];
        int iRow = indexRowU[get++];
        if (value) {
          elementU[put] = value;
          int n = numberInRow[iRow];
          int start = startRowU[iRow];
          int j;
          for (j = start; j < start + n; j++) {
            if (indexColumnU[j] == iColumn) {
              convertRowToColumnU[j] = put;
              break;
            }
          }
          assert(j < start + n);
          indexRowU[put++] = iRow;
        } else {
          assert(!numberInRow[iRow]);
          numberInColumn[iColumn]--;
        }
      }
      //insert
      int n = numberInRow[iRow];
      int start = startRowU[iRow];
      int j;
      for (j = start; j < start + n; j++) {
        if (indexColumnU[j] == iColumn) {
          convertRowToColumnU[j] = put;
          break;
        }
      }
      assert(j < start + n);
      elementU[put] = value;
      indexRowU[put] = iRow;
      numberInColumn[iColumn]++;
      //add 4 for luck
      startColumnU[maximumColumnsExtra_] = std::min(static_cast< int >(put + 4), lengthAreaU_);
    } else {
      // no room
      put = -1;
    }
  } else {
    // just slot in
    put = startColumnU[iColumn] + numberInColumn[iColumn];
    int n = numberInRow[iRow];
    int start = startRowU[iRow];
    int j;
    for (j = start; j < start + n; j++) {
      if (indexColumnU[j] == iColumn) {
        convertRowToColumnU[j] = put;
        break;
      }
    }
    assert(j < start + n);
    elementU[put] = value;
    indexRowU[put] = iRow;
    numberInColumn[iColumn]++;
  }
  return put;
}
/* Replaces one Row in basis,
   At present assumes just a singleton on row is in basis
   returns 0=OK, 1=Probably OK, 2=singular, 3 no space */
int CoinFactorization::replaceRow(int whichRow, int iNumberInRow,
  const int indicesColumn[], const double elements[])
{
  if (!iNumberInRow)
    return 0;
  int next = nextRowArray_[whichRow];
  int *COIN_RESTRICT numberInRow = numberInRowArray_;
#ifndef NDEBUG
  const int *numberInColumn = numberInColumnArray_;
#endif
  int numberNow = numberInRow[whichRow];
  const int *startRowU = startRowUArray_;
  CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;
  int start = startRowU[whichRow];
  CoinFactorizationDouble *COIN_RESTRICT elementU = elementUArray_;
  int *COIN_RESTRICT convertRowToColumnU = convertRowToColumnUArray_;
  if (numberNow && numberNow < 100) {
    int ind[100];
    CoinMemcpyN(indexColumnUArray_ + start, numberNow, ind);
    int i;
    for (i = 0; i < iNumberInRow; i++) {
      int jColumn = indicesColumn[i];
      int k;
      for (k = 0; k < numberNow; k++) {
        if (ind[k] == jColumn) {
          ind[k] = -1;
          break;
        }
      }
      if (k == numberNow) {
        printf("new column %d not in current\n", jColumn);
      } else {
        k = convertRowToColumnU[start + k];
        CoinFactorizationDouble oldValue = elementU[k];
        CoinFactorizationDouble newValue = elements[i] * pivotRegion[jColumn];
        if (fabs(oldValue - newValue) > 1.0e-3)
          printf("column %d, old value %g new %g (el %g, piv %g)\n", jColumn, oldValue,
            newValue, elements[i], pivotRegion[jColumn]);
      }
    }
    for (i = 0; i < numberNow; i++) {
      if (ind[i] >= 0)
        printf("current column %d not in new\n", ind[i]);
    }
    assert(numberNow == iNumberInRow);
  }
  assert(numberInColumn[whichRow] == 0);
  assert(pivotRegion[whichRow] == 1.0);
  int space;
  int returnCode = 0;

  space = startRowU[next] - (start + iNumberInRow);
  if (space < 0) {
    if (!getRowSpaceIterate(whichRow, iNumberInRow))
      returnCode = 3;
  }
  //return 0;
  if (!returnCode) {
    int *COIN_RESTRICT indexColumnU = indexColumnUArray_;
    numberInRow[whichRow] = iNumberInRow;
    start = startRowU[whichRow];
    int i;
    for (i = 0; i < iNumberInRow; i++) {
      int iColumn = indicesColumn[i];
      indexColumnU[start + i] = iColumn;
      assert(iColumn > whichRow);
      CoinFactorizationDouble value = elements[i] * pivotRegion[iColumn];
      int iWhere = getColumnSpaceIterate(iColumn, value, whichRow);
      if (iWhere >= 0) {
        convertRowToColumnU[start + i] = iWhere;
      } else {
        returnCode = 3;
        break;
      }
    }
  }
  return returnCode;
}
// Takes out all entries for given rows
void CoinFactorization::emptyRows(int numberToEmpty, const int which[])
{
  int i;
  int *COIN_RESTRICT delRow = new int[maximumRowsExtra_];
  int *COIN_RESTRICT indexRowU = indexRowUArray_;
#ifndef NDEBUG
  CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;
#endif
  for (i = 0; i < maximumRowsExtra_; i++)
    delRow[i] = 0;
  int *COIN_RESTRICT numberInRow = numberInRowArray_;
  int *COIN_RESTRICT numberInColumn = numberInColumnArray_;
  CoinFactorizationDouble *COIN_RESTRICT elementU = elementUArray_;
  const int *startColumnU = startColumnUArray_;
  for (i = 0; i < numberToEmpty; i++) {
    int iRow = which[i];
    delRow[iRow] = 1;
    assert(numberInColumn[iRow] == 0);
    assert(pivotRegion[iRow] == 1.0);
    numberInRow[iRow] = 0;
  }
  for (i = 0; i < numberU_; i++) {
    int k;
    int j = startColumnU[i];
    for (k = startColumnU[i]; k < startColumnU[i] + numberInColumn[i]; k++) {
      int iRow = indexRowU[k];
      if (!delRow[iRow]) {
        indexRowU[j] = indexRowU[k];
        elementU[j++] = elementU[k];
      }
    }
    numberInColumn[i] = j - startColumnU[i];
  }
  delete[] delRow;
  //space for cross reference
  int *COIN_RESTRICT convertRowToColumn = convertRowToColumnUArray_;
  int j = 0;
  int *COIN_RESTRICT startRow = startRowUArray_;

  int iRow;
  for (iRow = 0; iRow < numberRows_; iRow++) {
    startRow[iRow] = j;
    j += numberInRow[iRow];
  }
  factorElements_ = j;

  CoinZeroN(numberInRow, numberRows_);

  int *COIN_RESTRICT indexColumnU = indexColumnUArray_;
  for (i = 0; i < numberRows_; i++) {
    int start = startColumnU[i];
    int end = start + numberInColumn[i];

    int j;
    for (j = start; j < end; j++) {
      int iRow = indexRowU[j];
      int iLook = numberInRow[iRow];

      numberInRow[iRow] = iLook + 1;
      int k = startRow[iRow] + iLook;

      indexColumnU[k] = i;
      convertRowToColumn[k] = j;
    }
  }
}
// Updates part of column PFI (FTRAN)
void CoinFactorization::updateColumnPFI(CoinIndexedVector *regionSparse) const
{
  double *COIN_RESTRICT region = regionSparse->denseVector();
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  double tolerance = zeroTolerance_;
  const int *startColumn = startColumnUArray_ + numberRows_;
  const int *indexRow = indexRowUArray_;
  const CoinFactorizationDouble *COIN_RESTRICT element = elementUArray_;
  int numberNonZero = regionSparse->getNumElements();
  int i;
#ifdef COIN_DEBUG
  for (i = 0; i < numberNonZero; i++) {
    int iRow = regionIndex[i];
    assert(iRow >= 0 && iRow < numberRows_);
    assert(region[iRow]);
  }
#endif
  const CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_ + numberRows_;
  const int *pivotColumn = pivotColumnArray_ + numberRows_;

  for (i = 0; i < numberPivots_; i++) {
    int pivotRow = pivotColumn[i];
    CoinFactorizationDouble pivotValue = region[pivotRow];
    if (pivotValue) {
      if (fabs(pivotValue) > tolerance) {
        for (int j = startColumn[i]; j < startColumn[i + 1]; j++) {
          int iRow = indexRow[j];
          CoinFactorizationDouble oldValue = region[iRow];
          CoinFactorizationDouble value = oldValue - pivotValue * element[j];
          if (!oldValue) {
            if (fabs(value) > tolerance) {
              region[iRow] = value;
              regionIndex[numberNonZero++] = iRow;
            }
          } else {
            if (fabs(value) > tolerance) {
              region[iRow] = value;
            } else {
              region[iRow] = COIN_INDEXED_REALLY_TINY_ELEMENT;
            }
          }
        }
        pivotValue *= pivotRegion[i];
        region[pivotRow] = pivotValue;
      } else if (pivotValue) {
        region[pivotRow] = COIN_INDEXED_REALLY_TINY_ELEMENT;
      }
    }
  }
  regionSparse->setNumElements(numberNonZero);
}
// Updates part of column transpose PFI (BTRAN),

void CoinFactorization::updateColumnTransposePFI(CoinIndexedVector *regionSparse) const
{
  double *COIN_RESTRICT region = regionSparse->denseVector();
  int numberNonZero = regionSparse->getNumElements();
  int *COIN_RESTRICT index = regionSparse->getIndices();
  int i;
#ifdef COIN_DEBUG
  for (i = 0; i < numberNonZero; i++) {
    int iRow = index[i];
    assert(iRow >= 0 && iRow < numberRows_);
    assert(region[iRow]);
  }
#endif
  const int *pivotColumn = pivotColumnArray_ + numberRows_;
  const CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_ + numberRows_;
  double tolerance = zeroTolerance_;

  const int *startColumn = startColumnUArray_ + numberRows_;
  const int *indexRow = indexRowUArray_;
  const CoinFactorizationDouble *COIN_RESTRICT element = elementUArray_;

  for (i = numberPivots_ - 1; i >= 0; i--) {
    int pivotRow = pivotColumn[i];
    CoinFactorizationDouble pivotValue = region[pivotRow] * pivotRegion[i];
    for (int j = startColumn[i]; j < startColumn[i + 1]; j++) {
      int iRow = indexRow[j];
      CoinFactorizationDouble value = element[j];
      pivotValue -= value * region[iRow];
    }
    //pivotValue *= pivotRegion[i];
    if (fabs(pivotValue) > tolerance) {
      if (!region[pivotRow])
        index[numberNonZero++] = pivotRow;
      region[pivotRow] = pivotValue;
    } else {
      if (region[pivotRow])
        region[pivotRow] = COIN_INDEXED_REALLY_TINY_ELEMENT;
    }
  }
  //set counts
  regionSparse->setNumElements(numberNonZero);
}
/* Replaces one Column to basis for PFI
   returns 0=OK, 1=Probably OK, 2=singular, 3=no room
   Not sure what checkBeforeModifying means for PFI.
*/
int CoinFactorization::replaceColumnPFI(CoinIndexedVector *regionSparse,
  int pivotRow,
  double alpha)
{
  int *COIN_RESTRICT startColumn = startColumnUArray_ + numberRows_;
  int *COIN_RESTRICT indexRow = indexRowUArray_;
  CoinFactorizationDouble *COIN_RESTRICT element = elementUArray_;
  CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_ + numberRows_;
  // This has incoming column
  const double *region = regionSparse->denseVector();
  const int *index = regionSparse->getIndices();
  int numberNonZero = regionSparse->getNumElements();

  int iColumn = numberPivots_;

  if (!iColumn)
    startColumn[0] = startColumn[maximumColumnsExtra_];
  int start = startColumn[iColumn];

  //return at once if too many iterations
  if (numberPivots_ >= maximumPivots_) {
    return 5;
  }
  if (lengthAreaU_ - (start + numberNonZero) < 0) {
    return 3;
  }

  int i;
  if (numberPivots_) {
    if (fabs(alpha) < 1.0e-5) {
      if (fabs(alpha) < 1.0e-7)
        return 2;
      else
        return 1;
    }
  } else {
    if (fabs(alpha) < 1.0e-8)
      return 2;
  }
  CoinFactorizationDouble pivotValue = 1.0 / alpha;
  pivotRegion[iColumn] = pivotValue;
  double tolerance = zeroTolerance_;
  const int *pivotColumn = pivotColumnArray_;
  // Operations done before permute back
  if (regionSparse->packedMode()) {
    for (i = 0; i < numberNonZero; i++) {
      int iRow = index[i];
      if (iRow != pivotRow) {
        if (fabs(region[i]) > tolerance) {
          indexRow[start] = pivotColumn[iRow];
          element[start++] = region[i] * pivotValue;
        }
      }
    }
  } else {
    for (i = 0; i < numberNonZero; i++) {
      int iRow = index[i];
      if (iRow != pivotRow) {
        if (fabs(region[iRow]) > tolerance) {
          indexRow[start] = pivotColumn[iRow];
          element[start++] = region[iRow] * pivotValue;
        }
      }
    }
  }
  numberPivots_++;
  numberNonZero = start - startColumn[iColumn];
  startColumn[numberPivots_] = start;
  totalElements_ += numberNonZero;
  int *COIN_RESTRICT pivotColumn2 = pivotColumn_.array() + numberRows_;
  pivotColumn2[iColumn] = pivotColumn[pivotRow];
  return 0;
}
//  =
CoinFactorization &CoinFactorization::operator=(const CoinFactorization &other)
{
  if (this != &other) {
    gutsOfDestructor(2);
    gutsOfInitialize(3);
    persistenceFlag_ = other.persistenceFlag_;
    gutsOfCopy(other);
    setupPointers();
  }
  return *this;
}
void CoinFactorization::gutsOfCopy(const CoinFactorization &other)
{
  int lengthU = other.lengthAreaU_ + EXTRA_U_SPACE;
  elementU_.allocate(other.elementU_, lengthU * CoinSizeofAsInt(CoinFactorizationDouble));
  indexRowU_.allocate(other.indexRowU_, lengthU * CoinSizeofAsInt(int));
  elementL_.allocate(other.elementL_, other.lengthAreaL_ * CoinSizeofAsInt(CoinFactorizationDouble));
  indexRowL_.allocate(other.indexRowL_, other.lengthAreaL_ * CoinSizeofAsInt(int));
  startColumnL_.allocate(other.startColumnL_, (other.numberRows_ + 1) * CoinSizeofAsInt(int));
  int extraSpace;
  if (other.numberInColumnPlus_.array()) {
    extraSpace = other.maximumPivots_ + 1 + other.maximumColumnsExtra_ + 1;
  } else {
    extraSpace = other.maximumPivots_ + 1;
  }
  startColumnR_.allocate(other.startColumnR_, extraSpace * CoinSizeofAsInt(int));
  pivotRegion_.allocate(other.pivotRegion_, (other.maximumRowsExtra_ + 1) * CoinSizeofAsInt(CoinFactorizationDouble));
  permuteBack_.allocate(other.permuteBack_, (other.maximumRowsExtra_ + 1) * CoinSizeofAsInt(int));
  permute_.allocate(other.permute_, (other.maximumRowsExtra_ + 1) * CoinSizeofAsInt(int));
  pivotColumnBack_.allocate(other.pivotColumnBack_, (other.maximumRowsExtra_ + 1) * CoinSizeofAsInt(int));
  firstCount_.allocate(other.firstCount_, (other.maximumRowsExtra_ + 1) * CoinSizeofAsInt(int));
  startColumnU_.allocate(other.startColumnU_, (other.maximumColumnsExtra_ + 1) * CoinSizeofAsInt(int));
  numberInColumn_.allocate(other.numberInColumn_, (other.maximumColumnsExtra_ + 1) * CoinSizeofAsInt(int));
  pivotColumn_.allocate(other.pivotColumn_, (other.maximumColumnsExtra_ + 1) * CoinSizeofAsInt(int));
  nextColumn_.allocate(other.nextColumn_, (other.maximumColumnsExtra_ + 1) * CoinSizeofAsInt(int));
  lastColumn_.allocate(other.lastColumn_, (other.maximumColumnsExtra_ + 1) * CoinSizeofAsInt(int));
  indexColumnU_.allocate(other.indexColumnU_, lengthU * CoinSizeofAsInt(int));
  nextRow_.allocate(other.nextRow_, (other.maximumRowsExtra_ + 1) * CoinSizeofAsInt(int));
  lastRow_.allocate(other.lastRow_, (other.maximumRowsExtra_ + 1) * CoinSizeofAsInt(int));
  const int *convertUOther = other.convertRowToColumnU_.array();
#if COIN_ONE_ETA_COPY
  if (convertUOther) {
#endif
    convertRowToColumnU_.allocate(other.convertRowToColumnU_, lengthU * CoinSizeofAsInt(int));
    startRowU_.allocate(other.startRowU_, (other.maximumRowsExtra_ + 1) * CoinSizeofAsInt(int));
    numberInRow_.allocate(other.numberInRow_, (other.maximumRowsExtra_ + 1) * CoinSizeofAsInt(int));
#if COIN_ONE_ETA_COPY
  }
#endif
  if (other.sparseThreshold_) {
    elementByRowL_.allocate(other.elementByRowL_, other.lengthAreaL_);
    indexColumnL_.allocate(other.indexColumnL_, other.lengthAreaL_);
    startRowL_.allocate(other.startRowL_, other.numberRows_ + 1);
  }
  numberTrials_ = other.numberTrials_;
  biggerDimension_ = other.biggerDimension_;
  relaxCheck_ = other.relaxCheck_;
  numberSlacks_ = other.numberSlacks_;
  numberU_ = other.numberU_;
  maximumU_ = other.maximumU_;
  lengthU_ = other.lengthU_;
  lengthAreaU_ = other.lengthAreaU_;
  numberL_ = other.numberL_;
  baseL_ = other.baseL_;
  lengthL_ = other.lengthL_;
  lengthAreaL_ = other.lengthAreaL_;
  numberR_ = other.numberR_;
  lengthR_ = other.lengthR_;
  lengthAreaR_ = other.lengthAreaR_;
  pivotTolerance_ = other.pivotTolerance_;
  zeroTolerance_ = other.zeroTolerance_;
#ifndef COIN_FAST_CODE
  slackValue_ = other.slackValue_;
#endif
  areaFactor_ = other.areaFactor_;
  numberRows_ = other.numberRows_;
  numberRowsExtra_ = other.numberRowsExtra_;
  maximumRowsExtra_ = other.maximumRowsExtra_;
  numberColumns_ = other.numberColumns_;
  numberColumnsExtra_ = other.numberColumnsExtra_;
  maximumColumnsExtra_ = other.maximumColumnsExtra_;
  maximumPivots_ = other.maximumPivots_;
  numberGoodU_ = other.numberGoodU_;
  numberGoodL_ = other.numberGoodL_;
  numberPivots_ = other.numberPivots_;
  messageLevel_ = other.messageLevel_;
  totalElements_ = other.totalElements_;
  factorElements_ = other.factorElements_;
  status_ = other.status_;
  doForrestTomlin_ = other.doForrestTomlin_;
  ftranCountInput_ = other.ftranCountInput_;
  ftranCountAfterL_ = other.ftranCountAfterL_;
  ftranCountAfterR_ = other.ftranCountAfterR_;
  ftranCountAfterU_ = other.ftranCountAfterU_;
  btranCountInput_ = other.btranCountInput_;
  btranCountAfterU_ = other.btranCountAfterU_;
  btranCountAfterR_ = other.btranCountAfterR_;
  btranCountAfterL_ = other.btranCountAfterL_;
  numberFtranCounts_ = other.numberFtranCounts_;
  numberBtranCounts_ = other.numberBtranCounts_;
  ftranAverageAfterL_ = other.ftranAverageAfterL_;
  ftranAverageAfterR_ = other.ftranAverageAfterR_;
  ftranAverageAfterU_ = other.ftranAverageAfterU_;
  btranAverageAfterU_ = other.btranAverageAfterU_;
  btranAverageAfterR_ = other.btranAverageAfterR_;
  btranAverageAfterL_ = other.btranAverageAfterL_;
  biasLU_ = other.biasLU_;
  sparseThreshold_ = other.sparseThreshold_;
  sparseThreshold2_ = other.sparseThreshold2_;
  int space = lengthAreaL_ - lengthL_;

  numberDense_ = other.numberDense_;
  denseThreshold_ = other.denseThreshold_;
  if (numberDense_) {
    denseArea_ = new double[numberDense_ * numberDense_];
    denseAreaAddress_ = denseArea_;
    CoinMemcpyN(other.denseAreaAddress_,
      numberDense_ * numberDense_, denseAreaAddress_);
    densePermute_ = new int[numberDense_];
    CoinMemcpyN(other.densePermute_,
      numberDense_, densePermute_);
  }

  lengthAreaR_ = space;
  elementR_ = elementL_.array() + lengthL_;
  indexRowR_ = indexRowL_.array() + lengthL_;
  workArea_ = other.workArea_;
  workArea2_ = other.workArea2_;
  //now copy everything
  //assuming numberRowsExtra==numberColumnsExtra
  if (numberRowsExtra_) {
    if (convertUOther) {
      CoinMemcpyN(other.startRowU_.array(), numberRowsExtra_ + 1, startRowU_.array());
      CoinMemcpyN(other.numberInRow_.array(), numberRowsExtra_ + 1, numberInRow_.array());
      startRowU_.array()[maximumRowsExtra_] = other.startRowU_.array()[maximumRowsExtra_];
    }
    CoinMemcpyN(other.pivotRegion_.array(), numberRowsExtra_, pivotRegion_.array());
    CoinMemcpyN(other.permuteBack_.array(), numberRowsExtra_ + 1, permuteBack_.array());
    CoinMemcpyN(other.permute_.array(), numberRowsExtra_ + 1, permute_.array());
    CoinMemcpyN(other.pivotColumnBack_.array(), numberRowsExtra_ + 1, pivotColumnBack_.array());
    CoinMemcpyN(other.firstCount_.array(), numberRowsExtra_ + 1, firstCount_.array());
    CoinMemcpyN(other.startColumnU_.array(), numberRowsExtra_ + 1, startColumnU_.array());
    CoinMemcpyN(other.numberInColumn_.array(), numberRowsExtra_ + 1, numberInColumn_.array());
    CoinMemcpyN(other.pivotColumn_.array(), numberRowsExtra_ + 1, pivotColumn_.array());
    CoinMemcpyN(other.nextColumn_.array(), numberRowsExtra_ + 1, nextColumn_.array());
    CoinMemcpyN(other.lastColumn_.array(), numberRowsExtra_ + 1, lastColumn_.array());
    CoinMemcpyN(other.startColumnR_.array(), numberRowsExtra_ - numberColumns_ + 1,
      startColumnR_.array());
    //extra one at end
    startColumnU_.array()[maximumColumnsExtra_] = other.startColumnU_.array()[maximumColumnsExtra_];
    nextColumn_.array()[maximumColumnsExtra_] = other.nextColumn_.array()[maximumColumnsExtra_];
    lastColumn_.array()[maximumColumnsExtra_] = other.lastColumn_.array()[maximumColumnsExtra_];
    CoinMemcpyN(other.nextRow_.array(), numberRowsExtra_ + 1, nextRow_.array());
    CoinMemcpyN(other.lastRow_.array(), numberRowsExtra_ + 1, lastRow_.array());
    nextRow_.array()[maximumRowsExtra_] = other.nextRow_.array()[maximumRowsExtra_];
    lastRow_.array()[maximumRowsExtra_] = other.lastRow_.array()[maximumRowsExtra_];
  }
  CoinMemcpyN(other.elementR_, lengthR_, elementR_);
  CoinMemcpyN(other.indexRowR_, lengthR_, indexRowR_);
  //row and column copies of U
  /* as elements of U may have been zeroed but column counts zero
     copy all elements */
  const int *startColumnU = startColumnU_.array();
  const int *numberInColumn = numberInColumn_.array();
#ifndef NDEBUG
  int maxU = 0;
  for (int iRow = 0; iRow < numberRowsExtra_; iRow++) {
    int start = startColumnU[iRow];
    int numberIn = numberInColumn[iRow];
    maxU = std::max(maxU, start + numberIn);
  }
  assert(maximumU_ >= maxU);
#endif
  CoinMemcpyN(other.elementU_.array(), maximumU_, elementU_.array());
#if COIN_ONE_ETA_COPY
  if (convertUOther) {
#endif
    const int *indexColumnUOther = other.indexColumnU_.array();
    int *COIN_RESTRICT convertU = convertRowToColumnU_.array();
    int *COIN_RESTRICT indexColumnU = indexColumnU_.array();
    const int *startRowU = startRowU_.array();
    const int *numberInRow = numberInRow_.array();
    for (int iRow = 0; iRow < numberRowsExtra_; iRow++) {
      //row
      int start = startRowU[iRow];
      int numberIn = numberInRow[iRow];

      CoinMemcpyN(indexColumnUOther + start, numberIn, indexColumnU + start);
      CoinMemcpyN(convertUOther + start, numberIn, convertU + start);
    }
#if COIN_ONE_ETA_COPY
  }
#endif
  const int *indexRowUOther = other.indexRowU_.array();
  int *COIN_RESTRICT indexRowU = indexRowU_.array();
  for (int iRow = 0; iRow < numberRowsExtra_; iRow++) {
    //column
    int start = startColumnU[iRow];
    int numberIn = numberInColumn[iRow];
    CoinMemcpyN(indexRowUOther + start, numberIn, indexRowU + start);
  }
  // L is contiguous
  if (numberRows_)
    CoinMemcpyN(other.startColumnL_.array(), numberRows_ + 1, startColumnL_.array());
  CoinMemcpyN(other.elementL_.array(), lengthL_, elementL_.array());
  CoinMemcpyN(other.indexRowL_.array(), lengthL_, indexRowL_.array());
  setupPointers();
  if (other.sparseThreshold_) {
    goSparse();
  }
}
// See if worth going sparse
void CoinFactorization::checkSparse()
{
  // See if worth going sparse and when
  if (numberFtranCounts_ > 100) {
    ftranCountInput_ = std::max(ftranCountInput_, 1.0);
    ftranAverageAfterL_ = std::max(ftranCountAfterL_ / ftranCountInput_, 1.0);
    ftranAverageAfterR_ = std::max(ftranCountAfterR_ / ftranCountAfterL_, 1.0);
    ftranAverageAfterU_ = std::max(ftranCountAfterU_ / ftranCountAfterR_, 1.0);
    if (btranCountInput_ && btranCountAfterU_ && btranCountAfterR_) {
      btranAverageAfterU_ = std::max(btranCountAfterU_ / btranCountInput_, 1.0);
      btranAverageAfterR_ = std::max(btranCountAfterR_ / btranCountAfterU_, 1.0);
      btranAverageAfterL_ = std::max(btranCountAfterL_ / btranCountAfterR_, 1.0);
    } else {
      // we have not done any useful btrans (values pass?)
      btranAverageAfterU_ = 1.0;
      btranAverageAfterR_ = 1.0;
      btranAverageAfterL_ = 1.0;
    }
  }
  // scale back

  ftranCountInput_ *= 0.8;
  ftranCountAfterL_ *= 0.8;
  ftranCountAfterR_ *= 0.8;
  ftranCountAfterU_ *= 0.8;
  btranCountInput_ *= 0.8;
  btranCountAfterU_ *= 0.8;
  btranCountAfterR_ *= 0.8;
  btranCountAfterL_ *= 0.8;
}
// Condition number - product of pivots after factorization
double
CoinFactorization::conditionNumber() const
{
  double condition = 1.0;
  const CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;
  for (int i = 0; i < numberRows_; i++) {
    condition *= pivotRegion[i];
  }
  condition = std::max(fabs(condition), 1.0e-50);
  return 1.0 / condition;
}
#ifdef ABC_USE_COIN_FACTORIZATION
/* Checks if can replace one Column to basis,
   returns update alpha
   Fills in region for use later
   partial update already in U */
double
CoinFactorization::checkReplacePart1(CoinIndexedVector *regionSparse,
  int pivotRow)
{
  assert(numberU_ <= numberRowsExtra_);
  int *COIN_RESTRICT startColumnU = startColumnUArray_;
  CoinFactorizationDouble *COIN_RESTRICT element;
  int *COIN_RESTRICT numberInRow = numberInRowArray_;
  int *COIN_RESTRICT numberInColumn = numberInColumnArray_;
  int realPivotRow;
  realPivotRow = pivotColumnArray_[pivotRow];
  //zeroed out region
  double *COIN_RESTRICT region = regionSparse->denseVector();

  element = elementUArray_;
  //get entries in row (pivot not stored)
  int numberNonZero = 0;
  int *COIN_RESTRICT indexColumn = indexColumnUArray_;
  int *COIN_RESTRICT convertRowToColumn = convertRowToColumnUArray_;
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  int *COIN_RESTRICT startRow = startRowUArray_;
  int start = 0;
  int end = 0;

  assert(convertRowToColumn);
  start = startRow[realPivotRow];
  end = start + numberInRow[realPivotRow];

  int smallestIndex = numberRowsExtra_;
  for (int i = start; i < end; i++) {
    int iColumn = indexColumn[i];
    smallestIndex = std::min(smallestIndex, iColumn);
    int j = convertRowToColumn[i];

    region[iColumn] = element[j];
    regionIndex[numberNonZero++] = iColumn;
  }
  //do BTRAN - finding first one to use
  regionSparse->setNumElements(numberNonZero);
  updateColumnTransposeU(regionSparse, smallestIndex);
  numberNonZero = regionSparse->getNumElements();
  CoinFactorizationDouble saveFromU = 0.0;

  int startU = startColumnU[numberColumnsExtra_];
  int *COIN_RESTRICT indexU = &indexRowUArray_[startU];
  CoinFactorizationDouble *COIN_RESTRICT elementU = &elementUArray_[startU];
  double tolerance = zeroTolerance_;
  int number = numberInColumn[numberColumnsExtra_];

  for (int i = 0; i < number; i++) {
    int iRow = indexU[i];
    if (fabs(elementU[i]) > tolerance) {
      if (iRow != realPivotRow) {
        saveFromU -= elementU[i] * region[iRow];
      } else {
        saveFromU += elementU[i];
      }
    }
  }
  return saveFromU;
}
/* Checks if can replace one Column to basis,
   returns update alpha
   Fills in region for use later
   partial update in vector */
double
CoinFactorization::checkReplacePart1(CoinIndexedVector *regionSparse,
  CoinIndexedVector *partialUpdate,
  int pivotRow)
{
  CoinFactorizationDouble *COIN_RESTRICT element;
  int *COIN_RESTRICT numberInRow = numberInRowArray_;
  int realPivotRow;
  realPivotRow = pivotColumnArray_[pivotRow];
  //zeroed out region
  double *COIN_RESTRICT region = regionSparse->denseVector();

  element = elementUArray_;
  //get entries in row (pivot not stored)
  int numberNonZero = 0;
  int *COIN_RESTRICT indexColumn = indexColumnUArray_;
  int *COIN_RESTRICT convertRowToColumn = convertRowToColumnUArray_;
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  int *COIN_RESTRICT startRow = startRowUArray_;
  int start = 0;
  int end = 0;

  assert(convertRowToColumn);
  start = startRow[realPivotRow];
  end = start + numberInRow[realPivotRow];

  int smallestIndex = numberRowsExtra_;
  for (int i = start; i < end; i++) {
    int iColumn = indexColumn[i];
    smallestIndex = std::min(smallestIndex, iColumn);
    int j = convertRowToColumn[i];

    region[iColumn] = element[j];
    regionIndex[numberNonZero++] = iColumn;
  }
  //do BTRAN - finding first one to use
  regionSparse->setNumElements(numberNonZero);
  updateColumnTransposeU(regionSparse, smallestIndex);
  numberNonZero = regionSparse->getNumElements();
  CoinFactorizationDouble saveFromU = 0.0;
  double tolerance = zeroTolerance_;
  int *COIN_RESTRICT indexU2 = partialUpdate->getIndices();
  CoinFactorizationDouble *COIN_RESTRICT elementU2 = partialUpdate->denseVector();
  int numberInColumnU2 = partialUpdate->getNumElements();
  for (int i = 0; i < numberInColumnU2; i++) {
    int iRow = indexU2[i];
    if (fabs(elementU2[i]) > tolerance) {
      if (iRow != realPivotRow) {
        saveFromU -= elementU2[i] * region[iRow];
      } else {
        saveFromU += elementU2[i];
      }
    }
  }
  if (lengthU_ + numberInColumnU2 >= lengthAreaU_) {
    //not enough room
    saveFromU = 0.0;
  }
  if (lengthR_ + numberNonZero >= lengthAreaR_) {
    //not enough room
    saveFromU = 0.0;
  }
  return saveFromU;
}
/* Checks if can replace one Column in basis,
   returns 0=OK, 1=Probably OK, 2=singular, 3=no room, 5 max pivots */
int CoinFactorization::checkReplacePart2(int pivotRow,
  double btranAlpha,
  double ftranAlpha,
  double ftAlpha,
  double acceptablePivot)
{
  //return at once if too many iterations
  if (numberColumnsExtra_ >= maximumColumnsExtra_) {
    return 5;
  }
  if (lengthR_ + numberRows_ >= lengthAreaR_) {
    //not enough room
    return 3;
  }
  pivotRow = pivotColumnArray_[pivotRow];
  const CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;
  CoinFactorizationDouble oldPivot = pivotRegion[pivotRow];
  // for accuracy check
  CoinFactorizationDouble pivotCheck = ftranAlpha / oldPivot;
  //check accuracy
  int status = checkPivot(ftAlpha, pivotCheck);
  if (status == 1 && !numberPivots_) {
    printf("check status not ok\n");
    status = 2;
  }
  return status;
}
/* Replaces one Column to basis,
   partial update already in U */
void CoinFactorization::replaceColumnPart3(CoinIndexedVector *regionSparse,
  int pivotRow,
  double alpha)
{
  assert(numberU_ <= numberRowsExtra_);
  int *COIN_RESTRICT startColumnU = startColumnUArray_;
  int *COIN_RESTRICT startColumn;
  int *COIN_RESTRICT indexRow;
  CoinFactorizationDouble *COIN_RESTRICT element;
  int *COIN_RESTRICT numberInRow = numberInRowArray_;
  int *COIN_RESTRICT numberInColumn = numberInColumnArray_;
  int *COIN_RESTRICT numberInColumnPlus = numberInColumnPlusArray_;
  int realPivotRow;
  realPivotRow = pivotColumnArray_[pivotRow];
  //Filled in region
  double *COIN_RESTRICT region = regionSparse->denseVector();

  element = elementUArray_;
  //take out old pivot column
  totalElements_ -= numberInColumn[realPivotRow];
  CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;
  pivotRegion[realPivotRow] = 0.0;

  int saveEnd = startColumnU[realPivotRow]
    + numberInColumn[realPivotRow];
#ifdef CLP_FACTORIZATION_INSTRUMENT
  currentTakeoutU += numberInColumn[realPivotRow];
  currentTakeoutU += numberInRow[realPivotRow];
#endif
  // not necessary at present - but take no chances for future
  numberInColumn[realPivotRow] = 0;
  //get entries in row (pivot not stored)
  int numberNonZero = 0;
  int *COIN_RESTRICT indexColumn = indexColumnUArray_;
  int *COIN_RESTRICT convertRowToColumn = convertRowToColumnUArray_;
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  int *COIN_RESTRICT startRow = startRowUArray_;
  int start = 0;
  int end = 0;
  start = startRow[realPivotRow];
  end = start + numberInRow[realPivotRow];

  for (int i = start; i < end; i++) {
    int j = convertRowToColumn[i];
    element[j] = 0.0;
  }
  numberNonZero = regionSparse->getNumElements();
  int startU = startColumnU[numberColumnsExtra_];
  int *COIN_RESTRICT indexU = &indexRowUArray_[startU];
  CoinFactorizationDouble *COIN_RESTRICT elementU = &elementUArray_[startU];
  // Now zero out column of U
  //take out old pivot column
  for (int i = startColumnU[realPivotRow]; i < saveEnd; i++) {
    element[i] = 0.0;
  }
  //zero out pivot Row (before or after?)
  //add to R
  startColumn = startColumnRArray_;
  indexRow = indexRowR_;
  element = elementR_;
  int l = lengthR_;
  int number = numberR_;

  startColumn[number] = l; //for luck and first time
  number++;
  startColumn[number] = l + numberNonZero;
  numberR_ = number;
  lengthR_ = l + numberNonZero;
  totalElements_ += numberNonZero;
  assert(lengthR_ < lengthAreaR_);
  for (int i = 0; i < numberNonZero; i++) {
    int iRow = regionIndex[i];
    indexRow[l] = iRow;
    element[l] = region[iRow];
    l++;
  }
  int *COIN_RESTRICT nextRow;
  int *COIN_RESTRICT lastRow;
  int next;
  int last;
  //take out row
  nextRow = nextRowArray_;
  lastRow = lastRowArray_;
  next = nextRow[realPivotRow];
  last = lastRow[realPivotRow];

  nextRow[last] = next;
  lastRow[next] = last;
  numberInRow[realPivotRow] = 0;
  //do permute
  permuteArray_[numberRowsExtra_] = realPivotRow;
  // and other way
  permuteBackArray_[realPivotRow] = numberRowsExtra_;
  permuteBackArray_[numberRowsExtra_] = -1;
  ;
  //and for safety
  permuteArray_[numberRowsExtra_ + 1] = 0;

  pivotColumnArray_[pivotRow] = numberRowsExtra_;
  pivotColumnBack()[numberRowsExtra_] = pivotRow;
  startColumn = startColumnU;
  indexRow = indexRowUArray_;
  element = elementUArray_;

  numberU_++;
  number = numberInColumn[numberColumnsExtra_];

  totalElements_ += number;
  lengthU_ += number;

  double saveFromU = 0.0;

  //put in pivot
  //add row counts

  for (int i = 0; i < number; i++) {
    int iRow = indexU[i];
    if (iRow != realPivotRow) {
      int next = nextRow[iRow];
      int iNumberInRow = numberInRow[iRow];
      int space;
      int put = startRow[iRow] + iNumberInRow;

      space = startRow[next] - put;
      if (space <= 0) {
        getRowSpaceIterate(iRow, iNumberInRow + 4);
        put = startRow[iRow] + iNumberInRow;
      }
      indexColumn[put] = numberColumnsExtra_;
      convertRowToColumn[put] = i + startU;
      numberInRow[iRow] = iNumberInRow + 1;
      saveFromU = saveFromU - elementU[i] * region[iRow];
    } else {
      //zero out and save
      saveFromU += elementU[i];
      elementU[i] = 0.0;
    }
  }
  //in at end
  last = lastRow[maximumRowsExtra_];
  nextRow[last] = numberRowsExtra_;
  lastRow[maximumRowsExtra_] = numberRowsExtra_;
  lastRow[numberRowsExtra_] = last;
  nextRow[numberRowsExtra_] = maximumRowsExtra_;
  startRow[numberRowsExtra_] = startRow[maximumRowsExtra_];
  numberInRow[numberRowsExtra_] = 0;
  //column in at beginning (as empty)
  int *COIN_RESTRICT nextColumn = nextColumnArray_;
  int *COIN_RESTRICT lastColumn = lastColumnArray_;
  next = nextColumn[maximumColumnsExtra_];
  lastColumn[next] = numberColumnsExtra_;
  nextColumn[maximumColumnsExtra_] = numberColumnsExtra_;
  nextColumn[numberColumnsExtra_] = next;
  lastColumn[numberColumnsExtra_] = maximumColumnsExtra_;

  CoinFactorizationDouble pivotValue = 1.0 / saveFromU;

  pivotRegion[numberRowsExtra_] = pivotValue;
  //modify by pivot
  for (int i = 0; i < number; i++) {
    elementU[i] *= pivotValue;
  }
  maximumU_ = std::max(maximumU_, startU + number);
  numberRowsExtra_++;
  numberColumnsExtra_++;
  numberGoodU_++;
  numberPivots_++;
  if (numberRowsExtra_ > numberRows_ + 50) {
    int extra = factorElements_ >> 1;

    if (numberRowsExtra_ > numberRows_ + 100 + numberRows_ / 500) {
      if (extra < 2 * numberRows_) {
        extra = 2 * numberRows_;
      }
    } else {
      if (extra < 5 * numberRows_) {
        extra = 5 * numberRows_;
      }
    }
    int added = totalElements_ - factorElements_;

    if (added > extra && added > (factorElements_) << 1
      && 3 * totalElements_ > 2 * (lengthAreaU_ + lengthAreaL_)) {
      //status = 3;
      //if ( messageLevel_ & 4 ) {
      std::cout << "Factorization has " << totalElements_
                << ", basis had " << factorElements_ << std::endl;
      //}
      abort();
    }
  }
  // we are going to put another copy of R in R
  CoinFactorizationDouble *COIN_RESTRICT elementR = elementR_ + lengthAreaR_;
  int *COIN_RESTRICT indexRowR = indexRowR_ + lengthAreaR_;
  int *COIN_RESTRICT startR = startColumnRArray_ + maximumPivots_ + 1;
  int pivotRowNew = numberRowsExtra_ - 1;
  for (int i = 0; i < numberNonZero; i++) {
    int iRow = regionIndex[i];
    assert(pivotRowNew > iRow);
    next = nextColumn[iRow];
    int space;
    if (next != maximumColumnsExtra_)
      space = startR[next] - startR[iRow];
    else
      space = lengthAreaR_ - startR[iRow];
    int numberInR = numberInColumnPlus[iRow];
    if (space > numberInR) {
      // there is space
      int put = startR[iRow] + numberInR;
      numberInColumnPlus[iRow] = numberInR + 1;
      indexRowR[put] = pivotRowNew;
      elementR[put] = region[iRow];
      //add 4 for luck
      if (next == maximumColumnsExtra_)
        startR[maximumColumnsExtra_] = std::min(static_cast< int >(put + 4), lengthAreaR_);
    } else {
      // no space - do we shuffle?
      if (!getColumnSpaceIterateR(iRow, region[iRow], pivotRowNew)) {
        // printf("Need more space for R\n");
        numberInColumnPlus_.conditionalDelete();
        numberInColumnPlusArray_ = NULL;
        regionSparse->clear();
        break;
      }
    }
    region[iRow] = 0.0;
  }
  regionSparse->setNumElements(0);
}
/* Replaces one Column to basis,
   partial update in vector */
void CoinFactorization::replaceColumnPart3(CoinIndexedVector *regionSparse,
  CoinIndexedVector *partialUpdate,
  int pivotRow,
  double alpha)
{
  abort();
}
// Makes a non-singular basis by replacing variables
void CoinFactorization::makeNonSingular(int *COIN_RESTRICT sequence)
{
  // Replace bad ones by correct slack
  int *COIN_RESTRICT workArea = new int[numberRows_];
  for (int i = 0; i < numberRows_; i++)
    workArea[i] = -1;
  const int *COIN_RESTRICT pivot = pivotColumn();
  const int *COIN_RESTRICT permute = nextRowArray_;
  for (int i = 0; i < numberGoodU_; i++) {
    int iOriginal = pivot[i];
    assert(iOriginal >= 0);
    workArea[iOriginal] = i;
  }
  int iRow = 0;
  for (int i = 0; i < numberRows_; i++) {
    if (workArea[i] == -1) {
      while (permute[iRow] >= 0)
        iRow++;
      assert(iRow < numberRows_);
      // Put slack in basis
      sequence[i] = iRow;
      iRow++;
    }
  }
  delete[] workArea;
}
/* returns empty fake vector carved out of existing
   later - maybe use associated arrays */
CoinIndexedVector *
CoinFactorization::fakeVector(CoinIndexedVector *vector) const
{
  int oldCapacity = vector->capacity();
  CoinIndexedVector *newVector = new CoinIndexedVector();
  int n = (numberRows_ + 3) & ~3;
  newVector->setCapacity(oldCapacity - n);
  vector->setCapacity(n);
  newVector->setDenseVector(vector->denseVector() + n);
  newVector->setIndexVector(vector->getIndices() + n + ((n + 3) >> 2));
  // take this out when calmer - and think of best way
  memset(vector->getIndices() + vector->capacity(), 0, vector->capacity());
  //vector->checkClean();
  //newVector->checkClear();
  return newVector;
}
void CoinFactorization::deleteFakeVector(CoinIndexedVector *vector,
  CoinIndexedVector *fakeVector) const
{
  int n = fakeVector->capacity() + vector->capacity();
  //fakeVector->checkClear();
  fakeVector->setCapacity(0);
  fakeVector->setDenseVector(NULL);
  fakeVector->setIndexVector(NULL);
  delete fakeVector;
  vector->setCapacity(n);
  //vector->checkClean();
}
/* Updates one column (FTRAN) from regionSparse2
   Tries to do FT update
   number returned is negative if no room
   regionSparse starts as zero and is zero at end.
   Note - if regionSparse2 packed on input - will be packed on output
   long regions
*/
int CoinFactorization::updateColumnFT(CoinIndexedVector &regionSparse)
{
  CoinIndexedVector *newVector = fakeVector(&regionSparse);
  int returnCode = updateColumnFT(newVector, &regionSparse);
  deleteFakeVector(&regionSparse, newVector);
  return returnCode;
}
int CoinFactorization::updateColumnFTPart1(CoinIndexedVector &regionSparse2X)
{
  CoinIndexedVector *regionSparse = fakeVector(&regionSparse2X);
  CoinIndexedVector *regionSparse2 = &regionSparse2X;
  //permute and move indices into index array
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  int numberNonZero = regionSparse2->getNumElements();
  const int *permute = permuteArray_;
  int *COIN_RESTRICT index = regionSparse2->getIndices();
  double *COIN_RESTRICT region = regionSparse->denseVector();
  double *COIN_RESTRICT array = regionSparse2->denseVector();
  int *COIN_RESTRICT startColumnU = startColumnUArray_;
  // see if room
  int iColumn = numberColumnsExtra_;

  startColumnU[iColumn] = startColumnU[maximumColumnsExtra_];
  int start = startColumnU[iColumn];
  int space = lengthAreaU_ - (start + numberRowsExtra_);
  bool doFT = space >= 0;
  if (doFT) {
    regionIndex = indexRowUArray_ + start;
  } else {
    startColumnU[maximumColumnsExtra_] = lengthAreaU_ + 1;
  }
  assert(!regionSparse2->packedMode());
  for (int j = 0; j < numberNonZero; j++) {
    int iRow = index[j];
    double value = array[iRow];
    array[iRow] = 0.0;
    iRow = permute[iRow];
    region[iRow] = value;
    regionIndex[j] = iRow;
  }
  regionSparse->setNumElements(numberNonZero);
  if (collectStatistics_) {
    numberFtranCounts_++;
    ftranCountInput_ += numberNonZero;
  }

  //  ******* L
  updateColumnL(regionSparse, regionIndex);
  if (collectStatistics_)
    ftranCountAfterL_ += regionSparse->getNumElements();
  //permute extra
  //row bits here
  if (doFT)
    updateColumnRFT(regionSparse, regionIndex);
  else
    updateColumnR(regionSparse);
  if (collectStatistics_)
    ftranCountAfterR_ += regionSparse->getNumElements();
  deleteFakeVector(&regionSparse2X, regionSparse);
  // will be negative if no room
  if (doFT)
    return 1;
  else
    return -1;
}
void CoinFactorization::updateColumnFTPart2(CoinIndexedVector &regionSparse2X)
{
  CoinIndexedVector *regionSparse = fakeVector(&regionSparse2X);
  CoinIndexedVector *regionSparse2 = &regionSparse2X;
  //  ******* U
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  updateColumnU(regionSparse, regionIndex);
  permuteBack(regionSparse, regionSparse2);
  deleteFakeVector(&regionSparse2X, regionSparse);
}
/* Updates one column (FTRAN) - long region
   Tries to do FT update
   puts partial update in vector */
void CoinFactorization::updateColumnFT(CoinIndexedVector &regionSparseFT,
  CoinIndexedVector &partialUpdate,
  int which)
{
  abort();
}
/* Updates one column (FTRAN) long region */
int CoinFactorization::updateColumn(CoinIndexedVector &regionSparse) const
{
  CoinIndexedVector *newVector = fakeVector(&regionSparse);
  updateColumn(newVector, &regionSparse);
  deleteFakeVector(&regionSparse, newVector);
  return regionSparse.getNumElements();
}
/* Updates one column (FTRAN) from regionFT
   Tries to do FT update
   number returned is negative if no room.
   Also updates regionOther - long region*/
int CoinFactorization::updateTwoColumnsFT(CoinIndexedVector &regionSparseFT,
  CoinIndexedVector &regionSparseOther)
{
  CoinIndexedVector *newVector = fakeVector(&regionSparseFT);
  int returnCode = updateTwoColumnsFT(newVector,
    &regionSparseFT, &regionSparseOther);
  deleteFakeVector(&regionSparseFT, newVector);
  return returnCode;
}
/* Updates one column (BTRAN) - long region*/
int CoinFactorization::updateColumnTranspose(CoinIndexedVector &regionSparse) const
{
  CoinIndexedVector *newVector = fakeVector(&regionSparse);
  updateColumnTranspose(newVector, &regionSparse);
  deleteFakeVector(&regionSparse, newVector);
  return regionSparse.getNumElements();
}
/* Updates one column (FTRAN) - long region */
void CoinFactorization::updateColumnCpu(CoinIndexedVector &regionSparse, int whichCpu) const
{
  abort();
}
/* Updates one column (BTRAN) - long region */
void CoinFactorization::updateColumnTransposeCpu(CoinIndexedVector &regionSparse, int whichCpu) const
{
  abort();
}
/* Updates one full column (FTRAN) - long region */
void CoinFactorization::updateFullColumn(CoinIndexedVector &regionSparse) const
{
  if (!regionSparse.getNumElements())
    regionSparse.scan(0, numberRows_);
  CoinIndexedVector *newVector = fakeVector(&regionSparse);
  updateColumn(newVector, &regionSparse);
  deleteFakeVector(&regionSparse, newVector);
}
/* Updates one full column (BTRAN) - long region */
void CoinFactorization::updateFullColumnTranspose(CoinIndexedVector &regionSparse) const
{
  if (!regionSparse.getNumElements())
    regionSparse.scan(0, numberRows_);
  CoinIndexedVector *newVector = fakeVector(&regionSparse);
  updateColumnTranspose(newVector, &regionSparse);
  deleteFakeVector(&regionSparse, newVector);
}
/* Updates one column for dual steepest edge weights (FTRAN) - long region */
void CoinFactorization::updateWeights(CoinIndexedVector &regionSparse) const
{
  abort();
}
// Update partial Ftran by R update
void CoinFactorization::updatePartialUpdate(CoinIndexedVector &partialUpdate)
{
  abort();
}
#endif
#ifdef COIN_DEVELOP
extern double ncall_DZ;
extern double nrow_DZ;
extern double nslack_DZ;
extern double nU_DZ;
extern double nnz_DZ;
extern double nDone_DZ;
extern double ncall_SZ;
extern double nrow_SZ;
extern double nslack_SZ;
extern double nU_SZ;
extern double nnz_SZ;
extern double nDone_SZ;
void print_fac_stats()
{
  double mult = ncall_DZ ? 1.0 / ncall_DZ : 1.0;
  printf("UDen called %g times, average rows %g, average slacks %g, average (U-S) %g average nnz in %g average ops %g\n",
    ncall_DZ, mult * nrow_DZ, mult * nslack_DZ, mult * (nU_DZ - nslack_DZ), mult * nnz_DZ, mult * nDone_DZ);
  ncall_DZ = 0.0;
  nrow_DZ = 0.0;
  nslack_DZ = 0.0;
  nU_DZ = 0.0;
  nnz_DZ = 0.0;
  nDone_DZ = 0.0;
  mult = ncall_SZ ? 1.0 / ncall_SZ : 1.0;
  printf("USpars called %g times, average rows %g, average slacks %g, average (U-S) %g average nnz in %g average ops %g\n",
    ncall_SZ, mult * nrow_SZ, mult * nslack_SZ, mult * (nU_SZ - nslack_SZ), mult * nnz_SZ, mult * nDone_SZ);
  ncall_SZ = 0.0;
  nrow_SZ = 0.0;
  nslack_SZ = 0.0;
  nU_SZ = 0.0;
  nnz_SZ = 0.0;
  nDone_SZ = 0.0;
}
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
