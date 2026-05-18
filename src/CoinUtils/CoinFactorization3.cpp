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
  const double *A, const int lda, const int *ipiv, double *B,
  const int ldb);
}
#elif COIN_FACTORIZATION_DENSE_CODE == 3
// Intel compiler
#include "mkl_lapacke.h"
#endif
// For semi-sparse
#define BITS_PER_CHECK 8
#define CHECK_SHIFT 3
typedef unsigned char CoinCheckZero;
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

//:class CoinFactorization.  Deals with Factorization and Updates

/* Updates one column (FTRAN) from region2 and permutes.
   region1 starts as zero
   Note - if regionSparse2 packed on input - will be packed on output
   - returns un-permuted result in region2 and region1 is zero */
int CoinFactorization::updateColumn(CoinIndexedVector *regionSparse,
  CoinIndexedVector *regionSparse2,
  bool noPermute)
  const
{
#ifdef CLP_FACTORIZATION_INSTRUMENT
  double startTimeX = CoinCpuTime();
#endif
  //permute and move indices into index array
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  int numberNonZero;
  const int *permute = permuteArray_;
  double *COIN_RESTRICT region = regionSparse->denseVector();

#ifndef CLP_FACTORIZATION
  if (!noPermute) {
#endif
    numberNonZero = regionSparse2->getNumElements();
    int *COIN_RESTRICT index = regionSparse2->getIndices();
    double *COIN_RESTRICT array = regionSparse2->denseVector();
#ifndef CLP_FACTORIZATION
    bool packed = regionSparse2->packedMode();
    if (packed) {
      for (int j = 0; j < numberNonZero; j++) {
        int iRow = index[j];
        double value = array[j];
        array[j] = 0.0;
        iRow = permute[iRow];
        region[iRow] = value;
        regionIndex[j] = iRow;
      }
    } else {
#endif
      for (int j = 0; j < numberNonZero; j++) {
        int iRow = index[j];
        double value = array[iRow];
        array[iRow] = 0.0;
        iRow = permute[iRow];
        region[iRow] = value;
        regionIndex[j] = iRow;
      }
#ifndef CLP_FACTORIZATION
    }
#endif
    regionSparse->setNumElements(numberNonZero);
#ifndef CLP_FACTORIZATION
  } else {
    numberNonZero = regionSparse->getNumElements();
  }
#endif
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
  updateColumnR(regionSparse);
  if (collectStatistics_)
    ftranCountAfterR_ += regionSparse->getNumElements();

  //update counts
  //  ******* U
  updateColumnU(regionSparse, regionIndex);
  if (!doForrestTomlin_) {
    // Do PFI after everything else
    updateColumnPFI(regionSparse);
  }
#ifdef CLP_FACTORIZATION_INSTRUMENT
  numberUpdate++;
  timeInUpdate += CoinCpuTime() - startTimeX;
  averageLengthR += lengthR_;
  averageLengthU += lengthU_;
  averageLengthL += lengthL_;
#endif
  if (!noPermute) {
    permuteBack(regionSparse, regionSparse2);
    return regionSparse2->getNumElements();
  } else {
    return regionSparse->getNumElements();
  }
}
// Permutes back at end of updateColumn
void CoinFactorization::permuteBack(CoinIndexedVector *regionSparse,
  CoinIndexedVector *outVector) const
{
  // permute back
  int oldNumber = regionSparse->getNumElements();
  const int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  double *COIN_RESTRICT region = regionSparse->denseVector();
  int *COIN_RESTRICT outIndex = outVector->getIndices();
  double *COIN_RESTRICT out = outVector->denseVector();
  const int *COIN_RESTRICT permuteBack = pivotColumnBack();
  int number = 0;
  if (outVector->packedMode()) {
    for (int j = 0; j < oldNumber; j++) {
      int iRow = regionIndex[j];
      double value = region[iRow];
      region[iRow] = 0.0;
      if (fabs(value) > zeroTolerance_) {
        iRow = permuteBack[iRow];
        outIndex[number] = iRow;
        out[number++] = value;
      }
    }
  } else {
    int j = 0;
    if ((oldNumber & 1) != 0) {
      int iRow = regionIndex[0];
      j++;
      double value = region[iRow];
      region[iRow] = 0.0;
      if (fabs(value) > zeroTolerance_) {
        iRow = permuteBack[iRow];
        outIndex[number++] = iRow;
        out[iRow] = value;
      }
    }
    for (; j < oldNumber; j += 2) {
      int iRow0 = regionIndex[j];
      int iRow1 = regionIndex[j + 1];
      double value0 = region[iRow0];
      bool good0 = fabs(value0) > zeroTolerance_;
      double value1 = region[iRow1];
      bool good1 = fabs(value1) > zeroTolerance_;
      region[iRow0] = 0.0;
      region[iRow1] = 0.0;
      if (good0) {
        iRow0 = permuteBack[iRow0];
        outIndex[number++] = iRow0;
        out[iRow0] = value0;
      }
      if (good1) {
        iRow1 = permuteBack[iRow1];
        outIndex[number++] = iRow1;
        out[iRow1] = value1;
      }
    }
  }
  outVector->setNumElements(number);
  regionSparse->setNumElements(0);
}
//  updateColumnL.  Updates part of column (FTRANL)
void CoinFactorization::updateColumnL(CoinIndexedVector *regionSparse,
  int *COIN_RESTRICT regionIndex) const
{
  if (numberL_) {
    int number = regionSparse->getNumElements();
    int goSparse;
    // Guess at number at end
    if (sparseThreshold_ > 0) {
      if (ftranAverageAfterL_) {
        double newNumber = number * ftranAverageAfterL_;
        if (newNumber < sparseThreshold_ && (numberL_ << 2) > newNumber)
          goSparse = 2;
        else if (newNumber < sparseThreshold2_ && (numberL_ << 1) > newNumber)
          goSparse = 1;
        else
          goSparse = 0;
      } else {
        if (number < sparseThreshold_ && (numberL_ << 2) > number)
          goSparse = 2;
        else
          goSparse = 0;
      }
    } else {
      goSparse = 0;
    }
    switch (goSparse) {
    case 0: // densish
      updateColumnLDensish(regionSparse, regionIndex);
      break;
    case 1: // middling
      updateColumnLSparsish(regionSparse, regionIndex);
      break;
    case 2: // sparse
      updateColumnLSparse(regionSparse, regionIndex);
      break;
    }
  }
#ifdef COIN_FACTORIZATION_DENSE_CODE
  if (numberDense_) {
    //take off list
    int lastSparse = numberRows_ - numberDense_;
    int number = regionSparse->getNumElements();
    double *COIN_RESTRICT region = regionSparse->denseVector();
    int i = 0;
    bool doDense = false;
    while (i < number) {
      int iRow = regionIndex[i];
      if (iRow >= lastSparse) {
        doDense = true;
        regionIndex[i] = regionIndex[--number];
      } else {
        i++;
      }
    }
    if (doDense) {
#if COIN_FACTORIZATION_DENSE_CODE == 1
      char trans = 'N';
      int ione = 1;
      int info;
      COINUTILS_LAPACK_FUNC(dgetrs,DGETRS)(&trans,&numberDense_,&ione,
                       denseAreaAddress_,&numberDense_,densePermute_,
                       region+lastSparse,&numberDense_,&info,1);
#elif COIN_FACTORIZATION_DENSE_CODE==2
      clapack_dgetrs(CblasColMajor,CblasNoTrans,numberDense_,1,
                     denseAreaAddress_,numberDense_,densePermute_,
                     region+lastSparse,numberDense_);
#elif COIN_FACTORIZATION_DENSE_CODE==3
      LAPACKE_dgetrs(LAPACK_COL_MAJOR,'N',numberDense_,1,
                     denseAreaAddress_,numberDense_,densePermute_,
                     region+lastSparse,numberDense_);
#endif
      for (int i = lastSparse; i < numberRows_; i++) {
        double value = region[i];
        if (value) {
          if (fabs(value) >= 1.0e-15)
            regionIndex[number++] = i;
          else
            region[i] = 0.0;
        }
      }
      regionSparse->setNumElements(number);
    }
  }
#endif
}
// Updates part of column (FTRANL) when densish
void CoinFactorization::updateColumnLDensish(CoinIndexedVector *regionSparse,
  int *COIN_RESTRICT regionIndex)
  const
{
  double *COIN_RESTRICT region = regionSparse->denseVector();
  int number = regionSparse->getNumElements();
  int numberNonZero;
  double tolerance = zeroTolerance_;

  numberNonZero = 0;

  const int *COIN_RESTRICT startColumn = startColumnLArray_;
  const int *COIN_RESTRICT indexRow = indexRowLArray_;
  const CoinFactorizationDouble *COIN_RESTRICT element = elementLArray_;
  int last = numberRows_;
  assert(last == baseL_ + numberL_);
#if COIN_FACTORIZATION_DENSE_CODE
  //can take out last bit of sparse L as empty
  last -= numberDense_;
#endif
  int smallestIndex = numberRowsExtra_;
  // do easy ones
  for (int k = 0; k < number; k++) {
    int iPivot = regionIndex[k];
    if (iPivot >= baseL_)
      smallestIndex = std::min(iPivot, smallestIndex);
    else
      regionIndex[numberNonZero++] = iPivot;
  }
  // now others
  for (int i = smallestIndex; i < last; i++) {
    CoinFactorizationDouble pivotValue = region[i];

    if (fabs(pivotValue) > tolerance) {
      int start = startColumn[i];
      int end = startColumn[i + 1];
      for (int j = start; j < end; j++) {
        int iRow = indexRow[j];
        CoinFactorizationDouble result = region[iRow];
        CoinFactorizationDouble value = element[j];

        region[iRow] = result - value * pivotValue;
      }
      regionIndex[numberNonZero++] = i;
    } else {
      region[i] = 0.0;
    }
  }
  // and dense
  for (int i = last; i < numberRows_; i++) {
    CoinFactorizationDouble pivotValue = region[i];
    if (fabs(pivotValue) > tolerance) {
      regionIndex[numberNonZero++] = i;
    } else {
      region[i] = 0.0;
    }
  }
  regionSparse->setNumElements(numberNonZero);
}
// Updates part of column (FTRANL) when sparsish
void CoinFactorization::updateColumnLSparsish(CoinIndexedVector *regionSparse,
  int *COIN_RESTRICT regionIndex)
  const
{
  double *COIN_RESTRICT region = regionSparse->denseVector();
  int number = regionSparse->getNumElements();
  int numberNonZero;
  double tolerance = zeroTolerance_;

  numberNonZero = 0;

  const int *startColumn = startColumnLArray_;
  const int *indexRow = indexRowLArray_;
  const CoinFactorizationDouble *COIN_RESTRICT element = elementLArray_;
  int last = numberRows_;
  assert(last == baseL_ + numberL_);
#if COIN_FACTORIZATION_DENSE_CODE
  //can take out last bit of sparse L as empty
  last -= numberDense_;
#endif
  // mark known to be zero
  int nInBig = sizeof(int) / sizeof(int);
#if ABOCA_LITE_FACTORIZATION == 0
#define sparseOffset 0
#else
  int sparseOffset = ((regionSparse->capacity() & 0x80000000) != 0) ? sparseOffset_ : 0;
  assert(!sparseOffset);
#endif
  CoinCheckZero *COIN_RESTRICT mark = reinterpret_cast< CoinCheckZero * >(sparseArray_ + (2 + nInBig) * maximumRowsExtra_ + sparseOffset);
  int smallestIndex = numberRowsExtra_;
  // do easy ones
  for (int k = 0; k < number; k++) {
    int iPivot = regionIndex[k];
    if (iPivot < baseL_) {
      regionIndex[numberNonZero++] = iPivot;
    } else {
      smallestIndex = std::min(iPivot, smallestIndex);
      int iWord = iPivot >> CHECK_SHIFT;
      int iBit = iPivot - (iWord << CHECK_SHIFT);
      if (mark[iWord]) {
        mark[iWord] = static_cast< CoinCheckZero >(mark[iWord] | (1 << iBit));
      } else {
        mark[iWord] = static_cast< CoinCheckZero >(1 << iBit);
      }
    }
  }
  // now others
  // First do up to convenient power of 2
  int jLast = (smallestIndex + BITS_PER_CHECK - 1) >> CHECK_SHIFT;
  jLast = std::min((jLast << CHECK_SHIFT), last);
  int i;
  for (i = smallestIndex; i < jLast; i++) {
    CoinFactorizationDouble pivotValue = region[i];
    int start = startColumn[i];
    int end = startColumn[i + 1];

    if (fabs(pivotValue) > tolerance) {
      for (int j = start; j < end; j++) {
        int iRow = indexRow[j];
        CoinFactorizationDouble result = region[iRow];
        CoinFactorizationDouble value = element[j];
        region[iRow] = result - value * pivotValue;
        int iWord = iRow >> CHECK_SHIFT;
        int iBit = iRow - (iWord << CHECK_SHIFT);
        if (mark[iWord]) {
          mark[iWord] = static_cast< CoinCheckZero >(mark[iWord] | (1 << iBit));
        } else {
          mark[iWord] = static_cast< CoinCheckZero >(1 << iBit);
        }
      }
      regionIndex[numberNonZero++] = i;
    } else {
      region[i] = 0.0;
    }
  }

  int kLast = last >> CHECK_SHIFT;
  if (jLast < last) {
    // now do in chunks
    for (int k = (jLast >> CHECK_SHIFT); k < kLast; k++) {
      unsigned int iMark = mark[k];
      if (iMark) {
        // something in chunk - do all (as imark may change)
        i = k << CHECK_SHIFT;
        int iLast = i + BITS_PER_CHECK;
        for (; i < iLast; i++) {
          CoinFactorizationDouble pivotValue = region[i];
          int start = startColumn[i];
          int end = startColumn[i + 1];

          if (fabs(pivotValue) > tolerance) {
            int j;
            for (j = start; j < end; j++) {
              int iRow = indexRow[j];
              CoinFactorizationDouble result = region[iRow];
              CoinFactorizationDouble value = element[j];
              region[iRow] = result - value * pivotValue;
              int iWord = iRow >> CHECK_SHIFT;
              int iBit = iRow - (iWord << CHECK_SHIFT);
              if (mark[iWord]) {
                mark[iWord] = static_cast< CoinCheckZero >(mark[iWord] | (1 << iBit));
              } else {
                mark[iWord] = static_cast< CoinCheckZero >(1 << iBit);
              }
            }
            regionIndex[numberNonZero++] = i;
          } else {
            region[i] = 0.0;
          }
        }
        mark[k] = 0; // zero out marked
      }
    }
    i = kLast << CHECK_SHIFT;
  }
  for (; i < last; i++) {
    CoinFactorizationDouble pivotValue = region[i];
    int start = startColumn[i];
    int end = startColumn[i + 1];

    if (fabs(pivotValue) > tolerance) {
      for (int j = start; j < end; j++) {
        int iRow = indexRow[j];
        CoinFactorizationDouble result = region[iRow];
        CoinFactorizationDouble value = element[j];
        region[iRow] = result - value * pivotValue;
      }
      regionIndex[numberNonZero++] = i;
    } else {
      region[i] = 0.0;
    }
  }
  // Now dense part
  for (; i < numberRows_; i++) {
    double pivotValue = region[i];
    if (fabs(pivotValue) > tolerance) {
      regionIndex[numberNonZero++] = i;
    } else {
      region[i] = 0.0;
    }
  }
  // zero out ones that might have been skipped
  mark[smallestIndex >> CHECK_SHIFT] = 0;
  int kkLast = (numberRows_ + BITS_PER_CHECK - 1) >> CHECK_SHIFT;
  CoinZeroN(mark + kLast, kkLast - kLast);
  regionSparse->setNumElements(numberNonZero);
}
// Updates part of column (FTRANL) when sparse
void CoinFactorization::updateColumnLSparse(CoinIndexedVector *regionSparse,
  int *COIN_RESTRICT regionIndex)
  const
{
  double *COIN_RESTRICT region = regionSparse->denseVector();
  int number = regionSparse->getNumElements();
  int numberNonZero;
  double tolerance = zeroTolerance_;

  numberNonZero = 0;

  const int *startColumn = startColumnLArray_;
  const int *indexRow = indexRowLArray_;
  const CoinFactorizationDouble *COIN_RESTRICT element = elementLArray_;
  // use sparse_ as temporary area
  // mark known to be zero
#if ABOCA_LITE_FACTORIZATION == 0
#define sparseOffset 0
#else
  int sparseOffset = ((regionSparse->capacity() & 0x80000000) != 0) ? sparseOffset_ : 0;
  assert(!sparseOffset);
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
  for (int k = 0; k < number; k++) {
    int kPivot = regionIndex[k];
    if (kPivot >= baseL_) {
      assert(kPivot < numberRowsExtra_);
      //if (kPivot>=numberRowsExtra_) abort();
      if (!mark[kPivot]) {
        stack[0] = kPivot;
        int j = startColumn[kPivot + 1] - 1;
        int nStack = 0;
        while (nStack >= 0) {
          /* take off stack */
          if (j >= startColumn[kPivot]) {
            int jPivot = indexRow[j--];
            assert(jPivot >= baseL_ && jPivot < numberRowsExtra_);
            //if (jPivot<baseL_||jPivot>=numberRowsExtra_) abort();
            /* put back on stack */
            next[nStack] = j;
            if (!mark[jPivot]) {
              /* and new one */
              kPivot = jPivot;
              j = startColumn[kPivot + 1] - 1;
              stack[++nStack] = kPivot;
              assert(kPivot < numberRowsExtra_);
              //if (kPivot>=numberRowsExtra_) abort();
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
              assert(kPivot < numberRowsExtra_);
              j = next[nStack];
            }
          }
        }
      }
    } else {
      // just put on list
      regionIndex[numberNonZero++] = kPivot;
    }
  }
  for (int i = nList - 1; i >= 0; i--) {
    int iPivot = list[i];
    mark[iPivot] = 0;
    CoinFactorizationDouble pivotValue = region[iPivot];
    if (fabs(pivotValue) > tolerance) {
      regionIndex[numberNonZero++] = iPivot;
      for (int j = startColumn[iPivot];
           j < startColumn[iPivot + 1]; j++) {
        int iRow = indexRow[j];
        CoinFactorizationDouble value = element[j];
        region[iRow] -= value * pivotValue;
      }
    } else {
      region[iPivot] = 0.0;
    }
  }
  regionSparse->setNumElements(numberNonZero);
}
/* Updates one column (FTRAN) from region2
   Tries to do FT update
   number returned is negative if no room.
   Also updates region3
   region1 starts as zero and is zero at end */
int CoinFactorization::updateTwoColumnsFT(CoinIndexedVector *regionSparse1,
  CoinIndexedVector *regionSparse2,
  CoinIndexedVector *regionSparse3,
  bool noPermuteRegion3)
{
#ifdef CLP_FACTORIZATION_INSTRUMENT
  double startTimeX = CoinCpuTime();
#endif
#if 1
  //#ifdef NDEBUG
  //#undef NDEBUG
  //#endif
  //#define COIN_DEBUG
#ifdef COIN_DEBUG
  regionSparse1->checkClean();
  CoinIndexedVector save2(*regionSparse2);
  CoinIndexedVector save3(*regionSparse3);
#endif
  CoinIndexedVector *regionFT;
  CoinIndexedVector *regionUpdate;
  int *COIN_RESTRICT regionIndex;
  int numberNonZero;
  const int *permute = permuteArray_;
  int *COIN_RESTRICT index;
  double *COIN_RESTRICT region;
  if (!noPermuteRegion3) {
    regionFT = regionSparse3;
    regionUpdate = regionSparse1;
    //permute and move indices into index array
    regionIndex = regionUpdate->getIndices();
    //int numberNonZero;
    region = regionUpdate->denseVector();

    numberNonZero = regionSparse3->getNumElements();
    int *COIN_RESTRICT index = regionSparse3->getIndices();
    double *COIN_RESTRICT array = regionSparse3->denseVector();
    assert(!regionSparse3->packedMode());
    for (int j = 0; j < numberNonZero; j++) {
      int iRow = index[j];
      double value = array[iRow];
      array[iRow] = 0.0;
      iRow = permute[iRow];
      region[iRow] = value;
      regionIndex[j] = iRow;
    }
    regionUpdate->setNumElements(numberNonZero);
  } else {
    regionFT = regionSparse1;
    regionUpdate = regionSparse3;
  }
  //permute and move indices into index array (in U)
  regionIndex = regionFT->getIndices();
  numberNonZero = regionSparse2->getNumElements();
  index = regionSparse2->getIndices();
  region = regionFT->denseVector();
  double *COIN_RESTRICT array = regionSparse2->denseVector();
  int *COIN_RESTRICT startColumnU = startColumnUArray_;
  int start = startColumnU[maximumColumnsExtra_];
  startColumnU[numberColumnsExtra_] = start;
  regionIndex = indexRowUArray_ + start;

#ifndef ABC_USE_COIN_FACTORIZATION
  assert(regionSparse2->packedMode());
#else
  if (regionSparse2->packedMode()) {
#endif
  for (int j = 0; j < numberNonZero; j++) {
    int iRow = index[j];
    double value = array[j];
    array[j] = 0.0;
    iRow = permute[iRow];
    region[iRow] = value;
    regionIndex[j] = iRow;
  }
#ifdef ABC_USE_COIN_FACTORIZATION
}
else
{
  // not packed
  for (int j = 0; j < numberNonZero; j++) {
    int iRow = index[j];
    double value = array[iRow];
    array[iRow] = 0.0;
    iRow = permute[iRow];
    region[iRow] = value;
    regionIndex[j] = iRow;
  }
}
#endif
regionFT->setNumElements(numberNonZero);
if (collectStatistics_) {
  numberFtranCounts_ += 2;
  ftranCountInput_ += regionFT->getNumElements() + regionUpdate->getNumElements();
}

//  ******* L
updateColumnL(regionFT, regionIndex);
updateColumnL(regionUpdate, regionUpdate->getIndices());
if (collectStatistics_)
  ftranCountAfterL_ += regionFT->getNumElements() + regionUpdate->getNumElements();
//permute extra
//row bits here
updateColumnRFT(regionFT, regionIndex);
updateColumnR(regionUpdate);
if (collectStatistics_)
  ftranCountAfterR_ += regionFT->getNumElements() + regionUpdate->getNumElements();
//  ******* U - see if densish
// Guess at number at end
int goSparse = 0;
if (sparseThreshold_ > 0) {
  int numberNonZero = (regionUpdate->getNumElements() + regionFT->getNumElements()) >> 1;
  if (ftranAverageAfterR_) {
    double newNumber = numberNonZero * ftranAverageAfterU_;
    if (newNumber < sparseThreshold_)
      goSparse = 2;
    else if (newNumber < sparseThreshold2_)
      goSparse = 1;
  } else {
    if (numberNonZero < sparseThreshold_)
      goSparse = 2;
  }
}
#ifndef COIN_FAST_CODE
assert(slackValue_ == -1.0);
#endif
if (!goSparse && numberRows_ < 1000) {
  double *COIN_RESTRICT arrayFT = regionFT->denseVector();
  int *COIN_RESTRICT indexFT = regionFT->getIndices();
  int numberNonZeroFT;
  double *COIN_RESTRICT arrayUpdate = regionUpdate->denseVector();
  int *COIN_RESTRICT indexUpdate = regionUpdate->getIndices();
  int numberNonZeroUpdate;
  updateTwoColumnsUDensish(numberNonZeroFT, arrayFT, indexFT,
    numberNonZeroUpdate, arrayUpdate, indexUpdate);
  regionFT->setNumElements(numberNonZeroFT);
  regionUpdate->setNumElements(numberNonZeroUpdate);
  if (collectStatistics_) {
    ftranCountAfterU_ += numberNonZeroFT + numberNonZeroUpdate;
#ifdef CLP_FACTORIZATION_INSTRUMENT
    scaledLengthDense += numberDense_ * numberNonZeroFT;
    scaledLengthDenseSquared += numberDense_ * numberDense_ * numberNonZeroFT;
    scaledLengthL += lengthL_ * numberNonZeroFT;
    scaledLengthR += lengthR_ * numberNonZeroFT;
    scaledLengthU += lengthU_ * numberNonZeroFT;
    scaledLengthDense += numberDense_ * numberNonZeroUpdate;
    scaledLengthDenseSquared += numberDense_ * numberDense_ * numberNonZeroUpdate;
    scaledLengthL += lengthL_ * numberNonZeroUpdate;
    scaledLengthR += lengthR_ * numberNonZeroUpdate;
    scaledLengthU += lengthU_ * numberNonZeroUpdate;
#endif
  }
} else {
  // sparse
  updateColumnU(regionFT, regionIndex);
  updateColumnU(regionUpdate, regionUpdate->getIndices());
}
permuteBack(regionFT, regionSparse2);
if (!noPermuteRegion3) {
  permuteBack(regionUpdate, regionSparse3);
}
#ifdef COIN_DEBUG
int n2 = regionSparse2->getNumElements();
regionSparse1->checkClean();
int n2a = updateColumnFT(regionSparse1, &save2);
assert(n2 == n2a);
{
  int j;
  double *COIN_RESTRICT regionA = save2.denseVector();
  int *COIN_RESTRICT indexA = save2.getIndices();
  double *COIN_RESTRICT regionB = regionSparse2->denseVector();
  int *COIN_RESTRICT indexB = regionSparse2->getIndices();
  for (j = 0; j < n2; j++) {
    int k = indexA[j];
    assert(k == indexB[j]);
    CoinFactorizationDouble value = regionA[j];
    assert(value == regionB[j]);
  }
}
updateColumn(&save3,
  &save3,
  noPermuteRegion3);
int n3 = regionSparse3->getNumElements();
assert(n3 == save3.getNumElements());
{
  int j;
  double *COIN_RESTRICT regionA = save3.denseVector();
  int *COIN_RESTRICT indexA = save3.getIndices();
  double *COIN_RESTRICT regionB = regionSparse3->denseVector();
  int *COIN_RESTRICT indexB = regionSparse3->getIndices();
  for (j = 0; j < n3; j++) {
    int k = indexA[j];
    assert(k == indexB[j]);
    CoinFactorizationDouble value = regionA[k];
    assert(value == regionB[k]);
  }
}
//*regionSparse2=save2;
//*regionSparse3=save3;
printf("REGION2 %d els\n", regionSparse2->getNumElements());
regionSparse2->print();
printf("REGION3 %d els\n", regionSparse3->getNumElements());
regionSparse3->print();
#endif
return regionSparse2->getNumElements();
#else
  int returnCode = updateColumnFT(regionSparse1,
    regionSparse2);
  assert(noPermuteRegion3);
  updateColumn(regionSparse3,
    regionSparse3,
    noPermuteRegion3);
#ifdef CLP_FACTORIZATION_INSTRUMENT
  numberUpdateTwoFT++;
  timeInUpdateTwoFT += CoinCpuTime() - startTimeX;
  averageLengthR += 2 * lengthR_;
  averageLengthU += 2 * lengthU_;
  averageLengthL += 2 * lengthL_;
#endif
  //printf("REGION2 %d els\n",regionSparse2->getNumElements());
  //regionSparse2->print();
  //printf("REGION3 %d els\n",regionSparse3->getNumElements());
  //regionSparse3->print();
  return returnCode;
#endif
}
// Updates part of 2 columns (FTRANU) real work
void CoinFactorization::updateTwoColumnsUDensish(
  int &numberNonZero1,
  double *COIN_RESTRICT region1,
  int *COIN_RESTRICT index1,
  int &numberNonZero2,
  double *COIN_RESTRICT region2,
  int *COIN_RESTRICT index2) const
{
  double tolerance = zeroTolerance_;
  const int *COIN_RESTRICT startColumn = startColumnUArray_;
  const int *COIN_RESTRICT indexRow = indexRowUArray_;
  const CoinFactorizationDouble *COIN_RESTRICT element = elementUArray_;
  int numberNonZeroA = 0;
  int numberNonZeroB = 0;
  const int *numberInColumn = numberInColumnArray_;
  const CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;

  for (int i = numberU_ - 1; i >= numberSlacks_; i--) {
    CoinFactorizationDouble pivotValue2 = region2[i];
    bool doTwo = fabs(pivotValue2) > tolerance;
    region2[i] = 0.0;
    CoinFactorizationDouble pivotValue1 = region1[i];
    bool doOne = fabs(pivotValue1) > tolerance;
    region1[i] = 0.0;
    if (doTwo) {
      int start = startColumn[i];
      const CoinFactorizationDouble *COIN_RESTRICT thisElement = element + start;
      const int *COIN_RESTRICT thisIndex = indexRow + start;
      if (!doOne) {
	int n = numberInColumn[i];
	int end = n & (~1);
	bool doExtra = n != end;
        // just region 2
        for (int j = 0; j< end; j +=2 ) {
          int iRow0 = thisIndex[j];
          int iRow1 = thisIndex[j+1];
          CoinFactorizationDouble regionValue20 = region2[iRow0];
          CoinFactorizationDouble regionValue21 = region2[iRow1];
          CoinFactorizationDouble value0 = thisElement[j];
          CoinFactorizationDouble value1 = thisElement[j+1];
          region2[iRow0] = regionValue20 - value0 * pivotValue2;
          region2[iRow1] = regionValue21 - value1 * pivotValue2;
        }
	if (doExtra) {
          int iRow = thisIndex[end];
          CoinFactorizationDouble regionValue2 = region2[iRow];
          CoinFactorizationDouble value = thisElement[end];
          region2[iRow] = regionValue2 - value * pivotValue2;
	}
        pivotValue2 *= pivotRegion[i];
        region2[i] = pivotValue2;
        index2[numberNonZeroB++] = i;
      } else {
	int n = numberInColumn[i];
        for (int j = 0; j< n ;j++) {
          int iRow = thisIndex[j];
          CoinFactorizationDouble value = thisElement[j];
          CoinFactorizationDouble regionValue1 = region1[iRow];
          CoinFactorizationDouble regionValue2 = region2[iRow];
          region1[iRow] = regionValue1 - value * pivotValue1;
          region2[iRow] = regionValue2 - value * pivotValue2;
        }
        pivotValue1 *= pivotRegion[i];
        pivotValue2 *= pivotRegion[i];
        region1[i] = pivotValue1;
        index1[numberNonZeroA++] = i;
        region2[i] = pivotValue2;
        index2[numberNonZeroB++] = i;
      }
    } else if (doOne) {
      int n = numberInColumn[i];
      int end = n & (~1);
      bool doExtra = n != end;
      int start = startColumn[i];
      const CoinFactorizationDouble *COIN_RESTRICT thisElement = element + start;
      const int *COIN_RESTRICT thisIndex = indexRow + start;
      // just region 1
      for (int j = 0; j< end; j +=2 ) {
	int iRow0 = thisIndex[j];
	int iRow1 = thisIndex[j+1];
	CoinFactorizationDouble regionValue10 = region1[iRow0];
	CoinFactorizationDouble regionValue11 = region1[iRow1];
	CoinFactorizationDouble value0 = thisElement[j];
	CoinFactorizationDouble value1 = thisElement[j+1];
	region1[iRow0] = regionValue10 - value0 * pivotValue1;
	region1[iRow1] = regionValue11 - value1 * pivotValue1;
      }
      if (doExtra) {
	int iRow = thisIndex[end];
	CoinFactorizationDouble regionValue1 = region1[iRow];
	CoinFactorizationDouble value = thisElement[end];
	region1[iRow] = regionValue1 - value * pivotValue1;
      }
      pivotValue1 *= pivotRegion[i];
      region1[i] = pivotValue1;
      index1[numberNonZeroA++] = i;
    }
  }
  // Slacks

  for (int i = numberSlacks_ - 1; i >= 0; i--) {
    double value2 = region2[i];
    double value1 = region1[i];
    bool value1NonZero = (value1 != 0.0);
    if (fabs(value2) > tolerance) {
      region2[i] = -value2;
      index2[numberNonZeroB++] = i;
    } else {
      region2[i] = 0.0;
    }
    if (value1NonZero) {
      index1[numberNonZeroA] = i;
      if (fabs(value1) > tolerance) {
        region1[i] = -value1;
        numberNonZeroA++;
      } else {
        region1[i] = 0.0;
      }
    }
  }
  numberNonZero1 = numberNonZeroA;
  numberNonZero2 = numberNonZeroB;
}
#ifdef COIN_FACTORIZATION_DIAGNOSE
static int numberTimesX = 0;
static int numberSparseX = 0;
double sumNumberSparseX = 0.0;
static int numberSparsishX = 0;
double sumNumberSparsishX = 0.0;
static int numberDensishX = 0;
double sumNumberDensishX = 0.0;
#endif
//  updateColumnU.  Updates part of column (FTRANU)
void CoinFactorization::updateColumnU(CoinIndexedVector *regionSparse,
  int *COIN_RESTRICT indexIn) const
{
  int numberNonZero = regionSparse->getNumElements();

  int goSparse;
  // Guess at number at end
  if (sparseThreshold_ > 0) {
    if (ftranAverageAfterR_) {
      double newNumber = numberNonZero * ftranAverageAfterU_;
      if (newNumber < sparseThreshold_)
        goSparse = 2;
      else if (newNumber < sparseThreshold2_)
        goSparse = 1;
      else
        goSparse = 0;
    } else {
      if (numberNonZero < sparseThreshold_)
        goSparse = 2;
      else
        goSparse = 0;
    }
  } else {
    goSparse = 0;
  }
#ifdef COIN_FACTORIZATION_DIAGNOSE
  numberTimesX++;
  if (!goSparse) {
    numberDensishX++;
    sumNumberDensishX += numberNonZero;
  } else if (goSparse == 1) {
    numberSparsishX++;
    sumNumberSparsishX += numberNonZero;
  } else {
    numberSparseX++;
    sumNumberSparseX += numberNonZero;
  }
  if ((numberTimesX % 1000) == 0) {
    double averageDensish = (numberDensishX) ? sumNumberDensishX / numberDensishX : 0.0;
    double averageSparsish = (numberSparsishX) ? sumNumberSparsishX / numberSparsishX : 0.0;
    double averageSparse = (numberSparseX) ? sumNumberSparseX / numberSparseX : 0.0;
    printf("sparsity D,ish,S (%d,%g) , (%d,%g) , (%d,%g) - ftranFactor %g\n",
      numberDensishX, averageDensish,
      numberSparsishX, averageSparsish,
      numberSparseX, averageSparse, ftranAverageAfterU_);
  }
#endif
  switch (goSparse) {
  case 0: // densish
  {
    double *COIN_RESTRICT region = regionSparse->denseVector();
    int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
    int numberNonZero = updateColumnUDensish(region, regionIndex);
    regionSparse->setNumElements(numberNonZero);
  } break;
  case 1: // middling
    updateColumnUSparsish(regionSparse, indexIn);
    break;
  case 2: // sparse
    updateColumnUSparse(regionSparse, indexIn);
    break;
  }
  if (collectStatistics_) {
    ftranCountAfterU_ += regionSparse->getNumElements();
#ifdef CLP_FACTORIZATION_INSTRUMENT
    int numberNonZero = regionSparse->getNumElements();
    scaledLengthDense += numberDense_ * numberNonZero;
    scaledLengthDenseSquared += numberDense_ * numberDense_ * numberNonZero;
    scaledLengthL += lengthL_ * numberNonZero;
    scaledLengthR += lengthR_ * numberNonZero;
    scaledLengthU += lengthU_ * numberNonZero;
#endif
  }
}
#ifdef COIN_DEVELOP
double ncall_DZ = 0.0;
double nrow_DZ = 0.0;
double nslack_DZ = 0.0;
double nU_DZ = 0.0;
double nnz_DZ = 0.0;
double nDone_DZ = 0.0;
#endif
// Updates part of column (FTRANU) real work
int CoinFactorization::updateColumnUDensish(double *COIN_RESTRICT region,
  int *COIN_RESTRICT regionIndex) const
{
  double tolerance = zeroTolerance_;
  const int *startColumn = startColumnUArray_;
  const int *indexRow = indexRowUArray_;
  const CoinFactorizationDouble *COIN_RESTRICT element = elementUArray_;
  int numberNonZero = 0;
  const int *numberInColumn = numberInColumnArray_;
  const CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;
#ifdef COIN_DEVELOP
  ncall_DZ++;
  nrow_DZ += numberRows_;
  nslack_DZ += numberSlacks_;
  nU_DZ += numberU_;
#endif

  for (int i = numberU_ - 1; i >= numberSlacks_; i--) {
    CoinFactorizationDouble pivotValue = region[i];
    if (pivotValue) {
#ifdef COIN_DEVELOP
      nnz_DZ++;
#endif
      region[i] = 0.0;
      if (fabs(pivotValue) > tolerance) {
        int start = startColumn[i];
        const CoinFactorizationDouble *COIN_RESTRICT thisElement = element + start;
        const int *thisIndex = indexRow + start;
#ifdef COIN_DEVELOP
        nDone_DZ += numberInColumn[i];
#endif
        for (int j = numberInColumn[i] - 1; j >= 0; j--) {
          int iRow = thisIndex[j];
          CoinFactorizationDouble regionValue = region[iRow];
          CoinFactorizationDouble value = thisElement[j];
          region[iRow] = regionValue - value * pivotValue;
        }
        pivotValue *= pivotRegion[i];
        region[i] = pivotValue;
        regionIndex[numberNonZero++] = i;
      }
    }
  }

  // now do slacks
#ifndef COIN_FAST_CODE
  if (slackValue_ == -1.0) {
#endif
  for (int i = numberSlacks_ - 1; i >= 0; i--) {
    double value = region[i];
    if (value) {
      region[i] = -value;
      regionIndex[numberNonZero] = i;
      if (fabs(value) > tolerance)
        numberNonZero++;
      else
        region[i] = 0.0;
    }
  }
#ifndef COIN_FAST_CODE
  } else {
    assert(slackValue_ == 1.0);
    for (int i = numberSlacks_ - 1; i >= 0; i--) {
      double value = region[i];
      double absValue = fabs(value);
      if (value) {
        region[i] = 0.0;
        if (absValue > tolerance) {
          region[i] = value;
          regionIndex[numberNonZero++] = i;
        }
      }
    }
  }
#endif
  return numberNonZero;
}
//  updateColumnU.  Updates part of column (FTRANU)
/*
  Since everything is in order I should be able to do a better job of
  marking stuff - think.  Also as L is static maybe I can do something
  better there (I know I could if I marked the depth of every element
  but that would lead to other inefficiencies.
*/
void CoinFactorization::updateColumnUSparse(CoinIndexedVector *regionSparse,
  int *COIN_RESTRICT indexIn) const
{
  int numberNonZero = regionSparse->getNumElements();
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  double *COIN_RESTRICT region = regionSparse->denseVector();
  double tolerance = zeroTolerance_;
  const int *startColumn = startColumnUArray_;
  const int *indexRow = indexRowUArray_;
  const CoinFactorizationDouble *COIN_RESTRICT element = elementUArray_;
  const CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;
  // use sparse_ as temporary area
  // mark known to be zero
#if ABOCA_LITE_FACTORIZATION == 0
#define sparseOffset 0
#else
  int sparseOffset = ((regionSparse->capacity() & 0x80000000) != 0) ? sparseOffset_ : 0;
  assert(!sparseOffset);
#endif
  int *COIN_RESTRICT stack = sparseArray_ + sparseOffset; /* pivot */
  int *COIN_RESTRICT list = stack + maximumRowsExtra_; /* final list */
  int *COIN_RESTRICT next = reinterpret_cast< int * >(list + maximumRowsExtra_); /* jnext */
  char *COIN_RESTRICT mark = reinterpret_cast< char * >(next + maximumRowsExtra_);
#ifdef COIN_DEBUG
  for (int i = 0; i < maximumRowsExtra_; i++) {
    assert(!mark[i]);
  }
#endif

  // move slacks to end of stack list
  int *COIN_RESTRICT putLast = stack + maximumRowsExtra_;
  int *COIN_RESTRICT put = putLast;

  const int *numberInColumn = numberInColumnArray_;
  int nList = 0;
  for (int i = 0; i < numberNonZero; i++) {
    int kPivot = indexIn[i];
    stack[0] = kPivot;
    int j = startColumn[kPivot] + numberInColumn[kPivot] - 1;
    int nStack = 1;
    next[0] = j;
    while (nStack) {
      /* take off stack */
      int kPivot = stack[--nStack];
      if (mark[kPivot] != 1) {
        j = next[nStack];
        if (j >= startColumn[kPivot]) {
          kPivot = indexRow[j--];
          /* put back on stack */
          next[nStack++] = j;
          if (!mark[kPivot]) {
            /* and new one */
            int numberIn = numberInColumn[kPivot];
            if (numberIn) {
              j = startColumn[kPivot] + numberIn - 1;
              stack[nStack] = kPivot;
              mark[kPivot] = 2;
              next[nStack++] = j;
            } else {
              // can do immediately
              /* finished so mark */
              mark[kPivot] = 1;
              if (kPivot >= numberSlacks_) {
                list[nList++] = kPivot;
              } else {
                // slack - put at end
                --put;
                *put = kPivot;
              }
            }
          }
        } else {
          /* finished so mark */
          mark[kPivot] = 1;
          if (kPivot >= numberSlacks_) {
            list[nList++] = kPivot;
          } else {
            // slack - put at end
            assert(!numberInColumn[kPivot]);
            --put;
            *put = kPivot;
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
    region[iPivot] = 0.0;
    if (fabs(pivotValue) > tolerance) {
      int start = startColumn[iPivot];
      int number = numberInColumn[iPivot];

      int j;
      for (j = start; j < start + number; j++) {
        CoinFactorizationDouble value = element[j];
        int iRow = indexRow[j];
        region[iRow] -= value * pivotValue;
      }
      pivotValue *= pivotRegion[iPivot];
      region[iPivot] = pivotValue;
      regionIndex[numberNonZero++] = iPivot;
    }
  }
  // slacks
#ifndef COIN_FAST_CODE
  if (slackValue_ == 1.0) {
    for (; put < putLast; put++) {
      int iPivot = *put;
      mark[iPivot] = 0;
      CoinFactorizationDouble pivotValue = region[iPivot];
      region[iPivot] = 0.0;
      if (fabs(pivotValue) > tolerance) {
        region[iPivot] = pivotValue;
        regionIndex[numberNonZero++] = iPivot;
      }
    }
  } else {
#endif
    for (; put < putLast; put++) {
      int iPivot = *put;
      mark[iPivot] = 0;
      CoinFactorizationDouble pivotValue = region[iPivot];
      region[iPivot] = 0.0;
      if (fabs(pivotValue) > tolerance) {
        region[iPivot] = -pivotValue;
        regionIndex[numberNonZero++] = iPivot;
      }
    }
#ifndef COIN_FAST_CODE
  }
#endif
  regionSparse->setNumElements(numberNonZero);
}
//  updateColumnU.  Updates part of column (FTRANU)
/*
  Since everything is in order I should be able to do a better job of
  marking stuff - think.  Also as L is static maybe I can do something
  better there (I know I could if I marked the depth of every element
  but that would lead to other inefficiencies.
*/
#ifdef COIN_DEVELOP
double ncall_SZ = 0.0;
double nrow_SZ = 0.0;
double nslack_SZ = 0.0;
double nU_SZ = 0.0;
double nnz_SZ = 0.0;
double nDone_SZ = 0.0;
#endif
void CoinFactorization::updateColumnUSparsish(CoinIndexedVector *regionSparse,
  int *COIN_RESTRICT indexIn) const
{
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  // mark known to be zero
#if ABOCA_LITE_FACTORIZATION == 0
#define sparseOffset 0
#else
  int sparseOffset = ((regionSparse->capacity() & 0x80000000) != 0) ? sparseOffset_ : 0;
  assert(!sparseOffset);
#endif
  int *COIN_RESTRICT stack = sparseArray_ + sparseOffset; /* pivot */
  int *COIN_RESTRICT list = stack + maximumRowsExtra_; /* final list */
  int *COIN_RESTRICT next = reinterpret_cast< int * >(list + maximumRowsExtra_); /* jnext */
  CoinCheckZero *COIN_RESTRICT mark = reinterpret_cast< CoinCheckZero * >(next + maximumRowsExtra_);
  const int *numberInColumn = numberInColumnArray_;
#ifdef COIN_DEBUG
  for (int i = 0; i < maximumRowsExtra_; i++) {
    assert(!mark[i]);
  }
#endif

  int nMarked = 0;
  int numberNonZero = regionSparse->getNumElements();
  double *COIN_RESTRICT region = regionSparse->denseVector();
  double tolerance = zeroTolerance_;
  const int *startColumn = startColumnUArray_;
  const int *indexRow = indexRowUArray_;
  const CoinFactorizationDouble *COIN_RESTRICT element = elementUArray_;
  const CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;
#ifdef COIN_DEVELOP
  ncall_SZ++;
  nrow_SZ += numberRows_;
  nslack_SZ += numberSlacks_;
  nU_SZ += numberU_;
#endif

  for (int ii = 0; ii < numberNonZero; ii++) {
    int iPivot = indexIn[ii];
    int iWord = iPivot >> CHECK_SHIFT;
    int iBit = iPivot - (iWord << CHECK_SHIFT);
    if (mark[iWord]) {
      mark[iWord] = static_cast< CoinCheckZero >(mark[iWord] | (1 << iBit));
    } else {
      mark[iWord] = static_cast< CoinCheckZero >(1 << iBit);
      stack[nMarked++] = iWord;
    }
  }
  numberNonZero = 0;
  // First do down to convenient power of 2
  int jLast = (numberU_ - 1) >> CHECK_SHIFT;
  jLast = std::max((jLast << CHECK_SHIFT), static_cast< int >(numberSlacks_));
  int i;
  for (i = numberU_ - 1; i >= jLast; i--) {
    CoinFactorizationDouble pivotValue = region[i];
    region[i] = 0.0;
    if (fabs(pivotValue) > tolerance) {
#ifdef COIN_DEVELOP
      nnz_SZ++;
#endif
      int start = startColumn[i];
      const CoinFactorizationDouble *COIN_RESTRICT thisElement = element + start;
      const int *thisIndex = indexRow + start;

#ifdef COIN_DEVELOP
      nDone_SZ += numberInColumn[i];
#endif
      for (int j = numberInColumn[i] - 1; j >= 0; j--) {
        int iRow0 = thisIndex[j];
        CoinFactorizationDouble regionValue0 = region[iRow0];
        CoinFactorizationDouble value0 = thisElement[j];
        int iWord = iRow0 >> CHECK_SHIFT;
        int iBit = iRow0 - (iWord << CHECK_SHIFT);
        if (mark[iWord]) {
          mark[iWord] = static_cast< CoinCheckZero >(mark[iWord] | (1 << iBit));
        } else {
          mark[iWord] = static_cast< CoinCheckZero >(1 << iBit);
          stack[nMarked++] = iWord;
        }
        region[iRow0] = regionValue0 - value0 * pivotValue;
      }
      pivotValue *= pivotRegion[i];
      region[i] = pivotValue;
      regionIndex[numberNonZero++] = i;
    }
  }
  int kLast = (numberSlacks_ + BITS_PER_CHECK - 1) >> CHECK_SHIFT;
  if (jLast > numberSlacks_) {
    // now do in chunks
    for (int k = (jLast >> CHECK_SHIFT) - 1; k >= kLast; k--) {
      unsigned int iMark = mark[k];
      if (iMark) {
        // something in chunk - do all (as imark may change)
        int iLast = k << CHECK_SHIFT;
        for (i = iLast + BITS_PER_CHECK - 1; i >= iLast; i--) {
          CoinFactorizationDouble pivotValue = region[i];
          if (pivotValue) {
#ifdef COIN_DEVELOP
            nnz_SZ++;
#endif
            region[i] = 0.0;
            if (fabs(pivotValue) > tolerance) {
              int start = startColumn[i];
              const CoinFactorizationDouble *COIN_RESTRICT thisElement = element + start;
              const int *thisIndex = indexRow + start;
#ifdef COIN_DEVELOP
              nDone_SZ += numberInColumn[i];
#endif
              for (int j = numberInColumn[i] - 1; j >= 0; j--) {
                int iRow0 = thisIndex[j];
                CoinFactorizationDouble regionValue0 = region[iRow0];
                CoinFactorizationDouble value0 = thisElement[j];
                int iWord = iRow0 >> CHECK_SHIFT;
                int iBit = iRow0 - (iWord << CHECK_SHIFT);
                if (mark[iWord]) {
                  mark[iWord] = static_cast< CoinCheckZero >(mark[iWord] | (1 << iBit));
                } else {
                  mark[iWord] = static_cast< CoinCheckZero >(1 << iBit);
                  stack[nMarked++] = iWord;
                }
                region[iRow0] = regionValue0 - value0 * pivotValue;
              }
              pivotValue *= pivotRegion[i];
              region[i] = pivotValue;
              regionIndex[numberNonZero++] = i;
            }
          }
        }
        mark[k] = 0;
      }
    }
    i = (kLast << CHECK_SHIFT) - 1;
  }
  for (; i >= numberSlacks_; i--) {
    CoinFactorizationDouble pivotValue = region[i];
    region[i] = 0.0;
    if (fabs(pivotValue) > tolerance) {
#ifdef COIN_DEVELOP
      nnz_SZ++;
#endif
      int start = startColumn[i];
      const CoinFactorizationDouble *COIN_RESTRICT thisElement = element + start;
      const int *thisIndex = indexRow + start;
#ifdef COIN_DEVELOP
      nDone_SZ += numberInColumn[i];
#endif
      for (int j = numberInColumn[i] - 1; j >= 0; j--) {
        int iRow0 = thisIndex[j];
        CoinFactorizationDouble regionValue0 = region[iRow0];
        CoinFactorizationDouble value0 = thisElement[j];
        int iWord = iRow0 >> CHECK_SHIFT;
        int iBit = iRow0 - (iWord << CHECK_SHIFT);
        if (mark[iWord]) {
          mark[iWord] = static_cast< CoinCheckZero >(mark[iWord] | (1 << iBit));
        } else {
          mark[iWord] = static_cast< CoinCheckZero >(1 << iBit);
          stack[nMarked++] = iWord;
        }
        region[iRow0] = regionValue0 - value0 * pivotValue;
      }
      pivotValue *= pivotRegion[i];
      region[i] = pivotValue;
      regionIndex[numberNonZero++] = i;
    }
  }

  if (numberSlacks_) {
    // now do slacks
#ifndef COIN_FAST_CODE
    double factor = slackValue_;
    if (factor == 1.0) {
      // First do down to convenient power of 2
      int jLast = (numberSlacks_ - 1) >> CHECK_SHIFT;
      jLast = jLast << CHECK_SHIFT;
      for (i = numberSlacks_ - 1; i >= jLast; i--) {
        double value = region[i];
        double absValue = fabs(value);
        if (value) {
          region[i] = 0.0;
          if (absValue > tolerance) {
            region[i] = value;
            regionIndex[numberNonZero++] = i;
          }
        }
      }
      mark[jLast] = 0;
      // now do in chunks
      for (int k = (jLast >> CHECK_SHIFT) - 1; k >= 0; k--) {
        unsigned int iMark = mark[k];
        if (iMark) {
          // something in chunk - do all (as imark may change)
          int iLast = k << CHECK_SHIFT;
          i = iLast + BITS_PER_CHECK - 1;
          for (; i >= iLast; i--) {
            double value = region[i];
            double absValue = fabs(value);
            if (value) {
              region[i] = 0.0;
              if (absValue > tolerance) {
                region[i] = value;
                regionIndex[numberNonZero++] = i;
              }
            }
          }
          mark[k] = 0;
        }
      }
    } else {
      assert(factor == -1.0);
#endif
      // First do down to convenient power of 2
      int jLast = (numberSlacks_ - 1) >> CHECK_SHIFT;
      jLast = jLast << CHECK_SHIFT;
      for (i = numberSlacks_ - 1; i >= jLast; i--) {
        double value = region[i];
        double absValue = fabs(value);
        if (value) {
          region[i] = 0.0;
          if (absValue > tolerance) {
            region[i] = -value;
            regionIndex[numberNonZero++] = i;
          }
        }
      }
      mark[jLast] = 0;
      // now do in chunks
      for (int k = (jLast >> CHECK_SHIFT) - 1; k >= 0; k--) {
        unsigned int iMark = mark[k];
        if (iMark) {
          // something in chunk - do all (as imark may change)
          int iLast = k << CHECK_SHIFT;
          i = iLast + BITS_PER_CHECK - 1;
          for (; i >= iLast; i--) {
            double value = region[i];
            double absValue = fabs(value);
            if (value) {
              region[i] = 0.0;
              if (absValue > tolerance) {
                region[i] = -value;
                regionIndex[numberNonZero++] = i;
              }
            }
          }
          mark[k] = 0;
        }
      }
#ifndef COIN_FAST_CODE
    }
#endif
  }
  regionSparse->setNumElements(numberNonZero);
  mark[(numberU_ - 1) >> CHECK_SHIFT] = 0;
  mark[numberSlacks_ >> CHECK_SHIFT] = 0;
  if (numberSlacks_)
    mark[(numberSlacks_ - 1) >> CHECK_SHIFT] = 0;
#ifdef COIN_DEBUG
  for (i = 0; i < maximumRowsExtra_; i++) {
    assert(!mark[i]);
  }
#endif
}
//  updateColumnR.  Updates part of column (FTRANR)
void CoinFactorization::updateColumnR(CoinIndexedVector *regionSparse) const
{
  double *COIN_RESTRICT region = regionSparse->denseVector();
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  int numberNonZero = regionSparse->getNumElements();

  if (!numberR_)
    return; //return if nothing to do
  double tolerance = zeroTolerance_;

  const int *startColumn = startColumnRArray_ - numberRows_;
  const int *indexRow = indexRowR_;
  const CoinFactorizationDouble *COIN_RESTRICT element = elementR_;
  const int *permute = permuteArray_;

  // Work out very dubious idea of what would be fastest
  int method = -1;
  // Size of R
  double sizeR = startColumnRArray_[numberR_];
  // Average
  double averageR = sizeR / (static_cast< double >(numberRowsExtra_));
  // weights (relative to actual work)
  double setMark = 0.1; // setting mark
  double test1 = 1.0; // starting ftran (without testPivot)
  double testPivot = 2.0; // Seeing if zero etc
  double startDot = 2.0; // For starting dot product version
  // For final scan
  double final = numberNonZero * 1.0;
  double methodTime[3];
  // For second type
  methodTime[1] = numberPivots_ * (testPivot + ((static_cast< double >(numberNonZero)) / (static_cast< double >(numberRows_)) * averageR));
  methodTime[1] += numberNonZero * (test1 + averageR);
  // For first type
  methodTime[0] = methodTime[1] + (numberNonZero + numberPivots_) * setMark;
  methodTime[1] += numberNonZero * final;
  // third
  methodTime[2] = sizeR + numberPivots_ * startDot + numberNonZero * final;
  // switch off if necessary
  if (!numberInColumnPlus_.array()) {
    methodTime[0] = 1.0e100;
    methodTime[1] = 1.0e100;
  } else if (!sparseArray_) {
    methodTime[0] = 1.0e100;
  }
  double best = 1.0e100;
  for (int i = 0; i < 3; i++) {
    if (methodTime[i] < best) {
      best = methodTime[i];
      method = i;
    }
  }
  assert(method >= 0);
  const int *numberInColumnPlus = numberInColumnPlusArray_;
  //if (method==1)
  //printf(" methods %g %g %g - chosen %d\n",methodTime[0],methodTime[1],methodTime[2],method);

  switch (method) {
  case 0:
#ifdef STACK
  {
    // use sparse_ as temporary area
    // mark known to be zero
#if ABOCA_LITE_FACTORIZATION == 0
#define sparseOffset 0
#else
      int sparseOffset = ((regionSparse->capacity() & 0x80000000) != 0) ? sparseOffset_ : 0;
      assert(!sparseOffset);
#endif
    int *COIN_RESTRICT stack = sparseArray_ + sparseOffset; /* pivot */
    int *COIN_RESTRICT list = stack + maximumRowsExtra_; /* final list */
    int *COIN_RESTRICT next = (int *)(list + maximumRowsExtra_); /* jnext */
    char *COIN_RESTRICT mark = (char *)(next + maximumRowsExtra_);
    // we have another copy of R in R
    const CoinFactorizationDouble *COIN_RESTRICT elementR = elementR_ + lengthAreaR_;
    const int *indexRowR = indexRowR_ + lengthAreaR_;
    const int *startR = startColumnRArray_ + maximumPivots_ + 1;
    int nList = 0;
    const int *permuteBack = permuteBackArray_;
    for (int k = 0; k < numberNonZero; k++) {
      int kPivot = regionIndex[k];
      if (!mark[kPivot]) {
        stack[0] = kPivot;
        int j = -10;
        next[0] = j;
        int nStack = 0;
        while (nStack >= 0) {
          /* take off stack */
          if (j >= startR[kPivot]) {
            int jPivot = indexRowR[j--];
            /* put back on stack */
            next[nStack] = j;
            if (!mark[jPivot]) {
              /* and new one */
              kPivot = jPivot;
              j = -10;
              stack[++nStack] = kPivot;
              mark[kPivot] = 1;
              next[nStack] = j;
            }
          } else if (j == -10) {
            // before first - see if followon
            int jPivot = permuteBack[kPivot];
            if (jPivot < numberRows_) {
              // no
              j = startR[kPivot] + numberInColumnPlus[kPivot] - 1;
              next[nStack] = j;
            } else {
              // add to list
              if (!mark[jPivot]) {
                /* and new one */
                kPivot = jPivot;
                j = -10;
                stack[++nStack] = kPivot;
                mark[kPivot] = 1;
                next[nStack] = j;
              } else {
                j = startR[kPivot] + numberInColumnPlus[kPivot] - 1;
                next[nStack] = j;
              }
            }
          } else {
            // finished
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
      CoinFactorizationDouble pivotValue;
      if (iPivot < numberRows_) {
        pivotValue = region[iPivot];
      } else {
        int before = permute[iPivot];
        pivotValue = region[iPivot] + region[before];
        region[before] = 0.0;
      }
      if (fabs(pivotValue) > tolerance) {
        region[iPivot] = pivotValue;
        int start = startR[iPivot];
        int number = numberInColumnPlus[iPivot];
        int end = start + number;
        int j;
        for (j = start; j < end; j++) {
          int iRow = indexRowR[j];
          CoinFactorizationDouble value = elementR[j];
          region[iRow] -= value * pivotValue;
        }
        regionIndex[numberNonZero++] = iPivot;
      } else {
        region[iPivot] = 0.0;
      }
    }
  }
#else
  {

    // use sparse_ as temporary area
    // mark known to be zero
#if ABOCA_LITE_FACTORIZATION == 0
#define sparseOffset 0
#else
    int sparseOffset = ((regionSparse->capacity() & 0x80000000) != 0) ? sparseOffset_ : 0;
    assert(!sparseOffset);
#endif
    int *COIN_RESTRICT stack = sparseArray_ + sparseOffset; /* pivot */
    int *COIN_RESTRICT list = stack + maximumRowsExtra_; /* final list */
    int *COIN_RESTRICT next = reinterpret_cast< int * >(list + maximumRowsExtra_); /* jnext */
    char *COIN_RESTRICT mark = reinterpret_cast< char * >(next + maximumRowsExtra_);
    // mark all rows which will be permuted
    for (int i = numberRows_; i < numberRowsExtra_; i++) {
      int iRow = permute[i];
      mark[iRow] = 1;
    }
    // we have another copy of R in R
    const CoinFactorizationDouble *COIN_RESTRICT elementR = elementR_ + lengthAreaR_;
    const int *indexRowR = indexRowR_ + lengthAreaR_;
    const int *startR = startColumnRArray_ + maximumPivots_ + 1;
    // For current list order does not matter as
    // only affects end
    int newNumber = 0;
    for (int i = 0; i < numberNonZero; i++) {
      int iRow = regionIndex[i];
      assert(region[iRow]);
      if (!mark[iRow])
        regionIndex[newNumber++] = iRow;
      int number = numberInColumnPlus[iRow];
      if (number) {
        CoinFactorizationDouble pivotValue = region[iRow];
        int start = startR[iRow];
        int end = start + number;
        for (int j = start; j < end; j++) {
          CoinFactorizationDouble value = elementR[j];
          int jRow = indexRowR[j];
          region[jRow] -= pivotValue * value;
        }
      }
    }
    numberNonZero = newNumber;
    for (int i = numberRows_; i < numberRowsExtra_; i++) {
      //move using permute_ (stored in inverse fashion)
      int iRow = permute[i];
      CoinFactorizationDouble pivotValue = region[iRow] + region[i];
      //zero out pre-permuted
      region[iRow] = 0.0;
      if (fabs(pivotValue) > tolerance) {
        region[i] = pivotValue;
        if (!mark[i])
          regionIndex[numberNonZero++] = i;
        int number = numberInColumnPlus[i];
        int start = startR[i];
        int end = start + number;
        for (int j = start; j < end; j++) {
          CoinFactorizationDouble value = elementR[j];
          int jRow = indexRowR[j];
          region[jRow] -= pivotValue * value;
        }
      } else {
        region[i] = 0.0;
      }
      mark[iRow] = 0;
    }
  }
#endif
  break;
  case 1: {
    // no sparse region
    // we have another copy of R in R
    const CoinFactorizationDouble *COIN_RESTRICT elementR = elementR_ + lengthAreaR_;
    const int *indexRowR = indexRowR_ + lengthAreaR_;
    const int *startR = startColumnRArray_ + maximumPivots_ + 1;
    // For current list order does not matter as
    // only affects end
    for (int i = 0; i < numberNonZero; i++) {
      int iRow = regionIndex[i];
      assert(region[iRow]);
      int number = numberInColumnPlus[iRow];
      if (number) {
        CoinFactorizationDouble pivotValue = region[iRow];
        int start = startR[iRow];
        int end = start + number;
        for (int j = start; j < end; j++) {
          CoinFactorizationDouble value = elementR[j];
          int jRow = indexRowR[j];
          region[jRow] -= pivotValue * value;
        }
      }
    }
    for (int i = numberRows_; i < numberRowsExtra_; i++) {
      //move using permute_ (stored in inverse fashion)
      int iRow = permute[i];
      CoinFactorizationDouble pivotValue = region[iRow] + region[i];
      //zero out pre-permuted
      region[iRow] = 0.0;
      if (fabs(pivotValue) > tolerance) {
        region[i] = pivotValue;
        regionIndex[numberNonZero++] = i;
        int number = numberInColumnPlus[i];
        int start = startR[i];
        int end = start + number;
        for (int j = start; j < end; j++) {
          CoinFactorizationDouble value = elementR[j];
          int jRow = indexRowR[j];
          region[jRow] -= pivotValue * value;
        }
      } else {
        region[i] = 0.0;
      }
    }
  } break;
  case 2: {
    int start = startColumn[numberRows_];
    for (int i = numberRows_; i < numberRowsExtra_; i++) {
      //move using permute_ (stored in inverse fashion)
      int end = startColumn[i + 1];
      int iRow = permute[i];
      CoinFactorizationDouble pivotValue = region[iRow];
      //zero out pre-permuted
      region[iRow] = 0.0;

      for (int j = start; j < end; j++) {
        CoinFactorizationDouble value = element[j];
        int jRow = indexRow[j];
        value *= region[jRow];
        pivotValue -= value;
      }
      start = end;
      if (fabs(pivotValue) > tolerance) {
        region[i] = pivotValue;
        regionIndex[numberNonZero++] = i;
      } else {
        region[i] = 0.0;
      }
    }
  } break;
  }
  if (method) {
    // pack down
    int n = numberNonZero;
    numberNonZero = 0;
    for (int i = 0; i < n; i++) {
      int indexValue = regionIndex[i];
      double value = region[indexValue];
      if (value)
        regionIndex[numberNonZero++] = indexValue;
    }
  }
  //set counts
  regionSparse->setNumElements(numberNonZero);
}
//  updateColumnR.  Updates part of column (FTRANR)
void CoinFactorization::updateColumnRFT(CoinIndexedVector *regionSparse,
  int *COIN_RESTRICT regionIndex)
{
  double *COIN_RESTRICT region = regionSparse->denseVector();
  //int *regionIndex = regionSparse->getIndices (  );
  int *COIN_RESTRICT startColumnU = startColumnUArray_;
  int numberNonZero = regionSparse->getNumElements();

  if (numberR_) {
    double tolerance = zeroTolerance_;

    const int *startColumn = startColumnRArray_ - numberRows_;
    const int *indexRow = indexRowR_;
    const CoinFactorizationDouble *COIN_RESTRICT element = elementR_;
    const int *permute = permuteArray_;

    // Work out very dubious idea of what would be fastest
    int method = -1;
    // Size of R
    double sizeR = startColumnRArray_[numberR_];
    // Average
    double averageR = sizeR / (static_cast< double >(numberRowsExtra_));
    // weights (relative to actual work)
    double setMark = 0.1; // setting mark
    double test1 = 1.0; // starting ftran (without testPivot)
    double testPivot = 2.0; // Seeing if zero etc
    double startDot = 2.0; // For starting dot product version
    // For final scan
    double final = numberNonZero * 1.0;
    double methodTime[3];
    // For second type
    methodTime[1] = numberPivots_ * (testPivot + ((static_cast< double >(numberNonZero)) / (static_cast< double >(numberRows_)) * averageR));
    methodTime[1] += numberNonZero * (test1 + averageR);
    // For first type
    methodTime[0] = methodTime[1] + (numberNonZero + numberPivots_) * setMark;
    methodTime[1] += numberNonZero * final;
    // third
    methodTime[2] = sizeR + numberPivots_ * startDot + numberNonZero * final;
    // switch off if necessary
    if (!numberInColumnPlus_.array()) {
      methodTime[0] = 1.0e100;
      methodTime[1] = 1.0e100;
    } else if (!sparseArray_) {
      methodTime[0] = 1.0e100;
    }
    const int *numberInColumnPlus = numberInColumnPlusArray_;
    int *COIN_RESTRICT numberInColumn = numberInColumnArray_;
    // adjust for final scan
    methodTime[1] += final;
    double best = 1.0e100;
    for (int i = 0; i < 3; i++) {
      if (methodTime[i] < best) {
        best = methodTime[i];
        method = i;
      }
    }
    assert(method >= 0);

    switch (method) {
    case 0: {
      // use sparse_ as temporary area
      // mark known to be zero
#if ABOCA_LITE_FACTORIZATION == 0
#define sparseOffset 0
#else
      int sparseOffset = ((regionSparse->capacity() & 0x80000000) != 0) ? sparseOffset_ : 0;
      assert(!sparseOffset);
#endif
      int *COIN_RESTRICT stack = sparseArray_ + sparseOffset; /* pivot */
      int *COIN_RESTRICT list = stack + maximumRowsExtra_; /* final list */
      int *COIN_RESTRICT next = reinterpret_cast< int * >(list + maximumRowsExtra_); /* jnext */
      char *COIN_RESTRICT mark = reinterpret_cast< char * >(next + maximumRowsExtra_);
      // mark all rows which will be permuted
      for (int i = numberRows_; i < numberRowsExtra_; i++) {
        int iRow = permute[i];
        mark[iRow] = 1;
      }
      // we have another copy of R in R
      const CoinFactorizationDouble *COIN_RESTRICT elementR = elementR_ + lengthAreaR_;
      const int *indexRowR = indexRowR_ + lengthAreaR_;
      const int *startR = startColumnRArray_ + maximumPivots_ + 1;
      //save in U
      //in at end
      int iColumn = numberColumnsExtra_;

      startColumnU[iColumn] = startColumnU[maximumColumnsExtra_];
      int start = startColumnU[iColumn];

      //int * putIndex = indexRowU_ + start;
      CoinFactorizationDouble *COIN_RESTRICT putElement = elementUArray_ + start;
      // For current list order does not matter as
      // only affects end
      int newNumber = 0;
      for (int i = 0; i < numberNonZero; i++) {
        int iRow = regionIndex[i];
        CoinFactorizationDouble pivotValue = region[iRow];
        assert(region[iRow]);
        if (!mark[iRow]) {
          //putIndex[newNumber]=iRow;
          putElement[newNumber] = pivotValue;
          ;
          regionIndex[newNumber++] = iRow;
        }
        int number = numberInColumnPlus[iRow];
        if (number) {
          int start = startR[iRow];
          int end = start + number;
          for (int j = start; j < end; j++) {
            CoinFactorizationDouble value = elementR[j];
            int jRow = indexRowR[j];
            region[jRow] -= pivotValue * value;
          }
        }
      }
      numberNonZero = newNumber;
      for (int i = numberRows_; i < numberRowsExtra_; i++) {
        //move using permute_ (stored in inverse fashion)
        int iRow = permute[i];
        CoinFactorizationDouble pivotValue = region[iRow] + region[i];
        //zero out pre-permuted
        region[iRow] = 0.0;
        if (fabs(pivotValue) > tolerance) {
          region[i] = pivotValue;
          if (!mark[i]) {
            //putIndex[numberNonZero]=i;
            putElement[numberNonZero] = pivotValue;
            ;
            regionIndex[numberNonZero++] = i;
          }
          int number = numberInColumnPlus[i];
          int start = startR[i];
          int end = start + number;
          for (int j = start; j < end; j++) {
            CoinFactorizationDouble value = elementR[j];
            int jRow = indexRowR[j];
            region[jRow] -= pivotValue * value;
          }
        } else {
          region[i] = 0.0;
        }
        mark[iRow] = 0;
      }
      numberInColumn[iColumn] = numberNonZero;
      startColumnU[maximumColumnsExtra_] = start + numberNonZero;
    } break;
    case 1: {
      // no sparse region
      // we have another copy of R in R
      const CoinFactorizationDouble *COIN_RESTRICT elementR = elementR_ + lengthAreaR_;
      const int *indexRowR = indexRowR_ + lengthAreaR_;
      const int *startR = startColumnRArray_ + maximumPivots_ + 1;
      // For current list order does not matter as
      // only affects end
      for (int i = 0; i < numberNonZero; i++) {
        int iRow = regionIndex[i];
        assert(region[iRow]);
        int number = numberInColumnPlus[iRow];
        if (number) {
          CoinFactorizationDouble pivotValue = region[iRow];
          int start = startR[iRow];
          int end = start + number;
          for (int j = start; j < end; j++) {
            CoinFactorizationDouble value = elementR[j];
            int jRow = indexRowR[j];
            region[jRow] -= pivotValue * value;
          }
        }
      }
      for (int i = numberRows_; i < numberRowsExtra_; i++) {
        //move using permute_ (stored in inverse fashion)
        int iRow = permute[i];
        CoinFactorizationDouble pivotValue = region[iRow] + region[i];
        //zero out pre-permuted
        region[iRow] = 0.0;
        if (fabs(pivotValue) > tolerance) {
          region[i] = pivotValue;
          regionIndex[numberNonZero++] = i;
          int number = numberInColumnPlus[i];
          int start = startR[i];
          int end = start + number;
          for (int j = start; j < end; j++) {
            CoinFactorizationDouble value = elementR[j];
            int jRow = indexRowR[j];
            region[jRow] -= pivotValue * value;
          }
        } else {
          region[i] = 0.0;
        }
      }
    } break;
    case 2: {
      int start = startColumn[numberRows_];
      for (int i = numberRows_; i < numberRowsExtra_; i++) {
        //move using permute_ (stored in inverse fashion)
        int end = startColumn[i + 1];
        int iRow = permute[i];
        CoinFactorizationDouble pivotValue = region[iRow];
        //zero out pre-permuted
        region[iRow] = 0.0;

        for (int j = start; j < end; j++) {
          CoinFactorizationDouble value = element[j];
          int jRow = indexRow[j];
          value *= region[jRow];
          pivotValue -= value;
        }
        start = end;
        if (fabs(pivotValue) > tolerance) {
          region[i] = pivotValue;
          regionIndex[numberNonZero++] = i;
        } else {
          region[i] = 0.0;
        }
      }
    } break;
    }
    if (method) {
      // pack down
      int n = numberNonZero;
      numberNonZero = 0;
      //save in U
      //in at end
      int iColumn = numberColumnsExtra_;

      assert(startColumnU[iColumn] == startColumnU[maximumColumnsExtra_]);
      int start = startColumnU[iColumn];

      int *COIN_RESTRICT putIndex = indexRowUArray_ + start;
      CoinFactorizationDouble *COIN_RESTRICT putElement = elementUArray_ + start;
      for (int i = 0; i < n; i++) {
        int indexValue = regionIndex[i];
        double value = region[indexValue];
        if (value) {
          putIndex[numberNonZero] = indexValue;
          putElement[numberNonZero] = value;
          regionIndex[numberNonZero++] = indexValue;
        }
      }
      numberInColumn[iColumn] = numberNonZero;
      startColumnU[maximumColumnsExtra_] = start + numberNonZero;
    }
    //set counts
    regionSparse->setNumElements(numberNonZero);
  } else {
    // No R but we still need to save column
    //save in U
    //in at end
    int *COIN_RESTRICT numberInColumn = numberInColumnArray_;
    numberNonZero = regionSparse->getNumElements();
    int iColumn = numberColumnsExtra_;

    assert(startColumnU[iColumn] == startColumnU[maximumColumnsExtra_]);
    int start = startColumnU[iColumn];
    numberInColumn[iColumn] = numberNonZero;
    startColumnU[maximumColumnsExtra_] = start + numberNonZero;

    int *COIN_RESTRICT putIndex = indexRowUArray_ + start;
    CoinFactorizationDouble *COIN_RESTRICT putElement = elementUArray_ + start;
    for (int i = 0; i < numberNonZero; i++) {
      int indexValue = regionIndex[i];
      double value = region[indexValue];
      putIndex[i] = indexValue;
      putElement[i] = value;
    }
  }
}
/* Updates one column (FTRAN) from region2 and permutes.
   region1 starts as zero
   Note - if regionSparse2 packed on input - will be packed on output
   - returns un-permuted result in region2 and region1 is zero */
int CoinFactorization::updateColumnFT(CoinIndexedVector *regionSparse,
  CoinIndexedVector *regionSparse2)
{
#ifdef CLP_FACTORIZATION_INSTRUMENT
  double startTimeX = CoinCpuTime();
#endif
  //permute and move indices into index array
  int *COIN_RESTRICT regionIndex = regionSparse->getIndices();
  int numberNonZero = regionSparse2->getNumElements();
  const int *permute = permuteArray_;
  int *COIN_RESTRICT index = regionSparse2->getIndices();
  double *COIN_RESTRICT region = regionSparse->denseVector();
  double *COIN_RESTRICT array = regionSparse2->denseVector();
  int *COIN_RESTRICT startColumnU = startColumnUArray_;
  bool doFT = doForrestTomlin_;
  // see if room
  if (doFT) {
    int iColumn = numberColumnsExtra_;

    startColumnU[iColumn] = startColumnU[maximumColumnsExtra_];
    int start = startColumnU[iColumn];
    int space = lengthAreaU_ - (start + numberRowsExtra_);
    doFT = space >= 0;
    if (doFT) {
      regionIndex = indexRowUArray_ + start;
    } else {
      startColumnU[maximumColumnsExtra_] = lengthAreaU_ + 1;
    }
  }

#ifndef CLP_FACTORIZATION
  bool packed = regionSparse2->packedMode();
  if (packed) {
#else
  assert(regionSparse2->packedMode());
#endif
    for (int j = 0; j < numberNonZero; j++) {
      int iRow = index[j];
      double value = array[j];
      array[j] = 0.0;
      iRow = permute[iRow];
      region[iRow] = value;
      regionIndex[j] = iRow;
    }
#ifndef CLP_FACTORIZATION
  } else {
    for (int j = 0; j < numberNonZero; j++) {
      int iRow = index[j];
      double value = array[iRow];
      array[iRow] = 0.0;
      iRow = permute[iRow];
      region[iRow] = value;
      regionIndex[j] = iRow;
    }
  }
#endif
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
  //  ******* U
  updateColumnU(regionSparse, regionIndex);
  if (!doForrestTomlin_) {
    // Do PFI after everything else
    updateColumnPFI(regionSparse);
  }
  permuteBack(regionSparse, regionSparse2);
#ifdef CLP_FACTORIZATION_INSTRUMENT
  numberUpdateFT++;
  timeInUpdateFT += CoinCpuTime() - startTimeX;
  averageLengthR += lengthR_;
  averageLengthU += lengthU_;
  averageLengthL += lengthL_;
#endif
  // will be negative if no room
  if (doFT)
    return regionSparse2->getNumElements();
  else
    return -regionSparse2->getNumElements();
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
