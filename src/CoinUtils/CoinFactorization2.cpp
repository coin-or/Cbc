// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include "CoinPragma.hpp"
#include "CoinUtilsConfig.h"

#include <atomic>
#include <cassert>
#include <cfloat>
#include <stdio.h>
#include "CoinFactorization.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinFinite.hpp"
#if COIN_FACTORIZATION_DENSE_CODE == 1
// using simple lapack interface

extern "C" 
{
  /** LAPACK Fortran subroutine DGETRF. */
  void COINUTILS_LAPACK_FUNC(dgetrf,DGETRF)(ipfint *m, ipfint *n,
                                       double *A, ipfint *ldA,
                                       ipfint *ipiv, ipfint *info);
}
#elif COIN_FACTORIZATION_DENSE_CODE == 2
// C interface
enum CBLAS_ORDER { CblasRowMajor = 101,
  CblasColMajor = 102 };
extern "C" {
int clapack_dgetrf(const enum CBLAS_ORDER Order, const int M, const int N, double *A, const int lda, int *ipiv);
}
#elif COIN_FACTORIZATION_DENSE_CODE == 3
// Intel compiler
#include "mkl_lapacke.h"
#endif
#ifndef NDEBUG
static std::atomic<int> counter1(0);
#endif
//  factorSparse.  Does sparse phase of factorization
//return code is <0 error, 0= finished
int CoinFactorization::factorSparse()
{
  int larger;

  if (numberRows_ < numberColumns_) {
    larger = numberColumns_;
  } else {
    larger = numberRows_;
  }
  int returnCode;
#define LARGELIMIT 65530
#define SMALL_SET 65531
#define SMALL_UNSET (SMALL_SET + 1)
#define LARGE_SET COIN_INT_MAX - 10
#define LARGE_UNSET (LARGE_SET + 1)
  if (larger < LARGELIMIT)
    returnCode = factorSparseSmall();
  else
    returnCode = factorSparseLarge();
  return returnCode;
}
//  factorSparse.  Does sparse phase of factorization
//return code is <0 error, 0= finished
int CoinFactorization::factorSparseSmall()
{
  int *COIN_RESTRICT indexRow = indexRowUArray_;
  int *COIN_RESTRICT indexColumn = indexColumnUArray_;
  CoinFactorizationDouble *COIN_RESTRICT element = elementUArray_;
  int count = 1;
  workArea_.conditionalNew(numberRows_);
  CoinFactorizationDouble *COIN_RESTRICT workArea = workArea_.array();
  workAreaArray_ = workArea;
#ifndef NDEBUG
  counter1++;
#endif
  // when to go dense
  int denseThreshold = abs(denseThreshold_);

  CoinZeroN(workArea, numberRows_);
  //get space for bit work area
  int workSize = 1000;
  workArea2_.conditionalNew(workSize);
  unsigned int *COIN_RESTRICT workArea2 = workArea2_.array();
  workArea2Array_ = workArea2;

  //set markRow so no rows updated
  unsigned short *markRow = reinterpret_cast< unsigned short * >(markRowArray_);
  CoinFillN(markRow, numberRows_, static_cast< unsigned short >(SMALL_UNSET));
  int status = 0;
  //do slacks first
  int *COIN_RESTRICT numberInRow = numberInRowArray_;
  int *COIN_RESTRICT numberInColumn = numberInColumnArray_;
  int *COIN_RESTRICT numberInColumnPlus = numberInColumnPlusArray_;
  int *COIN_RESTRICT startColumnU = startColumnUArray_;
  int *COIN_RESTRICT startColumnL = startColumnLArray_;
  if (biasLU_ < 3 && numberColumns_ == numberRows_) {
    int iPivotColumn;
    int *COIN_RESTRICT pivotColumn = pivotColumnArray_;
    int *COIN_RESTRICT nextRow = nextRowArray_;
    int *COIN_RESTRICT lastRow = lastRowArray_;
    for (iPivotColumn = 0; iPivotColumn < numberColumns_;
         iPivotColumn++) {
      if (numberInColumn[iPivotColumn] == 1) {
        int start = startColumnU[iPivotColumn];
        CoinFactorizationDouble value = element[start];
        if (value == slackValue_ && numberInColumnPlus[iPivotColumn] == 0) {
          // treat as slack
          int iRow = indexRow[start];
          // but only if row not marked
          if (numberInRow[iRow] > 0) {
            totalElements_ -= numberInRow[iRow];
            //take out this bit of indexColumnU
            int next = nextRow[iRow];
            int last = lastRow[iRow];

            nextRow[last] = next;
            lastRow[next] = last;
            nextRow[iRow] = numberGoodU_; //use for permute
            lastRow[iRow] = -2; //mark
            //modify linked list for pivots
            deleteLink(iRow);
            numberInRow[iRow] = -1;
            numberInColumn[iPivotColumn] = 0;
            numberGoodL_++;
            startColumnL[numberGoodL_] = 0;
            pivotColumn[numberGoodU_] = iPivotColumn;
            numberGoodU_++;
          }
        }
      }
    }
    // redo
    preProcess(4);
    CoinFillN(markRow, numberRows_, static_cast< unsigned short >(SMALL_UNSET));
  }
  numberSlacks_ = numberGoodU_;
  int *COIN_RESTRICT nextCount = nextCountArray_;
  int *COIN_RESTRICT firstCount = firstCountArray_;
  int *COIN_RESTRICT startRow = startRowUArray_;
  int *COIN_RESTRICT startColumn = startColumnU;
  //#define UGLY_COIN_FACTOR_CODING
#ifdef UGLY_COIN_FACTOR_CODING
  CoinFactorizationDouble *COIN_RESTRICT elementL = elementLArray_;
  int *COIN_RESTRICT indexRowL = indexRowLArray_;
  int *COIN_RESTRICT saveColumn = saveColumnArray_;
  int *COIN_RESTRICT nextRow = nextRowArray_;
  int *COIN_RESTRICT lastRow = lastRowArray_;
#endif
  double pivotTolerance = pivotTolerance_;
  int numberTrials = numberTrials_;
  int numberRows = numberRows_;
  // Put column singletons first - (if false)
  separateLinks(1, (biasLU_ > 1));
#ifndef NDEBUG
  int counter2 = 0;
#endif
  while (count <= biggerDimension_) {
#ifndef NDEBUG
    counter2++;
    int badRow = -1;
    if (counter1 == -1 && counter2 >= 0) {
      // check counts consistent
      for (int iCount = 1; iCount < numberRows_; iCount++) {
        int look = firstCount[iCount];
        while (look >= 0) {
          if (look < numberRows_) {
            int iRow = look;
            if (iRow == badRow)
              printf("row count for row %d is %d\n", iCount, iRow);
            if (numberInRow[iRow] != iCount) {
              printf("failed debug on %d entry to factorSparse and %d try\n",
                (int)counter1, counter2);
              printf("row %d - count %d number %d\n", iRow, iCount, numberInRow[iRow]);
              abort();
            }
            look = nextCount[look];
          } else {
            int iColumn = look - numberRows;
            if (numberInColumn[iColumn] != iCount) {
              printf("failed debug on %d entry to factorSparse and %d try\n",
                (int)counter1, counter2);
              printf("column %d - count %d number %d\n", iColumn, iCount, numberInColumn[iColumn]);
              abort();
            }
            look = nextCount[look];
          }
        }
      }
    }
#endif
    int minimumCount = COIN_INT_MAX;
    double minimumCost = COIN_DBL_MAX;

    int iPivotRow = -1;
    int iPivotColumn = -1;
    int pivotRowPosition = -1;
    int pivotColumnPosition = -1;
    int look = firstCount[count];
    int trials = 0;
    int *COIN_RESTRICT pivotColumn = pivotColumnArray_;

    if (count == 1 && firstCount[1] >= 0 && !biasLU_) {
      //do column singletons first to put more in U
      while (look >= 0) {
        if (look < numberRows_) {
          look = nextCount[look];
        } else {
          int iColumn = look - numberRows_;

          assert(numberInColumn[iColumn] == count);
          int start = startColumnU[iColumn];
          int iRow = indexRow[start];

          iPivotRow = iRow;
          pivotRowPosition = start;
          iPivotColumn = iColumn;
          assert(iPivotRow >= 0 && iPivotColumn >= 0);
          pivotColumnPosition = -1;
          look = -1;
          break;
        }
      } /* endwhile */
      if (iPivotRow < 0) {
        //back to singletons
        look = firstCount[1];
      }
    }
    while (look >= 0) {
      if (look < numberRows_) {
        int iRow = look;
#ifndef NDEBUG
        if (numberInRow[iRow] != count) {
          printf("failed on %d entry to factorSparse and %d try\n",
            (int)counter1, counter2);
          printf("row %d - count %d number %d\n", iRow, count, numberInRow[iRow]);
          abort();
        }
#endif
        look = nextCount[look];
        bool rejected = false;
        int start = startRow[iRow];
        int end = start + count;

        int i;
        for (i = start; i < end; i++) {
          int iColumn = indexColumn[i];
          assert(numberInColumn[iColumn] > 0);
          double cost = (count - 1) * numberInColumn[iColumn];

          if (cost < minimumCost) {
            int where = startColumn[iColumn];
            double minimumValue = element[where];

            minimumValue = fabs(minimumValue) * pivotTolerance;
            while (indexRow[where] != iRow) {
              where++;
            } /* endwhile */
            assert(where < startColumn[iColumn] + numberInColumn[iColumn]);
            CoinFactorizationDouble value = element[where];

            value = fabs(value);
            if (value >= minimumValue) {
              minimumCost = cost;
              minimumCount = numberInColumn[iColumn];
              iPivotRow = iRow;
              pivotRowPosition = -1;
              iPivotColumn = iColumn;
              assert(iPivotRow >= 0 && iPivotColumn >= 0);
              pivotColumnPosition = i;
              rejected = false;
              if (minimumCount < count) {
                look = -1;
                break;
              }
            } else if (iPivotRow == -1) {
              rejected = true;
            }
          }
        }
        trials++;
        if (trials >= numberTrials && iPivotRow >= 0) {
          look = -1;
          break;
        }
        if (rejected) {
          //take out for moment
          //eligible when row changes
          deleteLink(iRow);
          addLink(iRow, biggerDimension_ + 1);
        }
      } else {
        int iColumn = look - numberRows;

        assert(numberInColumn[iColumn] == count);
        look = nextCount[look];
        int start = startColumn[iColumn];
        int end = start + numberInColumn[iColumn];
        CoinFactorizationDouble minimumValue = element[start];

        minimumValue = fabs(minimumValue) * pivotTolerance;
        int i;
        for (i = start; i < end; i++) {
          CoinFactorizationDouble value = element[i];

          value = fabs(value);
          if (value >= minimumValue) {
            int iRow = indexRow[i];
            int nInRow = numberInRow[iRow];
            assert(nInRow > 0);
            double cost = (count - 1) * nInRow;

            if (cost < minimumCost) {
              minimumCost = cost;
              minimumCount = nInRow;
              iPivotRow = iRow;
              pivotRowPosition = i;
              iPivotColumn = iColumn;
              assert(iPivotRow >= 0 && iPivotColumn >= 0);
              pivotColumnPosition = -1;
              if (minimumCount <= count + 1) {
                look = -1;
                break;
              }
            }
          }
        }
        trials++;
        if (trials >= numberTrials && iPivotRow >= 0) {
          look = -1;
          break;
        }
      }
    } /* endwhile */
    if (iPivotRow >= 0) {
      assert(iPivotRow < numberRows_);
      int numberDoRow = numberInRow[iPivotRow] - 1;
      int numberDoColumn = numberInColumn[iPivotColumn] - 1;

      totalElements_ -= (numberDoRow + numberDoColumn + 1);
      if (numberDoColumn > 0) {
        if (numberDoRow > 0) {
          if (numberDoColumn > 1) {
            //  if (1) {
            //need to adjust more for cache and SMP
            //allow at least 4 extra
            int increment = numberDoColumn + 1 + 4;

            if (increment & 15) {
              increment = increment & (~15);
              increment += 16;
            }
            int increment2 =

              (increment + COINFACTORIZATION_BITS_PER_INT - 1) >> COINFACTORIZATION_SHIFT_PER_INT;
            int size = increment2 * numberDoRow;

            if (size > workSize) {
              workSize = size;
              workArea2_.conditionalNew(workSize);
              workArea2 = workArea2_.array();
	      workArea2Array_ = workArea2;
            }
            bool goodPivot;
#ifndef UGLY_COIN_FACTOR_CODING
            //branch out to best pivot routine
            goodPivot = pivot(iPivotRow, iPivotColumn,
              pivotRowPosition, pivotColumnPosition,
              workArea, workArea2,
              increment2, markRow,
              SMALL_SET);
#else
#define FAC_SET SMALL_SET
#include "CoinFactorization.hpp"
#undef FAC_SET
#undef UGLY_COIN_FACTOR_CODING
#endif
            if (!goodPivot) {
              status = -99;
              count = biggerDimension_ + 1;
              break;
            }
          } else {
            if (!pivotOneOtherRow(iPivotRow, iPivotColumn)) {
              status = -99;
              count = biggerDimension_ + 1;
              break;
            }
          }
        } else {
          assert(!numberDoRow);
          if (!pivotRowSingleton(iPivotRow, iPivotColumn)) {
            status = -99;
            count = biggerDimension_ + 1;
            break;
          }
        }
      } else {
        assert(!numberDoColumn);
        if (!pivotColumnSingleton(iPivotRow, iPivotColumn)) {
          status = -99;
          count = biggerDimension_ + 1;
          break;
        }
      }
      assert(nextRow_.array()[iPivotRow] == numberGoodU_);
      pivotColumn[numberGoodU_] = iPivotColumn;
      numberGoodU_++;
      // This should not need to be trapped here - but be safe
      if (numberGoodU_ == numberRows_)
        count = biggerDimension_ + 1;
#if COIN_DEBUG == 2
      checkConsistency();
#endif
      if (denseThreshold) {
        // see whether to go dense
        int leftRows = numberRows_ - numberGoodU_;
        double full = leftRows;
        full *= full;
        assert(full >= 0.0);
        double leftElements = totalElements_;
        //if (leftRows==100)
        //printf("at 100 %d elements\n",totalElements_);
        double ratio;
#ifndef COIN_DENSE_MULTIPLIER
#define COIN_DENSE_MULTIPLIER 1.0
#endif
#define COIN_DENSE_MULTIPLIER2 1
        if (leftRows > 2000) {
          ratio = 4.0;
#if COIN_DENSE_MULTIPLIER2 == 1
          ratio = 3.5;
#endif
        } else if (leftRows > 800) {
          ratio = 3.0;
#if COIN_DENSE_MULTIPLIER2 == 1
          ratio = 2.75;
#endif
        } else if (leftRows > 256) {
          ratio = 2.0;
        } else {
          ratio = 1.5;
        }
#if COIN_DENSE_MULTIPLIER2 > 10
        ratio = 10000;
#endif
        ratio *= COIN_DENSE_MULTIPLIER;
        if ((ratio * leftElements > full && leftRows > denseThreshold)) {
#define COIN_ALIGN_DENSE 2
#if COIN_ALIGN_DENSE == 2
          if ((leftRows & 7) == 0) {
#endif
            //return to do dense
            if (status != 0)
              break;
            int check = 0;
            for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
              if (numberInColumn[iColumn])
                check++;
            }
            if (check != leftRows && denseThreshold) {
              //printf("** mismatch %d columns left, %d rows\n",check,leftRows);
              denseThreshold = 0;
            } else {
              status = 2;
              if ((messageLevel_ & 4) != 0)
                std::cout << "      Went dense at " << leftRows << " rows " << totalElements_ << " " << full << " " << leftElements << std::endl;
              //if (!denseThreshold_)
              //denseThreshold_=-check; // say how many
              break;
            }
#if COIN_ALIGN_DENSE == 2
          }
#endif
        }
      }
      // start at 1 again
      count = 1;
    } else {
      //end of this - onto next
      count++;
    }
  } /* endwhile */
  workArea_.conditionalDelete();
  workArea2_.conditionalDelete();
  return status;
}

//:method factorDense.  Does dense phase of factorization
//return code is <0 error, 0= finished
int CoinFactorization::factorDense()
{
  int status = 0;
  numberDense_ = numberRows_ - numberGoodU_;
  if (sizeof(int) == 4 && numberDense_ >= 2 << 15) {
    abort();
  }
  int full;
  if (denseThreshold_ > 0 || true)
    full = numberDense_ * numberDense_;
  else
    full = -denseThreshold_ * numberDense_;
  totalElements_ = full;
#ifdef COIN_ALIGN_DENSE
  int newSize = full + 8 * numberDense_;
  newSize += (numberDense_ + 1) / static_cast<int>(sizeof(CoinFactorizationDouble) / sizeof(int));
  newSize += 2 * ((numberDense_ + 3) / static_cast<int>(sizeof(CoinFactorizationDouble) / sizeof(short)));
  newSize += ((numberRows_ + 3) / static_cast<int>(sizeof(CoinFactorizationDouble) / sizeof(short)));
  // so we can align on 256 byte
  newSize += 32;
  denseArea_ = new double[newSize];
  denseAreaAddress_ = denseArea_;
  CoinInt64 xx = reinterpret_cast< CoinInt64 >(denseAreaAddress_);
  int iBottom = static_cast< int >(xx & 63);
  int offset = (256 - iBottom) >> 3;
  denseAreaAddress_ += offset;
  CoinZeroN(denseArea_, newSize);
#else
  denseArea_ = new double[full];
  denseAreaAddress_ = denseArea_;
  CoinZeroN(denseArea_, full);
#endif
  densePermute_ = new int[numberDense_];
  int *COIN_RESTRICT indexRowU = indexRowUArray_;
  //mark row lookup using lastRow
  int i;
  int *COIN_RESTRICT nextRow = nextRowArray_;
  int *COIN_RESTRICT lastRow = lastRowArray_;
  int *COIN_RESTRICT numberInColumn = numberInColumnArray_;
  int *COIN_RESTRICT numberInColumnPlus = numberInColumnPlusArray_;
  for (i = 0; i < numberRows_; i++) {
    if (lastRow[i] >= 0)
      lastRow[i] = 0;
  }
  int *COIN_RESTRICT indexRow = indexRowUArray_;
  CoinFactorizationDouble *COIN_RESTRICT element = elementUArray_;
  int which = 0;
  for (i = 0; i < numberRows_; i++) {
    if (!lastRow[i]) {
      lastRow[i] = which;
      nextRow[i] = numberGoodU_ + which;
      densePermute_[which] = i;
      which++;
    }
  }
  //for L part
  int *COIN_RESTRICT startColumnL = startColumnLArray_;
  CoinFactorizationDouble *COIN_RESTRICT elementL = elementLArray_;
  int *COIN_RESTRICT indexRowL = indexRowLArray_;
  int endL = startColumnL[numberGoodL_];
  //take out of U
  double *COIN_RESTRICT column = denseAreaAddress_;
  int rowsDone = 0;
  int iColumn = 0;
  int *COIN_RESTRICT pivotColumn = pivotColumnArray_;
  CoinFactorizationDouble *COIN_RESTRICT pivotRegion = pivotRegionArray_;
  int *COIN_RESTRICT startColumnU = startColumnUArray_;
  for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
    if (numberInColumn[iColumn]) {
      //move
      int start = startColumnU[iColumn];
      int number = numberInColumn[iColumn];
      int end = start + number;
      for (int i = start; i < end; i++) {
        int iRow = indexRow[i];
        iRow = lastRow[iRow];
        assert(iRow >= 0 && iRow < numberDense_);
        column[iRow] = element[i];
      } /* endfor */
      column += numberDense_;
      while (lastRow[rowsDone] < 0) {
        rowsDone++;
      } /* endwhile */
      nextRow[rowsDone] = numberGoodU_;
      rowsDone++;
      startColumnL[numberGoodU_ + 1] = endL;
      numberInColumn[iColumn] = 0;
      pivotColumn[numberGoodU_] = iColumn;
      pivotRegion[numberGoodU_] = 1.0;
      numberGoodU_++;
    }
  }
#ifdef COIN_FACTORIZATION_DENSE_CODE
  if (denseThreshold_ /*>0*/) {
    assert(numberGoodU_ == numberRows_);
    numberGoodL_ = numberRows_;
    //now factorize
    //dgef(denseAreaAddress_,&numberDense_,&numberDense_,densePermute_);
#if COIN_FACTORIZATION_DENSE_CODE == 1
    int info;
    COINUTILS_LAPACK_FUNC(dgetrf,DGETRF)(&numberDense_,&numberDense_,
                     denseAreaAddress_,&numberDense_,densePermute_,&info);
    // need to check size of pivots
    if (info) {
      //printf("Dense singular\n");
      status = -1;
    }
#elif COIN_FACTORIZATION_DENSE_CODE == 2
    status = clapack_dgetrf(CblasColMajor, numberDense_, numberDense_,
      denseAreaAddress_, numberDense_, densePermute_);
#elif COIN_FACTORIZATION_DENSE_CODE == 3
    status = LAPACKE_dgetrf(LAPACK_COL_MAJOR, numberDense_, numberDense_,
      denseAreaAddress_, numberDense_, densePermute_);
#endif
    return status;
  }
#endif
  //abort();
  numberGoodU_ = numberRows_ - numberDense_;
  int base = numberGoodU_;
  int iDense;
  int numberToDo = -denseThreshold_;
  denseThreshold_ = 0;
  double tolerance = zeroTolerance_;
  tolerance = 1.0e-30;
  int *COIN_RESTRICT nextColumn = nextColumnArray_;
  const int *pivotColumnConst = pivotColumnArray_;
  // make sure we have enough space in L and U
  for (iDense = 0; iDense < numberToDo; iDense++) {
    //how much space have we got
    iColumn = pivotColumnConst[base + iDense];
    int next = nextColumn[iColumn];
    int numberInPivotColumn = iDense;
    int space = startColumnU[next]
      - startColumnU[iColumn]
      - numberInColumnPlus[next];
    //assume no zero elements
    if (numberInPivotColumn > space) {
      //getColumnSpace also moves fixed part
      if (!getColumnSpace(iColumn, numberInPivotColumn)) {
        return -99;
      }
    }
    // set so further moves will work
    numberInColumn[iColumn] = numberInPivotColumn;
  }
  // Fill in ?
  for (iColumn = numberGoodU_ + numberToDo; iColumn < numberRows_; iColumn++) {
    nextRow[iColumn] = iColumn;
    startColumnL[iColumn + 1] = endL;
    pivotRegion[iColumn] = 1.0;
  }
  if (lengthL_ + full * 0.5 > lengthAreaL_) {
    //need more memory
    if ((messageLevel_ & 4) != 0)
      std::cout << "more memory needed in middle of invert" << std::endl;
    return -99;
  }
  CoinFactorizationDouble *COIN_RESTRICT elementU = elementUArray_;
  for (iDense = 0; iDense < numberToDo; iDense++) {
    int iRow;
    int jDense;
    int pivotRow = -1;
    double *COIN_RESTRICT element = denseAreaAddress_ + iDense * numberDense_;
    CoinFactorizationDouble largest = 1.0e-12;
    for (iRow = iDense; iRow < numberDense_; iRow++) {
      if (fabs(element[iRow]) > largest) {
        largest = fabs(element[iRow]);
        pivotRow = iRow;
      }
    }
    if (pivotRow >= 0) {
      iColumn = pivotColumnConst[base + iDense];
      CoinFactorizationDouble pivotElement = element[pivotRow];
      // get original row
      int originalRow = densePermute_[pivotRow];
      // do nextRow
      nextRow[originalRow] = numberGoodU_;
      lastRow[originalRow] = -2; //mark
      // swap
      densePermute_[pivotRow] = densePermute_[iDense];
      densePermute_[iDense] = originalRow;
      for (jDense = iDense; jDense < numberDense_; jDense++) {
        CoinFactorizationDouble value = element[iDense];
        element[iDense] = element[pivotRow];
        element[pivotRow] = value;
        element += numberDense_;
      }
      CoinFactorizationDouble pivotMultiplier = 1.0 / pivotElement;
      //printf("pivotMultiplier %g\n",pivotMultiplier);
      pivotRegion[numberGoodU_] = pivotMultiplier;
      // Do L
      element = denseAreaAddress_ + iDense * numberDense_;
      int l = lengthL_;
      startColumnL[numberGoodL_] = l; //for luck and first time
      for (iRow = iDense + 1; iRow < numberDense_; iRow++) {
        CoinFactorizationDouble value = element[iRow] * pivotMultiplier;
        element[iRow] = value;
        if (fabs(value) > tolerance) {
          indexRowL[l] = densePermute_[iRow];
          elementL[l++] = value;
        }
      }
      numberGoodL_++;
      lengthL_ = l;
      startColumnL[numberGoodL_] = l;
      // update U column
      int start = startColumnU[iColumn];
      for (iRow = 0; iRow < iDense; iRow++) {
        if (fabs(element[iRow]) > tolerance) {
          indexRowU[start] = densePermute_[iRow];
          elementU[start++] = element[iRow];
        }
      }
      numberInColumn[iColumn] = 0;
      numberInColumnPlus[iColumn] += start - startColumnU[iColumn];
      startColumnU[iColumn] = start;
      // update other columns
      double *COIN_RESTRICT element2 = element + numberDense_;
      for (jDense = iDense + 1; jDense < numberToDo; jDense++) {
        CoinFactorizationDouble value = element2[iDense];
        for (iRow = iDense + 1; iRow < numberDense_; iRow++) {
          //double oldValue=element2[iRow];
          element2[iRow] -= value * element[iRow];
          //if (oldValue&&!element2[iRow]) {
          //printf("Updated element for column %d, row %d old %g",
          //   pivotColumnConst[base+jDense],densePermute_[iRow],oldValue);
          //printf(" new %g\n",element2[iRow]);
          //}
        }
        element2 += numberDense_;
      }
      numberGoodU_++;
    } else {
      return -1;
    }
  }
  // free area (could use L?)
  delete[] denseArea_;
  denseArea_ = NULL;
  // check if can use another array for densePermute_
  delete[] densePermute_;
  densePermute_ = NULL;
  numberDense_ = 0;
  return status;
}
// Separate out links with same row/column count
void CoinFactorization::separateLinks(int count, bool rowsFirst)
{
  int *COIN_RESTRICT nextCount = nextCountArray_;
  int *COIN_RESTRICT firstCount = firstCountArray_;
  int *COIN_RESTRICT lastCount = lastCountArray_;
  int next = firstCount[count];
  int firstRow = -1;
  int firstColumn = -1;
  int lastRow = -1;
  int lastColumn = -1;
  while (next >= 0) {
    int next2 = nextCount[next];
    if (next >= numberRows_) {
      nextCount[next] = -1;
      // Column
      if (firstColumn >= 0) {
        lastCount[next] = lastColumn;
        nextCount[lastColumn] = next;
      } else {
        lastCount[next] = -2 - count;
        firstColumn = next;
      }
      lastColumn = next;
    } else {
      // Row
      if (firstRow >= 0) {
        lastCount[next] = lastRow;
        nextCount[lastRow] = next;
      } else {
        lastCount[next] = -2 - count;
        firstRow = next;
      }
      lastRow = next;
    }
    next = next2;
  }
  if (rowsFirst && firstRow >= 0) {
    firstCount[count] = firstRow;
    nextCount[lastRow] = firstColumn;
    if (firstColumn >= 0)
      lastCount[firstColumn] = lastRow;
  } else if (firstRow < 0) {
    firstCount[count] = firstColumn;
  } else if (firstColumn >= 0) {
    firstCount[count] = firstColumn;
    nextCount[lastColumn] = firstRow;
    if (firstRow >= 0)
      lastCount[firstRow] = lastColumn;
  }
}

// Debug - save on file
#if COINUTILS_BIGINDEX_IS_INT
int CoinFactorization::saveFactorization(const char *file) const
{
  FILE *fp = fopen(file, "wb");
  if (fp) {
    // Save so we can pick up scalars
    const char *first = reinterpret_cast< const char * >(&pivotTolerance_);
    const char *last = reinterpret_cast< const char * >(&biasLU_);
    // increment
    last += sizeof(CoinBigIndex);
    if (fwrite(first, last - first, 1, fp) != 1)
      return 1;
    // Now arrays
    if (CoinToFile(elementU_.array(), lengthAreaU_, fp))
      return 1;
    if (CoinToFile(indexRowU_.array(), lengthAreaU_, fp))
      return 1;
    if (CoinToFile(indexColumnU_.array(), lengthAreaU_, fp))
      return 1;
    if (CoinToFile(convertRowToColumnU_.array(), lengthAreaU_, fp))
      return 1;
    if (CoinToFile(elementByRowL_.array(), lengthAreaL_, fp))
      return 1;
    if (CoinToFile(indexColumnL_.array(), lengthAreaL_, fp))
      return 1;
    if (CoinToFile(startRowL_.array(), numberRows_ + 1, fp))
      return 1;
    if (CoinToFile(elementL_.array(), lengthAreaL_, fp))
      return 1;
    if (CoinToFile(indexRowL_.array(), lengthAreaL_, fp))
      return 1;
    if (CoinToFile(startColumnL_.array(), numberRows_ + 1, fp))
      return 1;
    if (CoinToFile(markRow_.array(), numberRows_, fp))
      return 1;
    if (CoinToFile(saveColumn_.array(), numberColumns_, fp))
      return 1;
    if (CoinToFile(startColumnR_.array(), maximumPivots_ + 1, fp))
      return 1;
    if (CoinToFile(startRowU_.array(), maximumRowsExtra_ + 1, fp))
      return 1;
    if (CoinToFile(numberInRow_.array(), maximumRowsExtra_ + 1, fp))
      return 1;
    if (CoinToFile(nextRow_.array(), maximumRowsExtra_ + 1, fp))
      return 1;
    if (CoinToFile(lastRow_.array(), maximumRowsExtra_ + 1, fp))
      return 1;
    if (CoinToFile(pivotRegion_.array(), maximumRowsExtra_ + 1, fp))
      return 1;
    if (CoinToFile(permuteBack_.array(), maximumRowsExtra_ + 1, fp))
      return 1;
    if (CoinToFile(permute_.array(), maximumRowsExtra_ + 1, fp))
      return 1;
    if (CoinToFile(pivotColumnBack_.array(), maximumRowsExtra_ + 1, fp))
      return 1;
    if (CoinToFile(startColumnU_.array(), maximumColumnsExtra_ + 1, fp))
      return 1;
    if (CoinToFile(numberInColumn_.array(), maximumColumnsExtra_ + 1, fp))
      return 1;
    if (CoinToFile(numberInColumnPlus_.array(), maximumColumnsExtra_ + 1, fp))
      return 1;
    if (CoinToFile(firstCount_.array(), biggerDimension_ + 2, fp))
      return 1;
    if (CoinToFile(nextCount_.array(), numberRows_ + numberColumns_, fp))
      return 1;
    if (CoinToFile(lastCount_.array(), numberRows_ + numberColumns_, fp))
      return 1;
    if (CoinToFile(pivotRowL_.array(), numberRows_ + 1, fp))
      return 1;
    if (CoinToFile(pivotColumn_.array(), maximumColumnsExtra_ + 1, fp))
      return 1;
    if (CoinToFile(nextColumn_.array(), maximumColumnsExtra_ + 1, fp))
      return 1;
    if (CoinToFile(lastColumn_.array(), maximumColumnsExtra_ + 1, fp))
      return 1;
    if (CoinToFile(denseAreaAddress_, numberDense_ * numberDense_, fp))
      return 1;
    if (CoinToFile(densePermute_, numberDense_, fp))
      return 1;
    fclose(fp);
  }
  return 0;
}
// Debug - restore from file
int CoinFactorization::restoreFactorization(const char *file, bool factorIt)
{
  FILE *fp = fopen(file, "rb");
  if (fp) {
    // Get rid of current
    gutsOfDestructor();
    int newSize = 0; // for checking - should be same
    // Restore so we can pick up scalars
    char *first = reinterpret_cast< char * >(&pivotTolerance_);
    char *last = reinterpret_cast< char * >(&biasLU_);
    // increment
    last += sizeof(CoinBigIndex);
    if (fread(first, last - first, 1, fp) != 1)
      return 1;
    int space = lengthAreaL_ - lengthL_;
    // Now arrays
    CoinFactorizationDouble *elementU = elementU_.array();
    if (CoinFromFile(elementU, lengthAreaU_, fp, newSize) == 1)
      return 1;
    assert(newSize == lengthAreaU_);
    int *indexRowU = indexRowU_.array();
    if (CoinFromFile(indexRowU, lengthAreaU_, fp, newSize) == 1)
      return 1;
    assert(newSize == lengthAreaU_);
    int *indexColumnU = indexColumnU_.array();
    if (CoinFromFile(indexColumnU, lengthAreaU_, fp, newSize) == 1)
      return 1;
    assert(newSize == lengthAreaU_);
    int *convertRowToColumnU = convertRowToColumnU_.array();
    if (CoinFromFile(convertRowToColumnU, lengthAreaU_, fp, newSize) == 1)
      return 1;
    assert(newSize == lengthAreaU_ || (newSize == 0 && !convertRowToColumnU_.array()));
    CoinFactorizationDouble *elementByRowL = elementByRowL_.array();
    if (CoinFromFile(elementByRowL, lengthAreaL_, fp, newSize) == 1)
      return 1;
    assert(newSize == lengthAreaL_ || (newSize == 0 && !elementByRowL_.array()));
    int *indexColumnL = indexColumnL_.array();
    if (CoinFromFile(indexColumnL, lengthAreaL_, fp, newSize) == 1)
      return 1;
    assert(newSize == lengthAreaL_ || (newSize == 0 && !indexColumnL_.array()));
    int *startRowL = startRowL_.array();
    if (CoinFromFile(startRowL, numberRows_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == numberRows_ + 1 || (newSize == 0 && !startRowL_.array()));
    CoinFactorizationDouble *elementL = elementL_.array();
    if (CoinFromFile(elementL, lengthAreaL_, fp, newSize) == 1)
      return 1;
    assert(newSize == lengthAreaL_);
    int *indexRowL = indexRowL_.array();
    if (CoinFromFile(indexRowL, lengthAreaL_, fp, newSize) == 1)
      return 1;
    assert(newSize == lengthAreaL_);
    int *startColumnL = startColumnL_.array();
    if (CoinFromFile(startColumnL, numberRows_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == numberRows_ + 1);
    int *markRow = markRow_.array();
    if (CoinFromFile(markRow, numberRows_, fp, newSize) == 1)
      return 1;
    assert(newSize == numberRows_);
    int *saveColumn = saveColumn_.array();
    if (CoinFromFile(saveColumn, numberColumns_, fp, newSize) == 1)
      return 1;
    assert(newSize == numberColumns_);
    int *startColumnR = startColumnR_.array();
    if (CoinFromFile(startColumnR, maximumPivots_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == maximumPivots_ + 1 || (newSize == 0 && !startColumnR_.array()));
    int *startRowU = startRowU_.array();
    if (CoinFromFile(startRowU, maximumRowsExtra_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == maximumRowsExtra_ + 1 || (newSize == 0 && !startRowU_.array()));
    int *numberInRow = numberInRow_.array();
    if (CoinFromFile(numberInRow, maximumRowsExtra_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == maximumRowsExtra_ + 1);
    int *nextRow = nextRow_.array();
    if (CoinFromFile(nextRow, maximumRowsExtra_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == maximumRowsExtra_ + 1);
    int *lastRow = lastRow_.array();
    if (CoinFromFile(lastRow, maximumRowsExtra_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == maximumRowsExtra_ + 1);
    CoinFactorizationDouble *pivotRegion = pivotRegion_.array();
    if (CoinFromFile(pivotRegion, maximumRowsExtra_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == maximumRowsExtra_ + 1);
    int *permuteBack = permuteBack_.array();
    if (CoinFromFile(permuteBack, maximumRowsExtra_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == maximumRowsExtra_ + 1 || (newSize == 0 && !permuteBack_.array()));
    int *permute = permute_.array();
    if (CoinFromFile(permute, maximumRowsExtra_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == maximumRowsExtra_ + 1 || (newSize == 0 && !permute_.array()));
    int *pivotColumnBack = pivotColumnBack_.array();
    if (CoinFromFile(pivotColumnBack, maximumRowsExtra_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == maximumRowsExtra_ + 1 || (newSize == 0 && !pivotColumnBack_.array()));
    int *startColumnU = startColumnU_.array();
    if (CoinFromFile(startColumnU, maximumColumnsExtra_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == maximumColumnsExtra_ + 1);
    int *numberInColumn = numberInColumn_.array();
    if (CoinFromFile(numberInColumn, maximumColumnsExtra_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == maximumColumnsExtra_ + 1);
    int *numberInColumnPlus = numberInColumnPlus_.array();
    if (CoinFromFile(numberInColumnPlus, maximumColumnsExtra_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == maximumColumnsExtra_ + 1);
    int *firstCount = firstCount_.array();
    if (CoinFromFile(firstCount, biggerDimension_ + 2, fp, newSize) == 1)
      return 1;
    assert(newSize == biggerDimension_ + 2);
    int *nextCount = nextCount_.array();
    if (CoinFromFile(nextCount, numberRows_ + numberColumns_, fp, newSize) == 1)
      return 1;
    assert(newSize == numberRows_ + numberColumns_);
    int *lastCount = lastCount_.array();
    if (CoinFromFile(lastCount, numberRows_ + numberColumns_, fp, newSize) == 1)
      return 1;
    assert(newSize == numberRows_ + numberColumns_);
    int *pivotRowL = pivotRowL_.array();
    if (CoinFromFile(pivotRowL, numberRows_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == numberRows_ + 1);
    int *pivotColumn = pivotColumn_.array();
    if (CoinFromFile(pivotColumn, maximumColumnsExtra_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == maximumColumnsExtra_ + 1);
    int *nextColumn = nextColumn_.array();
    if (CoinFromFile(nextColumn, maximumColumnsExtra_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == maximumColumnsExtra_ + 1);
    int *lastColumn = lastColumn_.array();
    if (CoinFromFile(lastColumn, maximumColumnsExtra_ + 1, fp, newSize) == 1)
      return 1;
    assert(newSize == maximumColumnsExtra_ + 1);
    if (CoinFromFile(denseAreaAddress_, numberDense_ * numberDense_, fp, newSize) == 1)
      return 1;
    assert(newSize == numberDense_ * numberDense_);
    if (CoinFromFile(densePermute_, numberDense_, fp, newSize) == 1)
      return 1;
    assert(newSize == numberDense_);
    lengthAreaR_ = space;
    elementR_ = elementL_.array() + lengthL_;
    indexRowR_ = indexRowL_.array() + lengthL_;
    fclose(fp);
    if (factorIt) {
      if (biasLU_ >= 3 || numberRows_ != numberColumns_)
        preProcess(2);
      else
        preProcess(3); // no row copy
      factor();
    }
  }
  return 0;
}
#endif
//  factorSparse.  Does sparse phase of factorization
//return code is <0 error, 0= finished
int CoinFactorization::factorSparseLarge()
{
  int *COIN_RESTRICT indexRow = indexRowUArray_;
  int *COIN_RESTRICT indexColumn = indexColumnUArray_;
  CoinFactorizationDouble *COIN_RESTRICT element = elementUArray_;
  int count = 1;
  workArea_.conditionalNew(numberRows_);
  CoinFactorizationDouble *COIN_RESTRICT workArea = workArea_.array();
  workAreaArray_ = workArea;
#ifndef NDEBUG
  counter1++;
#endif
  // when to go dense
  int denseThreshold = abs(denseThreshold_);

  CoinZeroN(workArea, numberRows_);
  //get space for bit work area
  int workSize = 1000;
  workArea2_.conditionalNew(workSize);
  unsigned int *COIN_RESTRICT workArea2 = workArea2_.array();
  workArea2Array_ = workArea2;

  //set markRow so no rows updated
  int *COIN_RESTRICT markRow = markRowArray_;
  CoinFillN(markRow, numberRows_, COIN_INT_MAX - 10 + 1);
  int status = 0;
  //do slacks first
  int *COIN_RESTRICT numberInRow = numberInRowArray_;
  int *COIN_RESTRICT numberInColumn = numberInColumnArray_;
  int *COIN_RESTRICT numberInColumnPlus = numberInColumnPlusArray_;
  int *COIN_RESTRICT startColumnU = startColumnUArray_;
  int *COIN_RESTRICT startColumnL = startColumnLArray_;
  if (biasLU_ < 3 && numberColumns_ == numberRows_) {
    int iPivotColumn;
    int *COIN_RESTRICT pivotColumn = pivotColumnArray_;
    int *COIN_RESTRICT nextRow = nextRowArray_;
    int *COIN_RESTRICT lastRow = lastRowArray_;
    for (iPivotColumn = 0; iPivotColumn < numberColumns_;
         iPivotColumn++) {
      if (numberInColumn[iPivotColumn] == 1) {
        int start = startColumnU[iPivotColumn];
        CoinFactorizationDouble value = element[start];
        if (value == slackValue_ && numberInColumnPlus[iPivotColumn] == 0) {
          // treat as slack
          int iRow = indexRow[start];
          // but only if row not marked
          if (numberInRow[iRow] > 0) {
            totalElements_ -= numberInRow[iRow];
            //take out this bit of indexColumnU
            int next = nextRow[iRow];
            int last = lastRow[iRow];

            nextRow[last] = next;
            lastRow[next] = last;
            nextRow[iRow] = numberGoodU_; //use for permute
            lastRow[iRow] = -2; //mark
            //modify linked list for pivots
            deleteLink(iRow);
            numberInRow[iRow] = -1;
            numberInColumn[iPivotColumn] = 0;
            numberGoodL_++;
            startColumnL[numberGoodL_] = 0;
            pivotColumn[numberGoodU_] = iPivotColumn;
            numberGoodU_++;
          }
        }
      }
    }
    // redo
    preProcess(4);
    CoinFillN(markRow, numberRows_, COIN_INT_MAX - 10 + 1);
  }
  numberSlacks_ = numberGoodU_;
  int *COIN_RESTRICT nextCount = nextCountArray_;
  int *COIN_RESTRICT firstCount = firstCountArray_;
  int *COIN_RESTRICT startRow = startRowUArray_;
  int *COIN_RESTRICT startColumn = startColumnU;
  //double *elementL = elementLArray_;
  //int *indexRowL = indexRowLArray_;
  //int *saveColumn = saveColumnArray_;
  //int *nextRow = nextRowArray_;
  //int *lastRow = lastRowArray_ ;
  double pivotTolerance = pivotTolerance_;
  int numberTrials = numberTrials_;
  int numberRows = numberRows_;
  // Put column singletons first - (if false)
  separateLinks(1, (biasLU_ > 1));
#ifndef NDEBUG
  int counter2 = 0;
#endif
  while (count <= biggerDimension_) {
#ifndef NDEBUG
    counter2++;
    int badRow = -1;
    if (counter1 == -1 && counter2 >= 0) {
      // check counts consistent
      for (int iCount = 1; iCount < numberRows_; iCount++) {
        int look = firstCount[iCount];
        while (look >= 0) {
          if (look < numberRows_) {
            int iRow = look;
            if (iRow == badRow)
              printf("row count for row %d is %d\n", iCount, iRow);
            if (numberInRow[iRow] != iCount) {
              printf("failed debug on %d entry to factorSparse and %d try\n",
                (int)counter1, counter2);
              printf("row %d - count %d number %d\n", iRow, iCount, numberInRow[iRow]);
              abort();
            }
            look = nextCount[look];
          } else {
            int iColumn = look - numberRows;
            if (numberInColumn[iColumn] != iCount) {
              printf("failed debug on %d entry to factorSparse and %d try\n",
                (int)counter1, counter2);
              printf("column %d - count %d number %d\n", iColumn, iCount, numberInColumn[iColumn]);
              abort();
            }
            look = nextCount[look];
          }
        }
      }
    }
#endif
    int minimumCount = COIN_INT_MAX;
    double minimumCost = COIN_DBL_MAX;

    int iPivotRow = -1;
    int iPivotColumn = -1;
    int pivotRowPosition = -1;
    int pivotColumnPosition = -1;
    int look = firstCount[count];
    int trials = 0;
    int *COIN_RESTRICT pivotColumn = pivotColumnArray_;

    if (count == 1 && firstCount[1] >= 0 && !biasLU_) {
      //do column singletons first to put more in U
      while (look >= 0) {
        if (look < numberRows_) {
          look = nextCount[look];
        } else {
          int iColumn = look - numberRows_;

          assert(numberInColumn[iColumn] == count);
          int start = startColumnU[iColumn];
          int iRow = indexRow[start];

          iPivotRow = iRow;
          pivotRowPosition = start;
          iPivotColumn = iColumn;
          assert(iPivotRow >= 0 && iPivotColumn >= 0);
          pivotColumnPosition = -1;
          look = -1;
          break;
        }
      } /* endwhile */
      if (iPivotRow < 0) {
        //back to singletons
        look = firstCount[1];
      }
    }
    while (look >= 0) {
      if (look < numberRows_) {
        int iRow = look;
#ifndef NDEBUG
        if (numberInRow[iRow] != count) {
          printf("failed on %d entry to factorSparse and %d try\n",
            (int)counter1, counter2);
          printf("row %d - count %d number %d\n", iRow, count, numberInRow[iRow]);
          abort();
        }
#endif
        look = nextCount[look];
        bool rejected = false;
        int start = startRow[iRow];
        int end = start + count;

        int i;
        for (i = start; i < end; i++) {
          int iColumn = indexColumn[i];
          assert(numberInColumn[iColumn] > 0);
          double cost = (count - 1) * numberInColumn[iColumn];

          if (cost < minimumCost) {
            int where = startColumn[iColumn];
            CoinFactorizationDouble minimumValue = element[where];

            minimumValue = fabs(minimumValue) * pivotTolerance;
            while (indexRow[where] != iRow) {
              where++;
            } /* endwhile */
            assert(where < startColumn[iColumn] + numberInColumn[iColumn]);
            CoinFactorizationDouble value = element[where];

            value = fabs(value);
            if (value >= minimumValue) {
              minimumCost = cost;
              minimumCount = numberInColumn[iColumn];
              iPivotRow = iRow;
              pivotRowPosition = -1;
              iPivotColumn = iColumn;
              assert(iPivotRow >= 0 && iPivotColumn >= 0);
              pivotColumnPosition = i;
              rejected = false;
              if (minimumCount < count) {
                look = -1;
                break;
              }
            } else if (iPivotRow == -1) {
              rejected = true;
            }
          }
        }
        trials++;
        if (trials >= numberTrials && iPivotRow >= 0) {
          look = -1;
          break;
        }
        if (rejected) {
          //take out for moment
          //eligible when row changes
          deleteLink(iRow);
          addLink(iRow, biggerDimension_ + 1);
        }
      } else {
        int iColumn = look - numberRows;

        assert(numberInColumn[iColumn] == count);
        look = nextCount[look];
        int start = startColumn[iColumn];
        int end = start + numberInColumn[iColumn];
        CoinFactorizationDouble minimumValue = element[start];

        minimumValue = fabs(minimumValue) * pivotTolerance;
        int i;
        for (i = start; i < end; i++) {
          CoinFactorizationDouble value = element[i];

          value = fabs(value);
          if (value >= minimumValue) {
            int iRow = indexRow[i];
            int nInRow = numberInRow[iRow];
            assert(nInRow > 0);
            double cost = (count - 1) * nInRow;

            if (cost < minimumCost) {
              minimumCost = cost;
              minimumCount = nInRow;
              iPivotRow = iRow;
              pivotRowPosition = i;
              iPivotColumn = iColumn;
              assert(iPivotRow >= 0 && iPivotColumn >= 0);
              pivotColumnPosition = -1;
              if (minimumCount <= count + 1) {
                look = -1;
                break;
              }
            }
          }
        }
        trials++;
        if (trials >= numberTrials && iPivotRow >= 0) {
          look = -1;
          break;
        }
      }
    } /* endwhile */
    if (iPivotRow >= 0) {
      if (iPivotRow >= 0) {
        int numberDoRow = numberInRow[iPivotRow] - 1;
        int numberDoColumn = numberInColumn[iPivotColumn] - 1;

        totalElements_ -= (numberDoRow + numberDoColumn + 1);
        if (numberDoColumn > 0) {
          if (numberDoRow > 0) {
            if (numberDoColumn > 1) {
              //  if (1) {
              //need to adjust more for cache and SMP
              //allow at least 4 extra
              int increment = numberDoColumn + 1 + 4;

              if (increment & 15) {
                increment = increment & (~15);
                increment += 16;
              }
              int increment2 =

                (increment + COINFACTORIZATION_BITS_PER_INT - 1) >> COINFACTORIZATION_SHIFT_PER_INT;
              int size = increment2 * numberDoRow;

              if (size > workSize) {
                workSize = size;
                workArea2_.conditionalNew(workSize);
                workArea2 = workArea2_.array();
		workArea2Array_ = workArea2;
              }
              bool goodPivot;

              //might be able to do better by permuting
#ifndef UGLY_COIN_FACTOR_CODING
              //branch out to best pivot routine
              goodPivot = pivot(iPivotRow, iPivotColumn,
                pivotRowPosition, pivotColumnPosition,
                workArea, workArea2,
                increment2, markRow,
                LARGE_SET);
#else
#define FAC_SET LARGE_SET
#include "CoinFactorization.hpp"
#undef FAC_SET
#endif
              if (!goodPivot) {
                status = -99;
                count = biggerDimension_ + 1;
                break;
              }
            } else {
              if (!pivotOneOtherRow(iPivotRow, iPivotColumn)) {
                status = -99;
                count = biggerDimension_ + 1;
                break;
              }
            }
          } else {
            assert(!numberDoRow);
            if (!pivotRowSingleton(iPivotRow, iPivotColumn)) {
              status = -99;
              count = biggerDimension_ + 1;
              break;
            }
          }
        } else {
          assert(!numberDoColumn);
          if (!pivotColumnSingleton(iPivotRow, iPivotColumn)) {
            status = -99;
            count = biggerDimension_ + 1;
            break;
          }
        }
        assert(nextRow_.array()[iPivotRow] == numberGoodU_);
        pivotColumn[numberGoodU_] = iPivotColumn;
        numberGoodU_++;
        // This should not need to be trapped here - but be safe
        if (numberGoodU_ == numberRows_)
          count = biggerDimension_ + 1;
      }
#if COIN_DEBUG == 2
      checkConsistency();
#endif
      if (denseThreshold) {
        // see whether to go dense
        int leftRows = numberRows_ - numberGoodU_;
        double full = leftRows;
        full *= full;
        assert(full >= 0.0);
        double leftElements = totalElements_;
        //if (leftRows==100)
        //printf("at 100 %d elements\n",totalElements_);
        double ratio;
        if (leftRows > 2000) {
          ratio = 4.0;
#if COIN_DENSE_MULTIPLIER2 == 1
          ratio = 3.5;
#endif
        } else if (leftRows > 800) {
          ratio = 3.0;
#if COIN_DENSE_MULTIPLIER2 == 1
          ratio = 2.75;
#endif
        } else if (leftRows > 256) {
          ratio = 2.0;
        } else {
          ratio = 1.5;
        }
#if COIN_DENSE_MULTIPLIER2 > 10
        ratio = 10000;
#endif
        ratio *= COIN_DENSE_MULTIPLIER;
        if ((ratio * leftElements > full && leftRows > denseThreshold)) {
#if COIN_ALIGN_DENSE == 2
          if ((leftRows & 7) == 0) {
#endif
            //return to do dense
            if (status != 0)
              break;
            int check = 0;
            for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
              if (numberInColumn[iColumn])
                check++;
            }
            if (check != leftRows && denseThreshold_) {
              //printf("** mismatch %d columns left, %d rows\n",check,leftRows);
              denseThreshold = 0;
            } else {
              status = 2;
              if ((messageLevel_ & 4) != 0 && true)
                std::cout << "      Went dense at " << leftRows << " rows " << totalElements_ << " " << full << " " << leftElements << std::endl;
              if (!denseThreshold_)
                denseThreshold_ = -check; // say how many
              break;
            }
#if COIN_ALIGN_DENSE == 2
          }
#endif
        }
      }
      // start at 1 again
      count = 1;
    } else {
      //end of this - onto next
      count++;
    }
  } /* endwhile */
  workArea_.conditionalDelete();
  workArea2_.conditionalDelete();
  return status;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
