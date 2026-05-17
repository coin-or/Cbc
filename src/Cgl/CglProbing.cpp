// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <climits>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <cassert>
#include <iostream>
#define PROBING100 0
// #define PRINT_DEBUG
// #define CGL_DEBUG 1
// #undef NDEBUG
#ifndef FIXED_BOTH_WAYS
#define FIXED_BOTH_WAYS 1
#endif
#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinFinite.hpp"
#include "CoinMpsIO.hpp"
#include "CoinTime.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiRowCut.hpp"
#include "CglProbing.hpp"
// #define PROBING_EXTRA_STUFF true
#define PROBING_EXTRA_STUFF false
#define FIXED_ALLOWANCE 10
#define SIZE_ROW_MULT 4
#define SIZE_ROW_ADD 2000
#define MANY_TIGHTEN 1
typedef struct {
  double infeasibility;
  int sequence;
} double_int_pair;
class double_int_pair_compare {
public:
  bool operator()(double_int_pair x, double_int_pair y) const
  {
    return (x.infeasibility < y.infeasibility);
  }
};
// for hashing
typedef struct {
  int index, next;
} CoinHashLink;

class row_cut {
public:
  row_cut(int nRows, bool initialPass)
  {
    numberCuts_ = 0;
    size_ = 100;
    if (initialPass)
      size_ = size_ >> 1;
    if (size_ < 1000)
      hashSize_ = 4 * size_;
    else
      hashSize_ = 2 * size_;
    nRows_ = nRows;
    rowCut_ = new OsiRowCut2 *[size_];
    hash_ = new CoinHashLink[hashSize_];
    for (int i = 0; i < hashSize_; i++) {
      hash_[i].index = -1;
      hash_[i].next = -1;
    }
    numberCuts_ = 0;
    lastHash_ = -1;
  }
  ~row_cut()
  {
    for (int i = 0; i < numberCuts_; i++)
      delete rowCut_[i];
    delete[] rowCut_;
    delete[] hash_;
  }
  OsiRowCut2 *cut(int i) const
  {
    return rowCut_[i];
  }
  int numberCuts() const
  {
    return numberCuts_;
  }
  OsiRowCut2 **rowCut_;
  /// Hash table
  CoinHashLink *hash_;
  int size_;
  int hashSize_;
  int nRows_;
  int numberCuts_;
  int lastHash_;
  // Return 0 if added, 1 if not, -1 if not added because of space
  int addCutIfNotDuplicate(OsiRowCut &cut, int whichRow = -1)
  {
    if (numberCuts_ == size_) {
      size_ = 2 * size_ + 100;
      if (size_ < 1000)
        hashSize_ = 4 * size_;
      else
        hashSize_ = 2 * size_;
      OsiRowCut2 **temp = new OsiRowCut2 *[size_];
      delete[] hash_;
      hash_ = new CoinHashLink[hashSize_];
      for (int i = 0; i < hashSize_; i++) {
        hash_[i].index = -1;
        hash_[i].next = -1;
      }
      lastHash_ = -1;
      for (int i = 0; i < numberCuts_; i++) {
        temp[i] = rowCut_[i];
        int ipos = hashCut(*temp[i], hashSize_);
        int jpos = ipos;
        int found = -1;
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
            hash_[ipos].index = i;
          } else {
            while (true) {
              ++lastHash_;
              assert(lastHash_ < hashSize_);
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
      // Sort the input cut's row in place (avoid copy)
      cut.mutableRow().sortIncrIndex();
      int numberElements = cut.row().getNumElements();
      const int *newIndices = cut.row().getIndices();
      const double *newElements = cut.row().getElements();
      bool bad = false;
      for (int i = 0; i < numberElements; i++) {
        double value = fabs(newElements[i]);
        if (value < 1.0e-12 || value > 1.0e12)
          bad = true;
      }
      if (newLb < -1.0e12 && newLb > -1.0e40)
        bad = true;
      if (newUb > 1.0e12 && newUb < 1.0e40)
        bad = true;
      if (bad)
        return 1;
      // Hash directly from the cut (no copy needed)
      int ipos = hashCut(cut, hashSize_);
      int found = -1;
      int jpos = ipos;
      while (true) {
        int j1 = hash_[ipos].index;
        if (j1 >= 0) {
          if (!same(cut, *rowCut_[j1])) {
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
          hash_[ipos].index = numberCuts_;
        } else {
          while (true) {
            ++lastHash_;
            assert(lastHash_ < hashSize_);
            if (hash_[lastHash_].index == -1)
              break;
          }
          hash_[ipos].next = lastHash_;
          hash_[lastHash_].index = numberCuts_;
        }
        // Only allocate and copy when actually storing
        OsiRowCut2 *newCutPtr = new OsiRowCut2(whichRow);
        newCutPtr->setLb(newLb);
        newCutPtr->setUb(newUb);
        newCutPtr->setRow(std::move(cut.mutableRow()));
        rowCut_[numberCuts_++] = newCutPtr;
        return 0;
      } else {
        return 1;
      }
    } else {
      return -1;
    }
  }
  void addCuts(OsiCuts &cs, OsiRowCut **whichRow, int iPass)
  {
    int numberCuts = cs.sizeRowCuts();
    int i;
    if (numberCuts_ < nRows_) {
      if ((iPass & 1) == 1) {
        for (i = 0; i < numberCuts_; i++) {
          cs.insert(*rowCut_[i]);
          if (whichRow) {
            int iRow = rowCut_[i]->whichRow();
            if (iRow >= 0 && !whichRow[iRow])
              whichRow[iRow] = cs.rowCutPtr(numberCuts);
            ;
          }
          numberCuts++;
        }
      } else {
        for (i = numberCuts_ - 1; i >= 0; i--) {
          cs.insert(*rowCut_[i]);
          if (whichRow) {
            int iRow = rowCut_[i]->whichRow();
            if (iRow >= 0 && !whichRow[iRow])
              whichRow[iRow] = cs.rowCutPtr(numberCuts);
            ;
          }
          numberCuts++;
        }
      }
    } else {
      // just best
      double *effectiveness = new double[numberCuts_];
      int iCut = 0;
      for (i = 0; i < numberCuts_; i++) {
        double value = -rowCut_[i]->effectiveness();
        if (whichRow) {
          int iRow = rowCut_[i]->whichRow();
          if (iRow >= 0)
            value -= 1.0e10;
        }
        effectiveness[iCut++] = value;
      }
      std::sort(effectiveness, effectiveness + numberCuts_);
      double threshold = -1.0e20;
      if (iCut > nRows_)
        threshold = effectiveness[nRows_];
      for (i = 0; i < numberCuts_; i++) {
        if (rowCut_[i]->effectiveness() > threshold) {
          cs.insert(*rowCut_[i]);
          if (whichRow) {
            int iRow = rowCut_[i]->whichRow();
            if (iRow >= 0 && !whichRow[iRow])
              whichRow[iRow] = cs.rowCutPtr(numberCuts);
            ;
          }
          numberCuts++;
        }
      }
      delete[] effectiveness;
    }
    for (i = 0; i < numberCuts_; i++) {
      delete rowCut_[i];
      rowCut_[i] = 0;
    }
    numberCuts_ = 0;
  }
};
// Adds in cut to list
#ifdef CGL_DEBUG
// Checks bounds okay against debugger
static void checkBounds(const OsiRowCutDebugger *debugger, OsiColCut &cut)
{
  if (debugger) {
    // on optimal path
    const double *optimal = debugger->optimalSolution();
    int i;
    int nIndex;
    const double *values;
    const int *index;
    const CoinPackedVector &lbs = cut.lbs();
    values = lbs.getElements();
    nIndex = lbs.getNumElements();
    index = lbs.getIndices();
    for (i = 0; i < nIndex; i++) {
      double value = values[i];
      int iColumn = index[i];
      printf("%d optimal %g lower %g\n", iColumn, optimal[iColumn], value);
      assert(value <= optimal[iColumn] + 1.0e-5);
    }
    const CoinPackedVector &ubs = cut.ubs();
    values = ubs.getElements();
    nIndex = ubs.getNumElements();
    index = ubs.getIndices();
    for (i = 0; i < nIndex; i++) {
      double value = values[i];
      int iColumn = index[i];
      printf("%d optimal %g upper %g\n", iColumn, optimal[iColumn], value);
      assert(value >= optimal[iColumn] - 1.0e-5);
    }
  }
}
#endif
#define CGL_REASONABLE_INTEGER_BOUND 1.23456789e10
// This tightens column bounds (and can declare infeasibility)
// It may also declare rows to be redundant
int CglProbing::tighten(double *colLower, double *colUpper,
  const int *column, const double *rowElements,
  const CoinBigIndex *rowStart,
  const CoinBigIndex *rowStartPos, const int *rowLength,
  const CoinPackedMatrix *columnCopy,
  double *rowLower, double *rowUpper,
  int nRows, int nCols, const char *intVar, int maxpass,
  double tolerance)
{
  int i, j;
  CoinBigIndex k, krs, kre;
  int dolrows;
  int iflagu, iflagl;
  int /*ntotal=0,*/ nchange = 1, jpass = 0;
  double dmaxup, dmaxdown, dbound;
  int ninfeas = 0;
  // relax tolerance
  tolerance *= 10.0;
  // For clique stuff
  double *cliqueMin = NULL;
  double *cliqueMax = NULL;
  // And second best ones
  double *cliqueMin2 = NULL;
  double *cliqueMax2 = NULL;
  if (cliqueRowStart_ && numberRows_ && cliqueRowStart_[numberRows_]) {
    cliqueMin = new double[nCols];
    cliqueMax = new double[nCols];
    cliqueMin2 = new double[nCols];
    cliqueMax2 = new double[nCols];
  } else {
    // do without cliques and using sorted version
    assert(rowStartPos);
    int *rowsToLookAtA[2];
    int numberLookNext = 0;
    if (!lookColumn_) {
      lookColumn_ = new int[2 * nRows + nCols];
      lookRow_[0] = lookColumn_ + nCols;
      lookRow_[1] = lookRow_[0] + nRows;
      rowsToLookAtA[0] = lookColumn_ + nCols;
      rowsToLookAtA[1] = rowsToLookAtA[0] + nRows;
      markRow_ = new char[nRows + nCols];
      markColumn_ = markRow_ + nRows;
    } else {
      rowsToLookAtA[0] = lookColumn_ + nCols;
      rowsToLookAtA[1] = rowsToLookAtA[0] + nRows;
      numberLookNext = numberRows_;
      numberRows_ = 0;
    }
    char *markRow = markRow_;
    char *markCol = markColumn_;
    int *whichColumn = lookColumn_;
    memset(markRow, 0, nRows + nCols);
    double *minR;
    double *maxR;
    if (minR_) {
      minR = minR_;
      maxR = maxR_;
    } else {
      minR = new double[2 * nRows];
      maxR = minR + nRows;
    }
    int *rowsToLookAt0 = rowsToLookAtA[0];
    int iWhichLook = 0;
    const int *row = columnCopy->getIndices();
    const CoinBigIndex *columnStart = columnCopy->getVectorStarts();
    const int *columnLength = columnCopy->getVectorLengths();
    const double *columnElements = columnCopy->getElements();
    // check if already filled in
    if (!numberLookNext) {
      for (int i = 0; i < nRows; ++i) {
        if (rowLower[i] > -1.0e20 || rowUpper[i] < 1.0e20) {
          int iflagu = 0;
          int iflagl = 0;
          double dmaxup = 0.0;
          double dmaxdown = 0.0;
          CoinBigIndex krs = rowStart[i];
          CoinBigIndex krs2 = rowStartPos[i];
          CoinBigIndex kre = rowStart[i] + rowLength[i];

          /* ------------------------------------------------------------*/
          /* Compute L(i) and U(i) */
          /* ------------------------------------------------------------*/
          for (k = krs; k < krs2; ++k) {
            double value = rowElements[k];
            int j = column[k];
            if (colUpper[j] < 1.0e12)
              dmaxdown += colUpper[j] * value;
            else
              ++iflagl;
            if (colLower[j] > -1.0e12)
              dmaxup += colLower[j] * value;
            else
              ++iflagu;
          }
          for (k = krs2; k < kre; ++k) {
            double value = rowElements[k];
            int j = column[k];
            if (colUpper[j] < 1.0e12)
              dmaxup += colUpper[j] * value;
            else
              ++iflagu;
            if (colLower[j] > -1.0e12)
              dmaxdown += colLower[j] * value;
            else
              ++iflagl;
          }
          if (iflagu)
            dmaxup = 1.0e31;
          if (iflagl)
            dmaxdown = -1.0e31;
          rowsToLookAt0[numberLookNext++] = i;
          minR[i] = dmaxdown;
          maxR[i] = dmaxup;
        } else {
          // ignore row
        }
      }
    }
    while (nchange) {
      nchange = 0;
      if (jpass == maxpass)
        break;
      jpass++;
      dolrows = (jpass & 1) == 1;
      int *rowsToLookAt = rowsToLookAtA[iWhichLook];
      // for next time
      iWhichLook = 1 - iWhichLook;
      int *rowsToLookAtNext = rowsToLookAtA[iWhichLook];
      int numberLook = numberLookNext;
      numberLookNext = 0;
      int numberBoundsChanged = 0;
      for (int iLook = 0; iLook < numberLook; iLook++) {
        int i = rowsToLookAt[iLook];
        if (rowLower[i] < -1.0e20 && rowUpper[i] > 1.0e20)
          continue; // crept in
        int iflagu = 0;
        int iflagl = 0;
        double dmaxup = 0.0;
        double dmaxdown = 0.0;
        CoinBigIndex krs = rowStart[i];
        CoinBigIndex krs2 = rowStartPos[i];
        CoinBigIndex kre = rowStart[i] + rowLength[i];

        /* ------------------------------------------------------------*/
        /* Compute L(i) and U(i) */
        /* ------------------------------------------------------------*/
        for (k = krs; k < krs2; ++k) {
          double value = rowElements[k];
          int j = column[k];
          if (colUpper[j] < 1.0e12)
            dmaxdown += colUpper[j] * value;
          else
            ++iflagl;
          if (colLower[j] > -1.0e12)
            dmaxup += colLower[j] * value;
          else
            ++iflagu;
        }
        for (k = krs2; k < kre; ++k) {
          double value = rowElements[k];
          int j = column[k];
          if (colUpper[j] < 1.0e12)
            dmaxup += colUpper[j] * value;
          else
            ++iflagu;
          if (colLower[j] > -1.0e12)
            dmaxdown += colLower[j] * value;
          else
            ++iflagl;
        }
        if (iflagu)
          dmaxup = 1.0e31;
        if (iflagl)
          dmaxdown = -1.0e31;
        if (dmaxup <= rowUpper[i] + tolerance && dmaxdown >= rowLower[i] - tolerance) {
          /*
           * The sum of the column maxs is at most the row ub, and
           * the sum of the column mins is at least the row lb;
           * this row says nothing at all.
           * I suspect that this corresponds to
           * an implied column singleton in the paper (viii, on p. 325),
           * where the singleton in question is the row slack.
           */
          ++nchange;
          rowLower[i] = -COIN_DBL_MAX;
          rowUpper[i] = COIN_DBL_MAX;
          markRow[i] = 2; // say don't use
          continue;
        } else {
          if (dmaxup < rowLower[i] - tolerance || dmaxdown > rowUpper[i] + tolerance) {
            ninfeas++;
            break;
          }
          /*        Finite U(i) */
          /* -------------------------------------------------------------*/
          /* below is deliberate mistake (previously was by chance) */
          /*        never do both */
          if (iflagu == 0 && rowLower[i] > 0.0 && iflagl == 0 && rowUpper[i] < 1e15) {
            if (dolrows) {
              iflagu = 1;
            } else {
              iflagl = 1;
            }
          }
          if (iflagu == 0 && rowLower[i] > -1e15) {
            for (k = krs; k < kre; ++k) {
              double value = rowElements[k];
              j = column[k];
              if (value > 0.0) {
                if (colUpper[j] < 1.0e12) {
                  dbound = colUpper[j] + (rowLower[i] - dmaxup) / value;
                  if (dbound > colLower[j] + 1.0e-8) {
                    /* we can tighten the lower bound */
                    /* the paper mentions this as a possibility on p. 227 */
                    colLower[j] = dbound;
                    ++nchange;
                    if (!markCol[j]) {
                      markCol[j] = 1;
                      whichColumn[numberBoundsChanged++] = j;
                    }
                    // see what rows need redoing
                    for (CoinBigIndex k = columnStart[j];
                      k < columnStart[j] + columnLength[j]; k++) {
                      int iRow = row[k];
                      if (markRow[iRow] == 0) {
                        markRow[iRow] = 1;
                        rowsToLookAtNext[numberLookNext++] = iRow;
                      }
                    }
                    /* this may have fixed the variable */
                    /* I believe that this roughly corresponds to a
                     * forcing constraint in the paper (p. 226).
                     * If there is a forcing constraint (with respect
                     * to the original, untightened bounds), then in this
                     * loop we will go through all the columns and fix
                     * each of them to their implied bound, rather than
                     * determining that the row as a whole is forced
                     * and just fixing them without doing computation for
                     * each column (as in the paper).
                     * By doing it this way, we can tighten bounds and
                     * get forcing constraints for free.
                     */
                    /* --------------------------------------------------*/
                    /*                check if infeasible !!!!! */
                    /* --------------------------------------------------*/
                    if (colUpper[j] - colLower[j] < -100.0 * tolerance) {
                      ninfeas++;
                    }
                  }
                }
              } else {
                if (colLower[j] > -1.0e12) {
                  dbound = colLower[j] + (rowLower[i] - dmaxup) / value;
                  if (dbound < colUpper[j] - 1.0e-8) {
                    colUpper[j] = dbound;
                    ++nchange;
                    if (!markCol[j]) {
                      markCol[j] = 1;
                      whichColumn[numberBoundsChanged++] = j;
                    }
                    // see what rows need redoing
                    for (CoinBigIndex k = columnStart[j];
                      k < columnStart[j] + columnLength[j]; k++) {
                      int iRow = row[k];
                      if (markRow[iRow] == 0) {
                        markRow[iRow] = 1;
                        rowsToLookAtNext[numberLookNext++] = iRow;
                      }
                    }
                    /* --------------------------------------------------*/
                    /*                check if infeasible !!!!! */
                    /* --------------------------------------------------*/
                    if (colUpper[j] - colLower[j] < -100.0 * tolerance) {
                      ninfeas++;
                    }
                  }
                }
              }
            }
          }

          /* ----------------------------------------------------------------*/
          /*        Finite L(i) */
          /* ----------------------------------------------------------------*/
          if (iflagl == 0 && rowUpper[i] < 1e15) {
            for (k = krs; k < kre; ++k) {
              double value = rowElements[k];
              j = column[k];
              if (value < 0.0) {
                if (colUpper[j] < 1.0e12) {
                  dbound = colUpper[j] + (rowUpper[i] - dmaxdown) / value;
                  if (dbound > colLower[j] + 1.0e-8) {
                    colLower[j] = dbound;
                    ++nchange;
                    if (!markCol[j]) {
                      markCol[j] = 1;
                      whichColumn[numberBoundsChanged++] = j;
                    }
                    // see what rows need redoing
                    for (CoinBigIndex k = columnStart[j];
                      k < columnStart[j] + columnLength[j]; k++) {
                      int iRow = row[k];
                      if (markRow[iRow] == 0) {
                        markRow[iRow] = 1;
                        rowsToLookAtNext[numberLookNext++] = iRow;
                      }
                    }
                    /* --------------------------------------------------*/
                    /*                check if infeasible !!!!! */
                    /* --------------------------------------------------*/
                    if (colUpper[j] - colLower[j] < -100.0 * tolerance) {
                      ninfeas++;
                    }
                  }
                }
              } else {
                if (colLower[j] > -1.0e12) {
                  dbound = colLower[j] + (rowUpper[i] - dmaxdown) / value;
                  if (dbound < colUpper[j] - 1.0e-8) {
                    colUpper[j] = dbound;
                    ++nchange;
                    if (!markCol[j]) {
                      markCol[j] = 1;
                      whichColumn[numberBoundsChanged++] = j;
                    }
                    // see what rows need redoing
                    for (CoinBigIndex k = columnStart[j];
                      k < columnStart[j] + columnLength[j]; k++) {
                      int iRow = row[k];
                      if (markRow[iRow] == 0) {
                        markRow[iRow] = 1;
                        rowsToLookAtNext[numberLookNext++] = iRow;
                      }
                    }
                    /* --------------------------------------------------*/
                    /*                check if infeasible !!!!! */
                    /* --------------------------------------------------*/
                    if (colUpper[j] - colLower[j] < -100.0 * tolerance) {
                      ninfeas++;
                    }
                  }
                }
              }
            }
          }
        }
      }
      // reset markRow
      for (int j = 0; j < numberLookNext; j++) {
        markRow[rowsToLookAtNext[j]] = 0;
      }
      // std::sort(rowsToLookAtNext,rowsToLookAtNext+numberLookNext);
      for (int k = 0; k < numberBoundsChanged; k++) {
        int j = whichColumn[k];
        markCol[j] = 0;
        if (intVar[j]) {
          if (colUpper[j] - colLower[j] > 1.0e-8) {
            if (floor(colUpper[j] + 1.0e-4) < colUpper[j])
              nchange++;
            // clean up anyway
            colUpper[j] = floor(colUpper[j] + 1.0e-4);
            if (ceil(colLower[j] - 1.0e-4) > colLower[j])
              nchange++;
            // clean up anyway
            colLower[j] = ceil(colLower[j] - 1.0e-4);
            if (colUpper[j] < colLower[j]) {
              /*printf("infeasible\n");*/
              ninfeas++;
            }
          } else {
            // clean
            colUpper[j] = floor(colUpper[j] + 1.0e-4);
            colLower[j] = ceil(colLower[j] - 1.0e-4);
            if (colUpper[j] < colLower[j]) {
              /*printf("infeasible\n");*/
              ninfeas++;
            }
          }
        }
      }
      if (ninfeas)
        break;
    }
    // if (ninfeas) printf("%d passes\n",jpass);
    return (ninfeas);
  }

  while (nchange) {
    int ilbred = 0; /* bounds reduced */
    int iubred = 0; /* bounds reduced */
    int nrwdrp = 0; /* redundant rows */
    if (jpass == maxpass)
      break;
    jpass++;
    dolrows = (jpass & 1) == 1;

    for (i = 0; i < nRows; ++i) {
      bool cliqueChanges = false;
      if (rowLower[i] > -1.0e20 || rowUpper[i] < 1.0e20) {
        iflagu = 0;
        iflagl = 0;
        dmaxup = 0.0;
        dmaxdown = 0.0;
        krs = rowStart[i];
        kre = rowStart[i] + rowLength[i];

        /* ------------------------------------------------------------*/
        /* Compute L(i) and U(i) */
        /* ------------------------------------------------------------*/
        if (!cliqueMin || i >= numberRows_ || cliqueRowStart_[i] == cliqueRowStart_[i + 1]) {
          // without cliques
          for (k = krs; k < kre; ++k) {
            double value = rowElements[k];
            j = column[k];
            if (value > 0.0) {
              if (colUpper[j] < 1.0e12)
                dmaxup += colUpper[j] * value;
              else
                ++iflagu;
              if (colLower[j] > -1.0e12)
                dmaxdown += colLower[j] * value;
              else
                ++iflagl;
            } else if (value < 0.0) {
              if (colUpper[j] < 1.0e12)
                dmaxdown += colUpper[j] * value;
              else
                ++iflagl;
              if (colLower[j] > -1.0e12)
                dmaxup += colLower[j] * value;
              else
                ++iflagu;
            }
          }
          if (iflagu)
            dmaxup = 1.0e31;
          if (iflagl)
            dmaxdown = -1.0e31;
        } else {
          // with cliques
          int nClique = 0;
          CoinBigIndex bias = cliqueRowStart_[i] - krs;
          double dmaxup2 = 0.0;
          double dmaxdown2 = 0.0;
          double sumZeroFixes = 0.0;
          for (k = krs; k < kre; ++k) {
            double value = rowElements[k];
            j = column[k];
            int iClique = sequenceInCliqueEntry(cliqueRow_[k + bias]);
            bool oneFixes = oneFixesInCliqueEntry(cliqueRow_[k + bias]);
            if (iClique >= numberColumns_ || colUpper[j] == colLower[j]) {
              if (value > 0.0) {
                if (colUpper[j] >= 1.0e12) {
                  dmaxup = 1e31;
                  ++iflagu;
                } else {
                  dmaxup += colUpper[j] * value;
                }
                if (colLower[j] <= -1.0e12) {
                  dmaxdown = -1e31;
                  ++iflagl;
                } else {
                  dmaxdown += colLower[j] * value;
                }
              } else if (value < 0.0) {
                if (colUpper[j] >= 1.0e12) {
                  dmaxdown = -1e31;
                  ++iflagl;
                } else {
                  dmaxdown += colUpper[j] * value;
                }
                if (colLower[j] <= -1.0e12) {
                  dmaxup = 1e31;
                  ++iflagu;
                } else {
                  dmaxup += colLower[j] * value;
                }
              }
            } else {
              // clique may help
              if (iClique >= nClique) {
                // zero out
                for (int j = nClique; j <= iClique; j++) {
                  cliqueMin[j] = 0.0;
                  cliqueMax[j] = 0.0;
                  cliqueMin2[j] = 0.0;
                  cliqueMax2[j] = 0.0;
                }
                nClique = iClique + 1;
              }
              //  Update best and second best
              if (oneFixes) {
                if (value > 0.0) {
                  dmaxup2 += value;
                  cliqueMax2[iClique] = cliqueMax[iClique];
                  cliqueMax[iClique] = std::max(cliqueMax[iClique], value);
                } else if (value < 0.0) {
                  dmaxdown2 += value;
                  cliqueMin2[iClique] = cliqueMin[iClique];
                  cliqueMin[iClique] = std::min(cliqueMin[iClique], value);
                }
              } else {
                sumZeroFixes += value;
                if (value > 0.0) {
                  dmaxup2 += value;
                  cliqueMin2[iClique] = cliqueMin[iClique];
                  cliqueMin[iClique] = std::min(cliqueMin[iClique], -value);
                } else if (value < 0.0) {
                  dmaxdown2 += value;
                  cliqueMax2[iClique] = cliqueMax[iClique];
                  cliqueMax[iClique] = std::max(cliqueMax[iClique], -value);
                }
              }
            }
          }
          double dmaxup3 = dmaxup + sumZeroFixes;
          double dmaxdown3 = dmaxdown + sumZeroFixes;
          for (int iClique = 0; iClique < nClique; iClique++) {
            dmaxup3 += cliqueMax[iClique];
            dmaxdown3 += cliqueMin[iClique];
          }
          dmaxup += dmaxup2;
          dmaxdown += dmaxdown2;
          assert(dmaxup3 <= dmaxup + 1.0e-8);
          assert(dmaxdown3 >= dmaxdown - 1.0e-8);
          if (dmaxup3 < dmaxup - 1.0e-8 || dmaxdown3 > dmaxdown + 1.0e-8) {
            cliqueChanges = true;
            // printf("normal min/max %g , %g clique %g , %g\n",
            //      dmaxdown,dmaxup,dmaxdown3,dmaxup3);
            dmaxdown = dmaxdown3;
            dmaxup = dmaxup3;
          }
        }
        if (dmaxup <= rowUpper[i] + tolerance && dmaxdown >= rowLower[i] - tolerance) {
          /*
           * The sum of the column maxs is at most the row ub, and
           * the sum of the column mins is at least the row lb;
           * this row says nothing at all.
           * I suspect that this corresponds to
           * an implied column singleton in the paper (viii, on p. 325),
           * where the singleton in question is the row slack.
           */
          ++nrwdrp;
          rowLower[i] = -COIN_DBL_MAX;
          rowUpper[i] = COIN_DBL_MAX;
        } else {
          if (dmaxup < rowLower[i] - tolerance || dmaxdown > rowUpper[i] + tolerance) {
            ninfeas++;
            assert(!cliqueChanges);
            break;
          }
          /*        Finite U(i) */
          /* -------------------------------------------------------------*/
          /* below is deliberate mistake (previously was by chance) */
          /*        never do both */
          if (iflagu == 0 && rowLower[i] > 0.0 && iflagl == 0 && rowUpper[i] < 1e15) {
            if (dolrows) {
              iflagu = 1;
            } else {
              iflagl = 1;
            }
          }
          if (!cliqueChanges) {
            // without cliques
            if (iflagu == 0 && rowLower[i] > -1e15) {
              for (k = krs; k < kre; ++k) {
                double value = rowElements[k];
                j = column[k];
                if (value > 0.0) {
                  if (colUpper[j] < 1.0e12) {
                    dbound = colUpper[j] + (rowLower[i] - dmaxup) / value;
                    if (dbound > colLower[j] + 1.0e-8) {
                      /* we can tighten the lower bound */
                      /* the paper mentions this as a possibility on p. 227 */
                      colLower[j] = dbound;
                      ++ilbred;

                      /* this may have fixed the variable */
                      /* I believe that this roughly corresponds to a
                       * forcing constraint in the paper (p. 226).
                       * If there is a forcing constraint (with respect
                       * to the original, untightened bounds), then in this
                       * loop we will go through all the columns and fix
                       * each of them to their implied bound, rather than
                       * determining that the row as a whole is forced
                       * and just fixing them without doing computation for
                       * each column (as in the paper).
                       * By doing it this way, we can tighten bounds and
                       * get forcing constraints for free.
                       */
                      if (colUpper[j] - colLower[j] <= tolerance) {
                        /* --------------------------------------------------*/
                        /*                check if infeasible !!!!! */
                        /* --------------------------------------------------*/
                        if (colUpper[j] - colLower[j] < -100.0 * tolerance) {
                          ninfeas++;
                        }
                      }
                    }
                  }
                } else {
                  if (colLower[j] > -1.0e12) {
                    dbound = colLower[j] + (rowLower[i] - dmaxup) / value;
                    if (dbound < colUpper[j] - 1.0e-8) {
                      colUpper[j] = dbound;
                      ++iubred;
                      if (colUpper[j] - colLower[j] <= tolerance) {
                        /* --------------------------------------------------*/
                        /*                check if infeasible !!!!! */
                        /* --------------------------------------------------*/
                        if (colUpper[j] - colLower[j] < -100.0 * tolerance) {
                          ninfeas++;
                        }
                      }
                    }
                  }
                }
              }
            }

            /* ----------------------------------------------------------------*/
            /*        Finite L(i) */
            /* ----------------------------------------------------------------*/
            if (iflagl == 0 && rowUpper[i] < 1e15) {
              for (k = krs; k < kre; ++k) {
                double value = rowElements[k];
                j = column[k];
                if (value < 0.0) {
                  if (colUpper[j] < 1.0e12) {
                    dbound = colUpper[j] + (rowUpper[i] - dmaxdown) / value;
                    if (dbound > colLower[j] + 1.0e-8) {
                      colLower[j] = dbound;
                      ++ilbred;
                      if (!(colUpper[j] - colLower[j] > tolerance)) {
                        /* --------------------------------------------------*/
                        /*                check if infeasible !!!!! */
                        /* --------------------------------------------------*/
                        if (colUpper[j] - colLower[j] < -100.0 * tolerance) {
                          ninfeas++;
                        }
                      }
                    }
                  }
                } else {
                  if (colLower[j] > -1.0e12) {
                    dbound = colLower[j] + (rowUpper[i] - dmaxdown) / value;
                    if (dbound < colUpper[j] - 1.0e-8) {
                      colUpper[j] = dbound;
                      ++iubred;
                      if (!(colUpper[j] - colLower[j] > tolerance)) {
                        /* --------------------------------------------------*/
                        /*                check if infeasible !!!!! */
                        /* --------------------------------------------------*/
                        if (colUpper[j] - colLower[j] < -100.0 * tolerance) {
                          ninfeas++;
                        }
                      }
                    }
                  }
                }
              }
            }
          } else {
            // with cliques
            CoinBigIndex bias = cliqueRowStart_[i] - krs;
            if (iflagu == 0 && rowLower[i] > -1e15) {
              for (k = krs; k < kre; ++k) {
                double value = rowElements[k];
                j = column[k];
                int iClique = sequenceInCliqueEntry(cliqueRow_[k + bias]);
                // bool oneFixes = (cliqueRow_[k+bias].oneFixes!=0);
                if (iClique >= numberColumns_) {
                  if (value > 0.0) {
                    if (colUpper[j] < 1.0e12) {
                      dbound = colUpper[j] + (rowLower[i] - dmaxup) / value;
                      if (dbound > colLower[j] + 1.0e-8) {
                        /* we can tighten the lower bound */
                        /* the paper mentions this as a possibility on p. 227 */
                        colLower[j] = dbound;
                        ++ilbred;

                        /* this may have fixed the variable */
                        /* I believe that this roughly corresponds to a
                         * forcing constraint in the paper (p. 226).
                         * If there is a forcing constraint (with respect
                         * to the original, untightened bounds), then in this
                         * loop we will go through all the columns and fix
                         * each of them to their implied bound, rather than
                         * determining that the row as a whole is forced
                         * and just fixing them without doing computation for
                         * each column (as in the paper).
                         * By doing it this way, we can tighten bounds and
                         * get forcing constraints for free.
                         */
                        if (colUpper[j] - colLower[j] <= tolerance) {
                          /* --------------------------------------------------*/
                          /*                check if infeasible !!!!! */
                          /* --------------------------------------------------*/
                          if (colUpper[j] - colLower[j] < -100.0 * tolerance) {
                            ninfeas++;
                          }
                        }
#if 0
                      } else if (intVar[j]==1 && rowUpper[i]>1.0e20) {
                        // can we modify coefficient
                        if (dmaxdown+value>rowLower[i]+1.0e-8) {
                          assert (dmaxdown<rowLower[i]+1.0e-8);
                          double change = dmaxdown+value - rowLower[i];
                          double newValue = value - change;
                          if (newValue<1.0e-12)
                            newValue=0.0;
                          printf("Could change value from %g to %g\n",
                                 value,newValue);
                          // dmaxup -= change; 
                        }
#endif
                      }
                    }
                  } else {
                    if (colLower[j] > -1.0e12) {
                      dbound = colLower[j] + (rowLower[i] - dmaxup) / value;
                      if (dbound < colUpper[j] - 1.0e-8) {
                        colUpper[j] = dbound;
                        ++iubred;
                        if (colUpper[j] - colLower[j] <= tolerance) {
                          /* --------------------------------------------------*/
                          /*                check if infeasible !!!!! */
                          /* --------------------------------------------------*/
                          if (colUpper[j] - colLower[j] < -100.0 * tolerance) {
                            ninfeas++;
                          }
                        }
#if 0
                      } else if (intVar[j]==1 && rowUpper[i]>1.0e20) {
                        // can we modify coefficient
                        if (dmaxdown-value>rowLower[i]+1.0e-8) {
                          assert (dmaxdown<rowLower[i]+1.0e-8);
                          double change = dmaxdown-value-rowLower[i];
                          double newValue = value+change;
                          double newLower = rowLower[i]+change;
                          if (newValue<1.0e-12)
                            newValue=0.0;
                          printf("Could change value from %g to %g and lorow from %g to %g\n",
                                 value,newValue,rowLower[i],newLower);
                         // dmaxdown += change                          
                        }
#endif
                      }
                    }
                  }
                } else if (colUpper[j] > colLower[j]) {
                  // in clique
                  // adjustment
                  double dmaxup2 = dmaxup;
                  assert(cliqueMax[iClique] >= 0);
                  assert(cliqueMax2[iClique] >= 0);
                  /* get max up if at other bound
                     May not go down at all but will not go up */
                  if (fabs(value) == fabs(cliqueMax[iClique]))
                    dmaxup2 -= cliqueMax[iClique] - cliqueMax2[iClique];
                  if (dmaxup2 < rowLower[i] - 1.0e-8) {
                    /* --------------------------------------------------*/
                    /*                check if infeasible !!!!! */
                    /* --------------------------------------------------*/
                    if (dmaxup < rowLower[i] - 1.0e-8) {
                      ninfeas++;
                    } else {
                      if (value > 0.0) {
                        colLower[j] = 1.0;
                        ++ilbred;
                      } else {
                        colUpper[j] = 0.0;
                        ++iubred;
                      }
                    }
                  }
                }
              }
            }

            /* ----------------------------------------------------------------*/
            /*        Finite L(i) */
            /* ----------------------------------------------------------------*/
            if (iflagl == 0 && rowUpper[i] < 1e15) {
              for (k = krs; k < kre; ++k) {
                double value = rowElements[k];
                j = column[k];
                int iClique = sequenceInCliqueEntry(cliqueRow_[k + bias]);
                // bool oneFixes = (cliqueRow_[k+bias].oneFixes!=0);
                if (iClique >= numberColumns_) {
                  if (value < 0.0) {
                    if (colUpper[j] < 1.0e12) {
                      dbound = colUpper[j] + (rowUpper[i] - dmaxdown) / value;
                      if (dbound > colLower[j] + 1.0e-8) {
                        colLower[j] = dbound;
                        ++ilbred;
                        if (!(colUpper[j] - colLower[j] > tolerance)) {
                          /* --------------------------------------------------*/
                          /*                check if infeasible !!!!! */
                          /* --------------------------------------------------*/
                          if (colUpper[j] - colLower[j] < -100.0 * tolerance) {
                            ninfeas++;
                          }
                        }
#if 0
                      } else if (intVar[j]==1 && rowLower[i]<-1.0e20) {
                        // can we modify coefficient
                        if (dmaxup+value<rowUpper[i]-1.0e-8) {
                          assert (dmaxup>rowUpper[i]-1.0e-8);
                          double change = dmaxup+value - rowUpper[i];
                          double newValue = value - change;
                          if (newValue<1.0e-12)
                            newValue=0.0;
                          printf("Could change value from %g to %g b\n",
                                 value,newValue);
                          // dmaxdown -= change; 
                        }
#endif
                      }
                    }
                  } else {
                    if (colLower[j] > -1.0e12) {
                      dbound = colLower[j] + (rowUpper[i] - dmaxdown) / value;
                      if (dbound < colUpper[j] - 1.0e-8) {
                        colUpper[j] = dbound;
                        ++iubred;
                        if (!(colUpper[j] - colLower[j] > tolerance)) {
                          /* --------------------------------------------------*/
                          /*                check if infeasible !!!!! */
                          /* --------------------------------------------------*/
                          if (colUpper[j] - colLower[j] < -100.0 * tolerance) {
                            ninfeas++;
                          }
                        }
#if 0
                      } else if (intVar[j]==1 && rowLower[i]<-1.0e20) {
                        // can we modify coefficient
                        if (dmaxup-value<rowUpper[i]-1.0e-8) {
                          assert (dmaxup>rowUpper[i]-1.0e-8);
                          double change = dmaxup-value-rowUpper[i];
                          double newValue = value+change;
                          double newUpper = rowUpper[i]+change;
                          if (newValue<1.0e-12)
                            newValue=0.0;
                          printf("Could change value from %g to %g and uprow from %g to %g b\n",
                                 value,newValue,rowLower[i],newUpper);
                         // dmaxup += change                          
                        }
#endif
                      }
                    }
                  }
                } else if (colUpper[j] > colLower[j]) {
                  // in clique
                  // adjustment
                  double dmaxdown2 = dmaxdown;
                  assert(cliqueMin[iClique] <= 0);
                  assert(cliqueMin2[iClique] <= 0);
                  /* get max down if this is at other bound
                     May not go up at all but will not go down */
                  if (fabs(value) == fabs(cliqueMin[iClique]))
                    dmaxdown2 -= cliqueMin[iClique] - cliqueMin2[iClique];
                  if (dmaxdown2 > rowUpper[i] + 1.0e-8) {
                    /* --------------------------------------------------*/
                    /*                check if infeasible !!!!! */
                    /* --------------------------------------------------*/
                    if (dmaxdown > rowUpper[i] + 1.0e-8) {
                      ninfeas++;
                    } else {
                      if (value < 0.0) {
                        colLower[j] = 1.0;
                        ++ilbred;
                      } else {
                        colUpper[j] = 0.0;
                        ++iubred;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    for (j = 0; j < nCols; ++j) {
      if (intVar[j]) {
        if (colUpper[j] - colLower[j] > 1.0e-8) {
          if (floor(colUpper[j] + 1.0e-4) < colUpper[j])
            nchange++;
          // clean up anyway
          colUpper[j] = floor(colUpper[j] + 1.0e-4);
          if (ceil(colLower[j] - 1.0e-4) > colLower[j])
            nchange++;
          // clean up anyway
          colLower[j] = ceil(colLower[j] - 1.0e-4);
          if (colUpper[j] < colLower[j]) {
            /*printf("infeasible\n");*/
            ninfeas++;
          }
        }
      }
    }
    nchange = ilbred + iubred + nrwdrp;
    // ntotal += nchange;
    if (ninfeas)
      break;
  }
  delete[] cliqueMin;
  delete[] cliqueMax;
  delete[] cliqueMin2;
  delete[] cliqueMax2;
  return (ninfeas);
}
// This just sets minima and maxima on rows
void CglProbing::tighten2(double *colLower, double *colUpper,
  const int *column, const double *rowElements,
  const CoinBigIndex *rowStart,
  const int *rowLength,
  double *rowLower, double *rowUpper,
  double *minR, double *maxR, int *markR,
  int nRows)
{
  int i, j;
  CoinBigIndex k;
  CoinBigIndex krs, kre;
  int iflagu, iflagl;
  double dmaxup, dmaxdown;

  for (i = 0; i < nRows; ++i) {
    if (rowLower[i] > -1.0e20 || rowUpper[i] < 1.0e20) {
      iflagu = 0;
      iflagl = 0;
      dmaxup = 0.0;
      dmaxdown = 0.0;
      krs = rowStart[i];
      kre = rowStart[i] + rowLength[i];

      /* ------------------------------------------------------------*/
      /* Compute L(i) and U(i) */
      /* ------------------------------------------------------------*/
      for (k = krs; k < kre; ++k) {
        double value = rowElements[k];
        j = column[k];
        if (value > 0.0) {
          if (colUpper[j] < 1.0e12)
            dmaxup += colUpper[j] * value;
          else
            ++iflagu;
          if (colLower[j] > -1.0e12)
            dmaxdown += colLower[j] * value;
          else
            ++iflagl;
        } else if (value < 0.0) {
          if (colUpper[j] < 1.0e12)
            dmaxdown += colUpper[j] * value;
          else
            ++iflagl;
          if (colLower[j] > -1.0e12)
            dmaxup += colLower[j] * value;
          else
            ++iflagu;
        }
      }
      if (iflagu)
        maxR[i] = 1.0e60;
      else
        maxR[i] = dmaxup;
      if (iflagl)
        minR[i] = -1.0e60;
      else
        minR[i] = dmaxdown;
#if 0
      if (minR[i]<-1.0e10&&maxR[i]>1.0e10) {
	markR[i]=-2;
      } else {
#endif
      markR[i] = -1;
#if 0
      }
#endif
    } else {
      minR[i] = -1.0e60;
      maxR[i] = 1.0e60;
#if 0
      markR[i]=-2;
      abort();
#else
      markR[i] = -1;
#endif
    }
  }
}
#ifdef CGL_DEBUG
static int nPath = 0;
#endif
//-------------------------------------------------------------------
// Generate disaggregation cuts
//-------------------------------------------------------------------
void CglProbing::generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
  const CglTreeInfo info2)
{
  si.getColType(true);

#ifdef CGL_DEBUG
  const OsiRowCutDebugger *debugger = si.getRowCutDebugger();
  if (debugger && debugger->onOptimalPath(si)) {
    printf("On optimal path %d\n", nPath);
    nPath++;
    int nCols = si.getNumCols();
    int i;
    const double *solution = si.getColSolution();
    const double *lower = si.getColLower();
    const double *upper = si.getColUpper();
    const double *optimal = debugger->optimalSolution();
    const double *objective = si.getObjCoefficients();
    double objval1 = 0.0, objval2 = 0.0;
    for (i = 0; i < nCols; i++) {
#if CGL_DEBUG > 1
      printf("%d %g %g %g %g\n", i, lower[i], solution[i], upper[i], optimal[i]);
#endif
      objval1 += solution[i] * objective[i];
      objval2 += optimal[i] * objective[i];
      assert(optimal[i] >= lower[i] && optimal[i] <= upper[i]);
    }
    printf("current obj %g, integer %g\n", objval1, objval2);
  }
#endif
  int saveRowCuts = rowCuts_;
  if (rowCuts_ < 0) {
    if (info2.inTree)
      rowCuts_ = 4;
    else
      rowCuts_ = -rowCuts_;
  }
  int nRows = si.getNumRows();
  double *rowLower = new double[nRows + 1];
  double *rowUpper = new double[nRows + 1];

  int nCols = si.getNumCols();
  // Set size if not set
  if (!rowCopy_) {
    numberRows_ = nRows;
    numberColumns_ = nCols;
  }
  double *colLower = new double[nCols];
  double *colUpper = new double[nCols];

  CglTreeInfo info = info2;
  int ninfeas = gutsOfGenerateCuts(si, cs, rowLower, rowUpper, colLower, colUpper, &info);
  if (ninfeas) {
    // generate infeasible cut and return
    OsiRowCut rc;
    rc.setLb(COIN_DBL_MAX);
    rc.setUb(0.0);
    cs.insert(rc);
#ifdef CGL_DEBUG
    const OsiRowCutDebugger *debugger = si.getRowCutDebugger();
    if (debugger && debugger->onOptimalPath(si))
      assert(!debugger->invalidCut(rc));
#endif
  }
  delete[] rowLower;
  delete[] rowUpper;
  delete[] colLower;
  delete[] colUpper;
  delete[] colLower_;
  delete[] colUpper_;
  colLower_ = NULL;
  colUpper_ = NULL;
  rowCuts_ = saveRowCuts;
  return;
}
int CglProbing::generateCutsAndModify(const OsiSolverInterface &si,
  OsiCuts &cs,
  CglTreeInfo *info)
{
#ifdef CGL_DEBUG
  const OsiRowCutDebugger *debugger = si.getRowCutDebugger();
  if (debugger && debugger->onOptimalPath(si)) {
    printf("On optimal path %d\n", nPath);
    nPath++;
    int nCols = si.getNumCols();
    int i;
    const double *solution = si.getColSolution();
    const double *lower = si.getColLower();
    const double *upper = si.getColUpper();
    const double *optimal = debugger->optimalSolution();
    const double *objective = si.getObjCoefficients();
    double objval1 = 0.0, objval2 = 0.0;
    for (i = 0; i < nCols; i++) {
#if CGL_DEBUG > 1
      printf("%d %g %g %g %g\n", i, lower[i], solution[i], upper[i], optimal[i]);
#endif
      objval1 += solution[i] * objective[i];
      objval2 += optimal[i] * objective[i];
      assert(optimal[i] >= lower[i] - 1.0e-5 && optimal[i] <= upper[i] + 1.0e-5);
    }
    printf("current obj %g, integer %g\n", objval1, objval2);
  }
#endif
  int saveRowCuts = rowCuts_;
  if (rowCuts_ < 0) {
    if (info->inTree)
      rowCuts_ = 4;
    else
      rowCuts_ = -rowCuts_;
  }
  int saveMode = mode_;
  bool rowCliques = false;
  if (!(mode_ & 15)) {
    if (info->pass != 4 || info->inTree) {
      mode_ = 1;
    } else {
      saveMode = 1; // make sure do just once
      mode_ = 0;
      rowCliques = true;
    }
  }
  int nRows = si.getNumRows();
  double *rowLower = new double[nRows + 1];
  double *rowUpper = new double[nRows + 1];

  int nCols = si.getNumCols();
  double *colLower = new double[nCols];
  double *colUpper = new double[nCols];

  int ninfeas = gutsOfGenerateCuts(si, cs, rowLower, rowUpper, colLower, colUpper, info);
  if (ninfeas) {
    // generate infeasible cut and return
    OsiRowCut rc;
    rc.setLb(COIN_DBL_MAX);
    rc.setUb(0.0);
    cs.insert(rc);
#ifdef CGL_DEBUG
    const OsiRowCutDebugger *debugger = si.getRowCutDebugger();
    if (debugger && debugger->onOptimalPath(si))
      assert(!debugger->invalidCut(rc));
#endif
  }
  rowCuts_ = saveRowCuts;
  mode_ = saveMode;
  // move bounds so can be used by user
  if (mode_ == 3) {
    delete[] rowLower_;
    delete[] rowUpper_;
    rowLower_ = rowLower;
    rowUpper_ = rowUpper;
  } else {
    delete[] rowLower;
    delete[] rowUpper;
  }
  delete[] colLower_;
  delete[] colUpper_;
  colLower_ = colLower;
  colUpper_ = colUpper;
  // Setup information
  if (rowCliques && numberRows_ && numberColumns_)
    setupRowCliqueInformation(si);
  return ninfeas;
}
#if 0 // moved to CglPreProcess
bool analyze(const OsiSolverInterface * solverX, char * intVar,
             double * lower, double * upper)
{
  OsiSolverInterface * solver = solverX->clone();
  const double *objective = solver->getObjCoefficients() ;
  int numberColumns = solver->getNumCols() ;
  int numberRows = solver->getNumRows();
  double direction = solver->getObjSense();
  int iRow,iColumn;

  // Row copy
  CoinPackedMatrix matrixByRow(*solver->getMatrixByRow());
  const double * elementByRow = matrixByRow.getElements();
  const int * column = matrixByRow.getIndices();
  const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
  const int * rowLength = matrixByRow.getVectorLengths();

  // Column copy
  CoinPackedMatrix  matrixByCol(*solver->getMatrixByCol());
  const double * element = matrixByCol.getElements();
  const int * row = matrixByCol.getIndices();
  const CoinBigIndex * columnStart = matrixByCol.getVectorStarts();
  const int * columnLength = matrixByCol.getVectorLengths();

  const double * rowLower = solver->getRowLower();
  const double * rowUpper = solver->getRowUpper();

  char * ignore = new char [numberRows];
  int * which = new int[numberRows];
  double * changeRhs = new double[numberRows];
  memset(changeRhs,0,numberRows*sizeof(double));
  memset(ignore,0,numberRows);
  int numberChanged=0;
  bool finished=false;
  while (!finished) {
    int saveNumberChanged = numberChanged;
    for (iRow=0;iRow<numberRows;iRow++) {
      int numberContinuous=0;
      double value1=0.0,value2=0.0;
      bool allIntegerCoeff=true;
      double sumFixed=0.0;
      int jColumn1=-1,jColumn2=-1;
      for (CoinBigIndex j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
        int jColumn = column[j];
        double value = elementByRow[j];
        if (upper[jColumn] > lower[jColumn]+1.0e-8) {
          if (!intVar[jColumn]) {
            if (numberContinuous==0) {
              jColumn1=jColumn;
              value1=value;
            } else {
              jColumn2=jColumn;
              value2=value;
            }
            numberContinuous++;
          } else {
            if (fabs(value-floor(value+0.5))>1.0e-12)
              allIntegerCoeff=false;
          }
        } else {
          sumFixed += lower[jColumn]*value;
        }
      }
      double low = rowLower[iRow];
      if (low>-1.0e20) {
        low -= sumFixed;
        if (fabs(low-floor(low+0.5))>1.0e-12)
          allIntegerCoeff=false;
      }
      double up = rowUpper[iRow];
      if (up<1.0e20) {
        up -= sumFixed;
        if (fabs(up-floor(up+0.5))>1.0e-12)
          allIntegerCoeff=false;
      }
      if (!allIntegerCoeff)
        continue; // can't do
      if (numberContinuous==1) {
        // see if really integer
        // This does not allow for complicated cases
        if (low==up) {
          if (fabs(value1)>1.0e-3) {
            value1 = 1.0/value1;
            if (fabs(value1-floor(value1+0.5))<1.0e-12) {
              // integer
              numberChanged++;
              intVar[jColumn1]=77;
            }
          }
        } else {
          if (fabs(value1)>1.0e-3) {
            value1 = 1.0/value1;
            if (fabs(value1-floor(value1+0.5))<1.0e-12) {
              // This constraint will not stop it being integer
              ignore[iRow]=1;
            }
          }
        }
      } else if (numberContinuous==2) {
        if (low==up) {
          /* need general theory - for now just look at 2 cases -
             1 - +- 1 one in column and just costs i.e. matching objective
             2 - +- 1 two in column but feeds into G/L row which will try and minimize
             (take out 2 for now - until fixed)
          */
          if (fabs(value1)==1.0&&value1*value2==-1.0&&!lower[jColumn1]
              &&!lower[jColumn2]&&columnLength[jColumn1]==1&&columnLength[jColumn2]==1) {
            int n=0;
            CoinBigIndex i;
            double objChange=direction*(objective[jColumn1]+objective[jColumn2]);
            double bound = std::min(upper[jColumn1],upper[jColumn2]);
            bound = std::min(bound,1.0e20);
            for ( i=columnStart[jColumn1];i<columnStart[jColumn1]+columnLength[jColumn1];i++) {
              int jRow = row[i];
              double value = element[i];
              if (jRow!=iRow) {
                which[n++]=jRow;
                changeRhs[jRow]=value;
              }
            }
            for ( i=columnStart[jColumn2];i<columnStart[jColumn2]+columnLength[jColumn2];i++) {
              int jRow = row[i];
              double value = element[i];
              if (jRow!=iRow) {
                if (!changeRhs[jRow]) {
                  which[n++]=jRow;
                  changeRhs[jRow]=value;
                } else {
                  changeRhs[jRow]+=value;
                }
              }
            }
            if (objChange>=0.0) {
              // see if all rows OK
              bool good=true;
              for (i=0;i<n;i++) {
                int jRow = which[i];
                double value = changeRhs[jRow];
                if (value) {
                  value *= bound;
                  if (rowLength[jRow]==1) {
                    if (value>0.0) {
                      double rhs = rowLower[jRow];
                      if (rhs>0.0) {
                        double ratio =rhs/value;
                        if (fabs(ratio-floor(ratio+0.5))>1.0e-12)
                          good=false;
                      }
                    } else {
                      double rhs = rowUpper[jRow];
                      if (rhs<0.0) {
                        double ratio =rhs/value;
                        if (fabs(ratio-floor(ratio+0.5))>1.0e-12)
                          good=false;
                      }
                    }
                  } else if (rowLength[jRow]==2) {
                    if (value>0.0) {
                      if (rowLower[jRow]>-1.0e20)
                        good=false;
                    } else {
                      if (rowUpper[jRow]<1.0e20)
                        good=false;
                    }
                  } else {
                    good=false;
                  }
                }
              }
              if (good) {
                // both can be integer
                numberChanged++;
                intVar[jColumn1]=77;
                numberChanged++;
                intVar[jColumn2]=77;
              }
            }
            // clear
            for (i=0;i<n;i++) {
              changeRhs[which[i]]=0.0;
            }
          }
        }
      }
    }
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (upper[iColumn] > lower[iColumn]+1.0e-8&&!intVar[iColumn]) {
        double value;
        value = upper[iColumn];
        if (value<1.0e20&&fabs(value-floor(value+0.5))>1.0e-12) 
          continue;
        value = lower[iColumn];
        if (value>-1.0e20&&fabs(value-floor(value+0.5))>1.0e-12) 
          continue;
        bool integer=true;
        for (CoinBigIndex j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
          int iRow = row[j];
          if (!ignore[iRow]) {
            integer=false;
            break;
          }
        }
        if (integer) {
          // integer
          numberChanged++;
          intVar[iColumn]=77;
        }
      }
    }
    finished = numberChanged==saveNumberChanged;
  }
  bool feasible=true;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (intVar[iColumn]==77) {
      if (upper[iColumn]>1.0e20) {
        upper[iColumn] = 1.0e20;
      } else {
        upper[iColumn] = floor(upper[iColumn]+1.0e-5);
      }
      if (lower[iColumn]<-1.0e20) {
        lower[iColumn] = -1.0e20;
      } else {
        lower[iColumn] = ceil(lower[iColumn]-1.0e-5);
        if (lower[iColumn]>upper[iColumn])
          feasible=false;
      }
      if (lower[iColumn]==0.0&&upper[iColumn]==1.0)
        intVar[iColumn]=1;
      else if (lower[iColumn]==upper[iColumn])
        intVar[iColumn]=0;
      else
        intVar[iColumn]=2;
    }
  }
  delete [] which;
  delete [] changeRhs;
  delete [] ignore;
  //if (numberChanged)
  //printf("%d variables could be made integer\n",numberChanged);
  delete solver;
  return feasible;
}
#endif
int CglProbing::gutsOfGenerateCuts(const OsiSolverInterface &si,
  OsiCuts &cs,
  double *rowLower, double *rowUpper,
  double *colLower, double *colUpper,
  CglTreeInfo *info)
{
  // printf("PASS\n");
  //  Get basic problem information
  int nRows;

  CoinPackedMatrix *rowCopy = NULL;
  int numberRowCutsBefore = cs.sizeRowCuts();

  // get branch and bound cutoff
  double cutoff;
  bool cutoff_available = si.getDblParam(OsiDualObjectiveLimit, cutoff);
  if (!cutoff_available || usingObjective_ < 0) { // cut off isn't set or isn't valid
    cutoff = si.getInfinity();
  }
  cutoff *= si.getObjSense();
  if (fabs(cutoff) > 1.0e30)
    assert(cutoff > 1.0e30);
  int mode = mode_;
  int nCols = si.getNumCols();

  // get integer variables
  const char *intVarOriginal = si.getColType(true);
  char *intVar = CoinCopyOfArray(intVarOriginal, nCols);
  int i;
#ifndef NDEBUG
  int numberIntegers = 0;
#endif
  CoinMemcpyN(si.getColLower(), nCols, colLower);
  CoinMemcpyN(si.getColUpper(), nCols, colUpper);
  const double *colsol = si.getColSolution();
  // and put reasonable bounds on variables
  for (i = 0; i < nCols; i++) {
    colUpper[i] = std::min(colUpper[i], 1.0e29);
    colLower[i] = std::max(colLower[i], -1.0e29);
    if (intVar[i]) {
#ifndef NDEBUG
      numberIntegers++;
#endif
      if (intVar[i] == 2) {
        // make sure reasonable bounds
        if (colsol[i] < 1.0e10 && colUpper[i] > 1.0e12)
          colUpper[i] = CGL_REASONABLE_INTEGER_BOUND;
        if (colsol[i] > -1.0e10 && colLower[i] < -1.0e12)
          colLower[i] = -CGL_REASONABLE_INTEGER_BOUND;
      }
    }
  }
  if (!info->inTree && !info->pass) {
#if CBC_PREPROCESS_EXPERIMENT > 1
    int nRows = si.getNumRows();
    double *columnLo = new double[2 * nRows + 2 * nCols];
    double *columnUp = columnLo + nCols;
    double *rowLo = columnUp + nCols;
    double *rowUp = rowLo + nRows;
    // temp
    memcpy(columnLo, si.getColLower(), nCols * sizeof(double));
    memcpy(columnUp, si.getColUpper(), nCols * sizeof(double));
    memcpy(rowLo, si.getRowLower(), nRows * sizeof(double));
    memcpy(rowUp, si.getRowUpper(), nRows * sizeof(double));
    // Get a row copy in standard format
    CoinPackedMatrix rowCopy = *si.getMatrixByRow();
    int columnsTightened;
    int rowsTightened;
    int elementsChanged;
    si.tightPrimalBounds(rowLo, rowUp, columnLo, columnUp,
      columnsTightened, rowsTightened,
      elementsChanged, rowCopy);
    int nChanged = 0;
    CoinPackedVector lbs;
    CoinPackedVector ubs;
    for (int iColumn = 0; iColumn < nCols; iColumn++) {
      if (columnLo[iColumn] > colLower[iColumn]) {
        nChanged++;
        colLower[iColumn] = columnLo[iColumn];
        if (intVar[iColumn])
          lbs.insert(iColumn, columnLo[iColumn]);
      }
      if (columnUp[iColumn] < colUpper[iColumn]) {
        nChanged++;
        colUpper[iColumn] = columnUp[iColumn];
        if (intVar[iColumn])
          ubs.insert(iColumn, columnUp[iColumn]);
      }
    }
#ifdef LOTS_OF_PRINTING
    printf("%d bounds changed, %d lower int, %d upper int\n",
      nChanged, lbs.getNumElements(), ubs.getNumElements());
#endif
    if (lbs.getNumElements()) {
      OsiColCut cc;
      cc.setLbs(lbs);
      cs.insert(cc);
    }
    if (ubs.getNumElements()) {
      OsiColCut cc;
      cc.setUbs(ubs);
      cs.insert(cc);
    }
    if (rowCopy_) {
      // printf("Can't tighten row bounds\n");
    } else {
      nChanged = 0;
      for (int iRow = 0; iRow < nRows; iRow++) {
        if (rowLo[iRow] > rowLower[iRow]) {
          nChanged++;
          rowLower[iRow] = rowLo[iRow];
        }
        if (rowUp[iRow] < rowUpper[iRow]) {
          nChanged++;
          rowUpper[iRow] = rowUp[iRow];
        }
      }
#ifdef LOTS_OF_PRINTING
      printf("%d row bounds changed\n", nChanged);
#endif
    }
#endif
    // make more integer
    // feasible = analyze(&si,intVar,colLower,colUpper);
  }
  if (PROBING_EXTRA_STUFF) {
    // tighten bounds on djs
    // should be in CbcCutGenerator and check if basic
    const double *djs = si.getReducedCost();
    const double *colsol = si.getColSolution();
    double direction = si.getObjSense();
    double cutoff;
    bool cutoff_available = si.getDblParam(OsiDualObjectiveLimit, cutoff);
    if (!cutoff_available || usingObjective_ < 0) { // cut off isn't set or isn't valid
      cutoff = si.getInfinity();
    }
    cutoff *= direction;
    if (fabs(cutoff) > 1.0e30)
      assert(cutoff > 1.0e30);
    double current = si.getObjValue();
    current *= direction;
    double gap = std::max(cutoff - current, 1.0e-1);
    for (int i = 0; i < nCols; ++i) {
      double djValue = djs[i] * direction;
      if (colUpper[i] - colLower[i] > 1.0e-8) {
        if (colsol[i] < colLower[i] + primalTolerance_) {
          if (djValue > gap) {
            if (intVar[i]) {
              printf("why int %d not fixed at lb\n", i);
              colUpper[i] = colLower[i];
            } else {
              double newUpper = colLower[i] + gap / djValue;
              if (newUpper < colUpper[i]) {
                // printf("%d ub from %g to %g\n",i,colUpper[i],newUpper);
                colUpper[i] = std::max(newUpper, colLower[i] + 1.0e-5);
              }
            }
          }
        } else if (colsol[i] > colUpper[i] - primalTolerance_) {
          if (-djValue > gap) {
            if (intVar[i]) {
              printf("why int %d not fixed at ub\n", i);
              colLower[i] = colUpper[i];
            } else {
              double newLower = colUpper[i] + gap / djValue;
              if (newLower > colLower[i]) {
                // printf("%d lb from %g to %g\n",i,colLower[i],newLower);
                colLower[i] = std::min(newLower, colUpper[i] - 1.0e-5);
              }
            }
          }
        }
      }
    }
  }
  int ninfeas = 0;
  // Set up maxes
  int maxProbe = info->inTree ? maxProbe_ : maxProbeRoot_;
  int maxElements = info->inTree ? maxElements_ : maxElementsRoot_;
  // if (!info->inTree&&!info->pass)
  // maxElements=nCols;
  //  Get objective offset
  double offset;
  si.getDblParam(OsiObjOffset, offset);
#ifdef COIN_DEVELOP
  if (offset && !info->inTree && !info->pass)
    printf("CglProbing obj offset %g\n", offset);
#endif
  // see if using cached copy or not
  if (!rowCopy_) {
    // create from current
    nRows = si.getNumRows();

    // mode==0 is invalid if going from current matrix
    if (mode == 0)
      mode = 1;
    // add in objective if there is a cutoff
    if (cutoff < 1.0e30 && usingObjective_ > 0) {
      rowCopy = new CoinPackedMatrix(*si.getMatrixByRow(), 1, nCols, false);
    } else {
      rowCopy = new CoinPackedMatrix(*si.getMatrixByRow());
    }
    // add in objective if there is a cutoff
    if (cutoff < 1.0e30 && usingObjective_ > 0) {
      int *columns = new int[nCols];
      double *elements = new double[nCols];
      int n = 0;
      const double *objective = si.getObjCoefficients();
      bool maximize = (si.getObjSense() == -1);
      for (i = 0; i < nCols; i++) {
        if (objective[i]) {
          elements[n] = (maximize) ? -objective[i] : objective[i];
          columns[n++] = i;
        }
      }
      rowCopy->appendRow(n, columns, elements);
      delete[] columns;
      delete[] elements;
      CoinMemcpyN(si.getRowLower(), nRows, rowLower);
      CoinMemcpyN(si.getRowUpper(), nRows, rowUpper);
      rowLower[nRows] = -COIN_DBL_MAX;
      rowUpper[nRows] = cutoff + offset;
      nRows++;
    } else {
      CoinMemcpyN(si.getRowLower(), nRows, rowLower);
      CoinMemcpyN(si.getRowUpper(), nRows, rowUpper);
    }
  } else {
    // use snapshot
    nRows = numberRows_;
    assert(nCols == numberColumns_);

    rowCopy = new CoinPackedMatrix(*rowCopy_);
    assert(rowCopy_->getNumRows() == numberRows_);
    rowLower = new double[nRows];
    rowUpper = new double[nRows];
    CoinMemcpyN(rowLower_, nRows, rowLower);
    CoinMemcpyN(rowUpper_, nRows, rowUpper);
    if (usingObjective_ > 0) {
      rowLower[nRows - 1] = -COIN_DBL_MAX;
      rowUpper[nRows - 1] = cutoff + offset;
    }
  }
  CoinBigIndex *rowStartPos = NULL;
  int *realRows = NULL;
  {
    // Now take out rows with too many elements
    int *rowLength = rowCopy->getMutableVectorLengths();
    // #define OUTRUBBISH
    double *elements = rowCopy->getMutableElements();
    int *column = rowCopy->getMutableIndices();
    CoinBigIndex *rowStart = rowCopy->getMutableVectorStarts();
#ifdef OUTRUBBISH
    double large = 1.0e3;
#endif
    int nDelete = 0;
    int nKeep = 0;
    int *which = new int[nRows];
    CoinBigIndex nElements = rowCopy->getNumElements();
    int nTotalOut = 0;
    int nRealRows = si.getNumRows();
    // keep objective
    for (i = 0; i < nRealRows; i++) {
      if (rowLength[i] > maxElements || (rowLower[i] < -1.0e20 && rowUpper[i] > 1.0e20)) {
        nTotalOut += rowLength[i];
      }
    }
    // keep all if only a few dense
    if (nTotalOut * 10 < nElements)
      maxElements = nCols;
#ifdef OUTRUBBISH
    int nExtraDel = 0;
#endif
    for (i = 0; i < nRows; i++) {
      if ((rowLength[i] > maxElements && i < nRealRows) || (rowLower[i] < -1.0e20 && rowUpper[i] > 1.0e20)) {
        which[nDelete++] = i;
      } else {
#ifdef OUTRUBBISH
        // out all rows with infinite plus and minus
        int nPlus = rowUpper[i] > -large ? 0 : 1;
        int nMinus = rowLower[i] < large ? 0 : 1;
        CoinBigIndex start = rowStart[i];
        CoinBigIndex end = start + rowLength[i];
        for (CoinBigIndex j = start; j < end; j++) {
          int iColumn = column[j];
          if (colUpper[iColumn] > large) {
            if (elements[j] > 0.0)
              nPlus++;
            else
              nMinus++;
          }
          if (colLower[iColumn] < -large) {
            if (elements[j] < 0.0)
              nPlus++;
            else
              nMinus++;
          }
        }
        if (!nPlus || !nMinus) {
          rowLower[nKeep] = rowLower[i];
          rowUpper[nKeep] = rowUpper[i];
          nKeep++;
        } else {
          nExtraDel++;
          which[nDelete++] = i;
        }
#else
        if (info->strengthenRow && !info->pass && (rowLower[i] < -1.0e20 || rowUpper[i] > 1.0e20)) {
          int nPlus = 0;
          int nMinus = 0;
          for (CoinBigIndex j = rowStart[i]; j < rowStart[i] + rowLength[i]; j++) {
            int jColumn = column[j];
            if (intVar[jColumn] && colLower[jColumn] == 0.0 && colUpper[jColumn] == 1.0) {
              double value = elements[j];
              if (value > 0.0) {
                nPlus++;
              } else {
                nMinus++;
              }
            } else {
              nPlus = 2;
              nMinus = 2;
              break;
            }
          }
          double effectiveness = 0.0;
          if (nPlus == 1 && rowUpper[i] > 0.0 && rowUpper[i] < 1.0e10) {
            // can make element smaller
            for (CoinBigIndex j = rowStart[i]; j < rowStart[i] + rowLength[i]; j++) {
              double value = elements[j];
              if (value > 0.0) {
                elements[j] -= rowUpper[i];
                // printf("pass %d row %d plus el from %g to %g\n",info->pass,
                //      i,elements[j]+rowUpper[i],elements[j]);
              }
              effectiveness += fabs(elements[j]);
            }
            rowUpper[i] = 0.0;
          } else if (nMinus == 1 && rowLower[i] < 0.0 && rowLower[i] > -1.0e10) {
            // can make element smaller in magnitude
            for (CoinBigIndex j = rowStart[i]; j < rowStart[i] + rowLength[i]; j++) {
              double value = elements[j];
              if (value < 0.0) {
                elements[j] -= rowLower[i];
                // printf("pass %d row %d minus el from %g to %g\n",info->pass,
                //      i,elements[j]+rowLower[i],elements[j]);
              }
              effectiveness += fabs(elements[j]);
            }
            rowLower[i] = 0.0;
          }
          if (effectiveness) {
            OsiRowCut rc;
            rc.setLb(rowLower[i]);
            rc.setUb(rowUpper[i]);
            CoinBigIndex start = rowStart[i];
            int n = rowLength[i];
            rc.setRow(rowLength[i], column + start, elements + start, false);
            // but get rid of tinies
            CoinPackedVector &row = rc.mutableRow();
            double *elements = row.getElements();
            int n2 = 0;
            for (n2 = 0; n2 < n; n2++) {
              if (fabs(elements[n2]) < 1.0e-12)
                break;
            }
            if (n2 < n) {
              int *columns = row.getIndices();
              for (int i = n2 + 1; i < n; i++) {
                if (fabs(elements[i]) > 1.0e-12) {
                  elements[n2] = elements[i];
                  columns[n2++] = columns[i];
                }
              }
              row.truncate(n2);
            }
            rc.setEffectiveness(effectiveness);
            assert(!info->strengthenRow[i]);
            info->strengthenRow[i] = new OsiRowCut(std::move(rc));
          }
        }
        rowLower[nKeep] = rowLower[i];
        rowUpper[nKeep] = rowUpper[i];
        nKeep++;
#endif
      }
    }
    if (nDelete) {
#ifdef OUTRUBBISH
      if (nExtraDel) {
        printf("%d rows deleted (extra %d)\n", nDelete, nExtraDel);
      }
#else
#endif
      if (info->strengthenRow) {
        // Set up pointers to real rows
        realRows = new int[nRows];
        CoinZeroN(realRows, nRows);
        for (i = 0; i < nDelete; i++)
          realRows[which[i]] = -1;
        int k = 0;
        for (i = 0; i < nRows; i++) {
          if (!realRows[i]) {
            if (i < nRealRows)
              realRows[k++] = i; // keep
            else
              realRows[k++] = -1; // objective - discard
          }
        }
      }
      rowCopy->deleteRows(nDelete, which);
      nRows = nKeep;
    }
    delete[] which;
    if (!nRows) {
#ifdef COIN_DEVELOP
      printf("All rows too long for probing\n");
#endif
      // nothing left!!
      // delete stuff
      delete rowCopy;
      if (rowCopy_) {
        delete[] rowLower;
        delete[] rowUpper;
      }
      delete[] intVar;
      // and put back unreasonable bounds on integer variables
      const double *trueLower = si.getColLower();
      const double *trueUpper = si.getColUpper();
      for (i = 0; i < nCols; i++) {
        if (intVarOriginal[i] == 2) {
          if (colUpper[i] == CGL_REASONABLE_INTEGER_BOUND)
            colUpper[i] = trueUpper[i];
          if (colLower[i] == -CGL_REASONABLE_INTEGER_BOUND)
            colLower[i] = trueLower[i];
        }
      }
      delete[] realRows;
      return 0;
    }
    // Out elements for fixed columns and sort
    elements = rowCopy->getMutableElements();
    column = rowCopy->getMutableIndices();
    rowStart = rowCopy->getMutableVectorStarts();
    rowLength = rowCopy->getMutableVectorLengths();
#if 0
    int nFixed=0;
    for (i=0;i<nCols;i++) {
      if (colUpper[i]==colLower[i])
	nFixed++;
    }
    printf("%d columns fixed\n",nFixed);
#endif
    CoinBigIndex newSize = 0;
    int *column2 = new int[nCols];
    double *elements2 = new double[nCols];
    rowStartPos = new CoinBigIndex[nRows];
    for (i = 0; i < nRows; i++) {
      double offset = 0.0;
      CoinBigIndex start = rowStart[i];
      rowStart[i] = newSize;
      CoinBigIndex save = newSize;
      CoinBigIndex end = start + rowLength[i];
      int nOther = 0;
      for (CoinBigIndex j = start; j < end; j++) {
        int iColumn = column[j];
        if (colUpper[iColumn] > colLower[iColumn]) {
          double value = elements[j];
          if (value < 0.0) {
            elements[newSize] = value;
            column[newSize++] = iColumn;
          } else if (value > 0.0) {
            elements2[nOther] = value;
            column2[nOther++] = iColumn;
          }
        } else {
          offset += colUpper[iColumn] * elements[j];
        }
      }
      rowStartPos[i] = newSize;
      for (int k = 0; k < nOther; k++) {
        elements[newSize] = elements2[k];
        column[newSize++] = column2[k];
      }
      rowLength[i] = static_cast< int >(newSize - save);
      if (offset) {
        if (rowLower[i] > -1.0e20)
          rowLower[i] -= offset;
        if (rowUpper[i] < 1.0e20)
          rowUpper[i] -= offset;
      }
    }
    delete[] column2;
    delete[] elements2;
    rowStart[nRows] = newSize;
    rowCopy->setNumElements(newSize);
  }
  CoinPackedMatrix *columnCopy = new CoinPackedMatrix(*rowCopy, 0, 0, true);
  int nRowsSafe = std::min(nRows, si.getNumRows());
#ifdef CGL_DEBUG
  const OsiRowCutDebugger *debugger = si.getRowCutDebugger();
  if (debugger && !debugger->onOptimalPath(si))
    debugger = NULL;
#else
  const OsiRowCutDebugger *debugger = NULL;
#endif

  const int *column = rowCopy->getIndices();
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
  const int *rowLength = rowCopy->getVectorLengths();
  const double *rowElements = rowCopy->getElements();
  // Arrays so user can find out what happened
  if (!lookedAt_) {
    lookedAt_ = new int[nCols];
  }
  numberThisTime_ = 0;
  // Let us never add more than twice the number of rows worth of row cuts
  // Keep cuts out of cs until end so we can find duplicates quickly
  int nRowsFake = info->inTree ? nRowsSafe / 3 : nRowsSafe;
  if (!info->inTree && !info->pass)
    nRowsFake *= 5;
  row_cut rowCut(nRowsFake, !info->inTree);
  int *markR = new int[nRows];
  double *minR = new double[nRows];
  double *maxR = new double[nRows];
  minR_ = minR;
  maxR_ = maxR;
  if (nRows > nRowsSafe) {
    // objective added
    maxR[nRowsSafe] = COIN_DBL_MAX;
    minR[nRowsSafe] = -COIN_DBL_MAX;
  }
  if (mode) {
    ninfeas = tighten(colLower, colUpper, column, rowElements,
      rowStart, rowStartPos, rowLength,
      columnCopy, rowLower, rowUpper,
      nRows, nCols, intVar,
      ((info->options & 2048) == 0) ? 2 : 3,
      primalTolerance_);
    if (!ninfeas) {
      // create column cuts where integer bounds have changed
      {
        const double *lower = si.getColLower();
        const double *upper = si.getColUpper();
        const double *colsol = si.getColSolution();
        int numberChanged = 0, ifCut = 0;
        CoinPackedVector lbs;
        CoinPackedVector ubs;
        for (i = 0; i < nCols; ++i) {
          if (intVar[i]) {
            // check integer
            assert(colUpper[i] == floor(colUpper[i] + 0.5));
            assert(colLower[i] == floor(colLower[i] + 0.5));
            // colUpper[i] = std::min(upper[i],floor(colUpper[i]+1.0e-4));
            if (colUpper[i] < upper[i] - 1.0e-8) {
              if (colUpper[i] < colsol[i] - 1.0e-8)
                ifCut = 1;
              ubs.insert(i, colUpper[i]);
              numberChanged++;
            }
            // colLower[i] = std::max(lower[i],ceil(colLower[i]-1.0e-4));
            if (colLower[i] > lower[i] + 1.0e-8) {
              if (colLower[i] > colsol[i] + 1.0e-8)
                ifCut = 1;
              lbs.insert(i, colLower[i]);
              numberChanged++;
            }
          }
        }
        if (numberChanged) {
          OsiColCut cc;
          cc.setUbs(ubs);
          cc.setLbs(lbs);
          if (ifCut) {
            cc.setEffectiveness(100.0);
          } else {
            cc.setEffectiveness(1.0e-5);
          }
#ifdef CGL_DEBUG
          checkBounds(debugger, cc);
#endif
          cs.insert(cc);
        }
      }
      if (maxProbe > 0) {
        numberThisTime_ = 0;
#if 0
        // get min max etc for rows
        tighten2(colLower, colUpper, column, rowElements,
                 rowStart, rowLength, rowLower, rowUpper,
                 minR , maxR , markR, nRows);
#else
        for (int j = 0; j < nRows; j++)
          markR[j] = -1;
#endif
        // decide what to look at
        if (mode == 1 && (info->options & 32768) == 0) {
          const double *colsol = si.getColSolution();
          double_int_pair *array = new double_int_pair[nCols];
#ifdef ZEROFAULT
          std::memset(array, 0, sizeof(double_int_pair) * nCols);
#endif
          /*
            >=100 throw away long > define/100
            1 - bias against long always but use away
            2 - bias against long always but don't use away
            4 - throw does not apply in tree
            8 - throw does not  apply at root
            16 - throw only applies if fast
            32 - bias does not apply at root
           */
          // start defined
          double multiplier = -1.0;
          if (info->inTree || (info->pass & 1) != 0)
            multiplier = 1.0;
          for (i = 0; i < nCols; i++) {
            if (intVar[i] && colUpper[i] - colLower[i] > 1.0e-8) {
              double away = fabs(0.5 - (colsol[i] - floor(colsol[i])));
              if (away < 0.49999 || !info->inTree) {
                // array[numberThisTime_].infeasibility=away;
                array[numberThisTime_].infeasibility = away * multiplier;
                // array[numberThisTime_].infeasibility=-columnLength[i];
                array[numberThisTime_++].sequence = i;
              }
            }
          }
          // end defined - start if
          std::sort(array, array + numberThisTime_, double_int_pair_compare());
          // numberThisTime_=std::min(numberThisTime_,maxProbe);
          for (i = 0; i < numberThisTime_; i++) {
            lookedAt_[i] = array[i].sequence;
          }
          delete[] array;
        } else {
          for (i = 0; i < nCols; i++) {
            if (intVar[i] && colUpper[i] - colLower[i] > 1.0e-8) {
              lookedAt_[numberThisTime_++] = i;
            }
          }
        }
        // sort to be clean
        // std::sort(lookedAt_,lookedAt_+numberThisTime_);
        if (!numberCliques_) {
          int numberColumnCutsBefore = cs.sizeColCuts();
#if FIXED_BOTH_WAYS
#define FIXED_BOTH_MAX 2000
          int nCols = si.getNumCols();
          if ((info->options & 2048) == 2048 && totalTimesCalled_ == 0
            && nCols < FIXED_BOTH_MAX) {
            delete[] endFixStart_;
            assert(!cliqueStart_);
            endFixStart_ = new int[2 * nCols];
            memset(endFixStart_, 0, 2 * nCols * sizeof(int));
          }
#endif
          /// FAST
          if ((info->options & 16384) == 0) {
            ninfeas = probe(si, debugger, cs, colLower, colUpper, rowCopy, columnCopy,
              rowStartPos, realRows, rowLower, rowUpper,
              intVar, minR, maxR, markR,
              info);
          } else {
            ninfeas = probeFast(si, debugger, cs, colLower, colUpper, rowCopy, columnCopy,
              rowStartPos, realRows, rowLower, rowUpper,
              intVar, minR, maxR, markR,
              info);
          }
#if FIXED_BOTH_WAYS
#define CC_AT_END 1
          if ((info->options & 2048) != 0 || (CC_AT_END && (info->options & 16384) != 0)) {
            int nCols = si.getNumCols();
            int numberColumnCutsAfter = cs.sizeColCuts();
            int nLo2 = 0;
            int nUp2 = 0;
            if (numberColumnCutsAfter > numberColumnCutsBefore) {
              for (int i = numberColumnCutsBefore; i < numberColumnCutsAfter; i++) {
                const OsiColCut *thisCut = cs.colCutPtr(i);
                const CoinPackedVector &lbs = thisCut->lbs();
                const CoinPackedVector &ubs = thisCut->ubs();
                nLo2 += lbs.getNumElements();
                nUp2 += ubs.getNumElements();
              }
            }
            int nLo = 0;
            int nUp = 0;
            int nFix = 0;
            int nLoC = 0;
            int nUpC = 0;
            const double *lower = si.getColLower();
            const double *upper = si.getColUpper();
            // will have redundant bound changes
            CoinPackedVector lbs;
            CoinPackedVector ubs;
            for (int iColumn = 0; iColumn < nCols; iColumn++) {
              if (colLower[iColumn] > lower[iColumn] && intVar[iColumn]) {
                lbs.insert(iColumn, colLower[iColumn]);
              }
              if (colUpper[iColumn] < upper[iColumn] && intVar[iColumn]) {
                ubs.insert(iColumn, colUpper[iColumn]);
              }
            }
            if (lbs.getNumElements()) {
              OsiColCut cc;
              cc.setLbs(lbs);
              cs.insert(cc);
            }
            if (ubs.getNumElements()) {
              OsiColCut cc;
              cc.setUbs(ubs);
              cs.insert(cc);
            }
            for (int i = 0; i < nCols; i++) {
              if (lower[i] != colLower[i]) {
                assert(lower[i] < colLower[i]);
                nLo++;
                if (!intVar[i])
                  nLoC++;
                if (colLower[i] == colUpper[i])
                  nFix++;
              }
              if (upper[i] != colUpper[i]) {
                assert(upper[i] > colUpper[i]);
                nUp++;
                if (!intVar[i])
                  nUpC++;
              }
            }
#if FIXED_BOTH_WAYS > 1
            if (nLo + nUp)
              printf("www %d lower increased (%d cont, %d cc), %d upper decreased (%d, %d) - %d fixed\n",
                nLo, nLoC, nLo2, nUp, nUpC, nUp2, nFix);
#endif
          }
#endif
        } else {
          ninfeas = probeCliques(si, debugger, cs, colLower, colUpper, rowCopy, columnCopy,
            realRows, rowLower, rowUpper,
            intVar, minR, maxR, markR,
            info);
        }
      }
    }
  } else if (maxProbe > 0) {
    // global cuts from previous calculations
    // could check more thoroughly that integers are correct
    assert(numberIntegers == numberIntegers_);
    // make up list of new variables to look at
    numberThisTime_ = 0;
    const double *colsol = si.getColSolution();
    double_int_pair *array = new double_int_pair[nCols];
#ifdef ZEROFAULT
    std::memset(array, 0, sizeof(double_int_pair) * nCols);
#endif
    for (i = 0; i < number01Integers_; i++) {
      int j = cutVector_[i].sequence;
      if (!cutVector_[i].index && colUpper[j] - colLower[j] > 1.0e-8) {
        double away = fabs(0.5 - (colsol[j] - floor(colsol[j])));
        array[numberThisTime_].infeasibility = away;
        array[numberThisTime_++].sequence = i;
      }
    }
    std::sort(array, array + numberThisTime_, double_int_pair_compare());
    numberThisTime_ = std::min(numberThisTime_, maxProbe);
    for (i = 0; i < numberThisTime_; i++) {
      lookedAt_[i] = array[i].sequence;
    }
    // sort to be clean
    // std::sort(lookedAt_,lookedAt_+numberThisTime_);
    delete[] array;
    // get min max etc for rows
    tighten2(colLower, colUpper, column, rowElements,
      rowStart, rowLength, rowLower, rowUpper,
      minR, maxR, markR, nRows);
    OsiCuts csNew;
    // don't do cuts at all if 0 (i.e. we are just checking bounds)
    if (rowCuts_) {
#if 0
      // Only look at short rows
      for (i=0;i<nRows;i++) {
        if (rowLength[i]>maxElements)
          abort(); //markR[i]=-2;
      }
#endif
      ninfeas = probeCliques(si, debugger, csNew, colLower, colUpper, rowCopy, columnCopy,
        realRows, rowLower, rowUpper,
        intVar, minR, maxR, markR,
        info);
    }
    if (!ninfeas) {
      // go through row cuts
      int nCuts = csNew.sizeRowCuts();
      int iCut;
      // need space for backward lookup
      // just for ones being looked at
      int *backward = new int[2 * nCols];
      int *onList = backward + nCols;
      for (i = 0; i < nCols; i++) {
        backward[i] = -1;
        onList[i] = 0;
      }
      for (i = 0; i < number01Integers_; i++) {
        int j = cutVector_[i].sequence;
        backward[j] = i;
        onList[j] = 1;
      }
      // first do counts
      // we know initialized to zero
      for (iCut = 0; iCut < nCuts; iCut++) {
        OsiRowCut rcut;
        CoinPackedVector rpv;
        rcut = csNew.rowCut(iCut);
        rpv = rcut.row();
        assert(rpv.getNumElements() == 2);
        const int *indices = rpv.getIndices();
        double *elements = rpv.getElements();
        double lb = rcut.lb();
        // find out which integer
        int which = 0;
        i = backward[indices[0]];
        if (i < 0 || !onList[indices[0]]) {
          which = 1;
          i = backward[indices[1]];
          // Just possible variable was general integer but now 0-1
          if (!onList[indices[which]])
            continue;
        }
        int other = indices[1 - which];
        if (lb == -COIN_DBL_MAX) {
          if (!rcut.ub()) {
            // UB
            if (elements[which] < 0.0) {
              // assert (elements[1-which]>0.0);
              //  delta to 0 => x to 0.0
              cutVector_[i].length++;
            } else {
              if (elements[1 - which] < 0.0 && fabs(elements[which] / elements[1 - which] - colUpper[other]) < 1.0e-5) {
                // delta to 1 => x to upper bound
                cutVector_[i].length++;
              } else {
                if (onList[other]) {
                  double value0 = elements[0];
                  double value1 = elements[1];
                  if (value0 * value1 == -1.0) {
                    // can do something ?
                    int j = backward[other];
                    cutVector_[i].length++;
                    cutVector_[j].length++;
                  } else {
                    continue;
                  }
                }
              }
            }
          } else {
            if (onList[other]) {
              if (elements[0] == 1.0 && elements[1] == 1.0 && rcut.ub() == 1.0) {
                // can do something ?
                int j = backward[other];
                cutVector_[i].length++;
                cutVector_[j].length++;
              } else {
                continue;
              }
            }
          }
        } else {
          assert(rcut.ub() == DBL_MAX);
          if (!lb) {
            // LB
            if (elements[which] > 0.0) {
              // assert (elements[1-which]<0.0);
              //  delta to 0 => x to 0.0
              //  flip so same as UB
              cutVector_[i].length++;
            } else {
              if (elements[1 - which] < 0.0 && fabs(elements[which] / elements[1 - which] - colUpper[other]) < 1.0e-5) {
                // delta to 1 => x to upper bound
                cutVector_[i].length++;
              } else {
                if (onList[other]) {
                  double value0 = elements[0];
                  double value1 = elements[1];
                  if (value0 * value1 == -1.0) {
                    // can do something ?
                    int j = backward[other];
                    cutVector_[i].length++;
                    cutVector_[j].length++;
                  } else {
                    continue;
                  }
                }
              }
            }
          }
        }
      }
      // allocate space
      for (i = 0; i < number01Integers_; i++) {
        int j = cutVector_[i].sequence;
        if (onList[j] && !cutVector_[i].index) {
          disaggregation thisOne = cutVector_[i];
          cutVector_[i].index = new disaggregationAction[thisOne.length];
          cutVector_[i].length = 0;
        }
      }
      // now put in
      for (iCut = 0; iCut < nCuts; iCut++) {
        OsiRowCut rcut;
        CoinPackedVector rpv;
        int iput;
        rcut = csNew.rowCut(iCut);
        rpv = rcut.row();
        assert(rpv.getNumElements() == 2);
        const int *indices = rpv.getIndices();
        double *elements = rpv.getElements();
        double lb = rcut.lb();
        // find out which integer
        // find out which integer
        int which = 0;
        i = backward[indices[0]];
        if (i < 0 || !onList[indices[0]]) {
          which = 1;
          i = backward[indices[1]];
          // Just possible variable was general integer but now 0-1
          if (!onList[indices[which]])
            continue;
        }
        int other = indices[1 - which];
        int j = other ? backward[other] : -1;
        if (lb == -COIN_DBL_MAX) {
          if (!rcut.ub()) {
            // UB
            if (elements[which] < 0.0) {
              iput = cutVector_[i].length;
              if (j >= 0)
                setAffectedInDisaggregation(cutVector_[i].index[iput], j);
              else
                setAffectedInDisaggregation(cutVector_[i].index[iput], other);
              setWhenAtUBInDisaggregation(cutVector_[i].index[iput], false);
              setAffectedToUBInDisaggregation(cutVector_[i].index[iput], false);
              setZeroOneInDisaggregation(cutVector_[i].index[iput], onList[other] != 0);
              cutVector_[i].length++;
            } else {
              if (elements[1 - which] < 0.0 && fabs(elements[which] / elements[1 - which] - colUpper[other]) < 1.0e-5) {
                // delta to 1 => x to upper bound
                iput = cutVector_[i].length;
                if (j >= 0)
                  setAffectedInDisaggregation(cutVector_[i].index[iput], j);
                else
                  setAffectedInDisaggregation(cutVector_[i].index[iput], other);
                setWhenAtUBInDisaggregation(cutVector_[i].index[iput], true);
                setAffectedToUBInDisaggregation(cutVector_[i].index[iput], true);
                setZeroOneInDisaggregation(cutVector_[i].index[iput], onList[other] != 0);
                cutVector_[i].length++;
              } else {
                if (onList[other]) {
                  double value0 = elements[0];
                  double value1 = elements[1];
                  if (value0 * value1 == -1.0) {
                    // can do something ?
                    int j = backward[other];
                    assert(j >= 0);
                    // flip so value0 1.0
                    if (value1 == 1.0) {
                      j = i;
                      i = backward[other];
                      value1 = value0;
                      value0 = 1.0;
                    }
                    assert(value0 == 1.0);
                    assert(value1 == -1.0);
                    iput = cutVector_[i].length;
                    setAffectedInDisaggregation(cutVector_[i].index[iput], j);
                    setWhenAtUBInDisaggregation(cutVector_[i].index[iput], true);
                    setAffectedToUBInDisaggregation(cutVector_[i].index[iput], true);
                    setZeroOneInDisaggregation(cutVector_[i].index[iput], true);
                    cutVector_[i].length++;
                    iput = cutVector_[j].length;
                    setAffectedInDisaggregation(cutVector_[j].index[iput], i);
                    setWhenAtUBInDisaggregation(cutVector_[j].index[iput], false);
                    setAffectedToUBInDisaggregation(cutVector_[j].index[iput], false);
                    setZeroOneInDisaggregation(cutVector_[j].index[iput], true);
                    cutVector_[j].length++;
                  }
                }
              }
            }
          } else {
            if (onList[other]) {
              if (elements[0] == 1.0 && elements[1] == 1.0 && rcut.ub() == 1.0) {
                // can do something ?
                int j = backward[other];
                assert(j >= 0);
                iput = cutVector_[i].length;
                setAffectedInDisaggregation(cutVector_[i].index[iput], j);
                setWhenAtUBInDisaggregation(cutVector_[i].index[iput], true);
                setAffectedToUBInDisaggregation(cutVector_[i].index[iput], false);
                setZeroOneInDisaggregation(cutVector_[i].index[iput], true);
                cutVector_[i].length++;
                iput = cutVector_[j].length;
                setAffectedInDisaggregation(cutVector_[j].index[iput], i);
                setWhenAtUBInDisaggregation(cutVector_[j].index[iput], true);
                setAffectedToUBInDisaggregation(cutVector_[j].index[iput], false);
                setZeroOneInDisaggregation(cutVector_[j].index[iput], true);
                cutVector_[j].length++;
              } else {
#ifdef COIN_DEVELOP
                abort();
#endif
                continue;
              }
            }
          }
        } else {
          assert(rcut.ub() == DBL_MAX);
          if (!lb) {
            // LB
            if (elements[which] > 0.0) {
              iput = cutVector_[i].length;
              if (j >= 0)
                setAffectedInDisaggregation(cutVector_[i].index[iput], j);
              else
                setAffectedInDisaggregation(cutVector_[i].index[iput], other);
              setWhenAtUBInDisaggregation(cutVector_[i].index[iput], false);
              setAffectedToUBInDisaggregation(cutVector_[i].index[iput], false);
              setZeroOneInDisaggregation(cutVector_[i].index[iput], onList[other] != 0);
              cutVector_[i].length++;
            } else {
              if (elements[1 - which] < 0.0 && fabs(elements[which] / elements[1 - which] - colUpper[other]) < 1.0e-5) {
                iput = cutVector_[i].length;
                if (j >= 0)
                  setAffectedInDisaggregation(cutVector_[i].index[iput], j);
                else
                  setAffectedInDisaggregation(cutVector_[i].index[iput], other);
                setWhenAtUBInDisaggregation(cutVector_[i].index[iput], true);
                setAffectedToUBInDisaggregation(cutVector_[i].index[iput], true);
                setZeroOneInDisaggregation(cutVector_[i].index[iput], onList[other] != 0);
                cutVector_[i].length++;
              } else {
                if (onList[other]) {
                  double value0 = elements[0];
                  double value1 = elements[1];
                  if (value0 * value1 == -1.0) {
                    // can do something ?
                    int j = backward[other];
                    assert(j >= 0);
                    // flip so value0 -1.0
                    if (value1 == -1.0) {
                      j = i;
                      i = backward[other];
                      value1 = value0;
                      value0 = -1.0;
                    }
                    assert(value0 == -1.0);
                    assert(value1 == 1.0);
                    iput = cutVector_[i].length;
                    setAffectedInDisaggregation(cutVector_[i].index[iput], j);
                    setWhenAtUBInDisaggregation(cutVector_[i].index[iput], true);
                    setAffectedToUBInDisaggregation(cutVector_[i].index[iput], true);
                    setZeroOneInDisaggregation(cutVector_[i].index[iput], true);
                    cutVector_[i].length++;
                    iput = cutVector_[j].length;
                    setAffectedInDisaggregation(cutVector_[j].index[iput], i);
                    setWhenAtUBInDisaggregation(cutVector_[j].index[iput], false);
                    setAffectedToUBInDisaggregation(cutVector_[j].index[iput], false);
                    setZeroOneInDisaggregation(cutVector_[j].index[iput], true);
                    cutVector_[j].length++;
                  }
                }
              }
            }
          }
        }
      }
      delete[] backward;
      // Now sort and get rid of duplicates
      // could also see if any are cliques
      int longest = 0;
      for (i = 0; i < number01Integers_; i++)
        longest = std::max(longest, cutVector_[i].length);
      unsigned int *sortit = new unsigned int[longest];
      for (i = 0; i < number01Integers_; i++) {
        disaggregation &thisOne = cutVector_[i];
        int k;
        int number = thisOne.length;
        for (k = 0; k < number; k++) {
          int affected = affectedInDisaggregation(thisOne.index[k]);
          int zeroOne = zeroOneInDisaggregation(thisOne.index[k]) ? 1 : 0;
          int whenAtUB = whenAtUBInDisaggregation(thisOne.index[k]) ? 1 : 0;
          int affectedToUB = affectedToUBInDisaggregation(thisOne.index[k]) ? 1 : 0;
          sortit[k] = (affected << 3) | (zeroOne << 2) | (whenAtUB << 1) | affectedToUB;
        }
        std::sort(sortit, sortit + number);
        int affectedLast = 0xffffffff;
        int zeroOneLast = 0;
        int whenAtUBLast = 0;
        int affectedToUBLast = 0;
        int put = 0;
        for (k = 0; k < number; k++) {
          int affected = sortit[k] >> 3;
          int zeroOne = (sortit[k] & 4) >> 2;
          int whenAtUB = (sortit[k] & 2) >> 1;
          int affectedToUB = sortit[k] & 1;
          disaggregationAction action;
          action.affected = 0;
          setAffectedInDisaggregation(action, affected);
          setZeroOneInDisaggregation(action, zeroOne != 0);
          setWhenAtUBInDisaggregation(action, whenAtUB != 0);
          setAffectedToUBInDisaggregation(action, affectedToUB != 0);
          if (affected != affectedLast || zeroOne != zeroOneLast) {
            // new variable
            thisOne.index[put++] = action;
          } else if (whenAtUB != whenAtUBLast || affectedToUB != affectedToUBLast) {
            // new action - what can we discover
            thisOne.index[put++] = action;
            int j = cutVector_[i].sequence;
            int k = affected;
            if (zeroOne) {
              k = cutVector_[k].sequence;
              if (logLevel_ > 1)
                printf("For %d %d 0-1 pair", j, k);
            } else {
              if (logLevel_ > 1)
                printf("For %d %d pair", j, k);
            }
            if (logLevel_ > 1)
              printf(" old whenAtUB, affectedToUB %d %d, new whenAtUB, affectedToUB %d %d\n",
                whenAtUBLast, affectedToUBLast, whenAtUB, affectedToUB);
          }
          affectedLast = affected;
          zeroOneLast = zeroOne;
          whenAtUBLast = whenAtUB;
          affectedToUBLast = affectedToUB;
        }
        if (put < number) {
          // printf("%d reduced from %d to %d\n",i,number,put);
          thisOne.length = put;
        }
      }
      // And look at all where two 0-1 variables involved
      for (i = 0; i < number01Integers_; i++) {
        disaggregation &thisOne = cutVector_[i];
        int k;
        int number = thisOne.length;
        for (k = 0; k < number; k++) {
          int affected = affectedInDisaggregation(thisOne.index[k]);
          bool zeroOne = zeroOneInDisaggregation(thisOne.index[k]);
          if (zeroOne && static_cast< int >(affected) > i) {
            bool whenAtUB = whenAtUBInDisaggregation(thisOne.index[k]);
            bool affectedToUB = affectedToUBInDisaggregation(thisOne.index[k]);
            disaggregation otherOne = cutVector_[affected];
            int numberOther = otherOne.length;
            // Could do binary search if a lot
            int lastAction = -1;
            for (int j = 0; j < numberOther; j++) {
              if (affectedInDisaggregation(otherOne.index[j]) == i) {
                bool whenAtUBOther = whenAtUBInDisaggregation(otherOne.index[j]);
                bool affectedToUBOther = affectedToUBInDisaggregation(otherOne.index[j]);
                /* action -
                   0 -> x + y <=1 (1,1 impossible)
                   1 -> x - y <=0 (1,0 impossible)
                   2 -> -x + y <=0 (0,1 impossible)
                   3 -> -x -y <= -1 (0,0 impossible)
                  10 -> x == y
                  11 -> x + y == 1
                  20 -> x == 0
                  21 -> x == 1
                  22 -> y == 0
                  23 -> y == 1
                */
                int action = -1;
                if (whenAtUB) {
                  if (affectedToUB) {
                    // x -> 1 => y -> 1
                    if (whenAtUBOther) {
                      if (affectedToUBOther) {
                        // y -> 1 => x -> 1
                        action = 10; // x,y must be same
                      } else {
                        // y -> 1 => x -> 0
                        action = 20; // If x is 1 then contradiction
                      }
                    } else {
                      if (affectedToUBOther) {
                        // y -> 0 => x -> 1
                        action = 23; // if y is 0 then contradiction
                      } else {
                        // y -> 0 => x -> 0
                        action = 1; // x,y 1,0 impossible
                      }
                    }
                  } else {
                    // x -> 1 => y -> 0
                    if (whenAtUBOther) {
                      if (affectedToUBOther) {
                        // y -> 1 => x -> 1
                        action = 22; // If y is 1 then contradiction
                      } else {
                        // y -> 1 => x -> 0
                        action = 0;
                      }
                    } else {
                      if (affectedToUBOther) {
                        // y -> 0 => x -> 1
                        action = 11; // x,y with same values impossible
                      } else {
                        // y -> 0 => x -> 0
                        action = 20; // If x is 1 then contradiction
                      }
                    }
                  }
                } else {
                  if (affectedToUB) {
                    // x -> 0 => y -> 1
                    if (whenAtUBOther) {
                      if (affectedToUBOther) {
                        // y -> 1 => x -> 1
                        action = 21; // If x is 0 then contradiction
                      } else {
                        // y -> 1 => x -> 0
                        action = 11; // x,y must be different
                      }
                    } else {
                      if (affectedToUBOther) {
                        // y -> 0 => x -> 1
                        action = 3; // one of x,y must be 1
                      } else {
                        // y -> 0 => x -> 0
                        action = 23; // if y is 0 then contradiction
                      }
                    }
                  } else {
                    // x -> 0 => y -> 0
                    if (whenAtUBOther) {
                      if (affectedToUBOther) {
                        // y -> 1 => x -> 1
                        action = 2; // x,y 0,1 impossible
                      } else {
                        // y -> 1 => x -> 0
                        action = 22; // If y is 1 then contradiction
                      }
                    } else {
                      if (affectedToUBOther) {
                        // y -> 0 => x -> 1
                        action = 21; // if x is 0 then contradiction
                      } else {
                        // y -> 0 => x -> 0
                        action = 10; // x,y must be same
                      }
                    }
                  }
                }
                assert(action >= 0);
                if (action < 4) {
                  // clique - see if there
                  if (oneFixStart_) {
                    switch (action) {
                    case 0:
                      break;
                    case 1:
                      break;
                    case 2:
                      break;
                    case 3:
                      break;
                    }
                    // If not can we add or strengthen
                  }
                  // check last action
                  if (lastAction >= 0) {
                    if (logLevel_ > 1)
                      printf("XX lastAction %d, this %d\n", lastAction, action);
                  }
                } else if (action < 12) {
                  if (logLevel_ > 1)
                    printf("XX Could eliminate one of %d %d 0-1 variables %c\n", i, affected,
                      (lastAction >= 0) ? '*' : ' ');
                  if (info->strengthenRow) {
                    OsiRowCut rc;
                    int index[2];
                    double element[2];
                    index[0] = cutVector_[i].sequence;
                    element[0] = 1.0;
                    index[1] = cutVector_[affected].sequence;
                    if (action == 10) {
                      // 10 -> x == y
                      rc.setLb(0.0);
                      rc.setUb(0.0);
                      element[1] = -1.0;
                    } else {
                      // 11 -> x + y == 1
                      rc.setLb(1.0);
                      rc.setUb(1.0);
                      element[1] = 1.0;
                    }
                    rc.setRow(2, index, element, false);
                    cs.insert(rc);
                  }
                } else {
                  if (action < 22) {
                    if (logLevel_ > 1)
                      printf("XX Could fix a 0-1 variable %d\n", i);
                  } else {
                    if (logLevel_ > 1)
                      printf("XX Could fix a 0-1 variable %d\n", affected);
                  }
                }
                // printf("%d when %d forces %d to %d , %d when %d forces %d to %d\n",
                //      i,whenAtUB,affected,affectedToUB,
                //      affected, whenAtUBOther,i, affectedToUBOther);
              }
            }
          }
        }
      }
      delete[] sortit;
    }
    if (cutVector_) {
      // now see if any disaggregation cuts are violated
      for (i = 0; i < number01Integers_; i++) {
        int j = cutVector_[i].sequence;
        double solInt = colsol[j];
        double upper, solValue;
        int icol;
        int index[2];
        double element[2];
        if (colUpper[j] - colLower[j] > 1.0e-8) {
          double away = fabs(0.5 - (solInt - floor(solInt)));
          if (away < 0.4999999) {
            disaggregation thisOne = cutVector_[i];
            int k;
            OsiRowCut rc;
            for (k = 0; k < thisOne.length; k++) {
              icol = affectedInDisaggregation(thisOne.index[k]);
              if (zeroOneInDisaggregation(thisOne.index[k]))
                icol = cutVector_[icol].sequence;
              solValue = colsol[icol];
              upper = colUpper_[icol];
              double infeasibility = 0.0;
              if (!whenAtUBInDisaggregation(thisOne.index[k])) {
                if (!affectedToUBInDisaggregation(thisOne.index[k])) {
                  // delta -> 0 => x to lb (at present just 0)
                  infeasibility = solValue - upper * solInt;
                  if (infeasibility > 1.0e-3) {
                    rc.setLb(-COIN_DBL_MAX);
                    rc.setUb(0.0);
                    index[0] = icol;
                    element[0] = 1.0;
                    index[1] = j;
                    element[1] = -upper;
                  } else {
                    infeasibility = 0.0;
                  }
                } else {
                  // delta -> 0 => x to ub
                  abort();
                }
              } else {
                if (affectedToUBInDisaggregation(thisOne.index[k])) {
                  // delta -> 1 => x to ub (?)
                  icol = affectedInDisaggregation(thisOne.index[k]);
                  if (zeroOneInDisaggregation(thisOne.index[k]))
                    icol = cutVector_[icol].sequence;
                  solValue = colsol[icol];
                  upper = colUpper_[icol];
                  if (!colLower[icol]) {
                    infeasibility = upper * solInt - solValue;
                    if (infeasibility > 1.0e-3) {
                      rc.setLb(-COIN_DBL_MAX);
                      rc.setUb(0.0);
                      index[0] = icol;
                      element[0] = -1.0;
                      index[1] = j;
                      element[1] = upper;
                    } else {
                      infeasibility = 0.0;
                    }
                  } else {
                    assert(upper == colLower[icol]);
                    infeasibility = 0.0;
                  }
                } else {
                  // delta + delta2 <= 1
                  assert(zeroOneInDisaggregation(thisOne.index[k]));
                  // delta -> 1 => delta2 -> 0
                  icol = affectedInDisaggregation(thisOne.index[k]);
                  icol = cutVector_[icol].sequence;
                  // only do if icol > j
                  if (icol > j && colUpper[icol]) {
                    solValue = colsol[icol];
                    if (!colLower[icol]) {
                      infeasibility = solInt + solValue - 1.0;
                      if (infeasibility > 1.0e-3) {
                        rc.setLb(-COIN_DBL_MAX);
                        rc.setUb(1.0);
                        index[0] = icol;
                        element[0] = 1.0;
                        index[1] = j;
                        element[1] = 1.0;
                      } else {
                        infeasibility = 0.0;
                      }
                    } else {
                      assert(upper == colLower[icol]);
                      infeasibility = 0.0;
                    }
                  }
                }
              }
              if (infeasibility) {
                rc.setEffectiveness(infeasibility);
                rc.setRow(2, index, element, false);
                if (logLevel_ > 1)
                  printf("%g <= %g * x%d + %g * x%d <= %g\n",
                    rc.lb(), element[0], index[0], element[1], index[1], rc.ub());
#ifdef CGL_DEBUG
                if (debugger)
                  assert(!debugger->invalidCut(rc));
#endif
                rowCut.addCutIfNotDuplicate(rc);
              }
            }
          }
        }
      }
    }
  }
  delete[] markR;
  delete[] minR;
  delete[] maxR;
  minR_ = NULL;
  maxR_ = NULL;
  if (markRow_) {
    delete[] markRow_;
    delete[] lookColumn_;
    markRow_ = NULL;
    markColumn_ = NULL;
    lookColumn_ = NULL;
  }
  // Add in row cuts
  if (!ninfeas) {
    rowCut.addCuts(cs, info->strengthenRow, 0);
  }
  // delete stuff
  delete rowCopy;
  delete columnCopy;
  if (rowCopy_) {
    delete[] rowLower;
    delete[] rowUpper;
  }
  delete[] intVar;
  delete[] rowStartPos;
  delete[] realRows;
  // and put back unreasonable bounds on integer variables
  const double *trueLower = si.getColLower();
  const double *trueUpper = si.getColUpper();
  if (!ninfeas) {
    for (i = 0; i < nCols; i++) {
      if (intVarOriginal[i] == 2) {
        if (colUpper[i] == CGL_REASONABLE_INTEGER_BOUND)
          colUpper[i] = trueUpper[i];
        if (colLower[i] == -CGL_REASONABLE_INTEGER_BOUND)
          colLower[i] = trueLower[i];
      }
    }
  } else {
    memcpy(colLower, trueLower, nCols * sizeof(double));
    memcpy(colUpper, trueUpper, nCols * sizeof(double));
  }
  if (!info->inTree && ((info->options & 4) == 4 || ((info->options & 8) && !info->pass))) {
    int numberRowCutsAfter = cs.sizeRowCuts();
    for (int i = numberRowCutsBefore; i < numberRowCutsAfter; i++)
      cs.rowCutPtr(i)->setGloballyValid();
  }
  return ninfeas;
}
// Choose tighter bounds
static void cleanBounds(int way, int &nFix, int *index, double *element)
{
  CoinSort_2(index, index + nFix, element);
  int last = -1;
  int n = 0;
  for (int i = 0; i < nFix; i++) {
    if (index[i] > last) {
      index[n] = index[i];
      element[n++] = element[i];
      last = index[i];
    } else if (way > 0) {
      // tighter ub
      if (element[n - 1] > element[i])
        element[n - 1] = element[i];
    } else {
      // tighter lb
      if (element[n - 1] < element[i])
        element[n - 1] = element[i];
    }
  }
  nFix = n;
}
#if CHECK_MINR_MAXR > 0
static void writeMps(int j, int iway, const OsiSolverInterface *solver,
  double *colLower, double *colUpper,
  const CoinPackedMatrix *rowCopy,
  const CoinPackedMatrix *columnCopy,
  const double *rowLower, const double *rowUpper,
  const char *intVar)
{
  char name[30];
  sprintf(name, "/tmp/probeJ%dWay%d.mps", j, iway);
  CoinMpsIO writer;
  // writer.passInMessageHandler(handler_);
  //*writer.messagesPointer() = coinMessages();
  writer.setMpsData(*columnCopy, COIN_DBL_MAX,
    colLower, colUpper,
    solver->getObjCoefficients(),
    intVar, //(const char *)NULL /*integrality*/,
    rowLower, rowUpper, NULL, NULL);
  double objectiveOffset;
  solver->getDblParam(OsiObjOffset, objectiveOffset);
  writer.setObjectiveOffset(objectiveOffset);
  // set name
  writer.setProblemName(name);
  writer.writeMps(name);
  int numberRows = rowCopy->getNumRows();
  int numberColumns = rowCopy->getNumCols();
  // get matrix data pointers
  const int *column = rowCopy->getIndices();
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
  const int *rowLength = rowCopy->getVectorLengths();
  const double *element = rowCopy->getElements();
  int numberChanged = 1, iPass = 0;
  int numberInfeasible = 0;
  int totalTightened = 0;

  double tolerance = 1.0e-8;

  // Save column bounds CHANGE if works
  double *saveLower = new double[numberColumns];
  CoinMemcpyN(colLower, numberColumns, saveLower);
  double *saveUpper = new double[numberColumns];
  CoinMemcpyN(colUpper, numberColumns, saveUpper);

  int iRow, iColumn;
  const double large = 1.0e30;
  const double large2 = 1.0e40;
#define MAXPASS 1

  // Loop round seeing if we can tighten bounds
  // Would be faster to have a stack of possible rows
  // and we put altered rows back on stack
  int numberCheck = -1;
  while (numberChanged > numberCheck) {

    numberChanged = 0; // Bounds tightened this pass

    if (iPass == MAXPASS)
      break;
    iPass++;

    for (iRow = 0; iRow < numberRows; iRow++) {

      if (rowLower[iRow] > -large || rowUpper[iRow] < large) {

        // possible row
        int infiniteUpper = 0;
        int infiniteLower = 0;
        double maximumUp = 0.0;
        double maximumDown = 0.0;
        double newBound;
        CoinBigIndex rStart = rowStart[iRow];
        CoinBigIndex rEnd = rowStart[iRow] + rowLength[iRow];
        CoinBigIndex j;
        // Compute possible lower and upper ranges

        for (j = rStart; j < rEnd; ++j) {
          double value = element[j];
          iColumn = column[j];
          if (value > 0.0) {
            if (colUpper[iColumn] >= large) {
              ++infiniteUpper;
            } else {
              maximumUp += colUpper[iColumn] * value;
            }
            if (colLower[iColumn] <= -large) {
              ++infiniteLower;
            } else {
              maximumDown += colLower[iColumn] * value;
            }
          } else if (value < 0.0) {
            if (colUpper[iColumn] >= large) {
              ++infiniteLower;
            } else {
              maximumDown += colUpper[iColumn] * value;
            }
            if (colLower[iColumn] <= -large) {
              ++infiniteUpper;
            } else {
              maximumUp += colLower[iColumn] * value;
            }
          }
        }
        // Build in a margin of error
        maximumUp += 1.0e-8 * fabs(maximumUp);
        maximumDown -= 1.0e-8 * fabs(maximumDown);
        double maxUp = maximumUp + infiniteUpper * 1.0e31;
        double maxDown = maximumDown - infiniteLower * 1.0e31;
        if (maxUp <= rowUpper[iRow] + tolerance && maxDown >= rowLower[iRow] - tolerance) {

          // Row is redundant - make totally free
          // NO - can't do this for postsolve
          // rowLower[iRow]=-COIN_DBLMAX;
          // rowUpper[iRow]=COIN_DBL_MAX;
          // printf("Redundant row in presolveX %d\n",iRow);

        } else {
          if (maxUp < rowLower[iRow] - 100.0 * tolerance || maxDown > rowUpper[iRow] + 100.0 * tolerance) {
            // problem is infeasible - exit at once
            numberInfeasible++;
            break;
          }
          double lower = rowLower[iRow];
          double upper = rowUpper[iRow];
          for (j = rStart; j < rEnd; ++j) {
            double value = element[j];
            iColumn = column[j];
            double nowLower = colLower[iColumn];
            double nowUpper = colUpper[iColumn];
            if (value > 0.0) {
              // positive value
              if (lower > -large) {
                if (!infiniteUpper) {
                  assert(nowUpper < large2);
                  newBound = nowUpper + (lower - maximumUp) / value;
                  // relax if original was large
                  if (fabs(maximumUp) > 1.0e8)
                    newBound -= 1.0e-12 * fabs(maximumUp);
                } else if (infiniteUpper == 1 && nowUpper > large) {
                  newBound = (lower - maximumUp) / value;
                  // relax if original was large
                  if (fabs(maximumUp) > 1.0e8)
                    newBound -= 1.0e-12 * fabs(maximumUp);
                } else {
                  newBound = -COIN_DBL_MAX;
                }
                if (newBound > nowLower + 1.0e-12 && newBound > -large) {
                  // Tighten the lower bound
                  numberChanged++;
                  // check infeasible (relaxed)
                  if (nowUpper < newBound) {
                    if (nowUpper - newBound < -100.0 * tolerance)
                      numberInfeasible++;
                    else
                      newBound = nowUpper;
                  }
                  colLower[iColumn] = newBound;
                  // adjust
                  double now;
                  if (nowLower < -large) {
                    now = 0.0;
                    infiniteLower--;
                  } else {
                    now = nowLower;
                  }
                  maximumDown += (newBound - now) * value;
                  nowLower = newBound;
                }
              }
              if (upper < large) {
                if (!infiniteLower) {
                  assert(nowLower > -large2);
                  newBound = nowLower + (upper - maximumDown) / value;
                  // relax if original was large
                  if (fabs(maximumDown) > 1.0e8)
                    newBound += 1.0e-12 * fabs(maximumDown);
                } else if (infiniteLower == 1 && nowLower < -large) {
                  newBound = (upper - maximumDown) / value;
                  // relax if original was large
                  if (fabs(maximumDown) > 1.0e8)
                    newBound += 1.0e-12 * fabs(maximumDown);
                } else {
                  newBound = COIN_DBL_MAX;
                }
                if (newBound < nowUpper - 1.0e-12 && newBound < large) {
                  // Tighten the upper bound
                  numberChanged++;
                  // check infeasible (relaxed)
                  if (nowLower > newBound) {
                    if (newBound - nowLower < -100.0 * tolerance)
                      numberInfeasible++;
                    else
                      newBound = nowLower;
                  }
                  colUpper[iColumn] = newBound;
                  // adjust
                  double now;
                  if (nowUpper > large) {
                    now = 0.0;
                    infiniteUpper--;
                  } else {
                    now = nowUpper;
                  }
                  maximumUp += (newBound - now) * value;
                  nowUpper = newBound;
                }
              }
            } else {
              // negative value
              if (lower > -large) {
                if (!infiniteUpper) {
                  assert(nowLower < large2);
                  newBound = nowLower + (lower - maximumUp) / value;
                  // relax if original was large
                  if (fabs(maximumUp) > 1.0e8)
                    newBound += 1.0e-12 * fabs(maximumUp);
                } else if (infiniteUpper == 1 && nowLower < -large) {
                  newBound = (lower - maximumUp) / value;
                  // relax if original was large
                  if (fabs(maximumUp) > 1.0e8)
                    newBound += 1.0e-12 * fabs(maximumUp);
                } else {
                  newBound = COIN_DBL_MAX;
                }
                if (newBound < nowUpper - 1.0e-12 && newBound < large) {
                  // Tighten the upper bound
                  numberChanged++;
                  // check infeasible (relaxed)
                  if (nowLower > newBound) {
                    if (newBound - nowLower < -100.0 * tolerance)
                      numberInfeasible++;
                    else
                      newBound = nowLower;
                  }
                  colUpper[iColumn] = newBound;
                  // adjust
                  double now;
                  if (nowUpper > large) {
                    now = 0.0;
                    infiniteLower--;
                  } else {
                    now = nowUpper;
                  }
                  maximumDown += (newBound - now) * value;
                  nowUpper = newBound;
                }
              }
              if (upper < large) {
                if (!infiniteLower) {
                  assert(nowUpper < large2);
                  newBound = nowUpper + (upper - maximumDown) / value;
                  // relax if original was large
                  if (fabs(maximumDown) > 1.0e8)
                    newBound -= 1.0e-12 * fabs(maximumDown);
                } else if (infiniteLower == 1 && nowUpper > large) {
                  newBound = (upper - maximumDown) / value;
                  // relax if original was large
                  if (fabs(maximumDown) > 1.0e8)
                    newBound -= 1.0e-12 * fabs(maximumDown);
                } else {
                  newBound = -COIN_DBL_MAX;
                }
                if (newBound > nowLower + 1.0e-12 && newBound > -large) {
                  // Tighten the lower bound
                  numberChanged++;
                  // check infeasible (relaxed)
                  if (nowUpper < newBound) {
                    if (nowUpper - newBound < -100.0 * tolerance)
                      numberInfeasible++;
                    else
                      newBound = nowUpper;
                  }
                  colLower[iColumn] = newBound;
                  // adjust
                  double now;
                  if (nowLower < -large) {
                    now = 0.0;
                    infiniteUpper--;
                  } else {
                    now = nowLower;
                  }
                  maximumUp += (newBound - now) * value;
                  nowLower = newBound;
                }
              }
            }
          }
        }
      }
    }
    totalTightened += numberChanged;
    if (iPass == 1)
      numberCheck = numberChanged >> 4;
    if (numberInfeasible)
      break;
  }
  CoinMemcpyN(saveLower, numberColumns, colLower);
  CoinMemcpyN(saveUpper, numberColumns, colUpper);
  delete[] saveLower;
  delete[] saveUpper;
  if (numberInfeasible) {
    printf("Can change bound on %d\n", j);
    // return 1;
  }
}
#endif
// Does probing and adding cuts
int CglProbing::probe(const OsiSolverInterface &si,
  const OsiRowCutDebugger *
#ifdef CGL_DEBUG
    debugger
#endif
  ,
  OsiCuts &cs,
  double *COIN_RESTRICT colLower, double *COIN_RESTRICT colUpper,
  CoinPackedMatrix *rowCopy,
  CoinPackedMatrix *columnCopy,
  const CoinBigIndex *COIN_RESTRICT rowStartPos, const int *COIN_RESTRICT realRows,
  const double *COIN_RESTRICT rowLower, const double *COIN_RESTRICT rowUpper,
  const char *COIN_RESTRICT intVar, double *COIN_RESTRICT minR, double *COIN_RESTRICT maxR,
  int *COIN_RESTRICT markR,
  CglTreeInfo *info)
{
  int nRows = rowCopy->getNumRows();
  int nRowsSafe = std::min(nRows, si.getNumRows());
  int nCols = rowCopy->getNumCols();
  const double *COIN_RESTRICT currentColLower = si.getColLower();
  const double *COIN_RESTRICT currentColUpper = si.getColUpper();
  // Set up maxes
  int maxStack = info->inTree ? maxStack_ : maxStackRoot_;
  int maxPass = info->inTree ? maxPass_ : maxPassRoot_;
  if ((totalTimesCalled_ % 10) == -1) {
    int newMax = std::min(2 * maxStack, 50);
    maxStack = std::max(newMax, maxStack);
  }
  maxStack = std::min(maxStack, 90000);
  // need a way for user to ask for more
  maxStack = 80;
  // if ((info->options&2048)!=0&&!info->pass)
  // maxStack=200;
  double relaxedTolerance = 2.0 * primalTolerance_;
#define ONE_ARRAY
#ifdef ONE_ARRAY
  unsigned int DIratio = sizeof(double) / sizeof(int);
  assert(DIratio == 1 || DIratio == 2);
  int nSpace = 10 * nCols + 4 * nRows + 2 * maxStack;
  nSpace += (6 * nCols + nRows + DIratio - 1) >> (DIratio - 1);
  // #define UPDATE_MINR_MAXR
#ifdef UPDATE_MINR_MAXR
  // If works try and use existing arrays
  nSpace += 2 * maxStack;
  nSpace += maxStack >> (DIratio - 1);
#endif
#if FIXED_BOTH_WAYS
  nSpace += 2 * nCols;
  nSpace += nCols >> (DIratio - 1);
#endif
  double *COIN_RESTRICT colsol = new double[nSpace];
  double *COIN_RESTRICT djs = colsol + nCols;
  double *COIN_RESTRICT columnGap = djs + nCols;
  double *COIN_RESTRICT saveL = columnGap + nCols;
  double *COIN_RESTRICT saveU = saveL + 2 * nCols;
  double *COIN_RESTRICT saveMin = saveU + 2 * nCols;
  double *COIN_RESTRICT saveMax = saveMin + nRows;
  double *COIN_RESTRICT largestPositiveInRow = saveMax + nRows;
  double *COIN_RESTRICT largestNegativeInRow = largestPositiveInRow + nRows;
  double *COIN_RESTRICT element = largestNegativeInRow + nRows;
  double *COIN_RESTRICT lo0 = element + nCols;
  double *COIN_RESTRICT up0 = lo0 + maxStack;
  double *COIN_RESTRICT newLo = up0 + maxStack;
  double *COIN_RESTRICT newUp = newLo + nCols;
#if FIXED_BOTH_WAYS
  bool doFixedBothWays = (totalTimesCalled_ == 0 && endFixStart_
    && (info->options & 2048) != 0);
  if (doFixedBothWays) {
    for (int i = 0; i < nCols; i++) {
      newLo[i] = COIN_DBL_MAX;
      newUp[i] = -COIN_DBL_MAX;
    }
  }
  int nFixedBoth = 0;
#if FIXED_BOTH_WAYS > 2
  int nToBoundZ = 0;
  int nCutZ = 0;
#endif
#endif
#ifndef UPDATE_MINR_MAXR
  int *COIN_RESTRICT markC = reinterpret_cast< int * COIN_RESTRICT >(newUp + nCols);
#else
  // If works try and use existing arrays
  // If + then increase in lower, negative decrease in upper
  double *COIN_RESTRICT realChange = newUp + nCols;
  int *COIN_RESTRICT markC = reinterpret_cast< int * COIN_RESTRICT >(realChange + nCols);
#endif
  int *COIN_RESTRICT stackC = markC + nCols;
  int *COIN_RESTRICT stackR = stackC + 2 * nCols;
  int *COIN_RESTRICT index = stackR + nRows;
  int *COIN_RESTRICT stackC0 = index + nCols;
#if FIXED_BOTH_WAYS
  int *COIN_RESTRICT fixedBoth = stackC0 + nCols;
  // Really must clean up arrays
  // Make more efficient when two versions of tightenPrimalBounds
  double *COIN_RESTRICT saveFColLower = new double[2 * nCols + 2 * nRows];
  double *COIN_RESTRICT saveFColUpper = saveFColLower + nCols;
  double *COIN_RESTRICT saveFRowLower = saveFColUpper + nCols;
  double *COIN_RESTRICT saveFRowUpper = saveFRowLower + nRows;
  if (!doFixedBothWays)
    fixedBoth = NULL;
#elif MANY_TIGHTEN
  // Really must clean up arrays
  // Make more efficient when two versions of tightenPrimalBounds
  double *COIN_RESTRICT saveFColLower = NULL;
  double *COIN_RESTRICT saveFColUpper = NULL;
  double *COIN_RESTRICT saveFRowLower = NULL;
  double *COIN_RESTRICT saveFRowUpper = NULL;
  if ((info->options & 2048) != 0) {
    saveFColLower = new double[2 * nCols + 2 * nRows];
    saveFColUpper = saveFColLower + nCols;
    saveFRowLower = saveFColUpper + nCols;
    saveFRowUpper = saveFRowLower + nRows;
  }
#endif
#ifdef UPDATE_MINR_MAXR
  int *COIN_RESTRICT whichRealChange = stackC0 + maxStack;
#if FIXED_BOTH_WAYS
  if (doFixedBothWays)
    whichRealChange += 2 * nCols;
#endif
  assert(whichRealChange + 2 * maxStack < ((int *)(colsol + nSpace)));
  ;
#endif
#else
  double *COIN_RESTRICT colsol = new double[nCols];
  double *COIN_RESTRICT djs = new double[nCols];
  double *COIN_RESTRICT columnGap = new double[nCols];
  double *COIN_RESTRICT saveL = new double[2 * nCols];
  double *COIN_RESTRICT saveU = new double[2 * nCols];
  double *COIN_RESTRICT saveMin = new double[nRows];
  double *COIN_RESTRICT saveMax = new double[nRows];
  double *COIN_RESTRICT largestPositiveInRow = new double[nRows];
  double *COIN_RESTRICT largestNegativeInRow = new double[nRows];
  double *COIN_RESTRICT element = new double[nCols];
  double *COIN_RESTRICT lo0 = new double[maxStack];
  double *COIN_RESTRICT up0 = new double[maxStack];
  int *COIN_RESTRICT markC = new int[nCols];
  int *COIN_RESTRICT stackC = new int[2 * nCols];
  int *COIN_RESTRICT stackR = new int[nRows];
  int *COIN_RESTRICT index = new int[nCols];
  int *COIN_RESTRICT stackC0 = new int[maxStack];
#endif
#ifdef UPDATE_MINR_MAXR
  int numberRealChanges = 0;
#endif
  int *COIN_RESTRICT changedDown = endFixStart_;
  int *COIN_RESTRICT changedUp = changedDown ? changedDown + nCols : NULL;
  // Let us never add more than twice the number of rows worth of row cuts
  // Keep cuts out of cs until end so we can find duplicates quickly
#define PROBING4
#ifdef PROBING4
  int nRowsFake = info->inTree ? nRowsSafe / 3 : nRowsSafe * 10;
#else
  int nRowsFake = info->inTree ? nRowsSafe / 3 : nRowsSafe;
#endif
  if (!info->inTree && !info->pass)
    nRowsFake *= 10;
  bool justReplace = ((info->options & 64) != 0) && (realRows != NULL) || ((info->options & 16384) != 0);
  if (justReplace) {
    nRowsFake = nRows;
  }
  row_cut rowCut(nRowsFake, !info->inTree);
  row_cut *rowCutFake = NULL;
  if ((mode_ & 8) != 0) {
    // redo to be more compact
    rowCutFake = new row_cut(10000, !info->inTree);
  }
  totalTimesCalled_++;
  const int *COIN_RESTRICT column = rowCopy->getIndices();
  const CoinBigIndex *COIN_RESTRICT rowStart = rowCopy->getVectorStarts();
  const double *COIN_RESTRICT rowElements = rowCopy->getElements();
  const int *COIN_RESTRICT row = columnCopy->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = columnCopy->getVectorStarts();
  const int *COIN_RESTRICT columnLength = columnCopy->getVectorLengths();
  const double *COIN_RESTRICT columnElements = columnCopy->getElements();
#define MOVE_SINGLETONS
#ifdef MOVE_SINGLETONS
  const double *COIN_RESTRICT objective = si.getObjCoefficients();
  const int *COIN_RESTRICT columnLength2 = si.getMatrixByCol()->getVectorLengths();
  // see whether to do
  int moveSingletonRows = nRows;
  if (nRows > si.getNumRows()) {
    // skip objective constraint
    moveSingletonRows = nRowsSafe;
  } else if (nRows < si.getNumRows()) {
    // some rows missing
    moveSingletonRows = 0;
  } else {
    // check objective constraint and one row dropped
    double direction = si.getObjSense();
    double cutoff;
    bool cutoff_available = si.getDblParam(OsiDualObjectiveLimit, cutoff);
    if (!cutoff_available || usingObjective_ < 1)
      cutoff = si.getInfinity();
    cutoff *= direction;
    if (cutoff < 1.0e30) {
      // switch off
      moveSingletonRows = 0;
    }
  }
#endif
  bool anyColumnCuts = false;
  int ninfeas = 0;
  int rowCuts;
  double disaggEffectiveness;
  /* clean up djs and solution */
  CoinMemcpyN(si.getReducedCost(), nCols, djs);
  CoinMemcpyN(si.getColSolution(), nCols, colsol);
  disaggEffectiveness = 1.0e-3;
  rowCuts = rowCuts_;
  // CoinBigIndex *COIN_RESTRICT rowStartPos = new CoinBigIndex [nRows];
  const int *COIN_RESTRICT rowLength = rowCopy->getVectorLengths();
  for (int i = 0; i < nRows; i++) {
    assert(rowStart[i] + rowLength[i] == rowStart[i + 1]);
    CoinBigIndex kk;
#ifndef NDEBUG
    for (kk = rowStart[i]; kk < rowStart[i + 1]; kk++) {
      double value = rowElements[kk];
      if (value > 0.0)
        break;
    }
    assert(rowStartPos[i] == kk);
#endif
    double value;
    value = 0.0;
    for (kk = rowStart[i]; kk < rowStartPos[i]; kk++) {
      int iColumn = column[kk];
      double gap = std::min(1.0e100, colUpper[iColumn] - colLower[iColumn]);
      value = std::min(value, gap * rowElements[kk]);
    }
    largestNegativeInRow[i] = value;
    value = 0.0;
    for (; kk < rowStart[i + 1]; kk++) {
      int iColumn = column[kk];
      double gap = std::min(1.0e100, colUpper[iColumn] - colLower[iColumn]);
      value = std::max(value, gap * rowElements[kk]);
    }
    largestPositiveInRow[i] = value;
  }
  double direction = si.getObjSense();
  for (int i = 0; i < nCols; ++i) {
    double djValue = djs[i] * direction;
    double gap = colUpper[i] - colLower[i];
    if (gap > 1.0e-8) {
      if (colsol[i] < colLower[i] + primalTolerance_) {
        colsol[i] = colLower[i];
        djs[i] = std::max(0.0, djValue);
      } else if (colsol[i] > colUpper[i] - primalTolerance_) {
        colsol[i] = colUpper[i];
        djs[i] = std::min(0.0, djValue);
      } else {
        djs[i] = 0.0;
      }
    }
    columnGap[i] = gap - primalTolerance_;
  }

  int ipass = 0, nfixed = -1;

  double cutoff;
  bool cutoff_available = si.getDblParam(OsiDualObjectiveLimit, cutoff);
  if (!cutoff_available || usingObjective_ < 0) { // cut off isn't set or isn't valid
    cutoff = si.getInfinity();
  }
  cutoff *= direction;
  if (fabs(cutoff) > 1.0e30)
    assert(cutoff > 1.0e30);
  double current = si.getObjValue();
  current *= direction;
  /* for both way coding */
  int nstackC0 = -1;
  int nstackR, nstackC;
  // int nFix=0;
  for (int i = 0; i < nCols; i++) {
    if (colUpper[i] - colLower[i] < 1.0e-8) {
      markC[i] = 3;
      // nFix++;
    } else {
      markC[i] = 0;
      if (colUpper[i] > 1.0e10)
        markC[i] |= 8;
      if (colLower[i] < -1.0e10)
        markC[i] |= 4;
    }
  }
  // printf("PROBE %d fixed out of %d\n",nFix,nCols);
  double tolerance = 1.0e1 * primalTolerance_;
  // If we are going to replace coefficient then we don't need to be effective
  // double needEffectiveness = info->strengthenRow ? -1.0e10 : 1.0e-3;
  double needEffectiveness = info->strengthenRow ? 1.0e-8 : 1.0e-3;
  if (justReplace && (info->pass & 1) != 0)
    needEffectiveness = -1.0e10;
  if (PROBING_EXTRA_STUFF) {
    int nCut = 0;
    for (int iRow = 0; iRow < nRows; iRow++) {
      int numberInt = 0;
      int whichInt = -1;
      int numberNeg = 0;
      double sumFixed = 0.0;
      double intValue = 0.0;
      for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow + 1]; j++) {
        int jColumn = column[j];
        double value = rowElements[j];
        if (colUpper[jColumn] > colLower[jColumn] + 1.0e-8) {
          if (intVar[jColumn]) {
            numberInt++;
            whichInt = jColumn;
            intValue = value;
          } else if (value < 0) {
            numberNeg++;
          }
        } else {
          sumFixed += colLower[jColumn] * value;
        }
      }
      if (numberInt == 1 && numberNeg == 0 && intValue < 0.0 && !rowUpper[iRow] && rowLower[iRow] < -1.0e30 && !sumFixed) {
        double intSol = colsol[whichInt];
        for (CoinBigIndex j = rowStart[iRow]; j < rowStart[iRow + 1]; j++) {
          int jColumn = column[j];
          // double value = rowElements[j];
          if (colUpper[jColumn] > colLower[jColumn] + 1.0e-8) {
            if (!intVar[jColumn]) {
              if (colLower[jColumn] || colUpper[jColumn] > 1.0)
                continue;
              ;
              double upper = colUpper[jColumn];
              if (colsol[jColumn] > intSol * upper + 1.0e-4) {
                nCut++;
                {
                  OsiRowCut *rc = new OsiRowCut;
                  rc->setLb(-COIN_DBL_MAX);
                  rc->setUb(0.0);
                  rc->setEffectiveness(1.0e-5);
                  int index[2];
                  double element[2];
                  index[0] = jColumn;
                  index[1] = whichInt;
                  element[0] = 1.0;
                  element[1] = -upper;
                  rc->setRow(2, index, element, false);
                  cs.insert(rc);
                }
              }
            }
          }
        }
      }
    }
    if (nCut)
      printf("%d possible cuts\n", nCut);
  }
  bool saveFixingInfo = false;
#if PROBING100
  CglTreeProbingInfo *probingInfo = dynamic_cast< CglTreeProbingInfo * >(info);
  const int *COIN_RESTRICT backward = NULL;
  const int *COIN_RESTRICT integerVariable = NULL;
  const int *COIN_RESTRICT toZero = NULL;
  const int *COIN_RESTRICT toOne = NULL;
  const fixEntry *fixEntries = NULL;
#endif
  if (info->inTree) {
#if PROBING100
    backward = probingInfo->backward();
    integerVariable = probingInfo->integerVariable();
    toZero = probingInfo->toZero();
    toOne = probingInfo->toOne();
    fixEntries = probingInfo->fixEntries();
#endif
  } else {
    saveFixingInfo = (info->initializeFixing(&si) > 0);
  }
  const double probingDeadline = (maxSeconds_ > 0.0) ? (CoinWallclockTime() + maxSeconds_) : 0.0;
  while (ipass < maxPass && nfixed) {
    if (probingDeadline && CoinWallclockTime() > probingDeadline)
      break;
    int iLook;
    ipass++;
    // printf("pass %d\n",ipass);
    nfixed = 0;
    int justFix = (!info->inTree && !info->pass) ? -1 : 0;
    int maxProbe = info->inTree ? maxProbe_ : maxProbeRoot_;
    if (justFix < 0)
      maxProbe = numberThisTime_;
    if (maxProbe == 123) {
      // Try and be a bit intelligent
      maxProbe = 0;
      if (!info->inTree) {
        if (!info->pass || numberThisTime_ < -100) {
          maxProbe = numberThisTime_;
        } else {
          int cutDown = 4;
          int offset = info->pass % cutDown;
          int i;
          int k = 0;
          int kk = offset;
          for (i = 0; i < numberThisTime_; i++) {
            if (!kk) {
#define XXXXXX
#ifdef XXXXXX
              lookedAt_[maxProbe] = lookedAt_[i];
#endif
              maxProbe++;
              kk = cutDown - 1;
            } else {
              stackC[k++] = lookedAt_[i];
              kk--;
            }
          }
#ifdef XXXXXX
          memcpy(lookedAt_ + maxProbe, stackC, k * sizeof(int));
#endif
        }
      } else {
        // in tree
        if (numberThisTime_ < 200) {
          maxProbe = numberThisTime_;
        } else {
          int cutDown = std::max(numberThisTime_ / 100, 4);
          int offset = info->pass % cutDown;
          int i;
          int k = 0;
          int kk = offset;
          for (i = 0; i < numberThisTime_; i++) {
            if (!kk) {
#ifdef XXXXXX
              lookedAt_[maxProbe] = lookedAt_[i];
#endif
              maxProbe++;
              kk = cutDown - 1;
            } else {
              stackC[k++] = lookedAt_[i];
              kk--;
            }
          }
#ifdef XXXXXX
          memcpy(lookedAt_ + maxProbe, stackC, k * sizeof(int));
#endif
        }
      }
    }
    // printf("maxStack %d maxPass %d numberThisTime %d info pass %d\n",
    //    maxStack,maxPass,numberThisTime_,info->pass);
    for (iLook = 0; iLook < numberThisTime_; iLook++) {
      if (probingDeadline && (iLook & 127) == 0 && CoinWallclockTime() > probingDeadline) {
        ipass = maxPass;
        break;
      }
      double solval;
      double down;
      double up;
      int j = lookedAt_[iLook];
#ifdef UPDATE_MINR_MAXR
      // update
      bool feasible2 = true;
      for (int i = 0; i < numberRealChanges; i++) {
        int jcol = whichRealChange[i];
        double movement = realChange[i];
        CoinBigIndex k;
        for (k = columnStart[jcol];
          k < columnStart[jcol] + columnLength[jcol]; k++) {
          int irow = row[k];
          double value = columnElements[k];
          if (movement > 0.0) {
            /* up */
            if (value > 0.0) {
              /* up does not change - down does */
              if (minR[irow] > -1.0e10)
                minR[irow] += value * movement;
              if (minR[irow] > rowUpper[irow] + 1.0e-5) {
                feasible2 = false;
                break;
              }
            } else {
              /* down does not change - up does */
              if (maxR[irow] < 1.0e10)
                maxR[irow] += value * movement;
              if (maxR[irow] < rowLower[irow] - 1.0e-5) {
                feasible2 = false;
                break;
              }
            }
          } else {
            /* down */
            if (value < 0.0) {
              /* up does not change - down does */
              if (minR[irow] > -1.0e10)
                minR[irow] += value * movement;
              if (minR[irow] > rowUpper[irow] + 1.0e-5) {
                feasible2 = false;
                break;
              }
            } else {
              /* down does not change - up does */
              if (maxR[irow] < 1.0e10)
                maxR[irow] += value * movement;
              if (maxR[irow] < rowLower[irow] - 1.0e-5) {
                feasible2 = false;
                break;
              }
            }
          }
        }
      }
      if (!feasible2) {
        printf("NOT feasible\n");
      }
#if CHECK_MINR_MAXR
      // check
      {
        int i, j;
        CoinBigIndex k;
        CoinBigIndex krs, kre;
        int iflagu, iflagl;
        double dmaxup, dmaxdown;

        for (i = 0; i < nRows; ++i) {
          if (markR[i] != -1)
            printf("why for %d markR is %d (not -1)\n", i, markR[i]);
          if (rowLower[i] > -1.0e20 || rowUpper[i] < 1.0e20) {
            iflagu = 0;
            iflagl = 0;
            dmaxup = 0.0;
            dmaxdown = 0.0;
            krs = rowStart[i];
            kre = rowStart[i + 1]; //+rowLength[i];

            /* ------------------------------------------------------------*/
            /* Compute L(i) and U(i) */
            /* ------------------------------------------------------------*/
            for (k = krs; k < kre; ++k) {
              double value = rowElements[k];
              j = column[k];
              if (value > 0.0) {
                if (colUpper[j] < 1.0e12)
                  dmaxup += colUpper[j] * value;
                else
                  ++iflagu;
                if (colLower[j] > -1.0e12)
                  dmaxdown += colLower[j] * value;
                else
                  ++iflagl;
              } else if (value < 0.0) {
                if (colUpper[j] < 1.0e12)
                  dmaxdown += colUpper[j] * value;
                else
                  ++iflagl;
                if (colLower[j] > -1.0e12)
                  dmaxup += colLower[j] * value;
                else
                  ++iflagu;
              }
            }
            if (iflagu) {
              if (maxR[i] < 1.0e40)
                printf("Why for %d maxRi is %g not 1.0e60\n",
                  i, maxR[i]);
            } else {
              if (fabs(maxR[i] - dmaxup) > 1.0e-5)
                printf("Why for %d maxRi is %g not %g\n",
                  i, maxR[i], dmaxup);
            }
            if (iflagl) {
              if (minR[i] > -1.0e40)
                printf("Why for %d minRi is %g not -1.0e60\n",
                  i, minR[i]);
            } else {
              minR[i] = dmaxdown;
              if (fabs(minR[i] - dmaxdown) > 1.0e-5)
                printf("Why for %d minRi is %g not %g\n",
                  i, minR[i], dmaxdown);
              // markR[i]=-1;
            }
          } else {
            if (minR[i] > -1.0e40)
              printf("Why for %d minRi is %g not -1.0e60\n",
                i, minR[i]);
            if (maxR[i] < 1.0e40)
              printf("Why for %d maxRi is %g not 1.0e60\n",
                i, maxR[i]);
            // markR[i]=-1;
          }
        }
      }
#endif
      numberRealChanges = 0;
#endif
#if FIXED_BOTH_WAYS
      // clean newLo
      if (doFixedBothWays) {
        for (int i = 0; i < nFixedBoth; i++) {
          int icol = fixedBoth[i];
          newLo[icol] = COIN_DBL_MAX;
        }
        nFixedBoth = 0;
      }
#endif
      solval = colsol[j];
      down = floor(solval + tolerance);
      up = ceil(solval - tolerance);
      if (columnGap[j] < 0.0)
        markC[j] = 3;
      if ((markC[j] & 3) != 0 || !intVar[j])
        continue;
      double saveSolval = solval;
      if (solval >= colUpper[j] - tolerance || solval <= colLower[j] + tolerance || up == down) {
        if (solval <= colLower[j] + 2.0 * tolerance) {
          solval = colLower[j] + 1.0e-1;
          down = colLower[j];
          up = down + 1.0;
        } else if (solval >= colUpper[j] - 2.0 * tolerance) {
          solval = colUpper[j] - 1.0e-1;
          up = colUpper[j];
          down = up - 1;
        } else {
          // odd
          up = down + 1.0;
          solval = down + 1.0e-1;
        }
      }
      assert(up <= colUpper[j]);
      assert(down >= colLower[j]);
      assert(up > down);
      int istackC, iway, istackR;
      int way[] = { 1, 2, 1 };
      int feas[] = { 1, 2, 4 };
      int feasible = 0;
      int notFeasible;
      for (iway = 0; iway < 3; iway++) {
        // int fixThis=0;
        double objVal = current;
        int goingToTrueBound = 0;
        stackC[0] = j;
        markC[j] = way[iway];
        double solMovement;
        double movement;
        if (way[iway] == 1) {
          movement = down - colUpper[j];
          solMovement = down - colsol[j];
          assert(movement < -0.99999);
          if (fabs(down - colLower[j]) < 1.0e-7) {
            goingToTrueBound = 2;
            down = colLower[j];
          }
        } else {
          movement = up - colLower[j];
          solMovement = up - colsol[j];
          assert(movement > 0.99999);
          if (fabs(up - colUpper[j]) < 1.0e-7) {
            goingToTrueBound = 2;
            up = colUpper[j];
          }
        }
        if (goingToTrueBound && (colUpper[j] - colLower[j] > 1.5 || colLower[j]))
          goingToTrueBound = 1;
        // switch off disaggregation if not wanted
        if ((rowCuts & 1) == 0)
          goingToTrueBound = 0;
        bool canReplace = info->strengthenRow && (goingToTrueBound == 2);
#ifdef PRINT_DEBUG
        if (fabs(movement) > 1.01) {
          printf("big %d %g %g %g\n", j, colLower[j], solval, colUpper[j]);
        }
#endif
        if (solMovement * djs[j] > 0.0)
          objVal += solMovement * djs[j];
        nstackC = 1;
        nstackR = 0;
        saveL[0] = colLower[j];
        saveU[0] = colUpper[j];
        assert(saveU[0] > saveL[0]);
        notFeasible = 0;
        if (movement < 0.0) {
          colUpper[j] += movement;
          colUpper[j] = floor(colUpper[j] + 0.5);
          columnGap[j] = colUpper[j] - colLower[j] - primalTolerance_;
#ifdef PRINT_DEBUG
          printf("** Trying %d down to 0\n", j);
#endif
        } else {
          colLower[j] += movement;
          colLower[j] = floor(colLower[j] + 0.5);
          columnGap[j] = colUpper[j] - colLower[j] - primalTolerance_;
#ifdef PRINT_DEBUG
          printf("** Trying %d up to 1\n", j);
#endif
        }
        if (fabs(colUpper[j] - colLower[j]) < 1.0e-6)
          markC[j] = 3; // say fixed
        markC[j] &= ~12;
        if (colUpper[j] > 1.0e10)
          markC[j] |= 8;
        if (colLower[j] < -1.0e10)
          markC[j] |= 4;
        istackC = 0;
        /* update immediately */
        CoinBigIndex k;
        for (k = columnStart[j]; k < columnStart[j] + columnLength[j]; k++) {
          int irow = row[k];
          double value = columnElements[k];
          assert(markR[irow] != -2);
          if (markR[irow] == -1) {
            stackR[nstackR] = irow;
            markR[irow] = nstackR;
            saveMin[nstackR] = minR[irow];
            saveMax[nstackR] = maxR[irow];
            nstackR++;
#if 0
          } else if (markR[irow]==-2) {
            continue;
#endif
          }
          /* could check immediately if violation */
          if (movement > 0.0) {
            /* up */
            if (value > 0.0) {
              /* up does not change - down does */
              if (minR[irow] > -1.0e10)
                minR[irow] += value * movement;
              if (minR[irow] > rowUpper[irow] + 1.0e-5) {
                notFeasible = 1;
                istackC = 1;
                break;
              }
            } else {
              /* down does not change - up does */
              if (maxR[irow] < 1.0e10)
                maxR[irow] += value * movement;
              if (maxR[irow] < rowLower[irow] - 1.0e-5) {
                notFeasible = 1;
                istackC = 1;
                break;
              }
            }
          } else {
            /* down */
            if (value < 0.0) {
              /* up does not change - down does */
              if (minR[irow] > -1.0e10)
                minR[irow] += value * movement;
              if (minR[irow] > rowUpper[irow] + 1.0e-5) {
                notFeasible = 1;
                istackC = 1;
                break;
              }
            } else {
              /* down does not change - up does */
              if (maxR[irow] < 1.0e10)
                maxR[irow] += value * movement;
              if (maxR[irow] < rowLower[irow] - 1.0e-5) {
                notFeasible = 1;
                istackC = 1;
                break;
              }
            }
          }
        }
        while (istackC < nstackC && nstackC < maxStack) { // could be istackC<maxStack?
          int jway;
          int jcol = stackC[istackC];
          jway = markC[jcol];
          // If not first and fixed then skip
          if ((jway & 3) == 3 && istackC) {
            // istackC++;
            // continue;
            // printf("fixed %d on stack\n",jcol);
          }
#if PROBING100
          if (backward) {
            int jColumn = backward[jcol];
            if (jColumn >= 0) {
              int nFix = 0;
              // 0-1 see what else could be fixed
              if (jway == 1) {
                // fixed to 0
                int j;
                for (j = toZero_[jColumn]; j < toOne_[jColumn]; j++) {
                  int kColumn = fixEntry_[j].sequence;
                  kColumn = integerVariable_[kColumn];
                  bool fixToOne = fixEntry_[j].oneFixed;
                  if (fixToOne) {
                    if (colLower[kColumn] == 0.0) {
                      if (colUpper[kColumn] == 1.0) {
                        // See if on list
                        if (!(markC[kColumn] & 3)) {
                          if (nStackC < nCols) {
                            stackC[nstackC] = kColumn;
                            saveL[nstackC] = colLower[kColumn];
                            saveU[nstackC] = colUpper[kColumn];
                            assert(saveU[nstackC] > saveL[nstackC]);
                            assert(nstackC < nCols);
                            nstackC++;
                            markC[kColumn] |= 2;
                            nFix++;
                          }
                        } else if ((markC[kColumn] & 3) == 1) {
                          notFeasible = true;
                        }
                      } else {
                        // infeasible!
                        notFeasible = true;
                      }
                    }
                  } else {
                    if (colUpper[kColumn] == 1.0) {
                      if (colLower[kColumn] == 0.0) {
                        // See if on list
                        if (!(markC[kColumn] & 3)) {
                          if (nStackC < nCols) {
                            stackC[nstackC] = kColumn;
                            saveL[nstackC] = colLower[kColumn];
                            saveU[nstackC] = colUpper[kColumn];
                            assert(saveU[nstackC] > saveL[nstackC]);
                            assert(nstackC < nCols);
                            nstackC++;
                            markC[kColumn] |= 1;
                            nFix++;
                          }
                        } else if ((markC[kColumn] & 3) == 2) {
                          notFeasible = true;
                        }
                      } else {
                        // infeasible!
                        notFeasible = true;
                      }
                    }
                  }
                }
              } else if (jway == 2) {
                int j;
                for (j = toOne_[jColumn]; j < toZero_[jColumn + 1]; j++) {
                  int kColumn = fixEntry_[j].sequence;
                  kColumn = integerVariable_[kColumn];
                  bool fixToOne = fixEntry_[j].oneFixed;
                  if (fixToOne) {
                    if (colLower[kColumn] == 0.0) {
                      if (colUpper[kColumn] == 1.0) {
                        // See if on list
                        if (!(markC[kColumn] & 3)) {
                          if (nStackC < nCols) {
                            stackC[nstackC] = kColumn;
                            saveL[nstackC] = colLower[kColumn];
                            saveU[nstackC] = colUpper[kColumn];
                            assert(saveU[nstackC] > saveL[nstackC]);
                            assert(nstackC < nCols);
                            nstackC++;
                            markC[kColumn] |= 2;
                            nFix++;
                          }
                        } else if ((markC[kColumn] & 3) == 1) {
                          notFeasible = true;
                        }
                      } else {
                        // infeasible!
                        notFeasible = true;
                      }
                    }
                  } else {
                    if (colUpper[kColumn] == 1.0) {
                      if (colLower[kColumn] == 0.0) {
                        // See if on list
                        if (!(markC[kColumn] & 3)) {
                          if (nStackC < nCols) {
                            stackC[nstackC] = kColumn;
                            saveL[nstackC] = colLower[kColumn];
                            saveU[nstackC] = colUpper[kColumn];
                            assert(saveU[nstackC] > saveL[nstackC]);
                            assert(nstackC < nCols);
                            nstackC++;
                            markC[kColumn] |= 1;
                            nFix++;
                          }
                        } else if ((markC[kColumn] & 3) == 2) {
                          notFeasible = true;
                        }
                      } else {
                        // infeasible!
                        notFeasible = true;
                      }
                    }
                  }
                }
              }
            }
          }
#endif
          for (k = columnStart[jcol]; k < columnStart[jcol] + columnLength[jcol]; k++) {
            // break if found not feasible
            if (notFeasible)
              break;
            int irow = row[k];
            /* see if anything forced */
            CoinBigIndex rStart = rowStart[irow];
            CoinBigIndex rEnd = rowStartPos[irow];
            double rowUp = rowUpper[irow];
            double rowUp2 = 0.0;
            bool doRowUpN;
            bool doRowUpP;
            if (rowUp < 1.0e10) {
              doRowUpN = true;
              doRowUpP = true;
              rowUp2 = rowUp - minR[irow];
              if (rowUp2 < -tolerance) {
                notFeasible = true;
                break;
              } else {
                if (rowUp2 + largestNegativeInRow[irow] > 0)
                  doRowUpN = false;
                if (rowUp2 - largestPositiveInRow[irow] > 0)
                  doRowUpP = false;
              }
            } else {
              doRowUpN = false;
              doRowUpP = false;
              rowUp2 = COIN_DBL_MAX;
            }
            double rowLo = rowLower[irow];
            double rowLo2 = 0.0;
            bool doRowLoN;
            bool doRowLoP;
            if (rowLo > -1.0e10) {
              doRowLoN = true;
              doRowLoP = true;
              rowLo2 = rowLo - maxR[irow];
              if (rowLo2 > tolerance) {
                notFeasible = true;
                break;
              } else {
                if (rowLo2 - largestNegativeInRow[irow] < 0)
                  doRowLoN = false;
                if (rowLo2 + largestPositiveInRow[irow] < 0)
                  doRowLoP = false;
              }
            } else {
              doRowLoN = false;
              doRowLoP = false;
              rowLo2 = -COIN_DBL_MAX;
            }
            if (doRowUpN && doRowLoN) {
              // doRowUpN=doRowLoN=false;
              //  Start neg values loop
              for (CoinBigIndex kk = rStart; kk < rEnd; kk++) {
                int kcol = column[kk];
                int markIt = markC[kcol];
                if ((markIt & 3) != 3) {
                  double value2 = rowElements[kk];
                  if (colUpper[kcol] <= 1e10)
                    assert((markIt & 8) == 0);
                  else
                    assert((markIt & 8) != 0);
                  if (colLower[kcol] >= -1e10)
                    assert((markIt & 4) == 0);
                  else
                    assert((markIt & 4) != 0);
                  assert(value2 < 0.0);
                  double gap = columnGap[kcol] * value2;
                  bool doUp = (rowUp2 + gap < -1.0e-7);
                  bool doDown = (rowLo2 - gap > 1.0e-7);
                  if (doUp || doDown) {
                    double moveUp = 0.0;
                    double moveDown = 0.0;
                    double newUpper = -1.0;
                    double newLower = 1.0;
                    if (doUp && (markIt & (2 + 8)) == 0) {
                      double dbound = colUpper[kcol] + rowUp2 / value2;
                      if (intVar[kcol]) {
                        markIt |= 2;
                        newLower = ceil(dbound - relaxedTolerance);
                      } else {
                        newLower = dbound;
                        if (newLower + primalTolerance_ > colUpper[kcol] && newLower - primalTolerance_ <= colUpper[kcol]) {
                          newLower = colUpper[kcol];
                          markIt |= 2;
                          // markIt=3;
                        } else {
                          // avoid problems - fix later ?
                          markIt = 3;
                        }
                      }
                      moveUp = newLower - colLower[kcol];
                    }
                    if (doDown && (markIt & (1 + 4)) == 0) {
                      double dbound = colLower[kcol] + rowLo2 / value2;
                      if (intVar[kcol]) {
                        markIt |= 1;
                        newUpper = floor(dbound + relaxedTolerance);
                      } else {
                        newUpper = dbound;
                        if (newUpper - primalTolerance_ < colLower[kcol] && newUpper + primalTolerance_ >= colLower[kcol]) {
                          newUpper = colLower[kcol];
                          markIt |= 1;
                          // markIt=3;
                        } else {
                          // avoid problems - fix later ?
                          markIt = 3;
                        }
                      }
                      moveDown = newUpper - colUpper[kcol];
                    }
                    if (!moveUp && !moveDown)
                      continue;
                    bool onList = ((markC[kcol] & 3) != 0);
                    if (nstackC < 2 * maxStack) {
                      markC[kcol] = markIt;
                    }
                    if (moveUp && nstackC < 2 * maxStack) {
                      // fixThis++;
                      if (!onList) {
                        stackC[nstackC] = kcol;
                        saveL[nstackC] = colLower[kcol];
                        saveU[nstackC] = colUpper[kcol];
                        assert(saveU[nstackC] > saveL[nstackC]);
                        assert(nstackC < nCols);
                        nstackC++;
                        onList = true;
                      }
                      if (newLower > colsol[kcol]) {
                        if (djs[kcol] < 0.0) {
                          /* should be infeasible if integer */
                          if (intVar[kcol])
                            assert(newLower > colUpper[kcol] + primalTolerance_);
                        } else {
                          objVal += moveUp * djs[kcol];
                        }
                      }
                      if (intVar[kcol])
                        newLower = std::max(colLower[kcol], ceil(newLower - 1.0e-4));
                      colLower[kcol] = newLower;
                      columnGap[kcol] = colUpper[kcol] - newLower - primalTolerance_;
                      if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                        markC[kcol] = 3; // say fixed
                      }
                      markC[kcol] &= ~12;
                      if (colUpper[kcol] > 1.0e10)
                        markC[kcol] |= 8;
                      if (colLower[kcol] < -1.0e10)
                        markC[kcol] |= 4;
                      /* update immediately */
                      for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                        int krow = row[jj];
                        double value = columnElements[jj];
                        assert(markR[krow] != -2);
                        if (markR[krow] == -1) {
                          stackR[nstackR] = krow;
                          markR[krow] = nstackR;
                          saveMin[nstackR] = minR[krow];
                          saveMax[nstackR] = maxR[krow];
                          nstackR++;
                        }
                        /* could check immediately if violation */
                        /* up */
                        if (value > 0.0) {
                          /* up does not change - down does */
                          if (minR[krow] > -1.0e10)
                            minR[krow] += value * moveUp;
                          if (krow == irow)
                            rowUp2 = rowUp - minR[irow];
                          if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            break;
                          }
                        } else {
                          /* down does not change - up does */
                          if (maxR[krow] < 1.0e10)
                            maxR[krow] += value * moveUp;
                          if (krow == irow)
                            rowLo2 = rowLo - maxR[irow];
                          if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            break;
                          }
                        }
                      }
                    }
                    if (moveDown && nstackC < 2 * maxStack) {
                      // fixThis++;
                      if (!onList) {
                        stackC[nstackC] = kcol;
                        saveL[nstackC] = colLower[kcol];
                        saveU[nstackC] = colUpper[kcol];
                        assert(saveU[nstackC] > saveL[nstackC]);
                        assert(nstackC < nCols);
                        nstackC++;
                        onList = true;
                      }
                      if (newUpper < colsol[kcol]) {
                        if (djs[kcol] > 0.0) {
                          /* should be infeasible if integer */
                          if (intVar[kcol])
                            assert(colLower[kcol] > newUpper + primalTolerance_);
                        } else {
                          objVal += moveDown * djs[kcol];
                        }
                      }
                      if (intVar[kcol])
                        newUpper = std::min(colUpper[kcol], floor(newUpper + 1.0e-4));
                      colUpper[kcol] = newUpper;
                      columnGap[kcol] = newUpper - colLower[kcol] - primalTolerance_;
                      if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                        markC[kcol] = 3; // say fixed
                      }
                      markC[kcol] &= ~12;
                      if (colUpper[kcol] > 1.0e10)
                        markC[kcol] |= 8;
                      if (colLower[kcol] < -1.0e10)
                        markC[kcol] |= 4;
                      /* update immediately */
                      for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                        int krow = row[jj];
                        double value = columnElements[jj];
                        assert(markR[krow] != -2);
                        if (markR[krow] == -1) {
                          stackR[nstackR] = krow;
                          markR[krow] = nstackR;
                          saveMin[nstackR] = minR[krow];
                          saveMax[nstackR] = maxR[krow];
                          nstackR++;
                        }
                        /* could check immediately if violation */
                        /* down */
                        if (value < 0.0) {
                          /* up does not change - down does */
                          if (minR[krow] > -1.0e10)
                            minR[krow] += value * moveDown;
                          if (krow == irow)
                            rowUp2 = rowUp - minR[irow];
                          if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        } else {
                          /* down does not change - up does */
                          if (maxR[krow] < 1.0e10)
                            maxR[krow] += value * moveDown;
                          if (krow == irow)
                            rowLo2 = rowLo - maxR[irow];
                          if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        }
                      }
                    }
                    if (colLower[kcol] > colUpper[kcol] + primalTolerance_) {
                      notFeasible = 1;
                      ;
                      k = columnStart[jcol] + columnLength[jcol];
                      istackC = nstackC + 1;
                      break;
                    }
                  }
                }
              } // end big loop rStart->rPos
            } else if (doRowUpN) {
              // Start neg values loop
              for (CoinBigIndex kk = rStart; kk < rEnd; kk++) {
                int kcol = column[kk];
                int markIt = markC[kcol];
                if ((markIt & 3) != 3) {
                  double value2 = rowElements[kk];
                  double gap = columnGap[kcol] * value2;
                  if (!(rowUp2 + gap < -1.0e-7))
                    continue;
                  double moveUp = 0.0;
                  double newLower = 1.0;
                  if ((markIt & (2 + 8)) == 0) {
                    double dbound = colUpper[kcol] + rowUp2 / value2;
                    if (intVar[kcol]) {
                      markIt |= 2;
                      newLower = ceil(dbound - relaxedTolerance);
                    } else {
                      newLower = dbound;
                      if (newLower + primalTolerance_ > colUpper[kcol] && newLower - primalTolerance_ <= colUpper[kcol]) {
                        newLower = colUpper[kcol];
                        markIt |= 2;
                        // markIt=3;
                      } else {
                        // avoid problems - fix later ?
                        markIt = 3;
                      }
                    }
                    moveUp = newLower - colLower[kcol];
                    if (!moveUp)
                      continue;
                    bool onList = ((markC[kcol] & 3) != 0);
                    if (nstackC < 2 * maxStack) {
                      markC[kcol] = markIt;
                    }
                    if (moveUp && nstackC < 2 * maxStack) {
                      // fixThis++;
                      if (!onList) {
                        stackC[nstackC] = kcol;
                        saveL[nstackC] = colLower[kcol];
                        saveU[nstackC] = colUpper[kcol];
                        assert(saveU[nstackC] > saveL[nstackC]);
                        assert(nstackC < nCols);
                        nstackC++;
                        onList = true;
                      }
                      if (newLower > colsol[kcol]) {
                        if (djs[kcol] < 0.0) {
                          /* should be infeasible if integer */
                          if (intVar[kcol])
                            assert(newLower > colUpper[kcol] + primalTolerance_);
                        } else {
                          objVal += moveUp * djs[kcol];
                        }
                      }
                      if (intVar[kcol])
                        newLower = std::max(colLower[kcol], ceil(newLower - 1.0e-4));
                      colLower[kcol] = newLower;
                      columnGap[kcol] = colUpper[kcol] - newLower - primalTolerance_;
                      if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                        markC[kcol] = 3; // say fixed
                      }
                      markC[kcol] &= ~12;
                      if (colUpper[kcol] > 1.0e10)
                        markC[kcol] |= 8;
                      if (colLower[kcol] < -1.0e10)
                        markC[kcol] |= 4;
                      /* update immediately */
                      for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                        int krow = row[jj];
                        double value = columnElements[jj];
                        assert(markR[krow] != -2);
                        if (markR[krow] == -1) {
                          stackR[nstackR] = krow;
                          markR[krow] = nstackR;
                          saveMin[nstackR] = minR[krow];
                          saveMax[nstackR] = maxR[krow];
                          nstackR++;
                        }
                        /* could check immediately if violation */
                        /* up */
                        if (value > 0.0) {
                          /* up does not change - down does */
                          if (minR[krow] > -1.0e10)
                            minR[krow] += value * moveUp;
                          if (krow == irow)
                            rowUp2 = rowUp - minR[irow];
                          if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        } else {
                          /* down does not change - up does */
                          if (maxR[krow] < 1.0e10)
                            maxR[krow] += value * moveUp;
                          if (krow == irow)
                            rowLo2 = rowLo - maxR[irow];
                          if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        }
                      }
                    }
                    if (colLower[kcol] > colUpper[kcol] + primalTolerance_) {
                      notFeasible = 1;
                      ;
                      k = columnStart[jcol] + columnLength[jcol];
                      istackC = nstackC + 1;
                      break;
                    }
                  }
                }
              } // end big loop rStart->rPos
            } else if (doRowLoN) {
              // Start neg values loop
              for (CoinBigIndex kk = rStart; kk < rEnd; kk++) {
                int kcol = column[kk];
                if ((markC[kcol] & 3) != 3) {
                  double moveDown = 0.0;
                  double newUpper = -1.0;
                  double value2 = rowElements[kk];
                  int markIt = markC[kcol];
                  assert(value2 < 0.0);
                  double gap = columnGap[kcol] * value2;
                  bool doDown = (rowLo2 - gap > 1.0e-7);
                  if (doDown && (markIt & (1 + 4)) == 0) {
                    double dbound = colLower[kcol] + rowLo2 / value2;
                    if (intVar[kcol]) {
                      markIt |= 1;
                      newUpper = floor(dbound + relaxedTolerance);
                    } else {
                      newUpper = dbound;
                      if (newUpper - primalTolerance_ < colLower[kcol] && newUpper + primalTolerance_ >= colLower[kcol]) {
                        newUpper = colLower[kcol];
                        markIt |= 1;
                        // markIt=3;
                      } else {
                        // avoid problems - fix later ?
                        markIt = 3;
                      }
                    }
                    moveDown = newUpper - colUpper[kcol];
                    if (!moveDown)
                      continue;
                    bool onList = ((markC[kcol] & 3) != 0);
                    if (nstackC < 2 * maxStack) {
                      markC[kcol] = markIt;
                    }
                    if (moveDown && nstackC < 2 * maxStack) {
                      // fixThis++;
                      if (!onList) {
                        stackC[nstackC] = kcol;
                        saveL[nstackC] = colLower[kcol];
                        saveU[nstackC] = colUpper[kcol];
                        assert(saveU[nstackC] > saveL[nstackC]);
                        assert(nstackC < nCols);
                        nstackC++;
                        onList = true;
                      }
                      if (newUpper < colsol[kcol]) {
                        if (djs[kcol] > 0.0) {
                          /* should be infeasible if integer */
                          if (intVar[kcol])
                            assert(colLower[kcol] > newUpper + primalTolerance_);
                        } else {
                          objVal += moveDown * djs[kcol];
                        }
                      }
                      if (intVar[kcol])
                        newUpper = std::min(colUpper[kcol], floor(newUpper + 1.0e-4));
                      colUpper[kcol] = newUpper;
                      columnGap[kcol] = newUpper - colLower[kcol] - primalTolerance_;
                      if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                        markC[kcol] = 3; // say fixed
                      }
                      markC[kcol] &= ~12;
                      if (colUpper[kcol] > 1.0e10)
                        markC[kcol] |= 8;
                      if (colLower[kcol] < -1.0e10)
                        markC[kcol] |= 4;
                      /* update immediately */
                      for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                        int krow = row[jj];
                        double value = columnElements[jj];
                        assert(markR[krow] != -2);
                        if (markR[krow] == -1) {
                          stackR[nstackR] = krow;
                          markR[krow] = nstackR;
                          saveMin[nstackR] = minR[krow];
                          saveMax[nstackR] = maxR[krow];
                          nstackR++;
                        }
                        /* could check immediately if violation */
                        /* down */
                        if (value < 0.0) {
                          /* up does not change - down does */
                          if (minR[krow] > -1.0e10)
                            minR[krow] += value * moveDown;
                          if (krow == irow)
                            rowUp2 = rowUp - minR[irow];
                          if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        } else {
                          /* down does not change - up does */
                          if (maxR[krow] < 1.0e10)
                            maxR[krow] += value * moveDown;
                          if (krow == irow)
                            rowLo2 = rowLo - maxR[irow];
                          if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        }
                      }
                    }
                    if (colLower[kcol] > colUpper[kcol] + primalTolerance_) {
                      notFeasible = 1;
                      ;
                      k = columnStart[jcol] + columnLength[jcol];
                      istackC = nstackC + 1;
                      break;
                    }
                  }
                }
              } // end big loop rStart->rPos
            }
            rStart = rEnd;
            rEnd = rowStart[irow + 1];
            if (doRowUpP && doRowLoP) {
              // doRowUpP=doRowLoP=false;
              //  Start pos values loop
              for (CoinBigIndex kk = rStart; kk < rEnd; kk++) {
                int kcol = column[kk];
                int markIt = markC[kcol];
                if ((markIt & 3) != 3) {
                  double value2 = rowElements[kk];
                  assert(value2 > 0.0);
                  /* positive element */
                  double gap = columnGap[kcol] * value2;
                  bool doDown = (rowLo2 + gap > 1.0e-7);
                  bool doUp = (rowUp2 - gap < -1.0e-7);
                  if (doDown || doUp) {
                    double moveUp = 0.0;
                    double moveDown = 0.0;
                    double newUpper = -1.0;
                    double newLower = 1.0;
                    if (doDown && (markIt & (2 + 8)) == 0) {
                      double dbound = colUpper[kcol] + rowLo2 / value2;
                      if (intVar[kcol]) {
                        markIt |= 2;
                        newLower = ceil(dbound - relaxedTolerance);
                      } else {
                        newLower = dbound;
                        if (newLower + primalTolerance_ > colUpper[kcol] && newLower - primalTolerance_ <= colUpper[kcol]) {
                          newLower = colUpper[kcol];
                          markIt |= 2;
                          // markIt=3;
                        } else {
                          // avoid problems - fix later ?
                          markIt = 3;
                        }
                      }
                      moveUp = newLower - colLower[kcol];
                    }
                    if (doUp && (markIt & (1 + 4)) == 0) {
                      double dbound = colLower[kcol] + rowUp2 / value2;
                      if (intVar[kcol]) {
                        markIt |= 1;
                        newUpper = floor(dbound + relaxedTolerance);
                      } else {
                        newUpper = dbound;
                        if (newUpper - primalTolerance_ < colLower[kcol] && newUpper + primalTolerance_ >= colLower[kcol]) {
                          newUpper = colLower[kcol];
                          markIt |= 1;
                          // markIt=3;
                        } else {
                          // avoid problems - fix later ?
                          markIt = 3;
                        }
                      }
                      moveDown = newUpper - colUpper[kcol];
                    }
                    if (!moveUp && !moveDown)
                      continue;
                    bool onList = ((markC[kcol] & 3) != 0);
                    if (nstackC < 2 * maxStack) {
                      markC[kcol] = markIt;
                    }
                    if (moveUp && nstackC < 2 * maxStack) {
                      // fixThis++;
                      if (!onList) {
                        stackC[nstackC] = kcol;
                        saveL[nstackC] = colLower[kcol];
                        saveU[nstackC] = colUpper[kcol];
                        assert(saveU[nstackC] > saveL[nstackC]);
                        assert(nstackC < nCols);
                        nstackC++;
                        onList = true;
                      }
                      if (newLower > colsol[kcol]) {
                        if (djs[kcol] < 0.0) {
                          /* should be infeasible if integer */
                          if (intVar[kcol])
                            assert(newLower > colUpper[kcol] + primalTolerance_);
                        } else {
                          objVal += moveUp * djs[kcol];
                        }
                      }
                      if (intVar[kcol])
                        newLower = std::max(colLower[kcol], ceil(newLower - 1.0e-4));
                      colLower[kcol] = newLower;
                      columnGap[kcol] = colUpper[kcol] - newLower - primalTolerance_;
                      if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                        markC[kcol] = 3; // say fixed
                      }
                      markC[kcol] &= ~12;
                      if (colUpper[kcol] > 1.0e10)
                        markC[kcol] |= 8;
                      if (colLower[kcol] < -1.0e10)
                        markC[kcol] |= 4;
                      /* update immediately */
                      for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                        int krow = row[jj];
                        double value = columnElements[jj];
                        assert(markR[krow] != -2);
                        if (markR[krow] == -1) {
                          stackR[nstackR] = krow;
                          markR[krow] = nstackR;
                          saveMin[nstackR] = minR[krow];
                          saveMax[nstackR] = maxR[krow];
                          nstackR++;
                        }
                        /* could check immediately if violation */
                        /* up */
                        if (value > 0.0) {
                          /* up does not change - down does */
                          if (minR[krow] > -1.0e10)
                            minR[krow] += value * moveUp;
                          if (krow == irow)
                            rowUp2 = rowUp - minR[irow];
                          if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        } else {
                          /* down does not change - up does */
                          if (maxR[krow] < 1.0e10)
                            maxR[krow] += value * moveUp;
                          if (krow == irow)
                            rowLo2 = rowLo - maxR[irow];
                          if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        }
                      }
                    }
                    if (moveDown && nstackC < 2 * maxStack) {
                      // fixThis++;
                      if (!onList) {
                        stackC[nstackC] = kcol;
                        saveL[nstackC] = colLower[kcol];
                        saveU[nstackC] = colUpper[kcol];
                        assert(saveU[nstackC] > saveL[nstackC]);
                        assert(nstackC < nCols);
                        nstackC++;
                        onList = true;
                      }
                      if (newUpper < colsol[kcol]) {
                        if (djs[kcol] > 0.0) {
                          /* should be infeasible if integer */
                          if (intVar[kcol])
                            assert(colLower[kcol] > newUpper + primalTolerance_);
                        } else {
                          objVal += moveDown * djs[kcol];
                        }
                      }
                      if (intVar[kcol])
                        newUpper = std::min(colUpper[kcol], floor(newUpper + 1.0e-4));
                      colUpper[kcol] = newUpper;
                      columnGap[kcol] = newUpper - colLower[kcol] - primalTolerance_;
                      if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                        markC[kcol] = 3; // say fixed
                      }
                      markC[kcol] &= ~12;
                      if (colUpper[kcol] > 1.0e10)
                        markC[kcol] |= 8;
                      if (colLower[kcol] < -1.0e10)
                        markC[kcol] |= 4;
                      /* update immediately */
                      for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                        int krow = row[jj];
                        double value = columnElements[jj];
                        assert(markR[krow] != -2);
                        if (markR[krow] == -1) {
                          stackR[nstackR] = krow;
                          markR[krow] = nstackR;
                          saveMin[nstackR] = minR[krow];
                          saveMax[nstackR] = maxR[krow];
                          nstackR++;
                        }
                        /* could check immediately if violation */
                        /* down */
                        if (value < 0.0) {
                          /* up does not change - down does */
                          if (minR[krow] > -1.0e10)
                            minR[krow] += value * moveDown;
                          if (krow == irow)
                            rowUp2 = rowUp - minR[irow];
                          if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        } else {
                          /* down does not change - up does */
                          if (maxR[krow] < 1.0e10)
                            maxR[krow] += value * moveDown;
                          if (krow == irow)
                            rowLo2 = rowLo - maxR[irow];
                          if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        }
                      }
                    }
                    if (colLower[kcol] > colUpper[kcol] + primalTolerance_) {
                      notFeasible = 1;
                      ;
                      k = columnStart[jcol] + columnLength[jcol];
                      istackC = nstackC + 1;
                      break;
                    }
                  }
                }
              } // end big loop rPos->rEnd
            } else if (doRowUpP) {
              // Start pos values loop
              for (CoinBigIndex kk = rStart; kk < rEnd; kk++) {
                int kcol = column[kk];
                int markIt = markC[kcol];
                if ((markIt & 3) != 3) {
                  double value2 = rowElements[kk];
                  assert(value2 > 0.0);
                  /* positive element */
                  double gap = columnGap[kcol] * value2;
                  bool doUp = (rowUp2 - gap < -1.0e-7);
                  if (doUp && (markIt & (1 + 4)) == 0) {
                    double newUpper = -1.0;
                    double dbound = colLower[kcol] + rowUp2 / value2;
                    if (intVar[kcol]) {
                      markIt |= 1;
                      newUpper = floor(dbound + relaxedTolerance);
                    } else {
                      newUpper = dbound;
                      if (newUpper - primalTolerance_ < colLower[kcol] && newUpper + primalTolerance_ >= colLower[kcol]) {
                        newUpper = colLower[kcol];
                        markIt |= 1;
                        // markIt=3;
                      } else {
                        // avoid problems - fix later ?
                        markIt = 3;
                      }
                    }
                    double moveDown = newUpper - colUpper[kcol];
                    if (!moveDown)
                      continue;
                    bool onList = ((markC[kcol] & 3) != 0);
                    if (nstackC < 2 * maxStack) {
                      markC[kcol] = markIt;
                    }
                    if (nstackC < 2 * maxStack) {
                      // fixThis++;
                      if (!onList) {
                        stackC[nstackC] = kcol;
                        saveL[nstackC] = colLower[kcol];
                        saveU[nstackC] = colUpper[kcol];
                        assert(saveU[nstackC] > saveL[nstackC]);
                        assert(nstackC < nCols);
                        nstackC++;
                        onList = true;
                      }
                      if (newUpper < colsol[kcol]) {
                        if (djs[kcol] > 0.0) {
                          /* should be infeasible if integer */
                          if (intVar[kcol])
                            assert(colLower[kcol] > newUpper + primalTolerance_);
                        } else {
                          objVal += moveDown * djs[kcol];
                        }
                      }
                      if (intVar[kcol])
                        newUpper = std::min(colUpper[kcol], floor(newUpper + 1.0e-4));
                      colUpper[kcol] = newUpper;
                      columnGap[kcol] = newUpper - colLower[kcol] - primalTolerance_;
                      if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                        markC[kcol] = 3; // say fixed
                      }
                      markC[kcol] &= ~12;
                      if (colUpper[kcol] > 1.0e10)
                        markC[kcol] |= 8;
                      if (colLower[kcol] < -1.0e10)
                        markC[kcol] |= 4;
                      /* update immediately */
                      for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                        int krow = row[jj];
                        double value = columnElements[jj];
                        assert(markR[krow] != -2);
                        if (markR[krow] == -1) {
                          stackR[nstackR] = krow;
                          markR[krow] = nstackR;
                          saveMin[nstackR] = minR[krow];
                          saveMax[nstackR] = maxR[krow];
                          nstackR++;
                        }
                        /* could check immediately if violation */
                        /* down */
                        if (value < 0.0) {
                          /* up does not change - down does */
                          if (minR[krow] > -1.0e10)
                            minR[krow] += value * moveDown;
                          if (krow == irow)
                            rowUp2 = rowUp - minR[irow];
                          if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        } else {
                          /* down does not change - up does */
                          if (maxR[krow] < 1.0e10)
                            maxR[krow] += value * moveDown;
                          if (krow == irow)
                            rowLo2 = rowLo - maxR[irow];
                          if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        }
                      }
                    }
                    if (colLower[kcol] > colUpper[kcol] + primalTolerance_) {
                      notFeasible = 1;
                      ;
                      k = columnStart[jcol] + columnLength[jcol];
                      istackC = nstackC + 1;
                      break;
                    }
                  }
                }
              } // end big loop rPos->rEnd
            } else if (doRowLoP) {
              // Start pos values loop
              for (CoinBigIndex kk = rStart; kk < rEnd; kk++) {
                int kcol = column[kk];
                if ((markC[kcol] & 3) != 3) {
                  double value2 = rowElements[kk];
                  int markIt = markC[kcol];
                  assert(value2 > 0.0);
                  /* positive element */
                  double gap = columnGap[kcol] * value2;
                  bool doDown = (rowLo2 + gap > 1.0e-7);
                  if (doDown && (markIt & (2 + 8)) == 0) {
                    double newLower = 1.0;
                    double dbound = colUpper[kcol] + rowLo2 / value2;
                    if (intVar[kcol]) {
                      markIt |= 2;
                      newLower = ceil(dbound - relaxedTolerance);
                    } else {
                      newLower = dbound;
                      if (newLower + primalTolerance_ > colUpper[kcol] && newLower - primalTolerance_ <= colUpper[kcol]) {
                        newLower = colUpper[kcol];
                        markIt |= 2;
                        // markIt=3;
                      } else {
                        // avoid problems - fix later ?
                        markIt = 3;
                      }
                    }
                    double moveUp = newLower - colLower[kcol];
                    if (!moveUp)
                      continue;
                    bool onList = ((markC[kcol] & 3) != 0);
                    if (nstackC < 2 * maxStack) {
                      markC[kcol] = markIt;
                    }
                    if (nstackC < 2 * maxStack) {
                      // fixThis++;
                      if (!onList) {
                        stackC[nstackC] = kcol;
                        saveL[nstackC] = colLower[kcol];
                        saveU[nstackC] = colUpper[kcol];
                        assert(saveU[nstackC] > saveL[nstackC]);
                        assert(nstackC < nCols);
                        nstackC++;
                        onList = true;
                      }
                      if (newLower > colsol[kcol]) {
                        if (djs[kcol] < 0.0) {
                          /* should be infeasible if integer */
                          if (intVar[kcol])
                            assert(newLower > colUpper[kcol] + primalTolerance_);
                        } else {
                          objVal += moveUp * djs[kcol];
                        }
                      }
                      if (intVar[kcol])
                        newLower = std::max(colLower[kcol], ceil(newLower - 1.0e-4));
                      colLower[kcol] = newLower;
                      columnGap[kcol] = colUpper[kcol] - newLower - primalTolerance_;
                      if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                        markC[kcol] = 3; // say fixed
                      }
                      markC[kcol] &= ~12;
                      if (colUpper[kcol] > 1.0e10)
                        markC[kcol] |= 8;
                      if (colLower[kcol] < -1.0e10)
                        markC[kcol] |= 4;
                      /* update immediately */
                      for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                        int krow = row[jj];
                        double value = columnElements[jj];
                        assert(markR[krow] != -2);
                        if (markR[krow] == -1) {
                          stackR[nstackR] = krow;
                          markR[krow] = nstackR;
                          saveMin[nstackR] = minR[krow];
                          saveMax[nstackR] = maxR[krow];
                          nstackR++;
                        }
                        /* could check immediately if violation */
                        /* up */
                        if (value > 0.0) {
                          /* up does not change - down does */
                          if (minR[krow] > -1.0e10)
                            minR[krow] += value * moveUp;
                          if (krow == irow)
                            rowUp2 = rowUp - minR[irow];
                          if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        } else {
                          /* down does not change - up does */
                          if (maxR[krow] < 1.0e10)
                            maxR[krow] += value * moveUp;
                          if (krow == irow)
                            rowLo2 = rowLo - maxR[irow];
                          if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        }
                      }
                    }
                    if (colLower[kcol] > colUpper[kcol] + primalTolerance_) {
                      notFeasible = 1;
                      ;
                      k = columnStart[jcol] + columnLength[jcol];
                      istackC = nstackC + 1;
                      break;
                    }
                  }
                }
              } // end big loop rPos->rEnd
            }
          }
          istackC++;
        }
#if FIXED_BOTH_WAYS || MANY_TIGHTEN
        // the nCols<1000 should be on some switch
        if (!notFeasible && iway < 2 && (info->options & 2048) != 0
          && nCols < 1000) {
          // try tighten
          double *COIN_RESTRICT rLower = const_cast< double * COIN_RESTRICT >(rowLower);
          double *COIN_RESTRICT rUpper = const_cast< double * COIN_RESTRICT >(rowUpper);
          // must be able to do better
          memcpy(saveFColLower, colLower, nCols * sizeof(double));
          memcpy(saveFColUpper, colUpper, nCols * sizeof(double));
          memcpy(saveFRowLower, rLower, nRows * sizeof(double));
          memcpy(saveFRowUpper, rUpper, nRows * sizeof(double));
          // move stackR
          memcpy(lookRow_[0], stackR, nstackR * sizeof(int));
          numberRows_ = nstackR;
          int returnCode = tighten(colLower, colUpper,
            column, rowElements, rowStart,
            rowStartPos, rowLength, columnCopy, rLower, rUpper,
            nRows, nCols, intVar, 10, primalTolerance_);
          if (returnCode) {
            // printf("INFEAS for j %d and way %d\n",j,iway);
            notFeasible = 1;
          }
#if 0
	  // maybe should do more than nstackC
	  for (int i=1;i<nstackC;i++) {
	    int icol = stackC[i];
	    assert (saveFColLower[icol]<=colLower[icol]);
	    assert (saveFColUpper[icol]>=colUpper[icol]);
	    saveFColLower[icol]=colLower[icol];
	    saveFColUpper[icol]>=colUpper[icol];
	  }
#endif
          memcpy(colLower, saveFColLower, nCols * sizeof(double));
          memcpy(colUpper, saveFColUpper, nCols * sizeof(double));
          memcpy(rLower, saveFRowLower, nRows * sizeof(double));
          memcpy(rUpper, saveFRowUpper, nRows * sizeof(double));
        }
#endif
#if CHECK_MINR_MAXR
        // double check
        if (!notFeasible && iway != 2) { // not sure about 2
#if CHECK_MINR_MAXR > 0
          writeMps(j, iway, &si, colLower, colUpper, rowCopy, columnCopy,
            rowLower, rowUpper, intVar);
#endif
          int verbose = 0;
          int checkBefore = -1;
          if (CHECK_MINR_MAXR > 10000) {
            verbose = 2;
            checkBefore = CHECK_MINR_MAXR - 10000;
          } else if (CHECK_MINR_MAXR > 0) {
            verbose = 1;
            checkBefore = CHECK_MINR_MAXR;
          } else if (CHECK_MINR_MAXR < 0) {
            verbose = 0;
            checkBefore = 0;
          }
          if (checkBefore)
            printf("Checking column %d way %d\n", j, iway);
          double *COIN_RESTRICT maxRow2 = new double[2 * nRows];
          double *COIN_RESTRICT minRow2 = maxRow2 + nRows;
          memset(maxRow2, 0, 2 * nRows * sizeof(double));
          int nRows2 = nRows;
          for (int i = 0; i < nCols; i++) {
            double low = std::max(colLower[i], -1.0e20);
            double high = std::min(colUpper[i], 1.0e20);
            if (high < -1.020) {
              nRows2 = 0;
              break;
            }
            CoinBigIndex k;
            for (k = columnStart[i]; k < columnStart[i] + columnLength[i]; k++) {
              int irow = row[k];
              double value = columnElements[k];
              if (value < 0.0) {
                maxRow2[irow] += value * low;
                minRow2[irow] += value * high;
              } else {
                maxRow2[irow] += value * high;
                minRow2[irow] += value * low;
              }
            }
          }
          for (int i = 0; i < nRows2; i++) {
            bool odd = false;
            if (fabs(minR[i] - minRow2[i]) > 1.0e-5) {
              if (minR[i] > -1.0e5 || minRow2[i] > -1.0e5)
                odd = true;
            }
            if (fabs(maxR[i] - maxRow2[i]) > 1.0e-5) {
              if (maxR[i] < 1.0e5 || maxRow2[i] < 1.0e5)
                odd = true;
            }
            if (minR[i] == -1.0e60 && maxR[i] == 1.0e60)
              odd = false;
            if (odd) {
              printf("col %d row %d now %g <= %g was %g <= %g\n",
                j, i, minRow2[i], maxRow2[i], minR[i], maxR[i]);
            } else if (j == -4) {
              printf("row %d lo %g min %g max %g up %g\n",
                i, rowLower[i], minRow2[i], maxRow2[i], rowUpper[i]);
            }
          }
          for (int k = 0; k < nstackR; k++) {
            int krow = stackR[k];
            if (checkBefore)
              printf("row %d minR,maxR %g,%g l,u %g,%g\n",
                krow, minR[krow], maxR[krow], rowLower[krow], rowUpper[krow]);
            if (minR[krow] > rowUpper[krow] + 1.0e-3) {
              // printf("row %d lo %g min %g max %g up %g\n",
              //	   krow,rowLower[krow],minRow2[krow],maxRow2[krow],rowUpper[krow]);
              notFeasible = true;
            } else if (maxR[krow] < rowLower[krow] - 1.0e-3) {
              // printf("row %d lo %g min %g max %g up %g\n",
              //	   krow,rowLower[krow],minRow2[krow],maxRow2[krow],rowUpper[krow]);
              notFeasible = true;
            }
          }
          delete[] maxRow2;
          if (verbose == 1) {
            for (int istackC = 0; istackC < nstackC; istackC++) {
              int icol = stackC[istackC];
              printf("col %d current clo,cup %g,%g orig %g,%g\n",
                icol, colLower[icol], colUpper[icol],
                saveL[istackC], saveU[istackC]);
            }
          } else if (verbose == 2) {
            double loX[1000], upX[1000];
            assert(nCols < 1000);
            for (int i = 0; i < nCols; i++)
              loX[i] = COIN_DBL_MAX;
            for (int istackC = 0; istackC < nstackC; istackC++) {
              int icol = stackC[istackC];
              loX[icol] = saveL[istackC];
              upX[icol] = saveU[istackC];
            }
            for (int icol = 0; icol < nCols; icol++) {
              if (loX[icol] == COIN_DBL_MAX) {
                printf("icolj col %d j %d clo,cup %g,%g jway %d\n",
                  icol, j, colLower[icol], colUpper[icol], iway);
                printf("jicol j %d col %d clo,cup %g,%g jway %d\n",
                  j, icol, colLower[icol], colUpper[icol], iway);
              } else {
                printf("icolj col %d j %d current clo,cup %g,%g orig %g,%g jway %d\n",
                  icol, j, colLower[icol], colUpper[icol],
                  loX[icol], upX[icol], iway);
                printf("jicol j %d col %d current clo,cup %g,%g orig %g,%g jway %d\n",
                  j, icol, colLower[icol], colUpper[icol],
                  loX[icol], upX[icol], iway);
              }
            }
          }
        } else {
          bool verbose = false;
          verbose = (j < CHECK_MINR_MAXR);
          if (verbose)
            printf("Checking column %d way %d - INFEAS\n", j, iway);
        }
#endif
        if (!notFeasible) {
          if (objVal <= cutoff) {
            feasible |= feas[iway];
          } else {
#ifdef PRINT_DEBUG
            printf("not feasible on dj\n");
#endif
            notFeasible = 1;
            if (iway == 1 && feasible == 0) {
              /* not feasible at all */
              ninfeas = 1;
              j = nCols - 1;
              break;
            }
          }
          if (rowCutFake && !notFeasible && iway == 1) {
            assert(j == stackC[0]);
            // add to cliques
            for (istackC = 1; istackC < nstackC; istackC++) {
              int icol = stackC[istackC];
              // for now back to just 0-1
              if (!colUpper[icol] && !saveL[istackC] && saveU[istackC] == 1.0) {
                OsiRowCut rc;
                rc.setLb(-COIN_DBL_MAX);
                rc.setUb(1.0);
                index[0] = icol;
                element[0] = 1.0;
                index[1] = j;
                element[1] = 1.0;
                rc.setRow(2, index, element, false);
                rowCutFake->addCutIfNotDuplicate(rc);
                // printf("XX x%d + x%d <= 1\n",j,icol);
              }
            }
          }
          if (!notFeasible && saveFixingInfo) {
            // save fixing info
            assert(j == stackC[0]);
            int toValue = (way[iway] == 1) ? -1 : +1;
            for (istackC = 1; istackC < nstackC; istackC++) {
              int icol = stackC[istackC];
              // for now back to just 0-1
              if (colUpper[icol] - colLower[icol] < 1.0e-12 && !saveL[istackC] && saveU[istackC] == 1.0) {
                assert(saveL[istackC] == colLower[icol] || saveU[istackC] == colUpper[icol]);
                saveFixingInfo = info->fixes(j, toValue,
                  icol, colLower[icol] == saveL[istackC]);
              }
            }
          }
        } else if (iway == 1 && feasible == 0) {
          /* not feasible at all */
          ninfeas = 1;
          j = nCols - 1;
          iLook = numberThisTime_;
          ipass = maxPass;
          break;
        }
        if (notFeasible)
          goingToTrueBound = 0;
        if (iway == 2 || (iway == 1 && feasible == 2)) {
          /* keep */
          iway = 3;
          nfixed++;
          OsiColCut *cc = new OsiColCut;
          int /*nTot=0,*/ nFix = 0, nInt = 0;
          bool ifCut = false;
          for (istackC = 0; istackC < nstackC; istackC++) {
            int icol = stackC[istackC];
            if (intVar[icol]) {
              if (colUpper[icol] < currentColUpper[icol] - 1.0e-4) {
                element[nFix] = colUpper[icol];
                index[nFix++] = icol;
                nInt++;
                if (colsol[icol] > colUpper[icol] + primalTolerance_) {
                  ifCut = true;
                  anyColumnCuts = true;
                }
              }
            }
          }
          if (nFix) {
            // nTot=nFix;
            cleanBounds(1, nFix, index, element);
            cc->setUbs(nFix, index, element);
            nFix = 0;
          }
          for (istackC = 0; istackC < nstackC; istackC++) {
            int icol = stackC[istackC];
            if (intVar[icol]) {
              if (colLower[icol] > currentColLower[icol] + 1.0e-4) {
                element[nFix] = colLower[icol];
                index[nFix++] = icol;
                nInt++;
                if (colsol[icol] < colLower[icol] - primalTolerance_) {
                  ifCut = true;
                  anyColumnCuts = true;
                }
              }
            }
          }
          if (nFix) {
            // nTot+=nFix;
            cleanBounds(-1, nFix, index, element);
            cc->setLbs(nFix, index, element);
          }
          // could tighten continuous as well
          if (nInt) {
            if (ifCut) {
              cc->setEffectiveness(100.0);
            } else {
              cc->setEffectiveness(1.0e-5);
            }
#ifdef CGL_DEBUG
            checkBounds(debugger, *cc);
#endif
            cs.insert(cc);
          } else {
            delete cc;
          }
          for (istackC = 0; istackC < nstackC; istackC++) {
            int icol = stackC[istackC];
            if (colUpper[icol] - colLower[icol] > primalTolerance_) {
              markC[icol] &= ~3;
            } else {
              markC[icol] = 3;
            }
          }
          for (istackR = 0; istackR < nstackR; istackR++) {
            int irow = stackR[istackR];
            markR[irow] = -1;
          }
        } else {
          /* is it worth seeing if can increase coefficients
             or maybe better see if it is a cut */
          if (iway == 0) {
            nstackC0 = std::min(nstackC, maxStack);
            double solMove = saveSolval - down;
            double boundChange;
            if (notFeasible) {
              nstackC0 = 0;
            } else {
              for (istackC = 0; istackC < nstackC0; istackC++) {
                int icol = stackC[istackC];
                stackC0[istackC] = icol;
                lo0[istackC] = colLower[icol];
                up0[istackC] = colUpper[icol];
              }
            }
            /* restore all */
            assert(iway == 0);
            if (changedDown)
              changedDown[j] = nstackC;
            for (istackC = nstackC - 1; istackC >= 0; istackC--) {
              int icol = stackC[istackC];
              double oldU = saveU[istackC];
              double oldL = saveL[istackC];
#if FIXED_BOTH_WAYS
              if (doFixedBothWays && intVar[icol]
                && colLower[icol] != colUpper[icol]) {
                assert(newLo[icol] == COIN_DBL_MAX);
                fixedBoth[nFixedBoth++] = icol;
                newLo[icol] = colLower[icol];
                newUp[icol] = colUpper[icol];
              }
#endif
              if (goingToTrueBound == 2 && istackC && !justReplace) {
                // upper disaggregation cut would be
                // xval < upper + (old_upper-upper) (jval-down)
                boundChange = oldU - colUpper[icol];
                if (boundChange > 0.0 && oldU < 1.0e10 && (colsol[icol] > colUpper[icol] + boundChange * solMove + primalTolerance_)) {
                  // create cut
                  OsiRowCut rc;
                  rc.setLb(-COIN_DBL_MAX);
                  double newUb = colUpper[icol] - down * boundChange;
                  if (!intVar[icol] && newUb > 1.0e3)
                    newUb += 1.0e-7 * newUb;
                  rc.setUb(newUb);
                  index[0] = icol;
                  element[0] = 1.0;
                  index[1] = j;
                  element[1] = -boundChange;
                  // effectiveness is how far j moves
                  double newSol = (colsol[icol] - colUpper[icol]) / boundChange;
                  assert(newSol > solMove);
                  rc.setEffectiveness(newSol - solMove);
                  if (rc.effectiveness() > disaggEffectiveness) {
                    rc.setRow(2, index, element, false);
#ifdef CGL_DEBUG
                    if (debugger)
                      assert(!debugger->invalidCut(rc));
#endif
                    rowCut.addCutIfNotDuplicate(rc);
                  }
                }
                // lower disaggregation cut would be
                // xval > lower + (old_lower-lower) (jval-down)
                boundChange = oldL - colLower[icol];
                if (boundChange < 0.0 && oldL > -1.0e10 && (colsol[icol] < colLower[icol] + boundChange * solMove - primalTolerance_)) {
                  // create cut
                  OsiRowCut rc;
                  double newLb = colLower[icol] - down * boundChange;
                  if (!intVar[icol] && fabs(newLb) > 1.0e3)
                    newLb -= 1.0e-7 * fabs(newLb);
                  rc.setLb(newLb);
                  rc.setUb(COIN_DBL_MAX);
                  index[0] = icol;
                  element[0] = 1.0;
                  index[1] = j;
                  element[1] = -boundChange;
                  // effectiveness is how far j moves
                  double newSol = (colsol[icol] - colLower[icol]) / boundChange;
                  assert(newSol > solMove);
                  rc.setEffectiveness(newSol - solMove);
                  if (rc.effectiveness() > disaggEffectiveness) {
                    rc.setRow(2, index, element, false);
#ifdef CGL_DEBUG
                    if (debugger)
                      assert(!debugger->invalidCut(rc));
#endif
                    rowCut.addCutIfNotDuplicate(rc);
#if 0
                    printf("%d original bounds %g, %g new Lo %g sol= %g int %d sol= %g\n",icol,oldL,oldU,colLower[icol],colsol[icol], j, colsol[j]);
                    printf("-1.0 * x(%d) + %g * y(%d) <= %g\n",
                           icol,boundChange,j,rc.ub());
#endif
                  }
                }
              }
              colUpper[icol] = oldU;
              colLower[icol] = oldL;
              columnGap[icol] = oldU - oldL - primalTolerance_;
              markC[icol] = 0;
              if (oldU > 1.0e10)
                markC[icol] |= 8;
              if (oldL < -1.0e10)
                markC[icol] |= 4;
            }
            for (istackR = 0; istackR < nstackR; istackR++) {
              int irow = stackR[istackR];
              // switch off strengthening if not wanted
              if ((rowCuts & 2) != 0 && goingToTrueBound) {
                bool ifCut = anyColumnCuts;
                double gap = rowUpper[irow] - maxR[irow];
                double sum = 0.0;
                if (!ifCut && (gap > primalTolerance_ && gap < 1.0e8)) {
                  // see if the strengthened row is a cut
                  // also see if singletons can go to good objective
                  // Taken out as should be found elsewhere
                  // and has to be original column length
#ifdef MOVE_SINGLETONS
                  bool moveSingletons = (irow < moveSingletonRows);
#endif
                  for (CoinBigIndex kk = rowStart[irow]; kk < rowStart[irow + 1];
                    kk++) {
                    int iColumn = column[kk];
                    double value = rowElements[kk];
                    sum += value * colsol[iColumn];
#ifdef MOVE_SINGLETONS
                    if (moveSingletons && j != iColumn) {
                      if (colUpper[iColumn] > colLower[iColumn]) {
                        if (columnLength2[iColumn] != 1) {
                          moveSingletons = false;
                        }
                      }
                    }
#endif
                  }
#ifdef MOVE_SINGLETONS
                  if (moveSingletons) {
                    // can fix any with good costs
                    for (CoinBigIndex kk = rowStart[irow]; kk < rowStart[irow + 1];
                      kk++) {
                      int iColumn = column[kk];
                      if (j != iColumn) {
                        if (colUpper[iColumn] > colLower[iColumn]) {
                          if (columnLength2[iColumn] == 1) {
                            double value = rowElements[kk];
                            if (direction * objective[iColumn] * value < 0.0 && !(markC[iColumn] & 3)) {
                              // Fix
                              if (nstackC0 + 1 < maxStack) {
                                stackC0[nstackC0] = iColumn;
                                if (value > 0.0) {
                                  lo0[nstackC0] = colUpper[iColumn];
                                  up0[nstackC0] = colUpper[iColumn];
                                } else {
                                  lo0[nstackC0] = colLower[iColumn];
                                  up0[nstackC0] = colLower[iColumn];
                                }
                                nstackC0++;
                              }
                            }
                          }
                        }
                      }
                    }
                  }
#endif
                  if (sum - gap * colsol[j] > maxR[irow] + primalTolerance_ || (info->strengthenRow && rowLower[irow] < -1.0e20)) {
                    // can be a cut
                    // subtract gap from upper and integer coefficient
                    // saveU and saveL spare
                    int *COIN_RESTRICT index = reinterpret_cast< int * COIN_RESTRICT >(saveL);
                    double *COIN_RESTRICT element = saveU;
                    int n = 0;
                    bool coefficientExists = false;
                    double sum2 = 0.0;
                    for (CoinBigIndex kk = rowStart[irow]; kk < rowStart[irow + 1];
                      kk++) {
                      int kColumn = column[kk];
                      double el = rowElements[kk];
                      if (kColumn != j) {
                        index[n] = kColumn;
                        element[n++] = el;
                      } else {
                        el = el - gap;
                        if (fabs(el) > 1.0e-12) {
                          index[n] = kColumn;
                          element[n++] = el;
                        }
                        coefficientExists = true;
                      }
                      sum2 += colsol[kColumn] * el;
                    }
                    if (!coefficientExists) {
                      index[n] = j;
                      element[n++] = -gap;
                      sum2 -= colsol[j] * gap;
                    }
                    OsiRowCut rc;
                    rc.setLb(-COIN_DBL_MAX);
                    double ub = rowUpper[irow] - gap * (colLower[j] + 1.0);
                    rc.setUb(ub);
                    double effectiveness = sum2 - ub;
                    effectiveness = std::max(effectiveness,
                      (sum - gap * colsol[j]
                        - maxR[irow])
                        / gap);
                    if (!coefficientExists)
                      effectiveness = std::max(1.0e-7,
                        effectiveness);
                    rc.setEffectiveness(effectiveness);
                    // rc.setEffectiveness((sum-gap*colsol[j]-maxR[irow])/gap);
                    if (rc.effectiveness() > needEffectiveness) {
                      rc.setRow(n, index, element, false);
#ifdef CGL_DEBUG
                      if (debugger)
                        assert(!debugger->invalidCut(rc));
#endif
                      // If strengthenRow point to row
                      // if(info->strengthenRow)
                      // printf("a point to row %d\n",irow);
                      // #define STRENGTHEN_PRINT
#ifdef STRENGTHEN_PRINT
                      if (canReplace && rowLower[irow] < -1.0e20) {
                        printf("1Cut %g <= ", rc.lb());
                        int k;
                        // printf("original row %d - %g <= <= %g - j = %d\n",iow,rowLower[irow],rowUpper[irow],j);
                        // for (int kk=rowStart[irow];kk<rowStart[irow+1];kk++)
                        // printf("(%d,%g) ",column[kk],rowElements[kk]);
                        // printf("\n");
                        for (k = 0; k < n; k++) {
                          int iColumn = index[k];
                          printf("%g*", element[k]);
                          if (intVar[iColumn])
                            printf("i%d ", iColumn);
                          else
                            printf("x%d ", iColumn);
                        }
                        printf("<= %g\n", rc.ub());
                        printf("Row %g <= ", rowLower[irow]);
                        for (k = rowStart[irow]; k < rowStart[irow + 1]; k++) {
                          int iColumn = column[k];
                          printf("%g*", rowElements[k]);
                          if (intVar[iColumn])
                            printf("i%d ", iColumn);
                          else
                            printf("x%d ", iColumn);
                        }
                        printf("<= %g\n", rowUpper[irow]);
                      }
#endif
                      int realRow = (canReplace && rowLower[irow] < -1.0e20) ? irow : -1;
                      if (realRows && realRow >= 0)
                        realRow = realRows[realRow];
                      if (!justReplace) {
                        rowCut.addCutIfNotDuplicate(rc, realRow);
                      } else if (realRow >= 0) {
                        double effectiveness = 0.0;
                        for (int i = 0; i < n; i++)
                          effectiveness += fabs(element[i]);
                        if (!info->strengthenRow[realRow] || info->strengthenRow[realRow]->effectiveness() > effectiveness) {
                          delete info->strengthenRow[realRow];
                          rc.setEffectiveness(effectiveness);
                          info->strengthenRow[realRow] = new OsiRowCut(std::move(rc));
                        }
                      }
                    }
                  }
                }
                gap = minR[irow] - rowLower[irow];
                if (!ifCut && (gap > primalTolerance_ && gap < 1.0e8)) {
                  // see if the strengthened row is a cut
                  sum = 0.0;
                  // also see if singletons can go to good objective
#ifdef MOVE_SINGLETONS
                  bool moveSingletons = (irow < moveSingletonRows);
#endif
                  for (CoinBigIndex kk = rowStart[irow]; kk < rowStart[irow + 1];
                    kk++) {
                    int iColumn = column[kk];
                    double value = rowElements[kk];
                    sum += value * colsol[iColumn];
#ifdef MOVE_SINGLETONS
                    if (moveSingletons && j != iColumn) {
                      if (colUpper[iColumn] > colLower[iColumn]) {
                        if (columnLength2[iColumn] != 1) {
                          moveSingletons = false;
                        }
                      }
                    }
#endif
                  }
#ifdef MOVE_SINGLETONS
                  if (moveSingletons) {
                    // can fix any with good costs
                    for (CoinBigIndex kk = rowStart[irow]; kk < rowStart[irow + 1];
                      kk++) {
                      int iColumn = column[kk];
                      if (j != iColumn) {
                        if (colUpper[iColumn] > colLower[iColumn]) {
                          if (columnLength2[iColumn] == 1) {
                            double value = rowElements[kk];
                            if (direction * objective[iColumn] * value > 0.0 && !(markC[iColumn] & 3)) {
                              // Fix
                              if (nstackC0 + 1 < maxStack) {
                                stackC0[nstackC0] = iColumn;
                                if (value < 0.0) {
                                  lo0[nstackC0] = colUpper[iColumn];
                                  up0[nstackC0] = colUpper[iColumn];
                                } else {
                                  lo0[nstackC0] = colLower[iColumn];
                                  up0[nstackC0] = colLower[iColumn];
                                }
                                nstackC0++;
                              }
                            }
                          }
                        }
                      }
                    }
                  }
#endif
                  if (sum + gap * colsol[j] < minR[irow] - primalTolerance_ || (info->strengthenRow && rowUpper[irow] > 1.0e20)) {
                    // can be a cut
                    // add gap to lower and integer coefficient
                    // saveU and saveL spare
                    int *COIN_RESTRICT index = reinterpret_cast< int * COIN_RESTRICT >(saveL);
                    double *COIN_RESTRICT element = saveU;
                    int n = 0;
                    bool coefficientExists = false;
                    double sum2 = 0.0;
                    for (CoinBigIndex kk = rowStart[irow]; kk < rowStart[irow + 1];
                      kk++) {
                      int kColumn = column[kk];
                      double el = rowElements[kk];
                      if (kColumn != j) {
                        index[n] = kColumn;
                        element[n++] = el;
                      } else {
                        el = el + gap;
                        if (fabs(el) > 1.0e-12) {
                          index[n] = kColumn;
                          element[n++] = el;
                        }
                        coefficientExists = true;
                      }
                      sum2 += colsol[kColumn] * el;
                    }
                    if (!coefficientExists) {
                      index[n] = j;
                      element[n++] = gap;
                      sum2 += colsol[j] * gap;
                    }
                    OsiRowCut rc;
                    double lb = rowLower[irow] + gap * (colLower[j] + 1.0);
                    rc.setLb(lb);
                    rc.setUb(COIN_DBL_MAX);
                    // effectiveness
                    double effectiveness = lb - sum2;
                    effectiveness = std::max(effectiveness,
                      (minR[irow] - sum - gap * colsol[j]) / gap);
                    if (!coefficientExists)
                      effectiveness = std::max(1.0e-7,
                        effectiveness);
                    rc.setEffectiveness(effectiveness);
                    if (rc.effectiveness() > needEffectiveness) {
                      rc.setRow(n, index, element, false);
#ifdef CGL_DEBUG
                      if (debugger)
                        assert(!debugger->invalidCut(rc));
#endif
                      // if(info->strengthenRow)
                      // printf("b point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
                      if (canReplace && rowUpper[irow] > 1.0e20) {
                        printf("2Cut %g <= ", rc.lb());
                        int k;
                        for (k = 0; k < n; k++) {
                          int iColumn = index[k];
                          printf("%g*", element[k]);
                          if (intVar[iColumn])
                            printf("i%d ", iColumn);
                          else
                            printf("x%d ", iColumn);
                        }
                        printf("<= %g\n", rc.ub());
                        printf("Row %g <= ", rowLower[irow]);
                        for (k = rowStart[irow]; k < rowStart[irow + 1]; k++) {
                          int iColumn = column[k];
                          printf("%g*", rowElements[k]);
                          if (intVar[iColumn])
                            printf("i%d ", iColumn);
                          else
                            printf("x%d ", iColumn);
                        }
                        printf("<= %g\n", rowUpper[irow]);
                      }
#endif
                      int realRow = (canReplace && rowUpper[irow] > 1.0e20) ? irow : -1;
                      if (realRows && realRow >= 0)
                        realRow = realRows[realRow];
                      if (!justReplace) {
                        rowCut.addCutIfNotDuplicate(rc, realRow);
                      } else if (realRow >= 0) {
                        double effectiveness = 0.0;
                        for (int i = 0; i < n; i++)
                          effectiveness += fabs(element[i]);
                        if (!info->strengthenRow[realRow] || info->strengthenRow[realRow]->effectiveness() > effectiveness) {
                          delete info->strengthenRow[realRow];
                          rc.setEffectiveness(effectiveness);
                          info->strengthenRow[realRow] = new OsiRowCut(std::move(rc));
                        }
                      }
                    }
                  }
                }
              }
              minR[irow] = saveMin[istackR];
              maxR[irow] = saveMax[istackR];
              markR[irow] = -1;
            }
          } else {
            if (iway == 1 && feasible == 3) {
              if ((info->options & 2048) != 0 && false) {
                printf("way %d j %d ", iway, j);
                for (int i = 0; i < std::min(nCols, 10); i++) {
                  if (intVar[i] && colUpper[i] < 1.5)
                    printf("(%d,%d) ", (int)(colLower[i]), (int)(colUpper[i]));
                }
                printf("\n");
              }
              iway = 3;
#ifdef MOVE_SINGLETONS
              // look for singletons that can move (just at root)
              if ((rowCuts & 2) != 0 && goingToTrueBound && info->strengthenRow) {
                for (istackR = 0; istackR < nstackR; istackR++) {
                  int irow = stackR[istackR];
                  if (irow >= moveSingletonRows)
                    continue;
                  double gap = rowUpper[irow] - maxR[irow];
                  if (gap > primalTolerance_) {
                    // also see if singletons can go to good objective
                    bool moveSingletons = true;
                    for (CoinBigIndex kk = rowStart[irow]; kk < rowStart[irow + 1];
                      kk++) {
                      int iColumn = column[kk];
                      if (moveSingletons && j != iColumn) {
                        if (colUpper[iColumn] > colLower[iColumn]) {
                          if (columnLength2[iColumn] != 1) {
                            moveSingletons = false;
                          }
                        }
                      }
                    }
                    if (moveSingletons) {
                      // can fix any with good costs
                      for (CoinBigIndex kk = rowStart[irow]; kk < rowStart[irow + 1];
                        kk++) {
                        int iColumn = column[kk];
                        if (j != iColumn) {
                          if (colUpper[iColumn] > colLower[iColumn]) {
                            if (columnLength2[iColumn] == 1) {
                              double value = rowElements[kk];
                              if (direction * objective[iColumn] * value < 0.0 && !(markC[iColumn] & 3)) {
                                // Fix
                                stackC[nstackC] = iColumn;
                                saveL[nstackC] = colLower[iColumn];
                                saveU[nstackC] = colUpper[iColumn];
                                assert(saveU[nstackC] > saveL[nstackC]);
                                if (value > 0.0) {
                                  colLower[iColumn] = colUpper[iColumn];
                                } else {
                                  colUpper[iColumn] = colLower[iColumn];
                                }
                                columnGap[iColumn] = -primalTolerance_;
                                assert(nstackC < nCols);
                                nstackC++;
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                  gap = minR[irow] - rowLower[irow];
                  if (gap > primalTolerance_) {
                    // also see if singletons can go to good objective
                    bool moveSingletons = true;
                    for (CoinBigIndex kk = rowStart[irow]; kk < rowStart[irow + 1];
                      kk++) {
                      int iColumn = column[kk];
                      if (moveSingletons && j != iColumn) {
                        if (colUpper[iColumn] > colLower[iColumn]) {
                          if (columnLength2[iColumn] != 1) {
                            moveSingletons = false;
                          }
                        }
                      }
                    }
                    if (moveSingletons) {
                      // can fix any with good costs
                      for (CoinBigIndex kk = rowStart[irow]; kk < rowStart[irow + 1];
                        kk++) {
                        int iColumn = column[kk];
                        if (j != iColumn) {
                          if (colUpper[iColumn] > colLower[iColumn]) {
                            if (columnLength2[iColumn] == 1) {
                              double value = rowElements[kk];
                              if (direction * objective[iColumn] * value > 0.0 && !(markC[iColumn] & 3)) {
                                // Fix
                                stackC[nstackC] = iColumn;
                                saveL[nstackC] = colLower[iColumn];
                                saveU[nstackC] = colUpper[iColumn];
                                assert(saveU[nstackC] > saveL[nstackC]);
                                if (value < 0.0) {
                                  colLower[iColumn] = colUpper[iColumn];
                                } else {
                                  colUpper[iColumn] = colLower[iColumn];
                                }
                                columnGap[iColumn] = -primalTolerance_;
                                assert(nstackC < nCols);
                                nstackC++;
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
#endif
              /* point back to stack */
              for (istackC = nstackC - 1; istackC >= 0; istackC--) {
                int icol = stackC[istackC];
                markC[icol] = istackC + 100000;
              }
              OsiColCut *cc = new OsiColCut;
              int /*nTot=0,*/ nFix = 0, nInt = 0;
              bool ifCut = false;
              for (istackC = 1; istackC < nstackC0; istackC++) {
                int icol = stackC0[istackC];
                int istackC1 = markC[icol] - 100000;
                if (istackC1 >= 0) {
                  markC[icol] -= 100000;
                  if (std::min(lo0[istackC], colLower[icol]) > saveL[istackC1] + 1.0e-4) {
                    double saveBound = saveL[istackC1];
                    saveL[istackC1] = std::min(lo0[istackC], colLower[icol]);
#ifdef UPDATE_MINR_MAXR
                    if (saveBound < saveL[istackC1]) {
                      // tighten minR and maxR
                      assert(numberRealChanges < 2 * maxStack);
                      realChange[numberRealChanges] = saveL[istackC1] - saveBound;
                      whichRealChange[numberRealChanges++] = icol;
                    }
#endif
                    if (intVar[icol] /*||!info->inTree*/) {
                      element[nFix] = saveL[istackC1];
                      index[nFix++] = icol;
                      nInt++;
                      if (colsol[icol] < saveL[istackC1] - primalTolerance_)
                        ifCut = true;
                    }
                    nfixed++;
                  } else {
                    // try up
                    // markC[icol] += 100000;
                  }
                }
              }
              if (nFix) {
                // nTot=nFix;
                cleanBounds(-1, nFix, index, element);
                cc->setLbs(nFix, index, element);
                nFix = 0;
              }
              for (istackC = 1; istackC < nstackC0; istackC++) {
                int icol = stackC0[istackC];
                int istackC1 = markC[icol] - 100000;
                if (istackC1 >= 0) {
                  markC[icol] -= 100000;
                  if (std::max(up0[istackC], colUpper[icol]) < saveU[istackC1] - 1.0e-4) {
                    double saveBound = saveU[istackC1];
                    saveU[istackC1] = std::max(up0[istackC], colUpper[icol]);
#ifdef UPDATE_MINR_MAXR
                    if (saveBound < saveU[istackC1]) {
                      // tighten minR and maxR
                      assert(numberRealChanges < 2 * maxStack);
                      realChange[numberRealChanges] = saveBound - saveU[istackC1];
                      whichRealChange[numberRealChanges++] = icol;
                    }
#endif
                    if (intVar[icol] /*||!info->inTree*/) {
                      element[nFix] = saveU[istackC1];
                      index[nFix++] = icol;
                      nInt++;
                      if (colsol[icol] > saveU[istackC1] + primalTolerance_)
                        ifCut = true;
                    }
                    nfixed++;
                  } else if (!info->inTree && saveL[0] == 0.0 && saveU[0] == 1.0) {
                    // See if can do two cut
                    double upperWhenDown = up0[istackC];
                    double lowerWhenDown = lo0[istackC];
                    double upperWhenUp = colUpper[icol];
                    double lowerWhenUp = colLower[icol];
                    double upperOriginal = saveU[istackC1];
                    double lowerOriginal = saveL[istackC1];
                    if (upperWhenDown < lowerOriginal + 1.0e-12 && lowerWhenUp > upperOriginal - 1.0e-12) {
                      OsiRowCut *rc = new OsiRowCut;
                      rc->setLb(lowerOriginal);
                      rc->setUb(lowerOriginal);
                      rc->setEffectiveness(1.0e-5);
                      int index[2];
                      double element[2];
                      index[0] = j;
                      index[1] = icol;
                      element[0] = -(upperOriginal - lowerOriginal);
                      // If zero then - must have been fixed without noticing!
                      if (fabs(element[0]) > 1.0e-8) {
                        element[1] = 1.0;
                        rc->setRow(2, index, element, false);
                        cs.insert(rc);
                      } else {
                        delete rc;
                      }
                    } else if (upperWhenUp < lowerOriginal + 1.0e-12 && lowerWhenDown > upperOriginal - 1.0e-12) {
                      OsiRowCut *rc = new OsiRowCut;
                      rc->setLb(upperOriginal);
                      rc->setUb(upperOriginal);
                      rc->setEffectiveness(1.0e-5);
                      int index[2];
                      double element[2];
                      index[0] = j;
                      index[1] = icol;
                      element[0] = upperOriginal - lowerOriginal;
                      element[1] = 1.0;
                      rc->setRow(2, index, element, false);
                      cs.insert(rc);
                    }
                  }
                }
              }
              if (nFix) {
                // nTot+=nFix;
                cleanBounds(1, nFix, index, element);
                cc->setUbs(nFix, index, element);
              }
              // could tighten continuous as well
              if (nInt) {
                if (ifCut) {
                  cc->setEffectiveness(100.0);
                } else {
                  cc->setEffectiveness(1.0e-5);
                }
#ifdef CGL_DEBUG
                checkBounds(debugger, *cc);
#endif
                cs.insert(cc);
              } else {
                delete cc;
              }
            } else {
              goingToTrueBound = 0;
            }
            double solMove = up - saveSolval;
            double boundChange;
            /* restore all */
            assert(iway > 0);
            if (changedUp)
              changedUp[j] = nstackC;
#if FIXED_BOTH_WAYS
            // #define JUST_ONE_CC
#ifdef JUST_ONE_CC
            CoinPackedVector lbsBoth;
            CoinPackedVector ubsBoth;
#endif
#endif
            for (istackC = nstackC - 1; istackC > 0; istackC--) {
              int icol = stackC[istackC];
              double oldU = saveU[istackC];
              double oldL = saveL[istackC];
#if FIXED_BOTH_WAYS
              if (doFixedBothWays && newLo[icol] != COIN_DBL_MAX) {
                double loDown = newLo[icol];
                double upDown = newUp[icol];
                double loUp = colLower[icol];
                double upUp = colUpper[icol];
                double newLower = oldL;
                double newUpper = oldU;
#if FIXED_BOTH_WAYS > 2
                bool printit = false;
#endif
                if (loDown > oldL + 1.0e-5 && loUp > oldL + 1.0e-5) {
                  newLower = std::min(loDown, loUp);
                  if (intVar[icol])
                    newLower = ceil(newLower - 1.0e-5);
#if FIXED_BOTH_WAYS > 2
                  nToBoundZ++;
#endif
#ifdef JUST_ONE_CC
                  lbsBoth.insert(icol, newLower);
#endif
#if FIXED_BOTH_WAYS > 2
                  printit = true;
                  printf("x%d >= %g || ", icol, newLower);
#endif
                  if (upDown < oldU - 1.0e-5 && upUp < oldU - 1.0e-5) {
                    newUpper = std::max(upDown, upUp);
#if FIXED_BOTH_WAYS > 2
                    nToBoundZ++;
#endif
                    if (intVar[icol])
                      newUpper = floor(newUpper + 1.0e-5);
#if FIXED_BOTH_WAYS > 2
                    printf("<= %g || ", newUpper);
#endif
                  }
                } else if (upDown < oldU - 1.0e-5 && upUp < oldU - 1.0e-5) {
                  newUpper = std::max(upDown, upUp);
                  if (intVar[icol])
                    newUpper = floor(newUpper + 1.0e-5);
#if FIXED_BOTH_WAYS > 2
                  nToBoundZ++;
#endif
#ifdef JUST_ONE_CC
                  ubsBoth.insert(icol, newUpper);
#endif
#if FIXED_BOTH_WAYS > 2
                  printit = true;
                  printf("x%d <= %g || ", icol, newUpper);
#endif
                } else if (loDown == upDown && loUp == upUp) {
                  if (saveL[0] == 0.0 && saveU[0] == 1.0) {
#if FIXED_BOTH_WAYS > 1
                    nCutZ++;
                    printf("Cutzz x%d == %g x%d + %g || ", icol, (loUp - loDown),
                      j, loDown);
                    printit = true;
#endif
                    OsiRowCut *rc = new OsiRowCut;
                    rc->setLb(loDown);
                    rc->setUb(loDown);
                    rc->setEffectiveness(1.0e-5);
                    int index[2];
                    double element[2];
                    index[0] = j;
                    index[1] = icol;
                    element[0] = loDown - loUp;
                    element[1] = 1.0;
                    rc->setRow(2, index, element, false);
                    cs.insert(rc);
                  }
                }
#if FIXED_BOTH_WAYS > 2
                if (printit)
                  printf("Col %d bnds on %d down %g,%g up %g,%g orig %g,%g\n",
                    j, icol, loDown, upDown, loUp,
                    upUp, oldL, oldU);
#endif
#ifdef UPDATE_MINR_MAXR
                if (oldL < newLower) {
                  // tighten minR and maxR
                  assert(numberRealChanges < 2 * maxStack);
                  realChange[numberRealChanges] = newLower - oldL;
                  whichRealChange[numberRealChanges++] = icol;
                }
                if (oldU > newUpper) {
                  // tighten minR and maxR
                  assert(numberRealChanges < 2 * maxStack);
                  realChange[numberRealChanges] = newUpper - oldU;
                  whichRealChange[numberRealChanges++] = icol;
                }
#ifndef JUST_ONE_CC
                // temp to check
                if (oldL < newLower || oldU > newUpper) {
                  OsiColCut *cc = new OsiColCut;
                  CoinPackedVector lbs;
                  CoinPackedVector ubs;
                  lbs.insert(icol, newLower);
                  ubs.insert(icol, newUpper);
                  cc->setLbs(lbs);
                  cc->setUbs(ubs);
                  cs.insert(cc);
                }
#endif
#endif
                oldL = newLower;
                oldU = newUpper;
              }
#endif
              if (goingToTrueBound == 2 && !justReplace) {
                // upper disaggregation cut would be
                // xval < upper + (old_upper-upper) (up-jval)
                boundChange = oldU - colUpper[icol];
                if (boundChange > 0.0 && oldU < 1.0e10 && (colsol[icol] > colUpper[icol] + boundChange * solMove + primalTolerance_)) {
                  // create cut
                  OsiRowCut rc;
                  rc.setLb(-COIN_DBL_MAX);
                  double newUb = colUpper[icol] + up * boundChange;
                  if (!intVar[icol] && newUb > 1.0e3)
                    newUb += 1.0e-7 * newUb;
                  rc.setUb(newUb);
                  index[0] = icol;
                  element[0] = 1.0;
                  index[1] = j;
                  element[1] = +boundChange;
                  // effectiveness is how far j moves
                  double newSol = (colsol[icol] - colUpper[icol]) / boundChange;
                  assert(newSol > solMove);
                  rc.setEffectiveness(newSol - solMove);
                  if (rc.effectiveness() > disaggEffectiveness) {
                    rc.setRow(2, index, element, false);
#ifdef CGL_DEBUG
                    if (debugger)
                      assert(!debugger->invalidCut(rc));
#endif
                    rowCut.addCutIfNotDuplicate(rc);
                  }
                }
                // lower disaggregation cut would be
                // xval > lower + (old_lower-lower) (up-jval)
                boundChange = oldL - colLower[icol];
                if (boundChange < 0.0 && oldL > -1.0e10 && (colsol[icol] < colLower[icol] + boundChange * solMove - primalTolerance_)) {
                  // create cut
                  OsiRowCut rc;
                  double newLb = colLower[icol] + up * boundChange;
                  if (!intVar[icol] && fabs(newLb) > 1.0e3)
                    newLb -= 1.0e-7 * fabs(newLb);
                  rc.setLb(newLb);
                  rc.setUb(COIN_DBL_MAX);
                  index[0] = icol;
                  element[0] = 1.0;
                  index[1] = j;
                  element[1] = +boundChange;
                  // effectiveness is how far j moves
                  double newSol = (colsol[icol] - colLower[icol]) / boundChange;
                  assert(newSol > solMove);
                  rc.setEffectiveness(newSol - solMove);
                  if (rc.effectiveness() > disaggEffectiveness) {
                    rc.setRow(2, index, element, false);
#ifdef CGL_DEBUG
                    if (debugger)
                      assert(!debugger->invalidCut(rc));
#endif
                    rowCut.addCutIfNotDuplicate(rc);
                  }
                }
              }
              colUpper[icol] = oldU;
              colLower[icol] = oldL;
              columnGap[icol] = oldU - oldL - primalTolerance_;
              if (oldU > oldL + 1.0e-4) {
                markC[icol] = 0;
                if (oldU > 1.0e10)
                  markC[icol] |= 8;
                if (oldL < -1.0e10)
                  markC[icol] |= 4;
              } else {
                markC[icol] = 3;
              }
            }
#if FIXED_BOTH_WAYS
#ifdef JUST_ONE_CC
            if (lbsBoth.getNumElements() || ubsBoth.getNumElements()) {
              OsiColCut *cc = new OsiColCut;
              if (lbsBoth.getNumElements())
                cc->setLbs(lbsBoth);
              if (ubsBoth.getNumElements())
                cc->setUbs(ubsBoth);
              cs.insert(cc);
            }
#endif
#endif
            // first
            {
              double oldU = saveU[0];
              double oldL = saveL[0];
              colUpper[j] = oldU;
              colLower[j] = oldL;
              columnGap[j] = oldU - oldL - primalTolerance_;
              assert(oldU > oldL + 1.0e-4);
              markC[j] = 0;
              if (oldU > 1.0e10)
                markC[j] |= 8;
              if (oldL < -1.0e10)
                markC[j] |= 4;
            }
            for (istackR = 0; istackR < nstackR; istackR++) {
              int irow = stackR[istackR];
              // switch off strengthening if not wanted
              if ((rowCuts & 2) != 0 && goingToTrueBound) {
                bool canReplace = info->strengthenRow && (goingToTrueBound == 2);
                bool ifCut = anyColumnCuts;
                double gap = rowUpper[irow] - maxR[irow];
                double sum = 0.0;
                if (!ifCut && (gap > primalTolerance_ && gap < 1.0e8)) {
                  // see if the strengthened row is a cut
                  for (CoinBigIndex kk = rowStart[irow]; kk < rowStart[irow + 1];
                    kk++) {
                    sum += rowElements[kk] * colsol[column[kk]];
                  }
                  if (sum + gap * colsol[j] > rowUpper[irow] + primalTolerance_ || (canReplace && rowLower[irow] < -1.e20)) {
                    // can be a cut
                    // add gap to integer coefficient
                    // saveU and saveL spare
                    int *COIN_RESTRICT index = reinterpret_cast< int * COIN_RESTRICT >(saveL);
                    double *COIN_RESTRICT element = saveU;
                    int n = 0;
                    bool coefficientExists = false;
                    double sum2 = 0.0;
                    for (CoinBigIndex kk = rowStart[irow]; kk < rowStart[irow + 1];
                      kk++) {
                      int kColumn = column[kk];
                      double el = rowElements[kk];
                      if (kColumn != j) {
                        index[n] = kColumn;
                        element[n++] = el;
                      } else {
                        el = el + gap;
                        if (fabs(el) > 1.0e-12) {
                          index[n] = kColumn;
                          element[n++] = el;
                        }
                        coefficientExists = true;
                      }
                      sum2 += colsol[kColumn] * el;
                    }
                    if (!coefficientExists) {
                      index[n] = j;
                      element[n++] = gap;
                      sum2 += colsol[j] * gap;
                    }
                    OsiRowCut rc;
                    rc.setLb(-COIN_DBL_MAX);
                    double ub = rowUpper[irow] + gap * (colUpper[j] - 1.0);
                    rc.setUb(ub);
                    // effectiveness
                    double effectiveness = sum2 - ub;
                    effectiveness = std::max(effectiveness,
                      (sum + gap * colsol[j] - rowUpper[irow]) / gap);
                    if (!coefficientExists)
                      effectiveness = std::max(1.0e-7,
                        effectiveness);
                    rc.setEffectiveness(effectiveness);
                    if (rc.effectiveness() > needEffectiveness) {
                      rc.setRow(n, index, element, false);
#ifdef CGL_DEBUG
                      if (debugger)
                        assert(!debugger->invalidCut(rc));
#endif
                      // if(canReplace)
                      // printf("c point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
                      if (canReplace && rowLower[irow] < -1.0e20) {
                        printf("3Cut %g <= ", rc.lb());
                        int k;
                        for (k = 0; k < n; k++) {
                          int iColumn = index[k];
                          printf("%g*", element[k]);
                          if (intVar[iColumn])
                            printf("i%d ", iColumn);
                          else
                            printf("x%d ", iColumn);
                        }
                        printf("<= %g\n", rc.ub());
                        printf("Row %g <= ", rowLower[irow]);
                        for (k = rowStart[irow]; k < rowStart[irow + 1]; k++) {
                          int iColumn = column[k];
                          printf("%g*", rowElements[k]);
                          if (intVar[iColumn])
                            printf("i%d ", iColumn);
                          else
                            printf("x%d ", iColumn);
                        }
                        printf("<= %g\n", rowUpper[irow]);
                      }
#endif
                      int realRow = (canReplace && rowLower[irow] < -1.0e20) ? irow : -1;
                      if (realRows && realRow >= 0)
                        realRow = realRows[realRow];
                      if (!justReplace) {
                        rowCut.addCutIfNotDuplicate(rc, realRow);
                      } else if (realRow >= 0) {
                        double effectiveness = 0.0;
                        for (int i = 0; i < n; i++)
                          effectiveness += fabs(element[i]);
                        if (!info->strengthenRow[realRow] || info->strengthenRow[realRow]->effectiveness() > effectiveness) {
                          delete info->strengthenRow[realRow];
                          rc.setEffectiveness(effectiveness);
                          info->strengthenRow[realRow] = new OsiRowCut(std::move(rc));
                        }
                      }
                    }
                  }
                }
                gap = minR[irow] - rowLower[irow];
                if (!ifCut && (gap > primalTolerance_ && gap < 1.0e8)) {
                  // see if the strengthened row is a cut
                  if (!sum) {
                    for (CoinBigIndex kk = rowStart[irow]; kk < rowStart[irow + 1];
                      kk++) {
                      sum += rowElements[kk] * colsol[column[kk]];
                    }
                  }
                  if (sum - gap * colsol[j] < rowLower[irow] - primalTolerance_ || (canReplace && rowUpper[irow] > 1.0e20)) {
                    // can be a cut
                    // subtract gap from integer coefficient
                    // saveU and saveL spare
                    int *COIN_RESTRICT index = reinterpret_cast< int * COIN_RESTRICT >(saveL);
                    double *COIN_RESTRICT element = saveU;
                    int n = 0;
                    bool coefficientExists = false;
                    double sum2 = 0.0;
                    for (CoinBigIndex kk = rowStart[irow]; kk < rowStart[irow + 1];
                      kk++) {
                      int kColumn = column[kk];
                      double el = rowElements[kk];
                      if (kColumn != j) {
                        index[n] = kColumn;
                        element[n++] = el;
                      } else {
                        el = el - gap;
                        if (fabs(el) > 1.0e-12) {
                          index[n] = kColumn;
                          element[n++] = el;
                        }
                        coefficientExists = true;
                      }
                      sum2 += colsol[kColumn] * el;
                    }
                    if (!coefficientExists) {
                      index[n] = j;
                      element[n++] = -gap;
                      sum2 -= colsol[j] * gap;
                    }
                    OsiRowCut rc;
                    double lb = rowLower[irow] - gap * (colUpper[j] - 1);
                    rc.setLb(lb);
                    rc.setUb(COIN_DBL_MAX);
                    double effectiveness = lb - sum2;
                    effectiveness = std::max(effectiveness,
                      (rowLower[irow] - sum + gap * colsol[j]) / gap);
                    if (!coefficientExists)
                      effectiveness = std::max(1.0e-7,
                        effectiveness);
                    rc.setEffectiveness(effectiveness);
                    if (rc.effectiveness() > needEffectiveness) {
                      rc.setRow(n, index, element, false);
#ifdef CGL_DEBUG
                      if (debugger) {
                        // assert(!debugger->invalidCut(rc));
                        if (debugger->invalidCut(rc)) {
                          printf("bad cut\n");
                          abort();
                        }
                      }
#endif
                      // if(canReplace)
                      // printf("d point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
                      if (canReplace && rowUpper[irow] > 1.0e20) {
                        printf("4Cut %g <= ", rc.lb());
                        int k;
                        for (k = 0; k < n; k++) {
                          int iColumn = index[k];
                          printf("%g*", element[k]);
                          if (intVar[iColumn])
                            printf("i%d ", iColumn);
                          else
                            printf("x%d ", iColumn);
                        }
                        printf("<= %g\n", rc.ub());
                        printf("Row %g <= ", rowLower[irow]);
                        for (k = rowStart[irow]; k < rowStart[irow + 1]; k++) {
                          int iColumn = column[k];
                          printf("%g*", rowElements[k]);
                          if (intVar[iColumn])
                            printf("i%d ", iColumn);
                          else
                            printf("x%d ", iColumn);
                        }
                        printf("<= %g\n", rowUpper[irow]);
                      }
#endif
                      int realRow = (canReplace && rowUpper[irow] > 1.0e20) ? irow : -1;
                      if (realRows && realRow >= 0)
                        realRow = realRows[realRow];
                      if (!justReplace) {
                        rowCut.addCutIfNotDuplicate(rc, realRow);
                      } else if (realRow >= 0) {
                        double effectiveness = 0.0;
                        for (int i = 0; i < n; i++)
                          effectiveness += fabs(element[i]);
                        if (!info->strengthenRow[realRow] || info->strengthenRow[realRow]->effectiveness() > effectiveness) {
                          delete info->strengthenRow[realRow];
                          rc.setEffectiveness(effectiveness);
                          info->strengthenRow[realRow] = new OsiRowCut(std::move(rc));
                        }
                      }
                    }
                  }
                }
              }
              minR[irow] = saveMin[istackR];
              maxR[irow] = saveMax[istackR];
              markR[irow] = -1;
            }
          }
        }
      }
    }
  }
  if (!ninfeas && (info->strengthenRow || !rowCut.numberCuts()) && rowCuts) {
    // Try and find ALL big M's
    for (int i = 0; i < nRowsSafe; ++i) {
      if ((rowLower[i] > -1.0e20 || rowUpper[i] < 1.0e20) && (!info->strengthenRow || !info->strengthenRow[i])) {
        int iflagu = 0;
        int iflagl = 0;
        double dmaxup = 0.0;
        double dmaxdown = 0.0;
        CoinBigIndex krs = rowStart[i];
        CoinBigIndex kre = rowStart[i + 1];
        int kInt = -1;
        double rhsAdjustment = 0.0;
        // int nPosInt=0;
        // int nNegInt=0;
        double valueInteger = 0.0;
        // Find largest integer coefficient
        CoinBigIndex k;
        for (k = krs; k < kre; ++k) {
          int j = column[k];
          if (intVar[j]) {
            double value = rowElements[k];
            if (colUpper[j] > colLower[j] && !colLower[j] && fabs(value) > fabs(valueInteger)) {
              kInt = j;
              valueInteger = value;
            }
          }
        }
        if (kInt >= 0) {
          double upperBound = std::min(colUpper[kInt], static_cast< double >(COIN_INT_MAX));
          double upAdjust = 0.0;
          double downAdjust = 0.0;
          for (k = krs; k < kre; ++k) {
            double value = rowElements[k];
            int j = column[k];
            if (colUpper[j] == colLower[j]) {
              rhsAdjustment += colUpper[j] * value;
              continue;
            }
            //	    if (intVar[j]) {
            //	      if (value>0.0)
            //		nPosInt++;
            //	      else
            //		nNegInt++;
            //	    } else {
            //	      nPosInt = -nCols;
            //	    }
            if (j != kInt) {
              // treat as continuous
              if (value > 0.0) {
                if (colUpper[j] >= 1e15) {
                  dmaxup = 1e31;
                  ++iflagu;
                } else {
                  dmaxup += colUpper[j] * value;
                }
                if (colLower[j] <= -1e15) {
                  dmaxdown = -1e31;
                  ++iflagl;
                } else {
                  dmaxdown += colLower[j] * value;
                }
              } else if (value < 0.0) {
                if (colUpper[j] >= 1e15) {
                  dmaxdown = -1e31;
                  ++iflagl;
                } else {
                  dmaxdown += colUpper[j] * value;
                }
                if (colLower[j] <= -1e15) {
                  dmaxup = 1e31;
                  ++iflagu;
                } else {
                  dmaxup += colLower[j] * value;
                }
              }
            } else {
              // Chosen variable
              if (value > 0.0) {
                if (colUpper[j] >= 1e15) {
                  upAdjust = 1e31;
                } else {
                  upAdjust = colUpper[j] * value;
                }
                if (colLower[j] <= -1e15) {
                  downAdjust = -1e31;
                } else {
                  downAdjust = colLower[j] * value;
                }
              } else if (value < 0.0) {
                if (colUpper[j] >= 1e15) {
                  downAdjust = -1e31;
                } else {
                  downAdjust = colUpper[j] * value;
                }
                if (colLower[j] <= -1e15) {
                  upAdjust = 1e31;
                } else {
                  upAdjust = colLower[j] * value;
                }
              }
            }
          }
          dmaxup += rhsAdjustment;
          dmaxdown += rhsAdjustment;
          // end of row
          if (iflagu)
            dmaxup = 1.0e31;
          if (iflagl)
            dmaxdown = -1.0e31;
          // See if redundant
          if (dmaxdown + downAdjust > rowLower[i] - tolerance && dmaxup + upAdjust < rowUpper[i] + tolerance)
            continue;
          if (dmaxdown + valueInteger * upperBound > rowLower[i] && dmaxup + valueInteger * upperBound < rowUpper[i]) {
            // check to see if always feasible at 1 but not always at 0
            if (dmaxdown + valueInteger > rowLower[i] && dmaxup + valueInteger < rowUpper[i] && (dmaxdown < rowLower[i] - primalTolerance_ || dmaxup > rowUpper[i] + primalTolerance_)) {
              // can tighten (maybe)
              double saveValue = valueInteger;
              if (valueInteger > 0.0) {
                assert(dmaxdown < rowLower[i]);
                valueInteger = rowLower[i] - dmaxdown;
              } else {
                assert(dmaxup > rowUpper[i]);
                valueInteger = rowUpper[i] - dmaxup;
              }
              if (fabs(saveValue - valueInteger) > 1.0e-12) {
                // take
                OsiRowCut rc;
                rc.setLb(rowLower[i]);
                rc.setUb(rowUpper[i]);
                int n = 0;
                double sum = 0.0;
                for (CoinBigIndex kk = rowStart[i]; kk < rowStart[i + 1]; kk++) {
                  int j = column[kk];
                  if (j != kInt) {
                    sum += colsol[j] * rowElements[kk];
                    index[n] = j;
                    element[n++] = rowElements[kk];
                  } else {
                    sum += colsol[j] * valueInteger;
                    assert(rowElements[kk] * valueInteger >= 0.0);
#if 0
                    if (fabs(rowElements[kk])>1.01*fabs(valueInteger)) {
                      printf("row %d changing coefficient of %d from %g to %g\n",
                             i,kInt,rowElements[kk],valueInteger);
                    }
#endif
                    if (fabs(valueInteger) > 1.0e-12) {
                      index[n] = column[kk];
                      element[n++] = valueInteger;
                    }
                  }
                }
                double gap = 0.0;
                if (sum < rowLower[i])
                  gap = rowLower[i] - sum;
                else if (sum > rowUpper[i])
                  gap = sum - rowUpper[i];
                if (gap > 1.0e-4 || info->strengthenRow != NULL) {
                  gap += 1.0e5;
                  rc.setEffectiveness(gap);
                  rc.setRow(n, index, element, false);
#ifdef STRENGTHEN_PRINT
                  {
                    printf("1aCut %g <= ", rc.lb());
                    int irow = i;
                    int k;
                    for (k = 0; k < n; k++) {
                      int iColumn = index[k];
                      printf("%g*", element[k]);
                      if (intVar[iColumn])
                        printf("i%d ", iColumn);
                      else
                        printf("x%d ", iColumn);
                    }
                    printf("<= %g\n", rc.ub());
                    printf("Row %g <= ", rowLower[irow]);
                    for (k = rowStart[irow]; k < rowStart[irow + 1]; k++) {
                      int iColumn = column[k];
                      printf("%g*", rowElements[k]);
                      if (intVar[iColumn])
                        printf("i%d ", iColumn);
                      else
                        printf("x%d ", iColumn);
                    }
                    printf("<= %g\n", rowUpper[irow]);
                  }
#endif
                  int returnCode = rowCut.addCutIfNotDuplicate(rc, i);
                  if (returnCode < 0)
                    break; // out of space
                }
              }
            }
          }
        }
      }
    }
  }
#ifndef ONE_ARRAY
  delete[] stackC0;
  delete[] lo0;
  delete[] up0;
  delete[] columnGap;
  delete[] markC;
  delete[] stackC;
  delete[] stackR;
  delete[] saveL;
  delete[] saveU;
  delete[] saveMin;
  delete[] saveMax;
  delete[] index;
  delete[] element;
  delete[] djs;
  delete[] largestPositiveInRow;
  delete[] largestNegativeInRow;
#endif
#if FIXED_BOTH_WAYS || MANY_TIGHTEN
  delete[] saveFColLower;
#endif
  delete[] colsol;
  // Add in row cuts
  if (!ninfeas && (info->options & 16384) == 0) {
    if (!justReplace) {
      rowCut.addCuts(cs, info->strengthenRow, info->pass);
    } else {
      for (int i = 0; i < nRows; i++) {
        int realRow = realRows[i];
        if (realRow >= 0) {
          OsiRowCut *cut = info->strengthenRow[realRow];
          if (cut) {
#ifdef CLP_INVESTIGATE
            printf("Row %d, real row %d effectiveness %g\n", i, realRow, cut->effectiveness());
#endif
            cs.insert(cut);
          }
        }
      }
    }
  }
#if 0
  {
    int numberRowCutsAfter = cs.sizeRowCuts() ;
    int k ;
    for (k = 0;k<numberRowCutsAfter;k++) {
      OsiRowCut thisCut = cs.rowCut(k) ;
      printf("Cut %d is %g <=",k,thisCut.lb());
      int n=thisCut.row().getNumElements();
      const int *COIN_RESTRICT column = thisCut.row().getIndices();
      const double *COIN_RESTRICT element = thisCut.row().getElements();
      assert (n);
      for (int i=0;i<n;i++) {
	printf(" %g*x%d",element[i],column[i]);
      }
      printf(" <= %g\n",thisCut.ub());
    }
  }
#endif
  if (rowCutFake) {
    numberCliques_ = rowCutFake->numberCuts();
    CoinBigIndex numberElements = 2 * numberCliques_;
    cliqueEntry_ = new CliqueEntry[numberElements];
    int *COIN_RESTRICT column = reinterpret_cast< int * COIN_RESTRICT >(cliqueEntry_);
    numberElements = 0;
    for (int iCut = 0; iCut < numberCliques_; iCut++) {
      column[numberElements++] = rowCutFake->cut(iCut)->row().getIndices()[0];
      column[numberElements++] = rowCutFake->cut(iCut)->row().getIndices()[1];
    }
  }
  delete rowCutFake;
#if FIXED_BOTH_WAYS > 1
  if (nToBoundZ || nCutZ)
    printf("zz %d fixed to bounds and %d == cuts\n",
      nToBoundZ, nCutZ);
#endif
#if CHECK_MINR_MAXR > 0
  printf("end of run due to CHECK_MINR_MAXR in CglProbing\n");
  exit(99);
#endif
  return (ninfeas);
}
// Does probing and adding cuts
int CglProbing::probeCliques(const OsiSolverInterface &si,
  const OsiRowCutDebugger *
#ifdef CGL_DEBUG
    debugger
#endif
  ,
  OsiCuts &cs,
  double *colLower, double *colUpper,
  CoinPackedMatrix *rowCopy,
  CoinPackedMatrix *columnCopy, const int *realRows,
  double *rowLower, double *rowUpper,
  char *intVar, double *minR, double *maxR,
  int *markR,
  CglTreeInfo *info)
{
  // Set up maxes
  int maxStack = info->inTree ? maxStack_ : maxStackRoot_;
  int nRows = rowCopy->getNumRows();
  int nCols = rowCopy->getNumCols();
  double *colsol = new double[nCols];
  double *djs = new double[nCols];
  const double *currentColLower = si.getColLower();
  const double *currentColUpper = si.getColUpper();
  double *tempL = new double[nCols];
  double *tempU = new double[nCols];
  int *markC = new int[nCols];
  int *stackC = new int[2 * nCols];
  int *stackR = new int[nRows];
  double *saveL = new double[2 * nCols];
  double *saveU = new double[2 * nCols];
  double *saveMin = new double[nRows];
  double *saveMax = new double[nRows];
  double *element = new double[nCols];
  int *index = new int[nCols];
  // For trying to extend cliques
  int *cliqueStack = NULL;
  int *cliqueCount = NULL;
  int *to_01 = NULL;
  int *cliqueAdd = NULL;
  int maxCliqueAdded = nCols + nRows;
  int numberCliqueAdded = 0;
  int *cliqueAdd2 = NULL;
  if (!mode_) {
    cliqueAdd = new int[2 * maxCliqueAdded];
    cliqueAdd2 = cliqueAdd + maxCliqueAdded;
    to_01 = new int[nCols];
    cliqueStack = new int[numberCliques_];
    cliqueCount = new int[numberCliques_];
    int i;
    for (i = 0; i < numberCliques_; i++) {
      cliqueCount[i] = static_cast< int >(cliqueStart_[i + 1] - cliqueStart_[i]);
    }
    for (i = 0; i < nCols; i++)
      to_01[i] = -1;
    for (i = 0; i < number01Integers_; i++) {
      int j = cutVector_[i].sequence;
      to_01[j] = i;
    }
  }
  // Let us never add more than twice the number of rows worth of row cuts
  // Keep cuts out of cs until end so we can find duplicates quickly
  int nRowsFake = info->inTree ? nRows / 3 : nRows;
  row_cut rowCut(nRowsFake, !info->inTree);
  const int *column = rowCopy->getIndices();
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
  const int *rowLength = rowCopy->getVectorLengths();
  const double *rowElements = rowCopy->getElements();
  const int *row = columnCopy->getIndices();
  const CoinBigIndex *columnStart = columnCopy->getVectorStarts();
  const int *columnLength = columnCopy->getVectorLengths();
  const double *columnElements = columnCopy->getElements();
  double movement;
  int i, j;
  CoinBigIndex jj, k, kk;
  int kcol, krow;
  bool anyColumnCuts = false;
  double dbound, value, value2;
  int ninfeas = 0;
  int rowCuts;
  double disaggEffectiveness;
  if (mode_) {
    /* clean up djs and solution */
    CoinMemcpyN(si.getReducedCost(), nCols, djs);
    CoinMemcpyN(si.getColSolution(), nCols, colsol);
    disaggEffectiveness = 1.0e-3;
    rowCuts = rowCuts_;
  } else {
    // need to go from a neutral place
    memset(djs, 0, nCols * sizeof(double));
    CoinMemcpyN(si.getColSolution(), nCols, colsol);
    disaggEffectiveness = -1.0e10;
    if (rowCuts_ != 4)
      rowCuts = 1;
    else
      rowCuts = 4;
  }
  for (i = 0; i < nCols; ++i) {
    /* was if (intVar[i]) */
    if (1) {
      if (colUpper[i] - colLower[i] > 1.0e-8) {
        if (colsol[i] < colLower[i] + primalTolerance_) {
          colsol[i] = colLower[i];
          djs[i] = std::max(0.0, djs[i]);
        } else if (colsol[i] > colUpper[i] - primalTolerance_) {
          colsol[i] = colUpper[i];
          djs[i] = std::min(0.0, djs[i]);
        } else {
          djs[i] = 0.0;
        }
        /*if (fabs(djs[i])<1.0e-5)
          djs[i]=0.0;*/
      }
    }
  }

  int ipass = 0, nfixed = -1;

  double cutoff;
  bool cutoff_available = si.getDblParam(OsiDualObjectiveLimit, cutoff);
  if (!cutoff_available || usingObjective_ < 0) { // cut off isn't set or isn't valid
    cutoff = si.getInfinity();
  }
  cutoff *= si.getObjSense();
  if (fabs(cutoff) > 1.0e30)
    assert(cutoff > 1.0e30);
  double current = si.getObjValue();
  // make irrelevant if mode is 0
  if (!mode_)
    cutoff = COIN_DBL_MAX;
  /* for both way coding */
  int nstackC0 = -1;
  int *stackC0 = new int[maxStack];
  double *lo0 = new double[maxStack];
  double *up0 = new double[maxStack];
  int nstackR, nstackC;
  for (i = 0; i < nCols; i++) {
    if (colUpper[i] - colLower[i] < 1.0e-8) {
      markC[i] = 3;
    } else {
      markC[i] = 0;
    }
  }
  double tolerance = 1.0e1 * primalTolerance_;
  int maxPass = info->inTree ? maxPass_ : maxPassRoot_;
  // If we are going to replace coefficient then we don't need to be effective
  double needEffectiveness = info->strengthenRow ? -1.0e10 : 1.0e-3;
  const double probingDeadline = (maxSeconds_ > 0.0) ? (CoinWallclockTime() + maxSeconds_) : 0.0;
  while (ipass < maxPass && nfixed) {
    int iLook;
    if (probingDeadline && CoinWallclockTime() > probingDeadline)
      break;
    ipass++;
    nfixed = 0;
    for (iLook = 0; iLook < numberThisTime_; iLook++) {
      double solval;
      if (probingDeadline && (iLook & 127) == 0 && CoinWallclockTime() > probingDeadline) {
        ipass = maxPass;
        break;
      }
      double down;
      double up;
      j = lookedAt_[iLook];
      solval = colsol[j];
      down = floor(solval + tolerance);
      up = ceil(solval - tolerance);
      if (colUpper[j] - colLower[j] < 1.0e-8)
        markC[j] = 3;
      if (markC[j] || !intVar[j])
        continue;
      double saveSolval = solval;
      if (solval >= colUpper[j] - tolerance || solval <= colLower[j] + tolerance || up == down) {
        if (solval <= colLower[j] + 2.0 * tolerance) {
          solval = colLower[j] + 1.0e-1;
          down = colLower[j];
          up = down + 1.0;
        } else if (solval >= colUpper[j] - 2.0 * tolerance) {
          solval = colUpper[j] - 1.0e-1;
          up = colUpper[j];
          down = up - 1;
        } else {
          // odd
          up = down + 1.0;
          solval = down + 1.0e-1;
        }
      }
      assert(up <= colUpper[j]);
      assert(down >= colLower[j]);
      assert(up > down);
      if ((solval - down > 1.0e-6 && up - solval > 1.0e-6) || mode_ != 1) {
        int istackC, iway, istackR;
        int way[] = { 1, 2, 1 };
        int feas[] = { 1, 2, 4 };
        int feasible = 0;
        int notFeasible;
        for (iway = 0; iway < 3; iway++) {
          // int fixThis=0;
          double objVal = current;
          int goingToTrueBound = 0;
          stackC[0] = j;
          markC[j] = way[iway];
          double solMovement;
          if (way[iway] == 1) {
            movement = down - colUpper[j];
            solMovement = down - colsol[j];
            assert(movement < -0.99999);
            if (fabs(down - colLower[j]) < 1.0e-7) {
              goingToTrueBound = 2;
              down = colLower[j];
            }
          } else {
            movement = up - colLower[j];
            solMovement = up - colsol[j];
            assert(movement > 0.99999);
            if (fabs(up - colUpper[j]) < 1.0e-7) {
              goingToTrueBound = 2;
              up = colUpper[j];
            }
          }
          if (goingToTrueBound && (colUpper[j] - colLower[j] > 1.5 || colLower[j]))
            goingToTrueBound = 1;
          // switch off disaggregation if not wanted
          if ((rowCuts & 1) == 0)
            goingToTrueBound = 0;
#ifdef PRINT_DEBUG
          if (fabs(movement) > 1.01) {
            printf("big %d %g %g %g\n", j, colLower[j], solval, colUpper[j]);
          }
#endif
          if (solMovement * djs[j] > 0.0)
            objVal += solMovement * djs[j];
          nstackC = 1;
          nstackR = 0;
          saveL[0] = colLower[j];
          saveU[0] = colUpper[j];
          assert(saveU[0] > saveL[0]);
          notFeasible = 0;
          if (movement < 0.0) {
            colUpper[j] += movement;
            colUpper[j] = floor(colUpper[j] + 0.5);
#ifdef PRINT_DEBUG
            printf("** Trying %d down to 0\n", j);
#endif
          } else {
            colLower[j] += movement;
            colLower[j] = floor(colLower[j] + 0.5);
#ifdef PRINT_DEBUG
            printf("** Trying %d up to 1\n", j);
#endif
          }
          if (fabs(colUpper[j] - colLower[j]) < 1.0e-6)
            markC[j] = 3; // say fixed
          istackC = 0;
          /* update immediately */
          for (k = columnStart[j]; k < columnStart[j] + columnLength[j]; k++) {
            int irow = row[k];
            value = columnElements[k];
            assert(markR[irow] != -2);
            if (markR[irow] == -1) {
              stackR[nstackR] = irow;
              markR[irow] = nstackR;
              saveMin[nstackR] = minR[irow];
              saveMax[nstackR] = maxR[irow];
              nstackR++;
#if 0
	    } else if (markR[irow]==-2) {
	      continue;
#endif
            }
            /* could check immediately if violation */
            if (movement > 0.0) {
              /* up */
              if (value > 0.0) {
                /* up does not change - down does */
                if (minR[irow] > -1.0e10)
                  minR[irow] += value;
                if (minR[irow] > rowUpper[irow] + 1.0e-5) {
                  notFeasible = 1;
                  istackC = 1;
                  break;
                }
              } else {
                /* down does not change - up does */
                if (maxR[irow] < 1.0e10)
                  maxR[irow] += value;
                if (maxR[irow] < rowLower[irow] - 1.0e-5) {
                  notFeasible = 1;
                  istackC = 1;
                  break;
                }
              }
            } else {
              /* down */
              if (value < 0.0) {
                /* up does not change - down does */
                if (minR[irow] > -1.0e10)
                  minR[irow] -= value;
                if (minR[irow] > rowUpper[irow] + 1.0e-5) {
                  notFeasible = 1;
                  istackC = 1;
                  break;
                }
              } else {
                /* down does not change - up does */
                if (maxR[irow] < 1.0e10)
                  maxR[irow] -= value;
                if (maxR[irow] < rowLower[irow] - 1.0e-5) {
                  notFeasible = 1;
                  istackC = 1;
                  break;
                }
              }
            }
          }
          while (istackC < nstackC && nstackC < maxStack) {
            int jway;
            int jcol = stackC[istackC];
            jway = markC[jcol];
            // If not first and fixed then skip
            if (jway == 3 && istackC) {
              // istackC++;
              // continue;
              // printf("fixed %d on stack\n",jcol);
            }
            // Do cliques
            if (oneFixStart_ && oneFixStart_[jcol] >= 0) {
              int start;
              int end;
              if (colLower[jcol] > saveL[istackC]) {
                // going up
                start = oneFixStart_[jcol];
                end = zeroFixStart_[jcol];
              } else {
                assert(colUpper[jcol] < saveU[istackC]);
                // going down
                start = zeroFixStart_[jcol];
                end = endFixStart_[jcol];
              }
              for (int i = start; i < end; i++) {
                int iClique = whichClique_[i];
                for (CoinBigIndex k = cliqueStart_[iClique]; k < cliqueStart_[iClique + 1]; k++) {
                  int kcol = sequenceInCliqueEntry(cliqueEntry_[k]);
                  if (jcol == kcol)
                    continue;
                  int kway = oneFixesInCliqueEntry(cliqueEntry_[k]);
                  if (kcol != jcol) {
                    if (!markC[kcol]) {
                      // not on list yet
                      if (nstackC < 2 * maxStack) {
                        markC[kcol] = 3; // say fixed
                        // fixThis++;
                        stackC[nstackC] = kcol;
                        saveL[nstackC] = colLower[kcol];
                        saveU[nstackC] = colUpper[kcol];
                        assert(saveU[nstackC] > saveL[nstackC]);
                        nstackC++;
                        if (!kway) {
                          // going up
                          double solMovement = 1.0 - colsol[kcol];
                          if (solMovement > 0.0001) {
                            assert(djs[kcol] >= 0.0);
                            objVal += djs[kcol] * solMovement;
                          }
                          colLower[kcol] = 1.0;
                          /* update immediately */
                          for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                            krow = row[jj];
                            value = columnElements[jj];
                            assert(markR[krow] != -2);
                            if (markR[krow] == -1) {
                              stackR[nstackR] = krow;
                              markR[krow] = nstackR;
                              saveMin[nstackR] = minR[krow];
                              saveMax[nstackR] = maxR[krow];
                              nstackR++;
#if 0
                            } else if (markR[krow]==-2) {
                              continue;
#endif
                            }
                            /* could check immediately if violation */
                            /* up */
                            if (value > 0.0) {
                              /* up does not change - down does */
                              if (minR[krow] > -1.0e10)
                                minR[krow] += value;
                              if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                                colUpper[kcol] = -1.0e50; /* force infeasible */
                                break;
                              }
                            } else {
                              /* down does not change - up does */
                              if (maxR[krow] < 1.0e10)
                                maxR[krow] += value;
                              if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                                notFeasible = 1;
                                break;
                              }
                            }
                          }
                        } else {
                          // going down
                          double solMovement = 0.0 - colsol[kcol];
                          if (solMovement < -0.0001) {
                            assert(djs[kcol] <= 0.0);
                            objVal += djs[kcol] * solMovement;
                          }
                          colUpper[kcol] = 0.0;
                          /* update immediately */
                          for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                            krow = row[jj];
                            value = columnElements[jj];
                            assert(markR[krow] != -2);
                            if (markR[krow] == -1) {
                              stackR[nstackR] = krow;
                              markR[krow] = nstackR;
                              saveMin[nstackR] = minR[krow];
                              saveMax[nstackR] = maxR[krow];
                              nstackR++;
#if 0
                            } else if (markR[krow]==-2) {
                              continue;
#endif
                            }
                            /* could check immediately if violation */
                            /* down */
                            if (value < 0.0) {
                              /* up does not change - down does */
                              if (minR[krow] > -1.0e10)
                                minR[krow] -= value;
                              if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                                notFeasible = 1;
                                break;
                              }
                            } else {
                              /* down does not change - up does */
                              if (maxR[krow] < 1.0e10)
                                maxR[krow] -= value;
                              if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                                notFeasible = 1;
                                break;
                              }
                            }
                          }
                        }
                      }
                    } else if (markC[kcol] == 1) {
                      // marked as going to 0
                      assert(!colUpper[kcol]);
                      if (!kway) {
                        // contradiction
                        notFeasible = 1;
                        break;
                      }
                    } else if (markC[kcol] == 2) {
                      // marked as going to 1
                      assert(colLower[kcol]);
                      if (kway) {
                        // contradiction
                        notFeasible = 1;
                        break;
                      }
                    } else {
                      // marked as fixed
                      assert(markC[kcol] == 3);
                      int jkway;
                      if (colLower[kcol])
                        jkway = 1;
                      else
                        jkway = 0;
                      if (kway == jkway) {
                        // contradiction
                        notFeasible = 1;
                        break;
                      }
                    }
                  }
                }
                if (notFeasible)
                  break;
              }
              if (notFeasible)
                istackC = nstackC + 1;
            }
            for (k = columnStart[jcol]; k < columnStart[jcol] + columnLength[jcol]; k++) {
              // break if found not feasible
              if (notFeasible)
                break;
              int irow = row[k];
              /*value = columnElements[k];*/
              assert(markR[irow] != -2);
#if 0
	      if (markR[irow]!=-2) {
#endif
              /* see if anything forced */
              for (kk = rowStart[irow]; kk < rowStart[irow] + rowLength[irow]; kk++) {
                double moveUp = 0.0;
                double moveDown = 0.0;
                double newUpper = -1.0, newLower = 1.0;
                kcol = column[kk];
                bool onList = (markC[kcol] != 0);
                if (markC[kcol] != 3) {
                  value2 = rowElements[kk];
                  int markIt = markC[kcol];
                  if (value2 < 0.0) {
                    if (colUpper[kcol] < 1e10 && (markIt & 2) == 0 && rowUpper[irow] < 1.0e10) {
                      dbound = colUpper[kcol] + (rowUpper[irow] - minR[irow]) / value2;
                      if (dbound > colLower[kcol] + primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 2;
                          newLower = ceil(dbound - primalTolerance_);
                        } else {
                          newLower = dbound;
                          if (newLower + primalTolerance_ > colUpper[kcol] && newLower - primalTolerance_ <= colUpper[kcol]) {
                            newLower = colUpper[kcol];
                            markIt |= 2;
                            // markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt = 3;
                          }
                        }
                        moveUp = newLower - colLower[kcol];
                      }
                    }
                    if (colLower[kcol] > -1e10 && (markIt & 1) == 0 && rowLower[irow] > -1.0e10) {
                      dbound = colLower[kcol] + (rowLower[irow] - maxR[irow]) / value2;
                      if (dbound < colUpper[kcol] - primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 1;
                          newUpper = floor(dbound + primalTolerance_);
                        } else {
                          newUpper = dbound;
                          if (newUpper - primalTolerance_ < colLower[kcol] && newUpper + primalTolerance_ >= colLower[kcol]) {
                            newUpper = colLower[kcol];
                            markIt |= 1;
                            // markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt = 3;
                          }
                        }
                        moveDown = newUpper - colUpper[kcol];
                      }
                    }
                  } else {
                    /* positive element */
                    if (colUpper[kcol] < 1e10 && (markIt & 2) == 0 && rowLower[irow] > -1.0e10) {
                      dbound = colUpper[kcol] + (rowLower[irow] - maxR[irow]) / value2;
                      if (dbound > colLower[kcol] + primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 2;
                          newLower = ceil(dbound - primalTolerance_);
                        } else {
                          newLower = dbound;
                          if (newLower + primalTolerance_ > colUpper[kcol] && newLower - primalTolerance_ <= colUpper[kcol]) {
                            newLower = colUpper[kcol];
                            markIt |= 2;
                            // markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt = 3;
                          }
                        }
                        moveUp = newLower - colLower[kcol];
                      }
                    }
                    if (colLower[kcol] > -1e10 && (markIt & 1) == 0 && rowUpper[irow] < 1.0e10) {
                      dbound = colLower[kcol] + (rowUpper[irow] - minR[irow]) / value2;
                      if (dbound < colUpper[kcol] - primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 1;
                          newUpper = floor(dbound + primalTolerance_);
                        } else {
                          newUpper = dbound;
                          if (newUpper - primalTolerance_ < colLower[kcol] && newUpper + primalTolerance_ >= colLower[kcol]) {
                            newUpper = colLower[kcol];
                            markIt |= 1;
                            // markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt = 3;
                          }
                        }
                        moveDown = newUpper - colUpper[kcol];
                      }
                    }
                  }
                  if (nstackC < 2 * maxStack) {
                    markC[kcol] = markIt;
                  }
                  if (moveUp && nstackC < 2 * maxStack) {
                    // fixThis++;
#ifdef PRINT_DEBUG
                    printf("lower bound on %d increased from %g to %g by row %d %g %g\n", kcol, colLower[kcol], newLower, irow, rowLower[irow], rowUpper[irow]);
                    value = 0.0;
                    for (jj = rowStart[irow]; jj < rowStart[irow] + rowLength[irow]; jj++) {
                      int ii = column[jj];
                      if (colUpper[ii] - colLower[ii] > primalTolerance_) {
                        printf("(%d, %g) ", ii, rowElements[jj]);
                      } else {
                        value += rowElements[jj] * colLower[ii];
                      }
                    }
                    printf(" - fixed %g\n", value);
                    for (jj = rowStart[irow]; jj < rowStart[irow] + rowLength[irow]; jj++) {
                      int ii = column[jj];
                      if (colUpper[ii] - colLower[ii] < primalTolerance_) {
                        printf("(%d, %g, %g) ", ii, rowElements[jj], colLower[ii]);
                      }
                    }
                    printf("\n");
#endif
                    if (!onList) {
                      stackC[nstackC] = kcol;
                      saveL[nstackC] = colLower[kcol];
                      saveU[nstackC] = colUpper[kcol];
                      assert(saveU[nstackC] > saveL[nstackC]);
                      nstackC++;
                      onList = true;
                    }
                    if (newLower > colsol[kcol]) {
                      if (djs[kcol] < 0.0) {
                        /* should be infeasible if integer */
                        if (intVar[kcol])
                          assert(newLower > colUpper[kcol] + primalTolerance_);
                      } else {
                        objVal += moveUp * djs[kcol];
                      }
                    }
                    if (intVar[kcol])
                      newLower = std::max(colLower[kcol], ceil(newLower - 1.0e-4));
                    colLower[kcol] = newLower;
                    if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                      markC[kcol] = 3; // say fixed
                    }
                    /* update immediately */
                    for (jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                      krow = row[jj];
                      value = columnElements[jj];
                      assert(markR[krow] != -2);
                      if (markR[krow] == -1) {
                        stackR[nstackR] = krow;
                        markR[krow] = nstackR;
                        saveMin[nstackR] = minR[krow];
                        saveMax[nstackR] = maxR[krow];
                        nstackR++;
#if 0
			} else if (markR[krow]==-2) {
			  continue;
#endif
                      }
                      /* could check immediately if violation */
                      /* up */
                      if (value > 0.0) {
                        /* up does not change - down does */
                        if (minR[krow] > -1.0e10)
                          minR[krow] += value * moveUp;
                        if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                          colUpper[kcol] = -1.0e50; /* force infeasible */
                          break;
                        }
                      } else {
                        /* down does not change - up does */
                        if (maxR[krow] < 1.0e10)
                          maxR[krow] += value * moveUp;
                        if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                          colUpper[kcol] = -1.0e50; /* force infeasible */
                          break;
                        }
                      }
                    }
                  }
                  if (moveDown && nstackC < 2 * maxStack) {
                    // fixThis++;
#ifdef PRINT_DEBUG
                    printf("upper bound on %d decreased from %g to %g by row %d %g %g\n", kcol, colUpper[kcol], newUpper, irow, rowLower[irow], rowUpper[irow]);
                    value = 0.0;
                    for (jj = rowStart[irow]; jj < rowStart[irow] + rowLength[irow]; jj++) {
                      int ii = column[jj];
                      if (colUpper[ii] - colLower[ii] > primalTolerance_) {
                        printf("(%d, %g) ", ii, rowElements[jj]);
                      } else {
                        value += rowElements[jj] * colLower[ii];
                      }
                    }
                    printf(" - fixed %g\n", value);
                    for (jj = rowStart[irow]; jj < rowStart[irow] + rowLength[irow]; jj++) {
                      int ii = column[jj];
                      if (colUpper[ii] - colLower[ii] < primalTolerance_) {
                        printf("(%d, %g, %g) ", ii, rowElements[jj], colLower[ii]);
                      }
                    }
                    printf("\n");
#endif
                    if (!onList) {
                      stackC[nstackC] = kcol;
                      saveL[nstackC] = colLower[kcol];
                      saveU[nstackC] = colUpper[kcol];
                      assert(saveU[nstackC] > saveL[nstackC]);
                      nstackC++;
                      onList = true;
                    }
                    if (newUpper < colsol[kcol]) {
                      if (djs[kcol] > 0.0) {
                        /* should be infeasible if integer */
                        if (intVar[kcol])
                          assert(colLower[kcol] > newUpper + primalTolerance_);
                      } else {
                        objVal += moveDown * djs[kcol];
                      }
                    }
                    if (intVar[kcol])
                      newUpper = std::min(colUpper[kcol], floor(newUpper + 1.0e-4));
                    colUpper[kcol] = newUpper;
                    if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                      markC[kcol] = 3; // say fixed
                    }
                    /* update immediately */
                    for (jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                      krow = row[jj];
                      value = columnElements[jj];
                      assert(markR[krow] != -2);
                      if (markR[krow] == -1) {
                        stackR[nstackR] = krow;
                        markR[krow] = nstackR;
                        saveMin[nstackR] = minR[krow];
                        saveMax[nstackR] = maxR[krow];
                        nstackR++;
#if 0
			} else if (markR[krow]==-2) {
#endif
                        continue;
                      }
                      /* could check immediately if violation */
                      /* down */
                      if (value < 0.0) {
                        /* up does not change - down does */
                        if (minR[krow] > -1.0e10)
                          minR[krow] += value * moveDown;
                        if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                          colUpper[kcol] = -1.0e50; /* force infeasible */
                          break;
                        }
                      } else {
                        /* down does not change - up does */
                        if (maxR[krow] < 1.0e10)
                          maxR[krow] += value * moveDown;
                        if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                          colUpper[kcol] = -1.0e50; /* force infeasible */
                          break;
                        }
                      }
                    }
                  }
                  if (colLower[kcol] > colUpper[kcol] + primalTolerance_) {
                    notFeasible = 1;
                    ;
                    k = columnStart[jcol] + columnLength[jcol];
                    istackC = nstackC + 1;
#ifdef PRINT_DEBUG
                    printf("** not feasible this way\n");
#endif
                    break;
                  }
                }
              }
#if 0
	      }
#endif
            }
            istackC++;
          }
          if (!notFeasible) {
            if (objVal <= cutoff) {
              feasible |= feas[iway];
            } else {
#ifdef PRINT_DEBUG
              printf("not feasible on dj\n");
#endif
              notFeasible = 1;
              if (iway == 1 && feasible == 0) {
                /* not feasible at all */
                ninfeas = 1;
                j = nCols - 1;
                break;
              }
            }
          } else if (iway == 1 && feasible == 0) {
            /* not feasible at all */
            ninfeas = 1;
            j = nCols - 1;
            iLook = numberThisTime_;
            ipass = maxPass;
            break;
          }
          if (notFeasible)
            goingToTrueBound = 0;
          if (iway == 2 || (iway == 1 && feasible == 2)) {
            /* keep */
            iway = 3;
            nfixed++;
            if (mode_) {
              OsiColCut cc;
              int /*nTot=0,*/ nFix = 0, nInt = 0;
              bool ifCut = false;
              for (istackC = 0; istackC < nstackC; istackC++) {
                int icol = stackC[istackC];
                if (intVar[icol]) {
                  if (colUpper[icol] < currentColUpper[icol] - 1.0e-4) {
                    element[nFix] = colUpper[icol];
                    index[nFix++] = icol;
                    nInt++;
                    if (colsol[icol] > colUpper[icol] + primalTolerance_) {
                      ifCut = true;
                      anyColumnCuts = true;
                    }
                  }
                }
              }
              if (nFix) {
                // nTot=nFix;
                cleanBounds(1, nFix, index, element);
                cc.setUbs(nFix, index, element);
                nFix = 0;
              }
              for (istackC = 0; istackC < nstackC; istackC++) {
                int icol = stackC[istackC];
                if (intVar[icol]) {
                  if (colLower[icol] > currentColLower[icol] + 1.0e-4) {
                    element[nFix] = colLower[icol];
                    index[nFix++] = icol;
                    nInt++;
                    if (colsol[icol] < colLower[icol] - primalTolerance_) {
                      ifCut = true;
                      anyColumnCuts = true;
                    }
                  }
                }
              }
              if (nFix) {
                // nTot+=nFix;
                cleanBounds(-1, nFix, index, element);
                cc.setLbs(nFix, index, element);
              }
              // could tighten continuous as well
              if (nInt) {
                if (ifCut) {
                  cc.setEffectiveness(100.0);
                } else {
                  cc.setEffectiveness(1.0e-5);
                }
#ifdef CGL_DEBUG
                checkBounds(debugger, cc);
#endif
                cs.insert(cc);
              }
            }
            for (istackC = 0; istackC < nstackC; istackC++) {
              int icol = stackC[istackC];
              if (colUpper[icol] - colLower[icol] > primalTolerance_) {
                markC[icol] = 0;
              } else {
                markC[icol] = 3;
              }
            }
            for (istackR = 0; istackR < nstackR; istackR++) {
              int irow = stackR[istackR];
              markR[irow] = -1;
            }
          } else {
            /* is it worth seeing if can increase coefficients
               or maybe better see if it is a cut */
            if (iway == 0) {
              nstackC0 = std::min(nstackC, maxStack);
              double solMove = saveSolval - down;
              double boundChange;
              if (notFeasible) {
                nstackC0 = 0;
              } else {
                for (istackC = 0; istackC < nstackC0; istackC++) {
                  int icol = stackC[istackC];
                  stackC0[istackC] = icol;
                  lo0[istackC] = colLower[icol];
                  up0[istackC] = colUpper[icol];
                }
              }
              /* restore all */
              int nCliquesAffected = 0;
              assert(iway == 0);
              for (istackC = nstackC - 1; istackC >= 0; istackC--) {
                int icol = stackC[istackC];
                double oldU = saveU[istackC];
                double oldL = saveL[istackC];
                if (goingToTrueBound == 2 && istackC) {
                  // Work for extending cliques
                  if (!mode_ && numberCliques_) {
                    int i_01 = to_01[icol];
                    if (i_01 >= 0) {
                      int start;
                      int end;
                      if (colLower[icol]) {
                        // going up - but we want weak way
                        start = zeroFixStart_[icol];
                        end = endFixStart_[icol];
                      } else {
                        // going down - but we want weak way
                        start = oneFixStart_[icol];
                        end = zeroFixStart_[icol];
                      }
                      // if (end>start)
                      // printf("j %d, other %d is in %d cliques\n",
                      //      j,i_01,end-start);
                      for (int i = start; i < end; i++) {
                        int iClique = whichClique_[i];
                        int size = static_cast< int >(cliqueStart_[iClique + 1] - cliqueStart_[iClique]);
                        if (cliqueCount[iClique] == size) {
                          // first time
                          cliqueStack[nCliquesAffected++] = iClique;
                        }
                        // decrement counts
                        cliqueCount[iClique]--;
                      }
                    }
                  }
                  // upper disaggregation cut would be
                  // xval < upper + (old_upper-upper) (jval-down)
                  boundChange = oldU - colUpper[icol];
                  if (boundChange > 0.0 && oldU < 1.0e10 && (!mode_ || colsol[icol] > colUpper[icol] + boundChange * solMove + primalTolerance_)) {
                    // create cut
                    OsiRowCut rc;
                    rc.setLb(-COIN_DBL_MAX);
                    double newUb = colUpper[icol] - down * boundChange;
                    if (!intVar[icol] && newUb > 1.0e3)
                      newUb += 1.0e-7 * newUb;
                    rc.setUb(newUb);
                    index[0] = icol;
                    element[0] = 1.0;
                    index[1] = j;
                    element[1] = -boundChange;
                    // effectiveness is how far j moves
                    double newSol = (colsol[icol] - colUpper[icol]) / boundChange;
                    if (mode_)
                      assert(newSol > solMove);
                    rc.setEffectiveness(newSol - solMove);
                    if (rc.effectiveness() > disaggEffectiveness) {
                      rc.setRow(2, index, element, false);
#ifdef CGL_DEBUG
                      if (debugger)
                        assert(!debugger->invalidCut(rc));
#endif
                      rowCut.addCutIfNotDuplicate(rc);
                    }
                  }
                  // lower disaggregation cut would be
                  // xval > lower + (old_lower-lower) (jval-down)
                  boundChange = oldL - colLower[icol];
                  if (boundChange < 0.0 && oldL > -1.0e10 && (!mode_ || colsol[icol] < colLower[icol] + boundChange * solMove - primalTolerance_)) {
                    // create cut
                    OsiRowCut rc;
                    double newLb = colLower[icol] - down * boundChange;
                    if (!intVar[icol] && fabs(newLb) > 1.0e3)
                      newLb -= 1.0e-7 * fabs(newLb);
                    rc.setLb(newLb);
                    rc.setUb(COIN_DBL_MAX);
                    index[0] = icol;
                    element[0] = 1.0;
                    index[1] = j;
                    element[1] = -boundChange;
                    // effectiveness is how far j moves
                    double newSol = (colsol[icol] - colLower[icol]) / boundChange;
                    if (mode_)
                      assert(newSol > solMove);
                    rc.setEffectiveness(newSol - solMove);
                    if (rc.effectiveness() > disaggEffectiveness) {
                      rc.setRow(2, index, element, false);
#ifdef CGL_DEBUG
                      if (debugger)
                        assert(!debugger->invalidCut(rc));
#endif
                      rowCut.addCutIfNotDuplicate(rc);
#if 0
		      printf("%d original bounds %g, %g new Lo %g sol= %g int %d sol= %g\n",icol,oldL,oldU,colLower[icol],colsol[icol], j, colsol[j]);
		      printf("-1.0 * x(%d) + %g * y(%d) <= %g\n",
			     icol,boundChange,j,rc.ub());
#endif
                    }
                  }
                }
                colUpper[icol] = oldU;
                colLower[icol] = oldL;
                markC[icol] = 0;
              }
              if (nCliquesAffected) {
                for (int i = 0; i < nCliquesAffected; i++) {
                  int iClique = cliqueStack[i];
                  int size = cliqueCount[iClique];
                  // restore
                  cliqueCount[iClique] = static_cast< int >(cliqueStart_[iClique + 1] - cliqueStart_[iClique]);
                  if (!size) {
                    if (numberCliqueAdded < maxCliqueAdded) {
                      printf("Can add %d (going to 0) to clique %d (%d entries)\n",
                        j, iClique, cliqueCount[iClique]);
                      cliqueAdd[numberCliqueAdded] = j;
                      cliqueAdd2[numberCliqueAdded] = iClique;
                      numberCliqueAdded++;
                    }
                    if (logLevel_ > 1)
                      printf("** could extend clique by adding j!\n");
                  }
                }
              }
              for (istackR = 0; istackR < nstackR; istackR++) {
                int irow = stackR[istackR];
                // switch off strengthening if not wanted
                if ((rowCuts & 2) != 0 && goingToTrueBound) {
                  bool ifCut = anyColumnCuts;
                  double gap = rowUpper[irow] - maxR[irow];
                  double sum = 0.0;
                  if (!ifCut && (gap > primalTolerance_ && gap < 1.0e8)) {
                    // see if the strengthened row is a cut
                    for (kk = rowStart[irow]; kk < rowStart[irow] + rowLength[irow];
                      kk++) {
                      sum += rowElements[kk] * colsol[column[kk]];
                    }
                    if (sum - gap * colsol[j] > maxR[irow] + primalTolerance_ || (info->strengthenRow && rowLower[irow] < -1.0e20)) {
                      // can be a cut
                      // subtract gap from upper and integer coefficient
                      // saveU and saveL spare
                      int *index = reinterpret_cast< int * >(saveL);
                      double *element = saveU;
                      int n = 0;
                      bool coefficientExists = false;
                      for (kk = rowStart[irow]; kk < rowStart[irow] + rowLength[irow];
                        kk++) {
                        if (column[kk] != j) {
                          index[n] = column[kk];
                          element[n++] = rowElements[kk];
                        } else {
                          double value = rowElements[kk] - gap;
                          if (fabs(value) > 1.0e-12) {
                            index[n] = column[kk];
                            element[n++] = value;
                          }
                          coefficientExists = true;
                        }
                      }
                      if (!coefficientExists) {
                        index[n] = j;
                        element[n++] = -gap;
                      }
                      OsiRowCut rc;
                      rc.setLb(-COIN_DBL_MAX);
                      rc.setUb(rowUpper[irow] - gap * (colLower[j] + 1.0));
                      // effectiveness is how far j moves
                      rc.setEffectiveness((sum - gap * colsol[j] - maxR[irow]) / gap);
                      if (rc.effectiveness() > needEffectiveness) {
                        rc.setRow(n, index, element, false);
#ifdef CGL_DEBUG
                        if (debugger)
                          assert(!debugger->invalidCut(rc));
#endif
                        // If strengthenRow point to row
                        // if(info->strengthenRow)
                        // printf("a point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
                        if (rowLower[irow] < -1.0e20) {
                          printf("5Cut %g <= ", rc.lb());
                          int k;
                          for (k = 0; k < n; k++) {
                            int iColumn = index[k];
                            printf("%g*", element[k]);
                            if (intVar[iColumn])
                              printf("i%d ", iColumn);
                            else
                              printf("x%d ", iColumn);
                          }
                          printf("<= %g\n", rc.ub());
                          printf("Row %g <= ", rowLower[irow]);
                          for (k = rowStart[irow]; k < rowStart[irow] + rowLength[irow]; k++) {
                            int iColumn = column[k];
                            printf("%g*", rowElements[k]);
                            if (intVar[iColumn])
                              printf("i%d ", iColumn);
                            else
                              printf("x%d ", iColumn);
                          }
                          printf("<= %g\n", rowUpper[irow]);
                        }
#endif
                        int realRow = (rowLower[irow] < -1.0e20) ? irow : -1;
                        if (realRows && realRow > 0)
                          realRow = realRows[realRow];
                        rowCut.addCutIfNotDuplicate(rc, realRow);
                      }
                    }
                  }
                  gap = minR[irow] - rowLower[irow];
                  if (!ifCut && (gap > primalTolerance_ && gap < 1.0e8)) {
                    // see if the strengthened row is a cut
                    if (!sum) {
                      for (kk = rowStart[irow]; kk < rowStart[irow] + rowLength[irow];
                        kk++) {
                        sum += rowElements[kk] * colsol[column[kk]];
                      }
                    }
                    if (sum + gap * colsol[j] < minR[irow] + primalTolerance_ || (info->strengthenRow && rowUpper[irow] > 1.0e20)) {
                      // can be a cut
                      // add gap to lower and integer coefficient
                      // saveU and saveL spare
                      int *index = reinterpret_cast< int * >(saveL);
                      double *element = saveU;
                      int n = 0;
                      bool coefficientExists = false;
                      for (kk = rowStart[irow]; kk < rowStart[irow] + rowLength[irow];
                        kk++) {
                        if (column[kk] != j) {
                          index[n] = column[kk];
                          element[n++] = rowElements[kk];
                        } else {
                          double value = rowElements[kk] + gap;
                          if (fabs(value) > 1.0e-12) {
                            index[n] = column[kk];
                            element[n++] = value;
                          }
                          coefficientExists = true;
                        }
                      }
                      if (!coefficientExists) {
                        index[n] = j;
                        element[n++] = gap;
                      }
                      OsiRowCut rc;
                      rc.setLb(rowLower[irow] + gap * (colLower[j] + 1.0));
                      rc.setUb(COIN_DBL_MAX);
                      // effectiveness is how far j moves
                      rc.setEffectiveness((minR[irow] - sum - gap * colsol[j]) / gap);
                      if (rc.effectiveness() > needEffectiveness) {
                        rc.setRow(n, index, element, false);
#ifdef CGL_DEBUG
                        if (debugger)
                          assert(!debugger->invalidCut(rc));
#endif
                        // if(info->strengthenRow)
                        // printf("b point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
                        if (rowUpper[irow] > 1.0e20) {
                          printf("6Cut %g <= ", rc.lb());
                          int k;
                          for (k = 0; k < n; k++) {
                            int iColumn = index[k];
                            printf("%g*", element[k]);
                            if (intVar[iColumn])
                              printf("i%d ", iColumn);
                            else
                              printf("x%d ", iColumn);
                          }
                          printf("<= %g\n", rc.ub());
                          printf("Row %g <= ", rowLower[irow]);
                          for (k = rowStart[irow]; k < rowStart[irow] + rowLength[irow]; k++) {
                            int iColumn = column[k];
                            printf("%g*", rowElements[k]);
                            if (intVar[iColumn])
                              printf("i%d ", iColumn);
                            else
                              printf("x%d ", iColumn);
                          }
                          printf("<= %g\n", rowUpper[irow]);
                        }
#endif
                        int realRow = (rowUpper[irow] > 1.0e20) ? irow : -1;
                        if (realRows && realRow > 0)
                          realRow = realRows[realRow];
                        rowCut.addCutIfNotDuplicate(rc, realRow);
                      }
                    }
                  }
                }
                minR[irow] = saveMin[istackR];
                maxR[irow] = saveMax[istackR];
                markR[irow] = -1;
              }
            } else {
              if (iway == 1 && feasible == 3) {
                iway = 3;
                /* point back to stack */
                for (istackC = nstackC - 1; istackC >= 0; istackC--) {
                  int icol = stackC[istackC];
                  markC[icol] = istackC + 100000;
                }
                if (mode_) {
                  OsiColCut cc;
                  int /*nTot=0,*/ nFix = 0, nInt = 0;
                  bool ifCut = false;
                  for (istackC = 0; istackC < nstackC0; istackC++) {
                    int icol = stackC0[istackC];
                    int istackC1 = markC[icol] - 100000;
                    if (istackC1 >= 0) {
                      markC[icol] -= 100000;
                      if (std::min(lo0[istackC], colLower[icol]) > saveL[istackC1] + 1.0e-4) {
                        saveL[istackC1] = std::min(lo0[istackC], colLower[icol]);
                        if (intVar[icol]) {
                          element[nFix] = saveL[istackC1];
                          index[nFix++] = icol;
                          nInt++;
                          if (colsol[icol] < saveL[istackC1] - primalTolerance_)
                            ifCut = true;
                        }
                        nfixed++;
                      }
                    }
                  }
                  if (nFix) {
                    // nTot=nFix;
                    cleanBounds(-1, nFix, index, element);
                    cc.setLbs(nFix, index, element);
                    nFix = 0;
                  }
                  for (istackC = 0; istackC < nstackC0; istackC++) {
                    int icol = stackC0[istackC];
                    int istackC1 = markC[icol] - 100000;
                    if (istackC1 >= 0) {
                      markC[icol] -= 100000;
                      if (std::max(up0[istackC], colUpper[icol]) < saveU[istackC1] - 1.0e-4) {
                        saveU[istackC1] = std::max(up0[istackC], colUpper[icol]);
                        if (intVar[icol]) {
                          element[nFix] = saveU[istackC1];
                          index[nFix++] = icol;
                          nInt++;
                          if (colsol[icol] > saveU[istackC1] + primalTolerance_)
                            ifCut = true;
                        }
                        nfixed++;
                      }
                    }
                  }
                  if (nFix) {
                    // nTot+=nFix;
                    cleanBounds(1, nFix, index, element);
                    cc.setUbs(nFix, index, element);
                  }
                  // could tighten continuous as well
                  if (nInt) {
                    if (ifCut) {
                      cc.setEffectiveness(100.0);
                    } else {
                      cc.setEffectiveness(1.0e-5);
                    }
#ifdef CGL_DEBUG
                    checkBounds(debugger, cc);
#endif
                    cs.insert(cc);
                  }
                }
              } else {
                goingToTrueBound = 0;
              }
              double solMove = up - saveSolval;
              double boundChange;
              /* restore all */
              int nCliquesAffected = 0;
              for (istackC = nstackC - 1; istackC >= 0; istackC--) {
                int icol = stackC[istackC];
                double oldU = saveU[istackC];
                double oldL = saveL[istackC];
                if (goingToTrueBound == 2 && istackC) {
                  // Work for extending cliques
                  if (!mode_ && numberCliques_ && iway == 3) {
                    int i_01 = to_01[icol];
                    if (i_01 >= 0) {
                      int start;
                      int end;
                      if (colLower[icol]) {
                        // going up - but we want weak way
                        start = zeroFixStart_[icol];
                        end = endFixStart_[icol];
                      } else {
                        // going down - but we want weak way
                        start = oneFixStart_[icol];
                        end = zeroFixStart_[icol];
                      }
                      // if (end>start)
                      // printf("up j %d, other %d is in %d cliques\n",
                      //      j,i_01,end-start);
                      for (int i = start; i < end; i++) {
                        int iClique = whichClique_[i];
                        int size = static_cast< int >(cliqueStart_[iClique + 1] - cliqueStart_[iClique]);
                        if (cliqueCount[iClique] == size) {
                          // first time
                          cliqueStack[nCliquesAffected++] = iClique;
                        }
                        // decrement counts
                        cliqueCount[iClique]--;
                      }
                    }
                  }
                  // upper disaggregation cut would be
                  // xval < upper + (old_upper-upper) (up-jval)
                  boundChange = oldU - colUpper[icol];
                  if (boundChange > 0.0 && oldU < 1.0e10 && (!mode_ || colsol[icol] > colUpper[icol] + boundChange * solMove + primalTolerance_)) {
                    // create cut
                    OsiRowCut rc;
                    rc.setLb(-COIN_DBL_MAX);
                    double newUb = colUpper[icol] + up * boundChange;
                    if (!intVar[icol] && newUb > 1.0e3)
                      newUb += 1.0e-7 * newUb;
                    rc.setUb(newUb);
                    index[0] = icol;
                    element[0] = 1.0;
                    index[1] = j;
                    element[1] = +boundChange;
                    // effectiveness is how far j moves
                    double newSol = (colsol[icol] - colUpper[icol]) / boundChange;
                    if (mode_)
                      assert(newSol > solMove);
                    rc.setEffectiveness(newSol - solMove);
                    if (rc.effectiveness() > disaggEffectiveness) {
                      rc.setRow(2, index, element, false);
#ifdef CGL_DEBUG
                      if (debugger)
                        assert(!debugger->invalidCut(rc));
#endif
                      rowCut.addCutIfNotDuplicate(rc);
                    }
                  }
                  // lower disaggregation cut would be
                  // xval > lower + (old_lower-lower) (up-jval)
                  boundChange = oldL - colLower[icol];
                  if (boundChange < 0.0 && oldL > -1.0e10 && (!mode_ || colsol[icol] < colLower[icol] + boundChange * solMove - primalTolerance_)) {
                    // create cut
                    OsiRowCut rc;
                    double newLb = colLower[icol] + up * boundChange;
                    if (!intVar[icol] && fabs(newLb) > 1.0e3)
                      newLb -= 1.0e-7 * fabs(newLb);
                    rc.setLb(newLb);
                    rc.setUb(COIN_DBL_MAX);
                    index[0] = icol;
                    element[0] = 1.0;
                    index[1] = j;
                    element[1] = +boundChange;
                    // effectiveness is how far j moves
                    double newSol = (colsol[icol] - colLower[icol]) / boundChange;
                    if (mode_)
                      assert(newSol > solMove);
                    rc.setEffectiveness(newSol - solMove);
                    if (rc.effectiveness() > disaggEffectiveness) {
                      rc.setRow(2, index, element, false);
#ifdef CGL_DEBUG
                      if (debugger)
                        assert(!debugger->invalidCut(rc));
#endif
                      rowCut.addCutIfNotDuplicate(rc);
                    }
                  }
                }
                colUpper[icol] = oldU;
                colLower[icol] = oldL;
                if (oldU > oldL + 1.0e-4)
                  markC[icol] = 0;
                else
                  markC[icol] = 3;
              }
              if (nCliquesAffected) {
                for (int i = 0; i < nCliquesAffected; i++) {
                  int iClique = cliqueStack[i];
                  int size = cliqueCount[iClique];
                  // restore
                  cliqueCount[iClique] = static_cast< int >(cliqueStart_[iClique + 1] - cliqueStart_[iClique]);
                  if (!size) {
                    if (numberCliqueAdded < maxCliqueAdded) {
                      printf("Can add %d (going to 1) to clique %d (%d entries)\n",
                        j, iClique, cliqueCount[iClique]);
                      cliqueAdd[numberCliqueAdded] = j | 0x80000000;
                      cliqueAdd2[numberCliqueAdded] = iClique;
                      numberCliqueAdded++;
                    }
                    if (logLevel_ > 1)
                      printf("** could extend clique by adding j!\n");
                  }
                }
              }
              for (istackR = 0; istackR < nstackR; istackR++) {
                int irow = stackR[istackR];
                // switch off strengthening if not wanted
                if ((rowCuts & 2) != 0 && goingToTrueBound) {
                  bool ifCut = anyColumnCuts;
                  double gap = rowUpper[irow] - maxR[irow];
                  double sum = 0.0;
                  if (!ifCut && (gap > primalTolerance_ && gap < 1.0e8)) {
                    // see if the strengthened row is a cut
                    for (kk = rowStart[irow]; kk < rowStart[irow] + rowLength[irow];
                      kk++) {
                      sum += rowElements[kk] * colsol[column[kk]];
                    }
                    if (sum + gap * colsol[j] > rowUpper[irow] + primalTolerance_ || (info->strengthenRow && rowLower[irow] < -1.0e20)) {
                      // can be a cut
                      // add gap to integer coefficient
                      // saveU and saveL spare
                      int *index = reinterpret_cast< int * >(saveL);
                      double *element = saveU;
                      int n = 0;
                      bool coefficientExists = false;
                      for (kk = rowStart[irow]; kk < rowStart[irow] + rowLength[irow];
                        kk++) {
                        if (column[kk] != j) {
                          index[n] = column[kk];
                          element[n++] = rowElements[kk];
                        } else {
                          double value = rowElements[kk] + gap;
                          if (fabs(value) > 1.0e-12) {
                            index[n] = column[kk];
                            element[n++] = value;
                          }
                          coefficientExists = true;
                        }
                      }
                      if (!coefficientExists) {
                        index[n] = j;
                        element[n++] = gap;
                      }
                      OsiRowCut rc;
                      rc.setLb(-COIN_DBL_MAX);
                      rc.setUb(rowUpper[irow] + gap * (colUpper[j] - 1.0));
                      // effectiveness is how far j moves
                      rc.setEffectiveness((sum + gap * colsol[j] - rowUpper[irow]) / gap);
                      if (rc.effectiveness() > needEffectiveness) {
                        rc.setRow(n, index, element, false);
#ifdef CGL_DEBUG
                        if (debugger)
                          assert(!debugger->invalidCut(rc));
#endif
                        // if(info->strengthenRow)
                        // printf("c point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
                        if (rowLower[irow] < -1.0e20) {
                          printf("7Cut %g <= ", rc.lb());
                          int k;
                          for (k = 0; k < n; k++) {
                            int iColumn = index[k];
                            printf("%g*", element[k]);
                            if (intVar[iColumn])
                              printf("i%d ", iColumn);
                            else
                              printf("x%d ", iColumn);
                          }
                          printf("<= %g\n", rc.ub());
                          printf("Row %g <= ", rowLower[irow]);
                          for (k = rowStart[irow]; k < rowStart[irow] + rowLength[irow]; k++) {
                            int iColumn = column[k];
                            printf("%g*", rowElements[k]);
                            if (intVar[iColumn])
                              printf("i%d ", iColumn);
                            else
                              printf("x%d ", iColumn);
                          }
                          printf("<= %g\n", rowUpper[irow]);
                        }
#endif
                        int realRow = (rowLower[irow] < -1.0e20) ? irow : -1;
                        if (realRows && realRow > 0)
                          realRow = realRows[realRow];
                        rowCut.addCutIfNotDuplicate(rc, realRow);
                      }
                    }
                  }
                  gap = minR[irow] - rowLower[irow];
                  if (!ifCut && (gap > primalTolerance_ && gap < 1.0e8)) {
                    // see if the strengthened row is a cut
                    if (!sum) {
                      for (kk = rowStart[irow]; kk < rowStart[irow] + rowLength[irow];
                        kk++) {
                        sum += rowElements[kk] * colsol[column[kk]];
                      }
                    }
                    if (sum - gap * colsol[j] < rowLower[irow] + primalTolerance_ || (info->strengthenRow && rowUpper[irow] > 1.0e20)) {
                      // can be a cut
                      // subtract gap from integer coefficient
                      // saveU and saveL spare
                      int *index = reinterpret_cast< int * >(saveL);
                      double *element = saveU;
                      int n = 0;
                      bool coefficientExists = false;
                      for (kk = rowStart[irow]; kk < rowStart[irow] + rowLength[irow];
                        kk++) {
                        if (column[kk] != j) {
                          index[n] = column[kk];
                          element[n++] = rowElements[kk];
                        } else {
                          double value = rowElements[kk] - gap;
                          if (fabs(value) > 1.0e-12) {
                            index[n] = column[kk];
                            element[n++] = value;
                          }
                          coefficientExists = true;
                        }
                      }
                      if (!coefficientExists) {
                        index[n] = j;
                        element[n++] = -gap;
                      }
                      OsiRowCut rc;
                      rc.setLb(rowLower[irow] - gap * (colUpper[j] - 1));
                      rc.setUb(COIN_DBL_MAX);
                      // effectiveness is how far j moves
                      rc.setEffectiveness((rowLower[irow] - sum + gap * colsol[j]) / gap);
                      if (rc.effectiveness() > needEffectiveness) {
                        rc.setRow(n, index, element, false);
#ifdef CGL_DEBUG
                        if (debugger)
                          assert(!debugger->invalidCut(rc));
#endif
                        // if(info->strengthenRow)
                        // printf("d point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
                        if (rowUpper[irow] > 1.0e20) {
                          printf("8Cut %g <= ", rc.lb());
                          int k;
                          for (k = 0; k < n; k++) {
                            int iColumn = index[k];
                            printf("%g*", element[k]);
                            if (intVar[iColumn])
                              printf("i%d ", iColumn);
                            else
                              printf("x%d ", iColumn);
                          }
                          printf("<= %g\n", rc.ub());
                          printf("Row %g <= ", rowLower[irow]);
                          for (k = rowStart[irow]; k < rowStart[irow] + rowLength[irow]; k++) {
                            int iColumn = column[k];
                            printf("%g*", rowElements[k]);
                            if (intVar[iColumn])
                              printf("i%d ", iColumn);
                            else
                              printf("x%d ", iColumn);
                          }
                          printf("<= %g\n", rowUpper[irow]);
                        }
#endif
                        int realRow = (rowUpper[irow] > 1.0e20) ? irow : -1;
                        if (realRows && realRow > 0)
                          realRow = realRows[realRow];
                        rowCut.addCutIfNotDuplicate(rc, realRow);
                      }
                    }
                  }
                }
                minR[irow] = saveMin[istackR];
                maxR[irow] = saveMax[istackR];
                markR[irow] = -1;
              }
            }
          }
        }
      }
    }
  }
  if (numberCliqueAdded) {
    CoinSort_2(cliqueAdd2, cliqueAdd2 + numberCliqueAdded, cliqueAdd);
    // do cliqueStart and cliqueEntry
    CoinBigIndex numberEntries = cliqueStart_[numberCliques_];
    numberEntries += numberCliqueAdded;
    CliqueEntry *entry = new CliqueEntry[numberEntries];
    // use cliqueCount
    memset(cliqueCount, 0, numberCliques_ * sizeof(int));
    for (int i = 0; i < numberCliqueAdded; i++) {
      cliqueCount[cliqueAdd2[i]]++;
    }
    int nAdd2 = numberCliqueAdded;
    CoinBigIndex last = cliqueStart_[numberCliques_];
    CoinBigIndex put = last + numberCliqueAdded;
    cliqueStart_[numberCliques_] = put;
    for (int iClique = numberCliques_ - 1; iClique >= 0; iClique--) {
      // add in new
      for (int i = 0; i < cliqueCount[iClique]; i++) {
        --nAdd2;
        int iColumn = cliqueAdd[nAdd2] & 0x7fffffff;
        bool oneFixes = ((cliqueAdd[nAdd2] & 0x80000000) != 0);
        --put;
        setSequenceInCliqueEntry(entry[put], iColumn);
        setOneFixesInCliqueEntry(entry[put], oneFixes);
        assert(iColumn >= 0 && iColumn < numberColumns_);
      }
      // move existing
      CoinBigIndex start = cliqueStart_[iClique];
      for (CoinBigIndex i = last - 1; i >= start; i--)
        entry[--put] = cliqueEntry_[i];
      last = start;
      cliqueStart_[iClique] = put;
    }
    assert(!put);
    delete[] cliqueEntry_;
    cliqueEntry_ = entry;
    for (int j = 0; j < numberEntries; j++) {
      assert(sequenceInCliqueEntry(cliqueEntry_[j]) >= 0);
      assert(sequenceInCliqueEntry(cliqueEntry_[j]) < numberColumns_);
    }
    // now create other stuff
    int *which = cliqueAdd;
    delete[] whichClique_;
    assert(nCols == numberColumns_);
    whichClique_ = new int[numberEntries];
    for (int i = 0; i < numberColumns_; i++) {
      oneFixStart_[i] = 0;
      zeroFixStart_[i] = 0;
      endFixStart_[i] = 0;
    }
    // Now do column lists
    // First do counts
    for (int iClique = 0; iClique < numberCliques_; iClique++) {
      for (CoinBigIndex j = cliqueStart_[iClique]; j < cliqueStart_[iClique + 1]; j++) {
        int iColumn = sequenceInCliqueEntry(cliqueEntry_[j]);
        if (oneFixesInCliqueEntry(cliqueEntry_[j]))
          oneFixStart_[iColumn]++;
        else
          zeroFixStart_[iColumn]++;
      }
    }
    // now get starts and use which and end as counters
    numberEntries = 0;
    for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
      if (oneFixStart_[iColumn] >= 0) {
        int n1 = oneFixStart_[iColumn];
        int n2 = zeroFixStart_[iColumn];
        oneFixStart_[iColumn] = numberEntries;
        which[iColumn] = numberEntries;
        numberEntries += n1;
        zeroFixStart_[iColumn] = numberEntries;
        endFixStart_[iColumn] = numberEntries;
        numberEntries += n2;
      }
    }
    // now put in
    for (int iClique = 0; iClique < numberCliques_; iClique++) {
      for (CoinBigIndex j = cliqueStart_[iClique]; j < cliqueStart_[iClique + 1]; j++) {
        int iColumn = sequenceInCliqueEntry(cliqueEntry_[j]);
        if (oneFixesInCliqueEntry(cliqueEntry_[j])) {
          int put = which[iColumn];
          which[iColumn]++;
          whichClique_[put] = iClique;
        } else {
          int put = endFixStart_[iColumn];
          endFixStart_[iColumn]++;
          whichClique_[put] = iClique;
        }
      }
    }
  }
  delete[] cliqueStack;
  delete[] cliqueCount;
  delete[] to_01;
  delete[] cliqueAdd;
  delete[] stackC0;
  delete[] lo0;
  delete[] up0;
  delete[] tempL;
  delete[] tempU;
  delete[] markC;
  delete[] stackC;
  delete[] stackR;
  delete[] saveL;
  delete[] saveU;
  delete[] saveMin;
  delete[] saveMax;
  delete[] index;
  delete[] element;
  delete[] djs;
  delete[] colsol;
  // Add in row cuts
  if (!ninfeas) {
    rowCut.addCuts(cs, info->strengthenRow, 0);
  }
  return (ninfeas);
}
#if 0 // bad code
// Does probing and adding cuts for clique slacks
int 
CglProbing::probeSlacks( const OsiSolverInterface & si, 
                          const OsiRowCutDebugger *
#ifdef CGL_DEBUG
			 debugger
#endif
			 ,OsiCuts & cs, 
                          double * colLower, double * colUpper, CoinPackedMatrix *rowCopy,
			 CoinPackedMatrix *columnCopy,
                          double * rowLower, double * rowUpper,
                          char * intVar, double * minR, double * maxR,int * markR,
                          CglTreeInfo * info)
{
  if (!numberCliques_)
    return 0;
  // Set up maxes
  int maxProbe = info->inTree ? maxProbe_ : maxProbeRoot_;
  int maxStack = info->inTree ? maxStack_ : maxStackRoot_;
  int nRows=rowCopy->getNumRows();
  int nCols=rowCopy->getNumCols();
  double * colsol = new double[nCols];
  CoinMemcpyN( si.getColSolution(),nCols,colsol);
  int rowCuts=rowCuts_;
  double_int_pair * array = new double_int_pair [numberCliques_];
  // look at <= cliques
  int iClique;
  int nLook=0;
  for (iClique=0;iClique<numberCliques_;iClique++) {
    if (!cliqueType_[iClique].equality) {
      double sum=0.0;
      for (CoinBigIndex j=cliqueStart_[iClique];j<cliqueStart_[iClique+1];j++) {
        int iColumn = sequenceInCliqueEntry(cliqueEntry_[j]);
        double value = colsol[iColumn];
        if (oneFixesInCliqueEntry(cliqueEntry_[j]))
          sum += value;
        else
          sum -= value;
      }
      double away = fabs(0.5-(sum-floor(sum)));
      if (away<0.49999) {
        array[nLook].infeasibility=away;
        array[nLook++].sequence=iClique;
      }
    }
  }
  std::sort(array,array+nLook,double_int_pair_compare());
  nLook=std::min(nLook,maxProbe);
  const double * currentColLower = si.getColLower();
  const double * currentColUpper = si.getColUpper();
  double * tempL = new double [nCols];
  double * tempU = new double [nCols];
  int * markC = new int [nCols];
  int * stackC = new int [2*nCols];
  int * stackR = new int [nRows];
  double * saveL = new double [2*nCols];
  double * saveU = new double [2*nCols];
  double * saveMin = new double [nRows];
  double * saveMax = new double [nRows];
  double * element = new double[nCols];
  int * index = new int[nCols];
  // Let us never add more than twice the number of rows worth of row cuts
  // Keep cuts out of cs until end so we can find duplicates quickly
  int nRowsFake = info->inTree ? nRows/3 : nRows;
  row_cut rowCut(nRowsFake, !info->inTree);
  const int * column = rowCopy->getIndices();
  const CoinBigIndex * rowStart = rowCopy->getVectorStarts();
  const int * rowLength = rowCopy->getVectorLengths(); 
  const double * rowElements = rowCopy->getElements();
  const int * row = columnCopy->getIndices();
  const CoinBigIndex * columnStart = columnCopy->getVectorStarts();
  const int * columnLength = columnCopy->getVectorLengths(); 
  const double * columnElements = columnCopy->getElements();
  double movement;
  int i, j;
  CoinBigIndex jj, k, kk;
  int kcol,irow,krow;
  bool anyColumnCuts=false;
  double dbound, value, value2;
  int ninfeas=0;
  for (i = 0; i < nCols; ++i) {
    if (colUpper[i]-colLower[i]>1.0e-8) {
      if (colsol[i]<colLower[i]+primalTolerance_) {
        colsol[i]=colLower[i];
      } else if (colsol[i]>colUpper[i]-primalTolerance_) {
        colsol[i]=colUpper[i];
      }
    }
  }

  int ipass=0,nfixed=-1;

  /* for both way coding */
  int nstackC0=-1;
  int * stackC0 = new int[maxStack];
  double * lo0 = new double[maxStack];
  double * up0 = new double[maxStack];
  int nstackR,nstackC;
  for (i=0;i<nCols;i++) {
    if (colUpper[i]-colLower[i]<1.0e-8) {
      markC[i]=3;
    } else {
      markC[i]=0;
    }
  }
  double tolerance = 1.0e1*primalTolerance_;
  // If we are going to replace coefficient then we don't need to be effective
  int maxPass = info->inTree ? maxPass_ : maxPassRoot_;
  double needEffectiveness = info->strengthenRow ? -1.0e10 : 1.0e-3;
  while (ipass<maxPass&&nfixed) {
    int iLook;
    ipass++;
    nfixed=0;
    for (iLook=0;iLook<nLook;iLook++) {
      double solval;
      double down;
      double up;
      int iClique=array[iLook].sequence;
      solval=0.0;
      j=0; 
      for (jj=cliqueStart_[iClique];jj<cliqueStart_[iClique+1];jj++) {
        int iColumn = sequenceInCliqueEntry(cliqueEntry_[jj]);
        double value = colsol[iColumn];
        if (oneFixesInCliqueEntry(cliqueEntry_[jj]))
          solval += value;
        else
          solval -= value;
      }
      down = floor(solval+tolerance);
      up = ceil(solval-tolerance);
      int istackC,iway, istackR;
      int way[]={1,2,1};
      int feas[]={1,2,4};
      int feasible=0;
      int notFeasible;
      for (iway=0;iway<3;iway ++) {
        int fixThis=0;
        stackC[0]=j;
        markC[j]=way[iway];
        if (way[iway]==1) {
          movement=down-colUpper[j];
          assert(movement<-0.99999);
          down=colLower[j];
        } else {
          movement=up-colLower[j];
          assert(movement>0.99999);
          up=colUpper[j];
        }
        nstackC=1;
        nstackR=0;
        saveL[0]=colLower[j];
        saveU[0]=colUpper[j];
        assert (saveU[0]>saveL[0]);
        notFeasible=0;
        if (movement<0.0) {
          colUpper[j] += movement;
          colUpper[j] = floor(colUpper[j]+0.5);
#ifdef PRINT_DEBUG
          printf("** Trying %d down to 0\n",j);
#endif
        } else {
          colLower[j] += movement;
          colLower[j] = floor(colLower[j]+0.5);
#ifdef PRINT_DEBUG
          printf("** Trying %d up to 1\n",j);
#endif
        }
        if (fabs(colUpper[j]-colLower[j])<1.0e-6) 
          markC[j]=3; // say fixed
        istackC=0;
        /* update immediately */
        for (k=columnStart[j];k<columnStart[j]+columnLength[j];k++) {
          int irow = row[k];
          value = columnElements[k];
          if (markR[irow]==-1) {
            stackR[nstackR]=irow;
            markR[irow]=nstackR;
            saveMin[nstackR]=minR[irow];
            saveMax[nstackR]=maxR[irow];
            nstackR++;
          } else if (markR[irow]==-2) {
            continue;
          }
          /* could check immediately if violation */
          if (movement>0.0) {
            /* up */
            if (value>0.0) {
              /* up does not change - down does */
              if (minR[irow]>-1.0e10)
                minR[irow] += value;
              if (minR[irow]>rowUpper[irow]+1.0e-5) {
                notFeasible=1;
                istackC=1;
                break;
              }
            } else {
              /* down does not change - up does */
              if (maxR[irow]<1.0e10)
                maxR[irow] += value;
              if (maxR[irow]<rowLower[irow]-1.0e-5) {
                notFeasible=1;
                istackC=1;
                break;
              }
            }
          } else {
            /* down */
            if (value<0.0) {
              /* up does not change - down does */
              if (minR[irow]>-1.0e10)
                minR[irow] -= value;
              if (minR[irow]>rowUpper[irow]+1.0e-5) {
                notFeasible=1;
                istackC=1;
                break;
              }
            } else {
              /* down does not change - up does */
              if (maxR[irow]<1.0e10)
                maxR[irow] -= value;
              if (maxR[irow]<rowLower[irow]-1.0e-5) {
                notFeasible=1;
                istackC=1;
                break;
              }
            }
          }
        }
        while (istackC<nstackC&&nstackC<maxStack) {
          int jway;
          int jcol =stackC[istackC];
          jway=markC[jcol];
          // If not first and fixed then skip
          if (jway==3&&istackC) {
            //istackC++;
            //continue;
            //printf("fixed %d on stack\n",jcol);
          }
          // Do cliques
          if (oneFixStart_&&oneFixStart_[jcol]>=0) {
            int start;
            int end;
            if (colLower[jcol]>saveL[istackC]) {
              // going up
              start = oneFixStart_[jcol];
              end = zeroFixStart_[jcol];
            } else {
              assert (colUpper[jcol]<saveU[istackC]);
              // going down
              start = zeroFixStart_[jcol];
              end = endFixStart_[jcol];
            }
            for (int i=start;i<end;i++) {
              int iClique = whichClique_[i];
              for (CoinBigIndex k=cliqueStart_[iClique];k<cliqueStart_[iClique+1];k++) {
                int kcol = sequenceInCliqueEntry(cliqueEntry_[k]);
                if (jcol==kcol)
                  continue;
                int kway = oneFixesInCliqueEntry(cliqueEntry_[k]);
                if (kcol!=jcol) {
                  if (!markC[kcol]) {
                    // not on list yet
                    if (nstackC<2*maxStack) {
                      markC[kcol] = 3; // say fixed
                      fixThis++;
                      stackC[nstackC]=kcol;
                      saveL[nstackC]=colLower[kcol];
                      saveU[nstackC]=colUpper[kcol];
                      assert (saveU[nstackC]>saveL[nstackC]);
                      nstackC++;
                      if (!kway) {
                        // going up
                        colLower[kcol]=1.0;
                        /* update immediately */
                        for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                          krow = row[jj];
                          value = columnElements[jj];
                          if (markR[krow]==-1) {
                            stackR[nstackR]=krow;
                            markR[krow]=nstackR;
                            saveMin[nstackR]=minR[krow];
                            saveMax[nstackR]=maxR[krow];
                            nstackR++;
                          } else if (markR[krow]==-2) {
                            continue;
                          }
                          /* could check immediately if violation */
                          /* up */
                          if (value>0.0) {
                            /* up does not change - down does */
                            if (minR[krow]>-1.0e10)
                              minR[krow] += value;
                            if (minR[krow]>rowUpper[krow]+1.0e-5) {
                              colUpper[kcol]=-1.0e50; /* force infeasible */
                              break;
                            }
                          } else {
                            /* down does not change - up does */
                            if (maxR[krow]<1.0e10)
                              maxR[krow] += value;
                            if (maxR[krow]<rowLower[krow]-1.0e-5) {
                              notFeasible=1;
                              break;
                            }
                          }
                        }
                      } else {
                        // going down
                        colUpper[kcol]=0.0;
                        /* update immediately */
                        for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                          krow = row[jj];
                          value = columnElements[jj];
                          if (markR[krow]==-1) {
                            stackR[nstackR]=krow;
                            markR[krow]=nstackR;
                            saveMin[nstackR]=minR[krow];
                            saveMax[nstackR]=maxR[krow];
                            nstackR++;
                          } else if (markR[krow]==-2) {
                            continue;
                          }
                          /* could check immediately if violation */
                          /* down */
                          if (value<0.0) {
                            /* up does not change - down does */
                            if (minR[krow]>-1.0e10)
                              minR[krow] -= value;
                            if (minR[krow]>rowUpper[krow]+1.0e-5) {
                              notFeasible=1;
                              break;
                            }
                          } else {
                            /* down does not change - up does */
                            if (maxR[krow]<1.0e10)
                              maxR[krow] -= value;
                            if (maxR[krow]<rowLower[krow]-1.0e-5) {
                              notFeasible=1;
                              break;
                            }
                          }
                        }
                      }
                    }
                  } else if (markC[kcol]==1) {
                    // marked as going to 0
                    assert (!colUpper[kcol]);
                    if (!kway) {
                      // contradiction
                      notFeasible=1;
                      break;
                    }
                  } else if (markC[kcol]==2) {
                    // marked as going to 1
                    assert (colLower[kcol]);
                    if (kway) {
                      // contradiction
                      notFeasible=1;
                      break;
                    }
                  } else {
                    // marked as fixed
                    assert (markC[kcol]==3);
                    int jkway;
                    if (colLower[kcol])
                      jkway=1;
                    else
                      jkway=0;
                    if (kway==jkway) {
                      // contradiction
                      notFeasible=1;
                      break;
                    }
                  }
                }
              }
              if (notFeasible)
                break;
            }
            if (notFeasible)
              istackC=nstackC+1;
          }
          for (k=columnStart[jcol];k<columnStart[jcol]+columnLength[jcol];k++) {
            // break if found not feasible
            if (notFeasible)
              break;
            irow = row[k];
            /*value = columnElements[k];*/
            if (markR[irow]!=-2) {
              /* see if anything forced */
              for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];kk++) {
                double moveUp=0.0;
                double moveDown=0.0;
                double newUpper=-1.0,newLower=1.0;
                kcol=column[kk];
                bool onList = (markC[kcol]!=0);
                if (markC[kcol]!=3) {
                  value2=rowElements[kk];
                  int markIt=markC[kcol];
                  if (value2 < 0.0) {
                    if (colUpper[kcol] < 1e10 && (markIt&2)==0 &&
                        rowUpper[irow]<1.0e10) {
                      dbound = colUpper[kcol]+
                        (rowUpper[irow]-minR[irow])/value2;
                      if (dbound > colLower[kcol] + primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 2;
                          newLower = ceil(dbound-primalTolerance_);
                        } else {
                          newLower=dbound;
                          if (newLower+primalTolerance_>colUpper[kcol]&&
                              newLower-primalTolerance_<=colUpper[kcol]) {
                            newLower=colUpper[kcol];
                            markIt |= 2;
                            markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt=3;
                          }
                        }
                        moveUp = newLower-colLower[kcol];
                      }
                    }
                    if (colLower[kcol] > -1e10 && (markIt&1)==0 &&
                        rowLower[irow]>-1.0e10) {
                      dbound = colLower[kcol] + 
                        (rowLower[irow]-maxR[irow])/value2;
                      if (dbound < colUpper[kcol] - primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 1;
                          newUpper = floor(dbound+primalTolerance_);
                        } else {
                          newUpper=dbound;
                          if (newUpper-primalTolerance_<colLower[kcol]&&
                              newUpper+primalTolerance_>=colLower[kcol]) {
                            newUpper=colLower[kcol];
                            markIt |= 1;
                            markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt=3;
                          }
                        }
                        moveDown = newUpper-colUpper[kcol];
                      }
                    }
                  } else {
                    /* positive element */
                    if (colUpper[kcol] < 1e10 && (markIt&2)==0 &&
                        rowLower[irow]>-1.0e10) {
                      dbound = colUpper[kcol] + 
                        (rowLower[irow]-maxR[irow])/value2;
                      if (dbound > colLower[kcol] + primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 2;
                          newLower = ceil(dbound-primalTolerance_);
                        } else {
                          newLower=dbound;
                          if (newLower+primalTolerance_>colUpper[kcol]&&
                              newLower-primalTolerance_<=colUpper[kcol]) {
                            newLower=colUpper[kcol];
                            markIt |= 2;
                            markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt=3;
                          }
                        }
                        moveUp = newLower-colLower[kcol];
                      }
                    }
                    if (colLower[kcol] > -1e10 && (markIt&1)==0 &&
                        rowUpper[irow]<1.0e10) {
                      dbound = colLower[kcol] + 
                        (rowUpper[irow]-minR[irow])/value2;
                      if (dbound < colUpper[kcol] - primalTolerance_) {
                        if (intVar[kcol]) {
                          markIt |= 1;
                          newUpper = floor(dbound+primalTolerance_);
                        } else {
                          newUpper=dbound;
                          if (newUpper-primalTolerance_<colLower[kcol]&&
                              newUpper+primalTolerance_>=colLower[kcol]) {
                            newUpper=colLower[kcol];
                            markIt |= 1;
                            markIt=3;
                          } else {
                            // avoid problems - fix later ?
                            markIt=3;
                          }
                        }
                        moveDown = newUpper-colUpper[kcol];
                      }
                    }
                  }
                  if (nstackC<2*maxStack) {
                    markC[kcol] = markIt;
		  }
                  if (moveUp&&nstackC<2*maxStack) {
                    fixThis++;
#ifdef PRINT_DEBUG
                    printf("lower bound on %d increased from %g to %g by row %d %g %g\n",kcol,colLower[kcol],newLower,irow,rowLower[irow],rowUpper[irow]);
                    value=0.0;
                    for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
                      int ii=column[jj];
                      if (colUpper[ii]-colLower[ii]>primalTolerance_) {
                        printf("(%d, %g) ",ii,rowElements[jj]);
                      } else {
                        value += rowElements[jj]*colLower[ii];
                      }
                    }
                    printf(" - fixed %g\n",value);
                    for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
                      int ii=column[jj];
                      if (colUpper[ii]-colLower[ii]<primalTolerance_) {
                        printf("(%d, %g, %g) ",ii,rowElements[jj],colLower[ii]);
                      }
                    }
                    printf("\n");
#endif
                    if (!onList) {
                      stackC[nstackC]=kcol;
                      saveL[nstackC]=colLower[kcol];
                      saveU[nstackC]=colUpper[kcol];
                      assert (saveU[nstackC]>saveL[nstackC]);
                      nstackC++;
                      onList=true;
                    }
                    if (intVar[kcol])
                      newLower = std::max(colLower[kcol],ceil(newLower-1.0e-4));
                    colLower[kcol]=newLower;
                    if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
                      markC[kcol]=3; // say fixed
		    }
                    /* update immediately */
                    for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                      krow = row[jj];
                      value = columnElements[jj];
                      if (markR[krow]==-1) {
                        stackR[nstackR]=krow;
                        markR[krow]=nstackR;
                        saveMin[nstackR]=minR[krow];
                        saveMax[nstackR]=maxR[krow];
                        nstackR++;
                      } else if (markR[krow]==-2) {
                        continue;
                      }
                      /* could check immediately if violation */
                      /* up */
                      if (value>0.0) {
                        /* up does not change - down does */
                        if (minR[krow]>-1.0e10)
                          minR[krow] += value*moveUp;
                        if (minR[krow]>rowUpper[krow]+1.0e-5) {
                          colUpper[kcol]=-1.0e50; /* force infeasible */
                          break;
                        }
                      } else {
                        /* down does not change - up does */
                        if (maxR[krow]<1.0e10)
                          maxR[krow] += value*moveUp;
                        if (maxR[krow]<rowLower[krow]-1.0e-5) {
                          colUpper[kcol]=-1.0e50; /* force infeasible */
                          break;
                        }
                      }
                    }
                  }
                  if (moveDown&&nstackC<2*maxStack) {
                    fixThis++;
#ifdef PRINT_DEBUG
                    printf("upper bound on %d decreased from %g to %g by row %d %g %g\n",kcol,colUpper[kcol],newUpper,irow,rowLower[irow],rowUpper[irow]);
                    value=0.0;
                    for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
                      int ii=column[jj];
                      if (colUpper[ii]-colLower[ii]>primalTolerance_) {
                        printf("(%d, %g) ",ii,rowElements[jj]);
                      } else {
                        value += rowElements[jj]*colLower[ii];
                      }
                    }
                    printf(" - fixed %g\n",value);
                    for (jj=rowStart[irow];jj<rowStart[irow]+rowLength[irow];jj++) {
                      int ii=column[jj];
                      if (colUpper[ii]-colLower[ii]<primalTolerance_) {
                        printf("(%d, %g, %g) ",ii,rowElements[jj],colLower[ii]);
                      }
                    }
                    printf("\n");
#endif
                    if (!onList) {
                      stackC[nstackC]=kcol;
                      saveL[nstackC]=colLower[kcol];
                      saveU[nstackC]=colUpper[kcol];
                      assert (saveU[nstackC]>saveL[nstackC]);
                      nstackC++;
                      onList=true;
                    }
                    if (intVar[kcol])
                      newUpper = std::min(colUpper[kcol],floor(newUpper+1.0e-4));
                    colUpper[kcol]=newUpper;
                    if (fabs(colUpper[kcol]-colLower[kcol])<1.0e-6) {
                      markC[kcol]=3; // say fixed
		    }
                    /* update immediately */
                    for (jj=columnStart[kcol];jj<columnStart[kcol]+columnLength[kcol];jj++) {
                      krow = row[jj];
                      value = columnElements[jj];
                      if (markR[krow]==-1) {
                        stackR[nstackR]=krow;
                        markR[krow]=nstackR;
                        saveMin[nstackR]=minR[krow];
                        saveMax[nstackR]=maxR[krow];
                        nstackR++;
                      } else if (markR[krow]==-2) {
                        continue;
                      }
                      /* could check immediately if violation */
                      /* down */
                      if (value<0.0) {
                        /* up does not change - down does */
                        if (minR[krow]>-1.0e10)
                          minR[krow] += value*moveDown;
                        if (minR[krow]>rowUpper[krow]+1.0e-5) {
                          colUpper[kcol]=-1.0e50; /* force infeasible */
                          break;
                        }
                      } else {
                        /* down does not change - up does */
                        if (maxR[krow]<1.0e10)
                          maxR[krow] += value*moveDown;
                        if (maxR[krow]<rowLower[krow]-1.0e-5) {
                          colUpper[kcol]=-1.0e50; /* force infeasible */
                          break;
                        }
                      }
                    }
                  }
                  if (colLower[kcol]>colUpper[kcol]+primalTolerance_) {
                    notFeasible=1;;
                    k=columnStart[jcol]+columnLength[jcol];
                    istackC=nstackC+1;
#ifdef PRINT_DEBUG
                    printf("** not feasible this way\n");
#endif
                    break;
                  }
                }
              }
            }
          }
          istackC++;
        }
        if (!notFeasible) {
          feasible |= feas[iway];
        } else if (iway==1&&feasible==0) {
          /* not feasible at all */
          ninfeas=1;
          j=nCols-1;
          iLook=nLook;
          ipass=maxPass;
          break;
        }
        if (iway==2||(iway==1&&feasible==2)) {
          /* keep */
          iway=3;
          nfixed++;
          if (mode_) {
            OsiColCut cc;
            int nTot=0,nFix=0,nInt=0;
            bool ifCut=false;
            for (istackC=0;istackC<nstackC;istackC++) {
              int icol=stackC[istackC];
              if (intVar[icol]) {
                if (colUpper[icol]<currentColUpper[icol]-1.0e-4) {
                  element[nFix]=colUpper[icol];
                  index[nFix++]=icol;
                  nInt++;
                  if (colsol[icol]>colUpper[icol]+primalTolerance_) {
                    ifCut=true;
                    anyColumnCuts=true;
                  }
                }
              }
            }
            if (nFix) {
              nTot=nFix;
	      cleanBounds(1,nFix,index,element);
              cc.setUbs(nFix,index,element);
              nFix=0;
            }
            for (istackC=0;istackC<nstackC;istackC++) {
              int icol=stackC[istackC];
              if (intVar[icol]) {
                if (colLower[icol]>currentColLower[icol]+1.0e-4) {
                  element[nFix]=colLower[icol];
                  index[nFix++]=icol;
                  nInt++;
                  if (colsol[icol]<colLower[icol]-primalTolerance_) {
                    ifCut=true;
                    anyColumnCuts=true;
                  }
                }
              }
            }
            if (nFix) {
              nTot+=nFix;
	      cleanBounds(-1,nFix,index,element);
              cc.setLbs(nFix,index,element);
            }
            // could tighten continuous as well
            if (nInt) {
              if (ifCut) {
                cc.setEffectiveness(100.0);
              } else {
                cc.setEffectiveness(1.0e-5);
              }
#ifdef CGL_DEBUG
              checkBounds(debugger,cc);
#endif
              cs.insert(cc);
            }
          }
          for (istackC=0;istackC<nstackC;istackC++) {
            int icol=stackC[istackC];
            if (colUpper[icol]-colLower[icol]>primalTolerance_) {
              markC[icol]=0;
            } else {
              markC[icol]=3;
            }
          }
          for (istackR=0;istackR<nstackR;istackR++) {
            int irow=stackR[istackR];
            markR[irow]=-1;
          }
        } else {
          /* is it worth seeing if can increase coefficients
             or maybe better see if it is a cut */
          if (iway==0) {
            nstackC0=std::min(nstackC,maxStack);
            if (notFeasible) {
              nstackC0=0;
            } else {
              for (istackC=0;istackC<nstackC0;istackC++) {
                int icol=stackC[istackC];
                stackC0[istackC]=icol;
                lo0[istackC]=colLower[icol];
                up0[istackC]=colUpper[icol];
              }
            }
            /* restore all */
            assert (iway==0);
            for (istackC=nstackC-1;istackC>=0;istackC--) {
              int icol=stackC[istackC];
              double oldU=saveU[istackC];
              double oldL=saveL[istackC];
              colUpper[icol]=oldU;
              colLower[icol]=oldL;
              markC[icol]=0;
            }
            for (istackR=0;istackR<nstackR;istackR++) {
              int irow=stackR[istackR];
              // switch off strengthening if not wanted
              if ((rowCuts&2)!=0) {
                bool ifCut=anyColumnCuts;
                double gap = rowUpper[irow]-maxR[irow];
                double sum=0.0;
                if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
                  // see if the strengthened row is a cut
                  for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                       kk++) {
                    sum += rowElements[kk]*colsol[column[kk]];
                  }
                  if (sum-gap*colsol[j]>maxR[irow]+primalTolerance_||info->strengthenRow) {
                    // can be a cut
                    // subtract gap from upper and integer coefficient
                    // saveU and saveL spare
                    int * index = reinterpret_cast<int *>(saveL);
                    double * element = saveU;
                    int n=0;
                    bool coefficientExists=false;
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      if (column[kk]!=j) {
                        index[n]=column[kk];
                        element[n++]=rowElements[kk];
                      } else {
                        double value=rowElements[kk]-gap;
                        if (fabs(value)>1.0e-12) {
                          index[n]=column[kk];
                          element[n++]=value;
                        }
                        coefficientExists=true;
                      }
                    }
                    if (!coefficientExists) {
                      index[n]=j;
                      element[n++]=-gap;
                    }
                    OsiRowCut rc;
                    rc.setLb(-COIN_DBL_MAX);
                    rc.setUb(rowUpper[irow]-gap*(colLower[j]+1.0));
                    // effectiveness is how far j moves
                    rc.setEffectiveness((sum-gap*colsol[j]-maxR[irow])/gap);
                    if (rc.effectiveness()>needEffectiveness) {
                      rc.setRow(n,index,element,false);
#ifdef CGL_DEBUG
                      if (debugger) assert(!debugger->invalidCut(rc));
#endif
                      // If strengthenRow point to row
                      //if(info->strengthenRow)
                      //printf("a point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
		      {
			printf("9Cut %g <= ",rc.lb());
			int k;
			for ( k=0;k<n;k++) {
			  int iColumn = index[k];
			  printf("%g*",element[k]);
			  if (intVar[iColumn])
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rc.ub());
			printf("Row %g <= ",rowLower[irow]);
			for (k=rowStart[irow];k<rowStart[irow]+rowLength[irow];k++) {
			  int iColumn = column[k];
			  printf("%g*",rowElements[k]);
			  if (intVar[iColumn])
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rowUpper[irow]);
		      }
#endif
                      rowCut.addCutIfNotDuplicate(rc,irow);
                    }
                  }
                }
                gap = minR[irow]-rowLower[irow];
                if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
                  // see if the strengthened row is a cut
                  if (!sum) {
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      sum += rowElements[kk]*colsol[column[kk]];
                    }
                  }
                  if (sum+gap*colsol[j]<minR[irow]+primalTolerance_||info->strengthenRow) {
                    // can be a cut
                    // add gap to lower and integer coefficient
                    // saveU and saveL spare
                    int * index = reinterpret_cast<int *>(saveL);
                    double * element = saveU;
                    int n=0;
                    bool coefficientExists=false;
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      if (column[kk]!=j) {
                        index[n]=column[kk];
                        element[n++]=rowElements[kk];
                      } else {
                        double value=rowElements[kk]+gap;
                        if (fabs(value)>1.0e-12) {
                          index[n]=column[kk];
                          element[n++]=value;
                        }
                        coefficientExists=true;
                      }
                    }
                    if (!coefficientExists) {
                      index[n]=j;
                      element[n++]=gap;
                    }
                    OsiRowCut rc;
                    rc.setLb(rowLower[irow]+gap*(colLower[j]+1.0));
                    rc.setUb(COIN_DBL_MAX);
                    // effectiveness is how far j moves
                    rc.setEffectiveness((minR[irow]-sum-gap*colsol[j])/gap);
                    if (rc.effectiveness()>needEffectiveness) {
                      rc.setRow(n,index,element,false);
#ifdef CGL_DEBUG
                      if (debugger) assert(!debugger->invalidCut(rc));
#endif
                      //if(info->strengthenRow)
                      //printf("b point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
		      {
			printf("10Cut %g <= ",rc.lb());
			int k;
			for ( k=0;k<n;k++) {
			  int iColumn = index[k];
			  printf("%g*",element[k]);
			  if (intVar[iColumn])
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rc.ub());
			printf("Row %g <= ",rowLower[irow]);
			for (k=rowStart[irow];k<rowStart[irow]+rowLength[irow];k++) {
			  int iColumn = column[k];
			  printf("%g*",rowElements[k]);
			  if (intVar[iColumn])
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rowUpper[irow]);
		      }
#endif
                      rowCut.addCutIfNotDuplicate(rc,irow);
                    }
                  }
                }
              }
              minR[irow]=saveMin[istackR];
              maxR[irow]=saveMax[istackR];
              markR[irow]=-1;
            }
          } else {
            if (iway==1&&feasible==3) {
              iway=3;
              /* point back to stack */
              for (istackC=nstackC-1;istackC>=0;istackC--) {
                int icol=stackC[istackC];
                markC[icol]=istackC+100000;
              }
              if (mode_) {
                OsiColCut cc;
                int nTot=0,nFix=0,nInt=0;
                bool ifCut=false;
                for (istackC=0;istackC<nstackC0;istackC++) {
                  int icol=stackC0[istackC];
                  int istackC1=markC[icol]-100000;
                  if (istackC1>=0) {
		    markC[icol] -= 100000;
                    if (std::min(lo0[istackC],colLower[icol])>saveL[istackC1]+1.0e-4) {
                      saveL[istackC1]=std::min(lo0[istackC],colLower[icol]);
                      if (intVar[icol]) {
                        element[nFix]=saveL[istackC1];
                        index[nFix++]=icol;
                        nInt++;
                        if (colsol[icol]<saveL[istackC1]-primalTolerance_)
                          ifCut=true;
                      }
                      nfixed++;
                    }
                  }
                }
                if (nFix) {
                  nTot=nFix;
		  cleanBounds(-1,nFix,index,element);
                  cc.setLbs(nFix,index,element);
                  nFix=0;
                }
                for (istackC=0;istackC<nstackC0;istackC++) {
                  int icol=stackC0[istackC];
                  int istackC1=markC[icol]-100000;
                  if (istackC1>=0) {
		  markC[icol] -= 100000;
		  if (std::max(up0[istackC],colUpper[icol])<saveU[istackC1]-1.0e-4) {
                      saveU[istackC1]=std::max(up0[istackC],colUpper[icol]);
                      if (intVar[icol]) {
                        element[nFix]=saveU[istackC1];
                        index[nFix++]=icol;
                        nInt++;
                        if (colsol[icol]>saveU[istackC1]+primalTolerance_)
                          ifCut=true;
                      }
                      nfixed++;
                    }
                  }
                }
                if (nFix) {
                  nTot+=nFix;
		  cleanBounds(1,nFix,index,element);
                  cc.setUbs(nFix,index,element);
                }
                // could tighten continuous as well
                if (nInt) {
                  if (ifCut) {
                    cc.setEffectiveness(100.0);
                  } else {
                    cc.setEffectiveness(1.0e-5);
                  }
#ifdef CGL_DEBUG
                  checkBounds(debugger,cc);
#endif
                  cs.insert(cc);
                }
              }
            }
            /* restore all */
            for (istackC=nstackC-1;istackC>=0;istackC--) {
              int icol=stackC[istackC];
              double oldU=saveU[istackC];
              double oldL=saveL[istackC];
              colUpper[icol]=oldU;
              colLower[icol]=oldL;
              if (oldU>oldL+1.0e-4)
                markC[icol]=0;
              else
                markC[icol]=3;
            }
            for (istackR=0;istackR<nstackR;istackR++) {
              int irow=stackR[istackR];
              // switch off strengthening if not wanted
              if ((rowCuts&2)!=0) {
                bool ifCut=anyColumnCuts;
                double gap = rowUpper[irow]-maxR[irow];
                double sum=0.0;
                if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
                  // see if the strengthened row is a cut
                  for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                       kk++) {
                    sum += rowElements[kk]*colsol[column[kk]];
                  }
                  if (sum+gap*colsol[j]>rowUpper[irow]+primalTolerance_||info->strengthenRow) {
                    // can be a cut
                    // add gap to integer coefficient
                    // saveU and saveL spare
                    int * index = reinterpret_cast<int *>(saveL);
                    double * element = saveU;
                    int n=0;
                    bool coefficientExists=false;
                    for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
                         kk++) {
                      if (column[kk]!=j) {
                        index[n]=column[kk];
                        element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]+gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(-COIN_DBL_MAX);
		      rc.setUb(rowUpper[irow]+gap*(colUpper[j]-1.0));
		      // effectiveness is how far j moves
		      rc.setEffectiveness((sum+gap*colsol[j]-rowUpper[irow])/gap);
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element,false);
#ifdef CGL_DEBUG
			if (debugger) assert(!debugger->invalidCut(rc));
#endif
                        //if(info->strengthenRow)
                        //printf("c point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
		      {
			printf("11Cut %g <= ",rc.lb());
			int k;
			for ( k=0;k<n;k++) {
			  int iColumn = index[k];
			  printf("%g*",element[k]);
			  if (intVar[iColumn])
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rc.ub());
			printf("Row %g <= ",rowLower[irow]);
			for (k=rowStart[irow];k<rowStart[irow]+rowLength[irow];k++) {
			  int iColumn = column[k];
			  printf("%g*",rowElements[k]);
			  if (intVar[iColumn])
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rowUpper[irow]);
		      }
#endif
			rowCut.addCutIfNotDuplicate(rc,irow);
		      }
		    }
		  }
		  gap = minR[irow]-rowLower[irow];
		  if (!ifCut&&(gap>primalTolerance_&&gap<1.0e8)) {
		    // see if the strengthened row is a cut
		    if (!sum) {
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			sum += rowElements[kk]*colsol[column[kk]];
		      }
		    }
		    if (sum-gap*colsol[j]<rowLower[irow]+primalTolerance_||info->strengthenRow) {
		      // can be a cut
		      // subtract gap from integer coefficient
		      // saveU and saveL spare
		      int * index = reinterpret_cast<int *>(saveL);
		      double * element = saveU;
		      int n=0;
                      bool coefficientExists=false;
		      for (kk=rowStart[irow];kk<rowStart[irow]+rowLength[irow];
			   kk++) {
			if (column[kk]!=j) {
			  index[n]=column[kk];
			  element[n++]=rowElements[kk];
			} else {
			  double value=rowElements[kk]-gap;
			  if (fabs(value)>1.0e-12) {
			    index[n]=column[kk];
			    element[n++]=value;
			  }
			  coefficientExists=true;
			}
		      }
		      if (!coefficientExists) {
			index[n]=j;
			element[n++]=-gap;
		      }
		      OsiRowCut rc;
		      rc.setLb(rowLower[irow]-gap*(colUpper[j]-1));
		      rc.setUb(COIN_DBL_MAX);
		      // effectiveness is how far j moves
		      rc.setEffectiveness((rowLower[irow]-sum+gap*colsol[j])/gap);
		      if (rc.effectiveness()>needEffectiveness) {
			rc.setRow(n,index,element,false);
#ifdef CGL_DEBUG
			if (debugger) assert(!debugger->invalidCut(rc));
#endif
                        //if(info->strengthenRow)
                        //printf("d point to row %d\n",irow);
#ifdef STRENGTHEN_PRINT
		      {
			printf("12Cut %g <= ",rc.lb());
			int k;
			for ( k=0;k<n;k++) {
			  int iColumn = index[k];
			  printf("%g*",element[k]);
			  if (intVar[iColumn])
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rc.ub());
			printf("Row %g <= ",rowLower[irow]);
			for (k=rowStart[irow];k<rowStart[irow]+rowLength[irow];k++) {
			  int iColumn = column[k];
			  printf("%g*",rowElements[k]);
			  if (intVar[iColumn])
			    printf("i%d ",iColumn);
			  else
			    printf("x%d ",iColumn);
			}
			printf("<= %g\n",rowUpper[irow]);
		      }
#endif
			rowCut.addCutIfNotDuplicate(rc,irow);
		      }
		    }
		  }
		}
		minR[irow]=saveMin[istackR];
		maxR[irow]=saveMax[istackR];
		markR[irow]=-1;
	      }
	    }
	  }
        }
    }
  }
  delete [] stackC0;
  delete [] lo0;
  delete [] up0;
  delete [] tempL;
  delete [] tempU;
  delete [] markC;
  delete [] stackC;
  delete [] stackR;
  delete [] saveL;
  delete [] saveU;
  delete [] saveMin;
  delete [] saveMax;
  delete [] index;
  delete [] element;
  delete [] colsol;
  // Add in row cuts
  if (!ninfeas) {
    rowCut.addCuts(cs,info->strengthenRow,0);
  }
  delete [] array;
  abort();
  return (ninfeas);
}
#endif
// Create a copy of matrix which is to be used
// this is to speed up process and to give global cuts
// Can give an array with 1 set to select, 0 to ignore
// column bounds are tightened
// If array given then values of 1 will be set to 0 if redundant
int CglProbing::snapshot(const OsiSolverInterface &si,
  char *possible, bool withObjective)
{
  deleteSnapshot();
  // Get basic problem information

  numberColumns_ = si.getNumCols();
  numberRows_ = si.getNumRows();
  colLower_ = new double[numberColumns_];
  colUpper_ = new double[numberColumns_];
  CoinMemcpyN(si.getColLower(), numberColumns_, colLower_);
  CoinMemcpyN(si.getColUpper(), numberColumns_, colUpper_);
  rowLower_ = new double[numberRows_ + 1];
  rowUpper_ = new double[numberRows_ + 1];
  CoinMemcpyN(si.getRowLower(), numberRows_, rowLower_);
  CoinMemcpyN(si.getRowUpper(), numberRows_, rowUpper_);

  int i;
  if (possible) {
    for (i = 0; i < numberRows_; i++) {
      if (!possible[i]) {
        rowLower_[i] = -COIN_DBL_MAX;
        rowUpper_[i] = COIN_DBL_MAX;
      }
    }
  }

  // get integer variables
  const char *intVarOriginal = si.getColType(true);
  char *intVar = CoinCopyOfArray(intVarOriginal, numberColumns_);
  numberIntegers_ = 0;
  number01Integers_ = 0;
  for (i = 0; i < numberColumns_; i++) {
    if (intVar[i]) {
      numberIntegers_++;
      if (intVar[i] == 1) {
        number01Integers_++;
      }
    }
  }

  rowCopy_ = new CoinPackedMatrix(*si.getMatrixByRow());

  int *column = rowCopy_->getMutableIndices();
  const CoinBigIndex *rowStart = rowCopy_->getVectorStarts();
  const int *rowLength = rowCopy_->getVectorLengths();
  double *rowElements = rowCopy_->getMutableElements();
  // Put negative first
  int *column2 = new int[numberColumns_];
  double *elements2 = new double[numberColumns_];
  CoinBigIndex *rowStartPos = new CoinBigIndex[numberRows_];
  for (int i = 0; i < numberRows_; i++) {
    CoinBigIndex start = rowStart[i];
    CoinBigIndex end = start + rowLength[i];
    int nOther = 0;
    for (CoinBigIndex j = start; j < end; j++) {
      int iColumn = column[j];
      double value = rowElements[j];
      if (value < 0.0) {
        rowElements[start] = value;
        column[start++] = iColumn;
      } else {
        elements2[nOther] = value;
        column2[nOther++] = iColumn;
      }
    }
    rowStartPos[i] = start;
    for (int k = 0; k < nOther; k++) {
      rowElements[start] = elements2[k];
      column[start++] = column2[k];
    }
  }
  delete[] column2;
  delete[] elements2;

  int returnCode = 0;
  int ninfeas = tighten(colLower_, colUpper_, column, rowElements,
    rowStart, rowStartPos, rowLength,
    si.getMatrixByCol(),
    rowLower_, rowUpper_,
    numberRows_, numberColumns_, intVar, 5, primalTolerance_);
  delete[] rowStartPos;
  if (ninfeas) {
    // let someone else find out
    returnCode = 1;
  }
  /*
    QUESTION: If ninfeas > 1 (one or more variables infeasible), shouldn't we
              bail out here?
  */

  // do integer stuff for mode 0
  cutVector_ = new disaggregation[number01Integers_];
  memset(cutVector_, 0, number01Integers_ * sizeof(disaggregation));
  number01Integers_ = 0;
  for (i = 0; i < numberColumns_; i++) {
    if (intVar[i] == 1)
      cutVector_[number01Integers_++].sequence = i;
  }
  delete[] intVar;

  // now delete rows
  if (possible) {
    for (i = 0; i < numberRows_; i++) {
      if (rowLower_[i] < -1.0e30 && rowUpper_[i] > 1.0e30)
        possible[i] = 0;
    }
  }
  int *index = new int[numberRows_];
  int nDrop = 0, nKeep = 0;
  for (i = 0; i < numberRows_; i++) {
    if (rowLower_[i] < -1.0e30 && rowUpper_[i] > 1.0e30 && false) {
      index[nDrop++] = i;
    } else {
      rowLower_[nKeep] = rowLower_[i];
      rowUpper_[nKeep++] = rowUpper_[i];
    }
  }
  numberRows_ = nKeep;
  if (nDrop)
    rowCopy_->deleteRows(nDrop, index);
  delete[] index;
  if (withObjective) {
    // add in objective
    int *columns = new int[numberColumns_];
    double *elements = new double[numberColumns_];
    int n = 0;
    const double *objective = si.getObjCoefficients();
    bool maximize = (si.getObjSense() == -1);
    for (i = 0; i < numberColumns_; i++) {
      if (objective[i]) {
        elements[n] = (maximize) ? -objective[i] : objective[i];
        columns[n++] = i;
      }
    }
    rowCopy_->appendRow(n, columns, elements);
    delete[] columns;
    delete[] elements;
    numberRows_++;
  }
  // create column copy
  if (rowCopy_->getNumElements()) {
    columnCopy_ = new CoinPackedMatrix(*rowCopy_, 0, 0, true);
  } else {
    columnCopy_ = new CoinPackedMatrix();
  }
  // make sure big enough - in case too many rows dropped
  columnCopy_->setDimensions(numberRows_, numberColumns_);
  rowCopy_->setDimensions(numberRows_, numberColumns_);
  return returnCode;
}
// Delete snapshot
void CglProbing::deleteSnapshot()
{
  delete[] rowLower_;
  delete[] rowUpper_;
  delete[] colLower_;
  delete[] colUpper_;
  delete rowCopy_;
  delete columnCopy_;
  rowCopy_ = NULL;
  columnCopy_ = NULL;
  rowLower_ = NULL;
  rowUpper_ = NULL;
  colLower_ = NULL;
  colUpper_ = NULL;
  int i;
  for (i = 0; i < number01Integers_; i++) {
    delete[] cutVector_[i].index;
  }
  delete[] cutVector_;
  numberIntegers_ = 0;
  number01Integers_ = 0;
  cutVector_ = NULL;
}
// Mode stuff
void CglProbing::setMode(int mode)
{
  int mode2 = mode & 7;
  if (mode2 >= 0 && mode2 < 3) {
    // take off bottom bit
    mode_ &= ~15;
    mode_ |= mode;
  }
}
int CglProbing::getMode() const
{
  return mode_ & 15;
}
// Set maximum number of passes per node
void CglProbing::setMaxPass(int value)
{
  if (value > 0)
    maxPass_ = value;
}
// Get maximum number of passes per node
int CglProbing::getMaxPass() const
{
  return maxPass_;
}
// Set log level
void CglProbing::setLogLevel(int value)
{
  if (value >= 0)
    logLevel_ = value;
}
// Get log level
int CglProbing::getLogLevel() const
{
  return logLevel_;
}
// Set maximum number of unsatisfied variables to look at
void CglProbing::setMaxProbe(int value)
{
  if (value >= 0)
    maxProbe_ = value;
}
// Get maximum number of unsatisfied variables to look at
int CglProbing::getMaxProbe() const
{
  return maxProbe_;
}
// Set maximum number of variables to look at in one probe
void CglProbing::setMaxLook(int value)
{
  if (value >= 0)
    maxStack_ = value;
}
// Get maximum number of variables to look at in one probe
int CglProbing::getMaxLook() const
{
  return maxStack_;
}
// Set maximum number of elements in row for scan
void CglProbing::setMaxElements(int value)
{
  if (value > 0)
    maxElements_ = value;
}
// Get maximum number of elements in row for scan
int CglProbing::getMaxElements() const
{
  return maxElements_;
}
// Set maximum number of passes per node (root node)
void CglProbing::setMaxPassRoot(int value)
{
  if (value > 0)
    maxPassRoot_ = value;
}
// Get maximum number of passes per node (root node)
int CglProbing::getMaxPassRoot() const
{
  return maxPassRoot_;
}
// Set maximum number of unsatisfied variables to look at (root node)
void CglProbing::setMaxProbeRoot(int value)
{
  if (value > 0)
    maxProbeRoot_ = value;
}
// Get maximum number of unsatisfied variables to look at (root node)
int CglProbing::getMaxProbeRoot() const
{
  return maxProbeRoot_;
}
// Set maximum number of variables to look at in one probe (root node)
void CglProbing::setMaxLookRoot(int value)
{
  if (value > 0)
    maxStackRoot_ = value;
}
// Get maximum number of variables to look at in one probe (root node)
int CglProbing::getMaxLookRoot() const
{
  return maxStackRoot_;
}
// Set maximum number of elements in row for scan (root node)
void CglProbing::setMaxElementsRoot(int value)
{
  if (value > 0)
    maxElementsRoot_ = value;
}
// Get maximum number of elements in row for scan (root node)
int CglProbing::getMaxElementsRoot() const
{
  return maxElementsRoot_;
}
// Set whether to use objective
void CglProbing::setUsingObjective(int yesNo)
{
  usingObjective_ = yesNo;
}
// Get whether objective is being used
int CglProbing::getUsingObjective() const
{
  return usingObjective_;
}
// Decide whether to do row cuts
void CglProbing::setRowCuts(int type)
{
  if (type > -5 && type < 5)
    rowCuts_ = type;
}
// Returns row cuts generation type
int CglProbing::rowCuts() const
{
  return rowCuts_;
}
// Returns tight lower
const double *CglProbing::tightLower() const
{
  return colLower_;
}
// Returns tight upper
const double *CglProbing::tightUpper() const
{
  return colUpper_;
}
// Returns relaxed Row lower
const double *CglProbing::relaxedRowLower() const
{
  return rowLower_;
}
// Returns relaxed Row upper
const double *CglProbing::relaxedRowUpper() const
{
  return rowUpper_;
}

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
CglProbing::CglProbing()
  : CglCutGenerator()
  , primalTolerance_(1.1e-07)
  , mode_(1)
  , rowCuts_(1)
  , maxPass_(3)
  , logLevel_(0)
  , maxProbe_(100)
  , maxStack_(50)
  , maxElements_(1000)
  , maxPassRoot_(3)
  , maxProbeRoot_(100)
  , maxStackRoot_(50)
  , maxElementsRoot_(10000)
  , usingObjective_(0)
  , maxSeconds_(0.0)
{

  numberRows_ = 0;
  numberColumns_ = 0;
  rowCopy_ = NULL;
  columnCopy_ = NULL;
  rowLower_ = NULL;
  rowUpper_ = NULL;
  colLower_ = NULL;
  colUpper_ = NULL;
  minR_ = NULL;
  maxR_ = NULL;
  markRow_ = NULL;
  markColumn_ = NULL;
  lookRow_[0] = NULL;
  lookRow_[1] = NULL;
  lookColumn_ = NULL;
  numberIntegers_ = 0;
  number01Integers_ = 0;
  numberThisTime_ = 0;
  totalTimesCalled_ = 0;
  lookedAt_ = NULL;
  cutVector_ = NULL;
  numberCliques_ = 0;
  cliqueType_ = NULL;
  cliqueStart_ = NULL;
  cliqueEntry_ = NULL;
  oneFixStart_ = NULL;
  zeroFixStart_ = NULL;
  endFixStart_ = NULL;
  whichClique_ = NULL;
  cliqueRow_ = NULL;
  cliqueRowStart_ = NULL;
  tightenBounds_ = NULL;
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CglProbing::CglProbing(const CglProbing &rhs)
  : CglCutGenerator(rhs)
  , primalTolerance_(rhs.primalTolerance_)
  , mode_(rhs.mode_)
  , rowCuts_(rhs.rowCuts_)
  , maxPass_(rhs.maxPass_)
  , logLevel_(rhs.logLevel_)
  , maxProbe_(rhs.maxProbe_)
  , maxStack_(rhs.maxStack_)
  , maxElements_(rhs.maxElements_)
  , maxPassRoot_(rhs.maxPassRoot_)
  , maxProbeRoot_(rhs.maxProbeRoot_)
  , maxStackRoot_(rhs.maxStackRoot_)
  , maxElementsRoot_(rhs.maxElementsRoot_)
  , usingObjective_(rhs.usingObjective_)
  , maxSeconds_(rhs.maxSeconds_)
{
  numberRows_ = rhs.numberRows_;
  numberColumns_ = rhs.numberColumns_;
  numberCliques_ = rhs.numberCliques_;
  if (rhs.rowCopy_) {
    rowCopy_ = new CoinPackedMatrix(*(rhs.rowCopy_));
    columnCopy_ = new CoinPackedMatrix(*(rhs.columnCopy_));
    rowLower_ = new double[numberRows_];
    CoinMemcpyN(rhs.rowLower_, numberRows_, rowLower_);
    rowUpper_ = new double[numberRows_];
    CoinMemcpyN(rhs.rowUpper_, numberRows_, rowUpper_);
    colLower_ = new double[numberColumns_];
    CoinMemcpyN(rhs.colLower_, numberColumns_, colLower_);
    colUpper_ = new double[numberColumns_];
    CoinMemcpyN(rhs.colUpper_, numberColumns_, colUpper_);
    int i;
    numberIntegers_ = rhs.numberIntegers_;
    number01Integers_ = rhs.number01Integers_;
    cutVector_ = new disaggregation[number01Integers_];
    CoinMemcpyN(rhs.cutVector_, number01Integers_, cutVector_);
    for (i = 0; i < number01Integers_; i++) {
      if (cutVector_[i].index) {
        cutVector_[i].index = CoinCopyOfArray(rhs.cutVector_[i].index, cutVector_[i].length);
      }
    }
  } else {
    rowCopy_ = NULL;
    columnCopy_ = NULL;
    rowLower_ = NULL;
    rowUpper_ = NULL;
    colLower_ = NULL;
    colUpper_ = NULL;
    numberIntegers_ = 0;
    number01Integers_ = 0;
    cutVector_ = NULL;
  }
  numberThisTime_ = rhs.numberThisTime_;
  totalTimesCalled_ = rhs.totalTimesCalled_;
  if (numberColumns_)
    lookedAt_ = CoinCopyOfArray(rhs.lookedAt_, numberColumns_);
  else
    lookedAt_ = NULL;
  if (numberCliques_) {
    cliqueType_ = new CliqueType[numberCliques_];
    CoinMemcpyN(rhs.cliqueType_, numberCliques_, cliqueType_);
    cliqueStart_ = new CoinBigIndex[numberCliques_ + 1];
    CoinMemcpyN(rhs.cliqueStart_, (numberCliques_ + 1), cliqueStart_);
    CoinBigIndex n = cliqueStart_[numberCliques_];
    cliqueEntry_ = new CliqueEntry[n];
    CoinMemcpyN(rhs.cliqueEntry_, n, cliqueEntry_);
    oneFixStart_ = new int[numberColumns_];
    CoinMemcpyN(rhs.oneFixStart_, numberColumns_, oneFixStart_);
    zeroFixStart_ = new int[numberColumns_];
    CoinMemcpyN(rhs.zeroFixStart_, numberColumns_, zeroFixStart_);
    endFixStart_ = new int[numberColumns_];
    CoinMemcpyN(rhs.endFixStart_, numberColumns_, endFixStart_);
#ifndef NDEBUG
    int n2 = -1;
    for (int i = numberColumns_ - 1; i >= 0; i--) {
      if (oneFixStart_[i] >= 0) {
        n2 = endFixStart_[i];
        break;
      }
    }
    assert(n == n2);
#endif
    whichClique_ = new int[n];
    CoinMemcpyN(rhs.whichClique_, n, whichClique_);
    if (rhs.cliqueRowStart_) {
      cliqueRowStart_ = CoinCopyOfArray(rhs.cliqueRowStart_, numberRows_ + 1);
      n = cliqueRowStart_[numberRows_];
      cliqueRow_ = CoinCopyOfArray(rhs.cliqueRow_, n);
    } else {
      cliqueRow_ = NULL;
      cliqueRowStart_ = NULL;
    }
  } else {
    cliqueType_ = NULL;
    cliqueStart_ = NULL;
    cliqueEntry_ = NULL;
    oneFixStart_ = NULL;
    zeroFixStart_ = NULL;
    endFixStart_ = NULL;
    cliqueRow_ = NULL;
    cliqueRowStart_ = NULL;
    whichClique_ = NULL;
  }
  if (rhs.tightenBounds_) {
    assert(numberColumns_);
    tightenBounds_ = CoinCopyOfArray(rhs.tightenBounds_, numberColumns_);
  } else {
    tightenBounds_ = NULL;
  }
  minR_ = NULL;
  maxR_ = NULL;
  markRow_ = NULL;
  markColumn_ = NULL;
  lookRow_[0] = NULL;
  lookRow_[1] = NULL;
  lookColumn_ = NULL;
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglProbing::clone() const
{
  return new CglProbing(*this);
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
CglProbing::~CglProbing()
{
  // free memory
  delete[] rowLower_;
  delete[] rowUpper_;
  delete[] colLower_;
  delete[] colUpper_;
  delete rowCopy_;
  delete columnCopy_;
  delete[] lookedAt_;
  delete[] cliqueType_;
  delete[] cliqueStart_;
  delete[] cliqueEntry_;
  delete[] oneFixStart_;
  delete[] zeroFixStart_;
  delete[] endFixStart_;
  delete[] whichClique_;
  delete[] cliqueRow_;
  delete[] cliqueRowStart_;
  if (cutVector_) {
    for (int i = 0; i < number01Integers_; i++) {
      delete[] cutVector_[i].index;
    }
    delete[] cutVector_;
  }
  delete[] tightenBounds_;
  assert(minR_ == NULL);
  assert(maxR_ == NULL);
  assert(markRow_ == NULL);
  assert(markColumn_ == NULL);
  assert(lookColumn_ == NULL);
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CglProbing &
CglProbing::operator=(
  const CglProbing &rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    primalTolerance_ = rhs.primalTolerance_;
    numberRows_ = rhs.numberRows_;
    numberColumns_ = rhs.numberColumns_;
    delete[] rowLower_;
    delete[] rowUpper_;
    delete[] colLower_;
    delete[] colUpper_;
    delete rowCopy_;
    delete columnCopy_;
    delete[] lookedAt_;
    delete[] cliqueType_;
    delete[] cliqueStart_;
    delete[] cliqueEntry_;
    delete[] oneFixStart_;
    delete[] zeroFixStart_;
    delete[] endFixStart_;
    delete[] whichClique_;
    delete[] cliqueRow_;
    delete[] cliqueRowStart_;
    delete[] tightenBounds_;
    mode_ = rhs.mode_;
    rowCuts_ = rhs.rowCuts_;
    maxPass_ = rhs.maxPass_;
    logLevel_ = rhs.logLevel_;
    maxProbe_ = rhs.maxProbe_;
    maxStack_ = rhs.maxStack_;
    maxElements_ = rhs.maxElements_;
    maxPassRoot_ = rhs.maxPassRoot_;
    maxProbeRoot_ = rhs.maxProbeRoot_;
    maxStackRoot_ = rhs.maxStackRoot_;
    maxElementsRoot_ = rhs.maxElementsRoot_;
    usingObjective_ = rhs.usingObjective_;
    maxSeconds_ = rhs.maxSeconds_;
    numberCliques_ = rhs.numberCliques_;
    assert(minR_ == NULL);
    assert(maxR_ == NULL);
    assert(markRow_ == NULL);
    assert(markColumn_ == NULL);
    assert(lookColumn_ == NULL);
    if (rhs.rowCopy_) {
      rowCopy_ = new CoinPackedMatrix(*(rhs.rowCopy_));
      columnCopy_ = new CoinPackedMatrix(*(rhs.columnCopy_));
      rowLower_ = new double[numberRows_];
      CoinMemcpyN(rhs.rowLower_, numberRows_, rowLower_);
      rowUpper_ = new double[numberRows_];
      CoinMemcpyN(rhs.rowUpper_, numberRows_, rowUpper_);
      colLower_ = new double[numberColumns_];
      CoinMemcpyN(rhs.colLower_, numberColumns_, colLower_);
      colUpper_ = new double[numberColumns_];
      CoinMemcpyN(rhs.colUpper_, numberColumns_, colUpper_);
      int i;
      numberIntegers_ = rhs.numberIntegers_;
      number01Integers_ = rhs.number01Integers_;
      for (i = 0; i < number01Integers_; i++) {
        delete[] cutVector_[i].index;
      }
      delete[] cutVector_;
      cutVector_ = new disaggregation[number01Integers_];
      CoinMemcpyN(rhs.cutVector_, number01Integers_, cutVector_);
      for (i = 0; i < number01Integers_; i++) {
        if (cutVector_[i].index) {
          cutVector_[i].index = CoinCopyOfArray(rhs.cutVector_[i].index, cutVector_[i].length);
        }
      }
    } else {
      rowCopy_ = NULL;
      columnCopy_ = NULL;
      rowLower_ = NULL;
      rowUpper_ = NULL;
      colLower_ = NULL;
      colUpper_ = NULL;
      numberIntegers_ = 0;
      number01Integers_ = 0;
      cutVector_ = NULL;
    }
    numberThisTime_ = rhs.numberThisTime_;
    totalTimesCalled_ = rhs.totalTimesCalled_;
    if (numberColumns_)
      lookedAt_ = CoinCopyOfArray(rhs.lookedAt_, numberColumns_);
    else
      lookedAt_ = NULL;
    if (numberCliques_) {
      cliqueType_ = new CliqueType[numberCliques_];
      CoinMemcpyN(rhs.cliqueType_, numberCliques_, cliqueType_);
      cliqueStart_ = new CoinBigIndex[numberCliques_ + 1];
      CoinMemcpyN(rhs.cliqueStart_, (numberCliques_ + 1), cliqueStart_);
      CoinBigIndex n = cliqueStart_[numberCliques_];
      cliqueEntry_ = new CliqueEntry[n];
      CoinMemcpyN(rhs.cliqueEntry_, n, cliqueEntry_);
      oneFixStart_ = new int[numberColumns_];
      CoinMemcpyN(rhs.oneFixStart_, numberColumns_, oneFixStart_);
      zeroFixStart_ = new int[numberColumns_];
      CoinMemcpyN(rhs.zeroFixStart_, numberColumns_, zeroFixStart_);
      endFixStart_ = new int[numberColumns_];
      CoinMemcpyN(rhs.endFixStart_, numberColumns_, endFixStart_);
#ifndef NDEBUG
      int n2 = -1;
      for (int i = numberColumns_ - 1; i >= 0; i--) {
        if (oneFixStart_[i] >= 0) {
          n2 = endFixStart_[i];
          break;
        }
      }
      assert(n == n2);
#endif
      whichClique_ = new int[n];
      CoinMemcpyN(rhs.whichClique_, n, whichClique_);
      if (rhs.cliqueRowStart_) {
        cliqueRowStart_ = CoinCopyOfArray(rhs.cliqueRowStart_, numberRows_ + 1);
        n = cliqueRowStart_[numberRows_];
        cliqueRow_ = CoinCopyOfArray(rhs.cliqueRow_, n);
      } else {
        cliqueRow_ = NULL;
        cliqueRowStart_ = NULL;
      }
    } else {
      cliqueType_ = NULL;
      cliqueStart_ = NULL;
      cliqueEntry_ = NULL;
      oneFixStart_ = NULL;
      zeroFixStart_ = NULL;
      endFixStart_ = NULL;
      whichClique_ = NULL;
      cliqueRow_ = NULL;
      cliqueRowStart_ = NULL;
    }
    if (rhs.tightenBounds_) {
      assert(numberColumns_);
      tightenBounds_ = CoinCopyOfArray(rhs.tightenBounds_, numberColumns_);
    } else {
      tightenBounds_ = NULL;
    }
  }
  return *this;
}

/// This can be used to refresh any inforamtion
void CglProbing::refreshSolver(OsiSolverInterface *solver)
{
  if (rowCopy_) {
    // snapshot existed - redo
    snapshot(*solver, NULL);
  }
  // Get integer information
  solver->getColType(true);
}
/* Creates cliques for use by probing.
   Can also try and extend cliques as a result of probing (root node).
   Returns number of cliques found.
*/
int CglProbing::createCliques(OsiSolverInterface &si,
  int minimumSize, int maximumSize)
{
  // get rid of what is there
  deleteCliques();
  CoinPackedMatrix matrixByRow(*si.getMatrixByRow());
  int numberRows = si.getNumRows();
  if (!rowCopy_)
    numberRows_ = numberRows;
  numberColumns_ = si.getNumCols();
  const char *intVar = si.getColType();
  const double *lower = si.getColLower();
  const double *upper = si.getColUpper();

  numberCliques_ = 0;
  int numberEntries = 0;
  int numberIntegers = 0;
  int *lookup = new int[numberColumns_];
  int i;
  for (i = 0; i < numberColumns_; i++) {
    if (intVar[i] == 1)
      lookup[i] = numberIntegers++;
    else
      lookup[i] = -1;
  }

  int *which = new int[numberColumns_];
  int *whichRow = new int[numberRows];
  // Statistics
  int totalP1 = 0, totalM1 = 0;
  int numberBig = 0, totalBig = 0;
  int numberFixed = 0;

  // Row copy
  const double *elementByRow = matrixByRow.getElements();
  const int *column = matrixByRow.getIndices();
  const CoinBigIndex *rowStart = matrixByRow.getVectorStarts();
  const int *rowLength = matrixByRow.getVectorLengths();

  const double *rowLower = si.getRowLower();
  const double *rowUpper = si.getRowUpper();
  int iRow;
  for (iRow = 0; iRow < numberRows; iRow++) {
    int numberP1 = 0, numberM1 = 0;
    CoinBigIndex j;
    double upperValue = rowUpper[iRow];
    double lowerValue = rowLower[iRow];
    bool good = true;
    for (j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
      int iColumn = column[j];
      int iInteger = lookup[iColumn];
      if (upper[iColumn] - lower[iColumn] < 1.0e-8) {
        // fixed
        upperValue -= lower[iColumn] * elementByRow[j];
        lowerValue -= lower[iColumn] * elementByRow[j];
        continue;
      } else if (upper[iColumn] != 1.0 || lower[iColumn] != 0.0) {
        good = false;
        break;
      } else if (iInteger < 0) {
        good = false;
        break;
      }
      if (fabs(elementByRow[j]) != 1.0) {
        good = false;
        break;
      } else if (elementByRow[j] > 0.0) {
        which[numberP1++] = iColumn;
      } else {
        numberM1++;
        which[numberIntegers - numberM1] = iColumn;
      }
    }
    int iUpper = upperValue > INT_MAX ? INT_MAX : static_cast< int >(floor(upperValue + 1.0e-5));
    int iLower = lowerValue < INT_MIN ? INT_MIN : static_cast< int >(ceil(lowerValue - 1.0e-5));
    int state = 0;
    if (upperValue < 1.0e6) {
      if (iUpper == 1 - numberM1)
        state = 1;
      else if (iUpper == -numberM1)
        state = 2;
      else if (iUpper < -numberM1)
        state = 3;
    }
    if (!state && lowerValue > -1.0e6) {
      if (-iLower == 1 - numberP1)
        state = -1;
      else if (-iLower == -numberP1)
        state = -2;
      else if (-iLower < -numberP1)
        state = -3;
    }
    if (good && state) {
      if (abs(state) == 3) {
        // infeasible
        numberCliques_ = -99999;
        break;
      } else if (abs(state) == 2) {
        // we can fix all
        numberFixed += numberP1 + numberM1;
        if (state > 0) {
          // fix all +1 at 0, -1 at 1
          for (i = 0; i < numberP1; i++)
            si.setColUpper(which[i], 0.0);
          for (i = 0; i < numberM1; i++)
            si.setColLower(which[numberIntegers - i - 1],
              1.0);
        } else {
          // fix all +1 at 1, -1 at 0
          for (i = 0; i < numberP1; i++)
            si.setColLower(which[i], 1.0);
          for (i = 0; i < numberM1; i++)
            si.setColUpper(which[numberIntegers - i - 1],
              0.0);
        }
      } else {
        int length = numberP1 + numberM1;
        totalP1 += numberP1;
        totalM1 += numberM1;
        if (length >= minimumSize && length < maximumSize) {
          whichRow[numberCliques_++] = iRow;
          numberEntries += length;
        } else if (numberP1 + numberM1 >= maximumSize) {
          // too big
          numberBig++;
          totalBig += numberP1 + numberM1;
        }
      }
    }
  }
  if (numberCliques_ < 0) {
    if (logLevel_)
      printf("*** Problem infeasible\n");
  } else {
    if (numberCliques_) {
      if (logLevel_)
        printf("%d cliques of average size %g found, %d P1, %d M1\n",
          numberCliques_,
          (static_cast< double >(totalP1 + totalM1)) / (static_cast< double >(numberCliques_)),
          totalP1, totalM1);
    } else {
      if (logLevel_ > 1)
        printf("No cliques found\n");
    }
    if (numberBig) {
      if (logLevel_)
        printf("%d large cliques ( >= %d) found, total %d\n",
          numberBig, maximumSize, totalBig);
    }
    if (numberFixed) {
      if (logLevel_)
        printf("%d variables fixed\n", numberFixed);
    }
  }
  if (numberCliques_ > 0) {
    cliqueType_ = new CliqueType[numberCliques_];
    cliqueStart_ = new CoinBigIndex[numberCliques_ + 1];
    cliqueEntry_ = new CliqueEntry[numberEntries];
    oneFixStart_ = new int[numberColumns_];
    zeroFixStart_ = new int[numberColumns_];
    endFixStart_ = new int[numberColumns_];
    whichClique_ = new int[numberEntries];
    numberEntries = 0;
    cliqueStart_[0] = 0;
    for (i = 0; i < numberColumns_; i++) {
      oneFixStart_[i] = -1;
      zeroFixStart_[i] = -1;
      endFixStart_[i] = -1;
    }
    int iClique;
    // Possible some have been fixed
    int numberCliques = numberCliques_;
    numberCliques_ = 0;
    for (iClique = 0; iClique < numberCliques; iClique++) {
      int iRow = whichRow[iClique];
      whichRow[numberCliques_] = iRow;
      int numberP1 = 0, numberM1 = 0;
      CoinBigIndex j;
      double upperValue = rowUpper[iRow];
      double lowerValue = rowLower[iRow];
      for (j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
        int iColumn = column[j];
        if (upper[iColumn] - lower[iColumn] < 1.0e-8) {
          // fixed
          upperValue -= lower[iColumn] * elementByRow[j];
          lowerValue -= lower[iColumn] * elementByRow[j];
          continue;
        }
        if (elementByRow[j] > 0.0) {
          which[numberP1++] = iColumn;
        } else {
          numberM1++;
          which[numberIntegers - numberM1] = iColumn;
        }
      }
      int iUpper = upperValue > INT_MAX ? INT_MAX : static_cast< int >(floor(upperValue + 1.0e-5));
      int iLower = lowerValue < INT_MIN ? INT_MIN : static_cast< int >(ceil(lowerValue - 1.0e-5));
      int state = 0;
      if (upperValue < 1.0e6) {
        if (iUpper == 1 - numberM1)
          state = 1;
      }
      if (!state && lowerValue > -1.0e6) {
        state = -1;
      }
      if (abs(state) != 1)
        continue; // must have been fixed
      if (iLower == iUpper) {
        cliqueType_[numberCliques_].equality = 1;
      } else {
        cliqueType_[numberCliques_].equality = 0;
      }
      if (state > 0) {
        for (i = 0; i < numberP1; i++) {
          // 1 is strong branch
          int iColumn = which[i];
          setSequenceInCliqueEntry(cliqueEntry_[numberEntries], iColumn);
          setOneFixesInCliqueEntry(cliqueEntry_[numberEntries], true);
          numberEntries++;
          // zero counts
          oneFixStart_[iColumn] = 0;
          zeroFixStart_[iColumn] = 0;
        }
        for (i = 0; i < numberM1; i++) {
          // 0 is strong branch
          int iColumn = which[numberIntegers - i - 1];
          setSequenceInCliqueEntry(cliqueEntry_[numberEntries], iColumn);
          setOneFixesInCliqueEntry(cliqueEntry_[numberEntries], false);
          numberEntries++;
          // zero counts
          oneFixStart_[iColumn] = 0;
          zeroFixStart_[iColumn] = 0;
        }
      } else {
        for (i = 0; i < numberP1; i++) {
          // 0 is strong branch
          int iColumn = which[i];
          setSequenceInCliqueEntry(cliqueEntry_[numberEntries], iColumn);
          setOneFixesInCliqueEntry(cliqueEntry_[numberEntries], false);
          numberEntries++;
          // zero counts
          oneFixStart_[iColumn] = 0;
          zeroFixStart_[iColumn] = 0;
        }
        for (i = 0; i < numberM1; i++) {
          // 1 is strong branch
          int iColumn = which[numberIntegers - i - 1];
          setSequenceInCliqueEntry(cliqueEntry_[numberEntries], iColumn);
          setOneFixesInCliqueEntry(cliqueEntry_[numberEntries], true);
          numberEntries++;
          // zero counts
          oneFixStart_[iColumn] = 0;
          zeroFixStart_[iColumn] = 0;
        }
      }
      numberCliques_++;
      cliqueStart_[numberCliques_] = numberEntries;
    }
    // Now do column lists
    // First do counts
    for (iClique = 0; iClique < numberCliques_; iClique++) {
      for (CoinBigIndex j = cliqueStart_[iClique]; j < cliqueStart_[iClique + 1]; j++) {
        int iColumn = sequenceInCliqueEntry(cliqueEntry_[j]);
        if (oneFixesInCliqueEntry(cliqueEntry_[j]))
          oneFixStart_[iColumn]++;
        else
          zeroFixStart_[iColumn]++;
      }
    }
    // now get starts and use which and end as counters
    numberEntries = 0;
    for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
      if (oneFixStart_[iColumn] >= 0) {
        int n1 = oneFixStart_[iColumn];
        int n2 = zeroFixStart_[iColumn];
        oneFixStart_[iColumn] = numberEntries;
        which[iColumn] = numberEntries;
        numberEntries += n1;
        zeroFixStart_[iColumn] = numberEntries;
        endFixStart_[iColumn] = numberEntries;
        numberEntries += n2;
      }
    }
    // now put in
    for (iClique = 0; iClique < numberCliques_; iClique++) {
      for (CoinBigIndex j = cliqueStart_[iClique]; j < cliqueStart_[iClique + 1]; j++) {
        int iColumn = sequenceInCliqueEntry(cliqueEntry_[j]);
        if (oneFixesInCliqueEntry(cliqueEntry_[j])) {
          int put = which[iColumn];
          which[iColumn]++;
          whichClique_[put] = iClique;
        } else {
          int put = endFixStart_[iColumn];
          endFixStart_[iColumn]++;
          whichClique_[put] = iClique;
        }
      }
    }
  }
  delete[] which;
  delete[] whichRow;
  delete[] lookup;
  return numberCliques_;
}
// Delete all clique information
void CglProbing::deleteCliques()
{
  delete[] cliqueType_;
  delete[] cliqueStart_;
  delete[] cliqueEntry_;
  delete[] oneFixStart_;
  delete[] zeroFixStart_;
  delete[] endFixStart_;
  delete[] whichClique_;
  delete[] cliqueRow_;
  delete[] cliqueRowStart_;
  cliqueType_ = NULL;
  cliqueStart_ = NULL;
  cliqueEntry_ = NULL;
  oneFixStart_ = NULL;
  zeroFixStart_ = NULL;
  endFixStart_ = NULL;
  whichClique_ = NULL;
  cliqueRow_ = NULL;
  cliqueRowStart_ = NULL;
  numberCliques_ = 0;
}
/* Create a fake model by adding cliques
   if type&4 then delete rest of model first,
   if 1 then add proper cliques, 2 add fake cliques */
OsiSolverInterface *
CglProbing::cliqueModel(const OsiSolverInterface *model,
  int type)
{
  OsiSolverInterface *newModel = model->clone();
  if ((type & 4) != 0) {
    int numberRows = newModel->getNumRows();
    int *which = new int[numberRows];
    for (int i = 0; i < numberRows; i++)
      which[i] = i;
    newModel->deleteRows(numberRows, which);
    delete[] which;
  }
  type &= 3;
  if (type != 2) {
    CoinBigIndex numberElements = cliqueStart_[numberCliques_];
    int *column = new int[numberElements];
    double *element = new double[numberElements];
    double *rowLower = new double[numberCliques_];
    double *rowUpper = new double[numberCliques_];
    numberElements = 0;
    for (int iClique = 0; iClique < numberCliques_; iClique++) {
      int nMinus = 0;
      for (CoinBigIndex j = cliqueStart_[iClique]; j < cliqueStart_[iClique + 1]; j++) {
        int iColumn = sequenceInCliqueEntry(cliqueEntry_[j]);
        column[numberElements] = iColumn;
        double value;
        if (oneFixesInCliqueEntry(cliqueEntry_[j])) {
          value = 1.0;
        } else {
          value = -1.0;
          nMinus++;
        }
        element[numberElements++] = value;
      }
      rowLower[iClique] = -COIN_DBL_MAX;
      rowUpper[iClique] = 1 - nMinus;
    }
    newModel->addRows(numberCliques_, cliqueStart_, column, element,
      rowLower, rowUpper);
    delete[] column;
    delete[] element;
    delete[] rowLower;
    delete[] rowUpper;
  } else {
    int numberCuts = numberCliques_;
    CoinBigIndex numberElements = 2 * numberCuts;
    int *column = new int[numberElements];
    double *element = new double[numberElements];
    double *rowLower = new double[numberCuts];
    double *rowUpper = new double[numberCuts];
    CoinBigIndex *start = new CoinBigIndex[numberCuts + 1];
    start[0] = 0;
    numberElements = 0;
    int *entry = reinterpret_cast< int * >(cliqueEntry_);
    for (int iCut = 0; iCut < numberCuts; iCut++) {
      column[numberElements] = entry[numberElements];
      element[numberElements++] = 1.0;
      column[numberElements] = entry[numberElements];
      element[numberElements++] = 1.0;
      rowLower[iCut] = 0.0;
      rowUpper[iCut] = 1.0;
      start[iCut + 1] = numberElements;
    }
    newModel->addRows(numberCuts, start, column, element,
      rowLower, rowUpper);
    delete[] column;
    delete[] element;
    delete[] rowLower;
    delete[] rowUpper;
    delete[] start;
  }
  return newModel;
}
/*
  Returns true if may generate Row cuts in tree (rather than root node).
  Used so know if matrix will change in tree.  Really
  meant so column cut generators can still be active
  without worrying code.
  Default is true
*/
bool CglProbing::mayGenerateRowCutsInTree() const
{
  return rowCuts_ > 0;
}
// Sets up clique information for each row
void CglProbing::setupRowCliqueInformation(const OsiSolverInterface &si)
{
  if (!numberCliques_)
    return;
  CoinPackedMatrix *rowCopy;
  if (!rowCopy_) {
    // create from current
    numberRows_ = si.getNumRows();
    numberColumns_ = si.getNumCols();
    rowCopy = new CoinPackedMatrix(*si.getMatrixByRow());
  } else {
    rowCopy = rowCopy_;
    assert(numberRows_ <= si.getNumRows());
    assert(numberColumns_ == si.getNumCols());
  }
  assert(numberRows_ && numberColumns_);
  cliqueRowStart_ = new int[numberRows_ + 1];
  cliqueRowStart_[0] = 0;
  // Temporary array while building list
  CliqueEntry **array = new CliqueEntry *[numberRows_];
  // Which cliques in use
  int *which = new int[numberCliques_];
  int *count = new int[numberCliques_];
  int *back = new int[numberColumns_];
  CoinZeroN(count, numberCliques_);
  CoinFillN(back, numberColumns_, -1);
  const int *column = rowCopy->getIndices();
  const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
  const int *rowLength = rowCopy->getVectorLengths();
  const double *lower = si.getColLower();
  const double *upper = si.getColUpper();
  int iRow;
  for (iRow = 0; iRow < numberRows_; iRow++) {
    CoinBigIndex j;
    int numberFree = 0;
    int numberUsed = 0;
    for (j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
      int iColumn = column[j];
      if (upper[iColumn] > lower[iColumn]) {
        back[iColumn] = static_cast< int >(j - rowStart[iRow]);
        numberFree++;
        for (int k = oneFixStart_[iColumn]; k < endFixStart_[iColumn]; k++) {
          int iClique = whichClique_[k];
          if (!count[iClique]) {
            which[numberUsed++] = iClique;
          }
          count[iClique]++;
        }
      }
    }
    // find largest cliques
    bool finished = false;
    int numberInThis = 0;
    CliqueEntry *entries = NULL;
    array[iRow] = entries;
    while (!finished) {
      int largest = 1;
      int whichClique = -1;
      for (int i = 0; i < numberUsed; i++) {
        int iClique = which[i];
        if (count[iClique] > largest) {
          largest = count[iClique];
          whichClique = iClique;
        }
      }
      // Add in if >1 (but not if all as that means clique==row)
      if (whichClique >= 0 && largest < numberFree) {
        if (!numberInThis) {
          int length = rowLength[iRow];
          entries = new CliqueEntry[length];
          array[iRow] = entries;
          for (int i = 0; i < length; i++) {
            setOneFixesInCliqueEntry(entries[i], false);
            setSequenceInCliqueEntry(entries[i], numberColumns_ + 1);
          }
        }
        // put in (and take out all counts)
        for (j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
          int iColumn = column[j];
          if (upper[iColumn] > lower[iColumn]) {
            bool found = false;
            CoinBigIndex k;
            for (k = oneFixStart_[iColumn]; k < endFixStart_[iColumn]; k++) {
              int iClique = whichClique_[k];
              if (iClique == whichClique) {
                found = true;
                break;
              }
            }
            if (found) {
              for (k = oneFixStart_[iColumn]; k < endFixStart_[iColumn]; k++) {
                int iClique = whichClique_[k];
                count[iClique]--;
              }
              for (k = cliqueStart_[whichClique]; k < cliqueStart_[whichClique + 1]; k++) {
                if (sequenceInCliqueEntry(cliqueEntry_[k]) == iColumn) {
                  int iback = back[iColumn];
                  setSequenceInCliqueEntry(entries[iback], numberInThis);
                  setOneFixesInCliqueEntry(entries[iback],
                    oneFixesInCliqueEntry(cliqueEntry_[k]));
                  break;
                }
              }
            }
          }
        }
        numberInThis++;
      } else {
        finished = true;
      }
    }
    if (numberInThis)
      cliqueRowStart_[iRow + 1] = cliqueRowStart_[iRow] + rowLength[iRow];
    else
      cliqueRowStart_[iRow + 1] = cliqueRowStart_[iRow];
    for (int i = 0; i < numberUsed; i++) {
      int iClique = which[i];
      count[iClique] = 0;
    }
    for (j = rowStart[iRow]; j < rowStart[iRow] + rowLength[iRow]; j++) {
      int iColumn = column[j];
      back[iColumn] = -1;
    }
  }
  delete[] which;
  delete[] count;
  delete[] back;
  // Now put info in one array
  cliqueRow_ = new CliqueEntry[cliqueRowStart_[numberRows_]];
  for (iRow = 0; iRow < numberRows_; iRow++) {
    if (array[iRow]) {
      int start = cliqueRowStart_[iRow];
      CoinMemcpyN(array[iRow], rowLength[iRow], cliqueRow_ + start);
      delete[] array[iRow];
    }
  }
  delete[] array;
  if (rowCopy != rowCopy_)
    delete rowCopy;
}
// Mark variables to be tightened
void CglProbing::tightenThese(const OsiSolverInterface &solver, int number, const int *which)
{
  delete[] tightenBounds_;
  int numberColumns = solver.getNumCols();
  if (numberColumns_)
    assert(numberColumns_ == numberColumns);
  tightenBounds_ = new char[numberColumns];
  memset(tightenBounds_, 0, numberColumns);
  for (int i = 0; i < number; i++) {
    int k = which[i];
    if (k >= 0 && k < numberColumns)
      tightenBounds_[k] = 1;
  }
}
// Create C++ lines to get to current state
std::string
CglProbing::generateCpp(FILE *fp)
{
  CglProbing other;
  fprintf(fp, "0#include \"CglProbing.hpp\"\n");
  fprintf(fp, "3  CglProbing probing;\n");
  if (getMode() != other.getMode())
    fprintf(fp, "3  probing.setMode(%d);\n", getMode());
  else
    fprintf(fp, "4  probing.setMode(%d);\n", getMode());
  if (getMaxPass() != other.getMaxPass())
    fprintf(fp, "3  probing.setMaxPass(%d);\n", getMaxPass());
  else
    fprintf(fp, "4  probing.setMaxPass(%d);\n", getMaxPass());
  if (getLogLevel() != other.getLogLevel())
    fprintf(fp, "3  probing.setLogLevel(%d);\n", getLogLevel());
  else
    fprintf(fp, "4  probing.setLogLevel(%d);\n", getLogLevel());
  if (getMaxProbe() != other.getMaxProbe())
    fprintf(fp, "3  probing.setMaxProbe(%d);\n", getMaxProbe());
  else
    fprintf(fp, "4  probing.setMaxProbe(%d);\n", getMaxProbe());
  if (getMaxLook() != other.getMaxLook())
    fprintf(fp, "3  probing.setMaxLook(%d);\n", getMaxLook());
  else
    fprintf(fp, "4  probing.setMaxLook(%d);\n", getMaxLook());
  if (getMaxElements() != other.getMaxElements())
    fprintf(fp, "3  probing.setMaxElements(%d);\n", getMaxElements());
  else
    fprintf(fp, "4  probing.setMaxElements(%d);\n", getMaxElements());
  if (getMaxPassRoot() != other.getMaxPassRoot())
    fprintf(fp, "3  probing.setMaxPassRoot(%d);\n", getMaxPassRoot());
  else
    fprintf(fp, "4  probing.setMaxPassRoot(%d);\n", getMaxPassRoot());
  if (getMaxProbeRoot() != other.getMaxProbeRoot())
    fprintf(fp, "3  probing.setMaxProbeRoot(%d);\n", getMaxProbeRoot());
  else
    fprintf(fp, "4  probing.setMaxProbeRoot(%d);\n", getMaxProbeRoot());
  if (getMaxLookRoot() != other.getMaxLookRoot())
    fprintf(fp, "3  probing.setMaxLookRoot(%d);\n", getMaxLookRoot());
  else
    fprintf(fp, "4  probing.setMaxLookRoot(%d);\n", getMaxLookRoot());
  if (getMaxElementsRoot() != other.getMaxElementsRoot())
    fprintf(fp, "3  probing.setMaxElementsRoot(%d);\n", getMaxElementsRoot());
  else
    fprintf(fp, "4  probing.setMaxElementsRoot(%d);\n", getMaxElementsRoot());
  if (rowCuts() != other.rowCuts())
    fprintf(fp, "3  probing.setRowCuts(%d);\n", rowCuts());
  else
    fprintf(fp, "4  probing.setRowCuts(%d);\n", rowCuts());
  if (getUsingObjective() != other.getUsingObjective())
    fprintf(fp, "3  probing.setUsingObjective(%d);\n", getUsingObjective());
  else
    fprintf(fp, "4  probing.setUsingObjective(%d);\n", getUsingObjective());
  if (getAggressiveness() != other.getAggressiveness())
    fprintf(fp, "3  probing.setAggressiveness(%d);\n", getAggressiveness());
  else
    fprintf(fp, "4  probing.setAggressiveness(%d);\n", getAggressiveness());
  return "probing";
}
//-------------------------------------------------------------
void CglImplication::generateCuts(const OsiSolverInterface &si, OsiCuts &cs,
  const CglTreeInfo info)
{
  if (probingInfo_) {
    // int n1=cs.sizeRowCuts();
    probingInfo_->generateCuts(si, cs, info);
    // int n2=cs.sizeRowCuts();
    // if (n2>n1)
    // printf("added %d cuts\n",n2-n1);
  }
}

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
CglImplication::CglImplication()
  : CglCutGenerator()
  , probingInfo_(NULL)
{
  // nothing to do here
}
//-------------------------------------------------------------------
// Constructor with info
//-------------------------------------------------------------------
CglImplication::CglImplication(CglTreeProbingInfo *info)
  : CglCutGenerator()
  , probingInfo_(info)
{
  // nothing to do here
}
//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CglImplication::CglImplication(
  const CglImplication &source)
  : CglCutGenerator(source)
  , probingInfo_(source.probingInfo_)
{
  // Nothing to do here
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CglCutGenerator *
CglImplication::clone() const
{
  return new CglImplication(*this);
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
CglImplication::~CglImplication()
{
  // Nothing to do here
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CglImplication &
CglImplication::operator=(
  const CglImplication &rhs)
{
  if (this != &rhs) {
    CglCutGenerator::operator=(rhs);
    probingInfo_ = rhs.probingInfo_;
  }
  return *this;
}
// Create C++ lines to get to current state
std::string
CglImplication::generateCpp(FILE *fp)
{
  CglImplication other;
  fprintf(fp, "0#include \"CglImplication.hpp\"\n");
  fprintf(fp, "3  CglImplication implication;\n");
  return "implication";
}
// Start FAST
// Does probing but no cuts
int CglProbing::probeFast(const OsiSolverInterface &si,
  const OsiRowCutDebugger *
#ifdef CGL_DEBUG
    debugger
#endif
  ,
  OsiCuts &cs,
  double *COIN_RESTRICT colLower, double *COIN_RESTRICT colUpper,
  CoinPackedMatrix *rowCopy,
  CoinPackedMatrix *columnCopy,
  const CoinBigIndex *COIN_RESTRICT rowStartPos, const int *COIN_RESTRICT realRows,
  const double *COIN_RESTRICT rowLower, const double *COIN_RESTRICT rowUpper,
  const char *COIN_RESTRICT intVar, double *COIN_RESTRICT minR, double *COIN_RESTRICT maxR,
  int *COIN_RESTRICT markR,
  CglTreeInfo *info)
{
  int nRows = rowCopy->getNumRows();
  int nRowsSafe = std::min(nRows, si.getNumRows());
  int nCols = rowCopy->getNumCols();
  const double *COIN_RESTRICT currentColLower = si.getColLower();
  const double *COIN_RESTRICT currentColUpper = si.getColUpper();
  // Set up maxes
  int maxPass = maxPass_;
  // need a way for user to ask for more
  int maxStack = 80;
  double relaxedTolerance = 2.0 * primalTolerance_;
#define ONE_ARRAY
#ifdef ONE_ARRAY
  unsigned int DIratio = sizeof(double) / sizeof(int);
  assert(DIratio == 1 || DIratio == 2);
  int nSpace = 10 * nCols + 4 * nRows + 2 * maxStack;
  nSpace += (6 * nCols + nRows + DIratio - 1) >> (DIratio - 1);
  // #define UPDATE_MINR_MAXR
#undef UPDATE_MINR_MAXR
#ifdef UPDATE_MINR_MAXR
  // If works try and use existing arrays
  nSpace += 2 * maxStack;
  nSpace += maxStack >> (DIratio - 1);
#endif
#if FIXED_BOTH_WAYS
  nSpace += 2 * nCols;
  nSpace += nCols >> (DIratio - 1);
#endif
  double *COIN_RESTRICT colsol = new double[nSpace];
  double *COIN_RESTRICT djs = colsol + nCols;
  double *COIN_RESTRICT columnGap = djs + nCols;
  double *COIN_RESTRICT saveL = columnGap + nCols;
  double *COIN_RESTRICT saveU = saveL + 2 * nCols;
  double *COIN_RESTRICT saveMin = saveU + 2 * nCols;
  double *COIN_RESTRICT saveMax = saveMin + nRows;
  double *COIN_RESTRICT largestPositiveInRow = saveMax + nRows;
  double *COIN_RESTRICT largestNegativeInRow = largestPositiveInRow + nRows;
  double *COIN_RESTRICT element = largestNegativeInRow + nRows;
  double *COIN_RESTRICT lo0 = element + nCols;
  double *COIN_RESTRICT up0 = lo0 + maxStack;
  double *COIN_RESTRICT newLo = up0 + maxStack;
  double *COIN_RESTRICT newUp = newLo + nCols;
#if FIXED_BOTH_WAYS
  bool doFixedBothWays = (totalTimesCalled_ == 0 && endFixStart_
    && (info->options & 2048) != 0);
  if (doFixedBothWays) {
    for (int i = 0; i < nCols; i++) {
      newLo[i] = COIN_DBL_MAX;
      newUp[i] = -COIN_DBL_MAX;
    }
  }
  int nFixedBoth = 0;
#if FIXED_BOTH_WAYS > 2
  int nToBoundZ = 0;
  int nCutZ = 0;
#endif
#endif
#ifndef UPDATE_MINR_MAXR
  int *COIN_RESTRICT markC = reinterpret_cast< int * COIN_RESTRICT >(newUp + nCols);
#else
  // If works try and use existing arrays
  // If + then increase in lower, negative decrease in upper
  double *COIN_RESTRICT realChange = newUp + nCols;
  int *COIN_RESTRICT markC = reinterpret_cast< int * COIN_RESTRICT >(realChange + nCols);
#endif
  int *COIN_RESTRICT stackC = markC + nCols;
  int *COIN_RESTRICT stackR = stackC + 2 * nCols;
  int *COIN_RESTRICT index = stackR + nRows;
  int *COIN_RESTRICT stackC0 = index + nCols;
#if FIXED_BOTH_WAYS
  int *COIN_RESTRICT fixedBoth = stackC0 + nCols;
  // Really must clean up arrays
  // Make more efficient when two versions of tightenPrimalBounds
  double *COIN_RESTRICT saveFColLower = new double[2 * nCols + 2 * nRows];
  double *COIN_RESTRICT saveFColUpper = saveFColLower + nCols;
  double *COIN_RESTRICT saveFRowLower = saveFColUpper + nCols;
  double *COIN_RESTRICT saveFRowUpper = saveFRowLower + nRows;
  if (!doFixedBothWays)
    fixedBoth = NULL;
#endif
#ifdef UPDATE_MINR_MAXR
  int *COIN_RESTRICT whichRealChange = stackC0 + maxStack;
#if FIXED_BOTH_WAYS
  if (doFixedBothWays)
    whichRealChange += 2 * nCols;
#endif
  assert(whichRealChange + 2 * maxStack < ((int *)(colsol + nSpace)));
  ;
#endif
#else
  double *COIN_RESTRICT colsol = new double[nCols];
  double *COIN_RESTRICT djs = new double[nCols];
  double *COIN_RESTRICT columnGap = new double[nCols];
  double *COIN_RESTRICT saveL = new double[2 * nCols];
  double *COIN_RESTRICT saveU = new double[2 * nCols];
  double *COIN_RESTRICT saveMin = new double[nRows];
  double *COIN_RESTRICT saveMax = new double[nRows];
  double *COIN_RESTRICT largestPositiveInRow = new double[nRows];
  double *COIN_RESTRICT largestNegativeInRow = new double[nRows];
  double *COIN_RESTRICT element = new double[nCols];
  double *COIN_RESTRICT lo0 = new double[maxStack];
  double *COIN_RESTRICT up0 = new double[maxStack];
  int *COIN_RESTRICT markC = new int[nCols];
  int *COIN_RESTRICT stackC = new int[2 * nCols];
  int *COIN_RESTRICT stackR = new int[nRows];
  int *COIN_RESTRICT index = new int[nCols];
  int *COIN_RESTRICT stackC0 = new int[maxStack];
#endif
#ifdef UPDATE_MINR_MAXR
  int numberRealChanges = 0;
#endif
  int *COIN_RESTRICT changedDown = endFixStart_;
  int *COIN_RESTRICT changedUp = changedDown ? changedDown + nCols : NULL;
  totalTimesCalled_++;
  const int *COIN_RESTRICT column = rowCopy->getIndices();
  const CoinBigIndex *COIN_RESTRICT rowStart = rowCopy->getVectorStarts();
  const double *COIN_RESTRICT rowElements = rowCopy->getElements();
  const int *COIN_RESTRICT row = columnCopy->getIndices();
  const CoinBigIndex *COIN_RESTRICT columnStart = columnCopy->getVectorStarts();
  const int *COIN_RESTRICT columnLength = columnCopy->getVectorLengths();
  const double *COIN_RESTRICT columnElements = columnCopy->getElements();
  bool anyColumnCuts = false;
  int ninfeas = 0;
  /* clean up djs and solution */
  CoinMemcpyN(si.getReducedCost(), nCols, djs);
  CoinMemcpyN(si.getColSolution(), nCols, colsol);
  // CoinBigIndex *COIN_RESTRICT rowStartPos = new CoinBigIndex [nRows];
  const int *COIN_RESTRICT rowLength = rowCopy->getVectorLengths();
  for (int i = 0; i < nRows; i++) {
    assert(rowStart[i] + rowLength[i] == rowStart[i + 1]);
    CoinBigIndex kk;
#ifndef NDEBUG
    for (kk = rowStart[i]; kk < rowStart[i + 1]; kk++) {
      double value = rowElements[kk];
      if (value > 0.0)
        break;
    }
    assert(rowStartPos[i] == kk);
#endif
    double value;
    value = 0.0;
    for (kk = rowStart[i]; kk < rowStartPos[i]; kk++) {
      int iColumn = column[kk];
      double gap = std::min(1.0e100, colUpper[iColumn] - colLower[iColumn]);
      value = std::min(value, gap * rowElements[kk]);
    }
    largestNegativeInRow[i] = value;
    value = 0.0;
    for (; kk < rowStart[i + 1]; kk++) {
      int iColumn = column[kk];
      double gap = std::min(1.0e100, colUpper[iColumn] - colLower[iColumn]);
      value = std::max(value, gap * rowElements[kk]);
    }
    largestPositiveInRow[i] = value;
  }
  for (int i = 0; i < nCols; ++i) {
    double djValue = djs[i];
    double gap = colUpper[i] - colLower[i];
    if (gap > 1.0e-8) {
      if (colsol[i] < colLower[i] + primalTolerance_) {
        colsol[i] = colLower[i];
        djs[i] = std::max(0.0, djValue);
      } else if (colsol[i] > colUpper[i] - primalTolerance_) {
        colsol[i] = colUpper[i];
        djs[i] = std::min(0.0, djValue);
      } else {
        djs[i] = 0.0;
      }
    }
    columnGap[i] = gap - primalTolerance_;
  }

  int ipass = 0, nfixed = -1;

  double cutoff;
  bool cutoff_available = si.getDblParam(OsiDualObjectiveLimit, cutoff);
  if (!cutoff_available || usingObjective_ < 0) { // cut off isn't set or isn't valid
    cutoff = si.getInfinity();
  }
  if (fabs(cutoff) > 1.0e30)
    assert(cutoff > 1.0e30);
  double current = si.getObjValue();
  /* for both way coding */
  int nstackC0 = -1;
  int nstackR, nstackC;
  // int nFix=0;
  for (int i = 0; i < nCols; i++) {
    if (colUpper[i] - colLower[i] < 1.0e-8) {
      markC[i] = 3;
      // nFix++;
    } else {
      markC[i] = 0;
      if (colUpper[i] > 1.0e10)
        markC[i] |= 8;
      if (colLower[i] < -1.0e10)
        markC[i] |= 4;
    }
  }
  double tolerance = 1.0e1 * primalTolerance_;
  const double probingDeadline = (maxSeconds_ > 0.0) ? (CoinWallclockTime() + maxSeconds_) : 1.0e100;
  while (ipass < maxPass && nfixed) {
    int iLook;
    if (probingDeadline && CoinWallclockTime() > probingDeadline)
      break;
    ipass++;
    // printf("pass %d\n",ipass);
    nfixed = 0;
    int justFix = (!info->inTree && !info->pass) ? -1 : 0;
    int maxProbe = info->inTree ? maxProbe_ : maxProbeRoot_;
    if (justFix < 0)
      maxProbe = numberThisTime_;
    if (maxProbe == 123) {
      // Try and be a bit intelligent
      maxProbe = 0;
      if (!info->inTree) {
        if (!info->pass || numberThisTime_ < -100) {
          maxProbe = numberThisTime_;
        } else {
          int cutDown = 4;
          int offset = info->pass % cutDown;
          int i;
          int k = 0;
          int kk = offset;
          for (i = 0; i < numberThisTime_; i++) {
            if (!kk) {
#define XXXXXX
#ifdef XXXXXX
              lookedAt_[maxProbe] = lookedAt_[i];
#endif
              maxProbe++;
              kk = cutDown - 1;
            } else {
              stackC[k++] = lookedAt_[i];
              kk--;
            }
          }
#ifdef XXXXXX
          memcpy(lookedAt_ + maxProbe, stackC, k * sizeof(int));
#endif
        }
      } else {
        // in tree
        if (numberThisTime_ < 200) {
          maxProbe = numberThisTime_;
        } else {
          int cutDown = std::max(numberThisTime_ / 100, 4);
          int offset = info->pass % cutDown;
          int i;
          int k = 0;
          int kk = offset;
          for (i = 0; i < numberThisTime_; i++) {
            if (!kk) {
#ifdef XXXXXX
              lookedAt_[maxProbe] = lookedAt_[i];
#endif
              maxProbe++;
              kk = cutDown - 1;
            } else {
              stackC[k++] = lookedAt_[i];
              kk--;
            }
          }
#ifdef XXXXXX
          memcpy(lookedAt_ + maxProbe, stackC, k * sizeof(int));
#endif
        }
      }
    }
    // printf("maxStack %d maxPass %d numberThisTime %d info pass %d\n",
    //    maxStack,maxPass,numberThisTime_,info->pass);
    for (iLook = 0; iLook < numberThisTime_; iLook++) {
      if (probingDeadline && (iLook & 127) == 0 && CoinWallclockTime() > probingDeadline) {
        ipass = maxPass;
        break;
      }
      double solval;
      double down;
      double up;
      int j = lookedAt_[iLook];
#ifdef UPDATE_MINR_MAXR
      // update
      bool feasible2 = true;
      for (int i = 0; i < numberRealChanges; i++) {
        int jcol = whichRealChange[i];
        double movement = realChange[i];
        CoinBigIndex k;
        for (k = columnStart[jcol];
          k < columnStart[jcol] + columnLength[jcol]; k++) {
          int irow = row[k];
          double value = columnElements[k];
          if (movement > 0.0) {
            /* up */
            if (value > 0.0) {
              /* up does not change - down does */
              if (minR[irow] > -1.0e10)
                minR[irow] += value * movement;
              if (minR[irow] > rowUpper[irow] + 1.0e-5) {
                feasible2 = false;
                break;
              }
            } else {
              /* down does not change - up does */
              if (maxR[irow] < 1.0e10)
                maxR[irow] += value * movement;
              if (maxR[irow] < rowLower[irow] - 1.0e-5) {
                feasible2 = false;
                break;
              }
            }
          } else {
            /* down */
            if (value < 0.0) {
              /* up does not change - down does */
              if (minR[irow] > -1.0e10)
                minR[irow] += value * movement;
              if (minR[irow] > rowUpper[irow] + 1.0e-5) {
                feasible2 = false;
                break;
              }
            } else {
              /* down does not change - up does */
              if (maxR[irow] < 1.0e10)
                maxR[irow] += value * movement;
              if (maxR[irow] < rowLower[irow] - 1.0e-5) {
                feasible2 = false;
                break;
              }
            }
          }
        }
      }
      if (!feasible2) {
        printf("NOT feasible\n");
      }
      numberRealChanges = 0;
#endif
#if FIXED_BOTH_WAYS
      // clean newLo
      if (doFixedBothWays) {
        for (int i = 0; i < nFixedBoth; i++) {
          int icol = fixedBoth[i];
          newLo[icol] = COIN_DBL_MAX;
        }
        nFixedBoth = 0;
      }
#endif
      solval = colsol[j];
      down = floor(solval + tolerance);
      up = ceil(solval - tolerance);
      if (columnGap[j] < 0.0)
        markC[j] = 3;
      if ((markC[j] & 3) != 0 || !intVar[j])
        continue;
      if (solval >= colUpper[j] - tolerance || solval <= colLower[j] + tolerance || up == down) {
        if (solval <= colLower[j] + 2.0 * tolerance) {
          solval = colLower[j] + 1.0e-1;
          down = colLower[j];
          up = down + 1.0;
        } else if (solval >= colUpper[j] - 2.0 * tolerance) {
          solval = colUpper[j] - 1.0e-1;
          up = colUpper[j];
          down = up - 1;
        } else {
          // odd
          up = down + 1.0;
          solval = down + 1.0e-1;
        }
      }
      assert(up <= colUpper[j]);
      assert(down >= colLower[j]);
      assert(up > down);
      int istackC, iway, istackR;
      int way[] = { 1, 2, 1 };
      int feas[] = { 1, 2, 4 };
      int feasible = 0;
      int notFeasible;
      for (iway = 0; iway < 3; iway++) {
        // int fixThis=0;
        double objVal = current;
        stackC[0] = j;
        markC[j] = way[iway];
        double solMovement;
        double movement;
        if (way[iway] == 1) {
          movement = down - colUpper[j];
          solMovement = down - colsol[j];
          assert(movement < -0.99999);
          if (fabs(down - colLower[j]) < 1.0e-7) {
            down = colLower[j];
          }
        } else {
          movement = up - colLower[j];
          solMovement = up - colsol[j];
          assert(movement > 0.99999);
          if (fabs(up - colUpper[j]) < 1.0e-7) {
            up = colUpper[j];
          }
        }
#ifdef PRINT_DEBUG
        if (fabs(movement) > 1.01) {
          printf("big %d %g %g %g\n", j, colLower[j], solval, colUpper[j]);
        }
#endif
        if (solMovement * djs[j] > 0.0)
          objVal += solMovement * djs[j];
        nstackC = 1;
        nstackR = 0;
        saveL[0] = colLower[j];
        saveU[0] = colUpper[j];
        assert(saveU[0] > saveL[0]);
        notFeasible = 0;
        if (movement < 0.0) {
          colUpper[j] += movement;
          colUpper[j] = floor(colUpper[j] + 0.5);
          columnGap[j] = colUpper[j] - colLower[j] - primalTolerance_;
#ifdef PRINT_DEBUG
          printf("** Trying %d down to 0\n", j);
#endif
        } else {
          colLower[j] += movement;
          colLower[j] = floor(colLower[j] + 0.5);
          columnGap[j] = colUpper[j] - colLower[j] - primalTolerance_;
#ifdef PRINT_DEBUG
          printf("** Trying %d up to 1\n", j);
#endif
        }
        if (fabs(colUpper[j] - colLower[j]) < 1.0e-6)
          markC[j] = 3; // say fixed
        markC[j] &= ~12;
        if (colUpper[j] > 1.0e10)
          markC[j] |= 8;
        if (colLower[j] < -1.0e10)
          markC[j] |= 4;
        istackC = 0;
        /* update immediately */
        CoinBigIndex k;
        for (k = columnStart[j]; k < columnStart[j] + columnLength[j]; k++) {
          int irow = row[k];
          double value = columnElements[k];
          assert(markR[irow] != -2);
          if (markR[irow] == -1) {
            stackR[nstackR] = irow;
            markR[irow] = nstackR;
            saveMin[nstackR] = minR[irow];
            saveMax[nstackR] = maxR[irow];
            nstackR++;
          }
          /* could check immediately if violation */
          if (movement > 0.0) {
            /* up */
            if (value > 0.0) {
              /* up does not change - down does */
              if (minR[irow] > -1.0e10)
                minR[irow] += value * movement;
              if (minR[irow] > rowUpper[irow] + 1.0e-5) {
                notFeasible = 1;
                istackC = 1;
                break;
              }
            } else {
              /* down does not change - up does */
              if (maxR[irow] < 1.0e10)
                maxR[irow] += value * movement;
              if (maxR[irow] < rowLower[irow] - 1.0e-5) {
                notFeasible = 1;
                istackC = 1;
                break;
              }
            }
          } else {
            /* down */
            if (value < 0.0) {
              /* up does not change - down does */
              if (minR[irow] > -1.0e10)
                minR[irow] += value * movement;
              if (minR[irow] > rowUpper[irow] + 1.0e-5) {
                notFeasible = 1;
                istackC = 1;
                break;
              }
            } else {
              /* down does not change - up does */
              if (maxR[irow] < 1.0e10)
                maxR[irow] += value * movement;
              if (maxR[irow] < rowLower[irow] - 1.0e-5) {
                notFeasible = 1;
                istackC = 1;
                break;
              }
            }
          }
        }
        while (istackC < nstackC && nstackC < maxStack) { // could be istackC<maxStack?
          int jway;
          int jcol = stackC[istackC];
          jway = markC[jcol];
          for (k = columnStart[jcol]; k < columnStart[jcol] + columnLength[jcol]; k++) {
            // break if found not feasible
            if (notFeasible)
              break;
            int irow = row[k];
            /* see if anything forced */
            CoinBigIndex rStart = rowStart[irow];
            CoinBigIndex rEnd = rowStartPos[irow];
            double rowUp = rowUpper[irow];
            double rowUp2 = 0.0;
            bool doRowUpN;
            bool doRowUpP;
            if (rowUp < 1.0e10) {
              doRowUpN = true;
              doRowUpP = true;
              rowUp2 = rowUp - minR[irow];
              if (rowUp2 < -tolerance) {
                notFeasible = true;
                break;
              } else {
                if (rowUp2 + largestNegativeInRow[irow] > 0)
                  doRowUpN = false;
                if (rowUp2 - largestPositiveInRow[irow] > 0)
                  doRowUpP = false;
              }
            } else {
              doRowUpN = false;
              doRowUpP = false;
              rowUp2 = COIN_DBL_MAX;
            }
            double rowLo = rowLower[irow];
            double rowLo2 = 0.0;
            bool doRowLoN;
            bool doRowLoP;
            if (rowLo > -1.0e10) {
              doRowLoN = true;
              doRowLoP = true;
              rowLo2 = rowLo - maxR[irow];
              if (rowLo2 > tolerance) {
                notFeasible = true;
                break;
              } else {
                if (rowLo2 - largestNegativeInRow[irow] < 0)
                  doRowLoN = false;
                if (rowLo2 + largestPositiveInRow[irow] < 0)
                  doRowLoP = false;
              }
            } else {
              doRowLoN = false;
              doRowLoP = false;
              rowLo2 = -COIN_DBL_MAX;
            }
            if (doRowUpN && doRowLoN) {
              // doRowUpN=doRowLoN=false;
              //  Start neg values loop
              for (CoinBigIndex kk = rStart; kk < rEnd; kk++) {
                int kcol = column[kk];
                int markIt = markC[kcol];
                if ((markIt & 3) != 3) {
                  double value2 = rowElements[kk];
                  if (colUpper[kcol] <= 1e10)
                    assert((markIt & 8) == 0);
                  else
                    assert((markIt & 8) != 0);
                  if (colLower[kcol] >= -1e10)
                    assert((markIt & 4) == 0);
                  else
                    assert((markIt & 4) != 0);
                  assert(value2 < 0.0);
                  double gap = columnGap[kcol] * value2;
                  bool doUp = (rowUp2 + gap < -1.0e-7);
                  bool doDown = (rowLo2 - gap > 1.0e-7);
                  if (doUp || doDown) {
                    double moveUp = 0.0;
                    double moveDown = 0.0;
                    double newUpper = -1.0;
                    double newLower = 1.0;
                    if (doUp && (markIt & (2 + 8)) == 0) {
                      double dbound = colUpper[kcol] + rowUp2 / value2;
                      if (intVar[kcol]) {
                        markIt |= 2;
                        newLower = ceil(dbound - relaxedTolerance);
                      } else {
                        newLower = dbound;
                        if (newLower + primalTolerance_ > colUpper[kcol] && newLower - primalTolerance_ <= colUpper[kcol]) {
                          newLower = colUpper[kcol];
                          markIt |= 2;
                          // markIt=3;
                        } else {
                          // avoid problems - fix later ?
                          markIt = 3;
                        }
                      }
                      moveUp = newLower - colLower[kcol];
                    }
                    if (doDown && (markIt & (1 + 4)) == 0) {
                      double dbound = colLower[kcol] + rowLo2 / value2;
                      if (intVar[kcol]) {
                        markIt |= 1;
                        newUpper = floor(dbound + relaxedTolerance);
                      } else {
                        newUpper = dbound;
                        if (newUpper - primalTolerance_ < colLower[kcol] && newUpper + primalTolerance_ >= colLower[kcol]) {
                          newUpper = colLower[kcol];
                          markIt |= 1;
                          // markIt=3;
                        } else {
                          // avoid problems - fix later ?
                          markIt = 3;
                        }
                      }
                      moveDown = newUpper - colUpper[kcol];
                    }
                    if (!moveUp && !moveDown)
                      continue;
                    bool onList = ((markC[kcol] & 3) != 0);
                    if (nstackC < 2 * maxStack) {
                      markC[kcol] = markIt;
                    }
                    if (moveUp && nstackC < 2 * maxStack) {
                      // fixThis++;
                      if (!onList) {
                        stackC[nstackC] = kcol;
                        saveL[nstackC] = colLower[kcol];
                        saveU[nstackC] = colUpper[kcol];
                        assert(saveU[nstackC] > saveL[nstackC]);
                        assert(nstackC < nCols);
                        nstackC++;
                        onList = true;
                      }
                      if (newLower > colsol[kcol]) {
                        if (djs[kcol] < 0.0) {
                          /* should be infeasible if integer */
                          if (intVar[kcol])
                            assert(newLower > colUpper[kcol] + primalTolerance_);
                        } else {
                          objVal += moveUp * djs[kcol];
                        }
                      }
                      if (intVar[kcol])
                        newLower = std::max(colLower[kcol], ceil(newLower - 1.0e-4));
                      colLower[kcol] = newLower;
                      columnGap[kcol] = colUpper[kcol] - newLower - primalTolerance_;
                      if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                        markC[kcol] = 3; // say fixed
                      }
                      markC[kcol] &= ~12;
                      if (colUpper[kcol] > 1.0e10)
                        markC[kcol] |= 8;
                      if (colLower[kcol] < -1.0e10)
                        markC[kcol] |= 4;
                      /* update immediately */
                      for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                        int krow = row[jj];
                        double value = columnElements[jj];
                        assert(markR[krow] != -2);
                        if (markR[krow] == -1) {
                          stackR[nstackR] = krow;
                          markR[krow] = nstackR;
                          saveMin[nstackR] = minR[krow];
                          saveMax[nstackR] = maxR[krow];
                          nstackR++;
                        }
                        /* could check immediately if violation */
                        /* up */
                        if (value > 0.0) {
                          /* up does not change - down does */
                          if (minR[krow] > -1.0e10)
                            minR[krow] += value * moveUp;
                          if (krow == irow)
                            rowUp2 = rowUp - minR[irow];
                          if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            break;
                          }
                        } else {
                          /* down does not change - up does */
                          if (maxR[krow] < 1.0e10)
                            maxR[krow] += value * moveUp;
                          if (krow == irow)
                            rowLo2 = rowLo - maxR[irow];
                          if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            break;
                          }
                        }
                      }
                    }
                    if (moveDown && nstackC < 2 * maxStack) {
                      // fixThis++;
                      if (!onList) {
                        stackC[nstackC] = kcol;
                        saveL[nstackC] = colLower[kcol];
                        saveU[nstackC] = colUpper[kcol];
                        assert(saveU[nstackC] > saveL[nstackC]);
                        assert(nstackC < nCols);
                        nstackC++;
                        onList = true;
                      }
                      if (newUpper < colsol[kcol]) {
                        if (djs[kcol] > 0.0) {
                          /* should be infeasible if integer */
                          if (intVar[kcol])
                            assert(colLower[kcol] > newUpper + primalTolerance_);
                        } else {
                          objVal += moveDown * djs[kcol];
                        }
                      }
                      if (intVar[kcol])
                        newUpper = std::min(colUpper[kcol], floor(newUpper + 1.0e-4));
                      colUpper[kcol] = newUpper;
                      columnGap[kcol] = newUpper - colLower[kcol] - primalTolerance_;
                      if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                        markC[kcol] = 3; // say fixed
                      }
                      markC[kcol] &= ~12;
                      if (colUpper[kcol] > 1.0e10)
                        markC[kcol] |= 8;
                      if (colLower[kcol] < -1.0e10)
                        markC[kcol] |= 4;
                      /* update immediately */
                      for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                        int krow = row[jj];
                        double value = columnElements[jj];
                        assert(markR[krow] != -2);
                        if (markR[krow] == -1) {
                          stackR[nstackR] = krow;
                          markR[krow] = nstackR;
                          saveMin[nstackR] = minR[krow];
                          saveMax[nstackR] = maxR[krow];
                          nstackR++;
                        }
                        /* could check immediately if violation */
                        /* down */
                        if (value < 0.0) {
                          /* up does not change - down does */
                          if (minR[krow] > -1.0e10)
                            minR[krow] += value * moveDown;
                          if (krow == irow)
                            rowUp2 = rowUp - minR[irow];
                          if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        } else {
                          /* down does not change - up does */
                          if (maxR[krow] < 1.0e10)
                            maxR[krow] += value * moveDown;
                          if (krow == irow)
                            rowLo2 = rowLo - maxR[irow];
                          if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        }
                      }
                    }
                    if (colLower[kcol] > colUpper[kcol] + primalTolerance_) {
                      notFeasible = 1;
                      ;
                      k = columnStart[jcol] + columnLength[jcol];
                      istackC = nstackC + 1;
                      break;
                    }
                  }
                }
              } // end big loop rStart->rPos
            } else if (doRowUpN) {
              // Start neg values loop
              for (CoinBigIndex kk = rStart; kk < rEnd; kk++) {
                int kcol = column[kk];
                int markIt = markC[kcol];
                if ((markIt & 3) != 3) {
                  double value2 = rowElements[kk];
                  double gap = columnGap[kcol] * value2;
                  if (!(rowUp2 + gap < -1.0e-7))
                    continue;
                  double moveUp = 0.0;
                  double newLower = 1.0;
                  if ((markIt & (2 + 8)) == 0) {
                    double dbound = colUpper[kcol] + rowUp2 / value2;
                    if (intVar[kcol]) {
                      markIt |= 2;
                      newLower = ceil(dbound - relaxedTolerance);
                    } else {
                      newLower = dbound;
                      if (newLower + primalTolerance_ > colUpper[kcol] && newLower - primalTolerance_ <= colUpper[kcol]) {
                        newLower = colUpper[kcol];
                        markIt |= 2;
                        // markIt=3;
                      } else {
                        // avoid problems - fix later ?
                        markIt = 3;
                      }
                    }
                    moveUp = newLower - colLower[kcol];
                    if (!moveUp)
                      continue;
                    bool onList = ((markC[kcol] & 3) != 0);
                    if (nstackC < 2 * maxStack) {
                      markC[kcol] = markIt;
                    }
                    if (moveUp && nstackC < 2 * maxStack) {
                      // fixThis++;
                      if (!onList) {
                        stackC[nstackC] = kcol;
                        saveL[nstackC] = colLower[kcol];
                        saveU[nstackC] = colUpper[kcol];
                        assert(saveU[nstackC] > saveL[nstackC]);
                        assert(nstackC < nCols);
                        nstackC++;
                        onList = true;
                      }
                      if (newLower > colsol[kcol]) {
                        if (djs[kcol] < 0.0) {
                          /* should be infeasible if integer */
                          if (intVar[kcol])
                            assert(newLower > colUpper[kcol] + primalTolerance_);
                        } else {
                          objVal += moveUp * djs[kcol];
                        }
                      }
                      if (intVar[kcol])
                        newLower = std::max(colLower[kcol], ceil(newLower - 1.0e-4));
                      colLower[kcol] = newLower;
                      columnGap[kcol] = colUpper[kcol] - newLower - primalTolerance_;
                      if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                        markC[kcol] = 3; // say fixed
                      }
                      markC[kcol] &= ~12;
                      if (colUpper[kcol] > 1.0e10)
                        markC[kcol] |= 8;
                      if (colLower[kcol] < -1.0e10)
                        markC[kcol] |= 4;
                      /* update immediately */
                      for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                        int krow = row[jj];
                        double value = columnElements[jj];
                        assert(markR[krow] != -2);
                        if (markR[krow] == -1) {
                          stackR[nstackR] = krow;
                          markR[krow] = nstackR;
                          saveMin[nstackR] = minR[krow];
                          saveMax[nstackR] = maxR[krow];
                          nstackR++;
                        }
                        /* could check immediately if violation */
                        /* up */
                        if (value > 0.0) {
                          /* up does not change - down does */
                          if (minR[krow] > -1.0e10)
                            minR[krow] += value * moveUp;
                          if (krow == irow)
                            rowUp2 = rowUp - minR[irow];
                          if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        } else {
                          /* down does not change - up does */
                          if (maxR[krow] < 1.0e10)
                            maxR[krow] += value * moveUp;
                          if (krow == irow)
                            rowLo2 = rowLo - maxR[irow];
                          if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        }
                      }
                    }
                    if (colLower[kcol] > colUpper[kcol] + primalTolerance_) {
                      notFeasible = 1;
                      ;
                      k = columnStart[jcol] + columnLength[jcol];
                      istackC = nstackC + 1;
                      break;
                    }
                  }
                }
              } // end big loop rStart->rPos
            } else if (doRowLoN) {
              // Start neg values loop
              for (CoinBigIndex kk = rStart; kk < rEnd; kk++) {
                int kcol = column[kk];
                if ((markC[kcol] & 3) != 3) {
                  double moveDown = 0.0;
                  double newUpper = -1.0;
                  double value2 = rowElements[kk];
                  int markIt = markC[kcol];
                  assert(value2 < 0.0);
                  double gap = columnGap[kcol] * value2;
                  bool doDown = (rowLo2 - gap > 1.0e-7);
                  if (doDown && (markIt & (1 + 4)) == 0) {
                    double dbound = colLower[kcol] + rowLo2 / value2;
                    if (intVar[kcol]) {
                      markIt |= 1;
                      newUpper = floor(dbound + relaxedTolerance);
                    } else {
                      newUpper = dbound;
                      if (newUpper - primalTolerance_ < colLower[kcol] && newUpper + primalTolerance_ >= colLower[kcol]) {
                        newUpper = colLower[kcol];
                        markIt |= 1;
                        // markIt=3;
                      } else {
                        // avoid problems - fix later ?
                        markIt = 3;
                      }
                    }
                    moveDown = newUpper - colUpper[kcol];
                    if (!moveDown)
                      continue;
                    bool onList = ((markC[kcol] & 3) != 0);
                    if (nstackC < 2 * maxStack) {
                      markC[kcol] = markIt;
                    }
                    if (moveDown && nstackC < 2 * maxStack) {
                      // fixThis++;
                      if (!onList) {
                        stackC[nstackC] = kcol;
                        saveL[nstackC] = colLower[kcol];
                        saveU[nstackC] = colUpper[kcol];
                        assert(saveU[nstackC] > saveL[nstackC]);
                        assert(nstackC < nCols);
                        nstackC++;
                        onList = true;
                      }
                      if (newUpper < colsol[kcol]) {
                        if (djs[kcol] > 0.0) {
                          /* should be infeasible if integer */
                          if (intVar[kcol])
                            assert(colLower[kcol] > newUpper + primalTolerance_);
                        } else {
                          objVal += moveDown * djs[kcol];
                        }
                      }
                      if (intVar[kcol])
                        newUpper = std::min(colUpper[kcol], floor(newUpper + 1.0e-4));
                      colUpper[kcol] = newUpper;
                      columnGap[kcol] = newUpper - colLower[kcol] - primalTolerance_;
                      if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                        markC[kcol] = 3; // say fixed
                      }
                      markC[kcol] &= ~12;
                      if (colUpper[kcol] > 1.0e10)
                        markC[kcol] |= 8;
                      if (colLower[kcol] < -1.0e10)
                        markC[kcol] |= 4;
                      /* update immediately */
                      for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                        int krow = row[jj];
                        double value = columnElements[jj];
                        assert(markR[krow] != -2);
                        if (markR[krow] == -1) {
                          stackR[nstackR] = krow;
                          markR[krow] = nstackR;
                          saveMin[nstackR] = minR[krow];
                          saveMax[nstackR] = maxR[krow];
                          nstackR++;
                        }
                        /* could check immediately if violation */
                        /* down */
                        if (value < 0.0) {
                          /* up does not change - down does */
                          if (minR[krow] > -1.0e10)
                            minR[krow] += value * moveDown;
                          if (krow == irow)
                            rowUp2 = rowUp - minR[irow];
                          if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        } else {
                          /* down does not change - up does */
                          if (maxR[krow] < 1.0e10)
                            maxR[krow] += value * moveDown;
                          if (krow == irow)
                            rowLo2 = rowLo - maxR[irow];
                          if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        }
                      }
                    }
                    if (colLower[kcol] > colUpper[kcol] + primalTolerance_) {
                      notFeasible = 1;
                      ;
                      k = columnStart[jcol] + columnLength[jcol];
                      istackC = nstackC + 1;
                      break;
                    }
                  }
                }
              } // end big loop rStart->rPos
            }
            rStart = rEnd;
            rEnd = rowStart[irow + 1];
            if (doRowUpP && doRowLoP) {
              // doRowUpP=doRowLoP=false;
              //  Start pos values loop
              for (CoinBigIndex kk = rStart; kk < rEnd; kk++) {
                int kcol = column[kk];
                int markIt = markC[kcol];
                if ((markIt & 3) != 3) {
                  double value2 = rowElements[kk];
                  assert(value2 > 0.0);
                  /* positive element */
                  double gap = columnGap[kcol] * value2;
                  bool doDown = (rowLo2 + gap > 1.0e-7);
                  bool doUp = (rowUp2 - gap < -1.0e-7);
                  if (doDown || doUp) {
                    double moveUp = 0.0;
                    double moveDown = 0.0;
                    double newUpper = -1.0;
                    double newLower = 1.0;
                    if (doDown && (markIt & (2 + 8)) == 0) {
                      double dbound = colUpper[kcol] + rowLo2 / value2;
                      if (intVar[kcol]) {
                        markIt |= 2;
                        newLower = ceil(dbound - relaxedTolerance);
                      } else {
                        newLower = dbound;
                        if (newLower + primalTolerance_ > colUpper[kcol] && newLower - primalTolerance_ <= colUpper[kcol]) {
                          newLower = colUpper[kcol];
                          markIt |= 2;
                          // markIt=3;
                        } else {
                          // avoid problems - fix later ?
                          markIt = 3;
                        }
                      }
                      moveUp = newLower - colLower[kcol];
                    }
                    if (doUp && (markIt & (1 + 4)) == 0) {
                      double dbound = colLower[kcol] + rowUp2 / value2;
                      if (intVar[kcol]) {
                        markIt |= 1;
                        newUpper = floor(dbound + relaxedTolerance);
                      } else {
                        newUpper = dbound;
                        if (newUpper - primalTolerance_ < colLower[kcol] && newUpper + primalTolerance_ >= colLower[kcol]) {
                          newUpper = colLower[kcol];
                          markIt |= 1;
                          // markIt=3;
                        } else {
                          // avoid problems - fix later ?
                          markIt = 3;
                        }
                      }
                      moveDown = newUpper - colUpper[kcol];
                    }
                    if (!moveUp && !moveDown)
                      continue;
                    bool onList = ((markC[kcol] & 3) != 0);
                    if (nstackC < 2 * maxStack) {
                      markC[kcol] = markIt;
                    }
                    if (moveUp && nstackC < 2 * maxStack) {
                      // fixThis++;
                      if (!onList) {
                        stackC[nstackC] = kcol;
                        saveL[nstackC] = colLower[kcol];
                        saveU[nstackC] = colUpper[kcol];
                        assert(saveU[nstackC] > saveL[nstackC]);
                        assert(nstackC < nCols);
                        nstackC++;
                        onList = true;
                      }
                      if (newLower > colsol[kcol]) {
                        if (djs[kcol] < 0.0) {
                          /* should be infeasible if integer */
                          if (intVar[kcol])
                            assert(newLower > colUpper[kcol] + primalTolerance_);
                        } else {
                          objVal += moveUp * djs[kcol];
                        }
                      }
                      if (intVar[kcol])
                        newLower = std::max(colLower[kcol], ceil(newLower - 1.0e-4));
                      colLower[kcol] = newLower;
                      columnGap[kcol] = colUpper[kcol] - newLower - primalTolerance_;
                      if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                        markC[kcol] = 3; // say fixed
                      }
                      markC[kcol] &= ~12;
                      if (colUpper[kcol] > 1.0e10)
                        markC[kcol] |= 8;
                      if (colLower[kcol] < -1.0e10)
                        markC[kcol] |= 4;
                      /* update immediately */
                      for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                        int krow = row[jj];
                        double value = columnElements[jj];
                        assert(markR[krow] != -2);
                        if (markR[krow] == -1) {
                          stackR[nstackR] = krow;
                          markR[krow] = nstackR;
                          saveMin[nstackR] = minR[krow];
                          saveMax[nstackR] = maxR[krow];
                          nstackR++;
                        }
                        /* could check immediately if violation */
                        /* up */
                        if (value > 0.0) {
                          /* up does not change - down does */
                          if (minR[krow] > -1.0e10)
                            minR[krow] += value * moveUp;
                          if (krow == irow)
                            rowUp2 = rowUp - minR[irow];
                          if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        } else {
                          /* down does not change - up does */
                          if (maxR[krow] < 1.0e10)
                            maxR[krow] += value * moveUp;
                          if (krow == irow)
                            rowLo2 = rowLo - maxR[irow];
                          if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        }
                      }
                    }
                    if (moveDown && nstackC < 2 * maxStack) {
                      // fixThis++;
                      if (!onList) {
                        stackC[nstackC] = kcol;
                        saveL[nstackC] = colLower[kcol];
                        saveU[nstackC] = colUpper[kcol];
                        assert(saveU[nstackC] > saveL[nstackC]);
                        assert(nstackC < nCols);
                        nstackC++;
                        onList = true;
                      }
                      if (newUpper < colsol[kcol]) {
                        if (djs[kcol] > 0.0) {
                          /* should be infeasible if integer */
                          if (intVar[kcol])
                            assert(colLower[kcol] > newUpper + primalTolerance_);
                        } else {
                          objVal += moveDown * djs[kcol];
                        }
                      }
                      if (intVar[kcol])
                        newUpper = std::min(colUpper[kcol], floor(newUpper + 1.0e-4));
                      colUpper[kcol] = newUpper;
                      columnGap[kcol] = newUpper - colLower[kcol] - primalTolerance_;
                      if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                        markC[kcol] = 3; // say fixed
                      }
                      markC[kcol] &= ~12;
                      if (colUpper[kcol] > 1.0e10)
                        markC[kcol] |= 8;
                      if (colLower[kcol] < -1.0e10)
                        markC[kcol] |= 4;
                      /* update immediately */
                      for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                        int krow = row[jj];
                        double value = columnElements[jj];
                        assert(markR[krow] != -2);
                        if (markR[krow] == -1) {
                          stackR[nstackR] = krow;
                          markR[krow] = nstackR;
                          saveMin[nstackR] = minR[krow];
                          saveMax[nstackR] = maxR[krow];
                          nstackR++;
                        }
                        /* could check immediately if violation */
                        /* down */
                        if (value < 0.0) {
                          /* up does not change - down does */
                          if (minR[krow] > -1.0e10)
                            minR[krow] += value * moveDown;
                          if (krow == irow)
                            rowUp2 = rowUp - minR[irow];
                          if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        } else {
                          /* down does not change - up does */
                          if (maxR[krow] < 1.0e10)
                            maxR[krow] += value * moveDown;
                          if (krow == irow)
                            rowLo2 = rowLo - maxR[irow];
                          if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        }
                      }
                    }
                    if (colLower[kcol] > colUpper[kcol] + primalTolerance_) {
                      notFeasible = 1;
                      ;
                      k = columnStart[jcol] + columnLength[jcol];
                      istackC = nstackC + 1;
                      break;
                    }
                  }
                }
              } // end big loop rPos->rEnd
            } else if (doRowUpP) {
              // Start pos values loop
              for (CoinBigIndex kk = rStart; kk < rEnd; kk++) {
                int kcol = column[kk];
                int markIt = markC[kcol];
                if ((markIt & 3) != 3) {
                  double value2 = rowElements[kk];
                  assert(value2 > 0.0);
                  /* positive element */
                  double gap = columnGap[kcol] * value2;
                  bool doUp = (rowUp2 - gap < -1.0e-7);
                  if (doUp && (markIt & (1 + 4)) == 0) {
                    double newUpper = -1.0;
                    double dbound = colLower[kcol] + rowUp2 / value2;
                    if (intVar[kcol]) {
                      markIt |= 1;
                      newUpper = floor(dbound + relaxedTolerance);
                    } else {
                      newUpper = dbound;
                      if (newUpper - primalTolerance_ < colLower[kcol] && newUpper + primalTolerance_ >= colLower[kcol]) {
                        newUpper = colLower[kcol];
                        markIt |= 1;
                        // markIt=3;
                      } else {
                        // avoid problems - fix later ?
                        markIt = 3;
                      }
                    }
                    double moveDown = newUpper - colUpper[kcol];
                    if (!moveDown)
                      continue;
                    bool onList = ((markC[kcol] & 3) != 0);
                    if (nstackC < 2 * maxStack) {
                      markC[kcol] = markIt;
                    }
                    if (nstackC < 2 * maxStack) {
                      // fixThis++;
                      if (!onList) {
                        stackC[nstackC] = kcol;
                        saveL[nstackC] = colLower[kcol];
                        saveU[nstackC] = colUpper[kcol];
                        assert(saveU[nstackC] > saveL[nstackC]);
                        assert(nstackC < nCols);
                        nstackC++;
                        onList = true;
                      }
                      if (newUpper < colsol[kcol]) {
                        if (djs[kcol] > 0.0) {
                          /* should be infeasible if integer */
                          if (intVar[kcol])
                            assert(colLower[kcol] > newUpper + primalTolerance_);
                        } else {
                          objVal += moveDown * djs[kcol];
                        }
                      }
                      if (intVar[kcol])
                        newUpper = std::min(colUpper[kcol], floor(newUpper + 1.0e-4));
                      colUpper[kcol] = newUpper;
                      columnGap[kcol] = newUpper - colLower[kcol] - primalTolerance_;
                      if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                        markC[kcol] = 3; // say fixed
                      }
                      markC[kcol] &= ~12;
                      if (colUpper[kcol] > 1.0e10)
                        markC[kcol] |= 8;
                      if (colLower[kcol] < -1.0e10)
                        markC[kcol] |= 4;
                      /* update immediately */
                      for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                        int krow = row[jj];
                        double value = columnElements[jj];
                        assert(markR[krow] != -2);
                        if (markR[krow] == -1) {
                          stackR[nstackR] = krow;
                          markR[krow] = nstackR;
                          saveMin[nstackR] = minR[krow];
                          saveMax[nstackR] = maxR[krow];
                          nstackR++;
                        }
                        /* could check immediately if violation */
                        /* down */
                        if (value < 0.0) {
                          /* up does not change - down does */
                          if (minR[krow] > -1.0e10)
                            minR[krow] += value * moveDown;
                          if (krow == irow)
                            rowUp2 = rowUp - minR[irow];
                          if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        } else {
                          /* down does not change - up does */
                          if (maxR[krow] < 1.0e10)
                            maxR[krow] += value * moveDown;
                          if (krow == irow)
                            rowLo2 = rowLo - maxR[irow];
                          if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        }
                      }
                    }
                    if (colLower[kcol] > colUpper[kcol] + primalTolerance_) {
                      notFeasible = 1;
                      ;
                      k = columnStart[jcol] + columnLength[jcol];
                      istackC = nstackC + 1;
                      break;
                    }
                  }
                }
              } // end big loop rPos->rEnd
            } else if (doRowLoP) {
              // Start pos values loop
              for (CoinBigIndex kk = rStart; kk < rEnd; kk++) {
                int kcol = column[kk];
                if ((markC[kcol] & 3) != 3) {
                  double value2 = rowElements[kk];
                  int markIt = markC[kcol];
                  assert(value2 > 0.0);
                  /* positive element */
                  double gap = columnGap[kcol] * value2;
                  bool doDown = (rowLo2 + gap > 1.0e-7);
                  if (doDown && (markIt & (2 + 8)) == 0) {
                    double newLower = 1.0;
                    double dbound = colUpper[kcol] + rowLo2 / value2;
                    if (intVar[kcol]) {
                      markIt |= 2;
                      newLower = ceil(dbound - relaxedTolerance);
                    } else {
                      newLower = dbound;
                      if (newLower + primalTolerance_ > colUpper[kcol] && newLower - primalTolerance_ <= colUpper[kcol]) {
                        newLower = colUpper[kcol];
                        markIt |= 2;
                        // markIt=3;
                      } else {
                        // avoid problems - fix later ?
                        markIt = 3;
                      }
                    }
                    double moveUp = newLower - colLower[kcol];
                    if (!moveUp)
                      continue;
                    bool onList = ((markC[kcol] & 3) != 0);
                    if (nstackC < 2 * maxStack) {
                      markC[kcol] = markIt;
                    }
                    if (nstackC < 2 * maxStack) {
                      // fixThis++;
                      if (!onList) {
                        stackC[nstackC] = kcol;
                        saveL[nstackC] = colLower[kcol];
                        saveU[nstackC] = colUpper[kcol];
                        assert(saveU[nstackC] > saveL[nstackC]);
                        assert(nstackC < nCols);
                        nstackC++;
                        onList = true;
                      }
                      if (newLower > colsol[kcol]) {
                        if (djs[kcol] < 0.0) {
                          /* should be infeasible if integer */
                          if (intVar[kcol])
                            assert(newLower > colUpper[kcol] + primalTolerance_);
                        } else {
                          objVal += moveUp * djs[kcol];
                        }
                      }
                      if (intVar[kcol])
                        newLower = std::max(colLower[kcol], ceil(newLower - 1.0e-4));
                      colLower[kcol] = newLower;
                      columnGap[kcol] = colUpper[kcol] - newLower - primalTolerance_;
                      if (fabs(colUpper[kcol] - colLower[kcol]) < 1.0e-6) {
                        markC[kcol] = 3; // say fixed
                      }
                      markC[kcol] &= ~12;
                      if (colUpper[kcol] > 1.0e10)
                        markC[kcol] |= 8;
                      if (colLower[kcol] < -1.0e10)
                        markC[kcol] |= 4;
                      /* update immediately */
                      for (CoinBigIndex jj = columnStart[kcol]; jj < columnStart[kcol] + columnLength[kcol]; jj++) {
                        int krow = row[jj];
                        double value = columnElements[jj];
                        assert(markR[krow] != -2);
                        if (markR[krow] == -1) {
                          stackR[nstackR] = krow;
                          markR[krow] = nstackR;
                          saveMin[nstackR] = minR[krow];
                          saveMax[nstackR] = maxR[krow];
                          nstackR++;
                        }
                        /* could check immediately if violation */
                        /* up */
                        if (value > 0.0) {
                          /* up does not change - down does */
                          if (minR[krow] > -1.0e10)
                            minR[krow] += value * moveUp;
                          if (krow == irow)
                            rowUp2 = rowUp - minR[irow];
                          if (minR[krow] > rowUpper[krow] + 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        } else {
                          /* down does not change - up does */
                          if (maxR[krow] < 1.0e10)
                            maxR[krow] += value * moveUp;
                          if (krow == irow)
                            rowLo2 = rowLo - maxR[irow];
                          if (maxR[krow] < rowLower[krow] - 1.0e-5) {
                            colUpper[kcol] = -1.0e50; /* force infeasible */
                            columnGap[kcol] = -1.0e50;
                            break;
                          }
                        }
                      }
                    }
                    if (colLower[kcol] > colUpper[kcol] + primalTolerance_) {
                      notFeasible = 1;
                      ;
                      k = columnStart[jcol] + columnLength[jcol];
                      istackC = nstackC + 1;
                      break;
                    }
                  }
                }
              } // end big loop rPos->rEnd
            }
          }
          istackC++;
        }
#if FIXED_BOTH_WAYS
        // the nCols<1000 should be on some switch
        if (!notFeasible && iway < 2 && (info->options & 2048) != 0
          && nCols < 1000) {
          // try tighten
          double *COIN_RESTRICT rLower = const_cast< double * COIN_RESTRICT >(rowLower);
          double *COIN_RESTRICT rUpper = const_cast< double * COIN_RESTRICT >(rowUpper);
          // must be able to do better
          memcpy(saveFColLower, colLower, nCols * sizeof(double));
          memcpy(saveFColUpper, colUpper, nCols * sizeof(double));
          memcpy(saveFRowLower, rLower, nRows * sizeof(double));
          memcpy(saveFRowUpper, rUpper, nRows * sizeof(double));
          // move stackR
          memcpy(lookRow_[0], stackR, nstackR * sizeof(int));
          numberRows_ = nstackR;
          int returnCode = tighten(colLower, colUpper,
            column, rowElements, rowStart,
            rowStartPos, rowLength, columnCopy, rLower, rUpper,
            nRows, nCols, intVar, 10, primalTolerance_);
          if (returnCode) {
            // printf("INFEAS for j %d and way %d\n",j,iway);
            notFeasible = 1;
          }
#if 0
	  // maybe should do more than nstackC
	  for (int i=1;i<nstackC;i++) {
	    int icol = stackC[i];
	    assert (saveFColLower[icol]<=colLower[icol]);
	    assert (saveFColUpper[icol]>=colUpper[icol]);
	    saveFColLower[icol]=colLower[icol];
	    saveFColUpper[icol]>=colUpper[icol];
	  }
#endif
          memcpy(colLower, saveFColLower, nCols * sizeof(double));
          memcpy(colUpper, saveFColUpper, nCols * sizeof(double));
          memcpy(rLower, saveFRowLower, nRows * sizeof(double));
          memcpy(rUpper, saveFRowUpper, nRows * sizeof(double));
        }
#endif
        if (!notFeasible) {
          if (objVal <= cutoff) {
            feasible |= feas[iway];
          } else {
#ifdef PRINT_DEBUG
            printf("not feasible on dj\n");
#endif
            notFeasible = 1;
            if (iway == 1 && feasible == 0) {
              /* not feasible at all */
              ninfeas = 1;
              j = nCols - 1;
              break;
            }
          }
        } else if (iway == 1 && feasible == 0) {
          /* not feasible at all */
          ninfeas = 1;
          j = nCols - 1;
          iLook = numberThisTime_;
          ipass = maxPass;
          break;
        }
        if (iway == 2 || (iway == 1 && feasible == 2)) {
          /* keep */
          iway = 3;
          nfixed++;
          OsiColCut cc;
          int /*nTot=0,*/ nFix = 0, nInt = 0;
          for (istackC = 0; istackC < nstackC; istackC++) {
            int icol = stackC[istackC];
            if (intVar[icol]) {
              if (colUpper[icol] < currentColUpper[icol] - 1.0e-4) {
                element[nFix] = colUpper[icol];
                index[nFix++] = icol;
                nInt++;
              }
            }
          }
          if (nFix) {
            // nTot=nFix;
            cleanBounds(1, nFix, index, element);
            cc.setUbs(nFix, index, element);
            nFix = 0;
          }
          for (istackC = 0; istackC < nstackC; istackC++) {
            int icol = stackC[istackC];
            if (intVar[icol]) {
              if (colLower[icol] > currentColLower[icol] + 1.0e-4) {
                element[nFix] = colLower[icol];
                index[nFix++] = icol;
                nInt++;
              }
            }
          }
          if (nFix) {
            // nTot+=nFix;
            cleanBounds(-1, nFix, index, element);
            cc.setLbs(nFix, index, element);
          }
          // could tighten continuous as well
          if (nInt) {
#ifdef CGL_DEBUG
            checkBounds(debugger, cc);
#endif
#ifndef CC_AT_END
            cs.insert(cc);
#endif
          }
          for (istackC = 0; istackC < nstackC; istackC++) {
            int icol = stackC[istackC];
            if (colUpper[icol] - colLower[icol] > primalTolerance_) {
              markC[icol] &= ~3;
            } else {
              markC[icol] = 3;
            }
          }
          for (istackR = 0; istackR < nstackR; istackR++) {
            int irow = stackR[istackR];
            markR[irow] = -1;
          }
        } else {
          /* is it worth seeing if can increase coefficients
             or maybe better see if it is a cut */
          if (iway == 0) {
            nstackC0 = std::min(nstackC, maxStack);
            double boundChange;
            if (notFeasible) {
              nstackC0 = 0;
            } else {
              for (istackC = 0; istackC < nstackC0; istackC++) {
                int icol = stackC[istackC];
                stackC0[istackC] = icol;
                lo0[istackC] = colLower[icol];
                up0[istackC] = colUpper[icol];
              }
            }
            /* restore all */
            assert(iway == 0);
            if (changedDown)
              changedDown[j] = nstackC;
            for (istackC = nstackC - 1; istackC >= 0; istackC--) {
              int icol = stackC[istackC];
              double oldU = saveU[istackC];
              double oldL = saveL[istackC];
#if FIXED_BOTH_WAYS
              if (doFixedBothWays && intVar[icol]) {
                assert(newLo[icol] == COIN_DBL_MAX);
                fixedBoth[nFixedBoth++] = icol;
                newLo[icol] = colLower[icol];
                newUp[icol] = colUpper[icol];
              }
#endif
              colUpper[icol] = oldU;
              colLower[icol] = oldL;
              columnGap[icol] = oldU - oldL - primalTolerance_;
              markC[icol] = 0;
              if (oldU > 1.0e10)
                markC[icol] |= 8;
              if (oldL < -1.0e10)
                markC[icol] |= 4;
            }
            for (istackR = 0; istackR < nstackR; istackR++) {
              int irow = stackR[istackR];
              minR[irow] = saveMin[istackR];
              maxR[irow] = saveMax[istackR];
              markR[irow] = -1;
            }
          } else {
            if (iway == 1 && feasible == 3) {
              iway = 3;
              /* point back to stack */
              for (istackC = nstackC - 1; istackC >= 0; istackC--) {
                int icol = stackC[istackC];
                markC[icol] = istackC + 100000;
              }
              OsiColCut cc;
              int /*nTot=0,*/ nFix = 0, nInt = 0;
              for (istackC = 1; istackC < nstackC0; istackC++) {
                int icol = stackC0[istackC];
                int istackC1 = markC[icol] - 100000;
                if (istackC1 >= 0) {
                  markC[icol] -= 100000;
                  if (std::min(lo0[istackC], colLower[icol]) > saveL[istackC1] + 1.0e-4) {
                    double saveBound = saveL[istackC1];
                    saveL[istackC1] = std::min(lo0[istackC], colLower[icol]);
#ifdef UPDATE_MINR_MAXR
                    if (saveBound < saveL[istackC1]) {
                      // tighten minR and maxR
                      assert(numberRealChanges < 2 * maxStack);
                      realChange[numberRealChanges] = saveL[istackC1] - saveBound;
                      whichRealChange[numberRealChanges++] = icol;
                    }
#endif
                    if (intVar[icol] /*||!info->inTree*/) {
                      element[nFix] = saveL[istackC1];
                      index[nFix++] = icol;
                      nInt++;
                    }
                    nfixed++;
                  }
                }
              }
              if (nFix) {
                // nTot=nFix;
                cleanBounds(-1, nFix, index, element);
                cc.setLbs(nFix, index, element);
                nFix = 0;
              }
              for (istackC = 1; istackC < nstackC0; istackC++) {
                int icol = stackC0[istackC];
                int istackC1 = markC[icol] - 100000;
                if (istackC1 >= 0) {
                  markC[icol] -= 100000;
                  if (std::max(up0[istackC], colUpper[icol]) < saveU[istackC1] - 1.0e-4) {
                    double saveBound = saveU[istackC1];
                    saveU[istackC1] = std::max(up0[istackC], colUpper[icol]);
#ifdef UPDATE_MINR_MAXR
                    if (saveBound < saveU[istackC1]) {
                      // tighten minR and maxR
                      assert(numberRealChanges < 2 * maxStack);
                      realChange[numberRealChanges] = saveBound - saveU[istackC1];
                      whichRealChange[numberRealChanges++] = icol;
                    }
#endif
                    if (intVar[icol] /*||!info->inTree*/) {
                      element[nFix] = saveU[istackC1];
                      index[nFix++] = icol;
                      nInt++;
                    }
                    nfixed++;
                  }
                }
              }
              if (nFix) {
                // nTot+=nFix;
                cleanBounds(1, nFix, index, element);
                cc.setUbs(nFix, index, element);
              }
              // could tighten continuous as well
              if (nInt) {
#ifdef CGL_DEBUG
                checkBounds(debugger, cc);
#endif
#ifndef CC_AT_END
                cs.insert(cc);
#endif
              }
            }
            double boundChange;
            /* restore all */
            assert(iway > 0);
            if (changedUp)
              changedUp[j] = nstackC;
#if FIXED_BOTH_WAYS
            // #define JUST_ONE_CC
#ifdef JUST_ONE_CC
            CoinPackedVector lbsBoth;
            CoinPackedVector ubsBoth;
#endif
#endif
            for (istackC = nstackC - 1; istackC > 0; istackC--) {
              int icol = stackC[istackC];
              double oldU = saveU[istackC];
              double oldL = saveL[istackC];
#if FIXED_BOTH_WAYS
              if (doFixedBothWays && newLo[icol] != COIN_DBL_MAX) {
                double loDown = newLo[icol];
                double upDown = newUp[icol];
                double loUp = colLower[icol];
                double upUp = colUpper[icol];
                double newLower = oldL;
                double newUpper = oldU;
#if FIXED_BOTH_WAYS > 2
                bool printit = false;
#endif
                if (loDown > oldL + 1.0e-5 && loUp > oldL + 1.0e-5) {
                  newLower = std::min(loDown, loUp);
                  if (intVar[icol])
                    newLower = ceil(newLower - 1.0e-5);
#if FIXED_BOTH_WAYS > 2
                  nToBoundZ++;
#endif
#ifdef JUST_ONE_CC
                  lbsBoth.insert(icol, newLower);
#endif
                  if (upDown < oldU - 1.0e-5 && upUp < oldU - 1.0e-5) {
                    newUpper = std::max(upDown, upUp);
#if FIXED_BOTH_WAYS > 2
                    nToBoundZ++;
#endif
                    if (intVar[icol])
                      newUpper = floor(newUpper + 1.0e-5);
                  }
                } else if (upDown < oldU - 1.0e-5 && upUp < oldU - 1.0e-5) {
                  newUpper = std::max(upDown, upUp);
                  if (intVar[icol])
                    newUpper = floor(newUpper + 1.0e-5);
#if FIXED_BOTH_WAYS > 2
                  nToBoundZ++;
#endif
#ifdef JUST_ONE_CC
                  ubsBoth.insert(icol, newUpper);
#endif
                }
#if FIXED_BOTH_WAYS > 2
                if (printit)
                  printf("Col %d bnds on %d down %g,%g up %g,%g orig %g,%g\n",
                    j, icol, loDown, upDown, loUp,
                    upUp, oldL, oldU);
#endif
#ifdef UPDATE_MINR_MAXR
                if (oldL < newLower) {
                  // tighten minR and maxR
                  assert(numberRealChanges < 2 * maxStack);
                  realChange[numberRealChanges] = newLower - oldL;
                  whichRealChange[numberRealChanges++] = icol;
                }
                if (oldU > newUpper) {
                  // tighten minR and maxR
                  assert(numberRealChanges < 2 * maxStack);
                  realChange[numberRealChanges] = newUpper - oldU;
                  whichRealChange[numberRealChanges++] = icol;
                }
#ifndef JUST_ONE_CC
                // temp to check
#ifndef CC_AT_END
                if (oldL < newLower || oldU > newUpper) {
                  OsiColCut cc;
                  CoinPackedVector lbs;
                  CoinPackedVector ubs;
                  lbs.insert(icol, newLower);
                  ubs.insert(icol, newUpper);
                  cc.setLbs(lbs);
                  cc.setUbs(ubs);
                  cs.insert(cc);
                }
#endif
#endif
#endif
                oldL = newLower;
                oldU = newUpper;
              }
#endif
              colUpper[icol] = oldU;
              colLower[icol] = oldL;
              columnGap[icol] = oldU - oldL - primalTolerance_;
              if (oldU > oldL + 1.0e-4) {
                markC[icol] = 0;
                if (oldU > 1.0e10)
                  markC[icol] |= 8;
                if (oldL < -1.0e10)
                  markC[icol] |= 4;
              } else {
                markC[icol] = 3;
              }
            }
#if FIXED_BOTH_WAYS
#ifdef JUST_ONE_CC
#ifndef CC_AT_END
            if (lbsBoth.getNumElements() || ubsBoth.getNumElements()) {
              OsiColCut cc;
              if (lbsBoth.getNumElements())
                cc.setLbs(lbsBoth);
              if (ubsBoth.getNumElements())
                cc.setUbs(ubsBoth);
              cs.insert(cc);
            }
#endif
#endif
#endif
            // first
            {
              double oldU = saveU[0];
              double oldL = saveL[0];
              colUpper[j] = oldU;
              colLower[j] = oldL;
              columnGap[j] = oldU - oldL - primalTolerance_;
              assert(oldU > oldL + 1.0e-4);
              markC[j] = 0;
              if (oldU > 1.0e10)
                markC[j] |= 8;
              if (oldL < -1.0e10)
                markC[j] |= 4;
            }
            for (istackR = 0; istackR < nstackR; istackR++) {
              int irow = stackR[istackR];
              minR[irow] = saveMin[istackR];
              maxR[irow] = saveMax[istackR];
              markR[irow] = -1;
            }
          }
        }
      }
    }
  }
#if 0 // ndef NDEBUG
  if (!ninfeas) {
    // check all column cuts reflected in bounds
    double * lower = colsol;
    double * upper = djs;
    memcpy(lower,colLower,nCols*sizeof(double));
    memcpy(upper,colUpper,nCols*sizeof(double));
    int nCuts = cs.sizeColCuts();
    for (int i=0;i<nCuts;i++) {
      OsiColCut * cut = cs.colCutPtr(i);
      const CoinPackedVector lbs = cut->lbs();
      const CoinPackedVector ubs = cut->ubs();
      int ncols;
      const int * cols;
      const double * values;
      ncols = lbs.getNumElements();
      cols = lbs.getIndices();
      values = lbs.getElements();
      for (int j=0;j<ncols;j++) {
	int jColumn = cols[j];
	double value = values[j];
	lower[jColumn] = std::max(lower[jColumn],value);
      }
      ncols = ubs.getNumElements();
      cols = ubs.getIndices();
      values = ubs.getElements();
      for (int j=0;j<ncols;j++) {
	int jColumn = cols[j];
	double value = values[j];
	upper[jColumn] = std::min(upper[jColumn],value);
      }
    }
    // arrays should be same
    for (int i=0;i<nCols;i++) {
      assert (lower[i]==colLower[i]);
      assert (upper[i]==colUpper[i]);
    }
  }
#endif
#if 0 // def CC_AT_END
  if (!ninfeas) {
    // One big cc
    const double * originalLower = si.getColLower();
    const double * originalUpper = si.getColUpper();
    CoinPackedVector lbs;
    CoinPackedVector ubs;
    for (int j=0;j<nCols;j++) {
      // try without later
      if (intVar[j]) {
	if (colLower[j]>originalLower[j]) {
	  lbs.insert(j,colLower[j]);
	}
	if (colUpper[j]<originalUpper[j]) {
	  ubs.insert(j,colUpper[j]);
	}
      }
    }
    int nCuts = 0;
    OsiColCut cc;
    if (lbs.getNumElements()) {
      cc.setLbs(lbs);
      nCuts = lbs.getNumElements();
    }
    if (ubs.getNumElements()) {
      cc.setUbs(ubs);
      nCuts += ubs.getNumElements();
    }
    if (nCuts) {
      cs.insert(cc);
    }
  }
#endif
#ifndef ONE_ARRAY
  delete[] stackC0;
  delete[] lo0;
  delete[] up0;
  delete[] columnGap;
  delete[] markC;
  delete[] stackC;
  delete[] stackR;
  delete[] saveL;
  delete[] saveU;
  delete[] saveMin;
  delete[] saveMax;
  delete[] index;
  delete[] element;
  delete[] djs;
  delete[] largestPositiveInRow;
  delete[] largestNegativeInRow;
#endif
#if FIXED_BOTH_WAYS
  delete[] saveFColLower;
#endif
  delete[] colsol;
  return (ninfeas);
}
