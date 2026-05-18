// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cstdio>

#include "CoinPragma.hpp"
#include "CoinIndexedVector.hpp"
#include "CoinPackedVector.hpp"
#include "CoinHelperFunctions.hpp"

#include "ClpSimplex.hpp"
#include "ClpFactorization.hpp"
// at end to get min/max!
#include "ClpPlusMinusOneMatrix.hpp"
#include "ClpMessage.hpp"
#ifdef CLP_PLUS_ONE_MATRIX
static int oneitcount[13] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
static void oneit(int i)
{
  if (!oneitcount[i])
    printf("Plus ones for call %d\n", i);
  oneitcount[i]++;
  oneitcount[12]++;
  if ((oneitcount[12] % 1000) == 0) {
    printf("Plus counts");
    for (int j = 0; j < 12; j++)
      printf(" %d", oneitcount[j]);
    printf("\n");
  }
}
#endif
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpPlusMinusOneMatrix::ClpPlusMinusOneMatrix()
  : ClpMatrixBase()
{
  setType(12);
  matrix_ = NULL;
  startPositive_ = NULL;
  startNegative_ = NULL;
  lengths_ = NULL;
  indices_ = NULL;
  numberRows_ = 0;
  numberColumns_ = 0;
#ifdef CLP_PLUS_ONE_MATRIX
  otherFlags_ = 0;
#endif
  columnOrdered_ = true;
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpPlusMinusOneMatrix::ClpPlusMinusOneMatrix(const ClpPlusMinusOneMatrix &rhs)
  : ClpMatrixBase(rhs)
{
  matrix_ = NULL;
  startPositive_ = NULL;
  startNegative_ = NULL;
  lengths_ = NULL;
  indices_ = NULL;
  numberRows_ = rhs.numberRows_;
  numberColumns_ = rhs.numberColumns_;
#ifdef CLP_PLUS_ONE_MATRIX
  otherFlags_ = rhs.otherFlags_;
#endif
  columnOrdered_ = rhs.columnOrdered_;
  if (numberColumns_) {
    CoinBigIndex numberElements = rhs.startPositive_[numberColumns_];
    indices_ = new int[numberElements];
    CoinMemcpyN(rhs.indices_, numberElements, indices_);
    startPositive_ = new CoinBigIndex[numberColumns_ + 1];
    CoinMemcpyN(rhs.startPositive_, (numberColumns_ + 1), startPositive_);
    startNegative_ = new CoinBigIndex[numberColumns_];
    CoinMemcpyN(rhs.startNegative_, numberColumns_, startNegative_);
  }
  int numberRows = getNumRows();
  if (rhs.rhsOffset_ && numberRows) {
    rhsOffset_ = ClpCopyOfArray(rhs.rhsOffset_, numberRows);
  } else {
    rhsOffset_ = NULL;
  }
}
// Constructor from arrays
ClpPlusMinusOneMatrix::ClpPlusMinusOneMatrix(int numberRows, int numberColumns,
  bool columnOrdered, const int *indices,
  const CoinBigIndex *startPositive,
  const CoinBigIndex *startNegative)
  : ClpMatrixBase()
{
  setType(12);
  matrix_ = NULL;
  lengths_ = NULL;
  numberRows_ = numberRows;
  numberColumns_ = numberColumns;
  columnOrdered_ = columnOrdered;
#ifdef CLP_PLUS_ONE_MATRIX
  otherFlags_ = 0;
#endif
  int numberMajor = (columnOrdered_) ? numberColumns_ : numberRows_;
  CoinBigIndex numberElements = startPositive[numberMajor];
  startPositive_ = ClpCopyOfArray(startPositive, numberMajor + 1);
  startNegative_ = ClpCopyOfArray(startNegative, numberMajor);
  indices_ = ClpCopyOfArray(indices, numberElements);
  // Check valid
  checkValid(false);
}

ClpPlusMinusOneMatrix::ClpPlusMinusOneMatrix(const CoinPackedMatrix &rhs)
  : ClpMatrixBase()
{
  setType(12);
  matrix_ = NULL;
  startPositive_ = NULL;
  startNegative_ = NULL;
  lengths_ = NULL;
  indices_ = NULL;
  int iColumn;
  assert(rhs.isColOrdered());
  // get matrix data pointers
  const int *row = rhs.getIndices();
  const CoinBigIndex *columnStart = rhs.getVectorStarts();
  const int *columnLength = rhs.getVectorLengths();
  const double *elementByColumn = rhs.getElements();
  numberColumns_ = rhs.getNumCols();
  numberRows_ = -1;
  indices_ = new int[rhs.getNumElements()];
  startPositive_ = new CoinBigIndex[numberColumns_ + 1];
  startNegative_ = new CoinBigIndex[numberColumns_];
  int *temp = new int[rhs.getNumRows()];
  CoinBigIndex j = 0;
  CoinBigIndex numberGoodP = 0;
  CoinBigIndex numberGoodM = 0;
  CoinBigIndex numberBad = 0;
  for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
    CoinBigIndex k;
    int iNeg = 0;
    startPositive_[iColumn] = j;
    for (k = columnStart[iColumn]; k < columnStart[iColumn] + columnLength[iColumn];
         k++) {
      int iRow;
      if (fabs(elementByColumn[k] - 1.0) < 1.0e-10) {
        iRow = row[k];
        numberRows_ = std::max(numberRows_, iRow);
        indices_[j++] = iRow;
        numberGoodP++;
      } else if (fabs(elementByColumn[k] + 1.0) < 1.0e-10) {
        iRow = row[k];
        numberRows_ = std::max(numberRows_, iRow);
        temp[iNeg++] = iRow;
        numberGoodM++;
      } else {
        numberBad++;
      }
    }
    // move negative
    startNegative_[iColumn] = j;
    for (k = 0; k < iNeg; k++) {
      indices_[j++] = temp[k];
    }
  }
  startPositive_[numberColumns_] = j;
  delete[] temp;
  if (numberBad) {
    delete[] indices_;
    indices_ = NULL;
    numberRows_ = 0;
    numberColumns_ = 0;
    delete[] startPositive_;
    delete[] startNegative_;
    // Put in statistics
    startPositive_ = new CoinBigIndex[3];
    startPositive_[0] = numberGoodP;
    startPositive_[1] = numberGoodM;
    startPositive_[2] = numberBad;
    startNegative_ = NULL;
  } else {
    numberRows_++; //  correct
    // but number should be same as rhs
    assert(numberRows_ <= rhs.getNumRows());
    numberRows_ = rhs.getNumRows();
    columnOrdered_ = true;
  }
  // Check valid
  if (!numberBad)
    checkValid(false);
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
ClpPlusMinusOneMatrix::~ClpPlusMinusOneMatrix()
{
  delete matrix_;
  delete[] startPositive_;
  delete[] startNegative_;
  delete[] lengths_;
  delete[] indices_;
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpPlusMinusOneMatrix &
ClpPlusMinusOneMatrix::operator=(const ClpPlusMinusOneMatrix &rhs)
{
  if (this != &rhs) {
    ClpMatrixBase::operator=(rhs);
    delete matrix_;
    delete[] startPositive_;
    delete[] startNegative_;
    delete[] lengths_;
    delete[] indices_;
    matrix_ = NULL;
    startPositive_ = NULL;
    lengths_ = NULL;
    indices_ = NULL;
    numberRows_ = rhs.numberRows_;
    numberColumns_ = rhs.numberColumns_;
    columnOrdered_ = rhs.columnOrdered_;
#ifdef CLP_PLUS_ONE_MATRIX
    otherFlags_ = rhs.otherFlags_;
#endif
    if (numberColumns_) {
      CoinBigIndex numberElements = rhs.startPositive_[numberColumns_];
      indices_ = new int[numberElements];
      CoinMemcpyN(rhs.indices_, numberElements, indices_);
      startPositive_ = new CoinBigIndex[numberColumns_ + 1];
      CoinMemcpyN(rhs.startPositive_, (numberColumns_ + 1), startPositive_);
      startNegative_ = new CoinBigIndex[numberColumns_];
      CoinMemcpyN(rhs.startNegative_, numberColumns_, startNegative_);
    }
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpMatrixBase *ClpPlusMinusOneMatrix::clone() const
{
  return new ClpPlusMinusOneMatrix(*this);
}
/* Subset clone (without gaps).  Duplicates are allowed
   and order is as given */
ClpMatrixBase *
ClpPlusMinusOneMatrix::subsetClone(int numberRows, const int *whichRows,
  int numberColumns,
  const int *whichColumns) const
{
  return new ClpPlusMinusOneMatrix(*this, numberRows, whichRows,
    numberColumns, whichColumns);
}
/* Subset constructor (without gaps).  Duplicates are allowed
   and order is as given */
ClpPlusMinusOneMatrix::ClpPlusMinusOneMatrix(
  const ClpPlusMinusOneMatrix &rhs,
  int numberRows, const int *whichRow,
  int numberColumns, const int *whichColumn)
  : ClpMatrixBase(rhs)
{
  matrix_ = NULL;
  startPositive_ = NULL;
  startNegative_ = NULL;
  lengths_ = NULL;
  indices_ = NULL;
  numberRows_ = 0;
  numberColumns_ = 0;
  columnOrdered_ = rhs.columnOrdered_;
#ifdef CLP_PLUS_ONE_MATRIX
  otherFlags_ = rhs.otherFlags_;
#endif
  if (numberRows <= 0 || numberColumns <= 0) {
    startPositive_ = new CoinBigIndex[1];
    startPositive_[0] = 0;
  } else {
    numberColumns_ = numberColumns;
    numberRows_ = numberRows;
    const int *index1 = rhs.indices_;
    CoinBigIndex *startPositive1 = rhs.startPositive_;

    int numberMinor = (!columnOrdered_) ? numberColumns_ : numberRows_;
    int numberMajor = (columnOrdered_) ? numberColumns_ : numberRows_;
    int numberMinor1 = (!columnOrdered_) ? rhs.numberColumns_ : rhs.numberRows_;
    int numberMajor1 = (columnOrdered_) ? rhs.numberColumns_ : rhs.numberRows_;
    // Also swap incoming if not column ordered
    if (!columnOrdered_) {
      int temp1 = numberRows;
      numberRows = numberColumns;
      numberColumns = temp1;
      const int *temp2;
      temp2 = whichRow;
      whichRow = whichColumn;
      whichColumn = temp2;
    }
    // Throw exception if rhs empty
    if (numberMajor1 <= 0 || numberMinor1 <= 0)
      throw CoinError("empty rhs", "subset constructor", "ClpPlusMinusOneMatrix");
    // Array to say if an old row is in new copy
    int *newRow = new int[numberMinor1];
    int iRow;
    for (iRow = 0; iRow < numberMinor1; iRow++)
      newRow[iRow] = -1;
    // and array for duplicating rows
    int *duplicateRow = new int[numberMinor];
    int numberBad = 0;
    for (iRow = 0; iRow < numberMinor; iRow++) {
      duplicateRow[iRow] = -1;
      int kRow = whichRow[iRow];
      if (kRow >= 0 && kRow < numberMinor1) {
        if (newRow[kRow] < 0) {
          // first time
          newRow[kRow] = iRow;
        } else {
          // duplicate
          int lastRow = newRow[kRow];
          newRow[kRow] = iRow;
          duplicateRow[iRow] = lastRow;
        }
      } else {
        // bad row
        numberBad++;
      }
    }

    if (numberBad)
      throw CoinError("bad minor entries",
        "subset constructor", "ClpPlusMinusOneMatrix");
    // now get size and check columns
    CoinBigIndex size = 0;
    int iColumn;
    numberBad = 0;
    for (iColumn = 0; iColumn < numberMajor; iColumn++) {
      int kColumn = whichColumn[iColumn];
      if (kColumn >= 0 && kColumn < numberMajor1) {
        CoinBigIndex i;
        for (i = startPositive1[kColumn]; i < startPositive1[kColumn + 1]; i++) {
          int kRow = index1[i];
          kRow = newRow[kRow];
          while (kRow >= 0) {
            size++;
            kRow = duplicateRow[kRow];
          }
        }
      } else {
        // bad column
        numberBad++;
        printf("%d %d %d %d\n", iColumn, numberMajor, numberMajor1, kColumn);
      }
    }
    if (numberBad)
      throw CoinError("bad major entries",
        "subset constructor", "ClpPlusMinusOneMatrix");
    // now create arrays
    startPositive_ = new CoinBigIndex[numberMajor + 1];
    startNegative_ = new CoinBigIndex[numberMajor];
    indices_ = new int[size];
    // and fill them
    size = 0;
    startPositive_[0] = 0;
    CoinBigIndex *startNegative1 = rhs.startNegative_;
    for (iColumn = 0; iColumn < numberMajor; iColumn++) {
      int kColumn = whichColumn[iColumn];
      CoinBigIndex i;
      for (i = startPositive1[kColumn]; i < startNegative1[kColumn]; i++) {
        int kRow = index1[i];
        kRow = newRow[kRow];
        while (kRow >= 0) {
          indices_[size++] = kRow;
          kRow = duplicateRow[kRow];
        }
      }
      startNegative_[iColumn] = size;
      for (; i < startPositive1[kColumn + 1]; i++) {
        int kRow = index1[i];
        kRow = newRow[kRow];
        while (kRow >= 0) {
          indices_[size++] = kRow;
          kRow = duplicateRow[kRow];
        }
      }
      startPositive_[iColumn + 1] = size;
    }
    delete[] newRow;
    delete[] duplicateRow;
  }
  // Check valid
  checkValid(false);
}

/* Returns a new matrix in reverse order without gaps */
ClpMatrixBase *
ClpPlusMinusOneMatrix::reverseOrderedCopy() const
{
  int numberMinor = (!columnOrdered_) ? numberColumns_ : numberRows_;
  int numberMajor = (columnOrdered_) ? numberColumns_ : numberRows_;
  // count number in each row/column
  CoinBigIndex *tempP = new CoinBigIndex[numberMinor];
  CoinBigIndex *tempN = new CoinBigIndex[numberMinor];
  memset(tempP, 0, numberMinor * sizeof(CoinBigIndex));
  memset(tempN, 0, numberMinor * sizeof(CoinBigIndex));
  CoinBigIndex j = 0;
  int i;
  for (i = 0; i < numberMajor; i++) {
    for (; j < startNegative_[i]; j++) {
      int iRow = indices_[j];
      tempP[iRow]++;
    }
    for (; j < startPositive_[i + 1]; j++) {
      int iRow = indices_[j];
      tempN[iRow]++;
    }
  }
  int *newIndices = new int[startPositive_[numberMajor]];
  CoinBigIndex *newP = new CoinBigIndex[numberMinor + 1];
  CoinBigIndex *newN = new CoinBigIndex[numberMinor];
  int iRow;
  j = 0;
  // do starts
  for (iRow = 0; iRow < numberMinor; iRow++) {
    newP[iRow] = j;
    j += tempP[iRow];
    tempP[iRow] = newP[iRow];
    newN[iRow] = j;
    j += tempN[iRow];
    tempN[iRow] = newN[iRow];
  }
  newP[numberMinor] = j;
  j = 0;
  for (i = 0; i < numberMajor; i++) {
    for (; j < startNegative_[i]; j++) {
      int iRow = indices_[j];
      CoinBigIndex put = tempP[iRow];
      newIndices[put++] = i;
      tempP[iRow] = put;
    }
    for (; j < startPositive_[i + 1]; j++) {
      int iRow = indices_[j];
      CoinBigIndex put = tempN[iRow];
      newIndices[put++] = i;
      tempN[iRow] = put;
    }
  }
  delete[] tempP;
  delete[] tempN;
  ClpPlusMinusOneMatrix *newCopy = new ClpPlusMinusOneMatrix();
  newCopy->passInCopy(numberMinor, numberMajor,
    !columnOrdered_, newIndices, newP, newN);
  return newCopy;
}
//static bool doPlusOnes=true;
//unscaled versions
void ClpPlusMinusOneMatrix::times(double scalar,
  const double *COIN_RESTRICT x, double *COIN_RESTRICT y) const
{
  int numberMajor = (columnOrdered_) ? numberColumns_ : numberRows_;
  int i;
  CoinBigIndex j;
  assert(columnOrdered_);
#ifdef CLP_PLUS_ONE_MATRIX
  if ((otherFlags_ & 1) == 0 || !doPlusOnes) {
#endif
    for (i = 0; i < numberMajor; i++) {
      double value = scalar * x[i];
      if (value) {
        for (j = startPositive_[i]; j < startNegative_[i]; j++) {
          int iRow = indices_[j];
          y[iRow] += value;
        }
        for (; j < startPositive_[i + 1]; j++) {
          int iRow = indices_[j];
          y[iRow] -= value;
        }
      }
    }
#ifdef CLP_PLUS_ONE_MATRIX
  } else {
    // plus one
    oneit(0);
    for (i = 0; i < numberMajor; i++) {
      double value = scalar * x[i];
      if (value) {
        for (j = startPositive_[i]; j < startPositive_[i + 1]; j++) {
          int iRow = indices_[j];
          y[iRow] += value;
        }
      }
    }
  }
#endif
}
void ClpPlusMinusOneMatrix::transposeTimes(double scalar,
  const double *COIN_RESTRICT x, double *COIN_RESTRICT y) const
{
  int numberMajor = (columnOrdered_) ? numberColumns_ : numberRows_;
  int i;
  CoinBigIndex j = 0;
  assert(columnOrdered_);
#ifdef CLP_PLUS_ONE_MATRIX
  if ((otherFlags_ & 1) == 0 || !doPlusOnes) {
#endif
    for (i = 0; i < numberMajor; i++) {
      double value = 0.0;
      for (; j < startNegative_[i]; j++) {
        int iRow = indices_[j];
        value += x[iRow];
      }
      for (; j < startPositive_[i + 1]; j++) {
        int iRow = indices_[j];
        value -= x[iRow];
      }
      y[i] += scalar * value;
    }
#ifdef CLP_PLUS_ONE_MATRIX
  } else {
    // plus one
    oneit(1);
    for (i = 0; i < numberMajor; i++) {
      double value = 0.0;
      for (; j < startPositive_[i + 1]; j++) {
        int iRow = indices_[j];
        value += x[iRow];
      }
      y[i] += scalar * value;
    }
  }
#endif
}
void ClpPlusMinusOneMatrix::times(double scalar,
  const double *COIN_RESTRICT x, double *COIN_RESTRICT y,
  const double * /*rowScale*/,
  const double * /*columnScale*/) const
{
  // we know it is not scaled
  times(scalar, x, y);
}
void ClpPlusMinusOneMatrix::transposeTimes(double scalar,
  const double *COIN_RESTRICT x, double *COIN_RESTRICT y,
  const double * /*rowScale*/,
  const double * /*columnScale*/,
  double * /*spare*/) const
{
  // we know it is not scaled
  transposeTimes(scalar, x, y);
}
/* Return <code>x * A + y</code> in <code>z</code>.
	Squashes small elements and knows about ClpSimplex */
void ClpPlusMinusOneMatrix::transposeTimes(const ClpSimplex *model, double scalar,
  const CoinIndexedVector *rowArray,
  CoinIndexedVector *y,
  CoinIndexedVector *columnArray) const
{
  // we know it is not scaled
  columnArray->clear();
  double *COIN_RESTRICT pi = rowArray->denseVector();
  int numberNonZero = 0;
  int *COIN_RESTRICT index = columnArray->getIndices();
  double *COIN_RESTRICT array = columnArray->denseVector();
  int numberInRowArray = rowArray->getNumElements();
  // maybe I need one in OsiSimplex
  double zeroTolerance = model->zeroTolerance();
  int numberRows = model->numberRows();
  bool packed = rowArray->packedMode();
#ifndef NO_RTTI
  ClpPlusMinusOneMatrix *rowCopy = dynamic_cast< ClpPlusMinusOneMatrix * >(model->rowCopy());
#else
  ClpPlusMinusOneMatrix *rowCopy = static_cast< ClpPlusMinusOneMatrix * >(model->rowCopy());
#endif
  double factor = 0.3;
  // We may not want to do by row if there may be cache problems
  int numberColumns = model->numberColumns();
  // It would be nice to find L2 cache size - for moment 512K
  // Be slightly optimistic
  if (numberColumns * sizeof(double) > 1000000) {
    if (numberRows * 10 < numberColumns)
      factor = 0.1;
    else if (numberRows * 4 < numberColumns)
      factor = 0.15;
    else if (numberRows * 2 < numberColumns)
      factor = 0.2;
  }
  if (numberInRowArray > factor * numberRows || !rowCopy) {
    assert(!y->getNumElements());
    // do by column
    // Need to expand if packed mode
    int iColumn;
    CoinBigIndex j = 0;
    assert(columnOrdered_);
    if (packed) {
      // need to expand pi into y
      assert(y->capacity() >= numberRows);
      double *COIN_RESTRICT piOld = pi;
      pi = y->denseVector();
      const int *COIN_RESTRICT whichRow = rowArray->getIndices();
      int i;
      // modify pi so can collapse to one loop
      for (i = 0; i < numberInRowArray; i++) {
        int iRow = whichRow[i];
        pi[iRow] = scalar * piOld[i];
      }
#ifdef CLP_PLUS_ONE_MATRIX
      if ((otherFlags_ & 1) == 0 || !doPlusOnes) {
#endif
        for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
          double value = 0.0;
          for (; j < startNegative_[iColumn]; j++) {
            int iRow = indices_[j];
            value += pi[iRow];
          }
          for (; j < startPositive_[iColumn + 1]; j++) {
            int iRow = indices_[j];
            value -= pi[iRow];
          }
          if (fabs(value) > zeroTolerance) {
            array[numberNonZero] = value;
            index[numberNonZero++] = iColumn;
          }
        }
#ifdef CLP_PLUS_ONE_MATRIX
      } else {
        // plus one
        oneit(2);
        for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
          double value = 0.0;
          for (; j < startPositive_[iColumn + 1]; j++) {
            int iRow = indices_[j];
            value += pi[iRow];
          }
          if (fabs(value) > zeroTolerance) {
            array[numberNonZero] = value;
            index[numberNonZero++] = iColumn;
          }
        }
      }
#endif
      for (i = 0; i < numberInRowArray; i++) {
        int iRow = whichRow[i];
        pi[iRow] = 0.0;
      }
    } else {
      for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
        double value = 0.0;
        for (; j < startNegative_[iColumn]; j++) {
          int iRow = indices_[j];
          value += pi[iRow];
        }
        for (; j < startPositive_[iColumn + 1]; j++) {
          int iRow = indices_[j];
          value -= pi[iRow];
        }
        value *= scalar;
        if (fabs(value) > zeroTolerance) {
          index[numberNonZero++] = iColumn;
          array[iColumn] = value;
        }
      }
    }
    columnArray->setNumElements(numberNonZero);
  } else {
    // do by row
    rowCopy->transposeTimesByRow(model, scalar, rowArray, y, columnArray);
  }
}
/* Return <code>x * A + y</code> in <code>z</code>.
	Squashes small elements and knows about ClpSimplex */
void ClpPlusMinusOneMatrix::transposeTimesByRow(const ClpSimplex *model, double scalar,
  const CoinIndexedVector *rowArray,
  CoinIndexedVector *y,
  CoinIndexedVector *columnArray) const
{
  columnArray->clear();
  double *COIN_RESTRICT pi = rowArray->denseVector();
  int numberNonZero = 0;
  int *COIN_RESTRICT index = columnArray->getIndices();
  double *COIN_RESTRICT array = columnArray->denseVector();
  int numberInRowArray = rowArray->getNumElements();
  // maybe I need one in OsiSimplex
  double zeroTolerance = model->zeroTolerance();
  const int *COIN_RESTRICT column = indices_;
  const CoinBigIndex *COIN_RESTRICT startPositive = startPositive_;
  const CoinBigIndex *COIN_RESTRICT startNegative = startNegative_;
  const int *COIN_RESTRICT whichRow = rowArray->getIndices();
  bool packed = rowArray->packedMode();
  if (numberInRowArray > 2) {
    // do by rows
    int iRow;
    double *COIN_RESTRICT markVector = y->denseVector(); // probably empty .. but
    int numberOriginal = 0;
    int i;
    if (packed) {
      CoinBigIndex numberCovered = 0;
      int numberColumns = getNumCols();
      bool sparse = true;
      int target = 1 * numberColumns;
      for (int i = 0; i < numberInRowArray; i++) {
        int iRow = whichRow[i];
        numberCovered += startPositive[iRow + 1] - startPositive[iRow];
        if (numberCovered > target) {
          sparse = false;
          break;
        }
      }
      numberNonZero = 0;
      if (sparse) {
        // and set up mark as char array
        char *COIN_RESTRICT marked = reinterpret_cast< char * >(index + columnArray->capacity());
        double *COIN_RESTRICT array2 = y->denseVector();
#ifdef CLP_DEBUG
        int numberColumns = model->numberColumns();
        for (int i = 0; i < numberColumns; i++) {
          assert(!marked[i]);
          assert(!array2[i]);
        }
#endif
#ifdef CLP_PLUS_ONE_MATRIX
        if ((otherFlags_ & 1) == 0 || !doPlusOnes) {
#endif
          for (int i = 0; i < numberInRowArray; i++) {
            iRow = whichRow[i];
            double value = pi[i] * scalar;
            CoinBigIndex j;
            for (j = startPositive[iRow]; j < startNegative[iRow]; j++) {
              int iColumn = column[j];
              if (!marked[iColumn]) {
                marked[iColumn] = 1;
                index[numberNonZero++] = iColumn;
              }
              array2[iColumn] += value;
            }
            for (j = startNegative[iRow]; j < startPositive[iRow + 1]; j++) {
              int iColumn = column[j];
              if (!marked[iColumn]) {
                marked[iColumn] = 1;
                index[numberNonZero++] = iColumn;
              }
              array2[iColumn] -= value;
            }
          }
#ifdef CLP_PLUS_ONE_MATRIX
        } else {
          // plus one
          oneit(4);
          for (int i = 0; i < numberInRowArray; i++) {
            iRow = whichRow[i];
            double value = pi[i] * scalar;
            CoinBigIndex j;
            for (j = startPositive[iRow]; j < startPositive[iRow + 1]; j++) {
              int iColumn = column[j];
              if (!marked[iColumn]) {
                marked[iColumn] = 1;
                index[numberNonZero++] = iColumn;
              }
              array2[iColumn] += value;
            }
          }
        }
#endif
        // get rid of tiny values and zero out marked
        numberOriginal = numberNonZero;
        numberNonZero = 0;
        for (int i = 0; i < numberOriginal; i++) {
          int iColumn = index[i];
          if (marked[iColumn]) {
            double value = array2[iColumn];
            array2[iColumn] = 0.0;
            marked[iColumn] = 0;
            if (fabs(value) > zeroTolerance) {
              array[numberNonZero] = value;
              index[numberNonZero++] = iColumn;
            }
          }
        }
      } else {
        // not sparse
#ifdef CLP_PLUS_ONE_MATRIX
        if ((otherFlags_ & 1) == 0 || !doPlusOnes) {
#endif
          for (int i = 0; i < numberInRowArray; i++) {
            iRow = whichRow[i];
            double value = pi[i] * scalar;
            CoinBigIndex j;
            for (j = startPositive[iRow]; j < startNegative[iRow]; j++) {
              int iColumn = column[j];
              array[iColumn] += value;
            }
            for (j = startNegative[iRow]; j < startPositive[iRow + 1]; j++) {
              int iColumn = column[j];
              array[iColumn] -= value;
            }
          }
#ifdef CLP_PLUS_ONE_MATRIX
        } else {
          // plus one
          oneit(5);
          for (int i = 0; i < numberInRowArray; i++) {
            iRow = whichRow[i];
            double value = pi[i] * scalar;
            CoinBigIndex j;
            for (j = startPositive[iRow]; j < startPositive[iRow + 1]; j++) {
              int iColumn = column[j];
              array[iColumn] += value;
            }
          }
        }
#endif
        // get rid of tiny values and count
        for (int i = 0; i < numberColumns; i++) {
          double value = array[i];
          if (value) {
            array[i] = 0.0;
            if (fabs(value) > zeroTolerance) {
              array[numberNonZero] = value;
              index[numberNonZero++] = i;
            }
          }
        }
      }
    } else {
      numberNonZero = 0;
      // and set up mark as char array
      char *COIN_RESTRICT marked = reinterpret_cast< char * >(markVector);
      for (i = 0; i < numberOriginal; i++) {
        int iColumn = index[i];
        marked[iColumn] = 0;
      }
      for (i = 0; i < numberInRowArray; i++) {
        iRow = whichRow[i];
        double value = pi[iRow] * scalar;
        CoinBigIndex j;
        for (j = startPositive[iRow]; j < startNegative[iRow]; j++) {
          int iColumn = column[j];
          if (!marked[iColumn]) {
            marked[iColumn] = 1;
            index[numberNonZero++] = iColumn;
          }
          array[iColumn] += value;
        }
        for (j = startNegative[iRow]; j < startPositive[iRow + 1]; j++) {
          int iColumn = column[j];
          if (!marked[iColumn]) {
            marked[iColumn] = 1;
            index[numberNonZero++] = iColumn;
          }
          array[iColumn] -= value;
        }
      }
      // get rid of tiny values and zero out marked
      numberOriginal = numberNonZero;
      numberNonZero = 0;
      for (i = 0; i < numberOriginal; i++) {
        int iColumn = index[i];
        marked[iColumn] = 0;
        if (fabs(array[iColumn]) > zeroTolerance) {
          index[numberNonZero++] = iColumn;
        } else {
          array[iColumn] = 0.0;
        }
      }
    }
  } else if (numberInRowArray == 2) {
    /* do by rows when two rows (do longer first when not packed
             and shorter first if packed */
    int iRow0 = whichRow[0];
    int iRow1 = whichRow[1];
    CoinBigIndex j;
    if (packed) {
      double pi0 = pi[0];
      double pi1 = pi[1];
      if (startPositive[iRow0 + 1] - startPositive[iRow0] > startPositive[iRow1 + 1] - startPositive[iRow1]) {
        int temp = iRow0;
        iRow0 = iRow1;
        iRow1 = temp;
        pi0 = pi1;
        pi1 = pi[0];
      }
      // and set up mark as char array
      char *COIN_RESTRICT marked = reinterpret_cast< char * >(index + columnArray->capacity());
      int *COIN_RESTRICT lookup = y->getIndices();
      double value = pi0 * scalar;
      int numberOriginal;
#ifdef CLP_PLUS_ONE_MATRIX
      if ((otherFlags_ & 1) == 0 || !doPlusOnes) {
#endif
        for (j = startPositive[iRow0]; j < startNegative[iRow0]; j++) {
          int iColumn = column[j];
          array[numberNonZero] = value;
          marked[iColumn] = 1;
          lookup[iColumn] = numberNonZero;
          index[numberNonZero++] = iColumn;
        }
        for (j = startNegative[iRow0]; j < startPositive[iRow0 + 1]; j++) {
          int iColumn = column[j];
          array[numberNonZero] = -value;
          marked[iColumn] = 1;
          lookup[iColumn] = numberNonZero;
          index[numberNonZero++] = iColumn;
        }
        numberOriginal = numberNonZero;
        value = pi1 * scalar;
        for (j = startPositive[iRow1]; j < startNegative[iRow1]; j++) {
          int iColumn = column[j];
          if (marked[iColumn]) {
            int iLookup = lookup[iColumn];
            array[iLookup] += value;
          } else {
            if (fabs(value) > zeroTolerance) {
              array[numberNonZero] = value;
              index[numberNonZero++] = iColumn;
            }
          }
        }
        for (j = startNegative[iRow1]; j < startPositive[iRow1 + 1]; j++) {
          int iColumn = column[j];
          if (marked[iColumn]) {
            int iLookup = lookup[iColumn];
            array[iLookup] -= value;
          } else {
            if (fabs(value) > zeroTolerance) {
              array[numberNonZero] = -value;
              index[numberNonZero++] = iColumn;
            }
          }
        }
#ifdef CLP_PLUS_ONE_MATRIX
      } else {
        // plus one
        oneit(7);
        for (j = startPositive[iRow0]; j < startPositive[iRow0 + 1]; j++) {
          int iColumn = column[j];
          array[numberNonZero] = value;
          marked[iColumn] = 1;
          lookup[iColumn] = numberNonZero;
          index[numberNonZero++] = iColumn;
        }
        numberOriginal = numberNonZero;
        value = pi1 * scalar;
        for (j = startPositive[iRow1]; j < startPositive[iRow1 + 1]; j++) {
          int iColumn = column[j];
          if (marked[iColumn]) {
            int iLookup = lookup[iColumn];
            array[iLookup] += value;
          } else {
            if (fabs(value) > zeroTolerance) {
              array[numberNonZero] = value;
              index[numberNonZero++] = iColumn;
            }
          }
        }
      }
#endif
      // get rid of tiny values and zero out marked
      int nDelete = 0;
      for (j = 0; j < numberOriginal; j++) {
        int iColumn = index[j];
        marked[iColumn] = 0;
        if (fabs(array[j]) <= zeroTolerance)
          nDelete++;
      }
      if (nDelete) {
        numberOriginal = numberNonZero;
        numberNonZero = 0;
        for (j = 0; j < numberOriginal; j++) {
          int iColumn = index[j];
          double value = array[j];
          array[j] = 0.0;
          if (fabs(value) > zeroTolerance) {
            array[numberNonZero] = value;
            index[numberNonZero++] = iColumn;
          }
        }
      }
    } else {
      if (startPositive[iRow0 + 1] - startPositive[iRow0] < startPositive[iRow1 + 1] - startPositive[iRow1]) {
        int temp = iRow0;
        iRow0 = iRow1;
        iRow1 = temp;
      }
      int numberOriginal;
      int i;
      numberNonZero = 0;
      double value;
      value = pi[iRow0] * scalar;
      CoinBigIndex j;
      for (j = startPositive[iRow0]; j < startNegative[iRow0]; j++) {
        int iColumn = column[j];
        index[numberNonZero++] = iColumn;
        array[iColumn] = value;
      }
      for (j = startNegative[iRow0]; j < startPositive[iRow0 + 1]; j++) {
        int iColumn = column[j];
        index[numberNonZero++] = iColumn;
        array[iColumn] = -value;
      }
      value = pi[iRow1] * scalar;
      for (j = startPositive[iRow1]; j < startNegative[iRow1]; j++) {
        int iColumn = column[j];
        double value2 = array[iColumn];
        if (value2) {
          value2 += value;
        } else {
          value2 = value;
          index[numberNonZero++] = iColumn;
        }
        array[iColumn] = value2;
      }
      for (j = startNegative[iRow1]; j < startPositive[iRow1 + 1]; j++) {
        int iColumn = column[j];
        double value2 = array[iColumn];
        if (value2) {
          value2 -= value;
        } else {
          value2 = -value;
          index[numberNonZero++] = iColumn;
        }
        array[iColumn] = value2;
      }
      // get rid of tiny values and zero out marked
      numberOriginal = numberNonZero;
      numberNonZero = 0;
      for (i = 0; i < numberOriginal; i++) {
        int iColumn = index[i];
        if (fabs(array[iColumn]) > zeroTolerance) {
          index[numberNonZero++] = iColumn;
        } else {
          array[iColumn] = 0.0;
        }
      }
    }
  } else if (numberInRowArray == 1) {
    // Just one row
    int iRow = rowArray->getIndices()[0];
    numberNonZero = 0;
    double value;
    iRow = whichRow[0];
    CoinBigIndex j;
    if (packed) {
      value = pi[0] * scalar;
      if (fabs(value) > zeroTolerance) {
#ifdef CLP_PLUS_ONE_MATRIX
        if ((otherFlags_ & 1) == 0 || !doPlusOnes) {
#endif
          for (j = startPositive[iRow]; j < startNegative[iRow]; j++) {
            int iColumn = column[j];
            array[numberNonZero] = value;
            index[numberNonZero++] = iColumn;
          }
          for (j = startNegative[iRow]; j < startPositive[iRow + 1]; j++) {
            int iColumn = column[j];
            array[numberNonZero] = -value;
            index[numberNonZero++] = iColumn;
          }
#ifdef CLP_PLUS_ONE_MATRIX
        } else {
          // plus one
          oneit(9);
          for (j = startPositive[iRow]; j < startPositive[iRow + 1]; j++) {
            int iColumn = column[j];
            array[numberNonZero] = value;
            index[numberNonZero++] = iColumn;
          }
        }
#endif
      }
    } else {
      value = pi[iRow] * scalar;
      if (fabs(value) > zeroTolerance) {
        for (j = startPositive[iRow]; j < startNegative[iRow]; j++) {
          int iColumn = column[j];
          array[iColumn] = value;
          index[numberNonZero++] = iColumn;
        }
        for (j = startNegative[iRow]; j < startPositive[iRow + 1]; j++) {
          int iColumn = column[j];
          array[iColumn] = -value;
          index[numberNonZero++] = iColumn;
        }
      }
    }
  }
  columnArray->setNumElements(numberNonZero);
  if (packed)
    columnArray->setPacked();
  y->setNumElements(0);
}
/* Return <code>x *A in <code>z</code> but
   just for indices in y. */
void ClpPlusMinusOneMatrix::subsetTransposeTimes(const ClpSimplex *,
  const CoinIndexedVector *rowArray,
  const CoinIndexedVector *y,
  CoinIndexedVector *columnArray) const
{
  columnArray->clear();
  double *COIN_RESTRICT pi = rowArray->denseVector();
  double *COIN_RESTRICT array = columnArray->denseVector();
  int jColumn;
  int numberToDo = y->getNumElements();
  const int *COIN_RESTRICT which = y->getIndices();
  assert(!rowArray->packedMode());
  columnArray->setPacked();
#ifdef CLP_PLUS_ONE_MATRIX
  if ((otherFlags_ & 1) == 0 || !doPlusOnes) {
#endif
    for (jColumn = 0; jColumn < numberToDo; jColumn++) {
      int iColumn = which[jColumn];
      double value = 0.0;
      CoinBigIndex j = startPositive_[iColumn];
      for (; j < startNegative_[iColumn]; j++) {
        int iRow = indices_[j];
        value += pi[iRow];
      }
      for (; j < startPositive_[iColumn + 1]; j++) {
        int iRow = indices_[j];
        value -= pi[iRow];
      }
      array[jColumn] = value;
    }
#ifdef CLP_PLUS_ONE_MATRIX
  } else {
    // plus one
    oneit(11);
    for (jColumn = 0; jColumn < numberToDo; jColumn++) {
      int iColumn = which[jColumn];
      double value = 0.0;
      CoinBigIndex j = startPositive_[iColumn];
      for (; j < startPositive_[iColumn + 1]; j++) {
        int iRow = indices_[j];
        value += pi[iRow];
      }
      array[jColumn] = value;
    }
  }
#endif
}
/// returns number of elements in column part of basis,
int ClpPlusMinusOneMatrix::countBasis(const int *whichColumn,
  int &numberColumnBasic)
{
  int i;
  CoinBigIndex numberElements = 0;
  for (i = 0; i < numberColumnBasic; i++) {
    int iColumn = whichColumn[i];
    numberElements += startPositive_[iColumn + 1] - startPositive_[iColumn];
  }
  if (numberElements > COIN_INT_MAX) {
    printf("Factorization too large\n");
    abort();
  }
  return static_cast< int >(numberElements);
}
void ClpPlusMinusOneMatrix::fillBasis(ClpSimplex *,
  const int *whichColumn,
  int &numberColumnBasic,
  int *indexRowU, int *start,
  int *rowCount, int *columnCount,
  CoinFactorizationDouble *elementU)
{
  int i;
  CoinBigIndex numberElements = start[0];
  assert(columnOrdered_);
  for (i = 0; i < numberColumnBasic; i++) {
    int iColumn = whichColumn[i];
    CoinBigIndex j = startPositive_[iColumn];
    for (; j < startNegative_[iColumn]; j++) {
      int iRow = indices_[j];
      indexRowU[numberElements] = iRow;
      rowCount[iRow]++;
      elementU[numberElements++] = 1.0;
    }
    for (; j < startPositive_[iColumn + 1]; j++) {
      int iRow = indices_[j];
      indexRowU[numberElements] = iRow;
      rowCount[iRow]++;
      elementU[numberElements++] = -1.0;
    }
    start[i + 1] = static_cast< int >(numberElements);
    columnCount[i] = static_cast< int >(numberElements - start[i]);
  }
  if (numberElements > COIN_INT_MAX) {
    printf("Factorization too large\n");
    abort();
  }
}
/* Unpacks a column into an CoinIndexedvector
 */
void ClpPlusMinusOneMatrix::unpack(const ClpSimplex *,
  CoinIndexedVector *rowArray,
  int iColumn) const
{
  CoinBigIndex j = startPositive_[iColumn];
  for (; j < startNegative_[iColumn]; j++) {
    int iRow = indices_[j];
    rowArray->add(iRow, 1.0);
  }
  for (; j < startPositive_[iColumn + 1]; j++) {
    int iRow = indices_[j];
    rowArray->add(iRow, -1.0);
  }
}
/* Unpacks a column into an CoinIndexedvector
** in packed foramt
Note that model is NOT const.  Bounds and objective could
be modified if doing column generation (just for this variable) */
void ClpPlusMinusOneMatrix::unpackPacked(ClpSimplex *,
  CoinIndexedVector *rowArray,
  int iColumn) const
{
  int *COIN_RESTRICT index = rowArray->getIndices();
  double *COIN_RESTRICT array = rowArray->denseVector();
  int number = 0;
  CoinBigIndex j = startPositive_[iColumn];
  for (; j < startNegative_[iColumn]; j++) {
    int iRow = indices_[j];
    array[number] = 1.0;
    index[number++] = iRow;
  }
  for (; j < startPositive_[iColumn + 1]; j++) {
    int iRow = indices_[j];
    array[number] = -1.0;
    index[number++] = iRow;
  }
  rowArray->setNumElements(number);
  rowArray->setPackedMode(true);
}
/* Adds multiple of a column into an CoinIndexedvector
      You can use quickAdd to add to vector */
void ClpPlusMinusOneMatrix::add(const ClpSimplex *, CoinIndexedVector *rowArray,
  int iColumn, double multiplier) const
{
  CoinBigIndex j = startPositive_[iColumn];
  for (; j < startNegative_[iColumn]; j++) {
    int iRow = indices_[j];
    rowArray->quickAdd(iRow, multiplier);
  }
  for (; j < startPositive_[iColumn + 1]; j++) {
    int iRow = indices_[j];
    rowArray->quickAdd(iRow, -multiplier);
  }
}
/* Adds multiple of a column into an array */
void ClpPlusMinusOneMatrix::add(const ClpSimplex *, double *array,
  int iColumn, double multiplier) const
{
  CoinBigIndex j = startPositive_[iColumn];
  for (; j < startNegative_[iColumn]; j++) {
    int iRow = indices_[j];
    array[iRow] += multiplier;
  }
  for (; j < startPositive_[iColumn + 1]; j++) {
    int iRow = indices_[j];
    array[iRow] -= multiplier;
  }
}

// Return a complete CoinPackedMatrix
CoinPackedMatrix *
ClpPlusMinusOneMatrix::getPackedMatrix() const
{
  if (!matrix_) {
    int numberMinor = (!columnOrdered_) ? numberColumns_ : numberRows_;
    int numberMajor = (columnOrdered_) ? numberColumns_ : numberRows_;
    CoinBigIndex numberElements = startPositive_[numberMajor];
    double *elements = new double[numberElements];
    CoinBigIndex j = 0;
    int i;
    for (i = 0; i < numberMajor; i++) {
      for (; j < startNegative_[i]; j++) {
        elements[j] = 1.0;
      }
      for (; j < startPositive_[i + 1]; j++) {
        elements[j] = -1.0;
      }
    }
    matrix_ = new CoinPackedMatrix(columnOrdered_, numberMinor, numberMajor,
      getNumElements(),
      elements, indices_,
      startPositive_, getVectorLengths());
    delete[] elements;
    delete[] lengths_;
    lengths_ = NULL;
  }
  return matrix_;
}
/* A vector containing the elements in the packed matrix. Note that there
   might be gaps in this list, entries that do not belong to any
   major-dimension vector. To get the actual elements one should look at
   this vector together with vectorStarts and vectorLengths. */
const double *
ClpPlusMinusOneMatrix::getElements() const
{
  if (!matrix_)
    getPackedMatrix();
  return matrix_->getElements();
}

const CoinBigIndex *
ClpPlusMinusOneMatrix::getVectorStarts() const
{
  return startPositive_;
}
/* The lengths of the major-dimension vectors. */
const int *
ClpPlusMinusOneMatrix::getVectorLengths() const
{
  if (!lengths_) {
    int numberMajor = (columnOrdered_) ? numberColumns_ : numberRows_;
    lengths_ = new int[numberMajor];
    int i;
    for (i = 0; i < numberMajor; i++) {
      lengths_[i] = static_cast< int >(startPositive_[i + 1] - startPositive_[i]);
    }
  }
  return lengths_;
}
/* Delete the columns whose indices are listed in <code>indDel</code>. */
void ClpPlusMinusOneMatrix::deleteCols(const int numDel, const int *indDel)
{
  int iColumn;
  CoinBigIndex newSize = startPositive_[numberColumns_];
  ;
  int numberBad = 0;
  // Use array to make sure we can have duplicates
  int *which = new int[numberColumns_];
  memset(which, 0, numberColumns_ * sizeof(int));
  int nDuplicate = 0;
  for (iColumn = 0; iColumn < numDel; iColumn++) {
    int jColumn = indDel[iColumn];
    if (jColumn < 0 || jColumn >= numberColumns_) {
      numberBad++;
    } else {
      newSize -= startPositive_[jColumn + 1] - startPositive_[jColumn];
      if (which[jColumn])
        nDuplicate++;
      else
        which[jColumn] = 1;
    }
  }
  if (numberBad)
    throw CoinError("Indices out of range", "deleteCols", "ClpPlusMinusOneMatrix");
  int newNumber = numberColumns_ - numDel + nDuplicate;
  // Get rid of temporary arrays
  delete[] lengths_;
  lengths_ = NULL;
  delete matrix_;
  matrix_ = NULL;
  CoinBigIndex *newPositive = new CoinBigIndex[newNumber + 1];
  CoinBigIndex *newNegative = new CoinBigIndex[newNumber];
  int *newIndices = new int[newSize];
  newNumber = 0;
  newSize = 0;
  for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
    if (!which[iColumn]) {
      CoinBigIndex start, end;
      CoinBigIndex i;
      start = startPositive_[iColumn];
      end = startNegative_[iColumn];
      newPositive[newNumber] = newSize;
      for (i = start; i < end; i++)
        newIndices[newSize++] = indices_[i];
      start = startNegative_[iColumn];
      end = startPositive_[iColumn + 1];
      newNegative[newNumber++] = newSize;
      for (i = start; i < end; i++)
        newIndices[newSize++] = indices_[i];
    }
  }
  newPositive[newNumber] = newSize;
  delete[] which;
  delete[] startPositive_;
  startPositive_ = newPositive;
  delete[] startNegative_;
  startNegative_ = newNegative;
  delete[] indices_;
  indices_ = newIndices;
  numberColumns_ = newNumber;
}
/* Delete the rows whose indices are listed in <code>indDel</code>. */
void ClpPlusMinusOneMatrix::deleteRows(const int numDel, const int *indDel)
{
  int iRow;
  int numberBad = 0;
  // Use array to make sure we can have duplicates
  int *which = new int[numberRows_];
  memset(which, 0, numberRows_ * sizeof(int));
  int nDuplicate = 0;
  for (iRow = 0; iRow < numDel; iRow++) {
    int jRow = indDel[iRow];
    if (jRow < 0 || jRow >= numberRows_) {
      numberBad++;
    } else {
      if (which[jRow])
        nDuplicate++;
      else
        which[jRow] = 1;
    }
  }
  if (numberBad)
    throw CoinError("Indices out of range", "deleteRows", "ClpPlusMinusOneMatrix");
  CoinBigIndex iElement;
  CoinBigIndex numberElements = startPositive_[numberColumns_];
  CoinBigIndex newSize = 0;
  for (iElement = 0; iElement < numberElements; iElement++) {
    iRow = indices_[iElement];
    if (!which[iRow])
      newSize++;
  }
  int newNumber = numberRows_ - numDel + nDuplicate;
  // Get rid of temporary arrays
  delete[] lengths_;
  lengths_ = NULL;
  delete matrix_;
  matrix_ = NULL;
  // redo which
  int numberRows = 0;
  for (iRow = 0; iRow < numberRows_; iRow++) {
    if (which[iRow]) {
      which[iRow] = -1; // not wanted
    } else {
      which[iRow] = numberRows;
      numberRows++;
    }
  }
  int *newIndices = new int[newSize];
  newSize = 0;
  int iColumn;
  for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
    CoinBigIndex start, end;
    CoinBigIndex i;
    start = startPositive_[iColumn];
    end = startNegative_[iColumn];
    startPositive_[newNumber] = newSize;
    for (i = start; i < end; i++) {
      iRow = which[indices_[i]];
      if (iRow >= 0)
        newIndices[newSize++] = iRow;
    }
    start = startNegative_[iColumn];
    end = startPositive_[iColumn + 1];
    startNegative_[newNumber] = newSize;
    for (i = start; i < end; i++) {
      iRow = which[indices_[i]];
      if (iRow >= 0)
        newIndices[newSize++] = iRow;
    }
  }
  startPositive_[numberColumns_] = newSize;
  delete[] which;
  delete[] indices_;
  indices_ = newIndices;
  numberRows_ = newNumber;
}
bool ClpPlusMinusOneMatrix::isColOrdered() const
{
  return columnOrdered_;
}
/* Number of entries in the packed matrix. */
CoinBigIndex
ClpPlusMinusOneMatrix::getNumElements() const
{
  int numberMajor = (columnOrdered_) ? numberColumns_ : numberRows_;
  if (startPositive_)
    return startPositive_[numberMajor];
  else
    return 0;
}
// pass in copy (object takes ownership)
void ClpPlusMinusOneMatrix::passInCopy(int numberRows, int numberColumns,
  bool columnOrdered, int *indices,
  CoinBigIndex *startPositive, CoinBigIndex *startNegative)
{
  columnOrdered_ = columnOrdered;
#ifdef CLP_PLUS_ONE_MATRIX
  otherFlags_ = 0;
#endif
  startPositive_ = startPositive;
  startNegative_ = startNegative;
  indices_ = indices;
  numberRows_ = numberRows;
  numberColumns_ = numberColumns;
  // Check valid
  checkValid(false);
}
// Just checks matrix valid - will say if dimensions not quite right if detail
void ClpPlusMinusOneMatrix::checkValid(bool detail) const
{
  int maxIndex = -1;
  int minIndex = columnOrdered_ ? numberRows_ : numberColumns_;
  int number = !columnOrdered_ ? numberRows_ : numberColumns_;
  CoinBigIndex numberElements = getNumElements();
  CoinBigIndex last = -1;
  int bad = 0;
  // say if all +1
#ifdef CLP_PLUS_ONE_MATRIX
  otherFlags_ = 1;
#endif
  for (int i = 0; i < number; i++) {
    if (startPositive_[i] < last)
      bad++;
    else
      last = startPositive_[i];
#ifdef CLP_PLUS_ONE_MATRIX
    if (startNegative_[i] < startPositive_[i + 1])
      otherFlags_ &= ~1; // not all +1
#endif
    if (startNegative_[i] < last)
      bad++;
    else
      last = startNegative_[i];
  }
  if (startPositive_[number] < last)
    bad++;
  CoinAssertHint(!bad, "starts are not monotonic");
  (void)bad;
  for (CoinBigIndex cbi = 0; cbi < numberElements; cbi++) {
    maxIndex = std::max(indices_[cbi], maxIndex);
    minIndex = std::min(indices_[cbi], minIndex);
  }
  CoinAssert(maxIndex < (columnOrdered_ ? numberRows_ : numberColumns_));
  CoinAssert(minIndex >= 0);
  if (detail) {
    if (minIndex > 0 || maxIndex + 1 < (columnOrdered_ ? numberRows_ : numberColumns_))
      printf("Not full range of indices - %d to %d\n", minIndex, maxIndex);
  }
}
// Append Columns
void ClpPlusMinusOneMatrix::appendCols(int number, const CoinPackedVectorBase *const *columns)
{
  int iColumn;
  CoinBigIndex size = 0;
  int numberBad = 0;
  // say if all +1
#ifdef CLP_PLUS_ONE_MATRIX
  otherFlags_ = 0;
#endif
  for (iColumn = 0; iColumn < number; iColumn++) {
    int n = columns[iColumn]->getNumElements();
    const double *element = columns[iColumn]->getElements();
    size += n;
    int i;
    for (i = 0; i < n; i++) {
      if (fabs(element[i]) != 1.0)
        numberBad++;
    }
  }
  if (numberBad)
    throw CoinError("Not +- 1", "appendCols", "ClpPlusMinusOneMatrix");
  // Get rid of temporary arrays
  delete[] lengths_;
  lengths_ = NULL;
  delete matrix_;
  matrix_ = NULL;
  CoinBigIndex numberNow = startPositive_[numberColumns_];
  CoinBigIndex *temp;
  temp = new CoinBigIndex[numberColumns_ + 1 + number];
  CoinMemcpyN(startPositive_, (numberColumns_ + 1), temp);
  delete[] startPositive_;
  startPositive_ = temp;
  temp = new CoinBigIndex[numberColumns_ + number];
  CoinMemcpyN(startNegative_, numberColumns_, temp);
  delete[] startNegative_;
  startNegative_ = temp;
  int *temp2 = new int[numberNow + size];
  CoinMemcpyN(indices_, numberNow, temp2);
  delete[] indices_;
  indices_ = temp2;
  // now add
  size = numberNow;
  for (iColumn = 0; iColumn < number; iColumn++) {
    int n = columns[iColumn]->getNumElements();
    const int *row = columns[iColumn]->getIndices();
    const double *element = columns[iColumn]->getElements();
    int i;
    for (i = 0; i < n; i++) {
      if (element[i] == 1.0)
        indices_[size++] = row[i];
    }
    startNegative_[iColumn + numberColumns_] = size;
    for (i = 0; i < n; i++) {
      if (element[i] == -1.0)
        indices_[size++] = row[i];
    }
    startPositive_[iColumn + numberColumns_ + 1] = size;
  }

  numberColumns_ += number;
}
// Append Rows
void ClpPlusMinusOneMatrix::appendRows(int number, const CoinPackedVectorBase *const *rows)
{
  // Allocate arrays to use for counting
  int *countPositive = new int[numberColumns_ + 1];
  memset(countPositive, 0, numberColumns_ * sizeof(int));
  int *countNegative = new int[numberColumns_];
  memset(countNegative, 0, numberColumns_ * sizeof(int));
  int iRow;
  CoinBigIndex size = 0;
  int numberBad = 0;
  // say if all +1
#ifdef CLP_PLUS_ONE_MATRIX
  otherFlags_ = 0;
#endif
  for (iRow = 0; iRow < number; iRow++) {
    int n = rows[iRow]->getNumElements();
    const int *column = rows[iRow]->getIndices();
    const double *element = rows[iRow]->getElements();
    size += n;
    int i;
    for (i = 0; i < n; i++) {
      int iColumn = column[i];
      if (element[i] == 1.0)
        countPositive[iColumn]++;
      else if (element[i] == -1.0)
        countNegative[iColumn]++;
      else
        numberBad++;
    }
  }
  if (numberBad)
    throw CoinError("Not +- 1", "appendRows", "ClpPlusMinusOneMatrix");
  // Get rid of temporary arrays
  delete[] lengths_;
  lengths_ = NULL;
  delete matrix_;
  matrix_ = NULL;
  CoinBigIndex numberNow = startPositive_[numberColumns_];
  int *newIndices = new int[numberNow + size];
  // Update starts and turn counts into positions
  // also move current indices
  int iColumn;
  CoinBigIndex numberAdded = 0;
  for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
    int n, move;
    CoinBigIndex now;
    now = startPositive_[iColumn];
    move = static_cast< int >(startNegative_[iColumn] - now);
    n = countPositive[iColumn];
    startPositive_[iColumn] += numberAdded;
    CoinMemcpyN(indices_ + now, move, newIndices + startPositive_[iColumn]);
    countPositive[iColumn] = static_cast< int >(startNegative_[iColumn] + numberAdded);
    numberAdded += n;
    now = startNegative_[iColumn];
    move = static_cast< int >(startPositive_[iColumn + 1] - now);
    n = countNegative[iColumn];
    startNegative_[iColumn] += numberAdded;
    CoinMemcpyN(indices_ + now, move, newIndices + startNegative_[iColumn]);
    countNegative[iColumn] = static_cast< int >(startPositive_[iColumn + 1] + numberAdded);
    numberAdded += n;
  }
  delete[] indices_;
  indices_ = newIndices;
  startPositive_[numberColumns_] += numberAdded;
  // Now put in
  for (iRow = 0; iRow < number; iRow++) {
    int newRow = numberRows_ + iRow;
    int n = rows[iRow]->getNumElements();
    const int *column = rows[iRow]->getIndices();
    const double *element = rows[iRow]->getElements();
    int i;
    for (i = 0; i < n; i++) {
      int iColumn = column[i];
      int put;
      if (element[i] == 1.0) {
        put = countPositive[iColumn];
        countPositive[iColumn] = put + 1;
      } else {
        put = countNegative[iColumn];
        countNegative[iColumn] = put + 1;
      }
      indices_[put] = newRow;
    }
  }
  delete[] countPositive;
  delete[] countNegative;
  numberRows_ += number;
}
/* Returns largest and smallest elements of both signs.
   Largest refers to largest absolute value.
*/
void ClpPlusMinusOneMatrix::rangeOfElements(double &smallestNegative, double &largestNegative,
  double &smallestPositive, double &largestPositive)
{
  int iColumn;
  bool plusOne = false;
  bool minusOne = false;
  for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
    if (startNegative_[iColumn] > startPositive_[iColumn])
      plusOne = true;
    if (startPositive_[iColumn + 1] > startNegative_[iColumn])
      minusOne = true;
  }
  if (minusOne) {
    smallestNegative = -1.0;
    largestNegative = -1.0;
  } else {
    smallestNegative = 0.0;
    largestNegative = 0.0;
  }
  if (plusOne) {
    smallestPositive = 1.0;
    largestPositive = 1.0;
  } else {
    smallestPositive = 0.0;
    largestPositive = 0.0;
  }
}
// Says whether it can do partial pricing
bool ClpPlusMinusOneMatrix::canDoPartialPricing() const
{
  return true;
}
// Partial pricing
void ClpPlusMinusOneMatrix::partialPricing(ClpSimplex *model, double startFraction, double endFraction,
  int &bestSequence, int &numberWanted)
{
  numberWanted = currentWanted_;
  int start = static_cast< int >(startFraction * numberColumns_);
  int end = std::min(static_cast< int >(endFraction * numberColumns_ + 1), numberColumns_);
  CoinBigIndex j;
  double tolerance = model->currentDualTolerance();
  double *COIN_RESTRICT reducedCost = model->djRegion();
  const double *COIN_RESTRICT duals = model->dualRowSolution();
  const double *COIN_RESTRICT cost = model->costRegion();
  double bestDj;
  if (bestSequence >= 0)
    bestDj = fabs(reducedCost[bestSequence]);
  else
    bestDj = tolerance;
  int sequenceOut = model->sequenceOut();
  int saveSequence = bestSequence;
  int iSequence;
  for (iSequence = start; iSequence < end; iSequence++) {
    if (iSequence != sequenceOut) {
      double value;
      ClpSimplex::Status status = model->getStatus(iSequence);

      switch (status) {

      case ClpSimplex::basic:
      case ClpSimplex::isFixed:
        break;
      case ClpSimplex::isFree:
      case ClpSimplex::superBasic:
        value = cost[iSequence];
        j = startPositive_[iSequence];
        for (; j < startNegative_[iSequence]; j++) {
          int iRow = indices_[j];
          value -= duals[iRow];
        }
        for (; j < startPositive_[iSequence + 1]; j++) {
          int iRow = indices_[j];
          value += duals[iRow];
        }
        value = fabs(value);
        if (value > FREE_ACCEPT * tolerance) {
          numberWanted--;
          // we are going to bias towards free (but only if reasonable)
          value *= FREE_BIAS;
          if (value > bestDj) {
            // check flagged variable and correct dj
            if (!model->flagged(iSequence)) {
              bestDj = value;
              bestSequence = iSequence;
            } else {
              // just to make sure we don't exit before got something
              numberWanted++;
            }
          }
        }
        break;
      case ClpSimplex::atUpperBound:
        value = cost[iSequence];
        j = startPositive_[iSequence];
        for (; j < startNegative_[iSequence]; j++) {
          int iRow = indices_[j];
          value -= duals[iRow];
        }
        for (; j < startPositive_[iSequence + 1]; j++) {
          int iRow = indices_[j];
          value += duals[iRow];
        }
        if (value > tolerance) {
          numberWanted--;
          if (value > bestDj) {
            // check flagged variable and correct dj
            if (!model->flagged(iSequence)) {
              bestDj = value;
              bestSequence = iSequence;
            } else {
              // just to make sure we don't exit before got something
              numberWanted++;
            }
          }
        }
        break;
      case ClpSimplex::atLowerBound:
        value = cost[iSequence];
        j = startPositive_[iSequence];
        for (; j < startNegative_[iSequence]; j++) {
          int iRow = indices_[j];
          value -= duals[iRow];
        }
        for (; j < startPositive_[iSequence + 1]; j++) {
          int iRow = indices_[j];
          value += duals[iRow];
        }
        value = -value;
        if (value > tolerance) {
          numberWanted--;
          if (value > bestDj) {
            // check flagged variable and correct dj
            if (!model->flagged(iSequence)) {
              bestDj = value;
              bestSequence = iSequence;
            } else {
              // just to make sure we don't exit before got something
              numberWanted++;
            }
          }
        }
        break;
      }
    }
    if (!numberWanted)
      break;
  }
  if (bestSequence != saveSequence) {
    // recompute dj
    double value = cost[bestSequence];
    j = startPositive_[bestSequence];
    for (; j < startNegative_[bestSequence]; j++) {
      int iRow = indices_[j];
      value -= duals[iRow];
    }
    for (; j < startPositive_[bestSequence + 1]; j++) {
      int iRow = indices_[j];
      value += duals[iRow];
    }
    reducedCost[bestSequence] = value;
    savedBestSequence_ = bestSequence;
    savedBestDj_ = reducedCost[savedBestSequence_];
  }
  currentWanted_ = numberWanted;
}
// Allow any parts of a created CoinMatrix to be deleted
void ClpPlusMinusOneMatrix::releasePackedMatrix() const
{
  delete matrix_;
  delete[] lengths_;
  matrix_ = NULL;
  lengths_ = NULL;
  // say if all +1
#ifdef CLP_PLUS_ONE_MATRIX
  otherFlags_ = 0;
#endif
}
/* Returns true if can combine transposeTimes and subsetTransposeTimes
   and if it would be faster */
bool ClpPlusMinusOneMatrix::canCombine(const ClpSimplex *model,
  const CoinIndexedVector *pi) const
{
  int numberInRowArray = pi->getNumElements();
  int numberRows = model->numberRows();
  bool packed = pi->packedMode();
  // factor should be smaller if doing both with two pi vectors
  double factor = 0.27;
  // We may not want to do by row if there may be cache problems
  // It would be nice to find L2 cache size - for moment 512K
  // Be slightly optimistic
  if (numberColumns_ * sizeof(double) > 1000000) {
    if (numberRows * 10 < numberColumns_)
      factor *= 0.333333333;
    else if (numberRows * 4 < numberColumns_)
      factor *= 0.5;
    else if (numberRows * 2 < numberColumns_)
      factor *= 0.66666666667;
    //if (model->numberIterations()%50==0)
    //printf("%d nonzero\n",numberInRowArray);
  }
  // if not packed then bias a bit more towards by column
  if (!packed)
    factor *= 0.9;
  return (numberInRowArray > factor * numberRows || !model->rowCopy());
}
// These have to match ClpPrimalColumnSteepest version
#define reference(i) (((reference[i >> 5] >> (i & 31)) & 1) != 0)
/* Updates two arrays for steepest and does devex weights
   Returns nonzero if updates reduced cost and infeas -
   new infeas in dj1  - This does not at present*/
int ClpPlusMinusOneMatrix::transposeTimes2(const ClpSimplex *model,
  const CoinIndexedVector *pi1, CoinIndexedVector *dj1,
  const CoinIndexedVector *pi2,
  CoinIndexedVector *spare,
  double *, double *,
  double referenceIn, double devex,
  // Array for exact devex to say what is in reference framework
  unsigned int *COIN_RESTRICT reference,
  double *COIN_RESTRICT weights, double scaleFactor)
{
  // put row of tableau in dj1
  double *COIN_RESTRICT pi = pi1->denseVector();
  int numberNonZero = 0;
  int *COIN_RESTRICT index = dj1->getIndices();
  double *COIN_RESTRICT array = dj1->denseVector();
  int numberInRowArray = pi1->getNumElements();
  double zeroTolerance = model->zeroTolerance();
  bool packed = pi1->packedMode();
  // do by column
  int iColumn;
  assert(!spare->getNumElements());
  double *COIN_RESTRICT piWeight = pi2->denseVector();
  assert(!pi2->packedMode());
  bool killDjs = (scaleFactor == 0.0);
  if (!scaleFactor)
    scaleFactor = 1.0;
  // Note scale factor was -1.0
  if (packed) {
    // need to expand pi into y
    assert(spare->capacity() >= model->numberRows());
    double *COIN_RESTRICT piOld = pi;
    pi = spare->denseVector();
    const int *COIN_RESTRICT whichRow = pi1->getIndices();
    int i;
    // modify pi so can collapse to one loop
    for (i = 0; i < numberInRowArray; i++) {
      int iRow = whichRow[i];
      pi[iRow] = piOld[i];
    }
    CoinBigIndex j;
    for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
      ClpSimplex::Status status = model->getStatus(iColumn);
      if (status == ClpSimplex::basic || status == ClpSimplex::isFixed)
        continue;
      double value = 0.0;
      for (j = startPositive_[iColumn]; j < startNegative_[iColumn]; j++) {
        int iRow = indices_[j];
        value -= pi[iRow];
      }
      for (; j < startPositive_[iColumn + 1]; j++) {
        int iRow = indices_[j];
        value += pi[iRow];
      }
      if (fabs(value) > zeroTolerance) {
        // and do other array
        double modification = 0.0;
        for (j = startPositive_[iColumn]; j < startNegative_[iColumn]; j++) {
          int iRow = indices_[j];
          modification += piWeight[iRow];
        }
        for (; j < startPositive_[iColumn + 1]; j++) {
          int iRow = indices_[j];
          modification -= piWeight[iRow];
        }
        double thisWeight = weights[iColumn];
        double pivot = value * scaleFactor;
        double pivotSquared = pivot * pivot;
        thisWeight += pivotSquared * devex + pivot * modification;
        if (thisWeight < DEVEX_TRY_NORM) {
          if (referenceIn < 0.0) {
            // steepest
            thisWeight = std::max(DEVEX_TRY_NORM, DEVEX_ADD_ONE + pivotSquared);
          } else {
            // exact
            thisWeight = referenceIn * pivotSquared;
            if (reference(iColumn))
              thisWeight += 1.0;
            thisWeight = std::max(thisWeight, DEVEX_TRY_NORM);
          }
        }
        weights[iColumn] = thisWeight;
        if (!killDjs) {
          array[numberNonZero] = value;
          index[numberNonZero++] = iColumn;
        }
      }
    }
    // zero out
    for (i = 0; i < numberInRowArray; i++) {
      int iRow = whichRow[i];
      pi[iRow] = 0.0;
    }
  } else {
    CoinBigIndex j;
    for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
      ClpSimplex::Status status = model->getStatus(iColumn);
      if (status == ClpSimplex::basic || status == ClpSimplex::isFixed)
        continue;
      double value = 0.0;
      for (j = startPositive_[iColumn]; j < startNegative_[iColumn]; j++) {
        int iRow = indices_[j];
        value -= pi[iRow];
      }
      for (; j < startPositive_[iColumn + 1]; j++) {
        int iRow = indices_[j];
        value += pi[iRow];
      }
      if (fabs(value) > zeroTolerance) {
        // and do other array
        double modification = 0.0;
        for (j = startPositive_[iColumn]; j < startNegative_[iColumn]; j++) {
          int iRow = indices_[j];
          modification += piWeight[iRow];
        }
        for (; j < startPositive_[iColumn + 1]; j++) {
          int iRow = indices_[j];
          modification -= piWeight[iRow];
        }
        double thisWeight = weights[iColumn];
        double pivot = value * scaleFactor;
        double pivotSquared = pivot * pivot;
        thisWeight += pivotSquared * devex + pivot * modification;
        if (thisWeight < DEVEX_TRY_NORM) {
          if (referenceIn < 0.0) {
            // steepest
            thisWeight = std::max(DEVEX_TRY_NORM, DEVEX_ADD_ONE + pivotSquared);
          } else {
            // exact
            thisWeight = referenceIn * pivotSquared;
            if (reference(iColumn))
              thisWeight += 1.0;
            thisWeight = std::max(thisWeight, DEVEX_TRY_NORM);
          }
        }
        weights[iColumn] = thisWeight;
        if (!killDjs) {
          array[iColumn] = value;
          index[numberNonZero++] = iColumn;
        }
      }
    }
  }
  dj1->setNumElements(numberNonZero);
  spare->setNumElements(0);
  if (packed)
    dj1->setPackedMode(true);
  return 0;
}
// Updates second array for steepest and does devex weights
void ClpPlusMinusOneMatrix::subsetTimes2(const ClpSimplex *,
  CoinIndexedVector *dj1,
  const CoinIndexedVector *pi2, CoinIndexedVector *,
  double referenceIn, double devex,
  // Array for exact devex to say what is in reference framework
  unsigned int *COIN_RESTRICT reference,
  double *COIN_RESTRICT weights, double scaleFactor)
{
  int number = dj1->getNumElements();
  const int *COIN_RESTRICT index = dj1->getIndices();
  double *COIN_RESTRICT array = dj1->denseVector();
  assert(dj1->packedMode());

  double *COIN_RESTRICT piWeight = pi2->denseVector();
  bool killDjs = (scaleFactor == 0.0);
  if (!scaleFactor)
    scaleFactor = 1.0;
  for (int k = 0; k < number; k++) {
    int iColumn = index[k];
    double pivot = array[k] * scaleFactor;
    if (killDjs)
      array[k] = 0.0;
    // and do other array
    double modification = 0.0;
    CoinBigIndex j;
    for (j = startPositive_[iColumn]; j < startNegative_[iColumn]; j++) {
      int iRow = indices_[j];
      modification += piWeight[iRow];
    }
    for (; j < startPositive_[iColumn + 1]; j++) {
      int iRow = indices_[j];
      modification -= piWeight[iRow];
    }
    double thisWeight = weights[iColumn];
    double pivotSquared = pivot * pivot;
    thisWeight += pivotSquared * devex + pivot * modification;
    if (thisWeight < DEVEX_TRY_NORM) {
      if (referenceIn < 0.0) {
        // steepest
        thisWeight = std::max(DEVEX_TRY_NORM, DEVEX_ADD_ONE + pivotSquared);
      } else {
        // exact
        thisWeight = referenceIn * pivotSquared;
        if (reference(iColumn))
          thisWeight += 1.0;
        thisWeight = std::max(thisWeight, DEVEX_TRY_NORM);
      }
    }
    weights[iColumn] = thisWeight;
  }
}
/* Set the dimensions of the matrix. In effect, append new empty
   columns/rows to the matrix. A negative number for either dimension
   means that that dimension doesn't change. Otherwise the new dimensions
   MUST be at least as large as the current ones otherwise an exception
   is thrown. */
void ClpPlusMinusOneMatrix::setDimensions(int newnumrows, int newnumcols)
{
  if (newnumrows < 0)
    newnumrows = numberRows_;
  if (newnumrows < numberRows_)
    throw CoinError("Bad new rownum (less than current)",
      "setDimensions", "CoinPackedMatrix");

  if (newnumcols < 0)
    newnumcols = numberColumns_;
  if (newnumcols < numberColumns_)
    throw CoinError("Bad new colnum (less than current)",
      "setDimensions", "CoinPackedMatrix");

  int number = 0;
  int length = 0;
  if (columnOrdered_) {
    length = numberColumns_;
    numberColumns_ = newnumcols;
    number = numberColumns_;

  } else {
    length = numberRows_;
    numberRows_ = newnumrows;
    number = numberRows_;
  }
  if (number > length) {
    CoinBigIndex *temp;
    int i;
    CoinBigIndex end = startPositive_[length];
    temp = new CoinBigIndex[number + 1];
    CoinMemcpyN(startPositive_, (length + 1), temp);
    delete[] startPositive_;
    for (i = length + 1; i < number + 1; i++)
      temp[i] = end;
    startPositive_ = temp;
    temp = new CoinBigIndex[number];
    CoinMemcpyN(startNegative_, length, temp);
    delete[] startNegative_;
    for (i = length; i < number; i++)
      temp[i] = end;
    startNegative_ = temp;
  }
}
#ifndef SLIM_CLP
/* Append a set of rows/columns to the end of the matrix. Returns number of errors
   i.e. if any of the new rows/columns contain an index that's larger than the
   number of columns-1/rows-1 (if numberOther>0) or duplicates
   If 0 then rows, 1 if columns */
int ClpPlusMinusOneMatrix::appendMatrix(int number, int type,
  const CoinBigIndex *starts, const int *index,
  const double *element, int /*numberOther*/)
{
  int numberErrors = 0;
  // say if all +1
#ifdef CLP_PLUS_ONE_MATRIX
  otherFlags_ = 0;
#endif
  // make into CoinPackedVector
  CoinPackedVectorBase **vectors = new CoinPackedVectorBase *[number];
  int iVector;
  for (iVector = 0; iVector < number; iVector++) {
    CoinBigIndex iStart = starts[iVector];
    vectors[iVector] = new CoinPackedVector(static_cast< int >(starts[iVector + 1] - iStart),
      index + iStart, element + iStart);
  }
  if (type == 0) {
    // rows
    appendRows(number, vectors);
  } else {
    // columns
    appendCols(number, vectors);
  }
  for (iVector = 0; iVector < number; iVector++)
    delete vectors[iVector];
  delete[] vectors;
  return numberErrors;
}
#endif
#if CLP_POOL_MATRIX
//#############################################################################
// Constructors / Destructor / Assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
ClpPoolMatrix::ClpPoolMatrix()
  : ClpMatrixBase()
{
  setType(13);
  matrix_ = NULL;
  lengths_ = NULL;
  elements_ = NULL;
  columnStart_ = NULL;
  stuff_ = NULL;
  numberRows_ = 0;
  numberColumns_ = 0;
  numberDifferent_ = 0;
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
ClpPoolMatrix::ClpPoolMatrix(const ClpPoolMatrix &rhs)
  : ClpMatrixBase(rhs)
{
  setType(13);
  matrix_ = NULL;
  lengths_ = NULL;
  elements_ = NULL;
  columnStart_ = NULL;
  stuff_ = NULL;
  numberRows_ = rhs.numberRows_;
  numberColumns_ = rhs.numberColumns_;
  numberDifferent_ = rhs.numberDifferent_;
  if (numberColumns_) {
    columnStart_ = CoinCopyOfArray(rhs.columnStart_, numberColumns_ + 1);
    CoinBigIndex numberElements = columnStart_[numberColumns_];
    stuff_ = CoinCopyOfArray(rhs.stuff_, numberElements);
    elements_ = CoinCopyOfArray(rhs.elements_, numberDifferent_);
  }
}
static const int mmult[] = {
  262139, 259459, 256889, 254291, 251701, 249133, 246709, 244247
};
inline int hashit(double value)
{
  const unsigned char *chars = reinterpret_cast< unsigned char * >(&value);
  int n = 0;
  for (int j = 0; j < 8; j++)
    n += mmult[j] * chars[j];
  n = n % (1 << CLP_POOL_SIZE);
  return n;
}
ClpPoolMatrix::ClpPoolMatrix(const CoinPackedMatrix &rhs)
  : ClpMatrixBase()
{
  setType(13);
  matrix_ = NULL;
  lengths_ = NULL;
  numberDifferent_ = 0;
  assert(rhs.isColOrdered());
  // get matrix data pointers
  const int *row = rhs.getIndices();
  const CoinBigIndex *columnStart = rhs.getVectorStarts();
  const int *columnLength = rhs.getVectorLengths();
  const double *elementByColumn = rhs.getElements();
  numberColumns_ = rhs.getNumCols();
  numberRows_ = rhs.getNumRows();
  assert(numberRows_ < (1 << CLP_POOL_MATRIX));
  elements_ = NULL;
  columnStart_ = new CoinBigIndex[numberColumns_ + 1];
  stuff_ = new poolInfo[rhs.getNumElements()];
  int maxPool = 1 << CLP_POOL_SIZE;
  double *tempDifferent = new double[maxPool];
#define HASH 2
#if HASH > 1
  // for hashing
  typedef struct {
    int index, next;
  } CoinHashLink;
  CoinHashLink *hashThis = new CoinHashLink[4 * maxPool];
  for (int i = 0; i < 4 * maxPool; i++) {
    hashThis[i].index = -1;
    hashThis[i].next = -1;
  }
#endif
  int numberElements = 0;
  int hashDifferent = 0;
  for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
    CoinBigIndex k;
    columnStart_[iColumn] = numberElements;
    for (k = columnStart[iColumn]; k < columnStart[iColumn] + columnLength[iColumn];
         k++) {
      int iRow = row[k];
      double value = elementByColumn[k];
      poolInfo stuff;
      stuff.row_ = iRow;
#if HASH == 1 || HASH == 3
      int j;
      for (j = 0; j < numberDifferent_; j++) {
        if (value == tempDifferent[j])
          break;
      }
#endif
#if HASH == 2 || HASH == 3
      int ipos = hashit(value);
      int j1;
      while (true) {
        j1 = hashThis[ipos].index;
        if (j1 == -1) {
          hashThis[ipos].index = numberDifferent_;
          j1 = numberDifferent_;
          break;
        } else if (value == tempDifferent[j1]) {
          break;
        } else {
          int k = hashThis[ipos].next;
          if (k == -1) {
            j1 = numberDifferent_;
            while (true) {
              ++hashDifferent;
              if (hashThis[hashDifferent].index == -1) {
                break;
              }
            }
            hashThis[ipos].next = hashDifferent;
            hashThis[hashDifferent].index = numberDifferent_;
            break;
          } else {
            ipos = k;
            /* nothing worked - try it again */
          }
        }
      }
#if HASH == 2
      int j = j1;
#endif
#endif
#if HASH == 3
      assert(j == j1);
#endif
      if (j == numberDifferent_) {
        if (j == maxPool)
          break; // too many
        tempDifferent[j] = value;
        numberDifferent_++;
      }
      stuff.pool_ = j;
      stuff_[numberElements++] = stuff;
    }
  }
  columnStart_[numberColumns_] = numberElements;
#if HASH > 1
  delete[] hashThis;
#endif
  elements_ = new double[numberDifferent_];
  memcpy(elements_, tempDifferent, numberDifferent_ * sizeof(double));
  delete[] tempDifferent;
  if (numberDifferent_ == maxPool) {
    delete[] stuff_;
    delete[] elements_;
    delete[] columnStart_;
    elements_ = NULL;
    columnStart_ = NULL;
    stuff_ = NULL;
    numberRows_ = -1;
    numberColumns_ = -1;
    numberDifferent_ = -1;
  }
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
ClpPoolMatrix::~ClpPoolMatrix()
{
  delete matrix_;
  delete[] elements_;
  delete[] columnStart_;
  delete[] lengths_;
  delete[] stuff_;
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
ClpPoolMatrix &
ClpPoolMatrix::operator=(const ClpPoolMatrix &rhs)
{
  if (this != &rhs) {
    ClpMatrixBase::operator=(rhs);
    delete matrix_;
    delete[] elements_;
    delete[] columnStart_;
    delete[] lengths_;
    delete[] stuff_;
    matrix_ = NULL;
    lengths_ = NULL;
    numberRows_ = rhs.numberRows_;
    numberColumns_ = rhs.numberColumns_;
    if (numberColumns_) {
      columnStart_ = CoinCopyOfArray(rhs.columnStart_, numberColumns_ + 1);
      CoinBigIndex numberElements = columnStart_[numberColumns_];
      stuff_ = CoinCopyOfArray(rhs.stuff_, numberElements);
      elements_ = CoinCopyOfArray(rhs.elements_, numberDifferent_);
    }
  }
  return *this;
}
// Constructor from arrays - handing over ownership
ClpPoolMatrix::ClpPoolMatrix(int numberColumns, CoinBigIndex *columnStart,
  poolInfo *stuff, double *elements)
  : ClpMatrixBase()
{
  setType(13);
  matrix_ = NULL;
  lengths_ = NULL;
  numberColumns_ = numberColumns;
  numberDifferent_ = 0;
  numberRows_ = 0;
  columnStart_ = columnStart;
  stuff_ = stuff;
  elements_ = elements;
  for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
    CoinBigIndex k;
    for (k = columnStart_[iColumn]; k < columnStart_[iColumn + 1];
         k++) {
      int iRow = stuff_[k].row_;
      int iPool = stuff_[k].pool_;
      numberDifferent_ = std::max(numberDifferent_, iPool);
      numberRows_ = std::max(numberRows_, iRow);
    }
  }
  // adjust
  numberDifferent_++;
  numberRows_++;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
ClpMatrixBase *ClpPoolMatrix::clone() const
{
  return new ClpPoolMatrix(*this);
}
/* Subset clone (without gaps).  Duplicates are allowed
   and order is as given */
ClpMatrixBase *
ClpPoolMatrix::subsetClone(int numberRows, const int *whichRows,
  int numberColumns,
  const int *whichColumns) const
{
  return new ClpPoolMatrix(*this, numberRows, whichRows,
    numberColumns, whichColumns);
}
/* Subset constructor (without gaps).  Duplicates are allowed
   and order is as given */
ClpPoolMatrix::ClpPoolMatrix(
  const ClpPoolMatrix &rhs,
  int numberRows, const int *whichRow,
  int numberColumns, const int *whichColumn)
  : ClpMatrixBase(rhs)
{
  setType(13);
  matrix_ = NULL;
  lengths_ = NULL;
  elements_ = NULL;
  columnStart_ = NULL;
  stuff_ = NULL;
  numberRows_ = numberRows;
  numberColumns_ = numberColumns;
  numberDifferent_ = 0;
  if (!numberRows_ || !numberColumns_)
    return;
  numberDifferent_ = rhs.numberDifferent_;
  elements_ = CoinCopyOfArray(rhs.elements_, numberDifferent_);
  int numberRowsOther = rhs.numberRows_;
  int numberColumnsOther = rhs.numberColumns_;
  int *newRow = new int[numberRowsOther + numberRows_];
  int *duplicateRow = newRow + numberRowsOther;
  for (int iRow = 0; iRow < numberRowsOther; iRow++)
    newRow[iRow] = -1;
  // and array for duplicating rows
  int numberBad = 0;
  for (int iRow = 0; iRow < numberRows_; iRow++) {
    duplicateRow[iRow] = -1;
    int kRow = whichRow[iRow];
    if (kRow >= 0 && kRow < numberRowsOther) {
      if (newRow[kRow] < 0) {
        // first time
        newRow[kRow] = iRow;
      } else {
        // duplicate
        int lastRow = newRow[kRow];
        newRow[kRow] = iRow;
        duplicateRow[iRow] = lastRow;
      }
    } else {
      // bad row
      numberBad++;
    }
  }

  if (numberBad)
    throw CoinError("bad row entries",
      "subset constructor", "ClpPoolMatrix");
  // now get size and check columns
  CoinBigIndex size = 0;
  numberBad = 0;
  const CoinBigIndex *starts = rhs.columnStart_;
  const poolInfo *stuff = rhs.stuff_;
  for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
    int kColumn = whichColumn[iColumn];
    if (kColumn >= 0 && kColumn < numberColumnsOther) {
      CoinBigIndex i;
      for (i = starts[kColumn]; i < starts[kColumn + 1]; i++) {
        int kRow = stuff[i].row_;
        kRow = newRow[kRow];
        while (kRow >= 0) {
          size++;
          kRow = duplicateRow[kRow];
        }
      }
    } else {
      // bad column
      numberBad++;
    }
  }
  if (numberBad)
    throw CoinError("bad major entries",
      "subset constructor", "ClpPoolMatrix");
  // now create arrays
  stuff_ = new poolInfo[size];
  columnStart_ = new CoinBigIndex[numberColumns_ + 1];
  size = 0;
  for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
    int kColumn = whichColumn[iColumn];
    columnStart_[iColumn] = size;
    CoinBigIndex i;
    for (i = starts[kColumn]; i < starts[kColumn + 1]; i++) {
      int kRow = stuff[i].row_;
      int iPool = stuff[i].pool_;
      kRow = newRow[kRow];
      while (kRow >= 0) {
        stuff_[size].row_ = kRow;
        stuff_[size++].pool_ = iPool;
        kRow = duplicateRow[kRow];
      }
    }
  }
  columnStart_[numberColumns_] = size;
  delete[] newRow;
}
// Create matrix_
ClpPackedMatrix *
ClpPoolMatrix::createMatrix() const
{
  if (!matrix_) {
    CoinBigIndex numberElements = columnStart_[numberColumns_];
    double *elements = new double[numberElements];
    int *rows = new int[numberElements];
    bool lengthsExist = lengths_ != NULL;
    numberElements = 0;
    if (!lengthsExist)
      lengths_ = new int[numberColumns_];
    for (int i = 0; i < numberColumns_; i++) {
      lengths_[i] = static_cast< int >(columnStart_[i + 1] - columnStart_[i]);
      for (CoinBigIndex j = columnStart_[i]; j < columnStart_[i + 1]; j++) {
        elements[numberElements] = elements_[stuff_[j].pool_];
        rows[numberElements++] = stuff_[j].row_;
      }
    }
    CoinPackedMatrix *matrix = new CoinPackedMatrix(true, numberRows_, numberColumns_,
      numberElements,
      elements, rows,
      columnStart_, lengths_);
    delete[] elements;
    delete[] rows;
    matrix_ = new ClpPackedMatrix(matrix);
  }
  return matrix_;
}

/* Returns a new matrix in reverse order without gaps */
ClpMatrixBase *
ClpPoolMatrix::reverseOrderedCopy() const
{
  // define to say no row copy
#define BIG_MATRIX
#ifndef BIG_MATRIX
  printf("XcreateMatrix at file %s line %d\n", __FILE__, __LINE__);
  createMatrix();
  return matrix_->reverseOrderedCopy();
#else
  return NULL;
#endif
}
//unscaled versions
void ClpPoolMatrix::times(double scalar,
  const double *COIN_RESTRICT x, double *COIN_RESTRICT y) const
{
#if 0
  printf("CcreateMatrix at file %s line %d\n",__FILE__,__LINE__);
  createMatrix()->times(scalar,x,y);
#else
  for (int iColumn = 0; iColumn < numberColumns_; iColumn++) {
    CoinBigIndex j;
    double value = x[iColumn];
    if (value) {
      CoinBigIndex start = columnStart_[iColumn];
      CoinBigIndex end = columnStart_[iColumn + 1];
      value *= scalar;
      for (j = start; j < end; j++) {
        int iRow = stuff_[j].row_;
        y[iRow] += value * elements_[stuff_[j].pool_];
      }
    }
  }
#endif
}
void ClpPoolMatrix::transposeTimes(double scalar,
  const double *COIN_RESTRICT x, double *COIN_RESTRICT y) const
{
#if 0
  printf("CcreateMatrix at file %s line %d\n",__FILE__,__LINE__);
  createMatrix()->transposeTimes(scalar,x,y);
#else
  int iColumn;
  CoinBigIndex start = columnStart_[0];
  if (scalar == -1.0) {
    for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
      CoinBigIndex j;
      CoinBigIndex next = columnStart_[iColumn + 1];
      double value = 0.0;
      // scaled
      for (j = start; j < next; j++) {
        int jRow = stuff_[j].row_;
        value += x[jRow] * elements_[stuff_[j].pool_];
      }
      start = next;
      y[iColumn] -= value;
    }
  } else {
    for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
      CoinBigIndex j;
      CoinBigIndex next = columnStart_[iColumn + 1];
      double value = 0.0;
      // scaled
      for (j = start; j < next; j++) {
        int jRow = stuff_[j].row_;
        value += x[jRow] * elements_[stuff_[j].pool_];
      }
      start = next;
      y[iColumn] += value * scalar;
    }
  }
#endif
}
void ClpPoolMatrix::times(double scalar,
  const double *COIN_RESTRICT x, double *COIN_RESTRICT y,
  const double *rowScale,
  const double *columnScale) const
{
  printf("createMatrix at file %s line %d\n", __FILE__, __LINE__);
  createMatrix()->times(scalar, x, y);
}
void ClpPoolMatrix::transposeTimes(double scalar,
  const double *COIN_RESTRICT x, double *COIN_RESTRICT y,
  const double *rowScale,
  const double *columnScale,
  double *spare) const
{
#if 0
  printf("CcreateMatrix at file %s line %d\n",__FILE__,__LINE__);
  createMatrix()->transposeTimes(scalar,x,y,rowScale,columnScale,spare);
#else
  if (rowScale) {
    int iColumn;
    if (!spare) {
      CoinBigIndex start = columnStart_[0];
      if (scalar == -1.0) {
        for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
          CoinBigIndex j;
          CoinBigIndex next = columnStart_[iColumn + 1];
          double value = 0.0;
          // scaled
          for (j = start; j < next; j++) {
            int jRow = stuff_[j].row_;
            value += x[jRow] * elements_[stuff_[j].pool_] * rowScale[jRow];
          }
          start = next;
          y[iColumn] -= value * columnScale[iColumn];
        }
      } else {
        for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
          CoinBigIndex j;
          CoinBigIndex next = columnStart_[iColumn + 1];
          double value = 0.0;
          // scaled
          for (j = start; j < next; j++) {
            int jRow = stuff_[j].row_;
            value += x[jRow] * elements_[stuff_[j].pool_] * rowScale[jRow];
          }
          start = next;
          y[iColumn] += value * scalar * columnScale[iColumn];
        }
      }
    } else {
      // can use spare region
      int iRow;
      int numberRows = matrix_->getNumRows();
      for (iRow = 0; iRow < numberRows; iRow++) {
        double value = x[iRow];
        if (value)
          spare[iRow] = value * rowScale[iRow];
        else
          spare[iRow] = 0.0;
      }
      CoinBigIndex start = columnStart_[0];
      for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
        CoinBigIndex j;
        CoinBigIndex next = columnStart_[iColumn + 1];
        double value = 0.0;
        // scaled
        for (j = start; j < next; j++) {
          int jRow = stuff_[j].row_;
          value += spare[jRow] * elements_[stuff_[j].pool_];
        }
        start = next;
        y[iColumn] += value * scalar * columnScale[iColumn];
      }
    }
  } else {
    transposeTimes(scalar, x, y);
  }
#endif
}
/* Return <code>x * A + y</code> in <code>z</code>.
   Squashes small elements and knows about ClpSimplex */
void ClpPoolMatrix::transposeTimes(const ClpSimplex *model, double scalar,
  const CoinIndexedVector *rowArray,
  CoinIndexedVector *y,
  CoinIndexedVector *columnArray) const
{
#if 0
  bool packed = rowArray->packedMode();
  assert (packed);
  printf("CcreateMatrix at file %s line %d\n",__FILE__,__LINE__);
  createMatrix()->transposeTimes(model,scalar,rowArray,y,columnArray);
#else
  columnArray->clear();
  // do by column
  double *COIN_RESTRICT pi = rowArray->denseVector();
  int numberNonZero = 0;
  int *COIN_RESTRICT index = columnArray->getIndices();
  double *COIN_RESTRICT array = columnArray->denseVector();
  int numberInRowArray = rowArray->getNumElements();
  // maybe I need one in OsiSimplex
  double zeroTolerance = model->zeroTolerance();
  bool packed = rowArray->packedMode();
  assert(packed);
  // do by column
  const double *COIN_RESTRICT rowScale = model->rowScale();
  assert(!y->getNumElements());
  // need to expand pi into y
  assert(y->capacity() >= model->numberRows());
  double *COIN_RESTRICT piOld = pi;
  pi = y->denseVector();
  const int *COIN_RESTRICT whichRow = rowArray->getIndices();
  int i;
  // later combine phases
  if (!rowScale) {
    // modify pi so can collapse to one loop
    for (i = 0; i < numberInRowArray; i++) {
      int iRow = whichRow[i];
      pi[iRow] = scalar * piOld[i];
    }
    double value = 0.0;
    CoinBigIndex j;
    CoinBigIndex end = columnStart_[1];
    for (j = columnStart_[0]; j < end; j++) {
      int iRow = stuff_[j].row_;
      value += pi[iRow] * elements_[stuff_[j].pool_];
    }
    int iColumn;
    for (iColumn = 0; iColumn < numberColumns_ - 1; iColumn++) {
      CoinBigIndex start = end;
      end = columnStart_[iColumn + 2];
      if (fabs(value) > zeroTolerance) {
        array[numberNonZero] = value;
        index[numberNonZero++] = iColumn;
      }
      value = 0.0;
      for (j = start; j < end; j++) {
        int iRow = stuff_[j].row_;
        value += pi[iRow] * elements_[stuff_[j].pool_];
      }
    }
    if (fabs(value) > zeroTolerance) {
      array[numberNonZero] = value;
      index[numberNonZero++] = iColumn;
    }
  } else {
    // scaled
    // modify pi so can collapse to one loop
    for (i = 0; i < numberInRowArray; i++) {
      int iRow = whichRow[i];
      pi[iRow] = scalar * piOld[i] * rowScale[iRow];
    }
    const double *columnScale = model->columnScale();
    double value = 0.0;
    double scale = columnScale[0];
    CoinBigIndex j;
    CoinBigIndex end = columnStart_[1];
    for (j = columnStart_[0]; j < end; j++) {
      int iRow = stuff_[j].row_;
      value += pi[iRow] * elements_[stuff_[j].pool_];
    }
    int iColumn;
    for (iColumn = 0; iColumn < numberColumns_ - 1; iColumn++) {
      value *= scale;
      CoinBigIndex start = end;
      scale = columnScale[iColumn + 1];
      end = columnStart_[iColumn + 2];
      if (fabs(value) > zeroTolerance) {
        array[numberNonZero] = value;
        index[numberNonZero++] = iColumn;
      }
      value = 0.0;
      for (j = start; j < end; j++) {
        int iRow = stuff_[j].row_;
        value += pi[iRow] * elements_[stuff_[j].pool_];
      }
    }
    value *= scale;
    if (fabs(value) > zeroTolerance) {
      array[numberNonZero] = value;
      index[numberNonZero++] = iColumn;
    }
  }
  // zero out
  int numberRows = model->numberRows();
  if (numberInRowArray * 4 < numberRows) {
    for (i = 0; i < numberInRowArray; i++) {
      int iRow = whichRow[i];
      pi[iRow] = 0.0;
    }
  } else {
    CoinZeroN(pi, numberRows);
  }
  columnArray->setNumElements(numberNonZero);
  y->setNumElements(0);
  columnArray->setPackedMode(true);
#endif
}
/* Return <code>x * A + y</code> in <code>z</code>.
   Squashes small elements and knows about ClpSimplex */
void ClpPoolMatrix::transposeTimesByRow(const ClpSimplex *model, double scalar,
  const CoinIndexedVector *rowArray,
  CoinIndexedVector *y,
  CoinIndexedVector *columnArray) const
{
  printf("createMatrix at file %s line %d\n", __FILE__, __LINE__);
  createMatrix()->transposeTimesByRow(model, scalar, rowArray, y, columnArray);
}
/* Return <code>x *A in <code>z</code> but
   just for indices in y. */
void ClpPoolMatrix::subsetTransposeTimes(const ClpSimplex *model,
  const CoinIndexedVector *rowArray,
  const CoinIndexedVector *y,
  CoinIndexedVector *columnArray) const
{
#if 0
  printf("CcreateMatrix at file %s line %d\n",__FILE__,__LINE__);
  createMatrix()->subsetTransposeTimes(model,rowArray,y,columnArray);
#else
  columnArray->clear();
  double *COIN_RESTRICT pi = rowArray->denseVector();
  double *COIN_RESTRICT array = columnArray->denseVector();
  int jColumn;
  // get matrix data pointers
  const double *COIN_RESTRICT rowScale = model->rowScale();
  int numberToDo = y->getNumElements();
  const int *COIN_RESTRICT which = y->getIndices();
  assert(!rowArray->packedMode());
  columnArray->setPacked();
  if (!rowScale) {
    for (jColumn = 0; jColumn < numberToDo; jColumn++) {
      int iColumn = which[jColumn];
      double value = 0.0;
      CoinBigIndex start = columnStart_[iColumn];
      CoinBigIndex end = columnStart_[iColumn + 1];
      array[jColumn] = value;
      CoinBigIndex j;
      for (j = start; j < end; j++) {
        int iRow = stuff_[j].row_;
        value += pi[iRow] * elements_[stuff_[j].pool_];
      }
      array[jColumn] = value;
    }
  } else {
    // scaled
    const double *columnScale = model->columnScale();
    for (jColumn = 0; jColumn < numberToDo; jColumn++) {
      int iColumn = which[jColumn];
      double value = 0.0;
      CoinBigIndex start = columnStart_[iColumn];
      CoinBigIndex end = columnStart_[iColumn + 1];
      CoinBigIndex j;
      for (j = start; j < end; j++) {
        int iRow = stuff_[j].row_;
        value += pi[iRow] * elements_[stuff_[j].pool_] * rowScale[iRow];
      }
      array[jColumn] = value * columnScale[iColumn];
    }
  }
#endif
}
/// returns number of elements in column part of basis,
int ClpPoolMatrix::countBasis(const int *whichColumn,
  int &numberColumnBasic)
{
  int i;
  CoinBigIndex numberElements = 0;
  for (i = 0; i < numberColumnBasic; i++) {
    int iColumn = whichColumn[i];
    numberElements += columnStart_[iColumn + 1] - columnStart_[iColumn];
  }
  if (numberElements > COIN_INT_MAX) {
    printf("Factorization too large\n");
    abort();
  }
  return static_cast< int >(numberElements);
}
void ClpPoolMatrix::fillBasis(ClpSimplex *,
  const int *whichColumn,
  int &numberColumnBasic,
  int *indexRowU, int *start,
  int *rowCount, int *columnCount,
  CoinFactorizationDouble *elementU)
{
  CoinBigIndex numberElements = start[0];
  for (int i = 0; i < numberColumnBasic; i++) {
    int iColumn = whichColumn[i];
    CoinBigIndex j = columnStart_[iColumn];
    for (; j < columnStart_[iColumn + 1]; j++) {
      int iRow = stuff_[j].row_;
      indexRowU[numberElements] = iRow;
      rowCount[iRow]++;
      elementU[numberElements++] = elements_[stuff_[j].pool_];
      ;
    }
    start[i + 1] = static_cast< int >(numberElements);
    columnCount[i] = static_cast< int >(numberElements - start[i]);
  }
  if (numberElements > COIN_INT_MAX) {
    printf("Factorization too large\n");
    abort();
  }
}
/* Unpacks a column into an CoinIndexedvector
 */
void ClpPoolMatrix::unpack(const ClpSimplex *,
  CoinIndexedVector *rowArray,
  int iColumn) const
{
  CoinBigIndex j = columnStart_[iColumn];
  for (; j < columnStart_[iColumn + 1]; j++) {
    int iRow = stuff_[j].row_;
    rowArray->add(iRow, elements_[stuff_[j].pool_]);
  }
}
/* Unpacks a column into an CoinIndexedvector
** in packed foramt
Note that model is NOT const.  Bounds and objective could
be modified if doing column generation (just for this variable) */
void ClpPoolMatrix::unpackPacked(ClpSimplex *,
  CoinIndexedVector *rowArray,
  int iColumn) const
{
  int *COIN_RESTRICT index = rowArray->getIndices();
  double *COIN_RESTRICT array = rowArray->denseVector();
  int number = 0;
  CoinBigIndex j = columnStart_[iColumn];
  for (; j < columnStart_[iColumn + 1]; j++) {
    int iRow = stuff_[j].row_;
    array[number] = elements_[stuff_[j].pool_];
    index[number++] = iRow;
  }
  rowArray->setNumElements(number);
  rowArray->setPackedMode(true);
}
/* Adds multiple of a column into an CoinIndexedvector
   You can use quickAdd to add to vector */
void ClpPoolMatrix::add(const ClpSimplex *, CoinIndexedVector *rowArray,
  int iColumn, double multiplier) const
{
  CoinBigIndex j = columnStart_[iColumn];
  for (; j < columnStart_[iColumn + 1]; j++) {
    int iRow = stuff_[j].row_;
    rowArray->quickAdd(iRow, multiplier * elements_[stuff_[j].pool_]);
  }
}
/* Adds multiple of a column into an array */
void ClpPoolMatrix::add(const ClpSimplex *, double *array,
  int iColumn, double multiplier) const
{
  CoinBigIndex j = columnStart_[iColumn];
  for (; j < columnStart_[iColumn + 1]; j++) {
    int iRow = stuff_[j].row_;
    array[iRow] += multiplier * elements_[stuff_[j].pool_];
  }
}

// Return a complete CoinPackedMatrix
CoinPackedMatrix *
ClpPoolMatrix::getPackedMatrix() const
{
  createMatrix();
  return matrix_->matrix();
}
/* A vector containing the elements in the packed matrix. Note that there
   might be gaps in this list, entries that do not belong to any
   major-dimension vector. To get the actual elements one should look at
   this vector together with vectorStarts and vectorLengths. */
const double *
ClpPoolMatrix::getElements() const
{
  printf("createMatrix at file %s line %d\n", __FILE__, __LINE__);
  return createMatrix()->getElements();
}

/* A vector containing the minor indices of the elements in the packed
   matrix. Note that there might be gaps in this list, entries that do not
   belong to any major-dimension vector. To get the actual elements one
   should look at this vector together with vectorStarts and
   vectorLengths. */
const int *
ClpPoolMatrix::getIndices() const
{
  printf("createMatrix at file %s line %d\n", __FILE__, __LINE__);
  return createMatrix()->getIndices();
}
const CoinBigIndex *
ClpPoolMatrix::getVectorStarts() const
{
  printf("createMatrix at file %s line %d\n", __FILE__, __LINE__);
  return createMatrix()->getVectorStarts();
}
/* The lengths of the major-dimension vectors. */
const int *
ClpPoolMatrix::getVectorLengths() const
{
#if 0
  printf("XcreateMatrix at file %s line %d\n",__FILE__,__LINE__);
  return createMatrix()->getVectorLengths();
#else
  if (!lengths_) {
    lengths_ = new int[numberColumns_];
    for (int i = 0; i < numberColumns_; i++)
      lengths_[i] = static_cast< int >(columnStart_[i + 1] - columnStart_[i]);
  }
  return lengths_;
#endif
}
/* The length of a major-dimension vector. */
int ClpPoolMatrix::getVectorLength(int index) const
{
  return static_cast< int >(columnStart_[index + 1] - columnStart_[index]);
}
/* Delete the columns whose indices are listed in <code>indDel</code>. */
void ClpPoolMatrix::deleteCols(const int numDel, const int *indDel)
{
  int iColumn;
  CoinBigIndex newSize = columnStart_[numberColumns_];
  ;
  int numberBad = 0;
  // Use array to make sure we can have duplicates
  int *which = new int[numberColumns_];
  memset(which, 0, numberColumns_ * sizeof(int));
  int nDuplicate = 0;
  for (iColumn = 0; iColumn < numDel; iColumn++) {
    int jColumn = indDel[iColumn];
    if (jColumn < 0 || jColumn >= numberColumns_) {
      numberBad++;
    } else {
      newSize -= columnStart_[jColumn + 1] - columnStart_[jColumn];
      if (which[jColumn])
        nDuplicate++;
      else
        which[jColumn] = 1;
    }
  }
  if (numberBad)
    throw CoinError("Indices out of range", "deleteCols", "ClpPoolMatrix");
  int newNumber = numberColumns_ - numDel + nDuplicate;
  // Get rid of temporary arrays
  delete matrix_;
  matrix_ = NULL;
  CoinBigIndex *columnStart = new CoinBigIndex[newNumber + 1];
  poolInfo *stuff = new poolInfo[newSize];
  newNumber = 0;
  newSize = 0;
  for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
    if (!which[iColumn]) {
      CoinBigIndex start, end;
      CoinBigIndex i;
      start = columnStart_[iColumn];
      end = columnStart_[iColumn + 1];
      columnStart[newNumber] = newSize;
      for (i = start; i < end; i++)
        stuff[newSize++] = stuff_[i];
    }
  }
  columnStart[newNumber] = newSize;
  delete[] which;
  delete[] columnStart_;
  columnStart_ = columnStart;
  delete[] stuff_;
  stuff_ = stuff;
  numberColumns_ = newNumber;
}
/* Delete the rows whose indices are listed in <code>indDel</code>. */
void ClpPoolMatrix::deleteRows(const int numDel, const int *indDel)
{
  int iRow;
  int numberBad = 0;
  // Use array to make sure we can have duplicates
  int *which = new int[numberRows_];
  memset(which, 0, numberRows_ * sizeof(int));
  int nDuplicate = 0;
  for (iRow = 0; iRow < numDel; iRow++) {
    int jRow = indDel[iRow];
    if (jRow < 0 || jRow >= numberRows_) {
      numberBad++;
    } else {
      if (which[jRow])
        nDuplicate++;
      else
        which[jRow] = 1;
    }
  }
  if (numberBad)
    throw CoinError("Indices out of range", "deleteRows", "ClpPoolMatrix");
  CoinBigIndex iElement;
  CoinBigIndex numberElements = columnStart_[numberColumns_];
  CoinBigIndex newSize = 0;
  for (iElement = 0; iElement < numberElements; iElement++) {
    iRow = stuff_[iElement].row_;
    if (!which[iRow])
      newSize++;
  }
  int newNumber = numberRows_ - numDel + nDuplicate;
  // Get rid of temporary arrays
  delete matrix_;
  matrix_ = NULL;
  poolInfo *stuff = new poolInfo[newSize];
  newSize = 0;
  // redo which
  int numberRows = 0;
  for (iRow = 0; iRow < numberRows_; iRow++) {
    if (which[iRow]) {
      which[iRow] = -1; // not wanted
    } else {
      which[iRow] = numberRows;
      numberRows++;
    }
  }
  int iColumn;
  for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
    CoinBigIndex start, end;
    CoinBigIndex i;
    start = columnStart_[iColumn];
    end = columnStart_[iColumn + 1];
    columnStart_[newNumber] = newSize;
    for (i = start; i < end; i++) {
      poolInfo newStuff = stuff_[i];
      iRow = which[newStuff.row_];
      if (iRow >= 0) {
        newStuff.row_ = iRow;
        stuff[newSize++] = newStuff;
      }
    }
  }
  columnStart_[numberColumns_] = newSize;
  delete[] which;
  delete[] stuff_;
  stuff_ = stuff;
  numberRows_ = newNumber;
}
bool ClpPoolMatrix::isColOrdered() const
{
  return true;
}
/* Number of entries in the packed matrix. */
CoinBigIndex
ClpPoolMatrix::getNumElements() const
{
  if (numberColumns_)
    return columnStart_[numberColumns_];
  else
    return 0;
}
/* Returns largest and smallest elements of both signs.
   Largest refers to largest absolute value.
*/
void ClpPoolMatrix::rangeOfElements(double &smallestNegative, double &largestNegative,
  double &smallestPositive, double &largestPositive)
{
  smallestNegative = -COIN_DBL_MAX;
  largestNegative = 0.0;
  smallestPositive = COIN_DBL_MAX;
  largestPositive = 0.0;
  for (int i = 0; i < numberDifferent_; i++) {
    double value = elements_[i];
    if (value > 0.0) {
      smallestPositive = std::min(smallestPositive, value);
      largestPositive = std::max(largestPositive, value);
    } else if (value < 0.0) {
      smallestNegative = std::max(smallestNegative, value);
      largestNegative = std::min(largestNegative, value);
    }
  }
}
// Says whether it can do partial pricing
bool ClpPoolMatrix::canDoPartialPricing() const
{
  return true;
}
// Partial pricing
void ClpPoolMatrix::partialPricing(ClpSimplex *model, double startFraction, double endFraction,
  int &bestSequence, int &numberWanted)
{
  printf("createMatrix at file %s line %d\n", __FILE__, __LINE__);
  createMatrix()->partialPricing(model, startFraction,
    endFraction, bestSequence,
    numberWanted);
}
// Allow any parts of a created CoinMatrix to be deleted
void ClpPoolMatrix::releasePackedMatrix() const
{
  delete matrix_;
  matrix_ = NULL;
}
// These have to match ClpPrimalColumnSteepest version
#define reference(i) (((reference[i >> 5] >> (i & 31)) & 1) != 0)
/* Updates two arrays for steepest and does devex weights
   Returns nonzero if updates reduced cost and infeas -
   new infeas in dj1  - This does not at present*/
int ClpPoolMatrix::transposeTimes2(const ClpSimplex *model,
  const CoinIndexedVector *pi1, CoinIndexedVector *dj1,
  const CoinIndexedVector *pi2,
  CoinIndexedVector *spare,
  double *infeas, double *reducedCost,
  double referenceIn, double devex,
  // Array for exact devex to say what is in reference framework
  unsigned int *COIN_RESTRICT reference,
  double *COIN_RESTRICT weights, double scaleFactor)
{
#if 0
  printf("CcreateMatrix at file %s line %d\n",__FILE__,__LINE__);
  return createMatrix()->transposeTimes2(model,pi1,dj1,pi2,spare,
					 infeas,reducedCost,referenceIn,devex,
					 reference,weights,scaleFactor);
#else
  int returnCode = 0;
  // put row of tableau in dj1
  double *COIN_RESTRICT pi = pi1->denseVector();
  int numberNonZero = 0;
  int *COIN_RESTRICT index = dj1->getIndices();
  double *COIN_RESTRICT array = dj1->denseVector();
  int numberInRowArray = pi1->getNumElements();
  double zeroTolerance = model->zeroTolerance();
  double dualTolerance = model->currentDualTolerance();
  // we can't really trust infeasibilities if there is dual error
  // this coding has to mimic coding in checkDualSolution
  double error = std::min(1.0e-2, model->largestDualError());
  // allow tolerance at least slightly bigger than standard
  dualTolerance = dualTolerance + error;
  bool packed = pi1->packedMode();
  assert(packed);
  // do by column
  int iColumn;
  const double *COIN_RESTRICT rowScale = model->rowScale();
  assert(!spare->getNumElements());
  assert(numberColumns_ > 0);
  double *COIN_RESTRICT piWeight = pi2->denseVector();
  assert(!pi2->packedMode());
  bool killDjs = (scaleFactor == 0.0);
  if (!scaleFactor)
    scaleFactor = 1.0;
  // need to expand pi into y
  assert(spare->capacity() >= model->numberRows());
  double *COIN_RESTRICT piOld = pi;
  pi = spare->denseVector();
  const int *COIN_RESTRICT whichRow = pi1->getIndices();
  int i;
  if (!rowScale) {
    // modify pi so can collapse to one loop
    for (i = 0; i < numberInRowArray; i++) {
      int iRow = whichRow[i];
      pi[iRow] = piOld[i];
    }
    if (infeas)
      returnCode = 1;
    CoinBigIndex j;
    CoinBigIndex end = columnStart_[0];
    for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
      CoinBigIndex start = end;
      end = columnStart_[iColumn + 1];
      ClpSimplex::Status status = model->getStatus(iColumn);
      if (status == ClpSimplex::basic || status == ClpSimplex::isFixed)
        continue;
      double value = 0.0;
      for (j = start; j < end; j++) {
        int iRow = stuff_[j].row_;
        value -= pi[iRow] * elements_[stuff_[j].pool_];
      }
      if (fabs(value) > zeroTolerance) {
        // and do other array
        double modification = 0.0;
        for (j = start; j < end; j++) {
          int iRow = stuff_[j].row_;
          modification += piWeight[iRow] * elements_[stuff_[j].pool_];
        }
        double thisWeight = weights[iColumn];
        double pivot = value * scaleFactor;
        double pivotSquared = pivot * pivot;
        thisWeight += pivotSquared * devex + pivot * modification;
        if (thisWeight < DEVEX_TRY_NORM) {
          if (referenceIn < 0.0) {
            // steepest
            thisWeight = std::max(DEVEX_TRY_NORM, DEVEX_ADD_ONE + pivotSquared);
          } else {
            // exact
            thisWeight = referenceIn * pivotSquared;
            if (reference(iColumn))
              thisWeight += 1.0;
            thisWeight = std::max(thisWeight, DEVEX_TRY_NORM);
          }
        }
        weights[iColumn] = thisWeight;
        if (!killDjs) {
          value = reducedCost[iColumn] - value;
          reducedCost[iColumn] = value;
          // simplify status
          ClpSimplex::Status status = model->getStatus(iColumn);

          switch (status) {

          case ClpSimplex::basic:
          case ClpSimplex::isFixed:
            break;
          case ClpSimplex::isFree:
          case ClpSimplex::superBasic:
            if (fabs(value) > FREE_ACCEPT * dualTolerance) {
              // we are going to bias towards free (but only if reasonable)
              value *= FREE_BIAS;
              value *= value;
              // store square in list
              if (infeas[iColumn]) {
                infeas[iColumn] = value; // already there
              } else {
                array[numberNonZero] = value;
                index[numberNonZero++] = iColumn;
              }
            } else {
              array[numberNonZero] = 0.0;
              index[numberNonZero++] = iColumn;
            }
            break;
          case ClpSimplex::atUpperBound:
            if (value > dualTolerance) {
              value *= value;
              // store square in list
              if (infeas[iColumn]) {
                infeas[iColumn] = value; // already there
              } else {
                array[numberNonZero] = value;
                index[numberNonZero++] = iColumn;
              }
            } else {
              array[numberNonZero] = 0.0;
              index[numberNonZero++] = iColumn;
            }
            break;
          case ClpSimplex::atLowerBound:
            if (value < -dualTolerance) {
              value *= value;
              // store square in list
              if (infeas[iColumn]) {
                infeas[iColumn] = value; // already there
              } else {
                array[numberNonZero] = value;
                index[numberNonZero++] = iColumn;
              }
            } else {
              array[numberNonZero] = 0.0;
              index[numberNonZero++] = iColumn;
            }
          }
        }
      }
    }
  } else {
    // scaled
    // modify pi so can collapse to one loop
    for (i = 0; i < numberInRowArray; i++) {
      int iRow = whichRow[i];
      pi[iRow] = piOld[i] * rowScale[iRow];
    }
    // can also scale piWeight as not used again
    int numberWeight = pi2->getNumElements();
    const int *indexWeight = pi2->getIndices();
    for (i = 0; i < numberWeight; i++) {
      int iRow = indexWeight[i];
      piWeight[iRow] *= rowScale[iRow];
    }
    if (infeas)
      returnCode = 1;
    const double *COIN_RESTRICT columnScale = model->columnScale();
    CoinBigIndex j;
    CoinBigIndex end = columnStart_[0];
    for (iColumn = 0; iColumn < numberColumns_; iColumn++) {
      CoinBigIndex start = end;
      end = columnStart_[iColumn + 1];
      ClpSimplex::Status status = model->getStatus(iColumn);
      if (status == ClpSimplex::basic || status == ClpSimplex::isFixed)
        continue;
      double scale = columnScale[iColumn];
      double value = 0.0;
      for (j = start; j < end; j++) {
        int iRow = stuff_[j].row_;
        value -= pi[iRow] * elements_[stuff_[j].pool_];
      }
      value *= scale;
      if (fabs(value) > zeroTolerance) {
        double modification = 0.0;
        for (j = start; j < end; j++) {
          int iRow = stuff_[j].row_;
          modification += piWeight[iRow] * elements_[stuff_[j].pool_];
        }
        modification *= scale;
        double thisWeight = weights[iColumn];
        double pivot = value * scaleFactor;
        double pivotSquared = pivot * pivot;
        thisWeight += pivotSquared * devex + pivot * modification;
        if (thisWeight < DEVEX_TRY_NORM) {
          if (referenceIn < 0.0) {
            // steepest
            thisWeight = std::max(DEVEX_TRY_NORM, DEVEX_ADD_ONE + pivotSquared);
          } else {
            // exact
            thisWeight = referenceIn * pivotSquared;
            if (reference(iColumn))
              thisWeight += 1.0;
            thisWeight = std::max(thisWeight, DEVEX_TRY_NORM);
          }
        }
        weights[iColumn] = thisWeight;
        if (!killDjs) {
          value = reducedCost[iColumn] - value;
          reducedCost[iColumn] = value;
          // simplify status
          ClpSimplex::Status status = model->getStatus(iColumn);

          switch (status) {

          case ClpSimplex::basic:
          case ClpSimplex::isFixed:
            break;
          case ClpSimplex::isFree:
          case ClpSimplex::superBasic:
            if (fabs(value) > FREE_ACCEPT * dualTolerance) {
              // we are going to bias towards free (but only if reasonable)
              value *= FREE_BIAS;
              value *= value;
              // store square in list
              if (infeas[iColumn]) {
                infeas[iColumn] = value; // already there
              } else {
                array[numberNonZero] = value;
                index[numberNonZero++] = iColumn;
              }
            } else {
              array[numberNonZero] = 0.0;
              index[numberNonZero++] = iColumn;
            }
            break;
          case ClpSimplex::atUpperBound:
            if (value > dualTolerance) {
              value *= value;
              // store square in list
              if (infeas[iColumn]) {
                infeas[iColumn] = value; // already there
              } else {
                array[numberNonZero] = value;
                index[numberNonZero++] = iColumn;
              }
            } else {
              array[numberNonZero] = 0.0;
              index[numberNonZero++] = iColumn;
            }
            break;
          case ClpSimplex::atLowerBound:
            if (value < -dualTolerance) {
              value *= value;
              // store square in list
              if (infeas[iColumn]) {
                infeas[iColumn] = value; // already there
              } else {
                array[numberNonZero] = value;
                index[numberNonZero++] = iColumn;
              }
            } else {
              array[numberNonZero] = 0.0;
              index[numberNonZero++] = iColumn;
            }
          }
        }
      }
    }
  }
  // zero out
  int numberRows = model->numberRows();
  if (numberInRowArray * 4 < numberRows) {
    for (i = 0; i < numberInRowArray; i++) {
      int iRow = whichRow[i];
      pi[iRow] = 0.0;
    }
  } else {
    CoinZeroN(pi, numberRows);
  }
  dj1->setNumElements(numberNonZero);
  spare->setNumElements(0);
  if (packed)
    dj1->setPackedMode(true);
  return returnCode;
#endif
}
// Updates second array for steepest and does devex weights
void ClpPoolMatrix::subsetTimes2(const ClpSimplex *model,
  CoinIndexedVector *dj1,
  const CoinIndexedVector *pi2, CoinIndexedVector *dj2,
  double referenceIn, double devex,
  // Array for exact devex to say what is in reference framework
  unsigned int *COIN_RESTRICT reference,
  double *COIN_RESTRICT weights, double scaleFactor)
{
  printf("createMatrix at file %s line %d\n", __FILE__, __LINE__);
  createMatrix()->subsetTimes2(model, dj1, pi2, dj2,
    referenceIn, devex,
    reference, weights, scaleFactor);
}
/* Set the dimensions of the matrix. In effect, append new empty
   columns/rows to the matrix. A negative number for either dimension
   means that that dimension doesn't change. Otherwise the new dimensions
   MUST be at least as large as the current ones otherwise an exception
   is thrown. */
void ClpPoolMatrix::setDimensions(int newnumrows, int newnumcols)
{
  // at present abort
  abort();
}
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
