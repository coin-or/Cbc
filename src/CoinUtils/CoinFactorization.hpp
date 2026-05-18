// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/* 
   Authors
   
   John Forrest

 */
#ifndef CoinFactorization_H
#define CoinFactorization_H
#define EXTRA_U_SPACE 4
//#define COIN_ONE_ETA_COPY 100

#include <iostream>
#include <string>
#include <cassert>
#include <cstdio>
#include <cmath>
#include "CoinTypes.h"
#include "CoinIndexedVector.hpp"

class CoinPackedMatrix;
/** This deals with Factorization and Updates

    This class started with a parallel simplex code I was writing in the
    mid 90's.  The need for parallelism led to many complications and
    I have simplified as much as I could to get back to this.

    I was aiming at problems where I might get speed-up so I was looking at dense
    problems or ones with structure.  This led to permuting input and output
    vectors and to increasing the number of rows each rank-one update.  This is 
    still in as a minor overhead.

    I have also put in handling for hyper-sparsity.  I have taken out
    all outer loop unrolling, dense matrix handling and most of the
    book-keeping for slacks.  Also I always use FTRAN approach to updating
    even if factorization fairly dense.  All these could improve performance.

    I blame some of the coding peculiarities on the history of the code
    but mostly it is just because I can't do elegant code (or useful
    comments).

    I am assuming that 32 bits is enough for number of rows or columns, but CoinBigIndex
    may be redefined to get 64 bits.
 */

class COINUTILSLIB_EXPORT CoinFactorization {
  friend void CoinFactorizationUnitTest(const std::string &mpsDir);

public:
  /**@name Constructors and destructor and copy */
  //@{
  /// Default constructor
  CoinFactorization();
  /// Copy constructor
  CoinFactorization(const CoinFactorization &other);

  /// Destructor
  ~CoinFactorization();
  /// Delete all stuff (leaves as after CoinFactorization())
  void almostDestructor();
  /// Debug show object (shows one representation)
  void show_self() const;
#if COINUTILS_BIGINDEX_IS_INT
  /// Debug - save on file - 0 if no error
  int saveFactorization(const char *file) const;
  /** Debug - restore from file - 0 if no error on file.
      If factor true then factorizes as if called from ClpFactorization
  */
  int restoreFactorization(const char *file, bool factor = false);
#endif
  /// Debug - sort so can compare
  void sort() const;
  /// = copy
  CoinFactorization &operator=(const CoinFactorization &other);
  //@}

  /**@name Do factorization */
  //@{
  /** When part of LP - given by basic variables.
  Actually does factorization.
  Arrays passed in have non negative value to say basic.
  If status is okay, basic variables have pivot row - this is only needed
  If status is singular, then basic variables have pivot row
  and ones thrown out have -1
  returns 0 -okay, -1 singular, -2 too many in basis, -99 memory */
  int factorize(const CoinPackedMatrix &matrix,
    int rowIsBasic[], int columnIsBasic[],
    double areaFactor = 0.0);
  /** When given as triplets.
  Actually does factorization.  maximumL is guessed maximum size of L part of
  final factorization, maximumU of U part.  These are multiplied by
  areaFactor which can be computed by user or internally.  
  Arrays are copied in.  I could add flag to delete arrays to save a 
  bit of memory.
  If status okay, permutation has pivot rows - this is only needed
  If status is singular, then basic variables have pivot row
  and ones thrown out have -1
  returns 0 -okay, -1 singular, -99 memory */
  int factorize(int numberRows,
    int numberColumns,
    int numberElements,
    int maximumL,
    int maximumU,
    const int indicesRow[],
    const int indicesColumn[], const double elements[],
    int permutation[],
    double areaFactor = 0.0);
  /** Two part version for maximum flexibility
      This part creates arrays for user to fill.
      estimateNumberElements is safe estimate of number
      returns 0 -okay, -99 memory */
  int factorizePart1(int numberRows,
    int numberColumns,
    int estimateNumberElements,
    int *COIN_RESTRICT indicesRow[],
    int *COIN_RESTRICT indicesColumn[],
    CoinFactorizationDouble *COIN_RESTRICT elements[],
    double areaFactor = 0.0);
  /** This is part two of factorization
      Arrays belong to factorization and were returned by part 1
      If status okay, permutation has pivot rows - this is only needed
      If status is singular, then basic variables have pivot row
      and ones thrown out have -1
      returns 0 -okay, -1 singular, -99 memory */
  int factorizePart2(int permutation[], int exactNumberElements);
  /// Condition number - product of pivots after factorization
  double conditionNumber() const;

  //@}

  /**@name general stuff such as permutation or status */
  //@{
  /// Returns status
  inline int status() const
  {
    return status_;
  }
  /// Sets status
  inline void setStatus(int value)
  {
    status_ = value;
  }
  /// Returns number of pivots since factorization
  inline int pivots() const
  {
    return numberPivots_;
  }
  /// Sets number of pivots since factorization
  inline void setPivots(int value)
  {
    numberPivots_ = value;
  }
  /// Returns address of permute region
  inline int *permute() const
  {
    return permute_.array();
  }
  /// Returns address of pivotColumn region (also used for permuting)
  inline int *pivotColumn() const
  {
    return pivotColumn_.array();
  }
  /// Returns address of pivot region
  inline CoinFactorizationDouble *pivotRegion() const
  {
    return pivotRegion_.array();
  }
  /// Returns address of permuteBack region
  inline int *permuteBack() const
  {
    return permuteBack_.array();
  }
  /// Returns address of lastRow region
  inline int *lastRow() const
  {
    return lastRow_.array();
  }
  /** Returns address of pivotColumnBack region (also used for permuting)
      Now uses firstCount to save memory allocation */
  inline int *pivotColumnBack() const
  {
    //return firstCount_.array();
    return pivotColumnBack_.array();
  }
  /// Start of each row in L
  inline int *startRowL() const
  {
    return startRowL_.array();
  }

  /// Start of each column in L
  inline int *startColumnL() const
  {
    return startColumnL_.array();
  }

  /// Index of column in row for L
  inline int *indexColumnL() const
  {
    return indexColumnL_.array();
  }

  /// Row indices of L
  inline int *indexRowL() const
  {
    return indexRowL_.array();
  }

  /// Elements in L (row copy)
  inline CoinFactorizationDouble *elementByRowL() const
  {
    return elementByRowL_.array();
  }

  /// Number of Rows after iterating
  inline int numberRowsExtra() const
  {
    return numberRowsExtra_;
  }
  /// Set number of Rows after factorization
  inline void setNumberRows(int value)
  {
    numberRows_ = value;
  }
  /// Number of Rows after factorization
  inline int numberRows() const
  {
    return numberRows_;
  }
  /// Number in L
  inline int numberL() const
  {
    return numberL_;
  }

  /// Base of L
  inline int baseL() const
  {
    return baseL_;
  }
  /// Maximum of Rows after iterating
  inline int maximumRowsExtra() const
  {
    return maximumRowsExtra_;
  }
  /// Total number of columns in factorization
  inline int numberColumns() const
  {
    return numberColumns_;
  }
  /// Total number of elements in factorization
  inline int numberElements() const
  {
    return totalElements_;
  }
  /// Length of FT vector
  inline int numberForrestTomlin() const
  {
    return numberInColumn_.array()[numberColumnsExtra_];
  }
  /// Number of good columns in factorization
  inline int numberGoodColumns() const
  {
    return numberGoodU_;
  }
  /// Whether larger areas needed
  inline double areaFactor() const
  {
    return areaFactor_;
  }
  inline void areaFactor(double value)
  {
    areaFactor_ = value;
  }
  /// Returns areaFactor but adjusted for dense
  double adjustedAreaFactor() const;
  /// Allows change of pivot accuracy check 1.0 == none >1.0 relaxed
  inline void relaxAccuracyCheck(double value)
  {
    relaxCheck_ = value;
  }
  inline double getAccuracyCheck() const
  {
    return relaxCheck_;
  }
  /// Level of detail of messages
  inline int messageLevel() const
  {
    return messageLevel_;
  }
  void messageLevel(int value);
  /// Maximum number of pivots between factorizations
  inline int maximumPivots() const
  {
    return maximumPivots_;
  }
  void maximumPivots(int value);

  /// Gets dense threshold
  inline int denseThreshold() const
  {
    return denseThreshold_;
  }
  /// Sets dense threshold
  inline void setDenseThreshold(int value)
  {
    denseThreshold_ = value;
  }
  /// Pivot tolerance
  inline double pivotTolerance() const
  {
    return pivotTolerance_;
  }
  void pivotTolerance(double value);
  /// Zero tolerance
  inline double zeroTolerance() const
  {
    return zeroTolerance_;
  }
  void zeroTolerance(double value);
#ifndef COIN_FAST_CODE
  /// Whether slack value is +1 or -1
  inline double slackValue() const
  {
    return slackValue_;
  }
  void slackValue(double value);
#endif
  /// Returns maximum absolute value in factorization
  double maximumCoefficient() const;
  /// true if Forrest Tomlin update, false if PFI
  inline bool forrestTomlin() const
  {
    return doForrestTomlin_;
  }
  inline void setForrestTomlin(bool value)
  {
    doForrestTomlin_ = value;
  }
  /// True if FT update and space
  inline bool spaceForForrestTomlin() const
  {
    int start = startColumnU_.array()[maximumColumnsExtra_];
    int space = lengthAreaU_ - (start + numberRowsExtra_);
    return (space >= 0) && doForrestTomlin_;
  }
  //@}

  /**@name some simple stuff */
  //@{

  /// Returns number of dense rows
  inline int numberDense() const
  {
    return numberDense_;
  }

  /// Returns number in U area
  inline int numberElementsU() const
  {
    return lengthU_;
  }
  /// Setss number in U area
  inline void setNumberElementsU(int value)
  {
    lengthU_ = value;
  }
  /// Returns length of U area
  inline int lengthAreaU() const
  {
    return lengthAreaU_;
  }
  /// Returns number in L area
  inline int numberElementsL() const
  {
    return lengthL_;
  }
  /// Returns length of L area
  inline int lengthAreaL() const
  {
    return lengthAreaL_;
  }
  /// Returns number in R area
  inline int numberElementsR() const
  {
    return lengthR_;
  }
  /// Number of compressions done
  inline int numberCompressions() const
  {
    return numberCompressions_;
  }
  /// Number of entries in each row
  inline int *numberInRow() const
  {
    return numberInRow_.array();
  }
  /// Number of entries in each column
  inline int *numberInColumn() const
  {
    return numberInColumn_.array();
  }
  /// Elements of U
  inline CoinFactorizationDouble *elementU() const
  {
    return elementU_.array();
  }
  /// Row indices of U
  inline int *indexRowU() const
  {
    return indexRowU_.array();
  }
  /// Start of each column in U
  inline int *startColumnU() const
  {
    return startColumnU_.array();
  }
  /// Maximum number of Columns after iterating
  inline int maximumColumnsExtra()
  {
    return maximumColumnsExtra_;
  }
  /** L to U bias
      0 - U bias, 1 - some U bias, 2 some L bias, 3 L bias
  */
  inline int biasLU() const
  {
    return biasLU_;
  }
  inline void setBiasLU(int value)
  {
    biasLU_ = value;
  }
  /** Array persistence flag
      If 0 then as now (delete/new)
      1 then only do arrays if bigger needed
      2 as 1 but give a bit extra if bigger needed
  */
  inline int persistenceFlag() const
  {
    return persistenceFlag_;
  }
  void setPersistenceFlag(int value);
  //@}

  /**@name rank one updates which do exist */
  //@{

  /** Replaces one Column to basis,
   returns 0=OK, 1=Probably OK, 2=singular, 3=no room
      If checkBeforeModifying is true will do all accuracy checks
      before modifying factorization.  Whether to set this depends on
      speed considerations.  You could just do this on first iteration
      after factorization and thereafter re-factorize
   partial update already in U */
  int replaceColumn(CoinIndexedVector *regionSparse,
    int pivotRow,
    double pivotCheck,
    bool checkBeforeModifying = false,
    double acceptablePivot = 1.0e-8);
  /** Combines BtranU and delete elements
      If deleted is NULL then delete elements
      otherwise store where elements are
  */
  void replaceColumnU(CoinIndexedVector *regionSparse,
    int *COIN_RESTRICT deleted,
    int internalPivotRow);
#ifdef ABC_USE_COIN_FACTORIZATION
  /** returns empty fake vector carved out of existing
      later - maybe use associated arrays */
  CoinIndexedVector *fakeVector(CoinIndexedVector *vector,
    int already = 0) const;
  void deleteFakeVector(CoinIndexedVector *vector,
    CoinIndexedVector *fakeVector) const;
  /** Checks if can replace one Column to basis,
      returns update alpha
      Fills in region for use later
      partial update already in U */
  double checkReplacePart1(CoinIndexedVector *regionSparse,
    int pivotRow);
  /** Checks if can replace one Column to basis,
      returns update alpha
      Fills in region for use later
      partial update in vector */
  double checkReplacePart1(CoinIndexedVector *regionSparse,
    CoinIndexedVector *partialUpdate,
    int pivotRow);
  /** Checks if can replace one Column in basis,
      returns 0=OK, 1=Probably OK, 2=singular, 3=no room, 5 max pivots */
  int checkReplacePart2(int pivotRow,
    double btranAlpha,
    double ftranAlpha,
    double ftAlpha,
    double acceptablePivot = 1.0e-8);
  /** Replaces one Column to basis,
      partial update already in U */
  void replaceColumnPart3(CoinIndexedVector *regionSparse,
    int pivotRow,
    double alpha);
  /** Replaces one Column to basis,
      partial update in vector */
  void replaceColumnPart3(CoinIndexedVector *regionSparse,
    CoinIndexedVector *partialUpdate,
    int pivotRow,
    double alpha);
  /** Updates one column (FTRAN) from regionSparse2
      Tries to do FT update
      number returned is negative if no room
      regionSparse starts as zero and is zero at end.
      Note - if regionSparse2 packed on input - will be packed on output
      long regions
  */
  int updateColumnFT(CoinIndexedVector &regionSparse);
  int updateColumnFTPart1(CoinIndexedVector &regionSparse);
  void updateColumnFTPart2(CoinIndexedVector &regionSparse);
  /** Updates one column (FTRAN) - long region
      Tries to do FT update
      puts partial update in vector */
  void updateColumnFT(CoinIndexedVector &regionSparseFT,
    CoinIndexedVector &partialUpdate,
    int which);
  /** Updates one column (FTRAN) long region */
  int updateColumn(CoinIndexedVector &regionSparse) const;
  /** Updates one column (FTRAN) from regionFT
      Tries to do FT update
      number returned is negative if no room.
      Also updates regionOther - long region*/
  int updateTwoColumnsFT(CoinIndexedVector &regionSparseFT,
    CoinIndexedVector &regionSparseOther);
  /** Updates one column (BTRAN) - long region*/
  int updateColumnTranspose(CoinIndexedVector &regionSparse) const;
  /** Updates one column (FTRAN) - long region */
  void updateColumnCpu(CoinIndexedVector &regionSparse, int whichCpu) const;
  /** Updates one column (BTRAN) - long region */
  void updateColumnTransposeCpu(CoinIndexedVector &regionSparse, int whichCpu) const;
  /** Updates one full column (FTRAN) - long region */
  void updateFullColumn(CoinIndexedVector &regionSparse) const;
  /** Updates one full column (BTRAN) - long region */
  void updateFullColumnTranspose(CoinIndexedVector &regionSparse) const;
  /** Updates one column for dual steepest edge weights (FTRAN) - long region */
  void updateWeights(CoinIndexedVector &regionSparse) const;
  /// Returns true if wants tableauColumn in replaceColumn
  inline bool wantsTableauColumn() const
  {
    return false;
  }
  /// Pivot tolerance
  inline double minimumPivotTolerance() const
  {
    return pivotTolerance_;
  }
  inline void minimumPivotTolerance(double value)
  {
    pivotTolerance(value);
  }
  /// Says parallel
  inline void setParallelMode(int value)
  {
    parallelMode_ = value;
  }
  /// Sets solve mode
  inline void setSolveMode(int value)
  {
    parallelMode_ &= 3;
    parallelMode_ |= (value << 2);
  }
  /// Sets solve mode
  inline int solveMode() const
  {
    return parallelMode_ >> 2;
  }
  /// Update partial Ftran by R update
  void updatePartialUpdate(CoinIndexedVector &partialUpdate);
  /// Makes a non-singular basis by replacing variables
  void makeNonSingular(int *COIN_RESTRICT sequence);
#endif
#if ABOCA_LITE_FACTORIZATION
  /// Does btranU part of replaceColumn (skipping entries)
  void replaceColumn1(CoinIndexedVector *regionSparse, int pivotRow);
  /// Does replaceColumn - having already done btranU
  int replaceColumn2(CoinIndexedVector *regionSparse,
    int pivotRow,
    double pivotCheck);
#endif
  //@}

  /**@name various uses of factorization (return code number elements) 
   which user may want to know about */
  //@{
  /** Updates one column (FTRAN) from regionSparse2
      Tries to do FT update
      number returned is negative if no room
      regionSparse starts as zero and is zero at end.
      Note - if regionSparse2 packed on input - will be packed on output
  */
  int updateColumnFT(CoinIndexedVector *regionSparse,
    CoinIndexedVector *regionSparse2);
  /** This version has same effect as above with FTUpdate==false
      so number returned is always >=0 */
  int updateColumn(CoinIndexedVector *regionSparse,
    CoinIndexedVector *regionSparse2,
    bool noPermute = false) const;
  /** Updates one column (FTRAN) from region2
      Tries to do FT update
      number returned is negative if no room.
      Also updates region3
      region1 starts as zero and is zero at end */
  int updateTwoColumnsFT(CoinIndexedVector *regionSparse1,
    CoinIndexedVector *regionSparse2,
    CoinIndexedVector *regionSparse3,
    bool noPermuteRegion3 = false);
  /** Updates one column (BTRAN) from regionSparse2
      regionSparse starts as zero and is zero at end 
      Note - if regionSparse2 packed on input - will be packed on output
  */
  int updateColumnTranspose(CoinIndexedVector *regionSparse,
    CoinIndexedVector *regionSparse2) const;
  /// Part of twocolumnsTranspose
  void updateOneColumnTranspose(CoinIndexedVector *regionWork, int &statistics) const;
  /** Updates two columns (BTRAN) from regionSparse2 and 3
      regionSparse starts as zero and is zero at end 
      Note - if regionSparse2 packed on input - will be packed on output - same for 3
  */
  void updateTwoColumnsTranspose(CoinIndexedVector *regionSparse,
    CoinIndexedVector *regionSparse2,
    CoinIndexedVector *regionSparse3,
    int type) const;
  /** makes a row copy of L for speed and to allow very sparse problems */
  void goSparse();
  /**  get sparse threshold */
  inline int sparseThreshold() const
  {
    return sparseThreshold_;
  }
  /**  set sparse threshold */
  void sparseThreshold(int value);
  //@}
  /// *** Below this user may not want to know about

  /**@name various uses of factorization (return code number elements) 
   which user may not want to know about (left over from my LP code) */
  //@{
  /// Get rid of all memory
  inline void clearArrays()
  {
    gutsOfDestructor();
  }
  //@}

  /**@name various updates - none of which have been written! */
  //@{

  /** Adds given elements to Basis and updates factorization,
      can increase size of basis. Returns rank */
  int add(int numberElements,
    int indicesRow[],
    int indicesColumn[], double elements[]);

  /** Adds one Column to basis,
      can increase size of basis. Returns rank */
  int addColumn(int numberElements,
    int indicesRow[], double elements[]);

  /** Adds one Row to basis,
      can increase size of basis. Returns rank */
  int addRow(int numberElements,
    int indicesColumn[], double elements[]);

  /// Deletes one Column from basis, returns rank
  int deleteColumn(int Row);
  /// Deletes one Row from basis, returns rank
  int deleteRow(int Row);

  /** Replaces one Row in basis,
      At present assumes just a singleton on row is in basis
      returns 0=OK, 1=Probably OK, 2=singular, 3 no space */
  int replaceRow(int whichRow, int numberElements,
    const int indicesColumn[], const double elements[]);
  /// Takes out all entries for given rows
  void emptyRows(int numberToEmpty, const int which[]);
  //@}
  /**@name used by ClpFactorization */
  /// See if worth going sparse
  void checkSparse();
  /// For statistics
  inline bool collectStatistics() const
  {
    return true;
  }
  /// For statistics
  inline void setCollectStatistics(bool onOff) const
  {
  }
  /// The real work of constructors etc 0 just scalars, 1 bit normal
  void gutsOfDestructor(int type = 1);
  /// 1 bit - tolerances etc, 2 more, 4 dummy arrays
  void gutsOfInitialize(int type);
  void gutsOfCopy(const CoinFactorization &other);

  /// Reset all sparsity etc statistics
  void resetStatistics();

  //@}

  /**@name used by factorization */
  /// Gets space for a factorization, called by constructors
  void getAreas(int numberRows,
    int numberColumns,
    int maximumL,
    int maximumU);

  /** PreProcesses raw triplet data.
      state is 0 - triplets, 1 - some counts etc , 2 - .. */
  void preProcess(int state,
    int possibleDuplicates = -1);
  /// Does most of factorization
  int factor();

protected:
  /** Does sparse phase of factorization
      return code is <0 error, 0= finished */
  int factorSparse();
  /** Does sparse phase of factorization (for smaller problems)
      return code is <0 error, 0= finished */
  int factorSparseSmall();
  /** Does sparse phase of factorization (for larger problems)
      return code is <0 error, 0= finished */
  int factorSparseLarge();
  /** Does dense phase of factorization
      return code is <0 error, 0= finished */
  int factorDense();

  /// Pivots when just one other row so faster?
  bool pivotOneOtherRow(int pivotRow,
    int pivotColumn);
  /// Does one pivot on Row Singleton in factorization
  bool pivotRowSingleton(int pivotRow,
    int pivotColumn);
  /// Does one pivot on Column Singleton in factorization
  bool pivotColumnSingleton(int pivotRow,
    int pivotColumn);

  /** Gets space for one Column with given length,
   may have to do compression  (returns True if successful),
   also moves existing vector,
   extraNeeded is over and above present */
  bool getColumnSpace(int iColumn,
    int extraNeeded);

  /** Reorders U so contiguous and in order (if there is space)
      Returns true if it could */
  bool reorderU();
  /**  getColumnSpaceIterateR.  Gets space for one extra R element in Column
       may have to do compression  (returns true)
       also moves existing vector */
  bool getColumnSpaceIterateR(int iColumn, double value,
    int iRow);
  /**  getColumnSpaceIterate.  Gets space for one extra U element in Column
       may have to do compression  (returns true)
       also moves existing vector.
       Returns -1 if no memory or where element was put
       Used by replaceRow (turns off R version) */
  int getColumnSpaceIterate(int iColumn, double value,
    int iRow);
  /** Gets space for one Row with given length,
  may have to do compression  (returns True if successful),
  also moves existing vector */
  bool getRowSpace(int iRow, int extraNeeded);

  /** Gets space for one Row with given length while iterating,
  may have to do compression  (returns True if successful),
  also moves existing vector */
  bool getRowSpaceIterate(int iRow,
    int extraNeeded);
  /// Checks that row and column copies look OK
  void checkConsistency();
  /// Adds a link in chain of equal counts
  inline void addLink(int index, int count)
  {
    int *COIN_RESTRICT nextCount = nextCount_.array();
    int *COIN_RESTRICT firstCount = firstCount_.array();
    int *COIN_RESTRICT lastCount = lastCount_.array();
    int next = firstCount[count];
    lastCount[index] = -2 - count;
    if (next < 0) {
      //first with that count
      firstCount[count] = index;
      nextCount[index] = -1;
    } else {
      firstCount[count] = index;
      nextCount[index] = next;
      lastCount[next] = index;
    }
  }
  /// Deletes a link in chain of equal counts
  inline void deleteLink(int index)
  {
    int *COIN_RESTRICT nextCount = nextCount_.array();
    int *COIN_RESTRICT firstCount = firstCount_.array();
    int *COIN_RESTRICT lastCount = lastCount_.array();
    int next = nextCount[index];
    int last = lastCount[index];
    if (last >= 0) {
      nextCount[last] = next;
    } else {
      int count = -last - 2;

      firstCount[count] = next;
    }
    if (next >= 0) {
      lastCount[next] = last;
    }
    nextCount[index] = -2;
    lastCount[index] = -2;
    return;
  }
  /// Separate out links with same row/column count
  void separateLinks(int count, bool rowsFirst);
  /// Cleans up at end of factorization
  void cleanup();
  /// Sets up all array pointers
  void setupPointers();

  /// Updates part of column (FTRANL)
  void updateColumnL(CoinIndexedVector *region, int *COIN_RESTRICT indexIn) const;
  /// Updates part of column (FTRANL) when densish
  void updateColumnLDensish(CoinIndexedVector *region, int *COIN_RESTRICT indexIn) const;
  /// Updates part of column (FTRANL) when sparse
  void updateColumnLSparse(CoinIndexedVector *region, int *COIN_RESTRICT indexIn) const;
  /// Updates part of column (FTRANL) when sparsish
  void updateColumnLSparsish(CoinIndexedVector *region, int *COIN_RESTRICT indexIn) const;

  /// Updates part of column (FTRANR) without FT update
  void updateColumnR(CoinIndexedVector *region) const;
  /** Updates part of column (FTRANR) with FT update.
      Also stores update after L and R */
  void updateColumnRFT(CoinIndexedVector *region, int *COIN_RESTRICT indexIn);

  /// Updates part of column (FTRANU)
  void updateColumnU(CoinIndexedVector *region, int *COIN_RESTRICT indexIn) const;

  /// Updates part of column (FTRANU) when sparse
  void updateColumnUSparse(CoinIndexedVector *regionSparse,
    int *COIN_RESTRICT indexIn) const;
  /// Updates part of column (FTRANU) when sparsish
  void updateColumnUSparsish(CoinIndexedVector *regionSparse,
    int *COIN_RESTRICT indexIn) const;
  /// Updates part of column (FTRANU)
  int updateColumnUDensish(double *COIN_RESTRICT region,
    int *COIN_RESTRICT regionIndex) const;
  /// Updates part of 2 columns (FTRANU) real work
  void updateTwoColumnsUDensish(
    int &numberNonZero1,
    double *COIN_RESTRICT region1,
    int *COIN_RESTRICT index1,
    int &numberNonZero2,
    double *COIN_RESTRICT region2,
    int *COIN_RESTRICT index2) const;
  /// Updates part of column PFI (FTRAN) (after rest)
  void updateColumnPFI(CoinIndexedVector *regionSparse) const;
  /// Permutes back at end of updateColumn
  void permuteBack(CoinIndexedVector *regionSparse,
    CoinIndexedVector *outVector) const;

  /// Updates part of column transpose PFI (BTRAN) (before rest)
  void updateColumnTransposePFI(CoinIndexedVector *region) const;
  /** Updates part of column transpose (BTRANU),
      assumes index is sorted i.e. region is correct */
  void updateColumnTransposeU(CoinIndexedVector *region,
    int smallestIndex) const;
  /** Updates part of column transpose (BTRANU) when sparsish,
      assumes index is sorted i.e. region is correct */
  void updateColumnTransposeUSparsish(CoinIndexedVector *region,
    int smallestIndex) const;
  /** Updates part of column transpose (BTRANU) when densish,
      assumes index is sorted i.e. region is correct */
  void updateColumnTransposeUDensish(CoinIndexedVector *region,
    int smallestIndex) const;
  /** Updates part of column transpose (BTRANU) when sparse,
      assumes index is sorted i.e. region is correct */
  void updateColumnTransposeUSparse(CoinIndexedVector *region) const;
  /** Updates part of column transpose (BTRANU) by column
      assumes index is sorted i.e. region is correct */
  void updateColumnTransposeUByColumn(CoinIndexedVector *region,
    int smallestIndex) const;

  /// Updates part of column transpose (BTRANR)
  void updateColumnTransposeR(CoinIndexedVector *region) const;
  /// Updates part of column transpose (BTRANR) when dense
  void updateColumnTransposeRDensish(CoinIndexedVector *region) const;
  /// Updates part of column transpose (BTRANR) when sparse
  void updateColumnTransposeRSparse(CoinIndexedVector *region) const;

  /// Updates part of column transpose (BTRANL)
  void updateColumnTransposeL(CoinIndexedVector *region) const;
  /// Updates part of column transpose (BTRANL) when densish by column
  void updateColumnTransposeLDensish(CoinIndexedVector *region) const;
  /// Updates part of column transpose (BTRANL) when densish by row
  void updateColumnTransposeLByRow(CoinIndexedVector *region) const;
  /// Updates part of column transpose (BTRANL) when sparsish by row
  void updateColumnTransposeLSparsish(CoinIndexedVector *region) const;
  /// Updates part of column transpose (BTRANL) when sparse (by Row)
  void updateColumnTransposeLSparse(CoinIndexedVector *region) const;

public:
  /** Replaces one Column to basis for PFI
   returns 0=OK, 1=Probably OK, 2=singular, 3=no room.
   In this case region is not empty - it is incoming variable (updated)
  */
  int replaceColumnPFI(CoinIndexedVector *regionSparse,
    int pivotRow, double alpha);

protected:
  /** Returns accuracy status of replaceColumn
      returns 0=OK, 1=Probably OK, 2=singular */
  int checkPivot(double saveFromU, double oldPivot) const;
  /********************************* START LARGE TEMPLATE ********/
#ifdef INT_IS_8
#define COINFACTORIZATION_BITS_PER_INT 64
#define COINFACTORIZATION_SHIFT_PER_INT 6
#define COINFACTORIZATION_MASK_PER_INT 0x3f
#else
#define COINFACTORIZATION_BITS_PER_INT 32
#define COINFACTORIZATION_SHIFT_PER_INT 5
#define COINFACTORIZATION_MASK_PER_INT 0x1f
#endif
  template < class T >
  inline bool
  pivot(int pivotRow,
    int pivotColumn,
    int pivotRowPosition,
    int pivotColumnPosition,
    CoinFactorizationDouble *COIN_RESTRICT work,
    unsigned int *COIN_RESTRICT workArea2,
    int increment2,
    T *COIN_RESTRICT markRow,
    int largeInteger)
  {
    int *COIN_RESTRICT indexColumnU = indexColumnU_.array();
    int *COIN_RESTRICT startColumnU = startColumnU_.array();
    int *COIN_RESTRICT numberInColumn = numberInColumn_.array();
    CoinFactorizationDouble *COIN_RESTRICT elementU = elementU_.array();
    int *COIN_RESTRICT indexRowU = indexRowU_.array();
    int *COIN_RESTRICT startRowU = startRowU_.array();
    int *COIN_RESTRICT numberInRow = numberInRow_.array();
    CoinFactorizationDouble *COIN_RESTRICT elementL = elementL_.array();
    int *COIN_RESTRICT indexRowL = indexRowL_.array();
    int *COIN_RESTRICT saveColumn = saveColumn_.array();
    int *COIN_RESTRICT nextRow = nextRow_.array();
    int *COIN_RESTRICT lastRow = lastRow_.array();

    //store pivot columns (so can easily compress)
    int numberInPivotRow = numberInRow[pivotRow] - 1;
    int startColumn = startColumnU[pivotColumn];
    int numberInPivotColumn = numberInColumn[pivotColumn] - 1;
    int endColumn = startColumn + numberInPivotColumn + 1;
    int put = 0;
    int startRow = startRowU[pivotRow];
    int endRow = startRow + numberInPivotRow + 1;

    if (pivotColumnPosition < 0) {
      for (pivotColumnPosition = startRow; pivotColumnPosition < endRow; pivotColumnPosition++) {
        int iColumn = indexColumnU[pivotColumnPosition];
        if (iColumn != pivotColumn) {
          saveColumn[put++] = iColumn;
        } else {
          break;
        }
      }
    } else {
      for (int i = startRow; i < pivotColumnPosition; i++) {
        saveColumn[put++] = indexColumnU[i];
      }
    }
    assert(pivotColumnPosition < endRow);
    assert(indexColumnU[pivotColumnPosition] == pivotColumn);
    pivotColumnPosition++;
    for (; pivotColumnPosition < endRow; pivotColumnPosition++) {
      saveColumn[put++] = indexColumnU[pivotColumnPosition];
    }
    //take out this bit of indexColumnU
    int next = nextRow[pivotRow];
    int last = lastRow[pivotRow];

    nextRow[last] = next;
    lastRow[next] = last;
    nextRow[pivotRow] = numberGoodU_; //use for permute
    lastRow[pivotRow] = -2;
    numberInRow[pivotRow] = 0;
    //store column in L, compress in U and take column out
    int l = lengthL_;

    if (l + numberInPivotColumn > lengthAreaL_) {
      //need more memory
      if ((messageLevel_ & 4) != 0)
        printf("more memory needed in middle of invert\n");
      return false;
    }
    //l+=currentAreaL_->elementByColumn-elementL;
    int lSave = l;

    int *COIN_RESTRICT startColumnL = startColumnL_.array();
    startColumnL[numberGoodL_] = l; //for luck and first time
    numberGoodL_++;
    startColumnL[numberGoodL_] = l + numberInPivotColumn;
    lengthL_ += numberInPivotColumn;
    if (pivotRowPosition < 0) {
      for (pivotRowPosition = startColumn; pivotRowPosition < endColumn; pivotRowPosition++) {
        int iRow = indexRowU[pivotRowPosition];
        if (iRow != pivotRow) {
          indexRowL[l] = iRow;
          elementL[l] = elementU[pivotRowPosition];
          markRow[iRow] = static_cast< T >(l - lSave);
          l++;
          //take out of row list
          int start = startRowU[iRow];
          int end = start + numberInRow[iRow];
          int where = start;

          while (indexColumnU[where] != pivotColumn) {
            where++;
          } /* endwhile */
#if DEBUG_COIN
          if (where >= end) {
            abort();
          }
#endif
          indexColumnU[where] = indexColumnU[end - 1];
          numberInRow[iRow]--;
        } else {
          break;
        }
      }
    } else {
      int i;

      for (i = startColumn; i < pivotRowPosition; i++) {
        int iRow = indexRowU[i];

        markRow[iRow] = static_cast< T >(l - lSave);
        indexRowL[l] = iRow;
        elementL[l] = elementU[i];
        l++;
        //take out of row list
        int start = startRowU[iRow];
        int end = start + numberInRow[iRow];
        int where = start;

        while (indexColumnU[where] != pivotColumn) {
          where++;
        } /* endwhile */
#if DEBUG_COIN
        if (where >= end) {
          abort();
        }
#endif
        indexColumnU[where] = indexColumnU[end - 1];
        numberInRow[iRow]--;
        assert(numberInRow[iRow] >= 0);
      }
    }
    assert(pivotRowPosition < endColumn);
    assert(indexRowU[pivotRowPosition] == pivotRow);
    CoinFactorizationDouble pivotElement = elementU[pivotRowPosition];
    CoinFactorizationDouble pivotMultiplier = 1.0 / pivotElement;

    pivotRegion_.array()[numberGoodU_] = pivotMultiplier;
    pivotRowPosition++;
    for (; pivotRowPosition < endColumn; pivotRowPosition++) {
      int iRow = indexRowU[pivotRowPosition];

      markRow[iRow] = static_cast< T >(l - lSave);
      indexRowL[l] = iRow;
      elementL[l] = elementU[pivotRowPosition];
      l++;
      //take out of row list
      int start = startRowU[iRow];
      int end = start + numberInRow[iRow];
      int where = start;

      while (indexColumnU[where] != pivotColumn) {
        where++;
      } /* endwhile */
#if DEBUG_COIN
      if (where >= end) {
        abort();
      }
#endif
      indexColumnU[where] = indexColumnU[end - 1];
      numberInRow[iRow]--;
      assert(numberInRow[iRow] >= 0);
    }
    markRow[pivotRow] = static_cast< T >(largeInteger);
    //compress pivot column (move pivot to front including saved)
    numberInColumn[pivotColumn] = 0;
    //use end of L for temporary space
    int *COIN_RESTRICT indexL = &indexRowL[lSave];
    CoinFactorizationDouble *COIN_RESTRICT multipliersL = &elementL[lSave];

    //adjust
    int j;

    for (j = 0; j < numberInPivotColumn; j++) {
      multipliersL[j] *= pivotMultiplier;
    }
    //zero out fill
    int iErase;
    for (iErase = 0; iErase < increment2 * numberInPivotRow;
         iErase++) {
      workArea2[iErase] = 0;
    }
    int added = numberInPivotRow * numberInPivotColumn;
    unsigned int *COIN_RESTRICT temp2 = workArea2;
    int *COIN_RESTRICT nextColumn = nextColumn_.array();

    //pack down and move to work
    int jColumn;
    for (jColumn = 0; jColumn < numberInPivotRow; jColumn++) {
      int iColumn = saveColumn[jColumn];
      int startColumn = startColumnU[iColumn];
      int endColumn = startColumn + numberInColumn[iColumn];
      int iRow = indexRowU[startColumn];
      CoinFactorizationDouble value = elementU[startColumn];
      double largest;
      int put = startColumn;
      int positionLargest = -1;
      CoinFactorizationDouble thisPivotValue = 0.0;

      //compress column and find largest not updated
      bool checkLargest;
      int mark = markRow[iRow];

      if (mark == largeInteger + 1) {
        largest = fabs(value);
        positionLargest = put;
        put++;
        checkLargest = false;
      } else {
        //need to find largest
        largest = 0.0;
        checkLargest = true;
        if (mark != largeInteger) {
          //will be updated
          work[mark] = value;
          int word = mark >> COINFACTORIZATION_SHIFT_PER_INT;
          int bit = mark & COINFACTORIZATION_MASK_PER_INT;

          temp2[word] = temp2[word] | (1 << bit); //say already in counts
          added--;
        } else {
          thisPivotValue = value;
        }
      }
      int i;
      for (i = startColumn + 1; i < endColumn; i++) {
        iRow = indexRowU[i];
        value = elementU[i];
        int mark = markRow[iRow];

        if (mark == largeInteger + 1) {
          //keep
          indexRowU[put] = iRow;
          elementU[put] = value;
          if (checkLargest) {
            double absValue = fabs(value);

            if (absValue > largest) {
              largest = absValue;
              positionLargest = put;
            }
          }
          put++;
        } else if (mark != largeInteger) {
          //will be updated
          work[mark] = value;
          int word = mark >> COINFACTORIZATION_SHIFT_PER_INT;
          int bit = mark & COINFACTORIZATION_MASK_PER_INT;

          temp2[word] = temp2[word] | (1 << bit); //say already in counts
          added--;
        } else {
          thisPivotValue = value;
        }
      }
      //slot in pivot
      elementU[put] = elementU[startColumn];
      indexRowU[put] = indexRowU[startColumn];
      if (positionLargest == startColumn) {
        positionLargest = put; //follow if was largest
      }
      put++;
      elementU[startColumn] = thisPivotValue;
      indexRowU[startColumn] = pivotRow;
      //clean up counts
      startColumn++;
      numberInColumn[iColumn] = put - startColumn;
      int *COIN_RESTRICT numberInColumnPlus = numberInColumnPlus_.array();
      numberInColumnPlus[iColumn]++;
      startColumnU[iColumn]++;
      //how much space have we got
      int next = nextColumn[iColumn];
      int space;

      space = startColumnU[next] - put - numberInColumnPlus[next];
      //assume no zero elements
      if (numberInPivotColumn > space) {
        //getColumnSpace also moves fixed part
        if (!getColumnSpace(iColumn, numberInPivotColumn)) {
          return false;
        }
        //redo starts
        if (positionLargest >= 0)
          positionLargest = positionLargest + startColumnU[iColumn] - startColumn;
        startColumn = startColumnU[iColumn];
        put = startColumn + numberInColumn[iColumn];
      }
      double tolerance = zeroTolerance_;

      int *COIN_RESTRICT nextCount = nextCount_.array();
      for (j = 0; j < numberInPivotColumn; j++) {
        value = work[j] - thisPivotValue * multipliersL[j];
        double absValue = fabs(value);

        if (absValue > tolerance) {
          work[j] = 0.0;
          assert(put < lengthAreaU_);
          elementU[put] = value;
          indexRowU[put] = indexL[j];
          if (absValue > largest) {
            largest = absValue;
            positionLargest = put;
          }
          put++;
        } else {
          work[j] = 0.0;
          added--;
          int word = j >> COINFACTORIZATION_SHIFT_PER_INT;
          int bit = j & COINFACTORIZATION_MASK_PER_INT;

          if (temp2[word] & (1 << bit)) {
            //take out of row list
            iRow = indexL[j];
            int start = startRowU[iRow];
            int end = start + numberInRow[iRow];
            int where = start;

            while (indexColumnU[where] != iColumn) {
              where++;
            } /* endwhile */
#if DEBUG_COIN
            if (where >= end) {
              abort();
            }
#endif
            indexColumnU[where] = indexColumnU[end - 1];
            numberInRow[iRow]--;
          } else {
            //make sure won't be added
            int word = j >> COINFACTORIZATION_SHIFT_PER_INT;
            int bit = j & COINFACTORIZATION_MASK_PER_INT;

            temp2[word] = temp2[word] | (1 << bit); //say already in counts
          }
        }
      }
      numberInColumn[iColumn] = put - startColumn;
      //move largest
      if (positionLargest >= 0) {
        value = elementU[positionLargest];
        iRow = indexRowU[positionLargest];
        elementU[positionLargest] = elementU[startColumn];
        indexRowU[positionLargest] = indexRowU[startColumn];
        elementU[startColumn] = value;
        indexRowU[startColumn] = iRow;
      }
      //linked list for column
      if (nextCount[iColumn + numberRows_] != -2) {
        //modify linked list
        deleteLink(iColumn + numberRows_);
        addLink(iColumn + numberRows_, numberInColumn[iColumn]);
      }
      temp2 += increment2;
    }
    //get space for row list
    unsigned int *COIN_RESTRICT putBase = workArea2;
    int bigLoops = numberInPivotColumn >> COINFACTORIZATION_SHIFT_PER_INT;
    int i = 0;

    // do linked lists and update counts
    while (bigLoops) {
      bigLoops--;
      int bit;
      for (bit = 0; bit < COINFACTORIZATION_BITS_PER_INT; i++, bit++) {
        unsigned int *COIN_RESTRICT putThis = putBase;
        int iRow = indexL[i];

        //get space
        int number = 0;
        int jColumn;

        for (jColumn = 0; jColumn < numberInPivotRow; jColumn++) {
          unsigned int test = *putThis;

          putThis += increment2;
          test = 1 - ((test >> bit) & 1);
          number += test;
        }
        int next = nextRow[iRow];
        int space;

        space = startRowU[next] - startRowU[iRow];
        number += numberInRow[iRow];
        if (space < number) {
          if (!getRowSpace(iRow, number)) {
            return false;
          }
        }
        // now do
        putThis = putBase;
        next = nextRow[iRow];
        number = numberInRow[iRow];
        int end = startRowU[iRow] + number;
        int saveIndex = indexColumnU[startRowU[next]];

        //add in
        for (jColumn = 0; jColumn < numberInPivotRow; jColumn++) {
          unsigned int test = *putThis;

          putThis += increment2;
          test = 1 - ((test >> bit) & 1);
          indexColumnU[end] = saveColumn[jColumn];
          end += test;
        }
        //put back next one in case zapped
        indexColumnU[startRowU[next]] = saveIndex;
        markRow[iRow] = static_cast< T >(largeInteger + 1);
        number = end - startRowU[iRow];
        numberInRow[iRow] = number;
        deleteLink(iRow);
        addLink(iRow, number);
      }
      putBase++;
    } /* endwhile */
    int bit;

    for (bit = 0; i < numberInPivotColumn; i++, bit++) {
      unsigned int *COIN_RESTRICT putThis = putBase;
      int iRow = indexL[i];

      //get space
      int number = 0;
      int jColumn;

      for (jColumn = 0; jColumn < numberInPivotRow; jColumn++) {
        unsigned int test = *putThis;

        putThis += increment2;
        test = 1 - ((test >> bit) & 1);
        number += test;
      }
      int next = nextRow[iRow];
      int space;

      space = startRowU[next] - startRowU[iRow];
      number += numberInRow[iRow];
      if (space < number) {
        if (!getRowSpace(iRow, number)) {
          return false;
        }
      }
      // now do
      putThis = putBase;
      next = nextRow[iRow];
      number = numberInRow[iRow];
      int end = startRowU[iRow] + number;
      int saveIndex;

      saveIndex = indexColumnU[startRowU[next]];

      //add in
      for (jColumn = 0; jColumn < numberInPivotRow; jColumn++) {
        unsigned int test = *putThis;

        putThis += increment2;
        test = 1 - ((test >> bit) & 1);

        indexColumnU[end] = saveColumn[jColumn];
        end += test;
      }
      indexColumnU[startRowU[next]] = saveIndex;
      markRow[iRow] = static_cast< T >(largeInteger + 1);
      number = end - startRowU[iRow];
      numberInRow[iRow] = number;
      deleteLink(iRow);
      addLink(iRow, number);
    }
    markRow[pivotRow] = static_cast< T >(largeInteger + 1);
    //modify linked list for pivots
    deleteLink(pivotRow);
    deleteLink(pivotColumn + numberRows_);
    totalElements_ += added;
    return true;
  }

  /********************************* END LARGE TEMPLATE ********/
  //@}
  ////////////////// data //////////////////
protected:
  /**@name data */
  //@{
  /// Pivot tolerance
  double pivotTolerance_;
  /// Zero tolerance
  double zeroTolerance_;
#ifndef COIN_FAST_CODE
  /// Whether slack value is  +1 or -1
  double slackValue_;
#else
#ifndef slackValue_
#define slackValue_ -1.0
#endif
#endif
  /// How much to multiply areas by
  double areaFactor_;
  /// Relax check on accuracy in replaceColumn
  double relaxCheck_;
  /// Number of Rows in factorization
  int numberRows_;
  /// Number of Rows after iterating
  int numberRowsExtra_;
  /// Maximum number of Rows after iterating
  int maximumRowsExtra_;
  /// Number of Columns in factorization
  int numberColumns_;
  /// Number of Columns after iterating
  int numberColumnsExtra_;
  /// Maximum number of Columns after iterating
  int maximumColumnsExtra_;
  /// Number factorized in U (not row singletons)
  int numberGoodU_;
  /// Number factorized in L
  int numberGoodL_;
  /// Maximum number of pivots before factorization
  int maximumPivots_;
  /// Number pivots since last factorization
  int numberPivots_;
  /// Number of elements in U (to go)
  ///       or while iterating total overall
  int totalElements_;
  /// Number of elements after factorization
  int factorElements_;
  /// Pivot order for each Column
  CoinIntArrayWithLength pivotColumn_;
  /// Permutation vector for pivot row order
  CoinIntArrayWithLength permute_;
  /// DePermutation vector for pivot row order
  CoinIntArrayWithLength permuteBack_;
  /// Inverse Pivot order for each Column
  CoinIntArrayWithLength pivotColumnBack_;
  /// Status of factorization
  int status_;

  /** 0 - no increasing rows - no permutations,
   1 - no increasing rows but permutations 
   2 - increasing rows 
     - taken out as always 2 */
  //int increasingRows_;

  /// Number of trials before rejection
  int numberTrials_;
  /// Start of each Row as pointer
  CoinIntArrayWithLength startRowU_;

  /// Number in each Row
  CoinIntArrayWithLength numberInRow_;

  /// Number in each Column
  CoinIntArrayWithLength numberInColumn_;

  /// Number in each Column including pivoted
  CoinIntArrayWithLength numberInColumnPlus_;

  /** First Row/Column with count of k,
      can tell which by offset - Rows then Columns */
  CoinIntArrayWithLength firstCount_;

  /// Next Row/Column with count
  CoinIntArrayWithLength nextCount_;

  /// Previous Row/Column with count
  CoinIntArrayWithLength lastCount_;

  /// Next Column in memory order
  CoinIntArrayWithLength nextColumn_;

  /// Previous Column in memory order
  CoinIntArrayWithLength lastColumn_;

  /// Next Row in memory order
  CoinIntArrayWithLength nextRow_;

  /// Previous Row in memory order
  CoinIntArrayWithLength lastRow_;

  /// Columns left to do in a single pivot
  CoinIntArrayWithLength saveColumn_;

  /// Marks rows to be updated
  CoinIntArrayWithLength markRow_;

  /// Detail in messages
  int messageLevel_;

  /// Larger of row and column size
  int biggerDimension_;

  /// Base address for U (may change)
  CoinIntArrayWithLength indexColumnU_;

  /// Pivots for L
  CoinIntArrayWithLength pivotRowL_;

  /// Inverses of pivot values
  CoinFactorizationDoubleArrayWithLength pivotRegion_;

  /// Number of slacks at beginning of U
  int numberSlacks_;

  /// Number in U
  int numberU_;

  /// Maximum space used in U
  int maximumU_;

  /// Base of U is always 0
  //int baseU_;

  /// Length of U
  int lengthU_;

  /// Length of area reserved for U
  int lengthAreaU_;

  /// Elements of U
  CoinFactorizationDoubleArrayWithLength elementU_;

  /// Row indices of U
  CoinIntArrayWithLength indexRowU_;

  /// Start of each column in U
  CoinIntArrayWithLength startColumnU_;

  /// Converts rows to columns in U
  CoinIntArrayWithLength convertRowToColumnU_;

  /// Number in L
  int numberL_;

  /// Base of L
  int baseL_;

  /// Length of L
  int lengthL_;

  /// Length of area reserved for L
  int lengthAreaL_;

  /// Elements of L
  CoinFactorizationDoubleArrayWithLength elementL_;

  /// Row indices of L
  CoinIntArrayWithLength indexRowL_;

  /// Start of each column in L
  CoinIntArrayWithLength startColumnL_;

  /// true if Forrest Tomlin update, false if PFI
  bool doForrestTomlin_;

  /// Number in R
  int numberR_;

  /// Length of R stuff
  int lengthR_;

  /// length of area reserved for R
  int lengthAreaR_;

  /// Elements of R
  CoinFactorizationDouble *elementR_;

  /// Row indices for R
  int *indexRowR_;

  /// Start of columns for R
  CoinIntArrayWithLength startColumnR_;

  /// Dense area
  double *denseArea_;

  /// Dense area - actually used (for alignment etc)
  double *denseAreaAddress_;

  /// Dense permutation
  int *densePermute_;

  /// Number of dense rows
  int numberDense_;

  /// Dense threshold
  int denseThreshold_;

  /// First work area
  CoinFactorizationDoubleArrayWithLength workArea_;

  /// Second work area
  CoinUnsignedIntArrayWithLength workArea2_;

  /// Number of compressions done
  int numberCompressions_;

public:
  /// Below are all to collect
  mutable double ftranCountInput_;
  mutable double ftranCountAfterL_;
  mutable double ftranCountAfterR_;
  mutable double ftranCountAfterU_;
  mutable double btranCountInput_;
  mutable double btranCountAfterU_;
  mutable double btranCountAfterR_;
  mutable double btranCountAfterL_;

  /// We can roll over factorizations
  mutable int numberFtranCounts_;
  mutable int numberBtranCounts_;

  /// While these are average ratios collected over last period
  double ftranAverageAfterL_;
  double ftranAverageAfterR_;
  double ftranAverageAfterU_;
  double btranAverageAfterU_;
  double btranAverageAfterR_;
  double btranAverageAfterL_;

protected:
  /// For statistics
#define collectStatistics_ 1

  /// Below this use sparse technology - if 0 then no L row copy
  int sparseThreshold_;

  /// And one for "sparsish"
  int sparseThreshold2_;

  /// Start of each row in L
  CoinIntArrayWithLength startRowL_;

  /// Index of column in row for L
  CoinIntArrayWithLength indexColumnL_;

  /// Elements in L (row copy)
  CoinFactorizationDoubleArrayWithLength elementByRowL_;

  /// Sparse regions
  mutable CoinIntArrayWithLength sparse_;
#if ABOCA_LITE_FACTORIZATION
  /// Offset to second version of sparse
  int sparseOffset_;
#endif
  /** L to U bias
      0 - U bias, 1 - some U bias, 2 some L bias, 3 L bias
  */
  int biasLU_;
  /** Array persistence flag
      If 0 then as now (delete/new)
      1 then only do arrays if bigger needed
      2 as 1 but give a bit extra if bigger needed
  */
  int persistenceFlag_;
#ifdef ABC_USE_COIN_FACTORIZATION
  /// Says if parallel
  int parallelMode_;
#endif
  /// Arrays - for speed of accessing
  /// Pivot order for each Column
  int * pivotColumnArray_;
  /// Permutation vector for pivot row order
  int * permuteArray_;
  /// DePermutation vector for pivot row order
  int * permuteBackArray_;
  /// Inverse Pivot order for each Column
  int * pivotColumnBackArray_;
  /// Start of each Row as pointer
  int * startRowUArray_;
  /// Number in each Row
  int * numberInRowArray_;
  /// Number in each Column
  int * numberInColumnArray_;
  /// Number in each Column including pivoted
  int * numberInColumnPlusArray_;
  /** First Row/Column with count of k,
      can tell which by offset - Rows then Columns */
  int * firstCountArray_;
  /// Next Row/Column with count
  int * nextCountArray_;
  /// Previous Row/Column with count
  int * lastCountArray_;
  /// Next Column in memory order
  int * nextColumnArray_;
  /// Previous Column in memory order
  int * lastColumnArray_;
  /// Next Row in memory order
  int * nextRowArray_;
  /// Previous Row in memory order
  int * lastRowArray_;
  /// Columns left to do in a single pivot
  int * saveColumnArray_;
  /// Marks rows to be updated
  int * markRowArray_;
  /// Base address for U (may change)
  int * indexColumnUArray_;
  /// Pivots for L
  int * pivotRowLArray_;
  /// Inverses of pivot values
  CoinFactorizationDouble * pivotRegionArray_;
  /// Elements of U
  CoinFactorizationDouble * elementUArray_;
  /// Row indices of U
  int * indexRowUArray_;
  /// Start of each column in U
  int * startColumnUArray_;
  /// Converts rows to columns in U
  int * convertRowToColumnUArray_;
  /// Elements of L
  CoinFactorizationDouble * elementLArray_;
  /// Row indices of L
  int * indexRowLArray_;
  /// Start of each column in L
  int * startColumnLArray_;
  /// Start of columns for R
  int * startColumnRArray_;
  /// First work area
  CoinFactorizationDouble * workAreaArray_;
  /// Second work area
  unsigned int * workArea2Array_;
  /// Start of each row in L
  int * startRowLArray_;
  /// Index of column in row for L
  int * indexColumnLArray_;
  /// Elements in L (row copy)
  CoinFactorizationDouble * elementByRowLArray_;
  /// Sparse regions
  mutable int * sparseArray_;
  //@}
};
// Dense coding
#ifdef INTEL_COMPILER
#define COIN_FACTORIZATION_DENSE_CODE 3
#endif
#ifdef COINUTILS_HAS_LAPACK
#ifndef COIN_FACTORIZATION_DENSE_CODE
#define COIN_FACTORIZATION_DENSE_CODE 1
#endif
#endif
#ifdef COIN_FACTORIZATION_DENSE_CODE
/* Type of Fortran integer translated into C */
#ifndef ipfint
//typedef ipfint FORTRAN_INTEGER_TYPE ;
typedef int ipfint;
typedef const int cipfint;
#endif
#endif
#endif
#ifdef CLP_NO_DENSE
#undef COIN_FACTORIZATION_DENSE_CODE
#endif
// Extra for ugly include
#ifdef UGLY_COIN_FACTOR_CODING
#define FAC_UNSET (FAC_SET + 1)
{
  goodPivot = false;
  //store pivot columns (so can easily compress)
  int startColumnThis = startColumn[iPivotColumn];
  int endColumn = startColumnThis + numberDoColumn + 1;
  int put = 0;
  int startRowThis = startRow[iPivotRow];
  int endRow = startRowThis + numberDoRow + 1;
  if (pivotColumnPosition < 0) {
    for (pivotColumnPosition = startRowThis; pivotColumnPosition < endRow; pivotColumnPosition++) {
      int iColumn = indexColumn[pivotColumnPosition];
      if (iColumn != iPivotColumn) {
        saveColumn[put++] = iColumn;
      } else {
        break;
      }
    }
  } else {
    for (int i = startRowThis; i < pivotColumnPosition; i++) {
      saveColumn[put++] = indexColumn[i];
    }
  }
  assert(pivotColumnPosition < endRow);
  assert(indexColumn[pivotColumnPosition] == iPivotColumn);
  pivotColumnPosition++;
  for (; pivotColumnPosition < endRow; pivotColumnPosition++) {
    saveColumn[put++] = indexColumn[pivotColumnPosition];
  }
  //take out this bit of indexColumn
  int next = nextRow[iPivotRow];
  int last = lastRow[iPivotRow];

  nextRow[last] = next;
  lastRow[next] = last;
  nextRow[iPivotRow] = numberGoodU_; //use for permute
  lastRow[iPivotRow] = -2;
  numberInRow[iPivotRow] = 0;
  //store column in L, compress in U and take column out
  int l = lengthL_;
  // **** HORRID coding coming up but a goto seems best!
  {
    if (l + numberDoColumn > lengthAreaL_) {
      //need more memory
      if ((messageLevel_ & 4) != 0)
        printf("more memory needed in middle of invert\n");
      goto BAD_PIVOT;
    }
    //l+=currentAreaL_->elementByColumn-elementL;
    int lSave = l;

    int *startColumnL = startColumnL_.array();
    startColumnL[numberGoodL_] = l; //for luck and first time
    numberGoodL_++;
    startColumnL[numberGoodL_] = l + numberDoColumn;
    lengthL_ += numberDoColumn;
    if (pivotRowPosition < 0) {
      for (pivotRowPosition = startColumnThis; pivotRowPosition < endColumn; pivotRowPosition++) {
        int iRow = indexRow[pivotRowPosition];
        if (iRow != iPivotRow) {
          indexRowL[l] = iRow;
          elementL[l] = element[pivotRowPosition];
          markRow[iRow] = l - lSave;
          l++;
          //take out of row list
          int start = startRow[iRow];
          int end = start + numberInRow[iRow];
          int where = start;

          while (indexColumn[where] != iPivotColumn) {
            where++;
          } /* endwhile */
#if DEBUG_COIN
          if (where >= end) {
            abort();
          }
#endif
          indexColumn[where] = indexColumn[end - 1];
          numberInRow[iRow]--;
        } else {
          break;
        }
      }
    } else {
      int i;

      for (i = startColumnThis; i < pivotRowPosition; i++) {
        int iRow = indexRow[i];

        markRow[iRow] = l - lSave;
        indexRowL[l] = iRow;
        elementL[l] = element[i];
        l++;
        //take out of row list
        int start = startRow[iRow];
        int end = start + numberInRow[iRow];
        int where = start;

        while (indexColumn[where] != iPivotColumn) {
          where++;
        } /* endwhile */
#if DEBUG_COIN
        if (where >= end) {
          abort();
        }
#endif
        indexColumn[where] = indexColumn[end - 1];
        numberInRow[iRow]--;
        assert(numberInRow[iRow] >= 0);
      }
    }
    assert(pivotRowPosition < endColumn);
    assert(indexRow[pivotRowPosition] == iPivotRow);
    CoinFactorizationDouble pivotElement = element[pivotRowPosition];
    CoinFactorizationDouble pivotMultiplier = 1.0 / pivotElement;

    pivotRegion_.array()[numberGoodU_] = pivotMultiplier;
    pivotRowPosition++;
    for (; pivotRowPosition < endColumn; pivotRowPosition++) {
      int iRow = indexRow[pivotRowPosition];

      markRow[iRow] = l - lSave;
      indexRowL[l] = iRow;
      elementL[l] = element[pivotRowPosition];
      l++;
      //take out of row list
      int start = startRow[iRow];
      int end = start + numberInRow[iRow];
      int where = start;

      while (indexColumn[where] != iPivotColumn) {
        where++;
      } /* endwhile */
#if DEBUG_COIN
      if (where >= end) {
        abort();
      }
#endif
      indexColumn[where] = indexColumn[end - 1];
      numberInRow[iRow]--;
      assert(numberInRow[iRow] >= 0);
    }
    markRow[iPivotRow] = FAC_SET;
    //compress pivot column (move pivot to front including saved)
    numberInColumn[iPivotColumn] = 0;
    //use end of L for temporary space
    int *indexL = &indexRowL[lSave];
    CoinFactorizationDouble *COIN_RESTRICT multipliersL = &elementL[lSave];

    //adjust
    int j;

    for (j = 0; j < numberDoColumn; j++) {
      multipliersL[j] *= pivotMultiplier;
    }
    //zero out fill
    int iErase;
    for (iErase = 0; iErase < increment2 * numberDoRow;
         iErase++) {
      workArea2[iErase] = 0;
    }
    int added = numberDoRow * numberDoColumn;
    unsigned int *temp2 = workArea2;
    int *nextColumn = nextColumn_.array();

    //pack down and move to work
    int jColumn;
    for (jColumn = 0; jColumn < numberDoRow; jColumn++) {
      int iColumn = saveColumn[jColumn];
      int startColumnThis = startColumn[iColumn];
      int endColumn = startColumnThis + numberInColumn[iColumn];
      int iRow = indexRow[startColumnThis];
      CoinFactorizationDouble value = element[startColumnThis];
      double largest;
      int put = startColumnThis;
      int positionLargest = -1;
      CoinFactorizationDouble thisPivotValue = 0.0;

      //compress column and find largest not updated
      bool checkLargest;
      int mark = markRow[iRow];

      if (mark == FAC_UNSET) {
        largest = fabs(value);
        positionLargest = put;
        put++;
        checkLargest = false;
      } else {
        //need to find largest
        largest = 0.0;
        checkLargest = true;
        if (mark != FAC_SET) {
          //will be updated
          workArea[mark] = value;
          int word = mark >> COINFACTORIZATION_SHIFT_PER_INT;
          int bit = mark & COINFACTORIZATION_MASK_PER_INT;

          temp2[word] = temp2[word] | (1 << bit); //say already in counts
          added--;
        } else {
          thisPivotValue = value;
        }
      }
      int i;
      for (i = startColumnThis + 1; i < endColumn; i++) {
        iRow = indexRow[i];
        value = element[i];
        int mark = markRow[iRow];

        if (mark == FAC_UNSET) {
          //keep
          indexRow[put] = iRow;
          element[put] = value;
          if (checkLargest) {
            double absValue = fabs(value);

            if (absValue > largest) {
              largest = absValue;
              positionLargest = put;
            }
          }
          put++;
        } else if (mark != FAC_SET) {
          //will be updated
          workArea[mark] = value;
          int word = mark >> COINFACTORIZATION_SHIFT_PER_INT;
          int bit = mark & COINFACTORIZATION_MASK_PER_INT;

          temp2[word] = temp2[word] | (1 << bit); //say already in counts
          added--;
        } else {
          thisPivotValue = value;
        }
      }
      //slot in pivot
      element[put] = element[startColumnThis];
      indexRow[put] = indexRow[startColumnThis];
      if (positionLargest == startColumnThis) {
        positionLargest = put; //follow if was largest
      }
      put++;
      element[startColumnThis] = thisPivotValue;
      indexRow[startColumnThis] = iPivotRow;
      //clean up counts
      startColumnThis++;
      numberInColumn[iColumn] = put - startColumnThis;
      int *numberInColumnPlus = numberInColumnPlus_.array();
      numberInColumnPlus[iColumn]++;
      startColumn[iColumn]++;
      //how much space have we got
      int next = nextColumn[iColumn];
      int space;

      space = startColumn[next] - put - numberInColumnPlus[next];
      //assume no zero elements
      if (numberDoColumn > space) {
        //getColumnSpace also moves fixed part
        if (!getColumnSpace(iColumn, numberDoColumn)) {
          goto BAD_PIVOT;
        }
        //redo starts
        positionLargest = positionLargest + startColumn[iColumn] - startColumnThis;
        startColumnThis = startColumn[iColumn];
        put = startColumnThis + numberInColumn[iColumn];
      }
      double tolerance = zeroTolerance_;

      int *nextCount = nextCount_.array();
      for (j = 0; j < numberDoColumn; j++) {
        value = workArea[j] - thisPivotValue * multipliersL[j];
        double absValue = fabs(value);

        if (absValue > tolerance) {
          workArea[j] = 0.0;
          element[put] = value;
          indexRow[put] = indexL[j];
          if (absValue > largest) {
            largest = absValue;
            positionLargest = put;
          }
          put++;
        } else {
          workArea[j] = 0.0;
          added--;
          int word = j >> COINFACTORIZATION_SHIFT_PER_INT;
          int bit = j & COINFACTORIZATION_MASK_PER_INT;

          if (temp2[word] & (1 << bit)) {
            //take out of row list
            iRow = indexL[j];
            int start = startRow[iRow];
            int end = start + numberInRow[iRow];
            int where = start;

            while (indexColumn[where] != iColumn) {
              where++;
            } /* endwhile */
#if DEBUG_COIN
            if (where >= end) {
              abort();
            }
#endif
            indexColumn[where] = indexColumn[end - 1];
            numberInRow[iRow]--;
          } else {
            //make sure won't be added
            int word = j >> COINFACTORIZATION_SHIFT_PER_INT;
            int bit = j & COINFACTORIZATION_MASK_PER_INT;

            temp2[word] = temp2[word] | (1 << bit); //say already in counts
          }
        }
      }
      numberInColumn[iColumn] = put - startColumnThis;
      //move largest
      if (positionLargest >= 0) {
        value = element[positionLargest];
        iRow = indexRow[positionLargest];
        element[positionLargest] = element[startColumnThis];
        indexRow[positionLargest] = indexRow[startColumnThis];
        element[startColumnThis] = value;
        indexRow[startColumnThis] = iRow;
      }
      //linked list for column
      if (nextCount[iColumn + numberRows_] != -2) {
        //modify linked list
        deleteLink(iColumn + numberRows_);
        addLink(iColumn + numberRows_, numberInColumn[iColumn]);
      }
      temp2 += increment2;
    }
    //get space for row list
    unsigned int *putBase = workArea2;
    int bigLoops = numberDoColumn >> COINFACTORIZATION_SHIFT_PER_INT;
    int i = 0;

    // do linked lists and update counts
    while (bigLoops) {
      bigLoops--;
      int bit;
      for (bit = 0; bit < COINFACTORIZATION_BITS_PER_INT; i++, bit++) {
        unsigned int *putThis = putBase;
        int iRow = indexL[i];

        //get space
        int number = 0;
        int jColumn;

        for (jColumn = 0; jColumn < numberDoRow; jColumn++) {
          unsigned int test = *putThis;

          putThis += increment2;
          test = 1 - ((test >> bit) & 1);
          number += test;
        }
        int next = nextRow[iRow];
        int space;

        space = startRow[next] - startRow[iRow];
        number += numberInRow[iRow];
        if (space < number) {
          if (!getRowSpace(iRow, number)) {
            goto BAD_PIVOT;
          }
        }
        // now do
        putThis = putBase;
        next = nextRow[iRow];
        number = numberInRow[iRow];
        int end = startRow[iRow] + number;
        int saveIndex = indexColumn[startRow[next]];

        //add in
        for (jColumn = 0; jColumn < numberDoRow; jColumn++) {
          unsigned int test = *putThis;

          putThis += increment2;
          test = 1 - ((test >> bit) & 1);
          indexColumn[end] = saveColumn[jColumn];
          end += test;
        }
        //put back next one in case zapped
        indexColumn[startRow[next]] = saveIndex;
        markRow[iRow] = FAC_UNSET;
        number = end - startRow[iRow];
        numberInRow[iRow] = number;
        deleteLink(iRow);
        addLink(iRow, number);
      }
      putBase++;
    } /* endwhile */
    int bit;

    for (bit = 0; i < numberDoColumn; i++, bit++) {
      unsigned int *putThis = putBase;
      int iRow = indexL[i];

      //get space
      int number = 0;
      int jColumn;

      for (jColumn = 0; jColumn < numberDoRow; jColumn++) {
        unsigned int test = *putThis;

        putThis += increment2;
        test = 1 - ((test >> bit) & 1);
        number += test;
      }
      int next = nextRow[iRow];
      int space;

      space = startRow[next] - startRow[iRow];
      number += numberInRow[iRow];
      if (space < number) {
        if (!getRowSpace(iRow, number)) {
          goto BAD_PIVOT;
        }
      }
      // now do
      putThis = putBase;
      next = nextRow[iRow];
      number = numberInRow[iRow];
      int end = startRow[iRow] + number;
      int saveIndex;

      saveIndex = indexColumn[startRow[next]];

      //add in
      for (jColumn = 0; jColumn < numberDoRow; jColumn++) {
        unsigned int test = *putThis;

        putThis += increment2;
        test = 1 - ((test >> bit) & 1);

        indexColumn[end] = saveColumn[jColumn];
        end += test;
      }
      indexColumn[startRow[next]] = saveIndex;
      markRow[iRow] = FAC_UNSET;
      number = end - startRow[iRow];
      numberInRow[iRow] = number;
      deleteLink(iRow);
      addLink(iRow, number);
    }
    markRow[iPivotRow] = FAC_UNSET;
    //modify linked list for pivots
    deleteLink(iPivotRow);
    deleteLink(iPivotColumn + numberRows_);
    totalElements_ += added;
    goodPivot = true;
    // **** UGLY UGLY UGLY
  }
BAD_PIVOT:

  ;
}
#undef FAC_UNSET
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
