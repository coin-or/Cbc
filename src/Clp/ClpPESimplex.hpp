// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
/*
   Authors

   Jeremy Omer

   Last update: june 10, 2014

 */

#ifndef ClpPESimplex_H
#define ClpPESimplex_H

#include "ClpSimplex.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpPackedMatrix.hpp"
#include "CoinTime.hpp"

#include <iostream>
#include <fstream>

// #define PE_DEBUG 0

/** SHARED METHODS FOR USEFUL ALGEBRAIC OPERATIONS */

/** inner product between a coin vector and a pointer */
CLPLIB_EXPORT
double PEdot(CoinIndexedVector &v1, const double *v2);

/** inner product between two coin vectors
    call the function with the sparser vector first for efficiency */
CLPLIB_EXPORT
double PEdot(CoinIndexedVector &v1, CoinIndexedVector &v2);

/** compute the product x^T*[A I] for the indices "which" of [A I] */
CLPLIB_EXPORT
void PEtransposeTimesSubsetAll(ClpSimplex *model, int number, const int *which,
  const double *COIN_RESTRICT x, double *COIN_RESTRICT y,
  const double *COIN_RESTRICT rowScale,
  const double *COIN_RESTRICT columnScale);

/** BASE CLASS FOR THE IMPROVED SIMPLEX
*/
class CLPLIB_EXPORT ClpPESimplex {

public:
  /** Constructor */
  ClpPESimplex(ClpSimplex *model);

  /** Destructor */
  ~ClpPESimplex();

  /** BASIC GET METHODS */
public:
  inline int coPrimalDegenerates() { return coPrimalDegenerates_; }
  inline int coDualDegenerates() { return coDualDegenerates_; }
  inline int coCompatibleCols() { return coCompatibleCols_; }
  inline int coCompatibleRows() { return coCompatibleRows_; }

  inline bool isCompatibleCol(int sequence) { return isCompatibleCol_[sequence]; }
  inline bool isCompatibleRow(int row)
  {
    assert(row >= 0 && row < numberRows_);
    return isCompatibleRow_[row];
  }

  inline ClpSimplex *clpModel() { return model_; }
  // check seems to be same model - returns false if size changed
  bool checkSize();
  /** PUBLIC METHODS RELATED TO COMPATIBILITY */
public:
  /** Updates the set of variables that are not at their bounds */
  void updatePrimalDegenerates();

  /** Updates the set of dual degenerate variables */
  void updateDualDegenerates();

  /** Identify the primal compatible columns
        The input argument is a temporary array that is needed for the Clp's BTRAN */
  void identifyCompatibleCols(int number, const int *which,
    CoinIndexedVector *spareRow2,
    CoinIndexedVector *wPrimal);

  /** Identify the dual compatible rows */
  void identifyCompatibleRows(CoinIndexedVector *spare,
    CoinIndexedVector *wDual);

  /** Update the dual compatible rows */
  void updateCompatibleRows(int sequence);

  /** DEBUG AND DISPLAY METHODS */
public:
#if PE_DEBUG >= 1
  /** Print the set of variables within their bounds */
  void printPrimalDegenerates();

  /** Print the set of primal compatible variables */
  void printCompatibleCols();

  /** Check that a nonbasic variable is indeed compatible */
  bool checkCompatibilityCol(int sequence, CoinIndexedVector *spareRow2);
#endif

  /** Check that a basic row is indeed compatible */
  bool checkCompatibilityRow(int pivotRow);

  /** Tracking the degenerate iterations after compatible pivots */
  inline double lastObjectiveValue() { return lastObjectiveValue_; }
  inline void updateLastObjectiveValue() { lastObjectiveValue_ = model_->objectiveValue(); }
  inline bool isDegeneratePivot() { return fabs(model_->objectiveValue() - lastObjectiveValue_) < model_->dualTolerance(); }
  inline bool isLastPivotCompatible() { return isLastPivotCompatible_; }
  inline void isLastPivotCompatible(bool yesOrNo) { isLastPivotCompatible_ = yesOrNo; }

  /** Start and stop the timer, and print the total recorded time */
  inline void startTimer() { timeTmp_ = CoinCpuTime(); }
  inline void stopTimer() { timeCompatibility_ += CoinCpuTime() - timeTmp_; }
  void printTimer(std::ostream &out);
  inline double timeMultRandom() { return timeMultRandom_; }
  inline double timeLinearSystem() { return timeLinearSystem_; }
  inline double timeCompatibility() { return timeCompatibility_; }

  /** Update and return the number of degenerate pivots and variables */
  inline void addDegeneratePivot() { coDegeneratePivots_++; }
  inline int coDegeneratePivots() { return coDegeneratePivots_; }
  inline void addDegeneratePivotConsecutive() { coDegeneratePivotsConsecutive_++; }
  inline void resetDegeneratePivotsConsecutive() { coDegeneratePivotsConsecutive_ = 0; }
  inline int coDegeneratePivotsConsecutive() { return coDegeneratePivotsConsecutive_; }
  void updateDualDegeneratesAvg(int coPivots);
  inline double coDualDegeneratesAvg() { return coDualDegeneratesAvg_; }
  void updatePrimalDegeneratesAvg(int coPivots);
  inline double coPrimalDegeneratesAvg() { return coPrimalDegeneratesAvg_; }
  inline double coCompatibleRowsAvg() { return coCompatibleRowsAvg_; }
  void updateCompatibleRowsAvg(int coPivots);
  inline double coCompatibleColsAvg() { return coCompatibleColsAvg_; }
  void updateCompatibleColsAvg(int coPivots);
  inline int coCompatiblePivots() { return coCompatiblePivots_; }
  inline void addCompatiblePivot() { coCompatiblePivots_++; }
  inline int coDegenerateCompatiblePivots() { return coDegenerateCompatiblePivots_; }
  inline void addDegenerateCompatiblePivot() { coDegenerateCompatiblePivots_++; }

  /* Get and update the number of compatible pivots that were done because of the priority factor */
  inline void addPriorityPivot() { coPriorityPivots_++; }
  inline int coPriorityPivots() { return coPriorityPivots_; }
  inline int doStatistics() const
  {
    return doStatistics_;
  }
  inline void setDoStatistics(int value)
  {
    doStatistics_ = value;
  }

protected:
  /** Indices of the variables that were not at one of their bounds
        during the last update (non primal degenerate variables) */
  int coPrimalDegenerates_;
  int *primalDegenerates_;
  bool *isPrimalDegenerate_;

  /** Indices of the non basic variables with a zero reduced cost
        during the last update (ndual-degenerate variables) */
  int coDualDegenerates_;
  int *dualDegenerates_;
  bool *isDualDegenerate_;

  /** Table of booleans indicating whether each variable is primal
        compatible (true) or not (false) */
  int coCompatibleCols_;
  double *compatibilityCol_;
  bool *isCompatibleCol_;

  /** Table of booleans indicating whether each constraint is dual
        compatible (true) or not (false) */
  int coCompatibleRows_;
  double *compatibilityRow_;
  bool *isCompatibleRow_;

private:
  /** pointer to the original model that shall be solved */
  ClpSimplex *model_;

  /** tolerance used for the tests of degeneracy and compatibility (resp.) */
  double epsDegeneracy_;
  double epsCompatibility_;

  /** size of the original model */
  int numberRows_;
  int numberColumns_;

  /** w vectors that are used to identify the compatible columns and
    rows. The name w, refers to the notations of the articles on
    positive edge */
  // now passed in CoinIndexedVector *wPrimal_;
  // now passed in CoinIndexedVector *wDual_;

  /** temporary vectors that are used to store colulns of the constraint
        matrix or random numbers */
  // not usedCoinIndexedVector *tempColumn_;
  double *tempRandom_;

  /** number of degenerate pivots and variables */
  int coPrimalDegeneratesAvg_;
  int coDualDegeneratesAvg_;
  int coCompatibleColsAvg_;
  int coCompatibleRowsAvg_;
  int coUpdateDegenerates_;
  //int coIdentifyCompatibles_;
  int coDegeneratePivots_;
  int coCompatiblePivots_;
  int coDegenerateCompatiblePivots_;
  int coDegeneratePivotsConsecutive_;

  /** number of compatible pivots that were done because of the priority factor */
  int coPriorityPivots_;
  /// Do statistics
  int doStatistics_;

  /** tracking the degenerate iterations after compatible pivots */
  double lastObjectiveValue_;
  bool isLastPivotCompatible_;

  /** Timer attribute recording the additional time spent in
        identifying compatible variables */
  double timeCompatibility_;
  double timeMultRandom_;
  double timeLinearSystem_;
  double timeTmp_;
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
