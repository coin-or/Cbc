/* $Id$ */
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcFathomDynamicProgramming_H
#define CbcFathomDynamicProgramming_H

#include "CbcFathom.hpp"

//#############################################################################
/** FathomDynamicProgramming class.

    The idea is that after some branching the problem will be effectively smaller than
    the original problem and maybe there will be a more specialized technique which can completely
    fathom this branch quickly.

    This is a dynamic programming implementation which is very fast for some
    specialized problems.  It expects small integral rhs, an all integer problem
    and positive integral coefficients. At present it can not do general set covering
    problems just set partitioning.  It can find multiple optima for various rhs
    combinations.

    The main limiting factor is size of state space.  Each 1 rhs doubles the size of the problem.
    2 or 3 rhs quadruples, 4,5,6,7 by 8 etc.
 */

class CbcFathomDynamicProgramming : public CbcFathom {
public:
  // Default Constructor
  CbcFathomDynamicProgramming();

  // Constructor with model - assumed before cuts
  CbcFathomDynamicProgramming(CbcModel &model);
  // Copy constructor
  CbcFathomDynamicProgramming(const CbcFathomDynamicProgramming &rhs);

  virtual ~CbcFathomDynamicProgramming();

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel *model);

  /// Clone
  virtual CbcFathom *clone() const;

  /// Resets stuff if model changes
  virtual void resetModel(CbcModel *model);

  /** returns 0 if no fathoming attempted, 1 fully fathomed ,
        2 incomplete search, 3 incomplete search but treat as complete.
        If solution then newSolution will not be NULL and
        will be freed by CbcModel.  It is expected that the solution is better
        than best so far but CbcModel will double check.

        If returns 3 then of course there is no guarantee of global optimum
    */
  virtual int fathom(double *&newSolution);

  /// Maximum size allowed
  inline int maximumSize() const
  {
    return maximumSizeAllowed_;
  }
  inline void setMaximumSize(int value)
  {
    maximumSizeAllowed_ = value;
  }
  /// Returns type of algorithm and sets up arrays
  int checkPossible(int allowableSize = 0);
  // set algorithm
  inline void setAlgorithm(int value)
  {
    algorithm_ = value;
  }
  /** Tries a column
        returns true if was used in making any changes.
    */
  bool tryColumn(int numberElements, const int *rows,
    const double *coefficients, double cost,
    int upper = COIN_INT_MAX);
  /// Returns cost array
  inline const double *cost() const
  {
    return cost_;
  }
  /// Returns back array
  inline const int *back() const
  {
    return back_;
  }
  /// Gets bit pattern for target result
  inline int target() const
  {
    return target_;
  }
  /// Sets bit pattern for target result
  inline void setTarget(int value)
  {
    target_ = value;
  }

private:
  /// Does deleteions
  void gutsOfDelete();

  /** Adds one attempt of one column of type 0,
        returns true if was used in making any changes
    */
  bool addOneColumn0(int numberElements, const int *rows,
    double cost);
  /** Adds one attempt of one column of type 1,
        returns true if was used in making any changes.
        At present the user has to call it once for each possible value
    */
  bool addOneColumn1(int numberElements, const int *rows,
    const int *coefficients, double cost);
  /** Adds one attempt of one column of type 1,
        returns true if was used in making any changes.
        At present the user has to call it once for each possible value.
        This version is when there are enough 1 rhs to do faster
    */
  bool addOneColumn1A(int numberElements, const int *rows,
    const int *coefficients, double cost);
  /// Gets bit pattern from original column
  int bitPattern(int numberElements, const int *rows,
    const int *coefficients);
  /// Gets bit pattern from original column
  int bitPattern(int numberElements, const int *rows,
    const double *coefficients);
  /// Fills in original column (dense) from bit pattern - returning number nonzero
  int decodeBitPattern(int bitPattern, int *values, int numberRows);

protected:
  /// Size of states (power of 2 unless just one constraint)
  int size_;
  /** Type - 0 coefficients and rhs all 1,
        1 - coefficients > 1 or rhs > 1
    */
  int type_;
  /// Space for states
  double *cost_;
  /// Which state produced this cheapest one
  int *back_;
  /// Some rows may be satisified so we need a lookup
  int *lookup_;
  /// Space for sorted indices
  int *indices_;
  /// Number of active rows
  int numberActive_;
  /// Maximum size allowed
  int maximumSizeAllowed_;
  /// Start bit for each active row
  int *startBit_;
  /// Number bits for each active row
  int *numberBits_;
  /// Effective rhs
  int *rhs_;
  /// Space for sorted coefficients
  int *coefficients_;
  /// Target pattern
  int target_;
  /// Number of Non 1 rhs
  int numberNonOne_;
  /// Current bit pattern
  int bitPattern_;
  /// Current algorithm
  int algorithm_;

private:
  /// Illegal Assignment operator
  CbcFathomDynamicProgramming &operator=(const CbcFathomDynamicProgramming &rhs);
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
