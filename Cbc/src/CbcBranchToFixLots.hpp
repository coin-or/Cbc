// $Id$
// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/13/2009-- carved out of CbcBranchCut

#ifndef CbcBranchToFixLots_H
#define CbcBranchToFixLots_H

#include "CbcBranchCut.hpp"
#include "CbcBranchBase.hpp"
#include "OsiRowCut.hpp"
#include "CoinPackedMatrix.hpp"

/** Define a branch class that branches so that one way variables are fixed
    while the other way cuts off that solution.
    a) On reduced cost
    b) When enough ==1 or <=1 rows have been satisfied (not fixed - satisfied)
*/

class CbcBranchToFixLots : public CbcBranchCut {

public:
  // Default Constructor
  CbcBranchToFixLots();

  /** Useful constructor - passed reduced cost tolerance and fraction we would like fixed.
        Also depth level to do at.
        Also passed number of 1 rows which when clean triggers fix
        Always does if all 1 rows cleaned up and number>0 or if fraction columns reached
        Also whether to create branch if can't reach fraction.
    */
  CbcBranchToFixLots(CbcModel *model, double djTolerance,
    double fractionFixed, int depth,
    int numberClean = 0,
    const char *mark = NULL,
    bool alwaysCreate = false);

  // Copy constructor
  CbcBranchToFixLots(const CbcBranchToFixLots &);

  /// Clone
  virtual CbcObject *clone() const;

  // Assignment operator
  CbcBranchToFixLots &operator=(const CbcBranchToFixLots &rhs);

  // Destructor
  ~CbcBranchToFixLots();

  /** Does a lot of the work,
        Returns 0 if no good, 1 if dj, 2 if clean, 3 if both
	FIXME: should use enum or equivalent to make these numbers clearer.
    */
  int shallWe() const;

  /// Infeasibility for an integer variable - large is 0.5, but also can be infinity when known infeasible.
  virtual double infeasibility(const OsiBranchingInformation *info,
    int &preferredWay) const;
  /** \brief Return true if object can take part in normal heuristics
    */
  virtual bool canDoHeuristics() const
  {
    return true;
  }

  /// Creates a branching object
  virtual CbcBranchingObject *createCbcBranch(OsiSolverInterface *solver, const OsiBranchingInformation *info, int way);
  /// Redoes data when sequence numbers change
  virtual void redoSequenceEtc(CbcModel *model, int numberColumns, const int *originalColumns);

protected:
  /// data

  /// Reduced cost tolerance i.e. dj has to be >= this before fixed
  double djTolerance_;
  /// We only need to make sure this fraction fixed
  double fractionFixed_;
  /// Never fix ones marked here
  char *mark_;
  /// Matrix by row
  CoinPackedMatrix matrixByRow_;
  /// Do if depth multiple of this
  int depth_;
  /// number of ==1 rows which need to be clean
  int numberClean_;
  /// If true then always create branch
  bool alwaysCreate_;
};
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
