// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/10/2009-- carved out of CbcBranchActual

#ifndef CbcSubProblem_H
#define CbcSubProblem_H

#ifdef COIN_HAS_CLP
#include "ClpSimplex.hpp"
#include "ClpNode.hpp"

/** Defines a general subproblem
    Basis will be made more compact later
*/
class CoinWarmStartDiff;
class CbcSubProblem {

public:
  /// Default constructor
  CbcSubProblem();

  /// Constructor from model
  CbcSubProblem(const OsiSolverInterface *solver,
    const double *lowerBefore,
    const double *upperBefore,
    const unsigned char *status,
    int depth);

  /// Copy constructor
  CbcSubProblem(const CbcSubProblem &);

  /// Assignment operator
  CbcSubProblem &operator=(const CbcSubProblem &rhs);

  /// Destructor
  virtual ~CbcSubProblem();

  /// Take over
  void takeOver(CbcSubProblem &, bool cleanup);
  /// Apply subproblem (1=bounds, 2=basis, 3=both)
  void apply(OsiSolverInterface *model, int what = 3) const;

public:
  /// Value of objective
  double objectiveValue_;
  /// Sum of infeasibilities
  double sumInfeasibilities_;
  /// Branch value
  double branchValue_;
  /// Dj on branching variable at end
  double djValue_;
  /** Which variable (top bit if upper bound changing)
        next bit if changing on down branch only */
  int *variables_;
  /// New bound
  double *newBounds_;
  /// Status
  mutable CoinWarmStartBasis *status_;
  /// Depth
  int depth_;
  /// Number of Extra bound changes
  int numberChangedBounds_;
  /// Number of infeasibilities
  int numberInfeasibilities_;
  /** Status 1 bit going up on first, 2 bit set first branch infeasible on second, 4 bit redundant branch,
	bits after 256 give reason for stopping (just last node)
	0 - solution
	1 - infeasible
	2 - maximum depth
	>2 - error or max time or something
    */
  int problemStatus_;
  /// Variable branched on
  int branchVariable_;
};

#endif //COIN_HAS_CLP
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
