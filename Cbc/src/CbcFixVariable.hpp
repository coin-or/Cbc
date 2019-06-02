// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/10/2009-- carved out of CbcBranchActual

#ifndef CbcFixVariable_H
#define CbcFixVariable_H

#include "CbcBranchBase.hpp"
/** Class for consequent bounds.
    When a variable is branched on it normally interacts with other variables by
    means of equations.  There are cases where we want to step outside LP and do something
    more directly e.g. fix bounds.  This class is for that.

    A state of -9999 means at LB, +9999 means at UB,
    others mean if fixed to that value.

 */

class CbcFixVariable : public CbcConsequence {

public:
  // Default Constructor
  CbcFixVariable();

  // One useful Constructor
  CbcFixVariable(int numberStates, const int *states, const int *numberNewLower, const int **newLowerValue,
    const int **lowerColumn,
    const int *numberNewUpper, const int **newUpperValue,
    const int **upperColumn);

  // Copy constructor
  CbcFixVariable(const CbcFixVariable &rhs);

  // Assignment operator
  CbcFixVariable &operator=(const CbcFixVariable &rhs);

  /// Clone
  virtual CbcConsequence *clone() const;

  /// Destructor
  virtual ~CbcFixVariable();

  /** Apply to an LP solver.  Action depends on state
     */
  virtual void applyToSolver(OsiSolverInterface *solver, int state) const;

protected:
  /// Number of states
  int numberStates_;
  /// Values of integers for various states
  int *states_;
  /// Start of information for each state (setting new lower)
  int *startLower_;
  /// Start of information for each state (setting new upper)
  int *startUpper_;
  /// For each variable new bounds
  double *newBound_;
  /// Variable
  int *variable_;
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
