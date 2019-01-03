/* $Id$ */
// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcHeuristicDiveFractional_H
#define CbcHeuristicDiveFractional_H

#include "CbcHeuristicDive.hpp"

/** DiveFractional class
 */

class CbcHeuristicDiveFractional : public CbcHeuristicDive {
public:
  // Default Constructor
  CbcHeuristicDiveFractional();

  // Constructor with model - assumed before cuts
  CbcHeuristicDiveFractional(CbcModel &model);

  // Copy constructor
  CbcHeuristicDiveFractional(const CbcHeuristicDiveFractional &);

  // Destructor
  ~CbcHeuristicDiveFractional();

  /// Clone
  virtual CbcHeuristicDiveFractional *clone() const;

  /// Assignment operator
  CbcHeuristicDiveFractional &operator=(const CbcHeuristicDiveFractional &rhs);

  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *fp);

  /// Selects the next variable to branch on
  /** Returns true if all the fractional variables can be trivially
        rounded. Returns false, if there is at least one fractional variable
        that is not trivially roundable. In this case, the bestColumn
        returned will not be trivially roundable.
    */
  virtual bool selectVariableToBranch(OsiSolverInterface *solver,
    const double *newSolution,
    int &bestColumn,
    int &bestRound);
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
