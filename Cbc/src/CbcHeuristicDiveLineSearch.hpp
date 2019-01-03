/* $Id$ */
// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcHeuristicDiveLineSearch_H
#define CbcHeuristicDiveLineSearch_H

#include "CbcHeuristicDive.hpp"

/** DiveLineSearch class
 */

class CbcHeuristicDiveLineSearch : public CbcHeuristicDive {
public:
  // Default Constructor
  CbcHeuristicDiveLineSearch();

  // Constructor with model - assumed before cuts
  CbcHeuristicDiveLineSearch(CbcModel &model);

  // Copy constructor
  CbcHeuristicDiveLineSearch(const CbcHeuristicDiveLineSearch &);

  // Destructor
  ~CbcHeuristicDiveLineSearch();

  /// Clone
  virtual CbcHeuristicDiveLineSearch *clone() const;

  /// Assignment operator
  CbcHeuristicDiveLineSearch &operator=(const CbcHeuristicDiveLineSearch &rhs);

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
