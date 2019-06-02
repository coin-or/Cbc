/* $Id$ */
// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcHeuristicDiveGuided_H
#define CbcHeuristicDiveGuided_H

#include "CbcHeuristicDive.hpp"

/** DiveGuided class
 */

class CbcHeuristicDiveGuided : public CbcHeuristicDive {
public:
  // Default Constructor
  CbcHeuristicDiveGuided();

  // Constructor with model - assumed before cuts
  CbcHeuristicDiveGuided(CbcModel &model);

  // Copy constructor
  CbcHeuristicDiveGuided(const CbcHeuristicDiveGuided &);

  // Destructor
  ~CbcHeuristicDiveGuided();

  /// Clone
  virtual CbcHeuristicDiveGuided *clone() const;

  /// Assignment operator
  CbcHeuristicDiveGuided &operator=(const CbcHeuristicDiveGuided &rhs);

  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *fp);

  /// Tests if the heuristic can run
  virtual bool canHeuristicRun();

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
