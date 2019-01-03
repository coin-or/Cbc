// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/12/2009 carved from CbcBranchBase

#ifndef CbcConsequence_H
#define CbcConsequence_H

class OsiSolverInterface;

/** Abstract base class for consequent bounds.
    When a variable is branched on it normally interacts with other variables by
    means of equations.  There are cases where we want to step outside LP and do something
    more directly e.g. fix bounds.  This class is for that.

    At present it need not be virtual as only instance is CbcFixVariable, but ...

 */

class CbcConsequence {

public:
  // Default Constructor
  CbcConsequence();

  // Copy constructor
  CbcConsequence(const CbcConsequence &rhs);

  // Assignment operator
  CbcConsequence &operator=(const CbcConsequence &rhs);

  /// Clone
  virtual CbcConsequence *clone() const = 0;

  /// Destructor
  virtual ~CbcConsequence();

  /** Apply to an LP solver.  Action depends on state
     */
  virtual void applyToSolver(OsiSolverInterface *solver, int state) const = 0;

protected:
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
