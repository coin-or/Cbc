/* $Id$ */
// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcFeasibilityBase_H
#define CbcFeasibilityBase_H

//#############################################################################
/*  There are cases where the user wants to control how CBC sees the problems feasibility.
    The user may want to examine the problem and say :
    a) The default looks OK
    b) Pretend this problem is Integer feasible
    c) Pretend this problem is infeasible even though it looks feasible

    This simple class allows user to do that.

*/

class CbcModel;
class CbcFeasibilityBase {
public:
  // Default Constructor
  CbcFeasibilityBase() {}

  /**
       On input mode:
       0 - called after a solve but before any cuts
       -1 - called after strong branching
       Returns :
       0 - no opinion
       -1 pretend infeasible
       1 pretend integer solution
    */
  virtual int feasible(CbcModel *, int)
  {
    return 0;
  }

  virtual ~CbcFeasibilityBase() {}

  // Copy constructor
  CbcFeasibilityBase(const CbcFeasibilityBase &) {}

  // Assignment operator
  CbcFeasibilityBase &operator=(const CbcFeasibilityBase &)
  {
    return *this;
  }

  /// Clone
  virtual CbcFeasibilityBase *clone() const
  {
    return new CbcFeasibilityBase(*this);
  }
};
#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
