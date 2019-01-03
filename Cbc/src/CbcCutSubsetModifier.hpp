// $Id$
// Copyright (C) 2003, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

//Edwin 11/25/09 carved out of CbcCutGenerator

#ifndef CbcCutSubsetModifier_H
#define CbcCutSubsetModifier_H

#include "OsiSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "CglCutGenerator.hpp"
#include "CbcCutModifier.hpp"

class CbcModel;
class OsiRowCut;
class OsiRowCutDebugger;
/** Simple cut modifier base class

    In exotic circumstances - cuts may need to be modified
    a) strengthened - changed
    b) weakened - changed
    c) deleted - set to NULL
    d) unchanged

    initially get rid of cuts with variables >= k
    could weaken
*/

class CbcCutSubsetModifier : public CbcCutModifier {
public:
  /// Default Constructor
  CbcCutSubsetModifier();

  /// Useful Constructor
  CbcCutSubsetModifier(int firstOdd);

  // Copy constructor
  CbcCutSubsetModifier(const CbcCutSubsetModifier &);

  /// Destructor
  virtual ~CbcCutSubsetModifier();

  /// Assignment
  CbcCutSubsetModifier &operator=(const CbcCutSubsetModifier &rhs);
  /// Clone
  virtual CbcCutModifier *clone() const;

  /** Returns
        0 unchanged
        1 strengthened
        2 weakened
        3 deleted
    */
  virtual int modify(const OsiSolverInterface *solver, OsiRowCut &cut);
  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *) {}

protected:
  /// data
  /// First odd variable
  int firstOdd_;
};

#endif //CbcCutSubsetModifier_H

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
