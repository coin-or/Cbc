// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

//Edwin 11/25/09 carved out of CbcCompareActual

#ifndef CbcCompareEstimate_H
#define CbcCompareEstimate_H

//#############################################################################
/*  These are alternative strategies for node traversal.
    They can take data etc for fine tuning

    At present the node list is stored as a heap and the "test"
    comparison function returns true if node y is better than node x.

*/
#include "CbcNode.hpp"
#include "CbcCompareBase.hpp"
#include "CbcCompare.hpp"
class CbcModel;

/* This is when rounding is being done
*/
class CbcCompareEstimate : public CbcCompareBase {
public:
  // Default Constructor
  CbcCompareEstimate();
  ~CbcCompareEstimate();
  // Copy constructor
  CbcCompareEstimate(const CbcCompareEstimate &rhs);

  // Assignment operator
  CbcCompareEstimate &operator=(const CbcCompareEstimate &rhs);

  /// Clone
  virtual CbcCompareBase *clone() const;
  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *fp);

  virtual bool test(CbcNode *x, CbcNode *y);
};

#endif //CbcCompareEstimate_H

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
