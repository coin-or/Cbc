// $Id$
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

//Edwin 11/24/09 carved out of CbcCompareActual

#ifndef CbcCompareDepth_H
#define CbcCompareDepth_H

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
// This is default before first solution
class CbcCompareDepth : public CbcCompareBase {
public:
  // Default Constructor
  CbcCompareDepth();

  ~CbcCompareDepth();
  // Copy constructor
  CbcCompareDepth(const CbcCompareDepth &rhs);

  // Assignment operator
  CbcCompareDepth &operator=(const CbcCompareDepth &rhs);

  /// Clone
  virtual CbcCompareBase *clone() const;
  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *fp);

  // This returns true if the depth of node y is greater than depth of node x
  virtual bool test(CbcNode *x, CbcNode *y);
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
