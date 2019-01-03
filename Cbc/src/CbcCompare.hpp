/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcCompare_H
#define CbcCompare_H

class CbcCompareBase;

class CbcCompare {
public:
  CbcCompareBase *test_;
  // Default Constructor
  CbcCompare()
  {
    test_ = NULL;
  }

  virtual ~CbcCompare() {}

  bool operator()(CbcNode *x, CbcNode *y)
  {
    return test_->test(x, y);
  }
  bool compareNodes(CbcNode *x, CbcNode *y)
  {
    return test_->test(x, y);
  }
  /// This is alternate test function
  inline bool alternateTest(CbcNode *x, CbcNode *y)
  {
    return test_->alternateTest(x, y);
  }

  /// return comparison object
  inline CbcCompareBase *comparisonObject() const
  {
    return test_;
  }
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
