/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CbcCompareBase_H
#define CbcCompareBase_H

//#############################################################################
/*  These are alternative strategies for node traversal.
    They can take data etc for fine tuning

    At present the node list is stored as a heap and the "test"
    comparison function returns true if node y is better than node x.

    This is rather inflexible so if the comparison functions wants
    it can signal to use alternative criterion on a complete pass
    throgh tree.

*/
#include "CbcNode.hpp"
#include "CbcConfig.h"

class CbcModel;
class CbcTree;
class CbcCompareBase {
public:
  // Default Constructor
  CbcCompareBase()
  {
    test_ = NULL;
    threaded_ = false;
  }

  /*! \brief Reconsider behaviour after discovering a new solution.
    
      This allows any method to change its behaviour. It is called
      after each solution.

      The method should return true if changes are made which will
      alter the evaluation criteria applied to a node. (So that in
      cases where the search tree is sorted, it can be properly
      rebuilt.)
    */
  virtual bool newSolution(CbcModel *) { return (false); }

  /*! \brief Reconsider behaviour after discovering a new solution.
    
      This allows any method to change its behaviour. It is called
      after each solution.

      The method should return true if changes are made which will
      alter the evaluation criteria applied to a node. (So that in
      cases where the search tree is sorted, it can be properly
      rebuilt.)
    */
  virtual bool newSolution(CbcModel *,
    double,
    int) { return (false); }

  // This allows any method to change behavior as it is called
  // after every 1000 nodes.
  // Return true if want tree re-sorted
  virtual bool every1000Nodes(CbcModel *, int)
  {
    return false;
  }

  /** Returns true if wants code to do scan with alternate criterion
        NOTE - this is temporarily disabled
    */
  virtual bool fullScan() const
  {
    return false;
  }

  virtual ~CbcCompareBase() {}
  /// Create C++ lines to get to current state
  virtual void generateCpp(FILE *) {}

  // Copy constructor
  CbcCompareBase(const CbcCompareBase &rhs)
  {
    test_ = rhs.test_;
    threaded_ = rhs.threaded_;
  }

  // Assignment operator
  CbcCompareBase &operator=(const CbcCompareBase &rhs)
  {
    if (this != &rhs) {
      test_ = rhs.test_;
      threaded_ = rhs.threaded_;
    }
    return *this;
  }

  /// Clone
  virtual CbcCompareBase *clone() const
  {
    abort();
    return NULL;
  }

  /// This is test function
  virtual bool test(CbcNode *, CbcNode *)
  {
    return true;
  }

  /// This is alternate test function
  virtual bool alternateTest(CbcNode *x, CbcNode *y)
  {
    return test(x, y);
  }

  bool operator()(CbcNode *x, CbcNode *y)
  {
    return test(x, y);
  }
  /// Further test if everything else equal
  inline bool equalityTest(CbcNode *x, CbcNode *y) const
  {
    assert(x);
    assert(y);
    if (!threaded_) {
      CbcNodeInfo *infoX = x->nodeInfo();
      assert(infoX);
      int nodeNumberX = infoX->nodeNumber();
      CbcNodeInfo *infoY = y->nodeInfo();
      assert(infoY);
      int nodeNumberY = infoY->nodeNumber();
      assert(nodeNumberX != nodeNumberY);
      return (nodeNumberX > nodeNumberY);
    } else {
      assert(x->nodeNumber() != y->nodeNumber());
      return (x->nodeNumber() > y->nodeNumber());
    }
  }
  /// Say threaded
  inline void sayThreaded()
  {
    threaded_ = true;
  }

protected:
  CbcCompareBase *test_;
  // If not threaded we can use better way to break ties
  bool threaded_;
};

#endif

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
