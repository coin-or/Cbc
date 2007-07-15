// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
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
  CbcCompareBase () {test_=NULL;};

  // This allows any method to change behavior as it is called
  // after each solution
  virtual void newSolution(CbcModel * model) {};

  // This Also allows any method to change behavior as it is called
  // after each solution
  virtual void newSolution(CbcModel * model,
			   double objectiveAtContinuous,
			   int numberInfeasibilitiesAtContinuous) {};

  // This allows any method to change behavior as it is called
  // after every 1000 nodes.
  // Return true if want tree re-sorted
  virtual bool every1000Nodes(CbcModel * model,int numberNodes) {return false;};

  /** Returns true if wants code to do scan with alternate criterion
      NOTE - this is temporarily disabled
  */
  virtual bool fullScan() const { return false;};

  virtual ~CbcCompareBase() {};
  /// Create C++ lines to get to current state
  virtual void generateCpp( FILE * fp) {};

  // Copy constructor 
  CbcCompareBase ( const CbcCompareBase & rhs)
  {test_=rhs.test_;};
   
  // Assignment operator 
  CbcCompareBase & operator=( const CbcCompareBase& rhs)
  {  if (this!=&rhs) {test_=rhs.test_;}
  return *this;
  };

  /// Clone
  virtual CbcCompareBase * clone() const=0;

  /// This is test function
  virtual bool test (CbcNode * x, CbcNode * y) {return true;};

  /// This is alternate test function
  virtual bool alternateTest (CbcNode * x, CbcNode * y) {return test(x,y);};

  bool operator() (CbcNode * x, CbcNode * y) {
    return test(x,y);
  }
  /// Further test if everything else equal
  inline bool equalityTest (CbcNode * x, CbcNode * y) const
  {
    assert (x);
    assert (y);
#ifndef CBC_THREAD
    CbcNodeInfo * infoX = x->nodeInfo();
    assert (infoX);
    int nodeNumberX = infoX->nodeNumber();
    CbcNodeInfo * infoY = y->nodeInfo();
    assert (infoY);
    int nodeNumberY = infoY->nodeNumber();
    assert (nodeNumberX!=nodeNumberY);
    return (nodeNumberX>nodeNumberY);
#else
    // doesn't work if threaded
    assert (x!=y);
    return (x>y);
#endif
  };
protected:
  CbcCompareBase * test_;
};
class CbcCompare {
public:
  CbcCompareBase * test_;
  // Default Constructor 
  CbcCompare () {test_=NULL;};

  virtual ~CbcCompare() {};

  bool operator() (CbcNode * x, CbcNode * y) {
    return test_->test(x,y);
  }
  /// This is alternate test function
  inline bool alternateTest (CbcNode * x, CbcNode * y) {return test_->alternateTest(x,y);};

  /// return comparison object
  inline CbcCompareBase * comparisonObject() const
  { return test_;};
};
//#############################################################################
/*  These can be alternative strategies for choosing variables
    Any descendant can be passed in by setVariableChoice
*/

class CbcChooseVariable {
public:
  // Default Constructor 
  CbcChooseVariable () {};

  virtual ~CbcChooseVariable() {};
  /** If strong branching, then only those passed in (and movement is that length)
      .  If not strong
      branching then all passed in and ignore movement.
      Returns which one chosen (or -1 if none).  way should be +1
      if branching up, -1 if down */
  virtual int chosen (const CbcModel * model,int numberToLookAt,
		      const int * which, const double * downMovement,
		      const double * upMovement, const double * solution,
		      int & way, double & value)=0;

};
#endif
