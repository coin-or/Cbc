// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcCompareActual_H
#define CbcCompareActual_H


//#############################################################################
/*  These are alternative strategies for node traversal.  
    They can take data etc for fine tuning 

    At present the node list is stored as a heap and the "test"
    comparison function returns true if node y is better than node x.

*/
#include "CbcNode.hpp"
#include "CbcCompareBase.hpp"

// This is default before first solution
class CbcCompareDepth : public CbcCompareBase{
public:
  // Default Constructor 
  CbcCompareDepth () {test_=this;};

  ~CbcCompareDepth() {};
  // Copy constructor 
  CbcCompareDepth ( const CbcCompareDepth &rhs)
    : CbcCompareBase(rhs)
  {};
   
  // Assignment operator 
  CbcCompareDepth & operator=( const CbcCompareDepth& rhs)
  {  if (this!=&rhs) { CbcCompareBase::operator=(rhs);}
  return *this;
  };

  /// Clone
  virtual CbcCompareBase * clone() const
  { return new CbcCompareDepth (*this);
  }

  // This returns true if the depth of node y is greater than depth of node x
  virtual bool test (CbcNode * x, CbcNode * y) {
    return x->depth() < y->depth();
  }
};
class CbcCompareObjective  : public CbcCompareBase {
public:
  // Default Constructor 
  CbcCompareObjective () {test_=this;};

  virtual ~CbcCompareObjective() {};
  // Copy constructor 
  CbcCompareObjective ( const CbcCompareObjective &rhs)
    : CbcCompareBase(rhs)
  {};
   
  // Assignment operator 
  CbcCompareObjective & operator=( const CbcCompareObjective& rhs)
  {  if (this!=&rhs) { CbcCompareBase::operator=(rhs);
  }
  return *this;
  };

  /// Clone
  virtual CbcCompareBase * clone() const
  { return new CbcCompareObjective (*this);
  }


  /* This returns true if objective value of node y is less than
     objective value of node x */
  virtual bool test (CbcNode * x, CbcNode * y) {
    return x->objectiveValue() > y->objectiveValue();
  }
};
/* This is an example of a more complex rule with data
   It is default after first solution
   If weight is 0.0 then it is computed to hit first solution
   less 2%
*/
class CbcCompareDefault  : public CbcCompareBase {
public:
  // Weight for each infeasibility
  double weight_;
  // Default Constructor 
  CbcCompareDefault () : weight_(0.0) {test_=this;};
  // Constructor with weight
  CbcCompareDefault (double weight) : weight_(weight) {test_=this;};

  // Copy constructor 
  CbcCompareDefault ( const CbcCompareDefault &rhs)
    : CbcCompareBase(rhs)
  {};
   
  // Assignment operator 
  CbcCompareDefault & operator=( const CbcCompareDefault& rhs)
  {  if (this!=&rhs) { CbcCompareBase::operator=(rhs);
}
  return *this;
  };

  /// Clone
  virtual CbcCompareBase * clone() const
  { return new CbcCompareDefault (*this);
  }

  ~CbcCompareDefault() {};
  /* This returns true if weighted value of node y is less than
     weighted value of node x */
  virtual bool test (CbcNode * x, CbcNode * y) {
    return x->objectiveValue()+ weight_*x->numberUnsatisfied() > 
      y->objectiveValue() + weight_*y->numberUnsatisfied();
  }
  double getWeight() const
  { return weight_;};
  void setWeight(double weight)
  { weight_ = weight;};
};

/* This is when rounding is being done
*/
class CbcCompareEstimate  : public CbcCompareBase {
public:
  // Default Constructor 
  CbcCompareEstimate () {test_=this;};
  ~CbcCompareEstimate() {};
  // Copy constructor 
  CbcCompareEstimate ( const CbcCompareEstimate &rhs)
    : CbcCompareBase(rhs)
  {};
   
  // Assignment operator 
  CbcCompareEstimate & operator=( const CbcCompareEstimate& rhs)
  {  if (this!=&rhs) { CbcCompareBase::operator=(rhs);
}
  return *this;
  };

  /// Clone
  virtual CbcCompareBase * clone() const
  { return new CbcCompareEstimate (*this);
  }

  virtual bool test (CbcNode * x, CbcNode * y) {
    return x->guessedObjectiveValue() >  y->guessedObjectiveValue() ;
  }
};

#endif
