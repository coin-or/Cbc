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
class CbcModel;
// This is default before first solution
class CbcCompareDepth : public CbcCompareBase{
public:
  // Default Constructor 
  CbcCompareDepth () ;

  ~CbcCompareDepth();
  // Copy constructor 
  CbcCompareDepth ( const CbcCompareDepth &rhs);
   
  // Assignment operator 
  CbcCompareDepth & operator=( const CbcCompareDepth& rhs);

  /// Clone
  virtual CbcCompareBase * clone() const;

  // This returns true if the depth of node y is greater than depth of node x
  virtual bool test (CbcNode * x, CbcNode * y);
};
class CbcCompareObjective  : public CbcCompareBase {
public:
  // Default Constructor 
  CbcCompareObjective ();

  virtual ~CbcCompareObjective();
  // Copy constructor 
  CbcCompareObjective ( const CbcCompareObjective &rhs);
   
  // Assignment operator 
  CbcCompareObjective & operator=( const CbcCompareObjective& rhs);

  /// Clone
  virtual CbcCompareBase * clone() const;

  /* This returns true if objective value of node y is less than
     objective value of node x */
  virtual bool test (CbcNode * x, CbcNode * y);
};
/* This is an example of a more complex rule with data
   It is default after first solution
   If weight is 0.0 then it is computed to hit first solution
   less 2%
*/
class CbcCompareDefault  : public CbcCompareBase {
public:
  // Default Constructor 
  CbcCompareDefault () ;
  // Constructor with weight
  CbcCompareDefault (double weight);

  // Copy constructor 
  CbcCompareDefault ( const CbcCompareDefault &rhs);
   
  // Assignment operator 
  CbcCompareDefault & operator=( const CbcCompareDefault& rhs);

  /// Clone
  virtual CbcCompareBase * clone() const;

  ~CbcCompareDefault() ;
  /* This returns true if weighted value of node y is less than
     weighted value of node x */
  virtual bool test (CbcNode * x, CbcNode * y) ;
  // This allows method to change behavior as it is called
  // after each solution
  virtual void newSolution(CbcModel * model,
			   double objectiveAtContinuous,
			   int numberInfeasibilitiesAtContinuous) ;
  // This allows method to change behavior 
  // Return true if want tree re-sorted
  virtual bool every1000Nodes(CbcModel * model,int numberNodes);

  inline double getWeight() const
  { return weight_;};
  inline void setWeight(double weight)
  { weight_ = weight;};
protected:
  // Weight for each infeasibility
  double weight_;
  // Number of solutions
  int numberSolutions_;
  // Tree size (at last check)
  int treeSize_;
};

/* This is when rounding is being done
*/
class CbcCompareEstimate  : public CbcCompareBase {
public:
  // Default Constructor 
  CbcCompareEstimate () ;
  ~CbcCompareEstimate() ;
  // Copy constructor 
  CbcCompareEstimate ( const CbcCompareEstimate &rhs);
   
  // Assignment operator 
  CbcCompareEstimate & operator=( const CbcCompareEstimate& rhs);

  /// Clone
  virtual CbcCompareBase * clone() const;

  virtual bool test (CbcNode * x, CbcNode * y) ;
};

#endif
