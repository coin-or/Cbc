// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#ifndef CbcCompareUser_H
#define CbcCompareUser_H


//#############################################################################
/*  These are alternative strategies for node traversal.  
    They can take data etc for fine tuning 

    At present the node list is stored as a heap and the "test"
    comparison function returns true if node y is better than node x.

    This has been broken out into .hpp and .cpp because of some odd bug

*/
class CbcNode;
class CbcModel;
#include "CbcCompareBase.hpp"

/* Before first solution do depth first,
   then it is computed to hit first solution less 2%
*/
class CbcCompareUser  : public CbcCompareBase {
public:
  // Default Constructor 
  CbcCompareUser ();

  // Copy constructor 
  CbcCompareUser ( const CbcCompareUser &rhs);
   
  // Assignment operator 
  CbcCompareUser & operator=( const CbcCompareUser& rhs);

  /// Clone
  virtual CbcCompareBase * clone() const;

  virtual ~CbcCompareUser();

  /* 
     Return true if y better than x
     Node y is better than node x if y has fewer unsatisfied (greater depth on tie) or
     after solution weighted value of y is less than weighted value of x
  */
  virtual bool test (CbcNode * x, CbcNode * y);
  // This allows method to change behavior as it is called
  // after each solution
  virtual void newSolution(CbcModel * model,
			   double objectiveAtContinuous,
			   int numberInfeasibilitiesAtContinuous) ;
  // This allows method to change behavior 
  virtual void every1000Nodes(CbcModel * model, int numberNodes);
  inline void setWeight(double value)
  { weight_=value;};
  inline void setModel(CbcModel * model)
  { model_=model;};
private:
  // Weight for each infeasibility
  double weight_;
  // Number of solutions
  int numberSolutions_;
  // Pointer to model
  CbcModel * model_;
};

#endif
