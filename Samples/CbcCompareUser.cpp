// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.


//#############################################################################
/*  These are alternative strategies for node traversal.  
    They can take data etc for fine tuning 

    At present the node list is stored as a heap and the "test"
    comparison function returns true if node y is better than node x.

*/
#include "CbcModel.hpp"
#include "CbcNode.hpp"
#include "CbcCompareUser.hpp"

// Default Constructor 
CbcCompareUser::CbcCompareUser ()
  : weight_(-1.0), numberSolutions_(0), model_(NULL)
{
  test_=this;
}

CbcCompareUser::~CbcCompareUser()
{
}
// Copy constructor 
CbcCompareUser::CbcCompareUser ( const CbcCompareUser &rhs)
: CbcCompareBase(rhs)
{
  weight_=rhs.weight_;
  numberSolutions_=rhs.numberSolutions_;
  model_=rhs.model_;
}
   
// Assignment operator 
CbcCompareUser & 
CbcCompareUser::operator=( const CbcCompareUser& rhs)
{  
  if (this!=&rhs) { 
    CbcCompareBase::operator=(rhs);
    weight_=rhs.weight_;
    numberSolutions_=rhs.numberSolutions_;
    model_=rhs.model_;
  }
  return *this;
}

// Clone
CbcCompareBase * 
CbcCompareUser::clone() const
{ 
  return new CbcCompareUser (*this);
}

/* 
   Return true if y better than x
   Node y is better than node x if y has fewer unsatisfied (greater depth on tie) or
   after solution weighted value of y is less than weighted value of x
*/
bool 
CbcCompareUser::test (CbcNode * x, CbcNode * y) 
{
  if (weight_<0.0) {
    // before solution
    /* printf("x %d %d %g, y %d %d %g\n",
       x->numberUnsatisfied(),x->depth(),x->objectiveValue(),
       y->numberUnsatisfied(),y->depth(),y->objectiveValue()); */
    if (x->numberUnsatisfied() > y->numberUnsatisfied())
      return true;
    else if (x->numberUnsatisfied() < y->numberUnsatisfied())
      return false;
    else
      return x->depth() < y->depth();
  } else {
    // after solution
    return x->objectiveValue()+ weight_*x->numberUnsatisfied() > 
      y->objectiveValue() + weight_*y->numberUnsatisfied();
  }
}
// This allows method to change behavior as it is called
// after each solution
void 
CbcCompareUser::newSolution(CbcModel * model,
			    double objectiveAtContinuous,
			    int numberInfeasibilitiesAtContinuous) 
{
  if (!model_)
    model_=model;
  if (model->getSolutionCount()==model->getNumberHeuristicSolutions())
    return; // solution was got by rounding
  // set to get close to this solution
  double costPerInteger = 
    (model->getObjValue()-objectiveAtContinuous)/
    ((double) numberInfeasibilitiesAtContinuous);
  weight_ = 0.98*costPerInteger;
  numberSolutions_++;
  if (numberSolutions_>50)
    weight_ =0.0; // this searches on objective
  if (model->messageHandler()->logLevel()>1)
    printf("new weight %g\n",weight_);
}
// This allows method to change behavior 
bool 
CbcCompareUser::every1000Nodes(CbcModel * model, int numberNodes)
{
  if (numberNodes>10000)
    weight_ =0.0; // this searches on objective
  return numberNodes==11000; // resort if first time
}
