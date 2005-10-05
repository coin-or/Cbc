// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "CbcMessage.hpp"
#include "CbcModel.hpp"
#include "CbcTree.hpp"
#include "CbcCompareUser.hpp"
#include "CoinError.hpp"
#include "CoinHelperFunctions.hpp"

/** Default Constructor

*/
CbcCompareUser::CbcCompareUser ()
  : CbcCompareBase(),
    weight_(-1.0),
    saveWeight_(0.0),
    numberSolutions_(0),
    count_(0),
    treeSize_(0)
{
  test_=this;
}

// Constructor with weight
CbcCompareUser::CbcCompareUser (double weight) 
  : CbcCompareBase(),
    weight_(weight) ,
    saveWeight_(0.0),
    numberSolutions_(0),
    count_(0),
    treeSize_(0)
{
  test_=this;
}


// Copy constructor 
CbcCompareUser::CbcCompareUser ( const CbcCompareUser & rhs)
  :CbcCompareBase(rhs)

{
  weight_=rhs.weight_;
  saveWeight_ = rhs.saveWeight_;
  numberSolutions_=rhs.numberSolutions_;
  count_ = rhs.count_;
  treeSize_ = rhs.treeSize_;
}

// Clone
CbcCompareBase *
CbcCompareUser::clone() const
{
  return new CbcCompareUser(*this);
}

// Assignment operator 
CbcCompareUser & 
CbcCompareUser::operator=( const CbcCompareUser& rhs)
{
  if (this!=&rhs) {
    CbcCompareBase::operator=(rhs);
    weight_=rhs.weight_;
    saveWeight_ = rhs.saveWeight_;
    numberSolutions_=rhs.numberSolutions_;
    count_ = rhs.count_;
    treeSize_ = rhs.treeSize_;
  }
  return *this;
}

// Destructor 
CbcCompareUser::~CbcCompareUser ()
{
}

// Returns true if y better than x
bool 
CbcCompareUser::test (CbcNode * x, CbcNode * y)
{
  if (x) {
    if (y) {
      if (weight_==-1.0) {
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
        double weight = CoinMax(weight_,0.0);
        return x->objectiveValue()+ weight*x->numberUnsatisfied() > 
          y->objectiveValue() + weight*y->numberUnsatisfied();
        //return x->guessedObjectiveValue()>y->guessedObjectiveValue();
      }
    } else {
      return false;
    }
  } else {
    return true;
  }
}
// This allows method to change behavior as it is called
// after each solution
void 
CbcCompareUser::newSolution(CbcModel * model,
			       double objectiveAtContinuous,
			       int numberInfeasibilitiesAtContinuous) 
{
  // set to get close to this solution
  double costPerInteger = 
    (model->getObjValue()-objectiveAtContinuous)/
    ((double) numberInfeasibilitiesAtContinuous);
  weight_ = 0.95*costPerInteger;
  saveWeight_=weight_;
  if (model->getSolutionCount()==model->getNumberHeuristicSolutions())
    return; // solution was got by rounding
  numberSolutions_++;
  if (numberSolutions_>5)
    weight_ =0.0; // this searches on objective
}
// This allows method to change behavior 
bool 
CbcCompareUser::every1000Nodes(CbcModel * model, int numberNodes)
{
  if (numberNodes>10000)
    weight_ =0.0; // this searches on objective
  else if (numberNodes==1000&&weight_==-2.0)
    weight_=-1.0; // Go to depth first
  // get size of tree
  treeSize_ = model->tree()->size();
  if (treeSize_>10000) {
    // set weight to reduce size most of time
    if (treeSize_>20000)
      weight_=-1.0;
    else if ((numberNodes%4000)!=0)
      weight_=-1.0;
    else
      weight_=saveWeight_;
  }
  return numberNodes==11000; // resort if first time
}
// Returns true if wants code to do scan with alternate criterion
bool 
CbcCompareUser::fullScan() const
{
  count_++;
  if (weight_)
    return (count_%10)==0;
  else
    return false;
}
// This is alternate test function
bool 
CbcCompareUser::alternateTest (CbcNode * x, CbcNode * y)
{
  if (x) {
    if (y) {
      return x->objectiveValue() > 
        y->objectiveValue() ;
    } else {
      return false;
    }
  } else {
    return true;
  }
}
