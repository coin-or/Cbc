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
#include "CbcCompareActual.hpp"
#include "CoinError.hpp"


/** Default Constructor

*/
CbcCompareDepth::CbcCompareDepth ()
  : CbcCompareBase()
{
  test_=this;
}

// Copy constructor 
CbcCompareDepth::CbcCompareDepth ( const CbcCompareDepth & rhs)
  :CbcCompareBase(rhs)

{
}

// Clone
CbcCompareBase *
CbcCompareDepth::clone() const
{
  return new CbcCompareDepth(*this);
}

// Assignment operator 
CbcCompareDepth & 
CbcCompareDepth::operator=( const CbcCompareDepth& rhs)
{
  if (this!=&rhs) {
    CbcCompareBase::operator=(rhs);
  }
  return *this;
}

// Destructor 
CbcCompareDepth::~CbcCompareDepth ()
{
}

// Returns true if y better than x
bool 
CbcCompareDepth::test (CbcNode * x, CbcNode * y)
{
  return x->depth() < y->depth();
}

/** Default Constructor

*/
CbcCompareObjective::CbcCompareObjective ()
  : CbcCompareBase()
{
  test_=this;
}

// Copy constructor 
CbcCompareObjective::CbcCompareObjective ( const CbcCompareObjective & rhs)
  :CbcCompareBase(rhs)

{
}

// Clone
CbcCompareBase *
CbcCompareObjective::clone() const
{
  return new CbcCompareObjective(*this);
}

// Assignment operator 
CbcCompareObjective & 
CbcCompareObjective::operator=( const CbcCompareObjective& rhs)
{
  if (this!=&rhs) {
    CbcCompareBase::operator=(rhs);
  }
  return *this;
}

// Destructor 
CbcCompareObjective::~CbcCompareObjective ()
{
}

// Returns true if y better than x
bool 
CbcCompareObjective::test (CbcNode * x, CbcNode * y)
{
  return x->objectiveValue() > y->objectiveValue();
}

/** Default Constructor

*/
CbcCompareDefault::CbcCompareDefault ()
  : CbcCompareBase(),
    weight_(-1.0),
    numberSolutions_(0),
    treeSize_(0)
{
  test_=this;
}

// Constructor with weight
CbcCompareDefault::CbcCompareDefault (double weight) 
  : CbcCompareBase(),
    weight_(weight) ,
    numberSolutions_(0),
    treeSize_(0)
{
  test_=this;
}


// Copy constructor 
CbcCompareDefault::CbcCompareDefault ( const CbcCompareDefault & rhs)
  :CbcCompareBase(rhs)

{
  weight_=rhs.weight_;
  numberSolutions_=rhs.numberSolutions_;
  treeSize_ = rhs.treeSize_;
}

// Clone
CbcCompareBase *
CbcCompareDefault::clone() const
{
  return new CbcCompareDefault(*this);
}

// Assignment operator 
CbcCompareDefault & 
CbcCompareDefault::operator=( const CbcCompareDefault& rhs)
{
  if (this!=&rhs) {
    CbcCompareBase::operator=(rhs);
    weight_=rhs.weight_;
    numberSolutions_=rhs.numberSolutions_;
    treeSize_ = rhs.treeSize_;
  }
  return *this;
}

// Destructor 
CbcCompareDefault::~CbcCompareDefault ()
{
}

// Returns true if y better than x
bool 
CbcCompareDefault::test (CbcNode * x, CbcNode * y)
{
  if (weight_<0.0||treeSize_>100000) {
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
CbcCompareDefault::newSolution(CbcModel * model,
			       double objectiveAtContinuous,
			       int numberInfeasibilitiesAtContinuous) 
{
  if (model->getSolutionCount()==model->getNumberHeuristicSolutions())
    return; // solution was got by rounding
  // set to get close to this solution
  double costPerInteger = 
    (model->getObjValue()-objectiveAtContinuous)/
    ((double) numberInfeasibilitiesAtContinuous);
  weight_ = 0.98*costPerInteger;
  numberSolutions_++;
  if (numberSolutions_>5)
    weight_ =0.0; // this searches on objective
}
// This allows method to change behavior 
bool 
CbcCompareDefault::every1000Nodes(CbcModel * model, int numberNodes)
{
  if (numberNodes>10000)
    weight_ =0.0; // this searches on objective
  // get size of tree
  treeSize_ = model->tree()->size();
  return numberNodes==11000; // resort if first time
}

/** Default Constructor

*/
CbcCompareEstimate::CbcCompareEstimate ()
  : CbcCompareBase()
{
  test_=this;
}

// Copy constructor 
CbcCompareEstimate::CbcCompareEstimate ( const CbcCompareEstimate & rhs)
  :CbcCompareBase(rhs)

{
}

// Clone
CbcCompareBase *
CbcCompareEstimate::clone() const
{
  return new CbcCompareEstimate(*this);
}

// Assignment operator 
CbcCompareEstimate & 
CbcCompareEstimate::operator=( const CbcCompareEstimate& rhs)
{
  if (this!=&rhs) {
    CbcCompareBase::operator=(rhs);
  }
  return *this;
}

// Destructor 
CbcCompareEstimate::~CbcCompareEstimate ()
{
}

// Returns true if y better than x
bool 
CbcCompareEstimate::test (CbcNode * x, CbcNode * y)
{
  return x->guessedObjectiveValue() >  y->guessedObjectiveValue() ;
}

