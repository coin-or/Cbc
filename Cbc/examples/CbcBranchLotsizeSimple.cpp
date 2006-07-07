// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cmath>
#include <cfloat>

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchLotsizeSimple.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"
#include "CoinHelperFunctions.hpp"
/** Default Constructor

*/
CbcLotsizeSimple::CbcLotsizeSimple ()
  : CbcObject(),
    columnNumber_(-1),
    numberRanges_(0),
    largestGap_(0),
    bound_(NULL),
    range_(0)
{
}

/** Useful constructor

  Loads actual upper & lower bounds for the specified variable.
*/
CbcLotsizeSimple::CbcLotsizeSimple (CbcModel * model, 
				    int iColumn, int numberPoints,
                                    const double * points)
  : CbcObject(model)
{
  columnNumber_ = iColumn ;
  // and set id so can be used for branching
  id_=iColumn;
  // sort ranges (just for safety)
  int * sort = new int[numberPoints];
  double * weight = new double [numberPoints];
  int i;
  for (i=0;i<numberPoints;i++) {
    sort[i]=i;
    weight[i]=points[i];
  }
  CoinSort_2(weight,weight+numberPoints,sort);
  numberRanges_=1;
  largestGap_=0;
  bound_ = new double[numberPoints+1];
  bound_[0]=weight[0];
  for (i=1;i<numberPoints;i++) {
    if (weight[i]!=weight[i-1]) 
      bound_[numberRanges_++]=weight[i]; // only different values for points
  }
  // and for safety add an extra
  bound_[numberRanges_]=bound_[numberRanges_-1];
  // Find largest gap (for averaging infeasibility)
  for (i=1;i<numberRanges_;i++) 
    largestGap_ = CoinMax(largestGap_,bound_[i]-bound_[i-1]);
  delete [] sort;
  delete [] weight;
  range_=0;
}

// Copy constructor 
CbcLotsizeSimple::CbcLotsizeSimple ( const CbcLotsizeSimple & rhs)
  :CbcObject(rhs)

{
  columnNumber_ = rhs.columnNumber_;
  numberRanges_ = rhs.numberRanges_;
  range_ = rhs.range_;
  largestGap_ = rhs.largestGap_;
  bound_ = CoinCopyOfArray(rhs.bound_,numberRanges_+1);
}

// Clone
CbcObject *
CbcLotsizeSimple::clone() const
{
  return new CbcLotsizeSimple(*this);
}

// Assignment operator 
CbcLotsizeSimple & 
CbcLotsizeSimple::operator=( const CbcLotsizeSimple& rhs)
{
  if (this!=&rhs) {
    CbcObject::operator=(rhs);
    columnNumber_ = rhs.columnNumber_;
    numberRanges_ = rhs.numberRanges_;
    largestGap_ = rhs.largestGap_;
    delete [] bound_;
    range_ = rhs.range_;
    bound_ = CoinCopyOfArray(rhs.bound_,numberRanges_+1);
  }
  return *this;
}

// Destructor 
CbcLotsizeSimple::~CbcLotsizeSimple ()
{
  delete [] bound_;
}
/* Finds range of interest so value is feasible in range range_ or infeasible 
   between bound_[range_] and bound_[range_+1].  Returns true if feasible.
*/
bool 
CbcLotsizeSimple::findRange(double value) const
{
  assert (range_>=0&&range_<numberRanges_+1);
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  int iLo;
  int iHi;
  double infeasibility=0.0;
  // see where we need to search
  if (value<bound_[range_]-integerTolerance) {
    // lower than current range
    iLo=0;
    iHi=range_-1;
  } else if (value<bound_[range_]+integerTolerance) {
    // feasible as close enough
    return true;
  } else if (value<bound_[range_+1]-integerTolerance) {
    // in this range - so infeasible
    return false;
  } else {
    // higher than current range
    iLo=range_+1;
    iHi=numberRanges_-1;
  }
  // check whether at one end of possible range
  if (value>bound_[iLo]-integerTolerance&&value<bound_[iLo+1]+integerTolerance) {
    range_=iLo;
  } else if (value>bound_[iHi]-integerTolerance&&value<bound_[iHi+1]+integerTolerance) {
    range_=iHi;
  } else {
    // do binary search to find range
    range_ = (iLo+iHi)>>1;
    while (true) {
      if (value<bound_[range_]) {
        if (value>=bound_[range_-1]) {
          // found
          range_--;
          break;
        } else {
          iHi = range_;
        }
      } else {
        if (value<bound_[range_+1]) {
          // found
          break;
        } else {
          iLo = range_;
        }
      }
      range_ = (iLo+iHi)>>1;
    }
  }
  // got range - see if feasible and find infeasibility
  if (value-bound_[range_]<=bound_[range_+1]-value) {
    infeasibility = value-bound_[range_];
  } else {
    infeasibility = bound_[range_+1]-value;
    if (infeasibility<integerTolerance)
      range_++; // feasible above - so move to that point
  }
  return (infeasibility<integerTolerance);
}
/* Returns floor and ceiling
   hardly necessary with this simple case
 */
void 
CbcLotsizeSimple::floorCeiling(double & floorLotsizeSimple, double & ceilingLotsizeSimple, double value,
			 double tolerance) const
{
  bool feasible=findRange(value);
  floorLotsizeSimple=bound_[range_];
  ceilingLotsizeSimple=bound_[range_+1];
  // may be able to adjust (only needed for fancy branching)
  if (feasible&&fabs(value-floorLotsizeSimple)>fabs(value-ceilingLotsizeSimple)) {
    floorLotsizeSimple=bound_[range_+1];
    ceilingLotsizeSimple=bound_[range_+2];
  }
}
// Simple function to return a solution value strictly feasible
double 
CbcLotsizeSimple::cleanValue() const
{
  // Get standard stuff
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  // Get value and clean it
  double value = solution[columnNumber_];
  value = CoinMax(value, lower[columnNumber_]);
  value = CoinMin(value, upper[columnNumber_]);
  return value;
}

// Infeasibility - large is 0.5
double 
CbcLotsizeSimple::infeasibility(int & preferredWay) const
{
  // Get value and clean it
  double value = cleanValue();
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  assert (value>=bound_[0]-integerTolerance
          &&value<=bound_[numberRanges_-1]+integerTolerance);
  double infeasibility=0.0;
  bool feasible = findRange(value);
  if (!feasible) {
    if (value-bound_[range_]<bound_[range_+1]-value) {
      preferredWay=-1;
      infeasibility = value-bound_[range_];
    } else {
      preferredWay=1;
      infeasibility = bound_[range_+1]-value;
    }
    // normalize
    infeasibility /= largestGap_;
    return infeasibility;
  } else {
    // satisfied
    preferredWay=-1;
    return 0.0;
  }
}
/* This looks at solution and sets bounds to contain solution
   More precisely: it first forces the variable within the existing
   bounds, and then tightens the bounds to make sure the variable is feasible
*/
void 
CbcLotsizeSimple::feasibleRegion()
{
  // Get value and clean it
  double value = cleanValue();
  findRange(value);
  double nearest;
  nearest = bound_[range_];
  // fix to nearest
  OsiSolverInterface * solver = model_->solver();
  solver->setColLower(columnNumber_,nearest);
  solver->setColUpper(columnNumber_,nearest);
  // Scaling may have moved it a bit
  // Lotsizing variables could be a lot larger
#ifndef NDEBUG
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  assert (fabs(value-nearest)<=(100.0+10.0*fabs(nearest))*integerTolerance);
#endif
}

// Creates a branching object
CbcBranchingObject * 
CbcLotsizeSimple::createBranch(int way) 
{
  // Get value and clean it
  double value = cleanValue();
  assert (!findRange(value)); // If it was feasible we would not be here
  return new CbcLotsizeSimpleBranchingObject(model_,columnNumber_,way,
					     value,this);
}

// Default Constructor 
CbcLotsizeSimpleBranchingObject::CbcLotsizeSimpleBranchingObject()
  :CbcBranchingObject()
{
  down_[0] = 0.0;
  down_[1] = 0.0;
  up_[0] = 0.0;
  up_[1] = 0.0;
}

// Useful constructor
CbcLotsizeSimpleBranchingObject::CbcLotsizeSimpleBranchingObject (CbcModel * model, 
						      int variable, int way , double value,
						      const CbcLotsizeSimple * lotsize)
  :CbcBranchingObject(model,variable,way,value)
{
  int iColumn = lotsize->modelSequence();
  assert (variable==iColumn);
  down_[0] = model_->solver()->getColLower()[iColumn];
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  lotsize->floorCeiling(down_[1],up_[0],value,integerTolerance);
  up_[1] = model->getColUpper()[iColumn];
}

// Copy constructor 
CbcLotsizeSimpleBranchingObject::CbcLotsizeSimpleBranchingObject ( const CbcLotsizeSimpleBranchingObject & rhs) :CbcBranchingObject(rhs)
{
  down_[0] = rhs.down_[0];
  down_[1] = rhs.down_[1];
  up_[0] = rhs.up_[0];
  up_[1] = rhs.up_[1];
}

// Assignment operator 
CbcLotsizeSimpleBranchingObject & 
CbcLotsizeSimpleBranchingObject::operator=( const CbcLotsizeSimpleBranchingObject& rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
    down_[0] = rhs.down_[0];
    down_[1] = rhs.down_[1];
    up_[0] = rhs.up_[0];
    up_[1] = rhs.up_[1];
  }
  return *this;
}
// Clone
CbcBranchingObject * 
CbcLotsizeSimpleBranchingObject::clone() const
{ 
  return (new CbcLotsizeSimpleBranchingObject(*this));
}


// Destructor 
CbcLotsizeSimpleBranchingObject::~CbcLotsizeSimpleBranchingObject ()
{
}

/*
  Perform a branch by adjusting the bounds of the specified variable. Note
  that each arm of the branch advances the object to the next arm by
  advancing the value of way_.

  Providing new values for the variable's lower and upper bounds for each
  branching direction gives a little bit of additional flexibility and will
  be easily extensible to multi-way branching.
*/
double
CbcLotsizeSimpleBranchingObject::branch(bool normalBranch)
{
  numberBranchesLeft_--;
  int iColumn = variable_;
  if (way_<0) {
    model_->solver()->setColLower(iColumn,down_[0]);
    model_->solver()->setColUpper(iColumn,down_[1]);
    way_=1;
  } else {
    model_->solver()->setColLower(iColumn,up_[0]);
    model_->solver()->setColUpper(iColumn,up_[1]);
    way_=-1;	  // Swap direction
  }
  return 0.0;
}
