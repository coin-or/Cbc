// Copyright (C) 2002, International Business Machines
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
#include "CbcBranchLotsize.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"


/** Default Constructor

*/
CbcLotsize::CbcLotsize ()
  : CbcObject(),
    columnNumber_(-1),
    rangeType_(0),
    numberRanges_(0),
    largestGap_(0),
    bound_(NULL),
    range_(0)
{
}

/** Useful constructor

  Loads actual upper & lower bounds for the specified variable.
*/
CbcLotsize::CbcLotsize (CbcModel * model, 
				    int iColumn, int numberPoints,
			const double * points, bool range)
  : CbcObject(model)
{
  assert (numberPoints>0);
  columnNumber_ = iColumn ;
  // and set id so can be used for branching
  id_=iColumn;
  // sort ranges
  int * sort = new int[numberPoints];
  double * weight = new double [numberPoints];
  int i;
  if (range) {
    rangeType_=2;
  } else {
    rangeType_=1;
  }
  for (i=0;i<numberPoints;i++) {
    sort[i]=i;
    weight[i]=points[i*rangeType_];
  }
  CoinSort_2(weight,weight+numberPoints,sort);
  numberRanges_=1;
  largestGap_=0;
  if (rangeType_==1) {
    bound_ = new double[numberPoints+1];
    bound_[0]=weight[0];
    for (i=1;i<numberPoints;i++) {
      if (weight[i]!=weight[i-1]) 
	bound_[numberRanges_++]=weight[i];
    }
    // and for safety
    bound_[numberRanges_]=bound_[numberRanges_-1];
    for (i=1;i<numberRanges_;i++) {
      largestGap_ = CoinMax(largestGap_,bound_[i]-bound_[i-1]);
    }
  } else {
    bound_ = new double[2*numberPoints+2];
    bound_[0]=points[sort[0]*2];
    bound_[1]=points[sort[0]*2+1];
    double lo=bound_[0];
    double hi=bound_[1];
    assert (hi>=lo);
    for (i=1;i<numberPoints;i++) {
      double thisLo =points[sort[i]*2];
      double thisHi =points[sort[i]*2+1];
      assert (thisHi>=thisLo);
      if (thisLo>hi) {
	bound_[2*numberRanges_]=thisLo;
	bound_[2*numberRanges_+1]=thisHi;
	numberRanges_++;
	lo=thisLo;
	hi=thisHi;
      } else {
	//overlap
	hi=CoinMax(hi,thisHi);
	bound_[2*numberRanges_-1]=hi;
      }
    }
    // and for safety
    bound_[2*numberRanges_]=bound_[2*numberRanges_-2];
    bound_[2*numberRanges_+1]=bound_[2*numberRanges_-1];
    for (i=1;i<numberRanges_;i++) {
      largestGap_ = CoinMax(largestGap_,bound_[2*i]-bound_[2*i-1]);
    }
  }
  delete [] sort;
  delete [] weight;
  range_=0;
}

// Copy constructor 
CbcLotsize::CbcLotsize ( const CbcLotsize & rhs)
  :CbcObject(rhs)

{
  columnNumber_ = rhs.columnNumber_;
  rangeType_ = rhs.rangeType_;
  numberRanges_ = rhs.numberRanges_;
  range_ = rhs.range_;
  largestGap_ = rhs.largestGap_;
  if (numberRanges_) {
    assert (rangeType_>0&&rangeType_<3);
    bound_= new double [(numberRanges_+1)*rangeType_];
    memcpy(bound_,rhs.bound_,(numberRanges_+1)*rangeType_*sizeof(double));
  } else {
    bound_=NULL;
  }
}

// Clone
CbcObject *
CbcLotsize::clone() const
{
  return new CbcLotsize(*this);
}

// Assignment operator 
CbcLotsize & 
CbcLotsize::operator=( const CbcLotsize& rhs)
{
  if (this!=&rhs) {
    CbcObject::operator=(rhs);
    columnNumber_ = rhs.columnNumber_;
    rangeType_ = rhs.rangeType_;
    numberRanges_ = rhs.numberRanges_;
    largestGap_ = rhs.largestGap_;
    delete [] bound_;
    range_ = rhs.range_;
    if (numberRanges_) {
      assert (rangeType_>0&&rangeType_<3);
      bound_= new double [(numberRanges_+1)*rangeType_];
      memcpy(bound_,rhs.bound_,(numberRanges_+1)*rangeType_*sizeof(double));
    } else {
      bound_=NULL;
    }
  }
  return *this;
}

// Destructor 
CbcLotsize::~CbcLotsize ()
{
  delete [] bound_;
}
/* Finds range of interest so value is feasible in range range_ or infeasible 
   between hi[range_] and lo[range_+1].  Returns true if feasible.
*/
bool 
CbcLotsize::findRange(double value) const
{
  assert (range_>=0&&range_<numberRanges_+1);
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  int iLo;
  int iHi;
  double infeasibility=0.0;
  if (rangeType_==1) {
    if (value<bound_[range_]-integerTolerance) {
      iLo=0;
      iHi=range_-1;
    } else if (value<bound_[range_]+integerTolerance) {
      return true;
    } else if (value<bound_[range_+1]-integerTolerance) {
      return false;
    } else {
      iLo=range_+1;
      iHi=numberRanges_-1;
    }
    // check lo and hi
    bool found=false;
    if (value>bound_[iLo]-integerTolerance&&value<bound_[iLo+1]+integerTolerance) {
      range_=iLo;
      found=true;
    } else if (value>bound_[iHi]-integerTolerance&&value<bound_[iHi+1]+integerTolerance) {
      range_=iHi;
      found=true;
    } else {
      range_ = (iLo+iHi)>>1;
    }
    //points
    while (!found) {
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
    if (value-bound_[range_]<=bound_[range_+1]-value) {
      infeasibility = value-bound_[range_];
    } else {
      infeasibility = bound_[range_+1]-value;
      if (infeasibility<integerTolerance)
	range_++;
    }
    return (infeasibility<integerTolerance);
  } else {
    // ranges
    if (value<bound_[2*range_]-integerTolerance) {
      iLo=0;
      iHi=range_-1;
    } else if (value<bound_[2*range_+1]+integerTolerance) {
      return true;
    } else if (value<bound_[2*range_+2]-integerTolerance) {
      return false;
    } else {
      iLo=range_+1;
      iHi=numberRanges_-1;
    }
    // check lo and hi
    bool found=false;
    if (value>bound_[2*iLo]-integerTolerance&&value<bound_[2*iLo+2]-integerTolerance) {
      range_=iLo;
      found=true;
    } else if (value>=bound_[2*iHi]-integerTolerance) {
      range_=iHi;
      found=true;
    } else {
      range_ = (iLo+iHi)>>1;
    }
    //points
    while (!found) {
      if (value<bound_[2*range_]) {
	if (value>=bound_[2*range_-2]) {
	  // found
	  range_--;
	  break;
	} else {
	  iHi = range_;
	}
      } else {
	if (value<bound_[2*range_+2]) {
	  // found
	  break;
	} else {
	  iLo = range_;
	}
      }
      range_ = (iLo+iHi)>>1;
    }
    if (value>=bound_[2*range_]-integerTolerance&&value<=bound_[2*range_+1]+integerTolerance)
      infeasibility=0.0;
    else if (value-bound_[2*range_+1]<bound_[2*range_+2]-value) {
      infeasibility = value-bound_[2*range_+1];
    } else {
      infeasibility = bound_[2*range_+2]-value;
    }
    return (infeasibility<integerTolerance);
  }
}
/* Returns floor and ceiling
 */
void 
CbcLotsize::floorCeiling(double & floorLotsize, double & ceilingLotsize, double value,
			 double tolerance) const
{
  bool feasible=findRange(value);
  if (rangeType_==1) {
    floorLotsize=bound_[range_];
    ceilingLotsize=bound_[range_+1];
    // may be able to adjust
    if (feasible&&fabs(value-floorLotsize)>fabs(value-ceilingLotsize)) {
      floorLotsize=bound_[range_+1];
      ceilingLotsize=bound_[range_+2];
    }
  } else {
    // ranges
    assert (value>=bound_[2*range_+1]);
    floorLotsize=bound_[2*range_+1];
    ceilingLotsize=bound_[2*range_+2];
  }
}

// Infeasibility - large is 0.5
double 
CbcLotsize::infeasibility(int & preferredWay) const
{
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->currentSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  double value = solution[columnNumber_];
  value = CoinMax(value, lower[columnNumber_]);
  value = CoinMin(value, upper[columnNumber_]);
  /*printf("%d %g %g %g %g\n",columnNumber_,value,lower[columnNumber_],
    solution[columnNumber_],upper[columnNumber_]);*/
  assert (value>=bound_[0]&&value<=bound_[rangeType_*numberRanges_-1]);
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double infeasibility=0.0;
  bool feasible = findRange(value);
  if (!feasible) {
    if (rangeType_==1) {
      if (value-bound_[range_]<bound_[range_+1]-value) {
	preferredWay=-1;
	infeasibility = value-bound_[range_];
      } else {
	preferredWay=1;
	infeasibility = bound_[range_+1]-value;
      }
    } else {
      // ranges
      if (value-bound_[2*range_+1]<bound_[2*range_+2]-value) {
	preferredWay=-1;
	infeasibility = value-bound_[2*range_+1];
      } else {
	preferredWay=1;
	infeasibility = bound_[2*range_+2]-value;
      }
    }
  } else {
    // always satisfied
    preferredWay=-1;
  }
  if (infeasibility<integerTolerance)
    infeasibility=0.0;
  else
    infeasibility /= largestGap_;
  return infeasibility;
}
// This looks at solution and sets bounds to contain solution
/** More precisely: it first forces the variable within the existing
    bounds, and then tightens the bounds to make sure the variable is feasible
*/
void 
CbcLotsize::feasibleRegion()
{
  OsiSolverInterface * solver = model_->solver();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  const double * solution = model_->currentSolution();
  double value = solution[columnNumber_];
  value = CoinMax(value, lower[columnNumber_]);
  value = CoinMin(value, upper[columnNumber_]);
  findRange(value);
  double nearest;
  if (rangeType_==1) {
    nearest = bound_[range_];
    solver->setColLower(columnNumber_,nearest);
    solver->setColUpper(columnNumber_,nearest);
  } else {
    // ranges
    solver->setColLower(columnNumber_,bound_[2*range_]);
    solver->setColUpper(columnNumber_,bound_[2*range_+1]);
    if (value>bound_[2*range_+1]) 
      nearest=bound_[2*range_+1];
    else if (value<bound_[2*range_]) 
      nearest = bound_[2*range_];
    else
      nearest = value;
  }
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  // Scaling may have moved it a bit
  assert (fabs(value-nearest)<=100.0*integerTolerance);
}

// Creates a branching object
CbcBranchingObject * 
CbcLotsize::createBranch(int way) const
{
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->currentSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  double value = solution[columnNumber_];
  value = CoinMax(value, lower[columnNumber_]);
  value = CoinMin(value, upper[columnNumber_]);
  bool feasible = findRange(value);
  assert (!feasible);
  return new CbcLotsizeBranchingObject(model_,columnNumber_,way,
					     value,this);
}


/* Given valid solution (i.e. satisfied) and reduced costs etc
   returns a branching object which would give a new feasible
   point in direction reduced cost says would be cheaper.
   If no feasible point returns null
*/
CbcBranchingObject * 
CbcLotsize::preferredNewFeasible() const
{
  OsiSolverInterface * solver = model_->solver();
  double value = model_->currentSolution()[columnNumber_];

  bool feasible = findRange(value);
  assert (feasible);
  double dj = solver->getObjSense()*solver->getReducedCost()[columnNumber_];
  CbcLotsizeBranchingObject * object = NULL;
  double lo,up;
  if (dj>=0.0) {
    // can we go down
    if (range_) {
      // yes
      if (rangeType_==1) {
	lo = bound_[range_-1];
	up = bound_[range_-1];
      } else {
	lo = bound_[2*range_-2];
	up = bound_[2*range_-1];
      }
      object = new CbcLotsizeBranchingObject(model_,columnNumber_,-1,
					     lo,up);
    }
  } else {
    // can we go up
    if (range_<numberRanges_-1) {
      // yes
      if (rangeType_==1) {
	lo = bound_[range_+1];
	up = bound_[range_+1];
      } else {
	lo = bound_[2*range_+2];
	up = bound_[2*range_+3];
      }
      object = new CbcLotsizeBranchingObject(model_,columnNumber_,-1,
					     lo,up);
    }
  }
  return object;
}
  
/* Given valid solution (i.e. satisfied) and reduced costs etc
   returns a branching object which would give a new feasible
   point in direction opposite to one reduced cost says would be cheaper.
   If no feasible point returns null
*/
CbcBranchingObject * 
CbcLotsize::notPreferredNewFeasible() const 
{
  OsiSolverInterface * solver = model_->solver();
  double value = model_->currentSolution()[columnNumber_];

  double nearest = floor(value+0.5);
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  assert (fabs(value-nearest)<=integerTolerance);
  double dj = solver->getObjSense()*solver->getReducedCost()[columnNumber_];
  CbcLotsizeBranchingObject * object = NULL;
  double lo,up;
  if (dj<=0.0) {
    // can we go down
    if (range_) {
      // yes
      if (rangeType_==1) {
	lo = bound_[range_-1];
	up = bound_[range_-1];
      } else {
	lo = bound_[2*range_-2];
	up = bound_[2*range_-1];
      }
      object = new CbcLotsizeBranchingObject(model_,columnNumber_,-1,
					     lo,up);
    }
  } else {
    // can we go up
    if (range_<numberRanges_-1) {
      // yes
      if (rangeType_==1) {
	lo = bound_[range_+1];
	up = bound_[range_+1];
      } else {
	lo = bound_[2*range_+2];
	up = bound_[2*range_+3];
      }
      object = new CbcLotsizeBranchingObject(model_,columnNumber_,-1,
					     lo,up);
    }
  }
  return object;
}
  
/*
  Bounds may be tightened, so it may be good to be able to refresh the local
  copy of the original bounds.
 */
void 
CbcLotsize::resetBounds()
{
  //printf("resetBounds needs coding for CbcLotSize\n");
}


// Default Constructor 
CbcLotsizeBranchingObject::CbcLotsizeBranchingObject()
  :CbcBranchingObject()
{
  down_[0] = 0.0;
  down_[1] = 0.0;
  up_[0] = 0.0;
  up_[1] = 0.0;
}

// Useful constructor
CbcLotsizeBranchingObject::CbcLotsizeBranchingObject (CbcModel * model, 
						      int variable, int way , double value,
						      const CbcLotsize * lotsize)
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
// Useful constructor for fixing
CbcLotsizeBranchingObject::CbcLotsizeBranchingObject (CbcModel * model, 
						      int variable, int way,
						      double lowerValue, 
						      double upperValue)
  :CbcBranchingObject(model,variable,way,lowerValue)
{
  numberBranchesLeft_=1;
  down_[0] = lowerValue;
  down_[1] = upperValue;
  up_[0] = lowerValue;
  up_[1] = upperValue;
}
  

// Copy constructor 
CbcLotsizeBranchingObject::CbcLotsizeBranchingObject ( const CbcLotsizeBranchingObject & rhs) :CbcBranchingObject(rhs)
{
  down_[0] = rhs.down_[0];
  down_[1] = rhs.down_[1];
  up_[0] = rhs.up_[0];
  up_[1] = rhs.up_[1];
}

// Assignment operator 
CbcLotsizeBranchingObject & 
CbcLotsizeBranchingObject::operator=( const CbcLotsizeBranchingObject& rhs)
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
CbcBranchingObject * 
CbcLotsizeBranchingObject::clone() const
{ 
  return (new CbcLotsizeBranchingObject(*this));
}


// Destructor 
CbcLotsizeBranchingObject::~CbcLotsizeBranchingObject ()
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
CbcLotsizeBranchingObject::branch(bool normalBranch)
{
  if (model_->messageHandler()->logLevel()>2&&normalBranch)
    print(normalBranch);
  numberBranchesLeft_--;
  int iColumn = variable_;
  if (way_<0) {
#ifdef CBC_DEBUG
  { double olb,oub ;
    olb = model_->solver()->getColLower()[iColumn] ;
    oub = model_->solver()->getColUpper()[iColumn] ;
    printf("branching down on var %d: [%g,%g] => [%g,%g]\n",
	   iColumn,olb,oub,down_[0],down_[1]) ; }
#endif
    model_->solver()->setColLower(iColumn,down_[0]);
    model_->solver()->setColUpper(iColumn,down_[1]);
    way_=1;
  } else {
#ifdef CBC_DEBUG
  { double olb,oub ;
    olb = model_->solver()->getColLower()[iColumn] ;
    oub = model_->solver()->getColUpper()[iColumn] ;
    printf("branching up on var %d: [%g,%g] => [%g,%g]\n",
	   iColumn,olb,oub,up_[0],up_[1]) ; }
#endif
    model_->solver()->setColLower(iColumn,up_[0]);
    model_->solver()->setColUpper(iColumn,up_[1]);
    way_=-1;	  // Swap direction
  }
  return 0.0;
}
// Print
void
CbcLotsizeBranchingObject::print(bool normalBranch)
{
  int iColumn = variable_;
  if (way_<0) {
  { double olb,oub ;
    olb = model_->solver()->getColLower()[iColumn] ;
    oub = model_->solver()->getColUpper()[iColumn] ;
    printf("branching down on var %d: [%g,%g] => [%g,%g]\n",
	   iColumn,olb,oub,down_[0],down_[1]) ; }
  } else {
  { double olb,oub ;
    olb = model_->solver()->getColLower()[iColumn] ;
    oub = model_->solver()->getColUpper()[iColumn] ;
    printf("branching up on var %d: [%g,%g] => [%g,%g]\n",
	   iColumn,olb,oub,up_[0],up_[1]) ; }
  }
}
