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

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcBranchCut.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"


/** Default Constructor

*/
CbcBranchCut::CbcBranchCut ()
  : CbcObject()
{
}

/* Constructor so model can be passed up
*/ 
CbcBranchCut::CbcBranchCut (CbcModel * model)
  : CbcObject(model)
{
}
// Copy constructor 
CbcBranchCut::CbcBranchCut ( const CbcBranchCut & rhs)
  :CbcObject(rhs)

{
}

// Clone
CbcObject *
CbcBranchCut::clone() const
{
  return new CbcBranchCut(*this);
}

// Assignment operator 
CbcBranchCut & 
CbcBranchCut::operator=( const CbcBranchCut& rhs)
{
  return *this;
}

// Destructor 
CbcBranchCut::~CbcBranchCut ()
{
}

// Infeasibility - large is 0.5
double 
CbcBranchCut::infeasibility(int & preferredWay) const
{
  throw CoinError("Use of base class","infeasibility","CbcBranchCut");
  preferredWay=-1;
  return 0.0;
}

// This looks at solution and sets bounds to contain solution
/** More precisely: it first forces the variable within the existing
    bounds, and then tightens the bounds to fix the variable at the
    nearest integer value.
*/
void 
CbcBranchCut::feasibleRegion()
{
}

// Creates a branching object
CbcBranchingObject * 
CbcBranchCut::createBranch(int way) const
{
  throw CoinError("Use of base class","createBranch","CbcBranchCut");
  return new CbcCutBranchingObject();
}


/* Given valid solution (i.e. satisfied) and reduced costs etc
   returns a branching object which would give a new feasible
   point in direction reduced cost says would be cheaper.
   If no feasible point returns null
*/
CbcBranchingObject * 
CbcBranchCut::preferredNewFeasible() const
{
  throw CoinError("Use of base class","preferredNewFeasible","CbcBranchCut");
  return new CbcCutBranchingObject();
}
  
/* Given valid solution (i.e. satisfied) and reduced costs etc
   returns a branching object which would give a new feasible
   point in direction opposite to one reduced cost says would be cheaper.
   If no feasible point returns null
*/
CbcBranchingObject * 
CbcBranchCut::notPreferredNewFeasible() const 
{
  throw CoinError("Use of base class","notPreferredNewFeasible","CbcBranchCut");
  return new CbcCutBranchingObject();
}
  
/*
  Bounds may be tightened, so it may be good to be able to refresh the local
  copy of the original bounds.
 */
void 
CbcBranchCut::resetBounds()
{
}


// Default Constructor 
CbcCutBranchingObject::CbcCutBranchingObject()
  :CbcBranchingObject()
{
  down_=OsiRowCut();
  up_=OsiRowCut();
}

// Useful constructor
CbcCutBranchingObject::CbcCutBranchingObject (CbcModel * model, 
					      OsiRowCut & down,
					      OsiRowCut &up)
  :CbcBranchingObject(model,0,-1,0.0)
{
  down_ = down;
  up_ = up;
}

// Copy constructor 
CbcCutBranchingObject::CbcCutBranchingObject ( const CbcCutBranchingObject & rhs) :CbcBranchingObject(rhs)
{
  down_ = rhs.down_;
  up_ = rhs.up_;
}

// Assignment operator 
CbcCutBranchingObject & 
CbcCutBranchingObject::operator=( const CbcCutBranchingObject& rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
    down_ = rhs.down_;
    up_ = rhs.up_;
  }
  return *this;
}
CbcBranchingObject * 
CbcCutBranchingObject::clone() const
{ 
  return (new CbcCutBranchingObject(*this));
}


// Destructor 
CbcCutBranchingObject::~CbcCutBranchingObject ()
{
}

/*
  Perform a branch by adjusting bounds and/or adding a cut. Note
  that each arm of the branch advances the object to the next arm by
  advancing the value of way_.

  Returns change in guessed objective on next branch
*/
double
CbcCutBranchingObject::branch(bool normalBranch)
{
  if (model_->messageHandler()->logLevel()>2&&normalBranch)
    print(normalBranch);
  numberBranchesLeft_--;
  OsiRowCut * cut;
  if (way_<0) {
    cut = &down_;
    way_=1;
  } else {
    cut = &up_;
    way_=-1;	  // Swap direction
  }
  // See if cut just fixes variables
  double lb = cut->lb();
  double ub = cut->ub();
  int n=cut->row().getNumElements();
  const int * column = cut->row().getIndices();
  const double * element = cut->row().getElements();
  OsiSolverInterface * solver = model_->solver();
  const double * upper = solver->getColUpper();
  const double * lower = solver->getColLower();
  double low = 0.0;
  double high=0.0;
  for (int i=0;i<n;i++) {
    int iColumn = column[i];
    double value = element[i];
    if (value>0.0) {
      high += upper[iColumn]*value;
      low += lower[iColumn]*value;
    } else {
      high += lower[iColumn]*value;
      low += upper[iColumn]*value;
    }
  }
  // assume cut was cunningly constructed so we need not worry too much about tolerances
  if (low+1.0e-8>=ub) {
    // fix
    for (int i=0;i<n;i++) {
      int iColumn = column[i];
      double value = element[i];
      if (value>0.0) {
	solver->setColUpper(iColumn,lower[iColumn]);
      } else {
	solver->setColLower(iColumn,upper[iColumn]);
      }
    }
  } else if (high-1.0e-8<=lb) {
    // fix
    for (int i=0;i<n;i++) {
      int iColumn = column[i];
      double value = element[i];
      if (value>0.0) {
	solver->setColLower(iColumn,upper[iColumn]);
      } else {
	solver->setColUpper(iColumn,lower[iColumn]);
      }
    }
  } else {
    // leave as cut
    model_->setNextRowCut(cut);
  }
  return 0.0;
}
// Print what would happen  
void
CbcCutBranchingObject::print(bool normalBranch)
{
  OsiRowCut * cut;
  if (way_<0) {
    cut = &down_;
    printf("CbcCut would branch down");
    way_=1;
  } else {
    cut = &up_;
    printf("CbcCut would branch up");
    way_=-1;	  // Swap direction
  }
  double lb = cut->lb();
  double ub = cut->ub();
  int n=cut->row().getNumElements();
  const int * column = cut->row().getIndices();
  const double * element = cut->row().getElements();
  if (n>5) {
    printf(" - %d elements, lo=%g, up=%g\n",n,lb,ub);
  } else {
    printf(" - %g <=",lb);
    for (int i=0;i<n;i++) {
      int iColumn = column[i];
      double value = element[i];
      printf(" (%d,%g)",iColumn,value);
    }
    printf(" <= %g\n",ub);
  }
}

// Return true if branch should fix variables
bool 
CbcCutBranchingObject::boundBranch() const
{
  return false;
}

/** Default Constructor

  Equivalent to an unspecified binary variable.
*/
CbcBranchOnReducedCost::CbcBranchOnReducedCost ()
  : CbcBranchCut(),
    djTolerance_(COIN_DBL_MAX),
    fractionFixed_(1.0),
    depth_(-1),
    alwaysCreate_(false)
{
}

/* Useful constructor - passed reduced cost tolerance and fraction we would like fixed.
   Also depth level to do at.
   Also whether to create branch if can't reach fraction.
*/ 
CbcBranchOnReducedCost::CbcBranchOnReducedCost (CbcModel * model, double djTolerance,
						double fractionFixed, int depth,
						bool alwaysCreate)
  : CbcBranchCut(model)
{
  djTolerance_ = djTolerance;
  fractionFixed_ = fractionFixed;
  depth_ = depth;
  alwaysCreate_ = alwaysCreate;
}
// Copy constructor 
CbcBranchOnReducedCost::CbcBranchOnReducedCost ( const CbcBranchOnReducedCost & rhs)
  :CbcBranchCut(rhs)
{
  djTolerance_ = rhs.djTolerance_;
  fractionFixed_ = rhs.fractionFixed_;
  depth_ = rhs.depth_;
  alwaysCreate_ = rhs.alwaysCreate_;
}

// Clone
CbcObject *
CbcBranchOnReducedCost::clone() const
{
  return new CbcBranchOnReducedCost(*this);
}

// Assignment operator 
CbcBranchOnReducedCost & 
CbcBranchOnReducedCost::operator=( const CbcBranchOnReducedCost& rhs)
{
  if (this!=&rhs) {
    CbcBranchCut::operator=(rhs);
    djTolerance_ = rhs.djTolerance_;
    fractionFixed_ = rhs.fractionFixed_;
    depth_ = rhs.depth_;
    alwaysCreate_ = rhs.alwaysCreate_;
  }
  return *this;
}

// Destructor 
CbcBranchOnReducedCost::~CbcBranchOnReducedCost ()
{
}
// Creates a branching object
CbcBranchingObject * 
CbcBranchOnReducedCost::createBranch(int way) const
{
  // by default way must be -1
  assert (way==-1);
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->currentSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  const double * dj = solver->getReducedCost();
  int i;
  int numberIntegers = model_->numberIntegers();
  const int * integerVariable = model_->integerVariable();
  int nSort=0;
  int nFixed=0;
  int * sort = new int[numberIntegers];
  double * dsort = new double[numberIntegers];
  double tolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  // make smaller ?
  tolerance = CoinMin(1.0e-8,tolerance);
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integerVariable[i];
    if (upper[iColumn]>lower[iColumn]) {
      if(solution[iColumn]<lower[iColumn]+tolerance) {
	if (dj[iColumn]>djTolerance_) {
	  dsort[nSort]=-dj[iColumn];
	  sort[nSort++]=iColumn;
	}
      } else if (solution[iColumn]>upper[iColumn]-tolerance) {
	if (dj[iColumn]<-djTolerance_) {
	  dsort[nSort]=dj[iColumn];
	  sort[nSort++]=iColumn;
	}
      }
    } else {
      nFixed++;
    }
  }
  // sort
  CoinSort_2(dsort,dsort+nSort,sort);
  // How many fixed are we aiming at
  int wantedFixed = (int) ((double)numberIntegers*fractionFixed_);
  nSort= CoinMin(nSort,wantedFixed-nFixed);
  OsiRowCut down;
  down.setLb(-COIN_DBL_MAX);
  double rhs=0.0;
  for (i=0;i<nSort;i++) {
    int iColumn = sort[i];
    if(solution[iColumn]<lower[iColumn]+tolerance) {
      rhs += lower[iColumn];
      dsort[i]=1.0;
    } else {
      rhs -= upper[iColumn];
      dsort[i]=-1.0;
    }
  }
  down.setUb(rhs);
  down.setRow(nSort,sort,dsort);
  delete [] sort;
  delete [] dsort;
  // up is same - just with rhs changed
  OsiRowCut up = down;
  up.setLb(rhs +1.0);
  up.setUb(COIN_DBL_MAX);
  CbcCutBranchingObject * newObject = 
    new CbcCutBranchingObject(model_,down,up);
  return newObject;
}
// Infeasibility - large is 0.5
double 
CbcBranchOnReducedCost::infeasibility(int & preferredWay) const
{
  preferredWay=-1;
  CbcNode * node = model_->currentNode();
  int depth;
  if (node)
    depth=CoinMax(node->depth(),0);
  else
    return 0.0;
  if (depth_<0) {
    return 0.0;
  } else if (depth_>0) {
    if ((depth%depth_)!=0)
      return 0.0;
  }
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->currentSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  const double * dj = solver->getReducedCost();
  int i;
  int numberIntegers = model_->numberIntegers();
  const int * integerVariable = model_->integerVariable();
  int nSort=0;
  int nFixed=0;
  double tolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  // make smaller ?
  tolerance = CoinMin(1.0e-8,tolerance);
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integerVariable[i];
    if (upper[iColumn]>lower[iColumn]) {
      if(solution[iColumn]<lower[iColumn]+tolerance) {
	if (dj[iColumn]>djTolerance_) {
	  nSort++;
	}
      } else if (solution[iColumn]>upper[iColumn]-tolerance) {
	if (dj[iColumn]<-djTolerance_) {
	  nSort++;
	}
      }
    } else {
      nFixed++;
    }
  }
  // How many fixed are we aiming at
  int wantedFixed = (int) ((double)numberIntegers*fractionFixed_);
  if (nFixed+nSort<wantedFixed&&!alwaysCreate_) {
    return 0.0;
  } else if (nFixed<wantedFixed) {
    return 1.0e20; // force ?
  } else {
    return 0.0;
  }
}
