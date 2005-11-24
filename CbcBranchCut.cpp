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
/* Return true if branch created by object should fix variables
 */
bool 
CbcBranchCut::boundBranch() const 
{return false;}

// Creates a branching object
CbcBranchingObject * 
CbcBranchCut::createBranch(int way) 
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
  canFix_=false;
}

// Useful constructor
CbcCutBranchingObject::CbcCutBranchingObject (CbcModel * model, 
					      OsiRowCut & down,
					      OsiRowCut &up,
                                              bool canFix)
  :CbcBranchingObject(model,0,-1,0.0)
{
  down_ = down;
  up_ = up;
  canFix_ = canFix;
}

// Copy constructor 
CbcCutBranchingObject::CbcCutBranchingObject ( const CbcCutBranchingObject & rhs) :CbcBranchingObject(rhs)
{
  down_ = rhs.down_;
  up_ = rhs.up_;
  canFix_ = rhs.canFix_;
}

// Assignment operator 
CbcCutBranchingObject & 
CbcCutBranchingObject::operator=( const CbcCutBranchingObject& rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
    down_ = rhs.down_;
    up_ = rhs.up_;
    canFix_ = rhs.canFix_;
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
  if (low+1.0e-8>=ub&&canFix_) {
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
  } else if (high-1.0e-8<=lb&&canFix_) {
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
    model_->setNextRowCut(*cut);
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
  } else {
    cut = &up_;
    printf("CbcCut would branch up");
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
CbcBranchToFixLots::CbcBranchToFixLots ()
  : CbcBranchCut(),
    djTolerance_(COIN_DBL_MAX),
    fractionFixed_(1.0),
    mark_(NULL),
    depth_(-1),
    numberClean_(0),
    alwaysCreate_(false)
{
}

/* Useful constructor - passed reduced cost tolerance and fraction we would like fixed.
   Also depth level to do at.
   Also passed number of 1 rows which when clean triggers fix
   Always does if all 1 rows cleaned up and number>0 or if fraction columns reached
   Also whether to create branch if can't reach fraction.
*/ 
CbcBranchToFixLots::CbcBranchToFixLots (CbcModel * model, double djTolerance,
					double fractionFixed, int depth,
					int numberClean,
					const char * mark, bool alwaysCreate)
  : CbcBranchCut(model)
{
  djTolerance_ = djTolerance;
  fractionFixed_ = fractionFixed;
  if (mark) {
    int numberColumns = model->getNumCols();
    mark_ = new char[numberColumns];
    memcpy(mark_,mark,numberColumns);
  }
  depth_ = depth;
  assert (model);
  OsiSolverInterface * solver = model_->solver();
  matrixByRow_ = *solver->getMatrixByRow();
  numberClean_ = numberClean;
  alwaysCreate_ = alwaysCreate;
}
// Copy constructor 
CbcBranchToFixLots::CbcBranchToFixLots ( const CbcBranchToFixLots & rhs)
  :CbcBranchCut(rhs)
{
  djTolerance_ = rhs.djTolerance_;
  fractionFixed_ = rhs.fractionFixed_;
  int numberColumns = model_->getNumCols();
  mark_ = CoinCopyOfArray(rhs.mark_,numberColumns);
  matrixByRow_=rhs.matrixByRow_;
  depth_ = rhs.depth_;
  numberClean_ = rhs.numberClean_;
  alwaysCreate_ = rhs.alwaysCreate_;
}

// Clone
CbcObject *
CbcBranchToFixLots::clone() const
{
  return new CbcBranchToFixLots(*this);
}

// Assignment operator 
CbcBranchToFixLots & 
CbcBranchToFixLots::operator=( const CbcBranchToFixLots& rhs)
{
  if (this!=&rhs) {
    CbcBranchCut::operator=(rhs);
    djTolerance_ = rhs.djTolerance_;
    fractionFixed_ = rhs.fractionFixed_;
    int numberColumns = model_->getNumCols();
    delete [] mark_;
    mark_ = CoinCopyOfArray(rhs.mark_,numberColumns);
    matrixByRow_=rhs.matrixByRow_;
    depth_ = rhs.depth_;
    numberClean_ = rhs.numberClean_;
    alwaysCreate_ = rhs.alwaysCreate_;
  }
  return *this;
}

// Destructor 
CbcBranchToFixLots::~CbcBranchToFixLots ()
{
  delete [] mark_;
}
// Creates a branching object
CbcBranchingObject * 
CbcBranchToFixLots::createBranch(int way) 
{
  // by default way must be -1
  assert (way==-1);
  OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  const double * dj = solver->getReducedCost();
  int i;
  int numberIntegers = model_->numberIntegers();
  const int * integerVariable = model_->integerVariable();
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  // make smaller ?
  double tolerance = CoinMin(1.0e-8,integerTolerance);
  // How many fixed are we aiming at
  int wantedFixed = (int) ((double)numberIntegers*fractionFixed_);
  int nSort=0;
  int numberFixed=0;
  int numberColumns = solver->getNumCols();
  int * sort = new int[numberColumns];
  double * dsort = new double[numberColumns];
  int type = shallWe();
  assert (type);
  // Take clean first
  if (type==1) {
    for (i=0;i<numberIntegers;i++) {
      int iColumn = integerVariable[i];
      if (upper[iColumn]>lower[iColumn]) {
	if (!mark_||!mark_[iColumn]) {
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
	}
      } else {
	numberFixed++;
      }
    }
    // sort
    CoinSort_2(dsort,dsort+nSort,sort);
    nSort= CoinMin(nSort,wantedFixed-numberFixed);
  } else {
    int i;
    //const double * rowLower = solver->getRowLower();
    const double * rowUpper = solver->getRowUpper();
    // Row copy
    const double * elementByRow = matrixByRow_.getElements();
    const int * column = matrixByRow_.getIndices();
    const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
    const int * rowLength = matrixByRow_.getVectorLengths();
    const double * columnLower = solver->getColLower();
    const double * columnUpper = solver->getColUpper();
    const double * solution = solver->getColSolution();
    int numberColumns = solver->getNumCols();
    int numberRows = solver->getNumRows();
    for (i=0;i<numberColumns;i++) {
      sort[i]=i;
      if (columnLower[i]!=columnUpper[i]){
	dsort[i]=1.0e100;
      } else {
	dsort[i]=1.0e50;
	numberFixed++;
      }
    }
    for (i=0;i<numberRows;i++) {
      double rhsValue = rowUpper[i];
      bool oneRow=true;
      // check elements
      int numberUnsatisfied=0;
      for (int j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
	int iColumn = column[j];
	double value = elementByRow[j];
	double solValue = solution[iColumn];
	if (columnLower[iColumn]!=columnUpper[iColumn]) {
	  if (solValue<1.0-integerTolerance&&solValue>integerTolerance)
	    numberUnsatisfied++;
	  if (value!=1.0) {
	    oneRow=false;
	    break;
	  }
	} else {
	  rhsValue -= value*floor(solValue+0.5);
	}
      }
      if (oneRow&&rhsValue<=1.0+tolerance) {
	if (!numberUnsatisfied) {
	  for (int j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
	    int iColumn = column[j];
	    if (dsort[iColumn]>1.0e50){
	      dsort[iColumn]=0;
	      nSort++;
	    }
	  }
	}
      }
    }
    // sort
    CoinSort_2(dsort,dsort+numberColumns,sort);
  }
  OsiRowCut down;
  down.setLb(-COIN_DBL_MAX);
  double rhs=0.0;
  for (i=0;i<nSort;i++) {
    int iColumn = sort[i];
    if(solution[iColumn]<lower[iColumn]+tolerance) {
      rhs += lower[iColumn];
      dsort[i]=1.0;
      assert (!lower[iColumn]);
    } else {
      assert (solution[iColumn]>upper[iColumn]-tolerance);
      rhs -= upper[iColumn];
      dsort[i]=-1.0;
      //printf("%d at ub of %g\n",iColumn,upper[iColumn]);
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
  // Say can fix one way
  CbcCutBranchingObject * newObject = 
    new CbcCutBranchingObject(model_,down,up,true);
  if (model_->messageHandler()->logLevel()>1)
    printf("creating cut in CbcBranchCut\n");
  return newObject;
}
/* Does a lot of the work,
   Returns 0 if no good, 1 if dj, 2 if clean, 3 if both
*/
int 
CbcBranchToFixLots::shallWe() const
{
  int returnCode=0;
  OsiSolverInterface * solver = model_->solver();
  int numberRows = matrixByRow_.getNumRows();
  //if (numberRows!=solver->getNumRows())
  //return 0;
  const double * solution = model_->testSolution();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  const double * dj = solver->getReducedCost();
  int i;
  int numberIntegers = model_->numberIntegers();
  const int * integerVariable = model_->integerVariable();
  double integerTolerance = 
    model_->getDblParam(CbcModel::CbcIntegerTolerance);
  // make smaller ?
  double tolerance = CoinMin(1.0e-8,integerTolerance);
  // How many fixed are we aiming at
  int wantedFixed = (int) ((double)numberIntegers*fractionFixed_);
  if (djTolerance_<1.0e10) {
    int nSort=0;
    int numberFixed=0;
    for (i=0;i<numberIntegers;i++) {
      int iColumn = integerVariable[i];
      if (upper[iColumn]>lower[iColumn]) {
	if (!mark_||!mark_[iColumn]) {
	  if(solution[iColumn]<lower[iColumn]+tolerance) {
	    if (dj[iColumn]>djTolerance_) {
	      nSort++;
	    }
	  } else if (solution[iColumn]>upper[iColumn]-tolerance) {
	    if (dj[iColumn]<-djTolerance_) {
	      nSort++;
	    }
	  }
	}
      } else {
	numberFixed++;
      }
    }
    if (numberFixed+nSort<wantedFixed&&!alwaysCreate_) {
      returnCode = 0;
    } else if (numberFixed<wantedFixed) {
      returnCode = 1;
    } else {
      returnCode = 0;
    }
  }
  if (numberClean_) {
    // see how many rows clean
    int i;
    //const double * rowLower = solver->getRowLower();
    const double * rowUpper = solver->getRowUpper();
    // Row copy
    const double * elementByRow = matrixByRow_.getElements();
    const int * column = matrixByRow_.getIndices();
    const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
    const int * rowLength = matrixByRow_.getVectorLengths();
    const double * columnLower = solver->getColLower();
    const double * columnUpper = solver->getColUpper();
    const double * solution = solver->getColSolution();
    int numberClean=0;
    bool someToDoYet=false;
    int numberColumns = solver->getNumCols();
    char * mark = new char[numberColumns];
    int numberFixed=0;
    for (i=0;i<numberColumns;i++) {
      if (columnLower[i]!=columnUpper[i]){
	mark[i]=0;
      } else {
	mark[i]=1;
	numberFixed++;
      }
    }
    int numberNewFixed=0;
    for (i=0;i<numberRows;i++) {
      double rhsValue = rowUpper[i];
      bool oneRow=true;
      // check elements
      int numberUnsatisfied=0;
      for (int j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
	int iColumn = column[j];
	double value = elementByRow[j];
	double solValue = solution[iColumn];
	if (columnLower[iColumn]!=columnUpper[iColumn]) {
	  if (solValue<1.0-integerTolerance&&solValue>integerTolerance)
	    numberUnsatisfied++;
	  if (value!=1.0) {
	    oneRow=false;
	    break;
	  }
	} else {
	  rhsValue -= value*floor(solValue+0.5);
	}
      }
      if (oneRow&&rhsValue<=1.0+tolerance) {
	if (numberUnsatisfied) {
	  someToDoYet=true;
	} else {
	  numberClean++;
	  for (int j=rowStart[i];j<rowStart[i]+rowLength[i];j++) {
	    int iColumn = column[j];
	    if (columnLower[iColumn]!=columnUpper[iColumn]&&!mark[iColumn]){
	      mark[iColumn]=1;
	      numberNewFixed++;
	    }
	  }
	}
      }
    }
    delete [] mark;
    //printf("%d clean, %d old fixed, %d new fixed\n",
    //   numberClean,numberFixed,numberNewFixed);
    if (someToDoYet&&numberClean<numberClean_
	&&numberNewFixed+numberFixed<wantedFixed) {
    } else if (numberFixed<wantedFixed) {
      returnCode |= 2;
    } else {
    }
  }
  return returnCode;
}
// Infeasibility - large is 0.5
double 
CbcBranchToFixLots::infeasibility(int & preferredWay) const
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
  if (!shallWe())
    return 0.0;
  else
    return 1.0e20;
}

/** Default Constructor
*/
CbcBranchAllDifferent::CbcBranchAllDifferent ()
  : CbcBranchCut(),
    numberInSet_(0),
    which_(NULL)
{
}

/* Useful constructor - passed set of variables
*/ 
CbcBranchAllDifferent::CbcBranchAllDifferent (CbcModel * model, int numberInSet,
                                              const int * members)
  : CbcBranchCut(model)
{
  numberInSet_=numberInSet;
  which_ = CoinCopyOfArray(members,numberInSet_);
}
// Copy constructor 
CbcBranchAllDifferent::CbcBranchAllDifferent ( const CbcBranchAllDifferent & rhs)
  :CbcBranchCut(rhs)
{
  numberInSet_=rhs.numberInSet_;
  which_ = CoinCopyOfArray(rhs.which_,numberInSet_);
}

// Clone
CbcObject *
CbcBranchAllDifferent::clone() const
{
  return new CbcBranchAllDifferent(*this);
}

// Assignment operator 
CbcBranchAllDifferent & 
CbcBranchAllDifferent::operator=( const CbcBranchAllDifferent& rhs)
{
  if (this!=&rhs) {
    CbcBranchCut::operator=(rhs);
    delete [] which_;
    numberInSet_=rhs.numberInSet_;
    which_ = CoinCopyOfArray(rhs.which_,numberInSet_);
  }
  return *this;
}

// Destructor 
CbcBranchAllDifferent::~CbcBranchAllDifferent ()
{
  delete [] which_;
}
// Creates a branching object
CbcBranchingObject * 
CbcBranchAllDifferent::createBranch(int way) 
{
  // by default way must be -1
  assert (way==-1);
  const double * solution = model_->testSolution();
  double * values = new double[numberInSet_];
  int * which = new int[numberInSet_];
  int i;
  for (i=0;i<numberInSet_;i++) {
    int iColumn = which_[i];
    values[i]=solution[iColumn];
    which[i]=iColumn;
  }
  CoinSort_2(values,values+numberInSet_,which);
  double last = -1.0;
  double closest=1.0;
  int worst=-1;
  for (i=0;i<numberInSet_;i++) {
    if (values[i]-last<closest) {
      closest=values[i]-last;
      worst=i-1;
    }
    last=values[i];
  }
  assert (closest<=0.99999);
  OsiRowCut down;
  down.setLb(-COIN_DBL_MAX);
  down.setUb(-1.0);
  int pair[2];
  double elements[]={1.0,-1.0};
  pair[0]=which[worst];
  pair[1]=which[worst+1];
  delete [] values;
  delete [] which;
  down.setRow(2,pair,elements);
  // up is same - just with rhs changed
  OsiRowCut up = down;
  up.setLb(1.0);
  up.setUb(COIN_DBL_MAX);
  // Say is not a fix type branch
  CbcCutBranchingObject * newObject = 
    new CbcCutBranchingObject(model_,down,up,false);
  if (model_->messageHandler()->logLevel()>1)
    printf("creating cut in CbcBranchCut\n");
  return newObject;
}
// Infeasibility - large is 0.5
double 
CbcBranchAllDifferent::infeasibility(int & preferredWay) const
{
  preferredWay=-1;
  //OsiSolverInterface * solver = model_->solver();
  const double * solution = model_->testSolution();
  //const double * lower = solver->getColLower();
  //const double * upper = solver->getColUpper();
  double * values = new double[numberInSet_];
  int i;
  for (i=0;i<numberInSet_;i++) {
    int iColumn = which_[i];
    values[i]=solution[iColumn];
  }
  std::sort(values,values+numberInSet_);
  double last = -1.0;
  double closest=1.0;
  for (i=0;i<numberInSet_;i++) {
    if (values[i]-last<closest) {
      closest=values[i]-last;
    }
    last=values[i];
  }
  delete [] values;
  if (closest>0.99999)
    return 0.0;
  else
    return 0.5*(1.0-closest);
}
