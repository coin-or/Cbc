// Copyright (C) 2006, International Business Machines
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
#include "OsiBranchLink.hpp"
#include "CoinError.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinWarmStartBasis.hpp"

// Default Constructor 
OsiOldLink::OsiOldLink ()
  : OsiSOS(),
    numberLinks_(0)
{
}

// Useful constructor (which are indices)
OsiOldLink::OsiOldLink (const OsiSolverInterface * solver,  int numberMembers,
	   int numberLinks, int first , const double * weights, int identifier)
  : OsiSOS(),
    numberLinks_(numberLinks)
{
  numberMembers_ = numberMembers;
  members_ = NULL;
  sosType_ = 1;
  if (numberMembers_) {
    weights_ = new double[numberMembers_];
    members_ = new int[numberMembers_*numberLinks_];
    if (weights) {
      memcpy(weights_,weights,numberMembers_*sizeof(double));
    } else {
      for (int i=0;i<numberMembers_;i++)
        weights_[i]=i;
    }
    // weights must be increasing
    int i;
    double last=-COIN_DBL_MAX;
    for (i=0;i<numberMembers_;i++) {
      assert (weights_[i]>last+1.0e-12);
      last=weights_[i];
    }
    for (i=0;i<numberMembers_*numberLinks_;i++) {
      members_[i]=first+i;
    }
  } else {
    weights_ = NULL;
  }
}

// Useful constructor (which are indices)
OsiOldLink::OsiOldLink (const OsiSolverInterface * solver,  int numberMembers,
	   int numberLinks, int sosType, const int * which , const double * weights, int identifier)
  : OsiSOS(),
    numberLinks_(numberLinks)
{
  numberMembers_ = numberMembers;
  members_ = NULL;
  sosType_ = 1;
  if (numberMembers_) {
    weights_ = new double[numberMembers_];
    members_ = new int[numberMembers_*numberLinks_];
    if (weights) {
      memcpy(weights_,weights,numberMembers_*sizeof(double));
    } else {
      for (int i=0;i<numberMembers_;i++)
        weights_[i]=i;
    }
    // weights must be increasing
    int i;
    double last=-COIN_DBL_MAX;
    for (i=0;i<numberMembers_;i++) {
      assert (weights_[i]>last+1.0e-12);
      last=weights_[i];
    }
    for (i=0;i<numberMembers_*numberLinks_;i++) {
      members_[i]= which[i];
    }
  } else {
    weights_ = NULL;
  }
}

// Copy constructor 
OsiOldLink::OsiOldLink ( const OsiOldLink & rhs)
  :OsiSOS(rhs)
{
  numberLinks_ = rhs.numberLinks_;
  if (numberMembers_) {
    delete [] members_;
    members_ = CoinCopyOfArray(rhs.members_,numberMembers_*numberLinks_);
  }
}

// Clone
OsiObject *
OsiOldLink::clone() const
{
  return new OsiOldLink(*this);
}

// Assignment operator 
OsiOldLink & 
OsiOldLink::operator=( const OsiOldLink& rhs)
{
  if (this!=&rhs) {
    OsiSOS::operator=(rhs);
    delete [] members_;
    numberLinks_ = rhs.numberLinks_;
    if (numberMembers_) {
      members_ = CoinCopyOfArray(rhs.members_,numberMembers_*numberLinks_);
    } else {
      members_ = NULL;
    }
  }
  return *this;
}

// Destructor 
OsiOldLink::~OsiOldLink ()
{
}

// Infeasibility - large is 0.5
double 
OsiOldLink::infeasibility(const OsiBranchingInformation * info,int & whichWay) const
{
  int j;
  int firstNonZero=-1;
  int lastNonZero = -1;
  const double * solution = info->solution_;
  //const double * lower = info->lower_;
  const double * upper = info->upper_;
  double integerTolerance = info->integerTolerance_;
  double weight = 0.0;
  double sum =0.0;

  // check bounds etc
  double lastWeight=-1.0e100;
  int base=0;
  for (j=0;j<numberMembers_;j++) {
    for (int k=0;k<numberLinks_;k++) {
      int iColumn = members_[base+k];
      if (lastWeight>=weights_[j]-1.0e-7)
        throw CoinError("Weights too close together in OsiLink","infeasibility","OsiLink");
      lastWeight = weights_[j];
      double value = CoinMax(0.0,solution[iColumn]);
      sum += value;
      if (value>integerTolerance&&upper[iColumn]) {
        // Possibly due to scaling a fixed variable might slip through
        if (value>upper[iColumn]+1.0e-8) {
#ifdef OSI_DEBUG
	  printf("** Variable %d (%d) has value %g and upper bound of %g\n",
		 iColumn,j,value,upper[iColumn]);
#endif
        } 
	value = CoinMin(value,upper[iColumn]);
        weight += weights_[j]*value;
        if (firstNonZero<0)
          firstNonZero=j;
        lastNonZero=j;
      }
    }
    base += numberLinks_;
  }
  double valueInfeasibility;
  whichWay=1;
  whichWay_=1;
  if (lastNonZero-firstNonZero>=sosType_) {
    // find where to branch
    assert (sum>0.0);
    weight /= sum;
    valueInfeasibility = lastNonZero-firstNonZero+1;
    valueInfeasibility *= 0.5/((double) numberMembers_);
    //#define DISTANCE
#ifdef DISTANCE
    assert (sosType_==1); // code up
    /* may still be satisfied.
       For LOS type 2 we might wish to move coding around
       and keep initial info in model_ for speed
    */
    int iWhere;
    bool possible=false;
    for (iWhere=firstNonZero;iWhere<=lastNonZero;iWhere++) {
      if (fabs(weight-weights_[iWhere])<1.0e-8) {
	possible=true;
	break;
      }
    }
    if (possible) {
      // One could move some of this (+ arrays) into model_
      const CoinPackedMatrix * matrix = solver->getMatrixByCol();
      const double * element = matrix->getMutableElements();
      const int * row = matrix->getIndices();
      const CoinBigIndex * columnStart = matrix->getVectorStarts();
      const int * columnLength = matrix->getVectorLengths();
      const double * rowSolution = solver->getRowActivity();
      const double * rowLower = solver->getRowLower();
      const double * rowUpper = solver->getRowUpper();
      int numberRows = matrix->getNumRows();
      double * array = new double [numberRows];
      CoinZeroN(array,numberRows);
      int * which = new int [numberRows];
      int n=0;
      int base=numberLinks_*firstNonZero;
      for (j=firstNonZero;j<=lastNonZero;j++) {
	for (int k=0;k<numberLinks_;k++) {
	  int iColumn = members_[base+k];
	  double value = CoinMax(0.0,solution[iColumn]);
	  if (value>integerTolerance&&upper[iColumn]) {
	    value = CoinMin(value,upper[iColumn]);
	    for (int j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
	      int iRow = row[j];
	      double a = array[iRow];
	      if (a) {
		a += value*element[j];
		if (!a)
		  a = 1.0e-100;
	      } else {
		which[n++]=iRow;
		a=value*element[j];
		assert (a);
	      }
	      array[iRow]=a;
	    }
	  }
	}
	base += numberLinks_;
      }
      base=numberLinks_*iWhere;
      for (int k=0;k<numberLinks_;k++) {
	int iColumn = members_[base+k];
	const double value = 1.0;
	for (int j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  double a = array[iRow];
	  if (a) {
	    a -= value*element[j];
	    if (!a)
	      a = 1.0e-100;
	  } else {
	    which[n++]=iRow;
	    a=-value*element[j];
	    assert (a);
	  }
	  array[iRow]=a;
	}
      }
      for (j=0;j<n;j++) {
	int iRow = which[j];
	// moving to point will increase row solution by this
	double distance = array[iRow];
	if (distance>1.0e-8) {
	  if (distance+rowSolution[iRow]>rowUpper[iRow]+1.0e-8) {
	    possible=false;
	    break;
	  }
	} else if (distance<-1.0e-8) {
	  if (distance+rowSolution[iRow]<rowLower[iRow]-1.0e-8) {
	    possible=false;
	    break;
	  } 
	}
      }
      for (j=0;j<n;j++)
	array[which[j]]=0.0;
      delete [] array;
      delete [] which;
      if (possible) {
	valueInfeasibility=0.0;
	printf("possible %d %d %d\n",firstNonZero,lastNonZero,iWhere);
      }
    }
#endif
  } else {
    valueInfeasibility = 0.0; // satisfied
  }
  infeasibility_=valueInfeasibility;
  otherInfeasibility_=1.0-valueInfeasibility;
  return valueInfeasibility;
}

// This looks at solution and sets bounds to contain solution
double
OsiOldLink::feasibleRegion(OsiSolverInterface * solver, const OsiBranchingInformation * info) const
{
  int j;
  int firstNonZero=-1;
  int lastNonZero = -1;
  const double * solution = info->solution_;
  const double * upper = info->upper_;
  double integerTolerance = info->integerTolerance_;
  double weight = 0.0;
  double sum =0.0;

  int base=0;
  for (j=0;j<numberMembers_;j++) {
    for (int k=0;k<numberLinks_;k++) {
      int iColumn = members_[base+k];
      double value = CoinMax(0.0,solution[iColumn]);
      sum += value;
      if (value>integerTolerance&&upper[iColumn]) {
        weight += weights_[j]*value;
        if (firstNonZero<0)
          firstNonZero=j;
        lastNonZero=j;
      }
    }
    base += numberLinks_;
  }
#ifdef DISTANCE
  if (lastNonZero-firstNonZero>sosType_-1) {
    /* may still be satisfied.
       For LOS type 2 we might wish to move coding around
       and keep initial info in model_ for speed
    */
    int iWhere;
    bool possible=false;
    for (iWhere=firstNonZero;iWhere<=lastNonZero;iWhere++) {
      if (fabs(weight-weights_[iWhere])<1.0e-8) {
	possible=true;
	break;
      }
    }
    if (possible) {
      // One could move some of this (+ arrays) into model_
      const CoinPackedMatrix * matrix = solver->getMatrixByCol();
      const double * element = matrix->getMutableElements();
      const int * row = matrix->getIndices();
      const CoinBigIndex * columnStart = matrix->getVectorStarts();
      const int * columnLength = matrix->getVectorLengths();
      const double * rowSolution = solver->getRowActivity();
      const double * rowLower = solver->getRowLower();
      const double * rowUpper = solver->getRowUpper();
      int numberRows = matrix->getNumRows();
      double * array = new double [numberRows];
      CoinZeroN(array,numberRows);
      int * which = new int [numberRows];
      int n=0;
      int base=numberLinks_*firstNonZero;
      for (j=firstNonZero;j<=lastNonZero;j++) {
	for (int k=0;k<numberLinks_;k++) {
	  int iColumn = members_[base+k];
	  double value = CoinMax(0.0,solution[iColumn]);
	  if (value>integerTolerance&&upper[iColumn]) {
	    value = CoinMin(value,upper[iColumn]);
	    for (int j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
	      int iRow = row[j];
	      double a = array[iRow];
	      if (a) {
		a += value*element[j];
		if (!a)
		  a = 1.0e-100;
	      } else {
		which[n++]=iRow;
		a=value*element[j];
		assert (a);
	      }
	      array[iRow]=a;
	    }
	  }
	}
	base += numberLinks_;
      }
      base=numberLinks_*iWhere;
      for (int k=0;k<numberLinks_;k++) {
	int iColumn = members_[base+k];
	const double value = 1.0;
	for (int j=columnStart[iColumn];j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  double a = array[iRow];
	  if (a) {
	    a -= value*element[j];
	    if (!a)
	      a = 1.0e-100;
	  } else {
	    which[n++]=iRow;
	    a=-value*element[j];
	    assert (a);
	  }
	  array[iRow]=a;
	}
      }
      for (j=0;j<n;j++) {
	int iRow = which[j];
	// moving to point will increase row solution by this
	double distance = array[iRow];
	if (distance>1.0e-8) {
	  if (distance+rowSolution[iRow]>rowUpper[iRow]+1.0e-8) {
	    possible=false;
	    break;
	  }
	} else if (distance<-1.0e-8) {
	  if (distance+rowSolution[iRow]<rowLower[iRow]-1.0e-8) {
	    possible=false;
	    break;
	  } 
	}
      }
      for (j=0;j<n;j++)
	array[which[j]]=0.0;
      delete [] array;
      delete [] which;
      if (possible) {
	printf("possible feas region %d %d %d\n",firstNonZero,lastNonZero,iWhere);
	firstNonZero=iWhere;
	lastNonZero=iWhere;
      }
    }
  }
#else
  assert (lastNonZero-firstNonZero<sosType_) ;
#endif
  base=0;
  for (j=0;j<firstNonZero;j++) {
    for (int k=0;k<numberLinks_;k++) {
      int iColumn = members_[base+k];
      solver->setColUpper(iColumn,0.0);
    }
    base += numberLinks_;
  }
  // skip
  base += numberLinks_;
  for (j=lastNonZero+1;j<numberMembers_;j++) {
    for (int k=0;k<numberLinks_;k++) {
      int iColumn = members_[base+k];
      solver->setColUpper(iColumn,0.0);
    }
    base += numberLinks_;
  }
  // go to coding as in OsiSOS
  abort();
  return -1.0;
}

// Redoes data when sequence numbers change
void 
OsiOldLink::resetSequenceEtc(int numberColumns, const int * originalColumns)
{
  int n2=0;
  for (int j=0;j<numberMembers_*numberLinks_;j++) {
    int iColumn = members_[j];
    int i;
    for (i=0;i<numberColumns;i++) {
      if (originalColumns[i]==iColumn)
        break;
    }
    if (i<numberColumns) {
      members_[n2]=i;
      weights_[n2++]=weights_[j];
    }
  }
  if (n2<numberMembers_) {
    printf("** SOS number of members reduced from %d to %d!\n",numberMembers_,n2/numberLinks_);
    numberMembers_=n2/numberLinks_;
  }
}

// Creates a branching object
OsiBranchingObject * 
OsiOldLink::createBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) const
{
  int j;
  const double * solution = info->solution_;
  double tolerance = info->primalTolerance_;
  const double * upper = info->upper_;
  int firstNonFixed=-1;
  int lastNonFixed=-1;
  int firstNonZero=-1;
  int lastNonZero = -1;
  double weight = 0.0;
  double sum =0.0;
  int base=0;
  for (j=0;j<numberMembers_;j++) {
    for (int k=0;k<numberLinks_;k++) {
      int iColumn = members_[base+k];
      if (upper[iColumn]) {
        double value = CoinMax(0.0,solution[iColumn]);
        sum += value;
        if (firstNonFixed<0)
          firstNonFixed=j;
        lastNonFixed=j;
        if (value>tolerance) {
          weight += weights_[j]*value;
          if (firstNonZero<0)
            firstNonZero=j;
          lastNonZero=j;
        }
      }
    }
    base += numberLinks_;
  }
  assert (lastNonZero-firstNonZero>=sosType_) ;
  // find where to branch
  assert (sum>0.0);
  weight /= sum;
  int iWhere;
  double separator=0.0;
  for (iWhere=firstNonZero;iWhere<lastNonZero;iWhere++) 
    if (weight<weights_[iWhere+1])
      break;
  if (sosType_==1) {
    // SOS 1
    separator = 0.5 *(weights_[iWhere]+weights_[iWhere+1]);
  } else {
    // SOS 2
    if (iWhere==firstNonFixed)
      iWhere++;;
    if (iWhere==lastNonFixed-1)
      iWhere = lastNonFixed-2;
    separator = weights_[iWhere+1];
  }
  // create object
  OsiBranchingObject * branch;
  branch = new OsiOldLinkBranchingObject(solver,this,way,separator);
  return branch;
}
OsiOldLinkBranchingObject::OsiOldLinkBranchingObject()
  :OsiSOSBranchingObject()
{
}

// Useful constructor
OsiOldLinkBranchingObject::OsiOldLinkBranchingObject (OsiSolverInterface * solver,
					      const OsiOldLink * set,
					      int way ,
					      double separator)
  :OsiSOSBranchingObject(solver,set,way,separator)
{
}

// Copy constructor 
OsiOldLinkBranchingObject::OsiOldLinkBranchingObject ( const OsiOldLinkBranchingObject & rhs) :OsiSOSBranchingObject(rhs)
{
}

// Assignment operator 
OsiOldLinkBranchingObject & 
OsiOldLinkBranchingObject::operator=( const OsiOldLinkBranchingObject& rhs)
{
  if (this != &rhs) {
    OsiSOSBranchingObject::operator=(rhs);
  }
  return *this;
}
OsiBranchingObject * 
OsiOldLinkBranchingObject::clone() const
{ 
  return (new OsiOldLinkBranchingObject(*this));
}


// Destructor 
OsiOldLinkBranchingObject::~OsiOldLinkBranchingObject ()
{
}
double
OsiOldLinkBranchingObject::branch(OsiSolverInterface * solver)
{
  const OsiOldLink * set =
    dynamic_cast <const OsiOldLink *>(originalObject_) ;
  assert (set);
  int way = (!branchIndex_) ? (2*firstBranch_-1) : -(2*firstBranch_-1);
  branchIndex_++;
  int numberMembers = set->numberMembers();
  const int * which = set->members();
  const double * weights = set->weights();
  int numberLinks = set->numberLinks();
  //const double * lower = info->lower_;
  //const double * upper = solver->getColUpper();
  // *** for way - up means fix all those in down section
  if (way<0) {
    int i;
    for ( i=0;i<numberMembers;i++) {
      if (weights[i] > value_)
	break;
    }
    assert (i<numberMembers);
    int base=i*numberLinks;;
    for (;i<numberMembers;i++) {
      for (int k=0;k<numberLinks;k++) {
        int iColumn = which[base+k];
        solver->setColUpper(iColumn,0.0);
      }
      base += numberLinks;
    }
  } else {
    int i;
    int base=0;
    for ( i=0;i<numberMembers;i++) { 
      if (weights[i] >= value_) {
	break;
      } else {
        for (int k=0;k<numberLinks;k++) {
	  int iColumn = which[base+k];
          solver->setColUpper(iColumn,0.0);
        }
        base += numberLinks;
      }
    }
    assert (i<numberMembers);
  }
  return 0.0;
}
// Print what would happen  
void
OsiOldLinkBranchingObject::print(const OsiSolverInterface * solver)
{
  const OsiOldLink * set =
    dynamic_cast <const OsiOldLink *>(originalObject_) ;
  assert (set);
  int way = (!branchIndex_) ? (2*firstBranch_-1) : -(2*firstBranch_-1);
  int numberMembers = set->numberMembers();
  int numberLinks = set->numberLinks();
  const double * weights = set->weights();
  const int * which = set->members();
  const double * upper = solver->getColUpper();
  int first=numberMembers;
  int last=-1;
  int numberFixed=0;
  int numberOther=0;
  int i;
  int base=0;
  for ( i=0;i<numberMembers;i++) {
    for (int k=0;k<numberLinks;k++) {
      int iColumn = which[base+k];
      double bound = upper[iColumn];
      if (bound) {
        first = CoinMin(first,i);
        last = CoinMax(last,i);
      }
    }
    base += numberLinks;
  }
  // *** for way - up means fix all those in down section
  base=0;
  if (way<0) {
    printf("SOS Down");
    for ( i=0;i<numberMembers;i++) {
      if (weights[i] > value_) 
	break;
      for (int k=0;k<numberLinks;k++) {
        int iColumn = which[base+k];
        double bound = upper[iColumn];
	if (bound)
          numberOther++;
      }
      base += numberLinks;
    }
    assert (i<numberMembers);
    for (;i<numberMembers;i++) {
      for (int k=0;k<numberLinks;k++) {
        int iColumn = which[base+k];
        double bound = upper[iColumn];
        if (bound)
          numberFixed++;
      }
      base += numberLinks;
    }
  } else {
    printf("SOS Up");
    for ( i=0;i<numberMembers;i++) {
      if (weights[i] >= value_)
	break;
      for (int k=0;k<numberLinks;k++) {
        int iColumn = which[base+k];
        double bound = upper[iColumn];
	if (bound)
          numberFixed++;
      }
      base += numberLinks;
    }
    assert (i<numberMembers);
    for (;i<numberMembers;i++) {
      for (int k=0;k<numberLinks;k++) {
        int iColumn = which[base+k];
        double bound = upper[iColumn];
        if (bound)
          numberOther++;
      }
      base += numberLinks;
    }
  }
  assert ((numberFixed%numberLinks)==0);
  assert ((numberOther%numberLinks)==0);
  printf(" - at %g, free range %d (%g) => %d (%g), %d would be fixed, %d other way\n",
	 value_,first,weights[first],last,weights[last],numberFixed/numberLinks,
         numberOther/numberLinks);
}
// Default Constructor 
OsiBiLinear::OsiBiLinear ()
  : OsiObject2(),
    coefficient_(0.0),
    xMeshSize_(0.0),
    yMeshSize_(0.0),
    xSatisfied_(1.0e-6),
    ySatisfied_(1.0e-6),
    xySatisfied_(1.0e-6),
    xyBranchValue_(0.0),
    xColumn_(-1),
    yColumn_(-1),
    firstLambda_(-1),
    branchingStrategy_(0),
    xRow_(-1),
    yRow_(-1),
    xyRow_(-1),
    convexity_(-1),
    chosen_(-1)
{
}

// Useful constructor
OsiBiLinear::OsiBiLinear (OsiSolverInterface * solver, int xColumn,
			  int yColumn, int xyRow, double coefficient,
			  double xMesh, double yMesh,
			  int numberExistingObjects,const OsiObject ** objects )
  : OsiObject2(),
    coefficient_(coefficient),
    xMeshSize_(xMesh),
    yMeshSize_(yMesh),
    xSatisfied_(1.0e-6),
    ySatisfied_(1.0e-6),
    xySatisfied_(1.0e-6),
    xyBranchValue_(0.0),
    xColumn_(xColumn),
    yColumn_(yColumn),
    firstLambda_(-1),
    branchingStrategy_(0),
    xRow_(-1),
    yRow_(-1),
    xyRow_(xyRow),
    convexity_(-1),
    chosen_(-1)
{
  double columnLower[4];
  double columnUpper[4];
  double objective[4];
  double rowLower[3];
  double rowUpper[3];
  CoinBigIndex starts[5];
  int index[16];
  double element[16];
  int i;
  starts[0]=0;
  // rows
  int numberRows = solver->getNumRows();
  // convexity
  rowLower[0]=1.0;
  rowUpper[0]=1.0;
  convexity_ = numberRows;
  starts[1]=0;
  // x
  rowLower[1]=0.0;
  rowUpper[1]=0.0;
  index[0]=xColumn_;
  element[0]=-1.0;
  xRow_ = numberRows+1;
  starts[2]=1;
  int nAdd=2;
  if (xColumn!=yColumn) {
    rowLower[2]=0.0;
    rowUpper[2]=0.0;
    index[1]=yColumn;
    element[1]=-1.0;
    nAdd=3;
    yRow_ = numberRows+2;
    starts[3]=2;
  } else {
    yRow_=-1;
    branchingStrategy_=1;
  }
  solver->addRows(nAdd,starts,index,element,rowLower,rowUpper);
  int n=0;
  // order is LxLy, LxUy, UxLy and UxUy
  firstLambda_ = solver->getNumCols();
  // bit sloppy as theoretically could be infeasible but otherwise need to do more work
  double xB[2];
  double yB[2];
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  xB[0]=lower[xColumn_];
  xB[1]=upper[xColumn_];
  // adjust
  double distance;
  int steps;
  distance = xB[1]-xB[0];
  steps = (int) ((distance+0.5*xMeshSize_)/xMeshSize_);
  distance = xB[0]+xMeshSize_*steps;
  if (fabs(xB[1]-distance)>1.0e-9) {
    printf("bad x mesh %g %g %g -> %g\n",xB[0],xMeshSize_,xB[1],distance);
    xB[1]=distance;
    solver->setColUpper(xColumn_,distance);
  }
  yB[0]=lower[yColumn_];
  yB[1]=upper[yColumn_];
  distance = yB[1]-yB[0];
  steps = (int) ((distance+0.5*yMeshSize_)/yMeshSize_);
  distance = yB[0]+yMeshSize_*steps;
  if (fabs(yB[1]-distance)>1.0e-9) {
    printf("bad y mesh %g %g %g -> %g\n",yB[0],yMeshSize_,yB[1],distance);
    yB[1]=distance;
    solver->setColUpper(yColumn_,distance);
  }
  for (i=0;i<4;i++) {
    double x = (i<2) ? xB[0] : xB[1];
    double y = ((i&1)==0) ? yB[0] : yB[1];
    columnLower[i]=0.0;
    columnUpper[i]=COIN_DBL_MAX;
    objective[i]=0.0;
    double value;
    // xy
    value=coefficient_*x*y;
    if (fabs(value)<1.0e-12)
      value = 1.0e-12;
    element[n]=value;
    index[n++]=xyRow_;
    // convexity
    value=1.0;
    element[n]=value;
    index[n++]=0+numberRows;
    // x
    value=coefficient_*x;
    if (fabs(value)<1.0e-12)
      value = 1.0e-12;
    element[n]=value;
    index[n++]=1+numberRows;
    if (xColumn_!=yColumn_) {
      // y
      value=coefficient_*y;
      if (fabs(value)<1.0e-12)
      value = 1.0e-12;
      element[n]=value;
      index[n++]=2+numberRows;
    }
    starts[i+1]=n;
  }
  solver->addCols(4,starts,index,element,columnLower,columnUpper,objective);
  // At least one has to have a mesh
  if (!xMeshSize_&&(!yMeshSize_||yRow_<0)) {
    printf("one of x and y must have a mesh size\n");
    abort();
  } else if (yRow_>=0) {
    if (!xMeshSize_)
      branchingStrategy_ = 2;
    else if (!yMeshSize_)
      branchingStrategy_ = 1;
  }
  // Now add constraints to link in x and or y to existing ones.
  bool xDone=false;
  bool yDone=false;
  // order is LxLy, LxUy, UxLy and UxUy
  for (i=numberExistingObjects-1;i>=0;i--) {
    const OsiObject * obj = objects[i];
    const OsiBiLinear * obj2 =
      dynamic_cast <const OsiBiLinear *>(obj) ;
    if (obj2) {
      if (xColumn_==obj2->xColumn_&&!xDone) {
	// make sure y equal
	double rhs=0.0;
	CoinBigIndex starts[2];
	int index[4];
	double element[4]= {1.0,1.0,-1.0,-1.0};
	starts[0]=0;
	starts[1]=4;
	index[0]=firstLambda_+0;
	index[1]=firstLambda_+1;
	index[2]=obj2->firstLambda_+0;
	index[3]=obj2->firstLambda_+1;
	solver->addRows(1,starts,index,element,&rhs,&rhs);
	xDone=true;
      }
      if (yColumn_==obj2->yColumn_&&yRow_>=0&&!yDone) {
	// make sure x equal
	double rhs=0.0;
	CoinBigIndex starts[2];
	int index[4];
	double element[4]= {1.0,1.0,-1.0,-1.0};
	starts[0]=0;
	starts[1]=4;
	index[0]=firstLambda_+0;
	index[1]=firstLambda_+2;
	index[2]=obj2->firstLambda_+0;
	index[3]=obj2->firstLambda_+2;
	solver->addRows(1,starts,index,element,&rhs,&rhs);
	yDone=true;
      }
    }
  }
}

// Copy constructor 
OsiBiLinear::OsiBiLinear ( const OsiBiLinear & rhs)
  :OsiObject2(rhs),
   coefficient_(rhs.coefficient_),
   xMeshSize_(rhs.xMeshSize_),
   yMeshSize_(rhs.yMeshSize_),
   xSatisfied_(rhs.xSatisfied_),
   ySatisfied_(rhs.ySatisfied_),
   xySatisfied_(rhs.xySatisfied_),
   xyBranchValue_(rhs.xyBranchValue_),
   xColumn_(rhs.xColumn_),
   yColumn_(rhs.yColumn_),
   firstLambda_(rhs.firstLambda_),
   branchingStrategy_(rhs.branchingStrategy_),
   xRow_(rhs.xRow_),
   yRow_(rhs.yRow_),
   xyRow_(rhs.xyRow_),
   convexity_(rhs.convexity_),
   chosen_(rhs.chosen_)
{
}

// Clone
OsiObject *
OsiBiLinear::clone() const
{
  return new OsiBiLinear(*this);
}

// Assignment operator 
OsiBiLinear & 
OsiBiLinear::operator=( const OsiBiLinear& rhs)
{
  if (this!=&rhs) {
    OsiObject2::operator=(rhs);
    coefficient_ = rhs.coefficient_;
    xMeshSize_ = rhs.xMeshSize_;
    yMeshSize_ = rhs.yMeshSize_;
    xSatisfied_ = rhs.xSatisfied_;
    ySatisfied_ = rhs.ySatisfied_;
    xySatisfied_ = rhs.xySatisfied_;
    xyBranchValue_ = rhs.xyBranchValue_;
    xColumn_ = rhs.xColumn_;
    yColumn_ = rhs.yColumn_;
    firstLambda_ = rhs.firstLambda_;
    branchingStrategy_ = rhs.branchingStrategy_;
    xRow_ = rhs.xRow_;
    yRow_ = rhs.yRow_;
    xyRow_ = rhs.xyRow_;
    convexity_ = rhs.convexity_;
    chosen_ = rhs.chosen_;
  }
  return *this;
}

// Destructor 
OsiBiLinear::~OsiBiLinear ()
{
}

// Infeasibility - large is 0.5
double 
OsiBiLinear::infeasibility(const OsiBranchingInformation * info,int & whichWay) const
{
  // order is LxLy, LxUy, UxLy and UxUy
  double xB[2];
  double yB[2];
  xB[0]=info->lower_[xColumn_];
  xB[1]=info->upper_[xColumn_];
  yB[0]=info->lower_[yColumn_];
  yB[1]=info->upper_[yColumn_];
  double x = info->solution_[xColumn_];
  x = CoinMax(x,xB[0]);
  x = CoinMin(x,xB[1]);
  double y = info->solution_[yColumn_];
  y = CoinMax(y,yB[0]);
  y = CoinMin(y,yB[1]);
  int j;
#ifndef NDEBUG
  double xLambda = 0.0;
  double yLambda = 0.0;
  for (j=0;j<4;j++) {
    int iX = j>>1;
    int iY = j&1;
    xLambda += xB[iX]*info->solution_[firstLambda_+j];
    yLambda += yB[iY]*info->solution_[firstLambda_+j];
  }
  assert (fabs(x-xLambda)<1.0e-4);
  assert (fabs(y-yLambda)<1.0e-4);
#endif
  // If x or y not satisfied then branch on that
  double distance;
  int steps;
  bool xSatisfied;
  double xNew;
  distance = x-xB[0];
  if (xMeshSize_) {
    steps = (int) ((distance+0.5*xMeshSize_)/xMeshSize_);
    xNew = xB[0]+steps*xMeshSize_;
    assert (xNew<=xB[1]+1.0e-5);
    xSatisfied =  (fabs(xNew-x)<xSatisfied_);
  } else {
    xSatisfied=true;
  }
  bool ySatisfied;
  double yNew;
  distance = y-yB[0];
  if (yMeshSize_) {
    steps = (int) ((distance+0.5*yMeshSize_)/yMeshSize_);
    yNew = yB[0]+steps*yMeshSize_;
    assert (yNew<=yB[1]+1.0e-5);
    ySatisfied =  (fabs(yNew-y)<ySatisfied_)||!y;
  } else {
    ySatisfied=true;
  }
  /* There are several possibilities
     1 - one or both are unsatisfied and branching strategy tells us what to do
     2 - both are unsatisfied and branching strategy is 0
     3 - both are satisfied but xy is not
         3a one has bounds within satisfied_ - other does not
	 (or neither have but branching strategy tells us what to do)
	 3b neither do - and branching strategy does not tell us
	 3c both do - treat as feasible knowing another copy of object will fix
     4 - both are satisfied and xy is satisfied - as 3c
  */
  chosen_=-1;
  xyBranchValue_=COIN_DBL_MAX;
  whichWay_=0;
  if ( !xSatisfied) {
    if (!ySatisfied) {
      if (branchingStrategy_==0) {
	// If pseudo shadow prices then see what would happen
	if (info->defaultDual_>=0.0) {
	  // need coding here
	  if (fabs(x-xNew)>fabs(y-yNew)) {
	    chosen_=0;
	    xyBranchValue_=x;
	  } else {
	    chosen_=1;
	    xyBranchValue_=y;
	  }
	} else {
	  if (fabs(x-xNew)>fabs(y-yNew)) {
	    chosen_=0;
	    xyBranchValue_=x;
	  } else {
	    chosen_=1;
	    xyBranchValue_=y;
	  }
	}
      } else if (branchingStrategy_==1) {
	chosen_=0;
	xyBranchValue_=x;
      } else {
	chosen_=1;
	xyBranchValue_=y;
      }
    } else {
      // y satisfied
      chosen_=0;
      xyBranchValue_=x;
    }
  } else {
    // x satisfied
    if (!ySatisfied) {
      chosen_=1;
      xyBranchValue_=y;
    } else {
      /*
	3 - both are satisfied but xy is not
	 3a one has bounds within satisfied_ - other does not
	 (or neither have but branching strategy tells us what to do)
	 3b neither do - and branching strategy does not tell us
	 3c both do - treat as feasible knowing another copy of object will fix
        4 - both are satisfied and xy is satisfied - as 3c
      */
      double xyTrue = x*y;
      double xyLambda = 0.0;
      for (j=0;j<4;j++) {
	int iX = j>>1;
	int iY = j&1;
	xyLambda += xB[iX]*yB[iY]*info->solution_[firstLambda_+j];
      }
      if (fabs(xyLambda-xyTrue)<xySatisfied_) {
	// satisfied
      } else {
	if (xB[1]-xB[0]>=xSatisfied_&&xMeshSize_) {
	  if (yB[1]-yB[0]>=ySatisfied_&&yMeshSize_) {
	    if (branchingStrategy_==0) {
	      // If pseudo shadow prices then see what would happen
	      if (info->defaultDual_>=0.0) {
		// need coding here
		if (xB[1]-xB[0]>yB[1]-yB[0]) {
		  chosen_=0;
		  xyBranchValue_=0.5*(xB[0]+xB[1]);
		} else {
		  chosen_=1;
		  xyBranchValue_=0.5*(yB[0]+yB[1]);
		}
	      } else {
		if (xB[1]-xB[0]>yB[1]-yB[0]) {
		  chosen_=0;
		  xyBranchValue_=0.5*(xB[0]+xB[1]);
		} else {
		  chosen_=1;
		  xyBranchValue_=0.5*(yB[0]+yB[1]);
		}
	      }
	    } else if (branchingStrategy_==1) {
	      chosen_=0;
	      xyBranchValue_=0.5*(xB[0]+xB[1]);
	    } else {
	      chosen_=1;
	      xyBranchValue_=0.5*(yB[0]+yB[1]);
	    }
	  } else {
	    // y satisfied
	    chosen_=0;
	    xyBranchValue_=0.5*(xB[0]+xB[1]);
	  }
	} else if (yB[1]-yB[0]>=ySatisfied_&&yMeshSize_) {
	  chosen_=1;
	  xyBranchValue_=0.5*(yB[0]+yB[1]);
	} else {
	  // treat as satisfied
	}
      }
    }
  }
  if (chosen_==-1) {
    infeasibility_=0.0;
  } else if (chosen_==0) {
    infeasibility_ = CoinMax(fabs(xyBranchValue_-x),1.0e-12);
    assert (xyBranchValue_>=info->lower_[xColumn_]&&xyBranchValue_<=info->upper_[xColumn_]);
  } else {
    infeasibility_ = CoinMax(fabs(xyBranchValue_-y),1.0e-12);
    assert (xyBranchValue_>=info->lower_[yColumn_]&&xyBranchValue_<=info->upper_[yColumn_]);
  }
  if (info->defaultDual_<0.0) {
    // not using pseudo shadow prices
    otherInfeasibility_ = 1.0-infeasibility_;
  } else {
    abort();
  }
  whichWay=whichWay_;
  return infeasibility_;
}

// This looks at solution and sets bounds to contain solution
double
OsiBiLinear::feasibleRegion(OsiSolverInterface * solver, const OsiBranchingInformation * info) const
{
  // order is LxLy, LxUy, UxLy and UxUy
  double xB[2];
  double yB[2];
  xB[0]=info->lower_[xColumn_];
  xB[1]=info->upper_[xColumn_];
  yB[0]=info->lower_[yColumn_];
  yB[1]=info->upper_[yColumn_];
  double x = info->solution_[xColumn_];
  double y = info->solution_[yColumn_];
  int j;
#ifndef NDEBUG
  double xLambda = 0.0;
  double yLambda = 0.0;
  for (j=0;j<4;j++) {
    int iX = j>>1;
    int iY = j&1;
    xLambda += xB[iX]*info->solution_[firstLambda_+j];
    yLambda += yB[iY]*info->solution_[firstLambda_+j];
  }
  if (fabs(x-xLambda)>1.0e-4||
      fabs(y-yLambda)>1.0e-4)
    printf("feasibleregion x %d given %g lambda %g y %d given %g lambda %g\n",
	   xColumn_,x,xLambda,
	   yColumn_,y,yLambda);
#endif
  double infeasibility=0.0;
  double distance;
  int steps;
  double xNew=x;
  distance = x-xB[0];
  if (xMeshSize_) {
    steps = (int) ((distance+0.5*xMeshSize_)/xMeshSize_);
    xNew = xB[0]+steps*xMeshSize_;
    assert (xNew<=xB[1]+1.0e-5);
    infeasibility +=  fabs(xNew-x);
    solver->setColLower(xColumn_,xNew);
    solver->setColUpper(xColumn_,xNew);
  }
  double yNew=y;
  distance = y-yB[0];
  if (yMeshSize_) {
    steps = (int) ((distance+0.5*yMeshSize_)/yMeshSize_);
    yNew = yB[0]+steps*yMeshSize_;
    assert (yNew<=yB[1]+1.0e-5);
    infeasibility +=  fabs(yNew-y);
    solver->setColLower(yColumn_,yNew);
    solver->setColUpper(yColumn_,yNew);
  }
  double xyTrue = xNew*yNew;
  double xyLambda = 0.0;
  for (j=0;j<4;j++) {
    int iX = j>>1;
    int iY = j&1;
    xyLambda += xB[iX]*yB[iY]*info->solution_[firstLambda_+j];
  }
  infeasibility += fabs(xyTrue-xyLambda);
  return infeasibility;
}

// Redoes data when sequence numbers change
void 
OsiBiLinear::resetSequenceEtc(int numberColumns, const int * originalColumns)
{
  int i;
  for (i=0;i<numberColumns;i++) {
    if (originalColumns[i]==firstLambda_)
      break;
  }
  if (i<numberColumns) {
    firstLambda_ = i;
    for (int j=0;j<4;j++) {
      assert (originalColumns[j+i]-firstLambda_==j);
    }
  } else {
    printf("lost set\n");
    abort();
  }
  // rows will be out anyway
  abort();
}

// Creates a branching object
OsiBranchingObject * 
OsiBiLinear::createBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) const
{
  // create object
  OsiBranchingObject * branch;
  assert (chosen_==0||chosen_==1);
  if (chosen_==0) 
    assert (xyBranchValue_>=info->lower_[xColumn_]&&xyBranchValue_<=info->upper_[xColumn_]);
  else
    assert (xyBranchValue_>=info->lower_[yColumn_]&&xyBranchValue_<=info->upper_[yColumn_]);
  branch = new OsiBiLinearBranchingObject(solver,this,way,xyBranchValue_,chosen_);
  return branch;
}
// Does work of branching
void 
OsiBiLinear::newBounds(OsiSolverInterface * solver, int way, short xOrY, double separator) const
{
  int iColumn;
  double mesh;
  if (xOrY==0) {
    iColumn=xColumn_;
    mesh=xMeshSize_;
  } else {
    iColumn=yColumn_;
    mesh=yMeshSize_;
  }
  double lower = solver->getColLower()[iColumn];
  double distance;
  int steps;
  double zNew=separator;
  distance = separator-lower;
  assert (mesh);
  steps = (int) ((distance+0.5*mesh)/mesh);
  zNew = lower+steps*mesh;
  assert (zNew<=solver->getColUpper()[iColumn]+1.0e-5);
#ifndef NDEBUG
    double oldUpper = solver->getColUpper()[iColumn] ;
    double oldLower = solver->getColLower()[iColumn] ;
#endif
  if (way<0) {
    if (zNew>separator)
      zNew -= mesh;
#ifndef NDEBUG
    double oldUpper = solver->getColUpper()[iColumn] ;
    assert (oldUpper>zNew-1.0e-8);
    if (oldUpper<zNew+1.0e-8)
      printf("null change on columnUpper %d - bounds %g,%g\n",iColumn,oldLower,oldUpper);
#endif
    solver->setColUpper(iColumn,zNew);
  } else {
    if (zNew<separator)
      zNew += mesh;
#ifndef NDEBUG
    double oldLower = solver->getColLower()[iColumn] ;
    assert (oldLower<zNew+1.0e-8);
    if (oldLower>zNew-1.0e-8)
      printf("null change on columnLower %d - bounds %g,%g\n",iColumn,oldLower,oldUpper);
#endif
    solver->setColLower(iColumn,zNew);
  }
#if 0
  // always free up lambda
  for (int i=firstLambda_;i<firstLambda_+4;i++) {
    solver->setColLower(i,0.0);
    solver->setColUpper(i,COIN_DBL_MAX);
  }
#endif
}
// Updates coefficients
void 
OsiBiLinear::updateCoefficients(const double * lower, const double * upper,
				CoinPackedMatrix * matrix, CoinWarmStartBasis * basis) const
{
  double * element = matrix->getMutableElements();
  const int * row = matrix->getIndices();
  const CoinBigIndex * columnStart = matrix->getVectorStarts();
  //const int * columnLength = matrix->getVectorLengths();
  // order is LxLy, LxUy, UxLy and UxUy
  double xB[2];
  double yB[2];
  xB[0]=lower[xColumn_];
  xB[1]=upper[xColumn_];
  yB[0]=lower[yColumn_];
  yB[1]=upper[yColumn_];
  //printf("x %d (%g,%g) y %d (%g,%g)\n",
  // xColumn_,xB[0],xB[1],
  // yColumn_,yB[0],yB[1]);
  CoinWarmStartBasis::Status status[4];
  int numStruct = basis->getNumStructural()-firstLambda_;
  for (int j=0;j<4;j++) {
    status[j]=(j<numStruct) ? basis->getStructStatus(j+firstLambda_) : CoinWarmStartBasis::atLowerBound;
    int iX = j>>1;
    double x = xB[iX];
    int iY = j&1;
    double y = yB[iY];
    CoinBigIndex k = columnStart[j+firstLambda_];
    double value;
    // xy
    value=coefficient_*x*y;
    assert (row[k]==xyRow_);
#if BI_PRINT > 1
    printf("j %d xy (%d,%d) coeff from %g to %g\n",j,xColumn_,yColumn_,element[k],value);
#endif
    element[k++]=value;
    // convexity
    assert (row[k]==convexity_);
    k++;
    // x
    value=coefficient_*x;
#if BI_PRINT > 1
    printf("j %d x (%d) coeff from %g to %g\n",j,xColumn_,element[k],value);
#endif
    assert (row[k]==xRow_);
    element[k++]=value;
    if (yRow_>=0) {
      // y
      value=coefficient_*y;
#if BI_PRINT > 1
      printf("j %d y (%d) coeff from %g to %g\n",j,yColumn_,element[k],value);
#endif
      assert (row[k]==yRow_);
      element[k++]=value;
    }
  }
  
  if (xB[0]==xB[1]) {
    if (yB[0]==yB[1]) {
      // only one basic
      bool first=true;
      for (int j=0;j<4;j++) {
	if (status[j]==CoinWarmStartBasis::basic) {
	  if (first) {
	    first=false;
	  } else {
	    basis->setStructStatus(j+firstLambda_,CoinWarmStartBasis::atLowerBound);
#if BI_PRINT
	    printf("zapping %d (x=%d,y=%d)\n",j,xColumn_,yColumn_);
#endif
	  }
	}
      }
    } else {
      if (status[0]==CoinWarmStartBasis::basic&&
	  status[2]==CoinWarmStartBasis::basic) {
	basis->setStructStatus(2+firstLambda_,CoinWarmStartBasis::atLowerBound);
#if BI_PRINT
	printf("zapping %d (x=%d,y=%d)\n",2,xColumn_,yColumn_);
#endif
      }
      if (status[1]==CoinWarmStartBasis::basic&&
	  status[3]==CoinWarmStartBasis::basic) {
	basis->setStructStatus(3+firstLambda_,CoinWarmStartBasis::atLowerBound);
#if BI_PRINT
	printf("zapping %d (x=%d,y=%d)\n",3,xColumn_,yColumn_);
#endif
      }
    }
  } else if (yB[0]==yB[1]) {
    if (status[0]==CoinWarmStartBasis::basic&&
	status[1]==CoinWarmStartBasis::basic) {
      basis->setStructStatus(1+firstLambda_,CoinWarmStartBasis::atLowerBound);
#if BI_PRINT
      printf("zapping %d (x=%d,y=%d)\n",1,xColumn_,yColumn_);
#endif
    }
    if (status[2]==CoinWarmStartBasis::basic&&
	status[3]==CoinWarmStartBasis::basic) {
      basis->setStructStatus(3+firstLambda_,CoinWarmStartBasis::atLowerBound);
#if BI_PRINT
      printf("zapping %d (x=%d,y=%d)\n",3,xColumn_,yColumn_);
#endif
    }
  }
}
// This does NOT set mutable stuff
double 
OsiBiLinear::checkInfeasibility(const OsiBranchingInformation * info) const
{
  int way;
  double saveInfeasibility = infeasibility_;
  int saveWhichWay = whichWay_;
  double saveXyBranchValue = xyBranchValue_;
  short saveChosen = chosen_;
  double value = infeasibility(info,way);
  infeasibility_ = saveInfeasibility;
  whichWay_ = saveWhichWay;
  xyBranchValue_ = saveXyBranchValue;
  chosen_ = saveChosen;
  return value;
}
OsiBiLinearBranchingObject::OsiBiLinearBranchingObject()
  :OsiTwoWayBranchingObject(),
   chosen_(0)
{
}

// Useful constructor
OsiBiLinearBranchingObject::OsiBiLinearBranchingObject (OsiSolverInterface * solver,
							const OsiBiLinear * set,
							int way ,
							double separator,
							int chosen)
  :OsiTwoWayBranchingObject(solver,set,way,separator),
   chosen_(chosen)
{
  assert (chosen_>=0&&chosen_<2);
}

// Copy constructor 
OsiBiLinearBranchingObject::OsiBiLinearBranchingObject ( const OsiBiLinearBranchingObject & rhs) 
  :OsiTwoWayBranchingObject(rhs),
   chosen_(rhs.chosen_)
{
}

// Assignment operator 
OsiBiLinearBranchingObject & 
OsiBiLinearBranchingObject::operator=( const OsiBiLinearBranchingObject& rhs)
{
  if (this != &rhs) {
    OsiTwoWayBranchingObject::operator=(rhs);
    chosen_ = rhs.chosen_;
  }
  return *this;
}
OsiBranchingObject * 
OsiBiLinearBranchingObject::clone() const
{ 
  return (new OsiBiLinearBranchingObject(*this));
}


// Destructor 
OsiBiLinearBranchingObject::~OsiBiLinearBranchingObject ()
{
}
double
OsiBiLinearBranchingObject::branch(OsiSolverInterface * solver)
{
  const OsiBiLinear * set =
    dynamic_cast <const OsiBiLinear *>(originalObject_) ;
  assert (set);
  int way = (!branchIndex_) ? (2*firstBranch_-1) : -(2*firstBranch_-1);
  branchIndex_++;
  set->newBounds(solver, way, chosen_, value_);
  return 0.0;
}
// Print what would happen  
void
OsiBiLinearBranchingObject::print(const OsiSolverInterface * solver)
{
  const OsiBiLinear * set =
    dynamic_cast <const OsiBiLinear *>(originalObject_) ;
  assert (set);
  int way = (!branchIndex_) ? (2*firstBranch_-1) : -(2*firstBranch_-1);
  int iColumn = (chosen_==1) ? set->xColumn() : set->yColumn();
  printf("OsiBiLinear would branch %s on %c variable %d from value %g\n",
	 (way<0) ? "down" : "up",
	 (chosen_==0) ? 'X' : 'Y', iColumn, value_);
}
/** Default Constructor

  Equivalent to an unspecified binary variable.
*/
OsiSimpleFixedInteger::OsiSimpleFixedInteger ()
  : OsiSimpleInteger()
{
}

/** Useful constructor

  Loads actual upper & lower bounds for the specified variable.
*/
OsiSimpleFixedInteger::OsiSimpleFixedInteger (const OsiSolverInterface * solver, int iColumn)
  : OsiSimpleInteger(solver,iColumn)
{
}

  
// Useful constructor - passed solver index and original bounds
OsiSimpleFixedInteger::OsiSimpleFixedInteger ( int iColumn, double lower, double upper)
  : OsiSimpleInteger(iColumn,lower,upper)
{
}

// Useful constructor - passed simple integer
OsiSimpleFixedInteger::OsiSimpleFixedInteger ( const OsiSimpleInteger &rhs)
  : OsiSimpleInteger(rhs)
{
}

// Copy constructor 
OsiSimpleFixedInteger::OsiSimpleFixedInteger ( const OsiSimpleFixedInteger & rhs)
  :OsiSimpleInteger(rhs)

{
}

// Clone
OsiObject *
OsiSimpleFixedInteger::clone() const
{
  return new OsiSimpleFixedInteger(*this);
}

// Assignment operator 
OsiSimpleFixedInteger & 
OsiSimpleFixedInteger::operator=( const OsiSimpleFixedInteger& rhs)
{
  if (this!=&rhs) {
    OsiSimpleInteger::operator=(rhs);
  }
  return *this;
}

// Destructor 
OsiSimpleFixedInteger::~OsiSimpleFixedInteger ()
{
}
// Infeasibility - large is 0.5
double 
OsiSimpleFixedInteger::infeasibility(const OsiBranchingInformation * info, int & whichWay) const
{
  double value = info->solution_[columnNumber_];
  value = CoinMax(value, info->lower_[columnNumber_]);
  value = CoinMin(value, info->upper_[columnNumber_]);
  double nearest = floor(value+(1.0-0.5));
  if (nearest>value) { 
    whichWay=1;
  } else {
    whichWay=0;
  }
  infeasibility_ = fabs(value-nearest);
  bool satisfied=false;
  if (infeasibility_<=info->integerTolerance_) {
    otherInfeasibility_ = 1.0;
    satisfied=true;
    if (info->lower_[columnNumber_]!=info->upper_[columnNumber_])
      infeasibility_ = 1.0e-5;
    else
      infeasibility_ = 0.0;
  } else if (info->defaultDual_<0.0) {
    otherInfeasibility_ = 1.0-infeasibility_;
  } else {
    const double * pi = info->pi_;
    const double * activity = info->rowActivity_;
    const double * lower = info->rowLower_;
    const double * upper = info->rowUpper_;
    const double * element = info->elementByColumn_;
    const int * row = info->row_;
    const CoinBigIndex * columnStart = info->columnStart_;
    const int * columnLength = info->columnLength_;
    double direction = info->direction_;
    double downMovement = value - floor(value);
    double upMovement = 1.0-downMovement;
    double valueP = info->objective_[columnNumber_]*direction;
    CoinBigIndex start = columnStart[columnNumber_];
    CoinBigIndex end = start + columnLength[columnNumber_];
    double upEstimate = 0.0;
    double downEstimate = 0.0;
    if (valueP>0.0)
      upEstimate = valueP*upMovement;
    else
      downEstimate -= valueP*downMovement;
    double tolerance = info->primalTolerance_;
    for (CoinBigIndex j=start;j<end;j++) {
      int iRow = row[j];
      if (lower[iRow]<-1.0e20) 
	assert (pi[iRow]<=1.0e-4);
      if (upper[iRow]>1.0e20) 
	assert (pi[iRow]>=-1.0e-4);
      valueP = pi[iRow]*direction;
      double el2 = element[j];
      double value2 = valueP*el2;
      double u=0.0;
      double d=0.0;
      if (value2>0.0)
	u = value2;
      else
	d = -value2;
      // if up makes infeasible then make at least default
      double newUp = activity[iRow] + upMovement*el2;
      if (newUp>upper[iRow]+tolerance||newUp<lower[iRow]-tolerance)
	u = CoinMax(u,info->defaultDual_);
      upEstimate += u*upMovement;
      // if down makes infeasible then make at least default
      double newDown = activity[iRow] - downMovement*el2;
      if (newDown>upper[iRow]+tolerance||newDown<lower[iRow]-tolerance)
	d = CoinMax(d,info->defaultDual_);
      downEstimate += d*downMovement;
    }
    if (downEstimate>=upEstimate) {
      infeasibility_ = CoinMax(1.0e-12,upEstimate);
      otherInfeasibility_ = CoinMax(1.0e-12,downEstimate);
      whichWay = 1;
    } else {
      infeasibility_ = CoinMax(1.0e-12,downEstimate);
      otherInfeasibility_ = CoinMax(1.0e-12,upEstimate);
      whichWay = 0;
    }
  }
  if (preferredWay_>=0&&!satisfied)
    whichWay = preferredWay_;
  whichWay_=whichWay;
  return infeasibility_;
}
// Creates a branching object
OsiBranchingObject * 
OsiSimpleFixedInteger::createBranch(OsiSolverInterface * solver, const OsiBranchingInformation * info, int way) const 
{
  double value = info->solution_[columnNumber_];
  value = CoinMax(value, info->lower_[columnNumber_]);
  value = CoinMin(value, info->upper_[columnNumber_]);
  assert (info->upper_[columnNumber_]>info->lower_[columnNumber_]);
  double nearest = floor(value+0.5);
  double integerTolerance = info->integerTolerance_;
  if (fabs(value-nearest)<integerTolerance) {
    // adjust value
    if (nearest!=info->upper_[columnNumber_])
      value = nearest+2.0*integerTolerance;
    else
      value = nearest-2.0*integerTolerance;
  }
  OsiBranchingObject * branch = new OsiIntegerBranchingObject(solver,this,way,
					     value);
  return branch;
}

