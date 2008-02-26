// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CbcConfig.h"

#include <string>
//#define CBC_DEBUG 1
//#define CHECK_CUT_COUNTS
//#define CHECK_NODE
//#define CBC_CHECK_BASIS
#define CBC_WEAK_STRONG
#include <cassert>
#include <cfloat>
#define CUTS
#include "OsiSolverInterface.hpp"
#include "OsiChooseVariable.hpp"
#include "OsiAuxInfo.hpp"
#include "OsiSolverBranch.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinTime.hpp"
#include "CbcModel.hpp"
#include "CbcNode.hpp"
#include "CbcStatistics.hpp"
#include "CbcStrategy.hpp"
#include "CbcBranchActual.hpp"
#include "CbcBranchDynamic.hpp"
#include "OsiRowCut.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiCuts.hpp"
#include "CbcCountRowCut.hpp"
#include "CbcFeasibilityBase.hpp"
#include "CbcMessage.hpp"
#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#include "ClpSimplexOther.hpp"
#endif
using namespace std;
#include "CglCutGenerator.hpp"
// Default Constructor 
CbcNodeInfo::CbcNodeInfo ()
  :
  numberPointingToThis_(0),
  parent_(NULL),
  owner_(NULL),
  numberCuts_(0),
  nodeNumber_(0),
  cuts_(NULL),
  numberRows_(0),
  numberBranchesLeft_(0),
  active_(7)
{
#ifdef CHECK_NODE
  printf("CbcNodeInfo %x Constructor\n",this);
#endif
}
// Constructor given parent
CbcNodeInfo::CbcNodeInfo (CbcNodeInfo * parent)
  :
  numberPointingToThis_(2),
  parent_(parent),
  owner_(NULL),
  numberCuts_(0),
  nodeNumber_(0),
  cuts_(NULL),
  numberRows_(0),
  numberBranchesLeft_(2),
  active_(7)
{
#ifdef CHECK_NODE
  printf("CbcNodeInfo %x Constructor from parent %x\n",this,parent_);
#endif
  if (parent_) {
    numberRows_ = parent_->numberRows_+parent_->numberCuts_;
    //parent_->increment();
  }
}
// Copy Constructor 
CbcNodeInfo::CbcNodeInfo (const CbcNodeInfo & rhs)
  :
  numberPointingToThis_(rhs.numberPointingToThis_),
  parent_(rhs.parent_),
  owner_(rhs.owner_),
  numberCuts_(rhs.numberCuts_),
  nodeNumber_(rhs.nodeNumber_),
  cuts_(NULL),
  numberRows_(rhs.numberRows_),
  numberBranchesLeft_(rhs.numberBranchesLeft_),
  active_(rhs.active_)
{
#ifdef CHECK_NODE
  printf("CbcNodeInfo %x Copy constructor\n",this);
#endif
  if (numberCuts_) {
    cuts_ = new CbcCountRowCut * [numberCuts_];
    int n=0;
    for (int i=0;i<numberCuts_;i++) {
      CbcCountRowCut * thisCut = rhs.cuts_[i];
      if (thisCut) {
	// I think this is correct - new one should take priority
	thisCut->setInfo(this,n);
	thisCut->increment(numberBranchesLeft_); 
	cuts_[n++] = thisCut;
      }
    }
    numberCuts_=n;
  }
}
// Constructor given parent and owner
CbcNodeInfo::CbcNodeInfo (CbcNodeInfo * parent, CbcNode * owner)
  :
  numberPointingToThis_(2),
  parent_(parent),
  owner_(owner),
  numberCuts_(0),
  nodeNumber_(0),
  cuts_(NULL),
  numberRows_(0),
  numberBranchesLeft_(2),
  active_(7)
{
#ifdef CHECK_NODE
  printf("CbcNodeInfo %x Constructor from parent %x\n",this,parent_);
#endif
  if (parent_) {
    numberRows_ = parent_->numberRows_+parent_->numberCuts_;
  }
}

/**
  Take care to detach from the owning CbcNode and decrement the reference
  count in the parent.  If this is the last nodeInfo object pointing to the
  parent, make a recursive call to delete the parent.
*/
CbcNodeInfo::~CbcNodeInfo()
{
#ifdef CHECK_NODE
  printf("CbcNodeInfo %x Destructor parent %x\n",this,parent_);
#endif

  assert(!numberPointingToThis_);
  // But there may be some left (max nodes?)
  for (int i=0;i<numberCuts_;i++) {
    if (cuts_[i]) {
#ifndef GLOBAL_CUTS_JUST_POINTERS
      delete cuts_[i];
#else
      if (cuts_[i]->globallyValidAsInteger()!=2)
	delete cuts_[i];
#endif
    }
  }
  delete [] cuts_;
  if (owner_) 
    owner_->nullNodeInfo();
  if (parent_) {
    int numberLinks = parent_->decrement();
    if (!numberLinks) delete parent_;
  }
}


//#define ALLCUTS
void
CbcNodeInfo::decrementCuts(int change)
{
  int i;
  // get rid of all remaining if negative
  int changeThis;
  if (change<0)
    changeThis = numberBranchesLeft_;
  else
    changeThis = change;
 // decrement cut counts
  for (i=0;i<numberCuts_;i++) {
    if (cuts_[i]) {
      int number = cuts_[i]->decrement(changeThis);
      if (!number) {
	//printf("info %x del cut %d %x\n",this,i,cuts_[i]);
#ifndef GLOBAL_CUTS_JUST_POINTERS
	delete cuts_[i];
#else
	if (cuts_[i]->globallyValidAsInteger()!=2)
	  delete cuts_[i];
#endif
	cuts_[i]=NULL;
      }
    }
  }
}
void
CbcNodeInfo::incrementCuts(int change)
{
  int i;
  assert (change>0);
  // increment cut counts
  for (i=0;i<numberCuts_;i++) {
    if (cuts_[i]) 
      cuts_[i]->increment(change);
  }
}
void
CbcNodeInfo::decrementParentCuts(CbcModel * model,int change)
{
  if (parent_) {
    // get rid of all remaining if negative
    int changeThis;
    if (change<0)
      changeThis = numberBranchesLeft_;
    else
      changeThis = change;
    int i;
    // Get over-estimate of space needed for basis
    CoinWarmStartBasis & dummy = model->workingBasis();
    dummy.setSize(0,numberRows_+numberCuts_);
    buildRowBasis(dummy);
    /* everything is zero (i.e. free) so we can use to see
       if latest basis */
    CbcNodeInfo * thisInfo = parent_;
    while (thisInfo) 
      thisInfo = thisInfo->buildRowBasis(dummy);
    // decrement cut counts
    thisInfo = parent_;
    int numberRows=numberRows_;
    while (thisInfo) {
      for (i=thisInfo->numberCuts_-1;i>=0;i--) {
	CoinWarmStartBasis::Status status = dummy.getArtifStatus(--numberRows);
#ifdef ALLCUTS
	status = CoinWarmStartBasis::isFree;
#endif
	if (thisInfo->cuts_[i]) {
	  int number=1;
	  if (status!=CoinWarmStartBasis::basic) {
	    // tight - drop 1 or 2
	    if (change<0)
	      number = thisInfo->cuts_[i]->decrement(changeThis);
	    else
	      number = thisInfo->cuts_[i]->decrement(change);
	  }
	  if (!number) {
#ifndef GLOBAL_CUTS_JUST_POINTERS
	    delete thisInfo->cuts_[i];
#else
	    if (thisInfo->cuts_[i]->globallyValidAsInteger()!=2)
	      delete thisInfo->cuts_[i];
#endif
	    thisInfo->cuts_[i]=NULL;
	  }
	}
      }
      thisInfo = thisInfo->parent_;
    }
  }
}
#if 0
void
CbcNodeInfo::incrementParentCuts(CbcModel * model, int change)
{
  if (parent_) {
    int i;
    // Get over-estimate of space needed for basis
    CoinWarmStartBasis & dummy = model->workingBasis();
    dummy.setSize(0,numberRows_+numberCuts_);
    /* everything is zero (i.e. free) so we can use to see
       if latest basis */
    buildRowBasis(dummy);
    CbcNodeInfo * thisInfo = parent_;
    while (thisInfo) 
      thisInfo = thisInfo->buildRowBasis(dummy);
    // increment cut counts
    thisInfo = parent_;
    int numberRows=numberRows_;
    while (thisInfo) {
      for (i=thisInfo->numberCuts_-1;i>=0;i--) {
	CoinWarmStartBasis::Status status = dummy.getArtifStatus(--numberRows);
#ifdef ALLCUTS
	status = CoinWarmStartBasis::isFree;
#endif
	if (thisInfo->cuts_[i]&&status!=CoinWarmStartBasis::basic) {
	  thisInfo->cuts_[i]->increment(change);
	}
      }
      thisInfo = thisInfo->parent_;
    }
  }
}
#endif
/*
  Append cuts to the cuts_ array in a nodeInfo. The initial reference count
  is set to numberToBranchOn, which will normally be the number of arms
  defined for the CbcBranchingObject attached to the CbcNode that owns this
  CbcNodeInfo.
*/
void
CbcNodeInfo::addCuts (OsiCuts & cuts, int numberToBranchOn,
		      int * whichGenerator)
{
  int numberCuts = cuts.sizeRowCuts();
  if (numberCuts) {
    int i;
    if (!numberCuts_) {
      cuts_ = new CbcCountRowCut * [numberCuts];
    } else {
      CbcCountRowCut ** temp = new CbcCountRowCut * [numberCuts+numberCuts_];
      memcpy(temp,cuts_,numberCuts_*sizeof(CbcCountRowCut *));
      delete [] cuts_;
      cuts_ = temp;
    }
    for (i=0;i<numberCuts;i++) {
      CbcCountRowCut * thisCut = new CbcCountRowCut(*cuts.rowCutPtr(i),
						    this,numberCuts_);
      thisCut->increment(numberToBranchOn); 
      cuts_[numberCuts_++] = thisCut;
#ifdef CBC_DEBUG
#if CBC_DEBUG>1
      int n=thisCut->row().getNumElements();
      printf("Cut %d has %d entries, rhs %g %g =>",i,n,thisCut->lb(),
	     thisCut->ub());
      int j;
      const int * index = thisCut->row().getIndices();
      const double * element = thisCut->row().getElements();
      for (j=0;j<n;j++) {
	printf(" (%d,%g)",index[j],element[j]);
	assert(fabs(element[j])>1.00e-12);
      }
      printf("\n");
#else
      int n=thisCut->row().getNumElements();
      int j;
      const double * element = thisCut->row().getElements();
      for (j=0;j<n;j++) {
	assert(fabs(element[j])>1.00e-12);
      }
#endif
#endif
    }
  }
}

void
CbcNodeInfo::addCuts(int numberCuts, CbcCountRowCut ** cut, 
		     int numberToBranchOn)
{
  if (numberCuts) {
    int i;
    if (!numberCuts_) {
      cuts_ = new CbcCountRowCut * [numberCuts];
    } else {
      CbcCountRowCut ** temp = new CbcCountRowCut * [numberCuts+numberCuts_];
      memcpy(temp,cuts_,numberCuts_*sizeof(CbcCountRowCut *));
      delete [] cuts_;
      cuts_ = temp;
    }
    for (i=0;i<numberCuts;i++) {
      CbcCountRowCut * thisCut = cut[i];
      thisCut->setInfo(this,numberCuts_);
      //printf("info %x cut %d %x\n",this,i,thisCut);
      thisCut->increment(numberToBranchOn); 
      cuts_[numberCuts_++] = thisCut;
#ifdef CBC_DEBUG
      int n=thisCut->row().getNumElements();
#if CBC_DEBUG>1
      printf("Cut %d has %d entries, rhs %g %g =>",i,n,thisCut->lb(),
	     thisCut->ub());
#endif
      int j;
#if CBC_DEBUG>1
      const int * index = thisCut->row().getIndices();
#endif
      const double * element = thisCut->row().getElements();
      for (j=0;j<n;j++) {
#if CBC_DEBUG>1
	printf(" (%d,%g)",index[j],element[j]);
#endif
	assert(fabs(element[j])>1.00e-12);
      }
      printf("\n");
#endif
    }
  }
}

// delete cuts
void
CbcNodeInfo::deleteCuts(int numberToDelete, CbcCountRowCut ** cuts)
{
  int i;
  int j;
  int last=-1;
  for (i=0;i<numberToDelete;i++) {
    CbcCountRowCut * next = cuts[i];
    for (j=last+1;j<numberCuts_;j++) {
      if (next==cuts_[j])
	break;
    }
    if (j==numberCuts_) {
      // start from beginning
      for (j=0;j<last;j++) {
	if (next==cuts_[j])
	  break;
      }
      assert(j<last);
    }
    last=j;
    int number = cuts_[j]->decrement();
    if (!number) {
#ifndef GLOBAL_CUTS_JUST_POINTERS
      delete cuts_[j];
#else
      if (cuts_[j]->globallyValidAsInteger()!=2)
	delete cuts_[j];
#endif
    }
    cuts_[j]=NULL;
  }
  j=0;
  for (i=0;i<numberCuts_;i++) {
    if (cuts_[i])
      cuts_[j++]=cuts_[i];
  }
  numberCuts_ = j;
}

// delete cuts
void
CbcNodeInfo::deleteCuts(int numberToDelete, int * which)
{
  int i;
  for (i=0;i<numberToDelete;i++) {
    int iCut=which[i];
    int number = cuts_[iCut]->decrement();
    if (!number) {
#ifndef GLOBAL_CUTS_JUST_POINTERS
      delete cuts_[iCut];
#else
      if (cuts_[iCut]->globallyValidAsInteger()!=2)
	delete cuts_[iCut];
#endif
    }
    cuts_[iCut]=NULL;
  }
  int j=0;
  for (i=0;i<numberCuts_;i++) {
    if (cuts_[i])
      cuts_[j++]=cuts_[i];
  }
  numberCuts_ = j;
}

// Really delete a cut
void 
CbcNodeInfo::deleteCut(int whichOne)
{
  assert(whichOne<numberCuts_);
  cuts_[whichOne]=NULL;
}
/* Deactivate node information.
   1 - bounds
   2 - cuts
   4 - basis!
*/
void 
CbcNodeInfo::deactivate(int mode)
{
  active_ &= (~mode);
}

CbcFullNodeInfo::CbcFullNodeInfo() :
  CbcNodeInfo(),
  basis_(),
  numberIntegers_(0),
  lower_(NULL),
  upper_(NULL)
{
}
CbcFullNodeInfo::CbcFullNodeInfo(CbcModel * model,
				 int numberRowsAtContinuous) :
  CbcNodeInfo()
{
  OsiSolverInterface * solver = model->solver();
  numberRows_ = numberRowsAtContinuous;
  numberIntegers_ = model->numberIntegers();
  int numberColumns = model->getNumCols();
  lower_ = new double [numberColumns];
  upper_ = new double [numberColumns];
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  int i;

  for (i=0;i<numberColumns;i++) {
    lower_[i]=lower[i];
    upper_[i]=upper[i];
  }

  basis_ =  dynamic_cast<CoinWarmStartBasis*>(solver->getWarmStart());
}

CbcFullNodeInfo::CbcFullNodeInfo(const CbcFullNodeInfo & rhs) :
  CbcNodeInfo(rhs)
{  
  basis_= dynamic_cast<CoinWarmStartBasis *>(rhs.basis_->clone()) ;
  numberIntegers_=rhs.numberIntegers_;
  lower_=NULL;
  upper_=NULL;
  if (rhs.lower_!=NULL) {
    int numberColumns = basis_->getNumStructural();
    lower_ = new double [numberColumns];
    upper_ = new double [numberColumns];
    assert (upper_!=NULL);
    memcpy(lower_,rhs.lower_,numberColumns*sizeof(double));
    memcpy(upper_,rhs.upper_,numberColumns*sizeof(double));
  }
}

CbcNodeInfo * 
CbcFullNodeInfo::clone() const
{ 
  return (new CbcFullNodeInfo(*this));
}

CbcFullNodeInfo::~CbcFullNodeInfo ()
{
  delete basis_ ;
  delete [] lower_;
  delete [] upper_;
}

/*
   The basis supplied as a parameter is deleted and replaced with a new basis
   appropriate for the node, and lower and upper bounds on variables are
   reset according to the stored bounds arrays. Any cuts associated with this
   node are added to the list in addCuts, but not actually added to the
   constraint system in the model.

   Why pass in a basis at all? The short answer is ``We need the parameter to
   pass out a basis, so might as well use it to pass in the size.''
   
   A longer answer is that in practice we take a memory allocation hit up in
   addCuts1 (the only place applyToModel is called) when we setSize() the
   basis that's passed in. It's immediately tossed here in favour of a clone
   of the basis attached to this nodeInfo. This can probably be fixed, given
   a bit of thought.
*/

void CbcFullNodeInfo::applyToModel (CbcModel *model,
				    CoinWarmStartBasis *&basis,
				    CbcCountRowCut **addCuts,
				    int &currentNumberCuts) const 

{ OsiSolverInterface *solver = model->solver() ;

  // branch - do bounds
  assert (active_==7||active_==15);
  int i;
  solver->setColLower(lower_);
  solver->setColUpper(upper_);
  int numberColumns = model->getNumCols();
  // move basis - but make sure size stays
  // for bon-min - should not be needed int numberRows = model->getNumRows();
  int numberRows=basis->getNumArtificial();
  delete basis ;
  if (basis_) {
    basis = dynamic_cast<CoinWarmStartBasis *>(basis_->clone()) ;
    basis->resize(numberRows,numberColumns);
  } else {
    // We have a solver without a basis
    basis=NULL;
  }
  for (i=0;i<numberCuts_;i++) 
    addCuts[currentNumberCuts+i]= cuts_[i];
  currentNumberCuts += numberCuts_;
  assert(!parent_);
  return ;
}
// Just apply bounds to one variable (1=>infeasible)
int 
CbcFullNodeInfo::applyBounds(int iColumn, double & lower, double & upper,int force) 
{
  if ((force&&1)==0) {
    if (lower>lower_[iColumn])
      printf("%d odd lower going from %g to %g\n",iColumn,lower,lower_[iColumn]);
    lower = lower_[iColumn];
  } else {
    lower_[iColumn]=lower;
  }
  if ((force&&2)==0) {
    if (upper<upper_[iColumn])
      printf("%d odd upper going from %g to %g\n",iColumn,upper,upper_[iColumn]);
    upper = upper_[iColumn];
  } else {
    upper_[iColumn]=upper;
  }
  return (upper_[iColumn]>=lower_[iColumn]) ? 0 : 1;
}

/* Builds up row basis backwards (until original model).
   Returns NULL or previous one to apply .
   Depends on Free being 0 and impossible for cuts
*/
CbcNodeInfo * 
CbcFullNodeInfo::buildRowBasis(CoinWarmStartBasis & basis ) const 
{
  const unsigned int * saved = 
    (const unsigned int *) basis_->getArtificialStatus();
  unsigned int * now = 
    (unsigned int *) basis.getArtificialStatus();
  int number=basis_->getNumArtificial()>>4;;
  int i;
  for (i=0;i<number;i++) { 
    if (!now[i])
      now[i] = saved[i];
  }
  return NULL;
}


// Default constructor
CbcPartialNodeInfo::CbcPartialNodeInfo()

  : CbcNodeInfo(),
    basisDiff_(NULL),
    variables_(NULL),
    newBounds_(NULL),
    numberChangedBounds_(0)

{ /* this space intentionally left blank */ }

// Constructor from current state 
CbcPartialNodeInfo::CbcPartialNodeInfo (CbcNodeInfo *parent, CbcNode *owner,
					int numberChangedBounds,
					const int *variables,
					const double *boundChanges,
					const CoinWarmStartDiff *basisDiff)
 : CbcNodeInfo(parent,owner)
{
  basisDiff_ = basisDiff->clone() ;

  numberChangedBounds_ = numberChangedBounds;
  int size = numberChangedBounds_*(sizeof(double)+sizeof(int));
  char * temp = new char [size];
  newBounds_ = (double *) temp;
  variables_ = (int *) (newBounds_+numberChangedBounds_);

  int i ;
  for (i=0;i<numberChangedBounds_;i++) {
    variables_[i]=variables[i];
    newBounds_[i]=boundChanges[i];
  }
}

CbcPartialNodeInfo::CbcPartialNodeInfo (const CbcPartialNodeInfo & rhs)

  : CbcNodeInfo(rhs.parent_)

{ basisDiff_ = rhs.basisDiff_->clone() ;

  numberChangedBounds_ = rhs.numberChangedBounds_;
  int size = numberChangedBounds_*(sizeof(double)+sizeof(int));
  char * temp = new char [size];
  newBounds_ = (double *) temp;
  variables_ = (int *) (newBounds_+numberChangedBounds_);

  int i ;
  for (i=0;i<numberChangedBounds_;i++) {
    variables_[i]=rhs.variables_[i];
    newBounds_[i]=rhs.newBounds_[i];
  }
}

CbcNodeInfo * 
CbcPartialNodeInfo::clone() const
{ 
  return (new CbcPartialNodeInfo(*this));
}


CbcPartialNodeInfo::~CbcPartialNodeInfo ()
{
  delete basisDiff_ ;
  delete [] newBounds_;
}


/**
   The basis supplied as a parameter is incrementally modified, and lower and
   upper bounds on variables in the model are incrementally modified. Any
   cuts associated with this node are added to the list in addCuts.
*/

void CbcPartialNodeInfo::applyToModel (CbcModel *model,
				       CoinWarmStartBasis *&basis,
				       CbcCountRowCut **addCuts,
				       int &currentNumberCuts) const 

{ OsiSolverInterface *solver = model->solver();
  if ((active_&4)!=0) {
    basis->applyDiff(basisDiff_) ;
  }

  // branch - do bounds
  int i;
  if ((active_&1)!=0) {
    for (i=0;i<numberChangedBounds_;i++) {
      int variable = variables_[i];
      int k = variable&0x3fffffff;
      if ((variable&0x80000000)==0) {
	// lower bound changing
	//#define CBC_PRINT2
#ifdef CBC_PRINT2
	if(solver->getColLower()[k]!=newBounds_[i])
	  printf("lower change for column %d - from %g to %g\n",k,solver->getColLower()[k],newBounds_[i]);
#endif
#ifndef NDEBUG
	if ((variable&0x40000000)==0&&false) {
	  double oldValue = solver->getColLower()[k];
	  assert (newBounds_[i]>oldValue-1.0e-8);
	  if (newBounds_[i]<oldValue+1.0e-8)
	    printf("bad null lower change for column %d - bound %g\n",k,oldValue);
	}
#endif
	solver->setColLower(k,newBounds_[i]);
      } else {
	// upper bound changing
#ifdef CBC_PRINT2
	if(solver->getColUpper()[k]!=newBounds_[i])
	  printf("upper change for column %d - from %g to %g\n",k,solver->getColUpper()[k],newBounds_[i]);
#endif
#ifndef NDEBUG
	if ((variable&0x40000000)==0&&false) {
	  double oldValue = solver->getColUpper()[k];
	  assert (newBounds_[i]<oldValue+1.0e-8);
	  if (newBounds_[i]>oldValue-1.0e-8)
	    printf("bad null upper change for column %d - bound %g\n",k,oldValue);
	}
#endif
	solver->setColUpper(k,newBounds_[i]);
      }
    }
  }
  if ((active_&2)!=0) {
    for (i=0;i<numberCuts_;i++) {
      addCuts[currentNumberCuts+i]= cuts_[i];
      if (cuts_[i]&&model->messageHandler()->logLevel()>4) {
	cuts_[i]->print();
      }
    }
    
    currentNumberCuts += numberCuts_;
  }
  return ;
}
// Just apply bounds to one variable (1=>infeasible)
int
CbcPartialNodeInfo::applyBounds(int iColumn, double & lower, double & upper,int force) 
{
  // branch - do bounds
  int i;
  int found=0;
  double newLower = -COIN_DBL_MAX;
  double newUpper = COIN_DBL_MAX;
  for (i=0;i<numberChangedBounds_;i++) {
    int variable = variables_[i];
    int k = variable&0x3fffffff;
    if (k==iColumn) {
      if ((variable&0x80000000)==0) {
	// lower bound changing
	found |= 1;
	newLower = CoinMax(newLower,newBounds_[i]);
	if ((force&1)==0) {
	  if (lower>newBounds_[i])
	    printf("%d odd lower going from %g to %g\n",iColumn,lower,newBounds_[i]);
	  lower = newBounds_[i];
	} else {
	  newBounds_[i]=lower;
	  variables_[i] |= 0x40000000; // say can go odd way
	}
      } else {
	// upper bound changing
	found |= 2;
	newUpper = CoinMin(newUpper,newBounds_[i]);
	if ((force&2)==0) {
	  if (upper<newBounds_[i])
	    printf("%d odd upper going from %g to %g\n",iColumn,upper,newBounds_[i]);
	  upper = newBounds_[i];
	} else {
	  newBounds_[i]=upper;
	  variables_[i] |= 0x40000000; // say can go odd way
	}
      }
    }
  }
  newLower = CoinMax(newLower,lower);
  newUpper = CoinMin(newUpper,upper);
  int nAdd=0;
  if ((force&2)!=0&&(found&2)==0) {
    // need to add new upper
    nAdd++;
  }
  if ((force&1)!=0&&(found&1)==0) {
    // need to add new lower
    nAdd++;
  }
  if (nAdd) { 
    int size = (numberChangedBounds_+nAdd)*(sizeof(double)+sizeof(int));
    char * temp = new char [size];
    double * newBounds = (double *) temp;
    int * variables = (int *) (newBounds+numberChangedBounds_+nAdd);

    int i ;
    for (i=0;i<numberChangedBounds_;i++) {
      variables[i]=variables_[i];
      newBounds[i]=newBounds_[i];
    }
    delete [] newBounds_;
    newBounds_ = newBounds;
    variables_ = variables;
    if ((force&2)!=0&&(found&2)==0) {
      // need to add new upper
      int variable = iColumn | 0x80000000;
      variables_[numberChangedBounds_]=variable;
      newBounds_[numberChangedBounds_++]=newUpper;
    }
    if ((force&1)!=0&&(found&1)==0) {
      // need to add new lower
      int variable = iColumn;
      variables_[numberChangedBounds_]=variable;
      newBounds_[numberChangedBounds_++]=newLower;
    }
  }
  
  return (newUpper>=newLower) ? 0 : 1;
}

/* Builds up row basis backwards (until original model).
   Returns NULL or previous one to apply .
   Depends on Free being 0 and impossible for cuts
*/

CbcNodeInfo * 
CbcPartialNodeInfo::buildRowBasis(CoinWarmStartBasis & basis ) const 

{ basis.applyDiff(basisDiff_) ;

  return parent_ ; }

CbcNode::CbcNode() :
  nodeInfo_(NULL),
  objectiveValue_(1.0e100),
  guessedObjectiveValue_(1.0e100),
  sumInfeasibilities_(0.0),
  branch_(NULL),
  depth_(-1),
  numberUnsatisfied_(0),
  nodeNumber_(-1),
  state_(0)
{
#ifdef CHECK_NODE
  printf("CbcNode %x Constructor\n",this);
#endif
}
// Print
void 
CbcNode::print() const
{
  printf("number %d obj %g depth %d sumun %g nunsat %d state %d\n",
	 nodeNumber_,objectiveValue_,depth_,sumInfeasibilities_,numberUnsatisfied_,state_);
}
CbcNode::CbcNode(CbcModel * model,
		 CbcNode * lastNode) :
  nodeInfo_(NULL),
  objectiveValue_(1.0e100),
  guessedObjectiveValue_(1.0e100),
  sumInfeasibilities_(0.0),
  branch_(NULL),
  depth_(-1),
  numberUnsatisfied_(0),
  nodeNumber_(-1),
  state_(0)
{
#ifdef CHECK_NODE
  printf("CbcNode %x Constructor from model\n",this);
#endif
  model->setObjectiveValue(this,lastNode);

  if (lastNode) {
    if (lastNode->nodeInfo_) {
       lastNode->nodeInfo_->increment();
    }
  }
  nodeNumber_= model->getNodeCount();
}

#define CBC_NEW_CREATEINFO
#ifdef CBC_NEW_CREATEINFO

/*
  New createInfo, with basis manipulation hidden inside mergeBasis. Allows
  solvers to override and carry over all information from one basis to
  another.
*/

void
CbcNode::createInfo (CbcModel *model,
		     CbcNode *lastNode,
		     const CoinWarmStartBasis *lastws,
		     const double *lastLower, const double *lastUpper,
		     int numberOldActiveCuts, int numberNewCuts)

{ OsiSolverInterface *solver = model->solver();
  CbcStrategy *strategy = model->strategy();
/*
  The root --- no parent. Create full basis and bounds information.
*/
  if (!lastNode)
  { 
    if (!strategy)
      nodeInfo_=new CbcFullNodeInfo(model,solver->getNumRows());
    else
      nodeInfo_ = strategy->fullNodeInfo(model,solver->getNumRows());
  } else {
/*
  Not the root. Create an edit from the parent's basis & bound information.
  This is not quite as straightforward as it seems. We need to reintroduce
  cuts we may have dropped out of the basis, in the correct position, because
  this whole process is strictly positional. Start by grabbing the current
  basis.
*/
    bool mustDeleteBasis;
    const CoinWarmStartBasis *ws =
      dynamic_cast<const CoinWarmStartBasis*>(solver->getPointerToWarmStart(mustDeleteBasis));
    assert(ws!=NULL); // make sure not volume
    //int numberArtificials = lastws->getNumArtificial();
    int numberColumns = solver->getNumCols();
    int numberRowsAtContinuous = model->numberRowsAtContinuous();
    int currentNumberCuts = model->currentNumberCuts();
#   ifdef CBC_CHECK_BASIS
    std::cout
      << "Before expansion: orig " << numberRowsAtContinuous
      << ", old " << numberOldActiveCuts
      << ", new " << numberNewCuts
      << ", current " << currentNumberCuts << "." << std::endl ;
    ws->print();
#   endif
/*
  Clone the basis and resize it to hold the structural constraints, plus
  all the cuts: old cuts, both active and inactive (currentNumberCuts),
  and new cuts (numberNewCuts). This will become the expanded basis.
*/
    CoinWarmStartBasis *expanded = 
      dynamic_cast<CoinWarmStartBasis *>(ws->clone()) ;
    int iCompact = numberRowsAtContinuous+numberOldActiveCuts+numberNewCuts ;
    // int nPartial = numberRowsAtContinuous+currentNumberCuts;
    int iFull = numberRowsAtContinuous+currentNumberCuts+numberNewCuts;
    // int maxBasisLength = ((iFull+15)>>4)+((numberColumns+15)>>4);
    // printf("l %d full %d\n",maxBasisLength,iFull);
    expanded->resize(iFull,numberColumns);
#   ifdef CBC_CHECK_BASIS
    std::cout
      << "\tFull basis " << iFull << " rows, "
      << numberColumns << " columns; compact "
      << iCompact << " rows." << std::endl ;
#   endif
/*
  Now flesh out the expanded basis. The clone already has the
  correct status information for the variables and for the structural
  (numberRowsAtContinuous) constraints. Any indices beyond nPartial must be
  cuts created while processing this node --- they can be copied en bloc
  into the correct position in the expanded basis. The space reserved for
  xferRows is a gross overestimate.
*/
    CoinWarmStartBasis::XferVec xferRows ;
    xferRows.reserve(iFull-numberRowsAtContinuous+1) ;
    if (numberNewCuts) {
      xferRows.push_back(
          CoinWarmStartBasis::XferEntry(iCompact-numberNewCuts,
					iFull-numberNewCuts,numberNewCuts)) ;
    }
/*
  From nPartial down, record the entries we want to copy from the current
  basis (the entries for the active cuts; non-zero in the list returned
  by addedCuts). Fill the expanded basis with entries showing a status of
  basic for the deactivated (loose) cuts.
*/
    CbcCountRowCut **cut = model->addedCuts();
    iFull -= (numberNewCuts+1) ;
    iCompact -= (numberNewCuts+1) ;
    int runLen = 0 ;
    CoinWarmStartBasis::XferEntry entry(-1,-1,-1) ;
    while (iFull >= numberRowsAtContinuous) { 
      for ( ; iFull >= numberRowsAtContinuous &&
	      cut[iFull-numberRowsAtContinuous] ; iFull--)
        runLen++ ;
      if (runLen) {
        iCompact -= runLen ;
        entry.first = iCompact+1 ;
        entry.second = iFull+1 ;
	entry.third = runLen ;
	runLen = 0 ;
        xferRows.push_back(entry) ;
      }
      for ( ; iFull >= numberRowsAtContinuous &&
      	      !cut[iFull-numberRowsAtContinuous] ; iFull--)
	expanded->setArtifStatus(iFull,CoinWarmStartBasis::basic);
    }
/*
  Finally, call mergeBasis to copy over entries from the current basis to
  the expanded basis. Since we cloned the expanded basis from the active basis
  and haven't changed the number of variables, only row status entries need
  to be copied.
*/
    expanded->mergeBasis(ws,&xferRows,0) ;

#ifdef CBC_CHECK_BASIS
    std::cout << "Expanded basis:" << std::endl ;
    expanded->print() ;
    std::cout << "Diffing against:" << std::endl ;
    lastws->print() ;
#endif    

/*
  Now that we have two bases in proper positional correspondence, creating
  the actual diff is dead easy.

  Note that we're going to compare the expanded basis here to the stripped
  basis (lastws) produced by addCuts. It doesn't affect the correctness (the
  diff process has no knowledge of the meaning of an entry) but it does
  mean that we'll always generate a whack of diff entries because the expanded
  basis is considerably larger than the stripped basis.
*/
    CoinWarmStartDiff *basisDiff = expanded->generateDiff(lastws) ;
/*
  Diff the bound vectors. It's assumed the number of structural variables
  is not changing. For branching objects that change bounds on integer
  variables, we should see at least one bound change as a consequence
  of applying the branch that generated this subproblem from its parent.
  This need not hold for other types of branching objects (hyperplane
  branches, for example).
*/
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();

    double *boundChanges = new double [2*numberColumns] ;
    int *variables = new int [2*numberColumns] ;
    int numberChangedBounds=0;
    
    int i;
    for (i=0;i<numberColumns;i++) {
      if (lower[i]!=lastLower[i]) {
	variables[numberChangedBounds]=i;
	boundChanges[numberChangedBounds++]=lower[i];
      }
      if (upper[i]!=lastUpper[i]) {
	variables[numberChangedBounds]=i|0x80000000;
	boundChanges[numberChangedBounds++]=upper[i];
      }
#ifdef CBC_DEBUG
      if (lower[i] != lastLower[i]) {
        std::cout
	  << "lower on " << i << " changed from "
	  << lastLower[i] << " to " << lower[i] << std::endl ;
      }
      if (upper[i] != lastUpper[i]) {
        std::cout
	  << "upper on " << i << " changed from "
	  << lastUpper[i] << " to " << upper[i] << std::endl ;
      }
#endif
    }
#ifdef CBC_DEBUG
    std::cout << numberChangedBounds << " changed bounds." << std::endl ;
#endif
    //if (lastNode->branchingObject()->boundBranch())
    //assert (numberChangedBounds);
/*
  Hand the lot over to the CbcPartialNodeInfo constructor, then clean up and
  return.
*/
    if (!strategy)
      nodeInfo_ =
        new CbcPartialNodeInfo(lastNode->nodeInfo_,this,numberChangedBounds,
                               variables,boundChanges,basisDiff) ;
    else
      nodeInfo_ =
        strategy->partialNodeInfo(model,lastNode->nodeInfo_,this,
				  numberChangedBounds,variables,boundChanges,
				  basisDiff) ;
    delete basisDiff ;
    delete [] boundChanges;
    delete [] variables;
    delete expanded ;
    if  (mustDeleteBasis)
      delete ws;
  }
  // Set node number
  nodeInfo_->setNodeNumber(model->getNodeCount2());
  state_ |= 2; // say active
}

#else	// CBC_NEW_CREATEINFO

/*
  Original createInfo, with bare manipulation of basis vectors. Fails if solver
  maintains additional information in basis.
*/

void
CbcNode::createInfo (CbcModel *model,
		     CbcNode *lastNode,
		     const CoinWarmStartBasis *lastws,
		     const double *lastLower, const double *lastUpper,
		     int numberOldActiveCuts,int numberNewCuts)
{ OsiSolverInterface * solver = model->solver();
 CbcStrategy * strategy = model->strategy();
/*
  The root --- no parent. Create full basis and bounds information.
*/
  if (!lastNode)
  { 
    if (!strategy)
      nodeInfo_=new CbcFullNodeInfo(model,solver->getNumRows());
    else
      nodeInfo_ = strategy->fullNodeInfo(model,solver->getNumRows());
  }
/*
  Not the root. Create an edit from the parent's basis & bound information.
  This is not quite as straightforward as it seems. We need to reintroduce
  cuts we may have dropped out of the basis, in the correct position, because
  this whole process is strictly positional. Start by grabbing the current
  basis.
*/
  else
  { 
    bool mustDeleteBasis;
    const CoinWarmStartBasis* ws =
      dynamic_cast<const CoinWarmStartBasis*>(solver->getPointerToWarmStart(mustDeleteBasis));
    assert(ws!=NULL); // make sure not volume
    //int numberArtificials = lastws->getNumArtificial();
    int numberColumns = solver->getNumCols();
    
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();

    int i;
/*
  Create a clone and resize it to hold all the structural constraints, plus
  all the cuts: old cuts, both active and inactive (currentNumberCuts), and
  new cuts (numberNewCuts).

  TODO: You'd think that the set of constraints (logicals) in the expanded
	basis should match the set represented in lastws. At least, that's
	what I thought. But at the point I first looked hard at this bit of
	code, it turned out that lastws was the stripped basis produced at
	the end of addCuts(), rather than the raw basis handed back by
	addCuts1(). The expanded basis here is equivalent to the raw basis of
	addCuts1(). I said ``whoa, that's not good, I must have introduced a
	bug'' and went back to John's code to see where I'd gone wrong.
	And discovered the same `error' in his code.

	After a bit of thought, my conclusion is that correctness is not
	affected by whether lastws is the stripped or raw basis. The diffs
	have no semantics --- just a set of changes that need to be made
	to convert lastws into expanded. I think the only effect is that we
	store a lot more diffs (everything in expanded that's not covered by
	the stripped basis). But I need to give this more thought. There
	may well be some subtle error cases.

	In the mean time, I've twiddled addCuts() to set lastws to the raw
	basis. Makes me (Lou) less nervous to compare apples to apples.
*/
    CoinWarmStartBasis *expanded = 
      dynamic_cast<CoinWarmStartBasis *>(ws->clone()) ;
    int numberRowsAtContinuous = model->numberRowsAtContinuous();
    int iFull = numberRowsAtContinuous+model->currentNumberCuts()+
      numberNewCuts;
    //int numberArtificialsNow = iFull;
    //int maxBasisLength = ((iFull+15)>>4)+((numberColumns+15)>>4);
    //printf("l %d full %d\n",maxBasisLength,iFull);
    if (expanded) 
      expanded->resize(iFull,numberColumns);
#ifdef CBC_CHECK_BASIS
    printf("Before expansion: orig %d, old %d, new %d, current %d\n",
	   numberRowsAtContinuous,numberOldActiveCuts,numberNewCuts,
	   model->currentNumberCuts()) ;
    ws->print();
#endif
/*
  Now fill in the expanded basis. Any indices beyond nPartial must
  be cuts created while processing this node --- they can be copied directly
  into the expanded basis. From nPartial down, pull the status of active cuts
  from ws, interleaving with a B entry for the deactivated (loose) cuts.
*/
    int numberDropped = model->currentNumberCuts()-numberOldActiveCuts;
    int iCompact=iFull-numberDropped;
    CbcCountRowCut ** cut = model->addedCuts();
    int nPartial = model->currentNumberCuts()+numberRowsAtContinuous;
    iFull--;
    for (;iFull>=nPartial;iFull--) {
      CoinWarmStartBasis::Status status = ws->getArtifStatus(--iCompact);
      //assert (status != CoinWarmStartBasis::basic); // may be permanent cut
      expanded->setArtifStatus(iFull,status);
    }
    for (;iFull>=numberRowsAtContinuous;iFull--) {
      if (cut[iFull-numberRowsAtContinuous]) {
	CoinWarmStartBasis::Status status = ws->getArtifStatus(--iCompact);
	// If no cut generator being used then we may have basic variables
	//if (model->getMaximumCutPasses()&&
	//  status == CoinWarmStartBasis::basic)
	//printf("cut basic\n");
	expanded->setArtifStatus(iFull,status);
      } else {
	expanded->setArtifStatus(iFull,CoinWarmStartBasis::basic);
      }
    }
#ifdef CBC_CHECK_BASIS
    printf("Expanded basis\n");
    expanded->print() ;
    printf("Diffing against\n") ;
    lastws->print() ;
#endif    
/*
  Now that we have two bases in proper positional correspondence, creating
  the actual diff is dead easy.
*/

    CoinWarmStartDiff *basisDiff = expanded->generateDiff(lastws) ;
/*
  Diff the bound vectors. It's assumed the number of structural variables is
  not changing. Assuming that branching objects all involve integer variables,
  we should see at least one bound change as a consequence of processing this
  subproblem. Different types of branching objects could break this assertion.
  Not true at all - we have not applied current branch - JJF.
*/
    double *boundChanges = new double [2*numberColumns] ;
    int *variables = new int [2*numberColumns] ;
    int numberChangedBounds=0;
    for (i=0;i<numberColumns;i++) {
      if (lower[i]!=lastLower[i]) {
	variables[numberChangedBounds]=i;
	boundChanges[numberChangedBounds++]=lower[i];
      }
      if (upper[i]!=lastUpper[i]) {
	variables[numberChangedBounds]=i|0x80000000;
	boundChanges[numberChangedBounds++]=upper[i];
      }
#ifdef CBC_DEBUG
      if (lower[i]!=lastLower[i])
	printf("lower on %d changed from %g to %g\n",
	       i,lastLower[i],lower[i]);
      if (upper[i]!=lastUpper[i])
	printf("upper on %d changed from %g to %g\n",
	       i,lastUpper[i],upper[i]);
#endif
    }
#ifdef CBC_DEBUG
    printf("%d changed bounds\n",numberChangedBounds) ;
#endif
    //if (lastNode->branchingObject()->boundBranch())
    //assert (numberChangedBounds);
/*
  Hand the lot over to the CbcPartialNodeInfo constructor, then clean up and
  return.
*/
    if (!strategy)
      nodeInfo_ =
        new CbcPartialNodeInfo(lastNode->nodeInfo_,this,numberChangedBounds,
                               variables,boundChanges,basisDiff) ;
    else
      nodeInfo_ = strategy->partialNodeInfo(model, lastNode->nodeInfo_,this,numberChangedBounds,
                               variables,boundChanges,basisDiff) ;
    delete basisDiff ;
    delete [] boundChanges;
    delete [] variables;
    delete expanded ;
    if  (mustDeleteBasis)
      delete ws;
  }
  // Set node number
  nodeInfo_->setNodeNumber(model->getNodeCount2());
  state_ |= 2; // say active
}

#endif	// CBC_NEW_CREATEINFO

/*
  The routine scans through the object list of the model looking for objects
  that indicate infeasibility. It tests each object using strong branching
  and selects the one with the least objective degradation.  A corresponding
  branching object is left attached to lastNode.

  If strong branching is disabled, a candidate object is chosen essentially
  at random (whatever object ends up in pos'n 0 of the candidate array).

  If a branching candidate is found to be monotone, bounds are set to fix the
  variable and the routine immediately returns (the caller is expected to
  reoptimize).

  If a branching candidate is found to result in infeasibility in both
  directions, the routine immediately returns an indication of infeasibility.

  Returns:  0	both branch directions are feasible
	   -1	branching variable is monotone
	   -2	infeasible

  Original comments:
    Here could go cuts etc etc
    For now just fix on objective from strong branching.
*/

int CbcNode::chooseBranch (CbcModel *model, CbcNode *lastNode,int numberPassesLeft)

{ if (lastNode)
    depth_ = lastNode->depth_+1;
  else
    depth_ = 0;
  delete branch_;
  branch_=NULL;
  OsiSolverInterface * solver = model->solver();
  double saveObjectiveValue = solver->getObjValue();
  double objectiveValue = CoinMax(solver->getObjSense()*saveObjectiveValue,objectiveValue_);
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  // See what user thinks
  int anyAction=model->problemFeasibility()->feasible(model,0);
  if (anyAction) {
    // will return -2 if infeasible , 0 if treat as integer
    return anyAction-1;
  }
  double integerTolerance = 
    model->getDblParam(CbcModel::CbcIntegerTolerance);
  // point to useful information
  OsiBranchingInformation usefulInfo = model->usefulInformation();
  // and modify
  usefulInfo.depth_=depth_;
  int i;
  bool beforeSolution = model->getSolutionCount()==0;
  int numberStrong=model->numberStrong();
  // switch off strong if hotstart
  if (model->hotstartSolution())
    numberStrong=0;
  int numberStrongDone=0;
  int numberUnfinished=0;
  int numberStrongInfeasible=0;
  int numberStrongIterations=0;
  int saveNumberStrong=numberStrong;
  int numberObjects = model->numberObjects();
  bool checkFeasibility = numberObjects>model->numberIntegers();
  int maximumStrong = CoinMax(CoinMin(numberStrong,numberObjects),1);
  int numberColumns = model->getNumCols();
  double * saveUpper = new double[numberColumns];
  double * saveLower = new double[numberColumns];

  // Save solution in case heuristics need good solution later
  
  double * saveSolution = new double[numberColumns];
  memcpy(saveSolution,solver->getColSolution(),numberColumns*sizeof(double));
  model->reserveCurrentSolution(saveSolution);
  /*
    Get a branching decision object. Use the default decision criteria unless
    the user has loaded a decision method into the model.
  */
  CbcBranchDecision *decision = model->branchingMethod();
  CbcDynamicPseudoCostBranchingObject * dynamicBranchingObject =
    dynamic_cast<CbcDynamicPseudoCostBranchingObject *>(decision);
  if (!decision||dynamicBranchingObject)
    decision = new CbcBranchDefaultDecision();
  decision->initialize(model);
  CbcStrongInfo * choice = new CbcStrongInfo[maximumStrong];
  for (i=0;i<numberColumns;i++) {
    saveLower[i] = lower[i];
    saveUpper[i] = upper[i];
  }
  // May go round twice if strong branching fixes all local candidates
  bool finished=false;
  double estimatedDegradation=0.0; 
  while(!finished) {
    finished=true;
    // Some objects may compute an estimate of best solution from here
    estimatedDegradation=0.0; 
    //int numberIntegerInfeasibilities=0; // without odd ones
    numberStrongDone=0;
    numberUnfinished=0;
    numberStrongInfeasible=0;
    numberStrongIterations=0;
    
    // We may go round this loop twice (only if we think we have solution)
    for (int iPass=0;iPass<2;iPass++) {
      
      // compute current state
      //int numberObjectInfeasibilities; // just odd ones
      //model->feasibleSolution(
      //                      numberIntegerInfeasibilities,
      //                      numberObjectInfeasibilities);
      const double * hotstartSolution = model->hotstartSolution();
      const int * hotstartPriorities = model->hotstartPriorities();
      
      // Some objects may compute an estimate of best solution from here
      estimatedDegradation=0.0; 
      numberUnsatisfied_ = 0;
      // initialize sum of "infeasibilities"
      sumInfeasibilities_ = 0.0;
      int bestPriority=COIN_INT_MAX;
      /*
        Scan for branching objects that indicate infeasibility. Choose the best
        maximumStrong candidates, using priority as the first criteria, then
        integer infeasibility.
        
        The algorithm is to fill the choice array with a set of good candidates (by
        infeasibility) with priority bestPriority.  Finding a candidate with
        priority better (less) than bestPriority flushes the choice array. (This
        serves as initialization when the first candidate is found.)
        
        A new candidate is added to choices only if its infeasibility exceeds the
        current max infeasibility (mostAway). When a candidate is added, it
        replaces the candidate with the smallest infeasibility (tracked by
        iSmallest).
      */
      int iSmallest = 0;
      double mostAway = 1.0e-100;
      for (i = 0 ; i < maximumStrong ; i++)
        choice[i].possibleBranch = NULL ;
      numberStrong=0;
      bool canDoOneHot=false;
      for (i=0;i<numberObjects;i++) {
        OsiObject * object = model->modifiableObject(i);
        int preferredWay;
        double infeasibility = object->infeasibility(&usefulInfo,preferredWay);
        int priorityLevel = object->priority();
        if (hotstartSolution) {
          // we are doing hot start
          const CbcSimpleInteger * thisOne = dynamic_cast <const CbcSimpleInteger *> (object);
          if (thisOne) {
            int iColumn = thisOne->columnNumber();
            bool canDoThisHot=true;
            double targetValue = hotstartSolution[iColumn];
            if (saveUpper[iColumn]>saveLower[iColumn]) {
              double value = saveSolution[iColumn];
              if (hotstartPriorities)
                priorityLevel=hotstartPriorities[iColumn]; 
              //double originalLower = thisOne->originalLower();
              //double originalUpper = thisOne->originalUpper();
              // switch off if not possible
              if (targetValue>=saveLower[iColumn]&&targetValue<=saveUpper[iColumn]) {
                /* priority outranks rest always if negative
                   otherwise can be downgraded if at correct level.
                   Infeasibility may be increased to choose 1.0 values first.
                   choose one near wanted value
                */
                if (fabs(value-targetValue)>integerTolerance) {
                  infeasibility = 1.0-fabs(value-targetValue);
                  if (targetValue==1.0)
                    infeasibility += 1.0;
                  if (value>targetValue) {
                    preferredWay=-1;
                  } else {
                    preferredWay=1;
                  }
                  priorityLevel = CoinAbs(priorityLevel);
                } else if (priorityLevel<0) {
                  priorityLevel = CoinAbs(priorityLevel);
                  if (targetValue==saveLower[iColumn]) {
                    infeasibility = integerTolerance+1.0e-12;
                    preferredWay=-1;
                  } else if (targetValue==saveUpper[iColumn]) {
                    infeasibility = integerTolerance+1.0e-12;
                    preferredWay=1;
                  } else {
                    // can't
                    priorityLevel += 10000000;
                    canDoThisHot=false;
                  }
                } else {
                  priorityLevel += 10000000;
                  canDoThisHot=false;
                }
              } else {
                // switch off if not possible
                canDoThisHot=false;
              }
              if (canDoThisHot)
                canDoOneHot=true;
            } else if (targetValue<saveLower[iColumn]||targetValue>saveUpper[iColumn]) {
            }
          } else {
            priorityLevel += 10000000;
          }
        }
        if (infeasibility) {
          // Increase estimated degradation to solution
          estimatedDegradation += CoinMin(object->upEstimate(),object->downEstimate());
          numberUnsatisfied_++;
	  sumInfeasibilities_ += infeasibility;
          // Better priority? Flush choices.
          if (priorityLevel<bestPriority) {
            int j;
            iSmallest=0;
            for (j=0;j<maximumStrong;j++) {
              choice[j].upMovement=0.0;
              delete choice[j].possibleBranch;
              choice[j].possibleBranch=NULL;
            }
            bestPriority = priorityLevel;
            mostAway=1.0e-100;
            numberStrong=0;
          } else if (priorityLevel>bestPriority) {
            continue;
          }
          // Check for suitability based on infeasibility.
          if (infeasibility>mostAway) {
            //add to list
            choice[iSmallest].upMovement=infeasibility;
            delete choice[iSmallest].possibleBranch;
	    CbcSimpleInteger * obj =
	      dynamic_cast <CbcSimpleInteger *>(object) ;
	    if (obj) {
	      choice[iSmallest].possibleBranch=obj->createBranch(solver,&usefulInfo,preferredWay);
	    } else {
	      CbcObject * obj =
		dynamic_cast <CbcObject *>(object) ;
	      assert (obj);
	      choice[iSmallest].possibleBranch=obj->createBranch(preferredWay);
	    }
            numberStrong = CoinMax(numberStrong,iSmallest+1);
            // Save which object it was
            choice[iSmallest].objectNumber=i;
            int j;
            iSmallest=-1;
            mostAway = 1.0e50;
            for (j=0;j<maximumStrong;j++) {
              if (choice[j].upMovement<mostAway) {
                mostAway=choice[j].upMovement;
                iSmallest=j;
              }
            }
          }
        }
      }
      if (!canDoOneHot&&hotstartSolution) {
        // switch off as not possible
        hotstartSolution=NULL;
        model->setHotstartSolution(NULL,NULL);
      }
      if (numberUnsatisfied_) {
        // some infeasibilities - go to next steps
        break;
      } else if (!iPass) {
        // looks like a solution - get paranoid
        bool roundAgain=false;
        // get basis
        CoinWarmStartBasis * ws = dynamic_cast<CoinWarmStartBasis*>(solver->getWarmStart());
        if (!ws)
          break;
        for (i=0;i<numberColumns;i++) {
          double value = saveSolution[i];
          if (value<lower[i]) {
            saveSolution[i]=lower[i];
            roundAgain=true;
            ws->setStructStatus(i,CoinWarmStartBasis::atLowerBound);
          } else if (value>upper[i]) {
            saveSolution[i]=upper[i];
            roundAgain=true;
            ws->setStructStatus(i,CoinWarmStartBasis::atUpperBound);
          } 
        }
        if (roundAgain&&saveNumberStrong) {
          // restore basis
          solver->setWarmStart(ws);
          delete ws;
          solver->resolve();
          memcpy(saveSolution,solver->getColSolution(),numberColumns*sizeof(double));
          model->reserveCurrentSolution(saveSolution);
          if (!solver->isProvenOptimal()) {
            // infeasible 
            anyAction=-2;
            break;
          }
        } else {
          delete ws;
          break;
        }
      }
    }
    /* Some solvers can do the strong branching calculations faster if
       they do them all at once.  At present only Clp does for ordinary
       integers but I think this coding would be easy to modify 
    */
    bool allNormal=true; // to say if we can do fast strong branching
    // Say which one will be best
    int bestChoice=0;
    double worstInfeasibility=0.0;
    for (i=0;i<numberStrong;i++) {
      choice[i].numIntInfeasUp = numberUnsatisfied_;
      choice[i].numIntInfeasDown = numberUnsatisfied_;
      choice[i].fix=0; // say not fixed
      if (!dynamic_cast <const CbcSimpleInteger *> (model->object(choice[i].objectNumber)))
        allNormal=false; // Something odd so lets skip clever fast branching
      if ( !model->object(choice[i].objectNumber)->boundBranch())
        numberStrong=0; // switch off
      if ( choice[i].possibleBranch->numberBranches()>2)
        numberStrong=0; // switch off
      // Do best choice in case switched off
      if (choice[i].upMovement>worstInfeasibility) {
        worstInfeasibility=choice[i].upMovement;
        bestChoice=i;
      }
    }
    // If we have hit max time don't do strong branching
    bool hitMaxTime = ( CoinCpuTime()-model->getDblParam(CbcModel::CbcStartSeconds) > 
                        model->getDblParam(CbcModel::CbcMaximumSeconds));
    // also give up if we are looping round too much
    if (hitMaxTime||numberPassesLeft<=0)
      numberStrong=0;
    /*
      Is strong branching enabled? If so, set up and do it. Otherwise, we'll
      fall through to simple branching.
      
      Setup for strong branching involves saving the current basis (for restoration
      afterwards) and setting up for hot starts.
    */
    if (numberStrong&&saveNumberStrong) {
      
      bool solveAll=false; // set true to say look at all even if some fixed (experiment)
      solveAll=true;
      // worth trying if too many times
      // Save basis
      CoinWarmStart * ws = solver->getWarmStart();
      // save limit
      int saveLimit;
      solver->getIntParam(OsiMaxNumIterationHotStart,saveLimit);
      if (beforeSolution&&saveLimit<100)
        solver->setIntParam(OsiMaxNumIterationHotStart,100); // go to end
#     ifdef COIN_HAS_CLP      
      /* If we are doing all strong branching in one go then we create new arrays
         to store information.  If clp NULL then doing old way.
         Going down -
         outputSolution[2*i] is final solution.
         outputStuff[2*i] is status (0 - finished, 1 infeas, other unknown
         outputStuff[2*i+numberStrong] is number iterations
         On entry newUpper[i] is new upper bound, on exit obj change
         Going up -
         outputSolution[2*i+1] is final solution.
         outputStuff[2*i+1] is status (0 - finished, 1 infeas, other unknown
         outputStuff[2*i+1+numberStrong] is number iterations
       On entry newLower[i] is new lower bound, on exit obj change
      */
      OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (solver);
      ClpSimplex * clp=NULL;
      double * newLower = NULL;
      double * newUpper = NULL;
      double ** outputSolution=NULL;
      int * outputStuff=NULL;
      // Go back to normal way if user wants it
      if (osiclp&&(osiclp->specialOptions()&16)!=0&&osiclp->specialOptions()>0)
        allNormal=false;
      if (osiclp&&!allNormal) {
        // say do fast
        int easy=1;
        osiclp->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,&easy) ;
      }
      if (osiclp&& allNormal) {
        clp = osiclp->getModelPtr();
        // Clp - do a different way
        newLower = new double[numberStrong];
        newUpper = new double[numberStrong];
        outputSolution = new double * [2*numberStrong];
        outputStuff = new int [4*numberStrong];
        int * which = new int[numberStrong];
        int startFinishOptions;
        int specialOptions = osiclp->specialOptions();
        int clpOptions = clp->specialOptions();
        int returnCode=0;
#define CRUNCH
#ifdef CRUNCH
        // Crunch down problem
        int numberRows = clp->numberRows();
        // Use dual region
        double * rhs = clp->dualRowSolution();
        int * whichRow = new int[3*numberRows];
        int * whichColumn = new int[2*numberColumns];
        int nBound;
        ClpSimplex * small = ((ClpSimplexOther *) clp)->crunch(rhs,whichRow,whichColumn,nBound,true);
        if (!small) {
          anyAction=-2;
          //printf("XXXX Inf by inspection\n");
          delete [] whichColumn;
          whichColumn=NULL;
          delete [] whichRow;
          whichRow=NULL;
          break;
        } else {
          clp = small;
        }
#else
        int saveLogLevel = clp->logLevel();
        int saveMaxIts = clp->maximumIterations();
#endif
        clp->setLogLevel(0);
        if((specialOptions&1)==0) {
          startFinishOptions=0;
          clp->setSpecialOptions(clpOptions|(64|1024));
        } else {
          startFinishOptions=1+2+4;
          //startFinishOptions=1+4; // for moment re-factorize
          if((specialOptions&4)==0) 
            clp->setSpecialOptions(clpOptions|(64|128|512|1024|4096));
          else
            clp->setSpecialOptions(clpOptions|(64|128|512|1024|2048|4096));
        }
        // User may want to clean up before strong branching
        if ((clp->specialOptions()&32)!=0) {
          clp->primal(1);
          if (clp->numberIterations())
            model->messageHandler()->message(CBC_ITERATE_STRONG,*model->messagesPointer())
              << clp->numberIterations()
              <<CoinMessageEol;
        }
        clp->setMaximumIterations(saveLimit);
#ifdef CRUNCH
        int * backColumn = whichColumn+numberColumns;
#endif
        for (i=0;i<numberStrong;i++) {
          int iObject = choice[i].objectNumber;
          const OsiObject * object = model->object(iObject);
          const CbcSimpleInteger * simple = dynamic_cast <const CbcSimpleInteger *> (object);
          int iSequence = simple->columnNumber();
          newLower[i]= ceil(saveSolution[iSequence]);
          newUpper[i]= floor(saveSolution[iSequence]);
#ifdef CRUNCH
          iSequence = backColumn[iSequence];
          assert (iSequence>=0);
#endif
          which[i]=iSequence;
          outputSolution[2*i]= new double [numberColumns];
          outputSolution[2*i+1]= new double [numberColumns];
        }
        //clp->writeMps("bad");
        returnCode=clp->strongBranching(numberStrong,which,
                                            newLower, newUpper,outputSolution,
                                            outputStuff,outputStuff+2*numberStrong,!solveAll,false,
                                            startFinishOptions);
#ifndef CRUNCH
        clp->setSpecialOptions(clpOptions); // restore
        clp->setMaximumIterations(saveMaxIts);
        clp->setLogLevel(saveLogLevel);
#endif
        if (returnCode==-2) {
          // bad factorization!!!
          // Doing normal way
          // Mark hot start
          solver->markHotStart();
          clp = NULL;
        } else {
#ifdef CRUNCH
          // extract solution
          //bool checkSol=true;
          for (i=0;i<numberStrong;i++) {
            int iObject = choice[i].objectNumber;
            const OsiObject * object = model->object(iObject);
            const CbcSimpleInteger * simple = dynamic_cast <const CbcSimpleInteger *> (object);
            int iSequence = simple->columnNumber();
            which[i]=iSequence;
            double * sol = outputSolution[2*i];
            double * sol2 = outputSolution[2*i+1];
            //bool x=true;
            //bool x2=true;
            for (int iColumn=numberColumns-1;iColumn>=0;iColumn--) {
              int jColumn = backColumn[iColumn];
              if (jColumn>=0) {
                sol[iColumn]=sol[jColumn];
                sol2[iColumn]=sol2[jColumn];
              } else {
                sol[iColumn]=saveSolution[iColumn];
                sol2[iColumn]=saveSolution[iColumn];
              }
            }
          }
#endif
        }
#ifdef CRUNCH
        delete [] whichColumn;
        delete [] whichRow;
        delete small;
#endif
        delete [] which;
      } else {
        // Doing normal way
        // Mark hot start
        solver->markHotStart();
      }
#     else	/* COIN_HAS_CLP */

      OsiSolverInterface *clp = NULL ;
      double **outputSolution = NULL ;
      int *outputStuff = NULL ;
      double * newLower = NULL ;
      double * newUpper = NULL ;

      solver->markHotStart();

#     endif	/* COIN_HAS_CLP */
      /*
        Open a loop to do the strong branching LPs. For each candidate variable,
        solve an LP with the variable forced down, then up. If a direction turns
        out to be infeasible or monotonic (i.e., over the dual objective cutoff),
        force the objective change to be big (1.0e100). If we determine the problem
        is infeasible, or find a monotone variable, escape the loop.
        
        TODO: The `restore bounds' part might be better encapsulated as an
	unbranch() method. Branching objects more exotic than simple integers
	or cliques might not restrict themselves to variable bounds.

        TODO: Virtuous solvers invalidate the current solution (or give bogus
	results :-) when the bounds are changed out from under them. So we
	need to do all the work associated with finding a new solution before
	restoring the bounds.
      */
      for (i = 0 ; i < numberStrong ; i++)
        { double objectiveChange ;
        double newObjectiveValue=1.0e100;
        // status is 0 finished, 1 infeasible and other
        int iStatus;
        /*
          Try the down direction first. (Specify the initial branching alternative as
          down with a call to way(-1). Each subsequent call to branch() performs the
          specified branch and advances the branch object state to the next branch
          alternative.)
        */
        if (!clp) {
          choice[i].possibleBranch->way(-1) ;
          choice[i].possibleBranch->branch() ;
          bool feasible=true;
          if (checkFeasibility) {
            // check branching did not make infeasible
            int iColumn;
            int numberColumns = solver->getNumCols();
            const double * columnLower = solver->getColLower();
            const double * columnUpper = solver->getColUpper();
            for (iColumn= 0;iColumn<numberColumns;iColumn++) {
              if (columnLower[iColumn]>columnUpper[iColumn]+1.0e-5)
                feasible=false;
            }
          }
          if (feasible) {
            solver->solveFromHotStart() ;
            numberStrongDone++;
            numberStrongIterations += solver->getIterationCount();
            /*
              We now have an estimate of objective degradation that we can use for strong
              branching. If we're over the cutoff, the variable is monotone up.
              If we actually made it to optimality, check for a solution, and if we have
              a good one, call setBestSolution to process it. Note that this may reduce the
              cutoff, so we check again to see if we can declare this variable monotone.
            */
            if (solver->isProvenOptimal())
              iStatus=0; // optimal
            else if (solver->isIterationLimitReached()
                     &&!solver->isDualObjectiveLimitReached())
              iStatus=2; // unknown 
            else
              iStatus=1; // infeasible
            newObjectiveValue = solver->getObjSense()*solver->getObjValue();
            choice[i].numItersDown = solver->getIterationCount();
          } else {
            iStatus=1; // infeasible
            newObjectiveValue = 1.0e100;
            choice[i].numItersDown = 0;
          }
        } else {
          iStatus = outputStuff[2*i];
          choice[i].numItersDown = outputStuff[2*numberStrong+2*i];
          numberStrongDone++;
          numberStrongIterations += choice[i].numItersDown;
          newObjectiveValue = objectiveValue+newUpper[i];
          solver->setColSolution(outputSolution[2*i]);
        }
        objectiveChange = CoinMax(newObjectiveValue  - objectiveValue_,0.0);
        if (!iStatus) {
          choice[i].finishedDown = true ;
          if (newObjectiveValue>=model->getCutoff()) {
            objectiveChange = 1.0e100; // say infeasible
            numberStrongInfeasible++;
          } else {
            // See if integer solution
            if (model->feasibleSolution(choice[i].numIntInfeasDown,
                                        choice[i].numObjInfeasDown)
                &&model->problemFeasibility()->feasible(model,-1)>=0) {
              model->setBestSolution(CBC_STRONGSOL,
                                     newObjectiveValue,
                                     solver->getColSolution()) ;
              // only needed for odd solvers
              newObjectiveValue = solver->getObjSense()*solver->getObjValue();
              objectiveChange = CoinMax(newObjectiveValue-objectiveValue_,0.0) ;
              model->setLastHeuristic(NULL);
              model->incrementUsed(solver->getColSolution());
              if (newObjectiveValue >= model->getCutoff()) {	//  *new* cutoff
                objectiveChange = 1.0e100 ;
                numberStrongInfeasible++;
              }
            }
          }
        } else if (iStatus==1) {
          objectiveChange = 1.0e100 ;
          numberStrongInfeasible++;
        } else {
          // Can't say much as we did not finish
          choice[i].finishedDown = false ;
          numberUnfinished++;
        }
        choice[i].downMovement = objectiveChange ;
        
        // restore bounds
        if (!clp)
          { for (int j=0;j<numberColumns;j++) {
            if (saveLower[j] != lower[j])
              solver->setColLower(j,saveLower[j]);
            if (saveUpper[j] != upper[j])
              solver->setColUpper(j,saveUpper[j]);
          }
          }
        //printf("Down on %d, status is %d, obj %g its %d cost %g finished %d inf %d infobj %d\n",
        //     choice[i].objectNumber,iStatus,newObjectiveValue,choice[i].numItersDown,
        //     choice[i].downMovement,choice[i].finishedDown,choice[i].numIntInfeasDown,
        //     choice[i].numObjInfeasDown);
        
        // repeat the whole exercise, forcing the variable up
        if (!clp) {
          bool feasible=true;
          // If odd branching then maybe just one possibility
          if(choice[i].possibleBranch->numberBranchesLeft()>0) {
            choice[i].possibleBranch->branch();
            if (checkFeasibility) {
              // check branching did not make infeasible
              int iColumn;
              int numberColumns = solver->getNumCols();
              const double * columnLower = solver->getColLower();
              const double * columnUpper = solver->getColUpper();
              for (iColumn= 0;iColumn<numberColumns;iColumn++) {
                if (columnLower[iColumn]>columnUpper[iColumn]+1.0e-5)
                  feasible=false;
              }
            }
          } else {
            // second branch infeasible
            feasible=false;
          }
          if (feasible) {
            solver->solveFromHotStart() ;
            numberStrongDone++;
            numberStrongIterations += solver->getIterationCount();
            /*
              We now have an estimate of objective degradation that we can use for strong
              branching. If we're over the cutoff, the variable is monotone up.
              If we actually made it to optimality, check for a solution, and if we have
              a good one, call setBestSolution to process it. Note that this may reduce the
              cutoff, so we check again to see if we can declare this variable monotone.
            */
            if (solver->isProvenOptimal())
              iStatus=0; // optimal
            else if (solver->isIterationLimitReached()
                     &&!solver->isDualObjectiveLimitReached())
              iStatus=2; // unknown 
            else
              iStatus=1; // infeasible
            newObjectiveValue = solver->getObjSense()*solver->getObjValue();
            choice[i].numItersUp = solver->getIterationCount();
          } else {
            iStatus=1; // infeasible
            newObjectiveValue = 1.0e100;
            choice[i].numItersDown = 0;
          }
        } else {
          iStatus = outputStuff[2*i+1];
          choice[i].numItersUp = outputStuff[2*numberStrong+2*i+1];
          numberStrongDone++;
          numberStrongIterations += choice[i].numItersUp;
          newObjectiveValue = objectiveValue+newLower[i];
          solver->setColSolution(outputSolution[2*i+1]);
        }
        objectiveChange = CoinMax(newObjectiveValue  - objectiveValue_,0.0);
        if (!iStatus) {
          choice[i].finishedUp = true ;
          if (newObjectiveValue>=model->getCutoff()) {
            objectiveChange = 1.0e100; // say infeasible
            numberStrongInfeasible++;
          } else {
            // See if integer solution
            if (model->feasibleSolution(choice[i].numIntInfeasUp,
                                        choice[i].numObjInfeasUp)
                &&model->problemFeasibility()->feasible(model,-1)>=0) {
              model->setBestSolution(CBC_STRONGSOL,
                                     newObjectiveValue,
                                     solver->getColSolution()) ;
              // only needed for odd solvers
              newObjectiveValue = solver->getObjSense()*solver->getObjValue();
              objectiveChange = CoinMax(newObjectiveValue-objectiveValue_,0.0) ;
              model->setLastHeuristic(NULL);
              model->incrementUsed(solver->getColSolution());
              if (newObjectiveValue >= model->getCutoff()) {	//  *new* cutoff
                objectiveChange = 1.0e100 ;
                numberStrongInfeasible++;
              }
            }
          }
        } else if (iStatus==1) {
          objectiveChange = 1.0e100 ;
          numberStrongInfeasible++;
        } else {
          // Can't say much as we did not finish
          choice[i].finishedUp = false ;
          numberUnfinished++;
        }
        choice[i].upMovement = objectiveChange ;
        
        // restore bounds
        if (!clp)
          { for (int j=0;j<numberColumns;j++) {
            if (saveLower[j] != lower[j])
              solver->setColLower(j,saveLower[j]);
            if (saveUpper[j] != upper[j])
              solver->setColUpper(j,saveUpper[j]);
          }
          }
        
        //printf("Up on %d, status is %d, obj %g its %d cost %g finished %d inf %d infobj %d\n",
        //     choice[i].objectNumber,iStatus,newObjectiveValue,choice[i].numItersUp,
        //     choice[i].upMovement,choice[i].finishedUp,choice[i].numIntInfeasUp,
        //     choice[i].numObjInfeasUp);
        
        /*
          End of evaluation for this candidate variable. Possibilities are:
          * Both sides below cutoff; this variable is a candidate for branching.
          * Both sides infeasible or above the objective cutoff: no further action
          here. Break from the evaluation loop and assume the node will be purged
          by the caller.
          * One side below cutoff: Install the branch (i.e., fix the variable). Break
          from the evaluation loop and assume the node will be reoptimised by the
          caller.
        */
        // reset
        choice[i].possibleBranch->resetNumberBranchesLeft();
        if (choice[i].upMovement<1.0e100) {
          if(choice[i].downMovement<1.0e100) {
            // feasible - no action
          } else {
            // up feasible, down infeasible
            anyAction=-1;
            //printf("Down infeasible for choice %d sequence %d\n",i,
            // model->object(choice[i].objectNumber)->columnNumber());
            if (!solveAll) {
              choice[i].possibleBranch->way(1);
              choice[i].possibleBranch->branch();
              break;
            } else {
              choice[i].fix=1;
            }
          }
        } else {
          if(choice[i].downMovement<1.0e100) {
            // down feasible, up infeasible
            anyAction=-1;
            //printf("Up infeasible for choice %d sequence %d\n",i,
            // model->object(choice[i].objectNumber)->columnNumber());
            if (!solveAll) {
              choice[i].possibleBranch->way(-1);
              choice[i].possibleBranch->branch();
              break;
            } else {
              choice[i].fix=-1;
            }
          } else {
            // neither side feasible
            anyAction=-2;
            //printf("Both infeasible for choice %d sequence %d\n",i,
            // model->object(choice[i].objectNumber)->columnNumber());
            break;
          }
        }
        bool hitMaxTime = ( CoinCpuTime()-model->getDblParam(CbcModel::CbcStartSeconds) > 
                            model->getDblParam(CbcModel::CbcMaximumSeconds));
        if (hitMaxTime) {
          numberStrong=i+1;
          break;
        }
        }
      if (!clp) {
        // Delete the snapshot
        solver->unmarkHotStart();
      } else {
        delete [] newLower;
        delete [] newUpper;
        delete [] outputStuff;
        int i;
        for (i=0;i<2*numberStrong;i++)
          delete [] outputSolution[i];
        delete [] outputSolution;
      }
      solver->setIntParam(OsiMaxNumIterationHotStart,saveLimit);
      // restore basis
      solver->setWarmStart(ws);
      // Unless infeasible we will carry on
      // But we could fix anyway
      if (anyAction==-1&&solveAll) {
        // apply and take off
        for (i = 0 ; i < numberStrong ; i++) {
          if (choice[i].fix) {
            choice[i].possibleBranch->way(choice[i].fix) ;
            choice[i].possibleBranch->branch() ;
          }
        }
        bool feasible=true;
        if (checkFeasibility) {
          // check branching did not make infeasible
          int iColumn;
          int numberColumns = solver->getNumCols();
          const double * columnLower = solver->getColLower();
          const double * columnUpper = solver->getColUpper();
          for (iColumn= 0;iColumn<numberColumns;iColumn++) {
            if (columnLower[iColumn]>columnUpper[iColumn]+1.0e-5)
              feasible=false;
          }
        }
        if (feasible) {
          // can do quick optimality check
          int easy=2;
          solver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,&easy) ;
          solver->resolve() ;
          solver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,NULL) ;
          feasible = solver->isProvenOptimal();
        }
        if (feasible) {
          memcpy(saveSolution,solver->getColSolution(),numberColumns*sizeof(double));
          model->reserveCurrentSolution(saveSolution);
          memcpy(saveLower,solver->getColLower(),numberColumns*sizeof(double));
          memcpy(saveUpper,solver->getColUpper(),numberColumns*sizeof(double));
          // Clean up all candidates whih are fixed
          int numberLeft=0;
          for (i = 0 ; i < numberStrong ; i++) {
            CbcStrongInfo thisChoice = choice[i];
            choice[i].possibleBranch=NULL;
            const OsiObject * object = model->object(thisChoice.objectNumber);
            int preferredWay;
            double infeasibility = object->infeasibility(&usefulInfo,preferredWay);
            if (!infeasibility) {
              // take out
              delete thisChoice.possibleBranch;
            } else {
              choice[numberLeft++]=thisChoice;
            }
          }
          numberStrong=numberLeft;
          for (;i<maximumStrong;i++) {
            delete choice[i].possibleBranch;
            choice[i].possibleBranch=NULL;
          }
          // If all fixed then round again
          if (!numberLeft) {
            finished=false;
            numberStrong=0;
            saveNumberStrong=0;
            maximumStrong=1;
          } else {
            anyAction=0;
          }
          // If these two uncommented then different action
          anyAction=-1;
          finished=true;
          //printf("some fixed but continuing %d left\n",numberLeft);
        } else {
          anyAction=-2; // say infeasible
        }
      }
      delete ws;
      int numberNodes = model->getNodeCount();
      // update number of strong iterations etc
      model->incrementStrongInfo(numberStrongDone,numberStrongIterations,
                                 anyAction==-2 ? 0:numberStrongInfeasible,anyAction==-2);
      
      /*
        anyAction >= 0 indicates that strong branching didn't produce any monotone
        variables. Sift through the candidates for the best one.
        
        QUERY: Setting numberNodes looks to be a distributed noop. numberNodes is
        local to this code block. Perhaps should be numberNodes_ from model?
        Unclear what this calculation is doing.
      */
      if (anyAction>=0) {
        
        // get average cost per iteration and assume stopped ones
        // would stop after 50% more iterations at average cost??? !!! ???
        double averageCostPerIteration=0.0;
        double totalNumberIterations=1.0;
        int smallestNumberInfeasibilities=COIN_INT_MAX;
        for (i=0;i<numberStrong;i++) {
          totalNumberIterations += choice[i].numItersDown +
            choice[i].numItersUp ;
          averageCostPerIteration += choice[i].downMovement +
            choice[i].upMovement;
          smallestNumberInfeasibilities= 
            CoinMin(CoinMin(choice[i].numIntInfeasDown ,
                            choice[i].numIntInfeasUp ),
                    smallestNumberInfeasibilities);
        }
        //if (smallestNumberInfeasibilities>=numberIntegerInfeasibilities)
        //numberNodes=1000000; // switch off search for better solution
        numberNodes=1000000; // switch off anyway
        averageCostPerIteration /= totalNumberIterations;
        // all feasible - choose best bet
        
        // New method does all at once so it can be more sophisticated
        // in deciding how to balance actions.
        // But it does need arrays
        double * changeUp = new double [numberStrong];
        int * numberInfeasibilitiesUp = new int [numberStrong];
        double * changeDown = new double [numberStrong];
        int * numberInfeasibilitiesDown = new int [numberStrong];
        CbcBranchingObject ** objects = new CbcBranchingObject * [ numberStrong];
        for (i = 0 ; i < numberStrong ; i++) {
          int iColumn = choice[i].possibleBranch->variable() ;
          model->messageHandler()->message(CBC_STRONG,*model->messagesPointer())
            << i << iColumn
            <<choice[i].downMovement<<choice[i].numIntInfeasDown 
            <<choice[i].upMovement<<choice[i].numIntInfeasUp 
            <<choice[i].possibleBranch->value()
            <<CoinMessageEol;
          changeUp[i]=choice[i].upMovement;
          numberInfeasibilitiesUp[i] = choice[i].numIntInfeasUp;
          changeDown[i]=choice[i].downMovement;
          numberInfeasibilitiesDown[i] = choice[i].numIntInfeasDown;
          objects[i] = choice[i].possibleBranch;
        }
        int whichObject = decision->bestBranch(objects,numberStrong,numberUnsatisfied_,
                                               changeUp,numberInfeasibilitiesUp,
                                               changeDown,numberInfeasibilitiesDown,
                                               objectiveValue_);
        // move branching object and make sure it will not be deleted
        if (whichObject>=0) {
          branch_ = objects[whichObject];
          if (model->messageHandler()->logLevel()>3) 
            printf("Choosing column %d\n",choice[whichObject].possibleBranch->variable()) ;
          choice[whichObject].possibleBranch=NULL;
        }
        delete [] changeUp;
        delete [] numberInfeasibilitiesUp;
        delete [] changeDown;
        delete [] numberInfeasibilitiesDown;
        delete [] objects;
      }
#     ifdef COIN_HAS_CLP
      if (osiclp&&!allNormal) {
        // back to normal
        osiclp->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,NULL) ;
      }
#     endif  
    }
    /*
      Simple branching. Probably just one, but we may have got here
      because of an odd branch e.g. a cut
    */
    else {
      // not strong
      // C) create branching object
      branch_ = choice[bestChoice].possibleBranch;
      choice[bestChoice].possibleBranch=NULL;
    }
  }
  // Set guessed solution value
  guessedObjectiveValue_ = objectiveValue_+estimatedDegradation;
/*
  Cleanup, then we're outta here.
*/
  if (!model->branchingMethod()||dynamicBranchingObject)
    delete decision;
    
  for (i=0;i<maximumStrong;i++)
    delete choice[i].possibleBranch;
  delete [] choice;
  delete [] saveLower;
  delete [] saveUpper;
  
  // restore solution
  solver->setColSolution(saveSolution);
  delete [] saveSolution;
  return anyAction;
}

/*
  Version for dynamic pseudo costs.

  **** For now just return if anything odd
  later allow even if odd

  The routine scans through the object list of the model looking for objects
  that indicate infeasibility. It tests each object using strong branching
  and selects the one with the least objective degradation.  A corresponding
  branching object is left attached to lastNode.
  This version gives preference in evaluation to variables which
  have not been evaluated many times.  It also uses numberStrong
  to say give up if last few tries have not changed incumbent.
  See Achterberg, Koch and Martin.

  If strong branching is disabled, a candidate object is chosen essentially
  at random (whatever object ends up in pos'n 0 of the candidate array).

  If a branching candidate is found to be monotone, bounds are set to fix the
  variable and the routine immediately returns (the caller is expected to
  reoptimize).

  If a branching candidate is found to result in infeasibility in both
  directions, the routine immediately returns an indication of infeasibility.

  Returns:  0	both branch directions are feasible
	   -1	branching variable is monotone
	   -2	infeasible
           -3   Use another method

           For now just fix on objective from strong branching.
*/

int CbcNode::chooseDynamicBranch (CbcModel *model, CbcNode *lastNode,
                                  OsiSolverBranch * & branches,int numberPassesLeft)

{ if (lastNode)
    depth_ = lastNode->depth_+1;
  else
    depth_ = 0;
  // Go to other choose if hot start
  if (model->hotstartSolution()) 
    return -3;
  delete branch_;
  branch_=NULL;
  OsiSolverInterface * solver = model->solver();
  // get information on solver type
  const OsiAuxInfo * auxInfo = solver->getAuxiliaryInfo();
  const OsiBabSolver * auxiliaryInfo = dynamic_cast<const OsiBabSolver *> (auxInfo);
  if (!auxiliaryInfo) {
    // use one from CbcModel
    auxiliaryInfo = model->solverCharacteristics();
  }
  // point to useful information
  OsiBranchingInformation usefulInfo = model->usefulInformation();
  // and modify
  usefulInfo.depth_=depth_;
  assert (auxiliaryInfo);
  //assert(objectiveValue_ == solver->getObjSense()*solver->getObjValue());
  double cutoff =model->getCutoff();
  double distanceToCutoff=cutoff-objectiveValue_;
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  // See what user thinks
  int anyAction=model->problemFeasibility()->feasible(model,0);
  if (anyAction) {
    // will return -2 if infeasible , 0 if treat as integer
    return anyAction-1;
  }
  int i;
  int saveStateOfSearch = model->stateOfSearch();
  int numberStrong=model->numberStrong();
  if (!auxiliaryInfo->warmStart())
    numberStrong=0;
  // But make more likely to get out after some times
  int changeStrategy=numberStrong;
  double changeFactor=1.0;
  // Use minimum of this and one stored in objects
  //int numberBeforeTrust = model->numberBeforeTrust();
  int numberObjects = model->numberObjects();
  bool checkFeasibility = numberObjects>model->numberIntegers();
  // For now return if not simple
  if (checkFeasibility)
    return -3;
  // Return if doing hot start (in BAB sense)
  if (model->hotstartSolution()) 
    return -3;
  //#define RANGING
#ifdef RANGING
  // Pass number
  int kPass=0;
  int numberRows = solver->getNumRows();
#endif
  int numberColumns = model->getNumCols();
  double * saveUpper = new double[numberColumns];
  double * saveLower = new double[numberColumns];

  // Save solution in case heuristics need good solution later
  
  double * saveSolution = new double[numberColumns];
  memcpy(saveSolution,solver->getColSolution(),numberColumns*sizeof(double));
  model->reserveCurrentSolution(saveSolution);
  /*
    Get a branching decision object. Use the default dynamic decision criteria unless
    the user has loaded a decision method into the model.
  */
  CbcBranchDecision *decision = model->branchingMethod();
  if (!decision)
    decision = new CbcBranchDynamicDecision();
  int numberMini=0;
  int xPen=0;
  int xMark=0;
  for (i=0;i<numberColumns;i++) {
    saveLower[i] = lower[i];
    saveUpper[i] = upper[i];
  }
  // Get arrays to sort 
  double * sort = new double[numberObjects];
  int * whichObject = new int[numberObjects];
  int * objectMark = new int[2*numberObjects+1];
  // Arrays with movements
  double * upEstimate = new double[numberObjects];
  double * downEstimate = new double[numberObjects];
  CbcStrongInfo * fixObject = new CbcStrongInfo[numberObjects];
  double estimatedDegradation=0.0; 
  int numberNodes=model->getNodeCount();
  int saveLogLevel = model->logLevel();
  if ((numberNodes%500)==0&&false) {
    model->setLogLevel(6);
    // Get average up and down costs
    double averageUp=0.0;
    double averageDown=0.0;
    int numberUp=0;
    int numberDown=0;
    int i;
    for ( i=0;i<numberObjects;i++) {
      OsiObject * object = model->modifiableObject(i);
      CbcSimpleIntegerDynamicPseudoCost * dynamicObject =
	dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object) ;
      assert(dynamicObject);
      int  numberUp2=0;
      int numberDown2=0;
      double up=0.0;
      double down=0.0;
      if (dynamicObject->numberTimesUp()) {
	numberUp++;
	averageUp += dynamicObject->upDynamicPseudoCost();
	numberUp2 += dynamicObject->numberTimesUp();
	up = dynamicObject->upDynamicPseudoCost();
      }
      if (dynamicObject->numberTimesDown()) {
	numberDown++;
	averageDown += dynamicObject->downDynamicPseudoCost();
	numberDown2 += dynamicObject->numberTimesDown();
	down = dynamicObject->downDynamicPseudoCost();
      }
      if (numberUp2||numberDown2)
	printf("col %d - up %d times cost %g, - down %d times cost %g\n",
	       dynamicObject->columnNumber(),numberUp2,up,numberDown2,down);
    }
    if (numberUp) 
      averageUp /= (double) numberUp;
    else
      averageUp=1.0;
    if (numberDown) 
      averageDown /= (double) numberDown;
    else
      averageDown=1.0;
    printf("total - up %d vars average %g, - down %d vars average %g\n",
	   numberUp,averageUp,numberDown,averageDown);
  }
  int numberBeforeTrust = model->numberBeforeTrust();
  int numberPenalties = model->numberPenalties();
  if (numberBeforeTrust>=1000000) {
    numberBeforeTrust = numberBeforeTrust % 1000000;
    numberPenalties=0;
  } else if (numberBeforeTrust<0) {
    if (numberBeforeTrust==-1)
      numberPenalties=numberColumns;
    else if (numberBeforeTrust==-2)
      numberPenalties=0;
    numberBeforeTrust=0;
  }
  // May go round twice if strong branching fixes all local candidates
  bool finished=false;
  int numberToFix=0;
# ifdef COIN_HAS_CLP
  OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (solver);
  int saveClpOptions=0;
  if (osiclp) {
    // for faster hot start
    saveClpOptions = osiclp->specialOptions();
    osiclp->setSpecialOptions(saveClpOptions|8192);
  }
# else
  OsiSolverInterface *osiclp = 0 ;
# endif
  const CglTreeProbingInfo * probingInfo = model->probingInfo();
  int saveSearchStrategy2 = model->searchStrategy();
  if (saveSearchStrategy2<999) {
    // Get average up and down costs
    double averageUp=0.0;
    double averageDown=0.0;
    {
      int numberUp=0;
      int numberDown=0;
      int i;
      for ( i=0;i<numberObjects;i++) {
        OsiObject * object = model->modifiableObject(i);
        CbcSimpleIntegerDynamicPseudoCost * dynamicObject =
          dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object) ;
        assert(dynamicObject);
        if (dynamicObject->numberTimesUp()) {
          numberUp++;
          averageUp += dynamicObject->upDynamicPseudoCost();
        }
        if (dynamicObject->numberTimesDown()) {
          numberDown++;
          averageDown += dynamicObject->downDynamicPseudoCost();
        }
      }
      if (numberUp) 
        averageUp /= (double) numberUp;
      else
        averageUp=1.0;
      if (numberDown) 
        averageDown /= (double) numberDown;
      else
        averageDown=1.0;
      for ( i=0;i<numberObjects;i++) {
        OsiObject * object = model->modifiableObject(i);
        CbcSimpleIntegerDynamicPseudoCost * dynamicObject =
          dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object) ;
        assert(dynamicObject);
        if (!dynamicObject->numberTimesUp()) 
          dynamicObject->setUpDynamicPseudoCost(averageUp);
      if (!dynamicObject->numberTimesDown()) 
        dynamicObject->setDownDynamicPseudoCost(averageDown);
      }
    }
  } else if (saveSearchStrategy2<1999) {
    // pseudo shadow prices
    model->pseudoShadow(NULL,NULL);
  } else if (saveSearchStrategy2<2999) {
    // leave old ones
  } else if (saveSearchStrategy2<3999) {
    // pseudo shadow prices at root
    if (!numberNodes)
      model->pseudoShadow(NULL,NULL);
  } else {
    abort();
  }
  if (saveSearchStrategy2>=0)
    saveSearchStrategy2 = saveSearchStrategy2 % 1000;
  if (saveSearchStrategy2==999)
    saveSearchStrategy2=-1;
  int px[4]={-1,-1,-1,-1};
  int saveSearchStrategy = saveSearchStrategy2<99 ? saveSearchStrategy2 : saveSearchStrategy2-100;
  bool newWay = saveSearchStrategy2>98;
  int numberNotTrusted=0;
  int numberStrongDone=0;
  int numberUnfinished=0;
  int numberStrongInfeasible=0;
  int numberStrongIterations=0;
  // so we can save lots of news
  CbcStrongInfo choice;
  CbcDynamicPseudoCostBranchingObject * choiceObject = NULL;
  if (model->allDynamic()) {
    CbcSimpleIntegerDynamicPseudoCost * object = NULL;
    choiceObject=new CbcDynamicPseudoCostBranchingObject(model,0,-1,0.5,object);
  }
  choice.possibleBranch=choiceObject;
  while(!finished) {
    finished=true;
    decision->initialize(model);
    // Some objects may compute an estimate of best solution from here
    estimatedDegradation=0.0; 
    numberToFix=0;
    int numberIntegerInfeasibilities=0; // without odd ones
    int numberToDo=0;
    int iBestNot=-1;
    int iBestGot=-1;
    double best=0.0;
    numberNotTrusted=0;
    numberStrongDone=0;
    numberUnfinished=0;
    numberStrongInfeasible=0;
    numberStrongIterations=0;
    int * which = objectMark+numberObjects+1;
    int neededPenalties;
    int branchingMethod=-1;
    // We may go round this loop three times (only if we think we have solution)
    for (int iPass=0;iPass<3;iPass++) {
      
      // compute current state
      int numberObjectInfeasibilities; // just odd ones
      model->feasibleSolution(
                              numberIntegerInfeasibilities,
                              numberObjectInfeasibilities);
      
      // Some objects may compute an estimate of best solution from here
      estimatedDegradation=0.0; 
      numberUnsatisfied_ = 0;
      // initialize sum of "infeasibilities"
      sumInfeasibilities_ = 0.0;
      int bestPriority=COIN_INT_MAX;
      int number01 = 0;
      const fixEntry * entry = NULL;
      const int * toZero = NULL;
      const int * toOne = NULL;
      const int * backward = NULL;
      int numberUnsatisProbed=0;
      int numberUnsatisNotProbed=0; // 0-1
      if (probingInfo) {
	number01 = probingInfo->numberIntegers();
	entry = probingInfo->fixEntries();
	toZero = probingInfo->toZero();
	toOne = probingInfo->toOne();
	backward = probingInfo->backward();
	if (!toZero[number01]||number01<numberObjects||true) {
	  // no info
	  probingInfo=NULL;
	}
      }
      /*
        Scan for branching objects that indicate infeasibility. Choose candidates
        using priority as the first criteria, then integer infeasibility.
        
        The algorithm is to fill the array with a set of good candidates (by
        infeasibility) with priority bestPriority.  Finding a candidate with
        priority better (less) than bestPriority flushes the choice array. (This
        serves as initialization when the first candidate is found.)
        
      */
      numberToDo=0;
      neededPenalties=0;
      iBestNot=-1;
      double bestNot=0.0;
      iBestGot=-1;
      best=0.0;
      /* Problem type as set by user or found by analysis.  This will be extended
	 0 - not known
	 1 - Set partitioning <=
	 2 - Set partitioning ==
	 3 - Set covering
	 4 - all +- 1 or all +1 and odd
      */
      int problemType = model->problemType();
#define PRINT_STUFF -1
      for (i=0;i<numberObjects;i++) {
        OsiObject * object = model->modifiableObject(i);
        CbcSimpleIntegerDynamicPseudoCost * dynamicObject =
          dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object) ;
        assert(dynamicObject);
        int preferredWay;
        double infeasibility = object->infeasibility(&usefulInfo,preferredWay);
        int priorityLevel = object->priority();
#define ZERO_ONE 0
#define ZERO_FAKE 1.0e20;
#if ZERO_ONE==1
        // branch on 0-1 first (temp)
        if (fabs(saveSolution[dynamicObject->columnNumber()])<1.0)
          priorityLevel--;
#endif
#if ZERO_ONE==2
        if (fabs(saveSolution[dynamicObject->columnNumber()])<1.0)
          infeasibility *= ZERO_FAKE;
#endif
        if (infeasibility) {
	  // check branching method
	  if (branchingMethod!=dynamicObject->method()) {
	    if (branchingMethod==-1)
	      branchingMethod = dynamicObject->method();
	    else
	      branchingMethod = 100;
	  }
          int iColumn = dynamicObject->columnNumber();
          //double gap = saveUpper[iColumn]-saveLower[iColumn];
          // Give precedence to ones with gap of 1.0 
          //assert(gap>0.0);
          //infeasibility /= CoinMin(gap,100.0);
          if (!depth_&&false) {
            // try closest to 0.5
            double part =saveSolution[iColumn]-floor(saveSolution[iColumn]);
            infeasibility = fabs(0.5-part);
          }
          if (problemType>0&&problemType<4&&false) {
            // try closest to 0.5
            double part =saveSolution[iColumn]-floor(saveSolution[iColumn]);
            infeasibility = 0.5-fabs(0.5-part);
          }
	  if (probingInfo) {
	    int iSeq = backward[iColumn];
	    assert (iSeq>=0);
	    infeasibility = 1.0 + (toZero[iSeq+1]-toZero[iSeq])+ 
	      5.0*CoinMin(toOne[iSeq]-toZero[iSeq],toZero[iSeq+1]-toOne[iSeq]);
	    if (toZero[iSeq+1]>toZero[iSeq]) {
	      numberUnsatisProbed++;
	    } else {
	      numberUnsatisNotProbed++;
	    }
	  }
          bool gotDown=false;
          int numberThisDown = dynamicObject->numberTimesDown();
          if (numberThisDown>=numberBeforeTrust)
            gotDown=true;
          bool gotUp=false;
          int numberThisUp = dynamicObject->numberTimesUp();
          if (numberThisUp>=numberBeforeTrust)
            gotUp=true;
          if ((numberNodes%PRINT_STUFF)==0&&PRINT_STUFF>0)
            printf("%d down %d %g up %d %g - infeas %g\n",
                   i,numberThisDown,object->downEstimate(),numberThisUp,object->upEstimate(),
                   infeasibility);
          // Increase estimated degradation to solution
          estimatedDegradation += CoinMin(object->upEstimate(),object->downEstimate());
          downEstimate[i]=object->downEstimate();
          upEstimate[i]=object->upEstimate();
          numberUnsatisfied_++;
	  sumInfeasibilities_ += infeasibility;
          // Better priority? Flush choices.
          if (priorityLevel<bestPriority) {
            numberToDo=0;
            bestPriority = priorityLevel;
            iBestGot=-1;
            best=0.0;
            numberNotTrusted=0;
          } else if (priorityLevel>bestPriority) {
            continue;
          }
          if (!gotUp||!gotDown)
            numberNotTrusted++;
          // Check for suitability based on infeasibility.
          if ((gotDown&&gotUp)&&numberStrong>0) {
            sort[numberToDo]=-infeasibility;
            if (infeasibility>best) {
              best=infeasibility;
              iBestGot=numberToDo;
            }
          } else {
            objectMark[neededPenalties]=numberToDo;
            which[neededPenalties++]=dynamicObject->columnNumber();
            int iColumn = dynamicObject->columnNumber();
            double part =saveSolution[iColumn]-floor(saveSolution[iColumn]);
            sort[numberToDo]=-10.0*infeasibility;
            if (!(numberThisUp+numberThisDown))
              sort[numberToDo] *= 100.0; // make even more likely
            if (1.0-fabs(part-0.5)>bestNot) {
              iBestNot=numberToDo;
              bestNot = 1.0-fabs(part-0.5);
            }
          }
	  if (model->messageHandler()->logLevel()>3) { 
            int iColumn = dynamicObject->columnNumber();
	    printf("%d (%d) down %d %g up %d %g - infeas %g - sort %g solution %g\n",
		   i,iColumn,numberThisDown,object->downEstimate(),numberThisUp,object->upEstimate(),
		   infeasibility,sort[numberToDo],saveSolution[iColumn]);
	  }
          whichObject[numberToDo++]=i;
        } else {
          // for debug
          downEstimate[i]=-1.0;
          upEstimate[i]=-1.0;
        }
      }
      if (numberUnsatisfied_) {
	if (probingInfo&&false)
	  printf("nunsat %d, %d probed, %d other 0-1\n",numberUnsatisfied_,
		 numberUnsatisProbed,numberUnsatisNotProbed);
        // some infeasibilities - go to next steps
        break;
      } else if (!iPass) {
        // may just need resolve
        solver->resolve();
        memcpy(saveSolution,solver->getColSolution(),numberColumns*sizeof(double));
        model->reserveCurrentSolution(saveSolution);
        if (!solver->isProvenOptimal()) {
          // infeasible 
          anyAction=-2;
          break;
        }
      } else if (iPass==1) {
        // looks like a solution - get paranoid
        bool roundAgain=false;
        // get basis
        CoinWarmStartBasis * ws = dynamic_cast<CoinWarmStartBasis*>(solver->getWarmStart());
        if (!ws)
          break;
        double tolerance;
        solver->getDblParam(OsiPrimalTolerance,tolerance);
        for (i=0;i<numberColumns;i++) {
          double value = saveSolution[i];
          if (value<lower[i]-tolerance) {
            saveSolution[i]=lower[i];
            roundAgain=true;
            ws->setStructStatus(i,CoinWarmStartBasis::atLowerBound);
          } else if (value>upper[i]+tolerance) {
            saveSolution[i]=upper[i];
            roundAgain=true;
            ws->setStructStatus(i,CoinWarmStartBasis::atUpperBound);
          } 
        }
        if (roundAgain) {
          // restore basis
          solver->setWarmStart(ws);
          solver->setColSolution(saveSolution);
          delete ws;
          bool takeHint;
          OsiHintStrength strength;
          solver->getHintParam(OsiDoDualInResolve,takeHint,strength);
          solver->setHintParam(OsiDoDualInResolve,false,OsiHintDo) ;
          solver->resolve();
          solver->setHintParam(OsiDoDualInResolve,takeHint,strength) ;
          memcpy(saveSolution,solver->getColSolution(),numberColumns*sizeof(double));
          model->reserveCurrentSolution(saveSolution);
          if (!solver->isProvenOptimal()) {
            // infeasible 
            anyAction=-2;
            break;
          }
        } else {
          delete ws;
          break;
        }
      }
    }
    if (anyAction==-2) {
      break;
    }
    bool solveAll=false; // set true to say look at all even if some fixed (experiment)
    solveAll=true;
    // skip if solution
    if (!numberUnsatisfied_)
      break;
    //bool skipAll = (numberBeforeTrust>20&&numberNodes>20000&&numberNotTrusted==0);
    bool skipAll = numberNotTrusted==0||numberToDo==1;
    bool doneHotStart=false;
    int searchStrategy = saveSearchStrategy>=0 ? (saveSearchStrategy%10) : -1;
#ifndef CBC_WEAK_STRONG
    if (((numberNodes%20)==0&&searchStrategy!=2)||(model->specialOptions()&8)!=0)
      skipAll=false;
#endif
    if (!newWay) {
    // 10 up always use %10, 20 up as 10 and allow penalties
    // But adjust depending on ratio of iterations
    if (searchStrategy>0&&saveSearchStrategy<10) {
      if (numberBeforeTrust>=5&&numberBeforeTrust<=10) {
        if (searchStrategy!=2) {
	  if (depth_>5) {
	    int numberIterations = model->getIterationCount();
	    int numberStrongIterations = model->numberStrongIterations();
	    if (numberStrongIterations>numberIterations+10000) {
	      searchStrategy=2;
	      //skipAll=true;
	    } else if (numberStrongIterations*4+1000<numberIterations||depth_<5) {
	      searchStrategy=3;
	      skipAll=false;
	    }
	  } else {
	    searchStrategy=3;
	    skipAll=false;
	  }
        } else {
          //skipAll=true;
        }
      }
    }
    } else {
    // But adjust depending on ratio of iterations
    if (saveSearchStrategy<0) {
      // unset
      if ((numberNodes%20)==0||(model->specialOptions()&8)!=0) {
        // Do numberStrong
        searchStrategy=3;
      } else if (depth_<5) {
        // Do numberStrong
        searchStrategy=2;
      } else {
        int numberIterations = model->getIterationCount();
        int numberStrongIterations = model->numberStrongIterations();
        int numberRows = solver->getNumRows();
        if (numberStrongIterations>numberIterations+CoinMin(10000,10*numberRows)) {
          // off
          searchStrategy=0;
        } else if (numberStrongIterations*4+1000<numberIterations) {
          // Do numberStrong if not trusted
          searchStrategy=2;
        } else {
          searchStrategy=1;
        }
      }
    }
    if (searchStrategy<3&&(!numberNotTrusted||!searchStrategy))
      skipAll=true;
    else
      skipAll=false;
    }
    // worth trying if too many times
    // Save basis
    CoinWarmStart * ws = NULL;
    // save limit
    int saveLimit=0;
    solver->getIntParam(OsiMaxNumIterationHotStart,saveLimit);
    if (!skipAll) {
      ws = solver->getWarmStart();
      int limit=100;
#if 0
      int averageBranchIterations = model->getIterationCount()/(model->getNodeCount()+1);
      if (numberNodes)
        limit = CoinMin(CoinMax(limit,2*averageBranchIterations),500);
      else
        limit = 500;
#endif
      if ((!saveStateOfSearch||searchStrategy>3)&&saveLimit<limit&&saveLimit==100)
        solver->setIntParam(OsiMaxNumIterationHotStart,limit); 
    }
    // Say which one will be best
    int whichChoice=0;
    int bestChoice;
    if (iBestGot>=0)
      bestChoice=iBestGot;
    else
      bestChoice=iBestNot;
    assert (bestChoice>=0);
    // If we have hit max time don't do strong branching
    bool hitMaxTime = ( CoinCpuTime()-model->getDblParam(CbcModel::CbcStartSeconds) > 
                        model->getDblParam(CbcModel::CbcMaximumSeconds));
    // also give up if we are looping round too much
    if (hitMaxTime||numberPassesLeft<=0||(!numberNotTrusted&&false)||branchingMethod==11) {
      int iObject = whichObject[bestChoice];
      OsiObject * object = model->modifiableObject(iObject);
      int preferredWay;
      object->infeasibility(&usefulInfo,preferredWay);
      CbcSimpleInteger * obj =
	dynamic_cast <CbcSimpleInteger *>(object) ;
      if (obj) {
	branch_=obj->createBranch(solver,&usefulInfo,preferredWay);
      } else {
	CbcObject * obj =
	  dynamic_cast <CbcObject *>(object) ;
	assert (obj);
	branch_=obj->createBranch(preferredWay);
      }
      {
	CbcBranchingObject * branchObj =
	  dynamic_cast <CbcBranchingObject *>(branch_) ;
	assert (branchObj);
	branchObj->way(preferredWay);
      }
      delete ws;
      ws=NULL;
      break;
    } else {
      // say do fast
      int easy=1;
      if (!skipAll)
        solver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,&easy) ;
      int iDo;
#ifdef RANGING
      if ((skipAll&&numberBeforeTrust&&saveSearchStrategy<20)||saveSearchStrategy<10)
        numberPenalties=0;
      {
        // off penalties if too much
        double needed = neededPenalties;
        needed *= numberRows;
        if (needed>1.0e6&&numberNodes&&saveSearchStrategy<20) {
          numberPenalties=0;
          neededPenalties=0;
        }
      }
#     ifdef COIN_HAS_CLP
      if (osiclp&&numberPenalties&&neededPenalties) {
        xPen += neededPenalties;
        which--;
        which[0]=neededPenalties;
        osiclp->passInRanges(which);
        // Mark hot start and get ranges
        if (kPass) {
          // until can work out why solution can go funny
          int save = osiclp->specialOptions();
          osiclp->setSpecialOptions(save|256);
          solver->markHotStart();
          osiclp->setSpecialOptions(save);
        } else {
          solver->markHotStart();
        }
        assert (auxiliaryInfo->warmStart());
        doneHotStart=true;
        xMark++;
        kPass++;
        osiclp->passInRanges(NULL);
        const double * downCost=osiclp->upRange();
        const double * upCost=osiclp->downRange();
        //printf("numberTodo %d needed %d numberpenalties %d\n",numberToDo,neededPenalties,numberPenalties);
        double invTrust = 1.0/((double) numberBeforeTrust);
        for (int i=0;i<neededPenalties;i++) {
          int j = objectMark[i];
          int iObject = whichObject[j];
          OsiObject * object = model->modifiableObject(iObject);
          CbcSimpleIntegerDynamicPseudoCost * dynamicObject =
            dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object) ;
          int iSequence=dynamicObject->columnNumber();
          double value = saveSolution[iSequence];
          value -= floor(value);
          double upPenalty = CoinMin(upCost[i],1.0e110)*(1.0-value);
          double downPenalty = CoinMin(downCost[i],1.0e110)*value;
          if (!numberBeforeTrust) {
            // override
            downEstimate[iObject]=downPenalty;
            upEstimate[iObject]=upPenalty;
          } else {
            int numberThisDown = dynamicObject->numberTimesDown();
            if (numberThisDown<numberBeforeTrust) {
              double fraction = ((double) numberThisDown)*invTrust;
              downEstimate[iObject] = fraction*downEstimate[iObject]+(1.0-fraction)*downPenalty;
            }
            int numberThisUp = dynamicObject->numberTimesUp();
            if (numberThisUp<numberBeforeTrust) {
              double fraction = ((double) numberThisUp)*invTrust;
              upEstimate[iObject] = fraction*upEstimate[iObject]+(1.0-fraction)*upPenalty;
            }
          }
          sort[j] = - CoinMin(downEstimate[iObject],upEstimate[iObject]);
#ifdef CBC_WEAK_STRONG
          sort[j] -= 1.0e10; // make more likely to be chosen
#endif
          //if ((numberNodes%PRINT_STUFF)==0&&PRINT_STUFF>0)
          if (!numberNodes)
            printf("%d pen down ps %g -> %g up ps %g -> %g\n",
                   iObject,downCost[i],downPenalty,upCost[i],upPenalty);
        }
      } else
#     endif	/* COIN_HAS_CLP */
      {
        if (!skipAll) {
          // Mark hot start
          solver->markHotStart();
          doneHotStart=true;
          assert (auxiliaryInfo->warmStart());
          xMark++;
          //if (solver->isProvenPrimalInfeasible())
          //printf("**** Hot start says node infeasible\n");
        }
        // make sure best will be first
        if (iBestGot>=0)
          sort[iBestGot]=-1.0e120;
      }
#else		/* RANGING */
      if (!skipAll) {
        // Mark hot start
        doneHotStart=true;
        assert (auxiliaryInfo->warmStart());
        solver->markHotStart();
        xMark++;
      }
      // make sure best will be first
      if (iBestGot>=0)
        sort[iBestGot]=-COIN_DBL_MAX;
#endif		/* RANGING */
      // Actions 0 - exit for repeat, 1 resolve and try old choice,2 exit for continue
#define ACTION 0 
#if ACTION<2
      if (anyAction)
        numberToDo=0; // skip as we will be trying again
#endif
      // Sort 
      CoinSort_2(sort,sort+numberToDo,whichObject);
      // Change in objective opposite infeasible
      double worstFeasible=0.0;
      // Just first if strong off
      if (!numberStrong)
        numberToDo=CoinMin(numberToDo,1);
      iDo=0;
      int saveLimit2;
      solver->getIntParam(OsiMaxNumIterationHotStart,saveLimit2);
      bool doQuickly = false; // numberToDo>2*numberStrong;
      if (searchStrategy==2)
	doQuickly=true;
      //printf("todo %d, strong %d\n",numberToDo,numberStrong);
      int numberTest=numberNotTrusted>0 ? numberStrong : (numberStrong+1)/2;
      int numberTest2 = 2*numberStrong;
      //double distanceToCutoff2 = model->getCutoff()-objectiveValue_;
      if (!newWay) {
      if (searchStrategy==3) {
        // Previously decided we need strong
        doQuickly=false;
        numberTest = numberStrong;
        //numberTest2 = 1000000;
      }
      //if (searchStrategy<0||searchStrategy==1)
        //numberTest2 = 1000000;
#if 0
      if (numberBeforeTrust>20&&(numberNodes>20000||(numberNodes>200&&numberNotTrusted==0))) {
        if ((numberNodes%20)!=0) {
          numberTest=0;
          doQuickly=true;
        }
      }
#else
      // Try nearly always off
      if (searchStrategy<2) {
        if ((numberNodes%20)!=0) {
          if ((model->specialOptions()&8)==0) {
            numberTest=0;
            doQuickly=true;
          }
        } else {
          doQuickly=false;
          numberTest=2*numberStrong;
          skipAll=false;
        }
      } else if (searchStrategy!=3) {
        doQuickly=true;
        numberTest=numberStrong;
      }
#endif
      if (depth_<8&&numberStrong) {
        if (searchStrategy!=2) {
          doQuickly=false;
	  int numberRows = solver->getNumRows();
	  // whether to do this or not is important - think
	  if (numberRows<300||numberRows+numberColumns<2500) {
	    if (depth_<7)
	      numberStrong = CoinMin(3*numberStrong,numberToDo);
	    if (!depth_) 
	      numberStrong=CoinMin(6*numberStrong,numberToDo);
	  }
          numberTest=numberStrong;
          skipAll=false;
        }
        //model->setStateOfSearch(2); // use min min
      }
      // could adjust using average iterations per branch
      // double average = ((double)model->getIterationCount())/
      //((double) model->getNodeCount()+1.0);
      // if too many and big then just do 10 its
      if (!skipAll&&saveStateOfSearch) {
        //if (numberNotTrusted>3*numberStrong&&numberRows>250&&numberColumns>1000&&saveLimit==100)
          // off solver->setIntParam(OsiMaxNumIterationHotStart,10); 
      }
      // make negative for test
      distanceToCutoff = - distanceToCutoff;
      if (numberObjects>-100) {
        // larger 
        distanceToCutoff *= 100.0;
      }
        distanceToCutoff = -COIN_DBL_MAX;
      // Do at least 5 strong
      if (numberColumns<1000&&(depth_<15||numberNodes<1000000))
        numberTest = CoinMax(numberTest,5);
      if ((model->specialOptions()&8)==0) {
        if (skipAll) {
          numberTest=0;
          doQuickly=true;
        }
      } else {
        // do 5 as strong is fixing
        numberTest = CoinMax(numberTest,5);
      }
      } else {
      int numberTest=numberNotTrusted>0 ? numberStrong : (numberStrong+1)/2;
      int numberTest2 = 2*numberStrong;
      if (searchStrategy>=3) {
        // Previously decided we need strong
        doQuickly=false;
        if (depth_<7)
          numberStrong *=3;
        if (!depth_) 
          numberStrong=CoinMin(6*numberStrong,numberToDo);
        numberTest = numberStrong;
        numberTest2 *= 2;
      } else if (searchStrategy==2||(searchStrategy==1&&depth_<6)) {
        numberStrong *=2;
        if (!depth_) 
          numberStrong=CoinMin(2*numberStrong,numberToDo);
        numberTest = numberStrong;
      } else if (searchStrategy==1&&numberNotTrusted) {
        numberTest = numberStrong;
      } else {
        numberTest=0;
        skipAll=true;
      }
      distanceToCutoff=model->getCutoff()-objectiveValue_;
      // make negative for test
      distanceToCutoff = - distanceToCutoff;
      if (numberObjects>-100) {
        // larger 
        distanceToCutoff *= 100.0;
      }
      distanceToCutoff = -COIN_DBL_MAX;
      if (skipAll) {
        numberTest=0;
        doQuickly=true;
      }
      }
#if 0
      // temp - always switch off
      if (0) {
        int numberIterations = model->getIterationCount();
        int numberStrongIterations = model->numberStrongIterations();
        if (numberStrongIterations>numberIterations+10000&&depth_>=5) {
          skipAll=true;
          newWay=false;
          numberTest=0;
          doQuickly=true;
        }
      }
      // temp - always switch on
      if (0) {
        int numberIterations = model->getIterationCount();
        int numberStrongIterations = model->numberStrongIterations();
        if (2*numberStrongIterations<numberIterations||depth_<=5) {
          skipAll=false;
          newWay=false;
          numberTest=CoinMax(numberTest,numberStrong);
          doQuickly=false;
        }
      }
#endif
      px[0]=numberTest;
      px[1]=numberTest2;
      px[2]= doQuickly ? 1 : -1;
      px[3]=numberStrong;
      if (!newWay) {
	if (numberColumns>8*solver->getNumRows()&&false) {
	  printf("skipAll %c doQuickly %c numberTest %d numberTest2 %d numberNot %d\n",
		 skipAll ? 'Y' : 'N',doQuickly ? 'Y' : 'N',numberTest,numberTest2,numberNotTrusted);
	  numberTest = CoinMin(numberTest,model->numberStrong());
	  numberTest2 = CoinMin(numberTest2,model->numberStrong());
	  printf("new test,test2 %d %d\n",numberTest,numberTest2);
	}
      }
      //printf("skipAll %c doQuickly %c numberTest %d numberTest2 %d numberNot %d\n",
      //   skipAll ? 'Y' : 'N',doQuickly ? 'Y' : 'N',numberTest,numberTest2,numberNotTrusted);
      // See if we want mini tree
      bool wantMiniTree=false;
      if (model->sizeMiniTree()&&depth_>7&&saveStateOfSearch>0)
        wantMiniTree=true;
      numberMini=0;
      //if (skipAll&&numberTest==0&&doQuickly)
      //numberToDo = 1; // trust previous stuff
      bool couldChooseFirst = false ; //(skipAll&&numberTest==0&&doQuickly);
      //skipAll=false;
      for ( iDo=0;iDo<numberToDo;iDo++) {
        int iObject = whichObject[iDo];
        OsiObject * object = model->modifiableObject(iObject);
        CbcSimpleIntegerDynamicPseudoCost * dynamicObject =
          dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object) ;
        int iColumn = dynamicObject->columnNumber();
        int preferredWay;
	double infeasibility = object->infeasibility(&usefulInfo,preferredWay);
	// may have become feasible
	if (!infeasibility)
	  continue;
	CbcSimpleInteger * obj =
	  dynamic_cast <CbcSimpleInteger *>(object) ;
	if (obj) {
	  if (choiceObject) {
	    obj->fillCreateBranch(choiceObject,&usefulInfo,preferredWay);
	    choiceObject->setObject(dynamicObject);
	  } else {
	    choice.possibleBranch=obj->createBranch(solver,&usefulInfo,preferredWay);
	  }
	} else {
	  CbcObject * obj =
	    dynamic_cast <CbcObject *>(object) ;
	  assert (obj);
	  choice.possibleBranch=obj->createBranch(preferredWay);
	}
        // Save which object it was
        choice.objectNumber=iObject;
        choice.numIntInfeasUp = numberUnsatisfied_;
        choice.numIntInfeasDown = numberUnsatisfied_;
        choice.upMovement = upEstimate[iObject];
        choice.downMovement = downEstimate[iObject];
        assert (choice.upMovement>=0.0);
        assert (choice.downMovement>=0.0);
        choice.fix=0; // say not fixed
	double maxChange = 0.5*(choice.upMovement+choice.downMovement);
	maxChange = CoinMin(choice.upMovement,choice.downMovement);
	maxChange = CoinMax(choice.upMovement,choice.downMovement);
	if (searchStrategy==2)
	  maxChange = COIN_DBL_MAX;
	//maxChange *= 5.0;
	if (dynamicObject->method()==1)
	  maxChange *= 0.1; // probing
        // see if can skip strong branching
        int canSkip = choice.possibleBranch->fillStrongInfo(choice);
#if 0
        if (!newWay) {
	  if ((maxChange>distanceToCutoff2)&&(!doQuickly||(numberTest>0&&searchStrategy!=2)))
          canSkip=0;
        } else {
        if (skipAll)
          canSkip=1;
        else if (numberTest>0&&searchStrategy>=3)
          canSkip=0;
        }
        if (!numberBeforeTrust) {
          canSkip=1;
        }
        if (sort[iDo]<distanceToCutoff)
          canSkip=0;
        if (((numberTest2<=0&&numberTest<=0)||skipAll)&&sort[iDo]>distanceToCutoff) {
          canSkip=1; // always skip
          if (iDo>20) {
	    if (!choiceObject) {
	      delete choice.possibleBranch;
	      choice.possibleBranch=NULL;
	    }
            break; // give up anyway
          }
        }
#else
        if (((numberTest2<=0&&numberTest<=0)||skipAll)&&sort[iDo]>distanceToCutoff) {
          //canSkip=1; // always skip
          if (iDo>20) {
	    if (!choiceObject) {
	      delete choice.possibleBranch;
	      choice.possibleBranch=NULL;
	    }
            break; // give up anyway
          }
        }
#endif
        if (model->messageHandler()->logLevel()>3&&numberBeforeTrust) 
          dynamicObject->print(1,choice.possibleBranch->value());
        // was if (!canSkip)
        if (newWay)
        numberTest2--;
        if (!canSkip) {
          //#ifndef RANGING 
          if (!doneHotStart) {
            // Mark hot start
            doneHotStart=true;
            assert (auxiliaryInfo->warmStart());
            solver->markHotStart();
            xMark++;
          }
          //#endif
          assert (!couldChooseFirst);
          numberTest--;
          if (!newWay)
          numberTest2--;
          // just do a few
          //if (canSkip)
          //solver->setIntParam(OsiMaxNumIterationHotStart,10); 
          double objectiveChange ;
          double newObjectiveValue=1.0e100;
          int j;
          // status is 0 finished, 1 infeasible and other
          int iStatus;
	  if (0) {
	    CbcDynamicPseudoCostBranchingObject * cbcobj = dynamic_cast<CbcDynamicPseudoCostBranchingObject *> (choice.possibleBranch);
	    if (cbcobj) {
	      CbcSimpleIntegerDynamicPseudoCost * object = cbcobj->object();
	      printf("strong %d ",iDo);
	      object->print(1,0.5);
	    }
	  }
          /*
            Try the down direction first. (Specify the initial branching alternative as
            down with a call to way(-1). Each subsequent call to branch() performs the
            specified branch and advances the branch object state to the next branch
            alternative.)
          */
          choice.possibleBranch->way(-1) ;
#if NEW_UPDATE_OBJECT==0
          decision->saveBranchingObject( choice.possibleBranch);
#endif
          choice.possibleBranch->branch() ;
          solver->solveFromHotStart() ;
          bool needHotStartUpdate=false;
          numberStrongDone++;
          numberStrongIterations += solver->getIterationCount();
          /*
            We now have an estimate of objective degradation that we can use for strong
            branching. If we're over the cutoff, the variable is monotone up.
            If we actually made it to optimality, check for a solution, and if we have
            a good one, call setBestSolution to process it. Note that this may reduce the
            cutoff, so we check again to see if we can declare this variable monotone.
          */
          if (solver->isProvenOptimal())
            iStatus=0; // optimal
          else if (solver->isIterationLimitReached()
                   &&!solver->isDualObjectiveLimitReached())
            iStatus=2; // unknown 
          else
            iStatus=1; // infeasible
          newObjectiveValue = solver->getObjSense()*solver->getObjValue();
          choice.numItersDown = solver->getIterationCount();
          objectiveChange = CoinMax(newObjectiveValue  - objectiveValue_,0.0);
	  // Update branching information if wanted
#if NEW_UPDATE_OBJECT==0
          decision->updateInformation( solver,this);
#elif NEW_UPDATE_OBJECT<2
	  CbcBranchingObject * cbcobj = dynamic_cast<CbcBranchingObject *> (choice.possibleBranch);
	  if (cbcobj) {
	    CbcObject * object = cbcobj->object();
	    CbcObjectUpdateData update = object->createUpdateInformation(solver,this,cbcobj);
	    object->updateInformation(update);
	  } else {
	    decision->updateInformation( solver,this);
	  }
#else
	  CbcBranchingObject * cbcobj = dynamic_cast<CbcBranchingObject *> (choice.possibleBranch);
	  if (cbcobj) {
	    CbcObject * object = cbcobj->object();
	    CbcObjectUpdateData update = object->createUpdateInformation(solver,this,cbcobj);
	    update.objectNumber_ = choice.objectNumber;
	    model->addUpdateInformation(update);
	  } else {
	    decision->updateInformation( solver,this);
	  }
#endif
          if (!iStatus) {
            choice.finishedDown = true ;
            if (newObjectiveValue>=cutoff) {
              objectiveChange = 1.0e100; // say infeasible
              numberStrongInfeasible++;
            } else {
              // See if integer solution
              if (model->feasibleSolution(choice.numIntInfeasDown,
                                          choice.numObjInfeasDown)
                  &&model->problemFeasibility()->feasible(model,-1)>=0) {
                if (auxiliaryInfo->solutionAddsCuts()) {
                  needHotStartUpdate=true;
                  solver->unmarkHotStart();
                }
                model->setBestSolution(CBC_STRONGSOL,
                                       newObjectiveValue,
                                       solver->getColSolution()) ;
                if (needHotStartUpdate) {
                  solver->resolve();
                  newObjectiveValue = solver->getObjSense()*solver->getObjValue();
                  objectiveChange = CoinMax(newObjectiveValue  - objectiveValue_,0.0);
                  model->feasibleSolution(choice.numIntInfeasDown,
                                          choice.numObjInfeasDown);
                }
                model->setLastHeuristic(NULL);
                model->incrementUsed(solver->getColSolution());
                cutoff =model->getCutoff();
                if (newObjectiveValue >= cutoff) {	//  *new* cutoff
                  objectiveChange = 1.0e100 ;
                  numberStrongInfeasible++;
                }
              }
            }
          } else if (iStatus==1) {
            objectiveChange = 1.0e100 ;
            numberStrongInfeasible++;
          } else {
            // Can't say much as we did not finish
            choice.finishedDown = false ;
            numberUnfinished++;
          }
          choice.downMovement = objectiveChange ;
          
          // restore bounds
          for ( j=0;j<numberColumns;j++) {
            if (saveLower[j] != lower[j])
              solver->setColLower(j,saveLower[j]);
            if (saveUpper[j] != upper[j])
              solver->setColUpper(j,saveUpper[j]);
          }
 	  if(needHotStartUpdate) {
            needHotStartUpdate = false;
            solver->resolve();
            //we may again have an integer feasible solution
            int numberIntegerInfeasibilities;
            int numberObjectInfeasibilities;
            if (model->feasibleSolution(
                                        numberIntegerInfeasibilities,
                                        numberObjectInfeasibilities)) {
#ifdef BONMIN
              //In this case node has become integer feasible, let us exit the loop
              std::cout<<"Node has become integer feasible"<<std::endl;
              numberUnsatisfied_ = 0;
              break;
#endif
              double objValue = solver->getObjValue();
              model->setBestSolution(CBC_STRONGSOL,
                                     objValue,
                                     solver->getColSolution()) ;
              solver->resolve();
              cutoff =model->getCutoff();
            }
            solver->markHotStart();
          }
          //printf("Down on %d, status is %d, obj %g its %d cost %g finished %d inf %d infobj %d\n",
          //printf("Down on %d, status is %d, obj %g its %d cost %g finished %d inf %d infobj %d\n",
          //     choice.objectNumber,iStatus,newObjectiveValue,choice.numItersDown,
          //     choice.downMovement,choice.finishedDown,choice.numIntInfeasDown,
          //     choice.numObjInfeasDown);
          
          // repeat the whole exercise, forcing the variable up
#if NEW_UPDATE_OBJECT==0
          decision->saveBranchingObject( choice.possibleBranch);
#endif
          choice.possibleBranch->branch();
          solver->solveFromHotStart() ;
          numberStrongDone++;
          numberStrongIterations += solver->getIterationCount();
          /*
            We now have an estimate of objective degradation that we can use for strong
            branching. If we're over the cutoff, the variable is monotone up.
            If we actually made it to optimality, check for a solution, and if we have
            a good one, call setBestSolution to process it. Note that this may reduce the
            cutoff, so we check again to see if we can declare this variable monotone.
          */
          if (solver->isProvenOptimal())
            iStatus=0; // optimal
          else if (solver->isIterationLimitReached()
                   &&!solver->isDualObjectiveLimitReached())
            iStatus=2; // unknown 
          else
            iStatus=1; // infeasible
          newObjectiveValue = solver->getObjSense()*solver->getObjValue();
          choice.numItersUp = solver->getIterationCount();
          objectiveChange = CoinMax(newObjectiveValue  - objectiveValue_,0.0);
	  // Update branching information if wanted
#if NEW_UPDATE_OBJECT==0
          decision->updateInformation( solver,this);
#elif NEW_UPDATE_OBJECT<2
	  cbcobj = dynamic_cast<CbcBranchingObject *> (choice.possibleBranch);
	  if (cbcobj) {
	    CbcObject * object = cbcobj->object();
	    CbcObjectUpdateData update = object->createUpdateInformation(solver,this,cbcobj);
	    object->updateInformation(update);
	  } else {
	    decision->updateInformation( solver,this);
	  }
#else
	  cbcobj = dynamic_cast<CbcBranchingObject *> (choice.possibleBranch);
	  if (cbcobj) {
	    CbcObject * object = cbcobj->object();
	    CbcObjectUpdateData update = object->createUpdateInformation(solver,this,cbcobj);
	    update.objectNumber_ = choice.objectNumber;
	    model->addUpdateInformation(update);
	  } else {
	    decision->updateInformation( solver,this);
	  }
#endif
          if (!iStatus) {
            choice.finishedUp = true ;
            if (newObjectiveValue>=cutoff) {
              objectiveChange = 1.0e100; // say infeasible
              numberStrongInfeasible++;
            } else {
              // See if integer solution
              if (model->feasibleSolution(choice.numIntInfeasUp,
                                          choice.numObjInfeasUp)
                  &&model->problemFeasibility()->feasible(model,-1)>=0) {
#ifdef BONMIN
                std::cout<<"Node has become integer feasible"<<std::endl;
                numberUnsatisfied_ = 0;
                break;
#endif
                if (auxiliaryInfo->solutionAddsCuts()) {
                  needHotStartUpdate=true;
                  solver->unmarkHotStart();
                }
                model->setBestSolution(CBC_STRONGSOL,
                                       newObjectiveValue,
                                       solver->getColSolution()) ;
                if (needHotStartUpdate) {
                  solver->resolve();
                  newObjectiveValue = solver->getObjSense()*solver->getObjValue();
                  objectiveChange = CoinMax(newObjectiveValue  - objectiveValue_,0.0);
                  model->feasibleSolution(choice.numIntInfeasDown,
                                          choice.numObjInfeasDown);
                }
                model->setLastHeuristic(NULL);
                model->incrementUsed(solver->getColSolution());
                cutoff =model->getCutoff();
                if (newObjectiveValue >= cutoff) {	//  *new* cutoff
                  objectiveChange = 1.0e100 ;
                  numberStrongInfeasible++;
                }
              }
            }
          } else if (iStatus==1) {
            objectiveChange = 1.0e100 ;
            numberStrongInfeasible++;
          } else {
            // Can't say much as we did not finish
            choice.finishedUp = false ;
            numberUnfinished++;
          }
          choice.upMovement = objectiveChange ;
          
          // restore bounds
          for ( j=0;j<numberColumns;j++) {
            if (saveLower[j] != lower[j])
              solver->setColLower(j,saveLower[j]);
            if (saveUpper[j] != upper[j])
              solver->setColUpper(j,saveUpper[j]);
          }
 	  if(needHotStartUpdate) {
            needHotStartUpdate = false;
            solver->resolve();
            //we may again have an integer feasible solution
            int numberIntegerInfeasibilities;
            int numberObjectInfeasibilities;
            if (model->feasibleSolution(
                                        numberIntegerInfeasibilities,
                                        numberObjectInfeasibilities)) {
              double objValue = solver->getObjValue();
              model->setBestSolution(CBC_STRONGSOL,
                                     objValue,
                                     solver->getColSolution()) ;
              solver->resolve();
              cutoff =model->getCutoff();
            }
            solver->markHotStart();
          }
          
          //printf("Up on %d, status is %d, obj %g its %d cost %g finished %d inf %d infobj %d\n",
          //     choice.objectNumber,iStatus,newObjectiveValue,choice.numItersUp,
          //     choice.upMovement,choice.finishedUp,choice.numIntInfeasUp,
          //     choice.numObjInfeasUp);
        }
    
        solver->setIntParam(OsiMaxNumIterationHotStart,saveLimit2); 
        /*
          End of evaluation for this candidate variable. Possibilities are:
          * Both sides below cutoff; this variable is a candidate for branching.
          * Both sides infeasible or above the objective cutoff: no further action
          here. Break from the evaluation loop and assume the node will be purged
          by the caller.
          * One side below cutoff: Install the branch (i.e., fix the variable). Break
          from the evaluation loop and assume the node will be reoptimised by the
          caller.
        */
        // reset
        choice.possibleBranch->resetNumberBranchesLeft();
        if (choice.upMovement<1.0e100) {
          if(choice.downMovement<1.0e100) {
            // In case solution coming in was odd
            choice.upMovement = CoinMax(0.0,choice.upMovement);
            choice.downMovement = CoinMax(0.0,choice.downMovement);
            if (couldChooseFirst)
              printf("candidate %d up %g down %g sort %g\n",iDo,choice.upMovement,choice.downMovement,sort[iDo]);
#if ZERO_ONE==2
            // branch on 0-1 first (temp)
            if (fabs(choice.possibleBranch->value())<1.0) {
              choice.upMovement *= ZERO_FAKE;
              choice.downMovement *= ZERO_FAKE;
            }
#endif
            // feasible - see which best
            if (!canSkip) {
              if (iColumn==-46) {
                printf("sort %g downest %g upest %g ",sort[iDo],downEstimate[iObject],
                     upEstimate[iObject]);
                printf("downMove %g upMove %g value %g current pseudo %g %g\n",
                       choice.downMovement,choice.upMovement,choice.possibleBranch->value(),
                       dynamicObject->downDynamicPseudoCost(),dynamicObject->upDynamicPseudoCost());
              }
              if (model->messageHandler()->logLevel()>3) 
                printf("sort %g downest %g upest %g ",sort[iDo],downEstimate[iObject],
                     upEstimate[iObject]);
              model->messageHandler()->message(CBC_STRONG,*model->messagesPointer())
                << iObject << iColumn
                <<choice.downMovement<<choice.numIntInfeasDown 
                <<choice.upMovement<<choice.numIntInfeasUp 
                <<choice.possibleBranch->value()
                <<CoinMessageEol;
            }
            //if (!stateOfSearch)
            //choice.numIntInfeasDown=99999; // temp fudge
            if (wantMiniTree)
              decision->setBestCriterion(-1.0);
            double bestCriterion = -1.0;
            //double gap = saveUpper[iColumn]-saveLower[iColumn];
            // Give precedence to ones with gap of 1.0 
            //assert(gap>0.0);
            double factor = 1.0; //changeFactor/CoinMin(gap,100.0);
	    int betterWay;
	    {
	      CbcBranchingObject * branchObj =
		dynamic_cast <CbcBranchingObject *>(branch_) ;
	      if (branch_)
		assert (branchObj);
	      betterWay = decision->betterBranch(choice.possibleBranch,
						     branchObj,
						     choice.upMovement*factor,
						     choice.numIntInfeasUp ,
						     choice.downMovement*factor,
						     choice.numIntInfeasDown );
	    }
            if (wantMiniTree) {
              double criterion = decision->getBestCriterion();
              sort[numberMini]=-criterion;
              whichObject[numberMini++]=whichObject[iDo];
              assert (betterWay);
              if (criterion>bestCriterion) 
                bestCriterion=criterion;
              else
                betterWay=0;
            }
            if (iDo>=changeStrategy) {
              // make less likely
              changeStrategy+=numberStrong;
              changeFactor *= 0.9;
            }
            if (betterWay) {
              // C) create branching object
	      if (choiceObject) {
		delete branch_;
		branch_ = choice.possibleBranch->clone();
	      } else {
		delete branch_;
		branch_ = choice.possibleBranch;
		choice.possibleBranch=NULL;
	      }
	      {
		CbcBranchingObject * branchObj =
		  dynamic_cast <CbcBranchingObject *>(branch_) ;
		assert (branchObj);
		//branchObj->way(preferredWay);
		branchObj->way(betterWay);
	      }
              if (couldChooseFirst)
                printf("choosing %d way %d\n",iDo,betterWay);
              bestChoice = choice.objectNumber;
              whichChoice = iDo;
              if (numberStrong<=1) {
                delete ws;
                ws=NULL;
                break;
              }
            } else {
	      if (!choiceObject) {
		delete choice.possibleBranch;
		choice.possibleBranch=NULL;
	      }
              if (iDo>=2*numberStrong) {
                delete ws;
                ws=NULL;
                break;
              }
              if (!dynamicObject||dynamicObject->numberTimesUp()>1) {
                if (iDo-whichChoice>=numberStrong) {
		  if (!choiceObject) {
		    delete choice.possibleBranch;
		    choice.possibleBranch=NULL;
		  }
                  break; // give up
                }
              } else {
                if (iDo-whichChoice>=2*numberStrong) {
                  delete ws;
                  ws=NULL;
		  if (!choiceObject) {
		    delete choice.possibleBranch;
		    choice.possibleBranch=NULL;
		  }
                  break; // give up
                }
              }
            }
          } else {
            // up feasible, down infeasible
            anyAction=-1;
            worstFeasible = CoinMax(worstFeasible,choice.upMovement);
	    model->messageHandler()->message(CBC_STRONG,*model->messagesPointer())
	      << iObject << iColumn
	      <<choice.downMovement<<choice.numIntInfeasDown 
	      <<choice.upMovement<<choice.numIntInfeasUp 
	      <<choice.possibleBranch->value()
	      <<CoinMessageEol;
            //printf("Down infeasible for choice %d sequence %d\n",i,
            // model->object(choice.objectNumber)->columnNumber());
            if (!solveAll) {
              choice.possibleBranch->way(1);
              choice.possibleBranch->branch();
	      if (!choiceObject) {
		delete choice.possibleBranch;
		choice.possibleBranch=NULL;
	      }
              delete ws;
              ws=NULL;
              break;
            } else {
              choice.fix=1;
              fixObject[numberToFix++]=choice;
	      if (!choiceObject) {
		choice.possibleBranch=NULL;
	      } else {
		choiceObject = new CbcDynamicPseudoCostBranchingObject(*choiceObject);
		choice.possibleBranch=choiceObject;
	      }
#define FIXNOW
#ifdef FIXNOW
              double value = ceil(saveSolution[iColumn]);
              saveLower[iColumn]=value;
              solver->setColLower(iColumn,value);
              assert(doneHotStart);
              solver->unmarkHotStart();
	      solver->resolve();
              solver->markHotStart();
	      // may be infeasible (if other way stopped on iterations)
	      if (!solver->isProvenOptimal()) {
		// neither side feasible
		anyAction=-2;
		if (!choiceObject) {
		  delete choice.possibleBranch;
		  choice.possibleBranch=NULL;
		}
		//printf("Both infeasible for choice %d sequence %d\n",i,
		// model->object(choice.objectNumber)->columnNumber());
		delete ws;
		ws=NULL;
		break;
	      }
#endif
            }
          }
        } else {
          if(choice.downMovement<1.0e100) {
            // down feasible, up infeasible
            anyAction=-1;
            worstFeasible = CoinMax(worstFeasible,choice.downMovement);
	    model->messageHandler()->message(CBC_STRONG,*model->messagesPointer())
	      << iObject << iColumn
	      <<choice.downMovement<<choice.numIntInfeasDown 
	      <<choice.upMovement<<choice.numIntInfeasUp 
	      <<choice.possibleBranch->value()
	      <<CoinMessageEol;
            //printf("Up infeasible for choice %d sequence %d\n",i,
            // model->object(choice.objectNumber)->columnNumber());
            if (!solveAll) {
              choice.possibleBranch->way(-1);
              choice.possibleBranch->branch();
	      if (!choiceObject) {
		delete choice.possibleBranch;
		choice.possibleBranch=NULL;
	      }
              delete ws;
              ws=NULL;
              break;
            } else {
              choice.fix=-1;
              fixObject[numberToFix++]=choice;
	      if (!choiceObject) {
		choice.possibleBranch=NULL;
	      } else {
		choiceObject = new CbcDynamicPseudoCostBranchingObject(*choiceObject);
		choice.possibleBranch=choiceObject;
	      }
#ifdef FIXNOW
              double value = floor(saveSolution[iColumn]);
              saveUpper[iColumn]=value;
              solver->setColUpper(iColumn,value);
              assert(doneHotStart);
              solver->unmarkHotStart();
	      solver->resolve();
              solver->markHotStart();
	      // may be infeasible (if other way stopped on iterations)
	      if (!solver->isProvenOptimal()) {
		// neither side feasible
		anyAction=-2;
		if (!choiceObject) {
		  delete choice.possibleBranch;
		  choice.possibleBranch=NULL;
		}
		//printf("Both infeasible for choice %d sequence %d\n",i,
		// model->object(choice.objectNumber)->columnNumber());
		delete ws;
		ws=NULL;
		break;
	      }
#endif
            }
          } else {
            // neither side feasible
            anyAction=-2;
	    if (!choiceObject) {
	      delete choice.possibleBranch;
	      choice.possibleBranch=NULL;
	    }
            //printf("Both infeasible for choice %d sequence %d\n",i,
            // model->object(choice.objectNumber)->columnNumber());
            delete ws;
            ws=NULL;
            break;
          }
        }
        // Check max time
        hitMaxTime = ( CoinCpuTime()-model->getDblParam(CbcModel::CbcStartSeconds) > 
                       model->getDblParam(CbcModel::CbcMaximumSeconds));
        if (hitMaxTime) {
          // make sure rest are fast
          doQuickly=true;
          for ( int jDo=iDo+1;jDo<numberToDo;jDo++) {
            int iObject = whichObject[iDo];
            OsiObject * object = model->modifiableObject(iObject);
            CbcSimpleIntegerDynamicPseudoCost * dynamicObject =
              dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object) ;
            dynamicObject->setNumberBeforeTrust(0);
          }
          numberTest=0;
          distanceToCutoff=-COIN_DBL_MAX;
        }
	if (!choiceObject) {
	  delete choice.possibleBranch;
	}
      }
      double averageChange = model->sumChangeObjective()/((double) model->getNodeCount());
      if (depth_<10||worstFeasible>0.2*averageChange) 
        solveAll=false;
      if (model->messageHandler()->logLevel()>3||false) { 
        if (anyAction==-2) {
          printf("infeasible\n");
        } else if(anyAction==-1) {
          if (!solveAll)
            printf("%d fixed\n",numberToFix);
          else
            printf("%d fixed AND choosing %d iDo %d iChosenWhen %d numberToDo %d\n",numberToFix,bestChoice,
                   iDo,whichChoice,numberToDo);
        } else {
	  int iObject = whichObject[whichChoice];
	  OsiObject * object = model->modifiableObject(iObject);
	  CbcSimpleIntegerDynamicPseudoCost * dynamicObject =
	    dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object) ;
	  int iColumn = dynamicObject->columnNumber();
          printf("choosing %d (column %d) iChosenWhen %d numberToDo %d\n",bestChoice,
                 iColumn,whichChoice,numberToDo);
	}
      }
      if (doneHotStart) {
        // Delete the snapshot
        solver->unmarkHotStart();
        // back to normal
        solver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,NULL) ;
        // restore basis
        solver->setWarmStart(ws);
      }
      solver->setIntParam(OsiMaxNumIterationHotStart,saveLimit);
      // Unless infeasible we will carry on
      // But we could fix anyway
      if (numberToFix&&!hitMaxTime) {
        if (anyAction==-2) {
          // take off
          for (i = 0 ; i < numberToFix ; i++) {
            delete fixObject[i].possibleBranch;
          }
        } else {
          // apply and take off
          for (i = 0 ; i < numberToFix ; i++) {
#ifndef FIXNOW
            fixObject[i].possibleBranch->way(fixObject[i].fix) ;
            fixObject[i].possibleBranch->branch() ;
#endif
            delete fixObject[i].possibleBranch;
          }
          bool feasible=true;
#if ACTION <2
          if (solveAll) {
            // can do quick optimality check
            int easy=2;
            solver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,&easy) ;
            solver->resolve() ;
            solver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,NULL) ;
            feasible = solver->isProvenOptimal();
            if (feasible) {
              anyAction=0;
              numberMini=0;
              memcpy(saveSolution,solver->getColSolution(),numberColumns*sizeof(double));
              model->reserveCurrentSolution(saveSolution);
              memcpy(saveLower,solver->getColLower(),numberColumns*sizeof(double));
              memcpy(saveUpper,solver->getColUpper(),numberColumns*sizeof(double));
              model->setPointers(solver);
              // See if candidate still possible
              if (branch_) {
                const OsiObject * object = model->object(bestChoice);
                int preferredWay;
                double infeasibility = object->infeasibility(&usefulInfo,preferredWay);
                if (!infeasibility) {
                  // take out
                  delete branch_;
                  branch_=NULL;
                } else {
		  CbcBranchingObject * branchObj =
		    dynamic_cast <CbcBranchingObject *>(branch_) ;
		  assert (branchObj);
		  branchObj->way(preferredWay);
                }
              }
            } else {
              anyAction=-2;
              finished=true;
            }
          }
#endif
          // If  fixed then round again
          if (!branch_&&anyAction!=-2) {
            finished=false;
          }
          // If these in then different action
#if ACTION == 1
          if (!anyAction)
            anyAction=-1;
          finished=true;
#endif
        }
      }
      delete ws;
    }
  }
  if (model->messageHandler()->logLevel()>2) 
    printf("%d strong, %d iters, %d pen, %d mark, %d fixed, action %d nnott %d nt %d, %d dq %s ns %d\n",
         numberStrongDone,numberStrongIterations,xPen,xMark,
           numberToFix,anyAction,numberNotTrusted,px[0],px[1],px[2]>0 ? "y" : "n",px[3]);
  // update number of strong iterations etc
  model->incrementStrongInfo(numberStrongDone,numberStrongIterations,
                             anyAction==-2 ? 0:numberToFix,anyAction==-2);
  if (!newWay) {
  if (((model->searchStrategy()+1)%1000)==0) {
    if (solver->messageHandler()->logLevel()>1)
      printf("%d strong, %d iters, %d inf, %d not finished, %d not trusted\n",
             numberStrongDone,numberStrongIterations,numberStrongInfeasible,numberUnfinished,
             numberNotTrusted);
    // decide what to do
    int strategy=1;
    if (numberUnfinished*4>numberStrongDone&&numberStrongInfeasible*10<numberStrongDone) {
      strategy=2;
      if (model->logLevel()>1)
        printf("going to strategy 2\n");
    }
    if (numberNodes)
      strategy=1;  // should only happen after hot start
    if (model->searchStrategy()<999)
      model->setSearchStrategy(strategy);
  }
  }
  //if (numberToFix&&depth_<5)
  //printf("%d fixed by strong at depth %d\n",numberToFix,depth_);
  // Set guessed solution value
  guessedObjectiveValue_ = objectiveValue_+estimatedDegradation;
  
  // Get collection of branches if mini tree wanted
  if (anyAction==0&&numberMini&&numberMini>1) {
    // Sort 
    CoinSort_2(sort,sort+numberMini,whichObject);
    delete branch_;
    branch_=NULL;
    numberMini = CoinMin(numberMini,model->sizeMiniTree());
    anyAction=numberMini;
    branches = new OsiSolverBranch[numberMini];
    for (int iDo=0;iDo<numberMini;iDo++) {
      int iObject = whichObject[iDo];
      OsiObject * object = model->modifiableObject(iObject);
      CbcSimpleInteger * obj =
	dynamic_cast <CbcSimpleInteger *>(object) ;
      OsiSolverBranch * oneBranch;
      if (obj) {
	oneBranch = obj->solverBranch(solver,&usefulInfo);
      } else {
	CbcObject * obj =
	  dynamic_cast <CbcObject *>(object) ;
	assert (obj);
	oneBranch = obj->solverBranch();
      }
      branches[iDo]=*oneBranch;
      delete oneBranch;
    }
  }
/*
  Cleanup, then we're finished
*/
  if (!model->branchingMethod())
    delete decision;

  delete choiceObject;
  delete [] fixObject;
  delete [] sort;
  delete [] whichObject;
  delete [] objectMark;
  delete [] saveLower;
  delete [] saveUpper;
  delete [] upEstimate;
  delete [] downEstimate;
# ifdef COIN_HAS_CLP
  if (osiclp) 
    osiclp->setSpecialOptions(saveClpOptions);
# endif
  // restore solution
  solver->setColSolution(saveSolution);
  model->reserveCurrentSolution(saveSolution);
  delete [] saveSolution;
  model->setStateOfSearch(saveStateOfSearch);
  model->setLogLevel(saveLogLevel);
  return anyAction;
}
int CbcNode::analyze (CbcModel *model, double * results)
{
  int i;
  int numberIterationsAllowed = model->numberAnalyzeIterations();
  OsiSolverInterface * solver = model->solver();
  objectiveValue_ = solver->getObjSense()*solver->getObjValue();
  double cutoff =model->getCutoff();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  const double * dj = solver->getReducedCost();
  int numberObjects = model->numberObjects();
  int numberColumns = model->getNumCols();
  // Initialize arrays
  int numberIntegers = model->numberIntegers();
  int * back = new int[numberColumns];
  const int * integerVariable = model->integerVariable();
  for (i=0;i<numberColumns;i++) 
    back[i]=-1;
  // What results is
  double * newLower = results;
  double * objLower = newLower+numberIntegers;
  double * newUpper = objLower+numberIntegers;
  double * objUpper = newUpper+numberIntegers;
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integerVariable[i];
    back[iColumn]=i;
    newLower[i]=0.0;
    objLower[i]=-COIN_DBL_MAX;
    newUpper[i]=0.0;
    objUpper[i]=-COIN_DBL_MAX;
  }
  double * saveUpper = new double[numberColumns];
  double * saveLower = new double[numberColumns];
  int anyAction=0;
  // Save solution in case heuristics need good solution later
  
  double * saveSolution = new double[numberColumns];
  memcpy(saveSolution,solver->getColSolution(),numberColumns*sizeof(double));
  model->reserveCurrentSolution(saveSolution);
  for (i=0;i<numberColumns;i++) {
    saveLower[i] = lower[i];
    saveUpper[i] = upper[i];
  }
  // Get arrays to sort 
  double * sort = new double[numberObjects];
  int * whichObject = new int[numberObjects];
  int numberToFix=0;
  int numberToDo=0;
  double integerTolerance = 
    model->getDblParam(CbcModel::CbcIntegerTolerance);
  // point to useful information
  OsiBranchingInformation usefulInfo = model->usefulInformation();
  // and modify
  usefulInfo.depth_=depth_;
      
  // compute current state
  int numberObjectInfeasibilities; // just odd ones
  int numberIntegerInfeasibilities;
  model->feasibleSolution(
                          numberIntegerInfeasibilities,
                          numberObjectInfeasibilities);
# ifdef COIN_HAS_CLP
  OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (solver);
  int saveClpOptions=0;
  bool fastIterations = (model->specialOptions()&8)!=0;
  if (osiclp&&fastIterations) {
    // for faster hot start
    saveClpOptions = osiclp->specialOptions();
    osiclp->setSpecialOptions(saveClpOptions|8192);
  }
# else
  bool fastIterations = false ;
# endif
  /*
    Scan for branching objects that indicate infeasibility. Choose candidates
    using priority as the first criteria, then integer infeasibility.
    
    The algorithm is to fill the array with a set of good candidates (by
    infeasibility) with priority bestPriority.  Finding a candidate with
    priority better (less) than bestPriority flushes the choice array. (This
    serves as initialization when the first candidate is found.)
    
  */
  numberToDo=0;
  for (i=0;i<numberObjects;i++) {
    OsiObject * object = model->modifiableObject(i);
    CbcSimpleIntegerDynamicPseudoCost * dynamicObject =
      dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object) ;
    if(!dynamicObject)
      continue;
    int preferredWay;
    double infeasibility = object->infeasibility(&usefulInfo,preferredWay);
    int iColumn = dynamicObject->columnNumber();
    if (saveUpper[iColumn]==saveLower[iColumn])
      continue;
    if (infeasibility)
      sort[numberToDo]=-1.0e10-infeasibility;
    else
      sort[numberToDo]=-fabs(dj[iColumn]);
    whichObject[numberToDo++]=i;
  }
  // Save basis
  CoinWarmStart * ws = solver->getWarmStart();
  int saveLimit;
  solver->getIntParam(OsiMaxNumIterationHotStart,saveLimit);
  int targetIterations = CoinMax(500,numberIterationsAllowed/numberObjects);
  if (saveLimit<targetIterations)
    solver->setIntParam(OsiMaxNumIterationHotStart,targetIterations); 
  // Mark hot start
  solver->markHotStart();
  // Sort 
  CoinSort_2(sort,sort+numberToDo,whichObject);
  //double distanceToCutoff=model->getCutoff()-objectiveValue_;
  double * currentSolution = model->currentSolution();
  double objMin = 1.0e50;
  double objMax = -1.0e50;
  bool needResolve=false;
  int iDo;
  for (iDo=0;iDo<numberToDo;iDo++) {
    CbcStrongInfo choice;
    int iObject = whichObject[iDo];
    OsiObject * object = model->modifiableObject(iObject);
    CbcSimpleIntegerDynamicPseudoCost * dynamicObject =
      dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object) ;
    int iColumn = dynamicObject->columnNumber();
    int preferredWay;
    object->infeasibility(&usefulInfo,preferredWay);
    double value = currentSolution[iColumn];
    double nearest = floor(value+0.5);
    double lowerValue = floor(value);
    bool satisfied=false;
    if (fabs(value-nearest)<=integerTolerance||value<saveLower[iColumn]||value>saveUpper[iColumn]) {
      satisfied=true;
      double newValue;
      if (nearest<saveUpper[iColumn]) {
        newValue = nearest + 1.0001*integerTolerance;
        lowerValue = nearest;
      } else {
        newValue = nearest - 1.0001*integerTolerance;
        lowerValue = nearest-1;
      }
      currentSolution[iColumn]=newValue;
    }
    double upperValue = lowerValue+1.0;
    CbcSimpleInteger * obj =
      dynamic_cast <CbcSimpleInteger *>(object) ;
    if (obj) {
      choice.possibleBranch=obj->createBranch(solver,&usefulInfo,preferredWay);
    } else {
      CbcObject * obj =
	dynamic_cast <CbcObject *>(object) ;
      assert (obj);
      choice.possibleBranch=obj->createBranch(preferredWay);
    }
    currentSolution[iColumn]=value;
    // Save which object it was
    choice.objectNumber=iObject;
    choice.numIntInfeasUp = numberUnsatisfied_;
    choice.numIntInfeasDown = numberUnsatisfied_;
    choice.downMovement = 0.0;
    choice.upMovement = 0.0;
    choice.numItersDown = 0;
    choice.numItersUp = 0;
    choice.fix=0; // say not fixed
    double objectiveChange ;
    double newObjectiveValue=1.0e100;
    int j;
    // status is 0 finished, 1 infeasible and other
    int iStatus;
    /*
      Try the down direction first. (Specify the initial branching alternative as
      down with a call to way(-1). Each subsequent call to branch() performs the
      specified branch and advances the branch object state to the next branch
      alternative.)
    */
    choice.possibleBranch->way(-1) ;
    choice.possibleBranch->branch() ;
    if (fabs(value-lowerValue)>integerTolerance) {
      solver->solveFromHotStart() ;
      /*
        We now have an estimate of objective degradation that we can use for strong
        branching. If we're over the cutoff, the variable is monotone up.
        If we actually made it to optimality, check for a solution, and if we have
        a good one, call setBestSolution to process it. Note that this may reduce the
        cutoff, so we check again to see if we can declare this variable monotone.
      */
      if (solver->isProvenOptimal())
        iStatus=0; // optimal
      else if (solver->isIterationLimitReached()
               &&!solver->isDualObjectiveLimitReached())
        iStatus=2; // unknown 
      else
        iStatus=1; // infeasible
      newObjectiveValue = solver->getObjSense()*solver->getObjValue();
      choice.numItersDown = solver->getIterationCount();
      numberIterationsAllowed -= choice.numItersDown;
      objectiveChange = newObjectiveValue  - objectiveValue_;
      if (!iStatus) {
        choice.finishedDown = true ;
        if (newObjectiveValue>=cutoff) {
          objectiveChange = 1.0e100; // say infeasible
        } else {
          // See if integer solution
          if (model->feasibleSolution(choice.numIntInfeasDown,
                                      choice.numObjInfeasDown)
              &&model->problemFeasibility()->feasible(model,-1)>=0) {
            model->setBestSolution(CBC_STRONGSOL,
                                   newObjectiveValue,
                                   solver->getColSolution()) ;
            model->setLastHeuristic(NULL);
            model->incrementUsed(solver->getColSolution());
            cutoff =model->getCutoff();
            if (newObjectiveValue >= cutoff)	//  *new* cutoff
              objectiveChange = 1.0e100 ;
          }
        }
      } else if (iStatus==1) {
        objectiveChange = 1.0e100 ;
      } else {
        // Can't say much as we did not finish
        choice.finishedDown = false ;
      }
      choice.downMovement = objectiveChange ;
    }
    // restore bounds
    for ( j=0;j<numberColumns;j++) {
      if (saveLower[j] != lower[j])
        solver->setColLower(j,saveLower[j]);
      if (saveUpper[j] != upper[j])
        solver->setColUpper(j,saveUpper[j]);
    }
    // repeat the whole exercise, forcing the variable up
    choice.possibleBranch->branch();
    if (fabs(value-upperValue)>integerTolerance) {
      solver->solveFromHotStart() ;
      /*
        We now have an estimate of objective degradation that we can use for strong
        branching. If we're over the cutoff, the variable is monotone up.
        If we actually made it to optimality, check for a solution, and if we have
        a good one, call setBestSolution to process it. Note that this may reduce the
        cutoff, so we check again to see if we can declare this variable monotone.
      */
      if (solver->isProvenOptimal())
        iStatus=0; // optimal
      else if (solver->isIterationLimitReached()
               &&!solver->isDualObjectiveLimitReached())
        iStatus=2; // unknown 
      else
        iStatus=1; // infeasible
      newObjectiveValue = solver->getObjSense()*solver->getObjValue();
      choice.numItersUp = solver->getIterationCount();
      numberIterationsAllowed -= choice.numItersUp;
      objectiveChange = newObjectiveValue  - objectiveValue_;
      if (!iStatus) {
        choice.finishedUp = true ;
        if (newObjectiveValue>=cutoff) {
          objectiveChange = 1.0e100; // say infeasible
        } else {
          // See if integer solution
          if (model->feasibleSolution(choice.numIntInfeasUp,
                                      choice.numObjInfeasUp)
              &&model->problemFeasibility()->feasible(model,-1)>=0) {
            model->setBestSolution(CBC_STRONGSOL,
                                   newObjectiveValue,
                                   solver->getColSolution()) ;
            model->setLastHeuristic(NULL);
            model->incrementUsed(solver->getColSolution());
            cutoff =model->getCutoff();
            if (newObjectiveValue >= cutoff)	//  *new* cutoff
              objectiveChange = 1.0e100 ;
          }
        }
      } else if (iStatus==1) {
        objectiveChange = 1.0e100 ;
      } else {
        // Can't say much as we did not finish
        choice.finishedUp = false ;
      }
      choice.upMovement = objectiveChange ;
      
      // restore bounds
      for ( j=0;j<numberColumns;j++) {
        if (saveLower[j] != lower[j])
          solver->setColLower(j,saveLower[j]);
        if (saveUpper[j] != upper[j])
          solver->setColUpper(j,saveUpper[j]);
      }
    }
    // If objective goes above certain amount we can set bound
    int jInt = back[iColumn];
    newLower[jInt]=upperValue;
    if (choice.finishedDown)
      objLower[jInt]=choice.downMovement+objectiveValue_;
    else
      objLower[jInt]=objectiveValue_;
    newUpper[jInt]=lowerValue;
    if (choice.finishedUp)
      objUpper[jInt]=choice.upMovement+objectiveValue_;
    else
      objUpper[jInt]=objectiveValue_;
    objMin = CoinMin(CoinMin(objLower[jInt],objUpper[jInt]),objMin);
    /*
      End of evaluation for this candidate variable. Possibilities are:
      * Both sides below cutoff; this variable is a candidate for branching.
      * Both sides infeasible or above the objective cutoff: no further action
      here. Break from the evaluation loop and assume the node will be purged
      by the caller.
      * One side below cutoff: Install the branch (i.e., fix the variable). Break
      from the evaluation loop and assume the node will be reoptimised by the
      caller.
    */
    if (choice.upMovement<1.0e100) {
      if(choice.downMovement<1.0e100) {
        objMax = CoinMax(CoinMax(objLower[jInt],objUpper[jInt]),objMax);
        // In case solution coming in was odd
        choice.upMovement = CoinMax(0.0,choice.upMovement);
        choice.downMovement = CoinMax(0.0,choice.downMovement);
        // feasible -
        model->messageHandler()->message(CBC_STRONG,*model->messagesPointer())
          << iObject << iColumn
          <<choice.downMovement<<choice.numIntInfeasDown 
          <<choice.upMovement<<choice.numIntInfeasUp 
          <<value
          <<CoinMessageEol;
      } else {
        // up feasible, down infeasible
        anyAction=-1;
        if (!satisfied)
          needResolve=true;
        choice.fix=1;
        numberToFix++;
        saveLower[iColumn]=upperValue;
        solver->setColLower(iColumn,upperValue);
      }
    } else {
      if(choice.downMovement<1.0e100) {
        // down feasible, up infeasible
        anyAction=-1;
        if (!satisfied)
          needResolve=true;
        choice.fix=-1;
        numberToFix++;
        saveUpper[iColumn]=lowerValue;
        solver->setColUpper(iColumn,lowerValue);
      } else {
        // neither side feasible
        anyAction=-2;
        printf("Both infeasible for choice %d sequence %d\n",i,
               model->object(choice.objectNumber)->columnNumber());
        delete ws;
        ws=NULL;
        //solver->writeMps("bad");
        numberToFix=-1;
        delete choice.possibleBranch;
        choice.possibleBranch=NULL;
        break;
      }
    }
    delete choice.possibleBranch;
    if (numberIterationsAllowed<=0)
      break;
    //printf("obj %d, col %d, down %g up %g value %g\n",iObject,iColumn,
    //     choice.downMovement,choice.upMovement,value);
  }
  printf("Best possible solution %g, can fix more if solution of %g found - looked at %d variables in %d iterations\n",
         objMin,objMax,iDo,model->numberAnalyzeIterations()-numberIterationsAllowed);
  model->setNumberAnalyzeIterations(numberIterationsAllowed);
  // Delete the snapshot
  solver->unmarkHotStart();
  // back to normal
  solver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,NULL) ;
  solver->setIntParam(OsiMaxNumIterationHotStart,saveLimit);
  // restore basis
  solver->setWarmStart(ws);
  delete ws;
    
  delete [] sort;
  delete [] whichObject;
  delete [] saveLower;
  delete [] saveUpper;
  delete [] back;
  // restore solution
  solver->setColSolution(saveSolution);
# ifdef COIN_HAS_CLP
  if (osiclp) 
    osiclp->setSpecialOptions(saveClpOptions);
# endif
  model->reserveCurrentSolution(saveSolution);
  delete [] saveSolution;
  if (needResolve)
    solver->resolve();
  return numberToFix;
}


CbcNode::CbcNode(const CbcNode & rhs) 
{  
#ifdef CHECK_NODE
  printf("CbcNode %x Constructor from rhs %x\n",this,&rhs);
#endif
  if (rhs.nodeInfo_)
    nodeInfo_ = rhs.nodeInfo_->clone();
  else
    nodeInfo_=NULL;
  objectiveValue_=rhs.objectiveValue_;
  guessedObjectiveValue_ = rhs.guessedObjectiveValue_;
  sumInfeasibilities_ = rhs.sumInfeasibilities_;
  if (rhs.branch_)
    branch_=rhs.branch_->clone();
  else
    branch_=NULL;
  depth_ = rhs.depth_;
  numberUnsatisfied_ = rhs.numberUnsatisfied_;
  nodeNumber_ = rhs.nodeNumber_;
  state_ = rhs.state_;
  if (nodeInfo_)
    assert ((state_&2)!=0);
  else
    assert ((state_&2)==0);
}

CbcNode &
CbcNode::operator=(const CbcNode & rhs)
{
  if (this != &rhs) {
    delete nodeInfo_;
    if (rhs.nodeInfo_)
      nodeInfo_ = rhs.nodeInfo_->clone();
    else
      nodeInfo_ = NULL;
    objectiveValue_=rhs.objectiveValue_;
    guessedObjectiveValue_ = rhs.guessedObjectiveValue_;
    sumInfeasibilities_ = rhs.sumInfeasibilities_;
    if (rhs.branch_)
      branch_=rhs.branch_->clone();
    else
      branch_=NULL,
    depth_ = rhs.depth_;
    numberUnsatisfied_ = rhs.numberUnsatisfied_;
    nodeNumber_ = rhs.nodeNumber_;
    state_ = rhs.state_;
    if (nodeInfo_)
      assert ((state_&2)!=0);
    else
      assert ((state_&2)==0);
  }
  return *this;
}
CbcNode::~CbcNode ()
{
#ifdef CHECK_NODE
  if (nodeInfo_) {
    printf("CbcNode %x Destructor nodeInfo %x (%d)\n",
	 this,nodeInfo_,nodeInfo_->numberPointingToThis());
    //assert(nodeInfo_->numberPointingToThis()>=0);
  } else {
    printf("CbcNode %x Destructor nodeInfo %x (?)\n",
	 this,nodeInfo_);
  }
#endif
  if (nodeInfo_) {
    // was if (nodeInfo_&&(state_&2)!=0) {
    nodeInfo_->nullOwner();
    int numberToDelete=nodeInfo_->numberBranchesLeft();
    //    CbcNodeInfo * parent = nodeInfo_->parent();
    //assert (nodeInfo_->numberPointingToThis()>0);
    if (nodeInfo_->decrement(numberToDelete)==0||(state_&2)==0) {
      if ((state_&2)==0) 
	nodeInfo_->nullParent();
      delete nodeInfo_;
    } else {
      //printf("node %x nodeinfo %x parent %x\n",this,nodeInfo_,nodeInfo_->parent());
      // anyway decrement parent
      //if (parent)
      ///parent->decrement(1);
    }
  }
  delete branch_;
}
// Decrement  active cut counts 
void 
CbcNode::decrementCuts(int change)
{
  if (nodeInfo_)
    assert ((state_&2)!=0);
  else
    assert ((state_&2)==0);
  if(nodeInfo_) {
    nodeInfo_->decrementCuts(change);
  }
}
void 
CbcNode::decrementParentCuts(CbcModel * model, int change)
{
  if (nodeInfo_)
    assert ((state_&2)!=0);
  else
    assert ((state_&2)==0);
  if(nodeInfo_) {
    nodeInfo_->decrementParentCuts(model, change);
  }
}

/*
  Initialize reference counts (numberPointingToThis, numberBranchesLeft_)
  in the attached nodeInfo_.
*/
void
CbcNode::initializeInfo()
{
  assert(nodeInfo_ && branch_) ;
  nodeInfo_->initializeInfo(branch_->numberBranches());
  assert ((state_&2)!=0);
}
// Nulls out node info
void 
CbcNode::nullNodeInfo()
{
  nodeInfo_=NULL;
  // say not active
  state_ &= ~2;
}

int
CbcNode::branch(OsiSolverInterface * solver)
{
  double changeInGuessed;
  if (!solver)
    changeInGuessed=branch_->branch();
  else
    changeInGuessed=branch_->branch(solver);
  guessedObjectiveValue_+= changeInGuessed;
  //#define PRINTIT
#ifdef PRINTIT
  int numberLeft = nodeInfo_->numberBranchesLeft();
  CbcNodeInfo * parent = nodeInfo_->parent();
  int parentNodeNumber = -1;
  //CbcBranchingObject * object1 = branch_->object_;
  //OsiObject * object = object1->
  //int sequence = object->columnNumber);
  int id=-1;
  double value=0.0;
  if (branch_) {
    id = branch_->variable();
    value = branch_->value();
  }
  printf("id %d value %g objvalue %g\n",id,value,objectiveValue_);
  if (parent)
    parentNodeNumber = parent->nodeNumber();
  printf("Node number %d, %s, way %d, depth %d, parent node number %d\n",
	 nodeInfo_->nodeNumber(),(numberLeft==2) ? "leftBranch" : "rightBranch",
	 way(),depth_,parentNodeNumber);
#endif
  return nodeInfo_->branchedOn();
}
/* Active arm of the attached OsiBranchingObject.
  
   In the simplest instance, coded -1 for the down arm of the branch, +1 for
   the up arm. But see OsiBranchingObject::way() 
     Use nodeInfo--.numberBranchesLeft_ to see how active
*/
int 
CbcNode::way() const
{
  if (branch_) {
    CbcBranchingObject * obj =
      dynamic_cast <CbcBranchingObject *>(branch_) ;
    if (obj) {
      return obj->way();
    } else {
      OsiTwoWayBranchingObject * obj2 =
      dynamic_cast <OsiTwoWayBranchingObject *>(branch_) ;
      assert (obj2);
      return obj2->way();
    }
  } else {
    return 0;
  }
}
/* Create a branching object for the node

    The routine scans the object list of the model and selects a set of
    unsatisfied objects as candidates for branching. The candidates are
    evaluated, and an appropriate branch object is installed.

    The numberPassesLeft is decremented to stop fixing one variable each time
    and going on and on (e.g. for stock cutting, air crew scheduling)

    If evaluation determines that an object is monotone or infeasible,
    the routine returns immediately. In the case of a monotone object,
    the branch object has already been called to modify the model.

    Return value:
    <ul>
      <li>  0: A branching object has been installed
      <li> -1: A monotone object was discovered
      <li> -2: An infeasible object was discovered
    </ul>
    Branch state:
    <ul>
      <li> -1: start
      <li> -1: A monotone object was discovered
      <li> -2: An infeasible object was discovered
    </ul>
*/
int 
CbcNode::chooseOsiBranch (CbcModel * model,
			  CbcNode * lastNode,
			  OsiBranchingInformation * usefulInfo,
			  int branchState)
{
  int returnStatus=0;
  if (lastNode)
    depth_ = lastNode->depth_+1;
  else
    depth_ = 0;
  OsiSolverInterface * solver = model->solver();
  objectiveValue_ = solver->getObjValue()*solver->getObjSense();
  usefulInfo->objectiveValue_ = objectiveValue_;
  usefulInfo->depth_ = depth_;
  const double * saveInfoSol = usefulInfo->solution_;
  double * saveSolution = new double[solver->getNumCols()];
  memcpy(saveSolution,solver->getColSolution(),solver->getNumCols()*sizeof(double));
  usefulInfo->solution_ = saveSolution;
  OsiChooseVariable * choose = model->branchingMethod()->chooseMethod();
  int numberUnsatisfied=-1;
  if (branchState<0) {
    // initialize
    // initialize sum of "infeasibilities"
    sumInfeasibilities_ = 0.0;
    numberUnsatisfied = choose->setupList(usefulInfo,true);
    numberUnsatisfied_ = numberUnsatisfied;
    branchState=0;
    if (numberUnsatisfied_<0) {
      // infeasible
      delete [] saveSolution;
      return -2;
    }
  }
  // unset best
  int best=-1;
  choose->setBestObjectIndex(-1);
  if (numberUnsatisfied) {
    if (branchState>0||!choose->numberOnList()) {
      // we need to return at once - don't do strong branching or anything
      if (choose->numberOnList()||!choose->numberStrong()) {
	best = choose->candidates()[0];
	choose->setBestObjectIndex(best);
      } else {
	// nothing on list - need to try again - keep any solution
	numberUnsatisfied = choose->setupList(usefulInfo, false);
	numberUnsatisfied_ = numberUnsatisfied;
	if (numberUnsatisfied) {
	  best = choose->candidates()[0];
	  choose->setBestObjectIndex(best);
	}
      }
    } else {
      // carry on with strong branching or whatever
      int returnCode = choose->chooseVariable(solver, usefulInfo,true);
      // update number of strong iterations etc
      model->incrementStrongInfo(choose->numberStrongDone(),choose->numberStrongIterations(),
                                 returnCode==-1 ? 0:choose->numberStrongFixed(),returnCode==-1);
      if (returnCode>1) {
	// has fixed some
	returnStatus=-1;
      } else if (returnCode==-1) {
	// infeasible
	returnStatus=-2;
      } else if (returnCode==0) {
	// normal
	returnStatus=0;
	numberUnsatisfied=1;
      } else {
	// ones on list satisfied - double check
	numberUnsatisfied = choose->setupList(usefulInfo, false);
	numberUnsatisfied_ = numberUnsatisfied;
	if (numberUnsatisfied) {
	  best = choose->candidates()[0];
	  choose->setBestObjectIndex(best);
	}
      }
    }
  } 
  delete branch_;
  branch_ = NULL;
  guessedObjectiveValue_ = COIN_DBL_MAX;//objectiveValue_; // for now
  if (!returnStatus) {
    if (numberUnsatisfied) {
      // create branching object
      const OsiObject * obj = model->solver()->object(choose->bestObjectIndex());
      //const OsiSolverInterface * solver = usefulInfo->solver_;
      branch_ = obj->createBranch(model->solver(),usefulInfo,obj->whichWay());
    }
  }
  usefulInfo->solution_=saveInfoSol;
  delete [] saveSolution;
  // may have got solution
  if (choose->goodSolution()
      &&model->problemFeasibility()->feasible(model,-1)>=0) {
    // yes
    double objValue = choose->goodObjectiveValue();
    model->setBestSolution(CBC_STRONGSOL,
                                     objValue,
                                     choose->goodSolution()) ;
    model->setLastHeuristic(NULL);
    model->incrementUsed(choose->goodSolution());
    choose->clearGoodSolution();
  }
  return returnStatus;
}
