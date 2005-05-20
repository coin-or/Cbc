// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <string>
//#define CBC_DEBUG 1
//#define CHECK_CUT_COUNTS
//#define CHECK_NODE
#include <cassert>
#include <cfloat>
#define CUTS
#include "OsiSolverInterface.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinTime.hpp"
#include "CbcModel.hpp"
#include "CbcNode.hpp"
#include "CbcBranchActual.hpp"
#include "CbcBranchDynamic.hpp"
#include "OsiRowCut.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiCuts.hpp"
#include "CbcCountRowCut.hpp"
#include "CbcMessage.hpp"
#include "OsiClpSolverInterface.hpp"
#include "ClpSimplexOther.hpp"
using namespace std;
#include "CglCutGenerator.hpp"
// Default Constructor 
CbcNodeInfo::CbcNodeInfo ()
  :
  numberPointingToThis_(0),
  parent_(NULL),
  owner_(NULL),
  numberCuts_(0),
  cuts_(NULL),
  numberRows_(0),
  numberBranchesLeft_(0)
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
  cuts_(NULL),
  numberRows_(0),
  numberBranchesLeft_(2)
{
#ifdef CHECK_NODE
  printf("CbcNodeInfo %x Constructor from parent %x\n",this,parent_);
#endif
  if (parent_) {
    numberRows_ = parent_->numberRows_+parent_->numberCuts_;
  }
}
// Copy Constructor 
CbcNodeInfo::CbcNodeInfo (const CbcNodeInfo & rhs)
  :
  numberPointingToThis_(rhs.numberPointingToThis_),
  parent_(rhs.parent_),
  owner_(rhs.owner_),
  numberCuts_(rhs.numberCuts_),
  cuts_(NULL),
  numberRows_(rhs.numberRows_),
  numberBranchesLeft_(rhs.numberBranchesLeft_)
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
  cuts_(NULL),
  numberRows_(0),
  numberBranchesLeft_(2)
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
  // But they may be some left (max nodes?)
  for (int i=0;i<numberCuts_;i++) 
    delete cuts_[i];

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
	delete cuts_[i];
	cuts_[i]=NULL;
      }
    }
  }
}
void
CbcNodeInfo::decrementParentCuts(int change)
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
    CoinWarmStartBasis dummy;
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
	    delete thisInfo->cuts_[i];
	    thisInfo->cuts_[i]=NULL;
	  }
	}
      }
      thisInfo = thisInfo->parent_;
    }
  }
}

void
CbcNodeInfo::incrementParentCuts(int change)
{
  if (parent_) {
    int i;
    // Get over-estimate of space needed for basis
    CoinWarmStartBasis dummy;
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
    if (!number)
      delete cuts_[j];
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
    if (!number)
      delete cuts_[iCut];
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
  int i;
  int numberColumns = model->getNumCols();
  for (i=0;i<numberColumns;i++) {
    solver->setColBounds(i,lower_[i],upper_[i]);
  }
  // move basis - but make sure size stays
  int numberRows = basis->getNumArtificial();
  delete basis ;
  basis = dynamic_cast<CoinWarmStartBasis *>(basis_->clone()) ;
  basis->resize(numberRows,numberColumns);
  for (i=0;i<numberCuts_;i++) 
    addCuts[currentNumberCuts+i]= cuts_[i];
  currentNumberCuts += numberCuts_;
  assert(!parent_);
  return ;
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
  variables_ = new int [numberChangedBounds_];
  newBounds_ = new double [numberChangedBounds_];

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
  variables_ = new int [numberChangedBounds_];
  newBounds_ = new double [numberChangedBounds_];

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
  delete [] variables_;
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

  basis->applyDiff(basisDiff_) ;

  // branch - do bounds
  int i;
  for (i=0;i<numberChangedBounds_;i++) {
    int variable = variables_[i];
    if ((variable&0x80000000)==0) {
      // lower bound changing
      solver->setColLower(variable,newBounds_[i]);
    } else {
      // upper bound changing
      solver->setColUpper(variable&0x7fffffff,newBounds_[i]);
    }
  }
  for (i=0;i<numberCuts_;i++) 
    addCuts[currentNumberCuts+i]= cuts_[i];
  currentNumberCuts += numberCuts_;
  return ;
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
  branch_(NULL),
  depth_(-1),
  numberUnsatisfied_(0),
  nodeNumber_(-1)
{
#ifdef CHECK_NODE
  printf("CbcNode %x Constructor\n",this);
#endif
}

CbcNode::CbcNode(CbcModel * model,
		 CbcNode * lastNode) :
  nodeInfo_(NULL),
  objectiveValue_(1.0e100),
  guessedObjectiveValue_(1.0e100),
  branch_(NULL),
  depth_(-1),
  numberUnsatisfied_(0),
  nodeNumber_(model->getNodeCount())
{
#ifdef CHECK_NODE
  printf("CbcNode %x Constructor from model\n",this);
#endif
  OsiSolverInterface * solver = model->solver();
  objectiveValue_ = solver->getObjSense()*solver->getObjValue();

  if (lastNode)
    lastNode->nodeInfo_->increment();
}


void
CbcNode::createInfo (CbcModel *model,
		     CbcNode *lastNode,
		     const CoinWarmStartBasis *lastws,
		     const double *lastLower, const double *lastUpper,
		     int numberOldActiveCuts,int numberNewCuts)
{ OsiSolverInterface * solver = model->solver();
/*
  The root --- no parent. Create full basis and bounds information.
*/
  if (!lastNode)
  { nodeInfo_=new CbcFullNodeInfo(model,solver->getNumRows()); }
/*
  Not the root. Create an edit from the parent's basis & bound information.
  This is not quite as straightforward as it seems. We need to reintroduce
  cuts we may have dropped out of the basis, in the correct position, because
  this whole process is strictly positional. Start by grabbing the current
  basis.
*/
  else
  { const CoinWarmStartBasis* ws =
      dynamic_cast<const CoinWarmStartBasis*>(solver->getWarmStart());
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
    expanded->resize(iFull,numberColumns);
#ifdef FULL_DEBUG
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
#ifdef FULL_DEBUG
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
  Yes - cuts will break
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
    if (lastNode->branchingObject()->boundBranch())
      assert (numberChangedBounds);
/*
  Hand the lot over to the CbcPartialNodeInfo constructor, then clean up and
  return.
*/
    nodeInfo_ =
      new CbcPartialNodeInfo(lastNode->nodeInfo_,this,numberChangedBounds,
			     variables,boundChanges,basisDiff) ;
    delete basisDiff ;
    delete [] boundChanges;
    delete [] variables;
    delete expanded ;
    delete ws;
  }
}

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
  double objectiveValue = solver->getObjSense()*saveObjectiveValue;
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  int anyAction=0;
  double integerTolerance = 
    model->getDblParam(CbcModel::CbcIntegerTolerance);
  int i;
  bool beforeSolution = model->getSolutionCount()==0;
  int numberStrong=model->numberStrong();
  int saveNumberStrong=numberStrong;
  int numberObjects = model->numberObjects();
  bool checkFeasibility = numberObjects>model->numberIntegers();
  int maximumStrong = CoinMax(CoinMin(model->numberStrong(),numberObjects),1);
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
  if (!decision)
    decision = new CbcBranchDefaultDecision();

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
    int numberIntegerInfeasibilities=0; // without odd ones
    
    // We may go round this loop twice (only if we think we have solution)
    for (int iPass=0;iPass<2;iPass++) {
      
      // compute current state
      int numberObjectInfeasibilities; // just odd ones
      model->feasibleSolution(
                              numberIntegerInfeasibilities,
                              numberObjectInfeasibilities);
      // If forcePriority > 0 then we want best solution
      const double * bestSolution = NULL;
      int hotstartStrategy=model->getHotstartStrategy();
      if (hotstartStrategy>0) {
        bestSolution = model->bestSolution();
      }
      
      // Some objects may compute an estimate of best solution from here
      estimatedDegradation=0.0; 
      numberUnsatisfied_ = 0;
      int bestPriority=INT_MAX;
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
      for (i=0;i<numberObjects;i++) {
        CbcObject * object = model->modifiableObject(i);
        int preferredWay;
        double infeasibility = object->infeasibility(preferredWay);
        int priorityLevel = object->priority();
        if (bestSolution) {
          // we are doing hot start
          const CbcSimpleInteger * thisOne = dynamic_cast <const CbcSimpleInteger *> (object);
          if (thisOne) {
            int iColumn = thisOne->modelSequence();
            if (saveUpper[iColumn]>saveLower[iColumn]) {
              double value = saveSolution[iColumn];
              double targetValue = bestSolution[iColumn];
              //double originalLower = thisOne->originalLower();
              //double originalUpper = thisOne->originalUpper();
              // switch off if not possible
              if (targetValue>=saveLower[iColumn]&&targetValue<=saveUpper[iColumn]) {
                /* priority outranks rest always if hotstartStrategy >1
                   otherwise can be downgraded if at correct level.
                   Infeasibility may be increased by targetValue to choose 1.0 values first.
                */
                if (fabs(value-targetValue)>integerTolerance) {
                  if (value>targetValue) {
                    infeasibility += value;
                    preferredWay=-1;
                  } else {
                    infeasibility += targetValue;
                    preferredWay=1;
                  }
                } else if (hotstartStrategy>1) {
                  if (targetValue==saveLower[iColumn]) {
                    infeasibility += integerTolerance+1.0e-12;
                    preferredWay=-1;
                  } else if (targetValue==saveUpper[iColumn]) {
                    infeasibility += integerTolerance+1.0e-12;
                    preferredWay=1;
                  } else {
                    infeasibility += integerTolerance+1.0e-12;
                    preferredWay=1;
                  }
                } else {
                  priorityLevel += 10000000;
                }
              } else {
                // switch off if not possible
                bestSolution=NULL;
                model->setHotstartStrategy(0);
              }
            }
          }
        }
        if (infeasibility) {
          // Increase estimated degradation to solution
          estimatedDegradation += CoinMin(object->upEstimate(),object->downEstimate());
          numberUnsatisfied_++;
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
            choice[iSmallest].possibleBranch=object->createBranch(preferredWay);
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
        if (roundAgain) {
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
            model->messageHandler()->message(CBC_ITERATE_STRONG,model->messages())
              << clp->numberIterations()
              <<CoinMessageEol;
        }
        clp->setMaximumIterations(saveLimit);
#ifdef CRUNCH
        int * backColumn = whichColumn+numberColumns;
#endif
        for (i=0;i<numberStrong;i++) {
          int iObject = choice[i].objectNumber;
          const CbcObject * object = model->object(iObject);
          const CbcSimpleInteger * simple = dynamic_cast <const CbcSimpleInteger *> (object);
          int iSequence = simple->modelSequence();
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
            const CbcObject * object = model->object(iObject);
            const CbcSimpleInteger * simple = dynamic_cast <const CbcSimpleInteger *> (object);
            int iSequence = simple->modelSequence();
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
          objectiveChange = newObjectiveValue-objectiveValue ;
        } else {
          iStatus = outputStuff[2*i];
          choice[i].numItersDown = outputStuff[2*numberStrong+2*i];
          newObjectiveValue = objectiveValue+newUpper[i];
          solver->setColSolution(outputSolution[2*i]);
        }
        objectiveChange = newObjectiveValue  - objectiveValue;
        if (!iStatus) {
          choice[i].finishedDown = true ;
          if (newObjectiveValue>=model->getCutoff()) {
            objectiveChange = 1.0e100; // say infeasible
          } else {
            // See if integer solution
            if (model->feasibleSolution(choice[i].numIntInfeasDown,
                                        choice[i].numObjInfeasDown))
              { model->setBestSolution(CBC_STRONGSOL,
                                       newObjectiveValue,
                                       solver->getColSolution()) ;
              model->incrementUsed(solver->getColSolution());
              if (newObjectiveValue >= model->getCutoff())	//  *new* cutoff
                objectiveChange = 1.0e100 ;
              }
          }
        } else if (iStatus==1) {
          objectiveChange = 1.0e100 ;
        } else {
          // Can't say much as we did not finish
          choice[i].finishedDown = false ;
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
          choice[i].possibleBranch->branch();
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
          objectiveChange = newObjectiveValue-objectiveValue ;
        } else {
          iStatus = outputStuff[2*i+1];
          choice[i].numItersUp = outputStuff[2*numberStrong+2*i+1];
          newObjectiveValue = objectiveValue+newLower[i];
          solver->setColSolution(outputSolution[2*i+1]);
        }
        objectiveChange = newObjectiveValue  - objectiveValue;
        if (!iStatus) {
          choice[i].finishedUp = true ;
          if (newObjectiveValue>=model->getCutoff()) {
            objectiveChange = 1.0e100; // say infeasible
          } else {
            // See if integer solution
            if (model->feasibleSolution(choice[i].numIntInfeasUp,
                                        choice[i].numObjInfeasUp))
              { model->setBestSolution(CBC_STRONGSOL,
                                       newObjectiveValue,
                                       solver->getColSolution()) ;
              model->incrementUsed(solver->getColSolution());
              if (newObjectiveValue >= model->getCutoff())	//  *new* cutoff
                objectiveChange = 1.0e100 ;
              }
          }
        } else if (iStatus==1) {
          objectiveChange = 1.0e100 ;
        } else {
          // Can't say much as we did not finish
          choice[i].finishedUp = false ;
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
            const CbcObject * object = model->object(thisChoice.objectNumber);
            int preferredWay;
            double infeasibility = object->infeasibility(preferredWay);
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
      
      /*
        anyAction >= 0 indicates that strong branching didn't produce any monotone
        variables. Sift through the candidates for the best one.
        
        QUERY: Setting numberNodes looks to be a distributed noop. numberNodes is
        local to this code block. Perhaps should be numberNodes_ from model?
        Unclear what this calculation is doing.
      */
      if (anyAction>=0) {
        
        int numberNodes = model->getNodeCount();
        // get average cost per iteration and assume stopped ones
        // would stop after 50% more iterations at average cost??? !!! ???
        double averageCostPerIteration=0.0;
        double totalNumberIterations=1.0;
        int smallestNumberInfeasibilities=INT_MAX;
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
        if (smallestNumberInfeasibilities>=numberIntegerInfeasibilities)
          numberNodes=1000000; // switch off search for better solution
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
          model->messageHandler()->message(CBC_STRONG,model->messages())
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
                                               objectiveValue);
        // move branching object and make sure it will not be deleted
        if (whichObject>=0) {
          branch_ = objects[whichObject];
          choice[whichObject].possibleBranch=NULL;
        }
        delete [] changeUp;
        delete [] numberInfeasibilitiesUp;
        delete [] changeDown;
        delete [] numberInfeasibilitiesDown;
        delete [] objects;
      }
      if (osiclp&&!allNormal) {
        // back to normal
        osiclp->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,NULL) ;
      }
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
  objectiveValue_ = solver->getObjSense()*saveObjectiveValue;
  guessedObjectiveValue_ = objectiveValue_+estimatedDegradation;
/*
  Cleanup, then we're outta here.
*/
  if (!model->branchingMethod())
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

int CbcNode::chooseDynamicBranch (CbcModel *model, CbcNode *lastNode,int numberPassesLeft)

{ if (lastNode)
    depth_ = lastNode->depth_+1;
  else
    depth_ = 0;
  delete branch_;
  branch_=NULL;
  OsiSolverInterface * solver = model->solver();
  objectiveValue_ = solver->getObjSense()*solver->getObjValue();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  int anyAction=0;
  int i;
  int stateOfSearch = model->stateOfSearch();
  int numberStrong=model->numberStrong();
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
  int hotstartStrategy=model->getHotstartStrategy();
  if (hotstartStrategy>0) 
    return -3;
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

  for (i=0;i<numberColumns;i++) {
    saveLower[i] = lower[i];
    saveUpper[i] = upper[i];
  }
  // Get arrays to sort 
  double * sort = new double[numberObjects];
  int * whichObject = new int[numberObjects];
  int * whichColumn = new int[2*numberObjects+1];
  CbcStrongInfo * fixObject = new CbcStrongInfo[numberObjects];
  double estimatedDegradation=0.0; 
  // May go round twice if strong branching fixes all local candidates
  bool finished=false;
  while(!finished) {
    finished=true;
    decision->initialize(model);
    // Some objects may compute an estimate of best solution from here
    estimatedDegradation=0.0; 
    int numberToFix=0;
    int numberIntegerInfeasibilities=0; // without odd ones
    int numberToDo=0;
    double averageDown=0.0;
    int numberDown=0;
    double averageUp=0.0;
    int numberUp=0;
    int iBestNot=-1;
    int iBestGot=-1;
    double best=0.0;
    
    // We may go round this loop twice (only if we think we have solution)
    for (int iPass=0;iPass<2;iPass++) {
      
      // compute current state
      int numberObjectInfeasibilities; // just odd ones
      model->feasibleSolution(
                              numberIntegerInfeasibilities,
                              numberObjectInfeasibilities);
      
      // Some objects may compute an estimate of best solution from here
      estimatedDegradation=0.0; 
      numberUnsatisfied_ = 0;
      int bestPriority=INT_MAX;
      /*
        Scan for branching objects that indicate infeasibility. Choose candidates
        using priority as the first criteria, then integer infeasibility.
        
        The algorithm is to fill the array with a set of good candidates (by
        infeasibility) with priority bestPriority.  Finding a candidate with
        priority better (less) than bestPriority flushes the choice array. (This
        serves as initialization when the first candidate is found.)
        
      */
      numberToDo=0;
      averageDown=0.0;
      numberDown=0;
      averageUp=0.0;
      numberUp=0;
      iBestNot=-1;
      double bestNot=0.0;
      iBestGot=-1;
      best=0.0;
      for (i=0;i<numberObjects;i++) {
        CbcObject * object = model->modifiableObject(i);
        CbcSimpleIntegerDynamicPseudoCost * dynamicObject =
          dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object) ;
        assert(dynamicObject);
        int preferredWay;
        double infeasibility = object->infeasibility(preferredWay);
        int priorityLevel = object->priority();
        bool gotDown;
        if (dynamicObject->numberTimesDown()) {
          averageDown += dynamicObject->downDynamicPseudoCost();
          numberDown++;
          gotDown=true;
        } else {
          gotDown=false;
        }
        bool gotUp;
        if (dynamicObject->numberTimesUp()) {
          averageUp += dynamicObject->upDynamicPseudoCost();
          numberUp++;
          gotUp=true;
        } else {
          gotUp=false;
        }
        if (infeasibility) {
          // Increase estimated degradation to solution
          estimatedDegradation += CoinMin(object->upEstimate(),object->downEstimate());
          numberUnsatisfied_++;
          // Better priority? Flush choices.
          if (priorityLevel<bestPriority) {
            numberToDo=0;
            bestPriority = priorityLevel;
            iBestGot=-1;
            best=0.0;
          } else if (priorityLevel>bestPriority) {
            continue;
          }
          // Check for suitability based on infeasibility.
          if (gotDown||gotUp) {
            sort[numberToDo]=-infeasibility;
            if (infeasibility>best) {
              best=infeasibility;
              iBestGot=numberToDo;
            }
            whichColumn[numberToDo]=-1; // range info not wanted
          } else {
            int iColumn = dynamicObject->columnNumber();
            sort[numberToDo]=saveSolution[iColumn]-floor(saveSolution[iColumn]);
            if (fabs(sort[numberToDo])>bestNot) {
              iBestNot=numberToDo;
              bestNot = fabs(sort[numberToDo]);
            }
            whichColumn[numberToDo]=dynamicObject->columnNumber();
          }
          whichObject[numberToDo++]=i;
        }
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
        if (roundAgain) {
          // restore basis
          solver->setWarmStart(ws);
          solver->setColSolution(saveSolution);
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
    if (anyAction==-2) {
      break;
    }
    bool solveAll=false; // set true to say look at all even if some fixed (experiment)
    solveAll=true;
    // skip if solution
    if (!numberUnsatisfied_)
      break;
    // worth trying if too many times
    // Save basis
    CoinWarmStart * ws = solver->getWarmStart();
    // save limit
    int saveLimit;
    solver->getIntParam(OsiMaxNumIterationHotStart,saveLimit);
    if (!stateOfSearch&&saveLimit<100)
      solver->setIntParam(OsiMaxNumIterationHotStart,100); 
    
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
    if (hitMaxTime||numberPassesLeft<=0) {
      int iObject = whichObject[bestChoice];
      CbcObject * object = model->modifiableObject(iObject);
      int preferredWay;
      object->infeasibility(preferredWay);
      branch_=object->createBranch(preferredWay);
      branch_->way(preferredWay);
      break;
    } else {
      // say do fast
      int easy=1;
      solver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,&easy) ;
      if (numberDown)
        averageDown /= (double) numberDown;
      else
        averageDown=1.0;
      if (numberUp)
        averageUp /= (double) numberUp;
      else
        averageUp=1.0;
      double weight;
      if (!stateOfSearch)
        weight=0.5;
      else
        weight=0.8;
      double average = 0.5*(averageUp+averageDown);
      int iDo;
#define RANGING
#ifdef RANGING
      OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (solver);
      if (osiclp&&iBestNot>=0) {
        // just get those not touched and best and not trusted
        int needed=0;
        int * which = whichColumn+numberToDo+1;
        int i;
        for ( i=0;i<numberToDo;i++) {
          if (whichColumn[i]>=0) {
            which[needed++]=whichColumn[i];
          }
        }
        which--;
        which[0]=needed;
        assert (needed);
        osiclp->passInRanges(which);
        // Mark hot start and get ranges
        solver->markHotStart();
        osiclp->passInRanges(NULL);
        int put=0;
        int get=0;
        const double * downRange=osiclp->downRange();
        const double * upRange=osiclp->upRange();;
        int numberBeforeTrust = model->numberBeforeTrust();
        for ( i=0;i<numberToDo;i++) {
          int iObject = whichObject[i];
          CbcObject * object = model->modifiableObject(iObject);
          CbcSimpleIntegerDynamicPseudoCost * dynamicObject =
            dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object) ;
          double estimate = sort[i];
          bool keep=false;
          if (i==iBestGot) {
            // keep
            keep=true;
            estimate=-COIN_DBL_MAX;
          } else {
            int numberDown = dynamicObject->numberTimesDown();
            int numberUp = dynamicObject->numberTimesUp();
            if (numberDown<numberBeforeTrust||
                numberUp<numberBeforeTrust) {
              keep=true;
              // make positive
              estimate += best + 0.1*average*(2*numberBeforeTrust-numberUp-numberDown);
              if (whichColumn[i]>=0) {
                estimate = - (downRange[get]+upRange[get]);
                get++;
              }
            }
          }
          if (keep) {
            sort[put]=estimate;
            whichObject[put++]=iObject;
          }
        }
        numberToDo=put;
      } else {
        // Mark hot start
        solver->markHotStart();
        // make sure best will be first
        if (iBestGot>=0)
          sort[iBestGot]=-COIN_DBL_MAX;
      }
#else
      // Mark hot start
      solver->markHotStart();
      // make sure best will be first
      if (iBestGot>=0)
        sort[iBestGot]=-COIN_DBL_MAX;
#endif
      // Sort
      CoinSort_2(sort,sort+numberToDo,whichObject);
      iDo=0;
      int saveLimit2;
      solver->getIntParam(OsiMaxNumIterationHotStart,saveLimit2);
      for ( iDo=0;iDo<numberToDo;iDo++) {
        CbcStrongInfo choice;
        int iObject = whichObject[iDo];
        CbcObject * object = model->modifiableObject(iObject);
        CbcSimpleIntegerDynamicPseudoCost * dynamicObject =
          dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object) ;
        int preferredWay;
        object->infeasibility(preferredWay);
        choice.possibleBranch=object->createBranch(preferredWay);
        // Save which object it was
        choice.objectNumber=iObject;
        choice.numIntInfeasUp = numberUnsatisfied_;
        choice.numIntInfeasDown = numberUnsatisfied_;
        choice.fix=0; // say not fixed
        // see if can skip strong branching
        int canSkip = choice.possibleBranch->fillStrongInfo(choice);
        // For now always do
        //canSkip=false;
        if (model->messageHandler()->logLevel()>3) 
          dynamicObject->print(1,choice.possibleBranch->value());
        // was if (!canSkip)
        if (!canSkip) {
          // just do a few
          if (canSkip)
            solver->setIntParam(OsiMaxNumIterationHotStart,10); 
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
          decision->saveBranchingObject( choice.possibleBranch);
          choice.possibleBranch->branch() ;
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
          objectiveChange = newObjectiveValue  - objectiveValue_;
          decision->updateInformation( solver,this);
          if (!iStatus) {
            choice.finishedDown = true ;
            if (newObjectiveValue>=model->getCutoff()) {
              objectiveChange = 1.0e100; // say infeasible
            } else {
              // See if integer solution
              if (model->feasibleSolution(choice.numIntInfeasDown,
                                          choice.numObjInfeasDown)) {
                model->setBestSolution(CBC_STRONGSOL,
                                       newObjectiveValue,
                                       solver->getColSolution()) ;
                model->incrementUsed(solver->getColSolution());
                if (newObjectiveValue >= model->getCutoff())	//  *new* cutoff
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
          
          // restore bounds
          for ( j=0;j<numberColumns;j++) {
            if (saveLower[j] != lower[j])
              solver->setColLower(j,saveLower[j]);
            if (saveUpper[j] != upper[j])
              solver->setColUpper(j,saveUpper[j]);
          }
          //printf("Down on %d, status is %d, obj %g its %d cost %g finished %d inf %d infobj %d\n",
          //     choice.objectNumber,iStatus,newObjectiveValue,choice.numItersDown,
          //     choice.downMovement,choice.finishedDown,choice.numIntInfeasDown,
          //     choice.numObjInfeasDown);
          
          // repeat the whole exercise, forcing the variable up
          decision->saveBranchingObject( choice.possibleBranch);
          choice.possibleBranch->branch();
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
          objectiveChange = newObjectiveValue  - objectiveValue_;
          decision->updateInformation( solver,this);
          if (!iStatus) {
            choice.finishedUp = true ;
            if (newObjectiveValue>=model->getCutoff()) {
              objectiveChange = 1.0e100; // say infeasible
            } else {
              // See if integer solution
              if (model->feasibleSolution(choice.numIntInfeasUp,
                                          choice.numObjInfeasUp)) { 
                model->setBestSolution(CBC_STRONGSOL,
                                       newObjectiveValue,
                                       solver->getColSolution()) ;
                model->incrementUsed(solver->getColSolution());
                if (newObjectiveValue >= model->getCutoff())	//  *new* cutoff
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
          
          //printf("Up on %d, status is %d, obj %g its %d cost %g finished %d inf %d infobj %d\n",
          //     choice.objectNumber,iStatus,newObjectiveValue,choice.numItersUp,
          //     choice.upMovement,choice.finishedUp,choice.numIntInfeasUp,
          //     choice.numObjInfeasUp);
          hitMaxTime = ( CoinCpuTime()-model->getDblParam(CbcModel::CbcStartSeconds) > 
                         model->getDblParam(CbcModel::CbcMaximumSeconds));
          if (hitMaxTime) {
            break;
          }
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
        if (choice.upMovement<1.0e100) {
          if(choice.downMovement<1.0e100) {
            // feasible - see which best
            int iColumn =
              model->integerVariable()[choice.possibleBranch->variable()] ;
            model->messageHandler()->message(CBC_STRONG,model->messages())
              << iObject << iColumn
              <<choice.downMovement<<choice.numIntInfeasDown 
              <<choice.upMovement<<choice.numIntInfeasUp 
              <<choice.possibleBranch->value()
              <<CoinMessageEol;
            //if (!stateOfSearch)
            //choice.numIntInfeasDown=99999; // temp fudge
            int betterWay = decision->betterBranch(choice.possibleBranch,
                                                   branch_,
                                                   choice.upMovement*changeFactor,
                                                   choice.numIntInfeasUp ,
                                                   choice.downMovement*changeFactor,
                                                   choice.numIntInfeasDown );
            if (iDo>=changeStrategy) {
              // make less likely
              changeStrategy+=numberStrong;
              changeFactor *= 0.9;
            }
            if (betterWay) {
              delete branch_;
              // C) create branching object
              branch_ = choice.possibleBranch;
              choice.possibleBranch=NULL;
              branch_->way(betterWay);
              bestChoice = choice.objectNumber;
              whichChoice = iDo;
            } else {
              delete choice.possibleBranch;
              if (iDo>=2*numberStrong)
                break;
              if (!dynamicObject||dynamicObject->numberTimesUp()>1) {
                if (iDo-whichChoice>=numberStrong)
                  break; // give up
              } else {
                if (iDo-whichChoice>=2*numberStrong)
                  break; // give up
              }
            }
          } else {
            // up feasible, down infeasible
            anyAction=-1;
            //printf("Down infeasible for choice %d sequence %d\n",i,
            // model->object(choice.objectNumber)->columnNumber());
            if (!solveAll) {
              choice.possibleBranch->way(1);
              choice.possibleBranch->branch();
              delete choice.possibleBranch;
              break;
            } else {
              choice.fix=1;
              fixObject[numberToFix++]=choice;
            }
          }
        } else {
          if(choice.downMovement<1.0e100) {
            // down feasible, up infeasible
            anyAction=-1;
            //printf("Up infeasible for choice %d sequence %d\n",i,
            // model->object(choice.objectNumber)->columnNumber());
            if (!solveAll) {
              choice.possibleBranch->way(-1);
              choice.possibleBranch->branch();
              delete choice.possibleBranch;
              break;
            } else {
              choice.fix=-1;
              fixObject[numberToFix++]=choice;
            }
          } else {
            // neither side feasible
            anyAction=-2;
            delete choice.possibleBranch;
            //printf("Both infeasible for choice %d sequence %d\n",i,
            // model->object(choice.objectNumber)->columnNumber());
            break;
          }
        }
      }
      if (model->messageHandler()->logLevel()>3) { 
        if (anyAction==-2)
          printf("infeasible\n");
        else if(anyAction==-1)
          printf("%d fixed\n",numberToFix);
        else
          printf("choosing %d iDo %d numberToDo %d\n",bestChoice,whichChoice,numberToDo);
      }
      // Delete the snapshot
      solver->unmarkHotStart();
      // back to normal
      solver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,NULL) ;
      solver->setIntParam(OsiMaxNumIterationHotStart,saveLimit);
      // restore basis
      solver->setWarmStart(ws);
      // Unless infeasible we will carry on
      // But we could fix anyway
      if (numberToFix) {
        if (anyAction==-2) {
          // take off
          for (i = 0 ; i < numberToFix ; i++) {
            delete fixObject[i].possibleBranch;
          }
        } else {
          // apply and take off
          for (i = 0 ; i < numberToFix ; i++) {
            fixObject[i].possibleBranch->way(fixObject[i].fix) ;
            fixObject[i].possibleBranch->branch() ;
            delete fixObject[i].possibleBranch;
          }
          bool feasible=true;
          // can do quick optimality check
          int easy=2;
          solver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,&easy) ;
          solver->resolve() ;
          solver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,NULL) ;
          feasible = solver->isProvenOptimal();
          if (feasible) {
            anyAction=0;
            memcpy(saveSolution,solver->getColSolution(),numberColumns*sizeof(double));
            model->reserveCurrentSolution(saveSolution);
            memcpy(saveLower,solver->getColLower(),numberColumns*sizeof(double));
            memcpy(saveUpper,solver->getColUpper(),numberColumns*sizeof(double));
            // See if candidate still possible
            if (branch_) {
              const CbcObject * object = model->object(bestChoice);
              int preferredWay;
              double infeasibility = object->infeasibility(preferredWay);
              if (!infeasibility) {
                // take out
                delete branch_;
                branch_=NULL;
              } else {
                branch_->way(preferredWay);
              }
            }
          } else {
            anyAction=-2;
            finished=true;
          }
          // If  fixed then round again
          if (!branch_) {
            finished=false;
          }
          // If these in then different action
#if 0
          //if (!anyAction)
          //anyAction=-1;
          //finished=true;
#endif
        }
      }
      delete ws;
    }
  }
  
  // Set guessed solution value
  guessedObjectiveValue_ = objectiveValue_+estimatedDegradation;
/*
  Cleanup, then we're outta here.
*/
  if (!model->branchingMethod())
    delete decision;
    
  delete [] fixObject;
  delete [] sort;
  delete [] whichObject;
  delete [] whichColumn;
  delete [] saveLower;
  delete [] saveUpper;
  
  // restore solution
  solver->setColSolution(saveSolution);
  model->reserveCurrentSolution(saveSolution);
  delete [] saveSolution;
  return anyAction;
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
  if (rhs.branch_)
    branch_=rhs.branch_->clone();
  else
    branch_=NULL;
  depth_ = rhs.depth_;
  numberUnsatisfied_ = rhs.numberUnsatisfied_;
  nodeNumber_=rhs.nodeNumber_;
}

CbcNode &
CbcNode::operator=(const CbcNode & rhs)
{
  if (this != &rhs) {
    delete nodeInfo_;
    if (nodeInfo_)
      nodeInfo_ = rhs.nodeInfo_->clone();
    else
      nodeInfo_ = NULL;
    objectiveValue_=rhs.objectiveValue_;
    guessedObjectiveValue_ = rhs.guessedObjectiveValue_;
    if (rhs.branch_)
      branch_=rhs.branch_->clone();
    else
      branch_=NULL,
    depth_ = rhs.depth_;
    numberUnsatisfied_ = rhs.numberUnsatisfied_;
    nodeNumber_=rhs.nodeNumber_;
  }
  return *this;
}


CbcNode::~CbcNode ()
{
#ifdef CHECK_NODE
  if (nodeInfo_) 
    printf("CbcNode %x Destructor nodeInfo %x (%d)\n",
	 this,nodeInfo_,nodeInfo_->numberPointingToThis());
  else
    printf("CbcNode %x Destructor nodeInfo %x (?)\n",
	 this,nodeInfo_);
#endif
  if (nodeInfo_) {
    nodeInfo_->nullOwner();
    int numberToDelete=nodeInfo_->numberBranchesLeft();
    //    CbcNodeInfo * parent = nodeInfo_->parent();
    if (nodeInfo_->decrement(numberToDelete)==0) {
      delete nodeInfo_;
    } else {
      //printf("node %x nodeinfo %x parent %x\n",this,nodeInfo_,parent);
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
  if(nodeInfo_) {
    nodeInfo_->decrementCuts(change);
  }
}
void 
CbcNode::decrementParentCuts(int change)
{
  if(nodeInfo_) {
    nodeInfo_->decrementParentCuts(change);
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
}
// Nulls out node info
void 
CbcNode::nullNodeInfo()
{
  nodeInfo_=NULL;
}

int
CbcNode::branch()
{
  double changeInGuessed=branch_->branch(true);
  guessedObjectiveValue_+= changeInGuessed;
  return nodeInfo_->branchedOn();
}
