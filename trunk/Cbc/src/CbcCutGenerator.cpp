// Copyright (C) 2003, International Business Machines
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
#include "CbcCutGenerator.hpp"
#include "CbcBranchDynamic.hpp"
#include "CglProbing.hpp"
#include "CoinTime.hpp"

// Default Constructor 
CbcCutGenerator::CbcCutGenerator ()
  : model_(NULL),
    generator_(NULL),
    whenCutGenerator_(-1),
    whenCutGeneratorInSub_(-100),
    switchOffIfLessThan_(0),
    depthCutGenerator_(-1),
    depthCutGeneratorInSub_(-1),
    generatorName_(NULL),
    normal_(true),
    atSolution_(false),
    whenInfeasible_(false),
    mustCallAgain_(false),
    switchedOff_(false),
    timing_(false),
    timeInCutGenerator_(0.0),
    numberTimes_(0),
    numberCuts_(0),
    numberColumnCuts_(0),
    numberCutsActive_(0),
    numberCutsAtRoot_(0),
    numberActiveCutsAtRoot_(0)
{
}
// Normal constructor
CbcCutGenerator::CbcCutGenerator(CbcModel * model,CglCutGenerator * generator,
				 int howOften, const char * name,
				 bool normal, bool atSolution, 
				 bool infeasible, int howOftenInSub,
				 int whatDepth, int whatDepthInSub,
                                 int switchOffIfLessThan)
  : 
    depthCutGenerator_(whatDepth),
    depthCutGeneratorInSub_(whatDepthInSub),
    mustCallAgain_(false),
    switchedOff_(false),
    timing_(false),
    timeInCutGenerator_(0.0),
    numberTimes_(0),
    numberCuts_(0),
    numberColumnCuts_(0),
    numberCutsActive_(0),
    numberCutsAtRoot_(0),
    numberActiveCutsAtRoot_(0)
{
  model_ = model;
  generator_=generator->clone();
  generator_->refreshSolver(model_->solver());
  whenCutGenerator_=howOften;
  whenCutGeneratorInSub_ = howOftenInSub;
  switchOffIfLessThan_=switchOffIfLessThan;
  if (name)
    generatorName_=strdup(name);
  else
    generatorName_ = strdup("Unknown");
  normal_=normal;
  atSolution_=atSolution;
  whenInfeasible_=infeasible;
}

// Copy constructor 
CbcCutGenerator::CbcCutGenerator ( const CbcCutGenerator & rhs)
{
  model_ = rhs.model_;
  generator_=rhs.generator_->clone();
  //generator_->refreshSolver(model_->solver());
  whenCutGenerator_=rhs.whenCutGenerator_;
  whenCutGeneratorInSub_ = rhs.whenCutGeneratorInSub_;
  switchOffIfLessThan_ = rhs.switchOffIfLessThan_;
  depthCutGenerator_=rhs.depthCutGenerator_;
  depthCutGeneratorInSub_ = rhs.depthCutGeneratorInSub_;
  generatorName_=strdup(rhs.generatorName_);
  normal_=rhs.normal_;
  atSolution_=rhs.atSolution_;
  whenInfeasible_=rhs.whenInfeasible_;
  mustCallAgain_ = rhs.mustCallAgain_;
  switchedOff_ = rhs.switchedOff_;
  timing_ = rhs.timing_;
  timeInCutGenerator_ = rhs.timeInCutGenerator_;
  numberTimes_ = rhs.numberTimes_;
  numberCuts_ = rhs.numberCuts_;
  numberColumnCuts_ = rhs.numberColumnCuts_;
  numberCutsActive_ = rhs.numberCutsActive_;
  numberCutsAtRoot_  = rhs.numberCutsAtRoot_;
  numberActiveCutsAtRoot_ = rhs.numberActiveCutsAtRoot_;
}

// Assignment operator 
CbcCutGenerator & 
CbcCutGenerator::operator=( const CbcCutGenerator& rhs)
{
  if (this!=&rhs) {
    delete generator_;
    free(generatorName_);
    model_ = rhs.model_;
    generator_=rhs.generator_->clone();
    generator_->refreshSolver(model_->solver());
    whenCutGenerator_=rhs.whenCutGenerator_;
    whenCutGeneratorInSub_ = rhs.whenCutGeneratorInSub_;
    switchOffIfLessThan_ = rhs.switchOffIfLessThan_;
    depthCutGenerator_=rhs.depthCutGenerator_;
    depthCutGeneratorInSub_ = rhs.depthCutGeneratorInSub_;
    generatorName_=strdup(rhs.generatorName_);
    normal_=rhs.normal_;
    atSolution_=rhs.atSolution_;
    whenInfeasible_=rhs.whenInfeasible_;
    mustCallAgain_ = rhs.mustCallAgain_;
    switchedOff_ = rhs.switchedOff_;
    timing_ = rhs.timing_;
    timeInCutGenerator_ = rhs.timeInCutGenerator_;
    numberTimes_ = rhs.numberTimes_;
    numberCuts_ = rhs.numberCuts_;
    numberColumnCuts_ = rhs.numberColumnCuts_;
    numberCutsActive_ = rhs.numberCutsActive_;
    numberCutsAtRoot_  = rhs.numberCutsAtRoot_;
    numberActiveCutsAtRoot_ = rhs.numberActiveCutsAtRoot_;
  }
  return *this;
}

// Destructor 
CbcCutGenerator::~CbcCutGenerator ()
{
  free(generatorName_);
  delete generator_;
}

/* This is used to refresh any inforamtion.
   It also refreshes the solver in the cut generator
   in case generator wants to do some work 
*/
void 
CbcCutGenerator::refreshModel(CbcModel * model)
{
  model_=model;
  generator_->refreshSolver(model_->solver());
}
/* Generate cuts for the model data contained in si.
   The generated cuts are inserted into and returned in the
   collection of cuts cs.
*/
bool
CbcCutGenerator::generateCuts( OsiCuts & cs , bool fullScan, OsiSolverInterface * solver, CbcNode * node)
{
  int howOften = whenCutGenerator_;
  if (howOften==-100)
    return false;
  if (howOften>0)
    howOften = howOften % 1000000;
  else 
    howOften=1;
  if (!howOften)
    howOften=1;
  bool returnCode=false;
  //OsiSolverInterface * solver = model_->solver();
  int depth;
  if (node)
    depth=node->depth();
  else
    depth=0;
  int pass=model_->getCurrentPassNumber()-1;
  bool doThis=(model_->getNodeCount()%howOften)==0;
  if (depthCutGenerator_>0) {
    doThis = (depth % depthCutGenerator_) ==0;
    if (depth<depthCutGenerator_)
      doThis=true; // and also at top of tree
  }
  // But turn off if 100
  if (howOften==100)
    doThis=false;
  // Switch off if special setting
  if (whenCutGeneratorInSub_==-200) {
    fullScan=false;
    doThis=false;
  }
  if (fullScan||doThis) {
    double time1=0.0;
    if (timing_)
      time1 = CoinCpuTime();
    //#define CBC_DEBUG
    int numberRowCutsBefore = cs.sizeRowCuts() ;
    int cutsBefore = cs.sizeCuts();
    CglTreeInfo info;
    info.level = depth;
    info.pass = pass;
    info.formulation_rows = model_->numberRowsAtContinuous();
    info.inTree = node!=NULL;
    incrementNumberTimesEntered();
    CglProbing* generator =
      dynamic_cast<CglProbing*>(generator_);
    if (!generator) {
      // Pass across model information in case it could be useful
      //void * saveData = solver->getApplicationData();
      //solver->setApplicationData(model_);
      generator_->generateCuts(*solver,cs,info);
      //solver->setApplicationData(saveData);
    } else {
      // Probing - return tight column bounds
      CglTreeProbingInfo * info2 = model_->probingInfo();
      if (info2) {
	info2->level = depth;
	info2->pass = pass;
	info2->formulation_rows = model_->numberRowsAtContinuous();
	info2->inTree = node!=NULL;
	generator->generateCutsAndModify(*solver,cs,info2);
      } else {
	generator->generateCutsAndModify(*solver,cs,&info);
      }
      const double * tightLower = generator->tightLower();
      const double * lower = solver->getColLower();
      const double * tightUpper = generator->tightUpper();
      const double * upper = solver->getColUpper();
      const double * solution = solver->getColSolution();
      int j;
      int numberColumns = solver->getNumCols();
      double primalTolerance = 1.0e-8;
      const char * tightenBounds = generator->tightenBounds();
      if ((model_->getThreadMode()&2)==0) {
	for (j=0;j<numberColumns;j++) {
	  if (solver->isInteger(j)) {
	    if (tightUpper[j]<upper[j]) {
	      double nearest = floor(tightUpper[j]+0.5);
	      //assert (fabs(tightUpper[j]-nearest)<1.0e-5); may be infeasible
	      solver->setColUpper(j,nearest);
	      if (nearest<solution[j]-primalTolerance)
		returnCode=true;
	    }
	    if (tightLower[j]>lower[j]) {
	      double nearest = floor(tightLower[j]+0.5);
	      //assert (fabs(tightLower[j]-nearest)<1.0e-5); may be infeasible
	      solver->setColLower(j,nearest);
	      if (nearest>solution[j]+primalTolerance)
		returnCode=true;
	    }
	  } else {
	    if (upper[j]>lower[j]) {
	      if (tightUpper[j]==tightLower[j]) {
		// fix
		solver->setColLower(j,tightLower[j]);
		solver->setColUpper(j,tightUpper[j]);
		if (tightLower[j]>solution[j]+primalTolerance||
		    tightUpper[j]<solution[j]-primalTolerance)
		  returnCode=true;
	      } else if (tightenBounds&&tightenBounds[j]) {
		solver->setColLower(j,CoinMax(tightLower[j],lower[j]));
		solver->setColUpper(j,CoinMin(tightUpper[j],upper[j]));
		if (tightLower[j]>solution[j]+primalTolerance||
		    tightUpper[j]<solution[j]-primalTolerance)
		  returnCode=true;
	      }
	    }
	  }
	}
      } else {
	CoinPackedVector lbs;
	CoinPackedVector ubs;
	int numberChanged=0;
	bool ifCut=false;
	for (j=0;j<numberColumns;j++) {
	  if (solver->isInteger(j)) {
	    if (tightUpper[j]<upper[j]) {
	      double nearest = floor(tightUpper[j]+0.5);
	      //assert (fabs(tightUpper[j]-nearest)<1.0e-5); may be infeasible
	      ubs.insert(j,nearest);
	      numberChanged++;
	      if (nearest<solution[j]-primalTolerance)
		ifCut=true;
	    }
	    if (tightLower[j]>lower[j]) {
	      double nearest = floor(tightLower[j]+0.5);
	      //assert (fabs(tightLower[j]-nearest)<1.0e-5); may be infeasible
	      lbs.insert(j,nearest);
	      numberChanged++;
	      if (nearest>solution[j]+primalTolerance)
		ifCut=true;
	    }
	  } else {
	    if (upper[j]>lower[j]) {
	      if (tightUpper[j]==tightLower[j]) {
		// fix
		lbs.insert(j,tightLower[j]);
		ubs.insert(j,tightUpper[j]);
		if (tightLower[j]>solution[j]+primalTolerance||
		    tightUpper[j]<solution[j]-primalTolerance)
		  ifCut=true;
	      } else if (tightenBounds&&tightenBounds[j]) {
		lbs.insert(j,CoinMax(tightLower[j],lower[j]));
		ubs.insert(j,CoinMin(tightUpper[j],upper[j]));
		if (tightLower[j]>solution[j]+primalTolerance||
		    tightUpper[j]<solution[j]-primalTolerance)
		  ifCut=true;
	      }
	    }
	  }
	}
	if (numberChanged) {
	  OsiColCut cc;
	  cc.setUbs(ubs);
	  cc.setLbs(lbs);
	  if (ifCut) {
	    cc.setEffectiveness(100.0);
	  } else {
	    cc.setEffectiveness(1.0e-5);
	  }
	  cs.insert(cc);
	}
      }
      //if (!solver->basisIsAvailable()) 
      //returnCode=true;
#if 0
      // Pass across info to pseudocosts
      char * mark = new char[numberColumns];
      memset(mark,0,numberColumns);
      int nLook = generator->numberThisTime();
      const int * lookedAt = generator->lookedAt();
      const int * fixedDown = generator->fixedDown();
      const int * fixedUp = generator->fixedUp();
      for (j=0;j<nLook;j++) 
	mark[lookedAt[j]]=1;
      int numberObjects = model_->numberObjects();
      for (int i=0;i<numberObjects;i++) {
	CbcSimpleIntegerDynamicPseudoCost * obj1 =
	  dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(model_->modifiableObject(i)) ;
	if (obj1) {
	  int iColumn = obj1->columnNumber();
	  if (mark[iColumn]) 
	    obj1->setProbingInformation(fixedDown[iColumn],fixedUp[iColumn]);
	}
      }
      delete [] mark;
#endif
    }
    CbcCutModifier * modifier = model_->cutModifier();
    if (modifier) {
      int numberRowCutsAfter = cs.sizeRowCuts() ;
      int k ;
      int nOdd=0;
      //const OsiSolverInterface * solver = model_->solver();
      for (k = numberRowCutsAfter-1;k>=numberRowCutsBefore;k--) {
	OsiRowCut & thisCut = cs.rowCut(k) ;
	int returnCode = modifier->modify(solver,thisCut);
	if (returnCode) {
	  nOdd++;
	  if (returnCode==3)
	    cs.eraseRowCut(k);
	}
      }
      if (nOdd) 
	printf("Cut generator %s produced %d cuts of which %d were modified\n",
		 generatorName_,numberRowCutsAfter-numberRowCutsBefore,nOdd);
    }
    { 
      // make all row cuts without test for duplicate
      int numberRowCutsAfter = cs.sizeRowCuts() ;
      int k ;
      for (k = numberRowCutsBefore;k<numberRowCutsAfter;k++) {
	OsiRowCut * thisCut = cs.rowCutPtr(k) ;
	thisCut->mutableRow().setTestForDuplicateIndex(false);
      }
    }
#ifdef CBC_DEBUG
    {
      int numberRowCutsAfter = cs.sizeRowCuts() ;
      int k ;
      int nBad=0;
      for (k = numberRowCutsBefore;k<numberRowCutsAfter;k++) {
	OsiRowCut thisCut = cs.rowCut(k) ;
	if (thisCut.lb()>thisCut.ub()||
	    thisCut.lb()>1.0e8||
	    thisCut.ub()<-1.0e8)
	  printf("cut from %s has bounds %g and %g!\n",
		 generatorName_,thisCut.lb(),thisCut.ub());
	if (thisCut.lb()<=thisCut.ub()) {
	  /* check size of elements.
	     We can allow smaller but this helps debug generators as it
	     is unsafe to have small elements */
	  int n=thisCut.row().getNumElements();
	  const int * column = thisCut.row().getIndices();
	  const double * element = thisCut.row().getElements();
	  assert (n);
	  for (int i=0;i<n;i++) {
	    double value = element[i];
	    if (fabs(value)<=1.0e-12||fabs(value)>=1.0e20)
	      nBad++;
	  }
	}
	if (nBad) 
	  printf("Cut generator %s produced %d cuts of which %d had tiny or large elements\n",
		 generatorName_,numberRowCutsAfter-numberRowCutsBefore,nBad);
      }
    }
#endif
    if (timing_)
      timeInCutGenerator_ += CoinCpuTime()-time1;
    // switch off if first time and no good
    if (node==NULL&&!pass) {
      if (cs.sizeCuts()-cutsBefore<CoinAbs(switchOffIfLessThan_)) {
        whenCutGenerator_=-99;
        whenCutGeneratorInSub_ = -200;
      }
    }
  }
  return returnCode;
}
void 
CbcCutGenerator::setHowOften(int howOften) 
{
  
  if (howOften>=1000000) {
    // leave Probing every 10
    howOften = howOften % 1000000;
    CglProbing* generator =
      dynamic_cast<CglProbing*>(generator_);
    
    if (generator&&howOften>10) 
      howOften=10+1000000;
    else
      howOften += 1000000;
  }
  whenCutGenerator_ = howOften;
}
void 
CbcCutGenerator::setWhatDepth(int value) 
{
  depthCutGenerator_ = value;
}
void 
CbcCutGenerator::setWhatDepthInSub(int value) 
{
  depthCutGeneratorInSub_ = value;
}


// Default Constructor
CbcCutModifier::CbcCutModifier() 
{
}


// Destructor 
CbcCutModifier::~CbcCutModifier ()
{
}

// Copy constructor 
CbcCutModifier::CbcCutModifier ( const CbcCutModifier & rhs)
{
}

// Assignment operator 
CbcCutModifier & 
CbcCutModifier::operator=( const CbcCutModifier& rhs)
{
  if (this!=&rhs) {
  }
  return *this;
}

// Default Constructor 
CbcCutSubsetModifier::CbcCutSubsetModifier ()
  : CbcCutModifier(),
    firstOdd_(COIN_INT_MAX)
{
}

// Useful constructor 
CbcCutSubsetModifier::CbcCutSubsetModifier (int firstOdd)
  : CbcCutModifier()
{
  firstOdd_=firstOdd;
}

// Copy constructor 
CbcCutSubsetModifier::CbcCutSubsetModifier ( const CbcCutSubsetModifier & rhs)
  :CbcCutModifier(rhs)
{
  firstOdd_ = rhs.firstOdd_;
}

// Clone
CbcCutModifier *
CbcCutSubsetModifier::clone() const
{
  return new CbcCutSubsetModifier(*this);
}

// Assignment operator 
CbcCutSubsetModifier & 
CbcCutSubsetModifier::operator=( const CbcCutSubsetModifier& rhs)
{
  if (this!=&rhs) {
    CbcCutModifier::operator=(rhs);
    firstOdd_ = rhs.firstOdd_;
  }
  return *this;
}

// Destructor 
CbcCutSubsetModifier::~CbcCutSubsetModifier ()
{
}
/* Returns
   0 unchanged
   1 strengthened
   2 weakened
   3 deleted
*/
int 
CbcCutSubsetModifier::modify(const OsiSolverInterface * solver, OsiRowCut & cut) 
{
  int n=cut.row().getNumElements();
  if (!n)
    return 0;
  const int * column = cut.row().getIndices();
  //const double * element = cut.row().getElements();
  int returnCode=0;
  for (int i=0;i<n;i++) {
    if (column[i]>=firstOdd_) {
      returnCode=3;
      break;
    }
  }
  if (!returnCode) {
    const double * element = cut.row().getElements();
    printf("%g <= ",cut.lb());
    for (int i=0;i<n;i++) {
      printf("%g*x%d ",element[i],column[i]);
    }
    printf("<= %g\n",cut.ub());
  }
  //return 3;
  return returnCode;
}

