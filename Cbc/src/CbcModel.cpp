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
//#define CHECK_NODE_FULL
//#define NODE_LOG
//#define GLOBAL_CUTS_JUST_POINTERS
#ifndef CLP_FAST_CODE
#ifdef NDEBUG
#undef NDEBUG
#endif
#endif
#include <cassert>
#include <cmath>
#include <cfloat>

#ifdef COIN_HAS_CLP
// include Presolve from Clp
#include "ClpPresolve.hpp"
#include "OsiClpSolverInterface.hpp"
#endif

#include "CbcEventHandler.hpp"

#include "OsiSolverInterface.hpp"
#include "OsiAuxInfo.hpp"
#include "OsiSolverBranch.hpp"
#include "OsiChooseVariable.hpp"
#include "CoinWarmStartBasis.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinHelperFunctions.hpp"
#include "CbcBranchActual.hpp"
#include "CbcBranchDynamic.hpp"
#include "CbcHeuristic.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcModel.hpp"
#include "CbcTreeLocal.hpp"
#include "CbcStatistics.hpp"
#include "CbcStrategy.hpp"
#include "CbcMessage.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiCuts.hpp"
#include "CbcCountRowCut.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcFeasibilityBase.hpp"
#include "CbcFathom.hpp"
// include Probing
#include "CglProbing.hpp"
// include preprocessing
#include "CglPreProcess.hpp"
#include "CglDuplicateRow.hpp"
#include "CglStored.hpp"

#include "CoinTime.hpp"
#include "CoinMpsIO.hpp"

#include "CbcCompareActual.hpp"
#include "CbcTree.hpp"
//#define CBC_THREAD
#ifdef CBC_THREAD
#include <pthread.h>
#ifndef CLP_FAST_CODE
//#define CBC_THREAD_DEBUG 1
#endif
#ifdef CBC_THREAD_DEBUG 
#ifdef NDEBUG
#undef NDEBUG
#undef assert
#         define assert(expression) { 				   \
             if (!(expression)) {					   \
                throw CoinError(__STRING(expression), __PRETTY_FUNCTION__, \
                                "", __FILE__, __LINE__);		   \
             }								   \
          }
#endif
#endif
// To Pass across to doOneNode
typedef struct {
  CbcModel * baseModel;
  CbcModel * thisModel;
  CbcNode * node; // filled in every time
  CbcNode * createdNode; // filled in every time on return
  pthread_t threadIdOfBase;
  pthread_mutex_t * mutex; // for locking data
  pthread_mutex_t * mutex2; // for waking up threads
  pthread_cond_t * condition2; // for waking up thread
  int returnCode; // -1 available, 0 busy, 1 finished , 2??
  double timeLocked;
  double timeWaitingToLock;
  double timeWaitingToStart;
  double timeInThread;
  int numberTimesLocked;
  int numberTimesUnlocked;
  int numberTimesWaitingToStart;
  int saveStuff[2];
  struct timespec absTime;
  bool locked;
#if CBC_THREAD_DEBUG
  int threadNumber;
#endif
} threadStruct;
static void * doNodesThread(void * voidInfo);
static void * doCutsThread(void * voidInfo);
#endif
/* Various functions local to CbcModel.cpp */

namespace {

//-------------------------------------------------------------------
// Returns the greatest common denominator of two 
// positive integers, a and b, found using Euclid's algorithm 
//-------------------------------------------------------------------
static int gcd(int a, int b) 
{
  int remainder = -1;
  // make sure a<=b (will always remain so)
  if(a > b) {
    // Swap a and b
    int temp = a;
    a = b;
    b = temp;
  }
  // if zero then gcd is nonzero (zero may occur in rhs of packed)
  if (!a) {
    if (b) {
      return b;
    } else {
      printf("**** gcd given two zeros!!\n");
      abort();
    }
  }
  while (remainder) {
    remainder = b % a;
    b = a;
    a = remainder;
  }
  return b;
}



#ifdef CHECK_NODE_FULL

/*
  Routine to verify that tree linkage is correct. The invariant that is tested
  is

  reference count = (number of actual references) + (number of branches left)

  The routine builds a set of paired arrays, info and count, by traversing the
  tree. Each CbcNodeInfo is recorded in info, and the number of times it is
  referenced (via the parent field) is recorded in count. Then a final check is
  made to see if the numberPointingToThis_ field agrees.
*/

void verifyTreeNodes (const CbcTree * branchingTree, const CbcModel &model)

{if (model.getNodeCount()==661) return;  printf("*** CHECKING tree after %d nodes\n",model.getNodeCount()) ;
 
  int j ;
  int nNodes = branchingTree->size() ;
# define MAXINFO 1000
  int *count = new int [MAXINFO] ;
  CbcNodeInfo **info = new CbcNodeInfo*[MAXINFO] ;
  int nInfo = 0 ;
/*
  Collect all CbcNodeInfo objects in info, by starting from each live node and
  traversing back to the root. Nodes in the live set should have unexplored
  branches remaining.

  TODO: The `while (nodeInfo)' loop could be made to break on reaching a
	common ancester (nodeInfo is found in info[k]). Alternatively, the
	check could change to signal an error if nodeInfo is not found above a
	common ancestor.
*/
  for (j = 0 ; j < nNodes ; j++)
  { CbcNode *node = branchingTree->nodePointer(j) ;
  if (!node)
    continue;
    CbcNodeInfo *nodeInfo = node->nodeInfo() ; 
    int change = node->nodeInfo()->numberBranchesLeft() ;
    assert(change) ;
    while (nodeInfo)
    { int k ;
      for (k = 0 ; k < nInfo ; k++)
      { if (nodeInfo == info[k]) break ; }
      if (k == nInfo)
      { assert(nInfo < MAXINFO) ;
	nInfo++ ;
	info[k] = nodeInfo ;
	count[k] = 0 ; }
      nodeInfo = nodeInfo->parent() ; } }
/*
  Walk the info array. For each nodeInfo, look up its parent in info and
  increment the corresponding count.
*/
  for (j = 0 ; j < nInfo ; j++)
  { CbcNodeInfo *nodeInfo = info[j] ;
    nodeInfo = nodeInfo->parent() ;
    if (nodeInfo)
    { int k ;
      for (k = 0 ; k < nInfo ; k++)
      { if (nodeInfo == info[k]) break ; }
      assert (k < nInfo) ;
      count[k]++ ; } }
/*
  Walk the info array one more time and check that the invariant holds. The
  number of references (numberPointingToThis()) should equal the sum of the
  number of actual references (held in count[]) plus the number of potential
  references (unexplored branches, numberBranchesLeft()).
*/
  for (j = 0;j < nInfo;j++) {
    CbcNodeInfo * nodeInfo = info[j] ;
    if (nodeInfo) {
      int k ;
      for (k = 0;k < nInfo;k++)
	if (nodeInfo == info[k])
	  break ;
      printf("Nodeinfo %x - %d left, %d count\n",
	     nodeInfo,
	     nodeInfo->numberBranchesLeft(),
	     nodeInfo->numberPointingToThis()) ;
      assert(nodeInfo->numberPointingToThis() ==
	     count[k]+nodeInfo->numberBranchesLeft()) ; } }

  delete [] count ;
  delete [] info ;
  
  return ; }

#endif	/* CHECK_NODE_FULL */



#ifdef CHECK_CUT_COUNTS

/*
  Routine to verify that cut reference counts are correct.
*/
void verifyCutCounts (const CbcTree * branchingTree, CbcModel &model)

{ printf("*** CHECKING cuts after %d nodes\n",model.getNodeCount()) ;

  int j ;
  int nNodes = branchingTree->size() ;

/*
  cut.tempNumber_ exists for the purpose of doing this verification. Clear it
  in all cuts. We traverse the tree by starting from each live node and working
  back to the root. At each CbcNodeInfo, check for cuts.
*/
  for (j = 0 ; j < nNodes ; j++)
  { CbcNode *node = branchingTree->nodePointer(j) ;
    CbcNodeInfo * nodeInfo = node->nodeInfo() ;
    assert (node->nodeInfo()->numberBranchesLeft()) ;
    while (nodeInfo)
    { int k ;
      for (k = 0 ; k < nodeInfo->numberCuts() ; k++)
      { CbcCountRowCut *cut = nodeInfo->cuts()[k] ;
	if (cut) cut->tempNumber_ = 0; }
      nodeInfo = nodeInfo->parent() ; } }
/*
  Walk the live set again, this time collecting the list of cuts in use at each
  node. addCuts1 will collect the cuts in model.addedCuts_. Take into account
  that when we recreate the basis for a node, we compress out the slack cuts.
*/
  for (j = 0 ; j < nNodes ; j++)
  { CoinWarmStartBasis *debugws = model.getEmptyBasis() ;
    CbcNode *node = branchingTree->nodePointer(j) ;
    CbcNodeInfo *nodeInfo = node->nodeInfo(); 
    int change = node->nodeInfo()->numberBranchesLeft() ;
    printf("Node %d %x (info %x) var %d way %d obj %g",j,node,
	   node->nodeInfo(),node->columnNumber(),node->way(),
	   node->objectiveValue()) ;

    model.addCuts1(node,debugws) ;

    int i ;
    int numberRowsAtContinuous = model.numberRowsAtContinuous() ;
    CbcCountRowCut **addedCuts = model.addedCuts() ;
    for (i = 0 ; i < model.currentNumberCuts() ; i++)
    { CoinWarmStartBasis::Status status = 
	debugws->getArtifStatus(i+numberRowsAtContinuous) ;
      if (status != CoinWarmStartBasis::basic && addedCuts[i])
      { addedCuts[i]->tempNumber_ += change ; } }

    while (nodeInfo)
    { nodeInfo = nodeInfo->parent() ;
      if (nodeInfo) printf(" -> %x",nodeInfo); }
    printf("\n") ;
    delete debugws ; }
/*
  The moment of truth: We've tallied up the references by direct scan of the
  search tree. Check for agreement with the count in the cut.

  TODO: Rewrite to check and print mismatch only when tempNumber_ == 0?
*/
  for (j = 0 ; j < nNodes ; j++)
  { CbcNode *node = branchingTree->nodePointer(j) ;
    CbcNodeInfo *nodeInfo = node->nodeInfo(); 
    while (nodeInfo)
    { int k ;
      for (k = 0 ; k < nodeInfo->numberCuts() ; k++)
      { CbcCountRowCut *cut = nodeInfo->cuts()[k] ;
	if (cut && cut->tempNumber_ >= 0)
	{ if (cut->tempNumber_ != cut->numberPointingToThis()) 
	    printf("mismatch %x %d %x %d %d\n",nodeInfo,k,
		    cut,cut->tempNumber_,cut->numberPointingToThis()) ;
	  else
	    printf("   match %x %d %x %d %d\n", nodeInfo,k,
		   cut,cut->tempNumber_,cut->numberPointingToThis()) ;
	  cut->tempNumber_ = -1 ; } }
      nodeInfo = nodeInfo->parent() ; } }

  return ; }

#endif /* CHECK_CUT_COUNTS */


//#define CHECK_CUT_SIZE
#ifdef CHECK_CUT_SIZE

/*
  Routine to verify that cut reference counts are correct.
*/
void verifyCutSize (const CbcTree * branchingTree, CbcModel &model)
{ 

  int j ;
  int nNodes = branchingTree->size() ;
  int totalCuts=0;

/*
  cut.tempNumber_ exists for the purpose of doing this verification. Clear it
  in all cuts. We traverse the tree by starting from each live node and working
  back to the root. At each CbcNodeInfo, check for cuts.
*/
  for (j = 0 ; j < nNodes ; j++) {
    CbcNode *node = branchingTree->nodePointer(j) ;
    CbcNodeInfo * nodeInfo = node->nodeInfo() ;
    assert (node->nodeInfo()->numberBranchesLeft()) ;
    while (nodeInfo) {
      totalCuts += nodeInfo->numberCuts();
      nodeInfo = nodeInfo->parent() ;
    }
  }
  printf("*** CHECKING cuts (size) after %d nodes - %d cuts\n",model.getNodeCount(),totalCuts) ;
  return ;
}

#endif /* CHECK_CUT_SIZE */

}

 /* End unnamed namespace for CbcModel.cpp */


static double trueIncrement=0.0;
void 
CbcModel::analyzeObjective ()
/*
  Try to find a minimum change in the objective function. The first scan
  checks that there are no continuous variables with non-zero coefficients,
  and grabs the largest objective coefficient associated with an unfixed
  integer variable. The second scan attempts to scale up the objective
  coefficients to a point where they are sufficiently close to integer that
  we can pretend they are integer, and calculate a gcd over the coefficients
  of interest. This will be the minimum increment for the scaled coefficients.
  The final action is to scale the increment back for the original coefficients
  and install it, if it's better than the existing value.

  John's note: We could do better than this.

  John's second note - apologies for changing s to z
*/
{ const double *objective = getObjCoefficients() ;
  const double *lower = getColLower() ;
  const double *upper = getColUpper() ;
/*
  Take a first scan to see if there are unfixed continuous variables in the
  objective.  If so, the minimum objective change could be arbitrarily small.
  Also pick off the maximum coefficient of an unfixed integer variable.

  If the objective is found to contain only integer variables, set the
  fathoming discipline to strict.
*/
  double maximumCost = 0.0 ;
  trueIncrement=0.0;
  bool possibleMultiple = true ;
  int iColumn ;
  int numberColumns = getNumCols() ;
  for (iColumn = 0 ; iColumn < numberColumns ; iColumn++)
  { if (upper[iColumn] > lower[iColumn]+1.0e-8)
    { if (isInteger(iColumn)) 
	maximumCost = CoinMax(maximumCost,fabs(objective[iColumn])) ;
      else if (objective[iColumn]) 
	possibleMultiple = false ; } }
  setIntParam(CbcModel::CbcFathomDiscipline,possibleMultiple) ;
/*
  If a nontrivial increment is possible, try and figure it out. We're looking
  for gcd(c<j>) for all c<j> that are coefficients of unfixed integer
  variables. Since the c<j> might not be integers, try and inflate them
  sufficiently that they look like integers (and we'll deflate the gcd
  later).

  2520.0 is used as it is a nice multiple of 2,3,5,7
*/
    if (possibleMultiple&&maximumCost)
    { int increment = 0 ;
      double multiplier = 2520.0 ;
      while (10.0*multiplier*maximumCost < 1.0e8)
	multiplier *= 10.0 ;
    int bigIntegers = 0; // Count of large costs which are integer
    for (iColumn = 0 ; iColumn < numberColumns ; iColumn++) {
      if (upper[iColumn] > lower[iColumn]+1.0e-8) {
	if (isInteger(iColumn)&&objective[iColumn]) {
	  double value = fabs(objective[iColumn])*multiplier ;
	  if (value <2.1e9) {
	    int nearest = (int) floor(value+0.5) ;
	    if (fabs(value-floor(value+0.5)) > 1.0e-8)
	      { increment = 0 ;
	      break ; }
	    else if (!increment)
	      { increment = nearest ; }
	    else
	      { increment = gcd(increment,nearest) ; }
	  } else {
	    // large value - may still be multiple of 1.0
	    value = fabs(objective[iColumn]);
	    if (fabs(value-floor(value+0.5)) > 1.0e-8) {
	      increment=0;
	      break;
	    } else {
	      bigIntegers++;
	    }
	  }
	}
      }
    }

/*
  If the increment beats the current value for objective change, install it.
*/
      if (increment)
      { double value = increment ;
	double cutoff = getDblParam(CbcModel::CbcCutoffIncrement) ;
	if (bigIntegers) {
	  // allow for 1.0
	  increment = gcd(increment,(int) multiplier);
	  value = increment;
	}
	value /= multiplier ;
	trueIncrement=CoinMax(cutoff,value);;
	if (value*0.999 > cutoff)
	{ messageHandler()->message(CBC_INTEGERINCREMENT,
					  messages())
	    << value << CoinMessageEol ;
	  setDblParam(CbcModel::CbcCutoffIncrement,value*0.999) ; } } }

  return ; 
}


/**
  \todo
  Normally, it looks like we enter here from command dispatch in the main
  routine, after calling the solver for an initial solution
  (CbcModel::initialSolve, which simply calls the solver's initialSolve
  routine.) The first thing we do is call resolve. Presumably there are
  circumstances where this is nontrivial? There's also a call from
  CbcModel::originalModel (tied up with integer presolve), which should be
  checked.

*/

/*
  The overall flow can be divided into three stages:
    * Prep: Check that the lp relaxation remains feasible at the root. If so,
      do all the setup for B&C.
    * Process the root node: Generate cuts, apply heuristics, and in general do
      the best we can to resolve the problem without B&C.
    * Do B&C search until we hit a limit or exhaust the search tree.
  
  Keep in mind that in general there is no node in the search tree that
  corresponds to the active subproblem. The active subproblem is represented
  by the current state of the model,  of the solver, and of the constraint
  system held by the solver.
*/
void CbcModel::branchAndBound(int doStatistics) 

{
/*
  Capture a time stamp before we start.
*/
  dblParam_[CbcStartSeconds] = CoinCpuTime();
  strongInfo_[0]=0;
  strongInfo_[1]=0;
  strongInfo_[2]=0;
  numberStrongIterations_ = 0;
#ifndef NDEBUG
  {
    int i;
    int n = solver_->getNumCols();
    const double *lower = solver_->getColLower() ;
    const double *upper = solver_->getColUpper() ;
    for (i=0;i<n;i++) {
      assert (lower[i]<1.0e10);
      assert (upper[i]>-1.0e10);
    }
    n = solver_->getNumRows();
    lower = solver_->getRowLower() ;
    upper = solver_->getRowUpper() ;
    for (i=0;i<n;i++) {
      assert (lower[i]<1.0e10);
      assert (upper[i]>-1.0e10);
    }
  }
#endif
  // original solver (only set if pre-processing)
  OsiSolverInterface * originalSolver=NULL;
  int numberOriginalObjects=numberObjects_;
  OsiObject ** originalObject = NULL;
  // Set up strategies
#if 0
  std::string problemName ;
  solver_->getStrParam(OsiProbName,problemName) ;
  if (!strcmp(problemName.c_str(),"PP08A")) solver_->activateRowCutDebugger(problemName.c_str()) ;
#endif
  if (strategy_) {
    // May do preprocessing
    originalSolver = solver_;
    strategy_->setupOther(*this);
    if (strategy_->preProcessState()) {
      // pre-processing done
      if (strategy_->preProcessState()<0) {
        // infeasible
        handler_->message(CBC_INFEAS,messages_)<< CoinMessageEol ;
        status_ = 0 ;
        secondaryStatus_ = 1;
        originalContinuousObjective_ = COIN_DBL_MAX;
        return ; 
      } else if (numberObjects_&&object_) {
        numberOriginalObjects=numberObjects_;
        // redo sequence
        numberIntegers_=0;
        int numberColumns = getNumCols();
        int nOrig = originalSolver->getNumCols();
        CglPreProcess * process = strategy_->process();
        assert (process);
        const int * originalColumns = process->originalColumns();
        // allow for cliques etc
        nOrig = CoinMax(nOrig,originalColumns[numberColumns-1]+1);
        // try and redo debugger
        OsiRowCutDebugger * debugger = const_cast<OsiRowCutDebugger *> (solver_->getRowCutDebuggerAlways());
        if (debugger)
          debugger->redoSolution(numberColumns,originalColumns);
        originalObject = object_;
        // object number or -1
        int * temp = new int[nOrig];
        int iColumn;
        for (iColumn=0;iColumn<nOrig;iColumn++) 
          temp[iColumn]=-1;
        int iObject;
        int nNonInt=0;
        for (iObject=0;iObject<numberOriginalObjects;iObject++) {
          iColumn = originalObject[iObject]->columnNumber();
          if (iColumn<0) {
            nNonInt++;
          } else {
            temp[iColumn]=iObject;
          }
        }
        int numberNewIntegers=0;
        int numberOldIntegers=0;
        int numberOldOther=0;
        for (iColumn=0;iColumn<numberColumns;iColumn++) {
          int jColumn = originalColumns[iColumn];
          if (temp[jColumn]>=0) {
            int iObject= temp[jColumn];
            CbcSimpleInteger * obj =
              dynamic_cast <CbcSimpleInteger *>(originalObject[iObject]) ;
            if (obj) 
              numberOldIntegers++;
            else
              numberOldOther++;
          } else if (isInteger(iColumn)) {
            numberNewIntegers++;
          }
        }
        /*
          Allocate an array to hold the indices of the integer variables.
          Make a large enough array for all objects
        */
        numberObjects_= numberNewIntegers+numberOldIntegers+numberOldOther+nNonInt;
        object_ = new OsiObject * [numberObjects_];
	delete [] integerVariable_;
        integerVariable_ = new int [numberNewIntegers+numberOldIntegers];
        /*
          Walk the variables again, filling in the indices and creating objects for
          the integer variables. Initially, the objects hold the index and upper &
          lower bounds.
        */
        numberIntegers_=0;
        for (iColumn=0;iColumn<numberColumns;iColumn++) {
          int jColumn = originalColumns[iColumn];
          if (temp[jColumn]>=0) {
            int iObject= temp[jColumn];
            CbcSimpleInteger * obj =
              dynamic_cast <CbcSimpleInteger *>(originalObject[iObject]) ;
            if (obj) {
              object_[numberIntegers_] = originalObject[iObject]->clone();
              // redo ids etc
              object_[numberIntegers_]->resetSequenceEtc(numberColumns,originalColumns);
              integerVariable_[numberIntegers_++]=iColumn;
            }
          } else if (isInteger(iColumn)) {
            object_[numberIntegers_] =
              new CbcSimpleInteger(this,iColumn);
            integerVariable_[numberIntegers_++]=iColumn;
          }
        }
        numberObjects_=numberIntegers_;
        // Now append other column stuff
        for (iColumn=0;iColumn<numberColumns;iColumn++) {
          int jColumn = originalColumns[iColumn];
          if (temp[jColumn]>=0) {
            int iObject= temp[jColumn];
            CbcSimpleInteger * obj =
              dynamic_cast <CbcSimpleInteger *>(originalObject[iObject]) ;
            if (!obj) {
              object_[numberObjects_] = originalObject[iObject]->clone();
              // redo ids etc
	      CbcObject * obj =
              dynamic_cast <CbcObject *>(object_[numberObjects_]) ;
	      assert (obj);
              obj->redoSequenceEtc(this,numberColumns,originalColumns);
              numberObjects_++;
            }
          }
        }
        // now append non column stuff
        for (iObject=0;iObject<numberOriginalObjects;iObject++) {
          iColumn = originalObject[iObject]->columnNumber();
          if (iColumn<0) {
            object_[numberObjects_] = originalObject[iObject]->clone();
            // redo ids etc
	    CbcObject * obj =
              dynamic_cast <CbcObject *>(object_[numberObjects_]) ;
	    assert (obj);
	    obj->redoSequenceEtc(this,numberColumns,originalColumns);
            numberObjects_++;
          }
        }
        delete [] temp;
        if (!numberObjects_)
          handler_->message(CBC_NOINT,messages_) << CoinMessageEol ;
      } else {
        int numberColumns = getNumCols();
        CglPreProcess * process = strategy_->process();
        assert (process);
        const int * originalColumns = process->originalColumns();
        // try and redo debugger
        OsiRowCutDebugger * debugger = const_cast<OsiRowCutDebugger *> (solver_->getRowCutDebuggerAlways());
        if (debugger)
          debugger->redoSolution(numberColumns,originalColumns);
      }
    } else {
      //no preprocessing
      originalSolver=NULL;
    }
    strategy_->setupCutGenerators(*this);
    strategy_->setupHeuristics(*this);
    // Set strategy print level to models
    strategy_->setupPrinting(*this,handler_->logLevel());
  }
  eventHappened_=false;
  CbcEventHandler *eventHandler = getEventHandler() ;
  // set up for probing
  probingInfo_ = new CglTreeProbingInfo(solver_);

  // Try for dominated columns
  if ((specialOptions_&64)!=0) {
    CglDuplicateRow dupcuts(solver_);
    dupcuts.setMode(2);
    CglStored * storedCuts = dupcuts.outDuplicates(solver_);
    addCutGenerator(storedCuts,1,"StoredCuts from dominated");
  }
  if (!nodeCompare_)
    nodeCompare_=new CbcCompareDefault();;
  // See if hot start wanted
  CbcCompareBase * saveCompare = NULL;
  if (hotstartSolution_) {
    if (strategy_&&strategy_->preProcessState()>0) {
      CglPreProcess * process = strategy_->process();
      assert (process);
      int n = solver_->getNumCols();
      const int * originalColumns = process->originalColumns();
      // columns should be in order ... but
      double * tempS = new double[n];
      for (int i=0;i<n;i++) {
        int iColumn = originalColumns[i];
        tempS[i]=hotstartSolution_[iColumn];
      }
      delete [] hotstartSolution_;
      hotstartSolution_=tempS;
      if (hotstartPriorities_) {
        int * tempP = new int [n];
        for (int i=0;i<n;i++) {
          int iColumn = originalColumns[i];
          tempP[i]=hotstartPriorities_[iColumn];
        }
        delete [] hotstartPriorities_;
        hotstartPriorities_=tempP;
      }
    }
    saveCompare = nodeCompare_;
    // depth first
    nodeCompare_ = new CbcCompareDepth();
  }
  if (!problemFeasibility_)
    problemFeasibility_=new CbcFeasibilityBase();
# ifdef CBC_DEBUG
  std::string problemName ;
  solver_->getStrParam(OsiProbName,problemName) ;
  printf("Problem name - %s\n",problemName.c_str()) ;
  solver_->setHintParam(OsiDoReducePrint,false,OsiHintDo,0) ;
# endif
/*
  Assume we're done, and see if we're proven wrong.
*/
  status_ = 0 ;
  secondaryStatus_ = 0;
  phase_=0;
/*
  Scan the variables, noting the integer variables. Create an
  CbcSimpleInteger object for each integer variable.
*/
  findIntegers(false) ;
  // If dynamic pseudo costs then do
  if (numberBeforeTrust_)
    convertToDynamic();
  // Set up char array to say if integer
  delete [] integerInfo_;
  {
    int n = solver_->getNumCols();
    integerInfo_ = new char [n];
    for (int i=0;i<n;i++) {
      if (solver_->isInteger(i))
        integerInfo_[i]=1;
      else
        integerInfo_[i]=0;
    }
  }
  if (preferredWay_) {
    // set all unset ones
    for (int iObject = 0 ; iObject < numberObjects_ ; iObject++) {
      CbcObject * obj =
	dynamic_cast <CbcObject *>(object_[iObject]) ;
      if (obj&&!obj->preferredWay())
	obj->setPreferredWay(preferredWay_);
    }
  }  
/*
  Ensure that objects on the lists of OsiObjects, heuristics, and cut
  generators attached to this model all refer to this model.
*/
  synchronizeModel() ;
  if (!solverCharacteristics_) {
    OsiBabSolver * solverCharacteristics = dynamic_cast<OsiBabSolver *> (solver_->getAuxiliaryInfo());
    if (solverCharacteristics) {
      solverCharacteristics_ = solverCharacteristics;
    } else {
      // replace in solver
      OsiBabSolver defaultC;
      solver_->setAuxiliaryInfo(&defaultC);
      solverCharacteristics_ = dynamic_cast<OsiBabSolver *> (solver_->getAuxiliaryInfo());
    }
  }
  solverCharacteristics_->setSolver(solver_);
  // Set so we can tell we are in initial phase in resolve
  continuousObjective_ = -COIN_DBL_MAX ;
/*
  Solve the relaxation.

  Apparently there are circumstances where this will be non-trivial --- i.e.,
  we've done something since initialSolve that's trashed the solution to the
  continuous relaxation.
*/
  bool feasible;
  // If NLP then we assume already solved outside branchAndbound
  if (!solverCharacteristics_->solverType()||solverCharacteristics_->solverType()==4) {
    feasible=resolve(NULL,0) != 0 ;
  } else {
    // pick up given status
    feasible = (solver_->isProvenOptimal() &&
                !solver_->isDualObjectiveLimitReached()) ;
  }
  if (problemFeasibility_->feasible(this,0)<0) {
    feasible=false; // pretend infeasible
  }
/*
  If the linear relaxation of the root is infeasible, bail out now. Otherwise,
  continue with processing the root node.
*/
  if (!feasible) {
    status_ = 0 ;
    if (!solver_->isProvenDualInfeasible()) {
      handler_->message(CBC_INFEAS,messages_)<< CoinMessageEol ;
      secondaryStatus_ = 1;
    } else {
      handler_->message(CBC_UNBOUNDED,messages_)<< CoinMessageEol ;
      secondaryStatus_ = 7;
    }
    originalContinuousObjective_ = COIN_DBL_MAX;
    solverCharacteristics_ = NULL;
    return ;
  }
  // Convert to Osi if wanted
  bool useOsiBranching=false;
  //OsiBranchingInformation * persistentInfo = NULL;
  if (branchingMethod_&&branchingMethod_->chooseMethod()) {
    useOsiBranching=true;
    //persistentInfo = new OsiBranchingInformation(solver_);
    if (numberOriginalObjects) {
      for (int iObject = 0 ; iObject < numberObjects_ ; iObject++) {
	CbcObject * obj =
	  dynamic_cast <CbcObject *>(object_[iObject]) ;
	if (obj) {
	  CbcSimpleInteger * obj2 =
	    dynamic_cast <CbcSimpleInteger *>(obj) ;
	  if (obj2) {
	    // back to Osi land
	    object_[iObject]=obj2->osiObject();
	    delete obj;
	  } else {
	    OsiSimpleInteger * obj3 =
	      dynamic_cast <OsiSimpleInteger *>(obj) ;
	    if (!obj3) {
	      OsiSOS * obj4 =
		dynamic_cast <OsiSOS *>(obj) ;
	      if (!obj4) {
		CbcSOS * obj5 =
		  dynamic_cast <CbcSOS *>(obj) ;
		if (obj5) {
		  // back to Osi land
		  object_[iObject]=obj5->osiObject(solver_);
		} else {
		  printf("Code up CbcObject type in Osi land\n");
		  abort();
		}
	      }
	    }
	  }
	}
      }
      // and add to solver if none
      if (!solver_->numberObjects()) {
	solver_->addObjects(numberObjects_,object_);
      } else {
	if (solver_->numberObjects()!=numberOriginalObjects) {
	  printf("should have trapped that solver has objects before\n");
	  abort();
	}
      }
    } else {
      // do from solver
      deleteObjects(false);
      solver_->findIntegersAndSOS(false);
      numberObjects_=solver_->numberObjects();
      object_ = solver_->objects();
      ownObjects_ = false;
    }
    branchingMethod_->chooseMethod()->setSolver(solver_);
  }
  // take off heuristics if have to
  if (numberHeuristics_) {
    int numberOdd=0;
    for (int i=0;i<numberObjects_;i++) {
      if (!object_[i]->canDoHeuristics()) 
	numberOdd++;
    }
    if (numberOdd) {
      int k=0;
      for (int i=0;i<numberHeuristics_;i++) {
	if (!heuristic_[i]->canDealWithOdd())
	  delete heuristic_[i];
	else
	  heuristic_[k++]=heuristic_[i];
      }
      if (!k) {
	delete [] heuristic_;
	heuristic_=NULL;
      }
      numberHeuristics_=k;
      handler_->message(CBC_HEURISTICS_OFF,messages_)<< numberOdd<<CoinMessageEol ;
    }
  }
  // Save objective (just so user can access it)
  originalContinuousObjective_ = solver_->getObjValue();
  bestPossibleObjective_=originalContinuousObjective_;
  sumChangeObjective1_=0.0;
  sumChangeObjective2_=0.0;
/*
  OsiRowCutDebugger knows an optimal answer for a subset of MIP problems.
  Assuming it recognises the problem, when called upon it will check a cut to
  see if it cuts off the optimal answer.
*/
  // If debugger exists set specialOptions_ bit
  if (solver_->getRowCutDebuggerAlways())
    specialOptions_ |= 1;

# ifdef CBC_DEBUG
  if ((specialOptions_&1)==0)
    solver_->activateRowCutDebugger(problemName.c_str()) ;
  if (solver_->getRowCutDebuggerAlways())
    specialOptions_ |= 1;
# endif

/*
  Begin setup to process a feasible root node.
*/
  bestObjective_ = CoinMin(bestObjective_,1.0e50) ;
  numberSolutions_ = 0 ;
  stateOfSearch_=0;
  numberHeuristicSolutions_ = 0 ;
  // Everything is minimization
  { 
    // needed to sync cutoffs
    double value ;
    solver_->getDblParam(OsiDualObjectiveLimit,value) ;
    dblParam_[CbcCurrentCutoff]= value * solver_->getObjSense();
  }
  double cutoff=getCutoff() ;
  double direction = solver_->getObjSense() ;
  dblParam_[CbcOptimizationDirection]=direction;
  if (cutoff < 1.0e20&&direction<0.0)
    messageHandler()->message(CBC_CUTOFF_WARNING1,
				    messages())
				      << cutoff << -cutoff << CoinMessageEol ;
  if (cutoff > bestObjective_)
    cutoff = bestObjective_ ;
  setCutoff(cutoff) ;
/*
  We probably already have a current solution, but just in case ...
*/
  int numberColumns = getNumCols() ;
  if (!currentSolution_)
    currentSolution_ = new double[numberColumns] ;
  testSolution_ = currentSolution_;
  /* Tell solver we are in Branch and Cut
     Could use last parameter for subtle differences */
  solver_->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,NULL) ;
#ifdef COIN_HAS_CLP
 {
   OsiClpSolverInterface * clpSolver 
     = dynamic_cast<OsiClpSolverInterface *> (solver_);
   if (clpSolver) {
     //#define CLP_QUICK_OPTIONS
#ifdef CLP_QUICK_OPTIONS
     // Try and re-use regions   
     ClpSimplex * simplex = clpSolver->getModelPtr();
     simplex->setPersistenceFlag(1);
     clpSolver->deleteScaleFactors();
     int value=131072;
     clpSolver->setSpecialOptions(clpSolver->specialOptions()|value);
     if ((clpSolver->specialOptions()&value)!=0)
       simplex->setSpecialOptions(simplex->specialOptions()|value);
     //if (simplex->numberRows()<50)
     //simplex->setAlphaAccuracy(1.0);
     //clpSolver->setSpecialOptions((clpSolver->specialOptions()&~128)|65536);
#endif
     if ((specialOptions_&32)==0) {
       ClpSimplex * clpSimplex = clpSolver->getModelPtr();
       // take off names
       clpSimplex->dropNames();
     }
     // no crunch if mostly continuous
     int numberColumns = solver_->getNumCols();
     if (numberColumns>1000&&numberIntegers_*4<numberColumns)
       clpSolver->setSpecialOptions(clpSolver->specialOptions()&(~1));
   }
 }
#endif
/*
  Create a copy of the solver, thus capturing the original (root node)
  constraint system (aka the continuous system).
*/
  continuousSolver_ = solver_->clone() ;

  numberRowsAtContinuous_ = getNumRows() ;
/*
  Check the objective to see if we can deduce a nontrivial increment. If
  it's better than the current value for CbcCutoffIncrement, it'll be
  installed.
*/
  if(solverCharacteristics_->reducedCostsAccurate())
    analyzeObjective() ;
/*
  Set up for cut generation. addedCuts_ holds the cuts which are relevant for
  the active subproblem. whichGenerator will be used to record the generator
  that produced a given cut.
*/
  maximumWhich_ = 1000 ;
  delete [] whichGenerator_;
  whichGenerator_ = new int[maximumWhich_] ;
  memset(whichGenerator_,0,maximumWhich_*sizeof(int));
  maximumNumberCuts_ = 0 ;
  currentNumberCuts_ = 0 ;
  delete [] addedCuts_ ;
  addedCuts_ = NULL ;
/*
  Set up an empty heap and associated data structures to hold the live set
  (problems which require further exploration).
*/
  tree_->setComparison(*nodeCompare_) ;
/*
  Used to record the path from a node to the root of the search tree, so that
  we can then traverse from the root to the node when restoring a subproblem.
*/
  maximumDepth_ = 10 ;
  delete [] walkback_ ;
  walkback_ = new CbcNodeInfo * [maximumDepth_] ;
/*
  Used to generate bound edits for CbcPartialNodeInfo.
*/
  double * lowerBefore = new double [numberColumns] ;
  double * upperBefore = new double [numberColumns] ;
/*
  
  Generate cuts at the root node and reoptimise. solveWithCuts does the heavy
  lifting. It will iterate a generate/reoptimise loop (including reduced cost
  fixing) until no cuts are generated, the change in objective falls off,  or
  the limit on the number of rounds of cut generation is exceeded.

  At the end of all this, any cuts will be recorded in cuts and also
  installed in the solver's constraint system. We'll have reoptimised, and
  removed any slack cuts (numberOldActiveCuts_ and numberNewCuts_ have been
  adjusted accordingly).

  Tell cut generators they can be a bit more aggressive at root node

  TODO: Why don't we make a copy of the solution after solveWithCuts?
  TODO: If numberUnsatisfied == 0, don't we have a solution?
*/
  phase_=1;
  int iCutGenerator;
  for (iCutGenerator = 0;iCutGenerator<numberCutGenerators_;iCutGenerator++) {
    CglCutGenerator * generator = generator_[iCutGenerator]->generator();
    generator->setAggressiveness(generator->getAggressiveness()+100);
  }
  OsiCuts cuts ;
  int anyAction = -1 ;
  numberOldActiveCuts_ = 0 ;
  numberNewCuts_ = 0 ;
  // Array to mark solution
  delete [] usedInSolution_;
  usedInSolution_ = new int[numberColumns];
  CoinZeroN(usedInSolution_,numberColumns);
/*
  For printing totals and for CbcNode (numberNodes_)
*/
  numberIterations_ = 0 ;
  numberNodes_ = 0 ;
  numberNodes2_ = 0 ;
  maximumStatistics_=0;
  maximumDepthActual_=0;
  numberDJFixed_=0.0;
  // Do heuristics
  {
    double * newSolution = new double [numberColumns] ;
    double heuristicValue = getCutoff() ;
    int found = -1; // no solution found

    currentPassNumber_ = 1; // so root heuristics will run
    int i;
    for (i = 0;i<numberHeuristics_;i++) {
      // see if heuristic will do anything
      double saveValue = heuristicValue ;
      int ifSol = heuristic_[i]->solution(heuristicValue,
					  newSolution);
      if (ifSol>0) {
	// better solution found
	found = i ;
	incrementUsed(newSolution);
	// increment number of solutions so other heuristics can test
	numberSolutions_++;
	numberHeuristicSolutions_++;
      } else {
	heuristicValue = saveValue ;
      }
    }
    currentPassNumber_ = 0;
    /*
      Did any of the heuristics turn up a new solution? Record it before we free
      the vector.
    */
    if (found >= 0) { 
      // For compiler error on terra cluster!
      if (found<numberHeuristics_)
	lastHeuristic_ = heuristic_[found];
      else
	lastHeuristic_ = heuristic_[0];
      setBestSolution(CBC_ROUNDING,heuristicValue,newSolution) ;
      CbcTreeLocal * tree 
          = dynamic_cast<CbcTreeLocal *> (tree_);
      if (tree)
        tree->passInSolution(bestSolution_,heuristicValue);
      for (i = 0;i<numberHeuristics_;i++) {
	// delete FPump
	CbcHeuristicFPump * pump 
	  = dynamic_cast<CbcHeuristicFPump *> (heuristic_[i]);
	if (pump) {
	  delete pump;
	  numberHeuristics_ --;
	  for (int j=i;j<numberHeuristics_;j++)
	    heuristic_[j] = heuristic_[j+1];
	}
      }
    }
    delete [] newSolution ;
  } 
  statistics_ = NULL;
  // Do on switch
  if (doStatistics>0&&doStatistics<100) {
    maximumStatistics_=10000;
    statistics_ = new CbcStatistics * [maximumStatistics_];
    memset(statistics_,0,maximumStatistics_*sizeof(CbcStatistics *));
  }

  { int iObject ;
    int preferredWay ;
    int numberUnsatisfied = 0 ;
    memcpy(currentSolution_,solver_->getColSolution(),
	   numberColumns*sizeof(double)) ;
    // point to useful information
    OsiBranchingInformation usefulInfo=usefulInformation();

    for (iObject = 0 ; iObject < numberObjects_ ; iObject++)
    { double infeasibility =
	  object_[iObject]->infeasibility(&usefulInfo,preferredWay) ;
      if (infeasibility ) numberUnsatisfied++ ; }
    // replace solverType 
    if(solverCharacteristics_->tryCuts())  {

      if (numberUnsatisfied)   {
        feasible = solveWithCuts(cuts,maximumCutPassesAtRoot_,
                                 NULL);
      }	else if (solverCharacteristics_->solutionAddsCuts()||
                 solverCharacteristics_->alwaysTryCutsAtRootNode()) {
        // may generate cuts and turn the solution
        //to an infeasible one
        feasible = solveWithCuts(cuts, 1,
                                 NULL);
      }
    }
    // check extra info on feasibility
    if (!solverCharacteristics_->mipFeasible())
      feasible = false;
  }
  // make cut generators less aggressive
  for (iCutGenerator = 0;iCutGenerator<numberCutGenerators_;iCutGenerator++) {
    CglCutGenerator * generator = generator_[iCutGenerator]->generator();
    generator->setAggressiveness(generator->getAggressiveness()-100);
  }
  currentNumberCuts_ = numberNewCuts_ ;
  // See if can stop on gap
  stoppedOnGap_ = false ;
  bestPossibleObjective_ = solver_->getObjValue()*solver_->getObjSense();
  double testGap = CoinMax(dblParam_[CbcAllowableGap],
			   CoinMax(fabs(bestObjective_),fabs(bestPossibleObjective_))
			   *dblParam_[CbcAllowableFractionGap]);
  if (bestObjective_-bestPossibleObjective_ < testGap && getCutoffIncrement()>=0.0) {
    if (bestPossibleObjective_<getCutoff())
      stoppedOnGap_ = true ;
    feasible = false;
  }
/*
  We've taken the continuous relaxation as far as we can. Time to branch.
  The first order of business is to actually create a node. chooseBranch
  currently uses strong branching to evaluate branch object candidates,
  unless forced back to simple branching. If chooseBranch concludes that a
  branching candidate is monotone (anyAction == -1) or infeasible (anyAction
  == -2) when forced to integer values, it returns here immediately.

  Monotone variables trigger a call to resolve(). If the problem remains
  feasible, try again to choose a branching variable. At the end of the loop,
  resolved == true indicates that some variables were fixed.

  Loss of feasibility will result in the deletion of newNode.
*/

  bool resolved = false ;
  CbcNode *newNode = NULL ;
  numberFixedAtRoot_=0;
  numberFixedNow_=0;
  int numberIterationsAtContinuous = numberIterations_;
  //solverCharacteristics_->setSolver(solver_);
  if (feasible) {
    newNode = new CbcNode ;
    // Set objective value (not so obvious if NLP etc)
    setObjectiveValue(newNode,NULL);
    anyAction = -1 ;
    // To make depth available we may need a fake node
    CbcNode fakeNode;
    if (!currentNode_) {
      // Not true if sub trees assert (!numberNodes_);
      currentNode_=&fakeNode;
    }
    phase_=3;
    // only allow 1000 passes
    int numberPassesLeft=1000;
    // This is first crude step
    if (numberAnalyzeIterations_) {
      delete [] analyzeResults_;
      analyzeResults_ = new double [4*numberIntegers_];
      numberFixedAtRoot_=newNode->analyze(this,analyzeResults_);
      if (numberFixedAtRoot_>0) {
        printf("%d fixed by analysis\n",numberFixedAtRoot_);
        setPointers(solver_);
        numberFixedNow_ = numberFixedAtRoot_;
      } else if (numberFixedAtRoot_<0) {
        printf("analysis found to be infeasible\n");
        anyAction=-2;
        delete newNode ;
	newNode = NULL ;
	feasible = false ;
      }
    }
    OsiSolverBranch * branches = NULL;
    anyAction = chooseBranch(newNode, numberPassesLeft, NULL, cuts,resolved,
			     NULL,NULL,NULL,branches);
    if (anyAction == -2||newNode->objectiveValue() >= cutoff) {
      if (anyAction != -2) {
	// zap parent nodeInfo
#ifdef COIN_DEVELOP
	printf("zapping CbcNodeInfo %x\n",newNode->nodeInfo()->parent());
#endif
	if (newNode->nodeInfo())
	  newNode->nodeInfo()->nullParent();
      }
      delete newNode ;
      newNode = NULL ;
      feasible = false ;
    }
  }
/*
  At this point, the root subproblem is infeasible or fathomed by bound
  (newNode == NULL), or we're live with an objective value that satisfies the
  current objective cutoff.
*/
  assert (!newNode || newNode->objectiveValue() <= cutoff) ;
  // Save address of root node as we don't want to delete it
  // initialize for print out
  int lastDepth=0;
  int lastUnsatisfied=0;
  if (newNode)
    lastUnsatisfied=newNode->numberUnsatisfied();
/*
  The common case is that the lp relaxation is feasible but doesn't satisfy
  integrality (i.e., newNode->branchingObject(), indicating we've been able to
  select a branching variable). Remove any cuts that have gone slack due to
  forcing monotone variables. Then tack on an CbcFullNodeInfo object and full
  basis (via createInfo()) and stash the new cuts in the nodeInfo (via
  addCuts()). If, by some miracle, we have an integral solution at the root
  (newNode->branchingObject() is NULL), takeOffCuts() will ensure that the solver holds
  a valid solution for use by setBestSolution().
*/
  CoinWarmStartBasis *lastws = NULL ;
  if (feasible && newNode->branchingObject())
  { if (resolved)
    { takeOffCuts(cuts,false,NULL) ;
#     ifdef CHECK_CUT_COUNTS
      { printf("Number of rows after chooseBranch fix (root)"
	       "(active only) %d\n",
		numberRowsAtContinuous_+numberNewCuts_+numberOldActiveCuts_) ;
	const CoinWarmStartBasis* debugws =
	  dynamic_cast <const CoinWarmStartBasis*>(solver_->getWarmStart()) ;
	debugws->print() ;
	delete debugws ; }
#     endif
    }
  //newNode->createInfo(this,NULL,NULL,NULL,NULL,0,0) ;
    newNode->nodeInfo()->addCuts(cuts,
				 newNode->numberBranches(),whichGenerator_) ;
    if (lastws) delete lastws ;
    lastws = dynamic_cast<CoinWarmStartBasis*>(solver_->getWarmStart()) ;
  }
/*
  Continuous data to be used later
*/
  continuousObjective_ = solver_->getObjValue()*solver_->getObjSense();
  continuousInfeasibilities_ = 0 ;
  if (newNode)
  { continuousObjective_ = newNode->objectiveValue() ;
    delete [] continuousSolution_;
    continuousSolution_ = CoinCopyOfArray(solver_->getColSolution(),
                                             numberColumns);
    continuousInfeasibilities_ = newNode->numberUnsatisfied() ; }
/*
  Bound may have changed so reset in objects
*/
  { int i ;
    for (i = 0;i < numberObjects_;i++)
      object_[i]->resetBounds(solver_) ; }
/*
  Feasible? Then we should have either a live node prepped for future
  expansion (indicated by variable() >= 0), or (miracle of miracles) an
  integral solution at the root node.

  initializeInfo sets the reference counts in the nodeInfo object.  Since
  this node is still live, push it onto the heap that holds the live set.
*/
  double bestValue = 0.0 ;
  if (newNode) {
    bestValue = newNode->objectiveValue();
    if (newNode->branchingObject()) {
      newNode->initializeInfo() ;
      tree_->push(newNode) ;
      if (statistics_) {
        if (numberNodes2_==maximumStatistics_) {
          maximumStatistics_ = 2*maximumStatistics_;
          CbcStatistics ** temp = new CbcStatistics * [maximumStatistics_];
          memset(temp,0,maximumStatistics_*sizeof(CbcStatistics *));
          memcpy(temp,statistics_,numberNodes2_*sizeof(CbcStatistics *));
          delete [] statistics_;
          statistics_=temp;
        }
        assert (!statistics_[numberNodes2_]);
        statistics_[numberNodes2_]=new CbcStatistics(newNode);
      }
      numberNodes2_++;
#     ifdef CHECK_NODE
      printf("Node %x on tree\n",newNode) ;
#     endif
    } else {
      // continuous is integer
      double objectiveValue = newNode->objectiveValue();
      setBestSolution(CBC_SOLUTION,objectiveValue,
		      solver_->getColSolution()) ;
      delete newNode ;
      newNode = NULL ;
    }
  }

  if (printFrequency_ <= 0) {
    printFrequency_ = 1000 ;
    if (getNumCols() > 2000)
      printFrequency_ = 100 ;
  }
  /*
    It is possible that strong branching fixes one variable and then the code goes round
    again and again.  This can take too long.  So we need to warn user - just once.
  */
  numberLongStrong_=0;
  double totalTime = 0.0;
#ifdef CBC_THREAD
  CbcNode * createdNode=NULL;
  CbcModel ** threadModel = NULL;
  pthread_t * threadId = NULL;
  int * threadCount = NULL;
  pthread_mutex_t mutex;
  pthread_cond_t condition_main;
  pthread_mutex_t condition_mutex;
  pthread_mutex_t * mutex2 = NULL;
  pthread_cond_t * condition2 = NULL;
  threadStruct * threadInfo = NULL;
  bool locked=false;
  int threadStats[6];
  memset(threadStats,0,sizeof(threadStats));
  double timeWaiting=0.0;
  // For now just one model
  if (numberThreads_) {
    threadId = new pthread_t [numberThreads_];
    threadCount = new int [numberThreads_];
    CoinZeroN(threadCount,numberThreads_);
    pthread_mutex_init(&mutex,NULL);
    pthread_cond_init(&condition_main,NULL);
    pthread_mutex_init(&condition_mutex,NULL);
    threadModel = new CbcModel * [numberThreads_+1];
    threadInfo = new threadStruct [numberThreads_+1];
    mutex2 = new pthread_mutex_t [numberThreads_];
    condition2 = new pthread_cond_t [numberThreads_];
    // we don't want a strategy object
    CbcStrategy * saveStrategy = strategy_;
    strategy_ = NULL;
    for (int i=0;i<numberThreads_;i++) {
      pthread_mutex_init(mutex2+i,NULL);
      pthread_cond_init(condition2+i,NULL);
      threadId[i]=0;
      threadInfo[i].baseModel=this;
      threadModel[i]=new CbcModel(*this);
#ifdef COIN_HAS_CLP
      // Solver may need to know about model
      CbcModel * thisModel = threadModel[i];
      CbcOsiSolver * solver =
              dynamic_cast<CbcOsiSolver *>(thisModel->solver()) ;
      if (solver)
	solver->setCbcModel(thisModel);
#endif
      mutex_ = (void *) (threadInfo+i);
      threadModel[i]->moveToModel(this,-1);
      threadInfo[i].thisModel=threadModel[i];
      threadInfo[i].node=NULL;
      threadInfo[i].createdNode=NULL;
      threadInfo[i].threadIdOfBase=pthread_self();
      threadInfo[i].mutex=&mutex;
      threadInfo[i].mutex2=mutex2+i;
      threadInfo[i].condition2=condition2+i;
      threadInfo[i].returnCode=-1;
      threadInfo[i].timeLocked=0.0;
      threadInfo[i].timeWaitingToLock=0.0;
      threadInfo[i].timeWaitingToStart=0.0;
      threadInfo[i].timeInThread=0.0;
      threadInfo[i].numberTimesLocked=0;
      threadInfo[i].numberTimesUnlocked=0;
      threadInfo[i].numberTimesWaitingToStart=0;
      threadInfo[i].locked=false;
#if CBC_THREAD_DEBUG
      threadInfo[i].threadNumber=i+2;
#endif
      pthread_create(threadId+i,NULL,doNodesThread,threadInfo+i);
    }
    strategy_ = saveStrategy;
    // Do a partial one for base model
    threadInfo[numberThreads_].baseModel=this;
    threadModel[numberThreads_]=this;
    mutex_ = (void *) (threadInfo+numberThreads_);
    threadInfo[numberThreads_].node=NULL;
    threadInfo[numberThreads_].mutex=&mutex;
    threadInfo[numberThreads_].condition2=&condition_main;
    threadInfo[numberThreads_].mutex2=&condition_mutex;
    threadInfo[numberThreads_].timeLocked=0.0;
    threadInfo[numberThreads_].timeWaitingToLock=0.0;
    threadInfo[numberThreads_].numberTimesLocked=0;
    threadInfo[numberThreads_].numberTimesUnlocked=0;
    threadInfo[numberThreads_].locked=false;
#if CBC_THREAD_DEBUG
    threadInfo[numberThreads_].threadNumber=1;
#endif
  }
#endif
/*
  At last, the actual branch-and-cut search loop, which will iterate until
  the live set is empty or we hit some limit (integrality gap, time, node
  count, etc.). The overall flow is to rebuild a subproblem, reoptimise using
  solveWithCuts(), choose a branching pattern with chooseBranch(), and finally
  add the node to the live set.

  The first action is to winnow the live set to remove nodes which are worse
  than the current objective cutoff.
*/
  while (true) {
#ifdef CBC_THREAD
    if (!locked) {
      lockThread();
      locked=true;
    }
#endif
    if (tree_->empty()) {
#ifdef CBC_THREAD
      if (numberThreads_) {
#ifdef COIN_DEVELOP
	printf("empty\n");
#endif
	// may still be outstanding nodes
	int iThread;
	for (iThread=0;iThread<numberThreads_;iThread++) {
	  if (threadId[iThread]) {
	    if (threadInfo[iThread].returnCode==0) 
	      break;
	  }
	}
	if (iThread<numberThreads_) {
#ifdef COIN_DEVELOP
	  printf("waiting for thread %d code 0\n",iThread);
#endif
	  unlockThread();
	  locked = false;
	  pthread_cond_signal(threadInfo[iThread].condition2); // unlock in case
	  while (true) {
	    pthread_mutex_lock(&condition_mutex);
	    struct timespec absTime;
	    clock_gettime(CLOCK_REALTIME,&absTime);
	    double time = absTime.tv_sec+1.0e-9*absTime.tv_nsec;
	    absTime.tv_nsec += 1000000; // millisecond
	    if (absTime.tv_nsec>=1000000000) {
	      absTime.tv_nsec -= 1000000000;
	      absTime.tv_sec++;
	    }
	    pthread_cond_timedwait(&condition_main,&condition_mutex,&absTime);
	    clock_gettime(CLOCK_REALTIME,&absTime);
	    double time2 = absTime.tv_sec+1.0e-9*absTime.tv_nsec;
	    timeWaiting += time2-time;
	    pthread_mutex_unlock(&condition_mutex);
	    if (threadInfo[iThread].returnCode!=0) 
	      break;
	    pthread_cond_signal(threadInfo[iThread].condition2); // unlock
	  }
	  threadModel[iThread]->moveToModel(this,1);
	  assert (threadInfo[iThread].returnCode==1);
	  // say available
	  threadInfo[iThread].returnCode=-1;
	  threadStats[4]++;
#ifdef COIN_DEVELOP
	  printf("thread %d code now -1\n",iThread);
#endif
	  continue;
	} else {
#ifdef COIN_DEVELOP
	  printf("no threads at code 0 \n");
#endif
	  // now check if any have just finished
	  for (iThread=0;iThread<numberThreads_;iThread++) {
	    if (threadId[iThread]) {
	      if (threadInfo[iThread].returnCode==1) 
		break;
	    }
	  }
	  if (iThread<numberThreads_) {
	    unlockThread();
	    locked = false;
	    threadModel[iThread]->moveToModel(this,1);
	    assert (threadInfo[iThread].returnCode==1);
	    // say available
	    threadInfo[iThread].returnCode=-1;
	    threadStats[4]++;
#ifdef COIN_DEVELOP
	    printf("thread %d code now -1\n",iThread);
#endif
	    continue;
	  }
	}
	if (!tree_->empty()) {
#ifdef COIN_DEVELOP
	  printf("tree not empty!!!!!!\n");
#endif
	  continue;
	}
	for (iThread=0;iThread<numberThreads_;iThread++) {
	  if (threadId[iThread]) {
	    if (threadInfo[iThread].returnCode!=-1) { 
	      printf("bad end of tree\n");
	      abort();
	    }
	  }
	}
#ifdef COIN_DEVELOP
	printf("finished ************\n");
#endif
      }
      unlockThread();
      locked=false; // not needed as break
#endif
      break;
    }
#ifdef CBC_THREAD
    unlockThread();
    locked = false;
#endif
/*
  Check for abort on limits: node count, solution count, time, integrality gap.
*/
    totalTime = getCurrentSeconds() ;
    if (!(numberNodes_ < intParam_[CbcMaxNumNode] &&
	  numberSolutions_ < intParam_[CbcMaxNumSol] &&
	  totalTime < dblParam_[CbcMaximumSeconds] &&
	  !stoppedOnGap_&&!eventHappened_)) {
      // out of loop
      break;
    }
#ifdef BONMIN
    assert(!solverCharacteristics_->solutionAddsCuts() || solverCharacteristics_->mipFeasible());
#endif
    if (cutoff > getCutoff()) {
      double newCutoff = getCutoff();
      if (analyzeResults_) {
	// see if we could fix any (more)
	int n=0;
	double * newLower = analyzeResults_;
	double * objLower = newLower+numberIntegers_;
	double * newUpper = objLower+numberIntegers_;
	double * objUpper = newUpper+numberIntegers_;
	for (int i=0;i<numberIntegers_;i++) {
	  if (objLower[i]>newCutoff) {
	    n++;
	    if (objUpper[i]>newCutoff) {
	      newCutoff = -COIN_DBL_MAX;
	      break;
	    }
	  } else if (objUpper[i]>newCutoff) {
	    n++;
	  }
	}
	if (newCutoff==-COIN_DBL_MAX) {
	  printf("Root analysis says finished\n");
	} else if (n>numberFixedNow_) {
	  printf("%d more fixed by analysis - now %d\n",n-numberFixedNow_,n);
	  numberFixedNow_=n;
	}
      }
      if (eventHandler) {
	if (!eventHandler->event(CbcEventHandler::solution)) {
	  eventHappened_=true; // exit
	}
      }
      lockThread();
      // Do from deepest
      tree_->cleanTree(this, newCutoff,bestPossibleObjective_) ;
      nodeCompare_->newSolution(this) ;
      nodeCompare_->newSolution(this,continuousObjective_,
				continuousInfeasibilities_) ;
      tree_->setComparison(*nodeCompare_) ;
      if (tree_->empty()) {
	unlockThread();
	// For threads we need to check further
	//break; // finished
	continue;
      }
      unlockThread();
    }
    cutoff = getCutoff() ;
/*
    Periodic activities: Opportunities to
    + tweak the nodeCompare criteria,
    + check if we've closed the integrality gap enough to quit, 
    + print a summary line to let the user know we're working
*/
    if ((numberNodes_%1000) == 0) {
      lockThread();
      bool redoTree=nodeCompare_->every1000Nodes(this, numberNodes_) ;
#ifdef CHECK_CUT_SIZE
      verifyCutSize (tree_, *this);
#endif
      // redo tree if wanted
      if (redoTree)
	tree_->setComparison(*nodeCompare_) ;
      unlockThread();
    }
    if (saveCompare&&!hotstartSolution_) {
      // hotstart switched off
      delete nodeCompare_; // off depth first
      nodeCompare_=saveCompare;
      saveCompare=NULL;
      // redo tree
      lockThread();
      tree_->setComparison(*nodeCompare_) ;
      unlockThread();
    }
    if ((numberNodes_%printFrequency_) == 0) {
      lockThread();
      int j ;
      int nNodes = tree_->size() ;
      bestPossibleObjective_ = 1.0e100 ;
      for (j = 0;j < nNodes;j++) {
	CbcNode * node = tree_->nodePointer(j) ;
	if (node&&node->objectiveValue() < bestPossibleObjective_)
	  bestPossibleObjective_ = node->objectiveValue() ;
      }
      unlockThread();
      if (!intParam_[CbcPrinting]) {
	messageHandler()->message(CBC_STATUS,messages())
	  << numberNodes_<< nNodes<< bestObjective_<< bestPossibleObjective_
	  <<getCurrentSeconds()
	  << CoinMessageEol ;
      } else {
	messageHandler()->message(CBC_STATUS2,messages())
	  << numberNodes_<< nNodes<< bestObjective_<< bestPossibleObjective_
	  <<lastDepth<<lastUnsatisfied<<numberIterations_
	  <<getCurrentSeconds()
	  << CoinMessageEol ;
      }
      if (!eventHandler->event(CbcEventHandler::treeStatus)) {
	eventHappened_=true; // exit
      }
    }
    // If no solution but many nodes - signal change in strategy
    if (numberNodes_>2*numberObjects_+1000&&stateOfSearch_!=2)
      stateOfSearch_=3;
    // See if can stop on gap
    double testGap = CoinMax(dblParam_[CbcAllowableGap],
			     CoinMax(fabs(bestObjective_),fabs(bestPossibleObjective_))
			     *dblParam_[CbcAllowableFractionGap]);
    if (bestObjective_-bestPossibleObjective_ < testGap && getCutoffIncrement()>=0.0) {
      stoppedOnGap_ = true ;
    }
    
#   ifdef CHECK_NODE_FULL
    verifyTreeNodes(tree_,*this) ;
#   endif
#   ifdef CHECK_CUT_COUNTS
    verifyCutCounts(tree_,*this) ;
#   endif
/*
  Now we come to the meat of the loop. To create the active subproblem, we'll
  pop the most promising node in the live set, rebuild the subproblem it
  represents, and then execute the current arm of the branch to create the
  active subproblem.
*/
    CbcNode *node = tree_->bestNode(cutoff) ;
    // Possible one on tree worse than cutoff
    if (!node)
      continue;
#ifndef CBC_THREAD
    int currentNumberCuts = 0 ;
    currentNode_=node; // so can be accessed elsewhere
#ifdef CBC_DEBUG
    printf("%d unsat, way %d, obj %g est %g\n",
	   node->numberUnsatisfied(),node->way(),node->objectiveValue(),
	   node->guessedObjectiveValue());
#endif
#if NEW_UPDATE_OBJECT==0
    // Save clone in branching decision
    if(branchingMethod_)
      branchingMethod_->saveBranchingObject(node->modifiableBranchingObject());
#endif
    // Say not on optimal path
    bool onOptimalPath=false;
#   ifdef CHECK_NODE
    printf("Node %x popped from tree - %d left, %d count\n",node,
	   node->nodeInfo()->numberBranchesLeft(),
	   node->nodeInfo()->numberPointingToThis()) ;
    printf("\tdepth = %d, z =  %g, unsat = %d, var = %d.\n",
	   node->depth(),node->objectiveValue(),
	   node->numberUnsatisfied(),
	   node->columnNumber()) ;
#   endif
    lastDepth=node->depth();
    lastUnsatisfied=node->numberUnsatisfied();

/*
  Rebuild the subproblem for this node:	 Call addCuts() to adjust the model
  to recreate the subproblem for this node (set proper variable bounds, add
  cuts, create a basis).  This may result in the problem being fathomed by
  bound or infeasibility. Returns 1 if node is fathomed.
  Execute the current arm of the branch: If the problem survives, save the
  resulting variable bounds and call branch() to modify variable bounds
  according to the current arm of the branching object. If we're processing
  the final arm of the branching object, flag the node for removal from the
  live set.
*/
    CbcNodeInfo * nodeInfo = node->nodeInfo() ;
    newNode = NULL ;
    int branchesLeft=0;
    if (!addCuts(node,lastws,numberFixedNow_>numberFixedAtRoot_))
    { int i ;
      const double * lower = getColLower() ;
      const double * upper = getColUpper() ;
      for (i = 0 ; i < numberColumns ; i++)
      { lowerBefore[i]= lower[i] ;
	upperBefore[i]= upper[i] ; }
      if ((solverCharacteristics_->extraCharacteristics()&2)!=0) {
	solverCharacteristics_->setBeforeLower(lowerBefore);
	solverCharacteristics_->setBeforeUpper(upperBefore);
      }
      if (messageHandler()->logLevel()>2)
	node->modifiableBranchingObject()->print();
      if (!useOsiBranching) 
	branchesLeft = node->branch(NULL); // old way
      else
	branchesLeft = node->branch(solver_); // new way
      if (branchesLeft) {
        // set nodenumber correctly
        node->nodeInfo()->setNodeNumber(numberNodes2_);
        tree_->push(node) ;
        if (statistics_) {
          if (numberNodes2_==maximumStatistics_) {
            maximumStatistics_ = 2*maximumStatistics_;
            CbcStatistics ** temp = new CbcStatistics * [maximumStatistics_];
            memset(temp,0,maximumStatistics_*sizeof(CbcStatistics *));
            memcpy(temp,statistics_,numberNodes2_*sizeof(CbcStatistics *));
            delete [] statistics_;
            statistics_=temp;
          }
          assert (!statistics_[numberNodes2_]);
          statistics_[numberNodes2_]=new CbcStatistics(node);
        }
        numberNodes2_++;
	//nodeOnTree=true; // back on tree
	//deleteNode = false ;
#	ifdef CHECK_NODE
	printf("Node %x pushed back on tree - %d left, %d count\n",node,
	       nodeInfo->numberBranchesLeft(),
	       nodeInfo->numberPointingToThis()) ;
#	endif
      } else {
	//deleteNode = true ;
	if (!nodeInfo->numberBranchesLeft())
	  nodeInfo->allBranchesGone(); // can clean up
      }
      if ((specialOptions_&1)!=0) {
        /*
          This doesn't work as intended --- getRowCutDebugger will return null
          unless the current feasible solution region includes the optimal solution
          that RowCutDebugger knows. There's no way to tell inactive from off the
          optimal path.
        */
        const OsiRowCutDebugger *debugger = solver_->getRowCutDebugger() ;
        if (debugger) {
          onOptimalPath=true;
          printf("On optimal path\n") ;
        }
      }
      
/*
  Reoptimize, possibly generating cuts and/or using heuristics to find
  solutions.  Cut reference counts are unaffected unless we lose feasibility,
  in which case solveWithCuts() will make the adjustment.
*/
      phase_=2;
      cuts = OsiCuts() ;
      currentNumberCuts = solver_->getNumRows()-numberRowsAtContinuous_ ;
      int saveNumber = numberIterations_;
      if(solverCharacteristics_->solutionAddsCuts()) {
        int returnCode=resolve(node ? node->nodeInfo() : NULL,1);
        feasible = returnCode != 0;
        if (feasible) {
          int iObject ;
          int preferredWay ;
          int numberUnsatisfied = 0 ;
          memcpy(currentSolution_,solver_->getColSolution(),
                 numberColumns*sizeof(double)) ;
	  // point to useful information
	  OsiBranchingInformation usefulInfo=usefulInformation();
          
          for (iObject = 0 ; iObject < numberObjects_ ; iObject++) {
            double infeasibility =
              object_[iObject]->infeasibility(&usefulInfo,preferredWay) ;
            if (infeasibility ) numberUnsatisfied++ ;
          }
          if (returnCode>0) {
            if (numberUnsatisfied)   {
              feasible = solveWithCuts(cuts,maximumCutPasses_,node);
            } else {
              // may generate cuts and turn the solution
              //to an infeasible one
              feasible = solveWithCuts(cuts, 1,
                                       node);
#if 0
              currentNumberCuts_ = cuts.sizeRowCuts();
              if (currentNumberCuts_ >= maximumNumberCuts_) {
                maximumNumberCuts_ = currentNumberCuts;
                delete [] addedCuts_;
                addedCuts_ = new CbcCountRowCut * [maximumNumberCuts_];
              }
#endif
            }
          }
          // check extra info on feasibility
          if (!solverCharacteristics_->mipFeasible()) {
            feasible = false;
	    solverCharacteristics_->setMipBound(-COIN_DBL_MAX);
	  }
        }
      } else {
        // normal
	//int zzzzzz=0;
	//if (zzzzzz)
	//solver_->writeMps("before");
        feasible = solveWithCuts(cuts,maximumCutPasses_,node);
      }
      if ((specialOptions_&1)!=0&&onOptimalPath) {
        if (!solver_->getRowCutDebugger()) {
	  // dj fix did something???
	  solver_->writeMps("infeas2");
	  assert (solver_->getRowCutDebugger()) ;
	}
      }
      if (statistics_) {
        assert (numberNodes2_);
        assert (statistics_[numberNodes2_-1]);
        assert (statistics_[numberNodes2_-1]->node()==numberNodes2_-1);
        statistics_[numberNodes2_-1]->endOfBranch(numberIterations_-saveNumber,
                                               feasible ? solver_->getObjValue()
                                               : COIN_DBL_MAX);
      }
/*
  Are we still feasible? If so, create a node and do the work to attach a
  branching object, reoptimising as needed if chooseBranch() identifies
  monotone objects.

  Finally, attach a partial nodeInfo object and store away any cuts that we
  created back in solveWithCuts. addCuts() will initialise the reference
  counts for these new cuts.

  This next test can be problematic if we've discovered an
  alternate equivalent answer and subsequently fathom the solution
  known to the row cut debugger due to bounds.
*/
        if (onOptimalPath) {
	  bool objLim = solver_->isDualObjectiveLimitReached() ;
          if (!feasible && !objLim) {
            printf("infeas2\n");
            solver_->writeMps("infeas");
            CoinWarmStartBasis *slack =
              dynamic_cast<CoinWarmStartBasis *>(solver_->getEmptyWarmStart()) ;
            solver_->setWarmStart(slack);
            delete slack ;
            solver_->setHintParam(OsiDoReducePrint,false,OsiHintDo,0) ;
            solver_->initialSolve();
            assert (!solver_->isProvenOptimal());
          }
          assert (feasible || objLim);
        }
        bool checkingNode=false;
	if (feasible) {
          newNode = new CbcNode ;
          // Set objective value (not so obvious if NLP etc)
          setObjectiveValue(newNode,node);
	  anyAction =-1 ;
	  resolved = false ;
	  if (newNode->objectiveValue() >= getCutoff()) 
	    anyAction=-2;
          // only allow at most a few passes
          int numberPassesLeft=5;
          checkingNode=true;
	OsiSolverBranch * branches=NULL;
	// point to useful information
	anyAction = chooseBranch(newNode, numberPassesLeft,node, cuts,resolved,
				 lastws, lowerBefore, upperBefore, branches);
/*
  If we end up infeasible, we can delete the new node immediately. Since this
  node won't be needing the cuts we collected, decrement the reference counts.
  If we are feasible, then we'll be placing this node into the live set, so
  increment the reference count in the current (parent) nodeInfo.
*/
	if (anyAction == -2)
	  { delete newNode ;
	  newNode = NULL ;
          // say strong doing well
          if (checkingNode)
            setSpecialOptions(specialOptions_|8);
	  for (i = 0 ; i < currentNumberCuts_ ; i++)
	    { if (addedCuts_[i])
	      { if (!addedCuts_[i]->decrement(1))
		delete addedCuts_[i] ; } } }
	else
	  { nodeInfo->increment() ;
	  if ((numberNodes_%20)==0) {
	    // say strong not doing as well
	    setSpecialOptions(specialOptions_&~8);
	  }
	}
        }
/*
  At this point, there are three possibilities:
    * newNode is live and will require further branching to resolve
      (variable() >= 0). Increment the cut reference counts by
      numberBranches() to allow for use by children of this node, and
      decrement by 1 because we've executed one arm of the branch of our
      parent (consuming one reference). Before we push newNode onto the
      search tree, try for a heuristic solution.
    * We have a solution, in which case newNode is non-null but we have no
      branching variable. Decrement the cut counts and save the solution.
    * The node was found to be infeasible, in which case it's already been
      deleted, and newNode is null.
*/
        if (!eventHandler->event(CbcEventHandler::node)) {
          eventHappened_=true; // exit
        }
	assert (!newNode || newNode->objectiveValue() <= getCutoff()) ;
        if (statistics_) {
          assert (numberNodes2_);
          assert (statistics_[numberNodes2_-1]);
          assert (statistics_[numberNodes2_-1]->node()==numberNodes2_-1);
          if (newNode)
            statistics_[numberNodes2_-1]->updateInfeasibility(newNode->numberUnsatisfied());
          else
            statistics_[numberNodes2_-1]->sayInfeasible();
        }
	if (newNode) {
	  if (newNode->branchingObject() == NULL&&solverCharacteristics_->solverType()==4) {
	    // need to check if any cuts would do anything
	    OsiCuts theseCuts;
	    // reset probing info
	    if (probingInfo_)
	      probingInfo_->initializeFixing();
	    for (int i = 0;i<numberCutGenerators_;i++) {
	      bool generate = generator_[i]->normal();
	      // skip if not optimal and should be (maybe a cut generator has fixed variables)
	      if (generator_[i]->needsOptimalBasis()&&!solver_->basisIsAvailable())
		generate=false;
	      if (!generator_[i]->mustCallAgain())
		generate=false; // only special cuts
	      if (generate) {
		generator_[i]->generateCuts(theseCuts,true,solver_,NULL) ;
		int numberRowCutsAfter = theseCuts.sizeRowCuts() ;
		if (numberRowCutsAfter) {
		  // need dummy branch
		  newNode->setBranchingObject(new CbcDummyBranchingObject(this));
		  newNode->nodeInfo()->initializeInfo(1);
		  break;
		}
	      }
	    }
	  }
	  if (newNode->branchingObject())
	  { handler_->message(CBC_BRANCH,messages_)
	       << numberNodes_<< newNode->objectiveValue()
	       << newNode->numberUnsatisfied()<< newNode->depth()
	       << CoinMessageEol ;
	    // Increment cut counts (taking off current)
	    int numberLeft = newNode->numberBranches() ;
	    for (i = 0;i < currentNumberCuts_;i++)
	    { if (addedCuts_[i])
	      {
#		ifdef CHECK_CUT_COUNTS
		printf("Count on cut %x increased by %d\n",addedCuts_[i],
			numberLeft-1) ;
#		endif
		addedCuts_[i]->increment(numberLeft-1) ; } }

	    double estValue = newNode->guessedObjectiveValue() ;
	    int found = -1 ;
	    // no - overhead on small problems solver_->resolve() ;	// double check current optimal
	    // assert (!solver_->getIterationCount());
	    double * newSolution = new double [numberColumns] ;
	    double heurValue = getCutoff() ;
	    int iHeur ;
	    for (iHeur = 0 ; iHeur < numberHeuristics_ ; iHeur++)
	    { double saveValue = heurValue ;
	      int ifSol = heuristic_[iHeur]->solution(heurValue,newSolution) ;
	      if (ifSol > 0) {
                // new solution found
                found = iHeur ;
                incrementUsed(newSolution);
              }
	      else
	      if (ifSol < 0)	// just returning an estimate
	      { estValue = CoinMin(heurValue,estValue) ;
		heurValue = saveValue ; } }
	    if (found >= 0) {
              lastHeuristic_ = heuristic_[found];
              setBestSolution(CBC_ROUNDING,heurValue,newSolution) ;
            }
	    delete [] newSolution ;
	    newNode->setGuessedObjectiveValue(estValue) ;
	    tree_->push(newNode) ;
            if (statistics_) {
              if (numberNodes2_==maximumStatistics_) {
                maximumStatistics_ = 2*maximumStatistics_;
                CbcStatistics ** temp = new CbcStatistics * [maximumStatistics_];
                memset(temp,0,maximumStatistics_*sizeof(CbcStatistics *));
                memcpy(temp,statistics_,numberNodes2_*sizeof(CbcStatistics *));
                delete [] statistics_;
                statistics_=temp;
              }
              assert (!statistics_[numberNodes2_]);
              statistics_[numberNodes2_]=new CbcStatistics(newNode);
            }
            numberNodes2_++;
#	    ifdef CHECK_NODE
	    printf("Node %x pushed on tree c\n",newNode) ;
#	    endif
	  }
	  else
	  { 
            if(solverCharacteristics_ && //we may be in a non standard bab
               solverCharacteristics_->solutionAddsCuts()// we are in some kind of OA based bab.
               )
              {
                std::cerr<<"You should never get here"<<std::endl;
                throw CoinError("Nodes should not be fathomed on integer infeasibility in this setting",
                                "branchAndBound","CbcModel") ;
              }
            for (i = 0 ; i < currentNumberCuts_ ; i++)
	    { if (addedCuts_[i])
	      { if (!addedCuts_[i]->decrement(1))
		  delete addedCuts_[i] ; } }
	  double objectiveValue = newNode->objectiveValue();
	    setBestSolution(CBC_SOLUTION,objectiveValue,
			    solver_->getColSolution()) ;
            lastHeuristic_ = NULL;
            incrementUsed(solver_->getColSolution());
	    //assert(nodeInfo->numberPointingToThis() <= 2) ;
	    // avoid accidental pruning, if newNode was final branch arm
	    nodeInfo->increment();
	    delete newNode ;
	    nodeInfo->decrement() ; } }
/*
  This node has been completely expanded and can be removed from the live
  set.
*/
      if (branchesLeft)
      { 
      }
      else
      { 
	if (!nodeInfo->numberBranchesLeft())
	  nodeInfo->allBranchesGone(); // can clean up
	delete node ; }
    } else {
      // add cuts found to be infeasible (on bound)!
      abort();
      delete node;
    }
/*
  Delete cuts to get back to the original system.

  I'm thinking this is redundant --- the call to addCuts that conditions entry
  to this code block also performs this action.
*/
      int numberToDelete = getNumRows()-numberRowsAtContinuous_ ;
      if (numberToDelete)
      { int * delRows = new int[numberToDelete] ;
	int i ;
	for (i = 0 ; i < numberToDelete ; i++)
	{ delRows[i] = i+numberRowsAtContinuous_ ; }
	solver_->deleteRows(numberToDelete,delRows) ;
	delete [] delRows ; }
#else // end of not CBC_THREAD
      if (!numberThreads_) {
	doOneNode(this,node,createdNode);
      } else {
	threadStats[0]++;
	//need to think
	int iThread;
	// Start one off if any available
	for (iThread=0;iThread<numberThreads_;iThread++) {
	  if (threadInfo[iThread].returnCode==-1) 
	    break;
	}
	if (iThread<numberThreads_) {
	  threadInfo[iThread].node=node;
	  assert (threadInfo[iThread].returnCode==-1);
	  // say in use
	  threadInfo[iThread].returnCode=0;
	  threadModel[iThread]->moveToModel(this,0);
	  pthread_cond_signal(threadInfo[iThread].condition2); // unlock
	  threadCount[iThread]++;
	}
	lockThread();
	locked=true;
	// see if any finished
	for (iThread=0;iThread<numberThreads_;iThread++) {
	  if (threadInfo[iThread].returnCode>0) 
	    break;
	}
	unlockThread();
	locked=false;
	if (iThread<numberThreads_) {
	  threadModel[iThread]->moveToModel(this,1);
	  assert (threadInfo[iThread].returnCode==1);
	  // say available
	  threadInfo[iThread].returnCode=-1;
	  // carry on
	  threadStats[3]++;
	} else {
	  // Start one off if any available
	  for (iThread=0;iThread<numberThreads_;iThread++) {
	    if (threadInfo[iThread].returnCode==-1) 
	      break;
	  }
	  if (iThread<numberThreads_) {
	    lockThread();
	    locked=true;
	    // If any on tree get
	    if (!tree_->empty()) {
	      //node = tree_->bestNode(cutoff) ;
	      //assert (node);
	      threadStats[1]++;
	      continue; // ** get another node
	    }
	    unlockThread();
	    locked=false;
	  }
	  // wait (for debug could sleep and use test)
	  bool finished=false;
	  while (!finished) {
	    pthread_mutex_lock(&condition_mutex);
	    struct timespec absTime;
	    clock_gettime(CLOCK_REALTIME,&absTime);
	    double time = absTime.tv_sec+1.0e-9*absTime.tv_nsec;
	    absTime.tv_nsec += 1000000; // millisecond
	    if (absTime.tv_nsec>=1000000000) {
	      absTime.tv_nsec -= 1000000000;
	      absTime.tv_sec++;
	    }
	    pthread_cond_timedwait(&condition_main,&condition_mutex,&absTime);
	    clock_gettime(CLOCK_REALTIME,&absTime);
	    double time2 = absTime.tv_sec+1.0e-9*absTime.tv_nsec;
	    timeWaiting += time2-time;
	    pthread_mutex_unlock(&condition_mutex);
	    for (iThread=0;iThread<numberThreads_;iThread++) {
	      if (threadInfo[iThread].returnCode>0) {
		finished=true;
		break;
	      } else if (threadInfo[iThread].returnCode==0) {
		pthread_cond_signal(threadInfo[iThread].condition2); // unlock
	      }
	    }
	  }
	  assert (iThread<numberThreads_);
	  threadModel[iThread]->moveToModel(this,1);
	  node = threadInfo[iThread].node;
	  threadInfo[iThread].node=NULL;
	  assert (threadInfo[iThread].returnCode==1);
	  // say available
	  threadInfo[iThread].returnCode=-1;
	  // carry on
	  threadStats[2]++;
	}
      }
      //lastDepth=node->depth();
      //lastUnsatisfied=node->numberUnsatisfied();
#endif // end of CBC_THREAD
  }
#ifdef CBC_THREAD
  if (numberThreads_) {
    //printf("stats ");
    //for (unsigned int j=0;j<sizeof(threadStats)/sizeof(int);j++)
    //printf("%d ",threadStats[j]);
    //printf("\n");
    int i;
    // Seems to be bug in CoinCpu on Linux - does threads as well despite documentation
    double time=0.0;
    for (i=0;i<numberThreads_;i++) 
      time += threadInfo[i].timeInThread;
    bool goodTimer = time<(getCurrentSeconds());
    //bool stopped = (!(numberNodes_ < intParam_[CbcMaxNumNode] &&
    //	      numberSolutions_ < intParam_[CbcMaxNumSol] &&
    //	      totalTime < dblParam_[CbcMaximumSeconds] &&
    //	      !stoppedOnGap_&&!eventHappened_));
    for (i=0;i<numberThreads_;i++) {
      while (threadInfo[i].returnCode==0) {
	pthread_cond_signal(threadInfo[i].condition2); // unlock
	pthread_mutex_lock(&condition_mutex);
	struct timespec absTime;
	clock_gettime(CLOCK_REALTIME,&absTime);
	absTime.tv_nsec += 1000000; // millisecond
	if (absTime.tv_nsec>=1000000000) {
	  absTime.tv_nsec -= 1000000000;
	  absTime.tv_sec++;
	}
	pthread_cond_timedwait(&condition_main,&condition_mutex,&absTime);
	clock_gettime(CLOCK_REALTIME,&absTime);
	pthread_mutex_unlock(&condition_mutex);
      }
      pthread_cond_signal(threadInfo[i].condition2); // unlock
      pthread_mutex_lock(&condition_mutex); // not sure necessary but have had one hang on interrupt
      threadModel[i]->numberThreads_=0; // say exit
      threadInfo[i].returnCode=0;
      pthread_mutex_unlock(&condition_mutex);
      pthread_cond_signal(threadInfo[i].condition2); // unlock
      //if (!stopped)
	pthread_join(threadId[i],NULL);
	//else
	//pthread_kill(threadId[i]); // kill rather than try and synchronize
      threadModel[i]->moveToModel(this,2);
      pthread_mutex_destroy (threadInfo[i].mutex2);
      pthread_cond_destroy (threadInfo[i].condition2);
      assert (threadInfo[i].numberTimesLocked==threadInfo[i].numberTimesUnlocked);
      handler_->message(CBC_THREAD_STATS,messages_)
	<<"Thread";
      handler_->printing(true)
	<<i<<threadCount[i]<<threadInfo[i].timeWaitingToStart;
      handler_->printing(goodTimer)<<threadInfo[i].timeInThread;
      handler_->printing(false)<<0.0;
      handler_->printing(true)<<threadInfo[i].numberTimesLocked
	<<threadInfo[i].timeLocked<<threadInfo[i].timeWaitingToLock
	<<CoinMessageEol;
    }
    assert (threadInfo[numberThreads_].numberTimesLocked==threadInfo[numberThreads_].numberTimesUnlocked);
    handler_->message(CBC_THREAD_STATS,messages_)
      <<"Main thread";
    handler_->printing(false)<<0<<0<<0.0;
    handler_->printing(false)<<0.0;
    handler_->printing(true)<<timeWaiting;
    handler_->printing(true)<<threadInfo[numberThreads_].numberTimesLocked
      <<threadInfo[numberThreads_].timeLocked<<threadInfo[numberThreads_].timeWaitingToLock
      <<CoinMessageEol;
    pthread_mutex_destroy (&mutex);
    pthread_cond_destroy (&condition_main);
    pthread_mutex_destroy (&condition_mutex);
    // delete models (here in case some point to others)
    for (i=0;i<numberThreads_;i++) {
      delete threadModel[i];
    }
    delete [] mutex2;
    delete [] condition2;
    delete [] threadId;
    delete [] threadInfo;
    delete [] threadModel;
    delete [] threadCount;
    mutex_=NULL;
    // adjust time to allow for children on some systems
    dblParam_[CbcStartSeconds] -= CoinCpuTimeJustChildren();
  }
#endif
/*
  End of the non-abort actions. The next block of code is executed if we've
  aborted because we hit one of the limits. Clean up by deleting the live set
  and break out of the node processing loop. Note that on an abort, node may
  have been pushed back onto the tree for further processing, in which case
  it'll be deleted in cleanTree. We need to check.
*/
    if (!(numberNodes_ < intParam_[CbcMaxNumNode] &&
	numberSolutions_ < intParam_[CbcMaxNumSol] &&
	totalTime < dblParam_[CbcMaximumSeconds] &&
	!stoppedOnGap_&&!eventHappened_)) {
      if (tree_->size())
	tree_->cleanTree(this,-COIN_DBL_MAX,bestPossibleObjective_) ;
      delete nextRowCut_;
      if (stoppedOnGap_)
	{ messageHandler()->message(CBC_GAP,messages())
	  << bestObjective_-bestPossibleObjective_
	  << dblParam_[CbcAllowableGap]
	  << dblParam_[CbcAllowableFractionGap]*100.0
	  << CoinMessageEol ;
        secondaryStatus_ = 2;
	status_ = 0 ; }
	else
	  if (isNodeLimitReached())
	    { handler_->message(CBC_MAXNODES,messages_) << CoinMessageEol ;
	    secondaryStatus_ = 3;
	    status_ = 1 ; }
	  else
	if (totalTime >= dblParam_[CbcMaximumSeconds])
	  { handler_->message(CBC_MAXTIME,messages_) << CoinMessageEol ; 
	  secondaryStatus_ = 4;
	  status_ = 1 ; }
	else
	  if (eventHappened_)
	    { handler_->message(CBC_EVENT,messages_) << CoinMessageEol ; 
	    secondaryStatus_ = 5;
	    status_ = 5 ; }
	  else
	    { handler_->message(CBC_MAXSOLS,messages_) << CoinMessageEol ;
	    secondaryStatus_ = 6;
	    status_ = 1 ; }
    }
/*
  That's it, we've exhausted the search tree, or broken out of the loop because
  we hit some limit on evaluation.

  We may have got an intelligent tree so give it one more chance
*/
  // Tell solver we are not in Branch and Cut
  solver_->setHintParam(OsiDoInBranchAndCut,false,OsiHintDo,NULL) ;
  tree_->endSearch();
  //  If we did any sub trees - did we give up on any?
  if ( numberStoppedSubTrees_)
    status_=1;
  if (!status_) {
    // Set best possible unless stopped on gap
    if(secondaryStatus_ != 2)
      bestPossibleObjective_=bestObjective_;
    handler_->message(CBC_END_GOOD,messages_)
      << bestObjective_ << numberIterations_ << numberNodes_<<getCurrentSeconds()
      << CoinMessageEol ;
  } else {
    handler_->message(CBC_END,messages_)
      << bestObjective_ <<bestPossibleObjective_
      << numberIterations_ << numberNodes_<<getCurrentSeconds()
      << CoinMessageEol ;
  }
  if (numberStrongIterations_)
    handler_->message(CBC_STRONG_STATS,messages_)
      << strongInfo_[0] << numberStrongIterations_ << strongInfo_[2]
      << strongInfo_[1] << CoinMessageEol ;
  handler_->message(CBC_OTHER_STATS,messages_)
      << maximumDepthActual_ 
      << numberDJFixed_ << CoinMessageEol ;
  if (doStatistics==100) {
    for (int i=0;i<numberObjects_;i++) {
      CbcSimpleIntegerDynamicPseudoCost * obj =
	dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object_[i]) ;
      if (obj)
	obj->print();
    }
  }
  if (statistics_) {
    // report in some way
    int * lookup = new int[numberObjects_];
    int i;
    for (i=0;i<numberObjects_;i++) 
      lookup[i]=-1;
    bool goodIds=true;
    for (i=0;i<numberObjects_;i++) {
      int iColumn = object_[i]->columnNumber();
      if(iColumn>=0&&iColumn<numberColumns) {
        if (lookup[i]==-1) {
          lookup[i]=iColumn;
        } else {
          goodIds=false;
          break;
        }
      } else {
        goodIds=false;
        break;
      }
    }
    if (!goodIds) {
      delete [] lookup;
      lookup=NULL;
    }
    if (doStatistics==3) {
      printf("  node parent depth column   value                    obj      inf\n");
      for ( i=0;i<numberNodes2_;i++) {
        statistics_[i]->print(lookup);
      }
    }
    if (doStatistics>1) {
      // Find last solution
      int k;
      for (k=numberNodes2_-1;k>=0;k--) {
        if (statistics_[k]->endingObjective()!=COIN_DBL_MAX&&
            !statistics_[k]->endingInfeasibility())
          break;
      }
      if (k>=0) {
        int depth=statistics_[k]->depth();
        int * which = new int[depth+1];
        for (i=depth;i>=0;i--) {
          which[i]=k;
          k=statistics_[k]->parentNode();
        }
        printf("  node parent depth column   value                    obj      inf\n");
        for (i=0;i<=depth;i++) {
          statistics_[which[i]]->print(lookup);
        }
        delete [] which;
      }
    }
    // now summary
    int maxDepth=0;
    double averageSolutionDepth=0.0;
    int numberSolutions=0;
    double averageCutoffDepth=0.0;
    double averageSolvedDepth=0.0;
    int numberCutoff=0;
    int numberDown=0;
    int numberFirstDown=0;
    double averageInfDown=0.0;
    double averageObjDown=0.0;
    int numberCutoffDown=0;
    int numberUp=0;
    int numberFirstUp=0;
    double averageInfUp=0.0;
    double averageObjUp=0.0;
    int numberCutoffUp=0;
    double averageNumberIterations1=0.0;
    double averageValue=0.0;
    for ( i=0;i<numberNodes2_;i++) {
      int depth =  statistics_[i]->depth(); 
      int way =  statistics_[i]->way(); 
      double value = statistics_[i]->value(); 
      double startingObjective =  statistics_[i]->startingObjective(); 
      int startingInfeasibility = statistics_[i]->startingInfeasibility(); 
      double endingObjective = statistics_[i]->endingObjective(); 
      int endingInfeasibility = statistics_[i]->endingInfeasibility(); 
      maxDepth = CoinMax(depth,maxDepth);
      // Only for completed
      averageNumberIterations1 += statistics_[i]->numberIterations();
      averageValue += value;
      if (endingObjective!=COIN_DBL_MAX&&!endingInfeasibility) {
        numberSolutions++;
        averageSolutionDepth += depth;
      }
      if (endingObjective==COIN_DBL_MAX) {
        numberCutoff++;
        averageCutoffDepth += depth;
        if (way<0) {
          numberDown++;
          numberCutoffDown++;
          if (way==-1)
            numberFirstDown++;
        } else {
          numberUp++;
          numberCutoffUp++;
          if (way==1)
            numberFirstUp++;
        }
      } else {
        averageSolvedDepth += depth;
        if (way<0) {
          numberDown++;
          averageInfDown += startingInfeasibility-endingInfeasibility;
          averageObjDown += endingObjective-startingObjective;
          if (way==-1)
            numberFirstDown++;
        } else {
          numberUp++;
          averageInfUp += startingInfeasibility-endingInfeasibility;
          averageObjUp += endingObjective-startingObjective;
          if (way==1)
            numberFirstUp++;
        }
      }
    }
    // Now print
    if (numberSolutions)
      averageSolutionDepth /= (double) numberSolutions;
    int numberSolved = numberNodes2_-numberCutoff;
    double averageNumberIterations2=numberIterations_-averageNumberIterations1
      -numberIterationsAtContinuous;
    if(numberCutoff) {
      averageCutoffDepth /= (double) numberCutoff;
      averageNumberIterations2 /= (double) numberCutoff;
    }
    if (numberNodes2_) 
      averageValue /= (double) numberNodes2_;
    if (numberSolved) {
      averageNumberIterations1 /= (double) numberSolved;
      averageSolvedDepth /= (double) numberSolved;
    }
    printf("%d solution(s) were found (by branching) at an average depth of %g\n",
           numberSolutions,averageSolutionDepth);
    printf("average value of variable being branched on was %g\n",
           averageValue);
    printf("%d nodes were cutoff at an average depth of %g with iteration count of %g\n",
           numberCutoff,averageCutoffDepth,averageNumberIterations2);
    printf("%d nodes were solved at an average depth of %g with iteration count of %g\n",
           numberSolved,averageSolvedDepth,averageNumberIterations1);
    if (numberDown) {
      averageInfDown /= (double) numberDown;
      averageObjDown /= (double) numberDown;
    }
    printf("Down %d nodes (%d first, %d second) - %d cutoff, rest decrease numinf %g increase obj %g\n",
           numberDown,numberFirstDown,numberDown-numberFirstDown,numberCutoffDown,
           averageInfDown,averageObjDown);
    if (numberUp) {
      averageInfUp /= (double) numberUp;
      averageObjUp /= (double) numberUp;
    }
    printf("Up %d nodes (%d first, %d second) - %d cutoff, rest decrease numinf %g increase obj %g\n",
           numberUp,numberFirstUp,numberUp-numberFirstUp,numberCutoffUp,
           averageInfUp,averageObjUp);
    for ( i=0;i<numberNodes2_;i++) 
      delete statistics_[i];
    delete [] statistics_;
    statistics_=NULL;
    maximumStatistics_=0;
    delete [] lookup;
  }
/*
  If we think we have a solution, restore and confirm it with a call to
  setBestSolution().  We need to reset the cutoff value so as not to fathom
  the solution on bounds.  Note that calling setBestSolution( ..., true)
  leaves the continuousSolver_ bounds vectors fixed at the solution value.

  Running resolve() here is a failsafe --- setBestSolution has already
  reoptimised using the continuousSolver_. If for some reason we fail to
  prove optimality, run the problem again after instructing the solver to
  tell us more.

  If all looks good, replace solver_ with continuousSolver_, so that the
  outside world will be able to obtain information about the solution using
  public methods.
*/
  if (bestSolution_&&(solverCharacteristics_->solverType()<2||solverCharacteristics_->solverType()==4)) 
  { setCutoff(1.0e50) ; // As best solution should be worse than cutoff
    phase_=5;
    double increment = getDblParam(CbcModel::CbcCutoffIncrement) ;
    if ((specialOptions_&4)!=0)
      bestObjective_ += 100.0*increment+1.0e-3; // only set if we are going to solve
    setBestSolution(CBC_END_SOLUTION,bestObjective_,bestSolution_,true) ;
    continuousSolver_->resolve() ;
    if (!continuousSolver_->isProvenOptimal())
    { continuousSolver_->messageHandler()->setLogLevel(2) ;
      continuousSolver_->initialSolve() ; }
    delete solver_ ;
    // above deletes solverCharacteristics_
    solverCharacteristics_ = NULL;
    solver_ = continuousSolver_ ;
    setPointers(solver_);
    continuousSolver_ = NULL ; }
/*
  Clean up dangling objects. continuousSolver_ may already be toast.
*/
  delete lastws ;
  delete [] whichGenerator_ ;
  whichGenerator_=NULL;
  delete [] lowerBefore ;
  delete [] upperBefore ;
  delete [] walkback_ ;
  walkback_ = NULL ;
  delete [] addedCuts_ ;
  addedCuts_ = NULL ;
  //delete persistentInfo;
  // Get rid of characteristics
  solverCharacteristics_=NULL;
  if (continuousSolver_)
  { delete continuousSolver_ ;
    continuousSolver_ = NULL ; }
/*
  Destroy global cuts by replacing with an empty OsiCuts object.
*/
  globalCuts_= OsiCuts() ;
  if (!bestSolution_) {
    // make sure lp solver is infeasible
    int numberColumns = solver_->getNumCols();
    const double * columnLower = solver_->getColLower();
    int iColumn;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (solver_->isInteger(iColumn)) 
	solver_->setColUpper(iColumn,columnLower[iColumn]);
    }
    solver_->initialSolve();
  }
  if (strategy_&&strategy_->preProcessState()>0) {
    // undo preprocessing
    CglPreProcess * process = strategy_->process();
    assert (process);
    int n = originalSolver->getNumCols();
    if (bestSolution_) {
      delete [] bestSolution_;
      bestSolution_ = new double [n];
      process->postProcess(*solver_);
    }
    strategy_->deletePreProcess();
    // Solution now back in originalSolver
    delete solver_;
    solver_=originalSolver;
    if (bestSolution_)
      memcpy(bestSolution_,solver_->getColSolution(),n*sizeof(double));
    // put back original objects if there were any
    if (originalObject) {
      int iColumn;
      assert (ownObjects_);
      for (iColumn=0;iColumn<numberObjects_;iColumn++) 
        delete object_[iColumn];
      delete [] object_;
      numberObjects_ = numberOriginalObjects;
      object_=originalObject;
      delete [] integerVariable_;
      numberIntegers_=0;
      for (iColumn=0;iColumn<n;iColumn++) {
        if (solver_->isInteger(iColumn))
          numberIntegers_++;
      }
      integerVariable_ = new int[numberIntegers_];
      numberIntegers_=0;
      for (iColumn=0;iColumn<n;iColumn++) {
        if (solver_->isInteger(iColumn))
            integerVariable_[numberIntegers_++]=iColumn;
      }
    }
  }
#ifdef CLP_QUICK_OPTIONS
  {
    OsiClpSolverInterface * clpSolver 
      = dynamic_cast<OsiClpSolverInterface *> (solver_);
    if (clpSolver) {
      // Try and re-use regions   
      ClpSimplex * simplex = clpSolver->getModelPtr();
      simplex->setPersistenceFlag(0);
      clpSolver->deleteScaleFactors();
      clpSolver->setSpecialOptions(clpSolver->specialOptions()&(~131072));
      simplex->setSpecialOptions(simplex->specialOptions()&(~131072));
      simplex->setAlphaAccuracy(-1.0);
      //clpSolver->setSpecialOptions((clpSolver->specialOptions()&~128)|65536);
    }
  }
#endif
  return ;
 }


// Solve the initial LP relaxation 
void 
CbcModel::initialSolve() 
{
  assert (solver_);
  assert (!solverCharacteristics_);
  OsiBabSolver * solverCharacteristics = dynamic_cast<OsiBabSolver *> (solver_->getAuxiliaryInfo());
  if (solverCharacteristics) {
    solverCharacteristics_ = solverCharacteristics;
  } else {
    // replace in solver
    OsiBabSolver defaultC;
    solver_->setAuxiliaryInfo(&defaultC);
    solverCharacteristics_ = dynamic_cast<OsiBabSolver *> (solver_->getAuxiliaryInfo());
  }
  solverCharacteristics_->setSolver(solver_);
  solver_->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,NULL) ;
  solver_->initialSolve();
  solver_->setHintParam(OsiDoInBranchAndCut,false,OsiHintDo,NULL) ;
  // But set up so Jon Lee will be happy
  status_=-1;
  secondaryStatus_ = -1;
  originalContinuousObjective_ = solver_->getObjValue()*solver_->getObjSense();
  delete [] continuousSolution_;
  continuousSolution_ = CoinCopyOfArray(solver_->getColSolution(),
                                             solver_->getNumCols());
  setPointers(solver_);
  solverCharacteristics_ = NULL;
}

/*! \brief Get an empty basis object

  Return an empty CoinWarmStartBasis object with the requested capacity,
  appropriate for the current solver. The object is cloned from the object
  cached as emptyWarmStart_. If there is no cached object, the routine
  queries the solver for a warm start object, empties it, and caches the
  result.
*/

CoinWarmStartBasis *CbcModel::getEmptyBasis (int ns, int na) const

{ CoinWarmStartBasis *emptyBasis ;
/*
  Acquire an empty basis object, if we don't yet have one.
*/
  if (emptyWarmStart_ == 0)
  { if (solver_ == 0)
    { throw CoinError("Cannot construct basis without solver!",
		      "getEmptyBasis","CbcModel") ; }
    emptyBasis =
	dynamic_cast<CoinWarmStartBasis *>(solver_->getEmptyWarmStart()) ;
    if (emptyBasis == 0)
    { throw CoinError(
	"Solver does not appear to use a basis-oriented warm start.",
	"getEmptyBasis","CbcModel") ; }
    emptyBasis->setSize(0,0) ;
    emptyWarmStart_ = dynamic_cast<CoinWarmStart *>(emptyBasis) ; }
/*
  Clone the empty basis object, resize it as requested, and return.
*/
  emptyBasis = dynamic_cast<CoinWarmStartBasis *>(emptyWarmStart_->clone()) ;
  assert(emptyBasis) ;
  if (ns != 0 || na != 0) emptyBasis->setSize(ns,na) ;

  return (emptyBasis) ; }
    

/** Default Constructor

  Creates an empty model without an associated solver.
*/
CbcModel::CbcModel() 

:
  solver_(NULL),
  ourSolver_(true),
  continuousSolver_(NULL),
  referenceSolver_(NULL),
  defaultHandler_(true),
  emptyWarmStart_(NULL),
  bestObjective_(COIN_DBL_MAX),
  bestPossibleObjective_(COIN_DBL_MAX),
  sumChangeObjective1_(0.0),
  sumChangeObjective2_(0.0),
  bestSolution_(NULL),
  currentSolution_(NULL),
  testSolution_(NULL),
  minimumDrop_(1.0e-4),
  numberSolutions_(0),
  stateOfSearch_(0),
  hotstartSolution_(NULL),
  hotstartPriorities_(NULL),
  numberHeuristicSolutions_(0),
  numberNodes_(0),
  numberNodes2_(0),
  numberIterations_(0),
  status_(-1),
  secondaryStatus_(-1),
  numberIntegers_(0),
  numberRowsAtContinuous_(0),
  maximumNumberCuts_(0),
  phase_(0),
  currentNumberCuts_(0),
  maximumDepth_(0),
  walkback_(NULL),
  addedCuts_(NULL),
  nextRowCut_(NULL),
  currentNode_(NULL),
  integerVariable_(NULL),
  integerInfo_(NULL),
  continuousSolution_(NULL),
  usedInSolution_(NULL),
  specialOptions_(0),
  subTreeModel_(NULL),
  numberStoppedSubTrees_(0),
  mutex_(NULL),
  presolve_(0),
  numberStrong_(5),
  numberBeforeTrust_(0),
  numberPenalties_(20),
  penaltyScaleFactor_(3.0),
  numberAnalyzeIterations_(0),
  analyzeResults_(NULL),
  numberInfeasibleNodes_(0),
  problemType_(0),
  printFrequency_(0),
  numberCutGenerators_(0),
  generator_(NULL),
  virginGenerator_(NULL),
  numberHeuristics_(0),
  heuristic_(NULL),
  lastHeuristic_(NULL),
  eventHandler_(0),
  numberObjects_(0),
  object_(NULL),
  ownObjects_(true),
  originalColumns_(NULL),
  howOftenGlobalScan_(1),
  numberGlobalViolations_(0),
  continuousObjective_(COIN_DBL_MAX),
  originalContinuousObjective_(COIN_DBL_MAX),
  continuousInfeasibilities_(COIN_INT_MAX),
  maximumCutPassesAtRoot_(20),
  maximumCutPasses_(10),
  preferredWay_(0),
  currentPassNumber_(0),
  maximumWhich_(1000),
  whichGenerator_(NULL),
  maximumStatistics_(0),
  statistics_(NULL),
  maximumDepthActual_(0),
  numberDJFixed_(0.0),
  probingInfo_(NULL),
  numberFixedAtRoot_(0),
  numberFixedNow_(0),
  stoppedOnGap_(false),
  eventHappened_(false),
  numberLongStrong_(0),
  numberOldActiveCuts_(0),
  numberNewCuts_(0),
  sizeMiniTree_(0),
  searchStrategy_(-1),
  numberStrongIterations_(0),
  resolveAfterTakeOffCuts_(true),
#if NEW_UPDATE_OBJECT>1
  numberUpdateItems_(0),
  maximumNumberUpdateItems_(0),
  updateItems_(NULL),
#endif
  numberThreads_(0),
  threadMode_(0)
{
  memset(intParam_,0,sizeof(intParam_));
  intParam_[CbcMaxNumNode] = 2147483647;
  intParam_[CbcMaxNumSol] = 9999999;
  intParam_[CbcFathomDiscipline] = 0;

  dblParam_[CbcIntegerTolerance] = 1e-6;
  dblParam_[CbcInfeasibilityWeight] = 0.0;
  dblParam_[CbcCutoffIncrement] = 1e-5;
  dblParam_[CbcAllowableGap] = 1.0e-10;
  dblParam_[CbcAllowableFractionGap] = 0.0;
  dblParam_[CbcMaximumSeconds] = 1.0e100;
  dblParam_[CbcCurrentCutoff] = 1.0e100;
  dblParam_[CbcOptimizationDirection] = 1.0;
  dblParam_[CbcCurrentObjectiveValue] = 1.0e100;
  dblParam_[CbcCurrentMinimizationObjectiveValue] = 1.0e100;
  dblParam_[CbcStartSeconds] = 0.0;
  strongInfo_[0]=0;
  strongInfo_[1]=0;
  strongInfo_[2]=0;
  solverCharacteristics_ = NULL;
  nodeCompare_=new CbcCompareDefault();;
  problemFeasibility_=new CbcFeasibilityBase();
  tree_= new CbcTree();
  branchingMethod_=NULL;
  cutModifier_=NULL;
  strategy_=NULL;
  parentModel_=NULL;
  cbcColLower_ = NULL;
  cbcColUpper_ = NULL;
  cbcRowLower_ = NULL;
  cbcRowUpper_ = NULL;
  cbcColSolution_ = NULL;
  cbcRowPrice_ = NULL;
  cbcReducedCost_ = NULL;
  cbcRowActivity_ = NULL;
  appData_=NULL;
  handler_ = new CoinMessageHandler();
  handler_->setLogLevel(2);
  messages_ = CbcMessage();
  eventHandler_ = new CbcEventHandler() ;
}

/** Constructor from solver.

  Creates a model complete with a clone of the solver passed as a parameter.
*/

CbcModel::CbcModel(const OsiSolverInterface &rhs)
:
  continuousSolver_(NULL),
  referenceSolver_(NULL),
  defaultHandler_(true),
  emptyWarmStart_(NULL),
  bestObjective_(COIN_DBL_MAX),
  bestPossibleObjective_(COIN_DBL_MAX),
  sumChangeObjective1_(0.0),
  sumChangeObjective2_(0.0),
  minimumDrop_(1.0e-4),
  numberSolutions_(0),
  stateOfSearch_(0),
  hotstartSolution_(NULL),
  hotstartPriorities_(NULL),
  numberHeuristicSolutions_(0),
  numberNodes_(0),
  numberNodes2_(0),
  numberIterations_(0),
  status_(-1),
  secondaryStatus_(-1),
  numberRowsAtContinuous_(0),
  maximumNumberCuts_(0),
  phase_(0),
  currentNumberCuts_(0),
  maximumDepth_(0),
  walkback_(NULL),
  addedCuts_(NULL),
  nextRowCut_(NULL),
  currentNode_(NULL),
  integerInfo_(NULL),
  specialOptions_(0),
  subTreeModel_(NULL),
  numberStoppedSubTrees_(0),
  mutex_(NULL),
  presolve_(0),
  numberStrong_(5),
  numberBeforeTrust_(0),
  numberPenalties_(20),
  penaltyScaleFactor_(3.0),
  numberAnalyzeIterations_(0),
  analyzeResults_(NULL),
  numberInfeasibleNodes_(0),
  problemType_(0),
  printFrequency_(0),
  numberCutGenerators_(0),
  generator_(NULL),
  virginGenerator_(NULL),
  numberHeuristics_(0),
  heuristic_(NULL),
  lastHeuristic_(NULL),
  eventHandler_(0),
  numberObjects_(0),
  object_(NULL),
  ownObjects_(true),
  originalColumns_(NULL),
  howOftenGlobalScan_(1),
  numberGlobalViolations_(0),
  continuousObjective_(COIN_DBL_MAX),
  originalContinuousObjective_(COIN_DBL_MAX),
  continuousInfeasibilities_(COIN_INT_MAX),
  maximumCutPassesAtRoot_(20),
  maximumCutPasses_(10),
  preferredWay_(0),
  currentPassNumber_(0),
  maximumWhich_(1000),
  whichGenerator_(NULL),
  maximumStatistics_(0),
  statistics_(NULL),
  maximumDepthActual_(0),
  numberDJFixed_(0.0),
  probingInfo_(NULL),
  numberFixedAtRoot_(0),
  numberFixedNow_(0),
  stoppedOnGap_(false),
  eventHappened_(false),
  numberLongStrong_(0),
  numberOldActiveCuts_(0),
  numberNewCuts_(0),
  sizeMiniTree_(0),
  searchStrategy_(-1),
  numberStrongIterations_(0),
  resolveAfterTakeOffCuts_(true),
#if NEW_UPDATE_OBJECT>1
  numberUpdateItems_(0),
  maximumNumberUpdateItems_(0),
  updateItems_(NULL),
#endif
  numberThreads_(0),
  threadMode_(0)
{
  memset(intParam_,0,sizeof(intParam_));
  intParam_[CbcMaxNumNode] = 2147483647;
  intParam_[CbcMaxNumSol] = 9999999;
  intParam_[CbcFathomDiscipline] = 0;

  dblParam_[CbcIntegerTolerance] = 1e-6;
  dblParam_[CbcInfeasibilityWeight] = 0.0;
  dblParam_[CbcCutoffIncrement] = 1e-5;
  dblParam_[CbcAllowableGap] = 1.0e-10;
  dblParam_[CbcAllowableFractionGap] = 0.0;
  dblParam_[CbcMaximumSeconds] = 1.0e100;
  dblParam_[CbcCurrentCutoff] = 1.0e100;
  dblParam_[CbcOptimizationDirection] = 1.0;
  dblParam_[CbcCurrentObjectiveValue] = 1.0e100;
  dblParam_[CbcCurrentMinimizationObjectiveValue] = 1.0e100;
  dblParam_[CbcStartSeconds] = 0.0;
  strongInfo_[0]=0;
  strongInfo_[1]=0;
  strongInfo_[2]=0;
  solverCharacteristics_ = NULL;

  nodeCompare_=new CbcCompareDefault();;
  problemFeasibility_=new CbcFeasibilityBase();
  tree_= new CbcTree();
  branchingMethod_=NULL;
  cutModifier_=NULL;
  strategy_=NULL;
  parentModel_=NULL;
  appData_=NULL;
  handler_ = new CoinMessageHandler();
  handler_->setLogLevel(2);
  messages_ = CbcMessage();
  eventHandler_ = new CbcEventHandler() ;
  solver_ = rhs.clone();
  referenceSolver_ = solver_->clone();
  ourSolver_ = true ;
  cbcColLower_ = NULL;
  cbcColUpper_ = NULL;
  cbcRowLower_ = NULL;
  cbcRowUpper_ = NULL;
  cbcColSolution_ = NULL;
  cbcRowPrice_ = NULL;
  cbcReducedCost_ = NULL;
  cbcRowActivity_ = NULL;

  // Initialize solution and integer variable vectors
  bestSolution_ = NULL; // to say no solution found
  numberIntegers_=0;
  int numberColumns = solver_->getNumCols();
  int iColumn;
  if (numberColumns) {
    // Space for current solution
    currentSolution_ = new double[numberColumns];
    continuousSolution_ = new double[numberColumns];
    usedInSolution_ = new int[numberColumns];
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if( solver_->isInteger(iColumn)) 
	numberIntegers_++;
    }
  } else {
    // empty model
    currentSolution_=NULL;
    continuousSolution_=NULL;
    usedInSolution_=NULL;
  }
  testSolution_=currentSolution_;
  if (numberIntegers_) {
    integerVariable_ = new int [numberIntegers_];
    numberIntegers_=0;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if( solver_->isInteger(iColumn)) 
	integerVariable_[numberIntegers_++]=iColumn;
    }
  } else {
    integerVariable_ = NULL;
  }
}

/*
  Assign a solver to the model (model assumes ownership)

  The integer variable vector is initialized if it's not already present.
  If deleteSolver then current solver deleted (if model owned)

  Assuming ownership matches usage in OsiSolverInterface
  (cf. assignProblem, loadProblem).

  TODO: What to do about solver parameters? A simple copy likely won't do it,
	because the SI must push the settings into the underlying solver. In
	the context of switching solvers in cbc, this means that command line
	settings will get lost. Stash the command line somewhere and reread it
	here, maybe?
  
  TODO: More generally, how much state should be transferred from the old
	solver to the new solver? Best perhaps to see how usage develops.
	What's done here mimics the CbcModel(OsiSolverInterface) constructor.
*/
void
CbcModel::assignSolver(OsiSolverInterface *&solver, bool deleteSolver)

{
  // resize best solution if exists
  if (bestSolution_) {
    int nOld = solver_->getNumCols();
    int nNew = solver->getNumCols();
    if (nNew>nOld) {
      double * temp = new double[nNew];
      memcpy(temp,bestSolution_,nOld*sizeof(double));
      memset(temp+nOld,0,(nNew-nOld)*sizeof(double));
      delete [] bestSolution_;
      bestSolution_=temp;
    }
  }
  // Keep the current message level for solver (if solver exists)
  if (solver_)
    solver->messageHandler()->setLogLevel(solver_->messageHandler()->logLevel()) ;

  if (ourSolver_&&deleteSolver) delete solver_ ;
  solver_ = solver;
  solver = NULL ;
  ourSolver_ = true ;
/*
  Basis information is solver-specific.
*/
  if (emptyWarmStart_)
  { delete emptyWarmStart_  ;
    emptyWarmStart_ = 0 ; }
/*
  Initialize integer variable vector.
*/
  numberIntegers_=0;
  int numberColumns = solver_->getNumCols();
  int iColumn;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if( solver_->isInteger(iColumn)) 
      numberIntegers_++;
  }
  delete [] integerVariable_;
  if (numberIntegers_) {
    integerVariable_ = new int [numberIntegers_];
    numberIntegers_=0;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if( solver_->isInteger(iColumn)) 
	integerVariable_[numberIntegers_++]=iColumn;
    }
  } else {
    integerVariable_ = NULL;
  }

  return ;
}

// Copy constructor.

CbcModel::CbcModel(const CbcModel & rhs, bool noTree)
:
  continuousSolver_(NULL),
  referenceSolver_(NULL),
  defaultHandler_(rhs.defaultHandler_),
  emptyWarmStart_(NULL),
  bestObjective_(rhs.bestObjective_),
  bestPossibleObjective_(rhs.bestPossibleObjective_),
  sumChangeObjective1_(rhs.sumChangeObjective1_),
  sumChangeObjective2_(rhs.sumChangeObjective2_),
  minimumDrop_(rhs.minimumDrop_),
  numberSolutions_(rhs.numberSolutions_),
  stateOfSearch_(rhs.stateOfSearch_),
  numberHeuristicSolutions_(rhs.numberHeuristicSolutions_),
  numberNodes_(rhs.numberNodes_),
  numberNodes2_(rhs.numberNodes2_),
  numberIterations_(rhs.numberIterations_),
  status_(rhs.status_),
  secondaryStatus_(rhs.secondaryStatus_),
  specialOptions_(rhs.specialOptions_),
  subTreeModel_(rhs.subTreeModel_),
  numberStoppedSubTrees_(rhs.numberStoppedSubTrees_),
  mutex_(NULL),
  presolve_(rhs.presolve_),
  numberStrong_(rhs.numberStrong_),
  numberBeforeTrust_(rhs.numberBeforeTrust_),
  numberPenalties_(rhs.numberPenalties_),
  penaltyScaleFactor_(rhs.penaltyScaleFactor_),
  numberAnalyzeIterations_(rhs.numberAnalyzeIterations_),
  analyzeResults_(NULL),
  numberInfeasibleNodes_(rhs.numberInfeasibleNodes_),
  problemType_(rhs.problemType_),
  printFrequency_(rhs.printFrequency_),
  howOftenGlobalScan_(rhs.howOftenGlobalScan_),
  numberGlobalViolations_(rhs.numberGlobalViolations_),
  continuousObjective_(rhs.continuousObjective_),
  originalContinuousObjective_(rhs.originalContinuousObjective_),
  continuousInfeasibilities_(rhs.continuousInfeasibilities_),
  maximumCutPassesAtRoot_(rhs.maximumCutPassesAtRoot_),
  maximumCutPasses_( rhs.maximumCutPasses_),
  preferredWay_(rhs.preferredWay_),
  currentPassNumber_(rhs.currentPassNumber_),
  maximumWhich_(rhs.maximumWhich_),
  whichGenerator_(NULL),
  maximumStatistics_(0),
  statistics_(NULL),
  maximumDepthActual_(0),
  numberDJFixed_(0.0),
  probingInfo_(NULL),
  numberFixedAtRoot_(rhs.numberFixedAtRoot_),
  numberFixedNow_(rhs.numberFixedNow_),
  stoppedOnGap_(rhs.stoppedOnGap_),
  eventHappened_(rhs.eventHappened_),
  numberLongStrong_(rhs.numberLongStrong_),
  numberOldActiveCuts_(rhs.numberOldActiveCuts_),
  numberNewCuts_(rhs.numberNewCuts_),
  sizeMiniTree_(rhs.sizeMiniTree_),
  searchStrategy_(rhs.searchStrategy_),
  numberStrongIterations_(rhs.numberStrongIterations_),
  resolveAfterTakeOffCuts_(rhs.resolveAfterTakeOffCuts_),
#if NEW_UPDATE_OBJECT>1
  numberUpdateItems_(rhs.numberUpdateItems_),
  maximumNumberUpdateItems_(rhs.maximumNumberUpdateItems_),
  updateItems_(NULL),
#endif
  numberThreads_(rhs.numberThreads_),
  threadMode_(rhs.threadMode_)
{
  memcpy(intParam_,rhs.intParam_,sizeof(intParam_));
  memcpy(dblParam_,rhs.dblParam_,sizeof(dblParam_));
  strongInfo_[0]=rhs.strongInfo_[0];
  strongInfo_[1]=rhs.strongInfo_[1];
  strongInfo_[2]=rhs.strongInfo_[2];
  solverCharacteristics_ = NULL;
  if (rhs.emptyWarmStart_) emptyWarmStart_ = rhs.emptyWarmStart_->clone() ;
  if (defaultHandler_) {
    handler_ = new CoinMessageHandler();
    handler_->setLogLevel(2);
  } else {
    handler_ = rhs.handler_;
  }
  messageHandler()->setLogLevel(rhs.messageHandler()->logLevel());
  numberCutGenerators_ = rhs.numberCutGenerators_;
  if (numberCutGenerators_) {
    generator_ = new CbcCutGenerator * [numberCutGenerators_];
    virginGenerator_ = new CbcCutGenerator * [numberCutGenerators_];
    int i;
    for (i=0;i<numberCutGenerators_;i++) {
      generator_[i]=new CbcCutGenerator(*rhs.generator_[i]);
      virginGenerator_[i]=new CbcCutGenerator(*rhs.virginGenerator_[i]);
    }
  } else {
    generator_=NULL;
    virginGenerator_=NULL;
  }
  if (!noTree)
    globalCuts_ = rhs.globalCuts_;
  numberHeuristics_ = rhs.numberHeuristics_;
  if (numberHeuristics_) {
    heuristic_ = new CbcHeuristic * [numberHeuristics_];
    int i;
    for (i=0;i<numberHeuristics_;i++) {
      heuristic_[i]=rhs.heuristic_[i]->clone();
    }
  } else {
    heuristic_=NULL;
  }
  lastHeuristic_ = NULL;
  if (rhs.eventHandler_)
  { eventHandler_ = new CbcEventHandler(*rhs.eventHandler_) ; }
  else
  { eventHandler_ = NULL ; }
  ownObjects_ = rhs.ownObjects_;
  if (ownObjects_) {
    numberObjects_=rhs.numberObjects_;
    if (numberObjects_) {
      object_ = new OsiObject * [numberObjects_];
      int i;
      for (i=0;i<numberObjects_;i++) {
	object_[i]=(rhs.object_[i])->clone();
	CbcObject * obj = dynamic_cast <CbcObject *>(object_[i]) ;
	// Could be OsiObjects
	if (obj)
	  obj->setModel(this);
      }
    } else {
      object_=NULL;
    }
  } else {
    // assume will be redone
    numberObjects_=0;
    object_=NULL;
  }
  if (rhs.referenceSolver_)
    referenceSolver_ = rhs.referenceSolver_->clone();
  else
    referenceSolver_=NULL;
  if (!noTree||!rhs.continuousSolver_)
    solver_ = rhs.solver_->clone();
  else
    solver_ = rhs.continuousSolver_->clone();
  if (rhs.originalColumns_) {
    int numberColumns = solver_->getNumCols();
    originalColumns_= new int [numberColumns];
    memcpy(originalColumns_,rhs.originalColumns_,numberColumns*sizeof(int));
  } else {
    originalColumns_=NULL;
  }
#if NEW_UPDATE_OBJECT>1
  if (maximumNumberUpdateItems_) {
    updateItems_ = new CbcObjectUpdateData [maximumNumberUpdateItems_];
    for (int i=0;i<maximumNumberUpdateItems_;i++)
      updateItems_[i] = rhs.updateItems_[i];
  }
#endif
  if (maximumWhich_&&rhs.whichGenerator_)
    whichGenerator_ = CoinCopyOfArray(rhs.whichGenerator_,maximumWhich_);
  nodeCompare_=rhs.nodeCompare_->clone();
  problemFeasibility_=rhs.problemFeasibility_->clone();
  tree_= rhs.tree_->clone();
  if (rhs.branchingMethod_)
    branchingMethod_=rhs.branchingMethod_->clone();
  else
    branchingMethod_=NULL;
  if (rhs.cutModifier_)
    cutModifier_=rhs.cutModifier_->clone();
  else
    cutModifier_=NULL;
  cbcColLower_ = NULL;
  cbcColUpper_ = NULL;
  cbcRowLower_ = NULL;
  cbcRowUpper_ = NULL;
  cbcColSolution_ = NULL;
  cbcRowPrice_ = NULL;
  cbcReducedCost_ = NULL;
  cbcRowActivity_ = NULL;
  if (rhs.strategy_)
    strategy_=rhs.strategy_->clone();
  else
    strategy_=NULL;
  parentModel_=rhs.parentModel_;
  appData_=rhs.appData_;
  messages_ = rhs.messages_;
  ourSolver_ = true ;
  messageHandler()->setLogLevel(rhs.messageHandler()->logLevel());
  numberIntegers_=rhs.numberIntegers_;
  if (numberIntegers_) {
    integerVariable_ = new int [numberIntegers_];
    memcpy(integerVariable_,rhs.integerVariable_,numberIntegers_*sizeof(int));
    integerInfo_ = CoinCopyOfArray(rhs.integerInfo_,solver_->getNumCols());
  } else {
    integerVariable_ = NULL;
    integerInfo_=NULL;
  }
  if (rhs.hotstartSolution_) {
    int numberColumns = solver_->getNumCols();
    hotstartSolution_ = CoinCopyOfArray(rhs.hotstartSolution_,numberColumns);
    hotstartPriorities_ = CoinCopyOfArray(rhs.hotstartPriorities_,numberColumns);
  } else {
    hotstartSolution_ = NULL;
    hotstartPriorities_ =NULL;
  }
  if (rhs.bestSolution_&&!noTree) {
    int numberColumns = solver_->getNumCols();
    bestSolution_ = new double[numberColumns];
    memcpy(bestSolution_,rhs.bestSolution_,numberColumns*sizeof(double));
  } else {
    bestSolution_=NULL;
  }
  if (!noTree) {
    int numberColumns = solver_->getNumCols();
    currentSolution_ = CoinCopyOfArray(rhs.currentSolution_,numberColumns);
    continuousSolution_ = CoinCopyOfArray(rhs.continuousSolution_,numberColumns);
    usedInSolution_ = CoinCopyOfArray(rhs.usedInSolution_,numberColumns);
  } else {
    currentSolution_=NULL;
    continuousSolution_=NULL;
    usedInSolution_=NULL;
  }
  testSolution_=currentSolution_;
  numberRowsAtContinuous_ = rhs.numberRowsAtContinuous_;
  maximumNumberCuts_=rhs.maximumNumberCuts_;
  phase_ = rhs.phase_;
  currentNumberCuts_=rhs.currentNumberCuts_;
  maximumDepth_= rhs.maximumDepth_;
  if (noTree) {
    bestObjective_ = COIN_DBL_MAX;
    numberSolutions_ =0;
    stateOfSearch_= 0;
    numberHeuristicSolutions_=0;
    numberNodes_=0;
    numberNodes2_=0;
    numberIterations_=0;
    status_=0;
    subTreeModel_=NULL;
    numberStoppedSubTrees_=0;
    continuousObjective_=COIN_DBL_MAX;
    originalContinuousObjective_=COIN_DBL_MAX;
    continuousInfeasibilities_=COIN_INT_MAX;
    maximumNumberCuts_=0;
    tree_->cleanTree(this,-COIN_DBL_MAX,bestPossibleObjective_);
    bestPossibleObjective_ = COIN_DBL_MAX;
  }
  // These are only used as temporary arrays so need not be filled
  if (maximumNumberCuts_) {
    addedCuts_ = new CbcCountRowCut * [maximumNumberCuts_];
  } else {
    addedCuts_ = NULL;
  }
  nextRowCut_ = NULL;
  currentNode_ = NULL;
  if (maximumDepth_)
    walkback_ = new CbcNodeInfo * [maximumDepth_];
  else
    walkback_ = NULL;
  synchronizeModel();
}
  
// Assignment operator 
CbcModel & 
CbcModel::operator=(const CbcModel& rhs)
{
  if (this!=&rhs) {
    if (ourSolver_) {
      delete solver_;
      solver_=NULL;
    }
    gutsOfDestructor();
    if (defaultHandler_)
    { delete handler_;
      handler_ = NULL; }
    defaultHandler_ = rhs.defaultHandler_;
    if (defaultHandler_)
    { handler_ = new CoinMessageHandler();
      handler_->setLogLevel(2); }
    else
    { handler_ = rhs.handler_; }
    messages_ = rhs.messages_;
    messageHandler()->setLogLevel(rhs.messageHandler()->logLevel());
    if (rhs.solver_)
    { solver_ = rhs.solver_->clone() ; }
    else
    { solver_ = 0 ; }
    ourSolver_ = true ;
    delete continuousSolver_ ;
    if (rhs.continuousSolver_)
    { continuousSolver_ = rhs.continuousSolver_->clone() ; }
    else
    { continuousSolver_ = 0 ; }
    delete referenceSolver_;
    if (rhs.referenceSolver_)
    { referenceSolver_ = rhs.referenceSolver_->clone() ; }
    else
    { referenceSolver_ = NULL ; }

    delete emptyWarmStart_ ;
    if (rhs.emptyWarmStart_)
    { emptyWarmStart_ = rhs.emptyWarmStart_->clone() ; }
    else
    { emptyWarmStart_ = 0 ; }

    bestObjective_ = rhs.bestObjective_;
    bestPossibleObjective_=rhs.bestPossibleObjective_;
    sumChangeObjective1_=rhs.sumChangeObjective1_;
    sumChangeObjective2_=rhs.sumChangeObjective2_;
    delete [] bestSolution_;
    if (rhs.bestSolution_) {
      int numberColumns = rhs.getNumCols();
      bestSolution_ = new double[numberColumns];
      memcpy(bestSolution_,rhs.bestSolution_,numberColumns*sizeof(double));
    } else {
      bestSolution_=NULL;
    }
    int numberColumns = rhs.getNumCols();
    currentSolution_ = CoinCopyOfArray(rhs.currentSolution_,numberColumns);
    continuousSolution_ = CoinCopyOfArray(rhs.continuousSolution_,numberColumns);
    usedInSolution_ = CoinCopyOfArray(rhs.usedInSolution_,numberColumns);
    testSolution_=currentSolution_;
    minimumDrop_ = rhs.minimumDrop_;
    numberSolutions_=rhs.numberSolutions_;
    stateOfSearch_= rhs.stateOfSearch_;
    numberHeuristicSolutions_=rhs.numberHeuristicSolutions_;
    numberNodes_ = rhs.numberNodes_;
    numberNodes2_ = rhs.numberNodes2_;
    numberIterations_ = rhs.numberIterations_;
    status_ = rhs.status_;
    secondaryStatus_ = rhs.secondaryStatus_;
    specialOptions_ = rhs.specialOptions_;
    subTreeModel_ = rhs.subTreeModel_;
    numberStoppedSubTrees_ = rhs.numberStoppedSubTrees_;
    mutex_ = NULL;
    presolve_ = rhs.presolve_;
    numberStrong_ = rhs.numberStrong_;
    numberBeforeTrust_ = rhs.numberBeforeTrust_;
    numberPenalties_ = rhs.numberPenalties_;
    penaltyScaleFactor_ = rhs.penaltyScaleFactor_;
    numberAnalyzeIterations_ = rhs.numberAnalyzeIterations_;
    delete [] analyzeResults_;
    analyzeResults_ = NULL;
    numberInfeasibleNodes_ = rhs.numberInfeasibleNodes_;
    problemType_ = rhs.problemType_;
    printFrequency_ = rhs.printFrequency_;
    howOftenGlobalScan_=rhs.howOftenGlobalScan_;
    numberGlobalViolations_=rhs.numberGlobalViolations_;
    continuousObjective_=rhs.continuousObjective_;
    originalContinuousObjective_ = rhs.originalContinuousObjective_;
    continuousInfeasibilities_ = rhs.continuousInfeasibilities_;
    maximumCutPassesAtRoot_ = rhs.maximumCutPassesAtRoot_;
    maximumCutPasses_ = rhs.maximumCutPasses_;
    preferredWay_ = rhs.preferredWay_;
    currentPassNumber_ = rhs.currentPassNumber_;
    memcpy(intParam_,rhs.intParam_,sizeof(intParam_));
    memcpy(dblParam_,rhs.dblParam_,sizeof(dblParam_));
    globalCuts_ = rhs.globalCuts_;
    int i;
    for (i=0;i<numberCutGenerators_;i++) {
      delete generator_[i];
      delete virginGenerator_[i];
    }
    delete [] generator_;
    delete [] virginGenerator_;
    delete [] heuristic_;
    maximumWhich_ = rhs.maximumWhich_;
    delete [] whichGenerator_;
    whichGenerator_ = NULL;
    if (maximumWhich_&&rhs.whichGenerator_)
      whichGenerator_ = CoinCopyOfArray(rhs.whichGenerator_,maximumWhich_);
    for (i=0;i<maximumStatistics_;i++)
      delete statistics_[i];
    delete [] statistics_;
    maximumStatistics_ = 0;
    statistics_ = NULL;
    delete probingInfo_;
    probingInfo_=NULL;
    numberFixedAtRoot_ = rhs.numberFixedAtRoot_;
    numberFixedNow_ = rhs.numberFixedNow_;
    stoppedOnGap_ = rhs.stoppedOnGap_;
    eventHappened_ = rhs.eventHappened_;
    numberLongStrong_ = rhs.numberLongStrong_;
    numberOldActiveCuts_ = rhs.numberOldActiveCuts_;
    numberNewCuts_ = rhs.numberNewCuts_;
    resolveAfterTakeOffCuts_=rhs.resolveAfterTakeOffCuts_;
#if NEW_UPDATE_OBJECT>1
    numberUpdateItems_ = rhs.numberUpdateItems_;
    maximumNumberUpdateItems_ = rhs.maximumNumberUpdateItems_;
    delete [] updateItems_;
    if (maximumNumberUpdateItems_) {
      updateItems_ = new CbcObjectUpdateData [maximumNumberUpdateItems_];
      for (i=0;i<maximumNumberUpdateItems_;i++)
	updateItems_[i] = rhs.updateItems_[i];
    } else {
      updateItems_ = NULL;
    }
#endif
    numberThreads_ = rhs.numberThreads_;
    threadMode_ = rhs.threadMode_;
    sizeMiniTree_ = rhs.sizeMiniTree_;
    searchStrategy_ = rhs.searchStrategy_;
    numberStrongIterations_ = rhs.numberStrongIterations_;
    strongInfo_[0]=rhs.strongInfo_[0];
    strongInfo_[1]=rhs.strongInfo_[1];
    strongInfo_[2]=rhs.strongInfo_[2];
    solverCharacteristics_ = NULL;
    lastHeuristic_ = NULL;
    numberCutGenerators_ = rhs.numberCutGenerators_;
    if (numberCutGenerators_) {
      generator_ = new CbcCutGenerator * [numberCutGenerators_];
      virginGenerator_ = new CbcCutGenerator * [numberCutGenerators_];
      int i;
      for (i=0;i<numberCutGenerators_;i++) {
	generator_[i]=new CbcCutGenerator(*rhs.generator_[i]);
	virginGenerator_[i]=new CbcCutGenerator(*rhs.virginGenerator_[i]);
      }
    } else {
      generator_=NULL;
      virginGenerator_=NULL;
    }
    numberHeuristics_ = rhs.numberHeuristics_;
    if (numberHeuristics_) {
      heuristic_ = new CbcHeuristic * [numberHeuristics_];
      memcpy(heuristic_,rhs.heuristic_,
	     numberHeuristics_*sizeof(CbcHeuristic *));
    } else {
      heuristic_=NULL;
    }
    lastHeuristic_ = NULL;
    if (eventHandler_)
      delete eventHandler_ ;
    if (rhs.eventHandler_)
    { eventHandler_ = new CbcEventHandler(*rhs.eventHandler_) ; }
    else
    { eventHandler_ = NULL ; }
    if (ownObjects_) {
      for (i=0;i<numberObjects_;i++)
	delete object_[i];
      delete [] object_;
      numberObjects_=rhs.numberObjects_;
      if (numberObjects_) {
	object_ = new OsiObject * [numberObjects_];
	int i;
	for (i=0;i<numberObjects_;i++) 
	  object_[i]=(rhs.object_[i])->clone();
      } else {
	object_=NULL;
    }
    } else {
      // assume will be redone
      numberObjects_=0;
      object_=NULL;
    }
    delete [] originalColumns_;
    if (rhs.originalColumns_) {
      int numberColumns = rhs.getNumCols();
      originalColumns_= new int [numberColumns];
      memcpy(originalColumns_,rhs.originalColumns_,numberColumns*sizeof(int));
    } else {
      originalColumns_=NULL;
    }
    nodeCompare_=rhs.nodeCompare_->clone();
    problemFeasibility_=rhs.problemFeasibility_->clone();
    delete tree_;
    tree_= rhs.tree_->clone();
    if (rhs.branchingMethod_)
      branchingMethod_=rhs.branchingMethod_->clone();
    else
      branchingMethod_=NULL;
    if (rhs.cutModifier_)
      cutModifier_=rhs.cutModifier_->clone();
    else
      cutModifier_=NULL;
    delete strategy_;
    if (rhs.strategy_)
      strategy_=rhs.strategy_->clone();
    else
      strategy_=NULL;
    parentModel_=rhs.parentModel_;
    appData_=rhs.appData_;

    delete [] integerVariable_;
    numberIntegers_=rhs.numberIntegers_;
    if (numberIntegers_) {
      integerVariable_ = new int [numberIntegers_];
      memcpy(integerVariable_,rhs.integerVariable_,
	     numberIntegers_*sizeof(int));
      integerInfo_ = CoinCopyOfArray(rhs.integerInfo_,rhs.getNumCols());
    } else {
      integerVariable_ = NULL;
      integerInfo_=NULL;
    }
    if (rhs.hotstartSolution_) {
      int numberColumns = solver_->getNumCols();
      hotstartSolution_ = CoinCopyOfArray(rhs.hotstartSolution_,numberColumns);
      hotstartPriorities_ = CoinCopyOfArray(rhs.hotstartPriorities_,numberColumns);
    } else {
      hotstartSolution_ = NULL;
      hotstartPriorities_ =NULL;
    }
    numberRowsAtContinuous_ = rhs.numberRowsAtContinuous_;
    maximumNumberCuts_=rhs.maximumNumberCuts_;
    phase_ = rhs.phase_;
    currentNumberCuts_=rhs.currentNumberCuts_;
    maximumDepth_= rhs.maximumDepth_;
    delete [] addedCuts_;
    delete [] walkback_;
    // These are only used as temporary arrays so need not be filled
    if (maximumNumberCuts_) {
      addedCuts_ = new CbcCountRowCut * [maximumNumberCuts_];
    } else {
      addedCuts_ = NULL;
    }
    nextRowCut_ = NULL;
    currentNode_ = NULL;
    if (maximumDepth_)
      walkback_ = new CbcNodeInfo * [maximumDepth_];
    else
      walkback_ = NULL;
    synchronizeModel();
    cbcColLower_ = NULL;
    cbcColUpper_ = NULL;
    cbcRowLower_ = NULL;
    cbcRowUpper_ = NULL;
    cbcColSolution_ = NULL;
    cbcRowPrice_ = NULL;
    cbcReducedCost_ = NULL;
    cbcRowActivity_ = NULL;
  }
  return *this;
}
  
// Destructor 
CbcModel::~CbcModel ()
{
  if (defaultHandler_) {
    delete handler_;
    handler_ = NULL;
  }
  delete tree_;
  tree_=NULL;
  if (ourSolver_) {
    delete solver_;
    solver_ = NULL;
  }
  gutsOfDestructor();
  delete eventHandler_ ;
  eventHandler_ = NULL ;
}
// Clears out as much as possible (except solver)
void 
CbcModel::gutsOfDestructor()
{
  delete referenceSolver_;
  referenceSolver_=NULL;
  int i;
  for (i=0;i<numberCutGenerators_;i++) {
    delete generator_[i];
    delete virginGenerator_[i];
  }
  delete [] generator_;
  delete [] virginGenerator_;
  generator_=NULL;
  virginGenerator_=NULL;
  for (i=0;i<numberHeuristics_;i++)
    delete heuristic_[i];
  delete [] heuristic_;
  heuristic_=NULL;
  delete nodeCompare_;
  nodeCompare_=NULL;
  delete problemFeasibility_;
  problemFeasibility_=NULL;
  delete [] originalColumns_;
  originalColumns_=NULL;
  delete strategy_;
#if NEW_UPDATE_OBJECT>1
  delete [] updateItems_;
  updateItems_=NULL;
  numberUpdateItems_=0;
  maximumNumberUpdateItems_=0;
#endif
  gutsOfDestructor2();
}
// Clears out enough to reset CbcModel
void 
CbcModel::gutsOfDestructor2()
{
  delete [] integerInfo_;
  integerInfo_=NULL;
  delete [] integerVariable_;
  integerVariable_=NULL;
  int i;
  if (ownObjects_) {
    for (i=0;i<numberObjects_;i++)
      delete object_[i];
    delete [] object_;
  }
  ownObjects_=true;
  object_=NULL;
  numberIntegers_=0;
  numberObjects_=0;
  // Below here is whatever consensus is
  ourSolver_=true;
  delete branchingMethod_;
  branchingMethod_=NULL;
  delete cutModifier_;
  cutModifier_=NULL;
  resetModel();
}
// Clears out enough to reset CbcModel
void 
CbcModel::resetModel()
{
  delete emptyWarmStart_ ;
  emptyWarmStart_ =NULL;
  delete continuousSolver_;
  continuousSolver_=NULL;
  delete [] bestSolution_;
  bestSolution_=NULL;
  delete [] currentSolution_;
  currentSolution_=NULL;
  delete [] continuousSolution_;
  continuousSolution_=NULL;
  solverCharacteristics_=NULL;
  delete [] usedInSolution_;
  usedInSolution_ = NULL;
  testSolution_=NULL;
  lastHeuristic_ = NULL;
  delete [] addedCuts_;
  addedCuts_=NULL;
  nextRowCut_ = NULL;
  currentNode_ = NULL;
  delete [] walkback_;
  walkback_=NULL;
  delete [] whichGenerator_;
  whichGenerator_ = NULL;
  for (int i=0;i<maximumStatistics_;i++)
    delete statistics_[i];
  delete [] statistics_;
  statistics_=NULL;
  maximumDepthActual_ = 0;
  numberDJFixed_ =0.0;
  delete probingInfo_;
  probingInfo_ = NULL;
  maximumStatistics_=0;
  delete [] analyzeResults_;
  analyzeResults_=NULL;
  bestObjective_=COIN_DBL_MAX;
  bestPossibleObjective_=COIN_DBL_MAX;
  sumChangeObjective1_=0.0;
  sumChangeObjective2_=0.0;
  numberSolutions_=0;
  stateOfSearch_=0;
  delete [] hotstartSolution_;
  hotstartSolution_=NULL;
  delete [] hotstartPriorities_;
  hotstartPriorities_=NULL;
  numberHeuristicSolutions_=0;
  numberNodes_=0;
  numberNodes2_=0;
  numberIterations_=0;
  status_=-1;
  secondaryStatus_=-1;
  maximumNumberCuts_=0;
  phase_=0;
  currentNumberCuts_=0;
  maximumDepth_=0;
  nextRowCut_=NULL;
  currentNode_=NULL;
  // clear out tree
  if (tree_&&tree_->size())
    tree_->cleanTree(this, -1.0e100,bestPossibleObjective_) ;
  subTreeModel_=NULL;
  numberStoppedSubTrees_=0;
  numberInfeasibleNodes_=0;
  numberGlobalViolations_=0;
  continuousObjective_=0.0;
  originalContinuousObjective_=0.0;
  continuousInfeasibilities_=0;
  numberFixedAtRoot_=0;
  numberFixedNow_=0;
  stoppedOnGap_=false;
  eventHappened_=false;
  numberLongStrong_=0;
  numberOldActiveCuts_=0;
  numberNewCuts_=0;
  searchStrategy_=-1;
  numberStrongIterations_=0;
  // Parameters which need to be reset
  setCutoff(COIN_DBL_MAX);
  dblParam_[CbcCutoffIncrement] = 1e-5;
  dblParam_[CbcCurrentCutoff] = 1.0e100;
  dblParam_[CbcCurrentObjectiveValue] = 1.0e100;
  dblParam_[CbcCurrentMinimizationObjectiveValue] = 1.0e100;
}
// Move status, nodes etc etc across
void 
CbcModel::moveInfo(const CbcModel & rhs)
{
  bestObjective_ = rhs.bestObjective_;
  bestPossibleObjective_=rhs.bestPossibleObjective_;
  numberSolutions_=rhs.numberSolutions_;
  numberHeuristicSolutions_=rhs.numberHeuristicSolutions_;
  numberNodes_ = rhs.numberNodes_;
  numberNodes2_ = rhs.numberNodes2_;
  numberIterations_ = rhs.numberIterations_;
  status_ = rhs.status_;
  secondaryStatus_ = rhs.secondaryStatus_;
  numberStoppedSubTrees_ = rhs.numberStoppedSubTrees_;
  numberInfeasibleNodes_ = rhs.numberInfeasibleNodes_;
  continuousObjective_=rhs.continuousObjective_;
  originalContinuousObjective_ = rhs.originalContinuousObjective_;
  continuousInfeasibilities_ = rhs.continuousInfeasibilities_;
  numberFixedAtRoot_ = rhs.numberFixedAtRoot_;
  numberFixedNow_ = rhs.numberFixedNow_;
  stoppedOnGap_ = rhs.stoppedOnGap_;
  eventHappened_ = rhs.eventHappened_;
  numberLongStrong_ = rhs.numberLongStrong_;
  numberStrongIterations_ = rhs.numberStrongIterations_;
  strongInfo_[0]=rhs.strongInfo_[0];
  strongInfo_[1]=rhs.strongInfo_[1];
  strongInfo_[2]=rhs.strongInfo_[2];
  numberRowsAtContinuous_ = rhs.numberRowsAtContinuous_;
  maximumDepth_= rhs.maximumDepth_;
}
// Save a copy of the current solver so can be reset to
void 
CbcModel::saveReferenceSolver()
{
  delete referenceSolver_;
  referenceSolver_= solver_->clone();
}

// Uses a copy of reference solver to be current solver
void 
CbcModel::resetToReferenceSolver()
{
  delete solver_;
  solver_ = referenceSolver_->clone();
  // clear many things
  gutsOfDestructor2();
  // Reset cutoff
  // Solvers know about direction
  double direction = solver_->getObjSense();
  double value;
  solver_->getDblParam(OsiDualObjectiveLimit,value); 
  setCutoff(value*direction);
}

// Are there a numerical difficulties?
bool 
CbcModel::isAbandoned() const
{
  return status_ == 2;
}
// Is optimality proven?
bool 
CbcModel::isProvenOptimal() const
{
  if (!status_ && bestObjective_<1.0e30)
    return true;
  else
    return false;
}
// Is  infeasiblity proven (or none better than cutoff)?
bool 
CbcModel::isProvenInfeasible() const
{
  if (!status_ && bestObjective_>=1.0e30)
    return true;
  else
    return false;
}
// Was continuous solution unbounded
bool 
CbcModel::isContinuousUnbounded() const
{
  if (!status_ && secondaryStatus_==7)
    return true;
  else
    return false;
}
// Was continuous solution unbounded
bool 
CbcModel::isProvenDualInfeasible() const
{
  if (!status_ && secondaryStatus_==7)
    return true;
  else
    return false;
}
// Node limit reached?
bool 
CbcModel::isNodeLimitReached() const
{
  return numberNodes_ >= intParam_[CbcMaxNumNode];
}
// Time limit reached?
bool 
CbcModel::isSecondsLimitReached() const
{
  if (status_==1&&secondaryStatus_==4)
    return true;
  else
    return false;
}
// Solution limit reached?
bool 
CbcModel::isSolutionLimitReached() const
{
  return numberSolutions_ >= intParam_[CbcMaxNumSol];
}
// Set language
void 
CbcModel::newLanguage(CoinMessages::Language language)
{
  messages_ = CbcMessage(language);
}
void 
CbcModel::setNumberStrong(int number)
{
  if (number<0)
    numberStrong_=0;
   else
    numberStrong_=number;
}
void 
CbcModel::setNumberBeforeTrust(int number)
{
  if (number<-1) {
    numberBeforeTrust_=0;
  } else {
    numberBeforeTrust_=number;
    //numberStrong_ = CoinMax(numberStrong_,1);
  }
}
void 
CbcModel::setNumberPenalties(int number)
{
  if (number<=0) {
    numberPenalties_=0;
  } else {
    numberPenalties_=number;
  }
}
void 
CbcModel::setPenaltyScaleFactor(double value)
{
  if (value<=0) {
    penaltyScaleFactor_=3.0;
  } else {
    penaltyScaleFactor_=value;
  }
}
void 
CbcModel::setHowOftenGlobalScan(int number)
{
  if (number<-1)
    howOftenGlobalScan_=0;
   else
    howOftenGlobalScan_=number;
}

// Add one generator
void 
CbcModel::addCutGenerator(CglCutGenerator * generator,
			  int howOften, const char * name,
			  bool normal, bool atSolution,
			  bool whenInfeasible,int howOftenInSub,
			  int whatDepth, int whatDepthInSub)
{
  CbcCutGenerator ** temp = generator_;
  generator_ = new CbcCutGenerator * [numberCutGenerators_+1];
  memcpy(generator_,temp,numberCutGenerators_*sizeof(CbcCutGenerator *));
  delete[] temp ;
  generator_[numberCutGenerators_]= 
    new CbcCutGenerator(this,generator, howOften, name,
			normal,atSolution,whenInfeasible,howOftenInSub,
			whatDepth, whatDepthInSub);
  // and before any cahnges
  temp = virginGenerator_;
  virginGenerator_ = new CbcCutGenerator * [numberCutGenerators_+1];
  memcpy(virginGenerator_,temp,numberCutGenerators_*sizeof(CbcCutGenerator *));
  delete[] temp ;
  virginGenerator_[numberCutGenerators_++]= 
    new CbcCutGenerator(this,generator, howOften, name,
			normal,atSolution,whenInfeasible,howOftenInSub,
			whatDepth, whatDepthInSub);
							  
}
// Add one heuristic
void 
CbcModel::addHeuristic(CbcHeuristic * generator, const char *name)
{
  CbcHeuristic ** temp = heuristic_;
  heuristic_ = new CbcHeuristic * [numberHeuristics_+1];
  memcpy(heuristic_,temp,numberHeuristics_*sizeof(CbcHeuristic *));
  delete [] temp;
  heuristic_[numberHeuristics_]=generator->clone();
  if (name)
  { heuristic_[numberHeuristics_]->setHeuristicName(name) ; }
  numberHeuristics_++ ;
}

/*
  The last subproblem handled by the solver is not necessarily related to the
  one being recreated, so the first action is to remove all cuts from the
  constraint system.  Next, traverse the tree from node to the root to
  determine the basis size required for this subproblem and create an empty
  basis with the right capacity.  Finally, traverse the tree from root to
  node, adjusting bounds in the constraint system, adjusting the basis, and
  collecting the cuts that must be added to the constraint system.
  applyToModel does the heavy lifting.

  addCuts1 is used in contexts where all that's desired is the list of cuts:
  the node is already fathomed, and we're collecting cuts so that we can
  adjust reference counts as we prune nodes. Arguably the two functions
  should be separated. The culprit is applyToModel, which performs cut
  collection and model adjustment.

  Certainly in the contexts where all we need is a list of cuts, there's no
  point in passing in a valid basis --- an empty basis will do just fine.
*/
void CbcModel::addCuts1 (CbcNode * node, CoinWarmStartBasis *&lastws)
{ int i;
  int nNode=0;
  int numberColumns = getNumCols();
  CbcNodeInfo * nodeInfo = node->nodeInfo();

/*
  Remove all cuts from the constraint system.
  (original comment includes ``see note below for later efficiency'', but
  the reference isn't clear to me).
*/
  int currentNumberCuts = solver_->getNumRows()-numberRowsAtContinuous_;
  int *which = new int[currentNumberCuts];
  for (i = 0 ; i < currentNumberCuts ; i++)
    which[i] = i+numberRowsAtContinuous_;
  solver_->deleteRows(currentNumberCuts,which);
  delete [] which;
/*
  Accumulate the path from node to the root in walkback_, and accumulate a
  cut count in currentNumberCuts.

  original comment: when working then just unwind until where new node joins
  old node (for cuts?)
*/
  currentNumberCuts = 0;
  while (nodeInfo) {
    //printf("nNode = %d, nodeInfo = %x\n",nNode,nodeInfo);
    walkback_[nNode++]=nodeInfo;
    currentNumberCuts += nodeInfo->numberCuts() ;
    nodeInfo = nodeInfo->parent() ;
    if (nNode==maximumDepth_) {
      maximumDepth_ *= 2;
      CbcNodeInfo ** temp = new CbcNodeInfo * [maximumDepth_];
      for (i=0;i<nNode;i++) 
	temp[i] = walkback_[i];
      delete [] walkback_;
      walkback_ = temp;
    }
  }
/*
  Create an empty basis with sufficient capacity for the constraint system
  we'll construct: original system plus cuts. Make sure we have capacity to
  record those cuts in addedCuts_.

  The method of adjusting the basis at a FullNodeInfo object (the root, for
  example) is to use a copy constructor to duplicate the basis held in the
  nodeInfo, then resize it and return the new basis object. Guaranteed,
  lastws will point to a different basis when it returns. We pass in a basis
  because we need the parameter to return the allocated basis, and it's an
  easy way to pass in the size. But we take a hit for memory allocation.
*/
  currentNumberCuts_=currentNumberCuts;
  if (currentNumberCuts > maximumNumberCuts_) {
    maximumNumberCuts_ = currentNumberCuts;
    delete [] addedCuts_;
    addedCuts_ = new CbcCountRowCut * [maximumNumberCuts_];
  }
  lastws->setSize(numberColumns,numberRowsAtContinuous_+currentNumberCuts);
/*
  This last bit of code traverses the path collected in walkback_ from the
  root back to node. At the end of the loop,
   * lastws will be an appropriate basis for node;
   * variable bounds in the constraint system will be set to be correct for
     node; and
   * addedCuts_ will be set to a list of cuts that need to be added to the
     constraint system at node.
  applyToModel does all the heavy lifting.
*/
  currentNumberCuts=0;
  while (nNode) {
    --nNode;
    walkback_[nNode]->applyToModel(this,lastws,addedCuts_,currentNumberCuts);
  }
}

/*
  adjustCuts might be a better name: If the node is feasible, we sift through
  the cuts collected by addCuts1, add the ones that are tight and omit the
  ones that are loose. If the node is infeasible, we just adjust the
  reference counts to reflect that we're about to prune this node and its
  descendants.
*/
int CbcModel::addCuts (CbcNode *node, CoinWarmStartBasis *&lastws,bool canFix)
{
/*
  addCuts1 performs step 1 of restoring the subproblem at this node; see the
  comments there.
*/
  addCuts1(node,lastws);
  int i;
  int numberColumns = getNumCols();
  CbcNodeInfo * nodeInfo = node->nodeInfo();
  double cutoff = getCutoff() ;
  int currentNumberCuts=currentNumberCuts_;
  if (canFix) {
    bool feasible=true;
    const double *lower = solver_->getColLower() ;
    const double *upper = solver_->getColUpper() ;
    double * newLower = analyzeResults_;
    double * objLower = newLower+numberIntegers_;
    double * newUpper = objLower+numberIntegers_;
    double * objUpper = newUpper+numberIntegers_;
    int n=0;
    for (i=0;i<numberIntegers_;i++) {
      int iColumn = integerVariable_[i];
      bool changed=false;
      double lo = 0.0;
      double up = 0.0;
      if (objLower[i]>cutoff) {
        lo = lower[iColumn];
        up = upper[iColumn];
        if (lo<newLower[i]) {
          lo = newLower[i];
          solver_->setColLower(iColumn,lo);
          changed=true;
          n++;
        }
        if (objUpper[i]>cutoff) {
          if (up>newUpper[i]) {
            up = newUpper[i];
            solver_->setColUpper(iColumn,up);
            changed=true;
            n++;
          }
        }
      } else if (objUpper[i]>cutoff) {
        lo = lower[iColumn];
        up = upper[iColumn];
        if (up>newUpper[i]) {
          up = newUpper[i];
          solver_->setColUpper(iColumn,up);
          changed=true;
          n++;
        }
      }
      if (changed&&lo>up) {
        feasible=false;
        break;
      }
    }
    if (!feasible) {
      printf("analysis says node infeas\n");
      cutoff=-COIN_DBL_MAX;
    }
  }
/*
  If the node can't be fathomed by bound, reinstall tight cuts in the
  constraint system. Even if there are no cuts, we'll want to set the
  reconstructed basis in the solver.
*/
  if (node->objectiveValue() < cutoff||numberThreads_)
  { 
#   ifdef CBC_CHECK_BASIS
    printf("addCuts: expanded basis; rows %d+%d\n",
	   numberRowsAtContinuous_,currentNumberCuts);
    lastws->print();
#   endif
/*
  Adjust the basis and constraint system so that we retain only active cuts.
  There are three steps:
    1) Scan the basis. Sort the cuts into effective cuts to be kept and
       loose cuts to be dropped.
    2) Drop the loose cuts and resize the basis to fit.
    3) Install the tight cuts in the constraint system (applyRowCuts) and
       and install the basis (setWarmStart).
  Use of compressRows conveys we're compressing the basis and not just
  tweaking the artificialStatus_ array.
*/
    if (currentNumberCuts > 0) {
      int numberToAdd = 0;
      const OsiRowCut **addCuts;
      int numberToDrop = 0 ;
      int *cutsToDrop ;
      addCuts = new const OsiRowCut* [currentNumberCuts];
      cutsToDrop = new int[currentNumberCuts] ;
      assert (currentNumberCuts+numberRowsAtContinuous_<=lastws->getNumArtificial());
      for (i=0;i<currentNumberCuts;i++) {
	CoinWarmStartBasis::Status status = 
	  lastws->getArtifStatus(i+numberRowsAtContinuous_);
	if (addedCuts_[i] &&
	    (status != CoinWarmStartBasis::basic ||
	     addedCuts_[i]->effectiveness()==COIN_DBL_MAX)) {
#	  ifdef CHECK_CUT_COUNTS
	  printf("Using cut %d %x as row %d\n",i,addedCuts_[i],
		 numberRowsAtContinuous_+numberToAdd);
#	  endif
	  addCuts[numberToAdd++] = new OsiRowCut(*addedCuts_[i]);
	} else {
#	  ifdef CHECK_CUT_COUNTS
	  printf("Dropping cut %d %x\n",i,addedCuts_[i]);
#	  endif
	  addedCuts_[i]=NULL;
	  cutsToDrop[numberToDrop++] = numberRowsAtContinuous_+i ;
	}
      }
      int numberRowsNow=numberRowsAtContinuous_+numberToAdd;
      lastws->compressRows(numberToDrop,cutsToDrop) ;
      lastws->resize(numberRowsNow,numberColumns);
      solver_->applyRowCuts(numberToAdd,addCuts);
#     ifdef CBC_CHECK_BASIS
      printf("addCuts: stripped basis; rows %d + %d\n",
	     numberRowsAtContinuous_,numberToAdd);
      lastws->print();
#     endif
      for (i=0;i<numberToAdd;i++)
	delete addCuts[i];
      delete [] addCuts;
      delete [] cutsToDrop ;
    }
/*
  Set the basis in the solver.
*/
    solver_->setWarmStart(lastws);
#if 0
    if ((numberNodes_%printFrequency_)==0) {
      printf("Objective %g, depth %d, unsatisfied %d\n",
	     node->objectiveValue(),
	     node->depth(),node->numberUnsatisfied());
    }
#endif
/*
  Clean up and we're out of here.
*/
    numberNodes_++;
    return 0;
  } 
/*
  This node has been fathomed by bound as we try to revive it out of the live
  set. Adjust the cut reference counts to reflect that we no longer need to
  explore the remaining branch arms, hence they will no longer reference any
  cuts. Cuts whose reference count falls to zero are deleted.  
*/
  else
  { int i;
    if (currentNumberCuts) {
      lockThread();
      int numberLeft = nodeInfo->numberBranchesLeft();
      for (i = 0 ; i < currentNumberCuts ; i++)
	{ if (addedCuts_[i])
	  { if (!addedCuts_[i]->decrement(numberLeft))
	    { delete addedCuts_[i];
	    addedCuts_[i] = NULL; } } }
      unlockThread();
    }
    return 1 ; }
}


/*
  Perform reduced cost fixing on integer variables.

  The variables in question are already nonbasic at bound. We're just nailing
  down the current situation.
*/

int CbcModel::reducedCostFix ()

{
  if(!solverCharacteristics_->reducedCostsAccurate())
    return 0; //NLP
  double cutoff = getCutoff() ;
  double direction = solver_->getObjSense() ;
  double gap = cutoff - solver_->getObjValue()*direction ;
  double tolerance;
  solver_->getDblParam(OsiDualTolerance,tolerance) ;
  if (gap<=0.0)
    return 0;
  gap += 100.0*tolerance;
  double integerTolerance = getDblParam(CbcIntegerTolerance) ;

  const double *lower = solver_->getColLower() ;
  const double *upper = solver_->getColUpper() ;
  const double *solution = solver_->getColSolution() ;
  const double *reducedCost = solver_->getReducedCost() ;

  int numberFixed = 0 ;

# ifdef COIN_HAS_CLP
  OsiClpSolverInterface * clpSolver 
    = dynamic_cast<OsiClpSolverInterface *> (solver_);
  ClpSimplex * clpSimplex=NULL;
  if (clpSolver) 
    clpSimplex = clpSolver->getModelPtr();
# endif
  for (int i = 0 ; i < numberIntegers_ ; i++)
  { int iColumn = integerVariable_[i] ;
    double djValue = direction*reducedCost[iColumn] ;
    if (upper[iColumn]-lower[iColumn] > integerTolerance)
    { if (solution[iColumn] < lower[iColumn]+integerTolerance && djValue > gap)
      { solver_->setColUpper(iColumn,lower[iColumn]) ;
#ifdef COIN_HAS_CLP
      // may just have been fixed before
      if (clpSimplex)
        assert(clpSimplex->getColumnStatus(iColumn)==ClpSimplex::atLowerBound||
	       clpSimplex->getColumnStatus(iColumn)==ClpSimplex::isFixed);
#endif
      numberFixed++ ; }
      else
      if (solution[iColumn] > upper[iColumn]-integerTolerance && -djValue > gap)
      { solver_->setColLower(iColumn,upper[iColumn]) ;
#ifdef COIN_HAS_CLP
      // may just have been fixed before
      if (clpSimplex)
        assert(clpSimplex->getColumnStatus(iColumn)==ClpSimplex::atUpperBound||
	       clpSimplex->getColumnStatus(iColumn)==ClpSimplex::isFixed);
#endif
      numberFixed++ ; } } }
  numberDJFixed_ += numberFixed;
  
  return numberFixed; }

// Collect coding to replace whichGenerator
void
CbcModel::resizeWhichGenerator(int numberNow, int numberAfter)
{
  if (numberAfter > maximumWhich_) {
    maximumWhich_ = CoinMax(maximumWhich_*2+100,numberAfter) ;
    int * temp = new int[2*maximumWhich_] ;
    memcpy(temp,whichGenerator_,numberNow*sizeof(int)) ;
    delete [] whichGenerator_ ;
    whichGenerator_ = temp ;
    memset(whichGenerator_+numberNow,0,(maximumWhich_-numberNow)*sizeof(int));
  }
}

/** Solve the model using cuts

  This version takes off redundant cuts from node.
  Returns true if feasible.

  \todo
  Why do I need to resolve the problem? What has been done between the last
  relaxation and calling solveWithCuts?

  If numberTries == 0 then user did not want any cuts.
*/

bool 
CbcModel::solveWithCuts (OsiCuts &cuts, int numberTries, CbcNode *node)
/*
  Parameters:
    numberTries: (i) the maximum number of iterations for this round of cut
		     generation; if negative then we don't mind if drop is tiny.
    
    cuts:	(o) all cuts generated in this round of cut generation

    node: (i)     So we can update dynamic pseudo costs
*/
			

{
# ifdef COIN_HAS_CLP
  OsiClpSolverInterface * clpSolver 
    = dynamic_cast<OsiClpSolverInterface *> (solver_);
  int saveClpOptions=0;
  if (clpSolver) 
    saveClpOptions = clpSolver->specialOptions();
# endif
#ifdef CBC_THREAD
  CbcModel ** threadModel = NULL;
  pthread_t * threadId = NULL;
  pthread_cond_t condition_main;
  pthread_mutex_t condition_mutex;
  pthread_mutex_t * mutex2 = NULL;
  pthread_cond_t * condition2 = NULL;
  threadStruct * threadInfo = NULL;
  void * saveMutex = NULL;
  if (numberThreads_&&(threadMode_&2)!=0&&!numberNodes_) {
    threadId = new pthread_t [numberThreads_];
    pthread_cond_init(&condition_main,NULL);
    pthread_mutex_init(&condition_mutex,NULL);
    threadModel = new CbcModel * [numberThreads_];
    threadInfo = new threadStruct [numberThreads_+1];
    mutex2 = new pthread_mutex_t [numberThreads_];
    condition2 = new pthread_cond_t [numberThreads_];
    saveMutex = mutex_;
    for (int i=0;i<numberThreads_;i++) {
      pthread_mutex_init(mutex2+i,NULL);
      pthread_cond_init(condition2+i,NULL);
      threadId[i]=0;
      threadModel[i]=new CbcModel;
      threadModel[i]->generator_ = new CbcCutGenerator * [1];
      delete threadModel[i]->solver_;
      threadModel[i]->solver_=NULL;
      threadModel[i]->numberThreads_=numberThreads_;
      mutex_ = (void *) (threadInfo+i);
      threadInfo[i].thisModel=(CbcModel *) threadModel[i];
      threadInfo[i].baseModel=this;
      threadInfo[i].threadIdOfBase=pthread_self();
      threadInfo[i].mutex2=mutex2+i;
      threadInfo[i].condition2=condition2+i;
      threadInfo[i].returnCode=-1;
      pthread_create(threadId+i,NULL,doCutsThread,threadInfo+i);
    }
    // Do a partial one for base model
    threadInfo[numberThreads_].baseModel=this;
    mutex_ = (void *) (threadInfo+numberThreads_);
    threadInfo[numberThreads_].condition2=&condition_main;
    threadInfo[numberThreads_].mutex2=&condition_mutex;
  }
#endif
  bool feasible = true ;
  int lastNumberCuts = 0 ;
  double lastObjective = -1.0e100 ;
  int violated = 0 ;
  int numberRowsAtStart = solver_->getNumRows() ;
  //printf("solver had %d rows\n",numberRowsAtStart);
  int numberColumns = solver_->getNumCols() ;
  CoinBigIndex numberElementsAtStart = solver_->getNumElements();

  numberOldActiveCuts_ = numberRowsAtStart-numberRowsAtContinuous_ ;
  numberNewCuts_ = 0 ;

  bool onOptimalPath = false ;
  const OsiRowCutDebugger *debugger = NULL;
  if ((specialOptions_&1)!=0) {
    /*
      See OsiRowCutDebugger for details. In a nutshell, make sure that current
      variable values do not conflict with a known optimal solution. (Obviously
      this can be fooled when there are multiple solutions.)
    */
    debugger = solver_->getRowCutDebugger() ;
    if (debugger) 
      onOptimalPath = (debugger->onOptimalPath(*solver_)) ;
  }
  OsiCuts slackCuts;
/*
  Resolve the problem. If we've lost feasibility, might as well bail out right
  after the debug stuff. The resolve will also refresh cached copies of the
  solver solution (cbcColLower_, ...) held by CbcModel.
*/
  double objectiveValue = solver_->getObjValue()*solver_->getObjSense();
  if (node)
    objectiveValue= node->objectiveValue();
  int returnCode = resolve(node ? node->nodeInfo() : NULL,1);
#ifdef COIN_DEVELOP
  //if (!solver_->getIterationCount()&&solver_->isProvenOptimal())
  //printf("zero iterations on first solve of branch\n");
#endif
  if (node&&node->nodeInfo()&&!node->nodeInfo()->numberBranchesLeft())
    node->nodeInfo()->allBranchesGone(); // can clean up
  feasible = returnCode  != 0 ;
  if (returnCode<0)
    numberTries=0;
  if (problemFeasibility_->feasible(this,0)<0) {
    feasible=false; // pretend infeasible
  }
  
#if NEW_UPDATE_OBJECT==0
  // Update branching information if wanted
  if(node &&branchingMethod_)
    branchingMethod_->updateInformation(solver_,node);
#elif NEW_UPDATE_OBJECT<2
  // Update branching information if wanted
  if(node &&branchingMethod_) {
    OsiBranchingObject * bobj = node->modifiableBranchingObject();
    CbcBranchingObject * cbcobj = dynamic_cast<CbcBranchingObject *> (bobj);
    if (cbcobj) {
      CbcObject * object = cbcobj->object();
      CbcObjectUpdateData update = object->createUpdateInformation(solver_,node,cbcobj);
      object->updateInformation(update);
    } else {
      branchingMethod_->updateInformation(solver_,node);
    }
  }
#else
  // Update branching information if wanted
  if(node &&branchingMethod_) {
    OsiBranchingObject * bobj = node->modifiableBranchingObject();
    CbcBranchingObject * cbcobj = dynamic_cast<CbcBranchingObject *> (bobj);
    if (cbcobj) {
      CbcObject * object = cbcobj->object();
      CbcObjectUpdateData update = object->createUpdateInformation(solver_,node,cbcobj);
      // have to compute object number as not saved
      CbcSimpleInteger * simpleObject =
	  dynamic_cast <CbcSimpleInteger *>(object) ;
      int iObject;
      int iColumn = simpleObject->columnNumber();
      for (iObject = 0 ; iObject < numberObjects_ ; iObject++) {
	simpleObject =
	  dynamic_cast <CbcSimpleInteger *>(object_[iObject]) ;
	if (simpleObject->columnNumber()==iColumn)
	  break;
      }
      assert (iObject<numberObjects_);
      update.objectNumber_ = iObject;
      addUpdateInformation(update);
    } else {
      OsiIntegerBranchingObject * obj = dynamic_cast<OsiIntegerBranchingObject *> (bobj);
      if (obj) {
	const OsiObject * object = obj->originalObject();
	// have to compute object number as not saved
	int iObject;
	int iColumn = object->columnNumber();
	for (iObject = 0 ; iObject < numberObjects_ ; iObject++) {
	  if (object_[iObject]->columnNumber()==iColumn)
	    break;
	}
	assert (iObject<numberObjects_);
	int branch = obj->firstBranch();
	if (obj->branchIndex()==2)
	  branch = 1-branch;
	assert (branch==0||branch==1);
	double originalValue=node->objectiveValue();
	double objectiveValue = solver_->getObjValue()*solver_->getObjSense();
	double changeInObjective = CoinMax(0.0,objectiveValue-originalValue);
	int iStatus = (feasible) ? 0 : 0;
	double value = obj->value();
	if (branch)
	  value = ceil(value)-value;
	else
	  value = value -floor(value);
	branchingMethod_->chooseMethod()->updateInformation(iObject,branch,changeInObjective,
							    value,iStatus);
      }
    }
  }
#endif

#ifdef CBC_DEBUG
  if (feasible)
  { printf("Obj value %g (%s) %d rows\n",solver_->getObjValue(),
	   (solver_->isProvenOptimal())?"proven":"unproven",
	   solver_->getNumRows()) ; }
  
  else
  { printf("Infeasible %d rows\n",solver_->getNumRows()) ; }
#endif
  if ((specialOptions_&1)!=0) {
/*
  If the RowCutDebugger said we were compatible with the optimal solution,
  and now we're suddenly infeasible, we might be confused. Then again, we
  may have fathomed by bound, heading for a rediscovery of an optimal solution.
*/
    if (onOptimalPath && !solver_->isDualObjectiveLimitReached()) {
      if (!feasible) {
        solver_->writeMps("infeas");
        CoinWarmStartBasis *slack =
          dynamic_cast<CoinWarmStartBasis *>(solver_->getEmptyWarmStart()) ;
        solver_->setWarmStart(slack);
        delete slack ;
        solver_->setHintParam(OsiDoReducePrint,false,OsiHintDo,0) ;
        solver_->initialSolve();
      }
      assert(feasible) ;
    }
  }

  if (!feasible) {
    numberInfeasibleNodes_++;
# ifdef COIN_HAS_CLP
    if (clpSolver) 
    clpSolver->setSpecialOptions(saveClpOptions);
# endif
    return (false) ;
  }
  sumChangeObjective1_ += solver_->getObjValue()*solver_->getObjSense()
    - objectiveValue ;
  if ( getCurrentSeconds() > dblParam_[CbcMaximumSeconds] )
    numberTries=0; // exit
  //if ((numberNodes_%100)==0)
  //printf("XXa sum obj changed by %g\n",sumChangeObjective1_);
  objectiveValue = solver_->getObjValue()*solver_->getObjSense();
  // Return at once if numberTries zero
  if (!numberTries) {
    cuts=OsiCuts();
    numberNewCuts_=0;
# ifdef COIN_HAS_CLP
    if (clpSolver) 
    clpSolver->setSpecialOptions(saveClpOptions);
# endif
    return true;
  }
/*
  Do reduced cost fixing.
*/
  reducedCostFix() ;
/*
  Set up for at most numberTries rounds of cut generation. If numberTries is
  negative, we'll ignore the minimumDrop_ cutoff and keep generating cuts for
  the specified number of rounds.
*/
  double minimumDrop = minimumDrop_ ;
  if (numberTries<0)
  { numberTries = -numberTries ;
    minimumDrop = -1.0 ; }
/*
  Is it time to scan the cuts in order to remove redundant cuts? If so, set
  up to do it.
*/
# define SCANCUTS 100  
  int *countColumnCuts = NULL ;
  // Always accumulate row cut counts
  int * countRowCuts =new int[numberCutGenerators_+numberHeuristics_] ;
  memset(countRowCuts,0,
	 (numberCutGenerators_+numberHeuristics_)*sizeof(int)) ;
  bool fullScan = false ;
  if ((numberNodes_%SCANCUTS) == 0)
  { fullScan = true ;
    countColumnCuts = new int[numberCutGenerators_+numberHeuristics_] ;
    memset(countColumnCuts,0,
	   (numberCutGenerators_+numberHeuristics_)*sizeof(int)) ; }

  double direction = solver_->getObjSense() ;
  double startObjective = solver_->getObjValue()*direction ;

  currentPassNumber_ = 0 ;
  double primalTolerance = 1.0e-7 ;
  // We may need to keep going on
  bool keepGoing=false;
/*
  Begin cut generation loop. Cuts generated during each iteration are
  collected in theseCuts. The loop can be divided into four phases:
   1) Prep: Fix variables using reduced cost. In the first iteration only,
      consider scanning globalCuts_ and activating any applicable cuts.
   2) Cut Generation: Call each generator and heuristic registered in the
      generator_ and heuristic_ arrays. Newly generated global cuts are
      copied to globalCuts_ at this time.
   3) Cut Installation and Reoptimisation: Install column and row cuts in
      the solver. Copy row cuts to cuts (parameter). Reoptimise.
   4) Cut Purging: takeOffCuts() removes inactive cuts from the solver, and
      does the necessary bookkeeping in the model.
*/
  do
  { currentPassNumber_++ ;
    numberTries-- ;
    if (numberTries<0&&keepGoing) {
      // switch off all normal ones
      for (int i = 0;i<numberCutGenerators_;i++) {
        if (!generator_[i]->mustCallAgain())
          generator_[i]->setSwitchedOff(true);
      }
    }
    keepGoing=false;
    OsiCuts theseCuts ;
/*
  Scan previously generated global column and row cuts to see if any are
  useful.
*/
    int numberViolated=0;
    if (currentPassNumber_ == 1 && howOftenGlobalScan_ > 0 &&
	(numberNodes_%howOftenGlobalScan_) == 0)
    { int numberCuts = globalCuts_.sizeColCuts() ;
      int i;
      // possibly extend whichGenerator
      resizeWhichGenerator(numberViolated, numberViolated+numberCuts);
      for ( i = 0 ; i < numberCuts ; i++)
      { OsiColCut *thisCut = globalCuts_.colCutPtr(i) ;
	if (thisCut->violated(cbcColSolution_)>primalTolerance) {
	  printf("Global cut added - violation %g\n",
		 thisCut->violated(cbcColSolution_)) ;
	  whichGenerator_[numberViolated++]=-1;
#ifndef GLOBAL_CUTS_JUST_POINTERS
	  theseCuts.insert(*thisCut) ;
#else
	  theseCuts.insert(thisCut) ;
#endif
	}
      }
      numberCuts = globalCuts_.sizeRowCuts() ;
      // possibly extend whichGenerator
      resizeWhichGenerator(numberViolated, numberViolated+numberCuts);
      for ( i = 0;i<numberCuts;i++) {
	OsiRowCut * thisCut = globalCuts_.rowCutPtr(i) ;
	if (thisCut->violated(cbcColSolution_)>primalTolerance) {
	  //printf("Global cut added - violation %g\n",
	  // thisCut->violated(cbcColSolution_)) ;
	  whichGenerator_[numberViolated++]=-1;
#ifndef GLOBAL_CUTS_JUST_POINTERS
	  theseCuts.insert(*thisCut) ;
#else
	  theseCuts.insert(thisCut) ;
#endif
	}
      }
      numberGlobalViolations_+=numberViolated;
    }
/*
  Generate new cuts (global and/or local) and/or apply heuristics.  If
  CglProbing is used, then it should be first as it can fix continuous
  variables.

  At present, CglProbing is the only case where generateCuts will return
  true. generateCuts actually modifies variable bounds in the solver when
  CglProbing indicates that it can fix a variable. Reoptimisation is required
  to take full advantage.

  The need to resolve here should only happen after a heuristic solution.
  (Note default OSI implementation of optimalBasisIsAvailable always returns
  false.)
*/
    if (solverCharacteristics_->warmStart()&&
        !solver_->optimalBasisIsAvailable()) {
      //printf("XXXXYY no opt basis\n");
      resolve(node ? node->nodeInfo() : NULL,3);
    }
    if (nextRowCut_) {
      // branch was a cut - add it
      theseCuts.insert(*nextRowCut_);
      if (handler_->logLevel()>1)
        nextRowCut_->print();
      const OsiRowCut * cut=nextRowCut_;
      double lb = cut->lb();
      double ub = cut->ub();
      int n=cut->row().getNumElements();
      const int * column = cut->row().getIndices();
      const double * element = cut->row().getElements();
      double sum=0.0;
      for (int i=0;i<n;i++) {
	int iColumn = column[i];
	double value = element[i];
	//if (cbcColSolution_[iColumn]>1.0e-7)
	//printf("value of %d is %g\n",iColumn,cbcColSolution_[iColumn]);
	sum += value * cbcColSolution_[iColumn];
      }
      delete nextRowCut_;
      nextRowCut_=NULL;
      if (handler_->logLevel()>1)
	printf("applying branch cut, sum is %g, bounds %g %g\n",sum,lb,ub);
      // possibly extend whichGenerator
      resizeWhichGenerator(numberViolated, numberViolated+1);
      // set whichgenerator (also serves as marker to say don't delete0
      whichGenerator_[numberViolated++]=-2;
    }

    // reset probing info
    if (probingInfo_)
      probingInfo_->initializeFixing();
    int i;
    if ((threadMode_&2)==0||numberNodes_) {
      for (i = 0;i<numberCutGenerators_;i++) {
	int numberRowCutsBefore = theseCuts.sizeRowCuts() ;
	int numberColumnCutsBefore = theseCuts.sizeColCuts() ;
	bool generate = generator_[i]->normal();
	// skip if not optimal and should be (maybe a cut generator has fixed variables)
	if (generator_[i]->needsOptimalBasis()&&!solver_->basisIsAvailable())
	  generate=false;
	if (generator_[i]->switchedOff())
	  generate=false;;
	if (generate) {
	  bool mustResolve = 
	    generator_[i]->generateCuts(theseCuts,fullScan,solver_,node) ;
	  int numberRowCutsAfter = theseCuts.sizeRowCuts() ;
	  if(numberRowCutsBefore < numberRowCutsAfter &&
	     generator_[i]->mustCallAgain())
	    keepGoing=true; // say must go round
	  // Check last cut to see if infeasible
	  if(numberRowCutsBefore < numberRowCutsAfter) {
	    const OsiRowCut * thisCut = theseCuts.rowCutPtr(numberRowCutsAfter-1) ;
	    if (thisCut->lb()>thisCut->ub()) {
	      feasible = false; // sub-problem is infeasible
	      break;
	    }
	  }
#ifdef CBC_DEBUG
	  {
	    int k ;
	    for (k = numberRowCutsBefore;k<numberRowCutsAfter;k++) {
	      OsiRowCut thisCut = theseCuts.rowCut(k) ;
	      /* check size of elements.
		 We can allow smaller but this helps debug generators as it
		 is unsafe to have small elements */
	      int n=thisCut.row().getNumElements();
	      const int * column = thisCut.row().getIndices();
	      const double * element = thisCut.row().getElements();
	      //assert (n);
	      for (int i=0;i<n;i++) {
		double value = element[i];
		assert(fabs(value)>1.0e-12&&fabs(value)<1.0e20);
	      }
	    }
	  }
#endif
	  if (mustResolve) {
	    int returncode = resolve(node ? node->nodeInfo() : NULL,2);
	    feasible = returnCode  != 0 ;
	    if (returncode<0)
	      numberTries=0;
	    if ((specialOptions_&1)!=0) {
	      debugger = solver_->getRowCutDebugger() ;
	      if (debugger) 
		onOptimalPath = (debugger->onOptimalPath(*solver_)) ;
	      else
		onOptimalPath=false;
	      if (onOptimalPath && !solver_->isDualObjectiveLimitReached())
		assert(feasible) ;
	    }
	    if (!feasible)
	      break ;
	  }
	}
	int numberRowCutsAfter = theseCuts.sizeRowCuts() ;
	int numberColumnCutsAfter = theseCuts.sizeColCuts() ;
	
	if ((specialOptions_&1)!=0) {
	  if (onOptimalPath) {
	    int k ;
	    for (k = numberRowCutsBefore;k<numberRowCutsAfter;k++) {
	      OsiRowCut thisCut = theseCuts.rowCut(k) ;
	      if(debugger->invalidCut(thisCut)) {
		solver_->writeMps("badCut");
#ifdef NDEBUG
		printf("Cut generator %d (%s) produced invalid cut (%dth in this go)\n",
		       i,generator_[i]->cutGeneratorName(),k-numberRowCutsBefore);
		abort();
#endif
	      }
	      assert(!debugger->invalidCut(thisCut)) ;
	    }
	  }
	}
/*
  The cut generator has done its thing, and maybe it generated some
  cuts.  Do a bit of bookkeeping: load
  whichGenerator[i] with the index of the generator responsible for a cut,
  and place cuts flagged as global in the global cut pool for the model.

  lastNumberCuts is the sum of cuts added in previous iterations; it's the
  offset to the proper starting position in whichGenerator.
*/
	int numberBefore =
	  numberRowCutsBefore+numberColumnCutsBefore+lastNumberCuts ;
	int numberAfter =
	  numberRowCutsAfter+numberColumnCutsAfter+lastNumberCuts ;
	// possibly extend whichGenerator
	resizeWhichGenerator(numberBefore, numberAfter);
	int j ;
	if (fullScan) {
	  // counts
	  countColumnCuts[i] += numberColumnCutsAfter-numberColumnCutsBefore ;
	}
	countRowCuts[i] += numberRowCutsAfter-numberRowCutsBefore ;
	
	bool dodgyCuts=false;
	for (j = numberRowCutsBefore;j<numberRowCutsAfter;j++) {
	  const OsiRowCut * thisCut = theseCuts.rowCutPtr(j) ;
	  if (thisCut->lb()>1.0e10||thisCut->ub()<-1.0e10) {
	    dodgyCuts=true;
	    break;
	  }
	  whichGenerator_[numberBefore++] = i ;
	  if (thisCut->lb()>thisCut->ub())
	    violated=-2; // sub-problem is infeasible
	  if (thisCut->globallyValid()) {
	    // add to global list
	    OsiRowCut newCut(*thisCut);
	    newCut.setGloballyValid(2);
	    newCut.mutableRow().setTestForDuplicateIndex(false);
	    globalCuts_.insert(newCut) ;
	  }
	}
	if (dodgyCuts) {
	  for (int k=numberRowCutsAfter-1;k>=j;k--) {
	    const OsiRowCut * thisCut = theseCuts.rowCutPtr(k) ;
	    if (thisCut->lb()>thisCut->ub())
	      violated=-2; // sub-problem is infeasible
	    if (thisCut->lb()>1.0e10||thisCut->ub()<-1.0e10) 
	      theseCuts.eraseRowCut(k);
	  }
	  numberRowCutsAfter = theseCuts.sizeRowCuts() ;
	  for (;j<numberRowCutsAfter;j++) {
	    const OsiRowCut * thisCut = theseCuts.rowCutPtr(j) ;
	    whichGenerator_[numberBefore++] = i ;
	    if (thisCut->globallyValid()) {
	      // add to global list
	      OsiRowCut newCut(*thisCut);
	      newCut.setGloballyValid(2);
	      newCut.mutableRow().setTestForDuplicateIndex(false);
	      globalCuts_.insert(newCut) ;
	    }
	  }
	}
	for (j = numberColumnCutsBefore;j<numberColumnCutsAfter;j++) {
	  whichGenerator_[numberBefore++] = i ;
	  const OsiColCut * thisCut = theseCuts.colCutPtr(j) ;
	  if (thisCut->globallyValid()) {
	    // add to global list
	    OsiColCut newCut(*thisCut);
	    newCut.setGloballyValid(2);
	    globalCuts_.insert(newCut) ;
	  }
	}
      }
      // Add in any violated saved cuts
      if (!theseCuts.sizeRowCuts()&&!theseCuts.sizeColCuts()) {
	int numberOld = theseCuts.sizeRowCuts();
	int numberCuts = slackCuts.sizeRowCuts() ;
	int i;
	// possibly extend whichGenerator
	resizeWhichGenerator(numberOld, numberOld+numberCuts);
	for ( i = 0;i<numberCuts;i++) {
	  const OsiRowCut * thisCut = slackCuts.rowCutPtr(i) ;
	  if (thisCut->violated(cbcColSolution_)>100.0*primalTolerance) {
	    if (messageHandler()->logLevel()>2)
	      printf("Old cut added - violation %g\n",
		     thisCut->violated(cbcColSolution_)) ;
	    whichGenerator_[numberOld++]=-1;
	    theseCuts.insert(*thisCut) ;
	  }
	}
      }
    } else {
      // do cuts independently
      OsiCuts * eachCuts = new OsiCuts [numberCutGenerators_];;
#ifdef CBC_THREAD
      if (!threadModel) {
#endif
	// generate cuts
	for (i = 0;i<numberCutGenerators_;i++) {
	  bool generate = generator_[i]->normal();
	  // skip if not optimal and should be (maybe a cut generator has fixed variables)
	  if (generator_[i]->needsOptimalBasis()&&!solver_->basisIsAvailable())
	    generate=false;
	  if (generator_[i]->switchedOff())
	    generate=false;;
	  if (generate) 
	    generator_[i]->generateCuts(eachCuts[i],fullScan,solver_,node) ;
	}
#ifdef CBC_THREAD
      } else {
	for (i=0;i<numberThreads_;i++) {
	  // set solver here after cloning
	  threadModel[i]->solver_=solver_->clone();
	  threadModel[i]->numberNodes_ = (fullScan) ? 1 : 0;
	}
	// generate cuts
	for (i = 0;i<numberCutGenerators_;i++) {
	  bool generate = generator_[i]->normal();
	  // skip if not optimal and should be (maybe a cut generator has fixed variables)
	  if (generator_[i]->needsOptimalBasis()&&!solver_->basisIsAvailable())
	    generate=false;
	  if (generator_[i]->switchedOff())
	    generate=false;;
	  if (generate) { 
	    bool finished=false;
	    int iThread=-1;
	    // see if any available
	    for (iThread=0;iThread<numberThreads_;iThread++) {
	      if (threadInfo[iThread].returnCode) {
		finished=true;
		break;
	      } else if (threadInfo[iThread].returnCode==0) {
		pthread_cond_signal(threadInfo[iThread].condition2); // unlock
	      }
	    }
	    while (!finished) {
	      pthread_mutex_lock(&condition_mutex);
	      struct timespec absTime;
	      clock_gettime(CLOCK_REALTIME,&absTime);
	      absTime.tv_nsec += 1000000; // millisecond
	      if (absTime.tv_nsec>=1000000000) {
		absTime.tv_nsec -= 1000000000;
		absTime.tv_sec++;
	      }
	      pthread_cond_timedwait(&condition_main,&condition_mutex,&absTime);
	      pthread_mutex_unlock(&condition_mutex);
	      for (iThread=0;iThread<numberThreads_;iThread++) {
		if (threadInfo[iThread].returnCode>0) {
		  finished=true;
		  break;
		} else if (threadInfo[iThread].returnCode==0) {
		  pthread_cond_signal(threadInfo[iThread].condition2); // unlock
		}
	      }
	    }
	    assert (iThread<numberThreads_);
	    assert (threadInfo[iThread].returnCode);
	    threadModel[iThread]->generator_[0]=generator_[i];
	    threadModel[iThread]->object_ = (OsiObject **) (eachCuts+i);
	    // allow to start
	    threadInfo[iThread].returnCode=0;
	    pthread_cond_signal(threadInfo[iThread].condition2); // unlock
	  }
	}
	// wait
	for (int iThread=0;iThread<numberThreads_;iThread++) {
	  if (threadInfo[iThread].returnCode==0) {
	    bool finished=false;
	    pthread_cond_signal(threadInfo[iThread].condition2); // unlock
	    while (!finished) {
	      pthread_mutex_lock(&condition_mutex);
	      struct timespec absTime;
	      clock_gettime(CLOCK_REALTIME,&absTime);
	      absTime.tv_nsec += 1000000; // millisecond
	      if (absTime.tv_nsec>=1000000000) {
		absTime.tv_nsec -= 1000000000;
		absTime.tv_sec++;
	      }
	      pthread_cond_timedwait(&condition_main,&condition_mutex,&absTime);
	      pthread_mutex_unlock(&condition_mutex);
	      if (threadInfo[iThread].returnCode>0) {
		finished=true;
		break;
	      } else if (threadInfo[iThread].returnCode==0) {
		pthread_cond_signal(threadInfo[iThread].condition2); // unlock
	      }
	    }
	  }
	  assert (threadInfo[iThread].returnCode);
	  // say available
	  threadInfo[iThread].returnCode=-1;
	  delete threadModel[iThread]->solver_;
	  threadModel[iThread]->solver_=NULL;
	}
      }
#endif
      // Now put together
      for (i = 0;i<numberCutGenerators_;i++) {
	// add column cuts
	int numberColumnCutsBefore = theseCuts.sizeColCuts() ;
	int numberColumnCuts = eachCuts[i].sizeColCuts();
	int numberColumnCutsAfter = numberColumnCutsBefore
	  + numberColumnCuts;
	int j;
	for (j=0;j<numberColumnCuts;j++) {
	  theseCuts.insert(eachCuts[i].colCut(j));
	}
	int numberRowCutsBefore = theseCuts.sizeRowCuts() ;
	int numberRowCuts = eachCuts[i].sizeRowCuts();
	int numberRowCutsAfter = numberRowCutsBefore
	  + numberRowCuts;
	if (numberRowCuts) {
	  for (j=0;j<numberRowCuts;j++) {
	    const OsiRowCut * thisCut = eachCuts[i].rowCutPtr(j) ;
	    if (thisCut->lb()<=1.0e10&&thisCut->ub()>=-1.0e10) 
	      theseCuts.insert(eachCuts[i].rowCut(j));
	  }
	  if (generator_[i]->mustCallAgain())
	    keepGoing=true; // say must go round
	  // Check last cut to see if infeasible
	  const OsiRowCut * thisCut = theseCuts.rowCutPtr(numberRowCutsAfter-1) ;
	  if (thisCut->lb()>thisCut->ub()) {
	    feasible = false; // sub-problem is infeasible
	    break;
	  }
	}
#ifdef CBC_DEBUG
	{
	  int k ;
	  for (k = numberRowCutsBefore;k<numberRowCutsAfter;k++) {
	    OsiRowCut thisCut = theseCuts.rowCut(k) ;
	    /* check size of elements.
	       We can allow smaller but this helps debug generators as it
	       is unsafe to have small elements */
	    int n=thisCut.row().getNumElements();
	    const int * column = thisCut.row().getIndices();
	    const double * element = thisCut.row().getElements();
	    //assert (n);
	    for (int i=0;i<n;i++) {
	      double value = element[i];
	      assert(fabs(value)>1.0e-12&&fabs(value)<1.0e20);
	    }
	  }
	}
#endif
	if ((specialOptions_&1)!=0) {
	  if (onOptimalPath) {
	    int k ;
	    for (k = numberRowCutsBefore;k<numberRowCutsAfter;k++) {
	      OsiRowCut thisCut = theseCuts.rowCut(k) ;
	      if(debugger->invalidCut(thisCut)) {
		solver_->writeMps("badCut");
#ifdef NDEBUG
		printf("Cut generator %d (%s) produced invalid cut (%dth in this go)\n",
		       i,generator_[i]->cutGeneratorName(),k-numberRowCutsBefore);
		abort();
#endif
	      }
	      assert(!debugger->invalidCut(thisCut)) ;
	    }
	  }
	}
/*
  The cut generator has done its thing, and maybe it generated some
  cuts.  Do a bit of bookkeeping: load
  whichGenerator[i] with the index of the generator responsible for a cut,
  and place cuts flagged as global in the global cut pool for the model.

  lastNumberCuts is the sum of cuts added in previous iterations; it's the
  offset to the proper starting position in whichGenerator.
*/
	int numberBefore =
	  numberRowCutsBefore+numberColumnCutsBefore+lastNumberCuts ;
	int numberAfter =
	  numberRowCutsAfter+numberColumnCutsAfter+lastNumberCuts ;
	// possibly extend whichGenerator
	resizeWhichGenerator(numberBefore, numberAfter);
	if (fullScan) {
	  // counts
	  countColumnCuts[i] += numberColumnCutsAfter-numberColumnCutsBefore ;
	}
	countRowCuts[i] += numberRowCutsAfter-numberRowCutsBefore ;
	
	for (j = numberRowCutsBefore;j<numberRowCutsAfter;j++) {
	  whichGenerator_[numberBefore++] = i ;
	  const OsiRowCut * thisCut = theseCuts.rowCutPtr(j) ;
	  if (thisCut->lb()>thisCut->ub())
	    violated=-2; // sub-problem is infeasible
	  if (thisCut->globallyValid()) {
	    // add to global list
	    OsiRowCut newCut(*thisCut);
	    newCut.setGloballyValid(2);
	    newCut.mutableRow().setTestForDuplicateIndex(false);
	    globalCuts_.insert(newCut) ;
	  }
	}
	for (j = numberColumnCutsBefore;j<numberColumnCutsAfter;j++) {
	  whichGenerator_[numberBefore++] = i ;
	  const OsiColCut * thisCut = theseCuts.colCutPtr(j) ;
	  if (thisCut->globallyValid()) {
	    // add to global list
	    OsiColCut newCut(*thisCut);
	    newCut.setGloballyValid(2);
	    globalCuts_.insert(newCut) ;
	  }
	}
      }
      // Add in any violated saved cuts
      if (!theseCuts.sizeRowCuts()&&!theseCuts.sizeColCuts()) {
	int numberOld = theseCuts.sizeRowCuts();
	int numberCuts = slackCuts.sizeRowCuts() ;
	int i;
	// possibly extend whichGenerator
	resizeWhichGenerator(numberOld, numberOld+numberCuts);
	for ( i = 0;i<numberCuts;i++) {
	  const OsiRowCut * thisCut = slackCuts.rowCutPtr(i) ;
	  if (thisCut->violated(cbcColSolution_)>100.0*primalTolerance) {
	    if (messageHandler()->logLevel()>2)
	      printf("Old cut added - violation %g\n",
		     thisCut->violated(cbcColSolution_)) ;
	    whichGenerator_[numberOld++]=-1;
	    theseCuts.insert(*thisCut) ;
	  }
	}
      }
      delete [] eachCuts;
    }
    //if (!feasible)
    //break;
/*
  End of the loop to exercise each generator - try heuristics
  - unless at root node and first pass
*/
    if (numberNodes_||currentPassNumber_!=1) {
      double * newSolution = new double [numberColumns] ;
      double heuristicValue = getCutoff() ;
      int found = -1; // no solution found
      for (i = 0;i<numberHeuristics_;i++) {
	// see if heuristic will do anything
	double saveValue = heuristicValue ;
	int ifSol = 
	  heuristic_[i]->solution(heuristicValue,
				  newSolution,
				  theseCuts) ;
	if (ifSol>0) {
	  // better solution found
	  found = i ;
          incrementUsed(newSolution);
	} else if (ifSol<0) {
	  heuristicValue = saveValue ;
	}
      }
/*
  Did any of the heuristics turn up a new solution? Record it before we free
  the vector.
*/
      if (found >= 0) { 
	phase_=4;
	incrementUsed(newSolution);
	lastHeuristic_ = heuristic_[found];
	setBestSolution(CBC_ROUNDING,heuristicValue,newSolution) ;
	CbcTreeLocal * tree 
          = dynamic_cast<CbcTreeLocal *> (tree_);
	if (tree)
	  tree->passInSolution(bestSolution_,heuristicValue);
      }
      delete [] newSolution ;
    }

#if 0
    // switch on to get all cuts printed
    theseCuts.printCuts() ;
#endif
    int numberColumnCuts = theseCuts.sizeColCuts() ;
    int numberRowCuts = theseCuts.sizeRowCuts() ;
    if (violated>=0)
      violated = numberRowCuts + numberColumnCuts ;
/*
  Apply column cuts (aka bound tightening). This may be partially redundant
  for column cuts returned by CglProbing, as generateCuts installs bounds
  from CglProbing when it determines it can fix a variable.

  TODO: Looks like the use of violated has evolved. The value set above is
	completely ignored. All that's left is violated == -1 indicates some
	cut is violated, violated == -2 indicates infeasibility. Only
	infeasibility warrants exceptional action.

  TODO: Strikes me that this code will fail to detect infeasibility, because
	the breaks escape the inner loops but the outer loop keeps going.
	Infeasibility in an early cut will be overwritten if a later cut is
	merely violated.
*/
    if (numberColumnCuts) {

#ifdef CBC_DEBUG
      double * oldLower = new double [numberColumns] ;
      double * oldUpper = new double [numberColumns] ;
      memcpy(oldLower,cbcColLower_,numberColumns*sizeof(double)) ;
      memcpy(oldUpper,cbcColUpper_,numberColumns*sizeof(double)) ;
#endif

      double integerTolerance = getDblParam(CbcIntegerTolerance) ;
      for (int i = 0;i<numberColumnCuts;i++) {
	const OsiColCut * thisCut = theseCuts.colCutPtr(i) ;
	const CoinPackedVector & lbs = thisCut->lbs() ;
	const CoinPackedVector & ubs = thisCut->ubs() ;
	int j ;
	int n ;
	const int * which ;
	const double * values ;
	n = lbs.getNumElements() ;
	which = lbs.getIndices() ;
	values = lbs.getElements() ;
	for (j = 0;j<n;j++) {
	  int iColumn = which[j] ;
	  double value = cbcColSolution_[iColumn] ;
#if CBC_DEBUG>1
	  printf("%d %g %g %g %g\n",iColumn,oldLower[iColumn],
		 cbcColSolution_[iColumn],oldUpper[iColumn],values[j]) ;
#endif
	  solver_->setColLower(iColumn,values[j]) ;
	  if (value<values[j]-integerTolerance)
	    violated = -1 ;
	  if (values[j]>cbcColUpper_[iColumn]+integerTolerance) {
	    // infeasible
	    violated = -2 ;
	    break ;
	  }
	}
	n = ubs.getNumElements() ;
	which = ubs.getIndices() ;
	values = ubs.getElements() ;
	for (j = 0;j<n;j++) {
	  int iColumn = which[j] ;
	  double value = cbcColSolution_[iColumn] ;
#if CBC_DEBUG>1
	  printf("%d %g %g %g %g\n",iColumn,oldLower[iColumn],
		 cbcColSolution_[iColumn],oldUpper[iColumn],values[j]) ;
#endif
	  solver_->setColUpper(iColumn,values[j]) ;
	  if (value>values[j]+integerTolerance)
	    violated = -1 ;
	  if (values[j]<cbcColLower_[iColumn]-integerTolerance) {
	    // infeasible
	    violated = -2 ;
	    break ;
	  }
	}
      }
#ifdef CBC_DEBUG
      delete [] oldLower ;
      delete [] oldUpper ;
#endif
    }
/*
  End installation of column cuts. The break here escapes the numberTries
  loop.
*/
    if (violated == -2||!feasible) {
      // infeasible
      feasible = false ;
      violated = -2;
      if (!numberNodes_) 
	messageHandler()->message(CBC_INFEAS,
				  messages())
				    << CoinMessageEol ;
      break ;
    }
/*
  Now apply the row (constraint) cuts. This is a bit more work because we need
  to obtain and augment the current basis.

  TODO: Why do this work, if there are no row cuts? The current basis will do
	just fine.
*/
    int numberRowsNow = solver_->getNumRows() ;
    assert(numberRowsNow == numberRowsAtStart+lastNumberCuts) ;
    int numberToAdd = theseCuts.sizeRowCuts() ;
    numberNewCuts_ = lastNumberCuts+numberToAdd ;
/*
  Now actually add the row cuts and reoptimise.

  Install the cuts in the solver using applyRowCuts and
  augment the basis with the corresponding slack. We also add each row cut to
  the set of row cuts (cuts.insert()) supplied as a parameter. The new basis
  must be set with setWarmStart().

  TODO: Seems to me the original code could allocate addCuts with size 0, if
	numberRowCuts was 0 and numberColumnCuts was nonzero. That might
	explain the memory fault noted in the comment by AJK.  Unfortunately,
	just commenting out the delete[] results in massive memory leaks. Try
	a revision to separate the row cut case. Why do we need addCuts at
	all? A typing issue, apparently: OsiCut vs. OsiRowCut.

  TODO: It looks to me as if numberToAdd and numberRowCuts are identical at
	this point. Confirm & get rid of one of them.

  TODO: Any reason why the three loops can't be consolidated?
*/
    if (numberRowCuts > 0 || numberColumnCuts > 0)
    { if (numberToAdd > 0)
      { int i ;
	// Faster to add all at once
	const OsiRowCut ** addCuts = new const OsiRowCut * [numberToAdd] ;
	for (i = 0 ; i < numberToAdd ; i++)
	{ addCuts[i] = &theseCuts.rowCut(i) ; }
	solver_->applyRowCuts(numberToAdd,addCuts) ;
	// AJK this caused a memory fault on Win32
	// may be okay now with ** form
	delete [] addCuts ;
	for (i = 0 ; i < numberToAdd ; i++)
	{ cuts.insert(theseCuts.rowCut(i)) ; }
        CoinWarmStartBasis * basis = dynamic_cast<CoinWarmStartBasis*>(solver_->getWarmStart()) ;
        assert(basis != NULL); // make sure not volume
/* dylp bug

  Consistent size used by OsiDylp as sanity check. Implicit resize seen
  as an error. Hence this call to resize is necessary.
*/
        basis->resize(numberRowsAtStart+numberNewCuts_,numberColumns) ;
	for (i = 0 ; i < numberToAdd ; i++) 
	{ basis->setArtifStatus(numberRowsNow+i,
				 CoinWarmStartBasis::basic) ; }
	if (solver_->setWarmStart(basis) == false)
	{ throw CoinError("Fail setWarmStart() after cut installation.",
			  "solveWithCuts","CbcModel") ; }
        delete basis;
      }
      feasible = resolve(node ? node->nodeInfo() : NULL,2) ;
      if ( getCurrentSeconds() > dblParam_[CbcMaximumSeconds] )
        numberTries=0; // exit
#     ifdef CBC_DEBUG
      printf("Obj value after cuts %g %d rows\n",solver_->getObjValue(),
	      solver_->getNumRows()) ;
      if (onOptimalPath && !solver_->isDualObjectiveLimitReached())
	assert(feasible) ;
#     endif
    }
/*
  No cuts. Cut short the cut generation (numberTries) loop.
*/
    else
    { numberTries = 0 ; }
/*
  If the problem is still feasible, first, call takeOffCuts() to remove cuts
  that are now slack. takeOffCuts() will call the solver to reoptimise if
  that's needed to restore a valid solution.

  Next, see if we should quit due to diminishing returns:
    * we've tried three rounds of cut generation and we're getting
      insufficient improvement in the objective; or
    * we generated no cuts; or
    * the solver declared optimality with 0 iterations after we added the
      cuts generated in this round.
  If we decide to keep going, prep for the next iteration.

  It just seems more safe to tell takeOffCuts() to call resolve(), even if
  we're not continuing cut generation. Otherwise code executed between here
  and final disposition of the node will need to be careful not to access the
  lp solution. It can happen that we lose feasibility in takeOffCuts ---
  numerical jitters when the cutoff bound is epsilon less than the current
  best, and we're evaluating an alternative optimum.

  TODO: After successive rounds of code motion, there seems no need to
	distinguish between the two checks for aborting the cut generation
	loop. Confirm and clean up.
*/
    if (feasible)
    { int cutIterations = solver_->getIterationCount() ;
      if (numberOldActiveCuts_+numberNewCuts_) {
        OsiCuts * saveCuts = node ? NULL : &slackCuts;
        takeOffCuts(cuts,resolveAfterTakeOffCuts_,saveCuts) ;
        if (solver_->isDualObjectiveLimitReached()&&resolveAfterTakeOffCuts_)
          { feasible = false ;
#	ifdef CBC_DEBUG
          double z = solver_->getObjValue() ;
          double cut = getCutoff() ;
          printf("Lost feasibility by %g in takeOffCuts; z = %g, cutoff = %g\n",
                 z-cut,z,cut) ;
#	endif
          }
      }
      if (feasible) {
	numberRowsAtStart = numberOldActiveCuts_+numberRowsAtContinuous_ ;
        lastNumberCuts = numberNewCuts_ ;
        if (direction*solver_->getObjValue() < lastObjective+minimumDrop &&
            currentPassNumber_ >= 3 && !keepGoing)
          { numberTries = 0 ; }
        if (numberRowCuts+numberColumnCuts == 0 || cutIterations == 0)
          { break ; }
        if (numberTries > 0) {
	  reducedCostFix() ;
	  lastObjective = direction*solver_->getObjValue() ;
	}
      }
    }
/*
  We've lost feasibility --- this node won't be referencing the cuts we've
  been collecting, so decrement the reference counts.
*/
    if (!feasible)
    { int i ;
      if (currentNumberCuts_) {
	lockThread();
	for (i = 0;i<currentNumberCuts_;i++) {
	  // take off node
	  if (addedCuts_[i]) {
	    if (!addedCuts_[i]->decrement())
	      delete addedCuts_[i] ;
	    addedCuts_[i] = NULL ;
	  }
	}
      }
      numberTries = 0 ;
    }
  } while (numberTries>0||keepGoing) ;
  {
    // switch on
    for (int i = 0;i<numberCutGenerators_;i++) 
      generator_[i]->setSwitchedOff(false);
  }
  //check feasibility.
  //If solution seems to be integer feasible calling setBestSolution
  //will eventually add extra global cuts which we need to install at
  //the nodes


  if(feasible&&solverCharacteristics_->solutionAddsCuts())  //check integer feasibility
  {
    bool integerFeasible = true;
    const double * save = testSolution_;
    testSolution_ = solver_->getColSolution();
    // point to useful information
    OsiBranchingInformation usefulInfo=usefulInformation();
    for (int i=0;i<numberObjects_ && integerFeasible;i++)
    {
      int preferredWay;
      double infeasibility = object_[i]->infeasibility(&usefulInfo,preferredWay);
      if(infeasibility)
        integerFeasible = false;
    }
    testSolution_ = save;
    if(integerFeasible)//update
    {
      double objValue = solver_->getObjValue();
      int numberGlobalBefore = globalCuts_.sizeRowCuts();
      // SOLUTION2 so won't up cutoff or print message
      setBestSolution(CBC_SOLUTION2, objValue, 
                      solver_->getColSolution(),0);
      int numberGlobalAfter = globalCuts_.sizeRowCuts();
      int numberToAdd = numberGlobalAfter - numberGlobalBefore;
      if(numberToAdd > 0)
        //We have added some cuts say they are tight at that node
        //Basis and lp should already have been updated
      {
        feasible = (solver_->isProvenOptimal() &&
                    !solver_->isDualObjectiveLimitReached()) ;
        if(feasible)
        {
          int numberCuts = numberNewCuts_ = cuts.sizeRowCuts();
      // possibly extend whichGenerator
	  resizeWhichGenerator(numberCuts, numberToAdd+numberCuts);
	  
          for (int i = numberGlobalBefore ; i < numberGlobalAfter ; i++)
	    { 	     
	      whichGenerator_[numberNewCuts_++]=-1;
#ifndef GLOBAL_CUTS_JUST_POINTERS
	      cuts.insert(globalCuts_.rowCut(i)) ; 
#else
	      OsiRowCut * rowCutPointer = globalCuts_.rowCutPtr(i);
	      cuts.insert(rowCutPointer) ; 
#endif
	    }
	  numberNewCuts_ = lastNumberCuts+numberToAdd;
          //now take off the cuts which are not tight anymore
          takeOffCuts(cuts,resolveAfterTakeOffCuts_, NULL) ;
          if (solver_->isDualObjectiveLimitReached()&&resolveAfterTakeOffCuts_)
	    { feasible = false ;}
        }
        if(!feasible) //node will be fathomed
	  { 
          for (int i = 0;i<currentNumberCuts_;i++) 
	    {
            // take off node
            if (addedCuts_[i])
	      {
		if (!addedCuts_[i]->decrement())
		  delete addedCuts_[i] ;
		addedCuts_[i] = NULL ;
	      }
	    }
	}
      }
    }
  }
/*
  Reduced cost fix at end. Must also check feasible, in case we've popped out
  because a generator indicated we're infeasible.
*/
  if (feasible && solver_->isProvenOptimal())
    reducedCostFix();
  // If at root node do heuristics
  if (!numberNodes_) {
    // First see if any cuts are slack
    int numberRows = solver_->getNumRows();
    int numberAdded = numberRows-numberRowsAtContinuous_;
    if (numberAdded) {
      CoinWarmStartBasis * basis = dynamic_cast<CoinWarmStartBasis*>(solver_->getWarmStart()) ;
      assert(basis != NULL);
      int * added = new int[numberAdded];
      int nDelete = 0;
      for (int j=numberRowsAtContinuous_;j<numberRows;j++) {
        if (basis->getArtifStatus(j)==CoinWarmStartBasis::basic) {
          //printf("%d slack!\n",j);
          added[nDelete++]=j;
        }
      }
      if (nDelete) {
        solver_->deleteRows(nDelete,added);
      }
      delete [] added;
      delete basis ;
    }
    // mark so heuristics can tell
    int savePass=currentPassNumber_;
    currentPassNumber_=999999;
    double * newSolution = new double [numberColumns] ;
    double heuristicValue = getCutoff() ;
    int found = -1; // no solution found
    for (int i = 0;i<numberHeuristics_;i++) {
      // see if heuristic will do anything
      double saveValue = heuristicValue ;
      int ifSol = heuristic_[i]->solution(heuristicValue,
                                          newSolution);
      if (ifSol>0) {
        // better solution found
        found = i ;
        incrementUsed(newSolution);
      } else {
        heuristicValue = saveValue ;
      }
    }
    currentPassNumber_=savePass;
    if (found >= 0) { 
      phase_=4;
      incrementUsed(newSolution);
      lastHeuristic_ = heuristic_[found];
      setBestSolution(CBC_ROUNDING,heuristicValue,newSolution) ;
    }
    delete [] newSolution ;
  }
  // Up change due to cuts
  if (feasible)
    sumChangeObjective2_ += solver_->getObjValue()*solver_->getObjSense()
      - objectiveValue ;
  //if ((numberNodes_%100)==0)
  //printf("XXb sum obj changed by %g\n",sumChangeObjective2_);
/*
  End of cut generation loop.

  Now, consider if we want to disable or adjust the frequency of use for any
  of the cut generators. If the client specified a positive number for
  howOften, it will never change. If the original value was negative, it'll
  be converted to 1000000+|howOften|, and this value will be adjusted each
  time fullScan is true. Actual cut generation is performed every
  howOften%1000000 nodes; the 1000000 offset is just a convenient way to
  specify that the frequency is adjustable.

  During cut generation, we recorded the number of cuts produced by each
  generator for this node. For all cuts, whichGenerator records the generator
  that produced a cut.

  TODO: All this should probably be hidden in a method of the CbcCutGenerator
  class.

  TODO: Can the loop that scans over whichGenerator to accumulate per generator
	counts be replaced by values in countRowCuts and countColumnCuts?
*/
#ifdef NODE_LOG
  int fatherNum = (node == NULL) ? -1 : node->nodeNumber();
  double value =  (node == NULL) ? -1 : node->branchingObject()->value();
  string bigOne = (solver_->getIterationCount() > 30) ? "*******" : "";
  string way = (node == NULL) ? "" : (node->branchingObject()->way())==1 ? "Down" : "Up";
  std::cout<<"Node "<<numberNodes_<<", father "<<fatherNum<<", #iterations "<<solver_->getIterationCount()<<", sol value : "<<solver_->getObjValue()<<std::endl;
#endif
  if (fullScan&&numberCutGenerators_) {
    /* If cuts just at root node then it will probably be faster to
       update matrix and leave all in */
    int willBeCutsInTree=0;
    double thisObjective = solver_->getObjValue()*direction ;
    // get sizes
    int numberRowsAdded = solver_->getNumRows()-numberRowsAtStart;
    CoinBigIndex numberElementsAdded =  solver_->getNumElements()-numberElementsAtStart ;
    double densityOld = ((double) numberElementsAtStart)/((double) numberRowsAtStart);
    double densityNew = numberRowsAdded ? ((double) (numberElementsAdded))/((double) numberRowsAdded)
      : 0.0;
    if (!numberNodes_) {
      if (numberRowsAdded)
        handler_->message(CBC_CUTS_STATS,messages_)
          <<numberRowsAdded
          <<densityNew
          <<CoinMessageEol ;
      if (thisObjective-startObjective<1.0e-5&&numberElementsAdded>0.2*numberElementsAtStart)
        willBeCutsInTree=-1;
    }
    if ((numberRowsAdded>100+0.5*numberRowsAtStart
         ||numberElementsAdded>0.5*numberElementsAtStart)
        &&(densityNew>200.0&&numberRowsAdded>100&&densityNew>2.0*densityOld)) {
      // much bigger
      //if (thisObjective-startObjective<0.1*fabs(startObjective)+1.0e-5)
      //willBeCutsInTree=-1;
      //printf("Cuts will be taken off , %d rows added with density %g\n",
      //     numberRowsAdded,densityNew);
    }
    if (densityNew>100.0&&numberRowsAdded>2&&densityNew>2.0*densityOld) {
      //if (thisObjective-startObjective<0.1*fabs(startObjective)+1.0e-5)
      //willBeCutsInTree=-2;
      //printf("Density says no cuts ? , %d rows added with density %g\n",
      //     numberRowsAdded,densityNew);
    }
    // Root node or every so often - see what to turn off
    int i ;
    for (i = 0;i<numberCutGenerators_;i++) {
      int howOften = generator_[i]->howOften() ;
      if (howOften>-90) 
        willBeCutsInTree=0;
    }
    if (!numberNodes_)
      handler_->message(CBC_ROOT,messages_)
	<<numberNewCuts_
	<<startObjective<<thisObjective
	<<currentPassNumber_
	<<CoinMessageEol ;
    int * count = new int[numberCutGenerators_] ;
    memset(count,0,numberCutGenerators_*sizeof(int)) ;
    int numberActiveGenerators=0;
    for (i = 0;i<numberNewCuts_;i++) {
      int iGenerator = whichGenerator_[i];
      if (iGenerator>=0&&iGenerator<numberCutGenerators_)
	count[iGenerator]++ ;
    }
    double totalCuts = 0.0 ;
    //#define JUST_ACTIVE
    for (i = 0;i<numberCutGenerators_;i++) {
      if (countRowCuts[i]||countColumnCuts[i])
        numberActiveGenerators++;
#ifdef JUST_ACTIVE
      totalCuts += count[i] + 5.0*countColumnCuts[i] ;
#else
      totalCuts += countRowCuts[i] + 5.0*countColumnCuts[i] ;
#endif
    }
    int iProbing=-1;
    double small = (0.5* totalCuts) /
      ((double) numberActiveGenerators) ;
    for (i = 0;i<numberCutGenerators_;i++) {
      int howOften = generator_[i]->howOften() ;
      /*  Probing can be set to just do column cuts in treee.
	  But if doing good then leave as on */
      bool probingWasOnBut=false;
      CglProbing * probing = dynamic_cast<CglProbing*>(generator_[i]->generator());
      if (probing&&!numberNodes_) {
	iProbing = i;
	if (probing->rowCuts()==-3) {
	  probingWasOnBut=true;
	  howOften = -98;
	  probing->setRowCuts(3);
	}
      }
      if (willBeCutsInTree<0&&howOften==-98)
        howOften =-99;
      if (howOften==-98&&generator_[i]->switchOffIfLessThan()<0) {
        if (thisObjective-startObjective<0.005*fabs(startObjective)+1.0e-5)
          howOften=-99; // switch off
        if (thisObjective-startObjective<0.1*fabs(startObjective)+1.0e-5
            &&5*solver_->getNumRows()<solver_->getNumCols())
          howOften=-99; // switch off
      }
      if (howOften<-99)
	continue ;
      if (howOften<0||howOften >= 1000000) {
        if( !numberNodes_) {
          // If small number switch mostly off
#ifdef JUST_ACTIVE
          double thisCuts = count[i] + 5.0*countColumnCuts[i] ;
#else
          double thisCuts = countRowCuts[i] + 5.0*countColumnCuts[i] ;
#endif
          if (!thisCuts||howOften == -99) {
            if (howOften == -99||howOften == -98) 
              howOften = -100 ;
            else
              howOften = 1000000+SCANCUTS; // wait until next time
          } else if (thisCuts<small) {
	    if (howOften!=1&&!probingWasOnBut) {
	      if (generator_[i]->whatDepth()<0||howOften!=-1) {
		int k = (int) sqrt(small/thisCuts) ;
		if (howOften!=-98)
		  howOften = k+1000000 ;
		else
		  howOften=-100;
	      } else {
		howOften=1;
	      }
	    } else {
	      howOften=1;
	      // allow cuts
	      probingWasOnBut=false;
	    }
          } else {
            if (thisObjective-startObjective<0.1*fabs(startObjective)+1.0e-5&&generator_[i]->whatDepth()<0)
              generator_[i]->setWhatDepth(5);
            howOften = 1+1000000 ;
          }
        }
        // If cuts useless switch off
        if (numberNodes_>=100000&&sumChangeObjective1_>2.0e2*(sumChangeObjective2_+1.0e-12)) {
          howOften = 1000000+SCANCUTS; // wait until next time
          //printf("switch off cut %d due to lack of use\n",i);
        }
      }
      if (!numberNodes_) {
	if (probingWasOnBut&&howOften==-100) {
	  probing->setRowCuts(-3);
	  howOften=1;
	}
	if (howOften==1)
	  generator_[i]->setWhatDepth(1);
	
	if (howOften>=0&&generator_[i]->generator()->mayGenerateRowCutsInTree())
	  willBeCutsInTree=1;
      }
      generator_[i]->setHowOften(howOften) ;
      if (howOften>=1000000&&howOften<2000000&&0) {
        // Go to depth
        int bias=1;
        if (howOften==1+1000000)
          generator_[i]->setWhatDepth(bias+1);
        else if (howOften<=10+1000000)
          generator_[i]->setWhatDepth(bias+2);
        else
          generator_[i]->setWhatDepth(bias+1000);
      }
      int newFrequency = generator_[i]->howOften()%1000000 ;
      // increment cut counts
      generator_[i]->incrementNumberCutsInTotal(countRowCuts[i]);
      generator_[i]->incrementNumberCutsActive(count[i]);
      if (handler_->logLevel()>1||!numberNodes_) {
	handler_->message(CBC_GENERATOR,messages_)
	  <<i
	  <<generator_[i]->cutGeneratorName()
	  <<countRowCuts[i]<<count[i]
	  <<countColumnCuts[i];
        handler_->printing(!numberNodes_&&generator_[i]->timing())
          <<generator_[i]->timeInCutGenerator();
        handler_->message()
	  <<newFrequency
	  <<CoinMessageEol ;
      }
    } 
    delete [] count ;
    if( !numberNodes_) {
      // save statistics
      for (i = 0;i<numberCutGenerators_;i++) {
	generator_[i]->setNumberCutsAtRoot(generator_[i]->numberCutsInTotal());
	generator_[i]->setNumberActiveCutsAtRoot(generator_[i]->numberCutsActive());
      }
      // decide on pseudo cost strategy
      int howOften = iProbing>=0 ? generator_[iProbing]->howOften() : 0;
      if ((howOften %1000000)!=1) 
	howOften = 0;
      //if (howOften) {
      //CglProbing * probing = dynamic_cast<CglProbing*>(generator_[iProbing]->generator());
      //}
      howOften = 0;
      if (howOften) {
	printf("** method 1\n");
	//CglProbing * probing = dynamic_cast<CglProbing*>(generator_[iProbing]->generator());
	generator_[iProbing]->setWhatDepth(1);
	// could set no row cuts
	//if (thisObjective-startObjective<0.001*fabs(startObjective)+1.0e-5)
	// probing->setRowCuts(0);
	for (int i=0;i<numberObjects_;i++) {
	  CbcSimpleIntegerDynamicPseudoCost * obj =
	    dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object_[i]) ;
	  if (obj) 
	    obj->setMethod(1);
	}
      }
      if (willBeCutsInTree==-2)
        willBeCutsInTree=0;
      if( willBeCutsInTree<=0) {
        // Take off cuts
        cuts = OsiCuts();
        numberNewCuts_=0;
        if (!willBeCutsInTree) {
          // update size of problem
          numberRowsAtContinuous_ = solver_->getNumRows() ;
        } else {
          // take off cuts
          int numberRows = solver_->getNumRows();
          int numberAdded = numberRows-numberRowsAtContinuous_;
          if (numberAdded) {
            int * added = new int[numberAdded];
            for (int i=0;i<numberAdded;i++)
              added[i]=i+numberRowsAtContinuous_;
            solver_->deleteRows(numberAdded,added);
            delete [] added;
            // resolve so optimal
            resolve(solver_);
          }
        }
#ifdef COIN_HAS_CLP
        OsiClpSolverInterface * clpSolver 
          = dynamic_cast<OsiClpSolverInterface *> (solver_);
	if (clpSolver) {
          // Maybe solver might like to know only column bounds will change
          //int options = clpSolver->specialOptions();
          //clpSolver->setSpecialOptions(options|128);
          clpSolver->synchronizeModel();
	}
#endif
      } else {
#ifdef COIN_HAS_CLP
        OsiClpSolverInterface * clpSolver 
          = dynamic_cast<OsiClpSolverInterface *> (solver_);
	if (clpSolver) {
        // make sure factorization can't carry over
          int options = clpSolver->specialOptions();
          clpSolver->setSpecialOptions(options&(~8));
	}
#endif
      }
    }
  } else {
#ifdef COIN_HAS_CLP
    OsiClpSolverInterface * clpSolver 
      = dynamic_cast<OsiClpSolverInterface *> (solver_);
    if (clpSolver) {
      // Maybe solver might like to know only column bounds will change
      //int options = clpSolver->specialOptions();
      //clpSolver->setSpecialOptions(options|128);
      clpSolver->synchronizeModel();
    }
#endif
    if (numberCutGenerators_) {
      int i;
      // add to counts anyway
      for (i = 0;i<numberCutGenerators_;i++) 
	generator_[i]->incrementNumberCutsInTotal(countRowCuts[i]);
      // What if not feasible as cuts may have helped
      if (feasible) {
	for (i = 0;i<numberNewCuts_;i++) {
	  int iGenerator = whichGenerator_[i];
	  if (iGenerator>=0)
	    generator_[iGenerator]->incrementNumberCutsActive();
	}
      }
    }
  }

  delete [] countRowCuts ;
  delete [] countColumnCuts ;


#ifdef CHECK_CUT_COUNTS
  if (feasible)
  { 
    CoinWarmStartBasis * basis = dynamic_cast<CoinWarmStartBasis*>(solver_->getWarmStart()) ;
    printf("solveWithCuts: Number of rows at end (only active cuts) %d\n",
	   numberRowsAtContinuous_+numberNewCuts_+numberOldActiveCuts_) ;
    basis->print() ; delete basis;}
#endif
#ifdef CBC_DEBUG
  if (onOptimalPath && !solver_->isDualObjectiveLimitReached())
    assert(feasible) ;
#endif
# ifdef COIN_HAS_CLP
  if (clpSolver) 
    clpSolver->setSpecialOptions(saveClpOptions);
# endif
#ifdef CBC_THREAD
  if (threadModel) {
    // stop threads
    int i;
    for (i=0;i<numberThreads_;i++) {
      while (threadInfo[i].returnCode==0) {
	pthread_cond_signal(threadInfo[i].condition2); // unlock
	pthread_mutex_lock(&condition_mutex);
	struct timespec absTime;
	clock_gettime(CLOCK_REALTIME,&absTime);
	absTime.tv_nsec += 1000000; // millisecond
	if (absTime.tv_nsec>=1000000000) {
	  absTime.tv_nsec -= 1000000000;
	  absTime.tv_sec++;
	}
	pthread_cond_timedwait(&condition_main,&condition_mutex,&absTime);
	clock_gettime(CLOCK_REALTIME,&absTime);
	pthread_mutex_unlock(&condition_mutex);
      }
      threadModel[i]->numberThreads_=0; // say exit
      threadInfo[i].returnCode=0;
      pthread_cond_signal(threadInfo[i].condition2); // unlock
      pthread_join(threadId[i],NULL);
      pthread_mutex_destroy (threadInfo[i].mutex2);
      pthread_cond_destroy (threadInfo[i].condition2);
      threadModel[i]->generator_[0]=NULL;
      delete [] threadModel[i]->generator_;
      threadModel[i]->generator_=NULL;
      threadModel[i]->object_=NULL;
    }
    pthread_cond_destroy (&condition_main);
    pthread_mutex_destroy (&condition_mutex);
    // delete models and solvers
    for (i=0;i<numberThreads_;i++) {
      delete threadModel[i];
    }
    delete [] mutex2;
    delete [] condition2;
    delete [] threadId;
    delete [] threadInfo;
    delete [] threadModel;
    mutex_ = saveMutex;
  }
#endif

  return feasible ; }


/*
  Remove slack cuts. We obtain a basis and scan it. Cuts with basic slacks
  are purged. If any cuts are purged, resolve() is called to restore the
  solution held in the solver.	If resolve() pivots, there's the possibility
  that a slack may be pivoted in (trust me :-), so the process iterates.
  Setting allowResolve to false will suppress reoptimisation (but see note
  below).

  At the level of the solver's constraint system, loose cuts are really
  deleted.  There's an implicit assumption that deleteRows will also update
  the active basis in the solver.

  At the level of nodes and models, it's more complicated.

  New cuts exist only in the collection of cuts passed as a parameter. They
  are deleted from the collection and that's the end of them.

  Older cuts have made it into addedCuts_. Two separate actions are needed.
  The reference count for the CbcCountRowCut object is decremented. If this
  count falls to 0, the node which owns the cut is located, the reference to
  the cut is removed, and then the cut object is destroyed (courtesy of the
  CbcCountRowCut destructor). We also need to set the addedCuts_ entry to
  NULL. This is important so that when it comes time to generate basis edits
  we can tell this cut was dropped from the basis during processing of the
  node.

  NOTE: In general, it's necessary to call resolve() after purging slack
	cuts.  Deleting constraints constitutes a change in the problem, and
	an OSI is not required to maintain a valid solution when the problem
	is changed. But ... it's really useful to maintain the active basis,
	and the OSI is supposed to do that. (Yes, it's splitting hairs.) In
	some places, it's possible to know that the solution will never be
	consulted after this call, only the basis.  (E.g., this routine is
	called as a last act before generating info to place the node in the
	live set.) For such use, set allowResolve to false.
  
  TODO: No real harm would be done if we just ignored the rare occasion when
	the call to resolve() pivoted a slack back into the basis. It's a
	minor inefficiency, at worst. But it does break assertions which
	check that there are no loose cuts in the basis. It might be better
	to remove the assertions.
*/

void
CbcModel::takeOffCuts (OsiCuts &newCuts,
		       bool allowResolve, OsiCuts * saveCuts)

{ // int resolveIterations = 0 ;
  int firstOldCut = numberRowsAtContinuous_ ;
  int totalNumberCuts = numberNewCuts_+numberOldActiveCuts_ ;
  int *solverCutIndices = new int[totalNumberCuts] ;
  int *newCutIndices = new int[numberNewCuts_] ;
  const CoinWarmStartBasis* ws ;
  CoinWarmStartBasis::Status status ;
  bool needPurge = true ;
/*
  The outer loop allows repetition of purge in the event that reoptimisation
  changes the basis. To start an iteration, clear the deletion counts and grab
  the current basis.
*/
  while (needPurge)
  { int numberNewToDelete = 0 ;
    int numberOldToDelete = 0 ;
    int i ;
    ws = dynamic_cast<const CoinWarmStartBasis*>(solver_->getWarmStart()) ;
/*
  Scan the basis entries of the old cuts generated prior to this round of cut
  generation.  Loose cuts are `removed' by decrementing their reference count
  and setting the addedCuts_ entry to NULL. (If the reference count falls to
  0, they're really deleted.  See CbcModel and CbcCountRowCut doc'n for
  principles of cut handling.)
*/
    int oldCutIndex = 0 ;
    if (numberOldActiveCuts_) {
      lockThread();
      for (i = 0 ; i < numberOldActiveCuts_ ; i++)
	{ status = ws->getArtifStatus(i+firstOldCut) ;
	while (!addedCuts_[oldCutIndex]) oldCutIndex++ ;
	assert(oldCutIndex < currentNumberCuts_) ;
	// always leave if from nextRowCut_
	if (status == CoinWarmStartBasis::basic&&
	    addedCuts_[oldCutIndex]->effectiveness()!=COIN_DBL_MAX)
	  { solverCutIndices[numberOldToDelete++] = i+firstOldCut ;
	  if (saveCuts) {
	    // send to cut pool
	    OsiRowCut * slackCut = addedCuts_[oldCutIndex];
	    if (slackCut->effectiveness()!=-1.234) {
	      slackCut->setEffectiveness(-1.234);
	      saveCuts->insert(*slackCut);
	    }
	  }
	  if (addedCuts_[oldCutIndex]->decrement() == 0)
	    delete addedCuts_[oldCutIndex] ;
	  addedCuts_[oldCutIndex] = NULL ;
	  oldCutIndex++ ; }
	else
	  { oldCutIndex++ ; } }
      unlockThread();
    }
/*
  Scan the basis entries of the new cuts generated with this round of cut
  generation.  At this point, newCuts is the only record of the new cuts, so
  when we delete loose cuts from newCuts, they're really gone. newCuts is a
  vector, so it's most efficient to compress it (eraseRowCut) from back to
  front.
*/
    int firstNewCut = firstOldCut+numberOldActiveCuts_ ;
    int k = 0 ;
    for (i = 0 ; i < numberNewCuts_ ; i++)
    { status = ws->getArtifStatus(i+firstNewCut) ;
      if (status == CoinWarmStartBasis::basic&&whichGenerator_[i]!=-2)
      { solverCutIndices[numberNewToDelete+numberOldToDelete] = i+firstNewCut ;
	newCutIndices[numberNewToDelete++] = i ; }
      else
      { // save which generator did it
	whichGenerator_[k++] = whichGenerator_[i] ; } }
    delete ws ;
    for (i = numberNewToDelete-1 ; i >= 0 ; i--)
    { int iCut = newCutIndices[i] ;
      if (saveCuts) {
        // send to cut pool
        OsiRowCut * slackCut = newCuts.rowCutPtr(iCut);
        if (slackCut->effectiveness()!=-1.234) {
          slackCut->setEffectiveness(-1.234);
          saveCuts->insert(*slackCut);
        }
      }
      newCuts.eraseRowCut(iCut) ; }
/*
  Did we delete anything? If so, delete the cuts from the constraint system
  held in the solver and reoptimise unless we're forbidden to do so. If the
  call to resolve() results in pivots, there's the possibility we again have
  basic slacks. Repeat the purging loop.
*/
    if (numberNewToDelete+numberOldToDelete > 0)
    { solver_->deleteRows(numberNewToDelete+numberOldToDelete,
			  solverCutIndices) ;
      numberNewCuts_ -= numberNewToDelete ;
      numberOldActiveCuts_ -= numberOldToDelete ;
#     ifdef CBC_DEBUG
      printf("takeOffCuts: purged %d+%d cuts\n", numberOldToDelete,
	     numberNewToDelete );
#     endif
      if (allowResolve)
      { 
	phase_=3;
        // can do quick optimality check
        int easy=2;
        solver_->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,&easy) ;
	resolve(solver_) ;
        setPointers(solver_);
        solver_->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,NULL) ;
	if (solver_->getIterationCount() == 0)
	{ needPurge = false ; }
#	ifdef CBC_DEBUG
	else
	  { printf( "Repeating purging loop. %d iters.\n",
		    solver_->getIterationCount()); }
#	endif
      }
      else
      { needPurge = false ; } }
    else
    { needPurge = false ; } }
/*
  Clean up and return.
*/
  delete [] solverCutIndices ;
  delete [] newCutIndices ;
}



int
CbcModel::resolve(CbcNodeInfo * parent, int whereFrom)
{
  // We may have deliberately added in violated cuts - check to avoid message
  int iRow;
  int numberRows = solver_->getNumRows();
  const double * rowLower = solver_->getRowLower();
  const double * rowUpper = solver_->getRowUpper();
  bool feasible=true;
  for (iRow= numberRowsAtContinuous_;iRow<numberRows;iRow++) {
    if (rowLower[iRow]>rowUpper[iRow]+1.0e-8)
      feasible=false;
  }
  // Can't happen if strong branching as would have been found before
  if (!numberStrong_&&numberObjects_>numberIntegers_) {
    int iColumn;
    int numberColumns = solver_->getNumCols();
    const double * columnLower = solver_->getColLower();
    const double * columnUpper = solver_->getColUpper();
    for (iColumn= 0;iColumn<numberColumns;iColumn++) {
      if (columnLower[iColumn]>columnUpper[iColumn]+1.0e-5)
        feasible=false;
    }
  }
#if 0
  {
    int iColumn;
    int numberColumns = solver_->getNumCols();
    const double * columnLower = solver_->getColLower();
    const double * columnUpper = solver_->getColUpper();
    const double * sol = solver_->getColSolution();
    bool feasible=true;
    double xx[37];
    memset(xx,0,sizeof(xx));
    xx[6]=xx[7]=xx[9]=xx[18]=xx[26]=xx[36]=1.0;
    for (iColumn= 0;iColumn<numberColumns;iColumn++) {
      if (solver_->isInteger(iColumn)) {
        //printf("yy %d at %g (bounds %g, %g)\n",iColumn,
        //     sol[iColumn],columnLower[iColumn],columnUpper[iColumn]);
        if (columnLower[iColumn]>xx[iColumn]||
            columnUpper[iColumn]<xx[iColumn])
          feasible=false;
        //solver_->setColLower(iColumn,CoinMin(xx[iColumn],columnUpper[iColumn]));
        //solver_->setColUpper(iColumn,CoinMax(xx[iColumn],columnLower[iColumn]));
      }
    }
    if (feasible)
      printf("On Feasible path\n");
  }
#endif
/*
  Reoptimize. Consider the possibility that we should fathom on bounds. But be
  careful --- where the objective takes on integral values, we may want to keep
  a solution where the objective is right on the cutoff.
*/
  if (feasible)
    {
      resolve(solver_) ;
      numberIterations_ += solver_->getIterationCount() ;
      feasible = (solver_->isProvenOptimal() &&
                  !solver_->isDualObjectiveLimitReached()) ;
    }
  if (0&&feasible) {
    const double * lb = solver_->getColLower();
    const double * ub = solver_->getColUpper();
    const double * x = solver_->getColSolution();
    const double * dj = solver_->getReducedCost();
    int numberColumns = solver_->getNumCols();
    for (int i=0;i<numberColumns;i++) {
      if (dj[i]>1.0e-4&&ub[i]-lb[i]>1.0e-4&&x[i]>lb[i]+1.0e-4)
        printf("error %d %g %g %g %g\n",i,dj[i],lb[i],x[i],ub[i]);
      if (dj[i]<-1.0e-4&&ub[i]-lb[i]>1.0e-4&&x[i]<ub[i]-1.0e-4)
        printf("error %d %g %g %g %g\n",i,dj[i],lb[i],x[i],ub[i]);
    }
  } 
  if (false&&!feasible&& continuousObjective_ <-1.0e30) {
    // at root node - double double check
    bool saveTakeHint;
    OsiHintStrength saveStrength;
    solver_->getHintParam(OsiDoDualInResolve,saveTakeHint,saveStrength);
    if (saveTakeHint||saveStrength==OsiHintIgnore) {
      solver_->setHintParam(OsiDoDualInResolve,false,OsiHintDo) ;
      resolve(solver_);
      solver_->setHintParam(OsiDoDualInResolve,saveTakeHint,saveStrength);
      numberIterations_ += solver_->getIterationCount() ;
      feasible = solver_->isProvenOptimal();
      //      solver_->writeMps("infeas");
    }
  }
  if (cutModifier_&&feasible&&!solverCharacteristics_->solutionAddsCuts()) {
    //double increment = getDblParam(CbcModel::CbcCutoffIncrement) ;
    double cutoff ;
    solver_->getDblParam(OsiDualObjectiveLimit,cutoff) ;
    double distance = fabs(cutoff-solver_->getObjValue());
    if (distance<10.0*trueIncrement) {
      double offset;
      solver_->getDblParam(OsiObjOffset,offset);
      double objFixedValue = -offset;
      double objValue=0.0;
      double direction = solver_->getObjSense();
      const double * solution = solver_->getColSolution();
      const double * objective = solver_->getObjCoefficients();
      const double * columnLower = solver_->getColLower();
      const double * columnUpper = solver_->getColUpper();
      int numberColumns = solver_->getNumCols();
      int increment = 0 ;
      double multiplier = 1.0/trueIncrement;
      int bigIntegers = 0; // Count of large costs which are integer
      for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	double value = solution[iColumn];
	// make sure clean
	value = CoinMin(value,columnUpper[iColumn]);
	value = CoinMax(value,columnLower[iColumn]);
	double cost = direction * objective[iColumn];
	if (cost) {
	  if (columnLower[iColumn]<columnUpper[iColumn]) {
	    objValue += value*cost;
	    value = fabs(cost)*multiplier ;
	    if (value <2.1e9) {
	      int nearest = (int) floor(value+0.5) ;
	      assert (fabs(value-floor(value+0.5)) < 1.0e-8);
	      if (!increment) 
		increment = nearest ; 
	      else
		increment = gcd(increment,nearest) ;
	    } else {
	      // large value - may still be multiple of 1.0
	      value = fabs(objective[iColumn]);
	      assert(fabs(value-floor(value+0.5)) < 1.0e-8);
	      bigIntegers++;
	    }
	  } else {
	    // fixed
	    objFixedValue += value*cost;
	  }
	}
      }
      if (increment) {
	double value = increment ;
	value /= multiplier ;
	if (value>trueIncrement) {
	  double x = objValue/value;
	  x = ceil(x-1.0e-5);
	  x *= value;
	  //printf("fixed %g, variable %g -> %g, sum %g - cutoff %g\n",
	  // objFixedValue,objValue,x,x+objFixedValue,cutoff);
	  x += objFixedValue;
	  if (x>cutoff + 1.0e-5*fabs(cutoff)+1.0e-5) {
	    //printf("Node cutoff\n");
	    feasible=false;
	  }
	} else {
	  value = trueIncrement;
	  double x = objValue/value;
	  x = ceil(x-1.0e-5);
	  x *= value;
	  x += objFixedValue;
	  if (x>cutoff + 1.0e-5*fabs(cutoff)+1.0e-5) {
	    //printf("Node cutoff\n");
	    feasible=false;
	  }
	}
      }
    }
  }

  setPointers(solver_);
  int returnStatus = feasible ? 1 : 0;
  if (strategy_) {
    // user can play clever tricks here
    int status = strategy_->status(this,parent,whereFrom);
    if (status>=0) {
      if (status==0)
        returnStatus = 1;
      else if(status==1)
        returnStatus=-1;
      else
        returnStatus=0;
    }
  }
  return returnStatus ;
}


/* Set up objects.  Only do ones whose length is in range.
   If makeEquality true then a new model may be returned if
   modifications had to be made, otherwise "this" is returned.

   Could use Probing at continuous to extend objects
*/
CbcModel * 
CbcModel::findCliques(bool makeEquality,
		      int atLeastThisMany, int lessThanThis, int defaultValue)
{
  // No objects are allowed to exist
  assert(numberObjects_==numberIntegers_||!numberObjects_);
  CoinPackedMatrix matrixByRow(*solver_->getMatrixByRow());
  int numberRows = solver_->getNumRows();
  int numberColumns = solver_->getNumCols();

  // We may want to add columns
  int numberSlacks=0;
  int * rows = new int[numberRows];
  double * element =new double[numberRows];

  int iRow;

  findIntegers(true);
  numberObjects_=numberIntegers_;

  int numberCliques=0;
  OsiObject ** object = new OsiObject * [numberRows];
  int * which = new int[numberIntegers_];
  char * type = new char[numberIntegers_];
  int * lookup = new int[numberColumns];
  int i;
  for (i=0;i<numberColumns;i++) 
    lookup[i]=-1;
  for (i=0;i<numberIntegers_;i++) 
    lookup[integerVariable_[i]]=i;

  // Statistics
  int totalP1=0,totalM1=0;
  int numberBig=0,totalBig=0;
  int numberFixed=0;

  // Row copy
  const double * elementByRow = matrixByRow.getElements();
  const int * column = matrixByRow.getIndices();
  const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
  const int * rowLength = matrixByRow.getVectorLengths();

  // Column lengths for slacks
  const int * columnLength = solver_->getMatrixByCol()->getVectorLengths();

  const double * lower = getColLower();
  const double * upper = getColUpper();
  const double * rowLower = getRowLower();
  const double * rowUpper = getRowUpper();

  for (iRow=0;iRow<numberRows;iRow++) {
    int numberP1=0, numberM1=0;
    int j;
    double upperValue=rowUpper[iRow];
    double lowerValue=rowLower[iRow];
    bool good=true;
    int slack = -1;
    for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
      int iColumn = column[j];
      int iInteger=lookup[iColumn];
      if (upper[iColumn]-lower[iColumn]<1.0e-8) {
	// fixed
	upperValue -= lower[iColumn]*elementByRow[j];
	lowerValue -= lower[iColumn]*elementByRow[j];
	continue;
      } else if (upper[iColumn]!=1.0||lower[iColumn]!=0.0) {
	good = false;
	break;
      } else if (iInteger<0) {
	good = false;
	break;
      } else {
	if (columnLength[iColumn]==1)
	  slack = iInteger;
      }
      if (fabs(elementByRow[j])!=1.0) {
	good=false;
	break;
      } else if (elementByRow[j]>0.0) {
	which[numberP1++]=iInteger;
      } else {
	numberM1++;
	which[numberIntegers_-numberM1]=iInteger;
      }
    }
    int iUpper = (int) floor(upperValue+1.0e-5);
    int iLower = (int) ceil(lowerValue-1.0e-5);
    int state=0;
    if (upperValue<1.0e6) {
      if (iUpper==1-numberM1)
	state=1;
      else if (iUpper==-numberM1)
	state=2;
      else if (iUpper<-numberM1)
	state=3;
    }
    if (!state&&lowerValue>-1.0e6) {
      if (-iLower==1-numberP1)
	state=-1;
      else if (-iLower==-numberP1)
	state=-2;
      else if (-iLower<-numberP1)
	state=-3;
    }
    if (good&&state) {
      if (abs(state)==3) {
	// infeasible
	numberObjects_=-1;
	break;
      } else if (abs(state)==2) {
	// we can fix all
	numberFixed += numberP1+numberM1;
	if (state>0) {
	  // fix all +1 at 0, -1 at 1
	  for (i=0;i<numberP1;i++)
	    solver_->setColUpper(integerVariable_[which[i]],0.0);
	  for (i=0;i<numberM1;i++)
	    solver_->setColLower(integerVariable_[which[numberIntegers_-i-1]],
				 1.0);
	} else {
	  // fix all +1 at 1, -1 at 0
	  for (i=0;i<numberP1;i++)
	    solver_->setColLower(integerVariable_[which[i]],1.0);
	  for (i=0;i<numberM1;i++)
	    solver_->setColUpper(integerVariable_[which[numberIntegers_-i-1]],
				 0.0);
	}
      } else {
	int length = numberP1+numberM1;
	if (length >= atLeastThisMany&&length<lessThanThis) {
	  // create object
	  bool addOne=false;
	  int objectType;
	  if (iLower==iUpper) {
	    objectType=1;
	  } else {
	    if (makeEquality) {
	      objectType=1;
	      element[numberSlacks]=state;
	      rows[numberSlacks++]=iRow;
	      addOne=true;
	    } else {
	      objectType=0;
	    }
	  }
	  if (state>0) {
	    totalP1 += numberP1;
	    totalM1 += numberM1;
	    for (i=0;i<numberP1;i++)
	      type[i]=1;
	    for (i=0;i<numberM1;i++) {
	      which[numberP1]=which[numberIntegers_-i-1];
	      type[numberP1++]=0;
	    }
	  } else {
	    totalP1 += numberM1;
	    totalM1 += numberP1;
	    for (i=0;i<numberP1;i++)
	      type[i]=0;
	    for (i=0;i<numberM1;i++) {
	      which[numberP1]=which[numberIntegers_-i-1];
	      type[numberP1++]=1;
	    }
	  }
	  if (addOne) {
	    // add in slack
	    which[numberP1]=numberIntegers_+numberSlacks-1;
	    slack = numberP1;
	    type[numberP1++]=1;
	  } else if (slack >= 0) {
	    for (i=0;i<numberP1;i++) {
	      if (which[i]==slack) {
		slack=i;
	      }
	    }
	  }
	  object[numberCliques] = new CbcClique(this,objectType,numberP1,
					      which,type,
					       1000000+numberCliques,slack);
          numberCliques++;
	} else if (numberP1+numberM1 >= lessThanThis) {
	  // too big
	  numberBig++;
	  totalBig += numberP1+numberM1;
	}
      }
    }
  }
  delete [] which;
  delete [] type;
  delete [] lookup;
#ifdef COIN_DEVELOP
  if (numberCliques<0) {
    printf("*** Problem infeasible\n");
  } else {
    if (numberCliques)
      printf("%d cliques of average size %g found, %d P1, %d M1\n",
	     numberCliques,
	     ((double)(totalP1+totalM1))/((double) numberCliques),
	     totalP1,totalM1);
    else
      printf("No cliques found\n");
    if (numberBig)
      printf("%d large cliques ( >= %d) found, total %d\n",
	     numberBig,lessThanThis,totalBig);
    if (numberFixed)
      printf("%d variables fixed\n",numberFixed);
  }
#endif
  if (numberCliques>0&&numberSlacks&&makeEquality) {
    printf("adding %d integer slacks\n",numberSlacks);
    // add variables to make equality rows
    int * temp = new int[numberIntegers_+numberSlacks];
    memcpy(temp,integerVariable_,numberIntegers_*sizeof(int));
    // Get new model
    CbcModel * newModel = new CbcModel(*this);
    OsiSolverInterface * newSolver = newModel->solver();
    for (i=0;i<numberSlacks;i++) {
      temp[i+numberIntegers_]=i+numberColumns;
      int iRow = rows[i];
      double value = element[i];
      double lowerValue = 0.0;
      double upperValue = 1.0;
      double objValue  = 0.0;
      CoinPackedVector column(1,&iRow,&value);
      newSolver->addCol(column,lowerValue,upperValue,objValue);
      // set integer
      newSolver->setInteger(numberColumns+i);
      if (value >0)
	newSolver->setRowLower(iRow,rowUpper[iRow]);
      else
	newSolver->setRowUpper(iRow,rowLower[iRow]);
    }
    // replace list of integers
    for (i=0;i<newModel->numberObjects_;i++)
      delete newModel->object_[i];
    newModel->numberObjects_ = 0;
    delete [] newModel->object_;
    newModel->object_=NULL;
    newModel->findIntegers(true); //Set up all integer objects
    for (i=0;i<numberIntegers_;i++) {
      newModel->modifiableObject(i)->setPriority(object_[i]->priority());
    }
    if (originalColumns_) {
      // old model had originalColumns
      delete [] newModel->originalColumns_;
      newModel->originalColumns_ = new int[numberColumns+numberSlacks];
      memcpy(newModel->originalColumns_,originalColumns_,numberColumns*sizeof(int));
      // mark as not in previous model
      for (i=numberColumns;i<numberColumns+numberSlacks;i++)
	newModel->originalColumns_[i]=-1;
    }
    delete [] rows;
    delete [] element;
    newModel->addObjects(numberCliques,object);
    assert (ownObjects_);
    for (;i<numberCliques;i++) 
      delete object[i];
    delete [] object;
    newModel->synchronizeModel();
    return newModel;
  } else {
    assert (ownObjects_);
    if (numberCliques>0) {
      addObjects(numberCliques,object);
      for (;i<numberCliques;i++) 
	delete object[i];
      synchronizeModel();
    }
    delete [] object;
    delete [] rows;
    delete [] element;
    return this;
  }
}
// Fill in useful estimates
void 
CbcModel::pseudoShadow(double * down, double * up)
{
  // Column copy of matrix
  const double * element = solver_->getMatrixByCol()->getElements();
  const int * row = solver_->getMatrixByCol()->getIndices();
  const CoinBigIndex * columnStart = solver_->getMatrixByCol()->getVectorStarts();
  const int * columnLength = solver_->getMatrixByCol()->getVectorLengths();
  const double *objective = solver_->getObjCoefficients() ;
  int numberColumns = solver_->getNumCols() ;
  double direction = solver_->getObjSense();
  int iColumn;
  const double * dual = cbcRowPrice_;
  down = new double[numberColumns];
  up = new double[numberColumns];
  double upSum=1.0e-20;
  double downSum = 1.0e-20;
  int numberIntegers=0;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    CoinBigIndex start = columnStart[iColumn];
    CoinBigIndex end = start + columnLength[iColumn];
    double upValue = 0.0;
    double downValue = 0.0;
    double value = direction*objective[iColumn];
    if (value) {
      if (value>0.0)
        upValue += value;
      else
        downValue -= value;
    }
    for (CoinBigIndex j=start;j<end;j++) {
      int iRow = row[j];
      value = -dual[iRow];
      if (value) {
        value *= element[j];
        if (value>0.0)
          upValue += value;
        else
          downValue -= value;
      }
    }
    // use dj if bigger
    double dj = cbcReducedCost_[iColumn];
    upValue = CoinMax(upValue,dj);
    downValue = CoinMax(downValue,-dj);
    up[iColumn]=upValue;
    down[iColumn]=downValue;
    if (solver_->isInteger(iColumn)) {
      if (!numberNodes_)
        printf("%d - dj %g up %g down %g cost %g\n",
               iColumn,dj,upValue,downValue,objective[iColumn]);
      upSum += upValue;
      downSum += downValue;
      numberIntegers++;
    }
  }
  if (numberIntegers) {
    double smallDown = 0.01*(downSum/((double) numberIntegers));
    double smallUp = 0.01*(upSum/((double) numberIntegers));
    for (int i=0;i<numberObjects_;i++) {
      CbcSimpleIntegerDynamicPseudoCost * obj1 =
        dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object_[i]) ;
      if (obj1) {
        iColumn = obj1->columnNumber();
        double upPseudoCost = obj1->upDynamicPseudoCost();
        double saveUp = upPseudoCost;
        upPseudoCost = CoinMax(upPseudoCost,smallUp);
        upPseudoCost = CoinMax(upPseudoCost,up[iColumn]);
        upPseudoCost = CoinMax(upPseudoCost,0.1*down[iColumn]);
        obj1->setUpDynamicPseudoCost(upPseudoCost);
        if (upPseudoCost>saveUp&&!numberNodes_)
          printf("For %d up went from %g to %g\n",
                 iColumn,saveUp,upPseudoCost);
        double downPseudoCost = obj1->downDynamicPseudoCost();
        double saveDown = downPseudoCost;
        downPseudoCost = CoinMax(downPseudoCost,smallDown);
        downPseudoCost = CoinMax(downPseudoCost,down[iColumn]);
        downPseudoCost = CoinMax(downPseudoCost,0.1*down[iColumn]);
        obj1->setDownDynamicPseudoCost(downPseudoCost);
        if (downPseudoCost>saveDown&&!numberNodes_)
          printf("For %d down went from %g to %g\n",
                 iColumn,saveDown,downPseudoCost);
      }
    }
  }
  delete [] down;
  delete [] up;
}

/*
  Set branching priorities.

  Setting integer priorities looks pretty robust; the call to findIntegers
  makes sure that SimpleInteger objects are in place. Setting priorities for
  other objects is entirely dependent on their existence, and the routine may
  quietly fail in several directions.
*/

void 
CbcModel::passInPriorities (const int * priorities,
			    bool ifObject)
{
  findIntegers(false);
  int i;
  if (priorities) {
    int i0=0;
    int i1=numberObjects_-1;
    if (ifObject) {
      for (i=numberIntegers_;i<numberObjects_;i++) {
        object_[i]->setPriority(priorities[i-numberIntegers_]);
      }
      i0=numberIntegers_;
    } else {
      for (i=0;i<numberIntegers_;i++) {
        object_[i]->setPriority(priorities[i]);
      }
      i1=numberIntegers_-1;
    }
    messageHandler()->message(CBC_PRIORITY,
                              messages())
                                << i0<<i1<<numberObjects_ << CoinMessageEol ;
  }
}

// Delete all object information
void 
CbcModel::deleteObjects(bool getIntegers)
{
  if (ownObjects_) {
    int i;
    for (i=0;i<numberObjects_;i++)
      delete object_[i];
    delete [] object_;
  }
  object_ = NULL;
  numberObjects_=0;
  if (getIntegers&&ownObjects_)
    findIntegers(true);
}

/*!
  Ensure all attached objects (OsiObjects, heuristics, and cut
  generators) point to this model.
*/
void CbcModel::synchronizeModel()
{
  int i;
  for (i=0;i<numberHeuristics_;i++) 
    heuristic_[i]->setModel(this);
  for (i=0;i<numberObjects_;i++) {
    CbcObject * obj =
      dynamic_cast <CbcObject *>(object_[i]) ;
    if (obj)
      obj->setModel(this);
  }
  for (i=0;i<numberCutGenerators_;i++)
    generator_[i]->refreshModel(this);
}

// Fill in integers and create objects

/**
  The routine first does a scan to count the number of integer variables.
  It then creates an array, integerVariable_, to store the indices of the
  integer variables, and an array of `objects', one for each variable.

  The scan is repeated, this time recording the index of each integer
  variable in integerVariable_, and creating an CbcSimpleInteger object that
  contains information about the integer variable. Initially, this is just
  the index and upper & lower bounds.

  \todo
  Note the assumption in cbc that the first numberIntegers_ objects are
  CbcSimpleInteger. In particular, the code which handles the startAgain
  case assumes that if the object_ array exists it can simply replace the first
  numberInteger_ objects. This is arguably unsafe.

  I am going to re-order if necessary
*/

void 
CbcModel::findIntegers(bool startAgain,int type)
{
  assert(solver_);
/*
  No need to do this if we have previous information, unless forced to start
  over.
*/
  if (numberIntegers_&&!startAgain&&object_)
    return;
/*
  Clear out the old integer variable list, then count the number of integer
  variables.
*/
  delete [] integerVariable_;
  integerVariable_ = NULL;
  numberIntegers_=0;
  int numberColumns = getNumCols();
  int iColumn;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (isInteger(iColumn)) 
      numberIntegers_++;
  }
  // Find out how many old non-integer objects there are
  int nObjects=0;
  OsiObject ** oldObject = object_;
  int iObject;
  for (iObject = 0;iObject<numberObjects_;iObject++) {
    CbcSimpleInteger * obj =
      dynamic_cast <CbcSimpleInteger *>(oldObject[iObject]) ;
    if (obj) 
      delete oldObject[iObject];
    else
      oldObject[nObjects++]=oldObject[iObject];
  }
  // See if there any SOS
#ifdef COIN_HAS_CLP
  if (!nObjects) {
    OsiClpSolverInterface * clpSolver 
      = dynamic_cast<OsiClpSolverInterface *> (solver_);
    if (clpSolver&&clpSolver->numberSOS()) {
      // deal with sos
      const CoinSet * setInfo = clpSolver->setInfo();
      int numberSOS = clpSolver->numberSOS();
      nObjects=0;
      delete [] oldObject;
      oldObject = new OsiObject * [numberSOS];
      for (int i=0;i<numberSOS;i++) {
	int type = setInfo[i].setType();
	int n=setInfo[i].numberEntries();
	const int * which = setInfo[i].which();
	const double * weights = setInfo[i].weights();
	oldObject[nObjects++] = new CbcSOS(this,n,which,weights,i,type);
      }
    }
  }
#endif
    
/*
  Found any? Allocate an array to hold the indices of the integer variables.
  Make a large enough array for all objects
*/
  delete [] integerVariable_;
  object_ = new OsiObject * [numberIntegers_+nObjects];
  numberObjects_=numberIntegers_+nObjects;
  integerVariable_ = new int [numberIntegers_];
/*
  Walk the variables again, filling in the indices and creating objects for
  the integer variables. Initially, the objects hold the index and upper &
  lower bounds.
*/
  numberIntegers_=0;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if(isInteger(iColumn)) {
      if (!type)
        object_[numberIntegers_] =
          new CbcSimpleInteger(this,iColumn);
      else if (type==1)
        object_[numberIntegers_] =
          new CbcSimpleIntegerPseudoCost(this,iColumn,0.3);
      integerVariable_[numberIntegers_++]=iColumn;
    }
  }
  // Now append other objects
  memcpy(object_+numberIntegers_,oldObject,nObjects*sizeof(OsiObject *));
  // Delete old array (just array)
  delete [] oldObject;
  
  if (!numberObjects_)
      handler_->message(CBC_NOINT,messages_) << CoinMessageEol ;
}
/* If numberBeforeTrust >0 then we are going to use CbcBranchDynamic.
   Scan and convert CbcSimpleInteger objects
*/
void 
CbcModel::convertToDynamic()
{
  int iObject;
  const double * cost = solver_->getObjCoefficients();
  bool allDynamic=true;
  for (iObject = 0;iObject<numberObjects_;iObject++) {
    CbcSimpleInteger * obj1 =
      dynamic_cast <CbcSimpleInteger *>(object_[iObject]) ;
    CbcSimpleIntegerPseudoCost * obj1a =
      dynamic_cast <CbcSimpleIntegerPseudoCost *>(object_[iObject]) ;
    CbcSimpleIntegerDynamicPseudoCost * obj2 =
      dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object_[iObject]) ;
    if (obj1&&!obj2) {
      // replace
      int iColumn = obj1->columnNumber();
      int priority = obj1->priority();
      int preferredWay = obj1->preferredWay();
      delete object_[iObject];
      double costValue = CoinMax(1.0e-5,fabs(cost[iColumn]));
      // treat as if will cost what it says up
      double upCost=costValue;
      // and balance at breakeven of 0.3
      double downCost=(0.7*upCost)/0.3;
      if (obj1a) {
        upCost=obj1a->upPseudoCost();
        downCost=obj1a->downPseudoCost();
      }
      CbcSimpleIntegerDynamicPseudoCost * newObject =
        new CbcSimpleIntegerDynamicPseudoCost(this,iColumn,1.0e0*downCost,1.0e0*upCost);
      newObject->setNumberBeforeTrust(numberBeforeTrust_);
      newObject->setPriority(priority);
      newObject->setPreferredWay(preferredWay);
      object_[iObject] = newObject;
    } else if (!obj2) {
      allDynamic=false;
    } else {
      // synchronize trust
      obj2->setNumberBeforeTrust(numberBeforeTrust_);
    }
  }
  if (branchingMethod_) {
    if ((branchingMethod_->whichMethod()&1)==0&&!branchingMethod_->chooseMethod()) {
      // Need a method which can do better
      delete branchingMethod_;
      branchingMethod_=NULL;
    }
  }
  if (!branchingMethod_&&allDynamic) {
    // create one
    branchingMethod_ = new CbcBranchDynamicDecision();
  }
}

/* Add in any object information (objects are cloned - owner can delete
   originals */
void 
CbcModel::addObjects(int numberObjects, CbcObject ** objects)
{
 // If integers but not enough objects fudge
  if (numberIntegers_>numberObjects_)
    findIntegers(true);
  /* But if incoming objects inherit from simple integer we just want
     to replace */
  int numberColumns = solver_->getNumCols();
  /** mark is -1 if not integer, >=0 if using existing simple integer and
      >=numberColumns if using new integer */
  int * mark = new int[numberColumns];
  int i;
  for (i=0;i<numberColumns;i++)
    mark[i]=-1;
  int newNumberObjects = numberObjects;
  int newIntegers=0;
  for (i=0;i<numberObjects;i++) { 
    CbcSimpleInteger * obj =
      dynamic_cast <CbcSimpleInteger *>(objects[i]) ;
    if (obj) {
      int iColumn = obj->columnNumber();
      mark[iColumn]=i+numberColumns;
      newIntegers++;
    }
  }
  // and existing
  for (i=0;i<numberObjects_;i++) { 
    CbcSimpleInteger * obj =
      dynamic_cast <CbcSimpleInteger *>(object_[i]) ;
    if (obj) {
      int iColumn = obj->columnNumber();
      if (mark[iColumn]<0) {
        newIntegers++;
        newNumberObjects++;
        mark[iColumn]=i;
      }
    }
  } 
  delete [] integerVariable_;
  integerVariable_=NULL;
#ifdef COIN_DEVELOP
  if (newIntegers!=numberIntegers_) 
    printf("changing number of integers from %d to %d\n",
           numberIntegers_,newIntegers);
#endif
  numberIntegers_ = newIntegers;
  integerVariable_ = new int [numberIntegers_];
  OsiObject ** temp  = new OsiObject * [newNumberObjects];
  // Put integers first
  newIntegers=0;
  numberIntegers_=0;
  for (i=0;i<numberColumns;i++) {
    int which = mark[i];
    if (which>=0) {
      if (!isInteger(i)) {
        newIntegers++;
        solver_->setInteger(i);
      }
      if (which<numberColumns) {
        temp[numberIntegers_]=object_[which];
        object_[which]=NULL;
      } else {
        temp[numberIntegers_]=objects[which-numberColumns]->clone();
      }
      integerVariable_[numberIntegers_++]=i;
    }
  }
#ifdef COIN_DEVELOP
  if (newIntegers)
    printf("%d variables were declared integer\n",newIntegers);
#endif
  int n=numberIntegers_;
  // Now rest of old
  for (i=0;i<numberObjects_;i++) { 
    if (object_[i]) {
      CbcSimpleInteger * obj =
        dynamic_cast <CbcSimpleInteger *>(object_[i]) ;
      if (obj) {
        delete object_[i];
      } else {
        temp[n++]=object_[i];
      }
    }
  }
  // and rest of new
  for (i=0;i<numberObjects;i++) { 
    CbcSimpleInteger * obj =
      dynamic_cast <CbcSimpleInteger *>(objects[i]) ;
    if (!obj) {
      temp[n]=objects[i]->clone();
      CbcObject * obj =
	dynamic_cast <CbcObject *>(temp[n]) ;
      if (obj)
	obj->setModel(this);
      n++;
    }
  }
  delete [] mark;
  assert (ownObjects_);
  delete [] object_;
  object_ = temp;
  assert (n==newNumberObjects);
  numberObjects_ = newNumberObjects;
}
/* Add in any object information (objects are cloned - owner can delete
   originals */
void 
CbcModel::addObjects(int numberObjects, OsiObject ** objects)
{
  // If integers but not enough objects fudge
  if (numberIntegers_>numberObjects_)
    findIntegers(true);
  /* But if incoming objects inherit from simple integer we just want
     to replace */
  int numberColumns = solver_->getNumCols();
  /** mark is -1 if not integer, >=0 if using existing simple integer and
      >=numberColumns if using new integer */
  int * mark = new int[numberColumns];
  int i;
  for (i=0;i<numberColumns;i++)
    mark[i]=-1;
  int newNumberObjects = numberObjects;
  int newIntegers=0;
  for (i=0;i<numberObjects;i++) { 
    CbcSimpleInteger * obj =
      dynamic_cast <CbcSimpleInteger *>(objects[i]) ;
    if (obj) {
      int iColumn = obj->columnNumber();
      mark[iColumn]=i+numberColumns;
      newIntegers++;
    } else {
      OsiSimpleInteger * obj2 =
      dynamic_cast <OsiSimpleInteger *>(objects[i]) ;
      if (obj2) {
	// Osi takes precedence
	int iColumn = obj2->columnNumber();
	mark[iColumn]=i+numberColumns;
	newIntegers++;
      }
    }
  }
  // and existing
  for (i=0;i<numberObjects_;i++) { 
    CbcSimpleInteger * obj =
      dynamic_cast <CbcSimpleInteger *>(object_[i]) ;
    if (obj) {
      int iColumn = obj->columnNumber();
      if (mark[iColumn]<0) {
        newIntegers++;
        newNumberObjects++;
        mark[iColumn]=i;
      }
    }
  } 
  delete [] integerVariable_;
  integerVariable_=NULL;
#ifdef COIN_DEVELOP
  if (newIntegers!=numberIntegers_) 
    printf("changing number of integers from %d to %d\n",
           numberIntegers_,newIntegers);
#endif
  numberIntegers_ = newIntegers;
  integerVariable_ = new int [numberIntegers_];
  OsiObject ** temp  = new OsiObject * [newNumberObjects];
  // Put integers first
  newIntegers=0;
  numberIntegers_=0;
  for (i=0;i<numberColumns;i++) {
    int which = mark[i];
    if (which>=0) {
      if (!isInteger(i)) {
        newIntegers++;
        solver_->setInteger(i);
      }
      if (which<numberColumns) {
        temp[numberIntegers_]=object_[which];
        object_[which]=NULL;
      } else {
        temp[numberIntegers_]=objects[which-numberColumns]->clone();
      }
      integerVariable_[numberIntegers_++]=i;
    }
  }
#ifdef COIN_DEVELOP
  if (newIntegers)
    printf("%d variables were declared integer\n",newIntegers);
#endif
  int n=numberIntegers_;
  // Now rest of old
  for (i=0;i<numberObjects_;i++) { 
    if (object_[i]) {
      CbcSimpleInteger * obj =
        dynamic_cast <CbcSimpleInteger *>(object_[i]) ;
      if (obj) {
        delete object_[i];
      } else {
	temp[n++]=object_[i];
      }
    }
  }
  // and rest of new
  for (i=0;i<numberObjects;i++) { 
    CbcSimpleInteger * obj =
      dynamic_cast <CbcSimpleInteger *>(objects[i]) ;
    OsiSimpleInteger * obj2 =
      dynamic_cast <OsiSimpleInteger *>(objects[i]) ;
    if (!obj&&!obj2) {
      temp[n]=objects[i]->clone();
      CbcObject * obj =
	dynamic_cast <CbcObject *>(temp[n]) ;
      if (obj)
	obj->setModel(this);
      n++;
    }
  }
  delete [] mark;
  assert (ownObjects_);
  delete [] object_;
  object_ = temp;
  assert (n==newNumberObjects);
  numberObjects_ = newNumberObjects;
}

/**
  This routine sets the objective cutoff value used for fathoming and
  determining monotonic variables.

  If the fathoming discipline is strict, a small tolerance is added to the
  new cutoff. This avoids problems due to roundoff when the target value
  is exact. The common example would be an IP with only integer variables in
  the objective. If the target is set to the exact value z of the optimum,
  it's possible to end up fathoming an ancestor of the solution because the
  solver returns z+epsilon.

  Determining if strict fathoming is needed is best done by analysis.
  In cbc, that's analyseObjective. The default is false.

  In cbc we always minimize so add epsilon
*/

void CbcModel::setCutoff (double value)

{
#if 0
  double tol = 0 ;
  int fathomStrict = getIntParam(CbcFathomDiscipline) ;
  if (fathomStrict == 1)
  { solver_->getDblParam(OsiDualTolerance,tol) ;
  tol = tol*(1+fabs(value)) ;
  
  value += tol ; }
#endif
  dblParam_[CbcCurrentCutoff]=value;
  if (solver_) {
    // Solvers know about direction
    double direction = solver_->getObjSense();
    solver_->setDblParam(OsiDualObjectiveLimit,value*direction);
  }
}



/*
  Call this to really test if a valid solution can be feasible. The cutoff is
  passed in as a parameter so that we don't need to worry here after swapping
  solvers.  The solution is assumed to be numberColumns in size.  If
  fixVariables is true then the bounds of the continuous solver are updated.
  The routine returns the objective value determined by reoptimizing from
  scratch. If the solution is rejected, this will be worse than the cutoff.

  TODO: There's an issue with getting the correct cutoff value: We update the
	cutoff in the regular solver, but not in continuousSolver_. But our only
	use for continuousSolver_ is verifying candidate solutions. Would it
	make sense to update the cutoff? Then we wouldn't need to step around
	isDualObjectiveLimitReached().
*/
double 
CbcModel::checkSolution (double cutoff, const double *solution,
			 bool fixVariables, double objectiveValue)

{
  if (!solverCharacteristics_->solutionAddsCuts()) {
    // Can trust solution
    int numberColumns = solver_->getNumCols();
    
    /*
      Grab the continuous solver (the pristine copy of the problem, made before
      starting to work on the root node). Save the bounds on the variables.
      Install the solution passed as a parameter, and copy it to the model's
      currentSolution_.
      
      TODO: This is a belt-and-suspenders approach. Once the code has settled
      a bit, we can cast a critical eye here.
    */
    OsiSolverInterface * saveSolver = solver_;
    if (continuousSolver_)
      solver_ = continuousSolver_;
    // move solution to continuous copy
    solver_->setColSolution(solution);
    // Put current solution in safe place
    // Point to current solution
    const double * save = testSolution_;
    // Safe as will be const inside infeasibility()
    testSolution_ = solver_->getColSolution();
    //memcpy(currentSolution_,solver_->getColSolution(),
    // numberColumns*sizeof(double));
    //solver_->messageHandler()->setLogLevel(4);
    
    // save original bounds
    double * saveUpper = new double[numberColumns];
    double * saveLower = new double[numberColumns];
    memcpy(saveUpper,getColUpper(),numberColumns*sizeof(double));
    memcpy(saveLower,getColLower(),numberColumns*sizeof(double));
    // point to useful information
    OsiBranchingInformation usefulInfo=usefulInformation();
    
    /*
      Run through the objects and use feasibleRegion() to set variable bounds
      so as to fix the variables specified in the objects at their value in this
      solution. Since the object list contains (at least) one object for every
      integer variable, this has the effect of fixing all integer variables.
    */
    int i;
    for (i=0;i<numberObjects_;i++)
      object_[i]->feasibleRegion(solver_,&usefulInfo);
    // We can switch off check
    if ((specialOptions_&4)==0) {
      if ((specialOptions_&2)==0&&solverCharacteristics_->warmStart()) {
        /*
          Remove any existing warm start information to be sure there is no
          residual influence on initialSolve().
        */
        CoinWarmStartBasis *slack =
          dynamic_cast<CoinWarmStartBasis *>(solver_->getEmptyWarmStart()) ;
        solver_->setWarmStart(slack);
        delete slack ;
      }
      // Give a hint to do dual
      bool saveTakeHint;
      OsiHintStrength saveStrength;
#ifndef NDEBUG
      bool gotHint = (solver_->getHintParam(OsiDoDualInInitial,saveTakeHint,saveStrength));
      assert (gotHint);
#else
      (solver_->getHintParam(OsiDoDualInInitial,saveTakeHint,saveStrength));
#endif
      solver_->setHintParam(OsiDoDualInInitial,true,OsiHintTry);
      solver_->initialSolve();
      if (!solver_->isProvenOptimal())
        { 
#ifdef COIN_DEVELOP
	  printf("checkSolution infeas! Retrying with primal.\n");
#endif
          //bool saveTakeHint;
          //OsiHintStrength saveStrength;
          //bool savePrintHint;
          //solver_->writeMps("infeas");
          //bool gotHint = (solver_->getHintParam(OsiDoReducePrint,savePrintHint,saveStrength));
          //gotHint = (solver_->getHintParam(OsiDoScale,saveTakeHint,saveStrength));
          //solver_->setHintParam(OsiDoScale,false,OsiHintTry);
          //solver_->setHintParam(OsiDoReducePrint,false,OsiHintTry) ;
          solver_->setHintParam(OsiDoDualInInitial,false,OsiHintTry);
          solver_->initialSolve();
          //solver_->setHintParam(OsiDoScale,saveTakeHint,saveStrength);
          //solver_->setHintParam(OsiDoReducePrint,savePrintHint,OsiHintTry) ;
	  // go from all slack now
	  specialOptions_ &= ~2;
	  if (!solver_->isProvenOptimal()) { 
	    CoinWarmStartBasis *slack =
	      dynamic_cast<CoinWarmStartBasis *>(solver_->getEmptyWarmStart()) ;
	    solver_->setWarmStart(slack);
	    delete slack ;
#ifdef COIN_DEVELOP
	    printf("checkSolution infeas! Retrying wihout basis and with primal.\n");
#endif
	    solver_->initialSolve();
#ifdef COIN_DEVELOP
	    if (!solver_->isProvenOptimal()) { 
	      printf("checkSolution still infeas!\n");
	    }
#endif
	  }
        }
      //assert(solver_->isProvenOptimal());
      solver_->setHintParam(OsiDoDualInInitial,saveTakeHint,saveStrength);
      objectiveValue = solver_->getObjValue()*solver_->getObjSense();
    }
    
    /*
      Check that the solution still beats the objective cutoff.
      
      If it passes, make a copy of the primal variable values and do some
      cleanup and checks:
      + Values of all variables are are within original bounds and values of
      all integer variables are within tolerance of integral.
      + There are no constraint violations.
      There really should be no need for the check against original bounds.
      Perhaps an opportunity for a sanity check?
    */
    if (objectiveValue>cutoff&&objectiveValue<cutoff+1.0e-8+1.0e-8*fabs(cutoff))
      cutoff = objectiveValue; // relax
    if ((solver_->isProvenOptimal()||(specialOptions_&4)!=0) && objectiveValue <= cutoff) { 
      double * solution = new double[numberColumns];
      memcpy(solution ,solver_->getColSolution(),numberColumns*sizeof(double)) ;
      
      int iColumn;
#ifndef NDEBUG
      double integerTolerance = getIntegerTolerance() ;
#endif
      for (iColumn = 0 ; iColumn < numberColumns ; iColumn++) {
        double value = solution[iColumn] ;
        value = CoinMax(value, saveLower[iColumn]) ;
        value = CoinMin(value, saveUpper[iColumn]) ;
        if (solver_->isInteger(iColumn)) 
          assert(fabs(value-solution[iColumn]) <= integerTolerance) ;
        solution[iColumn] = value ; 
      }
      if ((specialOptions_&16)==0) {
        const double * rowLower = solver_->getRowLower() ;
        const double * rowUpper = solver_->getRowUpper() ;
        int numberRows = solver_->getNumRows() ;
        double *rowActivity = new double[numberRows] ;
        memset(rowActivity,0,numberRows*sizeof(double)) ;
        solver_->getMatrixByCol()->times(solution,rowActivity) ;
        double primalTolerance ;
        solver_->getDblParam(OsiPrimalTolerance,primalTolerance) ;
        double largestInfeasibility =0.0;
        for (i=0 ; i < numberRows ; i++) {
          largestInfeasibility = CoinMax(largestInfeasibility,
                                         rowLower[i]-rowActivity[i]);
          largestInfeasibility = CoinMax(largestInfeasibility,
                                         rowActivity[i]-rowUpper[i]);
        }
        if (largestInfeasibility>100.0*primalTolerance) {
          handler_->message(CBC_NOTFEAS3, messages_)
            << largestInfeasibility << CoinMessageEol ;
          objectiveValue=1.0e50 ; 
        }
        delete [] rowActivity ;
      }
      delete [] solution;
    } else {
      objectiveValue=1.0e50 ; 
    }
    /*
      Regardless of what we think of the solution, we may need to restore the
      original bounds of the continuous solver. Unfortunately, const'ness
      prevents us from simply reversing the memcpy used to make these snapshots.
    */
    if (!fixVariables)
      { for (int iColumn = 0 ; iColumn < numberColumns ; iColumn++)
        { solver_->setColLower(iColumn,saveLower[iColumn]) ;
        solver_->setColUpper(iColumn,saveUpper[iColumn]) ; } }
    delete [] saveLower;
    delete [] saveUpper;
    
    /*
      Restore the usual solver.
    */
    solver_=saveSolver;
    testSolution_ = save;
    return objectiveValue;
  } else {
    // Outer approximation or similar
    //If this is true then the solution comes from the nlp we don't need to resolve the same nlp with ipopt
    //solverCharacteristics_->setSolver(solver_);
    bool solutionComesFromNlp = solverCharacteristics_->bestObjectiveValue()<cutoff;
    double objectiveValue;
    int numberColumns = solver_->getNumCols();
    double *saveLower = NULL;
    double * saveUpper = NULL;
    
    if(! solutionComesFromNlp)//Otherwise solution already comes from ipopt and cuts are known
      {
        if(fixVariables)//Will temporarily fix all integer valued var
          {
            // save original bounds
            saveUpper = new double[numberColumns];
            saveLower = new double[numberColumns];
            memcpy(saveUpper,solver_->getColUpper(),numberColumns*sizeof(double));
            memcpy(saveLower,solver_->getColLower(),numberColumns*sizeof(double));
            //in any case solution should be already loaded into solver_
            /*
              Run through the objects and use feasibleRegion() to set variable bounds
              so as to fix the variables specified in the objects at their value in this
              solution. Since the object list contains (at least) one object for every
              integer variable, this has the effect of fixing all integer variables.
            */
            const double * save = testSolution_;
            testSolution_ = solution;
	    // point to useful information
	    OsiBranchingInformation usefulInfo=usefulInformation();
            for (int i=0;i<numberObjects_;i++)
              object_[i]->feasibleRegion(solver_,&usefulInfo);
            testSolution_ = save;
            resolve(solver_);
          }
        
        /*
          Now step through the cut generators and see if any of them are flagged to
          run when a new solution is discovered. Only global cuts are useful. 
          (The solution being evaluated may not correspond to the current location in the
          search tree --- discovered by heuristic, for example.)
        */
        OsiCuts theseCuts;
        int i;
        int lastNumberCuts=0;
	// reset probing info
	if (probingInfo_)
	  probingInfo_->initializeFixing();
        for (i=0;i<numberCutGenerators_;i++) 
          {
            if (generator_[i]->atSolution()) 
              {
                generator_[i]->generateCuts(theseCuts,true,solver_,NULL);
                int numberCuts = theseCuts.sizeRowCuts();
                for (int j=lastNumberCuts;j<numberCuts;j++) 
                  {
                    const OsiRowCut * thisCut = theseCuts.rowCutPtr(j);
                    if (thisCut->globallyValid()) 
                      {
                        //           if ((specialOptions_&1)!=0) 
                        //           {
                        //             /* As these are global cuts -
                        //             a) Always get debugger object
                        // b) Not fatal error to cutoff optimal (if we have just got optimal)
                        // */
                        //             const OsiRowCutDebugger *debugger = solver_->getRowCutDebuggerAlways() ;
                        //             if (debugger) 
                        //             {
                        //               if(debugger->invalidCut(*thisCut))
                        //                 printf("ZZZZ Global cut - cuts off optimal solution!\n");
                        //             }
                        //           }
                        // add to global list
			OsiRowCut newCut(*thisCut);
			newCut.setGloballyValid(2);
			newCut.mutableRow().setTestForDuplicateIndex(false);
			globalCuts_.insert(newCut) ;
                        generator_[i]->incrementNumberCutsInTotal();
                      }
                  }
              }
          }
        //   int numberCuts = theseCuts.sizeColCuts();
        //   for (i=0;i<numberCuts;i++) {
        //     const OsiColCut * thisCut = theseCuts.colCutPtr(i);
        //     if (thisCut->globallyValid()) {
        //       // add to global list
        //       globalCuts_.insert(*thisCut);
        //     }
        //   }
        //have to retrieve the solution and its value from the nlp
      }
    double newObjectiveValue=cutoff;
    if(solverCharacteristics_->solution(newObjectiveValue,
                                        const_cast<double *> (solution),
                                        numberColumns))
      {
        objectiveValue = newObjectiveValue;
      }
    else
      {
        objectiveValue = 2e50;
      }
    if (!solutionComesFromNlp && fixVariables)
      { 
        for (int iColumn = 0 ; iColumn < numberColumns ; iColumn++)
          { 
            solver_->setColLower(iColumn,saveLower[iColumn]) ;
            solver_->setColUpper(iColumn,saveUpper[iColumn]) ;  
          }
        delete [] saveLower;
        delete [] saveUpper;
      }
    //If the variables were fixed the cutting plane procedure may have believed that the node could be fathomed
    //re-establish truth.- should do no harm for non nlp
    if(!solutionComesFromNlp && fixVariables)
      solverCharacteristics_->setMipBound(-COIN_DBL_MAX);
    return objectiveValue;
  }
}

/*
  Call this routine from anywhere when a solution is found. The solution
  vector is assumed to contain one value for each structural variable.

  The first action is to run checkSolution() to confirm the objective and
  feasibility. If this check causes the solution to be rejected, we're done.
  If fixVariables = true, the variable bounds held by the continuous solver
  will be left fixed to the values in the solution; otherwise they are
  restored to the original values.

  If the solution is accepted, install it as the best solution.

  The routine also contains a hook to run any cut generators that are flagged
  to run when a new solution is discovered. There's a potential hazard because
  the cut generators see the continuous solver >after< possible restoration of
  original bounds (which may well invalidate the solution).
*/

void
CbcModel::setBestSolution (CBC_Message how,
			   double & objectiveValue, const double *solution,
			   bool fixVariables)

{
  if (!solverCharacteristics_->solutionAddsCuts()) {
    // Can trust solution
    double cutoff = getCutoff();
    if (cutoff<1.0e30)
      cutoff = CoinMin(cutoff,bestObjective_) ;
    
    /*
      Double check the solution to catch pretenders.
    */
    double saveObjectiveValue = objectiveValue;
    objectiveValue =checkSolution(cutoff,solution,fixVariables,objectiveValue);
    if (objectiveValue>cutoff&&objectiveValue<cutoff+1.0e-8+1.0e-8*fabs(cutoff))
      cutoff = objectiveValue; // relax
    if (objectiveValue > cutoff||objectiveValue>1.0e30) {
      if (objectiveValue>1.0e30)
        handler_->message(CBC_NOTFEAS1, messages_) << CoinMessageEol ;
      else
        handler_->message(CBC_NOTFEAS2, messages_)
          << objectiveValue << cutoff << CoinMessageEol ; 
    } else if (objectiveValue<bestObjective_) {
      /*
        We have a winner. Install it as the new incumbent.
        Bump the objective cutoff value and solution counts. Give the user the
        good news.
      */
      bestObjective_ = objectiveValue;
      int numberColumns = solver_->getNumCols();
      if (!bestSolution_)
        bestSolution_ = new double[numberColumns];
      CoinCopyN(solution,numberColumns,bestSolution_);
      
      cutoff = bestObjective_-dblParam_[CbcCutoffIncrement];
      // But allow for rounding errors
      if (dblParam_[CbcCutoffIncrement] == 1e-5) {
#ifdef COIN_DEVELOP
	if (saveObjectiveValue+1.0e-7<bestObjective_)
	  printf("First try at solution had objective %.16g, rechecked as %.16g\n",
		 saveObjectiveValue,bestObjective_);
#endif
	saveObjectiveValue = CoinMax(saveObjectiveValue,bestObjective_-0.0000001*fabs(bestObjective_));
	cutoff = CoinMin(bestObjective_,saveObjectiveValue)-1.0e-5;
      }
      // This is not correct - that way cutoff can go up if maximization
      //double direction = solver_->getObjSense();
      //setCutoff(cutoff*direction);
      setCutoff(cutoff);
      
      if (how==CBC_ROUNDING)
        numberHeuristicSolutions_++;
      numberSolutions_++;
      if (numberHeuristicSolutions_==numberSolutions_) 
        stateOfSearch_ = 1;
      else 
        stateOfSearch_ = 2;

      if (how!=CBC_ROUNDING) {
	handler_->message(how,messages_)
	  <<bestObjective_<<numberIterations_
	  <<numberNodes_<<getCurrentSeconds()
	  <<CoinMessageEol;
      } else {
	assert (lastHeuristic_);
	const char * name = lastHeuristic_->heuristicName();
	handler_->message(CBC_ROUNDING,messages_)
	  <<bestObjective_
	  <<name
	  <<numberIterations_
	  <<numberNodes_<<getCurrentSeconds()
	  <<CoinMessageEol;
      }
      /*
        Now step through the cut generators and see if any of them are flagged to
        run when a new solution is discovered. Only global cuts are useful. (The
        solution being evaluated may not correspond to the current location in the
        search tree --- discovered by heuristic, for example.)
      */
      OsiCuts theseCuts;
      int i;
      int lastNumberCuts=0;
      // reset probing info
      if (probingInfo_)
	probingInfo_->initializeFixing();
      for (i=0;i<numberCutGenerators_;i++) {
        bool generate = generator_[i]->atSolution();
        // skip if not optimal and should be (maybe a cut generator has fixed variables)
        if (generator_[i]->needsOptimalBasis()&&!solver_->basisIsAvailable())
          generate=false;
        if (generate) {
          generator_[i]->generateCuts(theseCuts,true,solver_,NULL);
          int numberCuts = theseCuts.sizeRowCuts();
          for (int j=lastNumberCuts;j<numberCuts;j++) {
            const OsiRowCut * thisCut = theseCuts.rowCutPtr(j);
            if (thisCut->globallyValid()) {
              if ((specialOptions_&1)!=0) {
                /* As these are global cuts -
                   a) Always get debugger object
                   b) Not fatal error to cutoff optimal (if we have just got optimal)
                */
                const OsiRowCutDebugger *debugger = solver_->getRowCutDebuggerAlways() ;
                if (debugger) {
                  if(debugger->invalidCut(*thisCut))
                    printf("ZZZZ Global cut - cuts off optimal solution!\n");
                }
              }
              // add to global list
	      OsiRowCut newCut(*thisCut);
	      newCut.setGloballyValid(2);
	      newCut.mutableRow().setTestForDuplicateIndex(false);
	      globalCuts_.insert(newCut) ;
              generator_[i]->incrementNumberCutsInTotal();
            }
          }
        }
      }
      int numberCuts = theseCuts.sizeColCuts();
      for (i=0;i<numberCuts;i++) {
        const OsiColCut * thisCut = theseCuts.colCutPtr(i);
        if (thisCut->globallyValid()) {
          // add to global list
	  OsiColCut newCut(*thisCut);
	  newCut.setGloballyValid(2);
	  globalCuts_.insert(newCut) ;
        }
      }
    }
  } else {
    // Outer approximation or similar
    double cutoff = getCutoff() ;
    
    /*
      Double check the solution to catch pretenders.
    */
    
    int numberRowBefore = solver_->getNumRows();
    int numberColBefore = solver_->getNumCols();
    double *saveColSol=NULL;
    
    CoinWarmStart * saveWs=NULL;
    // if(how!=CBC_SOLUTION) return;
    if(how==CBC_ROUNDING)//We don't want to make any change to solver_
      //take a snapshot of current state
      {     
        //save solution
        saveColSol = new double[numberColBefore];
        CoinCopyN(solver_->getColSolution(), numberColBefore, saveColSol);
        //save warm start
        saveWs = solver_->getWarmStart();     
      }
    
    //run check solution this will eventually generate cuts
    //if in strongBranching or heuristic will do only one cut generation iteration
    // by fixing variables.
    fixVariables = fixVariables || (how==CBC_ROUNDING) || (how==CBC_STRONGSOL);
    double * candidate = new double[numberColBefore];
    CoinCopyN(solution, numberColBefore, candidate);
    objectiveValue = checkSolution(cutoff,candidate,fixVariables,objectiveValue);
    
    //If it was an heuristic solution we have to clean up the solver
    if (how==CBC_ROUNDING)
      {
        //delete the cuts
        int currentNumberRowCuts = solver_->getNumRows() - numberRowBefore;
        int currentNumberColCuts = solver_->getNumCols() - numberColBefore;
        if(CoinMax(currentNumberColCuts, currentNumberRowCuts)>0)
          {
            int *which = new int[CoinMax(currentNumberColCuts, currentNumberRowCuts)];
            if(currentNumberRowCuts)
              {
                for (int i = 0 ; i < currentNumberRowCuts ; i++)
                  which[i] = i + numberRowBefore;
                
                solver_->deleteRows(currentNumberRowCuts,which);
              }
            if(currentNumberColCuts)
              {
                for (int i = 0 ; i < currentNumberColCuts ; i++)
                  which[i] = i+numberColBefore;
                solver_->deleteCols(currentNumberColCuts,which);
              }
            delete [] which;
          }
        // Reset solution and warm start info
        solver_->setColSolution(saveColSol);
        solver_->setWarmStart(saveWs);
        delete [] saveColSol;
        delete saveWs;
      }
    
    if (objectiveValue > cutoff) {
      if (!solverCharacteristics_->solutionAddsCuts()) {
        if (objectiveValue>1.0e30)
          handler_->message(CBC_NOTFEAS1, messages_) << CoinMessageEol ;
        else
          handler_->message(CBC_NOTFEAS2, messages_)
            << objectiveValue << cutoff << CoinMessageEol ;
      }
    } else {
      /*
        We have a winner. Install it as the new incumbent.
        Bump the objective cutoff value and solution counts. Give the user the
        good news.
        NB - Not all of this if from solve with cuts
      */
      bestObjective_ = objectiveValue;
      int numberColumns = solver_->getNumCols();
      if (!bestSolution_)
        bestSolution_ = new double[numberColumns];
      CoinCopyN(candidate,numberColumns,bestSolution_);

      // don't update if from solveWithCuts
      if (how!=CBC_SOLUTION2) {
        if (how==CBC_ROUNDING)
          numberHeuristicSolutions_++;
        cutoff = bestObjective_-dblParam_[CbcCutoffIncrement];
        // This is not correct - that way cutoff can go up if maximization
        //double direction = solver_->getObjSense();
        //setCutoff(cutoff*direction);
        setCutoff(cutoff);
      
        numberSolutions_++;
        if (numberNodes_== 0 || numberHeuristicSolutions_==numberSolutions_) 
          stateOfSearch_ = 1;
        else 
          stateOfSearch_ = 2;
        
        handler_->message(how,messages_)
          <<bestObjective_<<numberIterations_
          <<numberNodes_<<getCurrentSeconds()
          <<CoinMessageEol;
      }
    }
    delete [] candidate;
  }
  return ;
}


/* Test the current solution for feasibility.

   Calculate the number of standard integer infeasibilities, then scan the
   remaining objects to see if any of them report infeasibilities.

   Currently (2003.08) the only object besides SimpleInteger is Clique, hence
   the comments about `odd ones' infeasibilities.
*/
bool 
CbcModel::feasibleSolution(int & numberIntegerInfeasibilities,
			int & numberObjectInfeasibilities) const
{
  int numberUnsatisfied=0;
  double sumUnsatisfied=0.0;
  int preferredWay;
  int j;
  // Point to current solution
  const double * save = testSolution_;
  // Safe as will be const inside infeasibility()
  testSolution_ = solver_->getColSolution();
  // Put current solution in safe place
  //memcpy(currentSolution_,solver_->getColSolution(),
  // solver_->getNumCols()*sizeof(double));
  // point to useful information
  OsiBranchingInformation usefulInfo=usefulInformation();
  for (j=0;j<numberIntegers_;j++) {
    const OsiObject * object = object_[j];
    double infeasibility = object->infeasibility(&usefulInfo,preferredWay);
    if (infeasibility) {
      assert (infeasibility>0);
      numberUnsatisfied++;
      sumUnsatisfied += infeasibility;
    }
  }
  numberIntegerInfeasibilities = numberUnsatisfied;
  for (;j<numberObjects_;j++) {
    const OsiObject * object = object_[j];
    double infeasibility = object->infeasibility(&usefulInfo,preferredWay);
    if (infeasibility) {
      assert (infeasibility>0);
      numberUnsatisfied++;
      sumUnsatisfied += infeasibility;
    }
  }
  // and restore
  testSolution_ = save;
  numberObjectInfeasibilities = numberUnsatisfied-numberIntegerInfeasibilities;
  return (!numberUnsatisfied);
}

/* For all vubs see if we can tighten bounds by solving Lp's
   type - 0 just vubs
   1 all (could be very slow)
   -1 just vubs where variable away from bound
   Returns false if not feasible
*/
bool 
CbcModel::tightenVubs(int type, bool allowMultipleBinary, double useCutoff)
{

  CoinPackedMatrix matrixByRow(*solver_->getMatrixByRow());
  int numberRows = solver_->getNumRows();
  int numberColumns = solver_->getNumCols();

  int iRow,iColumn;

  // Row copy
  //const double * elementByRow = matrixByRow.getElements();
  const int * column = matrixByRow.getIndices();
  const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
  const int * rowLength = matrixByRow.getVectorLengths();

  const double * colUpper = solver_->getColUpper();
  const double * colLower = solver_->getColLower();
  //const double * rowUpper = solver_->getRowUpper();
  //const double * rowLower = solver_->getRowLower();

  const double * objective = solver_->getObjCoefficients();
  //double direction = solver_->getObjSense();
  const double * colsol = solver_->getColSolution();

  int numberVub=0;
  int * continuous = new int[numberColumns];
  if (type >= 0) {
    double * sort = new double[numberColumns];
    for (iRow=0;iRow<numberRows;iRow++) {
      int j;
      int numberBinary=0;
      int numberUnsatisfiedBinary=0;
      int numberContinuous=0;
      int iCont=-1;
      double weight=1.0e30;
      for (j=rowStart[iRow];j<rowStart[iRow]+rowLength[iRow];j++) {
	int iColumn = column[j];
	if (colUpper[iColumn]-colLower[iColumn]>1.0e-8) {
	  if (solver_->isFreeBinary(iColumn)) {
	    numberBinary++;
	    /* For sort I make naive assumption:
	       x - a * delta <=0 or
	       -x + a * delta >= 0
	    */
	    if (colsol[iColumn]>colLower[iColumn]+1.0e-6&&
		colsol[iColumn]<colUpper[iColumn]-1.0e-6) {
	      numberUnsatisfiedBinary++;
	      weight = CoinMin(weight,fabs(objective[iColumn]));
	    }
	  } else {
	    numberContinuous++;
	    iCont=iColumn;
	  }
	}
      }
      if (numberContinuous==1&&numberBinary) {
	if (numberBinary==1||allowMultipleBinary) {
	  // treat as vub
	  if (!numberUnsatisfiedBinary)
	    weight=-1.0; // at end
	  sort[numberVub]=-weight;
	  continuous[numberVub++] = iCont;
	}
      }
    }
    if (type>0) {
      // take so many
      CoinSort_2(sort,sort+numberVub,continuous);
      numberVub = CoinMin(numberVub,type);
    }
    delete [] sort;
  } else {
    for (iColumn=0;iColumn<numberColumns;iColumn++) 
      continuous[iColumn]=iColumn;
    numberVub=numberColumns;
  }
  bool feasible = tightenVubs(numberVub,continuous,useCutoff);
  delete [] continuous;

  return feasible;
}
// This version is just handed a list of variables
bool 
CbcModel::tightenVubs(int numberSolves, const int * which,
		      double useCutoff)
{

  int numberColumns = solver_->getNumCols();

  int iColumn;
  
  OsiSolverInterface * solver = solver_;
  double saveCutoff = getCutoff() ;
  
  double * objective = new double[numberColumns];
  memcpy(objective,solver_->getObjCoefficients(),numberColumns*sizeof(double));
  double direction = solver_->getObjSense();
  
  // add in objective if there is a cutoff
  if (useCutoff<1.0e30) {
    // get new version of model
    solver = solver_->clone();
    CoinPackedVector newRow;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      solver->setObjCoeff(iColumn,0.0); // zero out in new model
      if (objective[iColumn]) 
	newRow.insert(iColumn,direction * objective[iColumn]);
      
    }
    solver->addRow(newRow,-COIN_DBL_MAX,useCutoff);
    // signal no objective
    delete [] objective;
    objective=NULL;
  }
  setCutoff(COIN_DBL_MAX);


  bool * vub = new bool [numberColumns];
  int iVub;

  // mark vub columns
  for (iColumn=0;iColumn<numberColumns;iColumn++) 
    vub[iColumn]=false;
  for (iVub=0;iVub<numberSolves;iVub++) 
    vub[which[iVub]]=true;
  OsiCuts cuts;
  // First tighten bounds anyway if CglProbing there
  CglProbing* generator = NULL;
  int iGen;
  // reset probing info
  if (probingInfo_)
    probingInfo_->initializeFixing();
  for (iGen=0;iGen<numberCutGenerators_;iGen++) {
    generator = dynamic_cast<CglProbing*>(generator_[iGen]->generator());
    if (generator)
      break;
  }
  int numberFixed=0;
  int numberTightened=0;
  int numberFixedByProbing=0;
  int numberTightenedByProbing=0;
  int printFrequency = (numberSolves+19)/20; // up to 20 messages
  int save[4]={0,0,0,0};
  if (generator) {
    // set to cheaper and then restore at end
    save[0]=generator->getMaxPass();
    save[1]=generator->getMaxProbe();
    save[2]=generator->getMaxLook();
    save[3]=generator->rowCuts();
    generator->setMaxPass(1);
    generator->setMaxProbe(10);
    generator->setMaxLook(50);
    generator->setRowCuts(0);
    
    // Probing - return tight column bounds
    CglTreeInfo info;
    generator->generateCutsAndModify(*solver,cuts,&info);
    const double * tightLower = generator->tightLower();
    const double * lower = solver->getColLower();
    const double * tightUpper = generator->tightUpper();
    const double * upper = solver->getColUpper();
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      double newUpper = tightUpper[iColumn];
      double newLower = tightLower[iColumn];
      if (newUpper<upper[iColumn]-1.0e-8*(fabs(upper[iColumn])+1)||
	  newLower>lower[iColumn]+1.0e-8*(fabs(lower[iColumn])+1)) {
	if (newUpper<newLower) {
	  fprintf(stderr,"Problem is infeasible\n");
	  return false;
	}
	if (newUpper==newLower) {
	  numberFixed++;
	  numberFixedByProbing++;
	  solver->setColLower(iColumn,newLower);
	  solver->setColUpper(iColumn,newUpper);
	  printf("Column %d, new bounds %g %g\n",iColumn,
		 newLower,newUpper);
	} else if (vub[iColumn]) {
	  numberTightened++;
	  numberTightenedByProbing++;
	  if (!solver->isInteger(iColumn)) {
	    // relax
	    newLower=CoinMax(lower[iColumn],
				    newLower
				    -1.0e-5*(fabs(lower[iColumn])+1));
	    newUpper=CoinMin(upper[iColumn],
				    newUpper
				    +1.0e-5*(fabs(upper[iColumn])+1));
	  }
	  solver->setColLower(iColumn,newLower);
	  solver->setColUpper(iColumn,newUpper);
	}
      }
    }
  }
  CoinWarmStart * ws = solver->getWarmStart();
  double * solution = new double [numberColumns];
  memcpy(solution,solver->getColSolution(),numberColumns*sizeof(double));
  for (iColumn=0;iColumn<numberColumns;iColumn++) 
    solver->setObjCoeff(iColumn,0.0);
  //solver->messageHandler()->setLogLevel(2);
  for (iVub=0;iVub<numberSolves;iVub++) {
    iColumn = which[iVub];
    int iTry;
    for (iTry=0;iTry<2;iTry++) {
      double saveUpper = solver->getColUpper()[iColumn];
      double saveLower = solver->getColLower()[iColumn];
      double value;
      if (iTry==1) {
	// try all way up
	solver->setObjCoeff(iColumn,-1.0);
      } else {
	// try all way down
	solver->setObjCoeff(iColumn,1.0);
      }
      solver->initialSolve();
      setPointers(continuousSolver_);
      value = solver->getColSolution()[iColumn];
      bool change=false;
      if (iTry==1) {
	if (value<saveUpper-1.0e-4) {
	  if (solver->isInteger(iColumn)) {
	    value = floor(value+0.00001);
	  } else {
	    // relax a bit
	    value=CoinMin(saveUpper,value+1.0e-5*(fabs(saveUpper)+1));
	  }
	  if (value-saveLower<1.0e-7) 
	    value = saveLower; // make sure exactly same
	  solver->setColUpper(iColumn,value);
	  saveUpper=value;
	  change=true;
	}
      } else {
	if (value>saveLower+1.0e-4) {
	  if (solver->isInteger(iColumn)) {
	    value = ceil(value-0.00001);
	  } else {
	    // relax a bit
	    value=CoinMax(saveLower,value-1.0e-5*(fabs(saveLower)+1));
	  }
	  if (saveUpper-value<1.0e-7) 
	    value = saveUpper; // make sure exactly same
	  solver->setColLower(iColumn,value);
	  saveLower=value;
	  change=true;
	}
      }
      solver->setObjCoeff(iColumn,0.0);
      if (change) {
	if (saveUpper==saveLower) 
	  numberFixed++;
	else
	  numberTightened++;
	int saveFixed=numberFixed;
	
	int jColumn;
	if (generator) {
	  // Probing - return tight column bounds
	  cuts = OsiCuts();
	  CglTreeInfo info;
	  generator->generateCutsAndModify(*solver,cuts,&info);
	  const double * tightLower = generator->tightLower();
	  const double * lower = solver->getColLower();
	  const double * tightUpper = generator->tightUpper();
	  const double * upper = solver->getColUpper();
	  for (jColumn=0;jColumn<numberColumns;jColumn++) {
	    double newUpper = tightUpper[jColumn];
	    double newLower = tightLower[jColumn];
	    if (newUpper<upper[jColumn]-1.0e-8*(fabs(upper[jColumn])+1)||
		newLower>lower[jColumn]+1.0e-8*(fabs(lower[jColumn])+1)) {
	      if (newUpper<newLower) {
		fprintf(stderr,"Problem is infeasible\n");
		return false;
	      }
	      if (newUpper==newLower) {
		numberFixed++;
		numberFixedByProbing++;
		solver->setColLower(jColumn,newLower);
		solver->setColUpper(jColumn,newUpper);
	      } else if (vub[jColumn]) {
		numberTightened++;
		numberTightenedByProbing++;
		if (!solver->isInteger(jColumn)) {
		  // relax
		  newLower=CoinMax(lower[jColumn],
			       newLower
			       -1.0e-5*(fabs(lower[jColumn])+1));
		  newUpper=CoinMin(upper[jColumn],
			       newUpper
			       +1.0e-5*(fabs(upper[jColumn])+1));
		}
		solver->setColLower(jColumn,newLower);
		solver->setColUpper(jColumn,newUpper);
	      }
	    }
	  }
	}
	if (numberFixed>saveFixed) {
	  // original solution may not be feasible
	  // go back to true costs to solve if exists
	  if (objective) {
	    for (jColumn=0;jColumn<numberColumns;jColumn++) 
	      solver->setObjCoeff(jColumn,objective[jColumn]);
	  }
	  solver->setColSolution(solution);
	  solver->setWarmStart(ws);
	  solver->resolve();
	  if (!solver->isProvenOptimal()) {
	    fprintf(stderr,"Problem is infeasible\n");
	    return false;
	  }
	  delete ws;
	  ws = solver->getWarmStart();
	  memcpy(solution,solver->getColSolution(),
		 numberColumns*sizeof(double));
	  for (jColumn=0;jColumn<numberColumns;jColumn++) 
	    solver->setObjCoeff(jColumn,0.0);
	}
      }
      solver->setColSolution(solution);
      solver->setWarmStart(ws);
    }
    if (iVub%printFrequency==0) 
      handler_->message(CBC_VUB_PASS,messages_)
	<<iVub+1<<numberFixed<<numberTightened
	<<CoinMessageEol;
  }
  handler_->message(CBC_VUB_END,messages_)
    <<numberFixed<<numberTightened
    <<CoinMessageEol;
  delete ws;
  delete [] solution;
  // go back to true costs to solve if exists
  if (objective) {
    for (iColumn=0;iColumn<numberColumns;iColumn++) 
      solver_->setObjCoeff(iColumn,objective[iColumn]);
    delete [] objective;
  }
  delete [] vub;
  if (generator) {
    /*printf("Probing fixed %d and tightened %d\n",
	   numberFixedByProbing,
	   numberTightenedByProbing);*/
    if (generator_[iGen]->howOften()==-1&&
	(numberFixedByProbing+numberTightenedByProbing)*5>
	(numberFixed+numberTightened))
      generator_[iGen]->setHowOften(1000000+1);
    generator->setMaxPass(save[0]);
    generator->setMaxProbe(save[1]);
    generator->setMaxLook(save[2]);
    generator->setRowCuts(save[3]);
  }

  if (solver!=solver_) {
    // move bounds across
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    const double * lowerOrig = solver_->getColLower();
    const double * upperOrig = solver_->getColUpper();
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      solver_->setColLower(iColumn,CoinMax(lower[iColumn],lowerOrig[iColumn]));
      solver_->setColUpper(iColumn,CoinMin(upper[iColumn],upperOrig[iColumn]));
    }
    delete solver;
  }
  setCutoff(saveCutoff);
  return true;
}
/* 
   Do Integer Presolve. Returns new model.
   I have to work out cleanest way of getting solution to
   original problem at end.  So this is very preliminary.
*/
CbcModel * 
CbcModel::integerPresolve(bool weak)
{
  status_ = 0;
  // solve LP
  //solver_->writeMps("bad");
  bool feasible = resolve(NULL,3);

  CbcModel * newModel = NULL;
  if (feasible) {

    // get a new model
    newModel = new CbcModel(*this);
    newModel->messageHandler()->setLogLevel(messageHandler()->logLevel());

    feasible = newModel->integerPresolveThisModel(solver_,weak);
  }
  if (!feasible) {
    handler_->message(CBC_INFEAS,messages_)
    <<CoinMessageEol;
    status_ = 0;
    secondaryStatus_ = 1;
    delete newModel;
    return NULL;
  } else {
    newModel->synchronizeModel(); // make sure everything that needs solver has it
    return newModel;
  }
}
/* 
   Do Integer Presolve - destroying current model
*/
bool 
CbcModel::integerPresolveThisModel(OsiSolverInterface * originalSolver,
				   bool weak)
{
  printf("DEPRECATED\n");
  status_ = 0;
  // solve LP
  bool feasible = resolve(NULL,3);

  bestObjective_=1.0e50;
  numberSolutions_=0;
  numberHeuristicSolutions_=0;
  double cutoff = getCutoff() ;
  double direction = solver_->getObjSense();
  if (cutoff < 1.0e20&&direction<0.0)
    messageHandler()->message(CBC_CUTOFF_WARNING1,
				    messages())
				      << cutoff << -cutoff << CoinMessageEol ;
  if (cutoff > bestObjective_)
    cutoff = bestObjective_ ;
  setCutoff(cutoff) ;
  int iColumn;
  int numberColumns = getNumCols();
  int originalNumberColumns = numberColumns;
  currentPassNumber_=0;
  synchronizeModel(); // make sure everything that needs solver has it
  if (!solverCharacteristics_) {
    OsiBabSolver * solverCharacteristics = dynamic_cast<OsiBabSolver *> (solver_->getAuxiliaryInfo());
    if (solverCharacteristics) {
      solverCharacteristics_ = solverCharacteristics;
    } else {
      // replace in solver
      OsiBabSolver defaultC;
      solver_->setAuxiliaryInfo(&defaultC);
      solverCharacteristics_ = dynamic_cast<OsiBabSolver *> (solver_->getAuxiliaryInfo());
    }
  }
  solverCharacteristics_->setSolver(solver_);
  // just point to solver_
  delete continuousSolver_;
  continuousSolver_ = solver_;
  // get a copy of original so we can fix bounds
  OsiSolverInterface * cleanModel = originalSolver->clone();
#ifdef CBC_DEBUG
  std::string problemName;
  cleanModel->getStrParam(OsiProbName,problemName);
  printf("Problem name - %s\n",problemName.c_str());
  cleanModel->activateRowCutDebugger(problemName.c_str());
  const OsiRowCutDebugger * debugger = cleanModel->getRowCutDebugger();
#endif

  // array which points from original columns to presolved
  int * original = new int[numberColumns];
  // arrays giving bounds - only ones found by probing 
  // rest will be found by presolve
  double * originalLower = new double[numberColumns];
  double * originalUpper = new double[numberColumns];
  {
    const double * lower = getColLower();
    const double * upper = getColUpper();
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      original[iColumn]=iColumn;
      originalLower[iColumn] = lower[iColumn];
      originalUpper[iColumn] = upper[iColumn];
    }
  }
  findIntegers(true);
  // save original integers
  int * originalIntegers = new int[numberIntegers_];
  int originalNumberIntegers = numberIntegers_;
  memcpy(originalIntegers,integerVariable_,numberIntegers_*sizeof(int));

  int todo=20;
  if (weak)
    todo=1;
  while (currentPassNumber_<todo) {
   
    currentPassNumber_++;
    numberSolutions_=0;
    // this will be set false to break out of loop with presolved problem
    bool doIntegerPresolve=(currentPassNumber_!=20);
    
    // Current number of free integer variables
    // Get increment in solutions
    {
      const double * objective = cleanModel->getObjCoefficients();
      const double * lower = cleanModel->getColLower();
      const double * upper = cleanModel->getColUpper();
      double maximumCost=0.0;
      bool possibleMultiple=true;
      int numberChanged=0;
      for (iColumn=0;iColumn<originalNumberColumns;iColumn++) {
	if (originalUpper[iColumn]>originalLower[iColumn]) {
	  if( cleanModel->isInteger(iColumn)) {
	    maximumCost = CoinMax(maximumCost,fabs(objective[iColumn]));
	  } else if (objective[iColumn]) {
	    possibleMultiple=false;
	  }
	}
	if (originalUpper[iColumn]<upper[iColumn]) {
#ifdef CBC_DEBUG
	  printf("Changing upper bound on %d from %g to %g\n",
		 iColumn,upper[iColumn],originalUpper[iColumn]);
#endif
	  cleanModel->setColUpper(iColumn,originalUpper[iColumn]);
	  numberChanged++;
	}
	if (originalLower[iColumn]>lower[iColumn]) {
#ifdef CBC_DEBUG
	  printf("Changing lower bound on %d from %g to %g\n",
		 iColumn,lower[iColumn],originalLower[iColumn]);
#endif
	  cleanModel->setColLower(iColumn,originalLower[iColumn]);
	  numberChanged++;
	}
      }
      // if first pass - always try
      if (currentPassNumber_==1)
	numberChanged += 1;
      if (possibleMultiple&&maximumCost) {
	int increment=0; 
	double multiplier = 2520.0;
	while (10.0*multiplier*maximumCost<1.0e8)
	  multiplier *= 10.0;
	for (int j =0;j<originalNumberIntegers;j++) {
          iColumn = originalIntegers[j];
	  if (originalUpper[iColumn]>originalLower[iColumn]) {
	    if(objective[iColumn]) {
	      double value = fabs(objective[iColumn])*multiplier;
	      int nearest = (int) floor(value+0.5);
	      if (fabs(value-floor(value+0.5))>1.0e-8||value>2.1e9) {
		increment=0;
		break; // no good
	      } else if (!increment) {
		// first
		increment=nearest;
	      } else {
		increment = gcd(increment,nearest);
	      }
	    }
	  }
	}
	if (increment) {
	  double value = increment;
	  value /= multiplier;
	  if (value*0.999>dblParam_[CbcCutoffIncrement]) {
	    messageHandler()->message(CBC_INTEGERINCREMENT,messages())
	      <<value
	      <<CoinMessageEol;
	    dblParam_[CbcCutoffIncrement]=value*0.999;
	  }
	}
      }
      if (!numberChanged) {
	doIntegerPresolve=false; // not doing any better
      }
    }
#ifdef CBC_DEBUG
    if (debugger) 
      assert(debugger->onOptimalPath(*cleanModel));
#endif
#ifdef COIN_HAS_CLP
    // do presolve - for now just clp but easy to get osi interface
    OsiClpSolverInterface * clpSolver 
      = dynamic_cast<OsiClpSolverInterface *> (cleanModel);
    if (clpSolver) {
      ClpSimplex * clp = clpSolver->getModelPtr();
      clp->messageHandler()->setLogLevel(cleanModel->messageHandler()->logLevel());
      ClpPresolve pinfo;
      //printf("integerPresolve - temp switch off doubletons\n");
      //pinfo.setPresolveActions(4);
      ClpSimplex * model2 = pinfo.presolvedModel(*clp,1.0e-8);
      if (!model2) {
	// presolve found to be infeasible
	feasible=false;
      } else {
	// update original array
	const int * originalColumns = pinfo.originalColumns();
	// just slot in new solver
	OsiClpSolverInterface * temp = new OsiClpSolverInterface(model2,true);
	numberColumns = temp->getNumCols();
	for (iColumn=0;iColumn<originalNumberColumns;iColumn++)
	  original[iColumn]=-1;
	for (iColumn=0;iColumn<numberColumns;iColumn++)
	  original[originalColumns[iColumn]]=iColumn;
	// copy parameters
	temp->copyParameters(*solver_);
	// and specialized ones
	temp->setSpecialOptions(clpSolver->specialOptions());
	delete solver_;
	solver_ = temp;
	setCutoff(cutoff);
	deleteObjects();
	if (!numberObjects_) {
	  // Nothing left
	  doIntegerPresolve=false;
	  weak=true;
	  break;
	}
	synchronizeModel(); // make sure everything that needs solver has it
	// just point to solver_
	continuousSolver_ = solver_;
	feasible=resolve(NULL,3);
	if (!feasible||!doIntegerPresolve||weak) break;
	// see if we can get solution by heuristics
	int found=-1;
	int iHeuristic;
	double * newSolution = new double [numberColumns];
	double heuristicValue=getCutoff();
	for (iHeuristic=0;iHeuristic<numberHeuristics_;iHeuristic++) {
	  double saveValue=heuristicValue;
	  int ifSol = heuristic_[iHeuristic]->solution(heuristicValue,
						       newSolution);
	  if (ifSol>0) {
	    // better solution found
	    found=iHeuristic;
            incrementUsed(newSolution);
	  } else if (ifSol<0) {
	    heuristicValue = saveValue;
	  }
	}
	if (found >= 0) {
	  // We probably already have a current solution, but just in case ...
	  int numberColumns = getNumCols() ;
	  if (!currentSolution_)
	    currentSolution_ = new double[numberColumns] ;
          testSolution_=currentSolution_;
	  // better solution save
          lastHeuristic_ = heuristic_[found];
	  setBestSolution(CBC_ROUNDING,heuristicValue,
			  newSolution);
	  // update cutoff
	  cutoff = getCutoff();
	}
	delete [] newSolution;
	// Space for type of cuts
	maximumWhich_=1000;
        delete [] whichGenerator_ ;
	whichGenerator_ = new int[maximumWhich_];
	// save number of rows
	numberRowsAtContinuous_ = getNumRows();
	maximumNumberCuts_=0;
	currentNumberCuts_=0;
	delete [] addedCuts_;
	addedCuts_ = NULL;
	
	// maximum depth for tree walkback
	maximumDepth_=10;
	delete [] walkback_;
	walkback_ = new CbcNodeInfo * [maximumDepth_];
	
	OsiCuts cuts;
	numberOldActiveCuts_=0;
	numberNewCuts_ = 0;
	feasible = solveWithCuts(cuts,maximumCutPassesAtRoot_,NULL);
	currentNumberCuts_=numberNewCuts_;
	delete [] whichGenerator_;
        whichGenerator_=NULL;
	delete [] walkback_;
	walkback_ = NULL;
	delete [] addedCuts_;
	addedCuts_=NULL;
	if (feasible) {
	  // fix anything in original which integer presolve fixed
	  // for now just integers
	  const double * lower = solver_->getColLower();
	  const double * upper = solver_->getColUpper();
	  int i;
	  for (i=0;i<originalNumberIntegers;i++) {
	    iColumn = originalIntegers[i];
	    int jColumn = original[iColumn];
	    if (jColumn >= 0) {
	      if (upper[jColumn]<originalUpper[iColumn]) 
		originalUpper[iColumn]	= upper[jColumn];
	      if (lower[jColumn]>originalLower[iColumn]) 
		originalLower[iColumn]	= lower[jColumn];
	    }
	  }
	}
      }
    }
#endif
    if (!feasible||!doIntegerPresolve) {
      break;
    }
  }
  //solver_->writeMps("xx");
  delete cleanModel;
  delete [] originalIntegers;
  numberColumns = getNumCols();
  delete [] originalColumns_;
  originalColumns_ = new int[numberColumns];
  numberColumns=0;
  for (iColumn=0;iColumn<originalNumberColumns;iColumn++) {
    int jColumn = original[iColumn];
    if (jColumn >= 0) 
      originalColumns_[numberColumns++]=iColumn;
  }
  delete [] original;
  delete [] originalLower;
  delete [] originalUpper;
  
  deleteObjects();
  synchronizeModel(); // make sure everything that needs solver has it
  continuousSolver_=NULL;
  currentNumberCuts_=0;
  return feasible;
}
// Put back information into original model - after integerpresolve 
void 
CbcModel::originalModel(CbcModel * presolvedModel,bool weak)
{
  solver_->copyParameters(*(presolvedModel->solver_));
  bestObjective_ = presolvedModel->bestObjective_;
  delete [] bestSolution_;
  findIntegers(true);
  if (presolvedModel->bestSolution_) {
    int numberColumns = getNumCols();
    int numberOtherColumns = presolvedModel->getNumCols();
    //bestSolution_ = new double[numberColumns];
    // set up map
    int * back = new int[numberColumns];
    int i;
    for (i=0;i<numberColumns;i++)
      back[i]=-1;
    for (i=0;i<numberOtherColumns;i++)
      back[presolvedModel->originalColumns_[i]]=i;
    int iColumn;
    // set ones in presolved model to values
    double * otherSolution = presolvedModel->bestSolution_;
    //const double * lower = getColLower();
    for (i=0;i<numberIntegers_;i++) {
      iColumn = integerVariable_[i];
      int jColumn = back[iColumn];
      //bestSolution_[iColumn]=lower[iColumn];
      if (jColumn >= 0) {
	double value=floor(otherSolution[jColumn]+0.5);
	solver_->setColLower(iColumn,value);
	solver_->setColUpper(iColumn,value);
	//bestSolution_[iColumn]=value;
      }
    }
    delete [] back;
#if 0
    // ** looks as if presolve needs more intelligence
    // do presolve - for now just clp but easy to get osi interface
    OsiClpSolverInterface * clpSolver 
      = dynamic_cast<OsiClpSolverInterface *> (solver_);
    assert (clpSolver);
    ClpSimplex * clp = clpSolver->getModelPtr();
    Presolve pinfo;
    ClpSimplex * model2 = pinfo.presolvedModel(*clp,1.0e-8);
    model2->primal(1);
    pinfo.postsolve(true);
    const double * solution = solver_->getColSolution();
    for (i=0;i<numberIntegers_;i++) {
      iColumn = integerVariable_[i];
      double value=floor(solution[iColumn]+0.5);
      solver_->setColLower(iColumn,value);
      solver_->setColUpper(iColumn,value);
    }
#else
    if (!weak) {
      // for now give up
      int save = numberCutGenerators_;
      numberCutGenerators_=0;
      bestObjective_=1.0e100;
      branchAndBound();
      numberCutGenerators_=save;
    }
#endif
    if (bestSolution_) {
      // solve problem
      resolve(NULL,3);
      // should be feasible
      if (!currentSolution_)
	currentSolution_ = new double[numberColumns] ;
      testSolution_ = currentSolution_;
#ifndef NDEBUG
      int numberIntegerInfeasibilities;
      int numberObjectInfeasibilities;
      assert(feasibleSolution(numberIntegerInfeasibilities,
			      numberObjectInfeasibilities));
#endif
    }
  } else {
    bestSolution_=NULL;
  }
  numberSolutions_=presolvedModel->numberSolutions_;
  numberHeuristicSolutions_=presolvedModel->numberHeuristicSolutions_;
  numberNodes_ = presolvedModel->numberNodes_;
  numberIterations_ = presolvedModel->numberIterations_;
  status_ = presolvedModel->status_;
  secondaryStatus_ = presolvedModel->secondaryStatus_;
  synchronizeModel();
} 
// Pass in Message handler (not deleted at end)
void 
CbcModel::passInMessageHandler(CoinMessageHandler * handler)
{
  if (defaultHandler_) {
    delete handler_;
    handler_ = NULL;
  }
  defaultHandler_=false;
  handler_=handler;
}
void 
CbcModel::passInTreeHandler(CbcTree & tree)
{
  delete tree_;
  tree_ = tree.clone();
}
// Make sure region there
void 
CbcModel::reserveCurrentSolution(const double * solution)
{
  int numberColumns = getNumCols() ;
  if (!currentSolution_)
    currentSolution_ = new double[numberColumns] ;
  testSolution_=currentSolution_;
  if (solution)
    memcpy(currentSolution_,solution,numberColumns*sizeof(double));
}
/* For passing in an CbcModel to do a sub Tree (with derived tree handlers).
   Passed in model must exist for duration of branch and bound
*/
void 
CbcModel::passInSubTreeModel(CbcModel & model)
{
  subTreeModel_=&model;
}
// For retrieving a copy of subtree model with given OsiSolver or NULL
CbcModel * 
CbcModel::subTreeModel(OsiSolverInterface * solver) const
{
  const CbcModel * subModel=subTreeModel_;
  if (!subModel)
    subModel=this;
  // Get new copy
  CbcModel * newModel = new CbcModel(*subModel);
  if (solver)
    newModel->assignSolver(solver);
  return newModel;
}
//#############################################################################
// Set/Get Application Data
// This is a pointer that the application can store into and retrieve
// from the solverInterface.
// This field is the application to optionally define and use.
//#############################################################################

void CbcModel::setApplicationData(void * appData)
{
  appData_ = appData;
}
//-----------------------------------------------------------------------------
void * CbcModel::getApplicationData() const
{
  return appData_;
}
/*  create a submodel from partially fixed problem

The method creates a new clean model with given bounds.
*/
CbcModel *  
CbcModel::cleanModel(const double * lower, const double * upper)
{
  OsiSolverInterface * solver = continuousSolver_->clone();

  int numberIntegers = numberIntegers_;
  const int * integerVariable = integerVariable_;
  
  int i;
  for (i=0;i<numberIntegers;i++) {
    int iColumn=integerVariable[i];
    const OsiObject * object = object_[i];
    const CbcSimpleInteger * integerObject = 
      dynamic_cast<const  CbcSimpleInteger *> (object);
    assert(integerObject);
    // get original bounds
    double originalLower = integerObject->originalLowerBound();
    double originalUpper = integerObject->originalUpperBound();
    solver->setColLower(iColumn,CoinMax(lower[iColumn],originalLower));
    solver->setColUpper(iColumn,CoinMin(upper[iColumn],originalUpper));
  }
  CbcModel * model = new CbcModel(*solver);
  // off some messages
  if (handler_->logLevel()<=1) {
    model->messagesPointer()->setDetailMessage(3,9);
    model->messagesPointer()->setDetailMessage(3,6);
    model->messagesPointer()->setDetailMessage(3,4);
    model->messagesPointer()->setDetailMessage(3,1);
    model->messagesPointer()->setDetailMessage(3,13);
    model->messagesPointer()->setDetailMessage(3,14);
    model->messagesPointer()->setDetailMessage(3,3007);
  }
  // Cuts
  for ( i = 0;i<numberCutGenerators_;i++) {
    int howOften = generator_[i]->howOftenInSub();
    if (howOften>-100) {
      CbcCutGenerator * generator = virginGenerator_[i];
      CglCutGenerator * cglGenerator = generator->generator();
      model->addCutGenerator(cglGenerator,howOften,
			      generator->cutGeneratorName(),
			      generator->normal(),
			      generator->atSolution(),
			      generator->whenInfeasible(),
			      -100, generator->whatDepthInSub(),-1);
    }
  }
  double cutoff = getCutoff();
  model->setCutoff(cutoff);
  return model;
}
/* Invoke the branch & cut algorithm on partially fixed problem
   
   The method uses a subModel created by cleanModel. The search 
   ends when the tree is exhausted or maximum nodes is reached.

   If better solution found then it is saved.
   
   Returns 0 if search completed and solution, 1 if not completed and solution,
   2 if completed and no solution, 3 if not completed and no solution.
   
   Normally okay to do subModel immediately followed by subBranchandBound
   (== other form of subBranchAndBound)
   but may need to get at model for advanced features.
   
   Deletes model
   
*/
  
int 
CbcModel::subBranchAndBound(CbcModel * model,
                            CbcModel * presolvedModel,
			    int maximumNodes)
{
  int i;
  double cutoff=model->getCutoff();
  CbcModel * model2;
  if (presolvedModel) 
    model2=presolvedModel;
  else
    model2=model;
  // Do complete search
  
  for (i=0;i<numberHeuristics_;i++) {
    model2->addHeuristic(heuristic_[i]);
    model2->heuristic(i)->resetModel(model2);
  }
  // Definition of node choice
  model2->setNodeComparison(nodeCompare_->clone());
  //model2->solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
  model2->messageHandler()->setLogLevel(CoinMax(0,handler_->logLevel()-1));
  //model2->solver()->messageHandler()->setLogLevel(2);
  model2->setMaximumCutPassesAtRoot(maximumCutPassesAtRoot_);
  model2->setPrintFrequency(50);
  model2->setIntParam(CbcModel::CbcMaxNumNode,maximumNodes);
  model2->branchAndBound();
  delete model2->nodeComparison();
  if (model2->getMinimizationObjValue()>cutoff) {
    // no good
    if (model!=model2)
      delete model2;
    delete model;
    return 2;
  }
  if (model!=model2) {
    // get back solution
    model->originalModel(model2,false);
    delete model2;
  }
  int status;
  if (model->getMinimizationObjValue()<cutoff&&model->bestSolution()) {
    double objValue = model->getObjValue();
    const double * solution = model->bestSolution();
    setBestSolution(CBC_TREE_SOL,objValue,solution);
    status = 0;
  } else {
    status=2;
  }
  if (model->status())
    status ++ ; // not finished search
  delete model;
  return status;
}
/* Invoke the branch & cut algorithm on partially fixed problem
   
The method creates a new model with given bounds, presolves it
then proceeds to explore the branch & cut search tree. The search 
ends when the tree is exhausted or maximum nodes is reached.
Returns 0 if search completed and solution, 1 if not completed and solution,
2 if completed and no solution, 3 if not completed and no solution.
*/
int 
CbcModel::subBranchAndBound(const double * lower, const double * upper,
			    int maximumNodes)
{
  OsiSolverInterface * solver = continuousSolver_->clone();

  int numberIntegers = numberIntegers_;
  const int * integerVariable = integerVariable_;
  
  int i;
  for (i=0;i<numberIntegers;i++) {
    int iColumn=integerVariable[i];
    const OsiObject * object = object_[i];
    const CbcSimpleInteger * integerObject = 
      dynamic_cast<const  CbcSimpleInteger *> (object);
    assert(integerObject);
    // get original bounds
    double originalLower = integerObject->originalLowerBound();
    double originalUpper = integerObject->originalUpperBound();
    solver->setColLower(iColumn,CoinMax(lower[iColumn],originalLower));
    solver->setColUpper(iColumn,CoinMin(upper[iColumn],originalUpper));
  }
  CbcModel model(*solver);
  // off some messages
  if (handler_->logLevel()<=1) {
    model.messagesPointer()->setDetailMessage(3,9);
    model.messagesPointer()->setDetailMessage(3,6);
    model.messagesPointer()->setDetailMessage(3,4);
    model.messagesPointer()->setDetailMessage(3,1);
    model.messagesPointer()->setDetailMessage(3,3007);
  }
  double cutoff = getCutoff();
  model.setCutoff(cutoff);
  // integer presolve
  CbcModel * model2 = model.integerPresolve(false);
  if (!model2||!model2->getNumRows()) {
    delete model2;
    delete solver;
    return 2;
  }
  if (handler_->logLevel()>1)
    printf("Reduced model has %d rows and %d columns\n",
	   model2->getNumRows(),model2->getNumCols());
  // Do complete search
  
  // Cuts
  for ( i = 0;i<numberCutGenerators_;i++) {
    int howOften = generator_[i]->howOftenInSub();
    if (howOften>-100) {
      CbcCutGenerator * generator = virginGenerator_[i];
      CglCutGenerator * cglGenerator = generator->generator();
      model2->addCutGenerator(cglGenerator,howOften,
			      generator->cutGeneratorName(),
			      generator->normal(),
			      generator->atSolution(),
			      generator->whenInfeasible(),
			      -100, generator->whatDepthInSub(),-1);
    }
  }
  for (i=0;i<numberHeuristics_;i++) {
    model2->addHeuristic(heuristic_[i]);
    model2->heuristic(i)->resetModel(model2);
  }
  // Definition of node choice
  model2->setNodeComparison(nodeCompare_->clone());
  //model2->solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
  model2->messageHandler()->setLogLevel(CoinMax(0,handler_->logLevel()-1));
  //model2->solver()->messageHandler()->setLogLevel(2);
  model2->setMaximumCutPassesAtRoot(maximumCutPassesAtRoot_);
  model2->setPrintFrequency(50);
  model2->setIntParam(CbcModel::CbcMaxNumNode,maximumNodes);
  model2->branchAndBound();
  delete model2->nodeComparison();
  if (model2->getMinimizationObjValue()>cutoff) {
    // no good
    delete model2;
    delete solver;
    return 2;
  }
  // get back solution
  model.originalModel(model2,false);
  delete model2;
  int status;
  if (model.getMinimizationObjValue()<cutoff&&model.bestSolution()) {
    double objValue = model.getObjValue();
    const double * solution = model.bestSolution();
    setBestSolution(CBC_TREE_SOL,objValue,solution);
    status = 0;
  } else {
    status=2;
  }
  if (model.status())
    status ++ ; // not finished search
  delete solver;
  return status;
}
// Set a pointer to a row cut which will be added instead of normal branching.
void 
CbcModel::setNextRowCut(const OsiRowCut & cut)
{ 
  nextRowCut_=new OsiRowCut(cut);
  nextRowCut_->setEffectiveness(COIN_DBL_MAX); // mark so will always stay
}
/* Process root node and return a strengthened model
   
The method assumes that initialSolve() has been called to solve the
LP relaxation. It processes the root node and then returns a pointer
to the strengthened model (or NULL if infeasible)
*/
OsiSolverInterface *  
CbcModel::strengthenedModel()
{
/*
  Switch off heuristics
*/
  int saveNumberHeuristics=numberHeuristics_;
  numberHeuristics_=0;
/*
  Scan the variables, noting the integer variables. Create an
  CbcSimpleInteger object for each integer variable.
*/
  findIntegers(false) ;
/*
  Ensure that objects on the lists of OsiObjects, heuristics, and cut
  generators attached to this model all refer to this model.
*/
  synchronizeModel() ;

  // Set so we can tell we are in initial phase in resolve
  continuousObjective_ = -COIN_DBL_MAX ;
/*
  Solve the relaxation.

  Apparently there are circumstances where this will be non-trivial --- i.e.,
  we've done something since initialSolve that's trashed the solution to the
  continuous relaxation.
*/
  bool feasible = resolve(NULL,0) != 0 ;
/*
  If the linear relaxation of the root is infeasible, bail out now. Otherwise,
  continue with processing the root node.
*/
  if (!feasible)
  { handler_->message(CBC_INFEAS,messages_)<< CoinMessageEol ;
    return NULL; }
  // Save objective (just so user can access it)
  originalContinuousObjective_ = solver_->getObjValue();

/*
  Begin setup to process a feasible root node.
*/
  bestObjective_ = CoinMin(bestObjective_,1.0e50) ;
  numberSolutions_ = 0 ;
  numberHeuristicSolutions_ = 0 ;
  // Everything is minimization
  double cutoff=getCutoff() ;
  double direction = solver_->getObjSense() ;
  if (cutoff < 1.0e20&&direction<0.0)
    messageHandler()->message(CBC_CUTOFF_WARNING1,
				    messages())
				      << cutoff << -cutoff << CoinMessageEol ;
  if (cutoff > bestObjective_)
    cutoff = bestObjective_ ;
  setCutoff(cutoff) ;
/*
  We probably already have a current solution, but just in case ...
*/
  int numberColumns = getNumCols() ;
  if (!currentSolution_)
    currentSolution_ = new double[numberColumns] ;
  testSolution_=currentSolution_;
/*
  Create a copy of the solver, thus capturing the original (root node)
  constraint system (aka the continuous system).
*/
  continuousSolver_ = solver_->clone() ;
  numberRowsAtContinuous_ = getNumRows() ;
/*
  Check the objective to see if we can deduce a nontrivial increment. If
  it's better than the current value for CbcCutoffIncrement, it'll be
  installed.
*/
  analyzeObjective() ;
/*
  Set up for cut generation. addedCuts_ holds the cuts which are relevant for
  the active subproblem. whichGenerator will be used to record the generator
  that produced a given cut.
*/
  maximumWhich_ = 1000 ;
  delete [] whichGenerator_ ;
  whichGenerator_ = new int[maximumWhich_] ;
  maximumNumberCuts_ = 0 ;
  currentNumberCuts_ = 0 ;
  delete [] addedCuts_ ;
  addedCuts_ = NULL ;
  /*  
  Generate cuts at the root node and reoptimise. solveWithCuts does the heavy
  lifting. It will iterate a generate/reoptimise loop (including reduced cost
  fixing) until no cuts are generated, the change in objective falls off,  or
  the limit on the number of rounds of cut generation is exceeded.

  At the end of all this, any cuts will be recorded in cuts and also
  installed in the solver's constraint system. We'll have reoptimised, and
  removed any slack cuts (numberOldActiveCuts_ and numberNewCuts_ have been
  adjusted accordingly).

  Tell cut generators they can be a bit more aggressive at root node

*/
  int iCutGenerator;
  for (iCutGenerator = 0;iCutGenerator<numberCutGenerators_;iCutGenerator++) {
    CglCutGenerator * generator = generator_[iCutGenerator]->generator();
    generator->setAggressiveness(generator->getAggressiveness()+100);
  }
  OsiCuts cuts ;
  numberOldActiveCuts_ = 0 ;
  numberNewCuts_ = 0 ;
  { int iObject ;
    int preferredWay ;
    int numberUnsatisfied = 0 ;
    memcpy(currentSolution_,solver_->getColSolution(),
	   numberColumns*sizeof(double)) ;

    // point to useful information
    OsiBranchingInformation usefulInfo=usefulInformation();
    for (iObject = 0 ; iObject < numberObjects_ ; iObject++)
    { double infeasibility =
	  object_[iObject]->infeasibility(&usefulInfo,preferredWay) ;
      if (infeasibility) numberUnsatisfied++ ; }
    if (numberUnsatisfied)
    { feasible = solveWithCuts(cuts,maximumCutPassesAtRoot_,
			       NULL) ; } }
/*
  We've taken the continuous relaxation as far as we can. 
*/

  OsiSolverInterface * newSolver=NULL;
  if (feasible) {
    // make copy of current solver
    newSolver = solver_->clone();
  }
/*
  Clean up dangling objects. continuousSolver_ may already be toast.
*/
  delete [] whichGenerator_ ;
  whichGenerator_ = NULL;
  delete [] walkback_ ;
  walkback_ = NULL ;
  delete [] addedCuts_ ;
  addedCuts_ = NULL ;
  if (continuousSolver_)
  { delete continuousSolver_ ;
    continuousSolver_ = NULL ; }
/*
  Destroy global cuts by replacing with an empty OsiCuts object.
*/
  globalCuts_= OsiCuts() ;
  numberHeuristics_ = saveNumberHeuristics;
  
  return newSolver; 
}
// Just update objectiveValue
void CbcModel::setBestObjectiveValue( double objectiveValue)
{
  bestObjective_=objectiveValue;
}
double 
CbcModel::getBestPossibleObjValue() const
{ 
  return CoinMin(bestPossibleObjective_,bestObjective_) * solver_->getObjSense() ;
}
// Make given rows (L or G) into global cuts and remove from lp
void 
CbcModel::makeGlobalCuts(int number,const int * which)
{
  const double * rowLower = solver_->getRowLower();
  const double * rowUpper = solver_->getRowUpper();

  int numberRows = solver_->getNumRows();

  // Row copy
  const double * elementByRow = solver_->getMatrixByRow()->getElements();
  const int * column = solver_->getMatrixByRow()->getIndices();
  const CoinBigIndex * rowStart = solver_->getMatrixByRow()->getVectorStarts();
  const int * rowLength = solver_->getMatrixByRow()->getVectorLengths();

  // Not all rows may be good so we need new array
  int * whichDelete = new int[numberRows];
  int nDelete=0;
  for (int i=0;i<number;i++) {
    int iRow = which[i];
    if (iRow>=0&&iRow<numberRows) {
      if (rowLower[iRow]<-1.0e20||rowUpper[iRow]>1.0e20) {
        whichDelete[nDelete++]=iRow;
        OsiRowCut  thisCut;
        thisCut.setLb(rowLower[iRow]);
        thisCut.setUb(rowUpper[iRow]);
        int start = rowStart[iRow];
        thisCut.setRow(rowLength[iRow],column+start,elementByRow+start,false);
	thisCut.setGloballyValid(2);
        globalCuts_.insert(thisCut) ;
      }
    }
  }
  if (nDelete)
    solver_->deleteRows(nDelete,whichDelete);
  delete [] whichDelete;
}
// Make given cut into a global cut
void 
CbcModel::makeGlobalCut(const OsiRowCut * cut)
{
  OsiRowCut newCut(*cut);
  newCut.setGloballyValidAsInteger(2);
  newCut.mutableRow().setTestForDuplicateIndex(false);
  globalCuts_.insert(newCut) ;
}
// Make given cut into a global cut
void 
CbcModel::makeGlobalCut(const OsiRowCut & cut)
{
  OsiRowCut newCut(cut);
  newCut.setGloballyValid(2);
  newCut.mutableRow().setTestForDuplicateIndex(false);
  globalCuts_.insert(newCut) ;
}
void 
CbcModel::setNodeComparison(CbcCompareBase * compare)
{ 
  delete nodeCompare_;
  nodeCompare_ = compare->clone();
}
void 
CbcModel::setNodeComparison(CbcCompareBase & compare)
{ 
  delete nodeCompare_;
  nodeCompare_ = compare.clone();
}
void 
CbcModel::setProblemFeasibility(CbcFeasibilityBase * feasibility)
{
  delete problemFeasibility_;
  problemFeasibility_ = feasibility->clone();
}
void 
CbcModel::setProblemFeasibility(CbcFeasibilityBase & feasibility)
{
  delete problemFeasibility_;
  problemFeasibility_ = feasibility.clone();
}
// Set the strategy. Clones
void 
CbcModel::setStrategy(CbcStrategy & strategy)
{
  delete strategy_;
  strategy_ = strategy.clone();
}
// Increases usedInSolution for nonzeros
void 
CbcModel::incrementUsed(const double * solution)
{
  // might as well mark all including continuous
  int numberColumns = solver_->getNumCols();
  for (int i=0;i<numberColumns;i++) {
    if (solution[i])
      usedInSolution_[i]++;
  }
}
// Are there numerical difficulties (for initialSolve) ?
bool 
CbcModel::isInitialSolveAbandoned() const 
{
  if (status_!=-1) {
    return false;
  } else {
    return solver_->isAbandoned();
  }
}
// Is optimality proven (for initialSolve) ?
bool 
CbcModel::isInitialSolveProvenOptimal() const 
{
  if (status_!=-1) {
    return originalContinuousObjective_<1.0e50;
  } else {
    return solver_->isProvenOptimal();
  }
}
// Is primal infeasiblity proven (for initialSolve) ?
bool 
CbcModel::isInitialSolveProvenPrimalInfeasible() const 
{
  if (status_!=-1) {
    if (status_==0&&secondaryStatus_==7)
      return false;
    else
      return originalContinuousObjective_>=1.0e50;
  } else {
    return solver_->isProvenPrimalInfeasible();
  }
}
// Is dual infeasiblity proven (for initialSolve) ?
bool 
CbcModel::isInitialSolveProvenDualInfeasible() const 
{
  if (status_!=-1) {
    if (status_==0&&secondaryStatus_==7)
      return true;
    else
      return false;
  } else {
    return solver_->isProvenDualInfeasible();
  }
}
// Set pointers for speed
void 
CbcModel::setPointers(const OsiSolverInterface * solver)
{
  /// Pointer to array[getNumCols()] (for speed) of column lower bounds
  cbcColLower_ = solver_->getColLower();
  /// Pointer to array[getNumCols()] (for speed) of column upper bounds
  cbcColUpper_ = solver_->getColUpper();
  /// Pointer to array[getNumRows()] (for speed) of row lower bounds
  cbcRowLower_ = solver_->getRowLower();
  /// Pointer to array[getNumRows()] (for speed) of row upper bounds
  cbcRowUpper_ = solver_->getRowUpper();
  /// Pointer to array[getNumCols()] (for speed) of primal solution vector
  cbcColSolution_ = solver_->getColSolution();
  /// Pointer to array[getNumRows()] (for speed) of dual prices
  cbcRowPrice_ = solver_->getRowPrice();
  /// Get a pointer to array[getNumCols()] (for speed) of reduced costs
  if(solverCharacteristics_&&solverCharacteristics_->reducedCostsAccurate())
    cbcReducedCost_ = solver_->getReducedCost();
  else
    cbcReducedCost_ = NULL;
  /// Pointer to array[getNumRows()] (for speed) of row activity levels.
  cbcRowActivity_ = solver_->getRowActivity();
  dblParam_[CbcCurrentObjectiveValue]=solver->getObjValue();
  dblParam_[CbcCurrentMinimizationObjectiveValue]=
    dblParam_[CbcCurrentObjectiveValue]* 
    dblParam_[CbcOptimizationDirection];
}

/*
  Delete any existing handler and create a clone of the one supplied.
*/
void CbcModel::passInEventHandler (const CbcEventHandler *eventHandler)
{
  delete eventHandler_;
  eventHandler_ = NULL ;
  if (eventHandler)
    eventHandler_ = eventHandler->clone();
}

/*
  CbcEventHandler* CbcModel::eventHandler is inlined in CbcModel.hpp.
*/

// Encapsulates solver resolve
int 
CbcModel::resolve(OsiSolverInterface * solver)
{
#ifdef COIN_HAS_CLP
  OsiClpSolverInterface * clpSolver 
    = dynamic_cast<OsiClpSolverInterface *> (solver);
  if (clpSolver) {
    /*bool takeHint;
    OsiHintStrength strength;
    bool gotHint = (clpSolver->getHintParam(OsiDoDualInResolve,takeHint,strength));
    assert (gotHint);
    int algorithm=-1;
    if (strength!=OsiHintIgnore)
      algorithm = takeHint ? -1 : 1;
      assert (algorithm==-1);*/
    //clpSolver->setHintParam(OsiDoDualInResolve,true,OsiHintTry);
    ClpSimplex * clpSimplex = clpSolver->getModelPtr();
    int save = clpSimplex->specialOptions();
    clpSimplex->setSpecialOptions(save|0x01000000); // say is Cbc (and in branch and bound)
    clpSolver->resolve();
    clpSimplex->setSpecialOptions(save);
  } else {
    solver->resolve();
  }
#else
  solver->resolve();
#endif
  return solver->isProvenOptimal() ? 1 : 0;
}

// Set log level
void 
CbcModel::setLogLevel(int value)
{ 
  handler_->setLogLevel(value);
  // Reduce print out in Osi
  if (solver_) {
    int oldLevel = solver_->messageHandler()->logLevel();
    if (value<oldLevel)
      solver_->messageHandler()->setLogLevel(value);
#ifdef COIN_HAS_CLP
    OsiClpSolverInterface * clpSolver 
      = dynamic_cast<OsiClpSolverInterface *> (solver_);
    if (clpSolver) {
      ClpSimplex * clpSimplex = clpSolver->getModelPtr();
      oldLevel = clpSimplex->logLevel();
      if (value<oldLevel)
        clpSimplex->setLogLevel(value);
    }
#else		// COIN_HAS_CLP
/*
  For generic OSI solvers, if the new log level is 0, try the
  DoReducePrint hint for emphasis.
*/
    if (value == 0) {
      solver_->setHintParam(OsiDoReducePrint,true,OsiHintDo) ;
    }
#endif		// COIN_HAS_CLP
  }
}

/* Pass in target solution and optional priorities.
   If priorities then >0 means only branch if incorrect
   while <0 means branch even if correct. +1 or -1 are
   highest priority */
void 
CbcModel::setHotstartSolution(const double * solution, const int * priorities)
{ 
  if (solution==NULL) {
    delete [] hotstartSolution_;
    hotstartSolution_=NULL;
    delete [] hotstartPriorities_;
    hotstartPriorities_=NULL;
  } else {
    int numberColumns = solver_->getNumCols();
    hotstartSolution_ = CoinCopyOfArray(solution,numberColumns);
    for (int i=0;i<numberColumns;i++) {
      if (solver_->isInteger(i)) 
        hotstartSolution_[i]=floor(hotstartSolution_[i]+0.5);
    }
    hotstartPriorities_ = CoinCopyOfArray(priorities,numberColumns);
  }
}
// Increment strong info
void 
CbcModel::incrementStrongInfo(int numberTimes, int numberIterations,
                           int numberFixed, bool ifInfeasible)
{
  strongInfo_[0] += numberTimes;
  numberStrongIterations_ += numberIterations;
  strongInfo_[1] += numberFixed;
  if (ifInfeasible) 
    strongInfo_[2] ++;
}
/* Set objective value in a node.  This is separated out so that
   odd solvers can use.  It may look at extra information in
   solverCharacteriscs_ and will also use bound from parent node
*/
void 
CbcModel::setObjectiveValue(CbcNode * thisNode, const CbcNode * parentNode) const
{
  double newObjValue = solver_->getObjSense()*solver_->getObjValue();
  // If odd solver take its bound
  if (solverCharacteristics_) {
    newObjValue = CoinMax(newObjValue,solverCharacteristics_->mipBound());
    // Reset bound anyway (no harm if not odd)
    solverCharacteristics_->setMipBound(-COIN_DBL_MAX);
  }
  // If not root then use max of this and parent
  if (parentNode)
    newObjValue = CoinMax(newObjValue,parentNode->objectiveValue());
  thisNode->setObjectiveValue(newObjValue);
}
// Current time since start of branchAndbound
double 
CbcModel::getCurrentSeconds() const {
  return CoinCpuTime()-getDblParam(CbcStartSeconds);
}
/* Encapsulates choosing a variable -
   anyAction -2, infeasible (-1 round again), 0 done
*/
int 
CbcModel::chooseBranch(CbcNode * newNode, int numberPassesLeft, 
		       CbcNode * oldNode, OsiCuts & cuts,
		       bool & resolved, CoinWarmStartBasis *lastws,
		       const double * lowerBefore,const double * upperBefore,
		       OsiSolverBranch * & branches)
{
  int anyAction =-1 ;
  resolved = false ;
  if (newNode->objectiveValue() >= getCutoff()) 
    anyAction=-2;
  branches=NULL;
  bool feasible=true;
  int branchingState=-1;
  currentNode_=newNode; // so can be used elsewhere
  while (anyAction == -1) {
    // Set objective value (not so obvious if NLP etc)
    setObjectiveValue(newNode,oldNode);
    //if (numberPassesLeft<=0)
    //branchingState=1;
    if (!branchingMethod_||!branchingMethod_->chooseMethod()) {
      if (numberBeforeTrust_==0 ) {
	anyAction = newNode->chooseBranch(this,oldNode,numberPassesLeft) ;
      } else {
	anyAction = newNode->chooseDynamicBranch(this,oldNode,branches,numberPassesLeft) ;
      if (anyAction==-3) 
	anyAction = newNode->chooseBranch(this,oldNode,numberPassesLeft) ; // dynamic did nothing
      }
    } else {
      OsiBranchingInformation usefulInfo=usefulInformation();
      anyAction = newNode->chooseOsiBranch(this,oldNode,&usefulInfo,branchingState) ;; // Osi method
      //branchingState=0;
    }
    if (solverCharacteristics_ && 
	solverCharacteristics_->solutionAddsCuts() && // we are in some OA based bab
	feasible && (newNode->numberUnsatisfied()==0) //solution has become integer feasible during strong branching
	) { 
      //in the present case we need to check here integer infeasibility if the node is not fathomed we will have to do the loop
      // again
      //std::cout<<solver_<<std::endl;
      resolve(solver_);
      double objval = solver_->getObjValue();
      lastHeuristic_ = NULL;
      setBestSolution(CBC_SOLUTION, objval,
		      solver_->getColSolution()) ;
      int easy=2;
      if (!solverCharacteristics_->mipFeasible())//did we prove that the node could be pruned?
	feasible = false;
      // Reset the bound now
      solverCharacteristics_->setMipBound(-COIN_DBL_MAX);
      
      
      solver_->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,&easy) ;
      feasible &= resolve(oldNode ? oldNode->nodeInfo() : NULL,11) != 0 ;
      solver_->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,NULL) ;
      resolved = true ;
      if (problemFeasibility_->feasible(this,0)<0) {
	feasible=false; // pretend infeasible
      }
      if(feasible)
	anyAction = -1;
      else
	anyAction = -2;
    }
    /*
      Yep, false positives for sure. And no easy way to distinguish honest
      infeasibility from `found a solution and tightened objective target.'
      
      if (onOptimalPath)
      assert (anyAction!=-2); // can be useful but gives false positives on strong
    */
    numberPassesLeft--;
    if (numberPassesLeft<=-1) {
      if (!numberLongStrong_&&!numberThreads_)
	messageHandler()->message(CBC_WARNING_STRONG,
				  messages()) << CoinMessageEol ;
      numberLongStrong_++;
    }
    if (anyAction == -1) {
      // can do quick optimality check
      int easy=2;
      solver_->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,&easy) ;
      feasible = resolve(oldNode ? oldNode->nodeInfo() : NULL,11) != 0 ;
      solver_->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo,NULL) ;
      resolved = true ;
      if (problemFeasibility_->feasible(this,0)<0) {
	feasible=false; // pretend infeasible
      }
      if (feasible) { 
	// Set objective value (not so obvious if NLP etc)
	setObjectiveValue(newNode,oldNode);
	reducedCostFix() ;
	if (newNode->objectiveValue() >= getCutoff()) 
	  anyAction=-2;
      } else { 
	anyAction = -2 ;
      }
    }
  }
  if (anyAction >= 0) {
    if (resolved) { 
/*
  Used to be that when the node was not fathomed (branching object present)
  the solution was not needed. But that's no longer the case --- heuristics
  are applied, and they may want the solution.
*/
      // bool needValidSolution = (newNode->branchingObject() == NULL) ;
      bool needValidSolution = true ;
      takeOffCuts(cuts,needValidSolution,NULL) ; 
#	      ifdef CHECK_CUT_COUNTS
      {
	printf("Number of rows after chooseBranch fix (node)"
	       "(active only) %d\n",
	       numberRowsAtContinuous_+numberNewCuts_+
	       numberOldActiveCuts_) ;
	const CoinWarmStartBasis* debugws =
	  dynamic_cast<const CoinWarmStartBasis*>
	  (solver_->getWarmStart()) ;
	debugws->print() ;
	delete debugws ; 
      }
#	      endif
    }
    newNode->createInfo(this,oldNode,lastws,lowerBefore,upperBefore,
			  numberOldActiveCuts_,numberNewCuts_) ;
    if (newNode->numberUnsatisfied()) {
      maximumDepthActual_ = CoinMax(maximumDepthActual_,newNode->depth());
      newNode->initializeInfo() ;
      if (cuts.sizeRowCuts()) {
	lockThread();
	newNode->nodeInfo()->addCuts(cuts,newNode->numberBranches(),
				     whichGenerator_) ;
	unlockThread();
      }
    }
  } else {
    anyAction = -2 ; 
    // Reset bound anyway (no harm if not odd)
    solverCharacteristics_->setMipBound(-COIN_DBL_MAX);
  }
  // May have slipped through i.e. anyAction == 0 and objective above cutoff
  // I think this will screw up cut reference counts if executed.
  // We executed addCuts just above. (lh)
  if ( anyAction >=0 ) {
    assert (newNode);
    if (newNode->objectiveValue() >= getCutoff()) {
      anyAction = -2; // say bad after all
      // zap parent nodeInfo
#ifdef COIN_DEVELOP
      printf("zapping3 CbcNodeInfo %x\n",newNode->nodeInfo()->parent());
#endif
      if (newNode->nodeInfo())
	newNode->nodeInfo()->nullParent();
    }
  }
  return anyAction;
}

/* 
   For advanced applications you may wish to modify the behavior of Cbc
   e.g. if the solver is a NLP solver then you may not have an exact
   optimum solution at each step.  Information could be built into
   OsiSolverInterface but this is an alternative so that that interface 
   does not have to be changed.  If something similar is useful to
   enough solvers then it could be migrated.
   You can also pass in by using solver->setAuxiliaryInfo.
   You should do that if solver is odd - if solver is normal simplex
   then use this
*/
void 
CbcModel::passInSolverCharacteristics(OsiBabSolver * solverCharacteristics)
{
  solverCharacteristics_ = solverCharacteristics;
}
/* preProcess problem - replacing solver
   If makeEquality true then <= cliques converted to ==.
   Presolve will be done numberPasses times.
   
   Returns NULL if infeasible
   
   If makeEquality is 1 add slacks to get cliques,
   if 2 add slacks to get sos (but only if looks plausible) and keep sos info
*/
CglPreProcess *
CbcModel::preProcess( int makeEquality, int numberPasses, int tuning)
{
  CglPreProcess * process = new CglPreProcess();
  // Default set of cut generators
  CglProbing generator1;
  generator1.setUsingObjective(true);
  generator1.setMaxPass(3);
  generator1.setMaxProbeRoot(solver_->getNumCols());
  generator1.setMaxElements(100);
  generator1.setMaxLookRoot(50);
  generator1.setRowCuts(3);
  // Add in generators
  process->addCutGenerator(&generator1);
  process->messageHandler()->setLogLevel(this->logLevel());
  /* model may not have created objects
     If none then create
  */
  if (!numberIntegers_||!numberObjects_) {
    this->findIntegers(true,1);
  }
  // Do SOS
  int i;
  int numberSOS2=0;
  for (i=0;i<numberObjects_;i++) {
    CbcSOS * objSOS =
      dynamic_cast <CbcSOS *>(object_[i]) ;
    if (objSOS) {
      int type = objSOS->sosType();
      if (type==2)
	numberSOS2++;
    }
  }
  if (numberSOS2) {
    // SOS
    int numberColumns = solver_->getNumCols();
    char * prohibited = new char[numberColumns];
    memset(prohibited,0,numberColumns);
    for (i=0;i<numberObjects_;i++) {
      CbcSOS * objSOS =
	dynamic_cast <CbcSOS *>(object_[i]) ;
      if (objSOS) {
	int type = objSOS->sosType();
	if (type==2) {
	  int n=objSOS->numberMembers();
	  const int * which = objSOS->members();
	  for (int j=0;j<n;j++) {
	    int iColumn = which[j];
	    prohibited[iColumn]=1;
	  }
	}
      }
    }
    process->passInProhibited(prohibited,numberColumns);
    delete [] prohibited;
  }
  // Tell solver we are not in Branch and Cut
  solver_->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo) ;
  OsiSolverInterface * newSolver = process->preProcessNonDefault(*solver_, makeEquality,
								numberPasses, tuning);
  // Tell solver we are not in Branch and Cut
  solver_->setHintParam(OsiDoInBranchAndCut,false,OsiHintDo) ;
  if (newSolver) {
    int numberOriginalObjects=numberObjects_;
    OsiSolverInterface * originalSolver = solver_;
    solver_=newSolver->clone(); // clone as process owns solver
    // redo sequence
    numberIntegers_=0;
    int numberColumns = solver_->getNumCols();
    int nOrig = originalSolver->getNumCols();
    const int * originalColumns = process->originalColumns();
    // allow for cliques etc
    nOrig = CoinMax(nOrig,originalColumns[numberColumns-1]+1);
    OsiObject ** originalObject = object_;
    // object number or -1
    int * temp = new int[nOrig];
    int iColumn;
    for (iColumn=0;iColumn<nOrig;iColumn++) 
      temp[iColumn]=-1;
    int iObject;
    numberObjects_=0;
    int nNonInt=0;
    for (iObject=0;iObject<numberOriginalObjects;iObject++) {
      iColumn = originalObject[iObject]->columnNumber();
      if (iColumn<0) {
	nNonInt++;
      } else {
	temp[iColumn]=iObject;
      }
    }
    int numberNewIntegers=0;
    int numberOldIntegers=0;
    int numberOldOther=0;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      int jColumn = originalColumns[iColumn];
      if (temp[jColumn]>=0) {
	int iObject= temp[jColumn];
	CbcSimpleInteger * obj =
	  dynamic_cast <CbcSimpleInteger *>(originalObject[iObject]) ;
	if (obj) 
	  numberOldIntegers++;
	else
	  numberOldOther++;
      } else if (isInteger(iColumn)) {
	numberNewIntegers++;
      }
    }
    /*
      Allocate an array to hold the indices of the integer variables.
      Make a large enough array for all objects
    */
    numberObjects_= numberNewIntegers+numberOldIntegers+numberOldOther+nNonInt;
    object_ = new OsiObject * [numberObjects_];
    delete [] integerVariable_;
    integerVariable_ = new int [numberNewIntegers+numberOldIntegers];
    /*
      Walk the variables again, filling in the indices and creating objects for
      the integer variables. Initially, the objects hold the index and upper &
      lower bounds.
    */
    numberIntegers_=0;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      int jColumn = originalColumns[iColumn];
      if (temp[jColumn]>=0) {
	int iObject= temp[jColumn];
	CbcSimpleInteger * obj =
	  dynamic_cast <CbcSimpleInteger *>(originalObject[iObject]) ;
	if (obj) {
	  object_[numberIntegers_] = originalObject[iObject]->clone();
	  // redo ids etc
	  object_[numberIntegers_]->resetSequenceEtc(numberColumns,originalColumns);
	  integerVariable_[numberIntegers_++]=iColumn;
	}
      } else if (isInteger(iColumn)) {
	object_[numberIntegers_] =
	  new CbcSimpleInteger(this,iColumn);
	integerVariable_[numberIntegers_++]=iColumn;
      }
    }
    numberObjects_=numberIntegers_;
    // Now append other column stuff
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      int jColumn = originalColumns[iColumn];
      if (temp[jColumn]>=0) {
	int iObject= temp[jColumn];
	CbcSimpleInteger * obj =
	  dynamic_cast <CbcSimpleInteger *>(originalObject[iObject]) ;
	if (!obj) {
	  object_[numberObjects_] = originalObject[iObject]->clone();
	  // redo ids etc
	  CbcObject * obj =
	    dynamic_cast <CbcObject *>(object_[numberObjects_]) ;
	  assert (obj);
	  obj->redoSequenceEtc(this,numberColumns,originalColumns);
	  numberObjects_++;
	}
      }
    }
    // now append non column stuff
    for (iObject=0;iObject<numberOriginalObjects;iObject++) {
      iColumn = originalObject[iObject]->columnNumber();
      if (iColumn<0) {
	object_[numberObjects_] = originalObject[iObject]->clone();
	// redo ids etc
	CbcObject * obj =
	  dynamic_cast <CbcObject *>(object_[numberObjects_]) ;
	assert (obj);
	obj->redoSequenceEtc(this,numberColumns,originalColumns);
	numberObjects_++;
      }
      delete originalObject[iObject];
    }
    delete [] originalObject;
    delete [] temp;
    if (!numberObjects_)
      handler_->message(CBC_NOINT,messages_) << CoinMessageEol ;
    return process;
  } else {
    // infeasible
    delete process;
    return NULL;
  }
    
}
// Generate an OsiBranchingInformation object
OsiBranchingInformation 
CbcModel::usefulInformation() const
{
  OsiBranchingInformation usefulInfo(solver_,normalSolver(),false);
  // and modify
  usefulInfo.solution_=testSolution_;
  usefulInfo.integerTolerance_= dblParam_[CbcIntegerTolerance] ;
  usefulInfo.hotstartSolution_=hotstartSolution_;
  usefulInfo.numberSolutions_=numberSolutions_;
  usefulInfo.numberBranchingSolutions_=numberSolutions_-numberHeuristicSolutions_;
  usefulInfo.depth_=-1;
  return usefulInfo;
}
/* Does postprocessing - original solver back.
   User has to delete process */
void 
CbcModel::postProcess(CglPreProcess * process)
{
  process->postProcess(*solver_);
  delete solver_;
  solver_ = process->originalModel();
}
void 
CbcModel::setBestSolution(const double * solution,int numberColumns,double objectiveValue)
{
  bestObjective_ = objectiveValue;
  int n = CoinMax(numberColumns,solver_->getNumCols());
  delete [] bestSolution_;
  bestSolution_ = new double [n];
  memset(bestSolution_,0,n*sizeof(double));
  memcpy(bestSolution_,solution,numberColumns*sizeof(double));
}
// Create C++ lines to get to current state
void 
CbcModel::generateCpp( FILE * fp,int options)
{
  // Do cut generators
  int i;
  for (i=0;i<numberCutGenerators_;i++) {
    CglCutGenerator * generator = generator_[i]->generator();
    std::string name = generator->generateCpp(fp);
    int howOften = generator_[i]->howOften();
    int howOftenInSub = generator_[i]->howOftenInSub();
    int whatDepth = generator_[i]->whatDepth();
    int whatDepthInSub = generator_[i]->whatDepthInSub();
    bool normal = generator_[i]->normal();
    bool atSolution = generator_[i]->atSolution();
    bool whenInfeasible = generator_[i]->whenInfeasible();
    bool timing = generator_[i]->timing();
    fprintf(fp,"3  cbcModel->addCutGenerator(&%s,%d,",
	    name.c_str(),howOften);
    // change name
    name[0]=toupper(name[0]);
    fprintf(fp,"\"%s\",%s,%s,%s,%d,%d,%d);\n",
	    name.c_str(),normal ? "true" : "false",
	    atSolution ? "true" : "false",
	    whenInfeasible ? "true" : "false",
	    howOftenInSub,whatDepth,whatDepthInSub);
    fprintf(fp,"3  cbcModel->cutGenerator(%d)->setTiming(%s);\n",
	    i,timing ? "true" : "false");
    fprintf(fp,"3  \n");
  }
  for (i=0;i<numberHeuristics_;i++) {
    CbcHeuristic * heuristic = heuristic_[i];
    heuristic->generateCpp(fp);
    fprintf(fp,"3  \n");
  }
  if (nodeCompare_)
    nodeCompare_->generateCpp(fp);
  tree_->generateCpp(fp);
  CbcModel defaultModel;
  CbcModel * other = &defaultModel;
  int iValue1, iValue2;
  double dValue1, dValue2;
  iValue1 = this->getMaximumNodes();
  iValue2 = other->getMaximumNodes();
  fprintf(fp,"%d  int save_getMaximumNodes = cbcModel->getMaximumNodes();\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setMaximumNodes(%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  cbcModel->setMaximumNodes(save_getMaximumNodes);\n",iValue1==iValue2 ? 7 : 6);
  iValue1 = this->getMaximumSolutions();
  iValue2 = other->getMaximumSolutions();
  fprintf(fp,"%d  int save_getMaximumSolutions = cbcModel->getMaximumSolutions();\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setMaximumSolutions(%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  cbcModel->setMaximumSolutions(save_getMaximumSolutions);\n",iValue1==iValue2 ? 7 : 6);
  iValue1 = this->numberStrong();
  iValue2 = other->numberStrong();
  fprintf(fp,"%d  int save_numberStrong = cbcModel->numberStrong();\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setNumberStrong(%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  cbcModel->setNumberStrong(save_numberStrong);\n",iValue1==iValue2 ? 7 : 6);
  iValue1 = this->numberBeforeTrust();
  iValue2 = other->numberBeforeTrust();
  fprintf(fp,"%d  int save_numberBeforeTrust = cbcModel->numberBeforeTrust();\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setNumberBeforeTrust(%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  cbcModel->setNumberBeforeTrust(save_numberBeforeTrust);\n",iValue1==iValue2 ? 7 : 6);
  iValue1 = this->numberPenalties();
  iValue2 = other->numberPenalties();
  fprintf(fp,"%d  int save_numberPenalties = cbcModel->numberPenalties();\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setNumberPenalties(%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  cbcModel->setNumberPenalties(save_numberPenalties);\n",iValue1==iValue2 ? 7 : 6);
  iValue1 = this->howOftenGlobalScan();
  iValue2 = other->howOftenGlobalScan();
  fprintf(fp,"%d  int save_howOftenGlobalScan = cbcModel->howOftenGlobalScan();\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setHowOftenGlobalScan(%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  cbcModel->setHowOftenGlobalScan(save_howOftenGlobalScan);\n",iValue1==iValue2 ? 7 : 6);
  iValue1 = this->printFrequency();
  iValue2 = other->printFrequency();
  fprintf(fp,"%d  int save_printFrequency = cbcModel->printFrequency();\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setPrintFrequency(%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  cbcModel->setPrintFrequency(save_printFrequency);\n",iValue1==iValue2 ? 7 : 6);
  iValue1 = this->getPrintingMode();
  iValue2 = other->getPrintingMode();
  fprintf(fp,"%d  int save_printingMode = cbcModel->getPrintingMode();\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setPrintingMode(%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  cbcModel->setPrintingMode(save_printingMode);\n",iValue1==iValue2 ? 7 : 6);
  iValue1 = this->searchStrategy();
  iValue2 = other->searchStrategy();
  fprintf(fp,"%d  int save_searchStrategy = cbcModel->searchStrategy();\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setSearchStrategy(%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  cbcModel->setSearchStrategy(save_searchStrategy);\n",iValue1==iValue2 ? 7 : 6);
  iValue1 = this->specialOptions();
  iValue2 = other->specialOptions();
  fprintf(fp,"%d  int save_cbcSpecialOptions = cbcModel->specialOptions();\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setSpecialOptions(%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  cbcModel->setSpecialOptions(save_cbcSpecialOptions);\n",iValue1==iValue2 ? 7 : 6);
  iValue1 = this->messageHandler()->logLevel();
  iValue2 = other->messageHandler()->logLevel();
  fprintf(fp,"%d  int save_cbcMessageLevel = cbcModel->messageHandler()->logLevel();\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->messageHandler()->setLogLevel(%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  cbcModel->messageHandler()->setLogLevel(save_cbcMessageLevel);\n",iValue1==iValue2 ? 7 : 6);
  iValue1 = this->getMaximumCutPassesAtRoot();
  iValue2 = other->getMaximumCutPassesAtRoot();
  fprintf(fp,"%d  int save_getMaximumCutPassesAtRoot = cbcModel->getMaximumCutPassesAtRoot();\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setMaximumCutPassesAtRoot(%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  cbcModel->setMaximumCutPassesAtRoot(save_getMaximumCutPassesAtRoot);\n",iValue1==iValue2 ? 7 : 6);
  iValue1 = this->getMaximumCutPasses();
  iValue2 = other->getMaximumCutPasses();
  fprintf(fp,"%d  int save_getMaximumCutPasses = cbcModel->getMaximumCutPasses();\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setMaximumCutPasses(%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  cbcModel->setMaximumCutPasses(save_getMaximumCutPasses);\n",iValue1==iValue2 ? 7 : 6);
  iValue1 = this->getPreferredWay();
  iValue2 = other->getPreferredWay();
  fprintf(fp,"%d  int save_getPreferredWay = cbcModel->getPreferredWay();\n",iValue1==iValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setPreferredWay(%d);\n",iValue1==iValue2 ? 4 : 3,iValue1);
  fprintf(fp,"%d  cbcModel->setPreferredWay(save_getPreferredWay);\n",iValue1==iValue2 ? 7 : 6);
  dValue1 = this->getMinimumDrop();
  dValue2 = other->getMinimumDrop();
  fprintf(fp,"%d  double save_getMinimumDrop = cbcModel->getMinimumDrop();\n",dValue1==dValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setMinimumDrop(%g);\n",dValue1==dValue2 ? 4 : 3,dValue1);
  fprintf(fp,"%d  cbcModel->setMinimumDrop(save_getMinimumDrop);\n",dValue1==dValue2 ? 7 : 6);
  dValue1 = this->getIntegerTolerance();
  dValue2 = other->getIntegerTolerance();
  fprintf(fp,"%d  double save_getIntegerTolerance = cbcModel->getIntegerTolerance();\n",dValue1==dValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setIntegerTolerance(%g);\n",dValue1==dValue2 ? 4 : 3,dValue1);
  fprintf(fp,"%d  cbcModel->setIntegerTolerance(save_getIntegerTolerance);\n",dValue1==dValue2 ? 7 : 6);
  dValue1 = this->getInfeasibilityWeight();
  dValue2 = other->getInfeasibilityWeight();
  fprintf(fp,"%d  double save_getInfeasibilityWeight = cbcModel->getInfeasibilityWeight();\n",dValue1==dValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setInfeasibilityWeight(%g);\n",dValue1==dValue2 ? 4 : 3,dValue1);
  fprintf(fp,"%d  cbcModel->setInfeasibilityWeight(save_getInfeasibilityWeight);\n",dValue1==dValue2 ? 7 : 6);
  dValue1 = this->getCutoffIncrement();
  dValue2 = other->getCutoffIncrement();
  fprintf(fp,"%d  double save_getCutoffIncrement = cbcModel->getCutoffIncrement();\n",dValue1==dValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setCutoffIncrement(%g);\n",dValue1==dValue2 ? 4 : 3,dValue1);
  fprintf(fp,"%d  cbcModel->setCutoffIncrement(save_getCutoffIncrement);\n",dValue1==dValue2 ? 7 : 6);
  dValue1 = this->getAllowableGap();
  dValue2 = other->getAllowableGap();
  fprintf(fp,"%d  double save_getAllowableGap = cbcModel->getAllowableGap();\n",dValue1==dValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setAllowableGap(%g);\n",dValue1==dValue2 ? 4 : 3,dValue1);
  fprintf(fp,"%d  cbcModel->setAllowableGap(save_getAllowableGap);\n",dValue1==dValue2 ? 7 : 6);
  dValue1 = this->getAllowableFractionGap();
  dValue2 = other->getAllowableFractionGap();
  fprintf(fp,"%d  double save_getAllowableFractionGap = cbcModel->getAllowableFractionGap();\n",dValue1==dValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setAllowableFractionGap(%g);\n",dValue1==dValue2 ? 4 : 3,dValue1);
  fprintf(fp,"%d  cbcModel->setAllowableFractionGap(save_getAllowableFractionGap);\n",dValue1==dValue2 ? 7 : 6);
  dValue1 = this->getMaximumSeconds();
  dValue2 = other->getMaximumSeconds();
  fprintf(fp,"%d  double save_cbcMaximumSeconds = cbcModel->getMaximumSeconds();\n",dValue1==dValue2 ? 2 : 1);
  fprintf(fp,"%d  cbcModel->setMaximumSeconds(%g);\n",dValue1==dValue2 ? 4 : 3,dValue1);
  fprintf(fp,"%d  cbcModel->setMaximumSeconds(save_cbcMaximumSeconds);\n",dValue1==dValue2 ? 7 : 6);
}
// So we can use osiObject or CbcObject during transition
void getIntegerInformation(const OsiObject * object, double & originalLower,
			   double & originalUpper) 
{
  const CbcSimpleInteger * integerObject = 
    dynamic_cast<const  CbcSimpleInteger *> (object);
  if (integerObject) {
    // get original bounds
    originalLower = integerObject->originalLowerBound();
    originalUpper = integerObject->originalUpperBound();
  } else {
    const OsiSimpleInteger * integerObject = 
      dynamic_cast<const  OsiSimpleInteger *> (object);
    assert (integerObject);
    // get original bounds
    originalLower = integerObject->originalLowerBound();
    originalUpper = integerObject->originalUpperBound();
  }
}
// Set the cut modifier method
void 
CbcModel::setCutModifier(CbcCutModifier * modifier)
{ 
  delete cutModifier_; 
  cutModifier_ = modifier->clone();
}
/* Set the cut modifier method
 */
void 
CbcModel::setCutModifier(CbcCutModifier & modifier)
{ 
  delete cutModifier_;
  cutModifier_ = modifier.clone();
}
#ifdef CBC_THREAD
/* Do one node - broken out for clarity?
   also for parallel (when baseModel!=this)
   Returns 1 if solution found
   node NULL on return if no branches left
   newNode NULL if no new node created
*/
int 
CbcModel::doOneNode(CbcModel * baseModel, CbcNode * & node, CbcNode * & newNode)
{
  int foundSolution=0;
  int currentNumberCuts = 0 ;
  currentNode_=node; // so can be accessed elsewhere
  double bestObjective = bestObjective_;
#ifdef CBC_DEBUG
  printf("%d unsat, way %d, obj %g est %g\n",
	 node->numberUnsatisfied(),node->way(),node->objectiveValue(),
	 node->guessedObjectiveValue());
#endif
#if 0
  if (node->objectiveValue()>getCutoff()) {
    printf("cutoff\n");
    assert (node->nodeInfo());
    node->nodeInfo()->decrement(node->nodeInfo()->numberBranchesLeft()) ;
    delete node;
    node=NULL;
    newNode=NULL;
    return 0;
  }
#endif
#if NEW_UPDATE_OBJECT==0
  // Save clone in branching decision
  if(branchingMethod_)
    branchingMethod_->saveBranchingObject(node->modifiableBranchingObject());
#elif NEW_UPDATE_OBJECT>1
  numberUpdateItems_=0;
#endif
  // Say not on optimal path
  bool onOptimalPath=false;
#   ifdef CHECK_NODE
  printf("Node %x popped from tree - %d left, %d count\n",node,
	 node->nodeInfo()->numberBranchesLeft(),
	 node->nodeInfo()->numberPointingToThis()) ;
  printf("\tdepth = %d, z =  %g, unsat = %d, var = %d.\n",
	 node->depth(),node->objectiveValue(),
	 node->numberUnsatisfied(),
	 node->columnNumber()) ;
#   endif
  
  /*
    Rebuild the subproblem for this node:	 Call addCuts() to adjust the model
    to recreate the subproblem for this node (set proper variable bounds, add
    cuts, create a basis).  This may result in the problem being fathomed by
    bound or infeasibility. Returns 1 if node is fathomed.
    Execute the current arm of the branch: If the problem survives, save the
    resulting variable bounds and call branch() to modify variable bounds
    according to the current arm of the branching object. If we're processing
    the final arm of the branching object, flag the node for removal from the
    live set.
  */
  /*
    Used to generate bound edits for CbcPartialNodeInfo.
  */
  int numberColumns = getNumCols() ;
  double * lowerBefore = new double [numberColumns] ;
  double * upperBefore = new double [numberColumns] ;
  newNode = NULL ;
  bool feasible=true;
  CoinWarmStartBasis *lastws = new CoinWarmStartBasis();
  lockThread();
  // point to genuine ones
  //int save1 = maximumNumberCuts_;
  int save2 = maximumDepth_;
  //maximumNumberCuts_ = baseModel->maximumNumberCuts_;
  //addedCuts_ = baseModel->addedCuts_;
  maximumDepth_ = baseModel->maximumDepth_;
  walkback_ = baseModel->walkback_;
  int retCode =addCuts(node,lastws,numberFixedNow_>numberFixedAtRoot_);
  //if (save1<maximumNumberCuts_) {
    // increased
    //baseModel->maximumNumberCuts_ = maximumNumberCuts_;
    //baseModel->addedCuts_ = addedCuts_;
  //}
  if (save2<maximumDepth_) {
    // increased
    baseModel->maximumDepth_ = maximumDepth_;
    baseModel->walkback_ = walkback_;
  }
  unlockThread();
  int branchesLeft=0;
  if (!retCode) {
    int i ;
    const double * lower = getColLower() ;
    const double * upper = getColUpper() ;
    for (i = 0 ; i < numberColumns ; i++) {
      lowerBefore[i]= lower[i] ;
      upperBefore[i]= upper[i] ;
    }
    if ((solverCharacteristics_->extraCharacteristics()&2)!=0) {
      solverCharacteristics_->setBeforeLower(lowerBefore);
      solverCharacteristics_->setBeforeUpper(upperBefore);
    }
    lockThread();
    if (messageHandler()->logLevel()>2)
      node->modifiableBranchingObject()->print();
    if (branchingMethod_&&branchingMethod_->chooseMethod()) {
      branchesLeft = node->branch(solver_); // new way
    } else {
      // old way so need to cheat
      OsiBranchingObject * branch2 = node->modifiableBranchingObject();
      CbcBranchingObject * branch = dynamic_cast <CbcBranchingObject *>(branch2) ;
      assert (branch);
      branch->setModel(this);
      branchesLeft = node->branch(NULL); // old way
      branch->setModel(baseModel);
    }
    if (mutex_) {
      assert (node->nodeInfo());
      node->nodeInfo()->increment() ;
    }
    unlockThread();
    if ((specialOptions_&1)!=0) {
      /*
	This doesn't work as intended --- getRowCutDebugger will return null
	unless the current feasible solution region includes the optimal solution
	that RowCutDebugger knows. There's no way to tell inactive from off the
	optimal path.
      */
      const OsiRowCutDebugger *debugger = solver_->getRowCutDebugger() ;
      if (debugger) {
	onOptimalPath=true;
	printf("On optimal path\n") ;
      }
    }
      
    /*
      Reoptimize, possibly generating cuts and/or using heuristics to find
      solutions.  Cut reference counts are unaffected unless we lose feasibility,
      in which case solveWithCuts() will make the adjustment.
    */
    phase_=2;
    OsiCuts cuts ;
    currentNumberCuts = solver_->getNumRows()-numberRowsAtContinuous_ ;
    int saveNumber = numberIterations_;
    if(solverCharacteristics_->solutionAddsCuts()) {
      int returnCode=resolve(node ? node->nodeInfo() : NULL,1);
      feasible = returnCode != 0;
      if (feasible) {
	int iObject ;
	int preferredWay ;
	int numberUnsatisfied = 0 ;
	memcpy(currentSolution_,solver_->getColSolution(),
	       numberColumns*sizeof(double)) ;
	// point to useful information
	OsiBranchingInformation usefulInfo=usefulInformation();
	
	for (iObject = 0 ; iObject < numberObjects_ ; iObject++) {
	  double infeasibility =
	    object_[iObject]->infeasibility(&usefulInfo,preferredWay) ;
	  if (infeasibility ) numberUnsatisfied++ ;
	}
	if (returnCode>0) {
	  if (numberUnsatisfied)   {
	    feasible = solveWithCuts(cuts,maximumCutPasses_,node);
	  } else {
	    // may generate cuts and turn the solution
	    //to an infeasible one
	    feasible = solveWithCuts(cuts, 1,
				     node);
	  }
	}
	// check extra info on feasibility
	if (!solverCharacteristics_->mipFeasible()) {
	  feasible = false;
	  solverCharacteristics_->setMipBound(-COIN_DBL_MAX);
	}
      }
    } else {
      // normal
      feasible = solveWithCuts(cuts,maximumCutPasses_,node);
    }
    if ((specialOptions_&1)!=0&&onOptimalPath) {
      if (!solver_->getRowCutDebugger()||!feasible) {
	// dj fix did something???
	solver_->writeMps("infeas2");
	assert (solver_->getRowCutDebugger()) ;
	assert (feasible);
      }
    }
    if (statistics_) {
      assert (numberNodes2_);
      assert (statistics_[numberNodes2_-1]);
      assert (statistics_[numberNodes2_-1]->node()==numberNodes2_-1);
      statistics_[numberNodes2_-1]->endOfBranch(numberIterations_-saveNumber,
						feasible ? solver_->getObjValue()
						: COIN_DBL_MAX);
    }
    /*
      Are we still feasible? If so, create a node and do the work to attach a
      branching object, reoptimising as needed if chooseBranch() identifies
      monotone objects.
      
      Finally, attach a partial nodeInfo object and store away any cuts that we
      created back in solveWithCuts. addCuts() will initialise the reference
      counts for these new cuts.
      
      This next test can be problematic if we've discovered an
      alternate equivalent answer and subsequently fathom the solution
      known to the row cut debugger due to bounds.
    */
    if (onOptimalPath) {
      bool objLim = solver_->isDualObjectiveLimitReached() ;
      if (!feasible && !objLim) {
	printf("infeas2\n");
	solver_->writeMps("infeas");
	CoinWarmStartBasis *slack =
	  dynamic_cast<CoinWarmStartBasis *>(solver_->getEmptyWarmStart()) ;
	solver_->setWarmStart(slack);
	delete slack ;
	solver_->setHintParam(OsiDoReducePrint,false,OsiHintDo,0) ;
	solver_->initialSolve();
	assert (!solver_->isProvenOptimal());
      }
      assert (feasible || objLim);
    }
    bool checkingNode=false;
    if (feasible) {
      newNode = new CbcNode ;
      // Set objective value (not so obvious if NLP etc)
      setObjectiveValue(newNode,node);
      int anyAction =-1 ;
      bool resolved = false ;
      if (newNode->objectiveValue() >= getCutoff()) 
	anyAction=-2;
      // only allow at most a few passes
      int numberPassesLeft=5;
      checkingNode=true;
      OsiSolverBranch * branches=NULL;
      // point to useful information
      anyAction = chooseBranch(newNode, numberPassesLeft,node, cuts,resolved,
			       lastws, lowerBefore, upperBefore, branches);
      /*
	If we end up infeasible, we can delete the new node immediately. Since this
	node won't be needing the cuts we collected, decrement the reference counts.
	If we are feasible, then we'll be placing this node into the live set, so
	increment the reference count in the current (parent) nodeInfo.
      */
      lockThread();
      if (anyAction == -2) {
	if (mutex_) {
	  assert (node->nodeInfo());
	  node->nodeInfo()->decrement() ;
	}
	delete newNode ;
	if (mutex_) {
	  assert (node->nodeInfo());
	  node->nodeInfo()->increment() ;
	}
	newNode = NULL ;
	// say strong doing well
	if (checkingNode)
	  setSpecialOptions(specialOptions_|8);
	for (i = 0 ; i < currentNumberCuts_ ; i++) {
	  if (addedCuts_[i]) {
	    if (!addedCuts_[i]->decrement(1)) {
	      delete addedCuts_[i] ;
	    }
	      addedCuts_[i]=NULL;
	      //}
	  }
	}
      }	else { 
	assert (node->nodeInfo());
	node->nodeInfo()->increment() ;
	if ((numberNodes_%20)==0) {
	  // say strong not doing as well
	  setSpecialOptions(specialOptions_&~8);
	}
      }
      unlockThread();
    }
    /*
      At this point, there are three possibilities:
      * newNode is live and will require further branching to resolve
      (variable() >= 0). Increment the cut reference counts by
      numberBranches() to allow for use by children of this node, and
      decrement by 1 because we've executed one arm of the branch of our
      parent (consuming one reference). Before we push newNode onto the
      search tree, try for a heuristic solution.
      * We have a solution, in which case newNode is non-null but we have no
      branching variable. Decrement the cut counts and save the solution.
      * The node was found to be infeasible, in which case it's already been
      deleted, and newNode is null.
    */
    if (!getEventHandler()->event(CbcEventHandler::node)) {
      eventHappened_=true; // exit
    }
    assert (!newNode || newNode->objectiveValue() <= getCutoff()) ;
    if (statistics_) {
      assert (numberNodes2_);
      assert (statistics_[numberNodes2_-1]);
      assert (statistics_[numberNodes2_-1]->node()==numberNodes2_-1);
      if (newNode)
	statistics_[numberNodes2_-1]->updateInfeasibility(newNode->numberUnsatisfied());
      else
	statistics_[numberNodes2_-1]->sayInfeasible();
    }
    lockThread();
#if NEW_UPDATE_OBJECT>1
      if (numberUpdateItems_&&!numberThreads_) {
	for (i=0;i<numberUpdateItems_;i++) {
	  CbcObjectUpdateData * update = updateItems_+i;
	  CbcObject * object = dynamic_cast<CbcObject *> (update->object_);
	  if (object) 
	    object->updateInformation(*update);
	}
	numberUpdateItems_=0;
      }
#endif
    if (newNode) {
      if (newNode->branchingObject() == NULL&&solverCharacteristics_->solverType()==4) {
	// need to check if any cuts would do anything
	OsiCuts theseCuts;
	// reset probing info
	if (probingInfo_)
	  probingInfo_->initializeFixing();
	for (int i = 0;i<numberCutGenerators_;i++) {
	  bool generate = generator_[i]->normal();
	  // skip if not optimal and should be (maybe a cut generator has fixed variables)
	  if (generator_[i]->needsOptimalBasis()&&!solver_->basisIsAvailable())
	    generate=false;
	  if (!generator_[i]->mustCallAgain())
	    generate=false; // only special cuts
	  if (generate) {
	    generator_[i]->generateCuts(theseCuts,true,solver_,NULL) ;
	    int numberRowCutsAfter = theseCuts.sizeRowCuts() ;
	    if (numberRowCutsAfter) {
	      // need dummy branch
	      newNode->setBranchingObject(new CbcDummyBranchingObject(this));
	      newNode->nodeInfo()->initializeInfo(1);
	      break;
	    }
	  }
	}
      }
      if (newNode->branchingObject()) {
	handler_->message(CBC_BRANCH,messages_)
	  << numberNodes_<< newNode->objectiveValue()
	  << newNode->numberUnsatisfied()<< newNode->depth()
	  << CoinMessageEol ;
	// Increment cut counts (taking off current)
	int numberLeft = newNode->numberBranches() ;
	for (i = 0;i < currentNumberCuts_;i++) {
	  if (addedCuts_[i]) {
#		ifdef CHECK_CUT_COUNTS
	    printf("Count on cut %x increased by %d\n",addedCuts_[i],
		   numberLeft-1) ;
#		endif
	    addedCuts_[i]->increment(numberLeft-1) ;
	  }
	}
	unlockThread();

	double estValue = newNode->guessedObjectiveValue() ;
	int found = -1 ;
	double * newSolution = new double [numberColumns] ;
	double heurValue = getCutoff() ;
	int iHeur ;
	for (iHeur = 0 ; iHeur < numberHeuristics_ ; iHeur++) {
	  double saveValue = heurValue ;
	  int ifSol = heuristic_[iHeur]->solution(heurValue,newSolution) ;
	  if (ifSol > 0) {
	    // new solution found
	    found = iHeur ;
	    lockThread();
	    baseModel->incrementUsed(newSolution);
	    unlockThread();
	  } else if (ifSol < 0)	{ // just returning an estimate 
	    estValue = CoinMin(heurValue,estValue) ;
	    heurValue = saveValue ;
	  }
	}
	if (found >= 0) {
	  lastHeuristic_ = heuristic_[found];
	  setBestSolution(CBC_ROUNDING,heurValue,newSolution) ;
	  foundSolution=1;
	}
	delete [] newSolution ;
	newNode->setGuessedObjectiveValue(estValue) ;
	lockThread();
#define PUSH_LATER
#ifdef PUSH_LATER
	if (!mutex_) // only if serial
#endif
	tree_->push(newNode) ;
	if (statistics_) {
	  if (numberNodes2_==maximumStatistics_) {
	    maximumStatistics_ = 2*maximumStatistics_;
	    CbcStatistics ** temp = new CbcStatistics * [maximumStatistics_];
	    memset(temp,0,maximumStatistics_*sizeof(CbcStatistics *));
	    memcpy(temp,statistics_,numberNodes2_*sizeof(CbcStatistics *));
	    delete [] statistics_;
	    statistics_=temp;
	  }
	  assert (!statistics_[numberNodes2_]);
	  statistics_[numberNodes2_]=new CbcStatistics(newNode);
	}
	numberNodes2_++;
#	    ifdef CHECK_NODE
	printf("Node %x pushed on tree c\n",newNode) ;
#	    endif
      } else { 
	if(solverCharacteristics_ && //we may be in a non standard bab
	   solverCharacteristics_->solutionAddsCuts()// we are in some kind of OA based bab.
	   ) {
	  
	  std::cerr<<"You should never get here"<<std::endl;
	  throw CoinError("Nodes should not be fathomed on integer infeasibility in this setting",
			  "branchAndBound","CbcModel") ;
	}
	for (i = 0 ; i < currentNumberCuts_ ; i++) { 
	  if (addedCuts_[i]) { 
	    if (!addedCuts_[i]->decrement(1)) {
	      delete addedCuts_[i] ;
	      addedCuts_[i]=NULL;
	    }
	  }
	}
	double objectiveValue = newNode->objectiveValue();
	lastHeuristic_ = NULL;
	// Just possible solver did not know about a solution from another thread!
	if (objectiveValue<getCutoff()) {
	  incrementUsed(solver_->getColSolution());
	  setBestSolution(CBC_SOLUTION,objectiveValue,
			  solver_->getColSolution()) ;
	  foundSolution=1;
	}
	//assert(nodeInfo->numberPointingToThis() <= 2) ;
	// avoid accidental pruning, if newNode was final branch arm
	node->nodeInfo()->increment();
	delete newNode ;
	newNode=NULL;
	node->nodeInfo()->decrement() ;
      }
    }
    if (branchesLeft) {
      // set nodenumber correctly
      if (node->nodeInfo())
	node->nodeInfo()->setNodeNumber(numberNodes2_);
#ifdef PUSH_LATER
	if (!mutex_) // only if serial
#endif
      tree_->push(node) ;
      if (statistics_) {
	if (numberNodes2_==maximumStatistics_) {
	  maximumStatistics_ = 2*maximumStatistics_;
	  CbcStatistics ** temp = new CbcStatistics * [maximumStatistics_];
	  memset(temp,0,maximumStatistics_*sizeof(CbcStatistics *));
	  memcpy(temp,statistics_,numberNodes2_*sizeof(CbcStatistics *));
	  delete [] statistics_;
	  statistics_=temp;
	}
	assert (!statistics_[numberNodes2_]);
	statistics_[numberNodes2_]=new CbcStatistics(node);
      }
      numberNodes2_++;
      //nodeOnTree=true; // back on tree
      //deleteNode = false ;
#	ifdef CHECK_NODE
      printf("Node %x pushed back on tree - %d left, %d count\n",node,
	     node->nodeInfo()->numberBranchesLeft(),
	     node->nodeInfo()->numberPointingToThis()) ;
#	endif
      if (mutex_) {
	assert (node->nodeInfo());
	node->nodeInfo()->decrement() ;
      }
    } else {
    /*
      This node has been completely expanded and can be removed from the live
      set.
    */
      if (mutex_) {
	assert (node->nodeInfo());
	node->nodeInfo()->decrement() ;
      }
      assert (node->nodeInfo());
      if (!node->nodeInfo()->numberBranchesLeft())
	node->nodeInfo()->allBranchesGone(); // can clean up
      delete node ;
      node=NULL;
    }
    unlockThread();
  } else {
    // add cuts found to be infeasible (on bound)!
    printf("found to be infeas! - branches left %d\n",node->nodeInfo()->numberBranchesLeft());
    //abort();
    assert (node->nodeInfo());
    if (!node->nodeInfo()->numberBranchesLeft())
      node->nodeInfo()->allBranchesGone(); // can clean up
    delete node;
    node=NULL;
  }
  /*
    Delete cuts to get back to the original system.
    
    I'm thinking this is redundant --- the call to addCuts that conditions entry
    to this code block also performs this action.
  */
  int numberToDelete = getNumRows()-numberRowsAtContinuous_ ;
  if (numberToDelete) {
    int * delRows = new int[numberToDelete] ;
    int i ;
    for (i = 0 ; i < numberToDelete ; i++) 
      delRows[i] = i+numberRowsAtContinuous_ ; 
    solver_->deleteRows(numberToDelete,delRows) ;
    delete [] delRows ;
  }
  delete lastws ;
  delete [] lowerBefore ;
  delete [] upperBefore ;
  if (bestObjective > bestObjective_)
    foundSolution=2;
  if (foundSolution) {
    lockThread();
    // might as well mark all including continuous
    int numberColumns = solver_->getNumCols();
    for (int i=0;i<numberColumns;i++) {
      baseModel->usedInSolution_[i] += usedInSolution_[i];
      usedInSolution_[i]=0;
    }
    baseModel->numberSolutions_ = numberSolutions_;
    if (bestObjective_ < baseModel->bestObjective_&&bestObjective_<baseModel->getCutoff()) {
      baseModel->bestObjective_ = bestObjective_ ;
      int numberColumns = solver_->getNumCols();
      if (!baseModel->bestSolution_)
        baseModel->bestSolution_ = new double[numberColumns];
      CoinCopyN(bestSolution_,numberColumns,baseModel->bestSolution_);
      baseModel->setCutoff(getCutoff());
    }
    unlockThread();
  }
  return foundSolution;
}
#if NEW_UPDATE_OBJECT>1
// Adds an update information object
void 
CbcModel::addUpdateInformation(const CbcObjectUpdateData & data)
{
  if (numberUpdateItems_==maximumNumberUpdateItems_) {
    maximumNumberUpdateItems_ += 10;
    CbcObjectUpdateData * temp = new CbcObjectUpdateData [maximumNumberUpdateItems_];
    for (int i=0;i<maximumNumberUpdateItems_-10;i++)
      temp[i] = updateItems_[i];
    delete [] updateItems_;
    updateItems_ = temp;
  }
  updateItems_[numberUpdateItems_++]=data;
}
#endif
/* Move/copy information from one model to another
   -1 - initial setup
   0 - from base model
   1 - to base model (and reset)
   2 - add in final statistics etc (and reset so can do clean destruction)
*/
void 
CbcModel::moveToModel(CbcModel * baseModel,int mode)
{
  if (mode==0) {
    setCutoff(baseModel->getCutoff());
    bestObjective_ = baseModel->bestObjective_;
    assert (!baseModel->globalCuts_.sizeRowCuts());
    numberSolutions_ = baseModel->numberSolutions_;
    stateOfSearch_ = baseModel->stateOfSearch_;
    numberNodes_ = baseModel->numberNodes_;
    numberIterations_ = baseModel->numberIterations_;
    numberFixedAtRoot_ = numberIterations_; // for statistics
    phase_ = baseModel->phase_;
    assert (!nextRowCut_);
    nodeCompare_ = baseModel->nodeCompare_;
    tree_ = baseModel->tree_;
    assert (!subTreeModel_);
    //branchingMethod_ = NULL; // need something but what
    numberOldActiveCuts_ = baseModel->numberOldActiveCuts_;
    cutModifier_ = NULL;
    assert (!analyzeResults_);
    threadStruct * stuff = (threadStruct *) mutex_;
    assert (stuff);
    //if (stuff)
    stuff->createdNode=NULL;
    // ?? searchStrategy_;
    searchStrategy_=baseModel->searchStrategy_;
    stuff->saveStuff[0]=searchStrategy_;
    stateOfSearch_=baseModel->stateOfSearch_;
    stuff->saveStuff[1]=stateOfSearch_;
#if NEW_UPDATE_OBJECT>1
    for (int iObject = 0 ; iObject < numberObjects_ ; iObject++) {
      CbcSimpleIntegerDynamicPseudoCost * dynamicObject =
	dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(object_[iObject]) ;
      if (dynamicObject) {
	CbcSimpleIntegerDynamicPseudoCost * baseObject =
	dynamic_cast <CbcSimpleIntegerDynamicPseudoCost *>(baseModel->object_[iObject]) ;
	assert (baseObject);
	dynamicObject->copySome(baseObject);
      }
    }
#endif
 } else if (mode==1) {
    lockThread();
    threadStruct * stuff = (threadStruct *) mutex_;
    assert (stuff);
    //stateOfSearch_
#if 1
    if(stuff->saveStuff[0]!=searchStrategy_) {
#ifdef COIN_DEVELOP
      printf("changing searchStrategy from %d to %d\n",
 	     baseModel->searchStrategy_,searchStrategy_);
#endif
      baseModel->searchStrategy_=searchStrategy_;
    }
    if(stuff->saveStuff[1]!=stateOfSearch_) {
#ifdef COIN_DEVELOP
      printf("changing stateOfSearch from %d to %d\n",
 	     baseModel->stateOfSearch_,stateOfSearch_);
#endif
      baseModel->stateOfSearch_=stateOfSearch_;
    }
#endif
#if NEW_UPDATE_OBJECT>1
    if (numberUpdateItems_) {
      for (int i=0;i<numberUpdateItems_;i++) {
	CbcObjectUpdateData * update = updateItems_+i;
	int objectNumber = update->objectNumber_;
	CbcObject * object = dynamic_cast<CbcObject *> (baseModel->object_[objectNumber]);
	if (object) 
	  object->updateInformation(*update);
      }
      numberUpdateItems_=0;
    }
#endif
    if (eventHappened_)
      baseModel->eventHappened_=true;
    baseModel->numberNodes_++;
    baseModel->numberIterations_ +=
      numberIterations_ - numberFixedAtRoot_; 
#ifdef PUSH_LATER
    if (stuff->node)
      baseModel->tree_->push(stuff->node);
    if (stuff->createdNode)
      baseModel->tree_->push(stuff->createdNode);
#endif
    unlockThread();
  } else if (mode==2) {
    baseModel->sumChangeObjective1_ += sumChangeObjective1_; 
    baseModel->sumChangeObjective2_ += sumChangeObjective2_; 
    //baseModel->numberIterations_ += numberIterations_;
    for (int iGenerator=0;iGenerator<numberCutGenerators_;iGenerator++) {
      CbcCutGenerator * generator = baseModel->generator_[iGenerator];
      CbcCutGenerator * generator2 = generator_[iGenerator];
      generator->incrementNumberTimesEntered(generator2->numberTimesEntered());
      generator->incrementNumberCutsInTotal(generator2->numberCutsInTotal());
      generator->incrementNumberCutsActive(generator2->numberCutsActive());
      generator->incrementTimeInCutGenerator(generator2->timeInCutGenerator());
    }
    nodeCompare_ = NULL;
    baseModel->maximumDepthActual_ = CoinMax(baseModel->maximumDepthActual_,maximumDepthActual_);
    baseModel->numberDJFixed_ += numberDJFixed_;
    baseModel->numberStrongIterations_ += numberStrongIterations_;
    int i;
    for (i=0;i<3;i++) 
      baseModel->strongInfo_[i] += strongInfo_[i];
    walkback_ = NULL;
    //addedCuts_ = NULL;
    tree_ = NULL;
    eventHandler_=NULL;
    delete solverCharacteristics_;
    solverCharacteristics_ = NULL;
    bool newMethod = (baseModel->branchingMethod_&&baseModel->branchingMethod_->chooseMethod());
    if (newMethod) {
      // new method - we were using base models
      numberObjects_=0;
      object_=NULL;
    }
  } else {
    delete eventHandler_;
    eventHandler_ = baseModel->eventHandler_;
    assert (!statistics_);
    delete [] walkback_;
    //delete [] addedCuts_;
    walkback_ = NULL;
    //addedCuts_ = NULL;
    assert(baseModel->solverCharacteristics_);
    solverCharacteristics_ = new OsiBabSolver (*baseModel->solverCharacteristics_);
    solverCharacteristics_->setSolver(solver_);
    delete tree_;
    tree_ = NULL;
    delete nodeCompare_;
    nodeCompare_ = NULL;
    continuousSolver_ = baseModel->continuousSolver_->clone();
    bool newMethod = (baseModel->branchingMethod_&&baseModel->branchingMethod_->chooseMethod());
    if (newMethod) {
      // new method uses solver - but point to base model
      // We may update an object in wrong order - shouldn't matter?
      numberObjects_=baseModel->numberObjects_;
      object_=baseModel->object_;
    }
    mutex_ = baseModel->mutex_;
    int i;
    for (i=0;i<numberHeuristics_;i++) {
      delete heuristic_[i];
      heuristic_[i] = baseModel->heuristic_[i]->clone();
      heuristic_[i]->setModelOnly(this);
    }
    for (i=0;i<numberCutGenerators_;i++) {
      delete generator_[i];
      generator_[i] = new CbcCutGenerator(*baseModel->generator_[i]);
      // refreshModel was overkill as thought too many rows
      generator_[i]->setModel(this);
    }
  }
}
#ifdef CBC_THREAD
static void * doNodesThread(void * voidInfo)
{
  threadStruct * stuff = (threadStruct *) voidInfo;
  pthread_mutex_t * mutex = stuff->mutex2;
  pthread_cond_t * condition = stuff->condition2;
  CbcModel * thisModel = stuff->thisModel;
  CbcModel * baseModel = stuff->baseModel;
  while (true) {
    pthread_mutex_lock (mutex);
    while (stuff->returnCode) {
      struct timespec absTime2;
      clock_gettime(CLOCK_REALTIME,&absTime2);
      double time2 = absTime2.tv_sec+1.0e-9*absTime2.tv_nsec;
      // timed wait as seems to hang on max nodes at times
      absTime2.tv_sec += 10;
      pthread_cond_timedwait(condition,mutex,&absTime2);
      clock_gettime(CLOCK_REALTIME,&stuff->absTime);
      double time = stuff->absTime.tv_sec+1.0e-9*stuff->absTime.tv_nsec;
      stuff->timeWaitingToStart+=time-time2;;
      stuff->numberTimesWaitingToStart++;
    }
    //printf("start node %x\n",stuff->node);
    int mode = thisModel->getNumberThreads();
    if (mode) {
      // normal
      double time2 = CoinCpuTime();
      assert (stuff->returnCode==0);
      assert (stuff->node->nodeInfo());
      thisModel->doOneNode(baseModel,stuff->node,stuff->createdNode);
      stuff->returnCode=1;
      //printf("end node %x\n",stuff->node);
      threadStruct * stuffMain = (threadStruct *) baseModel->mutex();
      //pthread_mutex_t * condition_mutex = stuffMain->mutex2;
      pthread_cond_t * condition_main = stuffMain->condition2;
      pthread_cond_signal(condition_main); // unlock
      pthread_mutex_unlock(mutex);
      stuff->timeInThread += CoinCpuTime()-time2;
    } else {
      // exit
      break;
    }
  }
  pthread_mutex_unlock(mutex);
  pthread_exit(NULL);
  return NULL;
}
static void * doCutsThread(void * voidInfo)
{
  threadStruct * stuff = (threadStruct *) voidInfo;
  pthread_mutex_t * mutex = stuff->mutex2;
  pthread_cond_t * condition = stuff->condition2;
  CbcModel * thisModel =  stuff->thisModel;
  CbcModel * baseModel = stuff->baseModel;
  while (true) {
    pthread_mutex_lock(mutex);
    while (stuff->returnCode) {
      pthread_cond_wait(condition,mutex);
    }
    //printf("start node %x\n",stuff->node);
    int mode = thisModel->getNumberThreads();
    if (mode) {
      // normal
      assert (stuff->returnCode==0);
      bool fullScan = thisModel->getNodeCount()>0;
      CbcCutGenerator * generator = thisModel->cutGenerator(0);
      OsiCuts * cuts = (OsiCuts *) thisModel->objects();
      OsiSolverInterface * thisSolver = thisModel->solver();
      generator->generateCuts(*cuts,fullScan,thisSolver,NULL);
      stuff->returnCode=1;
      //printf("end node %x\n",stuff->node);
      threadStruct * stuffMain = (threadStruct *) baseModel->mutex();
      //pthread_mutex_t * condition_mutex = stuffMain->mutex2;
      pthread_cond_t * condition_main = stuffMain->condition2;
      pthread_cond_signal(condition_main); // unlock
      pthread_mutex_unlock(mutex);
    } else {
      // exit
      break;
    }
  }
  pthread_mutex_unlock(mutex);
  pthread_exit(NULL);
  return NULL;
}
#endif
#endif
#if CBC_THREAD_DEBUG
static int threadFlag=0;
/*
  -1000 all quiet
  k locked
  k+n*1000 locked for n'th time
  -k unlocked when already unlocked
*/
static int lastAction=-1000;
#if CBC_THREAD_DEBUG>1
#define MAX_DEBUG 10000000
static int numberActions=0;
static int action[MAX_DEBUG];
#endif
#endif
#ifdef CBC_THREAD
/*
  Locks a thread if parallel so that stuff like cut pool
  can be updated and/or used.y
*/
void 
CbcModel::lockThread()
{
  threadStruct * stuff = (threadStruct *) mutex_;
  if (stuff) {
    if(!stuff->locked) {
      struct timespec absTime2;
      clock_gettime(CLOCK_REALTIME,&absTime2);
      double time2 = absTime2.tv_sec+1.0e-9*absTime2.tv_nsec;
      pthread_mutex_lock (stuff->mutex);
      stuff->locked=true;
      clock_gettime(CLOCK_REALTIME,&stuff->absTime);
      double time = stuff->absTime.tv_sec+1.0e-9*stuff->absTime.tv_nsec;
      stuff->timeWaitingToLock+=time-time2;;
      stuff->numberTimesLocked++;
#if CBC_THREAD_DEBUG
      if (lastAction==-1000) {
	lastAction = stuff->threadNumber;
      } else if (lastAction<=0) {
	printf ("thread %d was over unlocked - now locked by %d\n",
		lastAction,stuff->threadNumber);
	lastAction = stuff->threadNumber;
      } else if ((lastAction%1000)==stuff->threadNumber) {
	printf ("thread %d was already locked at level %d\n",
		stuff->threadNumber,lastAction);
	lastAction += 1000;
      } else {
	printf ("thread %d was locked - now lock wanted by %d\n",
		lastAction,stuff->threadNumber);
	abort();
      }
#if CBC_THREAD_DEBUG>1
      assert (numberActions<MAX_DEBUG);
      action[numberActions++] = lastAction;
#endif
      assert (threadFlag==0);
      threadFlag=1;
    } else {
      if (lastAction>0&&(lastAction%1000)==stuff->threadNumber) {
	//printf ("OKish - thread %d was already locked at level %d\n",
	//stuff->threadNumber,lastAction);
	lastAction += 1000;
      } else {
	printf ("lastAction %d - now second lock wanted by %d\n",
		lastAction,stuff->threadNumber);
	abort();
      }
#if CBC_THREAD_DEBUG>1
      assert (numberActions<MAX_DEBUG);
      action[numberActions++] = lastAction;
#endif
      assert (threadFlag==1);
#endif
    }
  }
}
/*
  Unlocks a thread if parallel 
*/
void 
CbcModel::unlockThread()
{
  threadStruct * stuff = (threadStruct *) mutex_;
  if (stuff) {
    if(stuff->locked) {
#if CBC_THREAD_DEBUG
      if (lastAction>0&&(lastAction%1000)==stuff->threadNumber) {
	lastAction = -1000;
      } else {
	printf ("lastAction %d - %d trying to unlock\n",
		lastAction,stuff->threadNumber);
	abort();
      }
#if CBC_THREAD_DEBUG>1
      assert (numberActions<MAX_DEBUG);
      action[numberActions++] = lastAction;
#endif
      assert (threadFlag==1);
      threadFlag=0;
#endif
      stuff->locked=false;
      pthread_mutex_unlock (stuff->mutex);
      struct timespec absTime2;
      clock_gettime(CLOCK_REALTIME,&absTime2);
      double time2 = absTime2.tv_sec+1.0e-9*absTime2.tv_nsec;
      double time = stuff->absTime.tv_sec+1.0e-9*stuff->absTime.tv_nsec;
      stuff->timeLocked+=time2-time;
      stuff->numberTimesUnlocked++;
#if CBC_THREAD_DEBUG
    } else {
      printf ("lastAction %d - %d trying to unlock when flag says unlocked!\n",
	      lastAction,stuff->threadNumber);
      lastAction = -1000;
#if CBC_THREAD_DEBUG>1
      assert (numberActions<MAX_DEBUG);
      action[numberActions++] = lastAction;
#endif
      assert (threadFlag==0);
#endif
    }
  }
}
#endif
// Returns true if locked
bool 
CbcModel::isLocked() const
{
#ifdef CBC_THREAD
  threadStruct * stuff = (threadStruct *) mutex_;
  if (stuff) {
    return (stuff->locked);
  } else {
    return true;
  }
#else
  return true;
#endif
}
