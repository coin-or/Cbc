// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

// Edwin 11/10/2009-- carved out of CbcBranchActual

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>
//#define CBC_DEBUG

#include "CoinTypes.h"
#include "OsiSolverInterface.hpp"
#include "OsiSolverBranch.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcGeneralDepth.hpp"
#include "CbcBranchActual.hpp"
#include "CbcHeuristicDive.hpp"
#include "CoinSort.hpp"
#include "CoinError.hpp"

#include "OsiClpSolverInterface.hpp"
#include "CoinWarmStartBasis.hpp"
#include "ClpNode.hpp"
#include "ClpFactorization.hpp"
#include "CbcBranchDynamic.hpp"
/* I (JJHF) am looking at the fathomMany option in Cbc.  This does a few 
   branches using Clp and tries to be much more efficient as it knows
   about Clp.  It will also lend itself well to deterministc parallel
   and that is what this update to master is all about.  As checked in
   this time it should give identical results as before but now allows
   fathomMany in deterministic parallel.  fathomMany is switched on by 
   -depth +n which varies size of mini tree.  Deterministic parallel
   is switched on by -thread 10k where k is number of threads
   (k > number of cpus can be good).

   I will be trying to tune further.  Ideas which seem possible but
   changed current behavior will be found in areas marked with
   DETERMINISTIC_TUNING e.g.
   // DETERMINISTIC_TUNING think harder or
   #if 0 // DETERMINISTIC_TUNING or just
   #if DETERMINISTIC_TUNING > n etc
*/
// Default Constructor
CbcGeneralDepth::CbcGeneralDepth()
  : CbcGeneral()
  , maximumDepth_(0)
  , maximumNodes_(0)
  , whichSolution_(-1)
  , numberNodes_(0)
  , nodeInfo_(NULL)
{
}

// Useful constructor (which are integer indices)
CbcGeneralDepth::CbcGeneralDepth(CbcModel *model, int maximumDepth)
  : CbcGeneral(model)
  , maximumDepth_(maximumDepth)
  , maximumNodes_(0)
  , whichSolution_(-1)
  , numberNodes_(0)
  , nodeInfo_(NULL)
{
  assert(maximumDepth_ < 10000000);
  if (maximumDepth_ == 0)
    maximumDepth_ = 1;
  // DETERMINISTIC_TUNING think harder
#define MAX_DEPTH 1000
#define MAX_NODES 200
  int realMaximumDepth = maximumDepth_ % MAX_DEPTH;
  if (maximumDepth_ > 0)
    maximumNodes_ = (1 << realMaximumDepth) + 1 + realMaximumDepth;
  else 
    maximumNodes_ = 1 + 1 - maximumDepth_;
  maximumNodes_ = std::min(maximumNodes_, 1 + abs(realMaximumDepth) + MAX_NODES);
  maximumNodes_ = std::max(maximumNodes_, 10);
  // special for 1
  if (maximumDepth_==1) {
    // DETERMINISTIC_TUNING think harder
    // amount will vary
    maximumNodes_ = 1 + 20 + MAX_NODES; 
  }
  {
    nodeInfo_ = new ClpNodeStuff();
    nodeInfo_->maximumNodes_ = maximumNodes_;
    ClpNodeStuff *info = nodeInfo_;
    // for reduced costs and duals
    info->solverOptions_ |= 7;
    if (maximumDepth_ > 0) {
      info->nDepth_ = realMaximumDepth;
    } else {
      info->nDepth_ = -maximumDepth_;
      info->solverOptions_ |= 32;
    }
    ClpNode **nodeInfo = new ClpNode *[maximumNodes_];
    for (int i = 0; i < maximumNodes_; i++)
      nodeInfo[i] = NULL;
    info->nodeInfo_ = nodeInfo;
  }
}

// Copy constructor
CbcGeneralDepth::CbcGeneralDepth(const CbcGeneralDepth &rhs)
  : CbcGeneral(rhs)
{
  maximumDepth_ = rhs.maximumDepth_;
  int realMaximumDepth = maximumDepth_ % MAX_DEPTH;
  maximumNodes_ = rhs.maximumNodes_;
  whichSolution_ = -1;
  numberNodes_ = 0;
  if (maximumNodes_) {
    assert(rhs.nodeInfo_);
    nodeInfo_ = new ClpNodeStuff(*rhs.nodeInfo_);
    nodeInfo_->maximumNodes_ = maximumNodes_;
    ClpNodeStuff *info = nodeInfo_;
    if (maximumDepth_ > 0) {
      info->nDepth_ = realMaximumDepth;
    } else {
      info->nDepth_ = -maximumDepth_;
      info->solverOptions_ |= 32;
    }
    if (!info->nodeInfo_) {
      ClpNode **nodeInfo = new ClpNode *[maximumNodes_];
      for (int i = 0; i < maximumNodes_; i++)
        nodeInfo[i] = NULL;
      info->nodeInfo_ = nodeInfo;
    }
  } else {
    nodeInfo_ = NULL;
  }
}

// Clone
CbcObject *
CbcGeneralDepth::clone() const
{
  return new CbcGeneralDepth(*this);
}

// Assignment operator
CbcGeneralDepth &
CbcGeneralDepth::operator=(const CbcGeneralDepth &rhs)
{
  if (this != &rhs) {
    CbcGeneral::operator=(rhs);
    delete nodeInfo_;
    maximumDepth_ = rhs.maximumDepth_;
    int realMaximumDepth = maximumDepth_ % MAX_DEPTH;
    maximumNodes_ = rhs.maximumNodes_;
    whichSolution_ = -1;
    numberNodes_ = 0;
    if (maximumDepth_) {
      assert(rhs.nodeInfo_);
      nodeInfo_ = new ClpNodeStuff(*rhs.nodeInfo_);
      nodeInfo_->maximumNodes_ = maximumNodes_;
    } else {
      nodeInfo_ = NULL;
    }
  }
  return *this;
}

// Destructor
CbcGeneralDepth::~CbcGeneralDepth()
{
  delete nodeInfo_;
}
// Infeasibility - large is 0.5
double
CbcGeneralDepth::infeasibility(const OsiBranchingInformation * /*info*/,
  int & /*preferredWay*/) const
{
  whichSolution_ = -1;
  if (maximumDepth_==1) {
    // DETERMINISTIC_TUNING think harder
    // vary
    cbc_node_count nfathoms = model_->getFathomCount();
    int depth =5 ;
    while (nfathoms && depth <20) {
      nfathoms /= 2;
      depth++;
    }
    int nnodes = (1 << depth) + 1 + depth;
    nnodes = std::min(nnodes, 1 + 20 + MAX_NODES);
    nnodes = std::max(nnodes, 10);
    nodeInfo_->maximumNodes_ = nnodes;
    nodeInfo_->nDepth_ = depth;
  }
  double returnValue;
  // should use genuine OsiBranchingInformation usefulInfo = model_->usefulInformation();
  // for now assume only called when correct
  //if (usefulInfo.depth_>=4&&!model_->parentModel()
  //     &&(usefulInfo.depth_%2)==0) {
  if (true) {
    OsiSolverInterface *solver = model_->solver();
    OsiClpSolverInterface *clpSolver
      = dynamic_cast< OsiClpSolverInterface * >(solver);
    if (clpSolver) {
      if ((model_->moreSpecialOptions() & 33554432) == 0) {
        ClpNodeStuff *info = nodeInfo_;
        info->integerTolerance_ = model_->getIntegerTolerance();
        info->integerIncrement_ = model_->getCutoffIncrement();
        info->numberBeforeTrust_ = model_->numberBeforeTrust();
        info->stateOfSearch_ = model_->stateOfSearch();
        // Compute "small" change in branch
        long int nBranches = model_->getNodeCount();
        if (nBranches) { 
          double average = model_->getDblParam(CbcModel::CbcSumChange) / static_cast< double >(nBranches);
          info->smallChange_ = std::max(average * 1.0e-5, model_->getDblParam(CbcModel::CbcSmallestChange));
          info->smallChange_ = std::max(info->smallChange_, 1.0e-8);
        } else {
          info->smallChange_ = 1.0e-8;
        }
        int numberIntegers = model_->numberIntegers();
        double *down = new double[numberIntegers];
        double *up = new double[numberIntegers];
        int *priority = new int[numberIntegers];
        int *numberDown = new int[numberIntegers];
        int *numberUp = new int[numberIntegers];
        int *numberDownInfeasible = new int[numberIntegers];
        int *numberUpInfeasible = new int[numberIntegers];
        model_->fillPseudoCosts(down, up, priority, numberDown, numberUp,
          numberDownInfeasible, numberUpInfeasible);
        info->fillPseudoCosts(down, up, priority, numberDown, numberUp,
          numberDownInfeasible,
          numberUpInfeasible, numberIntegers);
	// possible bug if 0 as can have mismatch on dense/nondense factorization
	// could check if number rows large enough
        info->presolveType_ = 1;
	ClpSimplex *simplex = clpSolver->getModelPtr();
	int numberRows = simplex->numberRows();
	if (simplex->factorization()->goDenseThreshold() < numberRows) {
	  const double *lower = simplex->columnLower();
	  const double *upper = simplex->columnUpper();
	  const int *integerVariable = model_->integerVariable();
	  int nFixed = 0;
	  for (int i=0;i<numberIntegers;i++) {
	    int iColumn = integerVariable[i];
	    if (upper[iColumn]==lower[iColumn])
	      nFixed++;
	  }
	  if (nFixed*2<numberIntegers) {
	    info->presolveType_=0;
	    for (int i=0;i<info->maximumNodes_;i++) {
	      if (info->nodeInfo_[i]) 
		info->nodeInfo_[i]->cleanUpForCrunch();
	    }
	  }
	}
        delete[] down;
        delete[] up;
	delete[] priority;
        delete[] numberDown;
        delete[] numberUp;
        delete[] numberDownInfeasible;
        delete[] numberUpInfeasible;
        bool takeHint;
        OsiHintStrength strength;
        solver->getHintParam(OsiDoReducePrint, takeHint, strength);
        int saveLevel = simplex->logLevel();
        if (strength != OsiHintIgnore && takeHint && saveLevel == 1)
          simplex->setLogLevel(0);
        clpSolver->setBasis();
	// DETERMINISTIC_TUNING think harder
	//info->presolveType_=1;
        whichSolution_ = simplex->fathomMany(info);
        //printf("FAT %d nodes, %d iterations\n",
        //info->numberNodesExplored_,info->numberIterations_);
        //printf("CbcBranch %d rows, %d columns\n",clpSolver->getNumRows(),
        //     clpSolver->getNumCols());
        model_->incrementExtra(info->numberNodesExplored_,
          info->numberIterations_);
        // update pseudo costs
        double smallest = 1.0e50;
        double largest = -1.0;
        OsiObject **objects = model_->objects();
#ifndef NDEBUG
        const int *integerVariable = model_->integerVariable();
#endif
        for (int i = 0; i < numberIntegers; i++) {
#ifndef NDEBUG
          CbcSimpleIntegerDynamicPseudoCost *obj = dynamic_cast< CbcSimpleIntegerDynamicPseudoCost * >(objects[i]);
          assert(obj && obj->columnNumber() == integerVariable[i]);
#else
          CbcSimpleIntegerDynamicPseudoCost *obj = static_cast< CbcSimpleIntegerDynamicPseudoCost * >(objects[i]);
#endif
          if (info->numberUp_[i] > 0) {
            if (info->downPseudo_[i] > largest)
              largest = info->downPseudo_[i];
            if (info->downPseudo_[i] < smallest)
              smallest = info->downPseudo_[i];
            if (info->upPseudo_[i] > largest)
              largest = info->upPseudo_[i];
            if (info->upPseudo_[i] < smallest)
              smallest = info->upPseudo_[i];
            obj->updateAfterMini(info->numberDown_[i],
              info->numberDownInfeasible_[i],
              info->downPseudo_[i],
              info->numberUp_[i],
              info->numberUpInfeasible_[i],
              info->upPseudo_[i]);
          }
        }
        //printf("range of costs %g to %g\n",smallest,largest);
        simplex->setLogLevel(saveLevel);
        numberNodes_ = info->nNodes_;
      } else {
        // Try diving
        // See if any diving heuristics set to do dive+save
        CbcHeuristicDive *dive = NULL;
        for (int i = 0; i < model_->numberHeuristics(); i++) {
          CbcHeuristicDive *possible = dynamic_cast< CbcHeuristicDive * >(model_->heuristic(i));
          if (possible && possible->maxSimplexIterations() == COIN_INT_MAX) {
            // if more than one then rotate later?
            //if (possible->canHeuristicRun()) {
            dive = possible;
            break;
          }
        }
        assert(dive); // otherwise moreSpecial should have been turned off
        CbcSubProblem **nodes = NULL;
        int branchState = dive->fathom(model_, numberNodes_, nodes);
        if (branchState) {
          printf("new solution\n");
          whichSolution_ = numberNodes_ - 1;
        } else {
          whichSolution_ = -1;
        }
#if 0
	    if (0) {
	      for (int iNode=0;iNode<numberNodes;iNode++) {
		//tree_->push(nodes[iNode]) ;
	      }
	      assert (node->nodeInfo());
	      if (node->nodeInfo()->numberBranchesLeft()) {
		tree_->push(node) ;
	      } else {
		node->setActive(false);
	      }
	    }
#endif
        //delete [] nodes;
        model_->setTemporaryPointer(reinterpret_cast< void * >(nodes));
        // end try diving
      }
      int numberDo = numberNodes_;
      if (numberDo > 0 || whichSolution_ >= 0) {
        returnValue = 0.5;
      } else {
        // no solution
        returnValue = COIN_DBL_MAX; // say infeasible
      }
    } else {
      returnValue = -1.0;
    }
  } else {
    returnValue = -1.0;
  }
  //nodeInfo_->maximumNodes_ = saveMaximumNodes;
  return returnValue;
}

// This looks at solution and sets bounds to contain solution
void CbcGeneralDepth::feasibleRegion()
{
  // Other stuff should have done this
}
// Redoes data when sequence numbers change
void CbcGeneralDepth::redoSequenceEtc(CbcModel * /*model*/,
  int /*numberColumns*/,
  const int * /*originalColumns*/)
{
}

//#define CHECK_PATH
#ifdef CHECK_PATH
extern const double *debuggerSolution_Z;
extern int numberColumns_Z;
extern int gotGoodNode_Z;
#endif
CbcBranchingObject *
CbcGeneralDepth::createCbcBranch(OsiSolverInterface *solver, const OsiBranchingInformation *info, int /*way*/)
{
  int numberDo = numberNodes_;
  if (whichSolution_ >= 0 && (model_->moreSpecialOptions() & 33554432) == 0) {
    numberDo--;
    if (numberDo <= 0)
      return NULL; // solution but complete search
  }
  assert(numberDo > 0);
  // create object
  CbcGeneralBranchingObject *branch = new CbcGeneralBranchingObject(model_);
  // skip solution
  branch->numberSubProblems_ = numberDo;
  // If parentBranch_ back in then will have to be 2*
  branch->numberSubLeft_ = numberDo;
  branch->setNumberBranches(numberDo);
  CbcSubProblem *sub = new CbcSubProblem[numberDo];
  int iProb = 0;
  branch->subProblems_ = sub;
  branch->numberRows_ = model_->solver()->getNumRows();
  int iNode;
  //OsiSolverInterface * solver = model_->solver();
  OsiClpSolverInterface *clpSolver
    = dynamic_cast< OsiClpSolverInterface * >(solver);
  assert(clpSolver);
  ClpSimplex *simplex = clpSolver->getModelPtr();
  int numberColumns = simplex->numberColumns();
  if ((model_->moreSpecialOptions() & 33554432) == 0) {
    double *lowerBefore = CoinCopyOfArray(simplex->getColLower(),
      numberColumns);
    double *upperBefore = CoinCopyOfArray(simplex->getColUpper(),
      numberColumns);
    ClpNodeStuff *info = nodeInfo_;
    double *weight = new double[numberNodes_];
    int *whichNode = new int[numberNodes_];
    // Sort
    for (iNode = 0; iNode < numberNodes_; iNode++) {
      if (iNode != whichSolution_) {
        double objectiveValue = info->nodeInfo_[iNode]->objectiveValue();
        double sumInfeasibilities = info->nodeInfo_[iNode]->sumInfeasibilities();
        int numberInfeasibilities = info->nodeInfo_[iNode]->numberInfeasibilities();
        double thisWeight = 0.0;
#if 1
        // just closest
        thisWeight = 1.0e9 * numberInfeasibilities;
        thisWeight += sumInfeasibilities;
        thisWeight += 1.0e-7 * objectiveValue;
        // Try estimate
        thisWeight = info->nodeInfo_[iNode]->estimatedSolution();
#else
        thisWeight = 1.0e-3 * numberInfeasibilities;
        thisWeight += 1.0e-5 * sumInfeasibilities;
        thisWeight += objectiveValue;
#endif
        whichNode[iProb] = iNode;
        weight[iProb++] = thisWeight;
      }
    }
    assert(iProb == numberDo);
    CoinSort_2(weight, weight + numberDo, whichNode);
    for (iProb = 0; iProb < numberDo; iProb++) {
      iNode = whichNode[iProb];
      ClpNode *node = info->nodeInfo_[iNode];
      // move bounds
      node->applyNode(simplex, 3);
      // create subproblem
      sub[iProb] = CbcSubProblem(clpSolver, lowerBefore, upperBefore,
        node->statusArray(), node->depth());
      sub[iProb].objectiveValue_ = node->objectiveValue();
      sub[iProb].sumInfeasibilities_ = node->sumInfeasibilities();
      sub[iProb].numberInfeasibilities_ = node->numberInfeasibilities();
#ifdef CHECK_PATH
      if (simplex->numberColumns() == numberColumns_Z) {
        bool onOptimal = true;
        const double *columnLower = simplex->columnLower();
        const double *columnUpper = simplex->columnUpper();
        for (int i = 0; i < numberColumns_Z; i++) {
          if (iNode == gotGoodNode_Z)
            printf("good %d %d %g %g\n", iNode, i, columnLower[i], columnUpper[i]);
          if (columnUpper[i] < debuggerSolution_Z[i] || columnLower[i] > debuggerSolution_Z[i] && simplex->isInteger(i)) {
            onOptimal = false;
            break;
          }
        }
        if (onOptimal) {
          printf("adding to node %x as %d - objs\n", this, iProb);
          for (int j = 0; j <= iProb; j++)
            printf("%d %g\n", j, sub[j].objectiveValue_);
        }
      }
#endif
    }
    delete[] weight;
    delete[] whichNode;
    const double *lower = solver->getColLower();
    const double *upper = solver->getColUpper();
    // restore bounds
    for (int j = 0; j < numberColumns; j++) {
      if (lowerBefore[j] != lower[j])
        solver->setColLower(j, lowerBefore[j]);
      if (upperBefore[j] != upper[j])
        solver->setColUpper(j, upperBefore[j]);
    }
    delete[] upperBefore;
    delete[] lowerBefore;
  } else {
    // from diving
    CbcSubProblem **nodes = reinterpret_cast< CbcSubProblem ** >(model_->temporaryPointer());
    assert(nodes);
    int adjustDepth = info->depth_;
    assert(numberDo);
    numberNodes_ = 0;
    for (iProb = 0; iProb < numberDo; iProb++) {
      if ((nodes[iProb]->problemStatus_ & 2) == 0) {
        // create subproblem (and swap way and/or make inactive)
        sub[numberNodes_].takeOver(*nodes[iProb], true);
        // but adjust depth
        sub[numberNodes_].depth_ += adjustDepth;
        numberNodes_++;
      }
      delete nodes[iProb];
    }
    branch->numberSubProblems_ = numberNodes_;
    branch->numberSubLeft_ = numberNodes_;
    branch->setNumberBranches(numberNodes_);
    if (!numberNodes_) {
      // infeasible
      delete branch;
      branch = NULL;
    }
    delete[] nodes;
  }
  return branch;
}

// Default Constructor
CbcGeneralBranchingObject::CbcGeneralBranchingObject()
  : CbcBranchingObject()
  , subProblems_(NULL)
  , node_(NULL)
  , numberSubProblems_(0)
  , numberSubLeft_(0)
  , whichNode_(-1)
  , numberRows_(0)
{
  //  printf("CbcGeneral %x default constructor\n",this);
}

// Useful constructor
CbcGeneralBranchingObject::CbcGeneralBranchingObject(CbcModel *model)
  : CbcBranchingObject(model, -1, -1, 0.5)
  , subProblems_(NULL)
  , node_(NULL)
  , numberSubProblems_(0)
  , numberSubLeft_(0)
  , whichNode_(-1)
  , numberRows_(0)
{
  //printf("CbcGeneral %x useful constructor\n",this);
}

// Copy constructor
CbcGeneralBranchingObject::CbcGeneralBranchingObject(const CbcGeneralBranchingObject &rhs)
  : CbcBranchingObject(rhs)
  , subProblems_(NULL)
  , node_(rhs.node_)
  , numberSubProblems_(rhs.numberSubProblems_)
  , numberSubLeft_(rhs.numberSubLeft_)
  , whichNode_(rhs.whichNode_)
  , numberRows_(rhs.numberRows_)
{
  abort();
  if (numberSubProblems_) {
    subProblems_ = new CbcSubProblem[numberSubProblems_];
    for (int i = 0; i < numberSubProblems_; i++)
      subProblems_[i] = rhs.subProblems_[i];
  }
}

// Assignment operator
CbcGeneralBranchingObject &
CbcGeneralBranchingObject::operator=(const CbcGeneralBranchingObject &rhs)
{
  if (this != &rhs) {
    abort();
    CbcBranchingObject::operator=(rhs);
    delete[] subProblems_;
    numberSubProblems_ = rhs.numberSubProblems_;
    numberSubLeft_ = rhs.numberSubLeft_;
    whichNode_ = rhs.whichNode_;
    numberRows_ = rhs.numberRows_;
    if (numberSubProblems_) {
      subProblems_ = new CbcSubProblem[numberSubProblems_];
      for (int i = 0; i < numberSubProblems_; i++)
        subProblems_[i] = rhs.subProblems_[i];
    } else {
      subProblems_ = NULL;
    }
    node_ = rhs.node_;
  }
  return *this;
}
CbcBranchingObject *
CbcGeneralBranchingObject::clone() const
{
  return (new CbcGeneralBranchingObject(*this));
}

// Destructor
CbcGeneralBranchingObject::~CbcGeneralBranchingObject()
{
  //printf("CbcGeneral %x destructor\n",this);
  delete[] subProblems_;
}
bool doingDoneBranch = false;
double
CbcGeneralBranchingObject::branch()
{
  double cutoff = model_->getCutoff();
  //printf("GenB %x whichNode %d numberLeft %d which %d\n",
  // this,whichNode_,numberBranchesLeft(),branchIndex());
  if (whichNode_ < 0) {
    assert(node_);
    bool applied = false;
    while (numberBranchesLeft()) {
      int which = branchIndex();
      decrementNumberBranchesLeft();
      CbcSubProblem *thisProb = subProblems_ + which;
      if (thisProb->objectiveValue_ < cutoff) {
        //printf("branch %x (sub %x) which now %d\n",this,
        //     subProblems_,which);
        OsiSolverInterface *solver = model_->solver();
        thisProb->apply(solver);
        OsiClpSolverInterface *clpSolver
          = dynamic_cast< OsiClpSolverInterface * >(solver);
        assert(clpSolver);
        // Move status to basis
        clpSolver->setWarmStart(NULL);
        //ClpSimplex * simplex = clpSolver->getModelPtr();
        node_->setObjectiveValue(thisProb->objectiveValue_);
        node_->setSumInfeasibilities(thisProb->sumInfeasibilities_);
        node_->setNumberUnsatisfied(thisProb->numberInfeasibilities_);
        applied = true;
        doingDoneBranch = true;
        break;
      } else if (numberBranchesLeft()) {
        node_->nodeInfo()->branchedOn();
      }
    }
    if (!applied) {
      // no good one
      node_->setObjectiveValue(cutoff + 1.0e20);
      node_->setSumInfeasibilities(1.0);
      node_->setNumberUnsatisfied(1);
      assert(whichNode_ < 0);
    }
  } else {
    decrementNumberBranchesLeft();
    CbcSubProblem *thisProb = subProblems_ + whichNode_;
    assert(thisProb->objectiveValue_ < cutoff);
    OsiSolverInterface *solver = model_->solver();
    thisProb->apply(solver);
    //OsiClpSolverInterface * clpSolver
    //= dynamic_cast<OsiClpSolverInterface *> (solver);
    //assert (clpSolver);
    // Move status to basis
    //clpSolver->setWarmStart(NULL);
  }
  return 0.0;
}
/* Double checks in case node can change its mind!
   Can change objective etc */
void CbcGeneralBranchingObject::checkIsCutoff(double cutoff)
{
  assert(node_);
  int first = branchIndex();
  int last = first + numberBranchesLeft();
  for (int which = first; which < last; which++) {
    CbcSubProblem *thisProb = subProblems_ + which;
    if (thisProb->objectiveValue_ < cutoff) {
      node_->setObjectiveValue(thisProb->objectiveValue_);
      node_->setSumInfeasibilities(thisProb->sumInfeasibilities_);
      node_->setNumberUnsatisfied(thisProb->numberInfeasibilities_);
      break;
    }
  }
}
// Print what would happen
void CbcGeneralBranchingObject::print()
{
  //printf("CbcGeneralObject has %d subproblems\n",numberSubProblems_);
}
// Fill in current objective etc
void CbcGeneralBranchingObject::state(double &objectiveValue,
  double &sumInfeasibilities,
  int &numberUnsatisfied, int which) const
{
  assert(which >= 0 && which < numberSubProblems_);
  const CbcSubProblem *thisProb = subProblems_ + which;
  objectiveValue = thisProb->objectiveValue_;
  sumInfeasibilities = thisProb->sumInfeasibilities_;
  numberUnsatisfied = thisProb->numberInfeasibilities_;
}
/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type.
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int CbcGeneralBranchingObject::compareOriginalObject(const CbcBranchingObject * /*brObj*/) const
{
  throw("must implement");
}

/** Compare the \c this with \c brObj. \c this and \c brObj must be os the
    same type and must have the same original object, but they may have
    different feasible regions.
    Return the appropriate CbcRangeCompare value (first argument being the
    sub/superset if that's the case). In case of overlap (and if \c
    replaceIfOverlap is true) replace the current branching object with one
    whose feasible region is the overlap.
*/
CbcRangeCompare
CbcGeneralBranchingObject::compareBranchingObject(const CbcBranchingObject * /*brObj*/, const bool /*replaceIfOverlap*/)
{
  throw("must implement");
}

// Default Constructor
CbcOneGeneralBranchingObject::CbcOneGeneralBranchingObject()
  : CbcBranchingObject()
  , object_(NULL)
  , whichOne_(-1)
{
  //printf("CbcOneGeneral %x default constructor\n",this);
}

// Useful constructor
CbcOneGeneralBranchingObject::CbcOneGeneralBranchingObject(CbcModel *model,
  CbcGeneralBranchingObject *object,
  int whichOne)
  : CbcBranchingObject(model, -1, -1, 0.5)
  , object_(object)
  , whichOne_(whichOne)
{
  //printf("CbcOneGeneral %x useful constructor object %x %d left\n",this,
  //	 object_,object_->numberSubLeft_);
  numberBranches_ = 1;
}

// Copy constructor
CbcOneGeneralBranchingObject::CbcOneGeneralBranchingObject(const CbcOneGeneralBranchingObject &rhs)
  : CbcBranchingObject(rhs)
  , object_(rhs.object_)
  , whichOne_(rhs.whichOne_)
{
}

// Assignment operator
CbcOneGeneralBranchingObject &
CbcOneGeneralBranchingObject::operator=(const CbcOneGeneralBranchingObject &rhs)
{
  if (this != &rhs) {
    CbcBranchingObject::operator=(rhs);
    object_ = rhs.object_;
    whichOne_ = rhs.whichOne_;
  }
  return *this;
}
CbcBranchingObject *
CbcOneGeneralBranchingObject::clone() const
{
  return (new CbcOneGeneralBranchingObject(*this));
}

// Destructor
CbcOneGeneralBranchingObject::~CbcOneGeneralBranchingObject()
{
  //printf("CbcOneGeneral %x destructor object %x %d left\n",this,
  // object_,object_->numberSubLeft_);
  assert(object_->numberSubLeft_ > 0 && object_->numberSubLeft_ < 1000000);
  if (!object_->decrementNumberLeft()) {
    // printf("CbcGeneral %x yy destructor\n",object_);
    delete object_;
  }
}
double
CbcOneGeneralBranchingObject::branch()
{
  assert(numberBranchesLeft());
  decrementNumberBranchesLeft();
  assert(!numberBranchesLeft());
  object_->setWhichNode(whichOne_);
  object_->setModel(model_);
  object_->branch();
  return 0.0;
}
/* Double checks in case node can change its mind!
   Can change objective etc */
void CbcOneGeneralBranchingObject::checkIsCutoff(double /*cutoff*/)
{
  assert(numberBranchesLeft());
}
// Print what would happen
void CbcOneGeneralBranchingObject::print()
{
  //printf("CbcOneGeneralObject has 1 subproblem\n");
}
/** Compare the original object of \c this with the original object of \c
    brObj. Assumes that there is an ordering of the original objects.
    This method should be invoked only if \c this and brObj are of the same
    type.
    Return negative/0/positive depending on whether \c this is
    smaller/same/larger than the argument.
*/
int CbcOneGeneralBranchingObject::compareOriginalObject(const CbcBranchingObject * /*brObj*/) const
{
  throw("must implement");
}

/** Compare the \c this with \c brObj. \c this and \c brObj must be os the
    same type and must have the same original object, but they may have
    different feasible regions.
    Return the appropriate CbcRangeCompare value (first argument being the
    sub/superset if that's the case). In case of overlap (and if \c
    replaceIfOverlap is true) replace the current branching object with one
    whose feasible region is the overlap.
*/
CbcRangeCompare
CbcOneGeneralBranchingObject::compareBranchingObject(const CbcBranchingObject * /*brObj*/, const bool /*replaceIfOverlap*/)
{
  throw("must implement");
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
