/* $Id$ */
// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if defined(_MSC_VER)
// Turn off compiler warning about long names
#pragma warning(disable : 4786)
#endif

#include "CbcConfig.h"

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>

//#define PRINT_DEBUG
#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#endif
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcHeuristic.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristicRINS.hpp"
#include "CbcEventHandler.hpp"
#include "CbcStrategy.hpp"
#include "CglPreProcess.hpp"
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "OsiAuxInfo.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiPresolve.hpp"
#include "CbcBranchActual.hpp"
#include "CbcCutGenerator.hpp"
#include "CoinMpsIO.hpp"
//==============================================================================

CbcHeuristicNode::CbcHeuristicNode(const CbcHeuristicNode &rhs)
{
  numObjects_ = rhs.numObjects_;
  brObj_ = new CbcBranchingObject *[numObjects_];
  for (int i = 0; i < numObjects_; ++i) {
    brObj_[i] = rhs.brObj_[i]->clone();
  }
}

void CbcHeuristicNodeList::gutsOfDelete()
{
  for (int i = (static_cast< int >(nodes_.size())) - 1; i >= 0; --i) {
    delete nodes_[i];
  }
}

void CbcHeuristicNodeList::gutsOfCopy(const CbcHeuristicNodeList &rhs)
{
  append(rhs);
}

CbcHeuristicNodeList::CbcHeuristicNodeList(const CbcHeuristicNodeList &rhs)
{
  gutsOfCopy(rhs);
}

CbcHeuristicNodeList &CbcHeuristicNodeList::operator=(const CbcHeuristicNodeList &rhs)
{
  if (this != &rhs) {
    gutsOfDelete();
    gutsOfCopy(rhs);
  }
  return *this;
}

CbcHeuristicNodeList::~CbcHeuristicNodeList()
{
  gutsOfDelete();
}

void CbcHeuristicNodeList::append(CbcHeuristicNode *&node)
{
  nodes_.push_back(node);
  node = NULL;
}

void CbcHeuristicNodeList::append(const CbcHeuristicNodeList &nodes)
{
  nodes_.reserve(nodes_.size() + nodes.size());
  for (int i = 0; i < nodes.size(); ++i) {
    CbcHeuristicNode *node = new CbcHeuristicNode(*nodes.node(i));
    append(node);
  }
}

//==============================================================================
#define DEFAULT_WHERE ((255 - 2 - 16) * (1 + 256))
// Default Constructor
CbcHeuristic::CbcHeuristic()
  : model_(NULL)
  , when_(2)
  , numberNodes_(200)
  , feasibilityPumpOptions_(-1)
  , fractionSmall_(1.0)
  , heuristicName_("Unknown")
  , howOften_(1)
  , decayFactor_(0.0)
  , switches_(0)
  , whereFrom_(DEFAULT_WHERE)
  , shallowDepth_(1)
  , howOftenShallow_(1)
  , numInvocationsInShallow_(0)
  , numInvocationsInDeep_(0)
  , lastRunDeep_(0)
  , numRuns_(0)
  , minDistanceToRun_(1)
  , runNodes_()
  , numCouldRun_(0)
  , numberSolutionsFound_(0)
  , numberNodesDone_(0)
  , inputSolution_(NULL)
{
  // As CbcHeuristic virtual need to modify .cpp if above change
}

// Constructor from model
CbcHeuristic::CbcHeuristic(CbcModel &model)
  : model_(&model)
  , when_(2)
  , numberNodes_(200)
  , feasibilityPumpOptions_(-1)
  , fractionSmall_(1.0)
  , heuristicName_("Unknown")
  , howOften_(1)
  , decayFactor_(0.0)
  , switches_(0)
  , whereFrom_(DEFAULT_WHERE)
  , shallowDepth_(1)
  , howOftenShallow_(1)
  , numInvocationsInShallow_(0)
  , numInvocationsInDeep_(0)
  , lastRunDeep_(0)
  , numRuns_(0)
  , minDistanceToRun_(1)
  , runNodes_()
  , numCouldRun_(0)
  , numberSolutionsFound_(0)
  , numberNodesDone_(0)
  , inputSolution_(NULL)
{
}

void CbcHeuristic::gutsOfCopy(const CbcHeuristic &rhs)
{
  model_ = rhs.model_;
  when_ = rhs.when_;
  numberNodes_ = rhs.numberNodes_;
  feasibilityPumpOptions_ = rhs.feasibilityPumpOptions_;
  fractionSmall_ = rhs.fractionSmall_;
  randomNumberGenerator_ = rhs.randomNumberGenerator_;
  heuristicName_ = rhs.heuristicName_;
  howOften_ = rhs.howOften_;
  decayFactor_ = rhs.decayFactor_;
  switches_ = rhs.switches_;
  whereFrom_ = rhs.whereFrom_;
  shallowDepth_ = rhs.shallowDepth_;
  howOftenShallow_ = rhs.howOftenShallow_;
  numInvocationsInShallow_ = rhs.numInvocationsInShallow_;
  numInvocationsInDeep_ = rhs.numInvocationsInDeep_;
  lastRunDeep_ = rhs.lastRunDeep_;
  numRuns_ = rhs.numRuns_;
  numCouldRun_ = rhs.numCouldRun_;
  minDistanceToRun_ = rhs.minDistanceToRun_;
  runNodes_ = rhs.runNodes_;
  numberSolutionsFound_ = rhs.numberSolutionsFound_;
  numberNodesDone_ = rhs.numberNodesDone_;
  if (rhs.inputSolution_) {
    int numberColumns = model_->getNumCols();
    setInputSolution(rhs.inputSolution_, rhs.inputSolution_[numberColumns]);
  }
}
// Copy constructor
CbcHeuristic::CbcHeuristic(const CbcHeuristic &rhs)
{
  inputSolution_ = NULL;
  gutsOfCopy(rhs);
}

// Assignment operator
CbcHeuristic &
CbcHeuristic::operator=(const CbcHeuristic &rhs)
{
  if (this != &rhs) {
    gutsOfDelete();
    gutsOfCopy(rhs);
  }
  return *this;
}

static
void CbcHeurDebugNodes(CbcModel *model_)
{
  CbcNode *node = model_->currentNode();
  CbcNodeInfo *nodeInfo = node->nodeInfo();
  std::cout << "===============================================================\n";
  while (nodeInfo) {
    const CbcNode *node = nodeInfo->owner();
    printf("nodeinfo: node %i\n", nodeInfo->nodeNumber());
    {
      const CbcIntegerBranchingObject *brPrint = dynamic_cast< const CbcIntegerBranchingObject * >(nodeInfo->parentBranch());
      if (!brPrint) {
        printf("    parentBranch: NULL\n");
      } else {
        const double *downBounds = brPrint->downBounds();
        const double *upBounds = brPrint->upBounds();
        int variable = brPrint->variable();
        int way = brPrint->way();
        printf("   parentBranch: var %i downBd [%i,%i] upBd [%i,%i] way %i\n",
          variable, static_cast< int >(downBounds[0]), static_cast< int >(downBounds[1]),
          static_cast< int >(upBounds[0]), static_cast< int >(upBounds[1]), way);
      }
    }
    if (!node) {
      printf("    owner: NULL\n");
    } else {
      printf("    owner: node %i depth %i onTree %i active %i",
        node->nodeNumber(), node->depth(), node->onTree(), node->active());
      const OsiBranchingObject *osibr = nodeInfo->owner()->branchingObject();
      const CbcBranchingObject *cbcbr = dynamic_cast< const CbcBranchingObject * >(osibr);
      const CbcIntegerBranchingObject *brPrint = dynamic_cast< const CbcIntegerBranchingObject * >(cbcbr);
      if (!brPrint) {
        printf("        ownerBranch: NULL\n");
      } else {
        const double *downBounds = brPrint->downBounds();
        const double *upBounds = brPrint->upBounds();
        int variable = brPrint->variable();
        int way = brPrint->way();
        printf("        ownerbranch: var %i downBd [%i,%i] upBd [%i,%i] way %i\n",
          variable, static_cast< int >(downBounds[0]), static_cast< int >(downBounds[1]),
          static_cast< int >(upBounds[0]), static_cast< int >(upBounds[1]), way);
      }
    }
    nodeInfo = nodeInfo->parent();
  }
}

void CbcHeuristic::debugNodes()
{
  CbcHeurDebugNodes(model_);
}

void CbcHeuristic::printDistanceToNodes()
{
  const CbcNode *currentNode = model_->currentNode();
  if (currentNode != NULL) {
    CbcHeuristicNode *nodeDesc = new CbcHeuristicNode(*model_);
    for (int i = runNodes_.size() - 1; i >= 0; --i) {
      nodeDesc->distance(runNodes_.node(i));
    }
    runNodes_.append(nodeDesc);
  }
}

bool CbcHeuristic::shouldHeurRun(int whereFrom)
{
  assert(whereFrom >= 0 && whereFrom < 16);
  // take off 8 (code - likes new solution)
  whereFrom &= 7;
  if ((whereFrom_ & (1 << whereFrom)) == 0)
    return false;
    // No longer used for original purpose - so use for ever run at all JJF
#ifndef JJF_ONE
  // Don't run if hot start or no rows!
  if (model_ && (model_->hotstartSolution() || !model_->getNumRows()))
    return false;
  else
    return true;
#else
#ifdef JJF_ZERO
  const CbcNode *currentNode = model_->currentNode();
  if (currentNode == NULL) {
    return false;
  }

  debugNodes();
  //   return false;

  const int depth = currentNode->depth();
#else
  int depth = model_->currentDepth();
#endif

  const int nodeCount = model_->getNodeCount(); // FIXME: check that this is
  // correct in parallel

  if (nodeCount == 0 || depth <= shallowDepth_) {
    // what to do when we are in the shallow part of the tree
    if (model_->getCurrentPassNumber() == 1) {
      // first time in the node...
      numInvocationsInShallow_ = 0;
    }
    ++numInvocationsInShallow_;
    // Very large howOftenShallow_ will give the original test:
    // (model_->getCurrentPassNumber() != 1)
    //    if ((numInvocationsInShallow_ % howOftenShallow_) != 1) {
    if ((numInvocationsInShallow_ % howOftenShallow_) != 0) {
      return false;
    }
    // LL: should we save these nodes in the list of nodes where the heur was
    // LL: run?
#ifndef JJF_ONE
    if (currentNode != NULL) {
      // Get where we are and create the appropriate CbcHeuristicNode object
      CbcHeuristicNode *nodeDesc = new CbcHeuristicNode(*model_);
      runNodes_.append(nodeDesc);
    }
#endif
  } else {
    // deeper in the tree
    if (model_->getCurrentPassNumber() == 1) {
      // first time in the node...
      ++numInvocationsInDeep_;
    }
    if (numInvocationsInDeep_ - lastRunDeep_ < howOften_) {
      return false;
    }
    if (model_->getCurrentPassNumber() > 1) {
      // Run the heuristic only when first entering the node.
      // LL: I don't think this is right. It should run just before strong
      // LL: branching, I believe.
      return false;
    }
    // Get where we are and create the appropriate CbcHeuristicNode object
    CbcHeuristicNode *nodeDesc = new CbcHeuristicNode(*model_);
    //#ifdef PRINT_DEBUG
#ifndef JJF_ONE
    const double minDistanceToRun = 1.5 * log((double)depth) / log((double)2);
#else
    const double minDistanceToRun = minDistanceToRun_;
#endif
#ifdef PRINT_DEBUG
    double minDistance = nodeDesc->minDistance(runNodes_);
    std::cout << "minDistance = " << minDistance
              << ", minDistanceToRun = " << minDistanceToRun << std::endl;
#endif
    if (nodeDesc->minDistanceIsSmall(runNodes_, minDistanceToRun)) {
      delete nodeDesc;
      return false;
    }
    runNodes_.append(nodeDesc);
    lastRunDeep_ = numInvocationsInDeep_;
    //    ++lastRunDeep_;
  }
  ++numRuns_;
  return true;
#endif
}

bool CbcHeuristic::shouldHeurRun_randomChoice()
{
  if (!when_)
    return false;
  int depth = model_->currentDepth();
  // when_ -999 is special marker to force to run
  if (depth != 0 && when_ != -999) {
    const double numerator = depth * depth;
    const double denominator = exp(depth * log(2.0));
    double probability = numerator / denominator;
    double randomNumber = randomNumberGenerator_.randomDouble();
    int when = when_ % 100;
    if (when > 2 && when < 8) {
      /* JJF adjustments
            3 only at root and if no solution
            4 only at root and if this heuristic has not got solution
            5 decay (but only if no solution)
            6 if depth <3 or decay
            7 run up to 2 times if solution found 4 otherwise
            */
      switch (when) {
      case 3:
      default:
        if (model_->bestSolution())
          probability = -1.0;
        break;
      case 4:
        if (numberSolutionsFound_)
          probability = -1.0;
        break;
      case 5:
        assert(decayFactor_);
        if (model_->bestSolution()) {
          probability = -1.0;
        } else if (numCouldRun_ > 1000) {
          decayFactor_ *= 0.99;
          probability *= decayFactor_;
        }
        break;
      case 6:
        if (depth >= 3) {
          if ((numCouldRun_ % howOften_) == 0 && numberSolutionsFound_ * howOften_ < numCouldRun_) {
            //#define COIN_DEVELOP
#ifdef COIN_DEVELOP
            int old = howOften_;
#endif
            howOften_ = CoinMin(CoinMax(static_cast< int >(howOften_ * 1.1), howOften_ + 1), 1000000);
#ifdef COIN_DEVELOP
            printf("Howoften changed from %d to %d for %s\n",
              old, howOften_, heuristicName_.c_str());
#endif
          }
          probability = 1.0 / howOften_;
          if (model_->bestSolution())
            probability *= 0.5;
        } else {
          probability = 1.1;
        }
        break;
      case 7:
        if ((model_->bestSolution() && numRuns_ >= 2) || numRuns_ >= 4)
          probability = -1.0;
        break;
      }
    }
    if (randomNumber > probability)
      return false;

    if (model_->getCurrentPassNumber() > 1)
      return false;
#ifdef COIN_DEVELOP
    printf("Running %s, random %g probability %g\n",
      heuristicName_.c_str(), randomNumber, probability);
#endif
  } else {
#ifdef COIN_DEVELOP
    printf("Running %s, depth %d when %d\n",
      heuristicName_.c_str(), depth, when_);
#endif
  }
  ++numRuns_;
  return true;
}

// Resets stuff if model changes
void CbcHeuristic::resetModel(CbcModel *model)
{
  model_ = model;
}
// Set seed
void CbcHeuristic::setSeed(int value)
{
  if (value == 0) {
    double time = fabs(CoinGetTimeOfDay());
    while (time >= COIN_INT_MAX)
      time *= 0.5;
    value = static_cast< int >(time);
    char printArray[100];
    sprintf(printArray, "using time of day seed was changed from %d to %d",
      randomNumberGenerator_.getSeed(), value);
    if (model_)
      model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
        << printArray
        << CoinMessageEol;
  }
  randomNumberGenerator_.setSeed(value);
}
// Get seed
int CbcHeuristic::getSeed() const
{
  return randomNumberGenerator_.getSeed();
}

// Create C++ lines to get to current state
void CbcHeuristic::generateCpp(FILE *fp, const char *heuristic)
{
  // hard coded as CbcHeuristic virtual
  if (when_ != 2)
    fprintf(fp, "3  %s.setWhen(%d);\n", heuristic, when_);
  else
    fprintf(fp, "4  %s.setWhen(%d);\n", heuristic, when_);
  if (numberNodes_ != 200)
    fprintf(fp, "3  %s.setNumberNodes(%d);\n", heuristic, numberNodes_);
  else
    fprintf(fp, "4  %s.setNumberNodes(%d);\n", heuristic, numberNodes_);
  if (feasibilityPumpOptions_ != -1)
    fprintf(fp, "3  %s.setFeasibilityPumpOptions(%d);\n", heuristic, feasibilityPumpOptions_);
  else
    fprintf(fp, "4  %s.setFeasibilityPumpOptions(%d);\n", heuristic, feasibilityPumpOptions_);
  if (fractionSmall_ != 1.0)
    fprintf(fp, "3  %s.setFractionSmall(%g);\n", heuristic, fractionSmall_);
  else
    fprintf(fp, "4  %s.setFractionSmall(%g);\n", heuristic, fractionSmall_);
  if (heuristicName_ != "Unknown")
    fprintf(fp, "3  %s.setHeuristicName(\"%s\");\n",
      heuristic, heuristicName_.c_str());
  else
    fprintf(fp, "4  %s.setHeuristicName(\"%s\");\n",
      heuristic, heuristicName_.c_str());
  if (decayFactor_ != 0.0)
    fprintf(fp, "3  %s.setDecayFactor(%g);\n", heuristic, decayFactor_);
  else
    fprintf(fp, "4  %s.setDecayFactor(%g);\n", heuristic, decayFactor_);
  if (switches_ != 0)
    fprintf(fp, "3  %s.setSwitches(%d);\n", heuristic, switches_);
  else
    fprintf(fp, "4  %s.setSwitches(%d);\n", heuristic, switches_);
  if (whereFrom_ != DEFAULT_WHERE)
    fprintf(fp, "3  %s.setWhereFrom(%d);\n", heuristic, whereFrom_);
  else
    fprintf(fp, "4  %s.setWhereFrom(%d);\n", heuristic, whereFrom_);
  if (shallowDepth_ != 1)
    fprintf(fp, "3  %s.setShallowDepth(%d);\n", heuristic, shallowDepth_);
  else
    fprintf(fp, "4  %s.setShallowDepth(%d);\n", heuristic, shallowDepth_);
  if (howOftenShallow_ != 1)
    fprintf(fp, "3  %s.setHowOftenShallow(%d);\n", heuristic, howOftenShallow_);
  else
    fprintf(fp, "4  %s.setHowOftenShallow(%d);\n", heuristic, howOftenShallow_);
  if (minDistanceToRun_ != 1)
    fprintf(fp, "3  %s.setMinDistanceToRun(%d);\n", heuristic, minDistanceToRun_);
  else
    fprintf(fp, "4  %s.setMinDistanceToRun(%d);\n", heuristic, minDistanceToRun_);
}
// Destructor
CbcHeuristic::~CbcHeuristic()
{
  delete[] inputSolution_;
}

// update model
void CbcHeuristic::setModel(CbcModel *model)
{
  model_ = model;
}
/* Clone but ..
   type 0 clone solver, 1 clone continuous solver
   Add 2 to say without integer variables which are at low priority
   Add 4 to say quite likely infeasible so give up easily.*/
OsiSolverInterface *
CbcHeuristic::cloneBut(int type)
{
  OsiSolverInterface *solver;
  if ((type & 1) == 0 || !model_->continuousSolver())
    solver = model_->solver()->clone();
  else
    solver = model_->continuousSolver()->clone();
#ifdef COIN_HAS_CLP
  OsiClpSolverInterface *clpSolver
    = dynamic_cast< OsiClpSolverInterface * >(solver);
#endif
  if ((type & 2) != 0) {
    int n = model_->numberObjects();
    int priority = model_->continuousPriority();
    if (priority < COIN_INT_MAX) {
      for (int i = 0; i < n; i++) {
        const OsiObject *obj = model_->object(i);
        const CbcSimpleInteger *thisOne = dynamic_cast< const CbcSimpleInteger * >(obj);
        if (thisOne) {
          int iColumn = thisOne->columnNumber();
          if (thisOne->priority() >= priority)
            solver->setContinuous(iColumn);
        }
      }
    }
#ifdef COIN_HAS_CLP
    if (clpSolver) {
      for (int i = 0; i < n; i++) {
        const OsiObject *obj = model_->object(i);
        const CbcSimpleInteger *thisOne = dynamic_cast< const CbcSimpleInteger * >(obj);
        if (thisOne) {
          int iColumn = thisOne->columnNumber();
          if (clpSolver->isOptionalInteger(iColumn))
            clpSolver->setContinuous(iColumn);
        }
      }
    }
#endif
  }
#ifdef COIN_HAS_CLP
  if ((type & 4) != 0 && clpSolver) {
    int options = clpSolver->getModelPtr()->moreSpecialOptions();
    clpSolver->getModelPtr()->setMoreSpecialOptions(options | 64);
  }
  if (clpSolver) {
    // take out zero cost integers which will be integer anyway
    const double *rowLower = clpSolver->getRowLower();
    const double *rowUpper = clpSolver->getRowUpper();
    const double *objective = clpSolver->getObjCoefficients();
    int numberRows = clpSolver->getNumRows();
    const CoinPackedMatrix *matrixByRow = clpSolver->getMatrixByRow();
    const double *elementByRow = matrixByRow->getElements();
    const int *column = matrixByRow->getIndices();
    const CoinBigIndex *rowStart = matrixByRow->getVectorStarts();
    const int *rowLength = matrixByRow->getVectorLengths();
    const CoinPackedMatrix *matrixByColumn = clpSolver->getMatrixByCol();
    //const double * element = matrixByColumn->getElements();
    //const int *row = matrixByColumn->getIndices();
    //const CoinBigIndex *columnStart = matrixByColumn->getVectorStarts();
    const int *columnLength = matrixByColumn->getVectorLengths();
    for (int i = 0; i < numberRows; i++) {
      if (rowLower[i] != floor(rowLower[i]) ||
	  rowUpper[i] != floor(rowUpper[i]))
	continue;
      int jColumn = -1;
      for (CoinBigIndex k = rowStart[i]; k < rowStart[i] + rowLength[i]; k++) {
	int iColumn = column[k];
	double value = elementByRow[k];
	if (!clpSolver->isInteger(iColumn) || floor(value) != value) {
	  jColumn = -2;
	  break;
	} else if (!objective[iColumn] && columnLength[iColumn] == 1) {
	  jColumn = iColumn;
	}
      }
      if (jColumn>=0)
	clpSolver->setContinuous(jColumn);
    }
  }
#endif
  return solver;
}
// Whether to exit at once on gap
bool CbcHeuristic::exitNow(double bestObjective) const
{
  if ((switches_ & 2048) != 0) {
    // exit may be forced - but unset for next time
    switches_ &= ~2048;
    if ((switches_ & 1024) != 0)
      return true;
  } else if ((switches_ & 1) == 0) {
    return false;
  }
  // See if can stop on gap
  OsiSolverInterface *solver = model_->solver();
  double bestPossibleObjective = solver->getObjValue() * solver->getObjSense();
  double absGap = CoinMax(model_->getAllowableGap(),
    model_->getHeuristicGap());
  double fracGap = CoinMax(model_->getAllowableFractionGap(),
    model_->getHeuristicFractionGap());
  double testGap = CoinMax(absGap, fracGap * CoinMax(fabs(bestObjective), fabs(bestPossibleObjective)));

  if (bestObjective - bestPossibleObjective < testGap
    && model_->getCutoffIncrement() >= 0.0) {
    return true;
  } else {
    return false;
  }
}
#ifdef HISTORY_STATISTICS
extern bool getHistoryStatistics_;
#endif
static double sizeRatio(int numberRowsNow, int numberColumnsNow,
  int numberRowsStart, int numberColumnsStart)
{
  double valueNow;
  if (numberRowsNow * 10 > numberColumnsNow || numberColumnsNow < 200) {
    valueNow = 2 * numberRowsNow + numberColumnsNow;
  } else {
    // long and thin - rows are more important
    if (numberRowsNow * 40 > numberColumnsNow)
      valueNow = 10 * numberRowsNow + numberColumnsNow;
    else
      valueNow = 200 * numberRowsNow + numberColumnsNow;
  }
  double valueStart;
  if (numberRowsStart * 10 > numberColumnsStart || numberColumnsStart < 200) {
    valueStart = 2 * numberRowsStart + numberColumnsStart;
  } else {
    // long and thin - rows are more important
    if (numberRowsStart * 40 > numberColumnsStart)
      valueStart = 10 * numberRowsStart + numberColumnsStart;
    else
      valueStart = 200 * numberRowsStart + numberColumnsStart;
  }
  //printf("sizeProblem Now %g, %d rows, %d columns\nsizeProblem Start %g, %d rows, %d columns\n",
  // valueNow,numberRowsNow,numberColumnsNow,
  // valueStart,numberRowsStart,numberColumnsStart);
  if (10 * numberRowsNow < 8 * numberRowsStart || 10 * numberColumnsNow < 7 * numberColumnsStart)
    return valueNow / valueStart;
  else if (10 * numberRowsNow < 9 * numberRowsStart)
    return 1.1 * (valueNow / valueStart);
  else if (numberRowsNow < numberRowsStart)
    return 1.5 * (valueNow / valueStart);
  else
    return 2.0 * (valueNow / valueStart);
}

//static int saveModel=0;
// Do mini branch and bound (return 1 if solution)
int CbcHeuristic::smallBranchAndBound(OsiSolverInterface *solver, int numberNodes,
  double *newSolution, double &newSolutionValue,
  double cutoff, std::string name) const
{
  CbcEventHandler *eventHandler = model_->getEventHandler();
  // Use this fraction
  double fractionSmall = fractionSmall_;
  int maximumSolutions = model_->getMaximumSolutions();
  int iterationMultiplier = 100;
  if (eventHandler) {
    typedef struct {
      double fractionSmall;
      double spareDouble[3];
      OsiSolverInterface *solver;
      void *sparePointer[2];
      int numberNodes;
      int maximumSolutions;
      int iterationMultiplier;
      int howOften;
      int spareInt[3];
    } SmallMod;
    SmallMod temp;
    temp.solver = solver;
    temp.fractionSmall = fractionSmall;
    temp.numberNodes = numberNodes;
    temp.iterationMultiplier = iterationMultiplier;
    temp.howOften = howOften_;
    temp.maximumSolutions = maximumSolutions;
    CbcEventHandler::CbcAction status = eventHandler->event(CbcEventHandler::smallBranchAndBound,
      &temp);
    if (status == CbcEventHandler::killSolution)
      return -1;
    if (status == CbcEventHandler::takeAction) {
      fractionSmall = temp.fractionSmall;
      numberNodes = temp.numberNodes;
      iterationMultiplier = temp.iterationMultiplier;
      howOften_ = temp.howOften;
      maximumSolutions = temp.maximumSolutions;
    }
  }
#if 0 
  if (saveModel || model_->getMaximumSolutions()==100) {
    printf("writing model\n");
    solver->writeMpsNative("before.mps", NULL, NULL, 2, 1);
  }
#endif
  // size before
  int shiftRows = 0;
  if (numberNodes < 0)
    shiftRows = solver->getNumRows() - numberNodes_;
  int numberRowsStart = solver->getNumRows() - shiftRows;
  int numberColumnsStart = solver->getNumCols();
#ifdef CLP_INVESTIGATE
  printf("%s has %d rows, %d columns\n",
    name.c_str(), solver->getNumRows(), solver->getNumCols());
#endif
  double before = 2 * numberRowsStart + numberColumnsStart;
  if (before > 40000.0) {
    // fairly large - be more conservative
    double multiplier = 1.0 - 0.3 * CoinMin(100000.0, before - 40000.0) / 100000.0;
    if (multiplier < 1.0) {
      fractionSmall *= multiplier;
#ifdef CLP_INVESTIGATE
      printf("changing fractionSmall from %g to %g for %s\n",
        fractionSmall_, fractionSmall, name.c_str());
#endif
    }
  }
#ifdef COIN_HAS_CLP
  OsiClpSolverInterface *clpSolver = dynamic_cast< OsiClpSolverInterface * >(solver);
  if (clpSolver && (clpSolver->specialOptions() & 65536) == 0) {
    // go faster stripes
    if (clpSolver->getNumRows() < 300 && clpSolver->getNumCols() < 500) {
      clpSolver->setupForRepeatedUse(2, 0);
    } else {
      clpSolver->setupForRepeatedUse(0, 0);
    }
    // Turn this off if you get problems
    // Used to be automatically set
    clpSolver->setSpecialOptions(clpSolver->specialOptions() | (128 + 64 - 128));
    ClpSimplex *lpSolver = clpSolver->getModelPtr();
    lpSolver->setSpecialOptions(lpSolver->specialOptions() | 0x01000000); // say is Cbc (and in branch and bound)
    lpSolver->setSpecialOptions(lpSolver->specialOptions() | (/*16384+*/ 4096 + 512 + 128));
  }
#endif
#ifdef HISTORY_STATISTICS
  getHistoryStatistics_ = false;
#endif
#ifdef COIN_DEVELOP
  int status = 0;
#endif
  int logLevel = model_->logLevel();
#define LEN_PRINT 250
  char generalPrint[LEN_PRINT];
  // Do presolve to see if possible
  int numberColumns = solver->getNumCols();
  char *reset = NULL;
  int returnCode = 1;
  int saveModelOptions = model_->specialOptions();
  //assert ((saveModelOptions&2048) == 0);
  model_->setSpecialOptions(saveModelOptions | 2048);
  if (fractionSmall < 1.0) {
    int saveLogLevel = solver->messageHandler()->logLevel();
    if (saveLogLevel == 1)
      solver->messageHandler()->setLogLevel(0);
    OsiPresolve *pinfo = new OsiPresolve();
    int presolveActions = 0;
    // Allow dual stuff on integers
    presolveActions = 1;
    // Do not allow all +1 to be tampered with
    //if (allPlusOnes)
    //presolveActions |= 2;
    // allow transfer of costs
    // presolveActions |= 4;
    pinfo->setPresolveActions(presolveActions);
    OsiSolverInterface *presolvedModel = pinfo->presolvedModel(*solver, 1.0e-8, true, 2);
    delete pinfo;
    // see if too big

    if (presolvedModel) {
      int afterRows = presolvedModel->getNumRows();
      int afterCols = presolvedModel->getNumCols();
      //#define COIN_DEVELOP
#ifdef COIN_DEVELOP_z
      if (numberNodes < 0) {
        solver->writeMpsNative("before.mps", NULL, NULL, 2, 1);
        presolvedModel->writeMpsNative("after1.mps", NULL, NULL, 2, 1);
      }
#endif
      delete presolvedModel;
      double ratio = sizeRatio(afterRows - shiftRows, afterCols,
        numberRowsStart, numberColumnsStart);
      double after = 2 * afterRows + afterCols;
      if (ratio > fractionSmall && after > 300 && numberNodes >= 0) {
        // Need code to try again to compress further using used
        const int *used = model_->usedInSolution();
        int maxUsed = 0;
        int iColumn;
        const double *lower = solver->getColLower();
        const double *upper = solver->getColUpper();
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (upper[iColumn] > lower[iColumn]) {
            if (solver->isBinary(iColumn))
              maxUsed = CoinMax(maxUsed, used[iColumn]);
          }
        }
        if (maxUsed) {
          reset = new char[numberColumns];
          int nFix = 0;
          for (iColumn = 0; iColumn < numberColumns; iColumn++) {
            reset[iColumn] = 0;
            if (upper[iColumn] > lower[iColumn]) {
              if (solver->isBinary(iColumn) && used[iColumn] == maxUsed) {
                bool setValue = true;
                if (maxUsed == 1) {
                  double randomNumber = randomNumberGenerator_.randomDouble();
                  if (randomNumber > 0.3)
                    setValue = false;
                }
                if (setValue) {
                  reset[iColumn] = 1;
                  solver->setColLower(iColumn, 1.0);
                  nFix++;
                }
              }
            }
          }
          pinfo = new OsiPresolve();
          presolveActions = 0;
          // Allow dual stuff on integers
          presolveActions = 1;
          // Do not allow all +1 to be tampered with
          //if (allPlusOnes)
          //presolveActions |= 2;
          // allow transfer of costs
          // presolveActions |= 4;
          pinfo->setPresolveActions(presolveActions);
          presolvedModel = pinfo->presolvedModel(*solver, 1.0e-8, true, 2);
          delete pinfo;
          if (presolvedModel) {
            // see if too big
            int afterRows2 = presolvedModel->getNumRows();
            int afterCols2 = presolvedModel->getNumCols();
            delete presolvedModel;
            double ratio = sizeRatio(afterRows2 - shiftRows, afterCols2,
              numberRowsStart, numberColumnsStart);
            double after = 2 * afterRows2 + afterCols2;
            if (ratio > fractionSmall && (after > 300 || numberNodes < 0)) {
              sprintf(generalPrint, "Full problem %d rows %d columns, reduced to %d rows %d columns - %d fixed gives %d, %d - still too large",
                solver->getNumRows(), solver->getNumCols(),
                afterRows, afterCols, nFix, afterRows2, afterCols2);
              // If much too big - give up
              if (ratio > 0.75)
                returnCode = -1;
            } else {
              sprintf(generalPrint, "Full problem %d rows %d columns, reduced to %d rows %d columns - %d fixed gives %d, %d - ok now",
                solver->getNumRows(), solver->getNumCols(),
                afterRows, afterCols, nFix, afterRows2, afterCols2);
            }
            model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
              << generalPrint
              << CoinMessageEol;
          } else {
            returnCode = 2; // infeasible
          }
        }
      } else if (ratio > fractionSmall && after > 300 && numberNodes >= 0) {
        returnCode = -1;
      }
    } else {
      returnCode = 2; // infeasible
    }
    solver->messageHandler()->setLogLevel(saveLogLevel);
  }
  if (returnCode == 2 || returnCode == -1) {
    model_->setSpecialOptions(saveModelOptions);
    delete[] reset;
#ifdef HISTORY_STATISTICS
    getHistoryStatistics_ = true;
#endif
    //printf("small no good\n");
    return returnCode;
  }
  // Reduce printout
  bool takeHint;
  OsiHintStrength strength;
  solver->getHintParam(OsiDoReducePrint, takeHint, strength);
  solver->setHintParam(OsiDoReducePrint, true, OsiHintTry);
  solver->setHintParam(OsiDoPresolveInInitial, false, OsiHintTry);
  double signedCutoff = cutoff * solver->getObjSense();
  solver->setDblParam(OsiDualObjectiveLimit, signedCutoff);
  solver->initialSolve();
  if (solver->isProvenOptimal()) {
    CglPreProcess process;
    OsiSolverInterface *solver2 = NULL;
    if ((model_->moreSpecialOptions() & 65536) != 0)
      process.setOptions(2 + 4 + 8 + 16); // no cuts
    else
      process.setOptions(16); // no complicated dupcol stuff
    /* Do not try and produce equality cliques and
	   do up to 2 passes (normally) 5 if restart */
    int numberPasses = 2;
    if ((model_->moreSpecialOptions2() & 16) != 0) {
      // quick
      process.setOptions(2 + 4 + 8 + 16); // no cuts
      numberPasses = 1;
    }
    if (numberNodes < 0) {
      numberPasses = 5;
      // Say some rows cuts
      int numberRows = solver->getNumRows();
      if (numberNodes_ < numberRows && true /* think */) {
        char *type = new char[numberRows];
        memset(type, 0, numberNodes_);
        memset(type + numberNodes_, 1, numberRows - numberNodes_);
        process.passInRowTypes(type, numberRows);
        delete[] type;
      }
    }
    if (logLevel <= 1)
      process.messageHandler()->setLogLevel(0);
    if (!solver->defaultHandler() && solver->messageHandler()->logLevel(0) != -1000)
      process.passInMessageHandler(solver->messageHandler());
#ifdef CGL_DEBUG
    /*
	  We're debugging. (specialOptions 1)
	*/
    if ((model_->specialOptions() & 1) != 0) {
      const OsiRowCutDebugger *debugger = solver->getRowCutDebugger();
      if (debugger) {
        process.setApplicationData(const_cast< double * >(debugger->optimalSolution()));
      }
    }
#endif
#ifdef COIN_HAS_CLP
    OsiClpSolverInterface *clpSolver = dynamic_cast< OsiClpSolverInterface * >(solver);
    // See if SOS
    if (clpSolver && clpSolver->numberSOS()) {
      // SOS
      int numberSOS = clpSolver->numberSOS();
      const CoinSet *setInfo = clpSolver->setInfo();
      int *sosStart = new int[numberSOS + 1];
      char *sosType = new char[numberSOS];
      int i;
      int nTotal = 0;
      sosStart[0] = 0;
      for (i = 0; i < numberSOS; i++) {
        int type = setInfo[i].setType();
        int n = setInfo[i].numberEntries();
        sosType[i] = static_cast< char >(type);
        nTotal += n;
        sosStart[i + 1] = nTotal;
      }
      int *sosIndices = new int[nTotal];
      double *sosReference = new double[nTotal];
      for (i = 0; i < numberSOS; i++) {
        int n = setInfo[i].numberEntries();
        const int *which = setInfo[i].which();
        const double *weights = setInfo[i].weights();
        int base = sosStart[i];
        for (int j = 0; j < n; j++) {
          int k = which[j];
          sosIndices[j + base] = k;
          sosReference[j + base] = weights ? weights[j] : static_cast< double >(j);
        }
      }
      int numberColumns = solver->getNumCols();
      char *prohibited = new char[numberColumns];
      memset(prohibited, 0, numberColumns);
      int n = sosStart[numberSOS];
      for (int i = 0; i < n; i++) {
        int iColumn = sosIndices[i];
        prohibited[iColumn] = 1;
      }
      delete[] sosIndices;
      delete[] sosReference;
      delete[] sosStart;
      delete[] sosType;
      process.passInProhibited(prohibited, numberColumns);
      delete[] prohibited;
    }
#endif
    solver2 = process.preProcessNonDefault(*solver, false,
      numberPasses);
    if (!solver2) {
      if (logLevel > 1)
        printf("Pre-processing says infeasible\n");
      returnCode = 2; // so will be infeasible
    } else {
#ifdef COIN_DEVELOP_z
      if (numberNodes < 0) {
        solver2->writeMpsNative("after2.mps", NULL, NULL, 2, 1);
      }
#endif
      // see if too big
      double ratio = sizeRatio(solver2->getNumRows() - shiftRows, solver2->getNumCols(),
        numberRowsStart, numberColumnsStart);
      double after = 2 * solver2->getNumRows() + solver2->getNumCols();
      if (ratio > fractionSmall && (after > 300 || numberNodes < 0)) {
        sprintf(generalPrint, "Full problem %d rows %d columns, reduced to %d rows %d columns - too large",
          solver->getNumRows(), solver->getNumCols(),
          solver2->getNumRows(), solver2->getNumCols());
        model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
          << generalPrint
          << CoinMessageEol;
        returnCode = -1;
        //printf("small no good2\n");
      } else {
        sprintf(generalPrint, "Full problem %d rows %d columns, reduced to %d rows %d columns",
          solver->getNumRows(), solver->getNumCols(),
          solver2->getNumRows(), solver2->getNumCols());
        model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
          << generalPrint
          << CoinMessageEol;
      }
#ifdef CGL_DEBUG
      if ((model_->specialOptions() & 1) != 0) {
        const OsiRowCutDebugger *debugger = solver2->getRowCutDebugger();
        if (debugger) {
          printf("On optimal path after preprocessing\n");
        }
      }
#endif
      if (returnCode == 1) {
        solver2->resolve();
        CbcModel model(*solver2);
        double startTime = model_->getDblParam(CbcModel::CbcStartSeconds);
        model.setDblParam(CbcModel::CbcStartSeconds, startTime);
        // move seed across
        model.randomNumberGenerator()->setSeed(model_->randomNumberGenerator()->getSeed());
#ifdef COIN_HAS_CLP
        // redo SOS
        OsiClpSolverInterface *clpSolver
          = dynamic_cast< OsiClpSolverInterface * >(model.solver());
        if (clpSolver && clpSolver->numberSOS()) {
          int numberColumns = clpSolver->getNumCols();
          const int *originalColumns = process.originalColumns();
          CoinSet *setInfo = const_cast< CoinSet * >(clpSolver->setInfo());
          int numberSOS = clpSolver->numberSOS();
          for (int iSOS = 0; iSOS < numberSOS; iSOS++) {
            //int type = setInfo[iSOS].setType();
            int n = setInfo[iSOS].numberEntries();
            int *which = setInfo[iSOS].modifiableWhich();
            double *weights = setInfo[iSOS].modifiableWeights();
            int n2 = 0;
            for (int j = 0; j < n; j++) {
              int iColumn = which[j];
              int i;
              for (i = 0; i < numberColumns; i++) {
                if (originalColumns[i] == iColumn)
                  break;
              }
              if (i < numberColumns) {
                which[n2] = i;
                weights[n2++] = weights[j];
              }
            }
            setInfo[iSOS].setNumberEntries(n2);
          }
        }
#endif
        if (numberNodes >= 0) {
          // normal
          model.setSpecialOptions(saveModelOptions | 2048);
          if (logLevel <= 1 && feasibilityPumpOptions_ != -3)
            model.setLogLevel(0);
          else
            model.setLogLevel(logLevel);
          // No small fathoming
          model.setFastNodeDepth(-1);
          model.setCutoff(signedCutoff);
          model.setStrongStrategy(0);
          // Don't do if original fraction > 1.0 and too large
          if (fractionSmall_ > 1.0 && fractionSmall_ < 1000000.0) {
            /* 1.4 means -1 nodes if >.4
			 2.4 means -1 nodes if >.5 and 0 otherwise
			 3.4 means -1 nodes if >.6 and 0 or 5
			 4.4 means -1 nodes if >.7 and 0, 5 or 10
		      */
            double fraction = fractionSmall_ - floor(fractionSmall_);
            if (ratio > fraction) {
              int type = static_cast< int >(floor(fractionSmall_ * 0.1));
              int over = static_cast< int >(ceil(ratio - fraction));
              int maxNodes[] = { -1, 0, 5, 10 };
              if (type > over)
                numberNodes = maxNodes[type - over];
              else
                numberNodes = -1;
            }
          }
          model.setMaximumNodes(numberNodes);
          model.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
          if ((saveModelOptions & 2048) == 0)
            model.setMoreSpecialOptions(model_->moreSpecialOptions());
          model.setMoreSpecialOptions2(model_->moreSpecialOptions2());
          // off conflict analysis
          model.setMoreSpecialOptions(model.moreSpecialOptions() & ~4194304);

          // Lightweight
          CbcStrategyDefaultSubTree strategy(model_, 1, 5, 1, 0);
          model.setStrategy(strategy);
          model.solver()->setIntParam(OsiMaxNumIterationHotStart, 10);
          model.setMaximumCutPassesAtRoot(CoinMin(20, CoinAbs(model_->getMaximumCutPassesAtRoot())));
          model.setMaximumCutPasses(CoinMin(10, model_->getMaximumCutPasses()));
          // Set best solution (even if bad for this submodel)
          if (model_->bestSolution()) {
            const double *bestSolution = model_->bestSolution();
            int numberColumns2 = model.solver()->getNumCols();
            double *bestSolution2 = new double[numberColumns2];
            const int *originalColumns = process.originalColumns();
            for (int iColumn = 0; iColumn < numberColumns2; iColumn++) {
              int jColumn = originalColumns[iColumn];
              bestSolution2[iColumn] = bestSolution[jColumn];
            }
            model.setBestSolution(bestSolution2, numberColumns2,
              1.0e50,
              false);
            model.setSolutionCount(1);
            maximumSolutions++;
            delete[] bestSolution2;
          }
        } else {
          // modify for event handler
          model.setSpecialOptions(saveModelOptions);
          model_->messageHandler()->message(CBC_RESTART, model_->messages())
            << solver2->getNumRows() << solver2->getNumCols()
            << CoinMessageEol;
          // going for full search and copy across more stuff
          model.gutsOfCopy(*model_, 2);
#ifdef CGL_DEBUG
          if ((model_->specialOptions() & 1) != 0) {
            const OsiRowCutDebugger *debugger = model.solver()->getRowCutDebugger();
            if (debugger) {
              printf("On optimal path BB\n");
            }
          }
#endif
          assert(!model_->heuristicModel());
          model_->setHeuristicModel(&model);
          for (int i = 0; i < model.numberCutGenerators(); i++) {
            CbcCutGenerator *generator = model.cutGenerator(i);
            CglGomory *gomory = dynamic_cast< CglGomory * >(generator->generator());
            if (gomory && gomory->originalSolver())
              gomory->passInOriginalSolver(model.solver());
            generator->setTiming(true);
            // Turn on if was turned on
            int iOften = model_->cutGenerator(i)->howOften();
#ifdef CLP_INVESTIGATE
            printf("Gen %d often %d %d\n",
              i, generator->howOften(),
              iOften);
#endif
            if (iOften > 0)
              generator->setHowOften(iOften % 1000000);
            if (model_->cutGenerator(i)->howOftenInSub() == -200)
              generator->setHowOften(-100);
          }
          model.setCutoff(signedCutoff);
          // make sure can't do nested search! but allow heuristics
          model.setSpecialOptions((model.specialOptions() & (~(512 + 2048))) | 1024);
          // but say we are doing full search
          model.setSpecialOptions(model.specialOptions() | 67108864);
          bool takeHint;
          OsiHintStrength strength;
          // Switch off printing if asked to
          model_->solver()->getHintParam(OsiDoReducePrint, takeHint, strength);
          model.solver()->setHintParam(OsiDoReducePrint, takeHint, strength);
          // no cut generators if none in parent
          CbcStrategyDefault
            strategy(model_->numberCutGenerators() ? 1 : -1,
              model_->numberStrong(),
              model_->numberBeforeTrust());
          // Set up pre-processing - no
          strategy.setupPreProcessing(0); // was (4);
          model.setStrategy(strategy);
          //model.solver()->writeMps("crunched");
          int numberCuts = process.cuts().sizeRowCuts();
          if (numberCuts) {
            // add in cuts
            CglStored cuts = process.cuts();
            model.addCutGenerator(&cuts, 1, "Stored from first");
            model.cutGenerator(model.numberCutGenerators() - 1)->setGlobalCuts(true);
          }
        }
        // Do search
        if (logLevel > 1)
          model_->messageHandler()->message(CBC_START_SUB, model_->messages())
            << name
            << model.getMaximumNodes()
            << CoinMessageEol;
        // probably faster to use a basis to get integer solutions
        model.setSpecialOptions(model.specialOptions() | 2);
#ifdef CBC_THREAD
        if (model_->getNumberThreads() > 0 && (model_->getThreadMode() & 4) != 0) {
          // See if at root node
          bool atRoot = model_->getNodeCount() == 0;
          int passNumber = model_->getCurrentPassNumber();
          if (atRoot && passNumber == 1)
            model.setNumberThreads(model_->getNumberThreads());
        }
#endif
        model.setParentModel(*model_);
        model.setMaximumSolutions(maximumSolutions);
        model.setOriginalColumns(process.originalColumns());
        model.setSearchStrategy(-1);
        // If no feasibility pump then insert a lightweight one
        if (feasibilityPumpOptions_ >= 0 || feasibilityPumpOptions_ == -2) {
          CbcHeuristicFPump *fpump = NULL;
          for (int i = 0; i < model.numberHeuristics(); i++) {
            CbcHeuristicFPump *pump = dynamic_cast< CbcHeuristicFPump * >(model.heuristic(i));
            if (pump) {
              fpump = pump;
              break;
            }
          }
          if (!fpump) {
            CbcHeuristicFPump heuristic4;
            // use any cutoff
            heuristic4.setFakeCutoff(0.5 * COIN_DBL_MAX);
            if (fractionSmall_ <= 1.0)
              heuristic4.setMaximumPasses(10);
            int pumpTune = feasibilityPumpOptions_;
            if (pumpTune == -2)
              pumpTune = 4; // proximity
            if (pumpTune > 0) {
              /*
                            >=10000000 for using obj
                            >=1000000 use as accumulate switch
                            >=1000 use index+1 as number of large loops
                            >=100 use 0.05 objvalue as increment
                            %100 == 10,20 etc for experimentation
                            1 == fix ints at bounds, 2 fix all integral ints, 3 and continuous at bounds
                            4 and static continuous, 5 as 3 but no internal integers
                            6 as 3 but all slack basis!
                            */
              double value = solver2->getObjSense() * solver2->getObjValue();
              int w = pumpTune / 10;
              int ix = w % 10;
              w /= 10;
              int c = w % 10;
              w /= 10;
              int r = w;
              int accumulate = r / 1000;
              r -= 1000 * accumulate;
              if (accumulate >= 10) {
                int which = accumulate / 10;
                accumulate -= 10 * which;
                which--;
                // weights and factors
                double weight[] = { 0.1, 0.1, 0.5, 0.5, 1.0, 1.0, 5.0, 5.0 };
                double factor[] = { 0.1, 0.5, 0.1, 0.5, 0.1, 0.5, 0.1, 0.5 };
                heuristic4.setInitialWeight(weight[which]);
                heuristic4.setWeightFactor(factor[which]);
              }
              // fake cutoff
              if (c) {
                double cutoff;
                solver2->getDblParam(OsiDualObjectiveLimit, cutoff);
                cutoff = CoinMin(cutoff, value + 0.1 * fabs(value) * c);
                heuristic4.setFakeCutoff(cutoff);
              }
              if (r) {
                // also set increment
                //double increment = (0.01*i+0.005)*(fabs(value)+1.0e-12);
                double increment = 0.0;
                heuristic4.setAbsoluteIncrement(increment);
                heuristic4.setAccumulate(accumulate);
                heuristic4.setMaximumRetries(r + 1);
              }
              pumpTune = pumpTune % 100;
              if (pumpTune == 6)
                pumpTune = 13;
              if (pumpTune != 13)
                pumpTune = pumpTune % 10;
              heuristic4.setWhen(pumpTune);
              if (ix) {
                heuristic4.setFeasibilityPumpOptions(ix * 10);
              }
            }
            model.addHeuristic(&heuristic4, "feasibility pump", 0);
          }
        } else if (feasibilityPumpOptions_ == -3) {
          // add all (except this)
          for (int i = 0; i < model_->numberHeuristics(); i++) {
            if (strcmp(heuristicName(), model_->heuristic(i)->heuristicName()))
              model.addHeuristic(model_->heuristic(i));
          }
        }
        // modify heuristics
        for (int i = 0; i < model.numberHeuristics(); i++) {
          // reset lastNode
          CbcHeuristicRINS *rins = dynamic_cast< CbcHeuristicRINS * >(model.heuristic(i));
          if (rins) {
            rins->setLastNode(-1000);
            rins->setSolutionCount(0);
          }
        }
        //printf("sol %x\n",inputSolution_);
#ifdef CGL_DEBUG
        if ((model_->specialOptions() & 1) != 0) {
          const OsiRowCutDebugger *debugger = model.solver()->getRowCutDebugger();
          if (debugger) {
            printf("On optimal path CC\n");
          }
        }
#endif
        if (inputSolution_) {
          // translate and add a serendipity heuristic
          int numberColumns = solver2->getNumCols();
          const int *which = process.originalColumns();
          OsiSolverInterface *solver3 = solver2->clone();
          for (int i = 0; i < numberColumns; i++) {
            if (isHeuristicInteger(solver3, i)) {
              int k = which[i];
              double value = inputSolution_[k];
              //if (value)
              //printf("orig col %d now %d val %g\n",
              //       k,i,value);
              solver3->setColLower(i, value);
              solver3->setColUpper(i, value);
            }
          }
          solver3->setDblParam(OsiDualObjectiveLimit, COIN_DBL_MAX);
          solver3->resolve();
          if (!solver3->isProvenOptimal()) {
            // Try just setting nonzeros
            OsiSolverInterface *solver4 = solver2->clone();
            for (int i = 0; i < numberColumns; i++) {
              if (isHeuristicInteger(solver4, i)) {
                int k = which[i];
                double value = floor(inputSolution_[k] + 0.5);
                if (value) {
                  solver3->setColLower(i, value);
                  solver3->setColUpper(i, value);
                }
              }
            }
            solver4->setDblParam(OsiDualObjectiveLimit, COIN_DBL_MAX);
            solver4->resolve();
            int nBad = -1;
            if (solver4->isProvenOptimal()) {
              nBad = 0;
              const double *solution = solver4->getColSolution();
              for (int i = 0; i < numberColumns; i++) {
                if (isHeuristicInteger(solver4, i)) {
                  double value = floor(solution[i] + 0.5);
                  if (fabs(value - solution[i]) > 1.0e-6)
                    nBad++;
                }
              }
            }
            if (nBad) {
              delete solver4;
            } else {
              delete solver3;
              solver3 = solver4;
            }
          }
          if (solver3->isProvenOptimal()) {
            // good
            CbcSerendipity heuristic(model);
            double value = solver3->getObjSense() * solver3->getObjValue();
            heuristic.setInputSolution(solver3->getColSolution(), value);
            value = value + 1.0e-7 * (1.0 + fabs(value));
            value *= solver3->getObjSense();
            model.setCutoff(value);
            model.addHeuristic(&heuristic, "Previous solution", 0);
            //printf("added seren\n");
          } else {
            double value = model_->getMinimizationObjValue();
            value = value + 1.0e-7 * (1.0 + fabs(value));
            value *= solver3->getObjSense();
            model.setCutoff(value);
            sprintf(generalPrint, "Unable to insert previous solution - using cutoff of %g",
		    trueObjValue(value));
            model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
              << generalPrint
              << CoinMessageEol;
#ifdef CLP_INVESTIGATE
            printf("NOT added seren\n");
            solver3->writeMps("bad_seren");
            solver->writeMps("orig_seren");
#endif
          }
          delete solver3;
        }
        if (model_->searchStrategy() == 2) {
          model.setNumberStrong(5);
          model.setNumberBeforeTrust(5);
        }
        if (model.getNumCols()) {
          if (numberNodes >= 0) {
            setCutAndHeuristicOptions(model);
            // not too many iterations
            model.setMaximumNumberIterations(iterationMultiplier * (numberNodes + 10));
            // Not fast stuff
            model.setFastNodeDepth(-1);
            //model.solver()->writeMps("before");
          } else if (model.fastNodeDepth() >= 1000000) {
            // already set
            model.setFastNodeDepth(model.fastNodeDepth() - 1000000);
          }
          model.setWhenCuts(999998);
#define ALWAYS_DUAL
#ifdef ALWAYS_DUAL
          OsiSolverInterface *solverD = model.solver();
          bool takeHint;
          OsiHintStrength strength;
          solverD->getHintParam(OsiDoDualInResolve, takeHint, strength);
          solverD->setHintParam(OsiDoDualInResolve, true, OsiHintDo);
#endif
          model.passInEventHandler(model_->getEventHandler());
          // say model_ is sitting there
          int saveOptions = model_->specialOptions();
          model_->setSpecialOptions(saveOptions | 1048576);
          // and switch off debugger
          model.setSpecialOptions(model.specialOptions() & (~1));
#if 0 //def COIN_HAS_CLP
		    OsiClpSolverInterface * clpSolver
		      = dynamic_cast<OsiClpSolverInterface *> (model.solver());
		    if (clpSolver)
		      clpSolver->zapDebugger();
#endif
#ifdef CONFLICT_CUTS
          if ((model_->moreSpecialOptions() & 4194304) != 0)
            model.zapGlobalCuts();
#endif
#ifdef CGL_DEBUG
          if ((model_->specialOptions() & 1) != 0) {
            const OsiRowCutDebugger *debugger = model.solver()->getRowCutDebugger();
            if (debugger) {
              printf("On optimal path DD\n");
            }
          }
#endif
          model.setPreProcess(&process);
          model.branchAndBound();
          model_->setHeuristicModel(NULL);
          model_->setSpecialOptions(saveOptions);
#ifdef ALWAYS_DUAL
          solverD = model.solver();
          solverD->setHintParam(OsiDoDualInResolve, takeHint, strength);
#endif
          numberNodesDone_ = model.getNodeCount();
#ifdef COIN_DEVELOP
          printf("sub branch %d nodes, %d iterations - max %d\n",
            model.getNodeCount(), model.getIterationCount(),
            100 * (numberNodes + 10));
#endif
          if (numberNodes < 0) {
            model_->incrementIterationCount(model.getIterationCount());
            model_->incrementNodeCount(model.getNodeCount());
            // update best solution (in case ctrl-c)
            // !!! not a good idea - think a bit harder
            //model_->setMinimizationObjValue(model.getMinimizationObjValue());
            for (int iGenerator = 0; iGenerator < model.numberCutGenerators(); iGenerator++) {
              CbcCutGenerator *generator = model.cutGenerator(iGenerator);
              sprintf(generalPrint,
                "%s was tried %d times and created %d cuts of which %d were active after adding rounds of cuts (%.3f seconds)",
                generator->cutGeneratorName(),
                generator->numberTimesEntered(),
                generator->numberCutsInTotal() + generator->numberColumnCuts(),
                generator->numberCutsActive(),
                generator->timeInCutGenerator());
              CglStored *stored = dynamic_cast< CglStored * >(generator->generator());
              if (stored && !generator->numberCutsInTotal())
                continue;
#ifndef CLP_INVESTIGATE
              CglImplication *implication = dynamic_cast< CglImplication * >(generator->generator());
              if (implication)
                continue;
#endif
              model_->messageHandler()->message(CBC_FPUMP1, model_->messages())
                << generalPrint
                << CoinMessageEol;
            }
          }
        } else {
          // empty model
          model.setMinimizationObjValue(model.solver()->getObjSense() * model.solver()->getObjValue());
        }
        if (logLevel > 1)
          model_->messageHandler()->message(CBC_END_SUB, model_->messages())
            << name
            << CoinMessageEol;
        if (model.getMinimizationObjValue() < CoinMin(cutoff, 1.0e30)) {
          // solution
          if (model.getNumCols())
            returnCode = model.isProvenOptimal() ? 3 : 1;
          else
            returnCode = 3;
            // post process
#ifdef COIN_HAS_CLP
          OsiClpSolverInterface *clpSolver = dynamic_cast< OsiClpSolverInterface * >(model.solver());
          if (clpSolver) {
            ClpSimplex *lpSolver = clpSolver->getModelPtr();
            lpSolver->setSpecialOptions(lpSolver->specialOptions() | 0x01000000); // say is Cbc (and in branch and bound)
          }
#endif
          //if (fractionSmall_ < 1000000.0)
          process.postProcess(*model.solver());
          if (solver->isProvenOptimal() && solver->getObjValue() * solver->getObjSense() < cutoff) {
            // Solution now back in solver
            int numberColumns = solver->getNumCols();
            memcpy(newSolution, solver->getColSolution(),
              numberColumns * sizeof(double));
            newSolutionValue = model.getMinimizationObjValue();
#ifdef COIN_HAS_CLP
            if (clpSolver) {
              if (clpSolver && clpSolver->numberSOS()) {
                // SOS
                int numberSOS = clpSolver->numberSOS();
                const CoinSet *setInfo = clpSolver->setInfo();
                int i;
                for (i = 0; i < numberSOS; i++) {
#ifndef NDEBUG
                  int type = setInfo[i].setType();
#endif
                  int n = setInfo[i].numberEntries();
                  const int *which = setInfo[i].which();
                  int first = -1;
                  int last = -1;
                  for (int j = 0; j < n; j++) {
                    int iColumn = which[j];
                    if (fabs(newSolution[iColumn]) > 1.0e-7) {
                      last = j;
                      if (first < 0)
                        first = j;
                    }
                  }
                  assert(last - first < type);
                  for (int j = 0; j < n; j++) {
                    if (j < first || j > last) {
                      int iColumn = which[j];
                      // do I want to fix??
                      solver->setColLower(iColumn, 0.0);
                      solver->setColUpper(iColumn, 0.0);
                      newSolution[iColumn] = 0.0;
                    }
                  }
                }
              }
            }
#endif
          } else {
            // odd - but no good
            returnCode = 0; // so will be infeasible
          }
        } else {
          // no good
          returnCode = model.isProvenInfeasible() ? 2 : 0; // so will be infeasible
        }
        int totalNumberIterations = model.getIterationCount() + process.numberIterationsPre() + process.numberIterationsPost();
        if (totalNumberIterations > 100 * (numberNodes + 10)
          && fractionSmall_ < 1000000.0) {
          // only allow smaller problems
          fractionSmall = fractionSmall_;
          fractionSmall_ *= 0.9;
#ifdef CLP_INVESTIGATE
          printf("changing fractionSmall from %g to %g for %s as %d iterations\n",
            fractionSmall, fractionSmall_, name.c_str(), totalNumberIterations);
#endif
        }
        if (model.status() == 5)
          model_->sayEventHappened();
#ifdef COIN_DEVELOP
        if (model.isProvenInfeasible())
          status = 1;
        else if (model.isProvenOptimal())
          status = 2;
#endif
      }
    }
  } else {
    returnCode = 2; // infeasible finished
    if (logLevel > 1) {
      printf("Infeasible on initial solve\n");
    }
  }
  model_->setSpecialOptions(saveModelOptions);
  model_->setLogLevel(logLevel);
  if (returnCode == 1 || returnCode == 2) {
    OsiSolverInterface *solverC = model_->continuousSolver();
    if (false && solverC) {
      const double *lower = solver->getColLower();
      const double *upper = solver->getColUpper();
      const double *lowerC = solverC->getColLower();
      const double *upperC = solverC->getColUpper();
      bool good = true;
      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
        if (isHeuristicInteger(solverC, iColumn)) {
          if (lower[iColumn] > lowerC[iColumn] && upper[iColumn] < upperC[iColumn]) {
            good = false;
            printf("CUT - can't add\n");
            break;
          }
        }
      }
      if (good) {
        double *cut = new double[numberColumns];
        int *which = new int[numberColumns];
        double rhs = -1.0;
        int n = 0;
        for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (isHeuristicInteger(solverC, iColumn)) {
            if (lower[iColumn] == upperC[iColumn]) {
              rhs += lower[iColumn];
              cut[n] = 1.0;
              which[n++] = iColumn;
            } else if (upper[iColumn] == lowerC[iColumn]) {
              rhs -= upper[iColumn];
              cut[n] = -1.0;
              which[n++] = iColumn;
            }
          }
        }
        printf("CUT has %d entries\n", n);
        OsiRowCut newCut;
        newCut.setLb(-COIN_DBL_MAX);
        newCut.setUb(rhs);
        newCut.setRow(n, which, cut, false);
        model_->makeGlobalCut(newCut);
        delete[] cut;
        delete[] which;
      }
    }
#ifdef COIN_DEVELOP
    if (status == 1)
      printf("heuristic could add cut because infeasible (%s)\n", heuristicName_.c_str());
    else if (status == 2)
      printf("heuristic could add cut because optimal (%s)\n", heuristicName_.c_str());
#endif
  }
  if (reset) {
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (reset[iColumn])
        solver->setColLower(iColumn, 0.0);
    }
    delete[] reset;
  }
#ifdef HISTORY_STATISTICS
  getHistoryStatistics_ = true;
#endif
  solver->setHintParam(OsiDoReducePrint, takeHint, strength);
  return returnCode;
}
// Set input solution
void CbcHeuristic::setInputSolution(const double *solution, double objValue)
{
  delete[] inputSolution_;
  inputSolution_ = NULL;
  if (model_ && solution) {
    int numberColumns = model_->getNumCols();
    inputSolution_ = new double[numberColumns + 1];
    memcpy(inputSolution_, solution, numberColumns * sizeof(double));
    inputSolution_[numberColumns] = objValue;
  }
}

//##############################################################################

inline int compare3BranchingObjects(const CbcBranchingObject *br0,
  const CbcBranchingObject *br1)
{
  const int t0 = br0->type();
  const int t1 = br1->type();
  if (t0 < t1) {
    return -1;
  }
  if (t0 > t1) {
    return 1;
  }
  return br0->compareOriginalObject(br1);
}

//==============================================================================

inline bool compareBranchingObjects(const CbcBranchingObject *br0,
  const CbcBranchingObject *br1)
{
  return compare3BranchingObjects(br0, br1) < 0;
}

//==============================================================================

void CbcHeuristicNode::gutsOfConstructor(CbcModel &model)
{
  //  CbcHeurDebugNodes(&model);
  CbcNode *node = model.currentNode();
  brObj_ = new CbcBranchingObject *[node->depth()];
  CbcNodeInfo *nodeInfo = node->nodeInfo();
  int cnt = 0;
  while (nodeInfo->parentBranch() != NULL) {
    const OsiBranchingObject *br = nodeInfo->parentBranch();
    const CbcBranchingObject *cbcbr = dynamic_cast< const CbcBranchingObject * >(br);
    if (!cbcbr) {
      throw CoinError("CbcHeuristicNode can be used only with CbcBranchingObjects.\n",
        "gutsOfConstructor",
        "CbcHeuristicNode",
        __FILE__, __LINE__);
    }
    brObj_[cnt] = cbcbr->clone();
    brObj_[cnt]->previousBranch();
    ++cnt;
    nodeInfo = nodeInfo->parent();
  }
  std::sort(brObj_, brObj_ + cnt, compareBranchingObjects);
  if (cnt <= 1) {
    numObjects_ = cnt;
  } else {
    numObjects_ = 0;
    CbcBranchingObject *br = NULL; // What should this be?
    for (int i = 1; i < cnt; ++i) {
      if (compare3BranchingObjects(brObj_[numObjects_], brObj_[i]) == 0) {
        int comp = brObj_[numObjects_]->compareBranchingObject(brObj_[i], br != 0);
        switch (comp) {
        case CbcRangeSame: // the same range
        case CbcRangeDisjoint: // disjoint decisions
          // should not happen! we are on a chain!
          abort();
        case CbcRangeSubset: // brObj_[numObjects_] is a subset of brObj_[i]
          delete brObj_[i];
          break;
        case CbcRangeSuperset: // brObj_[i] is a subset of brObj_[numObjects_]
          delete brObj_[numObjects_];
          brObj_[numObjects_] = brObj_[i];
          break;
        case CbcRangeOverlap: // overlap
          delete brObj_[i];
          delete brObj_[numObjects_];
          brObj_[numObjects_] = br;
          break;
        }
        continue;
      } else {
        brObj_[++numObjects_] = brObj_[i];
      }
    }
    ++numObjects_;
  }
}

//==============================================================================

CbcHeuristicNode::CbcHeuristicNode(CbcModel &model)
{
  gutsOfConstructor(model);
}

//==============================================================================

double
CbcHeuristicNode::distance(const CbcHeuristicNode *node) const
{

  const double disjointWeight = 1;
  const double overlapWeight = 0.4;
  const double subsetWeight = 0.2;
#ifdef COIN_DETAIL
  int countDisjointWeight = 0;
  int countOverlapWeight = 0;
  int countSubsetWeight = 0;
#endif
  int i = 0;
  int j = 0;
  double dist = 0.0;
#ifdef PRINT_DEBUG
  printf(" numObjects_ = %i, node->numObjects_ = %i\n",
    numObjects_, node->numObjects_);
#endif
  while (i < numObjects_ && j < node->numObjects_) {
    CbcBranchingObject *br0 = brObj_[i];
    const CbcBranchingObject *br1 = node->brObj_[j];
#ifdef PRINT_DEBUG
    const CbcIntegerBranchingObject *brPrint0 = dynamic_cast< const CbcIntegerBranchingObject * >(br0);
    const double *downBounds = brPrint0->downBounds();
    const double *upBounds = brPrint0->upBounds();
    int variable = brPrint0->variable();
    int way = brPrint0->way();
    printf("   br0: var %i downBd [%i,%i] upBd [%i,%i] way %i\n",
      variable, static_cast< int >(downBounds[0]), static_cast< int >(downBounds[1]),
      static_cast< int >(upBounds[0]), static_cast< int >(upBounds[1]), way);
    const CbcIntegerBranchingObject *brPrint1 = dynamic_cast< const CbcIntegerBranchingObject * >(br1);
    downBounds = brPrint1->downBounds();
    upBounds = brPrint1->upBounds();
    variable = brPrint1->variable();
    way = brPrint1->way();
    printf("   br1: var %i downBd [%i,%i] upBd [%i,%i] way %i\n",
      variable, static_cast< int >(downBounds[0]), static_cast< int >(downBounds[1]),
      static_cast< int >(upBounds[0]), static_cast< int >(upBounds[1]), way);
#endif
    const int brComp = compare3BranchingObjects(br0, br1);
    if (brComp < 0) {
      dist += subsetWeight;
#ifdef COIN_DETAIL
      countSubsetWeight++;
#endif
      ++i;
    } else if (brComp > 0) {
      dist += subsetWeight;
#ifdef COIN_DETAIL
      countSubsetWeight++;
#endif
      ++j;
    } else {
      const int comp = br0->compareBranchingObject(br1, false);
      switch (comp) {
      case CbcRangeSame:
        // do nothing
        break;
      case CbcRangeDisjoint: // disjoint decisions
        dist += disjointWeight;
#ifdef COIN_DETAIL
        countDisjointWeight++;
#endif
        break;
      case CbcRangeSubset: // subset one way or another
      case CbcRangeSuperset:
        dist += subsetWeight;
#ifdef COIN_DETAIL
        countSubsetWeight++;
#endif
        break;
      case CbcRangeOverlap: // overlap
        dist += overlapWeight;
#ifdef COIN_DETAIL
        countOverlapWeight++;
#endif
        break;
      }
      ++i;
      ++j;
    }
  }
  dist += subsetWeight * (numObjects_ - i + node->numObjects_ - j);
#ifdef COIN_DETAIL
  countSubsetWeight += (numObjects_ - i + node->numObjects_ - j);
#endif
  COIN_DETAIL_PRINT(printf("subset = %i, overlap = %i, disjoint = %i\n", countSubsetWeight,
    countOverlapWeight, countDisjointWeight));
  return dist;
}

//==============================================================================

CbcHeuristicNode::~CbcHeuristicNode()
{
  for (int i = 0; i < numObjects_; ++i) {
    delete brObj_[i];
  }
  delete[] brObj_;
}

//==============================================================================

double
CbcHeuristicNode::minDistance(const CbcHeuristicNodeList &nodeList) const
{
  double minDist = COIN_DBL_MAX;
  for (int i = nodeList.size() - 1; i >= 0; --i) {
    minDist = CoinMin(minDist, distance(nodeList.node(i)));
  }
  return minDist;
}

//==============================================================================

bool CbcHeuristicNode::minDistanceIsSmall(const CbcHeuristicNodeList &nodeList,
  const double threshold) const
{
  for (int i = nodeList.size() - 1; i >= 0; --i) {
    if (distance(nodeList.node(i)) >= threshold) {
      continue;
    } else {
      return true;
    }
  }
  return false;
}

//==============================================================================

double
CbcHeuristicNode::avgDistance(const CbcHeuristicNodeList &nodeList) const
{
  if (nodeList.size() == 0) {
    return COIN_DBL_MAX;
  }
  double sumDist = 0;
  for (int i = nodeList.size() - 1; i >= 0; --i) {
    sumDist += distance(nodeList.node(i));
  }
  return sumDist / nodeList.size();
}

//##############################################################################

// Default Constructor
CbcRounding::CbcRounding()
  : CbcHeuristic()
{
  // matrix and row copy will automatically be empty
  seed_ = 7654321;
  down_ = NULL;
  up_ = NULL;
  equal_ = NULL;
  //whereFrom_ |= 16*(1+256); // allow more often
}

// Constructor from model
CbcRounding::CbcRounding(CbcModel &model)
  : CbcHeuristic(model)
{
  // Get a copy of original matrix (and by row for rounding);
  assert(model.solver());
  if (model.solver()->getNumRows()) {
    matrix_ = *model.solver()->getMatrixByCol();
    matrixByRow_ = *model.solver()->getMatrixByRow();
    validate();
  }
  down_ = NULL;
  up_ = NULL;
  equal_ = NULL;
  seed_ = 7654321;
  //whereFrom_ |= 16*(1+256); // allow more often
}

// Destructor
CbcRounding::~CbcRounding()
{
  delete[] down_;
  delete[] up_;
  delete[] equal_;
}

// Clone
CbcHeuristic *
CbcRounding::clone() const
{
  return new CbcRounding(*this);
}
// Create C++ lines to get to current state
void CbcRounding::generateCpp(FILE *fp)
{
  CbcRounding other;
  fprintf(fp, "0#include \"CbcHeuristic.hpp\"\n");
  fprintf(fp, "3  CbcRounding rounding(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp, "rounding");
  if (seed_ != other.seed_)
    fprintf(fp, "3  rounding.setSeed(%d);\n", seed_);
  else
    fprintf(fp, "4  rounding.setSeed(%d);\n", seed_);
  fprintf(fp, "3  cbcModel->addHeuristic(&rounding);\n");
}
//#define NEW_ROUNDING
// Copy constructor
CbcRounding::CbcRounding(const CbcRounding &rhs)
  : CbcHeuristic(rhs)
  , matrix_(rhs.matrix_)
  , matrixByRow_(rhs.matrixByRow_)
  , seed_(rhs.seed_)
{
#ifdef NEW_ROUNDING
  int numberColumns = matrix_.getNumCols();
  down_ = CoinCopyOfArray(rhs.down_, numberColumns);
  up_ = CoinCopyOfArray(rhs.up_, numberColumns);
  equal_ = CoinCopyOfArray(rhs.equal_, numberColumns);
#else
  down_ = NULL;
  up_ = NULL;
  equal_ = NULL;
#endif
}

// Assignment operator
CbcRounding &
CbcRounding::operator=(const CbcRounding &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    matrix_ = rhs.matrix_;
    matrixByRow_ = rhs.matrixByRow_;
#ifdef NEW_ROUNDING
    delete[] down_;
    delete[] up_;
    delete[] equal_;
    int numberColumns = matrix_.getNumCols();
    down_ = CoinCopyOfArray(rhs.down_, numberColumns);
    up_ = CoinCopyOfArray(rhs.up_, numberColumns);
    equal_ = CoinCopyOfArray(rhs.equal_, numberColumns);
#else
    down_ = NULL;
    up_ = NULL;
    equal_ = NULL;
#endif
    seed_ = rhs.seed_;
  }
  return *this;
}

// Resets stuff if model changes
void CbcRounding::resetModel(CbcModel *model)
{
  model_ = model;
  // Get a copy of original matrix (and by row for rounding);
  assert(model_->solver());
  matrix_ = *model_->solver()->getMatrixByCol();
  matrixByRow_ = *model_->solver()->getMatrixByRow();
  validate();
}
/* Check whether the heuristic should run at all
   0 - before cuts at root node (or from doHeuristics)
   1 - during cuts at root
   2 - after root node cuts
   3 - after cuts at other nodes
   4 - during cuts at other nodes
   8 added if previous heuristic in loop found solution
*/
bool CbcRounding::shouldHeurRun(int whereFrom)
{
  if (whereFrom != 4) {
    return CbcHeuristic::shouldHeurRun(whereFrom);
  } else {
    numCouldRun_++;
    return shouldHeurRun_randomChoice();
  }
}
// See if rounding will give solution
// Sets value of solution
// Assumes rhs for original matrix still okay
// At present only works with integers
// Fix values if asked for
// Returns 1 if solution, 0 if not
int CbcRounding::solution(double &solutionValue,
  double *betterSolution)
{

  numCouldRun_++;
  // See if to do
  if (!when() || (when() % 10 == 1 && model_->phase() != 1) || (when() % 10 == 2 && (model_->phase() != 2 && model_->phase() != 3)))
    return 0; // switched off
  numRuns_++;
#ifdef HEURISTIC_INFORM
  printf("Entering heuristic %s - nRuns %d numCould %d when %d\n",
    heuristicName(), numRuns_, numCouldRun_, when_);
#endif
  OsiSolverInterface *solver = model_->solver();
  double direction = solver->getObjSense();
  double newSolutionValue = direction * solver->getObjValue();
  return solution(solutionValue, betterSolution, newSolutionValue);
}
// See if rounding will give solution
// Sets value of solution
// Assumes rhs for original matrix still okay
// At present only works with integers
// Fix values if asked for
// Returns 1 if solution, 0 if not
int CbcRounding::solution(double &solutionValue,
  double *betterSolution,
  double newSolutionValue)
{

  // See if to do
  if (!when() || (when() % 10 == 1 && model_->phase() != 1) || (when() % 10 == 2 && (model_->phase() != 2 && model_->phase() != 3)))
    return 0; // switched off
  OsiSolverInterface *solver = model_->solver();
  const double *lower = solver->getColLower();
  const double *upper = solver->getColUpper();
  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  const double *solution = solver->getColSolution();
  const double *objective = solver->getObjCoefficients();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance, primalTolerance);
  double useTolerance = primalTolerance;

  int numberRows = matrix_.getNumRows();
  assert(numberRows <= solver->getNumRows());
  if (numberRows == 0) {
    return 0;
  }
  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();
  int i;
  double direction = solver->getObjSense();
  //double newSolutionValue = direction*solver->getObjValue();
  int returnCode = 0;
  // Column copy
  const double *element = matrix_.getElements();
  const int *row = matrix_.getIndices();
  const CoinBigIndex *columnStart = matrix_.getVectorStarts();
  const int *columnLength = matrix_.getVectorLengths();
  // Row copy
  const double *elementByRow = matrixByRow_.getElements();
  const int *column = matrixByRow_.getIndices();
  const CoinBigIndex *rowStart = matrixByRow_.getVectorStarts();
  const int *rowLength = matrixByRow_.getVectorLengths();

  // Get solution array for heuristic solution
  int numberColumns = solver->getNumCols();
  double *newSolution = new double[numberColumns];
  memcpy(newSolution, solution, numberColumns * sizeof(double));

  double *rowActivity = new double[numberRows];
  memset(rowActivity, 0, numberRows * sizeof(double));
  for (i = 0; i < numberColumns; i++) {
    CoinBigIndex j;
    double value = newSolution[i];
    if (value < lower[i]) {
      value = lower[i];
      newSolution[i] = value;
    } else if (value > upper[i]) {
      value = upper[i];
      newSolution[i] = value;
    }
    if (value) {
      for (j = columnStart[i];
           j < columnStart[i] + columnLength[i]; j++) {
        int iRow = row[j];
        rowActivity[iRow] += value * element[j];
      }
    }
  }
  // check was feasible - if not adjust (cleaning may move)
  for (i = 0; i < numberRows; i++) {
    if (rowActivity[i] < rowLower[i]) {
      //assert (rowActivity[i]>rowLower[i]-1000.0*primalTolerance);
      rowActivity[i] = rowLower[i];
    } else if (rowActivity[i] > rowUpper[i]) {
      //assert (rowActivity[i]<rowUpper[i]+1000.0*primalTolerance);
      rowActivity[i] = rowUpper[i];
    }
  }
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    double value = newSolution[iColumn];
    //double thisTolerance = integerTolerance;
    if (fabs(floor(value + 0.5) - value) > integerTolerance) {
      double below = floor(value);
      double newValue = newSolution[iColumn];
      double cost = direction * objective[iColumn];
      double move;
      if (cost > 0.0) {
        // try up
        move = 1.0 - (value - below);
      } else if (cost < 0.0) {
        // try down
        move = below - value;
      } else {
        // won't be able to move unless we can grab another variable
        double randomNumber = randomNumberGenerator_.randomDouble();
        // which way?
        if (randomNumber < 0.5)
          move = below - value;
        else
          move = 1.0 - (value - below);
      }
      newValue += move;
      newSolution[iColumn] = newValue;
      newSolutionValue += move * cost;
      CoinBigIndex j;
      for (j = columnStart[iColumn];
           j < columnStart[iColumn] + columnLength[iColumn]; j++) {
        int iRow = row[j];
        rowActivity[iRow] += move * element[j];
      }
    }
  }

  // see if feasible - just using singletons
  for (i = 0; i < numberRows; i++) {
    double value = rowActivity[i];
    double thisInfeasibility = 0.0;
    if (value < rowLower[i] - primalTolerance)
      thisInfeasibility = value - rowLower[i];
    else if (value > rowUpper[i] + primalTolerance)
      thisInfeasibility = value - rowUpper[i];
    if (thisInfeasibility) {
      // See if there are any slacks I can use to fix up
      // maybe put in coding for multiple slacks?
      double bestCost = 1.0e50;
      CoinBigIndex k;
      int iBest = -1;
      double addCost = 0.0;
      double newValue = 0.0;
      double changeRowActivity = 0.0;
      double absInfeasibility = fabs(thisInfeasibility);
      for (k = rowStart[i]; k < rowStart[i] + rowLength[i]; k++) {
        int iColumn = column[k];
        // See if all elements help
        if (columnLength[iColumn] == 1) {
          double currentValue = newSolution[iColumn];
          double elementValue = elementByRow[k];
          double lowerValue = lower[iColumn];
          double upperValue = upper[iColumn];
          double gap = rowUpper[i] - rowLower[i];
          double absElement = fabs(elementValue);
          if (thisInfeasibility * elementValue > 0.0) {
            // we want to reduce
            if ((currentValue - lowerValue) * absElement >= absInfeasibility) {
              // possible - check if integer
              double distance = absInfeasibility / absElement;
              double thisCost = -direction * objective[iColumn] * distance;
              if (isHeuristicInteger(solver, iColumn)) {
                distance = ceil(distance - useTolerance);
                if (currentValue - distance >= lowerValue - useTolerance) {
                  if (absInfeasibility - distance * absElement < -gap - useTolerance)
                    thisCost = 1.0e100; // no good
                  else
                    thisCost = -direction * objective[iColumn] * distance;
                } else {
                  thisCost = 1.0e100; // no good
                }
              }
              if (thisCost < bestCost) {
                bestCost = thisCost;
                iBest = iColumn;
                addCost = thisCost;
                newValue = currentValue - distance;
                changeRowActivity = -distance * elementValue;
              }
            }
          } else {
            // we want to increase
            if ((upperValue - currentValue) * absElement >= absInfeasibility) {
              // possible - check if integer
              double distance = absInfeasibility / absElement;
              double thisCost = direction * objective[iColumn] * distance;
              if (isHeuristicInteger(solver, iColumn)) {
                distance = ceil(distance - 1.0e-7);
                assert(currentValue - distance <= upperValue + useTolerance);
                if (absInfeasibility - distance * absElement < -gap - useTolerance)
                  thisCost = 1.0e100; // no good
                else
                  thisCost = direction * objective[iColumn] * distance;
              }
              if (thisCost < bestCost) {
                bestCost = thisCost;
                iBest = iColumn;
                addCost = thisCost;
                newValue = currentValue + distance;
                changeRowActivity = distance * elementValue;
              }
            }
          }
        }
      }
      if (iBest >= 0) {
        /*printf("Infeasibility of %g on row %d cost %g\n",
                  thisInfeasibility,i,addCost);*/
        newSolution[iBest] = newValue;
        thisInfeasibility = 0.0;
        newSolutionValue += addCost;
        rowActivity[i] += changeRowActivity;
      }
    }
  }
  double penalty = 0.0;
  // integer variables may have wandered
  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    double value = newSolution[iColumn];
    if (fabs(floor(value + 0.5) - value) > integerTolerance)
      penalty += fabs(floor(value + 0.5) - value);
  }
  if (penalty) {
    // see if feasible using any
    // first continuous
    double penaltyChange = 0.0;
    int iColumn;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (isHeuristicInteger(solver, iColumn))
        continue;
      double currentValue = newSolution[iColumn];
      double lowerValue = lower[iColumn];
      double upperValue = upper[iColumn];
      CoinBigIndex j;
      int anyBadDown = 0;
      int anyBadUp = 0;
      double upImprovement = 0.0;
      double downImprovement = 0.0;
      for (j = columnStart[iColumn];
           j < columnStart[iColumn] + columnLength[iColumn]; j++) {
        int iRow = row[j];
        if (rowUpper[iRow] > rowLower[iRow]) {
          double value = element[j];
          if (rowActivity[iRow] > rowUpper[iRow] + primalTolerance) {
            // infeasible above
            downImprovement += value;
            upImprovement -= value;
            if (value > 0.0)
              anyBadUp++;
            else
              anyBadDown++;
          } else if (rowActivity[iRow] > rowUpper[iRow] - primalTolerance) {
            // feasible at ub
            if (value > 0.0) {
              upImprovement -= value;
              anyBadUp++;
            } else {
              downImprovement += value;
              anyBadDown++;
            }
          } else if (rowActivity[iRow] > rowLower[iRow] + primalTolerance) {
            // feasible in interior
          } else if (rowActivity[iRow] > rowLower[iRow] - primalTolerance) {
            // feasible at lb
            if (value < 0.0) {
              upImprovement += value;
              anyBadUp++;
            } else {
              downImprovement -= value;
              anyBadDown++;
            }
          } else {
            // infeasible below
            downImprovement -= value;
            upImprovement += value;
            if (value < 0.0)
              anyBadUp++;
            else
              anyBadDown++;
          }
        } else {
          // equality row
          double value = element[j];
          if (rowActivity[iRow] > rowUpper[iRow] + primalTolerance) {
            // infeasible above
            downImprovement += value;
            upImprovement -= value;
            if (value > 0.0)
              anyBadUp++;
            else
              anyBadDown++;
          } else if (rowActivity[iRow] < rowLower[iRow] - primalTolerance) {
            // infeasible below
            downImprovement -= value;
            upImprovement += value;
            if (value < 0.0)
              anyBadUp++;
            else
              anyBadDown++;
          } else {
            // feasible - no good
            anyBadUp = -1;
            anyBadDown = -1;
            break;
          }
        }
      }
      // could change tests for anyBad
      if (anyBadUp)
        upImprovement = 0.0;
      if (anyBadDown)
        downImprovement = 0.0;
      double way = 0.0;
      double improvement = 0.0;
      if (downImprovement > 0.0 && currentValue > lowerValue) {
        way = -1.0;
        improvement = downImprovement;
      } else if (upImprovement > 0.0 && currentValue < upperValue) {
        way = 1.0;
        improvement = upImprovement;
      }
      if (way) {
        // can improve
        double distance;
        if (way > 0.0)
          distance = upperValue - currentValue;
        else
          distance = currentValue - lowerValue;
        for (j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          double value = element[j] * way;
          if (rowActivity[iRow] > rowUpper[iRow] + primalTolerance) {
            // infeasible above
            assert(value < 0.0);
            double gap = rowActivity[iRow] - rowUpper[iRow];
            if (gap + value * distance < 0.0)
              distance = -gap / value;
          } else if (rowActivity[iRow] < rowLower[iRow] - primalTolerance) {
            // infeasible below
            assert(value > 0.0);
            double gap = rowActivity[iRow] - rowLower[iRow];
            if (gap + value * distance > 0.0)
              distance = -gap / value;
          } else {
            // feasible
            if (value > 0) {
              double gap = rowActivity[iRow] - rowUpper[iRow];
              if (gap + value * distance > 0.0)
                distance = -gap / value;
            } else {
              double gap = rowActivity[iRow] - rowLower[iRow];
              if (gap + value * distance < 0.0)
                distance = -gap / value;
            }
          }
        }
        //move
        penaltyChange += improvement * distance;
        distance *= way;
        newSolution[iColumn] += distance;
        newSolutionValue += direction * objective[iColumn] * distance;
        for (j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          double value = element[j];
          rowActivity[iRow] += distance * value;
        }
      }
    }
    // and now all if improving
    double lastChange = penaltyChange ? 1.0 : 0.0;
    int numberPasses = 0;
    while (lastChange > 1.0e-2 && numberPasses < 1000) {
      lastChange = 0;
      numberPasses++;
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        bool isInteger = isHeuristicInteger(solver, iColumn);
        double currentValue = newSolution[iColumn];
        double lowerValue = lower[iColumn];
        double upperValue = upper[iColumn];
        CoinBigIndex j;
        int anyBadDown = 0;
        int anyBadUp = 0;
        double upImprovement = 0.0;
        double downImprovement = 0.0;
        for (j = columnStart[iColumn];
             j < columnStart[iColumn] + columnLength[iColumn]; j++) {
          int iRow = row[j];
          double value = element[j];
          if (isInteger) {
            if (value > 0.0) {
              if (rowActivity[iRow] + value > rowUpper[iRow] + primalTolerance)
                anyBadUp++;
              if (rowActivity[iRow] - value < rowLower[iRow] - primalTolerance)
                anyBadDown++;
            } else {
              if (rowActivity[iRow] - value > rowUpper[iRow] + primalTolerance)
                anyBadDown++;
              if (rowActivity[iRow] + value < rowLower[iRow] - primalTolerance)
                anyBadUp++;
            }
          }
          if (rowUpper[iRow] > rowLower[iRow]) {
            if (rowActivity[iRow] > rowUpper[iRow] + primalTolerance) {
              // infeasible above
              downImprovement += value;
              upImprovement -= value;
              if (value > 0.0)
                anyBadUp++;
              else
                anyBadDown++;
            } else if (rowActivity[iRow] > rowUpper[iRow] - primalTolerance) {
              // feasible at ub
              if (value > 0.0) {
                upImprovement -= value;
                anyBadUp++;
              } else {
                downImprovement += value;
                anyBadDown++;
              }
            } else if (rowActivity[iRow] > rowLower[iRow] + primalTolerance) {
              // feasible in interior
            } else if (rowActivity[iRow] > rowLower[iRow] - primalTolerance) {
              // feasible at lb
              if (value < 0.0) {
                upImprovement += value;
                anyBadUp++;
              } else {
                downImprovement -= value;
                anyBadDown++;
              }
            } else {
              // infeasible below
              downImprovement -= value;
              upImprovement += value;
              if (value < 0.0)
                anyBadUp++;
              else
                anyBadDown++;
            }
          } else {
            // equality row
            if (rowActivity[iRow] > rowUpper[iRow] + primalTolerance) {
              // infeasible above
              downImprovement += value;
              upImprovement -= value;
              if (value > 0.0)
                anyBadUp++;
              else
                anyBadDown++;
            } else if (rowActivity[iRow] < rowLower[iRow] - primalTolerance) {
              // infeasible below
              downImprovement -= value;
              upImprovement += value;
              if (value < 0.0)
                anyBadUp++;
              else
                anyBadDown++;
            } else {
              // feasible - no good
              anyBadUp = -1;
              anyBadDown = -1;
              break;
            }
          }
        }
        // could change tests for anyBad
        if (anyBadUp)
          upImprovement = 0.0;
        if (anyBadDown)
          downImprovement = 0.0;
        double way = 0.0;
        double improvement = 0.0;
        if (downImprovement > 0.0 && currentValue > lowerValue) {
          way = -1.0;
          improvement = downImprovement;
          if (isInteger && currentValue < lowerValue + 0.99)
            continue; // no good
        } else if (upImprovement > 0.0 && currentValue < upperValue) {
          way = 1.0;
          improvement = upImprovement;
          if (isInteger && currentValue > upperValue - 0.99)
            continue; // no good
        }
        if (way) {
          // can improve
          double distance = COIN_DBL_MAX;
          for (j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            double value = element[j] * way;
            if (rowActivity[iRow] > rowUpper[iRow] + primalTolerance) {
              // infeasible above
              assert(value < 0.0);
              double gap = rowActivity[iRow] - rowUpper[iRow];
              if (gap + value * distance < 0.0) {
                // If integer then has to move by 1
                if (!isInteger)
                  distance = -gap / value;
                else
                  distance = CoinMax(-gap / value, 1.0);
              }
            } else if (rowActivity[iRow] < rowLower[iRow] - primalTolerance) {
              // infeasible below
              assert(value > 0.0);
              double gap = rowActivity[iRow] - rowLower[iRow];
              if (gap + value * distance > 0.0) {
                // If integer then has to move by 1
                if (!isInteger)
                  distance = -gap / value;
                else
                  distance = CoinMax(-gap / value, 1.0);
              }
            } else {
              // feasible
              if (value > 0) {
                double gap = rowActivity[iRow] - rowUpper[iRow];
                if (gap + value * distance > 0.0)
                  distance = -gap / value;
              } else {
                double gap = rowActivity[iRow] - rowLower[iRow];
                if (gap + value * distance < 0.0)
                  distance = -gap / value;
              }
            }
          }
          if (isInteger)
            distance = floor(distance + 1.05e-8);
          if (!distance) {
            // should never happen
            //printf("zero distance in CbcRounding - debug\n");
          }
          //move
          lastChange += improvement * distance;
          distance *= way;
          newSolution[iColumn] += distance;
          newSolutionValue += direction * objective[iColumn] * distance;
          for (j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            double value = element[j];
            rowActivity[iRow] += distance * value;
          }
        }
      }
      penaltyChange += lastChange;
    }
    penalty -= penaltyChange;
    if (penalty < 1.0e-5 * fabs(penaltyChange)) {
      // recompute
      penalty = 0.0;
      for (i = 0; i < numberRows; i++) {
        double value = rowActivity[i];
        if (value < rowLower[i] - primalTolerance)
          penalty += rowLower[i] - value;
        else if (value > rowUpper[i] + primalTolerance)
          penalty += value - rowUpper[i];
      }
    }
    // but integer variables may have wandered
    for (i = 0; i < numberIntegers; i++) {
      int iColumn = integerVariable[i];
      double value = newSolution[iColumn];
      if (fabs(floor(value + 0.5) - value) > integerTolerance)
	penalty += fabs(floor(value + 0.5) - value);
    }
  }

  // Could also set SOS (using random) and repeat
  if (!penalty) {
    // See if we can do better
    //seed_++;
    //CoinSeedRandom(seed_);
    // Random number between 0 and 1.
    double randomNumber = randomNumberGenerator_.randomDouble();
    int iPass;
    int start[2];
    int end[2];
    int iRandom = static_cast< int >(randomNumber * (static_cast< double >(numberIntegers)));
    start[0] = iRandom;
    end[0] = numberIntegers;
    start[1] = 0;
    end[1] = iRandom;
    for (iPass = 0; iPass < 2; iPass++) {
      int i;
      for (i = start[iPass]; i < end[iPass]; i++) {
        int iColumn = integerVariable[i];
#ifndef NDEBUG
        double value = newSolution[iColumn];
        assert(fabs(floor(value + 0.5) - value) <= integerTolerance);
#endif
        double cost = direction * objective[iColumn];
        double move = 0.0;
        if (cost > 0.0)
          move = -1.0;
        else if (cost < 0.0)
          move = 1.0;
        while (move) {
          bool good = true;
          double newValue = newSolution[iColumn] + move;
          if (newValue < lower[iColumn] - useTolerance || newValue > upper[iColumn] + useTolerance) {
            move = 0.0;
          } else {
            // see if we can move
            CoinBigIndex j;
            for (j = columnStart[iColumn];
                 j < columnStart[iColumn] + columnLength[iColumn]; j++) {
              int iRow = row[j];
              double newActivity = rowActivity[iRow] + move * element[j];
              if (newActivity < rowLower[iRow] - primalTolerance || newActivity > rowUpper[iRow] + primalTolerance) {
                good = false;
                break;
              }
            }
            if (good) {
              newSolution[iColumn] = newValue;
              newSolutionValue += move * cost;
              CoinBigIndex j;
              for (j = columnStart[iColumn];
                   j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                int iRow = row[j];
                rowActivity[iRow] += move * element[j];
              }
            } else {
              move = 0.0;
            }
          }
        }
      }
    }
    // Just in case of some stupidity
    double objOffset = 0.0;
    solver->getDblParam(OsiObjOffset, objOffset);
    newSolutionValue = -objOffset;
    for (i = 0; i < numberColumns; i++)
      newSolutionValue += objective[i] * newSolution[i];
    newSolutionValue *= direction;
    //printf("new solution value %g %g\n",newSolutionValue,solutionValue);
    if (newSolutionValue < solutionValue) {
      // paranoid check
      memset(rowActivity, 0, numberRows * sizeof(double));
      for (i = 0; i < numberColumns; i++) {
        CoinBigIndex j;
        double value = newSolution[i];
        if (value) {
          for (j = columnStart[i];
               j < columnStart[i] + columnLength[i]; j++) {
            int iRow = row[j];
            rowActivity[iRow] += value * element[j];
          }
        }
      }
      // check was approximately feasible
      bool feasible = true;
      for (i = 0; i < numberRows; i++) {
        if (rowActivity[i] < rowLower[i]) {
          if (rowActivity[i] < rowLower[i] - 1000.0 * primalTolerance)
            feasible = false;
        } else if (rowActivity[i] > rowUpper[i]) {
          if (rowActivity[i] > rowUpper[i] + 1000.0 * primalTolerance)
            feasible = false;
        }
      }
      if (feasible) {
        // new solution
        memcpy(betterSolution, newSolution, numberColumns * sizeof(double));
        solutionValue = newSolutionValue;
        //printf("** Solution of %g found by rounding\n",newSolutionValue);
        returnCode = 1;
      } else {
        // Can easily happen
        //printf("Debug CbcRounding giving bad solution\n");
      }
    }
  }
#ifdef NEW_ROUNDING
  if (!returnCode) {
#ifdef JJF_ZERO
    // back to starting point
    memcpy(newSolution, solution, numberColumns * sizeof(double));
    memset(rowActivity, 0, numberRows * sizeof(double));
    for (i = 0; i < numberColumns; i++) {
      int j;
      double value = newSolution[i];
      if (value < lower[i]) {
        value = lower[i];
        newSolution[i] = value;
      } else if (value > upper[i]) {
        value = upper[i];
        newSolution[i] = value;
      }
      if (value) {
        for (j = columnStart[i];
             j < columnStart[i] + columnLength[i]; j++) {
          int iRow = row[j];
          rowActivity[iRow] += value * element[j];
        }
      }
    }
    // check was feasible - if not adjust (cleaning may move)
    for (i = 0; i < numberRows; i++) {
      if (rowActivity[i] < rowLower[i]) {
        //assert (rowActivity[i]>rowLower[i]-1000.0*primalTolerance);
        rowActivity[i] = rowLower[i];
      } else if (rowActivity[i] > rowUpper[i]) {
        //assert (rowActivity[i]<rowUpper[i]+1000.0*primalTolerance);
        rowActivity[i] = rowUpper[i];
      }
    }
#endif
    int *candidate = new int[numberColumns];
    int nCandidate = 0;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      bool isInteger = isHeuristicInteger(solver, iColumn);
      if (isInteger) {
        double currentValue = newSolution[iColumn];
        if (fabs(currentValue - floor(currentValue + 0.5)) > 1.0e-8)
          candidate[nCandidate++] = iColumn;
      }
    }
    if (true) {
      // Rounding as in Berthold
      while (nCandidate) {
        double infeasibility = 1.0e-7;
        int iRow = -1;
        for (i = 0; i < numberRows; i++) {
          double value = 0.0;
          if (rowActivity[i] < rowLower[i]) {
            value = rowLower[i] - rowActivity[i];
          } else if (rowActivity[i] > rowUpper[i]) {
            value = rowActivity[i] - rowUpper[i];
          }
          if (value > infeasibility) {
            infeasibility = value;
            iRow = i;
          }
        }
        if (iRow >= 0) {
          // infeasible
        } else {
          // feasible
        }
      }
    } else {
      // Shifting as in Berthold
    }
    delete[] candidate;
  }
#endif
  delete[] newSolution;
  delete[] rowActivity;
  return returnCode;
}
// update model
void CbcRounding::setModel(CbcModel *model)
{
  model_ = model;
  // Get a copy of original matrix (and by row for rounding);
  assert(model_->solver());
  if (model_->solver()->getNumRows()) {
    matrix_ = *model_->solver()->getMatrixByCol();
    matrixByRow_ = *model_->solver()->getMatrixByRow();
    // make sure model okay for heuristic
    validate();
  }
}
// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
void CbcRounding::validate()
{
  if (model_ && (when() % 100) < 10) {
    if (model_->numberIntegers() != model_->numberObjects() && (model_->numberObjects() || (model_->specialOptions() & 1024) == 0)) {
      int numberOdd = 0;
      for (int i = 0; i < model_->numberObjects(); i++) {
        if (!model_->object(i)->canDoHeuristics())
          numberOdd++;
      }
      if (numberOdd)
        setWhen(0);
    }
  }
#ifdef NEW_ROUNDING
  int numberColumns = matrix_.getNumCols();
  down_ = new unsigned short[numberColumns];
  up_ = new unsigned short[numberColumns];
  equal_ = new unsigned short[numberColumns];
  // Column copy
  const double *element = matrix_.getElements();
  const int *row = matrix_.getIndices();
  const CoinBigIndex *columnStart = matrix_.getVectorStarts();
  const int *columnLength = matrix_.getVectorLengths();
  const double *rowLower = model.solver()->getRowLower();
  const double *rowUpper = model.solver()->getRowUpper();
  for (int i = 0; i < numberColumns; i++) {
    int down = 0;
    int up = 0;
    int equal = 0;
    if (columnLength[i] > 65535) {
      equal[0] = 65535;
      break; // unlikely to work
    }
    for (CoinBigIndex j = columnStart[i];
         j < columnStart[i] + columnLength[i]; j++) {
      int iRow = row[j];
      if (rowLower[iRow] > -1.0e20 && rowUpper[iRow] < 1.0e20) {
        equal++;
      } else if (element[j] > 0.0) {
        if (rowUpper[iRow] < 1.0e20)
          up++;
        else
          down--;
      } else {
        if (rowLower[iRow] > -1.0e20)
          up++;
        else
          down--;
      }
    }
    down_[i] = (unsigned short)down;
    up_[i] = (unsigned short)up;
    equal_[i] = (unsigned short)equal;
  }
#else
  down_ = NULL;
  up_ = NULL;
  equal_ = NULL;
#endif
}

// Default Constructor
CbcHeuristicPartial::CbcHeuristicPartial()
  : CbcHeuristic()
{
  fixPriority_ = 10000;
}

// Constructor from model
CbcHeuristicPartial::CbcHeuristicPartial(CbcModel &model, int fixPriority, int numberNodes)
  : CbcHeuristic(model)
{
  fixPriority_ = fixPriority;
  setNumberNodes(numberNodes);
  validate();
}

// Destructor
CbcHeuristicPartial::~CbcHeuristicPartial()
{
}

// Clone
CbcHeuristic *
CbcHeuristicPartial::clone() const
{
  return new CbcHeuristicPartial(*this);
}
// Create C++ lines to get to current state
void CbcHeuristicPartial::generateCpp(FILE *fp)
{
  CbcHeuristicPartial other;
  fprintf(fp, "0#include \"CbcHeuristic.hpp\"\n");
  fprintf(fp, "3  CbcHeuristicPartial partial(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp, "partial");
  if (fixPriority_ != other.fixPriority_)
    fprintf(fp, "3  partial.setFixPriority(%d);\n", fixPriority_);
  else
    fprintf(fp, "4  partial.setFixPriority(%d);\n", fixPriority_);
  fprintf(fp, "3  cbcModel->addHeuristic(&partial);\n");
}
//#define NEW_PARTIAL
// Copy constructor
CbcHeuristicPartial::CbcHeuristicPartial(const CbcHeuristicPartial &rhs)
  : CbcHeuristic(rhs)
  , fixPriority_(rhs.fixPriority_)
{
}

// Assignment operator
CbcHeuristicPartial &
CbcHeuristicPartial::operator=(const CbcHeuristicPartial &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    fixPriority_ = rhs.fixPriority_;
  }
  return *this;
}

// Resets stuff if model changes
void CbcHeuristicPartial::resetModel(CbcModel *model)
{
  model_ = model;
  // Get a copy of original matrix (and by row for partial);
  assert(model_->solver());
  validate();
}
// See if partial will give solution
// Sets value of solution
// Assumes rhs for original matrix still okay
// At present only works with integers
// Fix values if asked for
// Returns 1 if solution, 0 if not
int CbcHeuristicPartial::solution(double &solutionValue,
  double *betterSolution)
{
  // Return if already done
  if (fixPriority_ < 0)
    return 0; // switched off
#ifdef HEURISTIC_INFORM
  printf("Entering heuristic %s - nRuns %d numCould %d when %d\n",
    heuristicName(), numRuns_, numCouldRun_, when_);
#endif
  const double *hotstartSolution = model_->hotstartSolution();
  const int *hotstartPriorities = model_->hotstartPriorities();
  if (!hotstartSolution)
    return 0;
  OsiSolverInterface *solver = model_->solver();

  int numberIntegers = model_->numberIntegers();
  const int *integerVariable = model_->integerVariable();

  OsiSolverInterface *newSolver = model_->continuousSolver()->clone();
  const double *colLower = newSolver->getColLower();
  const double *colUpper = newSolver->getColUpper();

  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance, primalTolerance);

  int i;
  int numberFixed = 0;
  int returnCode = 0;

  for (i = 0; i < numberIntegers; i++) {
    int iColumn = integerVariable[i];
    if (abs(hotstartPriorities[iColumn]) <= fixPriority_) {
      double value = hotstartSolution[iColumn];
      double lower = colLower[iColumn];
      double upper = colUpper[iColumn];
      value = CoinMax(value, lower);
      value = CoinMin(value, upper);
      if (fabs(value - floor(value + 0.5)) < 1.0e-8) {
        value = floor(value + 0.5);
        newSolver->setColLower(iColumn, value);
        newSolver->setColUpper(iColumn, value);
        numberFixed++;
      }
    }
  }
  if (numberFixed > numberIntegers / 5 - 100000000) {
#ifdef COIN_DEVELOP
    printf("%d integers fixed\n", numberFixed);
#endif
    returnCode = smallBranchAndBound(newSolver, numberNodes_, betterSolution, solutionValue,
      model_->getCutoff(), "CbcHeuristicPartial");
    if (returnCode < 0)
      returnCode = 0; // returned on size
    //printf("return code %d",returnCode);
    if ((returnCode & 2) != 0) {
      // could add cut
      returnCode &= ~2;
      //printf("could add cut with %d elements (if all 0-1)\n",nFix);
    } else {
      //printf("\n");
    }
  }
  fixPriority_ = -1; // switch off

  delete newSolver;
  return returnCode;
}
// update model
void CbcHeuristicPartial::setModel(CbcModel *model)
{
  model_ = model;
  assert(model_->solver());
  // make sure model okay for heuristic
  validate();
}
// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
void CbcHeuristicPartial::validate()
{
  if (model_ && (when() % 100) < 10) {
    if (model_->numberIntegers() != model_->numberObjects())
      setWhen(0);
  }
}
bool CbcHeuristicPartial::shouldHeurRun(int /*whereFrom*/)
{
  return true;
}

// Default Constructor
CbcSerendipity::CbcSerendipity()
  : CbcHeuristic()
{
}

// Constructor from model
CbcSerendipity::CbcSerendipity(CbcModel &model)
  : CbcHeuristic(model)
{
}

// Destructor
CbcSerendipity::~CbcSerendipity()
{
}

// Clone
CbcHeuristic *
CbcSerendipity::clone() const
{
  return new CbcSerendipity(*this);
}
// Create C++ lines to get to current state
void CbcSerendipity::generateCpp(FILE *fp)
{
  fprintf(fp, "0#include \"CbcHeuristic.hpp\"\n");
  fprintf(fp, "3  CbcSerendipity serendipity(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp, "serendipity");
  fprintf(fp, "3  cbcModel->addHeuristic(&serendipity);\n");
}

// Copy constructor
CbcSerendipity::CbcSerendipity(const CbcSerendipity &rhs)
  : CbcHeuristic(rhs)
{
}

// Assignment operator
CbcSerendipity &
CbcSerendipity::operator=(const CbcSerendipity &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
  }
  return *this;
}

// Returns 1 if solution, 0 if not
int CbcSerendipity::solution(double &solutionValue,
  double *betterSolution)
{
  if (!model_)
    return 0;
#ifdef HEURISTIC_INFORM
  printf("Entering heuristic %s - nRuns %d numCould %d when %d\n",
    heuristicName(), numRuns_, numCouldRun_, when_);
#endif
  if (!inputSolution_) {
    // get information on solver type
    OsiAuxInfo *auxInfo = model_->solver()->getAuxiliaryInfo();
    OsiBabSolver *auxiliaryInfo = dynamic_cast< OsiBabSolver * >(auxInfo);
    if (auxiliaryInfo) {
      return auxiliaryInfo->solution(solutionValue, betterSolution, model_->solver()->getNumCols());
    } else {
      return 0;
    }
  } else {
    int numberColumns = model_->getNumCols();
    double value = inputSolution_[numberColumns];
    int returnCode = 0;
    if (value < solutionValue) {
      solutionValue = value;
      memcpy(betterSolution, inputSolution_, numberColumns * sizeof(double));
      returnCode = 1;
    }
    delete[] inputSolution_;
    inputSolution_ = NULL;
    model_ = NULL; // switch off
    return returnCode;
  }
}
// update model
void CbcSerendipity::setModel(CbcModel *model)
{
  model_ = model;
}
// Resets stuff if model changes
void CbcSerendipity::resetModel(CbcModel *model)
{
  model_ = model;
}

// Default Constructor
CbcHeuristicJustOne::CbcHeuristicJustOne()
  : CbcHeuristic()
  , probabilities_(NULL)
  , heuristic_(NULL)
  , numberHeuristics_(0)
{
}

// Constructor from model
CbcHeuristicJustOne::CbcHeuristicJustOne(CbcModel &model)
  : CbcHeuristic(model)
  , probabilities_(NULL)
  , heuristic_(NULL)
  , numberHeuristics_(0)
{
}

// Destructor
CbcHeuristicJustOne::~CbcHeuristicJustOne()
{
  for (int i = 0; i < numberHeuristics_; i++)
    delete heuristic_[i];
  delete[] heuristic_;
  delete[] probabilities_;
}

// Clone
CbcHeuristicJustOne *
CbcHeuristicJustOne::clone() const
{
  return new CbcHeuristicJustOne(*this);
}

// Create C++ lines to get to current state
void CbcHeuristicJustOne::generateCpp(FILE *fp)
{
  CbcHeuristicJustOne other;
  fprintf(fp, "0#include \"CbcHeuristicJustOne.hpp\"\n");
  fprintf(fp, "3  CbcHeuristicJustOne heuristicJustOne(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp, "heuristicJustOne");
  fprintf(fp, "3  cbcModel->addHeuristic(&heuristicJustOne);\n");
}

// Copy constructor
CbcHeuristicJustOne::CbcHeuristicJustOne(const CbcHeuristicJustOne &rhs)
  : CbcHeuristic(rhs)
  , probabilities_(NULL)
  , heuristic_(NULL)
  , numberHeuristics_(rhs.numberHeuristics_)
{
  if (numberHeuristics_) {
    probabilities_ = CoinCopyOfArray(rhs.probabilities_, numberHeuristics_);
    heuristic_ = new CbcHeuristic *[numberHeuristics_];
    for (int i = 0; i < numberHeuristics_; i++)
      heuristic_[i] = rhs.heuristic_[i]->clone();
  }
}

// Assignment operator
CbcHeuristicJustOne &
CbcHeuristicJustOne::operator=(const CbcHeuristicJustOne &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    for (int i = 0; i < numberHeuristics_; i++)
      delete heuristic_[i];
    delete[] heuristic_;
    delete[] probabilities_;
    probabilities_ = NULL;
    heuristic_ = NULL;
    numberHeuristics_ = rhs.numberHeuristics_;
    if (numberHeuristics_) {
      probabilities_ = CoinCopyOfArray(rhs.probabilities_, numberHeuristics_);
      heuristic_ = new CbcHeuristic *[numberHeuristics_];
      for (int i = 0; i < numberHeuristics_; i++)
        heuristic_[i] = rhs.heuristic_[i]->clone();
    }
  }
  return *this;
}
// Sets value of solution
// Returns 1 if solution, 0 if not
int CbcHeuristicJustOne::solution(double &solutionValue,
  double *betterSolution)
{
#ifdef DIVE_DEBUG
  std::cout << "solutionValue = " << solutionValue << std::endl;
#endif
  ++numCouldRun_;

  // test if the heuristic can run
  if (!shouldHeurRun_randomChoice() || !numberHeuristics_)
    return 0;
  double randomNumber = randomNumberGenerator_.randomDouble();
  int i;
  for (i = 0; i < numberHeuristics_; i++) {
    if (randomNumber < probabilities_[i])
      break;
  }
  assert(i < numberHeuristics_);
  int returnCode;
  //model_->unsetDivingHasRun();
#ifdef COIN_DEVELOP
  printf("JustOne running %s\n",
    heuristic_[i]->heuristicName());
#endif
  returnCode = heuristic_[i]->solution(solutionValue, betterSolution);
#ifdef COIN_DEVELOP
  if (returnCode)
    printf("JustOne running %s found solution\n",
      heuristic_[i]->heuristicName());
#endif
  return returnCode;
}
// Resets stuff if model changes
void CbcHeuristicJustOne::resetModel(CbcModel *model)
{
  CbcHeuristic::resetModel(model);
  for (int i = 0; i < numberHeuristics_; i++)
    heuristic_[i]->resetModel(model);
}
// update model (This is needed if cliques update matrix etc)
void CbcHeuristicJustOne::setModel(CbcModel *model)
{
  CbcHeuristic::setModel(model);
  for (int i = 0; i < numberHeuristics_; i++)
    heuristic_[i]->setModel(model);
}
// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
void CbcHeuristicJustOne::validate()
{
  CbcHeuristic::validate();
  for (int i = 0; i < numberHeuristics_; i++)
    heuristic_[i]->validate();
}
// Adds an heuristic with probability
void CbcHeuristicJustOne::addHeuristic(const CbcHeuristic *heuristic, double probability)
{
  CbcHeuristic *thisOne = heuristic->clone();
  thisOne->setWhen(-999);
  CbcHeuristic **tempH = CoinCopyOfArrayPartial(heuristic_, numberHeuristics_ + 1,
    numberHeuristics_);
  delete[] heuristic_;
  heuristic_ = tempH;
  heuristic_[numberHeuristics_] = thisOne;
  double *tempP = CoinCopyOfArrayPartial(probabilities_, numberHeuristics_ + 1,
    numberHeuristics_);
  delete[] probabilities_;
  probabilities_ = tempP;
  probabilities_[numberHeuristics_] = probability;
  numberHeuristics_++;
}
// Normalize probabilities
void CbcHeuristicJustOne::normalizeProbabilities()
{
  double sum = 0.0;
  for (int i = 0; i < numberHeuristics_; i++)
    sum += probabilities_[i];
  double multiplier = 1.0 / sum;
  sum = 0.0;
  for (int i = 0; i < numberHeuristics_; i++) {
    sum += probabilities_[i];
    probabilities_[i] = sum * multiplier;
  }
  assert(fabs(probabilities_[numberHeuristics_ - 1] - 1.0) < 1.0e-5);
  probabilities_[numberHeuristics_ - 1] = 1.000001;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
*/
