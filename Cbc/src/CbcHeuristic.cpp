// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
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
#include "CbcStrategy.hpp"
#include "CglPreProcess.hpp"
#include "OsiAuxInfo.hpp"
#include "OsiPresolve.hpp"
#include "CbcBranchActual.hpp"
#include "CbcCutGenerator.hpp"
//==============================================================================

CbcHeuristicNode::CbcHeuristicNode(const CbcHeuristicNode& rhs)
{
  numObjects_ = rhs.numObjects_;
  brObj_ = new CbcBranchingObject*[numObjects_];
  for (int i = 0; i < numObjects_; ++i) {
    brObj_[i] = rhs.brObj_[i]->clone();
  }
}

void
CbcHeuristicNodeList::gutsOfDelete()
{
  for (int i = nodes_.size() - 1; i >= 0; --i) {
    delete nodes_[i];
  }
}

void
CbcHeuristicNodeList::gutsOfCopy(const CbcHeuristicNodeList& rhs)
{
  append(rhs);
}

CbcHeuristicNodeList::CbcHeuristicNodeList(const CbcHeuristicNodeList& rhs)
{
  gutsOfCopy(rhs);
}

CbcHeuristicNodeList& CbcHeuristicNodeList::operator=
(const CbcHeuristicNodeList& rhs)
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

void
CbcHeuristicNodeList::append(CbcHeuristicNode*& node)
{
  nodes_.push_back(node);
  node = NULL;
}

void
CbcHeuristicNodeList::append(const CbcHeuristicNodeList& nodes)
{
  nodes_.reserve(nodes_.size() + nodes.size());
  for (int i = 0; i < nodes.size(); ++i) {
    CbcHeuristicNode* node = new CbcHeuristicNode(*nodes.node(i));
    append(node);
  }
}

//==============================================================================

// Default Constructor
CbcHeuristic::CbcHeuristic() :
  model_(NULL),
  when_(2),
  numberNodes_(200),
  feasibilityPumpOptions_(-1),
  fractionSmall_(1.0),
  heuristicName_("Unknown"),
  howOften_(1),
  decayFactor_(0.0),
  shallowDepth_(1),
  howOftenShallow_(1),
  numInvocationsInShallow_(0),
  numInvocationsInDeep_(0),
  lastRunDeep_(0),
  numRuns_(0),
  minDistanceToRun_(1),
  runNodes_(),
  numCouldRun_(0),
  numberSolutionsFound_(0),
  inputSolution_(NULL)
{
  // As CbcHeuristic virtual need to modify .cpp if above change
}

// Constructor from model
CbcHeuristic::CbcHeuristic(CbcModel & model) :
  model_(&model),
  when_(2),
  numberNodes_(200),
  feasibilityPumpOptions_(-1),
  fractionSmall_(1.0),
  heuristicName_("Unknown"),
  howOften_(1),
  decayFactor_(0.0),
  shallowDepth_(1),
  howOftenShallow_(1),
  numInvocationsInShallow_(0),
  numInvocationsInDeep_(0),
  lastRunDeep_(0),
  numRuns_(0),
  minDistanceToRun_(1),
  runNodes_(),
  numCouldRun_(0),
  numberSolutionsFound_(0),
  inputSolution_(NULL)
{}

void
CbcHeuristic::gutsOfCopy(const CbcHeuristic & rhs)
{
  model_ = rhs.model_;
  when_ = rhs.when_;
  numberNodes_ = rhs.numberNodes_;
  feasibilityPumpOptions_ = rhs.feasibilityPumpOptions_;
  fractionSmall_ = rhs.fractionSmall_;
  randomNumberGenerator_ = rhs.randomNumberGenerator_;
  heuristicName_ = rhs.heuristicName_;
  howOften_ = rhs.howOften_;
  decayFactor_ = rhs.howOften_;
  shallowDepth_= rhs.shallowDepth_;
  howOftenShallow_= rhs.howOftenShallow_;
  numInvocationsInShallow_ = rhs.numInvocationsInShallow_;
  numInvocationsInDeep_ = rhs.numInvocationsInDeep_;
  lastRunDeep_ = rhs.lastRunDeep_;
  numRuns_ = rhs.numRuns_;
  numCouldRun_ = rhs.numCouldRun_;
  minDistanceToRun_ = rhs.minDistanceToRun_;
  runNodes_ = rhs.runNodes_;
  numberSolutionsFound_ = rhs.numberSolutionsFound_;
  if (rhs.inputSolution_) {
    int numberColumns = model_->getNumCols();
    setInputSolution(rhs.inputSolution_,rhs.inputSolution_[numberColumns]);
  }
}
// Copy constructor 
CbcHeuristic::CbcHeuristic(const CbcHeuristic & rhs)
{
  inputSolution_ = NULL;
  gutsOfCopy(rhs);
}

// Assignment operator 
CbcHeuristic & 
CbcHeuristic::operator=( const CbcHeuristic& rhs)
{
  if (this!=&rhs) {
    gutsOfDelete();
    gutsOfCopy(rhs);
  }
  return *this;
}

void CbcHeurDebugNodes(CbcModel* model_)
{
  CbcNode* node = model_->currentNode();
  CbcNodeInfo* nodeInfo = node->nodeInfo();
  std::cout << "===============================================================\n";
  while (nodeInfo) {
    const CbcNode* node = nodeInfo->owner();
    printf("nodeinfo: node %i\n", nodeInfo->nodeNumber());
    {
      const CbcIntegerBranchingObject* brPrint =
	dynamic_cast<const CbcIntegerBranchingObject*>(nodeInfo->parentBranch());
      if (!brPrint) {
	printf("    parentBranch: NULL\n");
      } else {
	const double* downBounds = brPrint->downBounds();
	const double* upBounds = brPrint->upBounds();
	int variable = brPrint->variable();
	int way = brPrint->way();
	printf("   parentBranch: var %i downBd [%i,%i] upBd [%i,%i] way %i\n",
	       variable, (int)downBounds[0], (int)downBounds[1],
	       (int)upBounds[0], (int)upBounds[1], way);
      }
    }
    if (! node) {
      printf("    owner: NULL\n");
    } else {
      printf("    owner: node %i depth %i onTree %i active %i",
	     node->nodeNumber(), node->depth(), node->onTree(), node->active());
      const OsiBranchingObject* osibr =
	nodeInfo->owner()->branchingObject();
      const CbcBranchingObject* cbcbr =
	dynamic_cast<const CbcBranchingObject*>(osibr);
      const CbcIntegerBranchingObject* brPrint =
	dynamic_cast<const CbcIntegerBranchingObject*>(cbcbr);
      if (!brPrint) {
	printf("        ownerBranch: NULL\n");
      } else {
	const double* downBounds = brPrint->downBounds();
	const double* upBounds = brPrint->upBounds();
	int variable = brPrint->variable();
	int way = brPrint->way();
	printf("        ownerbranch: var %i downBd [%i,%i] upBd [%i,%i] way %i\n",
	       variable, (int)downBounds[0], (int)downBounds[1],
	       (int)upBounds[0], (int)upBounds[1], way);
      }
    }
    nodeInfo = nodeInfo->parent();
  }
}

void
CbcHeuristic::debugNodes()
{
  CbcHeurDebugNodes(model_);
}

void
CbcHeuristic::printDistanceToNodes()
{
  const CbcNode* currentNode = model_->currentNode();
  if (currentNode != NULL) {
    CbcHeuristicNode* nodeDesc = new CbcHeuristicNode(*model_);
    for (int i = runNodes_.size() - 1; i >= 0; --i) {
      nodeDesc->distance(runNodes_.node(i));
    }
    runNodes_.append(nodeDesc);
  }
}

bool
CbcHeuristic::shouldHeurRun()
{

#if 0
  const CbcNode* currentNode = model_->currentNode();
  if (currentNode == NULL) {
    return false;
  }

  debugNodes();
//   return false;

  const int depth = currentNode->depth();
#else
  int depth = 0;
  const CbcNode* currentNode = model_->currentNode();
  if (currentNode != NULL) {
    depth = currentNode->depth();
#ifdef PRINT_DEBUG
    debugNodes();
#endif
  }
#endif

  const int nodeCount = model_->getNodeCount();  // FIXME: check that this is
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
#if 1
    if (currentNode != NULL) {
      // Get where we are and create the appropriate CbcHeuristicNode object
      CbcHeuristicNode* nodeDesc = new CbcHeuristicNode(*model_);
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
    if (model_->getCurrentPassNumber() != 1) {
      // Run the heuristic only when first entering the node.
      // LL: I don't think this is right. It should run just before strong
      // LL: branching, I believe.
      return false;
    }
    // Get where we are and create the appropriate CbcHeuristicNode object
    CbcHeuristicNode* nodeDesc = new CbcHeuristicNode(*model_);
    //#ifdef PRINT_DEBUG
#if 1
    const double minDistanceToRun = 1.5 * log((double)depth) / log((double)2);
#else
    const double minDistanceToRun = minDistanceToRun_;
#endif
#ifdef PRINT_DEBUG
    double minDistance = nodeDesc->minDistance(runNodes_);
    std::cout<<"minDistance = "<<minDistance
	     <<", minDistanceToRun = "<<minDistanceToRun<<std::endl;
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
}

bool
CbcHeuristic::shouldHeurRun_randomChoice()
{
  int depth = 0;
  const CbcNode* currentNode = model_->currentNode();
  if (currentNode != NULL) {
    depth = currentNode->depth();
  }

  if(depth != 0) {
    const double numerator = depth * depth;
    const double denominator = exp(depth * log((double)2));
    double probability = numerator / denominator;
    double randomNumber = randomNumberGenerator_.randomDouble();
    if (when_>2&&when_<7) {
      /* JJF adjustments
	 3 only at root and if no solution
	 4 only at root and if this heuristic has not got solution
	 5 only at depth <4
	 6 decay
      */
      switch(when_) {
      case 3:
	if (model_->bestSolution())
	  probability=-1.0;
	break;
      case 4:
	if (numberSolutionsFound_)
	  probability=-1.0;
	break;
      case 5:
	if (depth>=4)
	  probability=-1.0;
	break;
      case 6:
	if (depth>=4) {
	  if (numberSolutionsFound_*howOften_<numCouldRun_)
	    howOften_ = (int) (howOften_*1.1);
	  probability = 1.0/howOften_;
	  if (model_->bestSolution())
	    probability *= 0.5;
	}
	break;
      }
    }
    if (randomNumber>probability)
      return false;
    
    if (model_->getCurrentPassNumber() > 1)
      return false;
  }
  ++numRuns_;
  return true;
}

// Resets stuff if model changes
void 
CbcHeuristic::resetModel(CbcModel * model)
{
  model_=model;
  }
// Set seed
void
CbcHeuristic::setSeed(int value)
{
  randomNumberGenerator_.setSeed(value);
}

// Create C++ lines to get to current state
void 
CbcHeuristic::generateCpp( FILE * fp, const char * heuristic) 
{
  // hard coded as CbcHeuristic virtual
  if (when_!=2)
    fprintf(fp,"3  %s.setWhen(%d);\n",heuristic,when_);
  else
    fprintf(fp,"4  %s.setWhen(%d);\n",heuristic,when_);
  if (numberNodes_!=200)
    fprintf(fp,"3  %s.setNumberNodes(%d);\n",heuristic,numberNodes_);
  else
    fprintf(fp,"4  %s.setNumberNodes(%d);\n",heuristic,numberNodes_);
  if (fractionSmall_!=1.0)
    fprintf(fp,"3  %s.setFractionSmall(%g);\n",heuristic,fractionSmall_);
  else
    fprintf(fp,"4  %s.setFractionSmall(%g);\n",heuristic,fractionSmall_);
  if (heuristicName_ != "Unknown")
    fprintf(fp,"3  %s.setHeuristicName(\"%s\");\n",
	    heuristic,heuristicName_.c_str()) ;
  else
    fprintf(fp,"4  %s.setHeuristicName(\"%s\");\n",
	    heuristic,heuristicName_.c_str()) ;
}
// Destructor 
CbcHeuristic::~CbcHeuristic ()
{
}

// update model
void CbcHeuristic::setModel(CbcModel * model)
{
  model_ = model;
}
#ifdef COIN_DEVELOP
extern bool getHistoryStatistics_;
#endif
// Do mini branch and bound (return 1 if solution)
int 
CbcHeuristic::smallBranchAndBound(OsiSolverInterface * solver,int numberNodes,
                                  double * newSolution, double & newSolutionValue,
                                  double cutoff, std::string name) const
{
  // size before
  double before = 2*solver->getNumRows()+solver->getNumCols();
  // Use this fraction
  double fractionSmall = fractionSmall_;
  if (before>80000.0) {
    // fairly large - be more conservative
    double multiplier = (0.7*200000.0)/CoinMin(before,200000.0);
    if (multiplier<1.0)
      fractionSmall *= multiplier;
#ifdef COIN_DEVELOP
    printf("changing fractionSmall from %g to %g for %s\n",
	   fractionSmall_,fractionSmall,name.c_str());
#endif
  }
#ifdef COIN_HAS_CLP
  OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (solver);
  if (osiclp&&(osiclp->specialOptions()&65536)==0) {
    // go faster stripes
    if (osiclp->getNumRows()<300&&osiclp->getNumCols()<500) {
      osiclp->setupForRepeatedUse(2,0);
    } else {
      osiclp->setupForRepeatedUse(0,0);
    }
    // Turn this off if you get problems
    // Used to be automatically set
    osiclp->setSpecialOptions(osiclp->specialOptions()|(128+64));
    ClpSimplex * lpSolver = osiclp->getModelPtr();
    lpSolver->setSpecialOptions(lpSolver->specialOptions()|0x01000000); // say is Cbc (and in branch and bound)
  }
#endif
#ifdef COIN_DEVELOP
  getHistoryStatistics_=false;
#endif
  int status=0;
  int logLevel = model_->logLevel();
#define LEN_PRINT 250
  char generalPrint[LEN_PRINT];
  // Do presolve to see if possible
  int numberColumns = solver->getNumCols();
  char * reset = NULL;
  int returnCode=1;
  {
    int saveLogLevel = solver->messageHandler()->logLevel();
    if (saveLogLevel==1)
      solver->messageHandler()->setLogLevel(0);
    OsiPresolve * pinfo = new OsiPresolve();
    int presolveActions=0;
    // Allow dual stuff on integers
    presolveActions=1;
    // Do not allow all +1 to be tampered with
    //if (allPlusOnes)
    //presolveActions |= 2;
    // allow transfer of costs
    // presolveActions |= 4;
    pinfo->setPresolveActions(presolveActions);
    OsiSolverInterface * presolvedModel = pinfo->presolvedModel(*solver,1.0e-8,true,2);
    delete pinfo;
    // see if too big
    if (presolvedModel) {
      int afterRows = presolvedModel->getNumRows();
      int afterCols = presolvedModel->getNumCols();
      delete presolvedModel;
      double after = 2*afterRows+afterCols;
      if (after>fractionSmall*before&&after>300&&numberNodes>=0) {
	// Need code to try again to compress further using used
	const int * used =  model_->usedInSolution();
	int maxUsed=0;
	int iColumn;
	const double * lower = solver->getColLower();
	const double * upper = solver->getColUpper();
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  if (upper[iColumn]>lower[iColumn]) {
	    if (solver->isBinary(iColumn))
	      maxUsed = CoinMax(maxUsed,used[iColumn]);
	  }
	}
	if (maxUsed) {
	  reset = new char [numberColumns];
	  int nFix=0;
	  for (iColumn=0;iColumn<numberColumns;iColumn++) {
	    reset[iColumn]=0;
	    if (upper[iColumn]>lower[iColumn]) {
	      if (solver->isBinary(iColumn)&&used[iColumn]==maxUsed) {
		bool setValue=true;
		if (maxUsed==1) {
		  double randomNumber = randomNumberGenerator_.randomDouble();
		  if (randomNumber>0.3)
		    setValue=false;
		}
		if (setValue) {
		  reset[iColumn]=1;
		  solver->setColLower(iColumn,1.0);
		  nFix++;
		}
	      }
	    }
	  }
	  pinfo = new OsiPresolve();
	  presolveActions=0;
	  // Allow dual stuff on integers
	  presolveActions=1;
	  // Do not allow all +1 to be tampered with
	  //if (allPlusOnes)
	  //presolveActions |= 2;
	  // allow transfer of costs
	  // presolveActions |= 4;
	  pinfo->setPresolveActions(presolveActions);
	  presolvedModel = pinfo->presolvedModel(*solver,1.0e-8,true,2);
	  delete pinfo;
	  if(presolvedModel) {
	    // see if too big
	    int afterRows2 = presolvedModel->getNumRows();
	    int afterCols2 = presolvedModel->getNumCols();
	    delete presolvedModel;
	    double after = 2*afterRows2+afterCols2;
	    if (after>fractionSmall*before&&(after>300||numberNodes<0)) {
	      sprintf(generalPrint,"Full problem %d rows %d columns, reduced to %d rows %d columns - %d fixed gives %d, %d - still too large",
		      solver->getNumRows(),solver->getNumCols(),
		      afterRows,afterCols,nFix,afterRows2,afterCols2);
	    } else {
	      sprintf(generalPrint,"Full problem %d rows %d columns, reduced to %d rows %d columns - %d fixed gives %d, %d - ok now",
		      solver->getNumRows(),solver->getNumCols(),
		      afterRows,afterCols,nFix,afterRows2,afterCols2);
	    }
	    model_->messageHandler()->message(CBC_FPUMP1,model_->messages())
	      << generalPrint
	      <<CoinMessageEol;
	  } else {
	    returnCode=-1; // infeasible
	  }
	}
      }
    } else {
      returnCode=-1; // infeasible
    }
    solver->messageHandler()->setLogLevel(saveLogLevel);
  }
  if (returnCode==-1) {
    delete [] reset;
#ifdef COIN_DEVELOP
    getHistoryStatistics_=true;
#endif
    return returnCode;
  }
  // Reduce printout
  bool takeHint;
  OsiHintStrength strength;
  solver->getHintParam(OsiDoReducePrint,takeHint,strength);
  solver->setHintParam(OsiDoReducePrint,true,OsiHintTry);
  solver->setHintParam(OsiDoPresolveInInitial,false,OsiHintTry);
  solver->setDblParam(OsiDualObjectiveLimit,cutoff*solver->getObjSense());
  solver->initialSolve();
  if (solver->isProvenOptimal()) {
    CglPreProcess process;
    /* Do not try and produce equality cliques and
       do up to 2 passes */
    if (logLevel<=1)
      process.messageHandler()->setLogLevel(0);
    OsiSolverInterface * solver2= process.preProcessNonDefault(*solver,false,2);
    if (!solver2) {
      if (logLevel>1)
        printf("Pre-processing says infeasible\n");
      returnCode=2; // so will be infeasible
    } else {
      // see if too big
      //double before = 2*solver->getNumRows()+solver->getNumCols();
      double after = 2*solver2->getNumRows()+solver2->getNumCols();
      if (after>fractionSmall*before&&(after>300||numberNodes<0)) {
	sprintf(generalPrint,"Full problem %d rows %d columns, reduced to %d rows %d columns - too large",
		solver->getNumRows(),solver->getNumCols(),
		solver2->getNumRows(),solver2->getNumCols());
	model_->messageHandler()->message(CBC_FPUMP1,model_->messages())
	  << generalPrint
	  <<CoinMessageEol;
	returnCode = -1;
      } else {
	sprintf(generalPrint,"Full problem %d rows %d columns, reduced to %d rows %d columns",
		solver->getNumRows(),solver->getNumCols(),
		solver2->getNumRows(),solver2->getNumCols());
	model_->messageHandler()->message(CBC_FPUMP1,model_->messages())
	  << generalPrint
	  <<CoinMessageEol;
      }
      if (returnCode==1) {
	solver2->resolve();
	CbcModel model(*solver2);
	if (numberNodes>=0) {
	  // normal
	  if (logLevel<=1)
	    model.setLogLevel(0);
	  else
	    model.setLogLevel(logLevel);
	  // No small fathoming
	  model.setFastNodeDepth(-1);
	  model.setCutoff(cutoff);
	  model.setMaximumNodes(numberNodes);
	  model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
	  // Lightweight
	  CbcStrategyDefaultSubTree strategy(model_,true,5,1,0);
	  model.setStrategy(strategy);
	  model.solver()->setIntParam(OsiMaxNumIterationHotStart,10);
	  model.setMaximumCutPassesAtRoot(CoinMin(20,model_->getMaximumCutPassesAtRoot()));
	} else {
	  model_->messageHandler()->message(CBC_RESTART,model_->messages())
	    <<solver2->getNumRows()<<solver2->getNumCols()
	    <<CoinMessageEol;
	  // going for full search and copy across more stuff
	  model.gutsOfCopy(*model_,2);
	  for (int i=0;i<model.numberCutGenerators();i++)
	    model.cutGenerator(i)->setTiming(true);
	  model.setCutoff(cutoff);
	  bool takeHint;
	  OsiHintStrength strength;
	  // Switch off printing if asked to
	  model_->solver()->getHintParam(OsiDoReducePrint,takeHint,strength);
	  model.solver()->setHintParam(OsiDoReducePrint,takeHint,strength);
	  CbcStrategyDefault strategy(true,model_->numberStrong(),
				      model_->numberBeforeTrust());
	  // Set up pre-processing - no 
	  strategy.setupPreProcessing(0); // was (4);
	  model.setStrategy(strategy);
	}
	if (inputSolution_) {
	  // translate and add a serendipity heuristic
	  int numberColumns = solver2->getNumCols();
	  const int * which = process.originalColumns();
	  OsiSolverInterface * solver3 = solver2->clone();
	  for (int i=0;i<numberColumns;i++) {
	    if (solver3->isInteger(i)) {
	      int k=which[i];
	      double value = inputSolution_[k];
	      solver3->setColLower(i,value);
	      solver3->setColUpper(i,value);
	    }
	  }
	  solver3->setDblParam(OsiDualObjectiveLimit,COIN_DBL_MAX);
	  solver3->resolve();
	  if (solver3->isProvenOptimal()) {
	    // good
	    CbcSerendipity heuristic(model);
	    double value = solver3->getObjSense()*solver3->getObjValue();
	    heuristic.setInputSolution(solver3->getColSolution(),value);
	    model.setCutoff(value+1.0e-7*(1.0+fabs(value)));
	    model.addHeuristic(&heuristic,"previous solution");
	  }
	  delete solver3;
	}
	// Do search
	if (logLevel>1)
	  model_->messageHandler()->message(CBC_START_SUB,model_->messages())
	    << name
	    << model.getMaximumNodes()
	    <<CoinMessageEol;
	// probably faster to use a basis to get integer solutions
	model.setSpecialOptions(2);
#ifdef CBC_THREAD
	if (model_->getNumberThreads()>0&&(model_->getThreadMode()&1)!=0) {
	  // See if at root node
	  bool atRoot = model_->getNodeCount()==0;
	  int passNumber = model_->getCurrentPassNumber();
	  if (atRoot&&passNumber==1)
	    model.setNumberThreads(model_->getNumberThreads());
	}
#endif
	model.setParentModel(*model_);
	model.setOriginalColumns(process.originalColumns());
	model.setSearchStrategy(-1);
	// If no feasibility pump then insert a lightweight one
	if (feasibilityPumpOptions_>=0) {
	  bool gotPump=false;
	  for (int i=0;i<model.numberHeuristics();i++) {
	    const CbcHeuristicFPump* pump =
	      dynamic_cast<const CbcHeuristicFPump*>(model_->heuristic(i));
	    if (pump) 
	      gotPump=true;
	  }
	  if (!gotPump) {
	    CbcHeuristicFPump heuristic4;
	    int pumpTune=feasibilityPumpOptions_;
	    if (pumpTune>0) {
	      /*
		>=10000000 for using obj
		>=1000000 use as accumulate switch
		>=1000 use index+1 as number of large loops
		>=100 use 0.05 objvalue as increment
		>=10 use +0.1 objvalue for cutoff (add)
		1 == fix ints at bounds, 2 fix all integral ints, 3 and continuous at bounds
		4 and static continuous, 5 as 3 but no internal integers
		6 as 3 but all slack basis!
	      */
	      double value = solver2->getObjSense()*solver2->getObjValue();
	      int w = pumpTune/10;
	      int c = w % 10;
	      w /= 10;
	      int i = w % 10;
	      w /= 10;
	      int r = w;
	      int accumulate = r/1000;
	      r -= 1000*accumulate;
	      if (accumulate>=10) {
		int which = accumulate/10;
		accumulate -= 10*which;
		which--;
		// weights and factors
		double weight[]={0.1,0.1,0.5,0.5,1.0,1.0,5.0,5.0};
		double factor[] = {0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5};
		heuristic4.setInitialWeight(weight[which]);
		heuristic4.setWeightFactor(factor[which]);
	      }
	      // fake cutoff
	      if (c) {
		double cutoff;
		solver2->getDblParam(OsiDualObjectiveLimit,cutoff);
		cutoff = CoinMin(cutoff,value + 0.1*fabs(value)*c);
		heuristic4.setFakeCutoff(cutoff);
	      }
	      if (i||r) {
		// also set increment
		//double increment = (0.01*i+0.005)*(fabs(value)+1.0e-12);
		double increment = 0.0;
		heuristic4.setAbsoluteIncrement(increment);
		heuristic4.setAccumulate(accumulate);
		heuristic4.setMaximumRetries(r+1);
	      }
	      pumpTune = pumpTune%100;
	      if (pumpTune==6)
		pumpTune =13;
	      heuristic4.setWhen(pumpTune+10);
	    }
	    model.addHeuristic(&heuristic4,"feasibility pump",0);
	  }
	}
	if (model_->searchStrategy()==2) {
	  model.setNumberStrong(5);
	  model.setNumberBeforeTrust(5);
	}
	if (model.getNumCols()) {
	  if (numberNodes>=0)
	    setCutAndHeuristicOptions(model);
	  model.branchAndBound();
	  if (numberNodes<0) {
	    model_->incrementIterationCount(model.getIterationCount());
	    model_->incrementNodeCount(model.getNodeCount());
	    for (int iGenerator=0;iGenerator<model.numberCutGenerators();iGenerator++) {
	      CbcCutGenerator * generator = model.cutGenerator(iGenerator);
	      printf("%s was tried %d times and created %d cuts of which %d were active after adding rounds of cuts",
		      generator->cutGeneratorName(),
		      generator->numberTimesEntered(),
		      generator->numberCutsInTotal()+
		      generator->numberColumnCuts(),
		      generator->numberCutsActive());
	      printf(" (%.3f seconds)\n)",generator->timeInCutGenerator());
	    }
	  }
	} else {
	  // empty model
	  model.setMinimizationObjValue(model.solver()->getObjSense()*model.solver()->getObjValue());
	}
	if (logLevel>1)
	  model_->messageHandler()->message(CBC_END_SUB,model_->messages())
	    << name
	    <<CoinMessageEol;
	if (model.getMinimizationObjValue()<CoinMin(cutoff,1.0e30)) {
	  // solution
	  if (model.getNumCols())
	    returnCode=model.isProvenOptimal() ? 3 : 1;
	  else
	    returnCode=3;
	  // post process
#ifdef COIN_HAS_CLP
	  OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (model.solver());
	  if (clpSolver) {
	    ClpSimplex * lpSolver = clpSolver->getModelPtr();
	    lpSolver->setSpecialOptions(lpSolver->specialOptions()|0x01000000); // say is Cbc (and in branch and bound)
	  }
#endif
	  process.postProcess(*model.solver());
	  if (solver->isProvenOptimal()) {
	    // Solution now back in solver
	    int numberColumns = solver->getNumCols();
	    memcpy(newSolution,solver->getColSolution(),
		   numberColumns*sizeof(double));
	    newSolutionValue = model.getMinimizationObjValue();
	  } else {
	    // odd - but no good
	    returnCode=0; // so will be infeasible
	  }
	} else {
        // no good
	  returnCode=model.isProvenInfeasible() ? 2 : 0; // so will be infeasible
	}
	if (model.status()==5)
	  returnCode=-2; // stop
	if (model.isProvenInfeasible())
	  status=1;
	else if (model.isProvenOptimal())
	  status=2;
      }
    }
  } else {
    returnCode=2; // infeasible finished
  }
  model_->setLogLevel(logLevel);
  if (reset) {
    for (int iColumn=0;iColumn<numberColumns;iColumn++) {
      if (reset[iColumn])
	solver->setColLower(iColumn,0.0);
    }
    delete [] reset;
  }
#ifdef COIN_DEVELOP
  getHistoryStatistics_=true;
  if (returnCode==1||returnCode==2) {
    if (status==1)
      printf("heuristic could add cut because infeasible (%s)\n",heuristicName_.c_str()); 
    else if (status==2)
      printf("heuristic could add cut because optimal (%s)\n",heuristicName_.c_str());
  } 
#endif
  solver->setHintParam(OsiDoReducePrint,takeHint,strength);
  return returnCode;
}
// Set input solution
void 
CbcHeuristic::setInputSolution(const double * solution, double objValue)
{
  delete [] inputSolution_;
  inputSolution_=NULL;
  if (model_&&solution) {
    int numberColumns = model_->getNumCols();
    inputSolution_ = new double [numberColumns+1];
    memcpy(inputSolution_,solution,numberColumns*sizeof(double));
    inputSolution_[numberColumns]=objValue;
  }
}

//##############################################################################

inline int compare3BranchingObjects(const CbcBranchingObject* br0,
				    const CbcBranchingObject* br1)
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

inline bool compareBranchingObjects(const CbcBranchingObject* br0,
				    const CbcBranchingObject* br1)
{
  return compare3BranchingObjects(br0, br1) < 0;
}

//==============================================================================

void
CbcHeuristicNode::gutsOfConstructor(CbcModel& model)
{
  //  CbcHeurDebugNodes(&model);
  CbcNode* node = model.currentNode();
  brObj_ = new CbcBranchingObject*[node->depth()];
  CbcNodeInfo* nodeInfo = node->nodeInfo();
  int cnt = 0;
  while (nodeInfo->parentBranch() != NULL) {
    const OsiBranchingObject* br = nodeInfo->parentBranch();
    const CbcBranchingObject* cbcbr = dynamic_cast<const CbcBranchingObject*>(br);
    if (! cbcbr) {
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
  std::sort(brObj_, brObj_+cnt, compareBranchingObjects);
  if (cnt <= 1) {
    numObjects_ = cnt;
  } else {
    numObjects_ = 0;
    CbcBranchingObject* br=NULL; // What should this be?
    for (int i = 1; i < cnt; ++i) {
      if (compare3BranchingObjects(brObj_[numObjects_], brObj_[i]) == 0) {
	int comp = brObj_[numObjects_]->compareBranchingObject(brObj_[i], &br);
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

CbcHeuristicNode::CbcHeuristicNode(CbcModel& model)
{
  gutsOfConstructor(model);
}

//==============================================================================

double
CbcHeuristicNode::distance(const CbcHeuristicNode* node) const 
{

  const double disjointWeight = 1;
  const double overlapWeight = 0.4;
  const double subsetWeight = 0.2;
  int countDisjointWeight = 0;
  int countOverlapWeight = 0;
  int countSubsetWeight = 0;
  int i = 0; 
  int j = 0;
  double dist = 0.0;
#ifdef PRINT_DEBUG
  printf(" numObjects_ = %i, node->numObjects_ = %i\n",
	 numObjects_, node->numObjects_);
#endif
  while( i < numObjects_ && j < node->numObjects_) {
    CbcBranchingObject* br0 = brObj_[i];
    const CbcBranchingObject* br1 = node->brObj_[j];
#ifdef PRINT_DEBUG
    const CbcIntegerBranchingObject* brPrint0 =
      dynamic_cast<const CbcIntegerBranchingObject*>(br0);
    const double* downBounds = brPrint0->downBounds();
    const double* upBounds = brPrint0->upBounds();
    int variable = brPrint0->variable();
    int way = brPrint0->way();
    printf("   br0: var %i downBd [%i,%i] upBd [%i,%i] way %i\n",
	   variable, (int)downBounds[0], (int)downBounds[1],
	   (int)upBounds[0], (int)upBounds[1], way);
    const CbcIntegerBranchingObject* brPrint1 =
      dynamic_cast<const CbcIntegerBranchingObject*>(br1);
    downBounds = brPrint1->downBounds();
    upBounds = brPrint1->upBounds();
    variable = brPrint1->variable();
    way = brPrint1->way();
    printf("   br1: var %i downBd [%i,%i] upBd [%i,%i] way %i\n",
	   variable, (int)downBounds[0], (int)downBounds[1],
	   (int)upBounds[0], (int)upBounds[1], way);
#endif
    const int brComp = compare3BranchingObjects(br0, br1);
    if (brComp < 0) {
      dist += subsetWeight;
      countSubsetWeight++;
      ++i;
    }
    else if (brComp > 0) {
      dist += subsetWeight;
      countSubsetWeight++;
      ++j;
    }
    else {
      const int comp = br0->compareBranchingObject(br1, false);
      switch (comp) {
      case CbcRangeSame:
	// do nothing
	break;
      case CbcRangeDisjoint: // disjoint decisions
	dist += disjointWeight;
	countDisjointWeight++;
	break;
      case CbcRangeSubset: // subset one way or another
      case CbcRangeSuperset:
	dist += subsetWeight;
	countSubsetWeight++;
	break;
      case CbcRangeOverlap: // overlap
	dist += overlapWeight;
	countOverlapWeight++;
	break;
      }
      ++i;
      ++j;
    }
  }
  dist += subsetWeight * (numObjects_ - i + node->numObjects_ - j);
  countSubsetWeight += (numObjects_ - i + node->numObjects_ - j);
  printf("subset = %i, overlap = %i, disjoint = %i\n", countSubsetWeight,
	 countOverlapWeight, countDisjointWeight);
  return dist;
}

//==============================================================================

CbcHeuristicNode::~CbcHeuristicNode()
{
  for (int i = 0; i < numObjects_; ++i) {
    delete brObj_[i];
  }
  delete [] brObj_;
}

//==============================================================================

double
CbcHeuristicNode::minDistance(const CbcHeuristicNodeList& nodeList) const
{
  double minDist = COIN_DBL_MAX;
  for (int i = nodeList.size() - 1; i >= 0; --i) {
    minDist = CoinMin(minDist, distance(nodeList.node(i)));
  }
  return minDist;
}

//==============================================================================

bool
CbcHeuristicNode::minDistanceIsSmall(const CbcHeuristicNodeList& nodeList,
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
CbcHeuristicNode::avgDistance(const CbcHeuristicNodeList& nodeList) const
{
  if (nodeList.size() == 0) {
    return COIN_DBL_MAX;
  }
  double sumDist = 0;
  for (int i = nodeList.size() - 1; i >= 0; --i) {
    sumDist += distance(nodeList.node(i));
  }
  return sumDist/nodeList.size();
}

//##############################################################################

// Default Constructor
CbcRounding::CbcRounding() 
  :CbcHeuristic()
{
  // matrix and row copy will automatically be empty
  seed_=1;
  down_ = NULL;
  up_ = NULL;
  equal_ = NULL;
}

// Constructor from model
CbcRounding::CbcRounding(CbcModel & model)
  :CbcHeuristic(model)
{
  // Get a copy of original matrix (and by row for rounding);
  assert(model.solver());
  matrix_ = *model.solver()->getMatrixByCol();
  matrixByRow_ = *model.solver()->getMatrixByRow();
  validate();
  seed_=1;
}

// Destructor 
CbcRounding::~CbcRounding ()
{
  delete [] down_;
  delete [] up_;
  delete [] equal_;
}

// Clone
CbcHeuristic *
CbcRounding::clone() const
{
  return new CbcRounding(*this);
}
// Create C++ lines to get to current state
void 
CbcRounding::generateCpp( FILE * fp) 
{
  CbcRounding other;
  fprintf(fp,"0#include \"CbcHeuristic.hpp\"\n");
  fprintf(fp,"3  CbcRounding rounding(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp,"rounding");
  if (seed_!=other.seed_)
    fprintf(fp,"3  rounding.setSeed(%d);\n",seed_);
  else
    fprintf(fp,"4  rounding.setSeed(%d);\n",seed_);
  fprintf(fp,"3  cbcModel->addHeuristic(&rounding);\n");
}
//#define NEW_ROUNDING
// Copy constructor 
CbcRounding::CbcRounding(const CbcRounding & rhs)
:
  CbcHeuristic(rhs),
  matrix_(rhs.matrix_),
  matrixByRow_(rhs.matrixByRow_),
  seed_(rhs.seed_)
{
#ifdef NEW_ROUNDING
  int numberColumns = matrix_.getNumCols();
  down_ = CoinCopyOfArray(rhs.down_,numberColumns);
  up_ = CoinCopyOfArray(rhs.up_,numberColumns);
  equal_ = CoinCopyOfArray(rhs.equal_,numberColumns);
#else
  down_ = NULL;
  up_ = NULL;
  equal_ = NULL;
#endif  
}

// Assignment operator 
CbcRounding & 
CbcRounding::operator=( const CbcRounding& rhs)
{
  if (this!=&rhs) {
    CbcHeuristic::operator=(rhs);
    matrix_ = rhs.matrix_;
    matrixByRow_ = rhs.matrixByRow_;
#ifdef NEW_ROUNDING
    delete [] down_;
    delete [] up_;
    delete [] equal_;
    int numberColumns = matrix_.getNumCols();
    down_ = CoinCopyOfArray(rhs.down_,numberColumns);
    up_ = CoinCopyOfArray(rhs.up_,numberColumns);
    equal_ = CoinCopyOfArray(rhs.equal_,numberColumns);
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
void 
CbcRounding::resetModel(CbcModel * model)
{
  model_=model;
  // Get a copy of original matrix (and by row for rounding);
  assert(model_->solver());
  matrix_ = *model_->solver()->getMatrixByCol();
  matrixByRow_ = *model_->solver()->getMatrixByRow();
  validate();
}
// See if rounding will give solution
// Sets value of solution
// Assumes rhs for original matrix still okay
// At present only works with integers 
// Fix values if asked for
// Returns 1 if solution, 0 if not
int
CbcRounding::solution(double & solutionValue,
		      double * betterSolution)
{

  // See if to do
  if (!when()||(when()%10==1&&model_->phase()!=1)||
      (when()%10==2&&(model_->phase()!=2&&model_->phase()!=3)))
    return 0; // switched off
  OsiSolverInterface * solver = model_->solver();
  double direction = solver->getObjSense();
  double newSolutionValue = direction*solver->getObjValue();
  return solution(solutionValue,betterSolution,newSolutionValue);
}
// See if rounding will give solution
// Sets value of solution
// Assumes rhs for original matrix still okay
// At present only works with integers 
// Fix values if asked for
// Returns 1 if solution, 0 if not
int
CbcRounding::solution(double & solutionValue,
		      double * betterSolution,
		      double newSolutionValue)
{

  // See if to do
  if (!when()||(when()%10==1&&model_->phase()!=1)||
      (when()%10==2&&(model_->phase()!=2&&model_->phase()!=3)))
    return 0; // switched off
  OsiSolverInterface * solver = model_->solver();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  const double * rowLower = solver->getRowLower();
  const double * rowUpper = solver->getRowUpper();
  const double * solution = solver->getColSolution();
  const double * objective = solver->getObjCoefficients();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance,primalTolerance);

  int numberRows = matrix_.getNumRows();
  assert (numberRows<=solver->getNumRows());
  int numberIntegers = model_->numberIntegers();
  const int * integerVariable = model_->integerVariable();
  int i;
  double direction = solver->getObjSense();
  //double newSolutionValue = direction*solver->getObjValue();
  int returnCode = 0;
  // Column copy
  const double * element = matrix_.getElements();
  const int * row = matrix_.getIndices();
  const CoinBigIndex * columnStart = matrix_.getVectorStarts();
  const int * columnLength = matrix_.getVectorLengths();
  // Row copy
  const double * elementByRow = matrixByRow_.getElements();
  const int * column = matrixByRow_.getIndices();
  const CoinBigIndex * rowStart = matrixByRow_.getVectorStarts();
  const int * rowLength = matrixByRow_.getVectorLengths();

  // Get solution array for heuristic solution
  int numberColumns = solver->getNumCols();
  double * newSolution = new double [numberColumns];
  memcpy(newSolution,solution,numberColumns*sizeof(double));

  double * rowActivity = new double[numberRows];
  memset(rowActivity,0,numberRows*sizeof(double));
  for (i=0;i<numberColumns;i++) {
    int j;
    double value = newSolution[i];
    if (value<lower[i]) {
      value=lower[i];
      newSolution[i]=value;
    } else if (value>upper[i]) {
      value=upper[i];
      newSolution[i]=value;
    }
    if (value) {
      for (j=columnStart[i];
	   j<columnStart[i]+columnLength[i];j++) {
	int iRow=row[j];
	rowActivity[iRow] += value*element[j];
      }
    }
  }
  // check was feasible - if not adjust (cleaning may move)
  for (i=0;i<numberRows;i++) {
    if(rowActivity[i]<rowLower[i]) {
      //assert (rowActivity[i]>rowLower[i]-1000.0*primalTolerance);
      rowActivity[i]=rowLower[i];
    } else if(rowActivity[i]>rowUpper[i]) {
      //assert (rowActivity[i]<rowUpper[i]+1000.0*primalTolerance);
      rowActivity[i]=rowUpper[i];
    }
  }
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integerVariable[i];
    double value=newSolution[iColumn];
    if (fabs(floor(value+0.5)-value)>integerTolerance) {
      double below = floor(value);
      double newValue=newSolution[iColumn];
      double cost = direction * objective[iColumn];
      double move;
      if (cost>0.0) {
	// try up
	move = 1.0 -(value-below);
      } else if (cost<0.0) {
	// try down
	move = below-value;
      } else {
	// won't be able to move unless we can grab another variable
        double randomNumber = randomNumberGenerator_.randomDouble();
	// which way?
        if (randomNumber<0.5) 
          move = below-value;
        else
          move = 1.0 -(value-below);
      }
      newValue += move;
      newSolution[iColumn] = newValue;
      newSolutionValue += move*cost;
      int j;
      for (j=columnStart[iColumn];
	   j<columnStart[iColumn]+columnLength[iColumn];j++) {
	int iRow = row[j];
	rowActivity[iRow] += move*element[j];
      }
    }
  }

  double penalty=0.0;
  const char * integerType = model_->integerType();
  // see if feasible - just using singletons
  for (i=0;i<numberRows;i++) {
    double value = rowActivity[i];
    double thisInfeasibility=0.0;
    if (value<rowLower[i]-primalTolerance)
      thisInfeasibility = value-rowLower[i];
    else if (value>rowUpper[i]+primalTolerance)
      thisInfeasibility = value-rowUpper[i];
    if (thisInfeasibility) {
      // See if there are any slacks I can use to fix up
      // maybe put in coding for multiple slacks?
      double bestCost = 1.0e50;
      int k;
      int iBest=-1;
      double addCost=0.0;
      double newValue=0.0;
      double changeRowActivity=0.0;
      double absInfeasibility = fabs(thisInfeasibility);
      for (k=rowStart[i];k<rowStart[i]+rowLength[i];k++) {
	int iColumn = column[k];
        // See if all elements help
	if (columnLength[iColumn]==1) {
	  double currentValue = newSolution[iColumn];
	  double elementValue = elementByRow[k];
	  double lowerValue = lower[iColumn];
	  double upperValue = upper[iColumn];
	  double gap = rowUpper[i]-rowLower[i];
	  double absElement=fabs(elementValue);
	  if (thisInfeasibility*elementValue>0.0) {
	    // we want to reduce
	    if ((currentValue-lowerValue)*absElement>=absInfeasibility) {
	      // possible - check if integer
	      double distance = absInfeasibility/absElement;
	      double thisCost = -direction*objective[iColumn]*distance;
	      if (integerType[iColumn]) {
		distance = ceil(distance-primalTolerance);
		if (currentValue-distance>=lowerValue-primalTolerance) {
		  if (absInfeasibility-distance*absElement< -gap-primalTolerance)
		    thisCost=1.0e100; // no good
		  else
		    thisCost = -direction*objective[iColumn]*distance;
		} else {
		  thisCost=1.0e100; // no good
		}
	      }
	      if (thisCost<bestCost) {
		bestCost=thisCost;
		iBest=iColumn;
		addCost = thisCost;
		newValue = currentValue-distance;
		changeRowActivity = -distance*elementValue;
	      }
	    }
	  } else {
	    // we want to increase
	    if ((upperValue-currentValue)*absElement>=absInfeasibility) {
	      // possible - check if integer
	      double distance = absInfeasibility/absElement;
	      double thisCost = direction*objective[iColumn]*distance;
	      if (integerType[iColumn]) {
		distance = ceil(distance-1.0e-7);
		assert (currentValue-distance<=upperValue+primalTolerance);
		if (absInfeasibility-distance*absElement< -gap-primalTolerance)
		  thisCost=1.0e100; // no good
		else
		  thisCost = direction*objective[iColumn]*distance;
	      }
	      if (thisCost<bestCost) {
		bestCost=thisCost;
		iBest=iColumn;
		addCost = thisCost;
		newValue = currentValue+distance;
		changeRowActivity = distance*elementValue;
	      }
	    }
	  }
	}
      }
      if (iBest>=0) {
	/*printf("Infeasibility of %g on row %d cost %g\n",
	  thisInfeasibility,i,addCost);*/
	newSolution[iBest]=newValue;
	thisInfeasibility=0.0;
	newSolutionValue += addCost;
	rowActivity[i] += changeRowActivity;
      }
      penalty += fabs(thisInfeasibility);
    }
  }
  if (penalty) {
    // see if feasible using any
    // first continuous
    double penaltyChange=0.0;
    int iColumn;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (integerType[iColumn])
        continue;
      double currentValue = newSolution[iColumn];
      double lowerValue = lower[iColumn];
      double upperValue = upper[iColumn];
      int j;
      int anyBadDown=0;
      int anyBadUp=0;
      double upImprovement=0.0;
      double downImprovement=0.0;
      for (j=columnStart[iColumn];
	   j<columnStart[iColumn]+columnLength[iColumn];j++) {
	int iRow = row[j];
        if (rowUpper[iRow]>rowLower[iRow]) {
          double value = element[j];
          if (rowActivity[iRow]>rowUpper[iRow]+primalTolerance) {
            // infeasible above
            downImprovement += value;
            upImprovement -= value;
            if (value>0.0) 
              anyBadUp++;
            else 
              anyBadDown++;
          } else if (rowActivity[iRow]>rowUpper[iRow]-primalTolerance) {
            // feasible at ub
            if (value>0.0) {
              upImprovement -= value;
              anyBadUp++;
            } else {
              downImprovement += value;
              anyBadDown++;
            }
          } else if (rowActivity[iRow]>rowLower[iRow]+primalTolerance) {
            // feasible in interior
          } else if (rowActivity[iRow]>rowLower[iRow]-primalTolerance) {
            // feasible at lb
            if (value<0.0) {
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
            if (value<0.0) 
              anyBadUp++;
            else 
              anyBadDown++;
          }
        } else {
          // equality row 
          double value = element[j];
          if (rowActivity[iRow]>rowUpper[iRow]+primalTolerance) {
            // infeasible above
            downImprovement += value;
            upImprovement -= value;
            if (value>0.0) 
              anyBadUp++;
            else 
              anyBadDown++;
          } else if (rowActivity[iRow]<rowLower[iRow]-primalTolerance) {
            // infeasible below
            downImprovement -= value;
            upImprovement += value;
            if (value<0.0) 
              anyBadUp++;
            else 
              anyBadDown++;
          } else {
            // feasible - no good
            anyBadUp=-1;
            anyBadDown=-1;
            break;
          }
        }
      }
      // could change tests for anyBad
      if (anyBadUp)
        upImprovement=0.0;
      if (anyBadDown)
        downImprovement=0.0;
      double way=0.0;
      double improvement=0.0;
      if (downImprovement>0.0&&currentValue>lowerValue) {
        way=-1.0;
        improvement = downImprovement;
      } else if (upImprovement>0.0&&currentValue<upperValue) {
        way=1.0;
        improvement = upImprovement;
      }
      if (way) {
        // can improve
        double distance;
        if (way>0.0)
          distance = upperValue-currentValue;
        else
          distance = currentValue-lowerValue;
        for (j=columnStart[iColumn];
             j<columnStart[iColumn]+columnLength[iColumn];j++) {
          int iRow = row[j];
          double value = element[j]*way;
          if (rowActivity[iRow]>rowUpper[iRow]+primalTolerance) {
            // infeasible above
            assert (value<0.0);
            double gap = rowActivity[iRow]-rowUpper[iRow];
            if (gap+value*distance<0.0) 
              distance = -gap/value;
          } else if (rowActivity[iRow]<rowLower[iRow]-primalTolerance) {
            // infeasible below
            assert (value>0.0);
            double gap = rowActivity[iRow]-rowLower[iRow];
            if (gap+value*distance>0.0) 
              distance = -gap/value;
          } else {
            // feasible
            if (value>0) {
              double gap = rowActivity[iRow]-rowUpper[iRow];
              if (gap+value*distance>0.0) 
              distance = -gap/value;
            } else {
              double gap = rowActivity[iRow]-rowLower[iRow];
              if (gap+value*distance<0.0) 
                distance = -gap/value;
            }
          }
        }
        //move
        penaltyChange += improvement*distance;
        distance *= way;
	newSolution[iColumn] += distance;
	newSolutionValue += direction*objective[iColumn]*distance;
        for (j=columnStart[iColumn];
             j<columnStart[iColumn]+columnLength[iColumn];j++) {
          int iRow = row[j];
          double value = element[j];
          rowActivity[iRow] += distance*value;
        }
      }
    }
    // and now all if improving
    double lastChange= penaltyChange ? 1.0 : 0.0;
    while (lastChange>1.0e-2) {
      lastChange=0;
      for (iColumn=0;iColumn<numberColumns;iColumn++) {
        bool isInteger = (integerType[iColumn]!=0);
        double currentValue = newSolution[iColumn];
        double lowerValue = lower[iColumn];
        double upperValue = upper[iColumn];
        int j;
        int anyBadDown=0;
        int anyBadUp=0;
        double upImprovement=0.0;
        double downImprovement=0.0;
        for (j=columnStart[iColumn];
             j<columnStart[iColumn]+columnLength[iColumn];j++) {
          int iRow = row[j];
          double value = element[j];
          if (isInteger) {
            if (value>0.0) {
              if (rowActivity[iRow]+value>rowUpper[iRow]+primalTolerance)
                anyBadUp++;
              if (rowActivity[iRow]-value<rowLower[iRow]-primalTolerance)
                anyBadDown++;
            } else {
              if (rowActivity[iRow]-value>rowUpper[iRow]+primalTolerance)
                anyBadDown++;
              if (rowActivity[iRow]+value<rowLower[iRow]-primalTolerance)
                anyBadUp++;
            }
          }
          if (rowUpper[iRow]>rowLower[iRow]) {
            if (rowActivity[iRow]>rowUpper[iRow]+primalTolerance) {
              // infeasible above
              downImprovement += value;
              upImprovement -= value;
              if (value>0.0) 
                anyBadUp++;
              else 
                anyBadDown++;
            } else if (rowActivity[iRow]>rowUpper[iRow]-primalTolerance) {
              // feasible at ub
              if (value>0.0) {
                upImprovement -= value;
                anyBadUp++;
              } else {
                downImprovement += value;
                anyBadDown++;
              }
            } else if (rowActivity[iRow]>rowLower[iRow]+primalTolerance) {
              // feasible in interior
            } else if (rowActivity[iRow]>rowLower[iRow]-primalTolerance) {
              // feasible at lb
              if (value<0.0) {
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
              if (value<0.0) 
                anyBadUp++;
              else 
                anyBadDown++;
            }
          } else {
            // equality row 
            if (rowActivity[iRow]>rowUpper[iRow]+primalTolerance) {
              // infeasible above
              downImprovement += value;
              upImprovement -= value;
              if (value>0.0) 
                anyBadUp++;
              else 
                anyBadDown++;
            } else if (rowActivity[iRow]<rowLower[iRow]-primalTolerance) {
              // infeasible below
              downImprovement -= value;
              upImprovement += value;
              if (value<0.0) 
                anyBadUp++;
              else 
                anyBadDown++;
            } else {
              // feasible - no good
              anyBadUp=-1;
              anyBadDown=-1;
              break;
            }
          }
        }
        // could change tests for anyBad
        if (anyBadUp)
          upImprovement=0.0;
        if (anyBadDown)
          downImprovement=0.0;
        double way=0.0;
        double improvement=0.0;
        if (downImprovement>0.0&&currentValue>lowerValue) {
          way=-1.0;
          improvement = downImprovement;
        } else if (upImprovement>0.0&&currentValue<upperValue) {
          way=1.0;
          improvement = upImprovement;
        }
        if (way) {
          // can improve
          double distance=COIN_DBL_MAX;
          for (j=columnStart[iColumn];
               j<columnStart[iColumn]+columnLength[iColumn];j++) {
            int iRow = row[j];
            double value = element[j]*way;
            if (rowActivity[iRow]>rowUpper[iRow]+primalTolerance) {
              // infeasible above
              assert (value<0.0);
              double gap = rowActivity[iRow]-rowUpper[iRow];
              if (gap+value*distance<0.0) {
                // If integer then has to move by 1
                if (!isInteger)
                  distance = -gap/value;
                else
                  distance = CoinMax(-gap/value,1.0);
              }
            } else if (rowActivity[iRow]<rowLower[iRow]-primalTolerance) {
              // infeasible below
              assert (value>0.0);
              double gap = rowActivity[iRow]-rowLower[iRow];
              if (gap+value*distance>0.0) {
                // If integer then has to move by 1
                if (!isInteger)
                  distance = -gap/value;
                else
                  distance = CoinMax(-gap/value,1.0);
              }
            } else {
              // feasible
              if (value>0) {
                double gap = rowActivity[iRow]-rowUpper[iRow];
                if (gap+value*distance>0.0) 
                  distance = -gap/value;
              } else {
                double gap = rowActivity[iRow]-rowLower[iRow];
                if (gap+value*distance<0.0) 
                  distance = -gap/value;
              }
            }
          }
          if (isInteger)
            distance = floor(distance+1.05e-8);
          if (!distance) {
            // should never happen
            //printf("zero distance in CbcRounding - debug\n");
          }
          //move
          lastChange += improvement*distance;
          distance *= way;
          newSolution[iColumn] += distance;
          newSolutionValue += direction*objective[iColumn]*distance;
          for (j=columnStart[iColumn];
               j<columnStart[iColumn]+columnLength[iColumn];j++) {
            int iRow = row[j];
            double value = element[j];
            rowActivity[iRow] += distance*value;
          }
        }
      }
      penaltyChange += lastChange;
    }
    penalty -= penaltyChange;
    if (penalty<1.0e-5*fabs(penaltyChange)) {
      // recompute
      penalty=0.0;
      for (i=0;i<numberRows;i++) {
        double value = rowActivity[i];
        if (value<rowLower[i]-primalTolerance)
          penalty += rowLower[i]-value;
        else if (value>rowUpper[i]+primalTolerance)
          penalty += value-rowUpper[i];
      }
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
    int iRandom = (int) (randomNumber*((double) numberIntegers));
    start[0]=iRandom;
    end[0]=numberIntegers;
    start[1]=0;
    end[1]=iRandom;
    for (iPass=0;iPass<2;iPass++) {
      int i;
      for (i=start[iPass];i<end[iPass];i++) {
	int iColumn = integerVariable[i];
#ifndef NDEBUG
	double value=newSolution[iColumn];
	assert (fabs(floor(value+0.5)-value)<integerTolerance);
#endif
	double cost = direction * objective[iColumn];
	double move=0.0;
	if (cost>0.0)
	  move = -1.0;
	else if (cost<0.0)
	  move=1.0;
	while (move) {
	  bool good=true;
	  double newValue=newSolution[iColumn]+move;
	  if (newValue<lower[iColumn]-primalTolerance||
	      newValue>upper[iColumn]+primalTolerance) {
	    move=0.0;
	  } else {
	    // see if we can move
	    int j;
	    for (j=columnStart[iColumn];
		 j<columnStart[iColumn]+columnLength[iColumn];j++) {
	      int iRow = row[j];
	      double newActivity = rowActivity[iRow] + move*element[j];
	      if (newActivity<rowLower[iRow]-primalTolerance||
		  newActivity>rowUpper[iRow]+primalTolerance) {
		good=false;
		break;
	      }
	    }
	    if (good) {
	      newSolution[iColumn] = newValue;
	      newSolutionValue += move*cost;
	      int j;
	      for (j=columnStart[iColumn];
		   j<columnStart[iColumn]+columnLength[iColumn];j++) {
		int iRow = row[j];
		rowActivity[iRow] += move*element[j];
	      }
	    } else {
	      move=0.0;
	    }
	  }
	}
      }
    }
    // Just in case of some stupidity
    double objOffset=0.0;
    solver->getDblParam(OsiObjOffset,objOffset);
    newSolutionValue = -objOffset;
    for ( i=0 ; i<numberColumns ; i++ )
      newSolutionValue += objective[i]*newSolution[i];
    newSolutionValue *= direction;
    //printf("new solution value %g %g\n",newSolutionValue,solutionValue);
    if (newSolutionValue<solutionValue) {
      // paranoid check
      memset(rowActivity,0,numberRows*sizeof(double));
      for (i=0;i<numberColumns;i++) {
	int j;
	double value = newSolution[i];
	if (value) {
	  for (j=columnStart[i];
	       j<columnStart[i]+columnLength[i];j++) {
	    int iRow=row[j];
	    rowActivity[iRow] += value*element[j];
	  }
	}
      }
      // check was approximately feasible
      bool feasible=true;
      for (i=0;i<numberRows;i++) {
	if(rowActivity[i]<rowLower[i]) {
	  if (rowActivity[i]<rowLower[i]-1000.0*primalTolerance)
	    feasible = false;
	} else if(rowActivity[i]>rowUpper[i]) {
	  if (rowActivity[i]>rowUpper[i]+1000.0*primalTolerance)
	    feasible = false;
	}
      }
      if (feasible) {
	// new solution
	memcpy(betterSolution,newSolution,numberColumns*sizeof(double));
	solutionValue = newSolutionValue;
	//printf("** Solution of %g found by rounding\n",newSolutionValue);
	returnCode=1;
      } else {
	// Can easily happen
	//printf("Debug CbcRounding giving bad solution\n");
      }
    }
  }
#ifdef NEW_ROUNDING
  if (!returnCode) {
#if 0
    // back to starting point
    memcpy(newSolution,solution,numberColumns*sizeof(double));
    memset(rowActivity,0,numberRows*sizeof(double));
    for (i=0;i<numberColumns;i++) {
      int j;
      double value = newSolution[i];
      if (value<lower[i]) {
	value=lower[i];
	newSolution[i]=value;
      } else if (value>upper[i]) {
	value=upper[i];
	newSolution[i]=value;
      }
      if (value) {
	for (j=columnStart[i];
	     j<columnStart[i]+columnLength[i];j++) {
	  int iRow=row[j];
	  rowActivity[iRow] += value*element[j];
	}
      }
    }
    // check was feasible - if not adjust (cleaning may move)
    for (i=0;i<numberRows;i++) {
      if(rowActivity[i]<rowLower[i]) {
	//assert (rowActivity[i]>rowLower[i]-1000.0*primalTolerance);
	rowActivity[i]=rowLower[i];
      } else if(rowActivity[i]>rowUpper[i]) {
	//assert (rowActivity[i]<rowUpper[i]+1000.0*primalTolerance);
	rowActivity[i]=rowUpper[i];
      }
    }
#endif
    int * candidate = new int [numberColumns];
    int nCandidate=0;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      bool isInteger = (integerType[iColumn]!=0);
      if (isInteger) {
	double currentValue = newSolution[iColumn];
	if (fabs(currentValue-floor(currentValue+0.5))>1.0e-8)
	  candidate[nCandidate++]=iColumn;
      }
    }
    if (true) {
      // Rounding as in Berthold
      while (nCandidate) {
	double infeasibility =1.0e-7;
	int iRow=-1;
	for (i=0;i<numberRows;i++) {
	  double value=0.0;
	  if(rowActivity[i]<rowLower[i]) {
	    value = rowLower[i]-rowActivity[i];
	  } else if(rowActivity[i]>rowUpper[i]) {
	    value = rowActivity[i]-rowUpper[i];
	  }
	  if (value>infeasibility) {
	    infeasibility = value;
	    iRow=i;
	  }
	}
	if (iRow>=0) {
	  // infeasible
	} else {
	  // feasible
	}
      }
    } else {
      // Shifting as in Berthold
    }
    delete [] candidate;
  }
#endif
  delete [] newSolution;
  delete [] rowActivity;
  return returnCode;
}
// update model
void CbcRounding::setModel(CbcModel * model)
{
  model_ = model;
  // Get a copy of original matrix (and by row for rounding);
  assert(model_->solver());
  matrix_ = *model_->solver()->getMatrixByCol();
  matrixByRow_ = *model_->solver()->getMatrixByRow();
  // make sure model okay for heuristic
  validate();
}
// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
void 
CbcRounding::validate() 
{
  if (model_&&when()<10) {
    if (model_->numberIntegers()!=
        model_->numberObjects()&&(model_->numberObjects()||
				  (model_->specialOptions()&1024)==0))
      setWhen(0);
  }
#ifdef NEW_ROUNDING
  int numberColumns = matrix_.getNumCols();
  down_ = new unsigned short [numberColumns];
  up_ = new unsigned short [numberColumns];
  equal_ = new unsigned short [numberColumns];
  // Column copy
  const double * element = matrix_.getElements();
  const int * row = matrix_.getIndices();
  const CoinBigIndex * columnStart = matrix_.getVectorStarts();
  const int * columnLength = matrix_.getVectorLengths();
  const double * rowLower = model.solver()->getRowLower();
  const double * rowUpper = model.solver()->getRowUpper();
  for (int i=0;i<numberColumns;i++) {
    int down=0;
    int up=0;
    int equal=0;
    if (columnLength[i]>65535) {
      equal[0]=65535; 
      break; // unlikely to work
    }
    for (CoinBigIndex j=columnStart[i];
	 j<columnStart[i]+columnLength[i];j++) {
      int iRow=row[j];
      if (rowLower[iRow]>-1.0e20&&rowUpper[iRow]<1.0e20) {
	equal++;
      } else if (element[j]>0.0) {
	if (rowUpper[iRow]<1.0e20)
	  up++;
	else
	  down--;
      } else {
	if (rowLower[iRow]>-1.0e20)
	  up++;
	else
	  down--;
      }
    }
    down_[i] = (unsigned short) down;
    up_[i] = (unsigned short) up;
    equal_[i] = (unsigned short) equal;
  }
#else
  down_ = NULL;
  up_ = NULL;
  equal_ = NULL;
#endif  
}

// Default Constructor
CbcHeuristicPartial::CbcHeuristicPartial() 
  :CbcHeuristic()
{
  fixPriority_ = 10000;
}

// Constructor from model
CbcHeuristicPartial::CbcHeuristicPartial(CbcModel & model, int fixPriority, int numberNodes)
  :CbcHeuristic(model)
{
  fixPriority_ = fixPriority;
  setNumberNodes(numberNodes);
  validate();
}

// Destructor 
CbcHeuristicPartial::~CbcHeuristicPartial ()
{
}

// Clone
CbcHeuristic *
CbcHeuristicPartial::clone() const
{
  return new CbcHeuristicPartial(*this);
}
// Create C++ lines to get to current state
void 
CbcHeuristicPartial::generateCpp( FILE * fp) 
{
  CbcHeuristicPartial other;
  fprintf(fp,"0#include \"CbcHeuristic.hpp\"\n");
  fprintf(fp,"3  CbcHeuristicPartial partial(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp,"partial");
  if (fixPriority_!=other.fixPriority_)
    fprintf(fp,"3  partial.setFixPriority(%d);\n",fixPriority_);
  else
    fprintf(fp,"4  partial.setFixPriority(%d);\n",fixPriority_);
  fprintf(fp,"3  cbcModel->addHeuristic(&partial);\n");
}
//#define NEW_PARTIAL
// Copy constructor 
CbcHeuristicPartial::CbcHeuristicPartial(const CbcHeuristicPartial & rhs)
:
  CbcHeuristic(rhs),
  fixPriority_(rhs.fixPriority_)
{
}

// Assignment operator 
CbcHeuristicPartial & 
CbcHeuristicPartial::operator=( const CbcHeuristicPartial& rhs)
{
  if (this!=&rhs) {
    CbcHeuristic::operator=(rhs);
    fixPriority_ = rhs.fixPriority_;
  }
  return *this;
}

// Resets stuff if model changes
void 
CbcHeuristicPartial::resetModel(CbcModel * model)
{
  model_=model;
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
int
CbcHeuristicPartial::solution(double & solutionValue,
		      double * betterSolution)
{
  // Return if already done
  if (fixPriority_<0)
    return 0; // switched off
  const double * hotstartSolution = model_->hotstartSolution();
  const int * hotstartPriorities = model_->hotstartPriorities();
  if (!hotstartSolution)
    return 0;
  OsiSolverInterface * solver = model_->solver();
  
  int numberIntegers = model_->numberIntegers();
  const int * integerVariable = model_->integerVariable();
  
  OsiSolverInterface * newSolver = model_->continuousSolver()->clone();
  const double * colLower = newSolver->getColLower();
  const double * colUpper = newSolver->getColUpper();

  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance,primalTolerance);
    
  int i;
  int numberFixed=0;
  int returnCode=0;

  for (i=0;i<numberIntegers;i++) {
    int iColumn=integerVariable[i];
    if (abs(hotstartPriorities[iColumn])<=fixPriority_) {
      double value = hotstartSolution[iColumn];
      double lower = colLower[iColumn];
      double upper = colUpper[iColumn];
      value = CoinMax(value,lower);
      value = CoinMin(value,upper);
      if (fabs(value-floor(value+0.5))<1.0e-8) {
	value = floor(value+0.5);
	newSolver->setColLower(iColumn,value);
	newSolver->setColUpper(iColumn,value);
	numberFixed++;
      }
    }
  }
  if (numberFixed>numberIntegers/5-100000000) {
#ifdef COIN_DEVELOP
    printf("%d integers fixed\n",numberFixed);
#endif
    returnCode = smallBranchAndBound(newSolver,numberNodes_,betterSolution,solutionValue,
				     model_->getCutoff(),"CbcHeuristicPartial");
    if (returnCode<0)
      returnCode=0; // returned on size
    //printf("return code %d",returnCode);
    if ((returnCode&2)!=0) {
      // could add cut
      returnCode &= ~2;
      //printf("could add cut with %d elements (if all 0-1)\n",nFix);
    } else {
      //printf("\n");
    }
  }
  fixPriority_=-1; // switch off
  
  delete newSolver;
  return returnCode;
}
// update model
void CbcHeuristicPartial::setModel(CbcModel * model)
{
  model_ = model;
  assert(model_->solver());
  // make sure model okay for heuristic
  validate();
}
// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
void 
CbcHeuristicPartial::validate() 
{
  if (model_&&when()<10) {
    if (model_->numberIntegers()!=
        model_->numberObjects())
      setWhen(0);
  }
}

// Default Constructor
CbcSerendipity::CbcSerendipity() 
  :CbcHeuristic()
{
}

// Constructor from model
CbcSerendipity::CbcSerendipity(CbcModel & model)
  :CbcHeuristic(model)
{
}

// Destructor 
CbcSerendipity::~CbcSerendipity ()
{
}

// Clone
CbcHeuristic *
CbcSerendipity::clone() const
{
  return new CbcSerendipity(*this);
}
// Create C++ lines to get to current state
void 
CbcSerendipity::generateCpp( FILE * fp) 
{
  fprintf(fp,"0#include \"CbcHeuristic.hpp\"\n");
  fprintf(fp,"3  CbcSerendipity serendipity(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp,"serendipity");
  fprintf(fp,"3  cbcModel->addHeuristic(&serendipity);\n");
}

// Copy constructor 
CbcSerendipity::CbcSerendipity(const CbcSerendipity & rhs)
:
  CbcHeuristic(rhs)
{
}

// Assignment operator 
CbcSerendipity & 
CbcSerendipity::operator=( const CbcSerendipity& rhs)
{
  if (this!=&rhs) {
    CbcHeuristic::operator=(rhs);
  }
  return *this;
}

// Returns 1 if solution, 0 if not
int
CbcSerendipity::solution(double & solutionValue,
			 double * betterSolution)
{
  if (!model_)
    return 0;
  if (!inputSolution_) {
    // get information on solver type
    OsiAuxInfo * auxInfo = model_->solver()->getAuxiliaryInfo();
    OsiBabSolver * auxiliaryInfo = dynamic_cast< OsiBabSolver *> (auxInfo);
    if (auxiliaryInfo) {
      return auxiliaryInfo->solution(solutionValue,betterSolution,model_->solver()->getNumCols());
    } else {
      return 0;
    }
  } else {
    int numberColumns = model_->getNumCols();
    double value =inputSolution_[numberColumns];
    int returnCode=0;
    if (value<solutionValue) {
      solutionValue = value;
      memcpy(betterSolution,inputSolution_,numberColumns*sizeof(double));
      returnCode=1;
    }
    delete [] inputSolution_;
    inputSolution_=NULL;
    model_ = NULL; // switch off
    return returnCode;
  }
}
// update model
void CbcSerendipity::setModel(CbcModel * model)
{
  model_ = model;
}
// Resets stuff if model changes
void 
CbcSerendipity::resetModel(CbcModel * model)
{
  model_ = model;
}
  
