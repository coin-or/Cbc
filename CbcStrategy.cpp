// Copyright (C) 2005, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cmath>
#include <cfloat>

#include "OsiSolverInterface.hpp"
#ifdef COIN_USE_CLP
#include "OsiClpSolverInterface.hpp"
#endif
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcStrategy.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcBranchActual.hpp"
#include "CbcNode.hpp"
#include "CoinWarmStart.hpp"
#include "CglPreProcess.hpp"
// Cuts

#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"

// Heuristics

#include "CbcHeuristic.hpp"

// Default Constructor
CbcStrategy::CbcStrategy() 
  :depth_(0),
   preProcessState_(0),
   process_(NULL)
{
}

// Destructor 
CbcStrategy::~CbcStrategy ()
{
  delete process_;
}
// Delete pre-processing object to save memory
void 
CbcStrategy::deletePreProcess()
{ 
  delete process_;
  process_=NULL;
}
// Return a new Full node information pointer (descendant of CbcFullNodeInfo)
CbcNodeInfo * 
CbcStrategy::fullNodeInfo(CbcModel * model,int numberRowsAtContinuous) const
{
  return new CbcFullNodeInfo(model,numberRowsAtContinuous);
}
// Return a new Partial node information pointer (descendant of CbcPartialNodeInfo)
CbcNodeInfo * 
CbcStrategy::partialNodeInfo(CbcModel * model, CbcNodeInfo * parent, CbcNode * owner,
                             int numberChangedBounds,const int * variables,
                             const double * boundChanges,
                             const CoinWarmStartDiff *basisDiff) const
{
  return new CbcPartialNodeInfo(parent, owner, numberChangedBounds, variables,
                            boundChanges,basisDiff);
}
/* After a CbcModel::resolve this can return a status
   -1 no effect
   0 treat as optimal
   1 as 0 but do not do any more resolves (i.e. no more cuts)
   2 treat as infeasible
*/
int
CbcStrategy::status(CbcModel * model, CbcNodeInfo * parent,int whereFrom)
{
  return -1;
}

// Default Constructor
CbcStrategyDefault::CbcStrategyDefault(bool cutsOnlyAtRoot,
                                       int numberStrong,
                                       int numberBeforeTrust,
                                       int printLevel)
  :CbcStrategy(),
   cutsOnlyAtRoot_(cutsOnlyAtRoot),
   numberStrong_(numberStrong),
   numberBeforeTrust_(numberBeforeTrust),
   printLevel_(printLevel),
   desiredPreProcess_(0),
   preProcessPasses_(0)
{
}


// Destructor 
CbcStrategyDefault::~CbcStrategyDefault ()
{
}

// Clone
CbcStrategy *
CbcStrategyDefault::clone() const
{
  return new CbcStrategyDefault(*this);
}

// Copy constructor 
CbcStrategyDefault::CbcStrategyDefault(const CbcStrategyDefault & rhs)
:
  CbcStrategy(rhs),
  cutsOnlyAtRoot_(rhs.cutsOnlyAtRoot_),
  numberStrong_(rhs.numberStrong_),
  numberBeforeTrust_(rhs.numberBeforeTrust_),
  printLevel_(rhs.printLevel_),
  desiredPreProcess_(rhs.desiredPreProcess_),
  preProcessPasses_(rhs.preProcessPasses_)
{
  setNested(rhs.getNested());
}

// Setup cut generators
void 
CbcStrategyDefault::setupCutGenerators(CbcModel & model)
{
  // Set up some cut generators and defaults
  // Probing first as gets tight bounds on continuous

  CglProbing generator1;
  generator1.setUsingObjective(true);
  generator1.setMaxPass(1);
  // Number of unsatisfied variables to look at
  generator1.setMaxProbe(10);
  // How far to follow the consequences
  generator1.setMaxLook(10);
  // Only look at rows with fewer than this number of elements
  generator1.setMaxElements(200);
  //generator1.setRowCuts(3);

  CglGomory generator2;
  // try larger limit
  generator2.setLimit(300);

  CglKnapsackCover generator3;

  //CglOddHole generator4;
  //generator4.setMinimumViolation(0.005);
  //generator4.setMinimumViolationPer(0.00002);
  // try larger limit
  //generator4.setMaximumEntries(200);

  CglClique generator5;
  generator5.setStarCliqueReport(false);
  generator5.setRowCliqueReport(false);

  CglMixedIntegerRounding2 mixedGen;
  CglFlowCover flowGen;
  
  // Add in generators
  int setting = cutsOnlyAtRoot_ ? -99 : -1;
  int numberGenerators = model.numberCutGenerators();
  int iGenerator;
  bool found;
  found=false;
  for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
    CglCutGenerator * generator = model.cutGenerator(iGenerator)->generator();
    CglProbing * cgl = dynamic_cast<CglProbing *>(generator);
    if (cgl) {
      found=true;
      break;
    }
  }
  if (!found)
    model.addCutGenerator(&generator1,setting,"Probing");
  found=false;
  for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
    CglCutGenerator * generator = model.cutGenerator(iGenerator)->generator();
    CglGomory * cgl = dynamic_cast<CglGomory *>(generator);
    if (cgl) {
      found=true;
      break;
    }
  }
  if (!found)
  model.addCutGenerator(&generator2,setting,"Gomory");
  found=false;
  for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
    CglCutGenerator * generator = model.cutGenerator(iGenerator)->generator();
    CglKnapsackCover * cgl = dynamic_cast<CglKnapsackCover *>(generator);
    if (cgl) {
      found=true;
      break;
    }
  }
  if (!found)
    model.addCutGenerator(&generator3,setting,"Knapsack");
  //model.addCutGenerator(&generator4,setting,"OddHole");
  found=false;
  for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
    CglCutGenerator * generator = model.cutGenerator(iGenerator)->generator();
    CglClique * cgl = dynamic_cast<CglClique *>(generator);
    if (cgl) {
      found=true;
      break;
    }
  }
  if (!found)
    model.addCutGenerator(&generator5,setting,"Clique");
  found=false;
  for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
    CglCutGenerator * generator = model.cutGenerator(iGenerator)->generator();
    CglFlowCover * cgl = dynamic_cast<CglFlowCover *>(generator);
    if (cgl) {
      found=true;
      break;
    }
  }
  if (!found)
    model.addCutGenerator(&flowGen,setting,"FlowCover");
  found=false;
  for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
    CglCutGenerator * generator = model.cutGenerator(iGenerator)->generator();
    CglMixedIntegerRounding2 * cgl = dynamic_cast<CglMixedIntegerRounding2 *>(generator);
    if (cgl) {
      found=true;
      break;
    }
  }
  if (!found)
    model.addCutGenerator(&mixedGen,setting,"MixedIntegerRounding2");
  // Say we want timings
  int newNumberGenerators = model.numberCutGenerators();
  for (iGenerator=numberGenerators;iGenerator<newNumberGenerators;iGenerator++) {
    CbcCutGenerator * generator = model.cutGenerator(iGenerator);
    generator->setTiming(true);
  }
  if (model.getNumCols()<500)
    model.setMaximumCutPassesAtRoot(-100); // always do 100 if possible
  else if (model.getNumCols()<5000)
    model.setMaximumCutPassesAtRoot(100); // use minimum drop
  else
    model.setMaximumCutPassesAtRoot(20);
}
// Setup heuristics
void 
CbcStrategyDefault::setupHeuristics(CbcModel & model)
{
  // Allow rounding heuristic

  CbcRounding heuristic1(model);
  int numberHeuristics = model.numberHeuristics();
  int iHeuristic;
  bool found;
  found=false;
  for (iHeuristic=0;iHeuristic<numberHeuristics;iHeuristic++) {
    CbcHeuristic * heuristic = model.heuristic(iHeuristic);
    CbcRounding * cgl = dynamic_cast<CbcRounding *>(heuristic);
    if (cgl) {
      found=true;
      break;
    }
  }
  if (!found)
    model.addHeuristic(&heuristic1);
}
// Do printing stuff
void 
CbcStrategyDefault::setupPrinting(CbcModel & model,int modelLogLevel)
{
  if (!modelLogLevel) {
    model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
    model.messageHandler()->setLogLevel(0);
    model.solver()->messageHandler()->setLogLevel(0);
  } else if (modelLogLevel==1) {
    model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
    model.messageHandler()->setLogLevel(1);
    model.solver()->messageHandler()->setLogLevel(0);
  } else {
    model.messageHandler()->setLogLevel(2);
    model.solver()->messageHandler()->setLogLevel(1);
    model.setPrintFrequency(50);
  }
}
// Other stuff e.g. strong branching
void 
CbcStrategyDefault::setupOther(CbcModel & model)
{
  // See if preprocessing wanted
  if (desiredPreProcess_) {
    delete process_;
    // solver_ should have been cloned outside
    CglPreProcess * process = new CglPreProcess();
    OsiSolverInterface * solver = model.solver();
    int logLevel = model.messageHandler()->logLevel();
#ifdef COIN_USE_CLP
    OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
    ClpSimplex * lpSolver=NULL;
    if (clpSolver) {
      if (clpSolver->messageHandler()->logLevel())
        clpSolver->messageHandler()->setLogLevel(1);
      if (logLevel>-1)
        clpSolver->messageHandler()->setLogLevel(CoinMin(logLevel,clpSolver->messageHandler()->logLevel()));
      lpSolver = clpSolver->getModelPtr();
      /// If user left factorization frequency then compute
      lpSolver->defaultFactorizationFrequency();
    }
#endif
    // Tell solver we are in Branch and Cut
    solver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo) ;
    // Default set of cut generators
    CglProbing generator1;
    generator1.setUsingObjective(true);
    generator1.setMaxPass(3);
    generator1.setMaxProbeRoot(solver->getNumCols());
    generator1.setMaxElements(100);
    generator1.setMaxLookRoot(50);
    generator1.setRowCuts(3);
    //generator1.messageHandler()->setLogLevel(logLevel);
    process->messageHandler()->setLogLevel(logLevel);
    // Add in generators
    process->addCutGenerator(&generator1);
    int translate[]={9999,0,2,3};
    OsiSolverInterface * solver2 = 
      process->preProcessNonDefault(*solver,
                                    translate[desiredPreProcess_],preProcessPasses_);
    // Tell solver we are not in Branch and Cut
    solver->setHintParam(OsiDoInBranchAndCut,false,OsiHintDo) ;
    if (solver2)
      solver2->setHintParam(OsiDoInBranchAndCut,false,OsiHintDo) ;
    bool feasible=true;
    if (!solver2) {
      feasible = false;
      //printf("Pre-processing says infeasible\n");
      delete process;
      preProcessState_=-1;
      process_=NULL;
    } else {
      // now tighten bounds
#ifdef COIN_USE_CLP
      if (clpSolver) {
        // model has changed
        solver = model.solver();
        OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
        ClpSimplex * lpSolver = clpSolver->getModelPtr();
        if (lpSolver->tightenPrimalBounds()==0) {
          lpSolver->dual();
        } else {
          feasible = false;
        }
      }
#endif
      if (feasible) {
        preProcessState_=1;
        process_=process;
        /* Note that original solver will be kept (with false)
           and that final solver will also be kept.
           This is for post-processing
        */
        OsiSolverInterface * solver3 = solver2->clone();
        model.assignSolver(solver3,false);
        if (process_->numberSOS()) {
          int numberSOS = process_->numberSOS();
          int numberIntegers = model.numberIntegers();
          /* model may not have created objects
             If none then create
          */
          if (!numberIntegers||!model.numberObjects()) {
            model.findIntegers(true);
            numberIntegers = model.numberIntegers();
          }
          CbcObject ** oldObjects = model.objects();
          // Do sets and priorities
          CbcObject ** objects = new CbcObject * [numberSOS];
          // set old objects to have low priority
          int numberOldObjects = model.numberObjects();
          int numberColumns = model.getNumCols();
          for (int iObj = 0;iObj<numberOldObjects;iObj++) {
            int oldPriority = oldObjects[iObj]->priority();
            oldObjects[iObj]->setPriority(numberColumns+oldPriority);
          }
          const int * starts = process_->startSOS();
          const int * which = process_->whichSOS();
          const int * type = process_->typeSOS();
          const double * weight = process_->weightSOS();
          int iSOS;
          for (iSOS =0;iSOS<numberSOS;iSOS++) {
            int iStart = starts[iSOS];
            int n=starts[iSOS+1]-iStart;
            objects[iSOS] = new CbcSOS(&model,n,which+iStart,weight+iStart,
                                       iSOS,type[iSOS]);
            // branch on long sets first
            objects[iSOS]->setPriority(numberColumns-n);
          }
          model.addObjects(numberSOS,objects);
          for (iSOS=0;iSOS<numberSOS;iSOS++)
            delete objects[iSOS];
          delete [] objects;
        }
      } else {
        //printf("Pre-processing says infeasible\n");
        delete process;
        preProcessState_=-1;
        process_=NULL;
      }
    }
  }
  model.setNumberStrong(numberStrong_);
  model.setNumberBeforeTrust(numberBeforeTrust_);
}
// Default Constructor
CbcStrategyDefaultSubTree::CbcStrategyDefaultSubTree(CbcModel * parent ,
                                                     bool cutsOnlyAtRoot,
                                       int numberStrong,
                                       int numberBeforeTrust,
                                       int printLevel)
  :CbcStrategy(),
   parentModel_(parent),
   cutsOnlyAtRoot_(cutsOnlyAtRoot),
   numberStrong_(numberStrong),
   numberBeforeTrust_(numberBeforeTrust),
   printLevel_(printLevel)
{
}


// Destructor 
CbcStrategyDefaultSubTree::~CbcStrategyDefaultSubTree ()
{
}

// Clone
CbcStrategy *
CbcStrategyDefaultSubTree::clone() const
{
  return new CbcStrategyDefaultSubTree(*this);
}

// Copy constructor 
CbcStrategyDefaultSubTree::CbcStrategyDefaultSubTree(const CbcStrategyDefaultSubTree & rhs)
:
  CbcStrategy(rhs),
  parentModel_(rhs.parentModel_),
  cutsOnlyAtRoot_(rhs.cutsOnlyAtRoot_),
  numberStrong_(rhs.numberStrong_),
  numberBeforeTrust_(rhs.numberBeforeTrust_),
  printLevel_(rhs.printLevel_)
{
  setNested(rhs.getNested());
}

// Setup cut generators
void 
CbcStrategyDefaultSubTree::setupCutGenerators(CbcModel & model)
{
  // Set up some cut generators and defaults
  // Probing first as gets tight bounds on continuous

  CglProbing generator1;
  generator1.setUsingObjective(true);
  generator1.setMaxPass(1);
  // Number of unsatisfied variables to look at
  generator1.setMaxProbe(10);
  // How far to follow the consequences
  generator1.setMaxLook(10);
  // Only look at rows with fewer than this number of elements
  generator1.setMaxElements(200);
  //generator1.setRowCuts(3);

  CglGomory generator2;
  // try larger limit
  generator2.setLimit(300);

  CglKnapsackCover generator3;

  //CglOddHole generator4;
  //generator4.setMinimumViolation(0.005);
  //generator4.setMinimumViolationPer(0.00002);
  // try larger limit
  //generator4.setMaximumEntries(200);

  CglClique generator5;
  generator5.setStarCliqueReport(false);
  generator5.setRowCliqueReport(false);

  CglMixedIntegerRounding2 mixedGen;
  CglFlowCover flowGen;
  
  // Add in generators
  int setting = cutsOnlyAtRoot_ ? -99 : -1;
  int numberGenerators = model.numberCutGenerators();
  int numberParentGenerators = parentModel_->numberCutGenerators();
  int iGenerator;
  bool found;
  found=false;
  int howOften=0;
  for (iGenerator=0;iGenerator<numberParentGenerators;iGenerator++) {
    CglCutGenerator * generator = parentModel_->cutGenerator(iGenerator)->generator();
    CglProbing * cgl = dynamic_cast<CglProbing *>(generator);
    if (cgl) {
      found=true;
      howOften = parentModel_->cutGenerator(iGenerator)->howOften();
      break;
    }
  }
  if (found&&howOften>=0) {
    found=false;
    for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
      CglCutGenerator * generator = model.cutGenerator(iGenerator)->generator();
      CglProbing * cgl = dynamic_cast<CglProbing *>(generator);
      if (cgl) {
        found=true;
        break;
      }
    }
    if (!found)
      model.addCutGenerator(&generator1,setting,"Probing");
  }
  found=false;
  for (iGenerator=0;iGenerator<numberParentGenerators;iGenerator++) {
    CglCutGenerator * generator = parentModel_->cutGenerator(iGenerator)->generator();
    CglGomory * cgl = dynamic_cast<CglGomory *>(generator);
    if (cgl) {
      found=true;
      howOften = parentModel_->cutGenerator(iGenerator)->howOften();
      break;
    }
  }
  if (found&&howOften>=0) {
    found=false;
    for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
      CglCutGenerator * generator = model.cutGenerator(iGenerator)->generator();
      CglGomory * cgl = dynamic_cast<CglGomory *>(generator);
      if (cgl) {
        found=true;
        break;
      }
    }
    if (!found)
      model.addCutGenerator(&generator2,setting,"Gomory");
  }
  found=false;
  for (iGenerator=0;iGenerator<numberParentGenerators;iGenerator++) {
    CglCutGenerator * generator = parentModel_->cutGenerator(iGenerator)->generator();
    CglKnapsackCover * cgl = dynamic_cast<CglKnapsackCover *>(generator);
    if (cgl) {
      found=true;
      howOften = parentModel_->cutGenerator(iGenerator)->howOften();
      break;
    }
  }
  if (found&&howOften>=0) {
    found=false;
    for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
      CglCutGenerator * generator = model.cutGenerator(iGenerator)->generator();
      CglKnapsackCover * cgl = dynamic_cast<CglKnapsackCover *>(generator);
      if (cgl) {
        found=true;
        break;
      }
    }
    if (!found)
      model.addCutGenerator(&generator3,setting,"Knapsack");
  }
  found=false;
  for (iGenerator=0;iGenerator<numberParentGenerators;iGenerator++) {
    CglCutGenerator * generator = parentModel_->cutGenerator(iGenerator)->generator();
    CglClique * cgl = dynamic_cast<CglClique *>(generator);
    if (cgl) {
      found=true;
      howOften = parentModel_->cutGenerator(iGenerator)->howOften();
      break;
    }
  }
  if (found&&howOften>=0) {
    found=false;
    for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
      CglCutGenerator * generator = model.cutGenerator(iGenerator)->generator();
      CglClique * cgl = dynamic_cast<CglClique *>(generator);
      if (cgl) {
        found=true;
        break;
      }
    }
    if (!found)
      model.addCutGenerator(&generator5,setting,"Clique");
  }
  found=false;
  for (iGenerator=0;iGenerator<numberParentGenerators;iGenerator++) {
    CglCutGenerator * generator = parentModel_->cutGenerator(iGenerator)->generator();
    CglFlowCover * cgl = dynamic_cast<CglFlowCover *>(generator);
    if (cgl) {
      found=true;
      howOften = parentModel_->cutGenerator(iGenerator)->howOften();
      break;
    }
  }
  if (found&&howOften>=0) {
    found=false;
    for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
      CglCutGenerator * generator = model.cutGenerator(iGenerator)->generator();
      CglFlowCover * cgl = dynamic_cast<CglFlowCover *>(generator);
      if (cgl) {
        found=true;
        break;
      }
    }
    if (!found)
      model.addCutGenerator(&flowGen,setting,"FlowCover");
    found=false;
  }
  for (iGenerator=0;iGenerator<numberParentGenerators;iGenerator++) {
    CglCutGenerator * generator = parentModel_->cutGenerator(iGenerator)->generator();
    CglMixedIntegerRounding2 * cgl = dynamic_cast<CglMixedIntegerRounding2 *>(generator);
    if (cgl) {
      found=true;
      howOften = parentModel_->cutGenerator(iGenerator)->howOften();
      break;
    }
  }
  if (found&&howOften>=0) {
    found=false;
    for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
      CglCutGenerator * generator = model.cutGenerator(iGenerator)->generator();
      CglMixedIntegerRounding2 * cgl = dynamic_cast<CglMixedIntegerRounding2 *>(generator);
      if (cgl) {
        found=true;
        break;
      }
    }
    if (!found)
      model.addCutGenerator(&mixedGen,setting,"MixedIntegerRounding2");
  }
#if 0
  // Say we want timings
  int newNumberGenerators = model.numberCutGenerators();
  for (iGenerator=numberGenerators;iGenerator<newNumberGenerators;iGenerator++) {
    CbcCutGenerator * generator = model.cutGenerator(iGenerator);
    generator->setTiming(true);
  }
#endif
  if (model.getNumCols()<500)
    model.setMaximumCutPassesAtRoot(-100); // always do 100 if possible
  else if (model.getNumCols()<5000)
    model.setMaximumCutPassesAtRoot(100); // use minimum drop
  else
    model.setMaximumCutPassesAtRoot(20);
}
// Setup heuristics
void 
CbcStrategyDefaultSubTree::setupHeuristics(CbcModel & model)
{
  // Allow rounding heuristic

  CbcRounding heuristic1(model);
  int numberHeuristics = model.numberHeuristics();
  int iHeuristic;
  bool found;
  found=false;
  for (iHeuristic=0;iHeuristic<numberHeuristics;iHeuristic++) {
    CbcHeuristic * heuristic = model.heuristic(iHeuristic);
    CbcRounding * cgl = dynamic_cast<CbcRounding *>(heuristic);
    if (cgl) {
      found=true;
      break;
    }
  }
  if (!found)
    model.addHeuristic(&heuristic1);
}
// Do printing stuff
void 
CbcStrategyDefaultSubTree::setupPrinting(CbcModel & model,int modelLogLevel)
{
  if (!modelLogLevel) {
    model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
    model.messageHandler()->setLogLevel(0);
    model.solver()->messageHandler()->setLogLevel(0);
  } else if (modelLogLevel==1) {
    model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
    model.messageHandler()->setLogLevel(1);
    model.solver()->messageHandler()->setLogLevel(0);
  } else {
    model.messageHandler()->setLogLevel(2);
    model.solver()->messageHandler()->setLogLevel(1);
    model.setPrintFrequency(50);
  }
}
// Other stuff e.g. strong branching
void 
CbcStrategyDefaultSubTree::setupOther(CbcModel & model)
{
  model.setNumberStrong(numberStrong_);
  model.setNumberBeforeTrust(numberBeforeTrust_);
}


  
