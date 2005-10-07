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
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcStrategy.hpp"
#include "CbcCutGenerator.hpp"
// Cuts

#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding.hpp"

// Heuristics

#include "CbcHeuristic.hpp"

// Default Constructor
CbcStrategy::CbcStrategy() 
  :depth_(0)
{
}

// Destructor 
CbcStrategy::~CbcStrategy ()
{
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
   printLevel_(printLevel)
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
  printLevel_(rhs.printLevel_)
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

  CglMixedIntegerRounding mixedGen;
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
    CglMixedIntegerRounding * cgl = dynamic_cast<CglMixedIntegerRounding *>(generator);
    if (cgl) {
      found=true;
      break;
    }
  }
  if (!found)
    model.addCutGenerator(&mixedGen,setting,"MixedIntegerRounding");
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

  CglMixedIntegerRounding mixedGen;
  CglFlowCover flowGen;
  
  // Add in generators
  int setting = cutsOnlyAtRoot_ ? -99 : -1;
  int numberGenerators = model.numberCutGenerators();
  int numberParentGenerators = parentModel_->numberCutGenerators();
  int iGenerator;
  bool found;
  found=false;
  for (iGenerator=0;iGenerator<numberParentGenerators;iGenerator++) {
    CglCutGenerator * generator = parentModel_->cutGenerator(iGenerator)->generator();
    CglProbing * cgl = dynamic_cast<CglProbing *>(generator);
    if (cgl) {
      found=true;
      break;
    }
  }
  if (found) {
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
      break;
    }
  }
  if (found) {
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
      break;
    }
  }
  if (found) {
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
      break;
    }
  }
  if (found) {
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
      break;
    }
  }
  if (found) {
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
    CglMixedIntegerRounding * cgl = dynamic_cast<CglMixedIntegerRounding *>(generator);
    if (cgl) {
      found=true;
      break;
    }
  }
  if (found) {
    found=false;
    for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
      CglCutGenerator * generator = model.cutGenerator(iGenerator)->generator();
      CglMixedIntegerRounding * cgl = dynamic_cast<CglMixedIntegerRounding *>(generator);
      if (cgl) {
        found=true;
        break;
      }
    }
    if (!found)
      model.addCutGenerator(&mixedGen,setting,"MixedIntegerRounding");
  }
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


  
