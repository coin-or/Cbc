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
#include "CglProbing.hpp"
#include "CoinTime.hpp"

// Default Constructor 
CbcCutGenerator::CbcCutGenerator ()
  : model_(NULL),
    generator_(NULL),
    whenCutGenerator_(-1),
    whenCutGeneratorInSub_(-100),
    depthCutGenerator_(-1),
    depthCutGeneratorInSub_(-1),
    generatorName_(NULL),
    normal_(true),
    atSolution_(false),
    whenInfeasible_(false),
    timing_(false),
    timeInCutGenerator_(0.0),
    numberTimes_(0),
    numberCuts_(0),
    numberCutsActive_(0)
{
}
// Normal constructor
CbcCutGenerator::CbcCutGenerator(CbcModel * model,CglCutGenerator * generator,
				 int howOften, const char * name,
				 bool normal, bool atSolution, 
				 bool infeasible, int howOftenInSub,
				 int whatDepth, int whatDepthInSub)
  : 
    depthCutGenerator_(whatDepth),
    depthCutGeneratorInSub_(whatDepthInSub),
    timing_(false),
    timeInCutGenerator_(0.0),
    numberTimes_(0),
    numberCuts_(0),
    numberCutsActive_(0)
{
  model_ = model;
  generator_=generator->clone();
  generator_->refreshSolver(model_->solver());
  whenCutGenerator_=howOften;
  whenCutGeneratorInSub_ = howOftenInSub;
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
  generator_->refreshSolver(model_->solver());
  whenCutGenerator_=rhs.whenCutGenerator_;
  whenCutGeneratorInSub_ = rhs.whenCutGeneratorInSub_;
  depthCutGenerator_=rhs.depthCutGenerator_;
  depthCutGeneratorInSub_ = rhs.depthCutGeneratorInSub_;
  generatorName_=strdup(rhs.generatorName_);
  normal_=rhs.normal_;
  atSolution_=rhs.atSolution_;
  whenInfeasible_=rhs.whenInfeasible_;
  timing_ = rhs.timing_;
  timeInCutGenerator_ = rhs.timeInCutGenerator_;
  numberTimes_ = rhs.numberTimes_;
  numberCuts_ = rhs.numberCuts_;
  numberCutsActive_ = rhs.numberCutsActive_;
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
    depthCutGenerator_=rhs.depthCutGenerator_;
    depthCutGeneratorInSub_ = rhs.depthCutGeneratorInSub_;
    generatorName_=strdup(rhs.generatorName_);
    normal_=rhs.normal_;
    atSolution_=rhs.atSolution_;
    whenInfeasible_=rhs.whenInfeasible_;
    timing_ = rhs.timing_;
    timeInCutGenerator_ = rhs.timeInCutGenerator_;
    numberTimes_ = rhs.numberTimes_;
    numberCuts_ = rhs.numberCuts_;
    numberCutsActive_ = rhs.numberCutsActive_;
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
CbcCutGenerator::generateCuts( OsiCuts & cs , bool fullScan, CbcNode * node)
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
  OsiSolverInterface * solver = model_->solver();
  int depth;
  if (node)
    depth=node->depth();
  else
    depth=0;
  int pass=model_->getCurrentPassNumber()-1;
  bool doThis=(model_->getNodeCount()%howOften)==0;
  if (depthCutGenerator_>0) 
    doThis = (depth % depthCutGenerator_) ==0;
  // But turn off if 100
  if (howOften==100)
    doThis=false;
  if (fullScan||doThis) {
    double time1=0.0;
    if (timing_)
      time1 = CoinCpuTime();
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
      //OsiSolverInterface * solver = model_->solver();
      //void * saveData = solver->getApplicationData();
      //solver->setApplicationData(model_);
      generator_->generateCuts(*solver,cs,info);
      //solver->setApplicationData(saveData);
    } else {
      // Probing - return tight column bounds
      generator->generateCutsAndModify(*solver,cs,info);
      const double * tightLower = generator->tightLower();
      const double * lower = solver->getColLower();
      const double * tightUpper = generator->tightUpper();
      const double * upper = solver->getColUpper();
      const double * solution = solver->getColSolution();
      int j;
      int numberColumns = solver->getNumCols();
      double primalTolerance = 1.0e-8;
      for (j=0;j<numberColumns;j++) {
	if (tightUpper[j]==tightLower[j]&&
	    upper[j]>lower[j]) {
	  // fix
	  solver->setColLower(j,tightLower[j]);
	  solver->setColUpper(j,tightUpper[j]);
	  if (tightLower[j]>solution[j]+primalTolerance||
	      tightUpper[j]<solution[j]-primalTolerance)
	    returnCode=true;
	}
      }
    }
    if (timing_)
      timeInCutGenerator_ += CoinCpuTime()-time1;
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
