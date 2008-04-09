// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cfloat>

#include "OsiSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcHeuristicRINS.hpp"
#include "CbcBranchActual.hpp"
#include "CbcStrategy.hpp"
#include "CglPreProcess.hpp"

// Default Constructor
CbcHeuristicRINS::CbcHeuristicRINS() 
  :CbcHeuristic()
{
  numberSolutions_=0;
  numberSuccesses_=0;
  numberTries_=0;
  howOften_=100;
  decayFactor_ = 0.5; 
  used_=NULL;
}

// Constructor with model - assumed before cuts

CbcHeuristicRINS::CbcHeuristicRINS(CbcModel & model)
  :CbcHeuristic(model)
{
  numberSolutions_=0;
  numberSuccesses_=0;
  numberTries_=0;
  howOften_=100;
  decayFactor_ = 0.5; 
  assert(model.solver());
  int numberColumns = model.solver()->getNumCols();
  used_ = new char[numberColumns];
  memset(used_,0,numberColumns);
}

// Destructor 
CbcHeuristicRINS::~CbcHeuristicRINS ()
{
  delete [] used_;
}

// Clone
CbcHeuristic *
CbcHeuristicRINS::clone() const
{
  return new CbcHeuristicRINS(*this);
}

// Assignment operator 
CbcHeuristicRINS & 
CbcHeuristicRINS::operator=( const CbcHeuristicRINS& rhs)
{
  if (this!=&rhs) {
    CbcHeuristic::operator=(rhs);
    numberSolutions_ = rhs.numberSolutions_;
    howOften_ = rhs.howOften_;
    decayFactor_ = rhs.decayFactor_;
    numberSuccesses_ = rhs.numberSuccesses_;
    numberTries_ = rhs.numberTries_;
    delete [] used_;
    if (model_&&rhs.used_) {
      int numberColumns = model_->solver()->getNumCols();
      used_ = new char[numberColumns];
      memcpy(used_,rhs.used_,numberColumns);
    } else {
      used_=NULL;
    }
  }
  return *this;
}

// Create C++ lines to get to current state
void 
CbcHeuristicRINS::generateCpp( FILE * fp) 
{
  CbcHeuristicRINS other;
  fprintf(fp,"0#include \"CbcHeuristicRINS.hpp\"\n");
  fprintf(fp,"3  CbcHeuristicRINS heuristicRINS(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp,"heuristicRINS");
  if (howOften_!=other.howOften_)
    fprintf(fp,"3  heuristicRINS.setHowOften(%d);\n",howOften_);
  else
    fprintf(fp,"4  heuristicRINS.setHowOften(%d);\n",howOften_);
  fprintf(fp,"3  cbcModel->addHeuristic(&heuristicRINS);\n");
}

// Copy constructor 
CbcHeuristicRINS::CbcHeuristicRINS(const CbcHeuristicRINS & rhs)
:
  CbcHeuristic(rhs),
  numberSolutions_(rhs.numberSolutions_),
  howOften_(rhs.howOften_),
  decayFactor_(rhs.decayFactor_),
  numberSuccesses_(rhs.numberSuccesses_),
  numberTries_(rhs.numberTries_)
{
  if (model_&&rhs.used_) {
    int numberColumns = model_->solver()->getNumCols();
    used_ = new char[numberColumns];
    memcpy(used_,rhs.used_,numberColumns);
  } else {
    used_=NULL;
  }
}
// Resets stuff if model changes
void 
CbcHeuristicRINS::resetModel(CbcModel * model)
{
  //CbcHeuristic::resetModel(model);
  delete [] used_;
  if (model_&&used_) {
    int numberColumns = model_->solver()->getNumCols();
    used_ = new char[numberColumns];
    memset(used_,0,numberColumns);
  } else {
    used_=NULL;
  }
}
/*
  First tries setting a variable to better value.  If feasible then
  tries setting others.  If not feasible then tries swaps
  Returns 1 if solution, 0 if not */
int
CbcHeuristicRINS::solution(double & solutionValue,
			 double * betterSolution)
{
  int returnCode=0;
  const double * bestSolution = model_->bestSolution();
  if (!bestSolution)
    return 0; // No solution found yet
  if (numberSolutions_<model_->getSolutionCount()) {
    // new solution - add info
    numberSolutions_=model_->getSolutionCount();

    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();
  
    int i;
    for (i=0;i<numberIntegers;i++) {
      int iColumn = integerVariable[i];
      const OsiObject * object = model_->object(i);
      // get original bounds
      double originalLower;
      double originalUpper;
      getIntegerInformation( object,originalLower, originalUpper); 
      double value=bestSolution[iColumn];
      if (value<originalLower) {
	value=originalLower;
      } else if (value>originalUpper) {
	value=originalUpper;
      }
      double nearest=floor(value+0.5);
      // if away from lower bound mark that fact
      if (nearest>originalLower) {
	used_[iColumn]=1;
      }
    }
  } 
  if ((model_->getNodeCount()%howOften_)==0&&model_->getCurrentPassNumber()==1) {
    OsiSolverInterface * solver = model_->solver();

    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();
  
    const double * currentSolution = solver->getColSolution();
    OsiSolverInterface * newSolver = model_->continuousSolver()->clone();
    const double * colLower = newSolver->getColLower();
    //const double * colUpper = newSolver->getColUpper();

    double primalTolerance;
    solver->getDblParam(OsiPrimalTolerance,primalTolerance);
    
    int i;
    int nFix=0;
    for (i=0;i<numberIntegers;i++) {
      int iColumn=integerVariable[i];
      const OsiObject * object = model_->object(i);
      // get original bounds
      double originalLower;
      double originalUpper;
      getIntegerInformation( object,originalLower, originalUpper); 
      double valueInt=bestSolution[iColumn];
      if (valueInt<originalLower) {
	valueInt=originalLower;
      } else if (valueInt>originalUpper) {
	valueInt=originalUpper;
      }
      if (fabs(currentSolution[iColumn]-valueInt)<10.0*primalTolerance) {
	double nearest=floor(valueInt+0.5);
	newSolver->setColLower(iColumn,nearest);
	newSolver->setColUpper(iColumn,colLower[iColumn]);
	nFix++;
      }
    }
    if (nFix>numberIntegers/5) {
      //printf("%d integers have same value\n",nFix);
      returnCode = smallBranchAndBound(newSolver,numberNodes_,betterSolution,solutionValue,
                                         model_->getCutoff(),"CbcHeuristicRINS");
      if (returnCode<0)
	returnCode=0; // returned on size
      if ((returnCode&1)!=0)
	numberSuccesses_++;
      //printf("return code %d",returnCode);
      if ((returnCode&2)!=0) {
	// could add cut
	returnCode &= ~2;
	//printf("could add cut with %d elements (if all 0-1)\n",nFix);
      } else {
	//printf("\n");
      }
      numberTries_++;
      if ((numberTries_%10)==0&&numberSuccesses_*3<numberTries_)
	howOften_ += (int) (howOften_*decayFactor_);
    }

    delete newSolver;
  }
  return returnCode;
}
// update model
void CbcHeuristicRINS::setModel(CbcModel * model)
{
  model_ = model;
  // Get a copy of original matrix
  assert(model_->solver());
  delete [] used_;
  int numberColumns = model->solver()->getNumCols();
  used_ = new char[numberColumns];
  memset(used_,0,numberColumns);
}
// Default Constructor
CbcHeuristicRENS::CbcHeuristicRENS() 
  :CbcHeuristic()
{
  numberTries_=0;
}

// Constructor with model - assumed before cuts

CbcHeuristicRENS::CbcHeuristicRENS(CbcModel & model)
  :CbcHeuristic(model)
{
  numberTries_=0;
}

// Destructor 
CbcHeuristicRENS::~CbcHeuristicRENS ()
{
}

// Clone
CbcHeuristic *
CbcHeuristicRENS::clone() const
{
  return new CbcHeuristicRENS(*this);
}

// Assignment operator 
CbcHeuristicRENS & 
CbcHeuristicRENS::operator=( const CbcHeuristicRENS& rhs)
{
  if (this!=&rhs) {
    CbcHeuristic::operator=(rhs);
    numberTries_ = rhs.numberTries_;
  }
  return *this;
}

// Copy constructor 
CbcHeuristicRENS::CbcHeuristicRENS(const CbcHeuristicRENS & rhs)
:
  CbcHeuristic(rhs),
  numberTries_(rhs.numberTries_)
{
}
// Resets stuff if model changes
void 
CbcHeuristicRENS::resetModel(CbcModel * model)
{
}
int
CbcHeuristicRENS::solution(double & solutionValue,
			 double * betterSolution)
{
  int returnCode=0;
  if (numberTries_)
    return 0; 
  numberTries_++;
  OsiSolverInterface * solver = model_->solver();
  
  int numberIntegers = model_->numberIntegers();
  const int * integerVariable = model_->integerVariable();
  
  const double * currentSolution = solver->getColSolution();
  OsiSolverInterface * newSolver = model_->continuousSolver()->clone();
  const double * colLower = newSolver->getColLower();
  const double * colUpper = newSolver->getColUpper();

  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance,primalTolerance);
    
  int i;
  int numberFixed=0;
  int numberTightened=0;

  for (i=0;i<numberIntegers;i++) {
    int iColumn=integerVariable[i];
    double value = currentSolution[iColumn];
    double lower = colLower[iColumn];
    double upper = colUpper[iColumn];
    value = CoinMax(value,lower);
    value = CoinMin(value,upper);
    if (fabs(value-floor(value+0.5))<1.0e-8) {
      value = floor(value+0.5);
      newSolver->setColLower(iColumn,value);
      newSolver->setColUpper(iColumn,value);
      numberFixed++;
    } else if (colUpper[iColumn]-colLower[iColumn]>=2.0) {
      numberTightened++;
      newSolver->setColLower(iColumn,floor(value));
      newSolver->setColUpper(iColumn,ceil(value));
    }
  }
  if (numberFixed>numberIntegers/5) {
#ifdef COIN_DEVELOP
    printf("%d integers fixed and %d tightened\n",numberFixed,numberTightened);
#endif
    returnCode = smallBranchAndBound(newSolver,numberNodes_,betterSolution,solutionValue,
				     model_->getCutoff(),"CbcHeuristicRENS");
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
  
  delete newSolver;
  return returnCode;
}
// update model
void CbcHeuristicRENS::setModel(CbcModel * model)
{
  model_ = model;
}

  
