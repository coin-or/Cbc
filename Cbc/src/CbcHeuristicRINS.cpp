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
  numCouldRun_++;
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
  if ((model_->getNodeCount()%howOften_)==0&&(model_->getCurrentPassNumber()==1||model_->getCurrentPassNumber()==999999)) {
    OsiSolverInterface * solver = model_->solver();

    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();
  
    const double * currentSolution = solver->getColSolution();
    OsiSolverInterface * newSolver = model_->continuousSolver()->clone();
    //const double * colLower = newSolver->getColLower();
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
	newSolver->setColUpper(iColumn,nearest);
	nFix++;
      }
    }
    if (nFix>numberIntegers/5) {
      //printf("%d integers have same value\n",nFix);
      returnCode = smallBranchAndBound(newSolver,numberNodes_,betterSolution,solutionValue,
                                         model_->getCutoff(),"CbcHeuristicRINS");
      if (returnCode<0)
	returnCode=0; // returned on size
      else
	numRuns_++;
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
	howOften_ += static_cast<int> (howOften_*decayFactor_);
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
  int numberAtBound=0;

  for (i=0;i<numberIntegers;i++) {
    int iColumn=integerVariable[i];
    double value = currentSolution[iColumn];
    double lower = colLower[iColumn];
    double upper = colUpper[iColumn];
    value = CoinMax(value,lower);
    value = CoinMin(value,upper);
    if (fabs(value-floor(value+0.5))<1.0e-8) {
      value = floor(value+0.5);
      if (value==lower||value==upper)
	numberAtBound++;
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
#ifdef COIN_DEVELOP
      if (!numberTightened&&numberFixed==numberAtBound)
	printf("could add cut with %d elements\n",numberFixed);
#endif
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

// Default Constructor
CbcHeuristicDINS::CbcHeuristicDINS() 
  :CbcHeuristic()
{
  numberSolutions_=0;
  numberSuccesses_=0;
  numberTries_=0;
  howOften_=100;
  decayFactor_ = 0.5; 
  maximumKeepSolutions_=5;
  numberKeptSolutions_=0;
  numberIntegers_=-1;
  localSpace_ = 10;
  values_=NULL;
}

// Constructor with model - assumed before cuts

CbcHeuristicDINS::CbcHeuristicDINS(CbcModel & model)
  :CbcHeuristic(model)
{
  numberSolutions_=0;
  numberSuccesses_=0;
  numberTries_=0;
  howOften_=100;
  decayFactor_ = 0.5; 
  assert(model.solver());
  maximumKeepSolutions_=5;
  numberKeptSolutions_=0;
  numberIntegers_=-1;
  localSpace_ = 10;
  values_=NULL;
}

// Destructor 
CbcHeuristicDINS::~CbcHeuristicDINS ()
{
  for (int i=0;i<numberKeptSolutions_;i++)
    delete [] values_[i];
  delete [] values_;
}

// Clone
CbcHeuristic *
CbcHeuristicDINS::clone() const
{
  return new CbcHeuristicDINS(*this);
}

// Assignment operator 
CbcHeuristicDINS & 
CbcHeuristicDINS::operator=( const CbcHeuristicDINS& rhs)
{
  if (this!=&rhs) {
    CbcHeuristic::operator=(rhs);
    numberSolutions_ = rhs.numberSolutions_;
    howOften_ = rhs.howOften_;
    decayFactor_ = rhs.decayFactor_;
    numberSuccesses_ = rhs.numberSuccesses_;
    numberTries_ = rhs.numberTries_;
    for (int i=0;i<numberKeptSolutions_;i++)
      delete [] values_[i];
    delete [] values_;
    maximumKeepSolutions_ = rhs.maximumKeepSolutions_;
    numberKeptSolutions_ = rhs.numberKeptSolutions_;
    numberIntegers_ = rhs.numberIntegers_;
    localSpace_ = rhs.localSpace_;
    if (model_&&rhs.values_) {
      assert (numberIntegers_>=0);
      values_ = new int * [maximumKeepSolutions_];
      for (int i=0;i<maximumKeepSolutions_;i++)
	values_[i]=CoinCopyOfArray(rhs.values_[i],numberIntegers_);
    } else {
      values_=NULL;
    }
  }
  return *this;
}

// Create C++ lines to get to current state
void 
CbcHeuristicDINS::generateCpp( FILE * fp) 
{
  CbcHeuristicDINS other;
  fprintf(fp,"0#include \"CbcHeuristicDINS.hpp\"\n");
  fprintf(fp,"3  CbcHeuristicDINS heuristicDINS(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp,"heuristicDINS");
  if (howOften_!=other.howOften_)
    fprintf(fp,"3  heuristicDINS.setHowOften(%d);\n",howOften_);
  else
    fprintf(fp,"4  heuristicDINS.setHowOften(%d);\n",howOften_);
  fprintf(fp,"3  cbcModel->addHeuristic(&heuristicDINS);\n");
}

// Copy constructor 
CbcHeuristicDINS::CbcHeuristicDINS(const CbcHeuristicDINS & rhs)
:
  CbcHeuristic(rhs),
  numberSolutions_(rhs.numberSolutions_),
  howOften_(rhs.howOften_),
  decayFactor_(rhs.decayFactor_),
  numberSuccesses_(rhs.numberSuccesses_),
  numberTries_(rhs.numberTries_),
  maximumKeepSolutions_(rhs.maximumKeepSolutions_),
  numberKeptSolutions_(rhs.numberKeptSolutions_),
  numberIntegers_(rhs.numberIntegers_),
  localSpace_(rhs.localSpace_)
{
  if (model_&&rhs.values_) {
    assert (numberIntegers_>=0);
    values_ = new int * [maximumKeepSolutions_];
    for (int i=0;i<maximumKeepSolutions_;i++)
      values_[i]=CoinCopyOfArray(rhs.values_[i],numberIntegers_);
  } else {
    values_=NULL;
  }
}
// Resets stuff if model changes
void 
CbcHeuristicDINS::resetModel(CbcModel * model)
{
  //CbcHeuristic::resetModel(model);
  for (int i=0;i<numberKeptSolutions_;i++)
    delete [] values_[i];
  delete [] values_;
  numberKeptSolutions_=0;
  numberIntegers_=-1;
  numberSolutions_=0;
  values_=NULL;
}
/*
  First tries setting a variable to better value.  If feasible then
  tries setting others.  If not feasible then tries swaps
  Returns 1 if solution, 0 if not */
int
CbcHeuristicDINS::solution(double & solutionValue,
			 double * betterSolution)
{
  numCouldRun_++;
  int returnCode=0;
  const double * bestSolution = model_->bestSolution();
  if (!bestSolution)
    return 0; // No solution found yet
  if (numberSolutions_<model_->getSolutionCount()) {
    // new solution - add info
    numberSolutions_=model_->getSolutionCount();

    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();
    if (numberIntegers_<0) {
      numberIntegers_ = numberIntegers;
      assert (!values_);
      values_ = new int * [maximumKeepSolutions_];
      for (int i=0;i<maximumKeepSolutions_;i++)
	values_[i]=NULL;
    } else {
      assert (numberIntegers==numberIntegers_);
    }
    // move solutions (0 will be most recent)
    {
      int * temp = values_[maximumKeepSolutions_-1];
      for (int i = maximumKeepSolutions_-1;i>0;i--)
	values_[i]=values_[i-1];
      if (!temp)
	temp = new int [numberIntegers_];
      values_[0]=temp;
    }
    int i;
    for (i=0;i<numberIntegers;i++) {
      int iColumn = integerVariable[i];
      double value=bestSolution[iColumn];
      double nearest=floor(value+0.5);
      values_[0][i]=static_cast<int> (nearest);
    }
    numberKeptSolutions_ = CoinMin(numberKeptSolutions_+1,maximumKeepSolutions_);
  } 
  int finalReturnCode=0;
  if (((model_->getNodeCount()%howOften_)==howOften_/2||!model_->getNodeCount())&&(model_->getCurrentPassNumber()==1||model_->getCurrentPassNumber()==999999)) {
    OsiSolverInterface * solver = model_->solver();

    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();
  
    const double * currentSolution = solver->getColSolution();
    int localSpace = localSpace_;
    // 0 means finished but no solution, 1 solution, 2 node limit
    int status=-1;
    double cutoff = model_->getCutoff();
    while(status) {
      status=0;
      OsiSolverInterface * newSolver = model_->continuousSolver()->clone();
      const double * colLower = solver->getColLower();
      const double * colUpper = solver->getColUpper();
      
      double primalTolerance;
      solver->getDblParam(OsiPrimalTolerance,primalTolerance);
      const double * continuousSolution = newSolver->getColSolution();
      // Space for added constraint
      double * element = new double [numberIntegers];
      int * column = new int [numberIntegers];
      int i;
      int nFix=0;
      int nCouldFix=0;
      int nCouldFix2=0;
      int nBound=0;
      int nEl=0;
      double bias=localSpace;
      int okSame = numberKeptSolutions_-1;
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
	int intValue = static_cast<int> (floor(valueInt+0.5));
	double currentValue = currentSolution[iColumn];
	double currentLower = colLower[iColumn];
	double currentUpper = colUpper[iColumn];
	if (fabs(valueInt-currentValue)>=0.5) {
	  // Re-bound
	  nBound++;
	  if (intValue>=currentValue) {
	    currentLower = CoinMax(currentLower,ceil(2*currentValue-intValue));
	    currentUpper = intValue;
	  } else {
	    currentLower = intValue;
	    currentUpper = CoinMin(currentUpper,floor(2*currentValue-intValue));
	  }
	  newSolver->setColLower(iColumn,currentLower);
	  newSolver->setColUpper(iColumn,currentUpper);
	} else {
	  // See if can fix
	  bool canFix=false;
	  double continuousValue = continuousSolution[iColumn];
	  if (fabs(currentValue-valueInt)<10.0*primalTolerance) {
	    if (currentUpper-currentLower>1.0) {
	      // General integer variable
	      canFix=true;
	    } else if(fabs(continuousValue-valueInt)<10.0*primalTolerance) {
	      int nSame=1;
	      //assert (intValue==values_[0][i]);
	      for (int k=1;k<numberKeptSolutions_;k++) {
		if (intValue==values_[k][i])
		  nSame++;
	      }
	      if (nSame>=okSame) {
		// can fix
		canFix=true;
	      } else {
		nCouldFix++;
	      }
	    } else {
	      nCouldFix2++;
	    }
	  }
	  if (canFix) {
	    newSolver->setColLower(iColumn,intValue);
	    newSolver->setColUpper(iColumn,intValue);
	    nFix++;
	  } else {
	    if (currentUpper-currentLower>1.0) {
	      // General integer variable
	      currentLower = floor(currentValue);
	      if (intValue>=currentLower&&intValue<=currentLower+1) {
		newSolver->setColLower(iColumn,currentLower);
		newSolver->setColUpper(iColumn,currentLower+1.0);
	      } else {
		// fix
		double value;
		if (intValue<currentLower)
		  value = currentLower;
		else
		  value = currentLower+1;
		newSolver->setColLower(iColumn,value);
		newSolver->setColUpper(iColumn,value);
		nFix++;
	      }
	    } else {
	      // 0-1 (ish)
	      column[nEl]=iColumn;
	      if (intValue==currentLower) {
		bias += currentLower;
		element[nEl++]=1.0;
	      } else if (intValue==currentUpper) {
		bias += currentUpper;
		element[nEl++]=-1.0;
	      } else {
		printf("bad DINS logic\n");
		abort();
	      }
	    }
	  }
	}
      }
      char generalPrint[200];
      sprintf(generalPrint,
	      "%d fixed, %d same as cont/int, %d same as int - %d bounded %d in cut\n",
	      nFix,nCouldFix,nCouldFix2,nBound,nEl);
      model_->messageHandler()->message(CBC_FPUMP2,model_->messages())
	<< generalPrint
	<<CoinMessageEol;
      if (nFix>numberIntegers/10) {
#if 0
	newSolver->initialSolve();
	printf("obj %g\n",newSolver->getObjValue());
	for (i=0;i<numberIntegers;i++) {
	  int iColumn=integerVariable[i];
	  printf("%d new bounds %g %g - solutions %g %g\n",
		 iColumn,newSolver->getColLower()[iColumn],
		 newSolver->getColUpper()[iColumn],
		 bestSolution[iColumn],
		 currentSolution[iColumn]);
	}
#endif
	if (nEl>0)
	  newSolver->addRow(nEl,column,element,-COIN_DBL_MAX,bias);
	//printf("%d integers have same value\n",nFix);
	returnCode = smallBranchAndBound(newSolver,numberNodes_,betterSolution,solutionValue,
                                         cutoff,"CbcHeuristicDINS");
	if (returnCode<0) {
	  returnCode=0; // returned on size
	  status=0;
	} else {
	  numRuns_++;
	  if ((returnCode&2)!=0) {
	    // could add cut as complete search
	    returnCode &= ~2;
	    if ((returnCode&1)!=0) {
	      numberSuccesses_++;
	      status=1;
	    } else {
	      // no solution
	      status=0;
	    }
	  } else {
	    if ((returnCode&1)!=0) {
	      numberSuccesses_++;
	      status=1;
	    } else {
	      // no solution but node limit
	      status=2;
	      if (nEl)
		localSpace -= 5;
	      else
		localSpace = -1;
	      if (localSpace<0)
		status=0;
	    }
	  }
	  if ((returnCode&1)!=0) {
	    cutoff = CoinMin(cutoff,solutionValue-model_->getCutoffIncrement());
	    finalReturnCode=1;
	  }
	}
      }
      delete [] element;
      delete [] column;
      delete newSolver;
    }
    numberTries_++;
    if ((numberTries_%10)==0&&numberSuccesses_*3<numberTries_)
      howOften_ += static_cast<int> (howOften_*decayFactor_);
  }
  return finalReturnCode;
}
// update model
void CbcHeuristicDINS::setModel(CbcModel * model)
{
  model_ = model;
  // Get a copy of original matrix
  assert(model_->solver());
  for (int i=0;i<numberKeptSolutions_;i++)
    delete [] values_[i];
  delete [] values_;
  numberKeptSolutions_=0;
  numberIntegers_=-1;
  numberSolutions_=0;
  values_=NULL;
}

  
