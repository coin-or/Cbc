// Copyright (C) 2004, International Business Machines
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
#include "CbcHeuristicFPump.hpp"
#include "CbcBranchActual.hpp"
#include "CoinHelperFunctions.hpp"
#include "CoinTime.hpp"


// Default Constructor
CbcHeuristicFPump::CbcHeuristicFPump() 
  :CbcHeuristic(),
   startTime_(0.0),
   maximumTime_(0.0),
   maximumPasses_(100),
   downValue_(0.5),
   roundExpensive_(false)
{
  setWhen(1);
}

// Constructor from model
CbcHeuristicFPump::CbcHeuristicFPump(CbcModel & model,
				     double downValue,bool roundExpensive)
  :CbcHeuristic(model),
   startTime_(0.0),
   maximumTime_(0.0),
   maximumPasses_(100),
   downValue_(downValue),
   roundExpensive_(roundExpensive)
{
  setWhen(1);
}

// Destructor 
CbcHeuristicFPump::~CbcHeuristicFPump ()
{
}

// Clone
CbcHeuristic *
CbcHeuristicFPump::clone() const
{
  return new CbcHeuristicFPump(*this);
}

// Copy constructor 
CbcHeuristicFPump::CbcHeuristicFPump(const CbcHeuristicFPump & rhs)
:
  CbcHeuristic(rhs),
  startTime_(rhs.startTime_),
  maximumTime_(rhs.maximumTime_),
  maximumPasses_(rhs.maximumPasses_),
  downValue_(rhs.downValue_),
  roundExpensive_(rhs.roundExpensive_)
{
  setWhen(rhs.when());
}
// Resets stuff if model changes
void 
CbcHeuristicFPump::resetModel(CbcModel * model)
{
}

/**************************BEGIN MAIN PROCEDURE ***********************************/

// See if feasibility pump will give better solution
// Sets value of solution
// Returns 1 if solution, 0 if not
int
CbcHeuristicFPump::solution(double & solutionValue,
			 double * betterSolution)
{
  if (!when()||(when()==1&&model_->phase()!=1))
    return 0; // switched off
  // See if at root node
  bool atRoot = model_->getNodeCount()==0;
  int passNumber = model_->getCurrentPassNumber();
  // just do once
  if (!atRoot||passNumber!=1)
    return 0;
  // probably a good idea
  if (model_->getSolutionCount()) return 0;
  // Clone solver - otherwise annoys root node computations
  OsiSolverInterface * solver = model_->solver()->clone();
  solver->resolve();
  const double * lower = solver->getColLower();
  const double * upper = solver->getColUpper();
  const double * solution = solver->getColSolution();
  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance,primalTolerance);
  
  int numberColumns = model_->getNumCols();
  int numberIntegers = model_->numberIntegers();
  const int * integerVariable = model_->integerVariable();

// 1. initially check 0-1
  int i,j;
  bool zeroOne=true;
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integerVariable[i];
    const CbcObject * object = model_->object(i);
    const CbcSimpleInteger * integerObject = 
      dynamic_cast<const  CbcSimpleInteger *> (object);
    assert(integerObject);
    if (upper[iColumn]-lower[iColumn]>1.000001) {
      zeroOne=false;
      break;
    }
  }
  if (!zeroOne) {
    delete solver;
    return 0;
  }

// 2. space for rounded solution
  double * newSolution = new double [numberColumns];
  // space for last rounded solutions
#define NUMBER_OLD 4
  double ** oldSolution = new double * [NUMBER_OLD];
  for (j=0;j<NUMBER_OLD;j++) {
    oldSolution[j]= new double[numberColumns];
    for (i=0;i<numberColumns;i++) oldSolution[j][i]=-COIN_DBL_MAX;
  }

// 3. Replace objective with an initial 0-valued objective
  double * saveObjective = new double [numberColumns];
  memcpy(saveObjective,solver->getObjCoefficients(),numberColumns*sizeof(double));
  for (i=0;i<numberColumns;i++) {
    solver->setObjCoeff(i,0.0);
  }
  bool finished=false;
  double direction = solver->getObjSense();
  int returnCode=0;
  bool takeHint;
  OsiHintStrength strength;
  solver->getHintParam(OsiDoDualInResolve,takeHint,strength);
  solver->setHintParam(OsiDoDualInResolve,false);
  solver->messageHandler()->setLogLevel(1);

// 4. Save objective offset so we can see progress
  double saveOffset;
  solver->getDblParam(OsiObjOffset,saveOffset);

// 5. MAIN WHILE LOOP
  int numberPasses=0;
  while (!finished) {
    if (numberPasses>=maximumPasses_) {
      break;
    }
    if (maximumTime_>0.0&&CoinCpuTime()>=startTime_+maximumTime_) break;
    numberPasses++;
    memcpy(newSolution,solution,numberColumns*sizeof(double));
    int flip;
    returnCode = rounds(newSolution,saveObjective,roundExpensive_,downValue_,&flip);
    if (returnCode) {
      // SOLUTION IS INTEGER
      // Put back correct objective
      printf("\n");
      for (i=0;i<numberColumns;i++)
        solver->setObjCoeff(i,saveObjective[i]);
      // solution - but may not be better
      // Compute using dot product
      double newSolutionValue = direction*solver->OsiSolverInterface::getObjValue();
      if (newSolutionValue<solutionValue) {
	memcpy(betterSolution,newSolution,numberColumns*sizeof(double));
      } else {
	returnCode=0;
      }      
      break;
    } else {
      // SOLUTION IS not INTEGER
      // 1. check for loop
      bool matched;
      for (int k = NUMBER_OLD-1; k > 0; k--) {
  	  double * b = oldSolution[k];
          matched = true;
          for (i = 0; i <numberIntegers; i++) {
	      int iColumn = integerVariable[i];
	      if (newSolution[iColumn]!=b[iColumn]) {
		matched=false;
		break;
	      }
	  }
	  if (matched) break;
      }
      if (matched || numberPasses%100 == 0) {
	 // perturbation
	 printf("Perturbation applied");
	 for (i=0;i<numberIntegers;i++) {
	     int iColumn = integerVariable[i];
	     double value = max(0.0,CoinDrand48()-0.3);
	     double difference = fabs(solution[iColumn]-newSolution[iColumn]);
	     if (difference+value>0.5) {
	        if (newSolution[iColumn]<lower[iColumn]+primalTolerance) newSolution[iColumn] += 1.0;
	     else if (newSolution[iColumn]>upper[iColumn]-primalTolerance) newSolution[iColumn] -= 1.0;
	          else abort();
	     }
	 }
      } else {
         for (j=NUMBER_OLD-1;j>0;j--) {
             for (i = 0; i < numberColumns; i++) oldSolution[j][i]=oldSolution[j-1][i];
	 }
         for (j = 0; j < numberColumns; j++) oldSolution[0][j] = newSolution[j];
      }

      // 2. update the objective function based on the new rounded solution
      double offset=0.0;
      for (i=0;i<numberIntegers;i++) {
	int iColumn = integerVariable[i];
	double costValue = 1.0;
	// deal with fixed variables (i.e., upper=lower)
	if (fabs(lower[iColumn]-upper[iColumn]) < primalTolerance) {
	   if (lower[iColumn] > 1. - primalTolerance) solver->setObjCoeff(iColumn,-costValue);
	   else                                       solver->setObjCoeff(iColumn,costValue);
	   continue;
	}
	if (newSolution[iColumn]<lower[iColumn]+primalTolerance) {
	  solver->setObjCoeff(iColumn,costValue);
	} else {
          if (newSolution[iColumn]>upper[iColumn]-primalTolerance) {
	    solver->setObjCoeff(iColumn,-costValue);
	  } else {
	    abort();
          }
	}
	offset += costValue*newSolution[iColumn];
      }
      solver->setDblParam(OsiObjOffset,-offset);
      solver->resolve();
      printf("\npass %3d: obj. %10.5lf --> ", numberPasses,solver->getObjValue());


    }
  } // END WHILE

  delete solver;
  delete [] newSolution;
  for ( j=0;j<NUMBER_OLD;j++) 
    delete [] oldSolution[j];
  delete [] oldSolution;
  delete [] saveObjective;
  return returnCode;
}

/**************************END MAIN PROCEDURE ***********************************/

// update model
void CbcHeuristicFPump::setModel(CbcModel * model)
{
  model_ = model;
}

/* Rounds solution - down if < downValue
   returns 1 if current is a feasible solution
*/
int 
CbcHeuristicFPump::rounds(double * solution,
			  const double * objective,
			  bool roundExpensive, double downValue, int *flip)
{
  OsiSolverInterface * solver = model_->solver();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double primalTolerance ;
  solver->getDblParam(OsiPrimalTolerance,primalTolerance) ;
  int numberIntegers = model_->numberIntegers();
  const int * integerVariable = model_->integerVariable();

  int i;

  static int iter = 0;
  int numberColumns = model_->getNumCols();
  // tmp contains the current obj coefficients 
  double * tmp = new double [numberColumns];
  memcpy(tmp,solver->getObjCoefficients(),numberColumns*sizeof(double));
  int flip_up = 0;
  int flip_down  = 0;
  double  v = CoinDrand48() * 20;
  int nn = 10 + (int) v;
  int nnv = 0;
  int * list = new int [nn];
  double * val = new double [nn];
  for (i = 0; i < nn; i++) val[i] = .001;

  // return rounded solution
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integerVariable[i];
    const CbcObject * object = model_->object(i);
    const CbcSimpleInteger * integerObject = 
      dynamic_cast<const  CbcSimpleInteger *> (object);
    assert(integerObject);
    double value=solution[iColumn];
    double round = floor(value+primalTolerance);
    if (value-round > .5) round += 1.;
    if (round < integerTolerance && tmp[iColumn] < -1. + integerTolerance) flip_down++;
    if (round > 1. - integerTolerance && tmp[iColumn] > 1. - integerTolerance) flip_up++;
    if (flip_up + flip_down == 0) { 
       for (int k = 0; k < nn; k++) {
           if (fabs(value-round) > val[k]) {
              nnv++;
              for (int j = nn-2; j >= k; j--) {
                  val[j+1] = val[j];
                  list[j+1] = list[j];
              } 
              val[k] = fabs(value-round);
              list[k] = iColumn;
              break;
           }
       }
    }
    solution[iColumn] = round;
  }

  if (nnv > nn) nnv = nn;
  if (iter != 0) printf("up = %5d , down = %5d", flip_up, flip_down); fflush(stdout);
  *flip = flip_up + flip_down;
  delete [] tmp;

  if (*flip == 0 && iter != 0) {
     printf(" -- rand = %4d (%4d) ", nnv, nn);
     for (i = 0; i < nnv; i++) solution[list[i]] = 1. - solution[list[i]];
     *flip = nnv;
  } else printf(" ");
  delete [] list; delete [] val;
  iter++;
    
  const double * rowLower = solver->getRowLower();
  const double * rowUpper = solver->getRowUpper();

  int numberRows = solver->getNumRows();
  // get row activities
  double * rowActivity = new double[numberRows];
  memset(rowActivity,0,numberRows*sizeof(double));
  solver->getMatrixByCol()->times(solution,rowActivity) ;
  double largestInfeasibility =0.0;
  for (i=0 ; i < numberRows ; i++) {
    largestInfeasibility = max(largestInfeasibility,
			       rowLower[i]-rowActivity[i]);
    largestInfeasibility = max(largestInfeasibility,
			       rowActivity[i]-rowUpper[i]);
  }
  return (largestInfeasibility>primalTolerance) ? 0 : 1;
}
// Set maximum Time (default off) - also sets starttime to current
void 
CbcHeuristicFPump::setMaximumTime(double value)
{
  startTime_=CoinCpuTime();
  maximumTime_=value;
}

  
