// Copyright (C) 2008, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CbcHeuristicDive.hpp"
#include "CbcStrategy.hpp"
#include  "CoinTime.hpp"

// Default Constructor
CbcHeuristicDive::CbcHeuristicDive() 
  :CbcHeuristic()
{
  // matrix and row copy will automatically be empty
  downLocks_ =NULL;
  upLocks_ = NULL;
  percentageToFix_ = 0.2;
  maxIterations_ = 100;
  maxTime_ = 60;
}

// Constructor from model
CbcHeuristicDive::CbcHeuristicDive(CbcModel & model)
  :CbcHeuristic(model)
{
  downLocks_ =NULL;
  upLocks_ = NULL;
  // Get a copy of original matrix
  assert(model.solver());
  // model may have empty matrix - wait until setModel
  const CoinPackedMatrix * matrix = model.solver()->getMatrixByCol();
  if (matrix) {
    matrix_ = *matrix;
    validate();
  }
  percentageToFix_ = 0.2;
  maxIterations_ = 100;
  maxTime_ = 60;
}

// Destructor 
CbcHeuristicDive::~CbcHeuristicDive ()
{
  delete [] downLocks_;
  delete [] upLocks_;
}

// Create C++ lines to get to current state
void 
CbcHeuristicDive::generateCpp( FILE * fp, const char * heuristic) 
{
  // hard coded as CbcHeuristic virtual
  CbcHeuristic::generateCpp(fp,heuristic);
  if (percentageToFix_!=0.2)
    fprintf(fp,"3  %s.setPercentageToFix(%.f);\n",heuristic,percentageToFix_);
  else
    fprintf(fp,"4  %s.setPercentageToFix(%.f);\n",heuristic,percentageToFix_);
  if (maxIterations_!=100)
    fprintf(fp,"3  %s.setMaxIterations(%d);\n",heuristic,maxIterations_);
  else
    fprintf(fp,"4  %s.setMaxIterations(%d);\n",heuristic,maxIterations_);
  if (maxTime_!=60)
    fprintf(fp,"3  %s.setMaxTime(%.2f);\n",heuristic,maxTime_);
  else
    fprintf(fp,"4  %s.setMaxTime(%.2f);\n",heuristic,maxTime_);
}

// Copy constructor 
CbcHeuristicDive::CbcHeuristicDive(const CbcHeuristicDive & rhs)
:
  CbcHeuristic(rhs),
  matrix_(rhs.matrix_),
  percentageToFix_(rhs.percentageToFix_),
  maxIterations_(rhs.maxIterations_),
  maxTime_(rhs.maxTime_)
{
  if (rhs.downLocks_) {
    int numberIntegers = model_->numberIntegers();
    downLocks_ = CoinCopyOfArray(rhs.downLocks_,numberIntegers);
    upLocks_ = CoinCopyOfArray(rhs.upLocks_,numberIntegers);
  } else {
    downLocks_ = NULL;
    upLocks_ = NULL;
  }
}

// Assignment operator 
CbcHeuristicDive & 
CbcHeuristicDive::operator=( const CbcHeuristicDive& rhs)
{
  if (this!=&rhs) {
    CbcHeuristic::operator=(rhs);
    matrix_ = rhs.matrix_;
    percentageToFix_ = rhs.percentageToFix_;
    maxIterations_ = rhs.maxIterations_;
    maxTime_ = rhs.maxTime_;
    delete [] downLocks_;
    delete [] upLocks_;
    if (rhs.downLocks_) {
      int numberIntegers = model_->numberIntegers();
      downLocks_ = CoinCopyOfArray(rhs.downLocks_,numberIntegers);
      upLocks_ = CoinCopyOfArray(rhs.upLocks_,numberIntegers);
    } else {
      downLocks_ = NULL;
      upLocks_ = NULL;
    }
  }
  return *this;
}

// Resets stuff if model changes
void 
CbcHeuristicDive::resetModel(CbcModel * model)
{
  model_=model;
  assert(model_->solver());
  // Get a copy of original matrix
  const CoinPackedMatrix * matrix = model_->solver()->getMatrixByCol();
  // model may have empty matrix - wait until setModel
  if (matrix) {
    matrix_ = *matrix;
    validate();
  }
}

// update model
void CbcHeuristicDive::setModel(CbcModel * model)
{
  model_ = model;
  assert(model_->solver());
  // Get a copy of original matrix
  const CoinPackedMatrix * matrix = model_->solver()->getMatrixByCol();
  if (matrix) {
    matrix_ = *matrix;
    // make sure model okay for heuristic
    validate();
  }
}

bool CbcHeuristicDive::canHeuristicRun()
{
  return shouldHeurRun_randomChoice();
}


// See if dive fractional will give better solution
// Sets value of solution
// Returns 1 if solution, 0 if not
int
CbcHeuristicDive::solution(double & solutionValue,
			   double * betterSolution)
{
  ++numCouldRun_;

  // test if the heuristic can run
  if(!canHeuristicRun())
    return 0;

#if 0
  // See if to do
  if (!when()||(when()%10==1&&model_->phase()!=1)||
      (when()%10==2&&(model_->phase()!=2&&model_->phase()!=3)))
    return 0; // switched off
#endif

  double time1 = CoinCpuTime();

  OsiSolverInterface * solver = model_->solver()->clone();
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
  double direction = solver->getObjSense(); // 1 for min, -1 for max
  double newSolutionValue = direction*solver->getObjValue();
  int returnCode = 0;
  // Column copy
  const double * element = matrix_.getElements();
  const int * row = matrix_.getIndices();
  const CoinBigIndex * columnStart = matrix_.getVectorStarts();
  const int * columnLength = matrix_.getVectorLengths();

  // Get solution array for heuristic solution
  int numberColumns = solver->getNumCols();
  double * newSolution = new double [numberColumns];
  memcpy(newSolution,solution,numberColumns*sizeof(double));

  // vectors to store the latest variables fixed at their bounds
  int* columnFixed = new int [numberIntegers];
  double* originalBound = new double [numberIntegers];
  bool * fixedAtLowerBound = new bool [numberIntegers];

  const int maxNumberAtBoundToFix = (int) floor(percentageToFix_ * numberIntegers);

  // count how many fractional variables
  int numberFractionalVariables = 0;
  for (int i=0; i<numberIntegers; i++) {
    int iColumn = integerVariable[i];
    double value=newSolution[iColumn];
    if (fabs(floor(value+0.5)-value)>integerTolerance) {
      numberFractionalVariables++;
    }
  }

  int iteration = 0;
  while(numberFractionalVariables) {
    iteration++;

    // select a fractional variable to bound
    int bestColumn;
    int bestRound; // -1 rounds down, +1 rounds up
    selectVariableToBranch(solver, newSolution, bestColumn, bestRound);

    bool canRoundSolution = true;
    if(bestColumn != -1)
      canRoundSolution = false;
      
    // fixed integer variables that are at there bounds
    int numberAtBoundFixed = 0;
    for (int i=0; i<numberIntegers; i++) {
      int iColumn = integerVariable[i];
      double value=newSolution[iColumn];
      if(fabs(floor(value+0.5)-value)<=integerTolerance && 
	 numberAtBoundFixed < maxNumberAtBoundToFix) {
	// fix the variable at one of its bounds
	if (fabs(lower[iColumn]-value)<=integerTolerance &&
	    lower[iColumn] != upper[iColumn]) {
	  columnFixed[numberAtBoundFixed] = iColumn;
	  originalBound[numberAtBoundFixed] = upper[iColumn];
	  fixedAtLowerBound[numberAtBoundFixed] = true;
	  solver->setColUpper(iColumn, lower[iColumn]);
	  numberAtBoundFixed++;
	}
	else if(fabs(upper[iColumn]-value)<=integerTolerance &&
	    lower[iColumn] != upper[iColumn]) {
	  columnFixed[numberAtBoundFixed] = iColumn;
	  originalBound[numberAtBoundFixed] = lower[iColumn];
	  fixedAtLowerBound[numberAtBoundFixed] = false;
	  solver->setColLower(iColumn, upper[iColumn]);
	  numberAtBoundFixed++;
	}
	if(numberAtBoundFixed == maxNumberAtBoundToFix)
	  break;
      }
    }

    if(canRoundSolution) {
      // Round all the fractional variables
      for (int i=0; i<numberIntegers; i++) {
	int iColumn = integerVariable[i];
	double value=newSolution[iColumn];
	if (fabs(floor(value+0.5)-value)>integerTolerance) {
	  if(downLocks_[i]==0 || upLocks_[i]==0) {
	    if(downLocks_[i]==0 && upLocks_[i]==0) {
	      if(direction * objective[iColumn] >= 0.0)
		newSolution[iColumn] = floor(value);
	      else
		newSolution[iColumn] = ceil(value);
	    }
	    else if(downLocks_[i]==0)
	      newSolution[iColumn] = floor(value);
	    else
	      newSolution[iColumn] = ceil(value);
	  }
	  else
	    break;
	}
      }
      break;
    }


    double originalBoundBestColumn;
    if(bestColumn >= 0) {
      if(bestRound < 0) {
	originalBoundBestColumn = upper[bestColumn];
	solver->setColUpper(bestColumn, floor(newSolution[bestColumn]));
      }
      else {
	originalBoundBestColumn = lower[bestColumn];
	solver->setColLower(bestColumn, ceil(newSolution[bestColumn]));
      }
    } else {
      break;
    }
    int originalBestRound = bestRound;
    while (1) {

      solver->resolve();

      if(!solver->isProvenOptimal()) {
	if(numberAtBoundFixed > 0) {
	  // Remove the bound fix for variables that were at bounds
	  for(int i=0; i<numberAtBoundFixed; i++) {
	    int iColFixed = columnFixed[i];
	    if(fixedAtLowerBound[i])
	      solver->setColUpper(iColFixed, originalBound[i]);
	    else
	      solver->setColLower(iColFixed, originalBound[i]);
	  }
	  numberAtBoundFixed = 0;
	}
	else if(bestRound == originalBestRound) {
	  bestRound *= (-1);
	  if(bestRound < 0) {
	    solver->setColLower(bestColumn, originalBoundBestColumn);
	    solver->setColUpper(bestColumn, floor(newSolution[bestColumn]));
	  }
	  else {
	    solver->setColLower(bestColumn, ceil(newSolution[bestColumn]));
	    solver->setColUpper(bestColumn, originalBoundBestColumn);
	  }
	}
	else
	  break;
      }
      else
	break;
    }

    if(!solver->isProvenOptimal())
      break;

    if(iteration > maxIterations_) {
      break;
    }

    if(CoinCpuTime()-time1 > maxTime_) {
      break;
    }

    memcpy(newSolution,solution,numberColumns*sizeof(double));
    numberFractionalVariables = 0;
    for (int i=0; i<numberIntegers; i++) {
      int iColumn = integerVariable[i];
      double value=newSolution[iColumn];
      if (fabs(floor(value+0.5)-value)>integerTolerance) {
	numberFractionalVariables++;
      }
    }

  }


  double * rowActivity = new double[numberRows];
  memset(rowActivity,0,numberRows*sizeof(double));

  // re-compute new solution value
  double objOffset=0.0;
  solver->getDblParam(OsiObjOffset,objOffset);
  newSolutionValue = -objOffset;
  for (int i=0 ; i<numberColumns ; i++ )
    newSolutionValue += objective[i]*newSolution[i];
  newSolutionValue *= direction;
    //printf("new solution value %g %g\n",newSolutionValue,solutionValue);
  if (newSolutionValue<solutionValue) {
    // paranoid check
    memset(rowActivity,0,numberRows*sizeof(double));
    for (int i=0;i<numberColumns;i++) {
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
    for (int i=0;i<numberRows;i++) {
      if(rowActivity[i]<rowLower[i]) {
	if (rowActivity[i]<rowLower[i]-1000.0*primalTolerance)
	  feasible = false;
      } else if(rowActivity[i]>rowUpper[i]) {
	if (rowActivity[i]>rowUpper[i]+1000.0*primalTolerance)
	  feasible = false;
      }
    }
    for (int i=0; i<numberIntegers; i++) {
      int iColumn = integerVariable[i];
      double value=newSolution[iColumn];
      if (fabs(floor(value+0.5)-value)>integerTolerance) {
	feasible = false;
	break;
      }
    }
    if (feasible) {
      // new solution
      memcpy(betterSolution,newSolution,numberColumns*sizeof(double));
      solutionValue = newSolutionValue;
      //printf("** Solution of %g found by CbcHeuristicDive\n",newSolutionValue);
      returnCode=1;
    } else {
      // Can easily happen
      //printf("Debug CbcHeuristicDive giving bad solution\n");
    }
  }

  delete [] newSolution;
  delete [] columnFixed;
  delete [] originalBound;
  delete [] fixedAtLowerBound;
  delete [] rowActivity;
  delete solver;
  return returnCode;
}

// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
void 
CbcHeuristicDive::validate() 
{
  if (model_&&when()<10) {
    if (model_->numberIntegers()!=
        model_->numberObjects())
      setWhen(0);
  }

  int numberIntegers = model_->numberIntegers();
  const int * integerVariable = model_->integerVariable();
  delete [] downLocks_;
  delete [] upLocks_;
  downLocks_ = new unsigned short [numberIntegers];
  upLocks_ = new unsigned short [numberIntegers];
  // Column copy
  const double * element = matrix_.getElements();
  const int * row = matrix_.getIndices();
  const CoinBigIndex * columnStart = matrix_.getVectorStarts();
  const int * columnLength = matrix_.getVectorLengths();
  const double * rowLower = model_->solver()->getRowLower();
  const double * rowUpper = model_->solver()->getRowUpper();
  for (int i=0;i<numberIntegers;i++) {
    int iColumn = integerVariable[i];
    int down=0;
    int up=0;
    if (columnLength[iColumn]>65535) {
      setWhen(0);
      break; // unlikely to work
    }
    for (CoinBigIndex j=columnStart[iColumn];
	 j<columnStart[iColumn]+columnLength[iColumn];j++) {
      int iRow=row[j];
      if (rowLower[iRow]>-1.0e20&&rowUpper[iRow]<1.0e20) {
	up++;
	down++;
      } else if (element[j]>0.0) {
	if (rowUpper[iRow]<1.0e20)
	  up++;
	else
	  down++;
      } else {
	if (rowLower[iRow]>-1.0e20)
	  up++;
	else
	  down++;
      }
    }
    downLocks_[i] = (unsigned short) down;
    upLocks_[i] = (unsigned short) up;
  }
}
