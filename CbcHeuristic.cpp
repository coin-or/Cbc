// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <cassert>
#include <cmath>
#include <cfloat>

#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcHeuristic.hpp"
#include "CbcStrategy.hpp"
#include "CglPreProcess.hpp"

// Default Constructor
CbcHeuristic::CbcHeuristic() 
  :model_(NULL),
   when_(2)
{
}

// Constructor from model
CbcHeuristic::CbcHeuristic(CbcModel & model)
:
  model_(&model),
  when_(2)
{
}
// Resets stuff if model changes
void 
CbcHeuristic::resetModel(CbcModel * model)
{
  model_=model;
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
// Do mini branch and bound (return 1 if solution)
int 
CbcHeuristic::smallBranchAndBound(OsiSolverInterface * solver,int numberNodes,
                                  double * newSolution, double & newSolutionValue,
                                  double cutoff, std::string name) const
{
  OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (solver);
  if (osiclp) {
    // go faster stripes
    if (osiclp->getNumRows()<300&&osiclp->getNumCols()<500) {
      osiclp->setupForRepeatedUse(2,0);
    } else {
      osiclp->setupForRepeatedUse(0,0);
    }
    // Turn this off if you get problems
    // Used to be automatically set
    osiclp->setSpecialOptions(osiclp->specialOptions()|(128+64));
  }
  // Reduce printout
  solver->setHintParam(OsiDoReducePrint,true,OsiHintTry);
  solver->setHintParam(OsiDoPresolveInInitial,false,OsiHintTry);
  solver->setDblParam(OsiDualObjectiveLimit,cutoff);
  solver->initialSolve();
  int returnCode=1;
  int logLevel = model_->logLevel();
  if (solver->isProvenOptimal()) {
    CglPreProcess process;
    /* Do not try and produce equality cliques and
       do up to 5 passes */
    if (logLevel<=1)
      process.messageHandler()->setLogLevel(0);
    OsiSolverInterface * solver2= process.preProcess(*solver);
    if (!solver2) {
      if (logLevel>1)
        printf("Pre-processing says infeasible\n");
      returnCode=0; // so will be infeasible
    } else {
      solver2->resolve();
      CbcModel model(*solver2);
      if (logLevel<=1)
        model.setLogLevel(0);
      else
        model.setLogLevel(logLevel);
      model.setCutoff(cutoff);
      model.setMaximumNodes(numberNodes);
      model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
      // Lightweight
      CbcStrategyDefaultSubTree strategy(model_,true,5,1,0);
      model.setStrategy(strategy);
      model.solver()->setIntParam(OsiMaxNumIterationHotStart,10);
      // Do search
      if (logLevel>1)
        model_->messageHandler()->message(CBC_START_SUB,model_->messages())
          << name
          << model.getMaximumNodes()
          <<CoinMessageEol;
      // probably faster to use a basis to get integer solutions
      model.setSpecialOptions(2);
      model.branchAndBound();
      if (logLevel>1)
        model_->messageHandler()->message(CBC_END_SUB,model_->messages())
          << name
          <<CoinMessageEol;
      if (model.getMinimizationObjValue()<cutoff) {
        // solution
        returnCode=1;
        // post process
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
        returnCode=0; // so will be infeasible
      }
    }
  } else {
    returnCode=0;
  }
  return returnCode;
}

// Default Constructor
CbcRounding::CbcRounding() 
  :CbcHeuristic()
{
  // matrix and row copy will automatically be empty
}

// Constructor from model
CbcRounding::CbcRounding(CbcModel & model)
  :CbcHeuristic(model)
{
  // Get a copy of original matrix (and by row for rounding);
  assert(model.solver());
  matrix_ = *model.solver()->getMatrixByCol();
  matrixByRow_ = *model.solver()->getMatrixByRow();
  seed_=1;
}

// Destructor 
CbcRounding::~CbcRounding ()
{
}

// Clone
CbcHeuristic *
CbcRounding::clone() const
{
  return new CbcRounding(*this);
}

// Copy constructor 
CbcRounding::CbcRounding(const CbcRounding & rhs)
:
  CbcHeuristic(rhs),
  matrix_(rhs.matrix_),
  matrixByRow_(rhs.matrixByRow_),
  seed_(rhs.seed_)
{
  setWhen(rhs.when());
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

  int numberIntegers = model_->numberIntegers();
  const int * integerVariable = model_->integerVariable();
  int i;
  double direction = solver->getObjSense();
  double newSolutionValue = direction*solver->getObjValue();
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
        double randomNumber = CoinDrand48();
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
        bool isInteger = integerType[iColumn];
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
    double randomNumber = CoinDrand48();
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
        model_->numberObjects())
      setWhen(0);
  }
}

  
