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
#include "CbcStrategy.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CoinSort.hpp"
#include "CglPreProcess.hpp"
// Default Constructor
CbcHeuristicGreedyCover::CbcHeuristicGreedyCover() 
  :CbcHeuristic()
{
  // matrix  will automatically be empty
  originalNumberRows_=0;
  algorithm_=0;
  numberTimes_=100;
}

// Constructor from model
CbcHeuristicGreedyCover::CbcHeuristicGreedyCover(CbcModel & model)
  :CbcHeuristic(model)
{
  // Get a copy of original matrix
  assert(model.solver());
  matrix_ = *model.solver()->getMatrixByCol();
  originalNumberRows_=model.solver()->getNumRows();
  algorithm_=0;
  numberTimes_=100;
}

// Destructor 
CbcHeuristicGreedyCover::~CbcHeuristicGreedyCover ()
{
}

// Clone
CbcHeuristic *
CbcHeuristicGreedyCover::clone() const
{
  return new CbcHeuristicGreedyCover(*this);
}

// Copy constructor 
CbcHeuristicGreedyCover::CbcHeuristicGreedyCover(const CbcHeuristicGreedyCover & rhs)
:
  CbcHeuristic(rhs),
  matrix_(rhs.matrix_),
  originalNumberRows_(rhs.originalNumberRows_),
  algorithm_(rhs.algorithm_),
  numberTimes_(rhs.numberTimes_)
{
  this->setWhen(rhs.when());
}

// Assignment operator 
CbcHeuristicGreedyCover & 
CbcHeuristicGreedyCover::operator=( const CbcHeuristicGreedyCover& rhs)
{
  if (this != &rhs) {
    setWhen(rhs.when());
    matrix_=rhs.matrix_;
    originalNumberRows_=rhs.originalNumberRows_;
    algorithm_=rhs.algorithm_;
    numberTimes_=rhs.numberTimes_;
  }
  return *this;
}
// Returns 1 if solution, 0 if not
int
CbcHeuristicGreedyCover::solution(double & solutionValue,
			 double * betterSolution)
{
  if (!model_)
    return 0;
  // See if to do
  if (!when()||(when()==1&&model_->phase()!=1))
    return 0; // switched off
  if (model_->getNodeCount()>numberTimes_)
    return 0;
  // See if at root node
  bool atRoot = model_->getNodeCount()==0;
  int passNumber = model_->getCurrentPassNumber();
  if (atRoot&&passNumber!=1)
    return 0;
  OsiSolverInterface * solver = model_->solver();
  const double * columnLower = solver->getColLower();
  const double * columnUpper = solver->getColUpper();
  // And original upper bounds in case we want to use them
  const double * originalUpper = model_->continuousSolver()->getColUpper();
  // But not if algorithm says so
  if ((algorithm_%10)==0)
    originalUpper = columnUpper;
  const double * rowLower = solver->getRowLower();
  const double * solution = solver->getColSolution();
  const double * objective = solver->getObjCoefficients();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance,primalTolerance);

  // This is number of rows when matrix was passed in
  int numberRows = originalNumberRows_;
  if (!numberRows)
    return 0; // switched off

  assert (numberRows==matrix_.getNumRows());
  int iRow, iColumn;
  double direction = solver->getObjSense();
  double offset;
  solver->getDblParam(OsiObjOffset,offset);
  double newSolutionValue = -offset;
  int returnCode = 0;

  // Column copy
  const double * element = matrix_.getElements();
  const int * row = matrix_.getIndices();
  const CoinBigIndex * columnStart = matrix_.getVectorStarts();
  const int * columnLength = matrix_.getVectorLengths();

  // Get solution array for heuristic solution
  int numberColumns = solver->getNumCols();
  double * newSolution = new double [numberColumns];
  double * rowActivity = new double[numberRows];
  memset(rowActivity,0,numberRows*sizeof(double));
  bool allOnes=true;
  // Get rounded down solution 
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    CoinBigIndex j;
    double value = solution[iColumn];
    if (solver->isInteger(iColumn)) {
      // Round down integer
      if (fabs(floor(value+0.5)-value)<integerTolerance) {
        value=floor(CoinMax(value+1.0e-3,columnLower[iColumn]));
      } else {
        value=CoinMax(floor(value),columnLower[iColumn]);
      }
    }
    // make sure clean
    value = CoinMin(value,columnUpper[iColumn]);
    value = CoinMax(value,columnLower[iColumn]);
    newSolution[iColumn]=value;
    double cost = direction * objective[iColumn];
    newSolutionValue += value*cost;
    for (j=columnStart[iColumn];
         j<columnStart[iColumn]+columnLength[iColumn];j++) {
      int iRow=row[j];
      rowActivity[iRow] += value*element[j];
      if (element[j]!=1.0)
        allOnes=false;
    }
  }
  // See if we round up
  bool roundup = ((algorithm_%100)!=0);
  if (roundup&&allOnes) {
    // Get rounded up solution 
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      CoinBigIndex j;
      double value = solution[iColumn];
      if (solver->isInteger(iColumn)) {
        // but round up if no activity
        if (roundup&&value>=0.499999&&!newSolution[iColumn]) {
          bool choose=true;
          for (j=columnStart[iColumn];
               j<columnStart[iColumn]+columnLength[iColumn];j++) {
            int iRow=row[j];
            if (rowActivity[iRow]) {
              choose=false;
              break;
            }
          }
          if (choose) {
            newSolution[iColumn]=1.0;
            double cost = direction * objective[iColumn];
            newSolutionValue += cost;
            for (j=columnStart[iColumn];
                 j<columnStart[iColumn]+columnLength[iColumn];j++) {
              int iRow=row[j];
              rowActivity[iRow] += 1.0;
            }
          }
        }
      }
    }
  }
  // Get initial list
  int * which = new int [numberColumns];
  for (iColumn=0;iColumn<numberColumns;iColumn++) 
    which[iColumn]=iColumn;
  int numberLook=numberColumns;
  // See if we want to perturb more
  double perturb = ((algorithm_%10)==0) ? 0.1 : 0.25;
  // Keep going round until a solution
  while (true) {
    // Get column with best ratio
    int bestColumn=-1;
    double bestRatio=COIN_DBL_MAX;
    double bestStepSize = 0.0;
    int newNumber=0;
    for (int jColumn=0;jColumn<numberLook;jColumn++) {
      int iColumn = which[jColumn];
      CoinBigIndex j;
      double value = newSolution[iColumn];
      double cost = direction * objective[iColumn];
      if (solver->isInteger(iColumn)) {
        // use current upper or original upper
        if (value+0.99<originalUpper[iColumn]) {
          double sum=0.0;
          int numberExact=0;
          for (j=columnStart[iColumn];
               j<columnStart[iColumn]+columnLength[iColumn];j++) {
            int iRow=row[j];
            double gap = rowLower[iRow]-rowActivity[iRow];
            double elementValue = allOnes ? 1.0 : element[j];
            if (gap>1.0e-7) {
              sum += CoinMin(elementValue,gap);
              if (fabs(elementValue-gap)<1.0e-7) 
                numberExact++;
	    }
	  }
          // could bias if exact
	  if (sum>0.0) {
            // add to next time
            which[newNumber++]=iColumn;
	    double ratio = (cost/sum)*(1.0+perturb*CoinDrand48());
            // If at root choose first
            if (atRoot)
              ratio = iColumn;
	    if (ratio<bestRatio) {
              bestRatio=ratio;
              bestColumn=iColumn;
              bestStepSize=1.0;
            }
	  }
	}
      } else {
        // continuous
        if (value<columnUpper[iColumn]) {
          // Go through twice - first to get step length
          double step=1.0e50;
          for (j=columnStart[iColumn];
               j<columnStart[iColumn]+columnLength[iColumn];j++) {
            int iRow=row[j];
            if (rowActivity[iRow]<rowLower[iRow]-1.0e-10&&
                element[j]*step+rowActivity[iRow]>=rowLower[iRow]) {
	      step = (rowLower[iRow]-rowActivity[iRow])/element[j];;
	    }
	  }
          // now ratio
          if (step<1.0e50) {
            // add to next time
            which[newNumber++]=iColumn;
            assert (step>0.0);
            double sum=0.0;
            for (j=columnStart[iColumn];
                 j<columnStart[iColumn]+columnLength[iColumn];j++) {
              int iRow=row[j];
              if (rowActivity[iRow]<rowLower[iRow]-1.0e-10&&
                element[j]*step+rowActivity[iRow]>=rowLower[iRow]) {
                sum += element[j];
              }
              assert (sum>0.0);
              double ratio = (cost/sum)*(1.0+perturb*CoinDrand48());
              if (ratio<bestRatio) {
                bestRatio=ratio;
                bestColumn=iColumn;
                bestStepSize=step;
              }
            }
          }
        }
      }
    }
    if (bestColumn<0)
      break; // we have finished
    // Increase chosen column
    newSolution[bestColumn] += bestStepSize;
    double cost = direction * objective[bestColumn];
    newSolutionValue += bestStepSize*cost;
    for (CoinBigIndex j=columnStart[bestColumn];
         j<columnStart[bestColumn]+columnLength[bestColumn];j++) {
      int iRow = row[j];
      rowActivity[iRow] += bestStepSize*element[j];
    }
  }
  delete [] which;
  if (newSolutionValue<solutionValue) {
    // check feasible
    memset(rowActivity,0,numberRows*sizeof(double));
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      CoinBigIndex j;
      double value = newSolution[iColumn];
      if (value) {
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow=row[j];
	  rowActivity[iRow] += value*element[j];
	}
      }
    }
    // check was approximately feasible
    bool feasible=true;
    for (iRow=0;iRow<numberRows;iRow++) {
      if(rowActivity[iRow]<rowLower[iRow]) {
	if (rowActivity[iRow]<rowLower[iRow]-10.0*primalTolerance)
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
      //printf("Debug CbcHeuristicGreedyCover giving bad solution\n");
    }
  }
  delete [] newSolution;
  delete [] rowActivity;
  return returnCode;
}
// update model
void CbcHeuristicGreedyCover::setModel(CbcModel * model)
{
#define SLOPPY
#ifndef SLOPPY
  model_ = model;
  assert(model_->solver());
  *this = CbcHeuristicGreedyCover(*model); 
#else
  if (model_&&model!=model_) {
    model_ = model;
    assert(model_->solver());
    *this = CbcHeuristicGreedyCover(*model); 
  }
#endif
  validate();
}
// Resets stuff if model changes
void 
CbcHeuristicGreedyCover::resetModel(CbcModel * model)
{
#ifndef SLOPPY
  model_ = model;
  assert(model_->solver());
  *this = CbcHeuristicGreedyCover(*model); 
#else
  // switch off
  model_ = NULL;
  matrix_ = CoinPackedMatrix();
#endif
}
// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
void 
CbcHeuristicGreedyCover::validate() 
{
  if (model_&&when()<10) {
    if (model_->numberIntegers()!=
        model_->numberObjects())
      setWhen(0);
    // Only works if costs positive, coefficients positive and all rows G
    OsiSolverInterface * solver = model_->solver();
    const double * columnLower = solver->getColLower();
    const double * rowUpper = solver->getRowUpper();
    const double * objective = solver->getObjCoefficients();
    double direction = solver->getObjSense();

    int numberRows = solver->getNumRows();
    // Column copy
    const double * element = matrix_.getElements();
    const CoinBigIndex * columnStart = matrix_.getVectorStarts();
    const int * columnLength = matrix_.getVectorLengths();
    bool good = true;
    for (int iRow=0;iRow<numberRows;iRow++) {
      if (rowUpper[iRow]<1.0e30)
        good = false;
    }
    int numberColumns = solver->getNumCols();
    for (int iColumn=0;iColumn<numberColumns;iColumn++) {
      if (objective[iColumn]*direction<0.0)
        good=false;
      if (columnLower[iColumn]<0.0)
        good=false;
      CoinBigIndex j;
      for (j=columnStart[iColumn];
           j<columnStart[iColumn]+columnLength[iColumn];j++) {
        if (element[j]<0.0)
          good=false;
      }
    }
    if (!good)
      setWhen(0); // switch off
  }
}
// Default Constructor
CbcHeuristicGreedyEquality::CbcHeuristicGreedyEquality() 
  :CbcHeuristic()
{
  // matrix  will automatically be empty
  fraction_=1.0; // no branch and bound
  originalNumberRows_=0;
  algorithm_=0;
  numberTimes_=100;
}

// Constructor from model
CbcHeuristicGreedyEquality::CbcHeuristicGreedyEquality(CbcModel & model)
  :CbcHeuristic(model)
{
  // Get a copy of original matrix
  assert(model.solver());
  matrix_ = *model.solver()->getMatrixByCol();
  fraction_=1.0; // no branch and bound
  originalNumberRows_=model.solver()->getNumRows();
  algorithm_=0;
  numberTimes_=100;
}

// Destructor 
CbcHeuristicGreedyEquality::~CbcHeuristicGreedyEquality ()
{
}

// Clone
CbcHeuristic *
CbcHeuristicGreedyEquality::clone() const
{
  return new CbcHeuristicGreedyEquality(*this);
}

// Copy constructor 
CbcHeuristicGreedyEquality::CbcHeuristicGreedyEquality(const CbcHeuristicGreedyEquality & rhs)
:
  CbcHeuristic(rhs),
  matrix_(rhs.matrix_),
  fraction_(rhs.fraction_),
  originalNumberRows_(rhs.originalNumberRows_),
  algorithm_(rhs.algorithm_),
  numberTimes_(rhs.numberTimes_)
{
  this->setWhen(rhs.when());
}

// Assignment operator 
CbcHeuristicGreedyEquality & 
CbcHeuristicGreedyEquality::operator=( const CbcHeuristicGreedyEquality& rhs)
{
  if (this != &rhs) {
    setWhen(rhs.when());
    matrix_=rhs.matrix_;
    fraction_ = rhs.fraction_;
    originalNumberRows_=rhs.originalNumberRows_;
    algorithm_=rhs.algorithm_;
    numberTimes_=rhs.numberTimes_;
  }
  return *this;
}
// Returns 1 if solution, 0 if not
int
CbcHeuristicGreedyEquality::solution(double & solutionValue,
			 double * betterSolution)
{
  if (!model_)
    return 0;
  // See if to do
  if (!when()||(when()==1&&model_->phase()!=1))
    return 0; // switched off
  if (model_->getNodeCount()>numberTimes_)
    return 0;
  // See if at root node
  bool atRoot = model_->getNodeCount()==0;
  int passNumber = model_->getCurrentPassNumber();
  if (atRoot&&passNumber!=1)
    return 0;
  OsiSolverInterface * solver = model_->solver();
  const double * columnLower = solver->getColLower();
  const double * columnUpper = solver->getColUpper();
  // And original upper bounds in case we want to use them
  const double * originalUpper = model_->continuousSolver()->getColUpper();
  // But not if algorithm says so
  if ((algorithm_%10)==0)
    originalUpper = columnUpper;
  const double * rowLower = solver->getRowLower();
  const double * rowUpper = solver->getRowUpper();
  const double * solution = solver->getColSolution();
  const double * objective = solver->getObjCoefficients();
  double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance,primalTolerance);

  // This is number of rows when matrix was passed in
  int numberRows = originalNumberRows_;
  if (!numberRows)
    return 0; // switched off

  assert (numberRows==matrix_.getNumRows());
  int iRow, iColumn;
  double direction = solver->getObjSense();
  double offset;
  solver->getDblParam(OsiObjOffset,offset);
  double newSolutionValue = -offset;
  int returnCode = 0;

  // Column copy
  const double * element = matrix_.getElements();
  const int * row = matrix_.getIndices();
  const CoinBigIndex * columnStart = matrix_.getVectorStarts();
  const int * columnLength = matrix_.getVectorLengths();

  // Get solution array for heuristic solution
  int numberColumns = solver->getNumCols();
  double * newSolution = new double [numberColumns];
  double * rowActivity = new double[numberRows];
  memset(rowActivity,0,numberRows*sizeof(double));
  double rhsNeeded=0;
  for (iRow=0;iRow<numberRows;iRow++) 
    rhsNeeded += rowUpper[iRow];
  rhsNeeded *= fraction_;
  bool allOnes=true;
  // Get rounded down solution 
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    CoinBigIndex j;
    double value = solution[iColumn];
    if (solver->isInteger(iColumn)) {
      // Round down integer
      if (fabs(floor(value+0.5)-value)<integerTolerance) {
        value=floor(CoinMax(value+1.0e-3,columnLower[iColumn]));
      } else {
        value=CoinMax(floor(value),columnLower[iColumn]);
      }
    }
    // make sure clean
    value = CoinMin(value,columnUpper[iColumn]);
    value = CoinMax(value,columnLower[iColumn]);
    newSolution[iColumn]=value;
    double cost = direction * objective[iColumn];
    newSolutionValue += value*cost;
    for (j=columnStart[iColumn];
         j<columnStart[iColumn]+columnLength[iColumn];j++) {
      int iRow=row[j];
      rowActivity[iRow] += value*element[j];
      rhsNeeded -= value*element[j];
      if (element[j]!=1.0)
        allOnes=false;
    }
  }
  // See if we round up
  bool roundup = ((algorithm_%100)!=0);
  if (roundup&&allOnes) {
    // Get rounded up solution 
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      CoinBigIndex j;
      double value = solution[iColumn];
      if (solver->isInteger(iColumn)) {
        // but round up if no activity
        if (roundup&&value>=0.6&&!newSolution[iColumn]) {
          bool choose=true;
          for (j=columnStart[iColumn];
               j<columnStart[iColumn]+columnLength[iColumn];j++) {
            int iRow=row[j];
            if (rowActivity[iRow]) {
              choose=false;
              break;
            }
          }
          if (choose) {
            newSolution[iColumn]=1.0;
            double cost = direction * objective[iColumn];
            newSolutionValue += cost;
            for (j=columnStart[iColumn];
                 j<columnStart[iColumn]+columnLength[iColumn];j++) {
              int iRow=row[j];
              rowActivity[iRow] += 1.0;
              rhsNeeded -= 1.0;
            }
          }
        }
      }
    }
  }
  // Get initial list
  int * which = new int [numberColumns];
  for (iColumn=0;iColumn<numberColumns;iColumn++) 
    which[iColumn]=iColumn;
  int numberLook=numberColumns;
  // See if we want to perturb more
  double perturb = ((algorithm_%10)==0) ? 0.1 : 0.25;
  // Keep going round until a solution
  while (true) {
    // Get column with best ratio
    int bestColumn=-1;
    double bestRatio=COIN_DBL_MAX;
    double bestStepSize = 0.0;
    int newNumber=0;
    for (int jColumn=0;jColumn<numberLook;jColumn++) {
      int iColumn = which[jColumn];
      CoinBigIndex j;
      double value = newSolution[iColumn];
      double cost = direction * objective[iColumn];
      if (solver->isInteger(iColumn)) {
        // use current upper or original upper
        if (value+0.9999<originalUpper[iColumn]) {
          double movement = 1.0;
          double sum=0.0;
          for (j=columnStart[iColumn];
               j<columnStart[iColumn]+columnLength[iColumn];j++) {
            int iRow=row[j];
            double gap = rowUpper[iRow]-rowActivity[iRow];
            double elementValue = allOnes ? 1.0 : element[j];
            sum += elementValue;
            if (movement*elementValue>gap) {
              movement = gap/elementValue;
            }
	  }
	  if (movement>0.999999) {
            // add to next time
            which[newNumber++]=iColumn;
	    double ratio = (cost/sum)*(1.0+perturb*CoinDrand48());
            // If at root
            if (atRoot) {
              if (fraction_==1.0)
                ratio = iColumn; // choose first
              else
                ratio = - solution[iColumn]; // choose largest
            }
	    if (ratio<bestRatio) {
              bestRatio=ratio;
              bestColumn=iColumn;
              bestStepSize=1.0;
            }
	  }
	}
      } else {
        // continuous
        if (value<columnUpper[iColumn]) {
          double movement=1.0e50;
          double sum=0.0;
          for (j=columnStart[iColumn];
               j<columnStart[iColumn]+columnLength[iColumn];j++) {
            int iRow=row[j];
            if (element[j]*movement+rowActivity[iRow]>rowUpper[iRow]) {
	      movement = (rowUpper[iRow]-rowActivity[iRow])/element[j];;
	    }
            sum += element[j];
	  }
          // now ratio
          if (movement>1.0e-7) {
            // add to next time
            which[newNumber++]=iColumn;
            double ratio = (cost/sum)*(1.0+perturb*CoinDrand48());
            if (ratio<bestRatio) {
              bestRatio=ratio;
              bestColumn=iColumn;
              bestStepSize=movement;
            }
          }
        }
      }
    }
    if (bestColumn<0)
      break; // we have finished
    // Increase chosen column
    newSolution[bestColumn] += bestStepSize;
    double cost = direction * objective[bestColumn];
    newSolutionValue += bestStepSize*cost;
    for (CoinBigIndex j=columnStart[bestColumn];
         j<columnStart[bestColumn]+columnLength[bestColumn];j++) {
      int iRow = row[j];
      rowActivity[iRow] += bestStepSize*element[j];
      rhsNeeded -= bestStepSize*element[j];
    }
    if (rhsNeeded<1.0e-8)
      break;
  }
  delete [] which;
  if (fraction_<1.0&&rhsNeeded<1.0e-8&&newSolutionValue<solutionValue) {
    // do branch and cut
    // fix all nonzero
    OsiSolverInterface * newSolver = model_->continuousSolver()->clone();
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      if (newSolver->isInteger(iColumn))
        newSolver->setColLower(iColumn,newSolution[iColumn]);
    }
    // Reduce printout
    newSolver->setHintParam(OsiDoReducePrint,true,OsiHintTry);
    newSolver->setHintParam(OsiDoPresolveInInitial,true,OsiHintTry);
    newSolver->setDblParam(OsiDualObjectiveLimit,solutionValue);
    newSolver->initialSolve();
    if (newSolver->isProvenOptimal()) {
      CglPreProcess process;
      /* Do not try and produce equality cliques and
         do up to 5 passes */
      OsiSolverInterface * solver2= process.preProcess(*newSolver);
      if (!solver2) {
        printf("Pre-processing says infeasible\n");
        rhsNeeded=1.0; // so will be infeasible
      } else {
        solver2->resolve();
        CbcModel model(*solver2);
        model.setCutoff(solutionValue);
        model.setMaximumNodes(200);
        model.solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
        CbcStrategyDefault strategy(true,5,5,0);
        model.setStrategy(strategy);
        // Lightweight
        model.setNumberStrong(5);
        model.setNumberBeforeTrust(1);
        model.solver()->setIntParam(OsiMaxNumIterationHotStart,10);
        // Do complete search
        model.branchAndBound();
        if (model.getMinimizationObjValue()<solutionValue) {
          // solution
          rhsNeeded=0.0;
          // post process
          process.postProcess(*model.solver());
          // Solution now back in newSolver
          memcpy(newSolution,newSolver->getColSolution(),
                 numberColumns*sizeof(double));
          newSolutionValue = model.getMinimizationObjValue();
        } else {
          // no good
          rhsNeeded=1.0; // so will be infeasible
        }
      }
    } else {
      rhsNeeded=1.0;
    }
    delete newSolver;
  }
  if (newSolutionValue<solutionValue&&rhsNeeded<1.0e-8) {
    // check feasible
    memset(rowActivity,0,numberRows*sizeof(double));
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      CoinBigIndex j;
      double value = newSolution[iColumn];
      if (value) {
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow=row[j];
	  rowActivity[iRow] += value*element[j];
	}
      }
    }
    // check was approximately feasible
    bool feasible=true;
    for (iRow=0;iRow<numberRows;iRow++) {
      if(rowActivity[iRow]<rowLower[iRow]) {
	if (rowActivity[iRow]<rowLower[iRow]-10.0*primalTolerance)
	  feasible = false;
      }
    }
    if (feasible) {
      // new solution
      memcpy(betterSolution,newSolution,numberColumns*sizeof(double));
      solutionValue = newSolutionValue;
      printf("** Solution of %g found by greedy\n",newSolutionValue);
      returnCode=1;
    }
  }
  delete [] newSolution;
  delete [] rowActivity;
  if (atRoot&&fraction_==1.0) {
    // try quick search
    fraction_=0.3;
    int newCode=this->solution(solutionValue,betterSolution);
    if (newCode)
      returnCode=1;
    fraction_=1.0;
  }
  return returnCode;
}
// update model
void CbcHeuristicGreedyEquality::setModel(CbcModel * model)
{
#define SLOPPY
#ifndef SLOPPY
  model_ = model;
  assert(model_->solver());
  *this = CbcHeuristicGreedyEquality(*model); 
#else
  if (model_&&model!=model_) {
    model_ = model;
    assert(model_->solver());
    *this = CbcHeuristicGreedyEquality(*model); 
  }
#endif
  validate();
}
// Resets stuff if model changes
void 
CbcHeuristicGreedyEquality::resetModel(CbcModel * model)
{
#ifndef SLOPPY
  model_ = model;
  assert(model_->solver());
  *this = CbcHeuristicGreedyEquality(*model); 
#else
  // switch off
  model_ = NULL;
  matrix_ = CoinPackedMatrix();
#endif
}
// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
void 
CbcHeuristicGreedyEquality::validate() 
{
  if (model_&&when()<10) {
    if (model_->numberIntegers()!=
        model_->numberObjects())
      setWhen(0);
    // Only works if costs positive, coefficients positive and all rows E or L
    // And if values are integer
    OsiSolverInterface * solver = model_->solver();
    const double * columnLower = solver->getColLower();
    const double * rowUpper = solver->getRowUpper();
    const double * rowLower = solver->getRowLower();
    const double * objective = solver->getObjCoefficients();
    double direction = solver->getObjSense();

    int numberRows = solver->getNumRows();
    // Column copy
    const double * element = matrix_.getElements();
    const CoinBigIndex * columnStart = matrix_.getVectorStarts();
    const int * columnLength = matrix_.getVectorLengths();
    bool good = true;
    for (int iRow=0;iRow<numberRows;iRow++) {
      if (rowUpper[iRow]>1.0e30)
        good = false;
      if (rowLower[iRow]>0.0&&rowLower[iRow]!=rowUpper[iRow])
        good = false;
      if (floor(rowUpper[iRow]+0.5)!=rowUpper[iRow])
        good = false;
    }
    int numberColumns = solver->getNumCols();
    for (int iColumn=0;iColumn<numberColumns;iColumn++) {
      if (objective[iColumn]*direction<0.0)
        good=false;
      if (columnLower[iColumn]<0.0)
        good=false;
      CoinBigIndex j;
      for (j=columnStart[iColumn];
           j<columnStart[iColumn]+columnLength[iColumn];j++) {
        if (element[j]<0.0)
          good=false;
        if (floor(element[j]+0.5)!=element[j])
          good = false;
      }
    }
    if (!good)
      setWhen(0); // switch off
  }
}

  
