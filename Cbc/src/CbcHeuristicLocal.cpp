// Copyright (C) 2002, International Business Machines
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
#include "CbcHeuristicLocal.hpp"
#include "CbcBranchActual.hpp"
#include "CbcStrategy.hpp"
#include "CglPreProcess.hpp"

// Default Constructor
CbcHeuristicLocal::CbcHeuristicLocal() 
  :CbcHeuristic()
{
  numberSolutions_=0;
  swap_=0;
  used_=NULL;
}

// Constructor with model - assumed before cuts

CbcHeuristicLocal::CbcHeuristicLocal(CbcModel & model)
  :CbcHeuristic(model)
{
  numberSolutions_=0;
  swap_=0;
  // Get a copy of original matrix
  assert(model.solver());
  matrix_ = *model.solver()->getMatrixByCol();
  int numberColumns = model.solver()->getNumCols();
  used_ = new char[numberColumns];
  memset(used_,0,numberColumns);
}

// Destructor 
CbcHeuristicLocal::~CbcHeuristicLocal ()
{
  delete [] used_;
}

// Clone
CbcHeuristic *
CbcHeuristicLocal::clone() const
{
  return new CbcHeuristicLocal(*this);
}
// Create C++ lines to get to current state
void 
CbcHeuristicLocal::generateCpp( FILE * fp) 
{
  CbcHeuristicLocal other;
  fprintf(fp,"0#include \"CbcHeuristicLocal.hpp\"\n");
  fprintf(fp,"3  CbcHeuristicLocal heuristicLocal(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp,"heuristicLocal");
  if (swap_!=other.swap_)
    fprintf(fp,"3  heuristicLocal.setSearchType(%d);\n",swap_);
  else
    fprintf(fp,"4  heuristicLocal.setSearchType(%d);\n",swap_);
  fprintf(fp,"3  cbcModel->addHeuristic(&heuristicLocal);\n");
}

// Copy constructor 
CbcHeuristicLocal::CbcHeuristicLocal(const CbcHeuristicLocal & rhs)
:
  CbcHeuristic(rhs),
  matrix_(rhs.matrix_),
  numberSolutions_(rhs.numberSolutions_),
  swap_(rhs.swap_)
{
  if (model_&&rhs.used_) {
    int numberColumns = model_->solver()->getNumCols();
    used_ = new char[numberColumns];
    memcpy(used_,rhs.used_,numberColumns);
  } else {
    used_=NULL;
  }
}

// Assignment operator 
CbcHeuristicLocal & 
CbcHeuristicLocal::operator=( const CbcHeuristicLocal& rhs)
{
  if (this!=&rhs) {
    CbcHeuristic::operator=(rhs);
    matrix_ = rhs.matrix_;
    numberSolutions_ = rhs.numberSolutions_;
    swap_ = rhs.swap_;
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

// Resets stuff if model changes
void 
CbcHeuristicLocal::resetModel(CbcModel * model)
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
// This version fixes stuff and does IP
int 
CbcHeuristicLocal::solutionFix(double & objectiveValue,
			    double * newSolution,
			    const int * keep)
{
  // See if to do
  if (!when()||(when()==1&&model_->phase()!=1))
    return 0; // switched off
  // Don't do if it was this heuristic which found solution!
  if (this==model_->lastHeuristic())
    return 0;
  OsiSolverInterface * newSolver = model_->continuousSolver()->clone();
  const double * colLower = newSolver->getColLower();
  //const double * colUpper = newSolver->getColUpper();

  int numberIntegers = model_->numberIntegers();
  const int * integerVariable = model_->integerVariable();
  
  int i;
  int nFix=0;
  for (i=0;i<numberIntegers;i++) {
    int iColumn=integerVariable[i];
    const OsiObject * object = model_->object(i);
    // get original bounds
    double originalLower;
    double originalUpper;
    getIntegerInformation( object,originalLower, originalUpper); 
    newSolver->setColLower(iColumn,CoinMax(colLower[iColumn],originalLower));
    if (!used_[iColumn]) {
      newSolver->setColUpper(iColumn,colLower[iColumn]);
      nFix++;
    }
  }
  int returnCode = smallBranchAndBound(newSolver,numberNodes_,newSolution,objectiveValue,
                                         objectiveValue,"CbcHeuristicLocal");
  if ((returnCode&2)!=0) {
    // could add cut
    returnCode &= ~2;
  }

  delete newSolver;
  return returnCode;
}
/*
  First tries setting a variable to better value.  If feasible then
  tries setting others.  If not feasible then tries swaps
  Returns 1 if solution, 0 if not */
int
CbcHeuristicLocal::solution(double & solutionValue,
			 double * betterSolution)
{

  if (numberSolutions_==model_->getSolutionCount())
    return 0;
  numberSolutions_=model_->getSolutionCount();
  if ((model_->getNumCols()>100000&&model_->getNumCols()>
       10*model_->getNumRows())||numberSolutions_<=1)
    return 0; // probably not worth it
  // worth trying

  OsiSolverInterface * solver = model_->solver();
  const double * rowLower = solver->getRowLower();
  const double * rowUpper = solver->getRowUpper();
  const double * solution = model_->bestSolution();
  if (!solution)
    return 0; // No solution found yet
  const double * objective = solver->getObjCoefficients();
  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance,primalTolerance);

  int numberRows = matrix_.getNumRows();

  int numberIntegers = model_->numberIntegers();
  const int * integerVariable = model_->integerVariable();
  
  int i;
  double direction = solver->getObjSense();
  double newSolutionValue = model_->getObjValue()*direction;
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

  // way is 1 if down possible, 2 if up possible, 3 if both possible
  char * way = new char[numberIntegers];
  // corrected costs
  double * cost = new double[numberIntegers];
  // for array to mark infeasible rows after iColumn branch
  char * mark = new char[numberRows];
  memset(mark,0,numberRows);
  // space to save values so we don't introduce rounding errors
  double * save = new double[numberRows];

  // clean solution
  for (i=0;i<numberIntegers;i++) {
    int iColumn = integerVariable[i];
    const OsiObject * object = model_->object(i);
    // get original bounds
    double originalLower;
    double originalUpper;
    getIntegerInformation( object,originalLower, originalUpper); 
    double value=newSolution[iColumn];
    if (value<originalLower) {
      value=originalLower;
      newSolution[iColumn]=value;
    } else if (value>originalUpper) {
      value=originalUpper;
      newSolution[iColumn]=value;
    }
    double nearest=floor(value+0.5);
    //assert(fabs(value-nearest)<10.0*primalTolerance);
    value=nearest;
    newSolution[iColumn]=nearest;
    // if away from lower bound mark that fact
    if (nearest>originalLower) {
      used_[iColumn]=1;
    }
    cost[i] = direction*objective[iColumn];
    int iway=0;
    
    if (value>originalLower+0.5) 
      iway = 1;
    if (value<originalUpper-0.5) 
      iway |= 2;
    way[i]=iway;
  }
  // get row activities
  double * rowActivity = new double[numberRows];
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
  // check was feasible - if not adjust (cleaning may move)
  // if very infeasible then give up
  bool tryHeuristic=true;
  for (i=0;i<numberRows;i++) {
    if(rowActivity[i]<rowLower[i]) {
      if (rowActivity[i]<rowLower[i]-10.0*primalTolerance)
	tryHeuristic=false;
      rowActivity[i]=rowLower[i];
    } else if(rowActivity[i]>rowUpper[i]) {
      if (rowActivity[i]<rowUpper[i]+10.0*primalTolerance)
	tryHeuristic=false;
      rowActivity[i]=rowUpper[i];
    }
  }
  // Switch off if may take too long
  if (model_->getNumCols()>10000&&model_->getNumCols()>
      10*model_->getNumRows())
    tryHeuristic=false;
  if (tryHeuristic) {
    
    // best change in objective
    double bestChange=0.0;
    
    for (i=0;i<numberIntegers;i++) {
      int iColumn = integerVariable[i];
      
      double objectiveCoefficient = cost[i];
      int k;
      int j;
      int goodK=-1;
      int wayK=-1,wayI=-1;
      if ((way[i]&1)!=0) {
	int numberInfeasible=0;
	// save row activities and adjust
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  save[iRow]=rowActivity[iRow];
	  rowActivity[iRow] -= element[j];
	  if(rowActivity[iRow]<rowLower[iRow]-primalTolerance||
	     rowActivity[iRow]>rowUpper[iRow]+primalTolerance) {
	    // mark row
	    mark[iRow]=1;
	    numberInfeasible++;
	  }
	}
	// try down
	for (k=i+1;k<numberIntegers;k++) {
	  if ((way[k]&1)!=0) {
	    // try down
	    if (-objectiveCoefficient-cost[k]<bestChange) {
	      // see if feasible down
	      bool good=true;
	      int numberMarked=0;
	      int kColumn = integerVariable[k];
	      for (j=columnStart[kColumn];
		   j<columnStart[kColumn]+columnLength[kColumn];j++) {
		int iRow = row[j];
		double newValue = rowActivity[iRow] - element[j];
		if(newValue<rowLower[iRow]-primalTolerance||
		   newValue>rowUpper[iRow]+primalTolerance) {
		  good=false;
		  break;
		} else if (mark[iRow]) {
		  // made feasible
		  numberMarked++;
		}
	      }
	      if (good&&numberMarked==numberInfeasible) {
		// better solution
		goodK=k;
		wayK=-1;
		wayI=-1;
		bestChange = -objectiveCoefficient-cost[k];
	      }
	    }
	  }
	  if ((way[k]&2)!=0) {
	    // try up
	    if (-objectiveCoefficient+cost[k]<bestChange) {
	      // see if feasible up
	      bool good=true;
	      int numberMarked=0;
	      int kColumn = integerVariable[k];
	      for (j=columnStart[kColumn];
		   j<columnStart[kColumn]+columnLength[kColumn];j++) {
		int iRow = row[j];
		double newValue = rowActivity[iRow] + element[j];
		if(newValue<rowLower[iRow]-primalTolerance||
		   newValue>rowUpper[iRow]+primalTolerance) {
		  good=false;
		  break;
		} else if (mark[iRow]) {
		  // made feasible
		  numberMarked++;
		}
	      }
	      if (good&&numberMarked==numberInfeasible) {
		// better solution
		goodK=k;
		wayK=1;
		wayI=-1;
		bestChange = -objectiveCoefficient+cost[k];
	      }
	    }
	  }
	}
	// restore row activities
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  rowActivity[iRow] = save[iRow];
	  mark[iRow]=0;
	}
      }
      if ((way[i]&2)!=0) {
	int numberInfeasible=0;
	// save row activities and adjust
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  save[iRow]=rowActivity[iRow];
	  rowActivity[iRow] += element[j];
	  if(rowActivity[iRow]<rowLower[iRow]-primalTolerance||
	     rowActivity[iRow]>rowUpper[iRow]+primalTolerance) {
	    // mark row
	    mark[iRow]=1;
	    numberInfeasible++;
	  }
	}
	// try up
	for (k=i+1;k<numberIntegers;k++) {
	  if ((way[k]&1)!=0) {
	    // try down
	    if (objectiveCoefficient-cost[k]<bestChange) {
	      // see if feasible down
	      bool good=true;
	      int numberMarked=0;
	      int kColumn = integerVariable[k];
	      for (j=columnStart[kColumn];
		   j<columnStart[kColumn]+columnLength[kColumn];j++) {
		int iRow = row[j];
		double newValue = rowActivity[iRow] - element[j];
		if(newValue<rowLower[iRow]-primalTolerance||
		   newValue>rowUpper[iRow]+primalTolerance) {
		  good=false;
		  break;
		} else if (mark[iRow]) {
		  // made feasible
		  numberMarked++;
		}
	      }
	      if (good&&numberMarked==numberInfeasible) {
		// better solution
		goodK=k;
		wayK=-1;
		wayI=1;
		bestChange = objectiveCoefficient-cost[k];
	      }
	    }
	  }
	  if ((way[k]&2)!=0) {
	    // try up
	    if (objectiveCoefficient+cost[k]<bestChange) {
	      // see if feasible up
	      bool good=true;
	      int numberMarked=0;
	      int kColumn = integerVariable[k];
	      for (j=columnStart[kColumn];
		   j<columnStart[kColumn]+columnLength[kColumn];j++) {
		int iRow = row[j];
		double newValue = rowActivity[iRow] + element[j];
		if(newValue<rowLower[iRow]-primalTolerance||
		   newValue>rowUpper[iRow]+primalTolerance) {
		  good=false;
		  break;
		} else if (mark[iRow]) {
		  // made feasible
		  numberMarked++;
		}
	      }
	      if (good&&numberMarked==numberInfeasible) {
		// better solution
		goodK=k;
		wayK=1;
		wayI=1;
		bestChange = objectiveCoefficient+cost[k];
	      }
	    }
	  }
	}
	// restore row activities
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  rowActivity[iRow] = save[iRow];
	  mark[iRow]=0;
	}
      }
      if (goodK>=0) {
	// we found something - update solution
	for (j=columnStart[iColumn];
	     j<columnStart[iColumn]+columnLength[iColumn];j++) {
	  int iRow = row[j];
	  rowActivity[iRow]  += wayI * element[j];
	}
	newSolution[iColumn] += wayI;
	int kColumn = integerVariable[goodK];
	for (j=columnStart[kColumn];
	     j<columnStart[kColumn]+columnLength[kColumn];j++) {
	  int iRow = row[j];
	  rowActivity[iRow]  += wayK * element[j];
	}
	newSolution[kColumn] += wayK;
	// See if k can go further ?
	const OsiObject * object = model_->object(goodK);
	// get original bounds
	double originalLower;
	double originalUpper;
	getIntegerInformation( object,originalLower, originalUpper); 
	
	double value=newSolution[kColumn];
	int iway=0;
	
	if (value>originalLower+0.5) 
	  iway = 1;
	if (value<originalUpper-0.5) 
	  iway |= 2;
	way[goodK]=iway;
      }
    }
    if (bestChange+newSolutionValue<solutionValue) {
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
      int numberBad=0;
      double sumBad=0.0;
      // check was approximately feasible
      for (i=0;i<numberRows;i++) {
	if(rowActivity[i]<rowLower[i]) {
          sumBad += rowLower[i]-rowActivity[i];
	  if (rowActivity[i]<rowLower[i]-10.0*primalTolerance)
            numberBad++;
	} else if(rowActivity[i]>rowUpper[i]) {
          sumBad += rowUpper[i]-rowActivity[i];
	  if (rowActivity[i]>rowUpper[i]+10.0*primalTolerance)
            numberBad++;
	}
      }
      if (!numberBad) {
        for (i=0;i<numberIntegers;i++) {
          int iColumn = integerVariable[i];
          const OsiObject * object = model_->object(i);
          // get original bounds
	  double originalLower;
	  double originalUpper;
	  getIntegerInformation( object,originalLower, originalUpper); 
          
          double value=newSolution[iColumn];
          // if away from lower bound mark that fact
          if (value>originalLower) {
            used_[iColumn]=1;
          }
        }
        // new solution
        memcpy(betterSolution,newSolution,numberColumns*sizeof(double));
        returnCode=1;
        solutionValue = newSolutionValue + bestChange;
      } else {
        // bad solution - should not happen so debug if see message
        printf("Local search got bad solution with %d infeasibilities summing to %g\n",
               numberBad,sumBad);
      }
    }
  }
  delete [] newSolution;
  delete [] rowActivity;
  delete [] way;
  delete [] cost;
  delete [] save;
  delete [] mark;
  if (numberSolutions_>1&&swap_==1) {
    // try merge
    int returnCode2=solutionFix( solutionValue, betterSolution,NULL);
    if (returnCode2)
      returnCode=1;
  }
  return returnCode;
}
// update model
void CbcHeuristicLocal::setModel(CbcModel * model)
{
  model_ = model;
  // Get a copy of original matrix
  assert(model_->solver());
  matrix_ = *model_->solver()->getMatrixByCol();
  delete [] used_;
  int numberColumns = model->solver()->getNumCols();
  used_ = new char[numberColumns];
  memset(used_,0,numberColumns);
}

  
