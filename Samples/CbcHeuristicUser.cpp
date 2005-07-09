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
#include "CbcHeuristicUser.hpp"
#include "CbcBranchActual.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcCompareUser.hpp"
// Cuts

#include "CglProbing.hpp"
// Default Constructor
CbcLocalSearch::CbcLocalSearch() 
  :CbcHeuristic()
{
  numberSolutions_=0;
  swap_=0;
  used_=NULL;
}

// Constructor with model - assumed before cuts

CbcLocalSearch::CbcLocalSearch(CbcModel & model)
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
CbcLocalSearch::~CbcLocalSearch ()
{
  delete [] used_;
}

// Clone
CbcHeuristic *
CbcLocalSearch::clone() const
{
  return new CbcLocalSearch(*this);
}

// Copy constructor 
CbcLocalSearch::CbcLocalSearch(const CbcLocalSearch & rhs)
:
  CbcHeuristic(rhs),
  matrix_(rhs.matrix_),
  numberSolutions_(rhs.numberSolutions_),
  swap_(rhs.swap_)
{
  setWhen(rhs.when());
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
CbcLocalSearch::resetModel(CbcModel * model)
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
CbcLocalSearch::solutionFix(double & objectiveValue,
			    double * newSolution,
			    const int * keep)
{
  // See if to do
  if (!when()||(when()==1&&model_->phase()!=1))
    return 0; // switched off
  OsiSolverInterface * solver = model_->continuousSolver()->clone();
  const double * colLower = solver->getColLower();
  //const double * colUpper = solver->getColUpper();

  int numberIntegers = model_->numberIntegers();
  const int * integerVariable = model_->integerVariable();
  
  int i;
  int nFix=0;
  for (i=0;i<numberIntegers;i++) {
    int iColumn=integerVariable[i];
    const CbcObject * object = model_->object(i);
    const CbcSimpleInteger * integerObject = 
      dynamic_cast<const  CbcSimpleInteger *> (object);
    assert(integerObject);
    // get original bounds
    double originalLower = integerObject->originalLowerBound();
    solver->setColLower(iColumn,CoinMax(colLower[iColumn],originalLower));
    if (!used_[iColumn]) {
      solver->setColUpper(iColumn,colLower[iColumn]);
      nFix++;
    }
  }
  CbcModel model(*solver);
  //model.solver()->writeMps("large");
  // integer presolve
  CbcModel * model2 = model.integerPresolve(false);
  if (!model2) {
    delete solver;
    return 0;
  }
  // Do complete search
  
  // Cuts
  // Probing first as gets tight bounds on continuous

  CglProbing generator1;
  generator1.setUsingObjective(true);
  generator1.setMaxPass(3);
  generator1.setMaxProbe(100);
  generator1.setMaxLook(50);
  generator1.setRowCuts(3);

  // Add in generators
  model2->addCutGenerator(&generator1,-1,"Probing");
  //model2->solver()->writeMps("small");
  // Definition of node choice
  CbcCompareUser compare;
  model2->setNodeComparison(compare);
  //model2->solver()->setHintParam(OsiDoReducePrint,true,OsiHintTry);
  model2->messageHandler()->setLogLevel(2);
  //CbcSolverUser * osiclp = dynamic_cast< Osi*> (model2->solver());
  //assert (osiclp);
  //ClpSimplex * clp = osiclp->getModelPtr();
  model2->solver()->messageHandler()->setLogLevel(2);
  model2->messageHandler()->setLogLevel(3);
  model2->setMaximumCutPassesAtRoot(-100); // always do 100 if possible
  model2->messageHandler()->setLogLevel(1);
  model2->setPrintFrequency(50);
  // ? bug if passed directly
  double cutoff = objectiveValue-1.0e-4;
  model2->setCutoff(cutoff);
  model2->setIntParam(CbcModel::CbcMaxNumNode,(1000*nFix)/20);
  model2->setIntParam(CbcModel::CbcMaxNumNode,500);
  model2->branchAndBound();
  // get back solution
  model.originalModel(model2,false);
  delete model2;

  // Get solution array for heuristic solution
  int numberColumns = solver->getNumCols();
  const double * solution = model.bestSolution();
  printf("obj %g old %g cutoff %g %g - fix %d\n",model.getObjValue(),
	 objectiveValue,model_->getObjValue(),model_->getCutoff(),
	 nFix);
  if (model.getObjValue()<objectiveValue) {
    objectiveValue=model.getObjValue();
    memcpy(newSolution,solution,numberColumns*sizeof(double));
    return 1;
  } else {
    return 0;
  }
}
/*
  First tries setting a variable to better value.  If feasible then
  tries setting others.  If not feasible then tries swaps
  Returns 1 if solution, 0 if not */
int
CbcLocalSearch::solution(double & solutionValue,
			 double * betterSolution)
{

  if (numberSolutions_==model_->getSolutionCount())
    return 0;
  if (model_->getNumCols()>1000&&model_->getNumCols()>
      10*model_->getNumRows())
    return 0; // probably not worth it
  // worth trying
  numberSolutions_=model_->getSolutionCount();

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
    const CbcObject * object = model_->object(i);
    const CbcSimpleInteger * integerObject = 
      dynamic_cast<const  CbcSimpleInteger *> (object);
    assert(integerObject);
    // get original bounds
    double originalLower = integerObject->originalLowerBound();
    double originalUpper = integerObject->originalUpperBound();

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
	const CbcObject * object = model_->object(goodK);
	const CbcSimpleInteger * integerObject = 
	  dynamic_cast<const  CbcSimpleInteger *> (object);
	// get original bounds
	double originalLower = integerObject->originalLowerBound();
	double originalUpper = integerObject->originalUpperBound();
	
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
          const CbcObject * object = model_->object(i);
          const CbcSimpleInteger * integerObject = 
            dynamic_cast<const  CbcSimpleInteger *> (object);
          // get original bounds
          double originalLower = integerObject->originalLowerBound();
          //double originalUpper = integerObject->originalUpperBound();
          
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
        if (bestChange>1.0e-12)
          printf("Local search heuristic improved solution by %g\n",
                 -bestChange);
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
void CbcLocalSearch::setModel(CbcModel * model)
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

  
