// Copyright (C) 2004, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <cassert>

#include "CoinTime.hpp"

#include "CoinHelperFunctions.hpp"
#include "CoinIndexedVector.hpp"
#include "ClpFactorization.hpp"
#include "ClpObjective.hpp"
#include "ClpSimplex.hpp"
#include "CbcSolverLongThin.hpp"
#include "CbcModel.hpp"
#include "ClpPresolve.hpp"
#include "CbcHeuristicUser.hpp"
#include "CbcBranchActual.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcCompareUser.hpp"
// Cuts

#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglOddHole.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding.hpp"
#include "CglTwomir.hpp"

static int timesBad_=0;
//#############################################################################
// Solve methods
//#############################################################################
void CbcSolverLongThin::initialSolve()
{
  modelPtr_->scaling(0);
  setBasis(basis_,modelPtr_);
  modelPtr_->dual();
  basis_ = getBasis(modelPtr_);
  assert(!modelPtr_->specialOptions());
  modelPtr_->setLogLevel(0);
}

//-----------------------------------------------------------------------------
void CbcSolverLongThin::resolve()
{
  if (nestedSearch_<1.0&&model_) {
    // problem may be small enough to do nested search
    const double * colLower = modelPtr_->getColLower();
    const double * colUpper = modelPtr_->getColUpper();

    int numberIntegers = model_->numberIntegers();
    const int * integerVariable = model_->integerVariable();
  
    int i;
    int nFix=0;
    for (i=0;i<numberIntegers;i++) {
      int iColumn=integerVariable[i];
      if (colLower[iColumn]==colUpper[iColumn])
	nFix++;
    }
    if (nFix>nestedSearch_*numberIntegers) {
      // Do nested search
     int returnCode= model_->subBranchAndBound(colLower,colUpper,50000);
     if (returnCode==0||returnCode==2) {
       modelPtr_->setProblemStatus(1);
       return;
     }
    }
  }
  if (count_<100||justCount_) {
    assert(!modelPtr_->specialOptions());
    modelPtr_->setSpecialOptions(64+128+512);
    setBasis(basis_,modelPtr_);
    //modelPtr_->setLogLevel(1);
    modelPtr_->dual(0,0);
    basis_ = getBasis(modelPtr_);
    modelPtr_->setSpecialOptions(0);
    if (modelPtr_->status()==0) {
      count_++;
      double * solution = modelPtr_->primalColumnSolution();
      int i;
      int numberColumns = modelPtr_->numberColumns();
      for (i=0;i<numberColumns;i++) {
	if (solution[i]>1.0e-6||modelPtr_->getStatus(i)==ClpSimplex::basic) {
	  node_[i]=CoinMax(count_,node_[i]);
	  howMany_[i]++;
	}
      }
    } else {
      if (!justCount_)
	printf("infeasible early on\n");
    }
  } else {
    // use counts
    int numberRows=modelPtr_->numberRows();
    int numberColumns = modelPtr_->numberColumns();
    int * whichRow = new int[numberRows];
    int * whichColumn = new int [numberColumns];
    int i;
    for (i=0;i<numberRows;i++) 
      whichRow[i]=i;
    const double * lower = modelPtr_->columnLower();
    const double * upper = modelPtr_->columnUpper();
    int n=0;
    for (i=0;i<numberColumns;i++) {
      if ((node_[i]>count_-memory_&&node_[i]>0&&upper[i])
	  ||modelPtr_->getStatus(i)!=ClpSimplex::atLowerBound
	  ||lower[i]>0.0)
	whichColumn[n++]=i;
    }
    setBasis(basis_,modelPtr_);
    ClpSimplex * temp = new ClpSimplex(modelPtr_,numberRows,whichRow,n,whichColumn);
    delete [] whichRow;
    temp->setSpecialOptions(128+512);
    temp->setDualObjectiveLimit(1.0e50);
    temp->dual();
    if (temp->status()) {
      // In some cases we know that it must be infeasible
      if (believeInfeasible_) {
	modelPtr_->setProblemStatus(1);
	printf("assuming infeasible!\n");
	//modelPtr_->writeMps("infeas.mps");
	//temp->writeMps("infeas2.mps");
	//abort();
	delete temp;
	delete [] whichColumn;
	return;
      }
    }
    double * solution = modelPtr_->primalColumnSolution();
    if (!temp->status()) {
      const double * solution2 = temp->primalColumnSolution();
      for (i=0;i<n;i++) {
	int iColumn = whichColumn[i];
	solution[iColumn]=solution2[i];
	modelPtr_->setStatus(iColumn,temp->getStatus(i));
      }
      memcpy(modelPtr_->statusArray()+numberColumns,temp->statusArray()+n,
	     numberRows);
      memcpy(modelPtr_->primalRowSolution(),temp->primalRowSolution(),
	     numberRows*sizeof(double));
      double * dual = modelPtr_->dualRowSolution();
      memcpy(dual,temp->dualRowSolution(),
	     numberRows*sizeof(double));
      // See if optimal
      double * dj = modelPtr_->dualColumnSolution();
      // get reduced cost for large problem
      // this assumes minimization
      memcpy(dj,modelPtr_->objective(),numberColumns*sizeof(double));
      modelPtr_->transposeTimes(-1.0,dual,dj);
      modelPtr_->setObjectiveValue(temp->objectiveValue());
      modelPtr_->setProblemStatus(0);
      int nBad=0;
      
      for (i=0;i<numberColumns;i++) {
	if (modelPtr_->getStatus(i)==ClpSimplex::atLowerBound
	    &&upper[i]>lower[i]&&dj[i]<-1.0e-5)
	  nBad++;
      }
      if (nBad) {
	//printf("%d bad\n",nBad);
	timesBad_++;
	modelPtr_->primal();
      }
    } else {
      // infeasible - do all
      modelPtr_->setSpecialOptions(64+128+512);
      setBasis(basis_,modelPtr_);
      //modelPtr_->setLogLevel(1);
      modelPtr_->dual(0,0);
      basis_ = getBasis(modelPtr_);
      modelPtr_->setSpecialOptions(0);
      if (modelPtr_->status()) {
	printf("really infeasible!\n");
	delete temp;
	delete [] whichColumn;
	return;
      } else {
	printf("initially infeasible\n");
      }
    }
    delete temp;
    delete [] whichColumn;
    basis_ = getBasis(modelPtr_);
    modelPtr_->setSpecialOptions(0);
    count_++;
    if ((count_%100)==0)
      printf("count %d, bad %d\n",count_,timesBad_);
    for (i=0;i<numberColumns;i++) {
      if (solution[i]>1.0e-6||modelPtr_->getStatus(i)==ClpSimplex::basic) {
	node_[i]=CoinMax(count_,node_[i]);
	howMany_[i]++;
      }
    }
    if (modelPtr_->objectiveValue()>=modelPtr_->dualObjectiveLimit())
      modelPtr_->setProblemStatus(1);
  }
}

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
CbcSolverLongThin::CbcSolverLongThin ()
  : OsiClpSolverInterface()
{
  node_=NULL;
  howMany_=NULL;
  count_=0;
  model_ = NULL;
  memory_=300;
  believeInfeasible_=false;
  nestedSearch_ = 1.0;
  justCount_=false;
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
OsiSolverInterface * 
CbcSolverLongThin::clone(bool CopyData) const
{
  if (CopyData) {
    return new CbcSolverLongThin(*this);
  } else {
    printf("warning CbcSolveUser clone with copyData false\n");
    return new CbcSolverLongThin();
  }
}


//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
CbcSolverLongThin::CbcSolverLongThin (
                  const CbcSolverLongThin & rhs)
  : OsiClpSolverInterface(rhs)
{
  model_ = rhs.model_;
  int numberColumns = modelPtr_->numberColumns();
  node_=CoinCopyOfArray(rhs.node_,numberColumns);
  howMany_=CoinCopyOfArray(rhs.howMany_,numberColumns);
  count_=rhs.count_;
  memory_=rhs.memory_;
  believeInfeasible_ = rhs.believeInfeasible_;
  nestedSearch_ = rhs.nestedSearch_;
  justCount_=rhs.justCount_;
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
CbcSolverLongThin::~CbcSolverLongThin ()
{
  delete [] node_;
  delete [] howMany_;
}

//-------------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
CbcSolverLongThin &
CbcSolverLongThin::operator=(const CbcSolverLongThin& rhs)
{
  if (this != &rhs) { 
    OsiClpSolverInterface::operator=(rhs);
    delete [] node_;
    delete [] howMany_;
    model_ = rhs.model_;
    int numberColumns = modelPtr_->numberColumns();
    node_=CoinCopyOfArray(rhs.node_,numberColumns);
    howMany_=CoinCopyOfArray(rhs.howMany_,numberColumns);
    count_=rhs.count_;
    memory_=rhs.memory_;
    believeInfeasible_ = rhs.believeInfeasible_;
    nestedSearch_ = rhs.nestedSearch_;
    justCount_=rhs.justCount_;
  }
  return *this;
}
//-------------------------------------------------------------------
// Real initializer
//-------------------------------------------------------------------
void
CbcSolverLongThin::initialize (CbcModel * model, const char * keep)
{
  model_=model;
  int numberColumns = modelPtr_->numberColumns();
  if (numberColumns) {
    node_ = new int[numberColumns];
    howMany_ = new int[numberColumns];
    for (int i=0;i<numberColumns;i++) {
      if (keep[i])
	node_[i]=INT_MAX;
      else
	node_[i]=0;
      howMany_[i]=0;
    }
  } else {
    node_=NULL;
    howMany_=NULL;
  }
}
