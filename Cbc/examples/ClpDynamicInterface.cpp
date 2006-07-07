// Copyright (C) 2006, International Business Machines
// Corporation and others.  All Rights Reserved.

#include <cassert>

#include "CoinHelperFunctions.hpp"
//#include "CoinIndexedVector.hpp"
//#include "ClpSimplex.hpp"
#include "ClpDynamicInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcModel.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CbcStrategy.hpp"
#include "OsiOpbdpSolve.hpp"
//#############################################################################
// Solve methods
//#############################################################################
void ClpDynamicInterface::initialSolve()
{
  OsiClpSolverInterface::initialSolve();
}

//-----------------------------------------------------------------------------
void ClpDynamicInterface::resolve()
{
  if (cbcModel_) {
    ClpDynamicInterface * solver 
      = dynamic_cast<ClpDynamicInterface *> (cbcModel_->continuousSolver());
    if (solver==this) {
      // use parent
      OsiClpSolverInterface::resolve();
      return;
    }
  }

  printf("%d proposals\n",proposals_.getNumCols());
  // build up model
  ClpSimplex * model = new ClpSimplex(*staticModel_);
  // Included are costed slacks (numberBlocks_ of them)
  int numberCommon = staticModel_->numberColumns();
  // make sure artificials can come in
  int j;
  for (j=numberCommon-numberArtificials_;j<numberCommon;j++) {
    model->setColumnUpper(j,COIN_DBL_MAX);
  }
  numberProposalsAdded_=0;
  setBasis(dynamicBasis_,model);
  int numberColumns = modelPtr_->numberColumns();
  const double * columnLower = modelPtr_->columnLower();
  const double * columnUpper = modelPtr_->columnUpper();
  // make sure bounds are good
  int * backward = backward_[numberBlocks_];
  double * columnLower2 = model->columnLower();
  double * columnUpper2 = model->columnUpper();
  int iColumn;
  for (iColumn=0;iColumn<numberCommon-numberArtificials_;iColumn++) {
    int jColumn = backward[iColumn];
    columnLower2[iColumn]=columnLower[jColumn];
    columnUpper2[iColumn]=columnUpper[jColumn];
  }
  double * solution = modelPtr_->primalColumnSolution();
  int numberProposals=proposals_.getNumCols();
  model->setLogLevel(0);
  // here solve and generate
  // save cutoff
  double cutoff = modelPtr_->dualObjectiveLimit();
  model->setDualObjectiveLimit(1.0e50);
  model->dual();
  // Big trouble if not feasible
  assert (!model->problemStatus());
  // add in existing valid ones
  {
    int * which = new int[numberProposals];
    for (int i=0;i<numberProposals;i++)
      which[i]=i;
    addProposals(model,numberProposals,which,false);
    delete [] which;
  }
  model->primal();
  // generate
  double lastObjective = model->objectiveValue();
  double sumDj=0.0;
  #define NPASS 5
  #define MAXPASS 100
  int maxPass=MAXPASS; // but no more than this
  for (int iPass=0;iPass<NPASS;iPass++) {
    maxPass--;
    double * saveObjective = new double[numberColumns];
    int numberRows2 = staticModel_->numberRows();
    int numberMasterRows = numberRows2-numberBlocks_;
    // we need big arrays for storing proposals
    int size = modelPtr_->numberColumns()+numberRows2+3;
    int * rows = new int[size];
    double * elements = new double[size];
    sumDj=0.0;
    double * pi = new double[numberRows2];
    for (int iBlock=0;iBlock<numberBlocks_;iBlock++) {
      CoinMemcpyN(model->dualRowSolution(),numberRows2,pi);;
      double * solution2 = model->primalColumnSolution();
      int numberBad=0;
      int i;
      for ( i=numberCommon-numberArtificials_;i<numberCommon;i++) {
        if (fabs(solution2[i])>1.0e-8) {
          numberBad++;
        }
      }
      int numberColumns2 = subProblem_[iBlock].numberColumns();
      double * objective = subProblem_[iBlock].objective();
      CoinMemcpyN(objective,numberColumns2,saveObjective);
      masterRow_[iBlock].transposeTimes(pi,objective);
      // and correct objective
      for ( i=0;i<numberColumns2;i++) {
        if (!numberBad)
          objective[i] = saveObjective[i] - objective[i];
        else
          objective[i] = - objective[i];
      }
      // Now solve IP by some means
      OsiClpSolverInterface solver (&subProblem_[iBlock],false);
      // fix variables and clean objective
      const int * backward = backward_[iBlock];
      int nInt=0;
      double * obj = solver.getModelPtr()->objective();
      double maxObj=0.0;
      for ( i=0;i<numberColumns2;i++) {
        maxObj = CoinMax(maxObj,fabs(obj[i]));
        int iColumn = backward[i];
        if (solver.isInteger(i)&&!columnLower[i]&&columnUpper[i]==1.0)
          nInt++;
        solver.setColLower(i,columnLower[iColumn]);
        solver.setColUpper(i,columnUpper[iColumn]);
      }
      if (maxObj>1.0e7) {
        double scaleFactor = 1.0e7/maxObj;
        for ( i=0;i<numberColumns2;i++) 
          obj[i] *= scaleFactor;
      }
      if (nInt==numberColumns2&&!iPass&&!cbcModel_->getNodeCount())
        printf("pure 0-1 problem\n");
#define ENUMERATE
#ifndef ENUMERATE
      CbcModel small(solver);
      // Normal type strategy
      CbcStrategyDefault strategy(true,5,5);
      small.setStrategy(strategy);
      small.setLogLevel(0);
      small.branchAndBound();
      if (small.getMinimizationObjValue()<1.0e50) {
        double * solution = small.bestSolution();
        // clean
        for (i=0;i<numberColumns2;i++) {
          if (solver.isInteger(i))
            solution[i]=floor(solution[i]+0.5);
        }
#else
        int rcode=solveOpbdp(&solver);
        if (rcode>0) {
          const double * solution = solver.getColSolution();
#endif
	double objValue =0.0;
	// proposal
        masterRow_[iBlock].times(solution,elements);
        for (i=0;i<numberColumns2;i++) {
          objValue += solution[i]*saveObjective[i];
        }
        // See if good dj and pack down
        int number=0;
        double dj = objValue;
        double smallest=1.0e100;
        double largest=0.0;
        for (i=0;i<numberMasterRows;i++) {
          double value = elements[i];
          if (fabs(value)>1.0e-15) {
            dj -= pi[i]*value;
            smallest = min(smallest,fabs(value));
            largest = max(largest,fabs(value));
            rows[number]=i;
            elements[number++]=value;
          }
        }
        // and convexity
        dj -= pi[numberMasterRows+iBlock];
        rows[number]=numberMasterRows+iBlock;
        elements[number++]=1.0;
        // if elements large then scale?
        //if (largest>1.0e8||smallest<1.0e-8)
        if (dj<0.0)
          sumDj += dj;
        if (dj<-1.0e-6) {
          printf("For subproblem %d smallest - %g, largest %g - dj %g\n",
                 iBlock,smallest,largest,dj);
          // take
          // objective
          int row = numberMasterRows+numberBlocks_;
          rows[number]=row;
          elements[number++]=objValue;
          // what went into it
          row++;
          for (i=0;i<numberColumns2;i++) {
            if (fabs(solution[i])>1.0e-15) {
              rows[number]=row+i;
              elements[number++]=solution[i];
            }
          }
          // add
          proposals_.appendCol(number,rows,elements);
          const double * solution2 = model->primalColumnSolution();
          // See if any artificials in
          bool artificialsIn=false;
          for (int j=numberCommon-numberArtificials_;j<numberCommon;j++) {
            if (fabs(solution2[j])>1.0e-8||model->getStatus(j)==ClpSimplex::basic) {
              artificialsIn=true;
              break;
            }
          }
          if ((maxPass>MAXPASS-5||artificialsIn)&&!cbcModel_->getNodeCount()) {
            // add in one at a time in case of symmetry
            int n=proposals_.getNumCols();
            int * which = new int[n-numberProposals];
            for (int i=numberProposals;i<n;i++)
              which[i-numberProposals]=i;
            addProposals(model,n-numberProposals,which,false);
            delete [] which;
            numberProposals=n;
            //model->setLogLevel(1);
            model->primal(1);
            assert (!model->problemStatus());
            double * solution2 = model->primalColumnSolution();
            // Take out any artificials if we can
            for (int j=numberCommon-numberArtificials_;j<numberCommon;j++) {
              if (fabs(solution2[j])<1.0e-8&&model->getStatus(j)==ClpSimplex::basic) {
                solution2[j]=0.0;
                model->setStatus(j,ClpSimplex::atLowerBound);
                model->setColumnUpper(j,0.0);
              }
            }
          }
        }
      }
      CoinMemcpyN(saveObjective,numberColumns2,objective);
    }
    delete [] pi;
    delete [] rows;
    delete [] elements;
    delete [] saveObjective;
    if ((sumDj>-1.0e-5||iPass==NPASS-1)) {
      // exit
      break;
    } else {
      int n=proposals_.getNumCols();
      int * which = new int[n-numberProposals];
      for (int i=numberProposals;i<n;i++)
        which[i-numberProposals]=i;
      addProposals(model,n-numberProposals,which,false);
      delete [] which;
      numberProposals=n;
      model->setLogLevel(1);
      // Take out any artificials if we can
      double * solution2 = model->primalColumnSolution();
      for (int j=numberCommon-numberArtificials_;j<numberCommon;j++) {
        if (fabs(solution2[j])<1.0e-8&&model->getStatus(j)==ClpSimplex::basic) {
          solution2[j]=0.0;
          model->setStatus(j,ClpSimplex::atLowerBound);
          model->setColumnUpper(j,0.0);
        }
      }
      model->primal(1);
      // Check if we have valid solution
      int status=setSolution(model);
      if (status<0)
        printf("*** we have a solution of %g\n",model->objectiveValue());
      if (lastObjective<model->objectiveValue()+1.0e-5&&cbcModel_->getNodeCount())
        iPass=NPASS-2; // just generate and exit
      else if (!cbcModel_->getNodeCount()&&status==1)
        iPass=0; // not feasible - make sure we keep going
      else if (!cbcModel_->getNodeCount()&&(sumDj<-1.0e-3&&maxPass>0))
        iPass=0; // not optimal - make sure we keep going
      lastObjective = model->objectiveValue();
    }
  }
  printf("DW Objective value is %g - best possible %g\n",model->objectiveValue(),
         model->objectiveValue()+sumDj);
  dynamicBasis_ = getBasis(model,proposalsAdded_,numberCommon);
  int numberRows = modelPtr_->numberRows();
  int problemStatus = model->problemStatus();
  CoinAssert (model->problemStatus()||model->objectiveValue()<1.0e50);
  if (model->objectiveValue()+sumDj>cutoff) {
    // say infeasible
    problemStatus=1;
  }
  modelPtr_->setObjectiveValue(model->objectiveValue()+sumDj);
  modelPtr_->setSumDualInfeasibilities(model->sumDualInfeasibilities());
  modelPtr_->setNumberDualInfeasibilities(model->numberDualInfeasibilities());
  modelPtr_->setSumPrimalInfeasibilities(model->sumPrimalInfeasibilities());
  modelPtr_->setNumberPrimalInfeasibilities(model->numberPrimalInfeasibilities());
  if (!problemStatus) {
    problemStatus=setSolution(model);
    if (problemStatus<0) {
      printf("*** We have a valid solution of %g\n",model->objectiveValue());
      problemStatus=0;
    }
  }
  //  If at root node try BAB
  if (!cbcModel_->getNodeCount()) {
    // build up model
    ClpSimplex * model = new ClpSimplex(*staticModel_);
    // Included are costed slacks (numberBlocks_ of them)
    int numberCommon = staticModel_->numberColumns();
    numberProposalsAdded_=0;
    const double * columnLower = modelPtr_->columnLower();
    const double * columnUpper = modelPtr_->columnUpper();
    // make sure bounds are good
    int * backward = backward_[numberBlocks_];
    double * columnLower2 = model->columnLower();
    double * columnUpper2 = model->columnUpper();
    int iColumn;
    for (iColumn=0;iColumn<numberCommon-numberArtificials_;iColumn++) {
      int jColumn = backward[iColumn];
      columnLower2[iColumn]=columnLower[jColumn];
      columnUpper2[iColumn]=columnUpper[jColumn];
    }
    int numberProposals=proposals_.getNumCols();
    // add in existing valid ones
    {
      int * which = new int[numberProposals];
      int i;
      for ( i=0;i<numberProposals;i++)
        which[i]=i;
      addProposals(model,numberProposals,which,false);
      // make integer
      for ( i=numberCommon;i<model->numberColumns();i++) {
        model->setColLower(i,0.0);
        model->setColUpper(i,1.0);
        model->setInteger(i);
      }
      delete [] which;
    }
    OsiClpSolverInterface * clpSolver1 = new OsiClpSolverInterface(model,false);
    CbcModel small(*clpSolver1);
    clpSolver1->releaseClp();
    delete clpSolver1;
    small.setLogLevel(0);
    CbcStrategyDefault strategy(true,5,5);
    small.setStrategy(strategy);
    small.branchAndBound();
    if (small.getMinimizationObjValue()<1.0e50) {
      // save solution round this
      double * saveSolution = CoinCopyOfArray(modelPtr_->primalColumnSolution(),
                                         modelPtr_->numberColumns());
      //const double * solution = small.bestSolution();
      OsiSolverInterface * solver1 = small.solver();
      OsiClpSolverInterface * clpSolver 
        = dynamic_cast<OsiClpSolverInterface *> (solver1);
      assert (clpSolver);
      ClpSimplex * clpSimplex = clpSolver->getModelPtr();
      int status=setSolution(clpSimplex);
      if (status<0)
        printf("*** we have a solution of %g\n",clpSimplex->objectiveValue());
      // copy back solution
      CoinMemcpyN(saveSolution,
                  modelPtr_->numberColumns(),
                  modelPtr_->primalColumnSolution());
      delete [] saveSolution ;
    }
    delete model;
  }
  modelPtr_->setProblemStatus(problemStatus);
  // Don't allow reduced cost fixing?
  CoinZeroN(modelPtr_->dualColumnSolution(),numberColumns);
  memset(modelPtr_->primalRowSolution(),0,numberRows*sizeof(double));
  modelPtr_->clpMatrix()->times(1.0,solution,modelPtr_->primalRowSolution());
  modelPtr_->setNumberIterations(model->numberIterations());
  delete model;
}

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
ClpDynamicInterface::ClpDynamicInterface ()
  : OsiClpSolverInterface()
{
  numberBlocks_=0;
  numberArtificials_ = 0;
  staticModel_=NULL;
  subProblem_=NULL;
  masterRow_=NULL;
  numberProposalsAdded_ = 0;
  maximumProposalsAdded_ =0;
  proposalsAdded_ = NULL;
  backward_=NULL;
  allGenerated_ = NULL;
  cbcModel_ = NULL;
  bestObjectiveValue_ =1.0e100;
  bestSolution_=NULL;
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
OsiSolverInterface * 
ClpDynamicInterface::clone(bool CopyData) const
{
  if (CopyData) {
    return new ClpDynamicInterface(*this);
  } else {
    printf("warning ClpDynamicInterface clone with copyData false\n");
    return new ClpDynamicInterface();
  }
}


//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
ClpDynamicInterface::ClpDynamicInterface (
                  const ClpDynamicInterface & rhs)
  : OsiClpSolverInterface(rhs)
{
  dynamicBasis_ = rhs.dynamicBasis_;
  numberBlocks_ = rhs.numberBlocks_;
  numberArtificials_ = rhs.numberArtificials_;
  numberProposalsAdded_ = rhs.numberProposalsAdded_;
  maximumProposalsAdded_ =rhs.maximumProposalsAdded_;
  bestObjectiveValue_ = rhs.bestObjectiveValue_;
  cbcModel_ = rhs.cbcModel_;
  if (rhs.staticModel_) {
    staticModel_=new ClpSimplex(*rhs.staticModel_);
    subProblem_=new ClpSimplex[numberBlocks_];
    masterRow_ = new CoinPackedMatrix[numberBlocks_];
    if (numberBlocks_) {
      backward_ = new int * [numberBlocks_+1];
      backward_[numberBlocks_]=CoinCopyOfArray(rhs.backward_[numberBlocks_],
                                               staticModel_->numberColumns()-numberArtificials_);
      for (int i=0;i<numberBlocks_;i++) {
        subProblem_[i]=rhs.subProblem_[i];
        masterRow_[i]=rhs.masterRow_[i];
        backward_[i]=CoinCopyOfArray(rhs.backward_[i],subProblem_[i].numberColumns());
      }
    } else {
      backward_=NULL;
    }
    proposalsAdded_ = CoinCopyOfArray(rhs.proposalsAdded_,maximumProposalsAdded_);;
    proposals_=rhs.proposals_;
    allGenerated_ = CoinCopyOfArray(rhs.allGenerated_,numberBlocks_);
    if (rhs.bestSolution_) {
      bestSolution_ = CoinCopyOfArray(rhs.bestSolution_,modelPtr_->getNumCols());
    } else {
      bestSolution_=NULL;
    }
  } else {
    staticModel_=NULL;
    subProblem_ = NULL;
    masterRow_=NULL;
    proposalsAdded_=NULL;
    backward_=NULL;
    allGenerated_ = NULL;
    bestSolution_=NULL;
  }
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
ClpDynamicInterface::~ClpDynamicInterface ()
{
  delete staticModel_;
  delete [] subProblem_;
  delete [] masterRow_;
  if (numberBlocks_) {
    for (int i=0;i<numberBlocks_+1;i++) {
      delete [] backward_[i];
    }
  }
  delete [] proposalsAdded_;
  delete [] backward_;
  delete [] allGenerated_;
  delete [] bestSolution_;
}

//-------------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
ClpDynamicInterface &
ClpDynamicInterface::operator=(const ClpDynamicInterface& rhs)
{
  if (this != &rhs) { 
    OsiClpSolverInterface::operator=(rhs);
    numberBlocks_ = rhs.numberBlocks_;
    numberArtificials_ = rhs.numberArtificials_;
    numberProposalsAdded_ = rhs.numberProposalsAdded_;
    maximumProposalsAdded_ =rhs.maximumProposalsAdded_;
    bestObjectiveValue_ = rhs.bestObjectiveValue_;
    delete staticModel_;
    delete [] subProblem_;
    delete [] masterRow_;
    if (numberBlocks_) {
      for (int i=0;i<numberBlocks_+1;i++) {
        delete [] backward_[i];
      }
    }
    delete [] proposalsAdded_;
    delete [] backward_;
    delete [] allGenerated_;
    delete [] bestSolution_;
    dynamicBasis_ = rhs.dynamicBasis_;
    proposalsAdded_ = CoinCopyOfArray(rhs.proposalsAdded_,maximumProposalsAdded_);;
    cbcModel_ = rhs.cbcModel_;
    if (rhs.staticModel_) {
      staticModel_=new ClpSimplex(*rhs.staticModel_);
      subProblem_=new ClpSimplex[numberBlocks_];
      masterRow_ = new CoinPackedMatrix[numberBlocks_];
      if (numberBlocks_) {
        backward_ = new int * [numberBlocks_+1];
        backward_[numberBlocks_]=CoinCopyOfArray(rhs.backward_[numberBlocks_],
                                                 staticModel_->numberColumns()-numberArtificials_);
        for (int i=0;i<numberBlocks_;i++) {
          subProblem_[i]=rhs.subProblem_[i];
          masterRow_[i]=rhs.masterRow_[i];
          backward_[i]=CoinCopyOfArray(rhs.backward_[i],subProblem_[i].numberColumns());
        }
      } else {
        backward_=NULL;
      }
      proposals_=rhs.proposals_;
      allGenerated_ = CoinCopyOfArray(rhs.allGenerated_,numberBlocks_);
      if (rhs.bestSolution_) {
        bestSolution_ = CoinCopyOfArray(rhs.bestSolution_,modelPtr_->getNumCols());
      } else {
        bestSolution_=NULL;
      }
    } else {
      staticModel_=NULL;
      subProblem_ = NULL;
      masterRow_=NULL;
      backward_=NULL;
      proposals_=CoinPackedMatrix();
      allGenerated_ = NULL;
      bestSolution_ = NULL;
    }
  }
  return *this;
}
//-------------------------------------------------------------------
// Real initializer
//-------------------------------------------------------------------
void
ClpDynamicInterface::initialize (int * rowBlock, CbcModel * cbcModel)
{
  // 
  delete staticModel_;
  delete [] subProblem_;
  delete [] masterRow_;
  cbcModel_ = cbcModel;
  int numberColumns = modelPtr_->numberColumns();
  int numberRows = modelPtr_->numberRows();
  CoinPackedMatrix * matrix = modelPtr_->matrix();
  int numberElements = matrix->getNumElements(); 
  // get row copy
  CoinPackedMatrix rowCopy = *matrix;
  rowCopy.reverseOrdering();
  const int * rowByColumn = matrix->getIndices();
  const int * columnLength = matrix->getVectorLengths();
  const CoinBigIndex * columnStart = matrix->getVectorStarts();
  //const double * elementByColumn = matrix->getElements();
  const int * column = rowCopy.getIndices();
  const int * rowLength = rowCopy.getVectorLengths();
  const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
  //const double * elementByRow = rowCopy.getElements();
  numberBlocks_ = 0;
  int * stack = new int [numberRows];
  // to say if column looked at
  int * columnBlock = new int[numberColumns];
  int iRow,iColumn;
  for (iColumn=0;iColumn<numberColumns;iColumn++)
    columnBlock[iColumn]=-2;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    int kstart = columnStart[iColumn];
    int kend = columnStart[iColumn]+columnLength[iColumn];
    if (columnBlock[iColumn]==-2) {
      // column not allocated
      int j;
      int nstack=0;
      for (j=kstart;j<kend;j++) {
	int iRow= rowByColumn[j];
	if (rowBlock[iRow]!=-1) {
	  assert(rowBlock[iRow]==-2);
	  rowBlock[iRow]=numberBlocks_; // mark
	  stack[nstack++] = iRow;
	}
      }
      if (nstack) {
	// new block - put all connected in
	numberBlocks_++;
	columnBlock[iColumn]=numberBlocks_-1;
	while (nstack) {
	  int iRow = stack[--nstack];
	  int k;
	  for (k=rowStart[iRow];k<rowStart[iRow]+rowLength[iRow];k++) {
	    int iColumn = column[k];
	    int kkstart = columnStart[iColumn];
	    int kkend = kkstart + columnLength[iColumn];
	    if (columnBlock[iColumn]==-2) {
	      columnBlock[iColumn]=numberBlocks_-1; // mark
	      // column not allocated
	      int jj;
	      for (jj=kkstart;jj<kkend;jj++) {
		int jRow= rowByColumn[jj];
		if (rowBlock[jRow]==-2) {
		  rowBlock[jRow]=numberBlocks_-1;
		  stack[nstack++]=jRow;
		}
	      }
	    } else {
	      assert (columnBlock[iColumn]==numberBlocks_-1);
	    }
	  }
	}
      } else {
	// Only in master
	columnBlock[iColumn]=-1;
      }
    }
  }
  printf("%d blocks found\n",numberBlocks_);
  delete [] stack;
  int i;
  subProblem_ = new ClpSimplex[numberBlocks_];
  masterRow_ = new CoinPackedMatrix[numberBlocks_];
  int * which = new int[numberColumns];
  int * whichRow = new int[numberRows];
  int * whichMaster = new int[numberRows];
  CoinBigIndex * starts = new int [numberColumns+1];
  int * row = new int[numberElements];
  double * element = new double[numberElements];
  // Master
  int numberRows2;
  int numberColumns2;
  numberRows2=0;
  numberColumns2=0;
  for (i=0;i<numberRows;i++) {
    if (rowBlock[i]==-1)
      whichMaster[numberRows2++]=i;
  }
  // and do backward pointers
  backward_ = new int * [numberBlocks_+1];
  int numberMasterRows = numberRows2;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (columnBlock[iColumn]==-1)
      which[numberColumns2++]=iColumn;
  }
  backward_[numberBlocks_]=CoinCopyOfArray(which,numberColumns2);
  staticModel_ = new ClpSimplex(modelPtr_,numberRows2,whichMaster,
                                numberColumns2,which,true,false,false);
  printf("Master has %d rows and %d columns\n",numberRows2,numberColumns2);
  // Add convexity rows
  for (i=0;i<numberBlocks_;i++) {
    staticModel_->addRow(0,NULL,NULL,1.0,1.0);
  }
  // Add costed slacks to make sure feasible
  // *** We could do better than this
  const double * rowLower = staticModel_->rowLower();
  const double * rowUpper = staticModel_->rowUpper();
  numberArtificials_ = 0;
  for (i=0;i<numberRows2;i++) {
    if (rowLower[i]>-1.0e20) {
      numberArtificials_ ++;
      double value = 1.0;
      staticModel_->addColumn(1,&i,&value,0.0,COIN_DBL_MAX,1.0e10);
    }
    if (rowUpper[i]<1.0e20) {
      numberArtificials_ ++;
      double value = -1.0;
      staticModel_->addColumn(1,&i,&value,0.0,COIN_DBL_MAX,1.0e10);
    }
  }
  for (i=0;i<numberBlocks_;i++) {
    int row=i+numberMasterRows;
    double value=1.0;
    staticModel_->addColumn(1,&row,&value,0.0,COIN_DBL_MAX,1.0e10);
  }
  numberArtificials_ += numberBlocks_;
  // Blocks
  for (i=0;i<numberBlocks_;i++) {
    numberRows2=0;
    numberColumns2=0;
    for (iRow=0;iRow<numberRows;iRow++) {
      if (rowBlock[iRow]==i)
        whichRow[numberRows2++]=iRow;
    }
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
    if (columnBlock[iColumn]==i)
      which[numberColumns2++]=iColumn;
    }
    backward_[i]=CoinCopyOfArray(which,numberColumns2);
    subProblem_[i] = ClpSimplex(modelPtr_,numberRows2,whichRow,
                                numberColumns2,which,true,false,false);
    printf("Subproblem %d has %d rows and %d columns\n",
           i,numberRows2,numberColumns2);
    masterRow_[i] = CoinPackedMatrix(*matrix,numberMasterRows,whichMaster,
                                     numberColumns2,which);
  }
  delete [] starts;
  delete [] row;
  delete [] element;
  delete [] which;
  delete [] whichRow;
  delete [] whichMaster;
  delete [] columnBlock;
  // Could see if we can generate all upfront
  allGenerated_ = new int[numberBlocks_];
  CoinZeroN(allGenerated_,numberBlocks_);
}
// Warm start
CoinWarmStartBasisDynamic
ClpDynamicInterface::getBasis(ClpSimplex * model, const int * whichColumns,
                              int numberCommon) const
{
  int iRow,iColumn;
  int numberRows = model->numberRows();
  int numberColumns = model->numberColumns();
  CoinWarmStartBasisDynamic basis;
  int * dynamic = new int [numberColumns];
  // compute number of columns
  int numberColumns2=numberCommon;
  const double * columnLower = model->columnLower();
  for (iColumn=numberCommon;iColumn<numberColumns;iColumn++) {
    int iStatus = model->getColumnStatus(iColumn);
    if (columnLower[iColumn]||iStatus==1) 
      numberColumns2++;
  }
  basis.setSize(numberColumns2,numberRows);
  if (model->statusExists()) {
    // Flip slacks
    int lookupA[]={0,1,3,2,0,2};
    for (iRow=0;iRow<numberRows;iRow++) {
      int iStatus = model->getRowStatus(iRow);
      iStatus = lookupA[iStatus];
      basis.setArtifStatus(iRow,(CoinWarmStartBasis::Status) iStatus);
    }
    int numberColumns2=0;
    int numberDynamic=0;
    const double * columnLower = model->columnLower();
    int lookupS[]={0,1,2,3,0,3};
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      int iStatus = model->getColumnStatus(iColumn);
      if (columnLower[iColumn]||iStatus==1||iColumn<numberCommon) {
        iStatus = lookupS[iStatus];
        basis.setStructStatus(numberColumns2,(CoinWarmStartBasis::Status) iStatus);
        if (iColumn>=numberCommon) 
          dynamic[numberDynamic++]=whichColumns[iColumn-numberCommon];
        numberColumns2++;
      }
    }
    basis.setDynamicVariables(numberDynamic,dynamic);
    delete [] dynamic;
    basis.setNumCommonVariables(numberCommon);
  }
  return basis;
}
// Sets up basis
void 
ClpDynamicInterface::setBasis ( const CoinWarmStartBasisDynamic & basis,
				  ClpSimplex * model)
{
  // Return if null basis
  if (!dynamicBasis_.getNumStructural())
    return;
  int numberCommon = basis.getNumCommonVariables();
  assert (numberCommon==model->numberColumns());
  int numberDynamic = basis.getNumDynamicVariables();
  const int * dynamic = basis.getDynamicVariables();
  // add in proposals
  addProposals(model, numberDynamic, dynamic,true);
  // Do rest of basis stuff
  OsiClpSolverInterface::setBasis(basis,model);
}
// Adds proposals to model
void 
ClpDynamicInterface::addProposals(ClpSimplex * model, int number,
                                  const int * which, bool addEvenIfFixed)
{
  const double * columnLower = modelPtr_->columnLower();
  const double * columnUpper = modelPtr_->columnUpper();
  double * cost = new double [number];
  double * lower = new double [number];
  double * upper = new double [number];
  int numberProposals = proposals_.getNumCols();
  // mark ones in
  char * mark = new char[numberProposals];
  memset(mark,0,numberProposals);
  int i;
  // re-allocate added if necessary
  if (number+numberProposalsAdded_>maximumProposalsAdded_) {
    int * temp = proposalsAdded_;
    maximumProposalsAdded_ = 
      CoinMax(2*maximumProposalsAdded_+100,number+numberProposalsAdded_);
    proposalsAdded_ = new int[maximumProposalsAdded_];
    CoinZeroN(proposalsAdded_+numberProposalsAdded_,
              maximumProposalsAdded_-numberProposalsAdded_);
    CoinMemcpyN(temp,numberProposalsAdded_,proposalsAdded_);
    delete [] temp;
  }
  for (i=0;i<numberProposalsAdded_;i++) 
    mark[proposalsAdded_[i]]=1;
  int numberRows2 = model->numberRows();
  int numberMasterRows = numberRows2-numberBlocks_;
  int objRow=numberRows2;
  int jColumn;
  // Get counts of variables fixed to 1
  int * numberAtOne = new int[numberBlocks_];
  for (int i=0;i<numberBlocks_;i++) {
    int * backward = backward_[i];
    int n=subProblem_[i].numberColumns();
    int count=0;
    for (int j=0;j<n;j++) {
      int iColumn = backward[j];
      if (columnLower[iColumn])
        count++;
    }
    numberAtOne[i]=count;
  }
  int numberElements = proposals_.getNumElements();
  CoinBigIndex * starts = new int [number+1];
  int * row = new int[numberElements];
  double * element = new double[numberElements];
  const int * rowByColumn = proposals_.getIndices();
  const int * columnLength = proposals_.getVectorLengths();
  const CoinBigIndex * columnStart = proposals_.getVectorStarts();
  const double * elementByColumn = proposals_.getElements();
  starts[0]=0;
  CoinBigIndex numberElements2=0;
  int numberColumns2=0;
  for (jColumn=0;jColumn<number;jColumn++) {
    int iColumn = which[jColumn];
    // skip if already in
    if (mark[iColumn])
      continue;
    double costValue=0.0;
    int iBlock=-1;
    int * backward = NULL;
    bool out=false;
    int numberOne=0;
    for (CoinBigIndex j=columnStart[iColumn];
         j<columnStart[iColumn]+columnLength[iColumn];j++) {
      int iRow = rowByColumn[j];
      if (iRow<objRow) {
        row[numberElements2]=iRow;
        element[numberElements2++]=elementByColumn[j];
        if (iRow>=numberMasterRows) {
          iBlock = iRow-numberMasterRows;
          backward = backward_[iBlock];
        }
      } else if (iRow==objRow) {
        costValue = elementByColumn[j];
      } else {
        assert(backward);
        if (iBlock==8&&proposals_.getNumCols()==1493)
          printf("block 8 row %d\n",iRow);
        int kColumn =backward[iRow-objRow-1];
        // only works for 0-1 at present
        if (!columnUpper[kColumn]) {
          out=true;
          break;
        } else if (columnLower[kColumn]) {
          numberOne++;
        }
      }
    }
    // If does not match state of problem - discard
    if (numberOne!=numberAtOne[iBlock])
      out=true;
    cost[numberColumns2]=costValue;
    lower[numberColumns2]=0.0;
    upper[numberColumns2]= out ? 0.0 : COIN_DBL_MAX;
    starts[numberColumns2+1]=numberElements2;
    // see if we take
    if (!out||addEvenIfFixed) {
      numberColumns2++;
      mark[iColumn]=1;
      proposalsAdded_[numberProposalsAdded_++]=iColumn;
    } else {
      // take off elements
      numberElements2=starts[numberColumns2];
    }
  }
  model->addColumns(numberColumns2,lower,upper,cost,
                   starts,row,element);
  delete [] starts;
  delete [] row;
  delete [] element;
  delete [] cost;
  delete [] lower;
  delete [] upper;
  delete [] numberAtOne;
  delete [] mark;
}
/* Creates modelPtr_ solution.
   Returns 1 if not feasible, -1 if integer solution
   0 otherwise */
int 
ClpDynamicInterface::setSolution(ClpSimplex * model)
{
  int numberColumns = modelPtr_->numberColumns();
  const double * columnLower = modelPtr_->columnLower();
  const double * columnUpper = modelPtr_->columnUpper();
  double * solution = modelPtr_->primalColumnSolution();
  int iColumn;
  for (iColumn=0;iColumn<numberColumns;iColumn++) {
    // make solution sensible
    solution[iColumn]=0.0;
  }
  const double * solution2 = model->primalColumnSolution();
  int problemStatus=0;
  int numberCommon = staticModel_->numberColumns();
  // copy back common (without costed slacks)
  int * backward = backward_[numberBlocks_];
  for (iColumn=0;iColumn<numberCommon-numberArtificials_;iColumn++) {
    solution[backward[iColumn]]= solution2[iColumn];
  }
  // infeasible if artificials in
  for (;iColumn<numberCommon;iColumn++) {
    if (fabs(solution2[iColumn])>1.0e-5)
      problemStatus=1;
  }
  const int * rowByColumn = proposals_.getIndices();
  const int * columnLength = proposals_.getVectorLengths();
  const CoinBigIndex * columnStart = proposals_.getVectorStarts();
  const double * elementByColumn = proposals_.getElements();
  int numberMasterRows = model->numberRows()-numberBlocks_;
  int objRow=model->numberRows();
  int numberColumns2=model->numberColumns()-numberCommon;
  for (int jColumn=0;jColumn<numberColumns2;jColumn++) {
    int iColumn = proposalsAdded_[jColumn];
    int iBlock=-1;
    int * backward = NULL;
    double value = solution2[jColumn+numberCommon];
    if (value) {
      for (CoinBigIndex j=columnStart[iColumn];
           j<columnStart[iColumn]+columnLength[iColumn];j++) {
        int iRow = rowByColumn[j];
        if (iRow<objRow&&iRow>=numberMasterRows) {
          iBlock = iRow-numberMasterRows;
          backward = backward_[iBlock];
        } else if (iRow>objRow) {
          assert(backward);
          int kColumn =backward[iRow-objRow-1];
          solution[kColumn] += value*elementByColumn[j];
        }
      }
    }
  }
  if (problemStatus==0) {
    problemStatus=-1;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      assert (solution[iColumn]>=columnLower[iColumn]-1.0e-7);
      assert (solution[iColumn]<=columnUpper[iColumn]+1.0e-7);
      if (modelPtr_->isInteger(iColumn)&&fabs(solution[iColumn]-floor(solution[iColumn]+0.5))>1.0e-7)
        problemStatus=0;
    }
    if (problemStatus&&model->objectiveValue()<bestObjectiveValue_) {
      delete [] bestSolution_;
      bestSolution_ = CoinCopyOfArray(solution,modelPtr_->getNumCols());
      bestObjectiveValue_ = model->objectiveValue();
      int numberRows = modelPtr_->numberRows();
      double * rowActivity = modelPtr_->primalRowSolution();
      CoinZeroN(rowActivity,numberRows);
      modelPtr_->times(1.0,bestSolution_,rowActivity);
      const double * rowLower = modelPtr_->rowLower();
      const double * rowUpper = modelPtr_->rowUpper();
      bool good=true;
      for (int i=0;i<numberRows;i++) {
        if(rowActivity[i]<rowLower[i]-1.0e-5) {
          printf("%d %g %g %g\n",i,rowLower[i],rowActivity[i],rowUpper[i]);
          good=false;
        } else if(rowActivity[i]>rowUpper[i]+1.0e-5) {
          printf("%d %g %g %g\n",i,rowLower[i],rowActivity[i],rowUpper[i]);
          good=false;
        }
      }
      assert (good);
    }
  }
  return problemStatus;
}
// Default Constructor
CbcHeuristicDynamic::CbcHeuristicDynamic() 
  :CbcHeuristic()
{
}

// Constructor from model
CbcHeuristicDynamic::CbcHeuristicDynamic(CbcModel & model)
  :CbcHeuristic(model)
{
}

// Destructor 
CbcHeuristicDynamic::~CbcHeuristicDynamic ()
{
}

// Clone
CbcHeuristic *
CbcHeuristicDynamic::clone() const
{
  return new CbcHeuristicDynamic(*this);
}

// Copy constructor 
CbcHeuristicDynamic::CbcHeuristicDynamic(const CbcHeuristicDynamic & rhs)
:
  CbcHeuristic(rhs)
{
}

// Returns 1 if solution, 0 if not
int
CbcHeuristicDynamic::solution(double & solutionValue,
			 double * betterSolution)
{
  if (!model_)
    return 0;
  ClpDynamicInterface * clpSolver 
    = dynamic_cast<ClpDynamicInterface *> (model_->solver());
  assert (clpSolver); 
  double newSolutionValue = clpSolver->bestObjectiveValue();
  const double * solution = clpSolver->bestSolution();
  if (newSolutionValue<solutionValue&&solution) {
    int numberColumns = clpSolver->getNumCols();
    // new solution
    memcpy(betterSolution,solution,numberColumns*sizeof(double));
    solutionValue = newSolutionValue;
    return 1;
  } else {
    return 0;
  }
}
// update model
void CbcHeuristicDynamic::setModel(CbcModel * model)
{
  model_ = model;
}
// Resets stuff if model changes
void 
CbcHeuristicDynamic::resetModel(CbcModel * model)
{
  model_ = model;
}
