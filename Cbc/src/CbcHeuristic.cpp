// Copyright (C) 2002, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include "CbcConfig.h"

#include <cassert>
#include <cmath>
#include <cfloat>

#ifdef COIN_HAS_CLP
#include "OsiClpSolverInterface.hpp"
#endif
#include "CbcModel.hpp"
#include "CbcMessage.hpp"
#include "CbcHeuristic.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcStrategy.hpp"
#include "CglPreProcess.hpp"
#include "OsiAuxInfo.hpp"
#include "OsiPresolve.hpp"

// Default Constructor
CbcHeuristic::CbcHeuristic() 
  :model_(NULL),
   when_(2),
   numberNodes_(200),
   feasibilityPumpOptions_(-1),
   fractionSmall_(1.0),
   heuristicName_("Unknown")
{
  // As CbcHeuristic virtual need to modify .cpp if above change
}

// Constructor from model
CbcHeuristic::CbcHeuristic(CbcModel & model)
:
  model_(&model),
  when_(2),
  numberNodes_(200),
  feasibilityPumpOptions_(-1),
  fractionSmall_(1.0),
  heuristicName_("Unknown")
{
  // As CbcHeuristic virtual need to modify .cpp if above change
}
// Copy constructor 
CbcHeuristic::CbcHeuristic(const CbcHeuristic & rhs)
:
  model_(rhs.model_),
  when_(rhs.when_),
  numberNodes_(rhs.numberNodes_),
  feasibilityPumpOptions_(rhs.feasibilityPumpOptions_),
  fractionSmall_(rhs.fractionSmall_),
  randomNumberGenerator_(rhs.randomNumberGenerator_),
  heuristicName_(rhs.heuristicName_)
{
}
// Assignment operator 
CbcHeuristic & 
CbcHeuristic::operator=( const CbcHeuristic& rhs)
{
  if (this!=&rhs) {
    model_ = rhs.model_;
    when_ = rhs.when_;
    numberNodes_ = rhs.numberNodes_;
    feasibilityPumpOptions_ = rhs.feasibilityPumpOptions_;
    fractionSmall_ = rhs.fractionSmall_;
    randomNumberGenerator_ = rhs.randomNumberGenerator_;
    heuristicName_ = rhs.heuristicName_ ;
  }
  return *this;
}

// Resets stuff if model changes
void 
CbcHeuristic::resetModel(CbcModel * model)
{
  model_=model;
}
// Set seed
void
CbcHeuristic::setSeed(int value)
{
  randomNumberGenerator_.setSeed(value);
}

// Create C++ lines to get to current state
void 
CbcHeuristic::generateCpp( FILE * fp, const char * heuristic) 
{
  // hard coded as CbcHeuristic virtual
  if (when_!=2)
    fprintf(fp,"3  %s.setWhen(%d);\n",heuristic,when_);
  else
    fprintf(fp,"4  %s.setWhen(%d);\n",heuristic,when_);
  if (numberNodes_!=200)
    fprintf(fp,"3  %s.setNumberNodes(%d);\n",heuristic,numberNodes_);
  else
    fprintf(fp,"4  %s.setNumberNodes(%d);\n",heuristic,numberNodes_);
  if (fractionSmall_!=1.0)
    fprintf(fp,"3  %s.setFractionSmall(%g);\n",heuristic,fractionSmall_);
  else
    fprintf(fp,"4  %s.setFractionSmall(%g);\n",heuristic,fractionSmall_);
  if (heuristicName_ != "Unknown")
    fprintf(fp,"3  %s.setHeuristicName(\"%s\");\n",
	    heuristic,heuristicName_.c_str()) ;
  else
    fprintf(fp,"4  %s.setHeuristicName(\"%s\");\n",
	    heuristic,heuristicName_.c_str()) ;
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
#ifdef COIN_DEVELOP
extern bool getHistoryStatistics_;
#endif
// Do mini branch and bound (return 1 if solution)
int 
CbcHeuristic::smallBranchAndBound(OsiSolverInterface * solver,int numberNodes,
                                  double * newSolution, double & newSolutionValue,
                                  double cutoff, std::string name) const
{
#ifdef COIN_HAS_CLP
  OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (solver);
  if (osiclp&&(osiclp->specialOptions()&65536)==0) {
    // go faster stripes
    if (osiclp->getNumRows()<300&&osiclp->getNumCols()<500) {
      osiclp->setupForRepeatedUse(2,0);
    } else {
      osiclp->setupForRepeatedUse(0,0);
    }
    // Turn this off if you get problems
    // Used to be automatically set
    osiclp->setSpecialOptions(osiclp->specialOptions()|(128+64));
    ClpSimplex * lpSolver = osiclp->getModelPtr();
    lpSolver->setSpecialOptions(lpSolver->specialOptions()|0x01000000); // say is Cbc (and in branch and bound)
  }
#endif
#ifdef COIN_DEVELOP
  getHistoryStatistics_=false;
#endif
  int status=0;
  int logLevel = model_->logLevel();
#define LEN_PRINT 250
  char generalPrint[LEN_PRINT];
  // Do presolve to see if possible
  int numberColumns = solver->getNumCols();
  char * reset = NULL;
  int returnCode=1;
  {
    int saveLogLevel = solver->messageHandler()->logLevel();
    if (saveLogLevel==1)
      solver->messageHandler()->setLogLevel(0);
    OsiPresolve * pinfo = new OsiPresolve();
    int presolveActions=0;
    // Allow dual stuff on integers
    presolveActions=1;
    // Do not allow all +1 to be tampered with
    //if (allPlusOnes)
    //presolveActions |= 2;
    // allow transfer of costs
    // presolveActions |= 4;
    pinfo->setPresolveActions(presolveActions);
    OsiSolverInterface * presolvedModel = pinfo->presolvedModel(*solver,1.0e-8,true,2);
    delete pinfo;
    // see if too big
    double before = 2*solver->getNumRows()+solver->getNumCols();
    if (presolvedModel) {
      int afterRows = presolvedModel->getNumRows();
      int afterCols = presolvedModel->getNumCols();
      delete presolvedModel;
      double after = 2*afterRows+afterCols;
      if (after>fractionSmall_*before&&after>300) {
	// Need code to try again to compress further using used
	const int * used =  model_->usedInSolution();
	int maxUsed=0;
	int iColumn;
	const double * lower = solver->getColLower();
	const double * upper = solver->getColUpper();
	for (iColumn=0;iColumn<numberColumns;iColumn++) {
	  if (upper[iColumn]>lower[iColumn]) {
	    if (solver->isBinary(iColumn))
	      maxUsed = CoinMax(maxUsed,used[iColumn]);
	  }
	}
	if (maxUsed) {
	  reset = new char [numberColumns];
	  int nFix=0;
	  for (iColumn=0;iColumn<numberColumns;iColumn++) {
	    reset[iColumn]=0;
	    if (upper[iColumn]>lower[iColumn]) {
	      if (solver->isBinary(iColumn)&&used[iColumn]==maxUsed) {
		bool setValue=true;
		if (maxUsed==1) {
		  double randomNumber = randomNumberGenerator_.randomDouble();
		  if (randomNumber>0.3)
		    setValue=false;
		}
		if (setValue) {
		  reset[iColumn]=1;
		  solver->setColLower(iColumn,1.0);
		  nFix++;
		}
	      }
	    }
	  }
	  pinfo = new OsiPresolve();
	  presolveActions=0;
	  // Allow dual stuff on integers
	  presolveActions=1;
	  // Do not allow all +1 to be tampered with
	  //if (allPlusOnes)
	  //presolveActions |= 2;
	  // allow transfer of costs
	  // presolveActions |= 4;
	  pinfo->setPresolveActions(presolveActions);
	  presolvedModel = pinfo->presolvedModel(*solver,1.0e-8,true,2);
	  delete pinfo;
	  if(presolvedModel) {
	    // see if too big
	    int afterRows2 = presolvedModel->getNumRows();
	    int afterCols2 = presolvedModel->getNumCols();
	    delete presolvedModel;
	    double after = 2*afterRows2+afterCols2;
	    if (after>fractionSmall_*before&&after>300) {
	      sprintf(generalPrint,"Full problem %d rows %d columns, reduced to %d rows %d columns - %d fixed gives %d, %d - still too large",
		      solver->getNumRows(),solver->getNumCols(),
		      afterRows,afterCols,nFix,afterRows2,afterCols2);
	    } else {
	      sprintf(generalPrint,"Full problem %d rows %d columns, reduced to %d rows %d columns - %d fixed gives %d, %d - ok now",
		      solver->getNumRows(),solver->getNumCols(),
		      afterRows,afterCols,nFix,afterRows2,afterCols2);
	    }
	    model_->messageHandler()->message(CBC_FPUMP1,model_->messages())
	      << generalPrint
	      <<CoinMessageEol;
	  } else {
	    returnCode=-1; // infeasible
	  }
	}
      }
    } else {
      returnCode=-1; // infeasible
    }
    solver->messageHandler()->setLogLevel(saveLogLevel);
  }
  if (returnCode==-1) {
    delete [] reset;
#ifdef COIN_DEVELOP
    getHistoryStatistics_=true;
#endif
    return returnCode;
  }
  // Reduce printout
  solver->setHintParam(OsiDoReducePrint,true,OsiHintTry);
  solver->setHintParam(OsiDoPresolveInInitial,false,OsiHintTry);
  solver->setDblParam(OsiDualObjectiveLimit,cutoff*solver->getObjSense());
  solver->initialSolve();
  if (solver->isProvenOptimal()) {
    CglPreProcess process;
    /* Do not try and produce equality cliques and
       do up to 2 passes */
    if (logLevel<=1)
      process.messageHandler()->setLogLevel(0);
    OsiSolverInterface * solver2= process.preProcessNonDefault(*solver,false,2);
    if (!solver2) {
      if (logLevel>1)
        printf("Pre-processing says infeasible\n");
      returnCode=2; // so will be infeasible
    } else {
      // see if too big
      double before = 2*solver->getNumRows()+solver->getNumCols();
      double after = 2*solver2->getNumRows()+solver2->getNumCols();
      if (after>fractionSmall_*before&&after>300) {
	sprintf(generalPrint,"Full problem %d rows %d columns, reduced to %d rows %d columns - too large",
		solver->getNumRows(),solver->getNumCols(),
		solver2->getNumRows(),solver2->getNumCols());
	model_->messageHandler()->message(CBC_FPUMP1,model_->messages())
	  << generalPrint
	  <<CoinMessageEol;
	returnCode = -1;
      } else {
	sprintf(generalPrint,"Full problem %d rows %d columns, reduced to %d rows %d columns",
		solver->getNumRows(),solver->getNumCols(),
		solver2->getNumRows(),solver2->getNumCols());
	model_->messageHandler()->message(CBC_FPUMP1,model_->messages())
	  << generalPrint
	  <<CoinMessageEol;
      }
      if (returnCode==1) {
	solver2->resolve();
	CbcModel model(*solver2);
	if (logLevel<=1)
	  model.setLogLevel(0);
	else
	  model.setLogLevel(logLevel);
	if (feasibilityPumpOptions_>=0) {
	  CbcHeuristicFPump heuristic4;
	  int pumpTune=feasibilityPumpOptions_;
	  if (pumpTune>0) {
	    /*
	      >=10000000 for using obj
	      >=1000000 use as accumulate switch
	      >=1000 use index+1 as number of large loops
	      >=100 use 0.05 objvalue as increment
	      >=10 use +0.1 objvalue for cutoff (add)
	      1 == fix ints at bounds, 2 fix all integral ints, 3 and continuous at bounds
	      4 and static continuous, 5 as 3 but no internal integers
	      6 as 3 but all slack basis!
	    */
	    double value = solver2->getObjSense()*solver2->getObjValue();
	    int w = pumpTune/10;
	    int c = w % 10;
	    w /= 10;
	    int i = w % 10;
	    w /= 10;
	    int r = w;
	    int accumulate = r/1000;
	    r -= 1000*accumulate;
	    if (accumulate>=10) {
	      int which = accumulate/10;
	      accumulate -= 10*which;
	      which--;
	      // weights and factors
	      double weight[]={0.1,0.1,0.5,0.5,1.0,1.0,5.0,5.0};
	      double factor[] = {0.1,0.5,0.1,0.5,0.1,0.5,0.1,0.5};
	      heuristic4.setInitialWeight(weight[which]);
	      heuristic4.setWeightFactor(factor[which]);
	    }
	    // fake cutoff
	    if (c) {
	      double cutoff;
	      solver2->getDblParam(OsiDualObjectiveLimit,cutoff);
	      cutoff = CoinMin(cutoff,value + 0.1*fabs(value)*c);
	      heuristic4.setFakeCutoff(cutoff);
	    }
	    if (i||r) {
	      // also set increment
	      //double increment = (0.01*i+0.005)*(fabs(value)+1.0e-12);
	      double increment = 0.0;
	      heuristic4.setAbsoluteIncrement(increment);
	      heuristic4.setAccumulate(accumulate);
	      heuristic4.setMaximumRetries(r+1);
	    }
	    pumpTune = pumpTune%100;
	    if (pumpTune==6)
	      pumpTune =13;
	    heuristic4.setWhen(pumpTune+10);
	  }
	  heuristic4.setHeuristicName("feasibility pump");
	  model.addHeuristic(&heuristic4);
	}
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
#ifdef CBC_THREAD
	if (model_->getNumberThreads()>0&&(model_->getThreadMode()&1)!=0) {
	  // See if at root node
	  bool atRoot = model_->getNodeCount()==0;
	  int passNumber = model_->getCurrentPassNumber();
	  if (atRoot&&passNumber==1)
	    model.setNumberThreads(model_->getNumberThreads());
	}
#endif
	model.setMaximumCutPassesAtRoot(CoinMin(20,model_->getMaximumCutPassesAtRoot()));
	model.setParentModel(*model_);
	model.setOriginalColumns(process.originalColumns());
	if (model.getNumCols()) {
	  setCutAndHeuristicOptions(model);
	  model.branchAndBound();
	} else {
	  // empty model
	  model.setMinimizationObjValue(model.solver()->getObjSense()*model.solver()->getObjValue());
	}
	if (logLevel>1)
	  model_->messageHandler()->message(CBC_END_SUB,model_->messages())
	    << name
	    <<CoinMessageEol;
	if (model.getMinimizationObjValue()<CoinMin(cutoff,1.0e30)) {
	  // solution
	  if (model.getNumCols())
	    returnCode=model.isProvenOptimal() ? 3 : 1;
	  else
	    returnCode=3;
	  // post process
#ifdef COIN_HAS_CLP
	  OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (model.solver());
	  if (clpSolver) {
	    ClpSimplex * lpSolver = clpSolver->getModelPtr();
	    lpSolver->setSpecialOptions(lpSolver->specialOptions()|0x01000000); // say is Cbc (and in branch and bound)
	  }
#endif
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
	  returnCode=model.isProvenInfeasible() ? 2 : 0; // so will be infeasible
	}
	if (model.status()==5)
	  returnCode=-2; // stop
	if (model.isProvenInfeasible())
	  status=1;
	else if (model.isProvenOptimal())
	  status=2;
      }
    }
  } else {
    returnCode=2; // infeasible finished
  }
  model_->setLogLevel(logLevel);
  if (reset) {
    for (int iColumn=0;iColumn<numberColumns;iColumn++) {
      if (reset[iColumn])
	solver->setColLower(iColumn,0.0);
    }
    delete [] reset;
  }
#ifdef COIN_DEVELOP
  getHistoryStatistics_=true;
  if (returnCode==1||returnCode==2) {
    if (status==1)
      printf("heuristic could add cut because infeasible (%s)\n",heuristicName_.c_str()); 
    else if (status==2)
      printf("heuristic could add cut because optimal (%s)\n",heuristicName_.c_str());
  } 
#endif
  return returnCode;
}

// Default Constructor
CbcRounding::CbcRounding() 
  :CbcHeuristic()
{
  // matrix and row copy will automatically be empty
  seed_=1;
  down_ = NULL;
  up_ = NULL;
  equal_ = NULL;
}

// Constructor from model
CbcRounding::CbcRounding(CbcModel & model)
  :CbcHeuristic(model)
{
  // Get a copy of original matrix (and by row for rounding);
  assert(model.solver());
  matrix_ = *model.solver()->getMatrixByCol();
  matrixByRow_ = *model.solver()->getMatrixByRow();
  validate();
  seed_=1;
}

// Destructor 
CbcRounding::~CbcRounding ()
{
  delete [] down_;
  delete [] up_;
  delete [] equal_;
}

// Clone
CbcHeuristic *
CbcRounding::clone() const
{
  return new CbcRounding(*this);
}
// Create C++ lines to get to current state
void 
CbcRounding::generateCpp( FILE * fp) 
{
  CbcRounding other;
  fprintf(fp,"0#include \"CbcHeuristic.hpp\"\n");
  fprintf(fp,"3  CbcRounding rounding(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp,"rounding");
  if (seed_!=other.seed_)
    fprintf(fp,"3  rounding.setSeed(%d);\n",seed_);
  else
    fprintf(fp,"4  rounding.setSeed(%d);\n",seed_);
  fprintf(fp,"3  cbcModel->addHeuristic(&rounding);\n");
}
//#define NEW_ROUNDING
// Copy constructor 
CbcRounding::CbcRounding(const CbcRounding & rhs)
:
  CbcHeuristic(rhs),
  matrix_(rhs.matrix_),
  matrixByRow_(rhs.matrixByRow_),
  seed_(rhs.seed_)
{
#ifdef NEW_ROUNDING
  int numberColumns = matrix_.getNumCols();
  down_ = CoinCopyOfArray(rhs.down_,numberColumns);
  up_ = CoinCopyOfArray(rhs.up_,numberColumns);
  equal_ = CoinCopyOfArray(rhs.equal_,numberColumns);
#else
  down_ = NULL;
  up_ = NULL;
  equal_ = NULL;
#endif  
}

// Assignment operator 
CbcRounding & 
CbcRounding::operator=( const CbcRounding& rhs)
{
  if (this!=&rhs) {
    CbcHeuristic::operator=(rhs);
    matrix_ = rhs.matrix_;
    matrixByRow_ = rhs.matrixByRow_;
#ifdef NEW_ROUNDING
    delete [] down_;
    delete [] up_;
    delete [] equal_;
    int numberColumns = matrix_.getNumCols();
    down_ = CoinCopyOfArray(rhs.down_,numberColumns);
    up_ = CoinCopyOfArray(rhs.up_,numberColumns);
    equal_ = CoinCopyOfArray(rhs.equal_,numberColumns);
#else
    down_ = NULL;
    up_ = NULL;
    equal_ = NULL;
#endif  
    seed_ = rhs.seed_;
  }
  return *this;
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
  validate();
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
  double direction = solver->getObjSense();
  double newSolutionValue = direction*solver->getObjValue();
  return solution(solutionValue,betterSolution,newSolutionValue);
}
// See if rounding will give solution
// Sets value of solution
// Assumes rhs for original matrix still okay
// At present only works with integers 
// Fix values if asked for
// Returns 1 if solution, 0 if not
int
CbcRounding::solution(double & solutionValue,
		      double * betterSolution,
		      double newSolutionValue)
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
  assert (numberRows<=solver->getNumRows());
  int numberIntegers = model_->numberIntegers();
  const int * integerVariable = model_->integerVariable();
  int i;
  double direction = solver->getObjSense();
  //double newSolutionValue = direction*solver->getObjValue();
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
        double randomNumber = randomNumberGenerator_.randomDouble();
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
        bool isInteger = (integerType[iColumn]!=0);
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
    double randomNumber = randomNumberGenerator_.randomDouble();
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
    // Just in case of some stupidity
    double objOffset=0.0;
    solver->getDblParam(OsiObjOffset,objOffset);
    newSolutionValue = -objOffset;
    for ( i=0 ; i<numberColumns ; i++ )
      newSolutionValue += objective[i]*newSolution[i];
    newSolutionValue *= direction;
    //printf("new solution value %g %g\n",newSolutionValue,solutionValue);
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
#ifdef NEW_ROUNDING
  if (!returnCode) {
#if 0
    // back to starting point
    memcpy(newSolution,solution,numberColumns*sizeof(double));
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
#endif
    int * candidate = new int [numberColumns];
    int nCandidate=0;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      bool isInteger = (integerType[iColumn]!=0);
      if (isInteger) {
	double currentValue = newSolution[iColumn];
	if (fabs(currentValue-floor(currentValue+0.5))>1.0e-8)
	  candidate[nCandidate++]=iColumn;
      }
    }
    if (true) {
      // Rounding as in Berthold
      while (nCandidate) {
	double infeasibility =1.0e-7;
	int iRow=-1;
	for (i=0;i<numberRows;i++) {
	  double value=0.0;
	  if(rowActivity[i]<rowLower[i]) {
	    value = rowLower[i]-rowActivity[i];
	  } else if(rowActivity[i]>rowUpper[i]) {
	    value = rowActivity[i]-rowUpper[i];
	  }
	  if (value>infeasibility) {
	    infeasibility = value;
	    iRow=i;
	  }
	}
	if (iRow>=0) {
	  // infeasible
	} else {
	  // feasible
	}
      }
    } else {
      // Shifting as in Berthold
    }
    delete [] candidate;
  }
#endif
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
#ifdef NEW_ROUNDING
  int numberColumns = matrix_.getNumCols();
  down_ = new unsigned short [numberColumns];
  up_ = new unsigned short [numberColumns];
  equal_ = new unsigned short [numberColumns];
  // Column copy
  const double * element = matrix_.getElements();
  const int * row = matrix_.getIndices();
  const CoinBigIndex * columnStart = matrix_.getVectorStarts();
  const int * columnLength = matrix_.getVectorLengths();
  const double * rowLower = model.solver()->getRowLower();
  const double * rowUpper = model.solver()->getRowUpper();
  for (int i=0;i<numberColumns;i++) {
    int down=0;
    int up=0;
    int equal=0;
    if (columnLength[i]>65535) {
      equal[0]=65535; 
      break; // unlikely to work
    }
    for (CoinBigIndex j=columnStart[i];
	 j<columnStart[i]+columnLength[i];j++) {
      int iRow=row[j];
      if (rowLower[iRow]>-1.0e20&&rowUpper[iRow]<1.0e20) {
	equal++;
      } else if (element[j]>0.0) {
	if (rowUpper[iRow]<1.0e20)
	  up++;
	else
	  down--;
      } else {
	if (rowLower[iRow]>-1.0e20)
	  up++;
	else
	  down--;
      }
    }
    down_[i] = (unsigned short) down;
    up_[i] = (unsigned short) up;
    equal_[i] = (unsigned short) equal;
  }
#else
  down_ = NULL;
  up_ = NULL;
  equal_ = NULL;
#endif  
}

// Default Constructor
CbcHeuristicPartial::CbcHeuristicPartial() 
  :CbcHeuristic()
{
  fixPriority_ = 10000;
}

// Constructor from model
CbcHeuristicPartial::CbcHeuristicPartial(CbcModel & model, int fixPriority, int numberNodes)
  :CbcHeuristic(model)
{
  fixPriority_ = fixPriority;
  setNumberNodes(numberNodes);
  validate();
}

// Destructor 
CbcHeuristicPartial::~CbcHeuristicPartial ()
{
}

// Clone
CbcHeuristic *
CbcHeuristicPartial::clone() const
{
  return new CbcHeuristicPartial(*this);
}
// Create C++ lines to get to current state
void 
CbcHeuristicPartial::generateCpp( FILE * fp) 
{
  CbcHeuristicPartial other;
  fprintf(fp,"0#include \"CbcHeuristic.hpp\"\n");
  fprintf(fp,"3  CbcHeuristicPartial partial(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp,"partial");
  if (fixPriority_!=other.fixPriority_)
    fprintf(fp,"3  partial.setFixPriority(%d);\n",fixPriority_);
  else
    fprintf(fp,"4  partial.setFixPriority(%d);\n",fixPriority_);
  fprintf(fp,"3  cbcModel->addHeuristic(&partial);\n");
}
//#define NEW_PARTIAL
// Copy constructor 
CbcHeuristicPartial::CbcHeuristicPartial(const CbcHeuristicPartial & rhs)
:
  CbcHeuristic(rhs),
  fixPriority_(rhs.fixPriority_)
{
}

// Assignment operator 
CbcHeuristicPartial & 
CbcHeuristicPartial::operator=( const CbcHeuristicPartial& rhs)
{
  if (this!=&rhs) {
    CbcHeuristic::operator=(rhs);
    fixPriority_ = rhs.fixPriority_;
  }
  return *this;
}

// Resets stuff if model changes
void 
CbcHeuristicPartial::resetModel(CbcModel * model)
{
  model_=model;
  // Get a copy of original matrix (and by row for partial);
  assert(model_->solver());
  validate();
}
// See if partial will give solution
// Sets value of solution
// Assumes rhs for original matrix still okay
// At present only works with integers 
// Fix values if asked for
// Returns 1 if solution, 0 if not
int
CbcHeuristicPartial::solution(double & solutionValue,
		      double * betterSolution)
{
  // Return if already done
  if (fixPriority_<0)
    return 0; // switched off
  const double * hotstartSolution = model_->hotstartSolution();
  const int * hotstartPriorities = model_->hotstartPriorities();
  if (!hotstartSolution)
    return 0;
  OsiSolverInterface * solver = model_->solver();
  
  int numberIntegers = model_->numberIntegers();
  const int * integerVariable = model_->integerVariable();
  
  OsiSolverInterface * newSolver = model_->continuousSolver()->clone();
  const double * colLower = newSolver->getColLower();
  const double * colUpper = newSolver->getColUpper();

  double primalTolerance;
  solver->getDblParam(OsiPrimalTolerance,primalTolerance);
    
  int i;
  int numberFixed=0;
  int returnCode=0;

  for (i=0;i<numberIntegers;i++) {
    int iColumn=integerVariable[i];
    if (abs(hotstartPriorities[iColumn])<=fixPriority_) {
      double value = hotstartSolution[iColumn];
      double lower = colLower[iColumn];
      double upper = colUpper[iColumn];
      value = CoinMax(value,lower);
      value = CoinMin(value,upper);
      if (fabs(value-floor(value+0.5))<1.0e-8) {
	value = floor(value+0.5);
	newSolver->setColLower(iColumn,value);
	newSolver->setColUpper(iColumn,value);
	numberFixed++;
      }
    }
  }
  if (numberFixed>numberIntegers/5-100000000) {
#ifdef COIN_DEVELOP
    printf("%d integers fixed\n",numberFixed);
#endif
    returnCode = smallBranchAndBound(newSolver,numberNodes_,betterSolution,solutionValue,
				     model_->getCutoff(),"CbcHeuristicPartial");
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
  fixPriority_=-1; // switch off
  
  delete newSolver;
  return returnCode;
}
// update model
void CbcHeuristicPartial::setModel(CbcModel * model)
{
  model_ = model;
  assert(model_->solver());
  // make sure model okay for heuristic
  validate();
}
// Validate model i.e. sets when_ to 0 if necessary (may be NULL)
void 
CbcHeuristicPartial::validate() 
{
  if (model_&&when()<10) {
    if (model_->numberIntegers()!=
        model_->numberObjects())
      setWhen(0);
  }
}

// Default Constructor
CbcSerendipity::CbcSerendipity() 
  :CbcHeuristic()
{
}

// Constructor from model
CbcSerendipity::CbcSerendipity(CbcModel & model)
  :CbcHeuristic(model)
{
}

// Destructor 
CbcSerendipity::~CbcSerendipity ()
{
}

// Clone
CbcHeuristic *
CbcSerendipity::clone() const
{
  return new CbcSerendipity(*this);
}
// Create C++ lines to get to current state
void 
CbcSerendipity::generateCpp( FILE * fp) 
{
  fprintf(fp,"0#include \"CbcHeuristic.hpp\"\n");
  fprintf(fp,"3  CbcSerendipity serendipity(*cbcModel);\n");
  CbcHeuristic::generateCpp(fp,"serendipity");
  fprintf(fp,"3  cbcModel->addHeuristic(&serendipity);\n");
}

// Copy constructor 
CbcSerendipity::CbcSerendipity(const CbcSerendipity & rhs)
:
  CbcHeuristic(rhs)
{
}

// Assignment operator 
CbcSerendipity & 
CbcSerendipity::operator=( const CbcSerendipity& rhs)
{
  if (this!=&rhs) {
    CbcHeuristic::operator=(rhs);
  }
  return *this;
}

// Returns 1 if solution, 0 if not
int
CbcSerendipity::solution(double & solutionValue,
			 double * betterSolution)
{
  if (!model_)
    return 0;
  // get information on solver type
  OsiAuxInfo * auxInfo = model_->solver()->getAuxiliaryInfo();
  OsiBabSolver * auxiliaryInfo = dynamic_cast< OsiBabSolver *> (auxInfo);
  if (auxiliaryInfo)
    return auxiliaryInfo->solution(solutionValue,betterSolution,model_->solver()->getNumCols());
  else
    return 0;
}
// update model
void CbcSerendipity::setModel(CbcModel * model)
{
  model_ = model;
}
// Resets stuff if model changes
void 
CbcSerendipity::resetModel(CbcModel * model)
{
  model_ = model;
}
  
