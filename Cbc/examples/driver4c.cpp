// Copyright (C) 2009, International Business Machines
// Corporation and others.  All Rights Reserved.
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif

#include <cassert>
#include <iomanip> 


#include "CbcModel.hpp"
#include "OsiClpSolverInterface.hpp"
//#include "CbcSolver.hpp"
//#include "CbcStrategy.hpp"
#include "CbcCutGenerator.hpp"
#include "CglProbing.hpp"
#include "CglGomory.hpp"
#include "CglKnapsackCover.hpp"
#include "CglClique.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglFlowCover.hpp"
#include "CglTwomir.hpp"
#include "CbcHeuristic.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcHeuristicDiveCoefficient.hpp"
#include "CbcHeuristicRINS.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcCompareActual.hpp"

#include  "CoinTime.hpp"


/************************************************************************

This main program shows how to try and obtain all solutions (within a certain range) to a problem.

It will finish saying no solutions found - but event handler will have
saved them.

************************************************************************/
/* Meaning of whereFrom:
   1 after initial solve by dualsimplex etc
   2 after preprocessing
   3 just before branchAndBound (so user can override)
   4 just after branchAndBound (before postprocessing)
   5 after postprocessing
*/
/* Meaning of model status is as normal
   status
      -1 before branchAndBound
      0 finished - check isProvenOptimal or isProvenInfeasible to see if solution found
      (or check value of best solution)
      1 stopped - on maxnodes, maxsols, maxtime
      2 difficulties so run was abandoned
      (5 event user programmed event occurred) 

      cbc secondary status of problem
        -1 unset (status_ will also be -1)
	0 search completed with solution
	1 linear relaxation not feasible (or worse than cutoff)
	2 stopped on gap
	3 stopped on nodes
	4 stopped on time
	5 stopped on user event
	6 stopped on solutions
	7 linear relaxation unbounded

   but initially check if status is 0 and secondary status is 1 -> infeasible
   or you can check solver status.
*/
/* Return non-zero to return quickly */   
static int callBack(CbcModel * model, int whereFrom)
{
  int returnCode=0;
  switch (whereFrom) {
  case 1:
  case 2:
    if (!model->status()&&model->secondaryStatus())
      returnCode=1;
    break;
  case 3:
    {
    }
    break;
  case 4:
    // If not good enough could skip postprocessing
    break;
  case 5:
    break;
  default:
    abort();
  }
  return returnCode;
}
#include "CbcEventHandler.hpp"
/** This is so user can trap events and do useful stuff.  

    CbcModel model_ is available as well as anything else you care 
    to pass in
*/

class MyEventHandler3 : public CbcEventHandler {
  
public:
  /**@name Overrides */
  //@{
  virtual CbcAction event(CbcEvent whichEvent);
  //@}

  /**@name sets and gets*/
  //@{
  /// Set maximum solutions
  inline void setMaximumSolutions(int value)
  { maximumSolutions_=value;}
  /// Set cutoff
  inline void setCutoff(double value)
  { cutoff_=value;}
  //@}

  /**@name Constructors, destructor etc*/
  //@{
  /** Default constructor. */
  MyEventHandler3();
  /// Constructor with pointer to model (redundant as setEventHandler does)
  MyEventHandler3(CbcModel * model);
  /** Destructor */
  virtual ~MyEventHandler3();
  /** The copy constructor. */
  MyEventHandler3(const MyEventHandler3 & rhs);
  /// Assignment
  MyEventHandler3& operator=(const MyEventHandler3 & rhs);
  /// Clone
  virtual CbcEventHandler * clone() const ;
  //@}
   
    
protected:
  // Save a solution if wanted
  void saveSolution(const char * whereFrom,
		    double objValue,
		    const double * solution);
  // returns true if duplicate
  bool sameSolution(double objValue,
		    const double * solution);
  // creates cut
  void createCut(const double * solution, int type);
  // data goes here
  // Best objective found
  double bestSoFar_;
  // Cutoff which user can set
  double cutoff_;
  // File for saving all solutions
  FILE * fp_;
  // Original model before preprocessing 
  CbcModel * originalModel_;
  // Saved solutions (last value is objective value)
  double ** solutions_;
  // Maximum number of saved solutions
  int maximumSolutions_;
  // Number of solutions found
  int numberSolutions_;

};
//-------------------------------------------------------------------
// Default Constructor 
//-------------------------------------------------------------------
MyEventHandler3::MyEventHandler3 () 
  : CbcEventHandler(),
    bestSoFar_(COIN_DBL_MAX),
    cutoff_(1.0e50),
    fp_(NULL),
    originalModel_(NULL),
    solutions_(NULL),
    maximumSolutions_(10),
    numberSolutions_(0)
{
}

//-------------------------------------------------------------------
// Copy constructor 
//-------------------------------------------------------------------
MyEventHandler3::MyEventHandler3 (const MyEventHandler3 & rhs) 
  : CbcEventHandler(rhs),
    bestSoFar_(rhs.bestSoFar_),
    cutoff_(rhs.cutoff_),
    fp_(NULL),
    originalModel_(rhs.originalModel_),
    solutions_(NULL),
    maximumSolutions_(rhs.maximumSolutions_),
    numberSolutions_(rhs.numberSolutions_)
{  
  if (rhs.solutions_) {
    int numberColumnsP=model_->getNumCols()+2;
    solutions_ = new double * [maximumSolutions_];
    for (int i=0;i<maximumSolutions_;i++) {
      solutions_[i]=CoinCopyOfArray(rhs.solutions_[i],numberColumnsP);
    }
  }
}

// Constructor with pointer to model
MyEventHandler3::MyEventHandler3(CbcModel * model)
  : CbcEventHandler(model),
    bestSoFar_(COIN_DBL_MAX),
    cutoff_(1.0e50),
    fp_(NULL),
    originalModel_(model),
    solutions_(NULL),
    maximumSolutions_(10),
    numberSolutions_(0)
{
}

//-------------------------------------------------------------------
// Destructor 
//-------------------------------------------------------------------
MyEventHandler3::~MyEventHandler3 ()
{
  if (solutions_) {
    for (int i=0;i<maximumSolutions_;i++) {
      delete [] solutions_[i];
    }
    delete [] solutions_;
  }
}

//----------------------------------------------------------------
// Assignment operator 
//-------------------------------------------------------------------
MyEventHandler3 &
MyEventHandler3::operator=(const MyEventHandler3& rhs)
{
  if (this != &rhs) {
    CbcEventHandler::operator=(rhs);
    originalModel_ = rhs.originalModel_;
    bestSoFar_ = rhs.bestSoFar_;
    cutoff_ = rhs.cutoff_;
    if (solutions_) {
      for (int i=0;i<maximumSolutions_;i++) {
	delete [] solutions_[i];
      }
      delete [] solutions_;
      solutions_ = NULL;
    }
    maximumSolutions_ = rhs.maximumSolutions_;
    if (rhs.solutions_) {
      int numberColumnsP=model_->getNumCols()+2;
      solutions_ = new double * [maximumSolutions_];
      for (int i=0;i<maximumSolutions_;i++) {
	solutions_[i]=CoinCopyOfArray(rhs.solutions_[i],numberColumnsP);
      }
    }
    numberSolutions_ = rhs.numberSolutions_;
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CbcEventHandler * MyEventHandler3::clone() const
{
  return new MyEventHandler3(*this);
}
// Save a solution if wanted
void 
MyEventHandler3::saveSolution(const char * whereFrom,
			      double objValue,
			      const double * solution)
{
  int numberColumns = model_->getNumCols();
  if (objValue<bestSoFar_) 
    bestSoFar_ = objValue;
  numberSolutions_++;
  if (!fp_) {
    fp_=fopen("sols.sol","w");
    assert (fp_);
  }
  fprintf(fp_,"%s %d %g ",whereFrom,numberSolutions_,objValue); 
  for (int i=0;i<numberColumns;i++) {
    if (fabs(solution[i])>1.0e-8) {
      fprintf(fp_,"(%d %g) ",i,solution[i]);
    }
  }
  fprintf(fp_,"\n");
  // Save if not bad
  if (!solutions_) {
    solutions_ = new double * [maximumSolutions_];
    for (int i=0;i<maximumSolutions_;i++)
      solutions_[i]=NULL;
  }
  int which=-1;
  double worst=-COIN_DBL_MAX;
  for (int i=0;i<maximumSolutions_;i++) {
    if (!solutions_[i]) {
      which=i;
      solutions_[i]= new double [numberColumns+2];
      worst=COIN_DBL_MAX;
      break;
    } else {
      double value = solutions_[i][numberColumns];
      if (value>worst) {
	worst=value;
	which=i;
      }
    }
  }
  OsiSolverInterface * solver=model_->solver();
  if (objValue<worst) {
    // save
    double * sol = solutions_[which];
    for (int i=0;i<numberColumns;i++) {
      double value=solution[i];
      if (solver->isInteger(i))
	value = floor(value+0.5);
      sol[i]=value;
    }
    sol[numberColumns]=objValue;
    sol[numberColumns+1]=numberSolutions_;
  }
}
// returns true if duplicate
bool 
MyEventHandler3::sameSolution(double objValue,
			      const double * solution)
{
  OsiSolverInterface * solver=model_->solver();
  int numberColumns = model_->getNumCols();
  double same=false;
  if (solutions_) {
    // Check for duplicates
    for (int j=0;j<maximumSolutions_;j++) {
      if (solutions_[j]) {
	double * sol = solutions_[j];
	if (fabs(sol[numberColumns]-objValue)<1.0e-5) {
	  bool thisSame=true;
	  for (int i=0;i<numberColumns;i++) {
	    if (fabs(sol[i]-solution[i])>1.0e-6&&
		solver->isInteger(i)) {
	      thisSame=false;
	      break;
	    }
	  }
	  if (thisSame) {
	    int k = static_cast<int>(sol[numberColumns+1]);
	    if (k!=numberSolutions_)
	      printf("*** Solutions %d and %d same\n",
		     numberSolutions_+1,k);
	    same=true;
	  }
	}
      }
    }
  }
  if (!same) {
    printf("value of solution is %g\n",objValue);
    for (int i=0;i<numberColumns;i++) {
      if (fabs(solution[i])>1.0e-8&&solver->isInteger(i)) {
	printf("%d %g\n",i,solution[i]);
      }
    }
  }
  return same;
}
CbcEventHandler::CbcAction 
MyEventHandler3::event(CbcEvent whichEvent)
{
  // If in sub tree carry on
  if (!model_->parentModel()&&maximumSolutions_>=0) {
    if (whichEvent==solution||whichEvent==heuristicSolution) {
      const double * bestSolution = model_->bestSolution();
      if (!bestSolution)
	return noAction;
      double objValue = model_->getObjValue();
      if (!sameSolution(objValue,bestSolution)) {
	saveSolution((whichEvent==solution) ? "solution" : "heuristicSolution",
		       objValue,bestSolution);
      }
      double newCutoff = CoinMin(bestSoFar_+0.1*fabs(bestSoFar_),cutoff_);
      // Set cutoff here as may have been changed
      if (fabs(model_->getCutoff()-newCutoff)>1.0e-3) {
	printf("Changing cutoff from %g to %g\n",
	       model_->getCutoff(),newCutoff);
	model_->setCutoff(newCutoff);
	// Also need to reset best objective value
	model_->setBestObjectiveValue(newCutoff+1.0e-1
				      +model_->getCutoffIncrement());
      }
    } else if (whichEvent==beforeSolution1||whichEvent==beforeSolution2) {
      double objValue = model_->getObjValue();
      if (objValue>model_->getCutoff())
	return noAction; // Must have slipped through
      /* If preprocessing was done solution will be to processed model
	 So would need to postprocess at end or maybe simpler to use
	 model_->originalColumns() */
      const double * bestSolution = model_->bestSolution();
      assert (bestSolution);
      // save as won't be done by 'solution'
      bool same = sameSolution(objValue,bestSolution);
      if (!same)
	saveSolution("killSolution",
		   objValue,bestSolution);
      // Set cutoff here as may have been changed
      double newCutoff = CoinMin(bestSoFar_+0.1*fabs(bestSoFar_),cutoff_);
      if (fabs(model_->getCutoff()-newCutoff)>1.0e-3) {
	printf("Changing cutoff from %g to %g\n",
	       model_->getCutoff(),newCutoff);
	model_->setCutoff(newCutoff);
	model_->setBestObjectiveValue(newCutoff+1.0e-1
				      +model_->getCutoffIncrement());
      }
      if (whichEvent==beforeSolution2) {
	return killSolution;
      } else {
	// add global cut
	createCut(bestSolution,0);
 	return addCuts;
      }
      return killSolution;
    } else if (whichEvent==endSearch) {
      printf("%d solutions found\n",numberSolutions_);
      double * obj = new double [maximumSolutions_];
      int * whichSolution = new int [maximumSolutions_];
      int nBest=0;
      int numberColumns=model_->getNumCols();
      for (int j=0;j<maximumSolutions_;j++) {
	if (solutions_[j]) {
	  double * sol = solutions_[j];
	  whichSolution[nBest]=j;
	  obj[nBest++]=sol[numberColumns];
	}
      }
      CoinSort_2(obj,obj+nBest,whichSolution);
      nBest = CoinMin(nBest,10);
      printf("Best %d solutions are :-\n",nBest);
      for (int j=0;j<nBest;j++)
	printf("Objective %g\n",obj[j]);
      fclose(fp_);
      /* Use one solution to solve original model (if preprocessed).
	 This would of course be tuned and may not be necessary
      */
      if (solutions_&&originalModel_&&model_->originalColumns()) {
	CbcModel model(*originalModel_);
	CbcEventHandler * handler=model.getEventHandler();
	MyEventHandler3	* myHandler = dynamic_cast<MyEventHandler3 *> (handler);
	assert (myHandler);
	myHandler->maximumSolutions_ = -1; // switch off handler
	OsiSolverInterface * solverOriginal = model.solver();
	OsiSolverInterface * solver = model_->solver();
	const int * originalColumns = model_->originalColumns();
	int numberColumns=solver->getNumCols();
	// use best
	int which=-1;
	double best=COIN_DBL_MAX;
	for (int i=0;i<maximumSolutions_;i++) {
	  if (solutions_[i]) {
	    double value = solutions_[i][numberColumns];
	    if (value<best) {
	      best=value;
	      which=i;
	    }
	  }
	}
	assert (which>=0);
	double * sol = solutions_[which];
	// Extra columns may have been added
	int nOriginal = solverOriginal->getNumCols();
	for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	  // Test for integer probably not needed
	  if (true||solver->isInteger(iColumn)) {
	    double value=sol[iColumn];
	    int kColumn=originalColumns[iColumn];
	    if (kColumn<nOriginal) {
	      solverOriginal->setColLower(kColumn,value);
	      solverOriginal->setColUpper(kColumn,value);
	    }
	  }
	}
	printf("======= start check Branch and Bound ======\n");
	model.branchAndBound();
	printf("======= end   check Branch and Bound ======\n");
      } else if (solutions_) {
	OsiSolverInterface * solver = model_->solver();
	for (int j=0;j<nBest;j++) {
	  double * sol = solutions_[whichSolution[j]];
	  printf("Solution %d has objective of %g\n",
		 j+1,sol[numberColumns]);
	  // For this example just print first 10 nonzero
	  int n=10;
	  for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	    if (solver->isInteger(iColumn)) {
	      double value=sol[iColumn];
	      if (fabs(value)>1.0e-7) {
		n--;
		if (n<0)
		  break; // stop printing
		printf("%d has value %g\n",iColumn,value);
	      }
	    }
	  }
	  if (n<0)
	    printf("........ etc .......\n");
	}
      }
      delete [] obj;
      delete [] whichSolution;
    }
  }
  return noAction; // carry on
}

// create cut
void
MyEventHandler3::createCut(const double * solution,int /*type*/)
{
  // Get information on continuous solver
  OsiSolverInterface * solver;
  //if (!type)
    solver = model_->continuousSolver();
    //else
    //solver = model_->solver();
  if (solver) {
    /* It is up to the user what sort of cuts to apply
       (also see FORCE_CUTS lower down)
       The classic one here, which just cuts off solution,
       does not do a lot - often the cuts are never applied!
    */
    int numberColumns=solver->getNumCols();
    const double * lower = solver->getColLower();
    const double * upper = solver->getColUpper();
    bool good=true;
    for (int iColumn=0;iColumn<numberColumns;iColumn++) {
      if (solver->isInteger(iColumn)) {
	double value=solution[iColumn];
	if (value>lower[iColumn]+0.9&&
	    value<upper[iColumn]-0.9) {
	  good=false;
	  printf("Can't add cut as general integers\n");
	  break;
	}
      }
    }
    if (good) {
      double * cut = new double [numberColumns];
      int * which = new int [numberColumns];
      /* It is up to the user what sort of cuts to apply
	 The classic row cut here, which just cuts off solution,
	 does not do a lot - often the cuts are never applied!
      */
      // Should be -1.0 for classic cut
      double rhs=-1.0;
      int n=0;
      double sum=0.0;
      for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	if (solver->isInteger(iColumn)) {
	  if (upper[iColumn]>lower[iColumn]) {
	    double value=solution[iColumn];
	    sum += value;
	    if (upper[iColumn]-value<0.1) {
	      rhs += upper[iColumn];
	      cut[n]=1.0;
	      which[n++]=iColumn;
	    } else {
	      rhs -= lower[iColumn];
	      cut[n]=-1.0;
	      which[n++]=iColumn;
	    }
	  }
	}
      }
      assert (sum>rhs+0.9);
      printf("CUT has %d entries\n",n);
      OsiRowCut newCut;
      newCut.setLb(-COIN_DBL_MAX);
      newCut.setUb(rhs);
      newCut.setRow(n,which,cut,false);
      delete [] which;
      delete [] cut;
      model_->makeGlobalCut(newCut);
    }
  }
}
// Optional function to do first solve to get best
static double firstSolve(OsiClpSolverInterface & solver1)
{
  double time1 = CoinCpuTime();
  //OsiClpSolverInterface * solver1 = solver->clone();
  CbcModel model(solver1);
  // Now do requested modifications
  CbcModel * cbcModel = & model;
  OsiSolverInterface * osiModel = model.solver();
  OsiClpSolverInterface * osiclpModel = dynamic_cast< OsiClpSolverInterface*> (osiModel);
  ClpSimplex * clpModel = osiclpModel->getModelPtr();

  // Set changed values

  CglProbing probing;
  probing.setMaxPass(1);
  probing.setMaxProbe(123);
  probing.setMaxLook(10);
  probing.setMaxElements(200);
  probing.setMaxPassRoot(1);
  probing.setMaxProbeRoot(123);
  probing.setMaxLookRoot(20);
  probing.setMaxElementsRoot(300);
  probing.setRowCuts(3);
  probing.setUsingObjective(1);
  cbcModel->addCutGenerator(&probing,-1,"Probing",true,false,false,-100,-1,-1);
  cbcModel->cutGenerator(0)->setTiming(true);
  
  CglGomory gomory;
  gomory.setLimitAtRoot(2000);
  gomory.setAwayAtRoot(0.005);
  cbcModel->addCutGenerator(&gomory,-98,"Gomory",true,false,false,-100,-1,-1);
  cbcModel->cutGenerator(1)->setTiming(true);
  
  CglKnapsackCover knapsackCover;
  cbcModel->addCutGenerator(&knapsackCover,-98,"KnapsackCover",true,false,false,-100,-1,-1);
  cbcModel->cutGenerator(2)->setTiming(true);
  
  CglClique clique;
  clique.setStarCliqueReport(false);
  clique.setRowCliqueReport(false);
  clique.setMinViolation(0.1);
  cbcModel->addCutGenerator(&clique,-98,"Clique",true,false,false,-100,-1,-1);
  cbcModel->cutGenerator(3)->setTiming(true);
  
  CglMixedIntegerRounding2 mixedIntegerRounding2;
  cbcModel->addCutGenerator(&mixedIntegerRounding2,-98,"MixedIntegerRounding2",true,false,false,-100,-1,-1);
  cbcModel->cutGenerator(4)->setTiming(true);
  
  CglFlowCover flowCover;
  cbcModel->addCutGenerator(&flowCover,-98,"FlowCover",true,false,false,-100,-1,-1);
  cbcModel->cutGenerator(5)->setTiming(true);
  
  CglTwomir twomir;
  twomir.setMaxElements(250);
  cbcModel->addCutGenerator(&twomir,-99,"Twomir",true,false,false,-100,-1,-1);
  cbcModel->cutGenerator(6)->setTiming(true);
  
  CbcRounding rounding(*cbcModel);
  rounding.setHeuristicName("rounding");
  cbcModel->addHeuristic(&rounding);
  
  CbcHeuristicGreedyCover heuristicGreedyCover(*cbcModel);
  heuristicGreedyCover.setHeuristicName("greedy cover");
  heuristicGreedyCover.setWhereFrom(1);
  cbcModel->addHeuristic(&heuristicGreedyCover);
  
  CbcHeuristicGreedyEquality heuristicGreedyEquality(*cbcModel);
  heuristicGreedyEquality.setHeuristicName("greedy equality");
  heuristicGreedyEquality.setWhereFrom(1);
  cbcModel->addHeuristic(&heuristicGreedyEquality);
  
  CbcHeuristicDiveCoefficient heuristicDiveCoefficient(*cbcModel);
  heuristicDiveCoefficient.setWhen(3);
  heuristicDiveCoefficient.setHeuristicName("DiveCoefficient");
  heuristicDiveCoefficient.setWhereFrom(493);
  cbcModel->addHeuristic(&heuristicDiveCoefficient);
  
  CbcHeuristicRINS heuristicRINS(*cbcModel);
  heuristicRINS.setFractionSmall(0.5);
  heuristicRINS.setHeuristicName("RINS");
  heuristicRINS.setDecayFactor(5);
  heuristicRINS.setWhereFrom(65289);
  cbcModel->addHeuristic(&heuristicRINS);
  
  CbcHeuristicLocal heuristicLocal(*cbcModel);
  heuristicLocal.setFractionSmall(0.5);
  heuristicLocal.setHeuristicName("combine solutions");
  heuristicLocal.setSearchType(1);
  cbcModel->addHeuristic(&heuristicLocal);
  
  CbcCompareDefault compare;
  cbcModel->setNodeComparison(compare);
  cbcModel->setPrintFrequency(100);
  cbcModel->setSpecialOptions(66050);
  cbcModel->messageHandler()->setLogLevel(1);
  cbcModel->setMaximumCutPasses(4);
  cbcModel->setMinimumDrop(3.83124e-05);
  clpModel->scaling(4);
  clpModel->setLogLevel(0);
  // For branchAndBound this may help
  clpModel->defaultFactorizationFrequency();
  clpModel->setDualBound(1.0001e+08);
  clpModel->setPerturbation(50);
  osiclpModel->setSpecialOptions(1057);
  osiclpModel->messageHandler()->setLogLevel(0);
  osiclpModel->setIntParam(OsiMaxNumIterationHotStart,100);
  osiclpModel->setHintParam(OsiDoReducePrint,true,OsiHintTry);

  // Solve

  cbcModel->initialSolve();
  if (clpModel->tightenPrimalBounds()!=0) {
    std::cout<<"Problem is infeasible - tightenPrimalBounds!"<<std::endl;
    exit(1);
  }
  clpModel->dual();  // clean up
  cbcModel->branchAndBound();
  // For best solution
  int numberColumns = solver1.getNumCols();
  if (cbcModel->getMinimizationObjValue()<1.0e50) {
    // put back in original solver
    solver1.setColSolution(cbcModel->bestSolution());
    const double * solution = solver1.getColSolution();

    // Get names from solver1 (as OsiSolverInterface may lose)
    std::vector<std::string> columnNames = *solver1.getModelPtr()->columnNames();
    
    int iColumn;
    std::cout<<std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14);
    
    std::cout<<"--------------------------------------"<<std::endl;
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      double value=solution[iColumn];
      if (fabs(value)>1.0e-7&&solver1.isInteger(iColumn)) 
	std::cout<<std::setw(6)<<iColumn<<" "
                 <<columnNames[iColumn]<<" "
                 <<value<<std::endl;
    }
    std::cout<<"--------------------------------------"<<std::endl;
  
    std::cout<<std::resetiosflags(std::ios::fixed|std::ios::showpoint|std::ios::scientific);
  }
  printf("First solve took %g seconds\n",
	 CoinCpuTime()-time1);
  return model.getMinimizationObjValue();
}
int main (int argc, const char *argv[])
{

  OsiClpSolverInterface solver1;
  // Read in model using argv[1]
  if (argc<2) {
    printf("First argument is file name\n");
    exit(1);
  }
  int numMpsReadErrors = solver1.readMps(argv[1],"");
  assert(numMpsReadErrors==0);
  // Optionally do first solve to get cutoff
  double cutoff=1.0e50;
  cutoff = firstSolve(solver1);
  // Set cutoff to whatever you want
  cutoff += 0.1*fabs(cutoff);
  printf("Looking for solutions better than %g\n",cutoff);
  // Pass to Cbc initialize defaults 
  CbcModel model(solver1);
  CbcMain0(model);
  // Event handler
  MyEventHandler3 eventHandler(&model);
  eventHandler.setCutoff(cutoff);
  eventHandler.setMaximumSolutions(300);
  model.passInEventHandler(&eventHandler);
  /* Now go into code for standalone solver
     Could copy arguments and add -quit at end to be safe
     but this will do
  */
  if (argc>2) {
    CbcMain1(argc-1,argv+1,model,callBack);
  } else {
    const char * argv2[]={"driver4c","-preprocess","off","-feas","off",
			  "-solve","-quit"};
    model.setCutoff(cutoff);
    CbcMain1(7,argv2,model,callBack);
  }
  return 0;
}    
