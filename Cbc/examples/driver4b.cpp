/* $Id: driver4b.cpp 1173 2009-06-04 09:44:10Z forrest $ */
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
#include "CbcCompareUser.hpp"
#include "CbcSolver.hpp"
#include "CbcStrategy.hpp"
#include "CbcCutGenerator.hpp"
#include "CglProbing.hpp"
#include "CglGomory.hpp"
#include "CglKnapsackCover.hpp"
#include "CbcHeuristic.hpp"

#include  "CoinTime.hpp"

//#############################################################################


/************************************************************************

This main program shows how to take advantage of the standalone cbc in your program,
while still making major modifications.
First it reads in an integer model from an mps file
Then it initializes the integer model with cbc defaults
Then it calls CbcMain1 passing all parameters apart from first but with callBack to modify stuff
It also shows how to use CbcStrategy to modify cut generators or heuristics
It can also be used to find many solutions (only coded for minimization and for 0-1 integers) (search for KILL_SOLUTION)
Finally it prints solution

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
      CbcCompareUser compare;
      model->setNodeComparison(compare);
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
  void createCut(const double * solution);
  // data goes here
  // Best objective found
  double bestSoFar_;
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
  int numberColumns = model_->getNumCols();
  double same=false;
  if (solutions_) {
    OsiSolverInterface * solver=model_->solver();
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
	    printf("*** Solutions %d and %d same\n",
		   numberSolutions_+1,static_cast<int>(sol[numberColumns+1]));
	    same=true;
	  }
	}
      }
    }
  }
  if (!same) {
    printf("value of solution is %g\n",objValue);
    for (int i=0;i<numberColumns;i++) {
      if (fabs(solution[i])>1.0e-8) {
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
#ifdef STOP_EARLY
      return stop; // say finished
#endif
      const double * bestSolution = model_->bestSolution();
      if (!bestSolution)
	return noAction;
      double objValue = model_->getObjValue();
      /* if KILL_SOLUTION not set then should be as driver4 was
	 if 1 then all solutions rejected and cuts added
	 if 2 then accepted and cuts added.  In practice this
	 might be decided by user.
	 
	 If FORCE_CUTS defined then the Global cuts will always
	 be added - even if not violated - see later on.
      */
#define KILL_SOLUTION 2
      //#define FORCE_CUTS
#ifndef KILL_SOLUTION
      // If preprocessing was done solution will be to processed model
      int numberColumns = model_->getNumCols();
      printf("value of solution is %g\n",objValue);
      for (int i=0;i<numberColumns;i++) {
	if (fabs(bestSolution[i])>1.0e-8) {
	  printf("%d %g\n",i,bestSolution[i]);
	}
      }
#else
      if (!sameSolution(objValue,bestSolution)) {
	saveSolution((whichEvent==solution) ? "solution" : "heuristicSolution",
		       objValue,bestSolution);
#if KILL_SOLUTION==2
	// Put in cut if we get here
	createCut(bestSolution);
#endif
      }
      double newCutoff = bestSoFar_+0.1*fabs(bestSoFar_);
      // Set cutoff here as may have been changed
      if (fabs(model_->getCutoff()-newCutoff)>1.0e-3) {
	printf("Changing cutoff from %g to %g\n",
	       model_->getCutoff(),newCutoff);
	model_->setCutoff(newCutoff);
	// Also need to reset best objective value
	model_->setBestObjectiveValue(newCutoff+1.0e-1
				      +model_->getCutoffIncrement());
      }
#endif
    } else if (whichEvent==beforeSolution) {
#ifdef KILL_SOLUTION
      double objValue = model_->getObjValue();
      if (objValue>model_->getCutoff())
	return noAction; // Must have slipped through
      /* If preprocessing was done solution will be to processed model
	 So would need to postprocess at end or maybe simpler to use
	 model_->originalColumns() */
      const double * bestSolution = model_->bestSolution();
      assert (bestSolution);
#if KILL_SOLUTION==1
      // save as won't be done by 'solution'
      bool same = sameSolution(objValue,bestSolution);
      if (!same)
	saveSolution("killSolution",
		   objValue,bestSolution);
      // Set cutoff here as may have been changed
      double newCutoff = bestSoFar_+0.1*fabs(bestSoFar_);
      if (fabs(model_->getCutoff()-newCutoff)>1.0e-3) {
	printf("Changing cutoff from %g to %g\n",
	       model_->getCutoff(),newCutoff);
	model_->setCutoff(newCutoff);
	model_->setBestObjectiveValue(newCutoff+1.0e-1
				      +model_->getCutoffIncrement());
      }
      if (same)
	return noAction;
#endif
      createCut(bestSolution);
#if KILL_SOLUTION==1
      return killSolution;
#endif
#endif
    } else if (whichEvent==endSearch) {
#ifdef KILL_SOLUTION
      printf("%d solutions found\n",numberSolutions_);
      fclose(fp_);
      /* Use one solution to solve original model (if preprocessed).
	 This would of course be tuned and may not be necessary
      */
      if (solutions_&&originalModel_) {
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
	for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	  // Test for integer probably not needed
	  if (true||solver->isInteger(iColumn)) {
	    double value=sol[iColumn];
	    int kColumn=originalColumns[iColumn];
	    solverOriginal->setColLower(kColumn,value);
	    solverOriginal->setColUpper(kColumn,value);
	  }
	}
	printf("======= start check Branch and Bound ======\n");
	model.branchAndBound();
	printf("======= end   check Branch and Bound ======\n");
      }
#endif
    }
  }
  return noAction; // carry on
}

// create cut
void
MyEventHandler3::createCut(const double * solution)
{
  /* TYPE_CUT 0 is off
     1 is row cut
     2 is column cut */
#define TYPE_CUT 1
#if TYPE_CUT
  // Get information on continuous solver
  OsiSolverInterface * solver = model_->continuousSolver();
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
#if TYPE_CUT==1
      double * cut = new double [numberColumns];
      int * which = new int [numberColumns];
      /* It is up to the user what sort of cuts to apply
	 The classic row cut here, which just cuts off solution,
	 does not do a lot - often the cuts are never applied!
      */
      // Should be -1.0 for classic cut
      double rhs=-6.0;
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
#else
      // Fix a variable other way
      int firstAtLb=-1;
      int firstAtUb=-1;
      for (int iColumn=0;iColumn<numberColumns;iColumn++) {
	if (solver->isInteger(iColumn)) {
	  if (upper[iColumn]>lower[iColumn]) {
	    double value=solution[iColumn];
	    if (upper[iColumn]-value<0.1) {
	      if (firstAtUb<0)
		firstAtUb=iColumn;
	    } else {
	      if (firstAtLb<0)
		firstAtLb=iColumn;
	    }
	  }
	}
      }
      OsiColCut newCut;
      CoinPackedVector bounds;
      if (firstAtUb>=0) {
	// fix to lb
	bounds.insert(firstAtUb,lower[firstAtUb]);
	newCut.setUbs(bounds);
      } else {
	// fix to ub
	bounds.insert(firstAtLb,upper[firstAtLb]);
	newCut.setLbs(bounds);
      }
#endif
#ifdef FORCE_CUTS
      /* The next statement  is very optional as it will
	 always add cut even if not violated.  If this is not done
	 you can get duplicates from heuristics and strong branching.
	 If you don't have it you may get some duplicates but it may be
	 much faster. On one problem I tried I got over 300 cuts each
	 of which was fully dense. 
	 This does not apply to column cuts, which might as well
	 always be forced on. 
	 
	 Even if this is done, you can still get duplicates
	 close together e.g. strong branching then heuristic
	 before cut generation entered.*/
      newCut.setEffectiveness(COIN_DBL_MAX);
#endif
      model_->makeGlobalCut(newCut);
    }
  }
#endif
}

/** Tuning class
 */

class CbcStrategyTuning : public CbcStrategy {
public:

  // Default Constructor 
  CbcStrategyTuning () {}

  // Copy constructor 
  CbcStrategyTuning ( const CbcStrategyTuning & rhs) : CbcStrategy(rhs) {}
   
  // Destructor 
  virtual ~CbcStrategyTuning () {}
  
  /// Clone
  virtual CbcStrategy * clone() const { return new CbcStrategyTuning(*this);}

  /// Setup cut generators
  virtual void setupCutGenerators(CbcModel & model);
  /// Setup heuristics
  virtual void setupHeuristics(CbcModel & model);
  /// Do printing stuff
  virtual void setupPrinting(CbcModel & ,int ) {}
  /// Other stuff e.g. strong branching
  virtual void setupOther(CbcModel & ) {}

protected:
  // Data
private:
  /// Illegal Assignment operator 
  CbcStrategyTuning & operator=(const CbcStrategyTuning& rhs);
};


// Setup cut generators
void 
CbcStrategyTuning::setupCutGenerators(CbcModel & model)
{
  int numberGenerators = model.numberCutGenerators();
  int iGenerator;
  for (iGenerator=0;iGenerator<numberGenerators;iGenerator++) {
    CbcCutGenerator * generator = model.cutGenerator(iGenerator);
    CglCutGenerator * cglGenerator = generator->generator();
    generator->generateTuning(stdout);
    /* Print options (in slightly garbled format) so user can change.
       Obviously in real use you would use generateCpp once and then
       cut and paste to modify what you need to modify.

       The output (which was designed for -cpp) gives the options you will 
       get if you do not modify the code.  So assuming you did not ask
       for more aggressive probing you will get a line like
       
       3  probing.setMaxElementsRoot(300);

       and you could modify that to 

       int numberColumns = model.getNumCols();
       probing->setMaxElementsRoot(100+numberColumns/2);
    */
    cglGenerator->generateCpp(stdout);
    CglProbing * probing = dynamic_cast<CglProbing *>(cglGenerator);
    if (probing) {
      // Could tune 
      continue;
    }
    CglGomory * gomory = dynamic_cast<CglGomory *>(cglGenerator);
    if (gomory) {
      // Could tune
      continue;
    }
    CglKnapsackCover * knapsackCover = dynamic_cast<CglKnapsackCover *>(cglGenerator);
    if (knapsackCover) {
      // Could tune
      continue;
    }
  }
}
// Setup heuristics
void 
CbcStrategyTuning::setupHeuristics(CbcModel & model)
{
  int numberHeuristics = model.numberHeuristics();
  int iHeuristic;
  for (iHeuristic=0;iHeuristic<numberHeuristics;iHeuristic++) {
    CbcHeuristic * heuristic = model.heuristic(iHeuristic);
    /* Print options (in slightly garbled format) so user can change.
       Obviously in real use you would use generateCpp once and then
       cut and paste to modify what you need to modify */
    heuristic->generateCpp(stdout);
    CbcRounding * rounding = dynamic_cast<CbcRounding *>(heuristic);
    if (rounding) {
      // Could tune
      rounding->setSeed(1234567);
      continue;
    }
  } 
}

int main (int argc, const char *argv[])
{

  OsiClpSolverInterface solver1;
  //#define USE_OSI_NAMES
#ifdef USE_OSI_NAMES
  // Say we are keeping names (a bit slower this way)
  solver1.setIntParam(OsiNameDiscipline,1);
#endif
  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName = "../../Data/Sample/p0033.mps";
  if (argc>=2) mpsFileName = argv[1];
  int numMpsReadErrors = solver1.readMps(mpsFileName.c_str(),"");
  assert(numMpsReadErrors==0);
  // Tell solver to return fast if presolve or initial solve infeasible
  solver1.getModelPtr()->setMoreSpecialOptions(3);

  /* Two ways of doing this depending on whether NEW_STYLE_SOLVER defined.
     So we need pointer to model.  Old way could use modelA. rather than model->
   */
#ifndef NEW_STYLE_SOLVER
#define NEW_STYLE_SOLVER 0
#endif
#if NEW_STYLE_SOLVER==0
  // Pass to Cbc initialize defaults 
  CbcModel modelA(solver1);
  CbcModel * model = &modelA;
  CbcMain0(modelA);
  // Event handler
  MyEventHandler3 eventHandler(model);
  eventHandler.setMaximumSolutions(300);
  model->passInEventHandler(&eventHandler);
  // Strategy so that we can tune generators and heuristics
  CbcStrategyTuning strategy;
  model->setStrategy(strategy);
  /* Now go into code for standalone solver
     Could copy arguments and add -quit at end to be safe
     but this will do
  */
  if (argc>2) {
    CbcMain1(argc-1,argv+1,modelA,callBack);
  } else {
    const char * argv2[]={"driver4","-solve","-quit"};
    CbcMain1(3,argv2,modelA,callBack);
  }
#else
  CbcSolver control(solver1);
  // initialize
  control.fillValuesInSolver();
  CbcModel * model = control.model();
  // Event handler
  MyEventHandler3 eventHandler(model);
  eventHandler.setMaximumSolutions(300);
  model->passInEventHandler(&eventHandler);
  // Strategy so that we can tune generators and heuristics
  CbcStrategyTuning strategy;
  model->setStrategy(strategy);
  control.solve (argc-1, argv+1, 1);
#endif
  // Solver was cloned so get current copy
  OsiSolverInterface * solver = model->solver();
  // Print solution if finished (could get from model->bestSolution() as well

  if (model->bestSolution()) {
    
    const double * solution = solver->getColSolution();
    
    int iColumn;
    int numberColumns = solver->getNumCols();
    std::cout<<std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14);
    
    std::cout<<"--------------------------------------"<<std::endl;
#ifdef USE_OSI_NAMES
    
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      double value=solution[iColumn];
      if (fabs(value)>1.0e-7&&solver->isInteger(iColumn)) 
	std::cout<<std::setw(6)<<iColumn<<" "<<std::setw(8)<<setiosflags(std::ios::left)<<solver->getColName(iColumn)
		 <<resetiosflags(std::ios::adjustfield)<<std::setw(14)<<" "<<value<<std::endl;
    }
#else
    // names may not be in current solver - use original
    
    for (iColumn=0;iColumn<numberColumns;iColumn++) {
      double value=solution[iColumn];
      if (fabs(value)>1.0e-7&&solver->isInteger(iColumn)) 
	std::cout<<std::setw(6)<<iColumn<<" "<<std::setw(8)<<setiosflags(std::ios::left)<<solver1.getModelPtr()->columnName(iColumn)
		 <<resetiosflags(std::ios::adjustfield)<<std::setw(14)<<" "<<value<<std::endl;
    }
#endif
    std::cout<<"--------------------------------------"<<std::endl;
  
    std::cout<<std::resetiosflags(std::ios::fixed|std::ios::showpoint|std::ios::scientific);
  } else {
    std::cout<<" No solution!"<<std::endl;
  }
  return 0;
}    
