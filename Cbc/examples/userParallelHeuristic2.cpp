// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip>

#include "CoinPragma.hpp"
#include "CbcModel.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcSolver.hpp"
#include "CbcHeuristic.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristicRENS.hpp"
#include "CbcHeuristicVND.hpp"

#include "CoinTime.hpp"
#define CBC_THREAD
// Standard build does not copy CbcThread.hpp to correct place
// use below if fixed
//#include "CbcThread.hpp"
// otherwise tweak this
#include "../../../../Cbc/Cbc/src/CbcThread.hpp" // see comments here
//#############################################################################

/************************************************************************

This main program shows how to take advantage of the standalone cbc in your program,
while trapping events
First it reads in an integer model from an mps file
Then it initializes the integer model with cbc defaults
Then it calls CbcMain1 passing all parameters apart from first but with callBack to modify stuff
Finally it could print solution

************************************************************************/
// Save useful CbcModel
static CbcModel * userModel = NULL;
// Restarts cause chaos - can use this instead of model_
static CbcModel * realModel =NULL;
class CbcHeuristicUser : public CbcHeuristic {
public:
  // Default Constructor
  CbcHeuristicUser();

  /* Constructor with model - assumed before cuts
    */
  CbcHeuristicUser(CbcModel &model);

  // Copy constructor
  CbcHeuristicUser(const CbcHeuristicUser &);

  // Destructor
  ~CbcHeuristicUser();

  /// Clone
  virtual CbcHeuristic *clone() const;

  /// Assignment operator
  CbcHeuristicUser &operator=(const CbcHeuristicUser &rhs);

  /// Resets stuff if model changes
  virtual void resetModel(CbcModel *model);

  /// update model (This is needed if cliques update matrix etc)
  virtual void setModel(CbcModel *model);

  using CbcHeuristic::solution;
  /** returns 0 if no solution, 1 if valid solution.
        Sets solution values if good, sets objective value (only if good)
  */
  virtual int solution(double &objectiveValue,
		       double *newSolution);

  virtual bool shouldHeurRun(int whereFrom);
protected:
  // Data
  // Number of solutions found by model
  int numberSolutions_;
  // -1 - off, 0 - startup, 1 - do heuristic
  int actions_;
  // Add any stuff needed by heuristic
  int k_;
  int nodePick_;
};
// Default Constructor
CbcHeuristicUser::CbcHeuristicUser()
  : CbcHeuristic()
{
  numberSolutions_ = 0;
  actions_= -1;
  k_ = 0;
  nodePick_ = 3;
}
// Constructor with model - assumed before cuts

CbcHeuristicUser::CbcHeuristicUser(CbcModel &model)
  : CbcHeuristic(model)
{
  numberSolutions_ = 0;
  actions_ = 0;
  k_ = 5;
  nodePick_ = 3;
}

// Destructor
CbcHeuristicUser::~CbcHeuristicUser()
{
}

// Clone
CbcHeuristic *
CbcHeuristicUser::clone() const
{
  return new CbcHeuristicUser(*this);
}

// Assignment operator
CbcHeuristicUser &
CbcHeuristicUser::operator=(const CbcHeuristicUser &rhs)
{
  if (this != &rhs) {
    CbcHeuristic::operator=(rhs);
    numberSolutions_ = rhs.numberSolutions_;
    actions_ = rhs.actions_;
    k_ = rhs.k_;
    nodePick_ = rhs.nodePick_;
  }
  return *this;
}

// Copy constructor
CbcHeuristicUser::CbcHeuristicUser(const CbcHeuristicUser &rhs)
  : CbcHeuristic(rhs)
{
  numberSolutions_ = rhs.numberSolutions_;
  actions_ = rhs.actions_;
  k_ = rhs.k_;
  nodePick_ = rhs.nodePick_;
}
// Resets stuff if model changes
void CbcHeuristicUser::resetModel(CbcModel *)
{
}
int CbcHeuristicUser::solution(double &solutionValue,
  double *betterSolution)
{
  if (actions_ <0 ) {
    return 0;
  } else if (actions_ == 0) {
    // see if wanted
    int numberThreads = model_->getNumberThreads();
    if (!numberThreads) {
      actions_ = 1;
    } else {
      // choose a thread in middle?? (don't want first)
      int wantedThread = numberThreads/2;
      /* This is really clumsy - may be better to modify Cbc a tiny bit
	 so that it is easy to pick up thread number */
      CbcThread * masterThread = model_->masterThread();
      if (!masterThread) 
	return 0;  // not going yet
      CbcModel * baseModel = masterThread->baseModel();
      CbcBaseModel * master = baseModel->master();
      CbcModel * wanted = master->model(wantedThread);
      if (model_ == wanted) {
	actions_ = 1;
      } else {
	actions_ = -1;
	return 0;
      }
    }
  }
  // if model restarts, some heuristics may not work
  if (model_->numberIntegers()!=realModel->numberIntegers()) {
    //printf("restartx? %d columns, %d rows, %d,%d ints\n",model_->getNumCols(),model_->getNumRows(),
    //	   model_->numberIntegers(),realModel->numberIntegers() );
    static int xxxxxx=0;
    if (!xxxxxx)
      printf("restart models not yet reset\n");
    xxxxxx++;
    return 0;
  }
  if (numberSolutions_== realModel->getSolutionCount())
    return 0;
  // but don't until in tree ?
  if (realModel->getNodeCount()<3)
    return 0;
  numberSolutions_ = realModel->getSolutionCount();
  if (!userModel) {
    // must have restarted model!
    return 0;
  }
  // Can either use LP solver or with root cuts
  // if LP solver then can use model_->
  OsiSolverInterface * solver = userModel->continuousSolver()->clone();
  //OsiSolverInterface * solver = userModel->solver()->clone();
  // Add local branching cut.
  int ones_sum = 0;
  CoinPackedVector local_branching_cut;
  const double * best_solution = realModel->bestSolution();
  int numberIntegers = realModel->numberIntegers();
  const int * integer_indices = realModel->integerVariable();
  
  for(int i=0;i<numberIntegers;i++) {
    int int_idx = integer_indices[i];
    double int_val = best_solution[int_idx];
    
    if(fabs(int_val - 1.0) <= realModel->getIntegerTolerance())
      {
	++ones_sum;
	local_branching_cut.insert(int_idx, -1.0);
      }
    else if(fabs(int_val) <= realModel->getIntegerTolerance())
      {
	local_branching_cut.insert(int_idx, 1.0);
      }
  }
  
  // k_ - neighbourhood size.
  solver->addRow(local_branching_cut, -COIN_DBL_MAX, k_ - ones_sum, "LocalBranchingCut"); 
  CbcModel model(*solver);
  // Changing settings and adding heuristics.
  model.setCutoffAsConstraint(true);
  model.setCutoff(realModel->getCutoff());
  model.setPreferredWay(1);
  model.setMaximumNodes(150);
  model.setLogLevel(0); // or 1 when starting off
  CbcHeuristicVND vnd(model);
  vnd.setHeuristicName("VND");
  vnd.setFractionSmall(0.4);
  vnd.setFeasibilityPumpOptions(1008003);
  int nodes[] = { -2, 50, 50, 50, 200, 1000, 10000 };
  vnd.setNumberNodes(nodes[nodePick_]);
  model.addHeuristic(&vnd);
  CbcHeuristicRENS rens(model);
  rens.setHeuristicName("RENSdj");
  rens.setFractionSmall(0.6 /*3.4*/);
  rens.setFeasibilityPumpOptions(3);
  rens.setNumberNodes(10);
  rens.setWhereFrom(4 * 256 + 4 * 1);
  rens.setWhen(2);
  rens.setRensType(1 + 16);
  model.addHeuristic(&rens);
  CbcHeuristicFPump fpump(model);
  fpump.setFractionSmall(0.5);
  fpump.setMaximumPasses(5);
  fpump.setFeasibilityPumpOptions(30);
  fpump.setWhen(13);
  fpump.setHeuristicName("feasibility pump");
  model.addHeuristic(&fpump);
  CbcRounding rounding(model);
  rounding.setHeuristicName("rounding");
  model.addHeuristic(&rounding);
  printf("starting user with cutoff of %g\n",realModel->getCutoff());
  model.branchAndBound();
  if (model.bestSolution()) {
    solutionValue = model.getMinimizationObjValue();
    // solution may be worse than one found now but ...
    printf("user_heuristic found solution of %g\n",solutionValue); // debug
    int numberColumns = solver->getNumCols();
    memcpy(betterSolution, model.bestSolution(),
	   numberColumns*sizeof(double));
    delete solver;
    return 1;
  } else {
    printf("user_failed\n"); // debug
    delete solver;
    return 0;
  }
}
// update model
void CbcHeuristicUser::setModel(CbcModel *model)
{
  model_ = model;
}
bool CbcHeuristicUser::shouldHeurRun(int whereFrom)
{
  switch (actions_) {
  case 0:
    return (model_->specialOptions()&2048)==0; //not in submodel
  case 1:
    return realModel->getSolutionCount()!=numberSolutions_;
  default:
    return false;
  }
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

  /**@name Constructors, destructor etc*/
  //@{
  /** Default constructor. */
  MyEventHandler3();
  /// Constructor with pointer to model (redundant as setEventHandler does)
  MyEventHandler3(CbcModel *model);
  /** Destructor */
  virtual ~MyEventHandler3();
  /** The copy constructor. */
  MyEventHandler3(const MyEventHandler3 &rhs);
  /// Assignment
  MyEventHandler3 &operator=(const MyEventHandler3 &rhs);
  /// Clone
  virtual CbcEventHandler *clone() const;
  /*! \brief Set model. */
  //@}

protected:
  // data goes here
};
//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
MyEventHandler3::MyEventHandler3()
  : CbcEventHandler()
{
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
MyEventHandler3::MyEventHandler3(const MyEventHandler3 &rhs)
  : CbcEventHandler(rhs)
{
}

// Constructor with pointer to model
MyEventHandler3::MyEventHandler3(CbcModel *model)
  : CbcEventHandler(model)
{
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
MyEventHandler3::~MyEventHandler3()
{
}

//----------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
MyEventHandler3 &
MyEventHandler3::operator=(const MyEventHandler3 &rhs)
{
  if (this != &rhs) {
    CbcEventHandler::operator=(rhs);
  }
  return *this;
}
//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
CbcEventHandler *MyEventHandler3::clone() const
{
  return new MyEventHandler3(*this);
}

CbcEventHandler::CbcAction
MyEventHandler3::event(CbcEvent whichEvent)
{
  // If in sub tree carry on (need to miss out when in heuristics)
  if ((model_->specialOptions()&2048)!=0) {
    return noAction; // carry on as in subproblem
  } else if (whichEvent == startUp) {
    // check heuristic not there (restart)
    int numberHeuristics = model_->numberHeuristics();
    int found;
    for (found=0;found<numberHeuristics;found++) {
      if (dynamic_cast<CbcHeuristicUser *>(model_->heuristic(found))) {
	// already there
	break;
      }
    }
    if (found==numberHeuristics) {
      // Add in dummy heuristic
      CbcHeuristicUser heuristicUser(*model_);
      heuristicUser.setHeuristicName("UserHeuristic");
      model_->addHeuristic(&heuristicUser);
    } else {
      model_->heuristic(found)->setModel(model_);
      delete userModel; 
      userModel = NULL;
    }
    realModel = model_;
    //printf("startup %p\n",model_);
  } else if (whichEvent == afterRootCuts) {
    // save model
    delete userModel; // in case restart
    userModel = new CbcModel(*model_);
    //printf("afterRootCuts %p\n",model_);
  } else if (whichEvent == endSearch) {
    //printf("endSearch %p\n",model_);
    // free model
    delete userModel;
    userModel = NULL;
  }
  return noAction; // carry on
}

static int dummyCallBack(CbcModel * /*model*/, int /*whereFrom*/)
{
  return 0;
}
int main(int argc, const char *argv[])
{

  OsiClpSolverInterface solver1;
#define USE_OSI_NAMES
#ifdef USE_OSI_NAMES
  // Say we are keeping names (a bit slower this way)
  solver1.setIntParam(OsiNameDiscipline, 1);
#endif
  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName;
  if (argc < 2) {
    fprintf(stderr, "Do not know where to find MPS file.\n");
    exit(1);
  } else {
    printf("Running"); 
    for (int i=1;i<argc;i++)
      printf("%s ",argv[i]);
    printf("\n");
    mpsFileName = argv[1];
  }
  int numMpsReadErrors = solver1.readMps(mpsFileName.c_str(), "");
  if (numMpsReadErrors != 0) {
    printf("%d errors reading MPS file\n", numMpsReadErrors);
    return numMpsReadErrors;
  }
  // Tell solver to return fast if presolve or initial solve infeasible
  solver1.getModelPtr()->setMoreSpecialOptions(3);

  // Pass to Cbc initialize defaults
  CbcModel modelA(solver1);
  CbcSolverUsefulData parameterData;
  CbcMain0(modelA,parameterData);
  // Event handler
  MyEventHandler3 eventHandler;
  modelA.passInEventHandler(&eventHandler);
  /* Now go into code for standalone solver
     Could copy arguments and add -quit at end to be safe
     but this will do
  */
  if (argc > 2) {
    CbcMain1(argc - 1, argv + 1, modelA,dummyCallBack,parameterData);
  } else {
    const char *argv2[] = { "user", "-solve", "-quit" };
    CbcMain1(3, argv2, modelA,dummyCallBack,parameterData);
  }
#ifdef PRINT_SOLUTION
  // Print solution if finished
  if (modelA.bestSolution()) {
    const double *solution = modelA.bestSolution();
    int numberColumns = modelA.getNumCols();
    std::cout << std::setiosflags(std::ios::fixed | std::ios::showpoint) << std::setw(14);

    std::cout << "--------------------------------------" << std::endl;
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = solution[iColumn];
      if (fabs(value) > 1.0e-7 && solver->isInteger(iColumn))
        std::cout << std::setw(6) << iColumn << " " << std::setw(8) << setiosflags(std::ios::left) << solver->getColName(iColumn)
                  << resetiosflags(std::ios::adjustfield) << std::setw(14) << " " << value << std::endl;
    }
    std::cout << "--------------------------------------" << std::endl;

    std::cout << std::resetiosflags(std::ios::fixed | std::ios::showpoint | std::ios::scientific);
  } else {
    std::cout << " No solution!" << std::endl;
  }
#else
  if (modelA.bestSolution()) {
    printf("Optimal solution of %g\n",modelA.getMinimizationObjValue());
  } else {
    printf(" No solution!\n");
  }
#endif
  return 0;
}
