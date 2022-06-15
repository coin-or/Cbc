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
  if (numberSolutions_== model_->getSolutionCount())
    return 0;
  // but don't until in tree ?
  if (model_->getNodeCount()<3)
    return 0;
  numberSolutions_ = model_->getSolutionCount();
  OsiSolverInterface * solver = model_->continuousSolver()->clone();
  // Add local branching cut.
  int ones_sum = 0;
  CoinPackedVector local_branching_cut;
  const double * best_solution = model_->bestSolution();
  int numberIntegers = model_->numberIntegers();
  const int * integer_indices = model_->integerVariable();
  
  for(int i=0;i<numberIntegers;i++) {
    int int_idx = integer_indices[i];
    double int_val = best_solution[int_idx];
    
    if(fabs(int_val - 1.0) <= model_->getIntegerTolerance())
      {
	++ones_sum;
	local_branching_cut.insert(int_idx, -1.0);
      }
    else if(fabs(int_val) <= model_->getIntegerTolerance())
      {
	local_branching_cut.insert(int_idx, 1.0);
      }
  }
  
  // k_ - neighbourhood size.
  solver->addRow(local_branching_cut, -COIN_DBL_MAX, k_ - ones_sum, "LocalBranchingCut"); 
  CbcModel model(*solver);
  // Changing settings and adding heuristics.
  model.setCutoffAsConstraint(true);
  model.setCutoff(model_->getCutoff());
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
    return model_->getSolutionCount()!=numberSolutions_;
  default:
    return false;
  }
}

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
static int callBack(CbcModel *model, int whereFrom)
{
  int returnCode = 0;
  switch (whereFrom) {
  case 1:
  case 2:
    if (!model->status() && model->secondaryStatus())
      returnCode = 1;
    break;
  case 3: {
    // Add in dummy heuristic
    CbcHeuristicUser heuristicUser(*model);
    heuristicUser.setHeuristicName("UserHeuristic");
    model->addHeuristic(&heuristicUser);
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
//#include "CbcEventHandler.hpp"
// May want to put back event handler - epecially if a few more events added ???

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
  //MyEventHandler3 eventHandler;
  //modelA.passInEventHandler(&eventHandler);
  /* Now go into code for standalone solver
     Could copy arguments and add -quit at end to be safe
     but this will do
  */
  if (argc > 2) {
    CbcMain1(argc - 1, argv + 1, modelA,callBack,parameterData);
  } else {
    const char *argv2[] = { "user", "-solve", "-quit" };
    CbcMain1(3, argv2, modelA,callBack,parameterData);
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
