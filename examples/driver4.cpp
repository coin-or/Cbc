// $Id$
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip>

#include "CoinPragma.hpp"
#include "CbcModel.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcSolver.hpp"
#include "CbcHeuristicDiveCoefficient.hpp"
#include "CbcHeuristicDiveFractional.hpp"
#include "CbcHeuristicDiveGuided.hpp"
#include "CbcHeuristicDiveVectorLength.hpp"
#include "CbcHeuristicDivePseudoCost.hpp"
#include "CbcHeuristicDiveLineSearch.hpp"

#include "CoinTime.hpp"

//#############################################################################

/************************************************************************

This main program shows how to take advantage of the standalone cbc in your program,
while still making major modifications.
First it reads in an integer model from an mps file
Then it initializes the integer model with cbc defaults
Then it calls CbcMain1 passing all parameters apart from first but with callBack to modify stuff
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
    // Add in some diving heuristics with different options
    CbcHeuristicDiveCoefficient heuristicDC(*model);
    heuristicDC.setHeuristicName("DiveCoefficient");
    // do if no solution
    heuristicDC.setWhen(3);
    // 150 passes and fix general integers
    heuristicDC.setMaxIterations(151);
    // make sure can do as many simplex iterations as wanted
    heuristicDC.setMaxSimplexIterations(COIN_INT_MAX);
    heuristicDC.setMaxSimplexIterationsAtRoot(COIN_INT_MAX);
    model->addHeuristic(&heuristicDC);
    CbcHeuristicDiveFractional heuristicDF(*model);
    heuristicDF.setHeuristicName("DiveFractional");
    // do if no solution
    heuristicDF.setWhen(3);
    // 150 passes and don't fix general integers
    heuristicDF.setMaxIterations(150);
    // make sure can do as many simplex iterations as wanted
    heuristicDF.setMaxSimplexIterations(COIN_INT_MAX);
    heuristicDF.setMaxSimplexIterationsAtRoot(COIN_INT_MAX);
    model->addHeuristic(&heuristicDF);
    CbcHeuristicDiveGuided heuristicDG(*model);
    heuristicDG.setHeuristicName("DiveGuided");
    // do if no solution
    heuristicDG.setWhen(3);
    // 200 passes and fix general integers
    heuristicDG.setMaxIterations(201);
    // make sure can do as many simplex iterations as wanted
    heuristicDG.setMaxSimplexIterations(COIN_INT_MAX);
    heuristicDG.setMaxSimplexIterationsAtRoot(COIN_INT_MAX);
    model->addHeuristic(&heuristicDG);
    CbcHeuristicDiveVectorLength heuristicDV(*model);
    heuristicDV.setHeuristicName("DiveVectorLength");
    // do if no solution
    heuristicDV.setWhen(3);
    // 150 passes and fix general integers
    heuristicDV.setMaxIterations(151);
    // make sure can do as many simplex iterations as wanted
    heuristicDV.setMaxSimplexIterations(COIN_INT_MAX);
    heuristicDV.setMaxSimplexIterationsAtRoot(COIN_INT_MAX);
    model->addHeuristic(&heuristicDV);
    // Second version!
    CbcHeuristicDiveVectorLength heuristicDV2(*model);
    heuristicDV2.setHeuristicName("DiveVectorLength");
    // do if no solution
    heuristicDV2.setWhen(3);
    // 300 passes and don't fix general integers
    heuristicDV2.setMaxIterations(300);
    // fix fewer
    heuristicDV2.setPercentageToFix(0.05);
    // make sure can do as many simplex iterations as wanted
    heuristicDV2.setMaxSimplexIterations(COIN_INT_MAX);
    heuristicDV2.setMaxSimplexIterationsAtRoot(COIN_INT_MAX);
    model->addHeuristic(&heuristicDV2);
    CbcHeuristicDivePseudoCost heuristicDP(*model);
    heuristicDP.setHeuristicName("DivePseudoCost");
    // do if no solution
    heuristicDP.setWhen(3);
    // 100 passes and don't fix general integers
    heuristicDP.setMaxIterations(100);
    // make sure can do as many simplex iterations as wanted
    heuristicDP.setMaxSimplexIterations(COIN_INT_MAX);
    heuristicDP.setMaxSimplexIterationsAtRoot(COIN_INT_MAX);
    model->addHeuristic(&heuristicDP);
    CbcHeuristicDiveLineSearch heuristicDL(*model);
    heuristicDL.setHeuristicName("DiveLineSearch");
    // do if no solution
    heuristicDL.setWhen(3);
    // 150 passes and fix general integers
    heuristicDL.setMaxIterations(151);
    // make sure can do as many simplex iterations as wanted
    heuristicDL.setMaxSimplexIterations(COIN_INT_MAX);
    heuristicDL.setMaxSimplexIterationsAtRoot(COIN_INT_MAX);
    model->addHeuristic(&heuristicDL);
    //CbcCompareUser compare;
    //model->setNodeComparison(compare);
  } break;
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
  // If in sub tree carry on
  if (!model_->parentModel()) {
    if (whichEvent == solution || whichEvent == heuristicSolution) {
#ifdef STOP_EARLY
      return stop; // say finished
#else
      // If preprocessing was done solution will be to processed model
      int numberColumns = model_->getNumCols();
      const double *bestSolution = model_->bestSolution();
      assert(bestSolution);
      printf("event handler received a solution with cost %g\n", model_->getObjValue());
      // printing solution contents
      //for (int i = 0; i < numberColumns; i++) {
      //  if (fabs(bestSolution[i]) > 1.0e-8)
      //    printf("%d %g\n", i, bestSolution[i]);
      //}
      return noAction; // carry on
#endif
    } else {
      return noAction; // carry on
    }
  } else {
    return noAction; // carry on
  }
}

int main(int argc, const char *argv[])
{

  OsiClpSolverInterface solver1;
  //#define USE_OSI_NAMES
#ifdef USE_OSI_NAMES
  // Say we are keeping names (a bit slower this way)
  solver1.setIntParam(OsiNameDiscipline, 1);
#endif
  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName;
#if defined(SAMPLEDIR)
  mpsFileName = SAMPLEDIR "/p0033.mps";
#else
  if (argc < 2) {
    fprintf(stderr, "Do not know where to find sample MPS files.\n");
    exit(1);
  }
#endif
  if (argc >= 2)
    mpsFileName = argv[1];
  int numMpsReadErrors = solver1.readMps(mpsFileName.c_str(), "");
  if (numMpsReadErrors != 0) {
    printf("%d errors reading MPS file\n", numMpsReadErrors);
    return numMpsReadErrors;
  }
  // Tell solver to return fast if presolve or initial solve infeasible
  solver1.getModelPtr()->setMoreSpecialOptions(3);

  /* Two ways of doing this depending on whether NEW_STYLE_SOLVER defined.
     So we need pointer to model.  Old way could use modelA. rather than model->
   */
  // Messy code below copied from CbcSolver.cpp
  // Pass to Cbc initialize defaults
  CbcModel modelA(solver1);
  CbcModel *model = &modelA;
  CbcMain0(modelA);
  // Event handler
  MyEventHandler3 eventHandler;
  model->passInEventHandler(&eventHandler);
  /* Now go into code for standalone solver
     Could copy arguments and add -quit at end to be safe
     but this will do
  */
  if (argc > 2) {
    CbcMain1(argc - 1, argv + 1, modelA, callBack);
  } else {
    const char *argv2[] = { "driver4", "-solve", "-quit" };
    CbcMain1(3, argv2, modelA, callBack);
  }
  // Solver was cloned so get current copy
  OsiSolverInterface *solver = model->solver();
  // Print solution if finished (could get from model->bestSolution() as well

  if (model->bestSolution()) {

    const double *solution = solver->getColSolution();

    int iColumn;
    int numberColumns = solver->getNumCols();
    std::cout << std::setiosflags(std::ios::fixed | std::ios::showpoint) << std::setw(14);

    std::cout << "--------------------------------------" << std::endl;
#ifdef USE_OSI_NAMES

    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = solution[iColumn];
      if (fabs(value) > 1.0e-7 && solver->isInteger(iColumn))
        std::cout << std::setw(6) << iColumn << " " << std::setw(8) << setiosflags(std::ios::left) << solver->getColName(iColumn)
                  << resetiosflags(std::ios::adjustfield) << std::setw(14) << " " << value << std::endl;
    }
#else
    // names may not be in current solver - use original

    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = solution[iColumn];
      if (fabs(value) > 1.0e-7 && solver->isInteger(iColumn))
        std::cout << std::setw(6) << iColumn << " " << std::setw(8) << setiosflags(std::ios::left) << solver1.getModelPtr()->columnName(iColumn)
                  << resetiosflags(std::ios::adjustfield) << std::setw(14) << " " << value << std::endl;
    }
#endif
    std::cout << "--------------------------------------" << std::endl;

    std::cout << std::resetiosflags(std::ios::fixed | std::ios::showpoint | std::ios::scientific);
  } else {
    std::cout << " No solution!" << std::endl;
  }
  return 0;
}
