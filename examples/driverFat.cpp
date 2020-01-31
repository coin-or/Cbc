// $Id: driverFat.cpp 1898 2013-04-09 18:06:04Z stefan $
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#include <cassert>
#include <iomanip>

#include "CoinPragma.hpp"
#include "CbcModel.hpp"
#include "CbcHeuristic.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcSolver.hpp"
#include "ClpSimplex.hpp"
#include "ClpFactorization.hpp"
#include "ClpPresolve.hpp"
#include "CoinSort.hpp"
#include "CoinHelperFunctions.hpp"

#include "CoinTime.hpp"
#include "CoinSignal.hpp"
/*
  This shows how to trap signals.
  This just traps ctrl-c and allows user to pause and then hit S or C
  In this simple version Stop may not be effective until a heuristic has exited
 */

static CbcModel *currentBranchModel = NULL;
extern "C" {
static void signal_handler(int whichSignal)
{
  int gotChar = 'X';
  while (toupper(gotChar) != 'S' && toupper(gotChar) != 'C') {
    // See what user wants to do
    fprintf(stderr, "Enter S to stop, C to continue:");
    gotChar = getchar();
  }
  if (currentBranchModel != NULL && toupper(gotChar) == 'S') {
    currentBranchModel->sayEventHappened(); // say why stopped
    if (currentBranchModel->heuristicModel())
      currentBranchModel->heuristicModel()->sayEventHappened();
  }
  return;
}
}
static CoinSighandler_t saveSignal = signal(SIGINT, signal_handler);
// Threshold below which use normal clp
static int rowsThreshold = -1;
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
    // set up signal trapping
    saveSignal = signal(SIGINT, signal_handler);
    currentBranchModel = model;
    /*******************************
        This tells code to be normal in heuristics i.e. smaller problems
        in practice you should probably make CoinMax(...,20000) or
        some such.
        You may wish to switch off strong branching and use priorities
        or something - as strong branching uses large model
      *******************************/
    rowsThreshold = model->getNumRows();
    // make heuristics do more nodes
    for (int i = 0; i < model->numberHeuristics(); i++) {
      CbcHeuristic *heuristic = model->heuristic(i);
      heuristic->setNumberNodes(5 * heuristic->numberNodes());
    }
    // could try doing feasibility after cuts?
    //model->setSpecialOptions(33554432|
    //		       model->specialOptions());
  } break;
  case 4: {
    // restore
    signal(SIGINT, saveSignal);
    currentBranchModel = NULL;
  }
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
#if 0
      int numberColumns = model_->getNumCols();
      const double * bestSolution = model_->bestSolution();
      assert (bestSolution);
      printf("value of solution is %g\n",model_->getObjValue());
      for (int i=0;i<numberColumns;i++) {
	if (fabs(bestSolution[i])>1.0e-8)
	  printf("%d %g\n",i,bestSolution[i]);
      }
#endif
      return noAction; // carry on
#endif
    } else {
      return noAction; // carry on
    }
  } else {
    return noAction; // carry on
  }
}
//#############################################################################

/**

    This is to allow the user to replace initialSolve and resolve
*/

class CbcSolverShortFat : public OsiClpSolverInterface {

public:
  //---------------------------------------------------------------------------
  /**@name Solve methods */
  //@{
  /// Solve initial LP relaxation
  virtual void initialSolve();

  /// Resolve an LP relaxation after problem modification
  virtual void resolve();

  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default Constructor
  CbcSolverShortFat();

  /// Clone
  virtual OsiSolverInterface *clone(bool CopyData = true) const;

  /// Copy constructor
  CbcSolverShortFat(const CbcSolverShortFat &);

  /// Assignment operator
  CbcSolverShortFat &operator=(const CbcSolverShortFat &rhs);

  /// Destructor
  virtual ~CbcSolverShortFat();

  //@}

  //---------------------------------------------------------------------------

private:
  /**@name Private member data */
  //@{
  //@}
};
static bool firstSolve = true;
void CbcSolverShortFat::initialSolve()
{
  ClpSimplex *model = getModelPtr();
  if (model->numberRows() < rowsThreshold) {
    OsiClpSolverInterface::initialSolve();
    return;
  }
#if LOGLEVEL > 1
  double time1 = CoinCpuTime();
#endif
  ClpPresolve pinfo;
#define INCREASE 6
#define PERTURB 5
#ifdef PERTURB
  bool externalPerturb = true;
#endif
  if (firstSolve) { // before preprocessing
    model = pinfo.presolvedModel(*model, 1.0e-8, false, 5, false);
  } else {
    externalPerturb = false;
    /* do initial factorization to get infeasibilities
       maybe fix all non basic
       initial fix checks if Osi already has basis - yes it has
       
     */
    setBasis(basis_, model);
  }
#define LOGLEVEL 0
#if LOGLEVEL < 3
  model->setLogLevel(0);
#endif
  int numberColumns = model->numberColumns();
  int originalNumberRows = model->numberRows();
  // change factorization frequency from 200
  model->setFactorizationFrequency(100 + model->numberRows() / 50);
  int numberIterations = 0;
#ifdef PERTURB
  double *objective = model->objective();
  double *saveObjective = NULL;
  if (externalPerturb) {
    saveObjective = CoinCopyOfArray(objective, numberColumns);
    double multiplier = 1.0;
    for (int i = 0; i < PERTURB; i++)
      multiplier *= 0.1;
    for (int i = 0; i < numberColumns; i++) {
      double value = CoinDrand48() * multiplier;
      // should be more sophisticated
      if (value > 1.0e-7) {
        if (objective[i] < 0.0)
          value = -value;
        objective[i] += value;
      }
    }
  }
#endif

  // We will need arrays to choose rows to add
  double *weight = new double[originalNumberRows];
  int *sort = new int[originalNumberRows];
  int numberSort = 0;
  char *take = new char[originalNumberRows];

  const double *rowLower = model->rowLower();
  const double *rowUpper = model->rowUpper();
  int iRow, iColumn;
  // Set up initial list
  numberSort = 0;
  int numberNonBasicRows = 0;
  for (iRow = 0; iRow < originalNumberRows; iRow++) {
    weight[iRow] = 1.123e50;
    if (model->getRowStatus(iRow) != ClpSimplex::basic) {
      sort[numberSort++] = iRow;
      weight[iRow] = -10.0;
      numberNonBasicRows++;
    } else if (rowLower[iRow] == rowUpper[iRow]) {
      sort[numberSort++] = iRow;
      weight[iRow] = 0.0;
    }
  }
  numberSort /= 2;
  numberSort = CoinMax(numberSort, numberNonBasicRows);
  // Just add this number of rows each time in small problem
  int smallNumberRows = 2 * numberColumns;
  smallNumberRows = CoinMin(smallNumberRows, originalNumberRows / 20);
  // and pad out with random rows
  double ratio = (static_cast< double >(smallNumberRows - numberSort)) / (static_cast< double >(originalNumberRows));
  bool primalInfeasible = false;
  if (firstSolve) {
    for (iRow = 0; iRow < originalNumberRows; iRow++) {
      if (weight[iRow] == 1.123e50 && CoinDrand48() < ratio)
        sort[numberSort++] = iRow;
    }
  }
  /* This is optional.
     The best thing to do is to miss out random rows and do a set which makes dual feasible.
     If that is not possible then make sure variables have bounds.
     
     One way that normally works is to automatically tighten bounds.
  */
  if (firstSolve) {
    // However for some we need to do anyway
    double *columnLower = model->columnLower();
    double *columnUpper = model->columnUpper();
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      columnLower[iColumn] = CoinMax(-1.0e6, columnLower[iColumn]);
      columnUpper[iColumn] = CoinMin(1.0e6, columnUpper[iColumn]);
    }
  }
  double *fullSolution = model->primalRowSolution();

  // Just do this number of passes
  int maxPass = 50;
  // And take out slack rows until this pass
  int takeOutPass = INCREASE;
  int iPass;

  const CoinBigIndex *start = model->clpMatrix()->getVectorStarts();
  const int *length = model->clpMatrix()->getVectorLengths();
  const int *row = model->clpMatrix()->getIndices();
  int *whichColumns = new int[numberColumns];
  for (int iRow = 0; iRow < originalNumberRows; iRow++)
    weight[iRow] = 0.0;
  for (int i = 0; i < numberSort; i++) {
    int iRow = sort[i];
    weight[iRow] = 1.0;
  }

  int numberSmallColumns = 0;
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    bool take = false;
    ;
    for (int j = start[iColumn]; j < start[iColumn] + length[iColumn]; j++) {
      int iRow = row[j];
      if (weight[iRow]) {
        take = true;
        break;
      }
    }
    if (take)
      whichColumns[numberSmallColumns++] = iColumn;
  }
#if LOGLEVEL > 1
  printf("%d rows, %d columns in initial problem\n", numberSort, numberSmallColumns);
#endif
  for (iPass = 0; iPass < maxPass; iPass++) {
#if LOGLEVEL > 2
    printf("Start of pass %d\n", iPass);
#endif
    // Cleaner this way
    std::sort(sort, sort + numberSort);
    // Create small problem
    ClpSimplex small(model, numberSort, sort, numberSmallColumns, whichColumns);
    small.setFactorizationFrequency(100 + numberSort / 200);
#ifndef PERTURB
    small.setPerturbation(50);
#else
    if (!externalPerturb)
      small.setPerturbation(50);
#endif
#if LOGLEVEL > 2
    small.setLogLevel(1);
#endif
    // A variation is to just do N iterations
    //if (iPass)
    //small.setMaximumIterations(100);
    // Solve
    //small.factorization()->messageLevel(8);
    small.dual();
    numberIterations += small.numberIterations();
    primalInfeasible = (small.status() == 1);
    if (primalInfeasible)
      break;
    bool dualInfeasible = (small.status() == 2);
    // move solution back
    double *solution = model->primalColumnSolution();
    const double *smallSolution = small.primalColumnSolution();
    for (int j = 0; j < numberSmallColumns; j++) {
      iColumn = whichColumns[j];
      solution[iColumn] = smallSolution[j];
      model->setColumnStatus(iColumn, small.getColumnStatus(j));
    }
    for (iRow = 0; iRow < numberSort; iRow++) {
      int kRow = sort[iRow];
      model->setRowStatus(kRow, small.getRowStatus(iRow));
    }
    // compute full solution
    memset(fullSolution, 0, originalNumberRows * sizeof(double));
    model->clpMatrix()->times(1.0, model->primalColumnSolution(), fullSolution);
    if (iPass != maxPass - 1) {
      // Mark row as not looked at
      for (iRow = 0; iRow < originalNumberRows; iRow++)
        weight[iRow] = 1.123e50;
      // Look at rows already in small problem
      int iSort;
      int numberDropped = 0;
      int numberKept = 0;
      int numberBinding = 0;
      int numberInfeasibilities = 0;
      double sumInfeasibilities = 0.0;
      for (iSort = 0; iSort < numberSort; iSort++) {
        iRow = sort[iSort];
        if (model->getRowStatus(iRow) == ClpSimplex::basic) {
          // Basic - we can get rid of if early on
          if (iPass < takeOutPass && !dualInfeasible) {
            // may have hit max iterations so check
            double infeasibility = CoinMax(fullSolution[iRow] - rowUpper[iRow],
              rowLower[iRow] - fullSolution[iRow]);
            weight[iRow] = -infeasibility;
            if (infeasibility > 1.0e-8) {
              numberInfeasibilities++;
              sumInfeasibilities += infeasibility;
            } else {
              weight[iRow] = 1.0;
              numberDropped++;
            }
          } else {
            // keep
            weight[iRow] = -1.0e40;
            numberKept++;
          }
        } else {
          // keep
          weight[iRow] = -1.0e50;
          numberKept++;
          numberBinding++;
        }
      }
      // Now rest
      for (iRow = 0; iRow < originalNumberRows; iRow++) {
        sort[iRow] = iRow;
        if (weight[iRow] == 1.123e50) {
          // not looked at yet
          double infeasibility = CoinMax(fullSolution[iRow] - rowUpper[iRow],
            rowLower[iRow] - fullSolution[iRow]);
          weight[iRow] = -infeasibility;
          if (infeasibility > 1.0e-8) {
            numberInfeasibilities++;
            sumInfeasibilities += infeasibility;
          }
        }
      }
      // sort
      CoinSort_2(weight, weight + originalNumberRows, sort);
      numberSort = CoinMin(originalNumberRows, smallNumberRows + numberKept);
      memset(take, 0, originalNumberRows);
      for (iRow = 0; iRow < numberSort; iRow++)
        take[sort[iRow]] = 1;
      numberSmallColumns = 0;
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        int n = 0;
        for (int j = start[iColumn]; j < start[iColumn] + length[iColumn]; j++) {
          int iRow = row[j];
          if (take[iRow])
            n++;
        }
        if (n)
          whichColumns[numberSmallColumns++] = iColumn;
      }
#if LOGLEVEL > 1
      printf("%d rows binding, %d rows kept, %d rows dropped - new size %d rows, %d columns\n",
        numberBinding, numberKept, numberDropped, numberSort, numberSmallColumns);
      printf("%d rows are infeasible - sum is %g\n", numberInfeasibilities,
        sumInfeasibilities);
#endif
      if (!numberInfeasibilities) {
#if LOGLEVEL > 1
        printf("Exiting as looks optimal\n");
#endif
        break;
      }
      numberInfeasibilities = 0;
      sumInfeasibilities = 0.0;
      for (iSort = 0; iSort < numberSort; iSort++) {
        if (weight[iSort] > -1.0e30 && weight[iSort] < -1.0e-8) {
          numberInfeasibilities++;
          sumInfeasibilities += -weight[iSort];
        }
      }
#if LOGLEVEL > 1
      printf("in small model %d rows are infeasible - sum is %g\n", numberInfeasibilities,
        sumInfeasibilities);
#endif
    }
  }
  delete[] weight;
  delete[] sort;
  delete[] whichColumns;
  delete[] take;
#ifdef PERTURB
  if (externalPerturb) {
    memcpy(objective, saveObjective, numberColumns * sizeof(double));
    delete[] saveObjective;
  }
  model->setPerturbation(50);
#endif
  if (!primalInfeasible) {
    model->primal(1);
    numberIterations += model->numberIterations();
  } else {
    model->setProblemStatus(1);
  }
  model->setNumberIterations(numberIterations);
  if (firstSolve) {
    pinfo.postsolve(true);
    model = getModelPtr();
    model->primal(1);
    firstSolve = false;
  }
  basis_ = getBasis(model);
#if LOGLEVEL > 1
  printf("solve took %g seconds and %d iterations\n", CoinCpuTime() - time1,
    numberIterations);
#endif
}

//-----------------------------------------------------------------------------
void CbcSolverShortFat::resolve()
{
  ClpSimplex *model = getModelPtr();
  if (model->numberRows() < rowsThreshold) {
    OsiClpSolverInterface::resolve();
  } else {
#if LOGLEVEL > 1
    printf("resolve\n");
#endif
    initialSolve();
  }
  return;
}

//#############################################################################
// Constructors, destructors clone and assignment
//#############################################################################

//-------------------------------------------------------------------
// Default Constructor
//-------------------------------------------------------------------
CbcSolverShortFat::CbcSolverShortFat()
  : OsiClpSolverInterface()
{
}

//-------------------------------------------------------------------
// Clone
//-------------------------------------------------------------------
OsiSolverInterface *
CbcSolverShortFat::clone(bool CopyData) const
{
  if (CopyData) {
    return new CbcSolverShortFat(*this);
  } else {
    printf("warning CbcSolveUser clone with copyData false\n");
    return new CbcSolverShortFat();
  }
}

//-------------------------------------------------------------------
// Copy constructor
//-------------------------------------------------------------------
CbcSolverShortFat::CbcSolverShortFat(
  const CbcSolverShortFat &rhs)
  : OsiClpSolverInterface(rhs)
{
}

//-------------------------------------------------------------------
// Destructor
//-------------------------------------------------------------------
CbcSolverShortFat::~CbcSolverShortFat()
{
}

//-------------------------------------------------------------------
// Assignment operator
//-------------------------------------------------------------------
CbcSolverShortFat &
CbcSolverShortFat::operator=(const CbcSolverShortFat &rhs)
{
  if (this != &rhs) {
    OsiClpSolverInterface::operator=(rhs);
  }
  return *this;
}
//-------------------------------------------------------------------

int main(int argc, const char *argv[])
{

  CbcSolverShortFat solver1;
  // Read in model using argv[1]
  // and assert that it is a clean model
  std::string mpsFileName;
  if (argc < 2) {
    fprintf(stderr, "Do not know where to find input file.\n");
    exit(1);
  }
  if (argc >= 2)
    mpsFileName = argv[1];
  int numMpsReadErrors;
  if (!strstr(mpsFileName.c_str(), ".lp"))
    numMpsReadErrors = solver1.readMps(mpsFileName.c_str(), "");
  else
    numMpsReadErrors = solver1.readLp(mpsFileName.c_str(), 1.0e-12);
  if (numMpsReadErrors != 0) {
    printf("%d errors reading MPS file\n", numMpsReadErrors);
    return numMpsReadErrors;
  }
  // Tell solver to return fast if presolve or initial solve infeasible
  solver1.getModelPtr()->setMoreSpecialOptions(3);

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
    const char *argv2[] = { "driverFat", "-solve", "-quit" };
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
    // names may not be in current solver - use original

    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      double value = solution[iColumn];
      if (fabs(value) > 1.0e-7 && solver->isInteger(iColumn))
        std::cout << std::setw(6) << iColumn << " " << std::setw(8) << setiosflags(std::ios::left) << solver1.getModelPtr()->columnName(iColumn)
                  << resetiosflags(std::ios::adjustfield) << std::setw(14) << " " << value << std::endl;
    }
    std::cout << "--------------------------------------" << std::endl;

    std::cout << std::resetiosflags(std::ios::fixed | std::ios::showpoint | std::ios::scientific);
  } else {
    std::cout << " No solution!" << std::endl;
  }
  return 0;
}
