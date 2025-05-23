// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file CbcSolver.cpp
    \brief Second level routines for the cbc stand-alone solver.
*/

// Need these up front to define symbols for other imports
#include "CoinUtilsConfig.h"
#include "ClpConfig.h"
#include "CbcConfig.h"

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sstream>
#include <iostream>

#if defined(NEW_DEBUG_AND_FILL) || defined(CLP_MALLOC_STATISTICS)
#include <exception>
#include <malloc.h>
#include <new>
#endif
#if defined(HAVE_SIGNAL_H) && defined(HAVE_EXECINFO_H)
#include <execinfo.h>
#include <signal.h>
void CbcCrashHandler(int sig);
#endif

#ifdef COINUTILS_HAS_GLPK
#include "glpk.h"
#endif
#ifdef DMALLOC
#include "dmalloc.h"
#endif

//###########################################################################
// COIN-OR headers (do we really need all these?)
// Should some be in CbcSolver.hpp?
//###########################################################################

#include "CoinHelperFunctions.hpp"
#include "CoinModel.hpp"
#include "CoinMpsIO.hpp"
#include "CoinParam.hpp"
#include "CoinPragma.hpp"
#include "CoinSignal.hpp"
#include "CoinWarmStartBasis.hpp"

#include "ClpDualRowDantzig.hpp"
#include "ClpDualRowSteepest.hpp"
#include "ClpFactorization.hpp"
#include "ClpLinearObjective.hpp"
#include "ClpNetworkMatrix.hpp"
#include "ClpPEDualRowDantzig.hpp"
#include "ClpPEDualRowSteepest.hpp"
#include "ClpPEPrimalColumnDantzig.hpp"
#include "ClpPEPrimalColumnSteepest.hpp"
#include "ClpPackedMatrix.hpp"
#include "ClpParamUtils.hpp"
#include "ClpPlusMinusOneMatrix.hpp"
#include "ClpPresolve.hpp"
#include "ClpPrimalColumnDantzig.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "ClpQuadraticObjective.hpp"
#include "ClpSimplex.hpp"
#include "ClpSimplexOther.hpp"
#include "ClpSolve.hpp"
#include "CoinTime.hpp"

#include "OsiAuxInfo.hpp"
#include "OsiChooseVariable.hpp"
#include "OsiClpSolverInterface.hpp"
#include "OsiColCut.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiSolverInterface.hpp"
/* CBC_OTHER_SOLVER == 1 is cplex. */
#if CBC_OTHER_SOLVER == 1
#ifndef CBC_HAS_OSICPX
#error "Configuration did not detect OsiCpx installation."
#else
#include "OsiCpxSolverInterface.hpp"
#endif
#endif

#include "CglBKClique.hpp"
#include "CglClique.hpp"
#include "CglCliqueStrengthening.hpp"
#include "CglCutGenerator.hpp"
#include "CglDuplicateRow.hpp"
#include "CglFlowCover.hpp"
#include "CglGMI.hpp"
#include "CglGomory.hpp"
#include "CglKnapsackCover.hpp"
#include "CglLandP.hpp"
#include "CglMessage.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglOddWheel.hpp"
#include "CglPreProcess.hpp"
#include "CglProbing.hpp"
#include "CglRedSplit.hpp"
#include "CglRedSplit2.hpp"
#include "CglResidualCapacity.hpp"
#include "CglStored.hpp"
#include "CglTwomir.hpp"
#include "CglZeroHalf.hpp"

#include "CbcBranchActual.hpp"
#include "CbcBranchCut.hpp"
#include "CbcBranchLotsize.hpp"
#include "CbcCompareActual.hpp"
#include "CbcCompareObjective.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcHeuristic.hpp"
#include "CbcHeuristicDiveCoefficient.hpp"
#include "CbcHeuristicDiveFractional.hpp"
#include "CbcHeuristicDiveGuided.hpp"
#include "CbcHeuristicDiveLineSearch.hpp"
#include "CbcHeuristicDivePseudoCost.hpp"
#include "CbcHeuristicDiveVectorLength.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcHeuristicPivotAndFix.hpp"
#include "CbcHeuristicRINS.hpp"
#include "CbcHeuristicRandRound.hpp"
#include "CbcMessage.hpp"
#include "CbcMipStartIO.hpp"
#include "CbcModel.hpp"
#include "CbcParam.hpp"
#include "CbcParamUtils.hpp"
#include "CbcSolver.hpp"
#include "CbcSolverAnalyze.hpp"
#include "CbcSolverExpandKnapsack.hpp"
#include "CbcSolverHeuristics.hpp"
#include "CbcStrategy.hpp"
#include "CbcTreeLocal.hpp"
//#define CBC_USE_OPENMP
#ifdef CBC_USE_OPENMP
#include "omp.h"
#endif

//#############################################################################
//#############################################################################

void printGeneralMessage(CbcModel &model, std::string message, int type)
{
   if (message.length()) {
      model.messageHandler()->message(type, model.messages())
         << message << CoinMessageEol;
   }
}

//#############################################################################
//#############################################################################

void printGeneralWarning(CbcModel &model, std::string message, int type)
{
   if (message.length()) {
      model.messageHandler()->message(type, model.messages())
         << message << CoinMessageEol;
   }
}

//###########################################################################
// Define symbols
//###########################################################################

#ifndef COIN_HAS_LINK
#define COIN_HAS_LINK
#endif
#ifdef COIN_HAS_LINK
#include "CbcLinked.hpp"
#endif

#ifndef CLP_OUTPUT_FORMAT
#define CLP_OUTPUT_FORMAT % 15.8g
#define CLP_INTEGER_OUTPUT_FORMAT % 15ld
#endif

#define CLP_QUOTE(s) CLP_STRING(s)
#define CLP_STRING(s) #s

#ifndef COINUTILS_HAS_GLPK
#define GLP_UNDEF 1
#define GLP_FEAS 2
#define GLP_INFEAS 3
#define GLP_NOFEAS 4
#define GLP_OPT 5
#endif

#ifndef CBC_QUIET
#define CBC_QUIET 0
#endif

#if defined(WSSMP_BARRIER) || defined(UFL_BARRIER) || defined(TAUCS_BARRIER)
#define FOREIGN_BARRIER
#endif

//#define IN_BRANCH_AND_BOUND (0x01000000|262144)
#define IN_BRANCH_AND_BOUND (0x01000000 | 262144 | 128 | 1024 | 2048)
//#define IN_BRANCH_AND_BOUND (0x01000000|262144|128)

//###########################################################################
// Global variables (can we get rid of these?)
//###########################################################################

static int initialPumpTune = -1;

#ifdef CGL_WRITEMPS
extern double *debugSolution;
extern int debugNumberColumns;
#endif

static CbcModel *currentBranchModel = NULL;

CglPreProcess *cbcPreProcessPointer = NULL;

static char cbcCrashAnnounced = 0;

#ifdef CLP_MALLOC_STATISTICS
static double malloc_times = 0.0;
static double malloc_total = 0.0;
static int malloc_amount[] = {0,    32,    128,   256,    1024,
                              4096, 16384, 65536, 262144, INT_MAX};
static int malloc_n = 10;
double malloc_counts[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
bool malloc_counts_on = true;
#endif

//###########################################################################
// Forward declarations (move to header?)
//###########################################################################

static void statistics(ClpSimplex *originalModel, ClpSimplex *model);
static bool maskMatches(const int *starts, char **masks, std::string &check);
static void generateCode(CbcModel *model, const char *fileName, int type,
                         int preProcess);
#ifdef CBC_HAS_NAUTY
// returns number of constraints added
static int nautiedConstraints(CbcModel &model, int maxPass);
#endif

// dummy fake main programs for UserClp and UserCbc
void fakeMain(ClpSimplex &model, OsiSolverInterface &osiSolver,
              CbcModel &babSolver);
void fakeMain2(ClpSimplex &model, OsiClpSolverInterface &osiSolver,
               int options);

int CbcClpUnitTest(const CbcModel &saveModel, const std::string &dirMiplib,
                   int testSwitch, const double *stuff, std::deque<std::string>,
                   int callBack(CbcModel *currentSolver, int whereFrom),
                   CbcParameters &parameters);

//###########################################################################
// Some debugging routines
//###########################################################################

//###########################################################################
//###########################################################################

#ifdef NEW_DEBUG_AND_FILL

void *operator new(size_t size) {
  void *p = malloc(size);
  char *xx = (char *)p;
  memset(xx, 0x20, size);
  return p;
}

//###########################################################################
//###########################################################################

void operator delete(void *p) throw() { free(p); }

#endif // end NEW_DEBUG_AND_FILL

//###########################################################################
//###########################################################################

#ifdef CLP_MALLOC_STATISTICS

#include "stolen_from_ekk_malloc.cpp"

void *operator new(size_t size) throw(std::bad_alloc) {
  malloc_times++;
  malloc_total += size;
  int i;
  for (i = 0; i < malloc_n; i++) {
    if ((int)size <= malloc_amount[i]) {
      malloc_counts[i]++;
      break;
    }
  }
#ifdef DEBUG_MALLOC
  void *p;
  if (malloc_counts_on)
    p = stolen_from_ekk_mallocBase(size);
  else
    p = malloc(size);
#else
  void *p = malloc(size);
#endif
  // char * xx = (char *) p;
  // memset(xx,0,size);
  // Initialize random seed
  // CoinSeedRandom(987654321);
  return p;
}

//###########################################################################
//###########################################################################

void operator delete(void *p) throw() {
#ifdef DEBUG_MALLOC
  if (malloc_counts_on)
    stolen_from_ekk_freeBase(p);
  else
    free(p);
#else
  free(p);
#endif
}

//###########################################################################
//###########################################################################

static void malloc_stats2() {
  double average = malloc_total / malloc_times;
  printf("count %g bytes %g - average %g\n", malloc_times, malloc_total,
         average);
  for (int i = 0; i < malloc_n; i++) {
    printf("%g ", malloc_counts[i]);
  }
  printf("\n");
  malloc_times = 0.0;
  malloc_total = 0.0;
  memset(malloc_counts, 0, sizeof(malloc_counts));
  // print results
}
#endif // CLP_MALLOC_STATISTICS

//###########################################################################
// Debug checks on special ordered sets.

// This is active only for debugging. The entire body of the routine becomes
// a noop when COIN_DEVELOP is not defined. To avoid compiler warnings, the
// formal parameters also need to go away.
//###########################################################################

#ifdef COIN_DEVELOP
void checkSOS(CbcModel *babModel, const OsiSolverInterface *solver)
#else
void checkSOS(CbcModel * /*babModel*/, const OsiSolverInterface * /*solver*/)
#endif
{
#ifdef COIN_DEVELOP
  if (!babModel->ownObjects())
    return;
#if COIN_DEVELOP > 2
  // const double *objective = solver->getObjCoefficients() ;
  const double *columnLower = solver->getColLower();
  const double *columnUpper = solver->getColUpper();
  const double *solution = solver->getColSolution();
  // int numberRows = solver->getNumRows();
  // double direction = solver->getObjSense();
  // int iRow,iColumn;
#endif

  // Row copy
  CoinPackedMatrix matrixByRow(*solver->getMatrixByRow());
  // const double * elementByRow = matrixByRow.getElements();
  // const int * column = matrixByRow.getIndices();
  // const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
  const int *rowLength = matrixByRow.getVectorLengths();

  // Column copy
  CoinPackedMatrix matrixByCol(*solver->getMatrixByCol());
  const double *element = matrixByCol.getElements();
  const int *row = matrixByCol.getIndices();
  const CoinBigIndex *columnStart = matrixByCol.getVectorStarts();
  const int *columnLength = matrixByCol.getVectorLengths();

  const double *rowLower = solver->getRowLower();
  const double *rowUpper = solver->getRowUpper();
  OsiObject **objects = babModel->objects();
  int numberObjects = babModel->numberObjects();
  int numberColumns = solver->getNumCols();
  for (int iObj = 0; iObj < numberObjects; iObj++) {
    CbcSOS *objSOS = dynamic_cast<CbcSOS *>(objects[iObj]);
    if (objSOS) {
      int n = objSOS->numberMembers();
      const int *which = objSOS->members();
#if COIN_DEVELOP > 2
      const double *weight = objSOS->weights();
#endif
      int type = objSOS->sosType();
      // convexity row?
      int iColumn;
      iColumn = which[0];
      int j;
      int convex = -1;
      for (j = columnStart[iColumn];
           j < columnStart[iColumn] + columnLength[iColumn]; j++) {
        int iRow = row[j];
        double value = element[j];
        if (rowLower[iRow] == 1.0 && rowUpper[iRow] == 1.0 && value == 1.0) {
          // possible
          if (rowLength[iRow] == n) {
            if (convex == -1)
              convex = iRow;
            else
              convex = -2;
          }
        }
      }
      printf("set %d of type %d has %d members - possible convexity row %d\n",
             iObj, type, n, convex);
      for (int i = 0; i < n; i++) {
        iColumn = which[i];
        // Column may have been added
        if (iColumn < numberColumns) {
          int convex2 = -1;
          for (j = columnStart[iColumn];
               j < columnStart[iColumn] + columnLength[iColumn]; j++) {
            int iRow = row[j];
            if (iRow == convex) {
              double value = element[j];
              if (value == 1.0) {
                convex2 = iRow;
              }
            }
          }
          if (convex2 < 0 && convex >= 0) {
            printf("odd convexity row\n");
            convex = -2;
          }
#if COIN_DEVELOP > 2
          printf("col %d has weight %g and value %g, bounds %g %g\n", iColumn,
                 weight[i], solution[iColumn], columnLower[iColumn],
                 columnUpper[iColumn]);
#endif
        }
      }
    }
  }
#endif // COIN_DEVELOP
}

//###########################################################################
// Signal handling
//###########################################################################

//###########################################################################
// Allow for interrupts
// But is this threadsafe? (so switched off by option)
//###########################################################################

extern "C" {
static void signal_handler(int whichSignal) {
  if (currentBranchModel != NULL) {
    currentBranchModel->sayEventHappened(); // say why stopped
    if (currentBranchModel->heuristicModel())
      currentBranchModel->heuristicModel()->sayEventHappened();
  }
  return;
}
}

//###########################################################################
// Handle crashes
//###########################################################################

#if defined(HAVE_SIGNAL_H) && defined(HAVE_EXECINFO_H)
void CbcCrashHandler(int sig) {
  char signame[256] = "";
  switch (sig) {
  case SIGILL:
    strcpy(signame, "SIGILL");
    break;
  case SIGSEGV:
    strcpy(signame, "SIGSEGV");
    break;
  case SIGABRT:
    strcpy(signame, "SIGABRT");
    break;
  }

  fflush(stderr);
  fflush(stdout);
  fprintf(
      stderr,
      "\n\nERROR while running Cbc. Signal %s caught. Getting stack trace.\n",
      signame);
  fflush(stderr);
  {
    char *st = getenv("RUNNING_TEST");
    if (st) {
      fprintf(stderr, "Error happened while running the \"%s\" test\n", st);
      fflush(stderr);
    }
  }

#define MAX_FRAMES 50
  void *array[MAX_FRAMES];
  size_t size;
  char **strings;
  size_t i;

  size = backtrace(array, MAX_FRAMES);
  strings = backtrace_symbols(array, size);

  for (i = 0; i < size; i++) {
    fprintf(stderr, "%s\n", strings[i]);
    fflush(stderr);
  }
  fprintf(stderr, "\n\n");
  fflush(stderr);

  free(strings);

  if (!cbcCrashAnnounced) {
    cbcCrashAnnounced = 1;
    abort();
  }
#undef MAX_FRAMES
}
#endif

//###########################################################################
// Some additional helper functions
//###########################################################################

static void putBackOtherSolutions(CbcModel *presolvedModel, CbcModel *model,
                                  CglPreProcess *preProcess) {
  int numberSolutions = presolvedModel->numberSavedSolutions();
  int numberColumns = presolvedModel->getNumCols();
  if (numberSolutions > 1) {
    model->deleteSolutions();
    double *bestSolution =
        CoinCopyOfArray(presolvedModel->bestSolution(), numberColumns);
    // double cutoff = presolvedModel->getCutoff();
    double objectiveValue = presolvedModel->getObjValue();
    // model->createSpaceForSavedSolutions(numberSolutions-1);
    for (int iSolution = numberSolutions - 1; iSolution >= 0; iSolution--) {
      presolvedModel->setCutoff(COIN_DBL_MAX);
      presolvedModel->solver()->setColSolution(
          presolvedModel->savedSolution(iSolution));
      // presolvedModel->savedSolutionObjective(iSolution));
      preProcess->postProcess(*presolvedModel->solver(), false);
      model->setBestSolution(preProcess->originalModel()->getColSolution(),
                             model->solver()->getNumCols(),
                             presolvedModel->savedSolutionObjective(iSolution));
    }
    presolvedModel->setBestObjectiveValue(objectiveValue);
    presolvedModel->solver()->setColSolution(bestSolution);
    // presolvedModel->setBestSolution(bestSolution,numberColumns,objectiveValue);
    delete [] bestSolution;
  }
}

//###########################################################################
//###########################################################################

// For when number of column is messed up e.g. BiLinear
static int numberPrintingColumns(const OsiSolverInterface *solver) {
#ifdef COIN_HAS_LINK
  const OsiSolverLink *linkSolver = dynamic_cast<const OsiSolverLink *>(solver);
  if (!linkSolver)
    return solver->getNumCols();
  return linkSolver->coinModel()->numberColumns();
#else
  return solver->getNumCols();
#endif
}

//###########################################################################
//###########################################################################

//###########################################################################
// CbcUser definitions
//###########################################################################

//###########################################################################
//###########################################################################

// User stuff (base class)
CbcUser::CbcUser() : coinModel_(NULL), userName_("null") {}

//###########################################################################
//###########################################################################

CbcUser::~CbcUser() { delete coinModel_; }

//###########################################################################
//###########################################################################

// Copy constructor
CbcUser::CbcUser(const CbcUser &rhs) {
  if (rhs.coinModel_)
    coinModel_ = new CoinModel(*rhs.coinModel_);
  else
    coinModel_ = NULL;
  userName_ = rhs.userName_;
}

//###########################################################################
//###########################################################################

// Assignment operator
CbcUser &CbcUser::operator=(const CbcUser &rhs) {
  if (this != &rhs) {
    if (rhs.coinModel_)
      coinModel_ = new CoinModel(*rhs.coinModel_);
    else
      coinModel_ = NULL;
    userName_ = rhs.userName_;
  }
  return *this;
}

//###########################################################################
//###########################################################################

//###########################################################################
// CbcStopNow definitions
//###########################################################################

//###########################################################################
//###########################################################################

CbcStopNow::CbcStopNow() {}

//###########################################################################
//###########################################################################

CbcStopNow::~CbcStopNow() {}

//###########################################################################
//###########################################################################

// Copy constructor
CbcStopNow::CbcStopNow(const CbcStopNow &) {}

//###########################################################################
//###########################################################################

// Assignment operator
CbcStopNow &CbcStopNow::operator=(const CbcStopNow &rhs) {
  if (this != &rhs) {
  }
  return *this;
}

//###########################################################################
//###########################################################################

// Clone
CbcStopNow *CbcStopNow::clone() const { return new CbcStopNow(*this); }

//###########################################################################
//###########################################################################

//###########################################################################
// CbcSsolver definitions
//###########################################################################

CbcSolver::CbcSolver()
    : babModel_(NULL), userFunction_(NULL), statusUserFunction_(NULL),
      originalSolver_(NULL), originalCoinModel_(NULL), cutGenerator_(NULL),
      numberUserFunctions_(0), numberCutGenerators_(0),
      startTime_(CoinCpuTime()), doMiplib_(false), noPrinting_(false),
      readMode_(1) {
  callBack_ = new CbcStopNow();
  //fillParameters();
}

//###########################################################################
//###########################################################################

CbcSolver::CbcSolver(const OsiClpSolverInterface &solver)
    : babModel_(NULL), userFunction_(NULL), statusUserFunction_(NULL),
      originalSolver_(NULL), originalCoinModel_(NULL), cutGenerator_(NULL),
      numberUserFunctions_(0), numberCutGenerators_(0),
      startTime_(CoinCpuTime()), doMiplib_(false), noPrinting_(false),
      readMode_(1) {
  callBack_ = new CbcStopNow();
  model_ = CbcModel(solver);
  //fillParameters();
}

//###########################################################################
//###########################################################################

CbcSolver::CbcSolver(const CbcModel &solver)
    : babModel_(NULL), userFunction_(NULL), statusUserFunction_(NULL),
      originalSolver_(NULL), originalCoinModel_(NULL), cutGenerator_(NULL),
      numberUserFunctions_(0), numberCutGenerators_(0),
      startTime_(CoinCpuTime()), doMiplib_(false), noPrinting_(false),
      readMode_(1) {
  callBack_ = new CbcStopNow();
  model_ = solver;
  //fillParameters();
}

//###########################################################################
//###########################################################################

CbcSolver::~CbcSolver() {
  int i;
  for (i = 0; i < numberUserFunctions_; i++)
    delete userFunction_[i];
  delete[] userFunction_;
  for (i = 0; i < numberCutGenerators_; i++)
    delete cutGenerator_[i];
  delete[] cutGenerator_;
  delete[] statusUserFunction_;
  delete originalSolver_;
  delete originalCoinModel_;
  delete babModel_;
  delete callBack_;
}

//###########################################################################
//###########################################################################

// Copy constructor
CbcSolver::CbcSolver(const CbcSolver &rhs)
    : model_(rhs.model_), babModel_(NULL), userFunction_(NULL),
      statusUserFunction_(NULL),
      cutGenerator_(new CglCutGenerator *[rhs.numberCutGenerators()]),
      numberUserFunctions_(rhs.numberUserFunctions_),
      numberCutGenerators_(rhs.numberCutGenerators()),
      startTime_(CoinCpuTime()), doMiplib_(rhs.doMiplib_),
      noPrinting_(rhs.noPrinting_), readMode_(rhs.readMode_) {
   //fillParameters();
  if (rhs.babModel_)
    babModel_ = new CbcModel(*rhs.babModel_);
  userFunction_ = new CbcUser *[numberUserFunctions_];
  int i;
  for (i = 0; i < numberUserFunctions_; i++)
    userFunction_[i] = rhs.userFunction_[i]->clone();
  this->parameters_ = rhs.parameters_;
  this->clpParameters_ = rhs.clpParameters_;
  for (i = 0; i < numberCutGenerators_; i++)
    cutGenerator_[i] = rhs.cutGenerator_[i]->clone();
  callBack_ = rhs.callBack_->clone();
  originalSolver_ = NULL;
  if (rhs.originalSolver_) {
    OsiSolverInterface *temp = rhs.originalSolver_->clone();
    originalSolver_ = dynamic_cast<OsiClpSolverInterface *>(temp);
    assert(originalSolver_);
  }
  originalCoinModel_ = NULL;
  if (rhs.originalCoinModel_)
    originalCoinModel_ = new CoinModel(*rhs.originalCoinModel_);
}

//###########################################################################
//###########################################################################

// Assignment operator
CbcSolver &CbcSolver::operator=(const CbcSolver &rhs) {
  if (this != &rhs) {
    int i;
    for (i = 0; i < numberUserFunctions_; i++)
      delete userFunction_[i];
    delete[] userFunction_;
    for (i = 0; i < numberCutGenerators_; i++)
      delete cutGenerator_[i];
    delete[] statusUserFunction_;
    delete originalSolver_;
    delete originalCoinModel_;
    statusUserFunction_ = NULL;
    delete babModel_;
    delete callBack_;
    numberUserFunctions_ = rhs.numberUserFunctions_;
    startTime_ = rhs.startTime_;
    this->parameters_ = rhs.parameters_;
    this->clpParameters_ = rhs.clpParameters_;
    for (i = 0; i < numberCutGenerators_; i++)
      cutGenerator_[i] = rhs.cutGenerator_[i]->clone();
    noPrinting_ = rhs.noPrinting_;
    readMode_ = rhs.readMode_;
    doMiplib_ = rhs.doMiplib_;
    model_ = rhs.model_;
    if (rhs.babModel_)
      babModel_ = new CbcModel(*rhs.babModel_);
    else
      babModel_ = NULL;
    userFunction_ = new CbcUser *[numberUserFunctions_];
    for (i = 0; i < numberUserFunctions_; i++)
      userFunction_[i] = rhs.userFunction_[i]->clone();
    callBack_ = rhs.callBack_->clone();
    originalSolver_ = NULL;
    if (rhs.originalSolver_) {
      OsiSolverInterface *temp = rhs.originalSolver_->clone();
      originalSolver_ = dynamic_cast<OsiClpSolverInterface *>(temp);
      assert(originalSolver_);
    }
    originalCoinModel_ = NULL;
    if (rhs.originalCoinModel_)
      originalCoinModel_ = new CoinModel(*rhs.originalCoinModel_);
  }
  return *this;
}

//###########################################################################
//###########################################################################

// Get int value
int CbcSolver::intValue(int code) {
  return parameters_[code]->intVal();
}

//###########################################################################
//###########################################################################

// Set int value
void CbcSolver::setIntValue(int code, int value) {
  parameters_[code]->setVal(value);
}

//###########################################################################
//###########################################################################

// Get double value
double CbcSolver::doubleValue(int code) {
  return parameters_[code]->dblVal();
}

//###########################################################################
//###########################################################################

// Set double value
void CbcSolver::setDoubleValue(int code, double value) {
  parameters_[code]->setVal(value);
}

//###########################################################################
//###########################################################################

// User function (NULL if no match)
CbcUser *CbcSolver::userFunction(const char *name) const {
  int i;
  for (i = 0; i < numberUserFunctions_; i++) {
    if (!strcmp(name, userFunction_[i]->name().c_str()))
      break;
  }
  if (i < numberUserFunctions_)
    return userFunction_[i];
  else
    return NULL;
}

//###########################################################################
//###########################################################################

// Add user function
void CbcSolver::addUserFunction(CbcUser *function) {
  CbcUser **temp = new CbcUser *[numberUserFunctions_ + 1];
  int i;
  for (i = 0; i < numberUserFunctions_; i++)
    temp[i] = userFunction_[i];
  delete[] userFunction_;
  userFunction_ = temp;
  userFunction_[numberUserFunctions_++] = function->clone();
  delete[] statusUserFunction_;
  statusUserFunction_ = NULL;
}

//###########################################################################
//###########################################################################

// Set user call back
void CbcSolver::setUserCallBack(CbcStopNow *function) {
  delete callBack_;
  callBack_ = function->clone();
}

//###########################################################################
//###########################################################################

// Copy of model on initial load (will contain output solutions)
void CbcSolver::setOriginalSolver(OsiClpSolverInterface *originalSolver) {
  delete originalSolver_;
  OsiSolverInterface *temp = originalSolver->clone();
  originalSolver_ = dynamic_cast<OsiClpSolverInterface *>(temp);
  assert(originalSolver_);
}

//###########################################################################
//###########################################################################

// Copy of model on initial load
void CbcSolver::setOriginalCoinModel(CoinModel *originalCoinModel) {
  delete originalCoinModel_;
  originalCoinModel_ = new CoinModel(*originalCoinModel);
}

//###########################################################################
//###########################################################################

// Add cut generator
void CbcSolver::addCutGenerator(CglCutGenerator *generator) {
  CglCutGenerator **temp = new CglCutGenerator *[numberCutGenerators_ + 1];
  int i;
  for (i = 0; i < numberCutGenerators_; i++)
    temp[i] = cutGenerator_[i];
  delete[] cutGenerator_;
  cutGenerator_ = temp;
  cutGenerator_[numberCutGenerators_++] = generator->clone();
}

//###########################################################################
//###########################################################################

//TODO Should we replace with the STL version
static bool ends_with(std::string const &value, std::string const &ending) {
  if (ending.size() > value.size())
    return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

//###########################################################################
// Wrappers for CbcMain0, CbcMain1. The various forms of callCbc will
// eventually resolve to a call to CbcMain0 followed by a call to callCbc1.
//###########################################################################

//###########################################################################
// This is the function for setting things up, default parameters, etc.
//###########################################################################

void CbcMain0(CbcModel &model, CbcParameters &parameters) {
#if defined(HAVE_SIGNAL_H) && defined(HAVE_EXECINFO_H)
  signal(SIGSEGV, CbcCrashHandler);
  signal(SIGABRT, CbcCrashHandler);
  signal(SIGFPE, CbcCrashHandler);

#endif

  ClpParameters &clpParameters = parameters.clpParameters();

#ifndef CBC_OTHER_SOLVER
  OsiClpSolverInterface *originalSolver =
      dynamic_cast<OsiClpSolverInterface *>(model.solver());
#elif CBC_OTHER_SOLVER == 1
  OsiCpxSolverInterface *originalSolver =
      dynamic_cast<OsiCpxSolverInterface *>(model.solver());
  // Dummy solvers
  OsiClpSolverInterface dummySolver;
  ClpSimplex *lpSolver = dummySolver.getModelPtr();
  OsiCpxSolverInterface *clpSolver = originalSolver;
#endif
  assert(originalSolver);
  CoinMessageHandler *generalMessageHandler = originalSolver->messageHandler();
  generalMessageHandler->setPrefix(true);
#ifndef CBC_OTHER_SOLVER
  OsiSolverInterface *solver = model.solver();
  OsiClpSolverInterface *clpSolver =
      dynamic_cast<OsiClpSolverInterface *>(solver);
  ClpSimplex *lpSolver = clpSolver->getModelPtr();
  lpSolver->setPerturbation(50);
  lpSolver->messageHandler()->setPrefix(false);
  clpParameters.setModel(lpSolver);
#endif
  // establishParams(numberParameters, parameters) ;
  int doIdiot = -1;
  int outputFormat = 2;
  int substitution = 3;
  int dualize = 3;
  int preSolve = 5;
  int doSprint = -1;
  int testOsiParameters = -1;
  clpParameters[ClpParam::DUALBOUND]->setVal(lpSolver->dualBound());
  clpParameters[ClpParam::DUALTOLERANCE]->setVal(lpSolver->dualTolerance());
  clpParameters[ClpParam::IDIOT]->setVal(doIdiot);
  clpParameters[ClpParam::PRESOLVETOLERANCE]->setVal(1.0e-8);
  clpParameters[ClpParam::MAXFACTOR]->setVal(lpSolver->factorizationFrequency());
  clpParameters[ClpParam::MAXITERATION]->setVal(lpSolver->maximumIterations());
  clpParameters[ClpParam::PRESOLVEPASS]->setVal(preSolve);
  clpParameters[ClpParam::PERTVALUE]->setVal(lpSolver->perturbation());
  clpParameters[ClpParam::PRIMALTOLERANCE]->setVal(lpSolver->primalTolerance());
  clpParameters[ClpParam::PRIMALWEIGHT]->setVal(lpSolver->infeasibilityCost());
  clpParameters[ClpParam::SPRINT]->setVal(doSprint);
  clpParameters[ClpParam::SUBSTITUTION]->setVal(substitution);
  clpParameters[ClpParam::DUALIZE]->setVal(dualize);
  parameters[CbcParam::OUTPUTFORMAT]->setVal(outputFormat);
  parameters[CbcParam::LOGLEVEL]->setVal(1);
  parameters[CbcParam::LPLOGLEVEL]->setVal(1);
  clpParameters[ClpParam::LOGLEVEL]->setType(CoinParam::paramInt); // as can be set
  clpSolver->messageHandler()->setLogLevel(1);
  lpSolver->setLogLevel(1);
  parameters[CbcParam::TIMELIMIT]->setVal(1.0e8);
  parameters[CbcParam::TESTOSI]->setVal(testOsiParameters);
  parameters[CbcParam::FPUMPTUNE]->setVal(1003);
  initialPumpTune = 1003;
#ifdef CBC_THREAD
  parameters[CbcParam::THREADS]->setVal(0);
#endif
  // Set up likely cut generators and defaults
  parameters[CbcParam::CLIQUECUTS]->setVal("ifmove");
  parameters[CbcParam::ODDWHEELCUTS]->setVal("ifmove");
  parameters[CbcParam::CLQSTRENGTHENING]->setVal("after");
  parameters[CbcParam::USECGRAPH]->setVal("on");
  parameters[CbcParam::AGGREGATEMIXED]->setVal(1);
  parameters[CbcParam::BKPIVOTINGSTRATEGY]->setVal(3);
  parameters[CbcParam::BKMAXCALLS]->setVal(1000);
  parameters[CbcParam::BKCLQEXTMETHOD]->setVal(4);
  parameters[CbcParam::ODDWEXTMETHOD]->setVal(2);
  parameters[CbcParam::PREPROCESS]->setVal("sos");
  parameters[CbcParam::MIPOPTIONS]->setVal(1057);
  parameters[CbcParam::CUTPASSINTREE]->setVal(10);
  parameters[CbcParam::MOREMIPOPTIONS]->setVal(-1);
  parameters[CbcParam::MAXHOTITS]->setVal(100);
  parameters[CbcParam::CUTSTRATEGY]->setVal("on");
  parameters[CbcParam::HEURISTICSTRATEGY]->setVal("on");
  parameters[CbcParam::NODESTRATEGY]->setVal("fewest");
  parameters[CbcParam::GOMORYCUTS]->setVal("ifmove");
  parameters[CbcParam::PROBINGCUTS]->setVal("ifmove");
  parameters[CbcParam::KNAPSACKCUTS]->setVal("ifmove");
  parameters[CbcParam::ZEROHALFCUTS]->setVal("ifmove");
  parameters[CbcParam::REDSPLITCUTS]->setVal("off");
  parameters[CbcParam::REDSPLIT2CUTS]->setVal("off");
  parameters[CbcParam::GMICUTS]->setVal("off");
  parameters[CbcParam::MIRCUTS]->setVal("ifmove");
  parameters[CbcParam::FLOWCUTS]->setVal("ifmove");
  parameters[CbcParam::TWOMIRCUTS]->setVal("ifmove");
  parameters[CbcParam::LANDPCUTS]->setVal("off");
  parameters[CbcParam::RESIDCAPCUTS]->setVal("off");
  parameters[CbcParam::ROUNDING]->setVal("on");
  parameters[CbcParam::FPUMP]->setVal("on");
  parameters[CbcParam::GREEDY]->setVal("on");
  parameters[CbcParam::COMBINE]->setVal("off");
  parameters[CbcParam::CROSSOVER]->setVal("off");
#ifdef CBC_HAS_NAUTY
#ifndef CBC_LIGHTWEIGHT_NAUTY
  parameters[CbcParam::ORBITAL]->setVal("on");
#else
  parameters[CbcParam::ORBITAL]->setVal("lightweight");
#endif
#endif
  //parameters[CbcParam::PIVOTANDCOMPLEMENT]->setVal("off");
  parameters[CbcParam::PIVOTANDFIX]->setVal("off");
  parameters[CbcParam::RANDROUND]->setVal("off");
  parameters[CbcParam::NAIVE]->setVal("off");
  parameters[CbcParam::RINS]->setVal("off");
  parameters[CbcParam::DINS]->setVal("off");
  parameters[CbcParam::RENS]->setVal("off");
  parameters[CbcParam::LOCALTREE]->setVal("off");
  parameters[CbcParam::BRANCHPRIORITY]->setVal("off");

  model.messageHandler()->setLogLevel(1);
  model.setNumberBeforeTrust(10);
  parameters[CbcParam::NUMBERBEFORE]->setVal(10);
  parameters[CbcParam::MAXNODES]->setVal(model.getMaximumNodes());
  model.setNumberStrong(5);
  parameters[CbcParam::STRONGBRANCHING]->setVal(model.numberStrong());
  parameters[CbcParam::INFEASIBILITYWEIGHT]->setVal(model.getDblParam(CbcModel::CbcInfeasibilityWeight));
  parameters[CbcParam::INTEGERTOLERANCE]->setVal(model.getDblParam(CbcModel::CbcIntegerTolerance));
  parameters[CbcParam::INCREMENT]->setVal(model.getDblParam(CbcModel::CbcCutoffIncrement));
}

//###########################################################################
//###########################################################################

//###########################################################################
// CbcMain 1
// Meaning of whereFrom:
//   1 after initial solve by dualsimplex etc
//   2 after preprocessing
//   3 just before branchAndBound (so user can override)
//   4 just after branchAndBound (before postprocessing)
//   5 after postprocessing
//   6 after a user called heuristic phase
//###########################################################################

int CbcMain1(std::deque<std::string> inputQueue, CbcModel &model,
             CbcParameters &parameters,
             int callBack(CbcModel *currentSolver, int whereFrom),
             ampl_info *info) {

  ClpParameters &clpParameters = parameters.clpParameters();

  std::ostringstream buffer;
  std::string field, message, fileName;
  FILE *fp;
#ifdef DEBUG_CBC_PYTHON
  // Probably being stupid but finding it difficult to debug from python
  static int startedPython=0;
  if (!startedPython) {
    printf("debug python:");
    getchar();
    startedPython=1;
  }
#endif

  double totalTime = parameters.getTotalTime();
  bool useSignalHandler = parameters.useSignalHandler();
  CbcModel &model_ = model;
  CglPreProcess *preProcessPointer = NULL;
  OsiSolverInterface *saveSolver = NULL;
  CglPreProcess process;
  // Save a copy of input for unit testing
  // Could also be used for friendly error messages?
  std::deque<std::string> saveInputQueue = inputQueue;
  // copy if we are reading from option file
  std::deque<std::string> partInputQueue;
  // Meaning 0 - start at very beginning
  // 1 start at beginning of preprocessing
  // 2 start at beginning of branch and bound
#ifndef CBC_USE_INITIAL_TIME
#define CBC_USE_INITIAL_TIME 1
#endif
#if CBC_USE_INITIAL_TIME == 0
  if (model_.useElapsedTime())
    model_.setDblParam(CbcModel::CbcStartSeconds, CoinGetTimeOfDay());
  else
    model_.setDblParam(CbcModel::CbcStartSeconds, CoinCpuTime());
#endif
  CbcModel *babModel_ = NULL;
#ifdef COINUTILS_HAS_GLPK
  glp_tran *coin_glp_tran = NULL;
  glp_prob *coin_glp_prob = NULL;
#endif
  int returnMode = 1;
  int statusUserFunction_[1];
  int numberUserFunctions_ = 1; // to allow for ampl
  // Statistics
  double statistics_seconds = 0.0, statistics_obj = 0.0;
  double statistics_sys_seconds = 0.0, statistics_elapsed_seconds = 0.0;
  CoinWallclockTime();
  double statistics_continuous = 0.0, statistics_tighter = 0.0;
  double statistics_cut_time = 0.0;
  int statistics_nodes = 0, statistics_iterations = 0;
  int statistics_nrows = 0, statistics_ncols = 0;
  int statistics_nprocessedrows = 0, statistics_nprocessedcols = 0;
  std::string statistics_result;
  int *statistics_number_cuts = NULL;
  const char **statistics_name_generators = NULL;
  int statistics_number_generators = 0;
  int currentBestSolution = 0;
  memset(statusUserFunction_, 0, numberUserFunctions_ * sizeof(int));
  /* Note
       This is meant as a stand-alone executable to do as much of coin as
     possible. It should only have one solver known to it.
    */
  CoinMessageHandler *generalMessageHandler = model_.messageHandler();
  generalMessageHandler->setPrefix(false);
  int numberLotSizing = 0;
  typedef struct {
    double low;
    double high;
    int column;
  } lotStruct;
  lotStruct *lotsize = NULL;
  typedef struct {lotStruct * lotsize;int numberLotSizing;} lotStruct2;
  lotStruct2 passSCtopreprocess;
#ifndef CBC_OTHER_SOLVER
  OsiClpSolverInterface *originalSolver =
      dynamic_cast<OsiClpSolverInterface *>(model_.solver());
  assert(originalSolver);
  // Move handler across if not default
  if (!originalSolver->defaultHandler() &&
      originalSolver->getModelPtr()->defaultHandler())
    originalSolver->getModelPtr()->passInMessageHandler(
        originalSolver->messageHandler());
  CoinMessages generalMessages = originalSolver->getModelPtr()->messages();
  if (originalSolver->getModelPtr()->logLevel() == 0) {
    generalMessageHandler->setLogLevel(0);
    parameters.disableWelcomePrinting();
  }
#elif CBC_OTHER_SOLVER == 1
  /*
    The only other solver that's ever been used is cplex, and the use is
    limited -- do the root with clp and all the cbc smarts, then give the
    problem over to cplex to finish. Although the defines can be read in some
    places to allow other options, nothing's been tested and success is
    unlikely.
  */
  OsiCpxSolverInterface *originalSolver =
      dynamic_cast<OsiCpxSolverInterface *>(model_.solver());
  assert(originalSolver);
  OsiClpSolverInterface dummySolver;
  OsiCpxSolverInterface *clpSolver = originalSolver;
  CoinMessages generalMessages = dummySolver.getModelPtr()->messages();
  parameters.enablePrinting();
#endif
  // Say not in integer
  int integerStatus = -1;
  // Say no resolve after cuts
  model_.setResolveAfterTakeOffCuts(false);
  double time0;
  double time0Elapsed = CoinGetTimeOfDay();
  {
    double time1 = CoinCpuTime(), time2;
    time0 = time1;
    double time1Elapsed = time0Elapsed;
    bool goodModel = (originalSolver->getNumCols()) ? true : false;

    // register signal handler
    // CoinSighandler_t saveSignal=signal(SIGINT,signal_handler);
#if CBC_QUIET < 2
    if (useSignalHandler)
      signal(SIGINT, signal_handler);
#endif
    // Set up all non-standard stuff
    int cutPass = -1234567;
    int cutPassInTree = -1234567;
    int tunePreProcess = 0;
    int testOsiParameters = -1;
    // 0 normal, 1 from ampl or MIQP etc (2 allows cuts)
    int complicatedInteger = 0;
    OsiSolverInterface *solver = model_.solver();
#ifndef CBC_OTHER_SOLVER
    OsiClpSolverInterface *clpSolver =
        dynamic_cast<OsiClpSolverInterface *>(solver);
    ClpSimplex *lpSolver = clpSolver->getModelPtr();
    if (parameters.noPrinting()) {
      lpSolver->setLogLevel(0);
    }
#else
    ClpSimplex *lpSolver = NULL;
#endif
    // For priorities etc
    int *priorities = NULL;
    int *branchDirection = NULL;
    double *pseudoDown = NULL;
    double *pseudoUp = NULL;
    double *solutionIn = NULL;
    int *prioritiesIn = NULL;
    std::vector<std::pair<std::string, double> > mipStart;
    std::vector<std::pair<std::string, double> > mipStartBefore;
    std::string mipStartFile = "";
    int numberSOS = 0;
    int *sosStart = NULL;
    int *sosIndices = NULL;
    char *sosType = NULL;
    double *sosReference = NULL;
    int *cut = NULL;
    int *sosPriority = NULL;
    CglStored storedAmpl;
    CoinModel *coinModel = NULL;
    CoinModel saveCoinModel;
    CoinModel saveTightenedModel;
    int *whichColumn = NULL;
    int *knapsackStart = NULL;
    int *knapsackRow = NULL;
    int numberKnapsack = 0;
    // default action on import
    int allowImportErrors = 0;
    int keepImportNames = 1;
    int doIdiot = -1;
    int outputFormat = 2;
    int slpValue = -1;
    int cppValue = -1;
    int printOptions = 0;
    int printMode = 0;
    int presolveOptions = 0;
    int substitution = 3;
    int dualize = 3;
    int doCrash = 0;
    int doVector = 0;
    int doSprint = -1;
    int doScaling = 4;

    if (info) {
       parameters[CbcParam::LOGLEVEL]->setVal(info->logLevel);
       goodModel = true;
       // FIXME Do we really need this? Seems to just mean that AMPL is being used
       statusUserFunction_[0] = 1;
       if (info->numberSos) {
          numberSOS = info->numberSos;
          sosStart = info->sosStart;
          sosIndices = info->sosIndices;
          sosType = info->sosType;
          sosReference = info->sosReference;
          sosPriority = info->sosPriority;
       }
       if (info->cut) {
          int numberRows = info->numberRows;
          int *whichRow = new int[numberRows];
          // Row copy
          const CoinPackedMatrix *matrixByRow = solver->getMatrixByRow();
          const double *elementByRow = matrixByRow->getElements();
          const int *column = matrixByRow->getIndices();
          const CoinBigIndex *rowStart = matrixByRow->getVectorStarts();
          const int *rowLength = matrixByRow->getVectorLengths();
          
          const double *rowLower = solver->getRowLower();
          const double *rowUpper = solver->getRowUpper();
          int nDelete = 0;
          for (int iRow = 0; iRow < numberRows; iRow++) {
             if (info->cut[iRow]) {
                whichRow[nDelete++] = iRow;
                int start = rowStart[iRow];
                storedAmpl.addCut(rowLower[iRow], rowUpper[iRow],
                                  rowLength[iRow], column + start,
                                  elementByRow + start);
             }
          }
          solver->deleteRows(nDelete, whichRow);
          delete[] whichRow;
       }
             solver->loadProblem(
                          info->numberColumns, info->numberRows,
                          reinterpret_cast<const CoinBigIndex *>(info->starts), info->rows,
                          info->elements, info->columnLower, info->columnUpper, info->objective,
                          info->rowLower, info->rowUpper);
      // take off cuts if ampl wants that
      if (info->cut && 0) {
         printf("AMPL CUTS OFF until global cuts fixed\n");
         info->cut = NULL;
      }

      //FIXME A big block was deleted which seemed as though it should never be entered.
      //      Check that.
      
      // If we had a solution use it
      if (info->primalSolution) {
         solver->setColSolution(info->primalSolution);
      }
      // status
      if (info->rowStatus) {
         unsigned char *statusArray = lpSolver->statusArray();
         int i;
         for (i = 0; i < info->numberColumns; i++)
            statusArray[i] = static_cast<unsigned char>(info->columnStatus[i]);
         statusArray += info->numberColumns;
         for (i = 0; i < info->numberRows; i++)
            statusArray[i] = static_cast<unsigned char>(info->rowStatus[i]);
         CoinWarmStartBasis *basis = lpSolver->getBasis();
         solver->setWarmStart(basis);
         delete basis;
      }
      freeArrays1(info);
      // modify objective if necessary
      solver->setObjSense(info->direction);
      solver->setDblParam(OsiObjOffset, -info->offset);
      if (info->offset) {
         buffer.str("");
         buffer << "Ampl objective offset is " << info->offset;
         printGeneralMessage(model, buffer.str());
      }
      // Set integer variables (unless nonlinear when set)
      if (!info->nonLinear) {
         for (int i = info->numberColumns - info->numberIntegers;
              i < info->numberColumns; i++)
            solver->setInteger(i);
      }
    }
    
    // set reasonable defaults
    int preSolve = 5;
    int preProcess = 4;
    bool useStrategy = false;
    bool preSolveFile = false;
    bool strongChanged = false;
    bool pumpChanged = false;

    double djFix = 1.0e100;
    double tightenFactor = 0.0;
    const char dirsep = CoinFindDirSeparator();
    double *debugValues = NULL;
    int numberDebugValues = -1;
    int basisHasValues = 0;
#ifndef CBC_OTHER_SOLVER
    double normalIncrement = model_.getCutoffIncrement();
    ;
#endif
    if (testOsiParameters >= 0) {
      // trying nonlinear - switch off some stuff
      preProcess = 0;
    }
    // Set up likely cut generators and defaults
    int nodeStrategy = 0;
    bool dominatedCuts = false;
    int doSOS = 1;
    int verbose = 0;
    CglGomory gomoryGen;
    // try larger limit
    gomoryGen.setLimitAtRoot(1000);
    gomoryGen.setLimit(50);
    // set default action (0=off,1=on,2=root,3=ifmove)
    int gomoryMode = CbcParameters::CGIfMove;
    assert (parameters[CbcParam::GOMORYCUTS]->modeVal()==gomoryMode);

    CglProbing probingGen;
    probingGen.setUsingObjective(1);
    probingGen.setMaxPass(1);
    probingGen.setMaxPassRoot(1);
    // Number of unsatisfied variables to look at
    probingGen.setMaxProbe(10);
    probingGen.setMaxProbeRoot(50);
    // How far to follow the consequences
    probingGen.setMaxLook(10);
    probingGen.setMaxLookRoot(50);
    probingGen.setMaxLookRoot(10);
    // Only look at rows with fewer than this number of elements
    probingGen.setMaxElements(200);
    probingGen.setMaxElementsRoot(300);
    probingGen.setRowCuts(3);
    // set default action (0=off,1=on,2=root,3=ifmove)
    int probingMode = CbcParameters::CGIfMove;
    assert (parameters[CbcParam::PROBINGCUTS]->modeVal()==probingMode);

    CglKnapsackCover knapsackGen;
    // knapsackGen.switchOnExpensive();
    // knapsackGen.setMaxInKnapsack(100);
    // set default action (0=off,1=on,2=root,3=ifmove)
    int knapsackMode = CbcParameters::CGIfMove;
    assert (parameters[CbcParam::KNAPSACKCUTS]->modeVal()==knapsackMode);

    CglRedSplit redsplitGen;
    // set default action (0=off,1=on,2=root,3=ifmove)
    // Off as seems to give some bad cuts
    int redsplitMode = CbcParameters::CGOff;
    assert (parameters[CbcParam::REDSPLITCUTS]->modeVal()==redsplitMode);

    CglRedSplit2 redsplit2Gen;
    // redsplit2Gen.setLimit(100);
    // set default action (0=off,1=on,2=root,3=ifmove)
    // Off
    int redsplit2Mode = CbcParameters::CGOff;
    assert (parameters[CbcParam::REDSPLIT2CUTS]->modeVal()==redsplit2Mode);

    CglGMI GMIGen; 
    int GMIMode = parameters[CbcParam::GMICUTS]->modeVal();

    std::string cgraphMode = "on";
    std::string clqstrMode = "after";
    CglBKClique bkCliqueGen;
    CoinBronKerbosch::PivotingStrategy bkPivotingStrategy = CoinBronKerbosch::PivotingStrategy::Weight;
    int maxCallsBK = 1000, bkClqExtMethod = 4;
    int cliqueMode = parameters[CbcParam::CLIQUECUTS]->modeVal();
    int oldCliqueMode = cliqueMode;
    CglOddWheel oddWheelGen;
    int oddWheelMode = CbcParameters::CGIfMove, oddWExtMethod = 2;
    assert (parameters[CbcParam::ODDWHEELCUTS]->modeVal()==oddWheelMode);

    // maxaggr,multiply,criterion(1-3)
    CglMixedIntegerRounding2 mixedGen(1, true, 1);
    // set default action (0=off,1=on,2=root,3=ifmove)
    int mixedMode = CbcParameters::CGIfMove;
    int mixedRoundStrategy = 1; 
    assert (parameters[CbcParam::MIRCUTS]->modeVal()==mixedMode);
    mixedGen.setDoPreproc(1); // safer (and better)

    CglFlowCover flowGen;
    // set default action (0=off,1=on,2=root,3=ifmove)
    int flowMode = CbcParameters::CGIfMove;
    assert (parameters[CbcParam::FLOWCUTS]->modeVal()==flowMode);

    CglTwomir twomirGen;
    twomirGen.setMaxElements(250);
    // set default action (0=off,1=on,2=root,3=ifmove)
    int twomirMode = CbcParameters::CGIfMove;
    assert (parameters[CbcParam::TWOMIRCUTS]->modeVal()==twomirMode);
#ifndef DEBUG_MALLOC
    CglLandP landpGen;
    landpGen.parameter().maximumCutLength = 2000;
    landpGen.validator().setMinViolation(1.0e-4);
#endif
    // set default action (0=off,1=on,2=root,3=ifmove)
    int landpMode = CbcParameters::CGOff;
    assert (parameters[CbcParam::LANDPCUTS]->modeVal()==landpMode);
    CglResidualCapacity residualCapacityGen;
    residualCapacityGen.setDoPreproc(1); // always preprocess
    // set default action (0=off,1=on,2=root,3=ifmove)
    int residualCapacityMode = CbcParameters::CGOff;
    assert (parameters[CbcParam::RESIDCAPCUTS]->modeVal()==residualCapacityMode);

    CglZeroHalf zerohalfGen;
    // zerohalfGen.switchOnExpensive();
    // set default action (0=off,1=on,2=root,3=ifmove)
    int zerohalfMode = CbcParameters::CGIfMove;
    assert (parameters[CbcParam::ZEROHALFCUTS]->modeVal()==zerohalfMode);

    // Stored cuts
    // bool storedCuts = false;

    int useCosts = 0;
    // don't use input solution
    int useSolution = -1;

    // total number of commands read
    int numberGoodCommands = 0;
    // Set false if user does anything advanced
    bool defaultSettings = true;

    // Hidden stuff for barrier
    int choleskyType = 0;
    int gamma = 0;
    int scaleBarrier = 0;
    int doKKT = 0;
    int crossover = 2; // do crossover unless quadratic
    bool biLinearProblem = false;
    // For names
    int lengthName = 0;
    std::vector<std::string> rowNames;
    std::vector<std::string> columnNames;
    // Default strategy stuff
    {
      // try changing tolerance at root
#define MORE_CUTS
#ifdef MORE_CUTS
      gomoryGen.setAwayAtRoot(0.005);
      twomirGen.setAwayAtRoot(0.005);
      twomirGen.setAway(0.01);
      // twomirGen.setMirScale(1,1);
      // twomirGen.setTwomirScale(1,1);
      // twomirGen.setAMax(2);
#else
      gomoryGen.setAwayAtRoot(0.01);
      twomirGen.setAwayAtRoot(0.01);
      twomirGen.setAway(0.01);
#endif
      // TODO: See what's going on with these free-standing parameters
      initialPumpTune = 1005043;
      tunePreProcess = 6;
      probingMode = CbcParameters::CGIfMove;
      parameters[CbcParam::DIVEOPT]->setVal(2);
      parameters[CbcParam::FPUMPITS]->setVal(30);
      parameters[CbcParam::FPUMPTUNE]->setVal(1005043);
      parameters[CbcParam::PROCESSTUNE]->setVal(7);
      parameters[CbcParam::DIVINGC]->setVal("on");
      parameters[CbcParam::RINS]->setVal("on");
      parameters[CbcParam::PROBINGCUTS]->setVal("ifmove");
      // parameters[iParam]->setVal("forceOnStrong");
      // probingMode = 8;
    }

    bool interactiveMode = false, canOpen;

    int status, iValue;
    double dValue;
    std::string prompt = "Cbc: ";

#if CBC_QUIET == 0
    if (parameters.printWelcome()) {
      buffer.str("");
      buffer << "Welcome to the CBC MILP Solver" << std::endl;
      if (strcmp(CBC_VERSION, "devel")) {
         buffer << "Version: " << CBC_VERSION  << std::endl;
      } else {
         buffer << "Version: Devel (unstable)" << std::endl;
      }
      buffer << "Build Date: " << __DATE__  << std::endl;
      printGeneralMessage(model_, buffer.str());
      // Print command line
      if (!inputQueue.empty()) {
        bool foundStrategy = false;
        buffer.str("");
        buffer << "command line -";
        for (int i = 0; i < inputQueue.size(); i++) {
           size_t found = inputQueue[i].find("strat"); 
           if (found != std::string::npos){
              foundStrategy = true;
           }
           buffer << " " << inputQueue[i];
        }
        if (!foundStrategy){
           buffer << " (default strategy 1)";
        }
        printGeneralMessage(model_, buffer.str());
      }
    }
#endif

    // If no arguments, just go into interactive mode
    if (!inputQueue.size()) {
      // see if AMPL
      if (!info) {
	interactiveMode = true;
	// let's give the sucker a hint
	std::cout
          << "Cbc takes input from arguments ( - switches to stdin)"
          << std::endl
          << "Enter ? for list of commands or help" << std::endl;
      } else {
	numberGoodCommands=1;
      }
    }
      
    while (1) {
      // Reset time
      time1 = CoinCpuTime();
      time1Elapsed = CoinGetTimeOfDay();

      field = CoinParamUtils::getNextField(inputQueue, interactiveMode, prompt);

      // exit if null or similar
      if (!field.length()) {
         if (numberGoodCommands == 1 && goodModel) {
            // we just had file name - do branch and bound
            field = "-solve";
         } else {
	   if (!partInputQueue.size())
	     break;
	   // switch back to original queue
	   inputQueue = partInputQueue;
	   std::deque<std::string> tempQueue;
	   partInputQueue = tempQueue;
	   continue;
         }
      }

      // These are some special case we deal with separately
      // (only in non-interactive mode)
      if (!interactiveMode){
         if (field == "-") {
            std::cout << "Switching to line mode" << std::endl;
            interactiveMode = true;
            while (!inputQueue.empty())
               inputQueue.pop_front();
	    field = CoinParamUtils::getNextField(inputQueue, interactiveMode, prompt);
         } else if (field[0] != '-') {
            // special dispensation - taken as -import name, put name back on queue
            inputQueue.push_front(field);
            field = "import";
         } else {
            if (field != "--") {
               // take off -
               field = field.substr(1);
	       if (field == "optionfile") {
		 if (!inputQueue.empty()) {
		   field = CoinParamUtils::getNextField(inputQueue, interactiveMode, prompt);
		   partInputQueue = inputQueue;
		   std::deque<std::string> tempQueue;
		   FILE * fp = fopen(field.c_str(),"r");
		   if (!fp) {
		     std::cout << "unable to open option file "
			       << field << std::endl;
		     std::deque<std::string> tempQueue;
		     partInputQueue = tempQueue;
		     continue;
		   }
		   /* format of file -
		      comments line starts with * or #
		      first character not blank
		      command is first - if no - then - added
		   */
		   char line[200];
		   std::cout << "Extra options - ";
		   while (fgets(line, 200, fp)) {
		     // skip comment
		     if (line[0]=='*'|| line[0]=='#')
		       continue;
		     int nchar = strlen(line);
		     if (nchar<2)
		       continue;
		     if (line[0]!='-') {
		       for (int i=nchar;i>=0;i--)
			 line[i+1] = line[i];
		       //memmove(line+1,line,nchar+1);
		       nchar++;
		       line[0]='-';
		     }
		     char *pos=line;
		     char *put=line;
		     while (true) {
		       while (*pos != '\n' && *pos != '\0'
			      && *pos != ' ' && *pos != '\t') 
			 pos++;
		       char save = *pos;
		       *pos = '\0';
		       if (strlen(put)) {
			 tempQueue.push_back(put);
			 std::cout << put << " ";
		       }
		       if (save == ' ' || save == '\t') {
			 pos++;
			 while (*pos == ' ' || *pos == '\t')
			   pos++;
			 put = pos;
		       } else {
			 break; // end of line
		       }
		     }
		   }
		   fclose(fp);
		   std::cout << std::endl;
   		   inputQueue = tempQueue;
		   continue;
		 } else if (partInputQueue.size()) {
		   inputQueue = partInputQueue;
		   continue;
		 }
	       }
            } else {
               // special dispensation - taken as -import --
               field = "import";
            }
         }
      }

      // find out if valid command
      // ? would it be better to do Cbc first ?
      int numberClpMatches(0), numberClpShortMatches(0), numberClpQuery(999);
      int clpParamCode = CoinParamUtils::lookupParam(field,
                                                     clpParameters.paramVec(),
                                                     &numberClpMatches,
                                                     &numberClpShortMatches,
                                                     &numberClpQuery);
      
      int numberCbcMatches(0), numberCbcShortMatches(0), numberCbcQuery(0);
      if (numberClpQuery>0)
	numberCbcQuery = 999; // so won't complain if no matches
      int cbcParamCode = CoinParamUtils::lookupParam(field,
                                                     parameters.paramVec(),
                                                     &numberCbcMatches,
                                                     &numberCbcShortMatches,
                                                     &numberCbcQuery);

      int numberMatches(numberClpMatches + numberCbcMatches);
      int numberShortMatches(numberClpShortMatches + numberCbcShortMatches);
      int numberQuery(numberClpQuery + numberCbcQuery);

      if (numberQuery > 0 || (numberShortMatches > 0 && !numberMatches)){
         continue;
      }
      if (!numberMatches) {
         if (!interactiveMode){
            std::cout << "Unrecognized parameter - " << field
                      << ", exiting..."
                      << std::endl;
            cbcParamCode = CbcParam::EXIT;
         } else {
            std::cout << "Unrecognized parameter - " << field
                      << " - enter valid command or end to exit"
                      << std::endl;
            continue;
	}
      }

      // Do some translation for backwards compatibility
      switch (cbcParamCode){
       case CbcParam::READMODEL_OLD:
         cbcParamCode = CbcParam::READMODEL;
         break;
       case CbcParam::WRITEGMPLSOL_OLD:
         cbcParamCode = CbcParam::WRITEGMPLSOL;
         break;
       case CbcParam::WRITEMODEL_OLD:
         cbcParamCode = CbcParam::WRITEMODEL;
         break;
       case CbcParam::WRITESOL_OLD:
         cbcParamCode = CbcParam::WRITESOL;
         break;
       case CbcParam::WRITESOLBINARY_OLD:
         cbcParamCode = CbcParam::WRITESOLBINARY;
         break;
       default:
         break;
      }
      
      CbcParam *cbcParam = parameters[CbcParam::INVALID];
      ClpParam *clpParam = clpParameters[ClpParam::INVALID];

      if (cbcParamCode >= 0){
         cbcParam = parameters[cbcParamCode];
         clpParamCode = ClpParam::INVALID;
      }else{
         clpParam = clpParameters[clpParamCode];
      }

#if 0
      // This logic is all captured in lookupParam
      if ((cbcParam->type() == CoinParam::paramInvalid &&
           clpParam->type() == CoinParam::paramInvalid) ||
          numberQuery) {
        if (!numberMatches) {
          std::cout << "No match for " << field << " - ? for list of commands"
                    << std::endl;
        } else if (numberMatches == 1) {
          if (!numberQuery) {
            std::cout << "Short match for " << field << " - completion: ";
            if (cbcParam->type() != CoinParam::paramInvalid) {
              std::cout << cbcParam->matchName() << std::endl;
            } else {
              std::cout << clpParam->matchName() << std::endl;
            }
          } else if (numberQuery) {
            if (cbcParam->type() != CoinParam::paramInvalid) {
              std::cout << cbcParam->matchName() << " : ";
              std::cout << cbcParam->shortHelp() << std::endl;
              if (numberQuery >= 2){
                cbcParam->printLongHelp();
              }
            } else {
              std::cout << clpParam->matchName() << " : ";
              std::cout << clpParam->shortHelp() << std::endl;
              if (numberQuery >= 2){
                clpParam->printLongHelp();
              }
            }
          }
        } else {
          if (!numberQuery)
            std::cout << "Multiple matches for " << field
                      << " - possible completions:" << std::endl;
          else
            std::cout << "Completions of " << field << ":" << std::endl;
          for (int iParam = CbcParam::FIRSTPARAM + 1;
               iParam < CbcParam::LASTPARAM; iParam++) {
            int match = parameters[iParam]->matches(field);
            if (match && parameters[iParam]->getDisplayPriority()) {
              std::cout << parameters[iParam]->matchName();
              if (numberQuery >= 2)
                std::cout << " : " << parameters[iParam]->shortHelp();
              std::cout << std::endl;
            }
          }
          for (int iParam = ClpParam::FIRSTPARAM + 1;
               iParam < ClpParam::LASTPARAM; iParam++) {
            int match = clpParameters[iParam]->matches(field);
            if (match && clpParameters[iParam]->getDisplayPriority()) {
              std::cout << clpParameters[iParam]->matchName();
              if (numberQuery >= 2)
                std::cout << " : " << clpParameters[iParam]->shortHelp();
              std::cout << std::endl;
            }
          }
        }
      } else {
#endif
        // found
        int status;
        numberGoodCommands++;
        if (cbcParamCode == CbcParam::BAB && goodModel) {
          if (clqstrMode == "before") { // performing clique strengthening
                                          // before initial solve
            CglCliqueStrengthening clqStr(model_.solver(),
                                          model_.messageHandler());
            clqStr.strengthenCliques(2);

            if (clqStr.constraintsExtended() + clqStr.constraintsDominated() >
                0) {
              model_.solver()->initialSolve();

              if (model_.solver()->isProvenPrimalInfeasible()) {
                 printGeneralMessage(model_,
                                     "Clique Strengthening says infeasible!");
              } else {
                 buffer.str("");
                 buffer << "After applying Clique Strengthening continuous "
                        << "objective value is " << model_.solver()->getObjValue();
                 printGeneralMessage(model_, buffer.str());
              }
            }
          }
#if CBC_USE_INITIAL_TIME == 1
          if (model_.useElapsedTime())
            model_.setDblParam(CbcModel::CbcStartSeconds, CoinGetTimeOfDay());
          else
            model_.setDblParam(CbcModel::CbcStartSeconds, CoinCpuTime());
#endif
          biLinearProblem = false;
          // check if any integers
#ifndef CBC_OTHER_SOLVER
          if (info && info->numberSos && doSOS && statusUserFunction_[0]) {
            // SOS
            numberSOS = info->numberSos;
          }
          lpSolver = clpSolver->getModelPtr();
          if (!lpSolver->integerInformation() && !numberSOS &&
              !clpSolver->numberSOS() && !model_.numberObjects() &&
              !clpSolver->numberObjects()) {
             clpParamCode = ClpParam::DUALSIMPLEX;
	     cbcParamCode = CbcParam::LASTACTIONPARAM; // acts like dummy
             clpParam = clpParameters[clpParamCode];
#ifdef CBC_MAXIMUM_BOUND
          } else {
            double *lower = lpSolver->columnLower();
            double *upper = lpSolver->columnUpper();
            int numberColumns = lpSolver->numberColumns();
            for (int i = 0; i < numberColumns; i++) {
              lower[i] = std::max(lower[i], -CBC_MAXIMUM_BOUND);
              upper[i] = std::min(upper[i], CBC_MAXIMUM_BOUND);
            }
#endif
          }
#endif
        }
        if (cbcParamCode == CbcParam::GENERALQUERY  ||
            cbcParamCode == CbcParam::FULLGENERALQUERY) {
          //TODO This should be made more transparent
          //bool evenHidden = false;
          int commandPrintLevel =
             parameters[CbcParam::COMMANDPRINTLEVEL]->modeVal();
          if (cbcParamCode == CbcParam::FULLGENERALQUERY) {
             verbose = 1;
             commandPrintLevel = 1;
          }	
#if 0
          if ((verbose & 8) != 0) {
            // even hidden
            // evenHidden = true;
            verbose &= ~8;
          }
          if (verbose < 4 && statusUserFunction_[0])
            verbose += 4;
#endif
          if (verbose) {
             std::cout << std::endl
                       << "Commands either invoke actions or set parameter values.\n"
                       << "When specifying multiple commands on one command line,\n"
                       << "parameter/action names should be prepended with a '-',\n"
                       << "followed by a value (some actions don't accept values as\n"
                       << "arguments). Specifying -stdin at anytime switches to stdin.\n"
                       << std::endl
                       << "In interactive mode, specify one command per line and\n"
                       << "don't prepend command names with '-'.\n"
                       << std::endl
                       << "Some actions take file names as arguments. If no file name\n"
                       << "is provided, then the previous name (or initial default)\n"
                       << "will be used.\n"
                       << std::endl
                       << "abcd? will list commands starting with 'abcd'.\n"
                       << "If there is only one match, a short explanation is given.\n"
                       << std::endl
                       << "abcd?? will list commands with explanations.\n"
                       << "If there is only one match, fuller help is given.\n"
                       << std::endl
                       << "abcd without value gives current value (for parameters).\n"
                       << "abcd 'value' sets value (for parameters)\n"
                       << std::endl
                       << "Commands are:" << std::endl;
          } else {
            std::cout << "Cbc options are set within AMPL with commands like:"
                      << std::endl
                      << std::endl;
            std::cout << "         option cbc_options \"cuts=root log=2 "
                         "feas=on slog=1\""
                      << std::endl
                      << std::endl;
            std::cout << "only maximize, dual, primal, help and quit are "
                         "recognized without ="
                      << std::endl;
          }
          // TODO Fix these!
          int maxAcross = 10;
          if ((verbose % 4) != 0){
             maxAcross = 1;
          }
          std::vector< std::string > types;
          types.push_back("Invalid parameters:");
          types.push_back("Action parameters:");
          types.push_back("Integer parameters:");
          types.push_back("Double parameters:");
          types.push_back("String parameters:");
	  types.push_back("Directory parameters:");
	  types.push_back("File parameters:");
          types.push_back("Keyword parameters:");
	  std::cout << std::endl << "#### Cbc Parameters ####" << std::endl;
	  for (int type = 1; type < 8; type++) {
	    int across = 0;
	    int lengthLine = 0;
	    bool first = true;
	    for (int iParam = CbcParam::FIRSTPARAM + 1;
                 iParam < CbcParam::LASTPARAM; iParam++) {
	      //TODO Print model parameters
              CbcParam *p = parameters[iParam];
	      if (p->type() != type ||  
                  p->getDisplayPriority() < commandPrintLevel){
                 continue;
              }
	      if (first) {
                std::cout << std::endl
                          << "*** " << types[type] << " ***"
                          << std::endl << std::endl;
		first = false;
	      }
#if 0
	      if ((verbose % 4) != 0) {
		std::cout << std::endl;
	      }
	      if ((verbose & 2) != 0){
		std::cout << std::endl;
	      }
              // TODO Fix AMPL mode stuff
              // but skip if not useful for ampl (and in ampl mode)
              if (verbose >= 4 && (p->whereUsed() & 4) == 0)
                 continue;
#endif
              int length = p->lengthMatchName() + 1;
              if (lengthLine + length > 80) {
                 std::cout << std::endl;
                 across = 0;
                 lengthLine = 0;
              }
              if (!across && (verbose & 2) != 0) {
		    std::cout << "Command ";
              }
              std::cout << p->matchName();
              lengthLine += length;
              across++;
              if (verbose) {
                 if ((verbose % 4) != 0) {
		    // put out description as well
		    if ((verbose & 1) != 0){
                       if (length < 8){
                          std::cout << "\t\t\t";
                       } else if (length < 16) {
                          std::cout << "\t\t";
                       } else {
                          std::cout << "\t";
                       }
                       std::cout << p->shortHelp();
                       std::cout << std::endl;
                    } else  if ((verbose & 2) != 0) {
                       std::cout << "---- description" << std::endl;
                       p->printLongHelp();
                       std::cout << "----" << std::endl << std::endl;
		    }
                 } 
                 across = 0;
                 lengthLine = 0;
              } else {
               std::cout << " ";
              }
              if (across == maxAcross) {
                 across = 0;
                 lengthLine = 0;
                 std::cout << std::endl;
              }
            }
            if (across){
               std::cout << std::endl;
            }
	  }
	  std::cout << std::endl << "#### Clp Parameters ####" << std::endl;
	  for (int type = 1;type < 8;type++) {
	    int across = 0;
	    int lengthLine = 0;
	    bool first = true;
	    for (int iParam = ClpParam::FIRSTPARAM+1;
		 iParam < ClpParam::LASTPARAM; iParam++) {
	      if (clpParameters[iParam]->type() != type)
		continue;
              ClpParam *p = clpParameters[iParam];
              if (p->type() != type ||  
                  p->getDisplayPriority() < commandPrintLevel){
                 continue;
              }
              if (first) {
                 std::cout << std::endl
                           << "*** " << types[type] << " ***"
                           << std::endl << std::endl;
                 first = false;
              }
              int length = p->lengthMatchName() + 1;
              if (lengthLine + length > 80) {
                 std::cout << std::endl;
                 across = 0;
                 lengthLine = 0;
              }
              if (!across && (verbose & 2) != 0){
                 std::cout << "Command ";
              }
              std::cout << p->matchName();
              lengthLine += length;
              across++;
              if (verbose) {
                 // put out description as well
                 if ((verbose & 1) != 0){
                    if (length < 8){
                       std::cout << "\t\t\t";
                    } else if (length < 16) {
                       std::cout << "\t\t";
                    } else {
                       std::cout << "\t";
                    }
                    std::cout << p->shortHelp();
                    std::cout << std::endl;
                 } else if ((verbose & 2) != 0) {
                    std::cout << "---- description" << std::endl;
                    p->printLongHelp();
                    std::cout << "----" << std::endl << std::endl;
                 }
                 across = 0;
                 lengthLine = 0;
              } else {
                 std::cout << " ";
              }
              if (across == maxAcross) {
                 across = 0;
                 lengthLine = 0;
                 std::cout << std::endl;
              }
            }
            if (across){
               std::cout << std::endl;
            }
          }
          continue;
        }

        if (cbcParam->type() == CoinParam::paramDbl) {
           if ((status = cbcParam->readValue(inputQueue, dValue, &message))){
              printGeneralMessage(model_, message);
              continue;
           }
           if (cbcParam->setVal(dValue, &message)){
              printGeneralMessage(model_, message);
              continue;
	   } else {
	     printGeneralMessage(model_, message);
           }
              // TODO: These should be moved to the push function
           switch (cbcParamCode) {
            case CbcParam::DJFIX:
              djFix = dValue;
#ifndef CBC_OTHER_SOLVER
              if (goodModel && djFix < 1.0e20) {
                 // do some fixing
                 clpSolver =
                    dynamic_cast<OsiClpSolverInterface *>(model_.solver());
                 clpSolver->initialSolve();
                 lpSolver = clpSolver->getModelPtr();
                 int numberColumns = lpSolver->numberColumns();
                 int i;
                 const char *type = lpSolver->integerInformation();
                 double *lower = lpSolver->columnLower();
                 double *upper = lpSolver->columnUpper();
                 double *solution = lpSolver->primalColumnSolution();
                 double *dj = lpSolver->dualColumnSolution();
                 int numberFixed = 0;
                 double dextra4 = parameters[CbcParam::DEXTRA4]->dblVal();
                 if (dextra4)
                    printf("Multiple for continuous dj fixing is %g\n",
                           dextra4);
                 for (i = 0; i < numberColumns; i++) {
                    double djValue = dj[i];
                    if (!type[i])
                       djValue *= dextra4;
                    if (type[i] || dextra4) {
                       double value = solution[i];
                       if (value < lower[i] + 1.0e-5 && djValue > djFix) {
                          solution[i] = lower[i];
                          upper[i] = lower[i];
                          numberFixed++;
                       } else if (value > upper[i] - 1.0e-5 &&
                                  djValue < -djFix) {
                          solution[i] = upper[i];
                          lower[i] = upper[i];
                          numberFixed++;
                       }
                    }
                 }
                 buffer.str("");
                 buffer << numberFixed << " columns fixed" << std::endl;
                 printGeneralMessage(model_, buffer.str());
              }
#endif
              break;
            case CbcParam::TIGHTENFACTOR:
              tightenFactor = dValue;
              if (!complicatedInteger)
                 defaultSettings = false; // user knows what she is doing
              break;
            default:
              break;
           }
        } else if (clpParam->type() == CoinParam::paramDbl) {
           if ((status = clpParam->readValue(inputQueue, dValue, &message))){
              printGeneralMessage(model_, message);
              continue;
           } 
           // We have a Clp parameter
           if (clpParam->setVal(dValue, &message)){
              printGeneralMessage(model_, message);
              continue;
           }      
        } else if (cbcParam->type() == CoinParam::paramInt){
           if ((status = cbcParam->readValue(inputQueue, iValue, &message))){
              printGeneralMessage(model_, message);
              continue;
           }
           if (cbcParamCode == CbcParam::LPLOGLEVEL){
              clpSolver->messageHandler()->setLogLevel(iValue); // as well
              clpParameters[ClpParam::LOGLEVEL]->setVal(iValue, &message);
              cbcParam->setVal(iValue, &message);
              printGeneralMessage(model_, message);
           }else{
              cbcParam->setVal(iValue, &message);
              printGeneralMessage(model_, message);
              if (cbcParamCode == CbcParam::CUTPASS){
                 cutPass = iValue;
              } else if (cbcParamCode== CbcParam::USESOLUTION) {
                 useSolution = iValue;
              } else if (cbcParamCode == CbcParam::CUTPASSINTREE) {
                 cutPassInTree = iValue;
              } else if (cbcParamCode == CbcParam::STRONGBRANCHING ||
                         cbcParamCode == CbcParam::NUMBERBEFORE) {
                 strongChanged = true;
              } else if (cbcParamCode == CbcParam::FPUMPTUNE || cbcParamCode == CbcParam::FPUMPTUNE2 ||
                         cbcParamCode == CbcParam::FPUMPITS) {
                 pumpChanged = true;
              } else if (cbcParamCode == CbcParam::AGGREGATEMIXED) {
                 mixedRoundStrategy = iValue;
              } else if (cbcParamCode == CbcParam::BKPIVOTINGSTRATEGY) {
                 bkPivotingStrategy = static_cast<CoinBronKerbosch::PivotingStrategy>(iValue);
              } else if (cbcParamCode == CbcParam::BKMAXCALLS) {
                 maxCallsBK = iValue;
              } else if (cbcParamCode == CbcParam::BKCLQEXTMETHOD) {
                 bkClqExtMethod = iValue;
              } else if (cbcParamCode == CbcParam::ODDWEXTMETHOD) {
                 oddWExtMethod = iValue;
              } else if (cbcParamCode == CbcParam::LOGLEVEL) {
                 model_.messageHandler()->setLogLevel(std::abs(iValue));
              } else if (cbcParamCode == CbcParam::MAXNODES) {
                 model_.setIntParam(CbcModel::CbcMaxNumNode, iValue);
              } else if (cbcParamCode == CbcParam::MAXSOLS) {
                 model_.setIntParam(CbcModel::CbcMaxNumSol, iValue);
              } else if (cbcParamCode == CbcParam::STRONGBRANCHING) {
                 model_.setNumberStrong(iValue);
              } else if (cbcParamCode == CbcParam::PROCESSTUNE) {
                 tunePreProcess = iValue;
	      } else if (cbcParamCode == CbcParam::PRINTOPTIONS) {
		printOptions = iValue;
              } else if (cbcParamCode == CbcParam::VERBOSE) {
                 verbose = iValue;
              } else if (cbcParamCode == CbcParam::EXPERIMENT && iValue < 10000) {
                 int addFlags = 0;
                 // switch on some later features if >999
                 if (iValue > 999) {
                    int switchValue = iValue / 1000;
                    iValue -= 1000 * switchValue;
                    parameters[CbcParam::EXPERIMENT]->setVal(0 /*value*/);
                    switch (switchValue) {
                     default:
                     case 4:
                       // hotstart 500, -200 cut passes
                       parameters[CbcParam::MAXHOTITS]->setVal(500, &message);
                       printGeneralMessage(model_, message);
                       parameters[CbcParam::CUTPASS]->setVal(-200, &message);
                       printGeneralMessage(model_, message);
                     case 3:
                       // multiple 4
                       parameters[CbcParam::MULTIPLEROOTS]->setVal(4, &message);
                       printGeneralMessage(model_, message);
                     case 2:
                       // rens plus all diving at root
                       parameters[CbcParam::DIVEOPT]->setVal(16, &message);
                       printGeneralMessage(model_, message);
                       model_.setNumberAnalyzeIterations(-iValue);
                       // -tune 7 zero,lagomory,gmi at root - probing on
                     case 1:
                       tunePreProcess = 7;
                       parameters[CbcParam::PROCESSTUNE]->setVal(7, &message);
                       printGeneralMessage(model_, message);
                       // message =
                       // parameters[CbcParam::MIPOPTIONS,
                       // parameters)]->setValWithMessage(1025); if
                       // (!noPrinting_&&message)
                       //    generalMessageHandler->message(CLP_GENERAL,
                       //    generalMessages)
                       //  << message << CoinMessageEol;
                       parameters[CbcParam::PROBINGCUTS]->setVal("on", &message);
                       probingMode = CbcParameters::CGOn;
                       printGeneralMessage(model_, message);
                       parameters[CbcParam::ZEROHALFCUTS]->setVal("root", &message);
                       printGeneralMessage(model_, message);
                       parameters[CbcParam::LAGOMORYCUTS]->setVal("root", &message);
                       printGeneralMessage(model_, message);
                       GMIMode = CbcParameters::CGRoot;
                       parameters[CbcParam::GMICUTS]->setVal("root", &message);
                       printGeneralMessage(model_, message);
                    }
                    iValue = 0;
                 }
                 if (iValue >= 10) {
                    addFlags = 1048576 * (iValue / 10);
                    iValue = iValue % 10;
                    parameters[CbcParam::EXPERIMENT]->setVal(iValue);
                 }
#ifndef CBC_EXPERIMENT_JJF
                 if (iValue == 1) {
                    // just experimental preprocessing and more restarts
                    tunePreProcess |= 8192;
                    model_.setSpecialOptions(model_.specialOptions() |
                                             (512 | 32768));
                 }
#else
                 experimentValue = iValue; // save
                 if (iValue > 0 && iValue < 5) {
                    // more restarts
                    // >1 go to end in strong branching
                    // >2 experimental preprocessing
                    // 4 try ranging before strong branching
                    model_.setSpecialOptions(model_.specialOptions() |
                                             (512 | 32768));
#ifndef CBC_OTHER_SOLVER
                    if (iValue > 1) {
                       OsiClpSolverInterface *osiclp =
                          dynamic_cast<OsiClpSolverInterface *>(model_.solver());
                       osiclp->setSpecialOptions(osiclp->specialOptions() & ~32);
                    }
#endif
                    if (iValue > 2)
                       tunePreProcess |= 8198; // was 8199;
                    if (iValue == 4) {
                       int more2 =
                          parameters[CbcParam::MOREMOREMIPOPTIONS]->intVal();
                       parameters[CbcParam::MOREMOREMIPOPTIONS]->setVal(more2 | 1048576);
                    }
                    iValue = 1;
                 }
#endif
                 if (iValue > 1) {
                    int values[] = {24003, 280003, 792003, 24003, 24003};
                    if (iValue >= 2 && iValue <= 3) {
                       // swap default diving
                       parameters[CbcParam::DIVINGC]->setVal("off");
                       parameters[CbcParam::DIVINGP]->setVal("on");
                    }
                    int extra4 = values[iValue - 1] + addFlags;
                    parameters[CbcParam::EXTRA4]->setVal(extra4);
                    buffer.str("");
                    buffer << "switching on global root cuts for gomory and "
                           << "knapsack" << std::endl;
                    //buffer << "using OSL factorization" << std::endl;
                    buffer << "extra options - -rens on -extra4 " << extra4
                           << " -passc 1000!" << std::endl;
                    printGeneralMessage(model_, buffer.str());
                    parameters[CbcParam::PROBINGCUTS]->setVal("forceOnStrong");
                    probingMode = CbcParameters::CGForceOnStrong;
                    parameters[CbcParam::GOMORYCUTS]->setVal("onGlobal");
                    gomoryMode = CbcParameters::CGOnGlobal;
                    parameters[CbcParam::KNAPSACKCUTS]->setVal("onGlobal");
                    knapsackMode = CbcParameters::CGOnGlobal;
                    clpParameters[ClpParam::FACTORIZATION]->setVal("osl");
                    //lpSolver->factorization()->forceOtherFactorization(3);
                    parameters[CbcParam::MAXHOTITS]->setVal(100);
                    parameters[CbcParam::CUTPASS]->setVal(1000);
                    cutPass = 1000;
                    parameters[CbcParam::RENS]->setVal("on");
                 }
              } else if (cbcParamCode == CbcParam::STRATEGY) {
                 if (iValue == 0) {
                    gomoryGen.setAwayAtRoot(0.05);
                    parameters[CbcParam::DIVEOPT]->setVal(-1);
                    parameters[CbcParam::FPUMPITS]->setVal(20);
                    parameters[CbcParam::FPUMPTUNE]->setVal(1003);
                    initialPumpTune = 1003;
                    parameters[CbcParam::PROCESSTUNE]->setVal(0);
                    tunePreProcess = 0;
                    parameters[CbcParam::RINS]->setVal("off");
                    // but not if cuts off
                    if (parameters[CbcParam::CUTSTRATEGY]->modeVal()){
                       parameters[CbcParam::PROBINGCUTS]->setVal("on");
                       probingMode = CbcParameters::CGOn;
                    } else {
                       parameters[CbcParam::PROBINGCUTS]->setVal("off");
                       probingMode = CbcParameters::CGOff;
                    }
                 }
              }
           }
        } else if (clpParam->type() == CoinParam::paramInt) {
           if ((status = clpParam->readValue(inputQueue, iValue, &message))){
              printGeneralMessage(model_, message);
              continue;
           } 
           if (clpParam->setVal(iValue, &message)){
              printGeneralMessage(model_, message);
              continue;
	   } else {
	     printGeneralMessage(model_, message);
           }
           if (clpParamCode == ClpParam::PRESOLVEPASS) {
              preSolve = iValue;
           } else if (clpParamCode == ClpParam::IDIOT) {
              doIdiot = iValue;
           } else if (clpParamCode == ClpParam::SPRINT) {
              doSprint = iValue;
           } else if (clpParamCode == ClpParam::OUTPUTFORMAT) {
              outputFormat = iValue;
           } else if (clpParamCode == ClpParam::SLPVALUE) {
              slpValue = iValue;
           } else if (clpParamCode == ClpParam::CPP) {
              cppValue = iValue;
           } else if (clpParamCode == ClpParam::PRESOLVEOPTIONS) {
              presolveOptions = iValue;
           } else if (clpParamCode == ClpParam::PRINTOPTIONS) {
              printOptions = iValue;
           } else if (clpParamCode == ClpParam::SUBSTITUTION) {
              substitution = iValue;
           } else if (clpParamCode == ClpParam::DUALIZE) {
              dualize = iValue;
           }
        } else if (cbcParam->type() == CoinParam::paramKwd) {
           if ((status = cbcParam->readValue(inputQueue, field, &message))){
              printGeneralMessage(model_, message);
	      continue;
           }
           if (cbcParam->setVal(field, &message)){
              printGeneralMessage(model_, message);
              continue;
           }
           int mode = cbcParam->modeVal();
           // TODO: this should be part of the push method
           switch (cbcParamCode) {
            case CbcParam::DIRECTION:
              if (mode == 0)
                 lpSolver->setOptimizationDirection(1);
              else if (mode == 1)
                 lpSolver->setOptimizationDirection(-1);
              else
                 lpSolver->setOptimizationDirection(0);
              break;
            case CbcParam::ERRORSALLOWED:
              allowImportErrors = mode;
              break;
            case CbcParam::INTPRINT:
              printMode = mode;
              break;
              // case CbcParam::CLP_PARAM_NOTUSED_ALGORITHM:
              // algorithm  = mode;
              // defaultSettings=false; // user knows what she is doing
              // abort();
              // break;
            case CbcParam::SOS:
              doSOS = mode;
              break;
            case CbcParam::CLQSTRENGTHENING:
              clqstrMode = field;
              break;
            case CbcParam::USECGRAPH:
              cgraphMode = field;
              break;
            case CbcParam::CLIQUECUTS:
              defaultSettings = false; // user knows what she is doing
              cliqueMode = mode;
	      oldCliqueMode = cliqueMode;
              break;
            case CbcParam::ODDWHEELCUTS:
              defaultSettings = false; // user knows what she is doing
              oddWheelMode = mode;
              break;
            case CbcParam::GOMORYCUTS:
              defaultSettings = false; // user knows what she is doing
              gomoryMode = mode;
              break;
            case CbcParam::PROBINGCUTS:
              defaultSettings = false; // user knows what she is doing
              probingMode = mode;
              break;
            case CbcParam::KNAPSACKCUTS:
              defaultSettings = false; // user knows what she is doing
              knapsackMode = mode;
              break;
            case CbcParam::REDSPLITCUTS:
              defaultSettings = false; // user knows what she is doing
              redsplitMode = mode;
              break;
            case CbcParam::REDSPLIT2CUTS:
              defaultSettings = false; // user knows what she is doing
              redsplit2Mode = mode;
              break;
            case CbcParam::GMICUTS:
              defaultSettings = false; // user knows what she is doing
              GMIMode = mode;
              break;
            case CbcParam::FLOWCUTS:
              defaultSettings = false; // user knows what she is doing
              flowMode = mode;
              break;
            case CbcParam::MIRCUTS:
              defaultSettings = false; // user knows what she is doing
              mixedMode = mode;
              break;
            case CbcParam::TWOMIRCUTS:
              defaultSettings = false; // user knows what she is doing
              twomirMode = mode;
              break;
            case CbcParam::LANDPCUTS:
              defaultSettings = false; // user knows what she is doing
              landpMode = mode;
              break;
            case CbcParam::RESIDCAPCUTS:
              defaultSettings = false; // user knows what she is doing
              residualCapacityMode = mode;
              break;
            case CbcParam::ZEROHALFCUTS:
              defaultSettings = false; // user knows what she is doing
              zerohalfMode = mode;
              break;
            case CbcParam::ROUNDING:
              defaultSettings = false; // user knows what she is doing
              break;
            case CbcParam::FPUMP:
              defaultSettings = false; // user knows what she is doing
              break;
            case CbcParam::RINS:
              break;
            case CbcParam::DINS:
              break;
            case CbcParam::RENS:
              break;
            case CbcParam::CUTSTRATEGY:
              gomoryMode = mode;
              probingMode = mode;
              knapsackMode = mode;
              cliqueMode = mode;
              flowMode = mode;
              mixedMode = mode;
              twomirMode = mode;
              zerohalfMode = mode;
              oddWheelMode = mode;
              parameters[CbcParam::GOMORYCUTS]->setVal(mode);
              parameters[CbcParam::PROBINGCUTS]->setVal(mode);
              parameters[CbcParam::KNAPSACKCUTS]->setVal(mode);
              parameters[CbcParam::CLIQUECUTS]->setVal(mode);
              parameters[CbcParam::FLOWCUTS]->setVal(mode);
              parameters[CbcParam::MIRCUTS]->setVal(mode);
              parameters[CbcParam::TWOMIRCUTS]->setVal(mode);
              parameters[CbcParam::ZEROHALFCUTS]->setVal(mode);
              parameters[CbcParam::ODDWHEELCUTS]->setVal(mode);
              if (!mode) {
                 // switch off clique strengthening
                 clqstrMode = "off";
		 cgraphMode = "off";
              }
              break;
            case CbcParam::HEURISTICSTRATEGY:
              parameters[CbcParam::ROUNDING]->setVal(mode);
              parameters[CbcParam::GREEDY]->setVal(mode);
              parameters[CbcParam::COMBINE]->setVal(mode);
              // parameters[CbcParam::LOCALTREE,numberParameters_,parameters)]->setVal(mode);
              parameters[CbcParam::FPUMP]->setVal(mode);
              parameters[CbcParam::DIVINGC]->setVal(mode);
              parameters[CbcParam::RINS]->setVal(mode);
              break;
            case CbcParam::GREEDY:
            case CbcParam::DIVINGS:
            case CbcParam::DIVINGC:
            case CbcParam::DIVINGF:
            case CbcParam::DIVINGG:
            case CbcParam::DIVINGL:
            case CbcParam::DIVINGP:
            case CbcParam::DIVINGV:
            case CbcParam::COMBINE:
            case CbcParam::PIVOTANDCOMPLEMENT:
            case CbcParam::PIVOTANDFIX:
            case CbcParam::RANDROUND:
            case CbcParam::LOCALTREE:
            case CbcParam::NAIVE:
            case CbcParam::CPX:
              defaultSettings = false; // user knows what she is doing
              break;
            case CbcParam::BRANCHPRIORITY:
              useCosts = mode;
              break;
            case CbcParam::NODESTRATEGY:
              nodeStrategy = mode;
              break;
            case CbcParam::PREPROCESS:
              preProcess = mode;
              break;
            case CbcParam::TIMEMODE:
              model_.setUseElapsedTime(mode != 0);
              break;
            default:
              // abort();
              break;
           }
        } else if (clpParam->type() == CoinParam::paramKwd) { 
           if ((status = clpParam->readValue(inputQueue, field, &message))){
              printGeneralMessage(model_, message);
              continue;
           }
           if (clpParam->setVal(field, &message)){
              printGeneralMessage(model_, message);
              continue;
           }
           int mode = clpParam->kwdToMode(field);
           // for now hard wired
           switch (clpParamCode) {
            case ClpParam::DUALPIVOT:
              if (mode == 0) {
                 ClpDualRowSteepest steep(3);
                 lpSolver->setDualRowPivotAlgorithm(steep);
              } else if (mode == 1) {
                 ClpDualRowDantzig dantzig;
                 // ClpDualRowSteepest dantzig(5);
                 lpSolver->setDualRowPivotAlgorithm(dantzig);
              } else if (mode == 2) {
                 // partial steep
                 ClpDualRowSteepest steep(2);
                 lpSolver->setDualRowPivotAlgorithm(steep);
              } else if (mode == 3) {
                 ClpDualRowSteepest steep;
                 lpSolver->setDualRowPivotAlgorithm(steep);
              } else if (mode == 4) {
                 // Positive edge steepest
                 ClpPEDualRowSteepest p(fabs(clpParameters[ClpParam::PSI]->dblVal()));
                 lpSolver->setDualRowPivotAlgorithm(p);
              } else if (mode == 5) {
                 // Positive edge Dantzig
                 ClpPEDualRowDantzig p(fabs(clpParameters[ClpParam::PSI]->dblVal()));
                 lpSolver->setDualRowPivotAlgorithm(p);
              }
              break;
            case ClpParam::PRIMALPIVOT:
              if (mode == 0) {
                 ClpPrimalColumnSteepest steep(3);
                 lpSolver->setPrimalColumnPivotAlgorithm(steep);
              } else if (mode == 1) {
                 ClpPrimalColumnSteepest steep(0);
                 lpSolver->setPrimalColumnPivotAlgorithm(steep);
              } else if (mode == 2) {
                 ClpPrimalColumnDantzig dantzig;
                 lpSolver->setPrimalColumnPivotAlgorithm(dantzig);
              } else if (mode == 3) {
                 ClpPrimalColumnSteepest steep(4);
                 lpSolver->setPrimalColumnPivotAlgorithm(steep);
              } else if (mode == 4) {
                 ClpPrimalColumnSteepest steep(1);
                 lpSolver->setPrimalColumnPivotAlgorithm(steep);
              } else if (mode == 5) {
                 ClpPrimalColumnSteepest steep(2);
                 lpSolver->setPrimalColumnPivotAlgorithm(steep);
              } else if (mode == 6) {
                 ClpPrimalColumnSteepest steep(10);
                 lpSolver->setPrimalColumnPivotAlgorithm(steep);
              } else if (mode == 7) {
                 // Positive edge steepest
                 ClpPEPrimalColumnSteepest p(fabs(clpParameters[ClpParam::PSI]->dblVal()));
                 lpSolver->setPrimalColumnPivotAlgorithm(p);
              } else if (mode == 8) {
                 // Positive edge Dantzig
                 ClpPEPrimalColumnDantzig p(fabs(clpParameters[ClpParam::PSI]->dblVal()));
                 lpSolver->setPrimalColumnPivotAlgorithm(p);
              }
              break;
            case ClpParam::SCALING:
              lpSolver->scaling(mode);
              solver->setHintParam(OsiDoScale, mode != 0, OsiHintTry);
              doScaling = mode;
              break;
            case ClpParam::AUTOSCALE:
              lpSolver->setAutomaticScaling(mode != 0);
              break;
            case ClpParam::SPARSEFACTOR:
              lpSolver->setSparseFactorization((1 - mode) != 0);
              break;
            case ClpParam::BIASLU:
              lpSolver->factorization()->setBiasLU(mode);
              break;
            case ClpParam::PERTURBATION:
              if (mode == 1) {
		if (lpSolver->perturbation()==100)
		  lpSolver->setPerturbation(50);
              } else {
                 lpSolver->setPerturbation(100);
		 clpParam->setVal(100);
	      }
              break;
            case ClpParam::KEEPNAMES:
              keepImportNames = mode;
              break;
            case ClpParam::PRESOLVE:
              if (mode == 0)
                 preSolve = 5;
              else if (mode == 1)
                 preSolve = 0;
              else if (mode == 2)
                 preSolve = 10;
              else
                 preSolveFile = true;
              break;
            case ClpParam::PFI:
              lpSolver->factorization()->setForrestTomlin(mode == 0);
              break;
            case ClpParam::FACTORIZATION:
              lpSolver->factorization()->forceOtherFactorization(mode);
              break;
            case ClpParam::CRASH:
              doCrash = mode;
              break;
            case ClpParam::VECTOR:
              doVector = mode;
              break;
            case ClpParam::MESSAGES:
              lpSolver->messageHandler()->setPrefix(mode != 0);
              break;
            case ClpParam::CHOLESKY:
              choleskyType = mode;
              break;
            case ClpParam::GAMMA:
              gamma = mode;
              break;
            case ClpParam::BARRIERSCALE:
              scaleBarrier = mode;
              break;
            case ClpParam::KKT:
              doKKT = mode;
              break;
            case ClpParam::CROSSOVER:
              crossover = mode;
              break;
            default:
              // abort();
              break;
           }
        } else if (cbcParam->type() == CoinParam::paramDir){
           if ((status = cbcParam->readValue(inputQueue, field, &message))){
              printGeneralMessage(model_, message);
              continue;
           }
           if (cbcParam->setVal(field, &message)){
              printGeneralMessage(model_, message);
               continue;
           }
        } else if (cbcParam->type() == CoinParam::paramFile){
           if ((status = cbcParam->readValue(inputQueue, field, &message))){
              printGeneralMessage(model_, message);
              continue;
           }
           if (cbcParam->setVal(field, &message)){
              printGeneralMessage(model_, message);
               continue;
           }
        } else  {
           // action
           if (cbcParamCode == CbcParam::EXIT ||
               cbcParamCode == CbcParam::STOP ||
               cbcParamCode == CbcParam::QUIT ||
               cbcParamCode == CbcParam::END) {
              if (statusUserFunction_[0]) {
                 if (info && (info->numberIntegers || info->numberBinary)) {
                    // integer
                 } else {
                    // linear
                 }
                 writeAmpl(info);
                 freeArrays2(info);
                 freeArgs(info);
              }
              break; // stop all
           }
           
           switch (clpParamCode) {
            case ClpParam::DUALSIMPLEX:
            case ClpParam::PRIMALSIMPLEX:
            case ClpParam::BARRIER:{
              if (!goodModel) {
               printGeneralWarning(model_, "** Current model not valid\n");
               continue;
              }
              // Say not in integer
              integerStatus = -1;
              double objScale = clpParameters[ClpParam::OBJSCALE2]->dblVal();
              // deal with positive edge
              double psi = clpParameters[ClpParam::PSI]->dblVal();
              if (psi > 0.0) {
                 ClpDualRowPivot *dualp = lpSolver->dualRowPivot();
                 ClpDualRowSteepest *d1 =
                    dynamic_cast<ClpDualRowSteepest *>(dualp);
                 ClpDualRowDantzig *d2 =
                    dynamic_cast<ClpDualRowDantzig *>(dualp);
                 if (d1) {
                    ClpPEDualRowSteepest p(psi, d1->mode());
                    lpSolver->setDualRowPivotAlgorithm(p);
                 } else if (d2) {
                    ClpPEDualRowDantzig p(psi);
                    lpSolver->setDualRowPivotAlgorithm(p);
                 }
                 ClpPrimalColumnPivot *primalp = lpSolver->primalColumnPivot();
                 ClpPrimalColumnSteepest *p1 =
                    dynamic_cast<ClpPrimalColumnSteepest *>(primalp);
                 ClpPrimalColumnDantzig *p2 =
                    dynamic_cast<ClpPrimalColumnDantzig *>(primalp);
                 if (p1) {
                    ClpPEPrimalColumnSteepest p(psi, p1->mode());
                    lpSolver->setPrimalColumnPivotAlgorithm(p);
                 } else if (p2) {
                    ClpPEPrimalColumnDantzig p(psi);
                    lpSolver->setPrimalColumnPivotAlgorithm(p);
                 }
              }
              if (objScale != 1.0) {
                 int iColumn;
                 int numberColumns = lpSolver->numberColumns();
                 double *dualColumnSolution = lpSolver->dualColumnSolution();
                 ClpObjective *obj = lpSolver->objectiveAsObject();
                 assert(dynamic_cast<ClpLinearObjective *>(obj));
                 double offset;
                 double *objective = obj->gradient(NULL, NULL, offset, true);
                 for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                    dualColumnSolution[iColumn] *= objScale;
                    objective[iColumn] *= objScale;
                    ;
                 }
                 int iRow;
                 int numberRows = lpSolver->numberRows();
                 double *dualRowSolution = lpSolver->dualRowSolution();
                 for (iRow = 0; iRow < numberRows; iRow++)
                    dualRowSolution[iRow] *= objScale;
                 lpSolver->setObjectiveOffset(objScale *
                                              lpSolver->objectiveOffset());
              }
              ClpSolve::SolveType method;
              ClpSolve::PresolveType presolveType;
              ClpSimplex *model2 = lpSolver;
              if (dualize) {
                 bool tryIt = true;
                 double fractionColumn = 1.0;
                 double fractionRow = 1.0;
                 if (dualize == 3) {
                    dualize = 1;
                    int numberColumns = lpSolver->numberColumns();
                    int numberRows = lpSolver->numberRows();
                    if (numberRows < 50000 || 5 * numberColumns > numberRows) {
                       tryIt = false;
                    } else {
                       fractionColumn = 0.1;
                       fractionRow = 0.1;
                    }
                 }
                 if (tryIt) {
                    model2 = static_cast<ClpSimplexOther *>(model2)->dualOfModel(
                                                                                 fractionRow, fractionColumn);
                    if (model2) {
                       buffer.str("");
                       buffer << "Dual of model has " << model2->numberRows()
                              << " rows and " << model2->numberColumns()
                              << " columns";
                       printGeneralMessage(model_, buffer.str());
                       model2->setOptimizationDirection(1.0);
                    } else {
                       model2 = lpSolver;
                       dualize = 0;
                    }
                 } else {
                    dualize = 0;
                 }
              }
              if (parameters.noPrinting()){
                 lpSolver->setLogLevel(0);
              }
              ClpSolve solveOptions;
              solveOptions.setPresolveActions(presolveOptions);
              solveOptions.setSubstitution(substitution);
              if (preSolve != 5 && preSolve) {
                 presolveType = ClpSolve::presolveNumber;
                 if (preSolve < 0) {
                    preSolve = -preSolve;
                    if (preSolve <= 100) {
                       presolveType = ClpSolve::presolveNumber;
                       buffer.str("");
                       buffer << "Doing " << preSolve
                              << " presolve passes - picking up non-costed slacks";
                       printGeneralMessage(model_, buffer.str());
                       solveOptions.setDoSingletonColumn(true);
                    } else {
                       preSolve -= 100;
                       presolveType = ClpSolve::presolveNumberCost;
                       buffer.str("");
                       buffer << "Doing " << preSolve
                              << " presolve passes - picking up non-costed slacks";
                       printGeneralMessage(model_, buffer.str());
                    }
                 }
              } else if (preSolve) {
                 presolveType = ClpSolve::presolveOn;
              } else {
                 presolveType = ClpSolve::presolveOff;
              }
              solveOptions.setPresolveType(presolveType, preSolve);
              if (clpParamCode == ClpParam::DUALSIMPLEX) {
                 method = ClpSolve::useDual;
              } else if (clpParamCode == ClpParam::PRIMALSIMPLEX) {
                 method = ClpSolve::usePrimalorSprint;
              } else {
                 method = ClpSolve::useBarrier;
                 if (crossover == 1) {
                    method = ClpSolve::useBarrierNoCross;
                 } else if (crossover == 2) {
                    ClpObjective *obj = lpSolver->objectiveAsObject();
                    if (obj->type() > 1) {
                       method = ClpSolve::useBarrierNoCross;
                       presolveType = ClpSolve::presolveOff;
                       solveOptions.setPresolveType(presolveType, preSolve);
                    }
                 }
              }
              solveOptions.setSolveType(method);
              if (preSolveFile)
                 presolveOptions |= 0x40000000;
              solveOptions.setSpecialOption(4, presolveOptions);
              solveOptions.setSpecialOption(5, printOptions);
              if (doVector) {
                 ClpMatrixBase *matrix = lpSolver->clpMatrix();
                 if (dynamic_cast<ClpPackedMatrix *>(matrix)) {
                    ClpPackedMatrix *clpMatrix =
                       dynamic_cast<ClpPackedMatrix *>(matrix);
                    clpMatrix->makeSpecialColumnCopy();
                 }
              }
              if (method == ClpSolve::useDual) {
                // dual
                if (doCrash)
                  solveOptions.setSpecialOption(0, 1, doCrash); // crash
                else if (doIdiot)
                  solveOptions.setSpecialOption(0, 2,
                                                doIdiot); // possible idiot
              } else if (method == ClpSolve::usePrimalorSprint) {
                // primal
                // if slp turn everything off
                if (slpValue > 0) {
                  doCrash = false;
                  doSprint = 0;
                  doIdiot = -1;
                  solveOptions.setSpecialOption(1, 10, slpValue); // slp
                  method = ClpSolve::usePrimal;
                }
                if (doCrash) {
                  solveOptions.setSpecialOption(1, 1, doCrash); // crash
                } else if (doSprint > 0) {
                  // sprint overrides idiot
                  solveOptions.setSpecialOption(1, 3, doSprint); // sprint
                } else if (doIdiot > 0) {
                  solveOptions.setSpecialOption(1, 2, doIdiot); // idiot
                } else if (slpValue <= 0) {
                  if (doIdiot == 0) {
                    if (doSprint == 0)
                      solveOptions.setSpecialOption(1, 4); // all slack
                    else
                      solveOptions.setSpecialOption(1,
                                                    9); // all slack or sprint
                  } else {
                    if (doSprint == 0)
                      solveOptions.setSpecialOption(1, 8); // all slack or idiot
                    else
                      solveOptions.setSpecialOption(1, 7); // initiative
                  }
                }
                if (basisHasValues == -1)
                  solveOptions.setSpecialOption(1, 11); // switch off values
              } else if (method == ClpSolve::useBarrier ||
                         method == ClpSolve::useBarrierNoCross) {
                int barrierOptions = choleskyType;
                if (scaleBarrier)
                  barrierOptions |= 8;
                if (doKKT)
                  barrierOptions |= 16;
                if (gamma)
                  barrierOptions |= 32 * gamma;
                if (crossover == 3)
                  barrierOptions |= 256; // try presolve in crossover
                solveOptions.setSpecialOption(4, barrierOptions);
              }
              model2->setMaximumSeconds(model_.getMaximumSeconds());
#ifdef COIN_HAS_LINK
              OsiSolverInterface *coinSolver = model_.solver();
              OsiSolverLink *linkSolver =
                  dynamic_cast<OsiSolverLink *>(coinSolver);
              if (!linkSolver) {
                model2->initialSolve(solveOptions);
              } else {
                // special solver
                int testOsiOptions = parameters[CbcParam::TESTOSI]->intVal();
                double *solution = NULL;
                if (testOsiOptions < 10) {
                  solution = linkSolver->nonlinearSLP(
                      slpValue > 0 ? slpValue : 20, 1.0e-5);
                } else if (testOsiOptions >= 10) {
                  CoinModel coinModel = *linkSolver->coinModel();
                  ClpSimplex *tempModel = approximateSolution(
                      coinModel, slpValue > 0 ? slpValue : 50, 1.0e-5, 0);
                  assert(tempModel);
                  solution = CoinCopyOfArray(tempModel->primalColumnSolution(),
                                             coinModel.numberColumns());
                  model2->setObjectiveValue(tempModel->objectiveValue());
                  model2->setProblemStatus(tempModel->problemStatus());
                  model2->setSecondaryStatus(tempModel->secondaryStatus());
                  delete tempModel;
                }
                if (solution) {
                  memcpy(model2->primalColumnSolution(), solution,
                         std::min(model2->numberColumns(),
                                 linkSolver->coinModel()->numberColumns()) *
                             sizeof(double));
                  delete[] solution;
                } else {
                  printf("No nonlinear solution\n");
                }
              }
#else
              model2->initialSolve(solveOptions);
#endif
              {
                // map states
                /* clp status
                                   -1 - unknown e.g. before solve or if
                   postSolve says not optimal 0 - optimal 1 - primal infeasible
                                   2 - dual infeasible
                                   3 - stopped on iterations or time
                                   4 - stopped due to errors
                                   5 - stopped by event handler (virtual int
                   ClpEventHandler::event()) */
                /* cbc status
                                   -1 before branchAndBound
                                   0 finished - check isProvenOptimal or
                   isProvenInfeasible to see if solution found (or check value
                   of best solution) 1 stopped - on maxnodes, maxsols, maxtime
                                   2 difficulties so run was abandoned
                                   (5 event user programmed event occurred) */
                /* clp secondary status of problem - may get extended
                                   0 - none
                                   1 - primal infeasible because dual limit
                   reached OR probably primal infeasible but can't prove it
                   (main status 4) 2 - scaled problem optimal - unscaled problem
                   has primal infeasibilities 3 - scaled problem optimal -
                   unscaled problem has dual infeasibilities 4 - scaled problem
                   optimal - unscaled problem has primal and dual
                   infeasibilities 5 - giving up in primal with flagged
                   variables 6 - failed due to empty problem check 7 - postSolve
                   says not optimal 8 - failed due to bad element check 9 -
                   status was 3 and stopped on time 100 up - translation of enum
                   from ClpEventHandler
                                */
                /* cbc secondary status of problem
                                   -1 unset (status_ will also be -1)
                                   0 search completed with solution
                                   1 linear relaxation not feasible (or worse
                   than cutoff) 2 stopped on gap 3 stopped on nodes 4 stopped on
                   time 5 stopped on user event 6 stopped on solutions 7 linear
                   relaxation unbounded 8 stopped on iterations limit
                                */
                int iStatus = model2->status();
                int iStatus2 = model2->secondaryStatus();
                if (iStatus == 0) {
                  iStatus2 = 0;
                  if (cbcParamCode == CbcParam::BAB) {
                    // set best solution in model as no integers
                    model_.setBestSolution(
                        model2->primalColumnSolution(), model2->numberColumns(),
                        model2->getObjValue() * model2->getObjSense());
		    // Cutoff
		    model_.solver()->setDblParam(OsiDualObjectiveLimit, 1.0e100);
                  }
                } else if (iStatus == 1) {
                  iStatus = 0;
                  iStatus2 = 1; // say infeasible
                } else if (iStatus == 2) {
                  iStatus = 0;
                  iStatus2 = 7; // say unbounded
                } else if (iStatus == 3) {
                  iStatus = 1;
                  if (iStatus2 == 9) // what does 9 mean ?????????????
                    iStatus2 = 4;
                  else
                    iStatus2 = 3; // Use nodes - as closer than solutions
                } else if (iStatus == 4) {
                  iStatus = 2; // difficulties
                  iStatus2 = 0;
                }
                model_.setProblemStatus(iStatus);
                model_.setSecondaryStatus(iStatus2);
                if ((iStatus == 2 || iStatus2 > 0) && parameters.noPrinting()) {
                  std::string statusName[] = {
                      "", "Stopped on ", "Run abandoned", "",
                      "", "User ctrl-c"};
                  std::string minor[] = {
                      "Optimal solution found",
                      "Linear relaxation infeasible",
                      "Optimal solution found (within gap tolerance)",
                      "node limit",
                      "time limit",
                      "user ctrl-c",
                      "solution limit",
                      "Linear relaxation unbounded",
                      "iterations limit",
                      "Problem proven infeasible"};
                  buffer.str("");
                  buffer << std::endl << "Result - " << statusName[iStatus].c_str()
                         << minor[iStatus2].c_str() << std::endl << std::endl;
                  buffer << "Enumerated nodes:           0" << std::endl;
                  buffer << "Total iterations:           0" << std::endl;
#if CBC_QUIET == 0
                  buffer << "Time (CPU seconds):         " << CoinCpuTime() - time0 << std::endl;
                  buffer << "Time (Wallclock Seconds):   " << CoinGetTimeOfDay() - time0Elapsed
                         << std::endl;
#endif
                  printGeneralMessage(model_, buffer.str());
                }
                // assert (lpSolver==clpSolver->getModelPtr());
                assert(clpSolver == model_.solver());
                clpSolver->setWarmStart(NULL);
                // and in babModel if exists
                if (babModel_) {
                  babModel_->setProblemStatus(iStatus);
                  babModel_->setSecondaryStatus(iStatus2);
                }
                int returnCode = callBack(&model, 1);
                if (returnCode) {
                  // exit if user wants
                  delete babModel_;
                  babModel_ = NULL;
                  return returnCode;
                }
              }
              basisHasValues = 1;
              if (dualize) {
                int returnCode =
                    static_cast<ClpSimplexOther *>(lpSolver)->restoreFromDual(
                        model2);
                if (model2->status() == 3)
                  returnCode = 0;
                delete model2;
                if (returnCode && dualize != 2)
                  lpSolver->primal(1);
                model2 = lpSolver;
              }
              if (statusUserFunction_[0]) {
                double value = model2->getObjValue();
                char buf[300];
                int pos = 0;
                int iStat = model2->status();
                if (iStat == 0) {
                  pos += sprintf(buf + pos, "optimal,");
                } else if (iStat == 1) {
                  // infeasible
                  pos += sprintf(buf + pos, "infeasible,");
                } else if (iStat == 2) {
                  // unbounded
                  pos += sprintf(buf + pos, "unbounded,");
                } else if (iStat == 3) {
                  pos += sprintf(buf + pos, "stopped on iterations or time,");
                } else if (iStat == 4) {
                  iStat = 7;
                  pos += sprintf(buf + pos, "stopped on difficulties,");
                } else if (iStat == 5) {
                  iStat = 3;
                  pos += sprintf(buf + pos, "stopped on ctrl-c,");
                } else if (iStat == 6) {
                  // bab infeasible
                  pos += sprintf(buf + pos, "integer infeasible,");
                  iStat = 1;
                } else {
                  pos += sprintf(buf + pos, "status unknown,");
                  iStat = 6;
                }
                if (info){
                   info->problemStatus = iStat;
                   info->objValue = value;
                   pos += sprintf(buf + pos, " objective %.*g", ampl_obj_prec(),
                                  value);
                   sprintf(buf + pos, "\n%d iterations",
                           model2->getIterationCount());
                   free(info->primalSolution);
                   int numberColumns = model2->numberColumns();
                   info->primalSolution = reinterpret_cast<double *>(
                       malloc(numberColumns * sizeof(double)));
                   CoinCopyN(model2->primalColumnSolution(), numberColumns,
                             info->primalSolution);
                   int numberRows = model2->numberRows();
                   free(info->dualSolution);
                   info->dualSolution = reinterpret_cast<double *>(
                      malloc(numberRows * sizeof(double)));
                   CoinCopyN(model2->dualRowSolution(), numberRows,
                             info->dualSolution);
                   CoinWarmStartBasis *basis = model2->getBasis();
                   free(info->rowStatus);
                   info->rowStatus =
                      reinterpret_cast<int *>(malloc(numberRows * sizeof(int)));
                   free(info->columnStatus);
                   info->columnStatus = reinterpret_cast<int *>(
                      malloc(numberColumns * sizeof(int)));
                   // Put basis in
                   int i;
                   // free,basic,ub,lb are 0,1,2,3
                   for (i = 0; i < numberRows; i++) {
                      CoinWarmStartBasis::Status status = basis->getArtifStatus(i);
                      info->rowStatus[i] = status;
                   }
                   for (i = 0; i < numberColumns; i++) {
                      CoinWarmStartBasis::Status status = basis->getStructStatus(i);
                      info->columnStatus[i] = status;
                   }
                   // put buffer into info
                   strcpy(info->buffer, buf);
                   delete basis;
                }
              }
            } break;
          case ClpParam::TIGHTEN:{
            if (!goodModel) {
              printGeneralMessage(model_, "** Current model not valid\n");
              continue;
            }
            int numberInfeasibilities = lpSolver->tightenPrimalBounds();
            if (numberInfeasibilities) {
               printGeneralMessage(model_,
                                   "** Analysis indicates model infeasible");
            }
          } break;
          case ClpParam::PLUSMINUS:{
            if (goodModel) {
               printGeneralWarning(model_, "** Current model not valid\n");
               continue;
            }
            ClpMatrixBase *saveMatrix = lpSolver->clpMatrix();
            ClpPackedMatrix *clpMatrix =
               dynamic_cast<ClpPackedMatrix *>(saveMatrix);
            if (clpMatrix) {
               ClpPlusMinusOneMatrix *newMatrix =
                  new ClpPlusMinusOneMatrix(*(clpMatrix->matrix()));
               if (newMatrix->getIndices()) {
                  lpSolver->replaceMatrix(newMatrix);
                  delete saveMatrix;
                  printGeneralMessage(model_,
                                      "Matrix converted to +- one matrix");
               } else {
                  printGeneralMessage(model_,
                                "Matrix can not be converted to +- 1 matrix");
               }
            } else {
               printGeneralMessage(model_, "Matrix not a ClpPackedMatrix");
            }
          }break;
          case ClpParam::OUTDUPROWS:
            dominatedCuts = true;
#ifdef JJF_ZERO
            if (goodModel) {
              int numberRows = clpSolver->getNumRows();
              // int nOut = outDupRow(clpSolver);
              CglDuplicateRow dupcuts(clpSolver);
              storedCuts = dupcuts.outDuplicates(clpSolver) != 0;
              int nOut = numberRows - clpSolver->getNumRows();
              if (nOut){
                 buffer.str("");
                 buffer << nOut << " rows eliminated";
                 printGeneralMessage(model_, buffer.str());
              }
            } else {
              printGeneralMessage(model_, "** Current model not valid\n");
            }
#endif
            break;
            case ClpParam::NETWORK:{
            if (!goodModel) {
               printGeneralWarning(model_, "** Current model not valid\n");
               continue;
            }
            ClpMatrixBase *saveMatrix = lpSolver->clpMatrix();
            ClpPackedMatrix *clpMatrix =
               dynamic_cast<ClpPackedMatrix *>(saveMatrix);
            if (clpMatrix) {
               ClpNetworkMatrix *newMatrix =
                  new ClpNetworkMatrix(*(clpMatrix->matrix()));
               if (newMatrix->getIndices()) {
                  lpSolver->replaceMatrix(newMatrix);
                  delete saveMatrix;
                  printGeneralMessage(model_,
                                      "Matrix converted to network matrix");
               } else {
                  printGeneralMessage(model_,
                          "Matrix can not be converted to network matrix");
               }
            } else {
               printGeneralMessage(model_, "Matrix not a ClpPackedMatrix");
            }
            }break;
          case ClpParam::BASISIN:{
            if (!goodModel){
               printGeneralWarning(model_, "** Current model not valid\n");
               continue;
            }
            clpParam->readValue(inputQueue, fileName, &message);
            CoinParamUtils::processFile(fileName,
                                parameters[CbcParam::DIRECTORY]->dirName(),
                                        &canOpen);
            if (!canOpen) {
               buffer.str("");
               buffer << "Unable to open file " << fileName.c_str();
               printGeneralMessage(model_, buffer.str());
               continue;
            }
            if (fileName == ""){
               fileName = clpParameters[ClpParam::BASISFILE]->fileName();
            }else{
               clpParameters[ClpParam::BASISFILE]->setFileName(fileName);
            }
#ifndef CBC_OTHER_SOLVER
            int values = lpSolver->readBasis(fileName.c_str());
            if (values == 0)
               basisHasValues = -1;
            else
               basisHasValues = 1;
            assert(lpSolver == clpSolver->getModelPtr());
            clpSolver->setWarmStart(NULL);
#endif
          } break;
          case ClpParam::BASISOUT:{
            if (!goodModel){
              printGeneralMessage(model_, "** Current model not valid\n");
              continue;
            }
            clpParam->readValue(inputQueue, fileName, &message);
            CoinParamUtils::processFile(fileName,
                                parameters[CbcParam::DIRECTORY]->dirName(),
                                        &canOpen);
            if (!canOpen) {
               buffer.str("");
               buffer << "Unable to open file " << fileName.c_str();
               printGeneralMessage(model_, buffer.str());
               continue;
            }
            if (fileName == ""){
               fileName = clpParameters[ClpParam::BASISFILE]->fileName();
            }else{
               clpParameters[ClpParam::BASISFILE]->setFileName(fileName);
            }
            ClpSimplex *model2 = lpSolver;
            model2->writeBasis(fileName.c_str(), outputFormat > 1,
                               outputFormat - 2);
            time2 = CoinCpuTime();
            totalTime += time2 - time1;
            time1 = time2;
          } break;
          case ClpParam::ALLSLACK:
            lpSolver->allSlackBasis(true);
            break;
          case ClpParam::NETLIB_DUAL:
          case ClpParam::NETLIB_EITHER:
          case ClpParam::NETLIB_BARRIER:
          case ClpParam::NETLIB_PRIMAL:
          case ClpParam::NETLIB_TUNE: {
            printf("unit test is now only from clp - does same thing\n");
            // return(22);
          } break;
          case ClpParam::FAKEBOUND:{
            if (!goodModel) {
               continue;
            }
            // get bound
            if ((status = clpParam->readValue(inputQueue, dValue, &message))){
               std::cout << "Must enter value for " << clpParam->name()
                         << std::endl;
               continue;
            }
            buffer.str("");
            buffer << "Setting " << clpParam->name() << " to DEBUG " << dValue
                   << std::endl;
            printGeneralMessage(model_, buffer.str());
            int iRow;
            int numberRows = lpSolver->numberRows();
            double *rowLower = lpSolver->rowLower();
            double *rowUpper = lpSolver->rowUpper();
            for (iRow = 0; iRow < numberRows; iRow++) {
               // leave free ones for now
               if (rowLower[iRow] > -1.0e20 || rowUpper[iRow] < 1.0e20) {
                  rowLower[iRow] = std::max(rowLower[iRow], -dValue);
                  rowUpper[iRow] = std::min(rowUpper[iRow], dValue);
               }
            }
            int iColumn;
            int numberColumns = lpSolver->numberColumns();
            double *columnLower = lpSolver->columnLower();
            double *columnUpper = lpSolver->columnUpper();
            for (iColumn = 0; iColumn < numberColumns; iColumn++) {
               // leave free ones for now
               if (columnLower[iColumn] > -1.0e20 ||
                   columnUpper[iColumn] < 1.0e20) {
                  columnLower[iColumn] =
                     std::max(columnLower[iColumn], -dValue);
                  columnUpper[iColumn] = std::min(columnUpper[iColumn], dValue);
               }
            }
          }break;
          case ClpParam::REALLY_SCALE:
            if (goodModel) {
              ClpSimplex newModel(*lpSolver, lpSolver->scalingFlag());
              printf("model really really scaled\n");
              *lpSolver = newModel;
            }
            break;
          case ClpParam::USERCLP:
#ifdef USER_HAS_FAKE_CLP
            // Replace the sample code by whatever you want
            if (goodModel) {
              // Way of using an existing piece of code
              OsiClpSolverInterface *clpSolver =
                  dynamic_cast<OsiClpSolverInterface *>(model_.solver());
              ClpSimplex *lpSolver = clpSolver->getModelPtr();
              // set time from integer model
              double timeToGo = model_.getMaximumSeconds();
              lpSolver->setMaximumSeconds(timeToGo);
              int extra1 = parameters[CbcParam::EXTRA1]->intVal();
              fakeMain2(*lpSolver, *clpSolver, extra1);
              lpSolver = clpSolver->getModelPtr();
              // My actual usage has objective only in clpSolver
              // double objectiveValue=clpSolver->getObjValue();
              // int iStat = lpSolver->status();
              // int iStat2 = lpSolver->secondaryStatus();
            }
#endif
            break;
            case ClpParam::PARAMETRICS:{
            if (!goodModel){
              printGeneralMessage(model_, "** Current model not valid\n");
            }
            cbcParam->readValue(inputQueue, fileName, &message);
            CoinParamUtils::processFile(fileName,
                                parameters[CbcParam::DIRECTORY]->dirName(),
                                        &canOpen);
            if (!canOpen) {
               buffer.str("");
               buffer << "Unable to open file " << fileName.c_str();
               printGeneralMessage(model_, buffer.str());
               continue;
            }
            if (fileName == ""){
               fileName = parameters[CbcParam::EXPORTFILE]->fileName();
            }else{
               parameters[CbcParam::EXPORTFILE]->setFileName(fileName);
            }
            static_cast<ClpSimplexOther *>(lpSolver)->parametrics(
                  fileName.c_str());
            time2 = CoinCpuTime();
            totalTime += time2 - time1;
            time1 = time2;
            } break;
          case ClpParam::GUESS:
            if (goodModel && model_.solver()) {
              OsiClpSolverInterface *clpSolver =
                  dynamic_cast<OsiClpSolverInterface *>(model_.solver());
              assert(clpSolver);
              lpSolver = clpSolver->getModelPtr();
              assert(lpSolver);
              ClpSimplexOther *model2 =
                  static_cast<ClpSimplexOther *>(lpSolver);
              //This is an ugly hack to make this compile, but we'll be getting rid of
              //this soon anyway.
              std::string input = model2->guess(0);
              if (input != ""){
                 while (!inputQueue.empty()){
                    inputQueue.pop_front();
                 }
                 std::istringstream inputStream(input);
                 CoinParamUtils::readFromStream(inputQueue, inputStream);
              }else{
                 printGeneralWarning(model,
                                     "** Guess unable to generate commands\n");
              }
            } else {
               printGeneralWarning(model, "** Guess needs a valid model\n");
            }
            break;
          default:
            // abort();
            break;
          }

          switch (cbcParamCode) {
           case CbcParam::SOLVECONTINUOUS:
             // Solve problem dual simplex
             inputQueue.push_front("-dualSimplex");
             continue;
           case CbcParam::STATISTICS:{
            if (!goodModel){
              printGeneralMessage(model_, "** Current model not valid\n");
              continue;
            }
            // If presolve on look at presolved
            bool deleteModel2 = false;
            ClpSimplex *model2 = lpSolver;
            if (preSolve) {
               ClpPresolve pinfo;
               int presolveOptions2 = presolveOptions & ~0x40000000;
               if ((presolveOptions2 & 0xffff) != 0)
                  pinfo.setPresolveActions(presolveOptions2);
               pinfo.setSubstitution(substitution);
               if ((printOptions & 1) != 0)
                  pinfo.statistics();
               double presolveTolerance =
                  clpParameters[ClpParam::PRESOLVETOLERANCE]->dblVal();
               model2 = pinfo.presolvedModel(*lpSolver, presolveTolerance,
                                             true, preSolve);
               if (model2) {
                  printf("Statistics for presolved model\n");
                  deleteModel2 = true;
               } else {
                  printf("Presolved model looks infeasible - will use "
                         "unpresolved\n");
                  model2 = lpSolver;
               }
            } else {
               printf("Statistics for unpresolved model\n");
               model2 = lpSolver;
            }
            statistics(lpSolver, model2);
            if (deleteModel2)
               delete model2;
           } break;
          case CbcParam::MIPLIB:
            // User can set options - main difference is lack of model and
            // CglPreProcess
            goodModel = true;
            parameters[CbcParam::MULTIPLEROOTS]->setVal(0);
            /*
                          Run branch-and-cut. First set a few options -- node
               comparison, scaling. Print elapsed time at the end.
                        */
          case CbcParam::BAB: {
#ifdef CBC_CLUMSY_CODING
	      parameters.setModel(&model_);
	      parameters.setGoodModel(true);
	      parameters.synchronizeModel();
#endif
              bool miplib = cbcParamCode == CbcParam::MIPLIB;
              int logLevel = parameters[CbcParam::LPLOGLEVEL]->intVal();
              int truncateColumns = COIN_INT_MAX;
              int truncateRows = -1;
              bool redoSOS = false;
              double *truncatedRhsLower = NULL;
              double *truncatedRhsUpper = NULL;
              int *newPriorities = NULL;
	      if (model_.solver()->getObjSense()==-1.0) {
		// If cutoff set flip
		if (model_.getCutoff()<1.0e30) {
		  double cutoff = model_.getCutoff();
		  model_.setCutoff(-cutoff);
		}
	      }
              // Reduce printout
              if (logLevel <= 1) {
                model_.solver()->setHintParam(OsiDoReducePrint, true,
                                              OsiHintTry);
              } else {
                model_.solver()->setHintParam(OsiDoReducePrint, false,
                                              OsiHintTry);
              }
              {
                OsiSolverInterface *solver = model_.solver();
#ifndef CBC_OTHER_SOLVER
                OsiClpSolverInterface *si =
                    dynamic_cast<OsiClpSolverInterface *>(solver);
                assert(si != NULL);
                si->getModelPtr()->scaling(doScaling);
                ClpSimplex *lpSolver = si->getModelPtr();
                // deal with positive edge
                double psi = clpParameters[ClpParam::PSI]->dblVal();
                if (psi > 0.0) {
                  ClpDualRowPivot *dualp = lpSolver->dualRowPivot();
                  ClpDualRowSteepest *d1 =
                      dynamic_cast<ClpDualRowSteepest *>(dualp);
                  ClpDualRowDantzig *d2 =
                      dynamic_cast<ClpDualRowDantzig *>(dualp);
                  if (d1) {
                    ClpPEDualRowSteepest p(psi, d1->mode());
                    lpSolver->setDualRowPivotAlgorithm(p);
                  } else if (d2) {
                    ClpPEDualRowDantzig p(psi);
                    lpSolver->setDualRowPivotAlgorithm(p);
                  }
                  ClpPrimalColumnPivot *primalp = lpSolver->primalColumnPivot();
                  ClpPrimalColumnSteepest *p1 =
                      dynamic_cast<ClpPrimalColumnSteepest *>(primalp);
                  ClpPrimalColumnDantzig *p2 =
                      dynamic_cast<ClpPrimalColumnDantzig *>(primalp);
                  if (p1) {
                    ClpPEPrimalColumnSteepest p(psi, p1->mode());
                    lpSolver->setPrimalColumnPivotAlgorithm(p);
                  } else if (p2) {
                    ClpPEPrimalColumnDantzig p(psi);
                    lpSolver->setPrimalColumnPivotAlgorithm(p);
                  }
                }
                if (doVector) {
                  ClpMatrixBase *matrix = lpSolver->clpMatrix();
                  if (dynamic_cast<ClpPackedMatrix *>(matrix)) {
                    ClpPackedMatrix *clpMatrix =
                        dynamic_cast<ClpPackedMatrix *>(matrix);
                    clpMatrix->makeSpecialColumnCopy();
                  }
                }
#elif CBC_OTHER_SOLVER == 1
                OsiCpxSolverInterface *si =
                    dynamic_cast<OsiCpxSolverInterface *>(solver);
                assert(si != NULL);
#endif
                statistics_nrows = si->getNumRows();
                statistics_ncols = si->getNumCols();
                statistics_nprocessedrows = si->getNumRows();
                statistics_nprocessedcols = si->getNumCols();
                // See if quadratic
#ifndef CBC_OTHER_SOLVER
#ifdef COIN_HAS_LINK
                if (!complicatedInteger) {
                  ClpQuadraticObjective *obj =
                      (dynamic_cast<ClpQuadraticObjective *>(
                          lpSolver->objectiveAsObject()));
                  if (obj) {
                    preProcess = 0;
                    int testOsiOptions = parameters[CbcParam::TESTOSI]->intVal();
                    parameters[CbcParam::TESTOSI]->setVal(std::max(0, testOsiOptions));
                    // create coin model
                    coinModel = lpSolver->createCoinModel();
                    assert(coinModel);
                    // load from coin model
                    OsiSolverLink solver1;
                    OsiSolverInterface *solver2 = solver1.clone();
                    model_.assignSolver(solver2, false);
                    OsiSolverLink *si =
                        dynamic_cast<OsiSolverLink *>(model_.solver());
                    assert(si != NULL);
                    si->setDefaultMeshSize(0.001);
                    // need some relative granularity
                    si->setDefaultBound(100.0);
                    double dextra3 = parameters[CbcParam::DEXTRA3]->dblVal();
                    if (dextra3)
                      si->setDefaultMeshSize(dextra3);
                    si->setDefaultBound(1000.0);
                    si->setIntegerPriority(1000);
                    si->setBiLinearPriority(10000);
                    biLinearProblem = true;
                    si->setSpecialOptions2(2 + 4 + 8);
                    CoinModel *model2 = coinModel;
                    si->load(*model2, true, parameters[CbcParam::LPLOGLEVEL]->intVal());
                    // redo
                    solver = model_.solver();
                    clpSolver = dynamic_cast<OsiClpSolverInterface *>(solver);
                    lpSolver = clpSolver->getModelPtr();
                    clpSolver->messageHandler()->setLogLevel(0);
                    testOsiParameters = 0;
                    complicatedInteger = 2; // allow cuts
                    OsiSolverInterface *coinSolver = model_.solver();
                    OsiSolverLink *linkSolver =
                        dynamic_cast<OsiSolverLink *>(coinSolver);
                    if (linkSolver->quadraticModel()) {
                      ClpSimplex *qp = linkSolver->quadraticModel();
                      // linkSolver->nonlinearSLP(std::max(slpValue,10),1.0e-5);
                      qp->nonlinearSLP(std::max(slpValue, 40), 1.0e-5);
                      qp->primal(1);
                      OsiSolverLinearizedQuadratic solver2(qp);
                      const double *solution = NULL;
                      // Reduce printout
                      solver2.setHintParam(OsiDoReducePrint, true, OsiHintTry);
                      CbcModel model2(solver2);
                      // Now do requested saves and modifications
                      CbcModel *cbcModel = &model2;
                      OsiSolverInterface *osiModel = model2.solver();
                      OsiClpSolverInterface *osiclpModel =
                          dynamic_cast<OsiClpSolverInterface *>(osiModel);
                      ClpSimplex *clpModel = osiclpModel->getModelPtr();

                      // Set changed values
                      int numCutGens = 0;

                      CglProbing probing;
                      probing.setMaxProbe(10);
                      probing.setMaxLook(10);
                      probing.setMaxElements(200);
                      probing.setMaxProbeRoot(50);
                      probing.setMaxLookRoot(10);
                      probing.setRowCuts(3);
                      probing.setUsingObjective(true);
                      cbcModel->addCutGenerator(&probing, -1, "Probing", true,
                                                false, false, -100, -1, -1);
                      cbcModel->cutGenerator(numCutGens++)->setTiming(true);

                      CglGomory gomory;
                      gomory.setLimitAtRoot(512);
                      cbcModel->addCutGenerator(&gomory, -98, "Gomory", true,
                                                false, false, -100, -1, -1);
                      cbcModel->cutGenerator(numCutGens++)->setTiming(true);

                      CglKnapsackCover knapsackCover;
                      cbcModel->addCutGenerator(&knapsackCover, -98,
                                                "KnapsackCover", true, false,
                                                false, -100, -1, -1);
                      cbcModel->cutGenerator(numCutGens++)->setTiming(true);

                      CglRedSplit redSplit;
                      cbcModel->addCutGenerator(&redSplit, -99, "RedSplit",
                                                true, false, false, -100, -1,
                                                -1);
                      cbcModel->cutGenerator(numCutGens++)->setTiming(true);

                      CglBKClique bkClique;
                      bkClique.setMaxCallsBK(1000);
                      bkClique.setExtendingMethod(4);
                      bkClique.setPivotingStrategy(CoinBronKerbosch::PivotingStrategy::Weight);
                      cbcModel->addCutGenerator(&bkClique, -98, "Clique", true,
                                                false, false, -100, -1, -1);
                      cbcModel->cutGenerator(numCutGens++)->setTiming(true);

                      CglOddWheel oddWheel;
                      oddWheel.setExtendingMethod(2);
                      cbcModel->addCutGenerator(&oddWheel, -98, "OddWheel",
                                                true, false, false, -100, -1,
                                                -1);
                      cbcModel->cutGenerator(numCutGens++)->setTiming(true);

                      CglMixedIntegerRounding2 mixedIntegerRounding2;
                      cbcModel->addCutGenerator(&mixedIntegerRounding2, -98,
                                                "MixedIntegerRounding2", true,
                                                false, false, -100, -1, -1);
                      cbcModel->cutGenerator(numCutGens++)->setTiming(true);

                      CglFlowCover flowCover;
                      cbcModel->addCutGenerator(&flowCover, -98, "FlowCover",
                                                true, false, false, -100, -1,
                                                -1);
                      cbcModel->cutGenerator(numCutGens++)->setTiming(true);

                      CglTwomir twomir;
                      twomir.setMaxElements(250);
                      cbcModel->addCutGenerator(&twomir, -99, "Twomir", true,
                                                false, false, -100, -1, -1);
                      cbcModel->cutGenerator(numCutGens++)->setTiming(true);
                      int heuristicOption =
                         parameters[CbcParam::HEURISTICSTRATEGY]->modeVal();
                      if (heuristicOption) {
                        CbcHeuristicFPump heuristicFPump(*cbcModel);
                        heuristicFPump.setWhen(13);
                        heuristicFPump.setMaximumPasses(20);
                        heuristicFPump.setMaximumRetries(7);
                        heuristicFPump.setHeuristicName("feasibility pump");
                        heuristicFPump.setInitialWeight(1);
                        heuristicFPump.setFractionSmall(0.6);
                        cbcModel->addHeuristic(&heuristicFPump);

                        CbcRounding rounding(*cbcModel);
                        rounding.setHeuristicName("rounding");
                        cbcModel->addHeuristic(&rounding);

                        CbcHeuristicLocal heuristicLocal(*cbcModel);
                        heuristicLocal.setHeuristicName("combine solutions");
                        heuristicLocal.setSearchType(1);
                        heuristicLocal.setFractionSmall(0.6);
                        cbcModel->addHeuristic(&heuristicLocal);

                        CbcHeuristicGreedyCover heuristicGreedyCover(*cbcModel);
                        heuristicGreedyCover.setHeuristicName("greedy cover");
                        cbcModel->addHeuristic(&heuristicGreedyCover);

                        CbcHeuristicGreedyEquality heuristicGreedyEquality(
                            *cbcModel);
                        heuristicGreedyEquality.setHeuristicName(
                            "greedy equality");
                        cbcModel->addHeuristic(&heuristicGreedyEquality);
#ifdef CBC_EXPERIMENT_JJF
#ifndef CBC_OTHER_SOLVER
                        CbcHeuristicRandRound heuristicRandRound(*cbcModel);
                        heuristicRandRound.setHeuristicName("random rounding");
                        cbcModel->addHeuristic(&heuristicRandRound);
#endif
#endif
                      }
                      CbcCompareDefault compare;
                      cbcModel->setNodeComparison(compare);
                      cbcModel->setNumberBeforeTrust(10);
                      cbcModel->setSpecialOptions(2);
                      cbcModel->messageHandler()->setLogLevel(1);
                      cbcModel->setMaximumCutPassesAtRoot(-100);
                      cbcModel->setMaximumCutPasses(1);
                      cbcModel->setMinimumDrop(0.05);
                      // For branchAndBound this may help
                      clpModel->defaultFactorizationFrequency();
                      clpModel->setDualBound(1.0001e+08);
                      clpModel->setPerturbation(50);
                      osiclpModel->setSpecialOptions(193);
                      osiclpModel->messageHandler()->setLogLevel(0);
                      osiclpModel->setIntParam(OsiMaxNumIterationHotStart, 100);
                      osiclpModel->setHintParam(OsiDoReducePrint, true,
                                                OsiHintTry);
                      // You can save some time by switching off message
                      // building
                      // clpModel->messagesPointer()->setDetailMessages(100,10000,(int
                      // *) NULL);

                      // Solve

                      cbcModel->initialSolve();
                      if (clpModel->tightenPrimalBounds() != 0) {
                        buffer.str("");
                        buffer << "Problem is infeasible - tightenPrimalBounds!";
                        printGeneralMessage(model_, buffer.str());
                        break;
                      }
                      clpModel->dual(); // clean up
                      cbcModel->initialSolve();
#ifdef CBC_THREAD
                      int numberThreads = parameters[CbcParam::THREADS]->intVal();
                      cbcModel->setNumberThreads(numberThreads % 100);
                      cbcModel->setThreadMode(std::min((numberThreads%1000) / 100, 7));
#endif
                      // setCutAndHeuristicOptions(*cbcModel);
                      cbcModel->branchAndBound();
                      OsiSolverLinearizedQuadratic *solver3 =
                          dynamic_cast<OsiSolverLinearizedQuadratic *>(
                              model2.solver());
                      assert(solver3);
                      solution = solver3->bestSolution();
                      double bestObjectiveValue = solver3->bestObjectiveValue();
                      linkSolver->setBestObjectiveValue(bestObjectiveValue);
                      if (solution) {
			// make sure array large enough
			int arraySize = solver3->getNumCols();
			int wanted = linkSolver->getNumCols();
			double * solution2 = new double[wanted];
			memset(solution2,0,wanted*sizeof(double));
			memcpy(solution2,solution,arraySize*sizeof(double));
                        linkSolver->setBestSolution(solution2, solver3->getNumCols());
			model_.setBestSolution(solution2,model_.getNumCols(),
					       bestObjectiveValue);
			model_.setCutoff(bestObjectiveValue+1.0e-4);
			delete [] solution2;
                      }
                      CbcHeuristicDynamic3 dynamic(model_);
                      dynamic.setHeuristicName("dynamic pass thru");
                      if (heuristicOption)
                        model_.addHeuristic(&dynamic);
                      // if convex
                      if ((linkSolver->specialOptions2() & 4) != 0 &&
                          solution) {
                        int numberColumns = coinModel->numberColumns();
                        assert(linkSolver->objectiveVariable() ==
                               numberColumns);
                        // add OA cut
                        double offset;
                        double *gradient = new double[numberColumns + 1];
                        memcpy(gradient,
                               qp->objectiveAsObject()->gradient(
                                   qp, solution, offset, true, 2),
                               numberColumns * sizeof(double));
                        double rhs = 0.0;
                        int *column = new int[numberColumns + 1];
                        int n = 0;
                        for (int i = 0; i < numberColumns; i++) {
                          double value = gradient[i];
                          if (fabs(value) > 1.0e-12) {
                            gradient[n] = value;
                            rhs += value * solution[i];
                            column[n++] = i;
                          }
                        }
                        gradient[n] = -1.0;
                        column[n++] = numberColumns;
                        storedAmpl.addCut(-COIN_DBL_MAX, offset + 1.0e-7, n,
                                          column, gradient);
                        linkSolver->addRow(n, column, gradient, -COIN_DBL_MAX,
                                           offset + 1.0e-7);
                        delete[] gradient;
                        delete[] column;
                      }
                      // could do three way branching round a) continuous b)
                      // best solution
                      // printf("obj %g\n", bestObjectiveValue);
                      linkSolver->initialSolve();
                    }
                  }
                }
#endif
#endif
                if (logLevel <= 1)
                  si->setHintParam(OsiDoReducePrint, true, OsiHintTry);
#ifndef CBC_OTHER_SOLVER
                si->setSpecialOptions(0x40000000);
#endif
              }
              if (!miplib) {
                if (!preSolve) {
                  model_.solver()->setHintParam(OsiDoPresolveInInitial, false,
                                                OsiHintTry);
                  model_.solver()->setHintParam(OsiDoPresolveInResolve, false,
                                                OsiHintTry);
                }
                double time1a = CoinCpuTime();
                OsiSolverInterface *solver = model_.solver();
#ifndef CBC_OTHER_SOLVER
                OsiClpSolverInterface *si =
                    dynamic_cast<OsiClpSolverInterface *>(solver);
                if (si)
                  si->setSpecialOptions(si->specialOptions() | 1024);
#endif
                model_.initialSolve();
#ifndef CBC_OTHER_SOLVER
                ClpSimplex *clpSolver = si->getModelPtr();
                int iStatus = clpSolver->status();
                int iStatus2 = clpSolver->secondaryStatus();
                if (iStatus == 0) {
                  iStatus2 = 0;
                } else if (iStatus == 1) {
                  iStatus = 0;
                  iStatus2 = 1; // say infeasible
                } else if (iStatus == 2) {
                  iStatus = 0;
                  iStatus2 = 7; // say unbounded
                } else if (iStatus == 3) {
                  iStatus = 1;
                  if (iStatus2 == 9)
                    iStatus2 = 4;
                  else
                    iStatus2 = 3; // Use nodes - as closer than solutions
                } else if (iStatus == 4) {
                  iStatus = 2; // difficulties
                  iStatus2 = 0;
                }
                model_.setProblemStatus(iStatus);
                model_.setSecondaryStatus(iStatus2);
                si->setWarmStart(NULL);
                int returnCode = 0;
                if (callBack != NULL)
                  returnCode = callBack(&model_, 1);
                if (returnCode) {
                  // exit if user wants
                  delete babModel_;
                  babModel_ = NULL;
                  return returnCode;
                }
                if (clpSolver->status() > 0) {
                  // and in babModel if exists
                  if (babModel_) {
                    babModel_->setProblemStatus(iStatus);
                    babModel_->setSecondaryStatus(iStatus2);
                  }
                  iStatus = clpSolver->status();
                  const char *msg[] = {"infeasible", "unbounded", "stopped",
                                       "difficulties", "other"};
                  buffer.str("");
                  buffer << "Problem is " << msg[iStatus - 1] << " - " <<  CoinCpuTime() - time1a
                         << " seconds";
                  printGeneralMessage(model_, buffer.str());
                  break;
                }
                clpSolver->setSpecialOptions(
                    clpSolver->specialOptions() |
                    IN_BRANCH_AND_BOUND); // say is Cbc (and in branch and
                                          // bound)
#elif CBC_OTHER_SOLVER == 1
#endif
                buffer.str("");
                buffer << "Continuous objective value is " << solver->getObjValue() << " - "
                       << CoinCpuTime() - time1a << " seconds";
                printGeneralMessage(model_, buffer.str());
                if (model_.getMaximumNodes() == -987654321) {
                  // See if No objective!
                  int numberColumns = clpSolver->getNumCols();
                  const double *obj = clpSolver->getObjCoefficients();
                  const double *lower = clpSolver->getColLower();
                  const double *upper = clpSolver->getColUpper();
                  int nObj = 0;
                  for (int i = 0; i < numberColumns; i++) {
                    if (upper[i] > lower[i] && obj[i])
                      nObj++;
                  }
                  if (!nObj) {
                    printf("************No objective!!\n");
                    model_.setMaximumSolutions(1);
                    // Column copy
                    CoinPackedMatrix matrixByCol(
                        *model_.solver()->getMatrixByCol());
                    // const double * element = matrixByCol.getElements();
                    // const int * row = matrixByCol.getIndices();
                    // const CoinBigIndex * columnStart =
                    // matrixByCol.getVectorStarts();
                    const int *columnLength = matrixByCol.getVectorLengths();
                    for (int i = 0; i < numberColumns; i++) {
                      double value = (CoinDrand48() + 0.5) * 10000;
                      value = 10;
                      value *= columnLength[i];
                      int iValue = static_cast<int>(value) / 10;
                      // iValue=1;
                      clpSolver->setObjCoeff(i, iValue);
                    }
                  }
                }
#ifndef CBC_OTHER_SOLVER
                if (!complicatedInteger && preProcess == 0 &&
                    clpSolver->tightenPrimalBounds(0.0, 0, true) != 0) {
                  printGeneralMessage(model_,
                              "Problem is infeasible - tightenPrimalBounds!");
                  model_.setProblemStatus(0);
                  model_.setSecondaryStatus(1);
                  // say infeasible for solution
                  integerStatus = 6;
                  // and in babModel if exists
                  if (babModel_) {
                    babModel_->setProblemStatus(0);
                    babModel_->setSecondaryStatus(1);
                  }
                  break;
                }
                if (clpSolver->dualBound() == 1.0e10) {
                  ClpSimplex temp = *clpSolver;
                  temp.setLogLevel(0);
                  temp.dual(0, 7);
                  // user did not set - so modify
                  // get largest scaled away from bound
                  double largest = 1.0e-12;
                  double largestScaled = 1.0e-12;
                  int numberRows = temp.numberRows();
                  const double *rowPrimal = temp.primalRowSolution();
                  const double *rowLower = temp.rowLower();
                  const double *rowUpper = temp.rowUpper();
                  const double *rowScale = temp.rowScale();
                  int iRow;
                  for (iRow = 0; iRow < numberRows; iRow++) {
                    double value = rowPrimal[iRow];
                    double above = value - rowLower[iRow];
                    double below = rowUpper[iRow] - value;
                    if (above < 1.0e12) {
                      largest = std::max(largest, above);
                    }
                    if (below < 1.0e12) {
                      largest = std::max(largest, below);
                    }
                    if (rowScale) {
                      double multiplier = rowScale[iRow];
                      above *= multiplier;
                      below *= multiplier;
                    }
                    if (above < 1.0e12) {
                      largestScaled = std::max(largestScaled, above);
                    }
                    if (below < 1.0e12) {
                      largestScaled = std::max(largestScaled, below);
                    }
                  }

                  int numberColumns = temp.numberColumns();
                  const double *columnPrimal = temp.primalColumnSolution();
                  const double *columnLower = temp.columnLower();
                  const double *columnUpper = temp.columnUpper();
                  const double *columnScale = temp.columnScale();
                  int iColumn;
                  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                    double value = columnPrimal[iColumn];
                    double above = value - columnLower[iColumn];
                    double below = columnUpper[iColumn] - value;
                    if (above < 1.0e12) {
                      largest = std::max(largest, above);
                    }
                    if (below < 1.0e12) {
                      largest = std::max(largest, below);
                    }
                    if (columnScale) {
                      double multiplier = 1.0 / columnScale[iColumn];
                      above *= multiplier;
                      below *= multiplier;
                    }
                    if (above < 1.0e12) {
                      largestScaled = std::max(largestScaled, above);
                    }
                    if (below < 1.0e12) {
                      largestScaled = std::max(largestScaled, below);
                    }
                  }
#ifdef COIN_DEVELOP
                  buffer.str("");
                  buffer << "Largest (scaled) away from bound "
                         << largestScaled << " unscaled " << largest
                         << std::endl;
                  printGeneralMessage(model_, buffer.str());
#endif
                  clpSolver->setDualBound(
                      std::max(1.0001e8, std::min(100.0 * largest, 1.00001e10)));
                }
                si->resolve(); // clean up
#endif
              }
              // If user made settings then use them
              if (!defaultSettings) {
                OsiSolverInterface *solver = model_.solver();
                if (!doScaling)
                  solver->setHintParam(OsiDoScale, false, OsiHintTry);
#ifndef CBC_OTHER_SOLVER
                OsiClpSolverInterface *si =
                    dynamic_cast<OsiClpSolverInterface *>(solver);
                assert(si != NULL);
                // get clp itself
                ClpSimplex *modelC = si->getModelPtr();
                // if (modelC->tightenPrimalBounds()!=0) {
                // std::cout<<"Problem is infeasible!"<<std::endl;
                // break;
                //}
                // bounds based on continuous
                if (tightenFactor && !complicatedInteger) {
                  if (modelC->tightenPrimalBounds(tightenFactor) != 0) {
                    printGeneralMessage(model_, "Problem is infeasible!");
                    model_.setProblemStatus(0);
                    model_.setSecondaryStatus(1);
                    // and in babModel if exists
                    if (babModel_) {
                      babModel_->setProblemStatus(0);
                      babModel_->setSecondaryStatus(1);
                    }
                    break;
                  }
                }
#endif
              }

              // See if we want preprocessing
              // FIXME: Why does code break if this is defined here instead of up top?
              //OsiSolverInterface *saveSolver = NULL;
              //CglPreProcess process;
              // Say integers in sync
              bool integersOK = true;
              delete babModel_;
              babModel_ = new CbcModel(model_);
#ifndef CBC_OTHER_SOLVER
              int numberChanged = 0;
              OsiSolverInterface *solver3 = clpSolver->clone();
              babModel_->assignSolver(solver3);
              OsiClpSolverInterface *clpSolver2 =
                  dynamic_cast<OsiClpSolverInterface *>(babModel_->solver());
              if (clpSolver2->messageHandler()->logLevel())
                clpSolver2->messageHandler()->setLogLevel(1);
              if (logLevel > -1)
                clpSolver2->messageHandler()->setLogLevel(logLevel);
              lpSolver = clpSolver2->getModelPtr();
              if (lpSolver->factorizationFrequency() == 200 && !miplib) {
                // User did not touch preset
                int numberRows = lpSolver->numberRows();
                const int cutoff1 = 10000;
                const int cutoff2 = 100000;
                const int base = 75;
                const int freq0 = 50;
                const int freq1 = 200;
                const int freq2 = 400;
                const int maximum = 1000;
                int frequency;
                if (numberRows < cutoff1)
                  frequency = base + numberRows / freq0;
                else if (numberRows < cutoff2)
                  frequency =
                      base + cutoff1 / freq0 + (numberRows - cutoff1) / freq1;
                else
                  frequency = base + cutoff1 / freq0 +
                              (cutoff2 - cutoff1) / freq1 +
                              (numberRows - cutoff2) / freq2;
                lpSolver->setFactorizationFrequency(
                    std::min(maximum, frequency));
              }
#elif CBC_OTHER_SOLVER == 1
              OsiSolverInterface *solver3 = model_.solver()->clone();
              babModel_->assignSolver(solver3);
#endif
              time2 = CoinCpuTime();
              totalTime += time2 - time1;
              // time1 = time2;
              double timeLeft = babModel_->getMaximumSeconds();
              int numberOriginalColumns = babModel_->solver()->getNumCols();
              if (preProcess == 7) {
                // use strategy instead
                preProcess = 0;
                useStrategy = true;
#ifdef COIN_HAS_LINK
                // empty out any cuts
                if (storedAmpl.sizeRowCuts()) {
                  printf(
                      "Emptying ampl stored cuts as internal preprocessing\n");
                  CglStored temp;
                  storedAmpl = temp;
                }
#endif
              }
              if (preProcess && cbcParamCode == CbcParam::BAB) {
                // see whether to switch off preprocessing
                // only allow SOS and integer
                OsiObject **objects = babModel_->objects();
                int numberObjects = babModel_->numberObjects();
                for (int iObj = 0; iObj < numberObjects; iObj++) {
                  CbcSOS *objSOS = dynamic_cast<CbcSOS *>(objects[iObj]);
                  CbcSimpleInteger *objSimpleInteger =
                      dynamic_cast<CbcSimpleInteger *>(objects[iObj]);
                  if (!objSimpleInteger && !objSOS) {
                    // find all integers anyway
                    babModel_->findIntegers(true);
                    preProcess = 0;
                    break;
                  }
                }
              }
              if (cbcParamCode == CbcParam::BAB) {
                if (preProcess == 0 && numberLotSizing) {
                  if (!babModel_->numberObjects()) {
                    /* model may not have created objects
                                       If none then create
                                    */
                    babModel_->findIntegers(true);
                  }
                  // Lotsizing
                  // int numberColumns = babModel_->solver()->getNumCols();
                  CbcObject **objects = new CbcObject *[numberLotSizing];
                  double points[] = {0.0, 0.0, 0.0, 0.0};
                  for (int i = 0; i < numberLotSizing; i++) {
                    int iColumn = lotsize[i].column;
                    points[2] = lotsize[i].low;
                    points[3] = lotsize[i].high;
                    objects[i] =
                        new CbcLotsize(&model_, iColumn, 2, points, true);
                  }
                  babModel_->addObjects(numberLotSizing, objects);
                  for (int i = 0; i < numberLotSizing; i++)
                    delete objects[i];
                  delete[] objects;
                }
                double limit;
                clpSolver->getDblParam(OsiDualObjectiveLimit, limit);
                if (clpSolver->getObjValue() * clpSolver->getObjSense() >=
                    limit * clpSolver->getObjSense())
                  preProcess = 0;
              }
              if (mipStartBefore.size()) {
                CbcModel tempModel = *babModel_;
                assert(babModel_->getNumCols() == model_.getNumCols());
                std::vector<std::string> colNames;
                for (int i = 0; (i < model_.solver()->getNumCols()); ++i)
                  colNames.push_back(model_.solver()->getColName(i));
                std::vector<double> x(model_.getNumCols(), 0.0);
                double obj;
                int status = CbcMipStartIO::computeCompleteSolution(
                    &tempModel, tempModel.solver(), colNames, mipStartBefore,
                    &x[0], obj, 0, tempModel.messageHandler(),
                    tempModel.messagesPointer());
                // set cutoff ( a trifle high)
                if (!status) {
                  double newCutoff =
                      std::min(babModel_->getCutoff(), obj + 1.0e-4);
                  babModel_->setBestSolution(&x[0], static_cast<int>(x.size()),
                                             obj, false);
                  babModel_->setCutoff(newCutoff);
                  babModel_->setSolutionCount(1);
                  model_.setBestSolution(&x[0], static_cast<int>(x.size()), obj,
                                         false);
                  model_.setCutoff(newCutoff);
                  model_.setSolutionCount(1);
                }
              }
	      // Do heuristics if asked for
	      //HEURSTART
 	      // Do heuristics if asked for
	      if (parameters[CbcParam::DOHEURISTIC]->modeVal()) {
#if CBC_USE_INITIAL_TIME == 1
		if (model_.useElapsedTime())
		  model_.setDblParam(CbcModel::CbcStartSeconds,
				     CoinGetTimeOfDay());
		else
		  model_.setDblParam(CbcModel::CbcStartSeconds, CoinCpuTime());
#endif
		int vubMode = parameters[CbcParam::VUBTRY]->intVal();
		if (vubMode != -1) {
		  // look at vubs
		  // extra1 is number of ints to leave free
		  // Just ones which affect >= extra3
		  int extra3 = parameters[CbcParam::EXTRA3]->intVal();
		  /* 2 is cost above which to fix if feasible
		     3 is fraction of integer variables fixed if
		     relaxing (0.97) 4 is fraction of all variables fixed if
		     relaxing (0.0)
		  */
		  double dextra[6];
		  int extra[5];
		  extra[1] = parameters[CbcParam::EXTRA1]->intVal();
		  int exp1 = parameters[CbcParam::EXPERIMENT]->intVal();
		  if (exp1 == 4 && extra[1] == -1)
		    extra[1] = 999998;
		  dextra[1] = parameters[CbcParam::FAKEINCREMENT]->dblVal();
		  dextra[2] = parameters[CbcParam::FAKECUTOFF]->dblVal();
		  dextra[3] = parameters[CbcParam::DEXTRA3]->dblVal();
		  dextra[4] = parameters[CbcParam::DEXTRA4]->dblVal();
		  dextra[5] = parameters[CbcParam::DEXTRA5]->dblVal();
		  if (!dextra[3])
		    dextra[3] = 0.97;
		  // OsiClpSolverInterface * newSolver =
		  fixVubs(model_, extra3, vubMode, generalMessageHandler,
			  debugValues, dextra, extra);
		  // assert (!newSolver);
		}
		// Actually do heuristics
		// may need to flip objective
		bool needFlip = model_.solver()->getObjSense() < 0.0;
		if (needFlip)
		  model_.flipModel();
		// if we do then - fix priorities in
		// clonebutmodel_.convertToDynamic();
		bool objectsExist = model_.objects() != NULL;
		if (!objectsExist) {
		  model_.findIntegers(false);
		  model_.convertToDynamic();
		}
		// set priorities etc
		if (priorities) {
		  OsiObject **objects = model_.objects();
		  int numberObjects = model_.numberObjects();
		  for (int iObj = 0; iObj < numberObjects; iObj++) {
		    CbcSimpleInteger *obj =
		      dynamic_cast<CbcSimpleInteger *>(objects[iObj]);
		    if (!obj)
		      continue;
		    int iColumn = obj->columnNumber();
		    if (branchDirection) {
		      obj->setPreferredWay(branchDirection[iColumn]);
		    }
		    if (priorities) {
		      int iPriority = priorities[iColumn];
		      if (iPriority > 0)
                        obj->setPriority(iPriority);
		    }
		    if (pseudoUp && pseudoUp[iColumn]) {
		      CbcSimpleIntegerPseudoCost *obj1a =
                        dynamic_cast<CbcSimpleIntegerPseudoCost *>(
								   objects[iObj]);
		      assert(obj1a);
		      if (pseudoDown[iColumn] > 0.0)
			obj1a->setDownPseudoCost(pseudoDown[iColumn]);
		      if (pseudoUp[iColumn] > 0.0)
			obj1a->setUpPseudoCost(pseudoUp[iColumn]);
		    }
		  }
		}
		doHeuristics(&model_, 2, parameters, parameters.noPrinting(),
			     initialPumpTune);
		if (!objectsExist) {
		  model_.deleteObjects(false);
		}
		if (needFlip)
		  model_.flipModel();
		if (model_.bestSolution()) {
		  model_.setProblemStatus(1);
		  model_.setSecondaryStatus(6);
		  if (statusUserFunction_[0] && info) {
		    double value = model_.getObjValue();
		    char buf[300];
		    int pos = 0;
		    pos += sprintf(buf + pos, "feasible,");
		    info->problemStatus = 0;
		    info->objValue = value;
		    pos += sprintf(buf + pos, " objective %.*g", ampl_obj_prec(),
				   value);
		    sprintf(buf + pos, "\n0 iterations");
		    free(info->primalSolution);
		    int numberColumns = lpSolver->numberColumns();
		    info->primalSolution = reinterpret_cast<double *>(
								      malloc(numberColumns * sizeof(double)));
		    CoinCopyN(model_.bestSolution(), numberColumns,
			      info->primalSolution);
		    int numberRows = lpSolver->numberRows();
		    free(info->dualSolution);
		    info->dualSolution = reinterpret_cast<double *>(
								    malloc(numberRows * sizeof(double)));
		    CoinZeroN(info->dualSolution, numberRows);
		    CoinWarmStartBasis *basis = lpSolver->getBasis();
		    free(info->rowStatus);
		    info->rowStatus =
                      reinterpret_cast<int *>(malloc(numberRows * sizeof(int)));
		    free(info->columnStatus);
		    info->columnStatus = reinterpret_cast<int *>(
								 malloc(numberColumns * sizeof(int)));
		    // Put basis in
		    int i;
		    // free,basic,ub,lb are 0,1,2,3
		    for (i = 0; i < numberRows; i++) {
		      CoinWarmStartBasis::Status status =
                        basis->getArtifStatus(i);
		      info->rowStatus[i] = status;
		    }
		    for (i = 0; i < numberColumns; i++) {
		      CoinWarmStartBasis::Status status =
                        basis->getStructStatus(i);
		      info->columnStatus[i] = status;
		    }
		    // put buffer into info
		    strcpy(info->buffer, buf);
		    delete basis;
		  }
		}
		int returnCode = callBack(&model, 6);
		if (returnCode) {
		  // exit if user wants
		  delete babModel_;
		  babModel_ = NULL;
		  return returnCode;
		}
	      }
	      //HEUREND
              bool hasTimePreproc = !babModel_->maximumSecondsReached();
              if (!hasTimePreproc)
                preProcess = 0;
              // See if we need to skip preprocessing
              if (preProcess) {
                int numberGenerators = model_.numberCutGenerators();
                int iGenerator;
                for (iGenerator = 0; iGenerator < numberGenerators;
                     iGenerator++) {
                  CglCutGenerator *generator =
                      babModel_->cutGenerator(iGenerator)->generator();

                  if (generator->needsOriginalModel())
                    break;
                }
                if (iGenerator < numberGenerators) {
                  preProcess = 0;
                  printGeneralMessage(model_,
                        "PreProcessing switched off due to lazy constraints");
                }
              }
              if (preProcess && cbcParamCode == CbcParam::BAB) {
                saveSolver = babModel_->solver()->clone();
                /* Do not try and produce equality cliques and
                                   do up to 10 passes */
                OsiSolverInterface *solver2;
                {
                  // Tell solver we are in Branch and Cut
                  saveSolver->setHintParam(OsiDoInBranchAndCut, true,
                                           OsiHintDo);
                  // Default set of cut generators
                  CglProbing generator1;
                  generator1.setUsingObjective(1);
                  generator1.setMaxPass(1);
                  generator1.setMaxPassRoot(1);
                  generator1.setMaxProbeRoot(
                      std::min(3000, saveSolver->getNumCols()));
                  generator1.setMaxElements(100);
                  generator1.setMaxElementsRoot(200);
                  generator1.setMaxLookRoot(50);
                  if (saveSolver->getNumCols() > 3000)
                    generator1.setMaxProbeRoot(123);
                  generator1.setRowCuts(3);
                  // switch off duplicate columns if we have a solution
                  if (model_.bestSolution() /*||debugValues*/)
                    tunePreProcess |= 4096;
                  // take off top
                  int tune2 = tunePreProcess % 10000;
                  if ((tune2 & (1 | 512)) != 0) {
                    // heavy probing
                    generator1.setMaxPassRoot(2);
#ifndef CBC_EXPERIMENT_JJF
                    generator1.setMaxElements(1000);
#else
                    if ((tune2 & 512) != 0)
                      generator1.setMaxElementsRoot(saveSolver->getNumCols());
                    else
                      generator1.setMaxElements(1000);
#endif
                    generator1.setMaxProbeRoot(saveSolver->getNumCols());
 		    if ((tune2 & 512) != 0) 
		      generator1.setMaxLookRoot(std::min(saveSolver->getNumCols(),1000));
		    else
		      generator1.setMaxLookRoot(std::min(saveSolver->getNumCols(),400));
                  }
                  if ((babModel_->specialOptions() & 65536) != 0)
                    process.setOptions(1);
                  // Add in generators
                  if ((model_.moreSpecialOptions() & 65536) == 0)
                    process.addCutGenerator(&generator1);
                  int translate[] = { 9999, 0, 0, -3, 2, 3, -2, 9999, 4, 5, 0 , -2};
                  process.passInMessageHandler(babModel_->messageHandler());
                  // process.messageHandler()->setLogLevel(babModel_->logLevel());
                  if (info && info->numberSos && doSOS && statusUserFunction_[0]) {
                    // SOS
                    numberSOS = info->numberSos;
                    sosStart = info->sosStart;
                    sosIndices = info->sosIndices;
                  }
                  if (numberSOS && doSOS) {
                    // SOS
                    int numberColumns = saveSolver->getNumCols();
                    char *prohibited = new char[numberColumns];
                    memset(prohibited, 0, numberColumns);
                    // worth looking to see if any members can be made integer

                    int numberRows = saveSolver->getNumRows();
                    const CoinPackedMatrix *matrixByCol =
                        saveSolver->getMatrixByCol();
                    const double *element = matrixByCol->getElements();
                    const int *row = matrixByCol->getIndices();
                    const CoinBigIndex *columnStart =
                        matrixByCol->getVectorStarts();
                    const int *columnLength = matrixByCol->getVectorLengths();
                    const double *columnLower = saveSolver->getColLower();
                    const double *columnUpper = saveSolver->getColUpper();
                    const double *rowLower = saveSolver->getRowLower();
                    const double *rowUpper = saveSolver->getRowUpper();
                    double *sameElement = new double[numberRows];
                    int *rowCount = new int[2 * numberRows];
                    int *rowUsed = rowCount + numberRows;
                    memset(sameElement, 0, numberRows * sizeof(double));
                    memset(rowCount, 0, numberRows * sizeof(int));
                    int numberInteresting1 = 0;
                    int numberInteresting2 = 0;
                    int numberChanged = 0;
                    for (int iSet = 0; iSet < numberSOS; iSet++) {
                      if (sosType[iSet] != 1) {
                        for (int i = sosStart[iSet]; i < sosStart[iSet + 1];
                             i++) {
                          numberInteresting2++;
                          int iColumn = sosIndices[i];
                          prohibited[iColumn] = 1;
                        }
                      } else {
                        int nUsed = 0;
                        for (int i = sosStart[iSet]; i < sosStart[iSet + 1];
                             i++) {
                          int iColumn = sosIndices[i];
                          for (CoinBigIndex j = columnStart[iColumn];
                               j < columnStart[iColumn] + columnLength[iColumn];
                               j++) {
                            int iRow = row[j];
                            double el = element[j];
                            if (rowCount[iRow]) {
                              if (el != sameElement[iRow])
                                sameElement[iRow] = 0.0;
                            } else {
                              sameElement[iRow] = el;
                              rowUsed[nUsed++] = iRow;
                            }
                            rowCount[iRow]++;
                          }
                        }
                        int nInSet = sosStart[iSet + 1] - sosStart[iSet];
                        double nonzeroValue = COIN_DBL_MAX;
                        for (int iUsed = 0; iUsed < nUsed; iUsed++) {
                          int iRow = rowUsed[iUsed];
                          if (rowCount[iRow] == nInSet && sameElement[iRow] &&
                              rowLower[iRow] == rowUpper[iRow]) {
                            // all entries must be 0.0 or xx
                            nonzeroValue = rowLower[iRow] / sameElement[iRow];
                          }
                          rowCount[iRow] = 0;
                          sameElement[iRow] = 0.0;
                        }
                        if (nonzeroValue != COIN_DBL_MAX) {
                          // could do scaling otherwise
                          if (fabs(nonzeroValue - 1.0) < 1.0e-8) {
                            for (int i = sosStart[iSet]; i < sosStart[iSet + 1];
                                 i++) {
                              int iColumn = sosIndices[i];
                              if (columnUpper[iColumn] < 0.0 ||
                                  columnLower[iColumn] > 1.0) {
                                printf("sos says infeasible\n");
                              }
                              if (!saveSolver->isInteger(iColumn)) {
                                numberChanged++;
                                saveSolver->setInteger(iColumn);
                              }
                              if (columnUpper[iColumn] < 1.0)
                                saveSolver->setColUpper(iColumn, 0.0);
                              else
                                saveSolver->setColUpper(iColumn, 1.0);
                              if (columnLower[iColumn] > 0.0)
                                saveSolver->setColLower(iColumn, 1.0);
                              else
                                saveSolver->setColLower(iColumn, 0.0);
#ifndef DO_LESS_PROHIBITED
                              prohibited[iColumn] = 1;
#endif
                            }
                          } else {
                            for (int i = sosStart[iSet]; i < sosStart[iSet + 1];
                                 i++) {
                              int iColumn = sosIndices[i];
#ifndef DO_LESS_PROHIBITED
                              prohibited[iColumn] = 1;
#endif
                            }
                          }
                        } else {
                          for (int i = sosStart[iSet]; i < sosStart[iSet + 1];
                               i++) {
                            int iColumn = sosIndices[i];
                            if (!saveSolver->isInteger(iColumn))
                              numberInteresting1++;
#ifdef DO_LESS_PROHIBITED
                            if (!saveSolver->isInteger(iColumn))
#endif
                              prohibited[iColumn] = 1;
                          }
                        }
                      }
                    }
                    if (numberChanged || numberInteresting1 || numberInteresting2) {
                      buffer.str("");
                      buffer << numberChanged << " variables in SOS1 sets made integer, "
                             << numberInteresting1 << " non integer in SOS1, "
                             << numberInteresting2 << " in SOS2" << std::endl;
                      printGeneralMessage(model_, buffer.str());
                    }
                    delete[] sameElement;
                    delete[] rowCount;
                    process.passInProhibited(prohibited, numberColumns);
                    delete[] prohibited;
                  }
                  if (0) {

                    // Special integers
                    int numberColumns = saveSolver->getNumCols();
                    char *prohibited = new char[numberColumns];
                    memset(prohibited, 0, numberColumns);
                    const CoinPackedMatrix *matrix =
                        saveSolver->getMatrixByCol();
                    const int *columnLength = matrix->getVectorLengths();
                    int numberProhibited = 0;
                    for (int iColumn = numberColumns - 1; iColumn >= 0;
                         iColumn--) {
                      if (!saveSolver->isInteger(iColumn) ||
                          columnLength[iColumn] > 1)
                        break;
                      numberProhibited++;
                      prohibited[iColumn] = 1;
                    }
                    if (numberProhibited) {
                      process.passInProhibited(prohibited, numberColumns);
                      printf("**** Treating last %d integers as special - give "
                             "high priority?\n",
                             numberProhibited);
                    }
                    delete[] prohibited;
                  }
                  if (!model_.numberObjects() && true) {
                    /* model may not have created objects
                                           If none then create
                                        */
                    model_.findIntegers(true);
                  }
                  // Lotsizing
                  if (numberLotSizing) {
                    int numberColumns = saveSolver->getNumCols();
                    int numberRows = saveSolver->getNumRows();
#if 0
		    // Create model which uses 0-1 variables
		    {
		      double * els = new double[8*numberLotSizing];
		      int * cols = new int[4*numberLotSizing];
		      CoinBigIndex * starts = new CoinBigIndex [2*numberLotSizing+1];
		      // add 0-1 variables
		      double * tup = els+4*numberLotSizing;
		      double * tlo = tup+2*numberLotSizing;
		      for (int i=0;i<numberLotSizing;i++) {
			tlo[i]=0.0;
			tup[i]=1.0;
			els[i]=0.0; // objective
		      }
		      OsiClpSolverInterface *si = dynamic_cast< OsiClpSolverInterface * >(saveSolver);
		      assert(si != NULL);
		      // get clp itself
		      ClpSimplex *model = si->getModelPtr();
		      model->tightenPrimalBounds(0.0,0);
		      model->addColumns(numberLotSizing,tlo,tup,els,NULL,NULL,NULL);
		      for (int i=0;i<numberLotSizing;i++) 
			model->setInteger(i+numberColumns);
		      starts[0]=0;
		      CoinBigIndex n=0;
		      const double * columnUpper = model->columnUpper();
		      for (int i=0;i<numberLotSizing;i++) {
			int iColumn = lotsize[i].column;
			tlo[n] = -COIN_DBL_MAX;
			tup[n]=0.0;
			els[2*n] = 1.0;
			// should be able to get correct value
			els[2*n+1] = -std::min(10000.0,columnUpper[iColumn]);
			cols[2*n] = iColumn;
			cols[2*n+1] = i+numberColumns;
			n++;
			starts[n]=2*n;
			tlo[n] = 0.0;
			tup[n]=COIN_DBL_MAX;
			els[2*n] = 1.0;
			els[2*n+1] = -lotsize[i].low;
			cols[2*n] = iColumn;
			cols[2*n+1] = i+numberColumns;
			n++;
			starts[n]=2*n;
		      }
		      model->addRows(n,tlo,tup,starts,cols,els);
		      model->writeMps("equivalentInteger.mps");
		      exit(1);
		    }
#endif
                    char *prohibited = new char[numberColumns+numberRows];
		    char *prohibitedRow = prohibited+numberColumns;
                    memset(prohibited, 0, numberColumns+numberRows);
                    const CoinPackedMatrix *matrix = saveSolver->getMatrixByCol();
                    //const double *element = matrix->getElements();
                    const int *row = matrix->getIndices();
                    const CoinBigIndex *columnStart = matrix->getVectorStarts();
                    const int *columnLength = matrix->getVectorLengths();
                    for (int i = 0; i < numberLotSizing; i++) {
                      int iColumn = lotsize[i].column;
                      prohibited[iColumn] = 2;
#if 0
		      // also leave all rows
		      for (CoinBigIndex j = columnStart[iColumn];
			   j < columnStart[iColumn] + columnLength[iColumn]; j++) {
			int iRow = row[j];
			prohibitedRow[iRow] = -2;
		      }
#endif
                    }
                    process.passInProhibited(prohibited, numberColumns);
                    process.passInRowTypes(prohibitedRow, numberRows);
		    passSCtopreprocess.lotsize=lotsize;
		    passSCtopreprocess.numberLotSizing=numberLotSizing;
		    process.setApplicationData(&passSCtopreprocess);
		    process.setOptions(128|256|process.options());
                    delete[] prohibited;
                  }
                  if (model_.numberObjects()) {
                    OsiObject **oldObjects = babModel_->objects();
                    int numberOldObjects = babModel_->numberObjects();
                    if (!numberOldObjects) {
                      oldObjects = model_.objects();
                      numberOldObjects = model_.numberObjects();
                    }
                    // SOS
                    int numberColumns = saveSolver->getNumCols();
                    char *prohibited = new char[numberColumns];
                    memset(prohibited, 0, numberColumns);
                    int numberProhibited = 0;
                    for (int iObj = 0; iObj < numberOldObjects; iObj++) {
                      CbcSOS *obj = dynamic_cast<CbcSOS *>(oldObjects[iObj]);
                      if (obj) {
                        int n = obj->numberMembers();
                        const int *which = obj->members();
                        for (int i = 0; i < n; i++) {
                          int iColumn = which[i];
#ifdef DO_LESS_PROHIBITED
                          if (!saveSolver->isInteger(iColumn))
#endif
                            prohibited[iColumn] = 1;
                          numberProhibited++;
                        }
                      }
                      CbcLotsize *obj2 =
                          dynamic_cast<CbcLotsize *>(oldObjects[iObj]);
                      if (obj2) {
                        int iColumn = obj2->columnNumber();
                        prohibited[iColumn] = 1;
                        numberProhibited++;
                      }
                    }
                    if (numberProhibited) 
                      process.passInProhibited(prohibited, numberColumns);
                    delete[] prohibited;
                  }
                  int numberPasses = 10;
#ifndef CBC_OTHER_SOLVER
                  if (doSprint > 0) {
                    // Sprint for primal solves
                    ClpSolve::SolveType method = ClpSolve::usePrimalorSprint;
                    ClpSolve::PresolveType presolveType = ClpSolve::presolveOff;
                    int numberPasses = 5;
                    int options[] = {0, 3, 0, 0, 0, 0};
                    int extraInfo[] = {-1, 20, -1, -1, -1, -1};
                    extraInfo[1] = doSprint;
                    int independentOptions[] = {0, 0, 3};
                    ClpSolve clpSolve(method, presolveType, numberPasses,
                                      options, extraInfo, independentOptions);
                    // say use in OsiClp
                    clpSolve.setSpecialOption(6, 1);
                    OsiClpSolverInterface *osiclp =
                        dynamic_cast<OsiClpSolverInterface *>(saveSolver);
                    osiclp->setSolveOptions(clpSolve);
                    osiclp->setHintParam(OsiDoDualInResolve, false);
                    // switch off row copy
                    osiclp->getModelPtr()->setSpecialOptions(
                        osiclp->getModelPtr()->specialOptions() | 256);
                    osiclp->getModelPtr()->setInfeasibilityCost(1.0e11);
                  }
#endif
#ifndef CBC_OTHER_SOLVER
                  {
                    OsiClpSolverInterface *osiclp =
                        dynamic_cast<OsiClpSolverInterface *>(saveSolver);
                    osiclp->setSpecialOptions(osiclp->specialOptions() | 1024);
		    osiclp->getModelPtr()->setClpScaledMatrix(NULL); // safer
                    int savePerturbation =
                       osiclp->getModelPtr()->perturbation();
                    //#define CBC_TEMP1
#ifdef CBC_TEMP1
                    if (savePerturbation == 50)
                      osiclp->getModelPtr()->setPerturbation(52); // try less
#endif
                    if ((model_.moreSpecialOptions() & 65536) != 0)
                      process.setOptions(2 + 4 + 8); // no cuts
                    cbcPreProcessPointer = &process;
                    preProcessPointer = &process; // threadsafe
                    int saveOptions =
                        osiclp->getModelPtr()->moreSpecialOptions();
                    if ((model_.specialOptions() & 16777216) != 0 &&
                        model_.getCutoff() > 1.0e30) {
                      osiclp->getModelPtr()->setMoreSpecialOptions(saveOptions |
                                                                   262144);
                    }
#if DEBUG_PREPROCESS > 1
                    if (debugValues) {
                      process.setApplicationData(
                          const_cast<double *>(debugValues));
                    }
#endif
                    if (parameters[CbcParam::DEBUGFILE]->fileName() ==
                        "unitTest") {
                      //This is probably wrong, will need to debug
                      babModel_->solver()->activateRowCutDebugger(inputQueue[0].c_str());
                      OsiRowCutDebugger *debugger =
                          babModel_->solver()->getRowCutDebuggerAlways();
                      numberDebugValues = babModel_->getNumCols();
                      debugValues = CoinCopyOfArray(debugger->optimalSolution(),
                                                    numberDebugValues);
                    }
                    redoSOS = true;
                    bool keepPPN = parameters[CbcParam::PREPROCNAMES]->modeVal();
#ifdef SAVE_NAUTY
                    keepPPN = 1;
#endif
                    process.setKeepColumnNames(keepPPN);
                    process.setTimeLimit(babModel_->getMaximumSeconds() -
                                             babModel_->getCurrentSeconds(),
                                         babModel_->useElapsedTime());
                    if (model.getKeepNamesPreproc())
                      process.setKeepColumnNames(true);
		    if (keepPPN)
		      babModel_->setKeepNamesPreproc(1);
		    setPreProcessingMode(saveSolver,1);
#if CBC_USE_PAPILO
		    extern void zapPapilo(int pOptions,CglPreProcess * process);
		    int pOptions = 0;
		    int tune2 = preProcess;
		    // Convert to minimize if papilo
		    bool maximize = false;
		    if (tune2>11) {
		      OsiClpSolverInterface * clpSolver =
			dynamic_cast<OsiClpSolverInterface *>(saveSolver);
		      if (clpSolver->getObjSense()==-1.0) {
			maximize = true;
			clpSolver->setObjSense(1.0);
			double objOffset;
			clpSolver->getDblParam(OsiObjOffset, objOffset);
			int numberColumns = clpSolver->getNumCols();
			double * objective = clpSolver->getModelPtr()->objective();
			for (int i=0;i<numberColumns;i++)
			  objective[i] = -objective[i];
			clpSolver->setDblParam(OsiObjOffset, -objOffset);
		      }
		      bool stopAfter = false;
		      if (tune2>15) {
			preProcess=10; // say stop
			tune2 -=4;
		      } else {
			preProcess=1;
		      }
#ifdef CBC_THREAD
		      pOptions = (tune2&1) != 0 ? 2 : 1;
#endif
		      if ((tune2&2) == 0)
			pOptions|= 8; // at beginning
		      else
			pOptions|= 16; // at end
		    }
		    zapPapilo(pOptions,&process);
#endif
                    solver2 = process.preProcessNonDefault(*saveSolver, translate[preProcess], numberPasses,
                      tunePreProcess);
 		    setPreProcessingMode(saveSolver,0);
#if CBC_USE_PAPILO
		    // Convert back
		    if (maximize) {
		      OsiClpSolverInterface * clpSolver =
			dynamic_cast<OsiClpSolverInterface *>(saveSolver);
		      double objOffset;
		      clpSolver->setObjSense(-1.0);
		      clpSolver->getDblParam(OsiObjOffset, objOffset);
		      int numberColumns = clpSolver->getNumCols();
		      double * objective = clpSolver->getModelPtr()->objective();
		      for (int i=0;i<numberColumns;i++)
			objective[i] = -objective[i];
		      clpSolver->setDblParam(OsiObjOffset, -objOffset);
		      if (solver2) {
			OsiClpSolverInterface * clpSolver =
			  dynamic_cast<OsiClpSolverInterface *>(solver2);
			double objOffset;
			clpSolver->setObjSense(-1.0);
			clpSolver->getDblParam(OsiObjOffset, objOffset);
			int numberColumns = clpSolver->getNumCols();
			double * objective = clpSolver->getModelPtr()->objective();
			for (int i=0;i<numberColumns;i++)
			  objective[i] = -objective[i];
			clpSolver->setDblParam(OsiObjOffset, -objOffset);
		      }
		    }
#endif
                    if (solver2) {
                      setPreProcessingMode(solver2, 0);
                      model_.setOriginalColumns(process.originalColumns(),
                                                solver2->getNumCols());

                      osiclp->getModelPtr()->setPerturbation(savePerturbation);
                      osiclp->getModelPtr()->setMoreSpecialOptions(saveOptions);
                      /* clean solvers - should be done in preProcess but
                         that doesn't know about Clp */
                      OsiClpSolverInterface *solver;
                      solver = dynamic_cast<OsiClpSolverInterface *>(solver2);
                      solver->getModelPtr()->cleanScalingEtc();
                      solver = dynamic_cast<OsiClpSolverInterface *>(
                          process.originalModel());
                      solver->getModelPtr()->cleanScalingEtc();
                      solver = dynamic_cast<OsiClpSolverInterface *>(
                          process.startModel());
                      if (solver)
                        solver->getModelPtr()->cleanScalingEtc();
                      int numberSolvers = process.numberSolvers();
                      if (numberSolvers == 99)
                        numberSolvers = 1; // really just 1
                      // some of these may be same
                      for (int i = 0; i < numberSolvers; i++) {
                        solver = dynamic_cast<OsiClpSolverInterface *>(
                            process.modelAtPass(i));
                        if (solver)
                          solver->getModelPtr()->cleanScalingEtc();
                        solver = dynamic_cast<OsiClpSolverInterface *>(
                            process.modifiedModel(i));
                        if (solver)
                          solver->getModelPtr()->cleanScalingEtc();
                      }
                    }
                  }
#elif CBC_OTHER_SOLVER == 1
                  cbcPreProcessPointer = &process;
                  preProcessPointer = &process; // threadsafe
                  redoSOS = true;
                  if (model.getKeepNamesPreproc())
                    process.setKeepColumnNames(true);
                  solver2 = process.preProcessNonDefault(
                      *saveSolver, translate[preProcess], numberPasses,
                      tunePreProcess);
#endif
                  integersOK = false; // We need to redo if CbcObjects exist
                  // Tell solver we are not in Branch and Cut
                  saveSolver->setHintParam(OsiDoInBranchAndCut, false,
                                           OsiHintDo);
                  if (solver2)
                    solver2->setHintParam(OsiDoInBranchAndCut, false,
                                          OsiHintDo);
                }
                if (info && !solver2 && statusUserFunction_[0]) {
                  // infeasible
                  info->problemStatus = 1;
                  info->objValue = 1.0e100;
                  sprintf(info->buffer,
                          "infeasible/unbounded by pre-processing");
                  info->primalSolution = NULL;
                  info->dualSolution = NULL;
                  break;
                }
                if (!solver2) {
                   printGeneralMessage(model_,
                                "Pre-processing says infeasible or unbounded");
                  // say infeasible for solution
                  integerStatus = 6;
                  delete saveSolver;
                  saveSolver = NULL;
                  model_.setProblemStatus(0);
                  model_.setSecondaryStatus(1);
                  babModel_->setProblemStatus(0);
                  babModel_->setSecondaryStatus(1);
                } else {
                    // printf("processed model has %d rows, %d columns and %d
                    // elements\n",
                    //     solver2->getNumRows(),solver2->getNumCols(),solver2->getNumElements());
                  statistics_nprocessedrows = solver2->getNumRows();
                  statistics_nprocessedcols = solver2->getNumCols();
                  model_.setProblemStatus(-1);
                  babModel_->setProblemStatus(-1);
                }
                int returnCode = 0;
                if (callBack != NULL)
                  returnCode = callBack(babModel_, 2);
                if (returnCode) {
                  // exit if user wants
                  delete babModel_;
                  babModel_ = NULL;
                  return returnCode;
                }
                if (!solver2)
                  break;
                if (model_.bestSolution()) {
                  // need to redo - in case no better found in BAB
                  // just get integer part right
                  const int *originalColumns = process.originalColumns();
                  int numberColumns =
                      std::min(solver2->getNumCols(), babModel_->getNumCols());
#if 0		  
                  double *bestSolution = babModel_->bestSolution();
                  const double *oldBestSolution = model_.bestSolution();
                  for (int i = 0; i < numberColumns; i++) {
                    int jColumn = originalColumns[i];
                    bestSolution[i] = oldBestSolution[jColumn];
                  }
#else
		  int numberColumnsB = babModel_->getNumCols();
                  int numberColumns2 =
                      std::max(solver2->getNumCols(), numberColumnsB);
                  double *bestSolution = new double [numberColumns2];
		  memset(bestSolution,0,numberColumns2*sizeof(double));
                  const double *oldBestSolution = model_.bestSolution();
                  for (int i = 0; i < numberColumns; i++) {
                    int jColumn = originalColumns[i];
		    if (jColumn<numberColumnsB)
		      bestSolution[i] = oldBestSolution[jColumn];
                  }
		  double obj=model_.getObjValue();
                  double newCutoff =
                      std::min(model_.getCutoff(), obj + 1.0e-4);
		  babModel_->setBestSolution(bestSolution,numberColumns,1.0e10,false);
                  babModel_->setCutoff(newCutoff);
		  delete [] bestSolution;
#endif
                }
                //solver2->resolve();
#ifdef CBC_NAMES_FOR_COMPARE
		{
		  OsiClpSolverInterface * solver =
		    dynamic_cast<OsiClpSolverInterface *>(solver2);
		  ClpSimplex * simplex = solver->getModelPtr();
		  int numberRows = simplex->numberRows();
		  int numberColumns = simplex->numberColumns();
		  int numberIntegers = 0;
		  for (int i=0;i<numberColumns;i++) {
		    if (simplex->isInteger(i))
		      numberIntegers++;
		  }
		  printf("CbC %s after preprocessing %d rows, %d columns, %d integers\n",
			 simplex->problemName().c_str(),numberRows,
			 numberColumns,numberIntegers);
		}
#endif
                if (preProcess == 2 || preProcess >= 10) {
                  // names are wrong - redo 
		  const int *originalColumns = process.originalColumns();
                  int numberColumns = solver2->getNumCols();
                  OsiSolverInterface *originalSolver = model.solver();
		  int numberOriginalColumns = originalSolver->getNumCols();
		  int oddColumn=1;
                  for (int i = 0; i < numberColumns; i++) {
                    int iColumn = originalColumns[i];
		    if (iColumn<numberOriginalColumns) {
		      std::string columnName = originalSolver->getColName(iColumn);
		      // check if odd name
		      if (columnName[0]=='N') {
			try {
                           int sequence = std::atoi(columnName.substr(1).c_str());
			  oddColumn = sequence+1;
			} catch (...) {
			}
		      }
		      solver2->setColName(i, columnName);
		    } else {
		      char name[15];
		      sprintf(name,"N%.7d",oddColumn);
		      oddColumn++;
		      solver2->setColName(i, name);
		    }
                  }
                  OsiClpSolverInterface *clpSolver2 = dynamic_cast< OsiClpSolverInterface * >(solver2);
                  ClpSimplex *lpSolver = clpSolver2->getModelPtr();
                  char name[100];
                  if (preProcess == 2) {
                    strcpy(name, "presolved.mps");
                  } else {
                    strcpy(name,
                           parameters[CbcParam::IMPORTFILE]->fileName().c_str());
                    char *dot = strstr(name, ".mps");
                    if (!dot)
                      dot = strstr(name, ".lp");
		    if (preProcess>=10&&!dot) {
		      // get some sort of name
		      dot = name+strlen(name);
		    }
                    if (dot) {
                      *dot = '\0';
                      int n = static_cast<int>(dot - name);
                      int i;
                      for (i = n - 1; i >= 0; i--) {
                        if (name[i] == '/')
                          break;
                      }
                      if (i >= 0)
                        memmove(name, name + i + 1, n);
                      strcat(name, "_preprocessed.mps");
                    } else {
                      strcpy(name, "preprocessed.mps");
                    }
                  }
                  lpSolver->writeMps(name, 0, 1, lpSolver->optimizationDirection());
                  printf("Preprocessed model (minimization) on %s - size %d %d \n",
			 name,lpSolver->getNumRows(),lpSolver->getNumCols());
                  if (preProcess >= 10) {
                    printf("user wanted to stop\n");
                    return 0; //exit(0);
                  }
                }
                {
                  // look at new integers
                  int numberOriginalColumns =
                      process.originalModel()->getNumCols();
                  const int *originalColumns = process.originalColumns();
                  OsiClpSolverInterface *osiclp2 =
                      dynamic_cast<OsiClpSolverInterface *>(solver2);
                  int numberColumns = osiclp2->getNumCols();
                  OsiClpSolverInterface *osiclp =
                      dynamic_cast<OsiClpSolverInterface *>(saveSolver);
                  for (int i = 0; i < numberColumns; i++) {
                    int iColumn = originalColumns[i];
                    if (iColumn < numberOriginalColumns) {
                      if (osiclp2->isInteger(i) && !osiclp->isInteger(iColumn))
                        osiclp2->setOptionalInteger(i); // say optional
                    }
                  }
                  // do lotsizing
                  if (numberLotSizing) {
                    const double *columnLower = solver2->getColLower();
                    const double *columnUpper = solver2->getColUpper();
                    CbcObject **objects = new CbcObject *[numberLotSizing];
                    double points[] = { 0.0, 0.0, 0.0, 0.0 };
                    int *back = new int[numberOriginalColumns];
                    for (int i = 0; i < numberOriginalColumns; i++)
                      back[i] = -1;
                    for (int i = 0; i < numberColumns; i++) {
                      int iColumn = originalColumns[i];
                      back[iColumn] = i;
                    }
		    int n = numberLotSizing;
		    numberLotSizing = 0;
                    for (int i = 0; i < n; i++) {
                      int iColumn = lotsize[i].column;
		      iColumn = back[iColumn];
		      if (iColumn>=0) {
			if (columnLower[iColumn]<lotsize[i].low&&
			    columnUpper[iColumn]>lotsize[i].low-1.0e-8) {
			  points[2] = lotsize[i].low;
			  points[3] = columnUpper[iColumn];
			  objects[numberLotSizing++] =
			    new CbcLotsize(babModel_, iColumn, 2,
					   points, true);
			} else {
			  //printf("SC %d new bounds %g,%g was %d low %g high %g\n",iColumn,
			  //	 columnLower[iColumn],columnUpper[iColumn],
			  //	 lotsize[i].column,lotsize[i].low,lotsize[i].high);
			  if (columnUpper[iColumn]<lotsize[i].low) {
			    if (columnLower[iColumn]) {
			      printf("Infeasible due to SC variables?!\n");
			    } else {
			      solver2->setColUpper(iColumn,0.0);
			    }
			  }
			}
		      }
                    }
		    delete [] back;
                    babModel_->addObjects(numberLotSizing, objects);
                    for (int i = 0; i < numberLotSizing; i++)
                      delete objects[i];
                    delete[] objects;
		    if (numberLotSizing < n) {
                       buffer.str("");
                       buffer << "Pre-processing reduced number of SC"
                              << " variables from " << n << " to "
                              << numberLotSizing;
                       printGeneralMessage(model_, buffer.str());
		    }
                  }
                  // redo existing SOS
                  if (osiclp->numberSOS()) {
                    redoSOS = false;
                    int *back = new int[numberOriginalColumns];
                    for (int i = 0; i < numberOriginalColumns; i++)
                      back[i] = -1;
                    for (int i = 0; i < numberColumns; i++) {
                      int iColumn = originalColumns[i];
                      back[iColumn] = i;
                    }
                    int numberSOSOld = osiclp->numberSOS();
                    int numberSOS = osiclp2->numberSOS();
                    assert(numberSOS == numberSOSOld);
                    CoinSet *setInfo =
                        const_cast<CoinSet *>(osiclp2->setInfo());
                    for (int i = 0; i < numberSOS; i++) {
                      // int type = setInfo[i].setType();
                      int n = setInfo[i].numberEntries();
                      int *which = setInfo[i].modifiableWhich();
#ifndef DO_LESS_PROHIBITED
                      for (int j = 0; j < n; j++) {
                        int iColumn = which[j];
                        iColumn = back[iColumn];
                        assert(iColumn >= 0);
                        which[j] = iColumn;
                      }
#else
                      double *weights = setInfo[i].modifiableWeights();
                      int n2 = 0;
                      for (int j = 0; j < n; j++) {
                        int iColumn = which[j];
                        iColumn = back[iColumn];
                        if (iColumn >= 0) {
                          which[n2] = iColumn;
                          weights[n2++] = weights[j];
                        }
                      }
                      setInfo[i].setNumberEntries(n2);
#endif
                    }
                    delete[] back;
                  }
                }
                // we have to keep solver2 so pass clone
                solver2 = solver2->clone();
                // see if extra variables wanted
                int threshold = parameters[CbcParam::EXTRAVARIABLES]->intVal();
                int more2 = parameters[CbcParam::MOREMOREMIPOPTIONS]->intVal();
                if (threshold || (more2 & (512 | 1024)) != 0) {
                  int numberColumns = solver2->getNumCols();
                  truncateRows = solver2->getNumRows();
                  bool modifiedModel = false;
                  int highPriority = 0;
                  /*
                                    normal - no priorities
                                    >10000 equal high priority
                                    >20000 higher priority for higher cost
                                  */
                  if (threshold > 10000) {
                    highPriority = threshold / 10000;
                    threshold -= 10000 * highPriority;
                  }
                  // If 1000 set then don't put obj on ne variables
                  bool moveObjective = true;
                  if (threshold > 1000) {
                    moveObjective = false;
                    threshold -= 1000;
                  }
                  const double *columnLower = solver2->getColLower();
                  const double *columnUpper = solver2->getColUpper();
                  const double *objective = solver2->getObjCoefficients();
                  int numberIntegers = 0;
                  int numberBinary = 0;
                  int numberTotalIntegers = 0;
                  double *obj = new double[numberColumns];
                  int *which = new int[numberColumns];
                  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
                    if (solver2->isInteger(iColumn)) {
                      numberTotalIntegers++;
                      if (columnUpper[iColumn] > columnLower[iColumn]) {
                        numberIntegers++;
                        if (columnLower[iColumn] == 0.0 &&
                            columnUpper[iColumn] == 1)
                          numberBinary++;
                      }
                    }
                  }
                  int numberSort = 0;
                  int numberZero = 0;
                  int numberZeroContinuous = 0;
                  int numberDifferentObj = 0;
                  int numberContinuous = 0;
                  for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
                    if (columnUpper[iColumn] > columnLower[iColumn]) {
                      if (solver2->isInteger(iColumn)) {
                        if (!objective[iColumn]) {
                          numberZero++;
                        } else {
                          obj[numberSort] = fabs(objective[iColumn]);
                          which[numberSort++] = iColumn;
                        }
                      } else if (objective[iColumn]) {
                        numberContinuous++;
                      } else {
                        numberZeroContinuous++;
                      }
                    }
                  }
                  CoinSort_2(obj, obj + numberSort, which);
                  double last = obj[0];
                  for (int jColumn = 1; jColumn < numberSort; jColumn++) {
                    if (fabs(obj[jColumn] - last) > 1.0e-12) {
                      numberDifferentObj++;
                      last = obj[jColumn];
                    }
                  }
                  numberDifferentObj++;
                  buffer.str("");
                  buffer << "Problem has " << numberIntegers << " integers (" << numberBinary
                         << " of which are binary) and " << numberColumns - numberIntegers
                         << " continuous";
                  printGeneralMessage(model_, buffer.str());
                  if (numberColumns > numberIntegers) {
                    buffer.str("");
                    buffer << numberContinuous << " continuous have nonzero objective, "
                           << numberZeroContinuous << " have zero objective";
                    printGeneralMessage(model_, buffer.str());
                  }
                  buffer.str("");
                  buffer << numberSort << " integer have nonzero objective, "
                         << numberZero << " have zero objective, " 
                         << numberDifferentObj << " different nonzero (taking abs)";
                  printGeneralMessage(model_, buffer.str());
                  if (numberDifferentObj <= threshold + (numberZero)
                          ? 1
                          : 0 && numberDifferentObj) {
                    int *backward = NULL;
                    if (highPriority) {
                      newPriorities =
                          new int[numberTotalIntegers + numberDifferentObj +
                                  numberColumns];
                      backward = newPriorities + numberTotalIntegers +
                                 numberDifferentObj;
                      numberTotalIntegers = 0;
                      for (int iColumn = 0; iColumn < numberColumns;
                           iColumn++) {
                        if (solver2->isInteger(iColumn)) {
                          backward[iColumn] = numberTotalIntegers;
                          newPriorities[numberTotalIntegers++] = 10000;
                        }
                      }
                    }
                    int iLast = 0;
                    double last = obj[0];
                    for (int jColumn = 1; jColumn < numberSort; jColumn++) {
                      if (fabs(obj[jColumn] - last) > 1.0e-12) {
                        buffer.str("");
                        buffer << jColumn - iLast << " variables have objective of "
                               << last;
                        printGeneralMessage(model_, buffer.str());
                        iLast = jColumn;
                        last = obj[jColumn];
                      }
                    }
                    buffer.str("");
                    buffer << numberSort - iLast << " variables have objective of "
                           << last;
                    printGeneralMessage(model_, buffer.str());
                    int spaceNeeded = numberSort + numberDifferentObj;
                    CoinBigIndex *columnAddDummy =
                        new CoinBigIndex[numberDifferentObj + 1];
                    int *columnAdd = new int[spaceNeeded];
                    double *elementAdd = new double[spaceNeeded];
                    CoinBigIndex *rowAdd =
                        new CoinBigIndex[numberDifferentObj + 1];
                    double *objectiveNew = new double[3 * numberDifferentObj];
                    double *lowerNew = objectiveNew + numberDifferentObj;
                    double *upperNew = lowerNew + numberDifferentObj;
                    memset(columnAddDummy, 0,
                           (numberDifferentObj + 1) * sizeof(CoinBigIndex));
                    iLast = 0;
                    last = obj[0];
                    numberDifferentObj = 0;
                    int priorityLevel = 9999;
                    int numberElements = 0;
                    rowAdd[0] = 0;
                    for (int jColumn = 1; jColumn < numberSort + 1; jColumn++) {
                      if (jColumn == numberSort ||
                          fabs(obj[jColumn] - last) > 1.0e-12) {
                        // not if just one
                        if (jColumn - iLast > 1) {
                          // do priority
                          if (highPriority == 1) {
                            newPriorities[numberTotalIntegers +
                                          numberDifferentObj] = 500;
                          } else if (highPriority == 2) {
                            newPriorities[numberTotalIntegers +
                                          numberDifferentObj] = priorityLevel;
                            priorityLevel--;
                          }
                          int iColumn = which[iLast];
                          objectiveNew[numberDifferentObj] = objective[iColumn];
                          double lower = 0.0;
                          double upper = 0.0;
                          for (int kColumn = iLast; kColumn < jColumn;
                               kColumn++) {
                            iColumn = which[kColumn];
                            if (moveObjective)
                              solver2->setObjCoeff(iColumn, 0.0);
                            double lowerValue = columnLower[iColumn];
                            double upperValue = columnUpper[iColumn];
                            double elementValue = -1.0;
                            if (objectiveNew[numberDifferentObj] *
                                    objective[iColumn] <
                                0.0) {
                              lowerValue = -columnUpper[iColumn];
                              upperValue = -columnLower[iColumn];
                              elementValue = 1.0;
                            }
                            if (!moveObjective)
                              objectiveNew[numberDifferentObj] = 0.0;
                            columnAdd[numberElements] = iColumn;
                            elementAdd[numberElements++] = elementValue;
                            if (lower != -COIN_DBL_MAX) {
                              if (lowerValue != -COIN_DBL_MAX)
                                lower += lowerValue;
                              else
                                lower = -COIN_DBL_MAX;
                            }
                            if (upper != COIN_DBL_MAX) {
                              if (upperValue != COIN_DBL_MAX)
                                upper += upperValue;
                              else
                                upper = COIN_DBL_MAX;
                            }
                          }
                          columnAdd[numberElements] =
                              numberColumns + numberDifferentObj;
                          elementAdd[numberElements++] = 1.0;
                          lowerNew[numberDifferentObj] = lower;
                          upperNew[numberDifferentObj] = upper;
                          numberDifferentObj++;
                          rowAdd[numberDifferentObj] = numberElements;
                        } else if (highPriority) {
                          // just one
                          // do priority
                          int iColumn = which[iLast];
                          int iInt = backward[iColumn];
                          if (highPriority == 1) {
                            newPriorities[iInt] = 500;
                          } else {
                            newPriorities[iInt] = priorityLevel;
                            priorityLevel--;
                          }
                        }
                        if (jColumn < numberSort) {
                          iLast = jColumn;
                          last = obj[jColumn];
                        }
                      }
                    }
                    if (numberDifferentObj) {
                      // add columns
                      solver2->addCols(numberDifferentObj, columnAddDummy, NULL,
                                       NULL, lowerNew, upperNew, objectiveNew);
                      // add constraints and make integer if all integer in
                      // group
                      OsiClpSolverInterface *clpSolver2 =
                          dynamic_cast<OsiClpSolverInterface *>(solver2);
                      for (int iObj = 0; iObj < numberDifferentObj; iObj++) {
                        lowerNew[iObj] = 0.0;
                        upperNew[iObj] = 0.0;
                        solver2->setInteger(numberColumns + iObj);
                        if (clpSolver2)
                          clpSolver2->setOptionalInteger(numberColumns + iObj);
                      }
                      solver2->addRows(numberDifferentObj, rowAdd, columnAdd,
                                       elementAdd, lowerNew, upperNew);
                      buffer.str("");
                      buffer << "Replacing model - "
                             << numberDifferentObj << " new variables";
                      modifiedModel = true;
                    }
                    delete[] columnAdd;
                    delete[] columnAddDummy;
                    delete[] elementAdd;
                    delete[] rowAdd;
                    delete[] objectiveNew;
                  }
                  delete[] which;
                  delete[] obj;
                  if ((more2 & (512 | 1024)) != 0) {
                    // try for row slacks etc
                    // later do row branching
                    int iRow, iColumn;
                    int numberColumns = solver2->getNumCols();
                    int numberRows = solver2->getNumRows();
                    int fudgeObjective = more2 & 512;
                    int addSlacks = more2 & 1024;
                    if (fudgeObjective) {
                      bool moveObj = false;
                      fudgeObjective = 0;
                      const double *objective = solver2->getObjCoefficients();
                      const double *columnLower = solver2->getColLower();
                      const double *columnUpper = solver2->getColUpper();
                      double *newValues = new double[numberColumns + 1];
                      int *newColumn = new int[numberColumns + 1];
                      bool allInteger = true;
                      int n = 0;
                      double newLower = 0.0;
                      double newUpper = 0.0;
                      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                        if (objective[iColumn]) {
                          if (!solver2->isInteger(iColumn)) {
                            allInteger = false;
                            break;
                          } else {
                            double value = objective[iColumn];
                            double nearest = floor(value + 0.5);
                            if (fabs(value - nearest) > 1.0e-8) {
                              allInteger = false;
                              break;
                            } else {
                              newValues[n] = nearest;
                              newColumn[n++] = iColumn;
                              if (nearest > 0.0) {
                                newLower +=
                                    std::max(columnLower[iColumn], -1.0e20) *
                                    nearest;
                                newUpper +=
                                    std::min(columnUpper[iColumn], 1.0e20) *
                                    nearest;
                              } else {
                                newUpper +=
                                    std::max(columnLower[iColumn], -1.0e20) *
                                    nearest;
                                newLower +=
                                    std::min(columnUpper[iColumn], 1.0e20) *
                                    nearest;
                              }
                            }
                          }
                        }
                      }
                      if (allInteger && n) {
                        fudgeObjective = n;
                        solver2->addCol(0, NULL, NULL, newLower, newUpper, 0.0,
                                        "obj_col");
                        solver2->setInteger(numberColumns);
                        newValues[n] = -1.0;
                        newColumn[n++] = numberColumns;
                        solver2->addRow(n, newColumn, newValues, 0.0, 0.0);
                        if (moveObj) {
                          memset(newValues, 0, numberColumns * sizeof(double));
                          newValues[numberColumns] = 1.0;
                          solver2->setObjective(newValues);
                        }
                        numberRows++;
                        numberColumns++;
                      }
                      delete[] newValues;
                      delete[] newColumn;
                    }
                    if (addSlacks) {
                      bool moveObj = false;
                      addSlacks = 0;
                      // get row copy
                      const CoinPackedMatrix *matrix =
                          solver2->getMatrixByRow();
                      const double *element = matrix->getElements();
                      const int *column = matrix->getIndices();
                      const CoinBigIndex *rowStart = matrix->getVectorStarts();
                      const int *rowLength = matrix->getVectorLengths();
                      const double *rowLower = solver2->getRowLower();
                      const double *rowUpper = solver2->getRowUpper();
                      const double *columnLower = solver2->getColLower();
                      const double *columnUpper = solver2->getColUpper();

                      // maximum space for additional columns
                      CoinBigIndex *newColumnStart =
                          new CoinBigIndex[numberRows + 1];
                      newColumnStart[0] = 0;
                      int *newRow = new int[numberRows];
                      double *newElement = new double[numberRows];
                      double *newObjective = new double[numberRows];
                      double *newColumnLower = new double[numberRows];
                      double *newColumnUpper = new double[numberRows];
                      double *oldObjective = CoinCopyOfArray(
                          solver2->getObjCoefficients(), numberColumns);
                      for (iRow = 0; iRow < numberRows; iRow++) {
                        if (rowLower[iRow] != rowUpper[iRow]) {
                          bool allInteger = true;
                          double newLower = 0.0;
                          double newUpper = 0.0;
                          double constantObjective = 0.0;
                          for (CoinBigIndex j = rowStart[iRow];
                               j < rowStart[iRow] + rowLength[iRow]; j++) {
                            int iColumn = column[j];
                            if (!solver2->isInteger(iColumn)) {
                              allInteger = false;
                              break;
                            } else {
                              double value = element[j];
                              double nearest = floor(value + 0.5);
                              if (fabs(value - nearest) > 1.0e-8) {
                                allInteger = false;
                                break;
                              } else {
                                if (!oldObjective[iColumn])
                                  constantObjective = COIN_DBL_MAX;
                                if (!constantObjective) {
                                  constantObjective =
                                      oldObjective[iColumn] / nearest;
                                } else if (constantObjective != COIN_DBL_MAX) {
                                  double newConstant =
                                      oldObjective[iColumn] / nearest;
                                  if (constantObjective > 0.0) {
                                    if (newConstant <= 0.0)
                                      constantObjective = COIN_DBL_MAX;
                                    else
                                      constantObjective = std::min(
                                          constantObjective, newConstant);
                                  } else {
                                    if (newConstant >= 0.0)
                                      constantObjective = COIN_DBL_MAX;
                                    else
                                      constantObjective = std::max(
                                          constantObjective, newConstant);
                                  }
                                }
                                if (nearest > 0.0) {
                                  newLower +=
                                      std::max(columnLower[iColumn], -1.0e20) *
                                      nearest;
                                  newUpper +=
                                      std::min(columnUpper[iColumn], 1.0e20) *
                                      nearest;
                                } else {
                                  newUpper +=
                                      std::max(columnLower[iColumn], -1.0e20) *
                                      nearest;
                                  newLower +=
                                      std::min(columnUpper[iColumn], 1.0e20) *
                                      nearest;
                                }
                              }
                            }
                          }
                          if (allInteger) {
                            newColumnStart[addSlacks + 1] = addSlacks + 1;
                            newRow[addSlacks] = iRow;
                            newElement[addSlacks] = -1.0;
                            newObjective[addSlacks] = 0.0;
                            if (moveObj && constantObjective != COIN_DBL_MAX) {
                              // move some of objective here if looks constant
                              newObjective[addSlacks] = constantObjective;
                              for (CoinBigIndex j = rowStart[iRow];
                                   j < rowStart[iRow] + rowLength[iRow]; j++) {
                                int iColumn = column[j];
                                double value = element[j];
                                double nearest = floor(value + 0.5);
                                oldObjective[iColumn] -=
                                    nearest * constantObjective;
                              }
                            }
                            newColumnLower[addSlacks] =
                                std::max(newLower, ceil(rowLower[iRow]));
                            ;
                            newColumnUpper[addSlacks] =
                                std::min(newUpper, floor(rowUpper[iRow]));
                            addSlacks++;
                          }
                        }
                      }
                      if (addSlacks) {
                        solver2->setObjective(oldObjective);
                        solver2->addCols(addSlacks, newColumnStart, newRow,
                                         newElement, newColumnLower,
                                         newColumnUpper, newObjective);
                        truncatedRhsLower =
                            CoinCopyOfArray(solver2->getRowLower(), numberRows);
                        truncatedRhsUpper =
                            CoinCopyOfArray(solver2->getRowUpper(), numberRows);
                        for (int j = 0; j < addSlacks; j++) {
                          int iRow = newRow[j];
                          solver2->setRowLower(iRow, 0.0);
                          solver2->setRowUpper(iRow, 0.0);
                          int iColumn = j + numberColumns;
                          solver2->setInteger(iColumn);
                          std::string name = solver2->getRowName(iRow);
                          name += "_int";
                          solver2->setColName(iColumn, name);
                        }
                      }
                    }
                    if (fudgeObjective || addSlacks) {
                      modifiedModel = true;
                      if (fudgeObjective && addSlacks) {
                        buffer.str("");
                        buffer << "Objective integer added with "
                               << fudgeObjective << " elements and "
                               <<  addSlacks << " Integer slacks added";
                      } else if (fudgeObjective) {
                        // just objective
                        buffer.str("");
                        buffer << "Objective integer added with "
                               << fudgeObjective << " elements",
                        more2 &= ~1024;
                      } else {
                        // just slacks
                        buffer.str("");
                        buffer << addSlacks << " Integer slacks added",
                        more2 &= ~512;
                      }
                    } else {
                      more2 &= ~(512 | 1024);
                    }
                    parameters[CbcParam::MOREMOREMIPOPTIONS]->setVal(more2);
                  }
                  if (modifiedModel) {
                    printGeneralMessage(model_, buffer.str());
                    truncateColumns = numberColumns;
                  }
                }
                babModel_->assignSolver(solver2);
                babModel_->setOriginalColumns(process.originalColumns(),
                                              truncateColumns);
                babModel_->initialSolve();
#if CBC_USE_INITIAL_TIME == 2
                // time starts from here?
                time1Elapsed = CoinGetTimeOfDay();
                time1 = CoinCpuTime();
                if (babModel_->useElapsedTime())
                  babModel_->setDblParam(CbcModel::CbcStartSeconds,
                                         CoinGetTimeOfDay());
                else
                  babModel_->setDblParam(CbcModel::CbcStartSeconds,
                                         CoinCpuTime());
                  // babModel_->setMaximumSeconds(timeLeft - (CoinCpuTime() -
                  // time2));
#endif
              }

              if (cgraphMode == "off") {
                // switch off new clique, odd wheel and clique strengthening
                cliqueMode = CbcParameters::CGOff;
                oddWheelMode = CbcParameters::CGOff;
		parameters[CbcParam::ODDWHEELCUTS]->setModeVal(CbcParameters::CGOff);
                clqstrMode = "off";
              } else if (cgraphMode == "clq") {
                // old style
                //CglClique clique;
                //clique.setStarCliqueReport(false);
                //clique.setRowCliqueReport(false);
                //clique.setMinViolation(0.05);
                //int translate[] = {-100, -1, -99, -98, 1, -1098};
                //babModel_->addCutGenerator(&clique, translate[cliqueMode],
		//                         "Clique");
                cliqueMode = CbcParameters::CGOff;
 		parameters[CbcParam::CLIQUECUTS]->setVal("off");
                oddWheelMode = CbcParameters::CGOff;
 		parameters[CbcParam::ODDWHEELCUTS]->setVal("off");
                clqstrMode = "off";
              }
              if (clqstrMode == "after") {
                CglCliqueStrengthening clqStr(babModel_->solver(),
                                              babModel_->messageHandler());
                clqStr.strengthenCliques(4);

                if (clqStr.constraintsExtended() +
                        clqStr.constraintsDominated() >
                    0) {
                  OsiSolverInterface *solver = babModel_->solver();
                  bool takeHint;
                  OsiHintStrength strength;
                  solver->getHintParam(OsiDoDualInResolve, takeHint, strength);
                  solver->setHintParam(OsiDoDualInResolve, false, OsiHintTry);
                  solver->resolve();
                  solver->setHintParam(OsiDoDualInResolve, takeHint, strength);

                  if (solver->isProvenPrimalInfeasible()) {
                     printGeneralMessage(model_,
                                   "Clique Strengthening says infeasible!");
                  } else {
                     buffer.str("");
                     buffer << "After applying Clique Strengthening continuous "
                            << "objective value is "
                            << solver->getObjValue();
                     printGeneralMessage(model_, buffer.str());
                  } // results of clique Strengthening stengthening
                }   // checking impact of clique Strengthening stengthening
              }     // clique Strengthening

              // now tighten bounds
              if (!miplib) {
#ifndef CBC_OTHER_SOLVER
                OsiClpSolverInterface *si =
                    dynamic_cast<OsiClpSolverInterface *>(babModel_->solver());
                assert(si != NULL);
                // get clp itself
                ClpSimplex *modelC = si->getModelPtr();
                // if (noPrinting_)
                // modelC->setLogLevel(0);
                if (!complicatedInteger && modelC->tightenPrimalBounds() != 0) {
                  printGeneralMessage(model_, "Problem is infeasible!");
                  model_.setProblemStatus(0);
                  model_.setSecondaryStatus(1);
                  // say infeasible for solution
                  integerStatus = 6;
                  delete saveSolver;
                  saveSolver = NULL;
                  // and in babModel_ if exists
                  if (babModel_) {
                    babModel_->setProblemStatus(0);
                    babModel_->setSecondaryStatus(1);
                  }
                  break;
                }
                si->resolve();
#elif CBC_OTHER_SOLVER == 1
#endif
              }
              if (debugValues) {
                // for debug
                std::string problemName;
                babModel_->solver()->getStrParam(OsiProbName, problemName);
                babModel_->solver()->activateRowCutDebugger(
                    problemName.c_str());
                twomirGen.probname_ = CoinStrdup(problemName.c_str());
                // checking seems odd
                // redsplitGen.set_given_optsol(babModel_->solver()->getRowCutDebuggerAlways()->optimalSolution(),
                //                         babModel_->getNumCols());
              }
              int testOsiOptions = parameters[CbcParam::TESTOSI]->intVal();
#ifndef JJF_ONE
              // If linked then see if expansion wanted
              {
                OsiSolverLink *solver3 =
                    dynamic_cast<OsiSolverLink *>(babModel_->solver());
                int options = parameters[CbcParam::MIPOPTIONS]->intVal() / 10000;
                if (solver3 || (options & 16) != 0) {
                  if (options) {
                    /*
                                          1 - force mini branch and bound
                                          2 - set priorities high on continuous
                                          4 - try adding OA cuts
                                          8 - try doing quadratic linearization
                                          16 - try expanding knapsacks
                                        */
                    if ((options & 16)) {
                       // FIXME As far as I can see, saveCoinModel is never set anywhere. The only
                       // place it should have been set is in the AMPL parsing, but that block
                       // was never entered as far as I can see. Check whether this block is
                       // needed.
                      int numberColumns = saveCoinModel.numberColumns();
                      int numberRows = saveCoinModel.numberRows();
                      whichColumn = new int[numberColumns];
                      knapsackStart = new int[numberRows + 1];
                      knapsackRow = new int[numberRows];
                      numberKnapsack = 10000;
                      int extra1 = parameters[CbcParam::EXTRA1]->intVal();
                      int extra2 = parameters[CbcParam::EXTRA2]->intVal();
                      int logLevel = parameters[CbcParam::LPLOGLEVEL]->intVal();
                      OsiSolverInterface *solver = expandKnapsack(
                          saveCoinModel, whichColumn, knapsackStart,
                          knapsackRow, numberKnapsack, storedAmpl, logLevel,
                          extra1, extra2, saveTightenedModel);
                      if (solver) {
#ifndef CBC_OTHER_SOLVER
                        clpSolver =
                            dynamic_cast<OsiClpSolverInterface *>(solver);
                        assert(clpSolver);
                        lpSolver = clpSolver->getModelPtr();
#endif
                        babModel_->assignSolver(solver);
                        testOsiOptions = 0;
                        // allow gomory
                        complicatedInteger = 0;
                        // Priorities already done
                        if (info){
                           free(info->priorities);
                           info->priorities = NULL;
                        }
                      } else {
                        numberKnapsack = 0;
                        delete[] whichColumn;
                        delete[] knapsackStart;
                        delete[] knapsackRow;
                        whichColumn = NULL;
                        knapsackStart = NULL;
                        knapsackRow = NULL;
                      }
                    }
                  }
                }
              }
#endif
              if (useCosts && testOsiOptions < 0) {
                int numberColumns = babModel_->getNumCols();
                int *sort = new int[numberColumns];
                double *dsort = new double[numberColumns];
                int *priority = new int[numberColumns];
                const double *objective = babModel_->getObjCoefficients();
                const double *lower = babModel_->getColLower();
                const double *upper = babModel_->getColUpper();
                const CoinPackedMatrix *matrix =
                    babModel_->solver()->getMatrixByCol();
                const int *columnLength = matrix->getVectorLengths();
                int iColumn;
                int n = 0;
                for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                  if (babModel_->isInteger(iColumn)) {
                    sort[n] = n;
                    if (useCosts == 1) {
                      dsort[n++] = -fabs(objective[iColumn]);
                    } else if (useCosts == 2) {
                      dsort[n++] = iColumn;
                    } else if (useCosts == 3) {
                      dsort[n++] = (upper[iColumn]-lower[iColumn])==1.0 ?1:2;
                    } else if (useCosts == 4) {
                      dsort[n++] = (upper[iColumn]-lower[iColumn])!=1.0 ?1:2;
                    } else if (useCosts == 5) {
                      dsort[n++] = -columnLength[iColumn];
                    } else if (useCosts == 6) {
                      dsort[n++] = (columnLength[iColumn] == 1) ? -1.0 : 0.0;
                    } else if (useCosts == 7) {
                      dsort[n++] = (objective[iColumn]) ? -1.0 : 0.0;
		    }
                  }
                }
                CoinSort_2(dsort, dsort + n, sort);
                int level = 999;
                double last = -1.0e100;
                for (int i = 0; i < n; i++) {
                  int iPut = sort[i];
                  if (dsort[i] != last) {
                    level++;
                    last = dsort[i];
                  }
                  priority[iPut] = level;
                }
                if (newPriorities) {
                  // get rid of
                  delete[] newPriorities;
                  newPriorities = NULL;
                }
                babModel_->passInPriorities(priority, false);
                integersOK = true;
                delete[] priority;
                delete[] sort;
                delete[] dsort;
              }
              // Set up heuristics
              doHeuristics(babModel_, ((!miplib) ? 1 : 10), parameters,
                           parameters.noPrinting(), initialPumpTune);
              if (!miplib) {
                if (parameters[CbcParam::LOCALTREE]->modeVal()) {
                  CbcTreeLocal localTree(babModel_, NULL, 10, 0, 0, 10000,
                                         2000);
                  babModel_->passInTreeHandler(localTree);
                }
              }
              if (cbcParamCode == CbcParam::MIPLIB) {
                if (babModel_->numberStrong() == 5 &&
                    babModel_->numberBeforeTrust() == 5)
                  babModel_->setNumberBeforeTrust(10);
              }
              int experimentFlag = parameters[CbcParam::EXPERIMENT]->intVal();
              int strategyFlag = parameters[CbcParam::STRATEGY]->modeVal();
              int bothFlags = std::max(std::min(experimentFlag, 1), strategyFlag);
              // add cut generators if wanted
              int switches[30] = {};
              int accuracyFlag[30] = {};
              char doAtEnd[30] = {};
	      int lagrangeanFlag =
		(parameters[CbcParam::MOREMOREMIPOPTIONS]->intVal()&(7*33554432))>>9;
#define ALL_LAGRANGEAN 1
              int numberGenerators = 0;
              std::map<int, int> translate;
              translate[CbcParameters::CGOff] = -100;
              translate[CbcParameters::CGOn] = -1;
              translate[CbcParameters::CGRoot] = -99;
              translate[CbcParameters::CGIfMove] = -98;
              translate[CbcParameters::CGForceOn] = 1;
              translate[CbcParameters::CGOnGlobal] = -1098;
              translate[CbcParameters::CGForceOnGlobal] = -999;
              translate[CbcParameters::CGForceOnBut] = 1;
              translate[CbcParameters::CGForceOnStrong] = 1;
              translate[CbcParameters::CGForceOnButStrong] = 1;
              translate[CbcParameters::CGStrongRoot] = -1;
              std::map<int, int> laTranslate;
	      laTranslate[CbcParameters::CGEndOnlyRoot] = 1;
	      laTranslate[CbcParameters::CGEndCleanRoot] = 2;
	      laTranslate[CbcParameters::CGEndBothRoot] = 3;
	      laTranslate[CbcParameters::CGEndOnly] = 4;
	      laTranslate[CbcParameters::CGEndClean] = 5;
	      laTranslate[CbcParameters::CGEndBoth] = 6;
	      laTranslate[CbcParameters::CGOnlyAsWell] = 7;
	      laTranslate[CbcParameters::CGOnlyAsWellRoot] = 13;
	      laTranslate[CbcParameters::CGCleanAsWell] = 8;
	      laTranslate[CbcParameters::CGCleanAsWellRoot] =14;
	      laTranslate[CbcParameters::CGBothAsWell] = 9;
	      laTranslate[CbcParameters::CGBothAsWellRoot] = 15;
	      laTranslate[CbcParameters::CGOnlyInstead] = 10;
	      laTranslate[CbcParameters::CGCleanInstead] = 11;
	      laTranslate[CbcParameters::CGBothInstead] = 12;
              int maximumSlowPasses = parameters[CbcParam::MAXSLOWCUTS]->intVal();
 	     assert (parameters[CbcParam::PROBINGCUTS]->modeVal()==probingMode);
              if (probingMode) {
                int numberColumns = babModel_->solver()->getNumCols();
                if (probingMode > CbcParameters::CGForceOnBut) {
                  probingGen.setMaxElements(numberColumns);
                  probingGen.setMaxElementsRoot(numberColumns);
                }
                probingGen.setMaxProbeRoot(std::min(2000, numberColumns));
                probingGen.setMaxProbeRoot(123);
                probingGen.setMaxProbe(123);
                probingGen.setMaxLookRoot(20);
                if (probingMode == CbcParameters::CGForceOnBut ||
		    probingMode == CbcParameters::CGForceOnButStrong)
                  probingGen.setRowCuts(-3); // strengthening etc just at root
                if (probingMode == CbcParameters::CGForceOnStrong ||
		    probingMode == CbcParameters::CGForceOnButStrong) {
                  // Number of unsatisfied variables to look at
                  probingGen.setMaxProbeRoot(numberColumns);
                  probingGen.setMaxProbe(numberColumns);
                  // How far to follow the consequences
                  probingGen.setMaxLook(50);
                  probingGen.setMaxLookRoot(50);
                }
                if (probingMode == CbcParameters::CGStrongRoot) {
                  probingGen.setMaxPassRoot(2);
                  probingGen.setMaxProbeRoot(numberColumns);
                  probingGen.setMaxLookRoot(numberColumns);
                }
                // If 5 then force on
                int iMode = translate[probingMode];
                if (probingMode == CbcParameters::CGOnGlobal)
                  iMode = 1;
                babModel_->addCutGenerator(&probingGen, iMode, "Probing");
                accuracyFlag[numberGenerators] = 5;
                switches[numberGenerators++] = 0;
              }
 	      assert (parameters[CbcParam::GOMORYCUTS]->modeVal()==gomoryMode);
              if (gomoryMode && (complicatedInteger != 1 ||
                                   (gomoryMode == 1 || gomoryMode >= 4))) {
                // try larger limit
                int numberColumns = babModel_->getNumCols();
                if (gomoryMode == CbcParameters::CGForceOnBut) {
                  gomoryMode = CbcParameters::CGForceOn;
                  gomoryGen.setLimitAtRoot(numberColumns);
                  gomoryGen.setLimit(numberColumns);
                } else if (gomoryMode == CbcParameters::CGForceOnStrong) {
                  gomoryMode = CbcParameters::CGIfMove;
                  gomoryGen.setLimitAtRoot(numberColumns);
                  gomoryGen.setLimit(200);
                } else if (gomoryMode == CbcParameters::CGForceOnButStrong) {
                  gomoryMode = CbcParameters::CGIfMove;
                  gomoryGen.setLimitAtRoot(500);
                  gomoryGen.setLimit(200);
                } else if (numberColumns > 5000) {
                  //#define MORE_CUTS2
#ifdef MORE_CUTS2
                  // try larger limit
                  gomoryGen.setLimitAtRoot(numberColumns);
                  gomoryGen.setLimit(200);
#else
                  gomoryGen.setLimitAtRoot(2000);
                  // gomoryGen.setLimit(200);
#endif
                } else {
#ifdef MORE_CUTS2
                  // try larger limit
                  gomoryGen.setLimitAtRoot(numberColumns);
                  gomoryGen.setLimit(200);
#endif
                }
                int cutLength = parameters[CbcParam::CUTLENGTH]->intVal();
                if (cutLength != -1) {
                  gomoryGen.setLimitAtRoot(cutLength);
                  if (cutLength < 10000000) {
                    gomoryGen.setLimit(cutLength);
                  } else {
                    gomoryGen.setLimit(cutLength % 10000000);
                  }
                }
                int laGomory = parameters[CbcParam::LAGOMORYCUTS]->modeVal();
                int gType = translate[gomoryMode];
                if (!laGomory) {
                  // Normal
                  babModel_->addCutGenerator(&gomoryGen,
                                             translate[gomoryMode], "Gomory");
                  accuracyFlag[numberGenerators] = 3;
                  switches[numberGenerators++] = lagrangeanFlag;
                } else {
		  laGomory = laTranslate[laGomory]-1;
                  int type = (laGomory % 3) + 1;
                  int when = laGomory / 3;
                  char atEnd = (when < 2) ? 1 : 0;
                  int gomoryTypeMajor = 10;
                  if (when != 3) {
                    // normal as well
                    babModel_->addCutGenerator(&gomoryGen, gType, "Gomory");
                    accuracyFlag[numberGenerators] = 3;
                    switches[numberGenerators++] = 0;
                    if (when == 2) {
                      gomoryTypeMajor = 20;
                    } else if (when == 4) {
                      gomoryTypeMajor = 20;
                      when = 0;
                    }
                  } else {
                    when--; // so on
                    gomoryTypeMajor = 20;
                  }
                  if (!when)
                    gType = -99; // root
                  gomoryGen.passInOriginalSolver(babModel_->solver());
                  if ((type & 1) != 0) {
                    // clean
                    gomoryGen.setGomoryType(gomoryTypeMajor + 1);
                    babModel_->addCutGenerator(&gomoryGen, gType, "GomoryL1");
                    accuracyFlag[numberGenerators] = 3;
                    doAtEnd[numberGenerators] = atEnd;
                    if (atEnd) {
                      babModel_->cutGenerator(numberGenerators)
                          ->setMaximumTries(99999999);
                      babModel_->cutGenerator(numberGenerators)->setHowOften(1);
                    }
                    switches[numberGenerators++] = 16384;
                  }
                  if ((type & 2) != 0) {
                    // simple
                    gomoryGen.setGomoryType(gomoryTypeMajor + 2);
                    babModel_->addCutGenerator(&gomoryGen, gType, "GomoryL2");
                    accuracyFlag[numberGenerators] = 3;
                    doAtEnd[numberGenerators] = atEnd;
                    if (atEnd) {
                      babModel_->cutGenerator(numberGenerators)
                          ->setMaximumTries(99999999);
                      babModel_->cutGenerator(numberGenerators)->setHowOften(1);
                    }
                    switches[numberGenerators++] = 32768;
                  }
                }
              }
#ifdef CLIQUE_ANALYSIS
              if (miplib && !storedAmpl.sizeRowCuts()) {
                printf("looking at probing\n");
                babModel_->addCutGenerator(&storedAmpl, 1, "Stored");
                accuracyFlag[numberGenerators] = 0;
                switches[numberGenerators++] = 0;
              }
#endif
 	      assert (parameters[CbcParam::KNAPSACKCUTS]->modeVal()==knapsackMode);
              if (knapsackMode) {
                babModel_->addCutGenerator(
                    &knapsackGen, translate[knapsackMode], "Knapsack");
                accuracyFlag[numberGenerators] = 1;
                switches[numberGenerators++] = -2;
              }
              if (redsplitMode && !complicatedInteger) {
                babModel_->addCutGenerator(&redsplitGen,
                                           translate[redsplitMode],
                                           "Reduce-and-split");
                accuracyFlag[numberGenerators] = 5;
                // slow ? - just do a few times
                if (redsplitMode != CbcParameters::CGOn) {
                  babModel_->cutGenerator(numberGenerators)
                      ->setMaximumTries(maximumSlowPasses);
                  babModel_->cutGenerator(numberGenerators)->setHowOften(10);
                }
                switches[numberGenerators++] = 1 | (ALL_LAGRANGEAN*lagrangeanFlag); 
              }
              if (redsplit2Mode && !complicatedInteger) {
                int maxLength = 256;
                if (redsplit2Mode > CbcParameters::CGRoot) {
                  redsplit2Mode -= 2;
                  maxLength = COIN_INT_MAX;
                }
                CglRedSplit2Param &parameters = redsplit2Gen.getParam();
                parameters.setMaxNonzeroesTab(maxLength);
                babModel_->addCutGenerator(&redsplit2Gen,
                                           translate[redsplit2Mode],
                                           "Reduce-and-split(2)");
                accuracyFlag[numberGenerators] = 5;
                // slow ? - just do a few times
                if (redsplit2Mode != CbcParameters::CGOn) {
                  babModel_->cutGenerator(numberGenerators)
                      ->setHowOften(maximumSlowPasses);
                  babModel_->cutGenerator(numberGenerators)
                      ->setMaximumTries(maximumSlowPasses);
                  babModel_->cutGenerator(numberGenerators)->setHowOften(5);
                }

                switches[numberGenerators++] = 1 | (ALL_LAGRANGEAN*lagrangeanFlag);
              }
              if (GMIMode && !complicatedInteger) {
                if (GMIMode > CbcParameters::CGOnGlobal) {
                  // long
                  GMIMode -= 5;
                  CglGMIParam &parameters = GMIGen.getParam();
                  parameters.setMaxSupportRel(1.0);
                }
                babModel_->addCutGenerator(&GMIGen, translate[GMIMode],
                                           "Gomory(2)");
                if (GMIMode == CbcParameters::CGOnGlobal) {
                  // just at end and root
                  GMIMode = CbcParameters::CGRoot;
                  doAtEnd[numberGenerators] = 1;
                  babModel_->cutGenerator(numberGenerators)
                      ->setMaximumTries(99999999);
                  babModel_->cutGenerator(numberGenerators)->setHowOften(1);
                }
                accuracyFlag[numberGenerators] = 5;
                switches[numberGenerators++] = 0 | (ALL_LAGRANGEAN*lagrangeanFlag);
              }
              if (cliqueMode) {
                bkCliqueGen.setMaxCallsBK(maxCallsBK);
                bkCliqueGen.setExtendingMethod(bkClqExtMethod);
                bkCliqueGen.setPivotingStrategy(bkPivotingStrategy);
                babModel_->addCutGenerator(&bkCliqueGen,
                                           translate[cliqueMode], "Clique");
                accuracyFlag[numberGenerators] = 0;
                switches[numberGenerators++] = 0;
	      } else if (cgraphMode == "clq") {
		// old style
		CglClique clique;
		clique.setStarCliqueReport(false);
		clique.setRowCliqueReport(false);
		clique.setMinViolation(0.05);
		oddWheelMode = 0;
		parameters[CbcParam::ODDWHEELCUTS]->setModeVal(CbcParameters::CGOff);
		babModel_->addCutGenerator(&clique, translate[oldCliqueMode],
					   "Clique");
                accuracyFlag[numberGenerators] = 0;
                switches[numberGenerators++] = 0;
              }
 	      assert (parameters[CbcParam::ODDWHEELCUTS]->modeVal()==oddWheelMode);
              if (oddWheelMode) {
                oddWheelGen.setExtendingMethod(oddWExtMethod);
                babModel_->addCutGenerator(
                    &oddWheelGen, translate[oddWheelMode], "OddWheel");
                accuracyFlag[numberGenerators] = 0;
                switches[numberGenerators++] = 0;
              }
 	      assert (parameters[CbcParam::MIRCUTS]->modeVal()==mixedMode);
              if (mixedMode) {
		if (mixedRoundStrategy != 1)
		  mixedGen.setMAXAGGR_(mixedRoundStrategy);
                babModel_->addCutGenerator(&mixedGen, translate[mixedMode],
                                           "MixedIntegerRounding2");
                accuracyFlag[numberGenerators] = 2;
                switches[numberGenerators++] = 0 | (ALL_LAGRANGEAN*lagrangeanFlag);
              }
 	      assert (parameters[CbcParam::FLOWCUTS]->modeVal()==flowMode);
              if (flowMode) {
                babModel_->addCutGenerator(&flowGen, translate[flowMode],
                                           "FlowCover");
                accuracyFlag[numberGenerators] = 2;
                switches[numberGenerators++] = 0 | (ALL_LAGRANGEAN*lagrangeanFlag);
              }
              if (twomirMode && (complicatedInteger != 1 ||
                                   (twomirMode == CbcParameters::CGOn ||
				    twomirMode >= CbcParameters::CGForceOn))) {
                // try larger limit
                int numberColumns = babModel_->getNumCols();
                if (twomirMode == CbcParameters::CGForceOnBut) {
                  twomirMode = CbcParameters::CGForceOn;
                  twomirGen.setMaxElements(numberColumns);
                } else if (numberColumns > 5000 && twomirMode == CbcParameters::CGForceOn) {
                  twomirGen.setMaxElements(2000);
                }
                int laTwomir = parameters[CbcParam::LATWOMIRCUTS]->modeVal();
                int twomirType = translate[twomirMode];
                if (!laTwomir) {
                  // Normal
                  babModel_->addCutGenerator(
                      &twomirGen, translate[twomirMode], "TwoMirCuts");
                  accuracyFlag[numberGenerators] = 4;
                  switches[numberGenerators++] = 1 | lagrangeanFlag;
                } else {
		  laTwomir = laTranslate[laTwomir]-1;
                  int type = (laTwomir % 3) + 1;
                  int when = laTwomir / 3;
                  char atEnd = (when < 2) ? 1 : 0;
                  int twomirTypeMajor = 10;
                  if (when < 3) {
                    // normal as well
                    babModel_->addCutGenerator(
                        &twomirGen, translate[twomirMode], "TwoMirCuts");
                    accuracyFlag[numberGenerators] = 4;
                    switches[numberGenerators++] = 1;
                    if (when == 2)
                      twomirTypeMajor = 10;
                  } else {
                    when--; // so on
                    twomirTypeMajor = 20;
                  }
                  if (!when)
                    twomirType = -99; // root
                  twomirGen.passInOriginalSolver(babModel_->solver());
                  if ((type & 1) != 0) {
                    // clean
                    twomirGen.setTwomirType(twomirTypeMajor + 1);
                    babModel_->addCutGenerator(&twomirGen, twomirType,
                                               "TwoMirCutsL1");
                    accuracyFlag[numberGenerators] = 4;
                    doAtEnd[numberGenerators] = atEnd;
                    switches[numberGenerators++] = (atEnd ? 0 : 1) | 16384;
                  }
                  if ((type & 2) != 0) {
                    // simple
                    twomirGen.setTwomirType(twomirTypeMajor + 2);
                    babModel_->addCutGenerator(&twomirGen, twomirType,
                                               "TwoMirCutsL2");
                    accuracyFlag[numberGenerators] = 4;
                    doAtEnd[numberGenerators] = atEnd;
                    switches[numberGenerators++] = (atEnd ? 0 : 1) | 32768;
                  }
                }
              }
#ifndef DEBUG_MALLOC
 	      assert (parameters[CbcParam::LANDPCUTS]->modeVal()==landpMode);
              if (landpMode) {
                if (landpMode == CbcParameters::CGOnGlobal) {
                  // allow longer
                  landpGen.parameter().maximumCutLength = 2000000;
                  landpMode = CbcParameters::CGIfMove;
                }
                babModel_->addCutGenerator(&landpGen, translate[landpMode],
                                           "LiftAndProject");
                accuracyFlag[numberGenerators] = 5;
                // slow ? - just do a few times
                if (landpMode != CbcParameters::CGOn) {
                  babModel_->cutGenerator(numberGenerators)
                      ->setMaximumTries(maximumSlowPasses);
                  babModel_->cutGenerator(numberGenerators)->setHowOften(10);
                }
                switches[numberGenerators++] = 1 | (ALL_LAGRANGEAN*lagrangeanFlag);
              }
#endif
 	      assert (parameters[CbcParam::RESIDCAPCUTS]->modeVal()==residualCapacityMode);
              if (residualCapacityMode) {
                babModel_->addCutGenerator(&residualCapacityGen,
                                           translate[residualCapacityMode],
                                           "ResidualCapacity");
                accuracyFlag[numberGenerators] = 5;
                switches[numberGenerators++] = 1 | (ALL_LAGRANGEAN*lagrangeanFlag);
              }
              if (zerohalfMode) {
                if (zerohalfMode > CbcParameters::CGForceOn) {
                  // zerohalfMode -=4;
                  zerohalfGen.setFlags(1);
                }
                babModel_->addCutGenerator(
                    &zerohalfGen, translate[zerohalfMode], "ZeroHalf");
                accuracyFlag[numberGenerators] = 5;
                babModel_->cutGenerator(numberGenerators)
                    ->setNeedsRefresh(true);
                switches[numberGenerators++] = 2; //| (ALL_LAGRANGEAN*lagrangeanFlag);
              }
              if (dominatedCuts)
                babModel_->setSpecialOptions(babModel_->specialOptions() | 64);
              // Say we want timings
              numberGenerators = babModel_->numberCutGenerators();
              int iGenerator;
              int cutDepth = parameters[CbcParam::CUTDEPTH]->intVal();
              for (iGenerator = 0; iGenerator < numberGenerators;
                   iGenerator++) {
                CbcCutGenerator *generator =
                    babModel_->cutGenerator(iGenerator);
                int howOften = generator->howOften();
		int iSwitch = switches[iGenerator];
		int iSwitch2, iSwitch3;
		if (iSwitch>=0) {
		  iSwitch2 = iSwitch&127;
		  iSwitch3 = iSwitch & ~16383;
		  // say if wants modified solver
		  generator->setSwitches(generator->switches()|iSwitch3);
		} else {
		  iSwitch2 = iSwitch;
		}
                if (howOften == -98 || howOften == -99 || generator->maximumTries() > 0)
                  generator->setSwitchOffIfLessThan(iSwitch2);
                // Use if any at root as more likely later and fairly cheap
                // if (switches[iGenerator]==-2)
                // generator->setWhetherToUse(true);
                generator->setInaccuracy(accuracyFlag[iGenerator]);
                if (doAtEnd[iGenerator]) {
                  generator->setWhetherCallAtEnd(true);
                  // generator->setMustCallAgain(true);
                }
                generator->setTiming(true);
                if (cutDepth >= 0)
                  generator->setWhatDepth(cutDepth);
              }
              // Could tune more
              if (!miplib) {
                double minimumDrop =
                    fabs(babModel_->solver()->getObjValue()) * 1.0e-5 + 1.0e-5;
                babModel_->setMinimumDrop(std::min(5.0e-2, minimumDrop));
                if (cutPass == -1234567) {
                  if (babModel_->getNumCols() < 500)
                    babModel_->setMaximumCutPassesAtRoot(
                        -100); // always do 100 if possible
                  else if (babModel_->getNumCols() < 5000)
                    babModel_->setMaximumCutPassesAtRoot(
                        100); // use minimum drop
                  else
                    babModel_->setMaximumCutPassesAtRoot(50);
                } else {
                  babModel_->setMaximumCutPassesAtRoot(cutPass);
                }
                if (cutPassInTree == -1234567)
                  babModel_->setMaximumCutPasses(4);
                else
                  babModel_->setMaximumCutPasses(cutPassInTree);
              } else if (cutPass != -1234567) {
                babModel_->setMaximumCutPassesAtRoot(cutPass);
              }
              // Do more strong branching if small
              // if (babModel_->getNumCols()<5000)
              // babModel_->setNumberStrong(20);
              // Switch off strong branching if wanted
              // if (babModel_->getNumCols()>10*babModel_->getNumRows())
              // babModel_->setNumberStrong(0);
              if (!parameters.noPrinting()) {
                int iLevel = parameters[CbcParam::LOGLEVEL]->intVal();
                if (iLevel < 0) {
                  if (iLevel > -10) {
                    babModel_->setPrintingMode(1);
                  } else {
                    babModel_->setPrintingMode(2);
                    iLevel += 10;
                    parameters[CbcParam::LOGLEVEL]->setVal(iLevel);
                  }
                  iLevel = -iLevel;
                }
                babModel_->messageHandler()->setLogLevel(iLevel);
                // if (babModel_->getNumCols() > 2000 || babModel_->getNumRows()
                // > 1500 || babModel_->messageHandler()->logLevel() > 1)
                //  babModel_->setPrintFrequency(100);
                babModel_->setPrintFrequency(1);
              }

              babModel_->solver()->setIntParam(
                  OsiMaxNumIterationHotStart,
                  parameters[CbcParam::MAXHOTITS]->intVal());
#ifndef CBC_OTHER_SOLVER
              OsiClpSolverInterface *osiclp =
                  dynamic_cast<OsiClpSolverInterface *>(babModel_->solver());
              // go faster stripes
              if ((osiclp->getNumRows() < 300 && osiclp->getNumCols() < 500)) {
                osiclp->setupForRepeatedUse(2, parameters[CbcParam::LPLOGLEVEL]->intVal());
                if (bothFlags >= 1) {
                  ClpSimplex *lp = osiclp->getModelPtr();
                  int specialOptions = lp->specialOptions();
                  lp->setSpecialOptions(specialOptions | (2048 + 4096));
                }
              } else {
                osiclp->setupForRepeatedUse(0, parameters[CbcParam::LPLOGLEVEL]->intVal());
              }
              if (bothFlags >= 2) {
                ClpSimplex *lp = osiclp->getModelPtr();
                int specialOptions = lp->specialOptions();
                lp->setSpecialOptions(specialOptions | (2048 + 4096));
              }
              double increment = babModel_->getCutoffIncrement();
              ;
              int *changed = NULL;
#ifdef UNSAFE_FOR_LAZY_CUTS
              if (!miplib && increment == normalIncrement) 
                changed = analyze(osiclp, numberChanged, increment, false,
                                  generalMessageHandler, parameters.noPrinting());
#endif
#elif CBC_OTHER_SOLVER == 1
              double increment = babModel_->getCutoffIncrement();
              ;
#endif
              if (debugValues) {
                int numberColumns = babModel_->solver()->getNumCols();
                if (numberDebugValues == numberColumns) {
                  // for debug
                  babModel_->solver()->activateRowCutDebugger(debugValues);
                } else {
                  int numberOriginalColumns = process.originalModel()->getNumCols();
#if DEBUG_PREPROCESS < 2
                  if (numberDebugValues <= numberOriginalColumns) {
                    const int *originalColumns = process.originalColumns();
                    double *newValues = new double[numberColumns];
                    // in case preprocess added columns!
                    // need to find values
                    OsiSolverInterface *siCopy = babModel_->solver()->clone();
                    for (int i = 0; i < numberColumns; i++) {
                      int jColumn = originalColumns[i];
                      if (jColumn < numberDebugValues && siCopy->isInteger(i)) {
                        // integer variable
                        double soln = floor(debugValues[jColumn] + 0.5);
                        // Set bounds to fix variable to its solution
                        siCopy->setColUpper(i, soln);
                        siCopy->setColLower(i, soln);
                      }
                    }
                    // All integers have been fixed at optimal value.
                    // Now solve to get continuous values
                    siCopy->setHintParam(OsiDoScale, false);
                    siCopy->initialSolve();
                    if (siCopy->isProvenOptimal()) {
                      memcpy(newValues, siCopy->getColSolution(),
                             numberColumns * sizeof(double));
                    } else {
                      printf("BAD debug file\n");
                      siCopy->writeMps("Bad");
                      exit(22);
                    }
                    delete siCopy;
                    // for debug
                    babModel_->solver()->activateRowCutDebugger(newValues);
                    delete[] newValues;
                  } else {
                    printf("debug file has incorrect number of columns\n");
                  }
#else
		  // for debug
		  babModel_->solver()->activateRowCutDebugger(debugSolution);
#endif
                }
              }
              babModel_->setCutoffIncrement(
                  std::max(babModel_->getCutoffIncrement(), increment));
              // Turn this off if you get problems
              // Used to be automatically set
              int mipOptions = parameters[CbcParam::MIPOPTIONS]->intVal() % 10000;
              if (mipOptions != (1057) && mipOptions != 1025) {
                buffer.str("");
                buffer << "mip options " << mipOptions;
                printGeneralMessage(model_, buffer.str());
              }
#ifndef CBC_OTHER_SOLVER
              osiclp->setSpecialOptions(mipOptions);
#elif CBC_OTHER_SOLVER == 1
#endif
              // probably faster to use a basis to get integer solutions
              babModel_->setSpecialOptions(babModel_->specialOptions() | 2);
              currentBranchModel = babModel_;
              // OsiSolverInterface * strengthenedModel=NULL;
              if (cbcParamCode == CbcParam::BAB ||
                  cbcParamCode == CbcParam::MIPLIB) {
                if (strategyFlag == 1) {
                  // try reduced model
                  babModel_->setSpecialOptions(babModel_->specialOptions() |
                                               512);
                }
                if (experimentFlag >= 5 || strategyFlag == 2) {
                  // try reduced model at root
                  babModel_->setSpecialOptions(babModel_->specialOptions() |
                                               32768);
                  if (experimentFlag >= 10000) {
                    // try reduced model at root with cuts
                    babModel_->setSpecialOptions(babModel_->specialOptions() |
                                                 512);
                  }
                }
                {
                  int depthMiniBab = parameters[CbcParam::DEPTHMINIBAB]->intVal();
                  if (depthMiniBab != -1)
                    babModel_->setFastNodeDepth(depthMiniBab);
                }
                int extra4 = parameters[CbcParam::EXTRA4]->intVal();
                if (extra4 >= 0) {
                  int strategy = extra4 % 10;
                  extra4 /= 10;
                  int method = extra4 % 100;
                  extra4 /= 100;
                  extra4 = strategy + method * 8 + extra4 * 1024;
                  babModel_->setMoreSpecialOptions(extra4);
                }
                int moreMipOptions = parameters[CbcParam::MOREMIPOPTIONS]->intVal();
                if (moreMipOptions >= 0) {
                  buffer.str("");
                  buffer << "more mip options " << moreMipOptions;
                  printGeneralMessage(model_, buffer.str());
#if 1
                  // some options may have been set already
                  // e.g. use elapsed time
                  babModel_->setMoreSpecialOptions(
                      moreMipOptions | babModel_->moreSpecialOptions());
#else
                  OsiClpSolverInterface *osiclp =
                      dynamic_cast<OsiClpSolverInterface *>(
                          babModel_->solver());
                  if (moreMipOptions == 10000) {
                    // test memory saving
                    moreMipOptions -= 10000;
                    ClpSimplex *lpSolver = osiclp->getModelPtr();
                    lpSolver->setPersistenceFlag(1);
                    // switch off row copy if few rows
                    if (lpSolver->numberRows() < 150)
                      lpSolver->setSpecialOptions(lpSolver->specialOptions() |
                                                  256);
                  }
                  if (moreMipOptions < 10000 && moreMipOptions) {
                    if (((moreMipOptions + 1) % 1000000) != 0)
                      babModel_->setSearchStrategy(moreMipOptions % 1000000);
                  } else if (moreMipOptions < 100000) {
                    // try reduced model
                    babModel_->setSpecialOptions(babModel_->specialOptions() |
                                                 512);
                  }
                  // go faster stripes
                  if (moreMipOptions >= 999999) {
                    if (osiclp) {
                      int save = osiclp->specialOptions();
                      osiclp->setupForRepeatedUse(2, 0);
                      osiclp->setSpecialOptions(save |
                                                osiclp->specialOptions());
                    }
                  }
#endif
                }
              }
              {
                int extra1 = parameters[CbcParam::EXTRA1]->intVal();
                if (extra1 != -1) {
                  if (extra1 < 0) {
                    if (extra1 == -7777)
                      extra1 = -1;
                    else if (extra1 == -8888)
                      extra1 = 1;
                    babModel_->setWhenCuts(-extra1);
                  } else if (extra1 < 19000) {
                    babModel_->setSearchStrategy(extra1);
                    printf("XXXXX searchStrategy %d\n", extra1);
                  } else {
                    int n = extra1 - 20000;
                    if (!n)
                      n--;
                    babModel_->setNumberAnalyzeIterations(n);
                    printf("XXXXX analyze %d\n", extra1);
                  }
                } else if (bothFlags >= 1) {
                  babModel_->setWhenCuts(999998);
                }
              }
              if (cbcParamCode == CbcParam::BAB) {
                if (info && statusUserFunction_[0]) {
                  priorities = info->priorities;
                  branchDirection = info->branchDirection;
                  pseudoDown = info->pseudoDown;
                  pseudoUp = info->pseudoUp;
                  solutionIn = info->primalSolution;
                  prioritiesIn = info->priorities;
                  if (info->numberSos && doSOS) {
                    // SOS
                    numberSOS = info->numberSos;
                    sosStart = info->sosStart;
                    sosIndices = info->sosIndices;
                    sosType = info->sosType;
                    sosReference = info->sosReference;
                    sosPriority = info->sosPriority;
                  }
                }
                const int *originalColumns =
                    preProcess ? process.originalColumns() : NULL;
                if (model.getMIPStart().size())
                  mipStart = model.getMIPStart();
                std::string testEmpty = mipStartFile.substr(0, 6);
                if ((mipStart.size() || testEmpty == "empty.") &&
                    !mipStartBefore.size() && babModel_->getNumCols()) {
                  std::vector<std::string> colNames;
                  if (preProcess) {
                    /* translating mipstart solution */
                    std::map<std::string, double> mipStartV;
                    for (size_t i = 0; (i < mipStart.size()); ++i)
                      mipStartV[mipStart[i].first] = mipStart[i].second;

                    std::vector<std::pair<std::string, double> > mipStart2;
                    for (int i = 0; (i < babModel_->solver()->getNumCols());
                         ++i) {
                      int iColumn = babModel_->originalColumns()[i];
                      if (iColumn >= 0 && iColumn < model.getNumCols()) {
                        std::string cname =
                            model_.solver()->getColName(iColumn);
                        colNames.push_back(cname);
                        babModel_->solver()->setColName(i, cname);
                        std::map<std::string, double>::const_iterator msIt =
                            mipStartV.find(cname);
                        if (msIt != mipStartV.end())
                          mipStart2.push_back(std::pair<std::string, double>(
                              cname, msIt->second));
                      } else {
                        // created variable
                        char newName[15];
                        sprintf(newName, "C%7.7d", i);
                        colNames.push_back(newName);
                      }
                    }
                    mipStart = mipStart2;
                  } else {
                    for (int i = 0; (i < babModel_->solver()->getNumCols());
                         ++i)
                      colNames.push_back(model_.solver()->getColName(i));
                  }
                  // printf("--- %s %d\n",
                  // babModel_->solver()->getColName(0).c_str(),
                  // babModel_->solver()->getColNames().size() ); printf("--
                  // SIZES of models %d %d %d\n", model_.getNumCols(),
                  // babModel_->solver()->getNumCols(),
                  // babModel_->solver()->getColNames().size() );
                  std::vector<double> x(babModel_->getNumCols(), 0.0);
                  double obj;
                  babModel_->findIntegers(true);
                  int extraActions = 0;
                  int lengthFileName = mipStartFile.size();
                  const char *checkEnd[6] = {".low",          ".high",
                                             ".lowcheap",     ".highcheap",
                                             ".lowexpensive", ".highexpensive"};
                  for (extraActions = 0; extraActions < 6; extraActions++) {
                    if (ends_with(mipStartFile,
                                  std::string(checkEnd[extraActions])))
                      break;
                  }
                  if (extraActions == 6)
                    extraActions = 0;
                  else
                    extraActions++;
                  int status = CbcMipStartIO::computeCompleteSolution(
                      babModel_, babModel_->solver(), colNames, mipStart, &x[0],
                      obj, extraActions, babModel_->messageHandler(),
                      babModel_->messagesPointer());
                  if (!status) {
                    // need to check more babModel_->setBestSolution( &x[0],
                    // static_cast<int>(x.size()), obj, false );
                    OsiBabSolver dummy;
                    babModel_->passInSolverCharacteristics(&dummy);
                    babModel_->createContinuousSolver();
		    /* This is a mess - mipStart does not really know about
		       SOS (but I am not going to try and fix that at present).
		       But at least we should make sure solution is valid */
		    int numberObjects = babModel_->numberObjects();
		    OsiObject **oldObjects = babModel_->objects();
		    int numberSOS = 0;
		    if (preProcess  && redoSOS) {
		      for (int iObj = 0; iObj < numberObjects; iObj++) {
			int iColumn = oldObjects[iObj]->columnNumber();
			if (iColumn < 0 || iColumn >= numberOriginalColumns) {
                          CbcSOS *objSOS = dynamic_cast< CbcSOS * >(oldObjects[iObj]);
                          if (objSOS) {
			    numberSOS++;
			  } else {
			    printf("not sure mipstart can deal with this\n");
			    abort();
			  }
			}
		      }
		    }
		    if (!numberSOS) {
		      // normal
		      babModel_->setBestSolution(CBC_ROUNDING,
						 obj, &x[0], 1);
		    } else {
		      // save and restore!
		      OsiObject **sosObjects = new OsiObject *[numberSOS];
		      int *newColumn = new int[numberOriginalColumns];
		      int i;
		      for (i = 0; i < numberOriginalColumns; i++)
			newColumn[i] = -1;
		      for (i = 0; i < babModel_->getNumCols(); i++)
			newColumn[originalColumns[i]] = i;
		      int nSOS = 0;
		      for (int iObj = 0; iObj < numberObjects; iObj++) {
			int iColumn = oldObjects[iObj]->columnNumber();
			if (iColumn < 0 || iColumn >= numberOriginalColumns) {
                          CbcSOS *objSOS = dynamic_cast< CbcSOS * >(oldObjects[iObj]);
                          if (objSOS) {
			    sosObjects[nSOS++] = objSOS->clone();
                            int n = objSOS->numberMembers();
                            int *which = objSOS->mutableMembers();
                            int nn = 0;
                            for (i = 0; i < n; i++) {
                              int iColumn = which[i];
                              int jColumn = newColumn[iColumn];
                              if (jColumn >= 0) {
                                which[nn++] = jColumn;
                              }
                            }
                            objSOS->setNumberMembers(nn);
                          }
			}
		      }
		      delete [] newColumn;
		      babModel_->setBestSolution(CBC_ROUNDING,
						 obj, &x[0], 1);
		      // restore!
		      nSOS = 0;
		      for (int iObj = 0; iObj < numberObjects; iObj++) {
			int iColumn = oldObjects[iObj]->columnNumber();
			if (iColumn < 0 || iColumn >= numberOriginalColumns) {
                          CbcSOS *objSOS = dynamic_cast< CbcSOS * >(oldObjects[iObj]);
                          if (objSOS) {
			    delete oldObjects[iObj];
			    oldObjects[iObj] = sosObjects[nSOS++];
                          }
			}
		      }
		      delete [] sosObjects;
		    }
		    /* But this is outside branchAndBound so needs to know 
		       about direction */
		    if (babModel_->getObjSense()==-1.0) {
		      babModel_->setCutoff(-obj);
		      babModel_->setMinimizationObjValue(-obj);
		    }
                    babModel_->clearContinuousSolver();
                    babModel_->passInSolverCharacteristics(NULL);
                    if (useSolution == 0)
                      babModel_->setHotstartSolution(&x[0]);
                  }
                }

                if (solutionIn && useSolution >= 0) {
                  if (!prioritiesIn) {
                    int n;
                    if (preProcess) {
                      int numberColumns = babModel_->getNumCols();
                      // extend arrays in case SOS
                      n = originalColumns[numberColumns - 1] + 1;
                    } else {
                      n = babModel_->getNumCols();
                    }
                    prioritiesIn =
                        reinterpret_cast<int *>(malloc(n * sizeof(int)));
                    for (int i = 0; i < n; i++)
                      prioritiesIn[i] = 100;
                  }
                  if (preProcess) {
                    int numberColumns = babModel_->getNumCols();
                    // extend arrays in case SOS
                    int n = originalColumns[numberColumns - 1] + 1;
                    int nSmaller = std::min(n, numberOriginalColumns);
                    double *solutionIn2 = new double[n];
                    int *prioritiesIn2 = new int[n];
                    int i;
                    for (i = 0; i < nSmaller; i++) {
                      solutionIn2[i] = solutionIn[i];
                      prioritiesIn2[i] = prioritiesIn[i];
                    }
                    for (; i < n; i++) {
                      solutionIn2[i] = 0.0;
                      prioritiesIn2[i] = 1000000;
                    }
#ifndef NDEBUG
                    int iLast = -1;
#endif
                    for (i = 0; i < numberColumns; i++) {
                      int iColumn = originalColumns[i];
#ifndef NDEBUG
                      assert(iColumn > iLast);
                      iLast = iColumn;
#endif
                      solutionIn2[i] = solutionIn2[iColumn];
                      if (prioritiesIn)
                        prioritiesIn2[i] = prioritiesIn2[iColumn];
                    }
                    if (useSolution)
                      babModel_->setHotstartSolution(solutionIn2,
                                                     prioritiesIn2);
                    else
                      babModel_->setBestSolution(solutionIn2, numberColumns,
                                                 COIN_DBL_MAX, true);
                    delete[] solutionIn2;
                    delete[] prioritiesIn2;
                  } else {
                    if (useSolution)
                      babModel_->setHotstartSolution(solutionIn, prioritiesIn);
                    else
                      babModel_->setBestSolution(solutionIn,
                                                 babModel_->getNumCols(),
                                                 COIN_DBL_MAX, true);
                  }
                }
                OsiSolverInterface *testOsiSolver =
                    (testOsiOptions >= 0) ? babModel_->solver() : NULL;
                if (!testOsiSolver) {
                  // *************************************************************
                  // CbcObjects
                  if (preProcess &&
                      (process.numberSOS() || babModel_->numberObjects())) {
                    int numberSOS = process.numberSOS();
                    int numberIntegers = babModel_->numberIntegers();
                    /* model may not have created objects
                                           If none then create
                                        */
                    if (!integersOK) {
                      int type = (pseudoUp) ? 1 : 0;
                      babModel_->findIntegers(true, type);
                      numberIntegers = babModel_->numberIntegers();
                    }
                    OsiObject **oldObjects = babModel_->objects();
                    // Do sets and priorities
                    OsiObject **objects = new OsiObject *[numberSOS];
                    // set old objects to have low priority
                    int numberOldObjects = babModel_->numberObjects();
                    int numberColumns = babModel_->getNumCols();
                    // backward pointer to new variables
                    // extend arrays in case SOS
                    assert(originalColumns);
                    int n = std::min(truncateColumns, numberColumns);
                    // allow for empty problem
                    n = (n) ? originalColumns[n - 1] + 1 : 0;
                    n = std::max(n,
                                std::max(numberColumns, numberOriginalColumns));
                    int *newColumn = new int[n];
                    int i;
                    for (i = 0; i < numberOriginalColumns; i++)
                      newColumn[i] = -1;
                    for (i = 0; i < std::min(truncateColumns, numberColumns);
                         i++)
                      newColumn[originalColumns[i]] = i;
                    int nMissing = 0;
                    for (int iObj = 0; iObj < numberOldObjects; iObj++) {
                      if (process.numberSOS())
                        oldObjects[iObj]->setPriority(numberColumns + 1);
                      int iColumn = oldObjects[iObj]->columnNumber();
                      if (iColumn < 0 || iColumn >= numberOriginalColumns) {
                        if (redoSOS) { // now done earlier??
                          CbcSOS *obj =
                              dynamic_cast<CbcSOS *>(oldObjects[iObj]);
                          if (obj) {
                            int n = obj->numberMembers();
                            int *which = obj->mutableMembers();
                            double *weights = obj->mutableWeights();
                            int nn = 0;
                            for (i = 0; i < n; i++) {
                              int iColumn = which[i];
                              int jColumn = newColumn[iColumn];
                              if (jColumn >= 0) {
                                which[nn] = jColumn;
                                weights[nn++] = weights[i];
                              } else {
                                nMissing++;
                              }
                            }
                            obj->setNumberMembers(nn);
                          }
                        }
                        continue;
                      }
                      if (originalColumns)
                        iColumn = originalColumns[iColumn];
                      if (branchDirection) {
                        CbcSimpleInteger *obj =
                            dynamic_cast<CbcSimpleInteger *>(oldObjects[iObj]);
                        if (obj) {
                          obj->setPreferredWay(branchDirection[iColumn]);
                        } else {
                          CbcObject *obj =
                              dynamic_cast<CbcObject *>(oldObjects[iObj]);
                          assert(obj);
                          obj->setPreferredWay(branchDirection[iColumn]);
                        }
                      }
                      if (pseudoUp) {
                        CbcSimpleIntegerPseudoCost *obj1a =
                            dynamic_cast<CbcSimpleIntegerPseudoCost *>(
                                oldObjects[iObj]);
                        assert(obj1a);
                        if (pseudoDown[iColumn] > 0.0)
                          obj1a->setDownPseudoCost(pseudoDown[iColumn]);
                        if (pseudoUp[iColumn] > 0.0)
                          obj1a->setUpPseudoCost(pseudoUp[iColumn]);
                      }
                    }
                    if (nMissing) {
#ifndef DO_LESS_PROHIBITED
                      buffer.str("");
                      buffer << nMissing << " SOS variables vanished due to preprocessing?"
                             << " - check validity?";
#else
                      buffer.str("");
                      buffer << nMissing << " SOS variables eliminated by preprocessing";
#endif
                      printGeneralMessage(model_, buffer.str());
                    }
                    delete[] newColumn;
                    const int *starts = process.startSOS();
                    const int *which = process.whichSOS();
                    const int *type = process.typeSOS();
                    const double *weight = process.weightSOS();
                    int iSOS;
                    for (iSOS = 0; iSOS < numberSOS; iSOS++) {
                      int iStart = starts[iSOS];
                      int n = starts[iSOS + 1] - iStart;
                      //#define MAKE_SOS_CLIQUES
#ifndef MAKE_SOS_CLIQUES
                      objects[iSOS] =
                          new CbcSOS(babModel_, n, which + iStart,
                                     weight + iStart, iSOS, type[iSOS]);
#else
                      objects[iSOS] = new CbcClique(
                          babModel_, 1, n, which + iStart, NULL, -iSOS - 1);
#endif
                      // branch on long sets first
                      objects[iSOS]->setPriority(numberColumns - n);
                    }
                    if (numberSOS)
                      babModel_->addObjects(numberSOS, objects);
                    for (iSOS = 0; iSOS < numberSOS; iSOS++)
                      delete objects[iSOS];
                    delete[] objects;
                  } else if (priorities || branchDirection || pseudoDown ||
                             pseudoUp || numberSOS) {
                    // do anyway for priorities etc
                    int numberIntegers = babModel_->numberIntegers();
                    /* model may not have created objects
                                           If none then create
                                        */
                    if (!numberIntegers || !babModel_->numberObjects()) {
                      int type = (pseudoUp) ? 1 : 0;
                      babModel_->findIntegers(true, type);
                    }
                    if (numberSOS) {
                      // Do sets and priorities
                      OsiObject **objects = new OsiObject *[numberSOS];
                      int iSOS;
                      if (originalColumns) {
                        // redo sequence numbers
                        int numberColumns = babModel_->getNumCols();
                        int nOld = originalColumns[numberColumns - 1] + 1;
                        int *back = new int[nOld];
                        int i;
                        for (i = 0; i < nOld; i++)
                          back[i] = -1;
                        for (i = 0; i < numberColumns; i++)
                          back[originalColumns[i]] = i;
                        // Really need better checks
                        int nMissing = 0;
#ifndef DO_LESS_PROHIBITED
                        int n = sosStart[numberSOS];
                        for (i = 0; i < n; i++) {
                          int iColumn = sosIndices[i];
                          int jColumn = back[iColumn];
                          if (jColumn >= 0)
                            sosIndices[i] = jColumn;
                          else
                            nMissing++;
                        }
#else
                        int startSet = 0;
                        int iPut = 0;
                        for (int j = 0; j < numberSOS; j++) {
                          for (i = startSet; i < sosStart[j + 1]; i++) {
                            int iColumn = sosIndices[i];
                            int jColumn = back[iColumn];
                            if (jColumn >= 0) {
                              sosReference[iPut] = sosReference[i];
                              sosIndices[iPut++] = jColumn;
                            } else {
                              nMissing++;
                            }
                          }
                          startSet = sosStart[j + 1];
                          sosStart[j + 1] = iPut;
                        }
#endif
                        delete[] back;
                        if (nMissing) {
#ifndef DO_LESS_PROHIBITED
                          buffer.str("");
                          buffer << nMissing << " SOS variables vanished due to preprocessing?"
                                 << " - check validity?";
#else
                          buffer.str("");
                          buffer << nMissing << " SOS variables eliminated by preprocessing";
#endif
                          printGeneralMessage(model_, buffer.str());
                        }
                      }
                      int sosPriorityOption =
                         parameters[CbcParam::SOSPRIORITIZE]->modeVal();
                      if (sosPriorityOption) {
                        const char *msg[4] = {
                            "high with equal priority",
                            "low with equal priority",
                            "high but with decreasing priority",
                            "low and decreasing priority"};
                        buffer.str("");
                        buffer << "Setting " << numberSOS << " SOS priorities "
                               << msg[sosPriorityOption - 1];
                        printGeneralMessage(model_, buffer.str());
                      }
                      for (iSOS = 0; iSOS < numberSOS; iSOS++) {
                        int iStart = sosStart[iSOS];
                        int n = sosStart[iSOS + 1] - iStart;
                        CbcSOS *sosObject = new CbcSOS(
                            babModel_, n, sosIndices + iStart,
                            sosReference + iStart, iSOS, sosType[iSOS]);
                        if (sosPriority) {
                          sosObject->setPriority(sosPriority[iSOS]);
                        } else if (sosPriorityOption) {
                          int priority = 10;
                          switch (sosPriorityOption) {
                          case 2:
                            priority = 100000;
                            break;
                          case 3:
                            // really should check <990 sets
                            priority = 10 + iSOS;
                            break;
                          case 4:
                            priority = 100000 + iSOS;
                            break;
                          }
                          sosObject->setPriority(priority);
                        } else if (!prioritiesIn) {
                          sosObject->setPriority(10); // rather than 1000
                        }
                        objects[iSOS] = sosObject;
                      }
                      // delete any existing SOS objects
                      int numberObjects = babModel_->numberObjects();
                      OsiObject **oldObjects = babModel_->objects();
                      int nNew = 0;
                      for (int i = 0; i < numberObjects; i++) {
                        OsiObject *objThis = oldObjects[i];
                        CbcSOS *obj1 = dynamic_cast<CbcSOS *>(objThis);
                        OsiSOS *obj2 = dynamic_cast<OsiSOS *>(objThis);
                        if (!obj1 && !obj2) {
                          oldObjects[nNew++] = objThis;
                        } else {
                          delete objThis;
                        }
                      }
                      babModel_->setNumberObjects(nNew);
                      babModel_->addObjects(numberSOS, objects);
                      for (iSOS = 0; iSOS < numberSOS; iSOS++)
                        delete objects[iSOS];
                      delete[] objects;
                    }
                  }
                  OsiObject **objects = babModel_->objects();
                  int numberObjects = babModel_->numberObjects();
                  for (int iObj = 0; iObj < numberObjects; iObj++) {
                    // skip sos
                    CbcSOS *objSOS = dynamic_cast<CbcSOS *>(objects[iObj]);
                    if (objSOS)
                      continue;
#ifdef MAKE_SOS_CLIQUES
                    // skip cliques
                    CbcClique *objClique =
                        dynamic_cast<CbcClique *>(objects[iObj]);
                    if (objClique)
                      continue;
#endif
                    int iColumn = objects[iObj]->columnNumber();
                    assert(iColumn >= 0);
                    if (originalColumns)
                      iColumn = originalColumns[iColumn];
                    if (branchDirection) {
                      CbcSimpleInteger *obj =
                          dynamic_cast<CbcSimpleInteger *>(objects[iObj]);
                      if (obj) {
                        obj->setPreferredWay(branchDirection[iColumn]);
                      } else {
                        CbcObject *obj =
                            dynamic_cast<CbcObject *>(objects[iObj]);
                        assert(obj);
                        obj->setPreferredWay(branchDirection[iColumn]);
                      }
                    }
                    if (priorities) {
                      int iPriority = priorities[iColumn];
                      if (iPriority > 0)
                        objects[iObj]->setPriority(iPriority);
                    }
                    if (pseudoUp && pseudoUp[iColumn]) {
                      CbcSimpleIntegerPseudoCost *obj1a =
                          dynamic_cast<CbcSimpleIntegerPseudoCost *>(
                              objects[iObj]);
                      assert(obj1a);
                      if (pseudoDown[iColumn] > 0.0)
                        obj1a->setDownPseudoCost(pseudoDown[iColumn]);
                      if (pseudoUp[iColumn] > 0.0)
                        obj1a->setUpPseudoCost(pseudoUp[iColumn]);
                    }
                  }
                  // *************************************************************
                } else {
                  // *************************************************************
                  // OsiObjects
                  // Find if none
                  int numberIntegers = testOsiSolver->getNumIntegers();
                  /* model may not have created objects
                                       If none then create
                                    */
                  if (!numberIntegers || !testOsiSolver->numberObjects()) {
                    // int type = (pseudoUp) ? 1 : 0;
                    testOsiSolver->findIntegers(false);
                    numberIntegers = testOsiSolver->getNumIntegers();
                  }
                  if (preProcess && process.numberSOS()) {
                    int numberSOS = process.numberSOS();
                    OsiObject **oldObjects = testOsiSolver->objects();
                    // Do sets and priorities
                    OsiObject **objects = new OsiObject *[numberSOS];
                    // set old objects to have low priority
                    int numberOldObjects = testOsiSolver->numberObjects();
                    int numberColumns = testOsiSolver->getNumCols();
                    for (int iObj = 0; iObj < numberOldObjects; iObj++) {
                      oldObjects[iObj]->setPriority(numberColumns + 1);
                      int iColumn = oldObjects[iObj]->columnNumber();
                      assert(iColumn >= 0);
                      if (iColumn >= numberOriginalColumns)
                        continue;
                      if (originalColumns)
                        iColumn = originalColumns[iColumn];
                      if (branchDirection) {
                        OsiSimpleInteger *obj =
                            dynamic_cast<OsiSimpleInteger *>(oldObjects[iObj]);
                        if (obj) {
                          obj->setPreferredWay(branchDirection[iColumn]);
                        } else {
                          OsiObject2 *obj =
                              dynamic_cast<OsiObject2 *>(oldObjects[iObj]);
                          if (obj)
                            obj->setPreferredWay(branchDirection[iColumn]);
                        }
                      }
                      if (pseudoUp) {
                        abort();
                      }
                    }
                    const int *starts = process.startSOS();
                    const int *which = process.whichSOS();
                    const int *type = process.typeSOS();
                    const double *weight = process.weightSOS();
                    int iSOS;
                    for (iSOS = 0; iSOS < numberSOS; iSOS++) {
                      int iStart = starts[iSOS];
                      int n = starts[iSOS + 1] - iStart;
                      objects[iSOS] =
                          new OsiSOS(testOsiSolver, n, which + iStart,
                                     weight + iStart, type[iSOS]);
                      // branch on long sets first
                      objects[iSOS]->setPriority(numberColumns - n);
                    }
                    testOsiSolver->addObjects(numberSOS, objects);
                    for (iSOS = 0; iSOS < numberSOS; iSOS++)
                      delete objects[iSOS];
                    delete[] objects;
                  } else if (priorities || branchDirection || pseudoDown ||
                             pseudoUp || numberSOS) {
                    if (numberSOS) {
                      // Do sets and priorities
                      OsiObject **objects = new OsiObject *[numberSOS];
                      int iSOS;
                      if (originalColumns) {
                        // redo sequence numbers
                        int numberColumns = testOsiSolver->getNumCols();
                        int nOld = originalColumns[numberColumns - 1] + 1;
                        int *back = new int[nOld];
                        int i;
                        for (i = 0; i < nOld; i++)
                          back[i] = -1;
                        for (i = 0; i < numberColumns; i++)
                          back[originalColumns[i]] = i;
                        // Really need better checks
                        int nMissing = 0;
                        int n = sosStart[numberSOS];
                        for (i = 0; i < n; i++) {
                          int iColumn = sosIndices[i];
                          int jColumn = back[iColumn];
                          if (jColumn >= 0)
                            sosIndices[i] = jColumn;
                          else
                            nMissing++;
                        }
                        delete[] back;
                        if (nMissing) {
                          buffer.str("");
                          buffer << nMissing << " SOS variables vanished due to preprocessing?"
                                 << " - check validity?";
                          printGeneralMessage(model_, buffer.str());
                        }
                      }
                      for (iSOS = 0; iSOS < numberSOS; iSOS++) {
                        int iStart = sosStart[iSOS];
                        int n = sosStart[iSOS + 1] - iStart;
                        objects[iSOS] =
                            new OsiSOS(testOsiSolver, n, sosIndices + iStart,
                                       sosReference + iStart, sosType[iSOS]);
                        if (sosPriority)
                          objects[iSOS]->setPriority(sosPriority[iSOS]);
                        else if (!prioritiesIn)
                          objects[iSOS]->setPriority(10); // rather than 1000
                      }
                      // delete any existing SOS objects
                      int numberObjects = testOsiSolver->numberObjects();
                      OsiObject **oldObjects = testOsiSolver->objects();
                      int nNew = 0;
                      for (int i = 0; i < numberObjects; i++) {
                        OsiObject *objThis = oldObjects[i];
                        OsiSOS *obj1 = dynamic_cast<OsiSOS *>(objThis);
                        OsiSOS *obj2 = dynamic_cast<OsiSOS *>(objThis);
                        if (!obj1 && !obj2) {
                          oldObjects[nNew++] = objThis;
                        } else {
                          delete objThis;
                        }
                      }
                      testOsiSolver->setNumberObjects(nNew);
                      testOsiSolver->addObjects(numberSOS, objects);
                      for (iSOS = 0; iSOS < numberSOS; iSOS++)
                        delete objects[iSOS];
                      delete[] objects;
                    }
                  }
                  OsiObject **objects = testOsiSolver->objects();
                  int numberObjects = testOsiSolver->numberObjects();
                  int logLevel = parameters[CbcParam::LPLOGLEVEL]->intVal();
                  for (int iObj = 0; iObj < numberObjects; iObj++) {
                    // skip sos
                    OsiSOS *objSOS = dynamic_cast<OsiSOS *>(objects[iObj]);
                    if (objSOS) {
                      if (logLevel > 2)
                        printf("Set %d is SOS - priority %d\n", iObj,
                               objSOS->priority());
                      continue;
                    }
                    int iColumn = objects[iObj]->columnNumber();
                    if (iColumn >= 0) {
                      if (originalColumns)
                        iColumn = originalColumns[iColumn];
                      if (branchDirection) {
                        OsiSimpleInteger *obj =
                            dynamic_cast<OsiSimpleInteger *>(objects[iObj]);
                        if (obj) {
                          obj->setPreferredWay(branchDirection[iColumn]);
                        } else {
                          OsiObject2 *obj =
                              dynamic_cast<OsiObject2 *>(objects[iObj]);
                          if (obj)
                            obj->setPreferredWay(branchDirection[iColumn]);
                        }
                      }
                      if (priorities) {
                        int iPriority = priorities[iColumn];
                        if (iPriority > 0)
                          objects[iObj]->setPriority(iPriority);
                      }
                      if (logLevel > 2)
                        printf("Obj %d is int? - priority %d\n", iObj,
                               objects[iObj]->priority());
                      if (pseudoUp && pseudoUp[iColumn]) {
                        abort();
                      }
                    }
                  }
                  // *************************************************************
                }
                int statistics = (printOptions > 0) ? printOptions : 0;
                if (!statusUserFunction_[0]) {
                  free(priorities);
                  priorities = NULL;
                  free(branchDirection);
                  branchDirection = NULL;
                  free(pseudoDown);
                  pseudoDown = NULL;
                  free(pseudoUp);
                  pseudoUp = NULL;
                  free(solutionIn);
                  solutionIn = NULL;
                  free(prioritiesIn);
                  prioritiesIn = NULL;
                  free(sosStart);
                  sosStart = NULL;
                  free(sosIndices);
                  sosIndices = NULL;
                  free(sosType);
                  sosType = NULL;
                  free(sosReference);
                  sosReference = NULL;
                  free(cut);
                  cut = NULL;
                  free(sosPriority);
                  sosPriority = NULL;
                }
                if (nodeStrategy) {
                  // change default
                  if (nodeStrategy > 2) {
                    // up or down
                    int way = (((nodeStrategy - 1) % 2) == 1) ? -1 : +1;
                    babModel_->setPreferredWay(way);
#ifdef JJF_ZERO
                    OsiObject **objects = babModel_->objects();
                    int numberObjects = babModel_->numberObjects();
                    for (int iObj = 0; iObj < numberObjects; iObj++) {
                      CbcObject *obj = dynamic_cast<CbcObject *>(objects[iObj]);
                      assert(obj);
                      obj->setPreferredWay(way);
                    }
#endif
                  }
                  if (nodeStrategy == 2 || nodeStrategy > 4) {
                    // depth
                    CbcCompareDefault compare;
                    compare.setWeight(-3.0);
                    babModel_->setNodeComparison(compare);
                  } else if (nodeStrategy == 0) {
                    // hybrid was default i.e. mixture of low depth and
                    // infeasibility
                  } else if (nodeStrategy == 1) {
                    // real fewest
                    CbcCompareDefault compare;
                    compare.setWeight(-2.0);
                    babModel_->setNodeComparison(compare);
                  }
                }
                if (cppValue >= 0) {
                  int prepro = useStrategy ? -1 : preProcess;
                  // generate code
                  fp = fopen("user_driver.cpp", "w");
                  if (fp) {
                    // generate enough to do BAB
                    babModel_->generateCpp(fp, 1);
                    OsiClpSolverInterface *osiclp =
                        dynamic_cast<OsiClpSolverInterface *>(
                            babModel_->solver());
                    // Make general so do factorization
                    int factor =
                        osiclp->getModelPtr()->factorizationFrequency();
                    osiclp->getModelPtr()->setFactorizationFrequency(200);
                    osiclp->generateCpp(fp);
                    osiclp->getModelPtr()->setFactorizationFrequency(factor);
                    // solveOptions.generateCpp(fp);
                    fclose(fp);
                    // now call generate code
                    generateCode(babModel_, "user_driver.cpp", cppValue,
                                 prepro);
                  } else {
                    std::cout << "Unable to open file user_driver.cpp"
                              << std::endl;
                  }
                }
                if (!babModel_->numberStrong() &&
                    babModel_->numberBeforeTrust() > 0)
                  babModel_->setNumberBeforeTrust(0);
                if (useStrategy) {
                  CbcStrategyDefault strategy(1, babModel_->numberStrong(),
                                              babModel_->numberBeforeTrust());
                  strategy.setupPreProcessing(1);
                  babModel_->setStrategy(strategy);
                }
                if (testOsiOptions >= 0) {
                  buffer.str("");
                  buffer << "Testing OsiObject options " << testOsiOptions;
                  printGeneralMessage(model_, buffer.str());
                  if (!numberSOS) {
                    babModel_->solver()->findIntegersAndSOS(false);
#ifdef COIN_HAS_LINK
                    // If linked then pass in model
                    OsiSolverLink *solver3 =
                        dynamic_cast<OsiSolverLink *>(babModel_->solver());
                    if (solver3) {
                      CbcHeuristicDynamic3 serendipity(*babModel_);
                      serendipity.setHeuristicName("linked");
                      int heuristicOption =
                         parameters[CbcParam::HEURISTICSTRATEGY]->modeVal();
                      if (heuristicOption)
                        babModel_->addHeuristic(&serendipity);
                      double dextra3 = parameters[CbcParam::DEXTRA3]->dblVal();
                      if (dextra3)
                        solver3->setMeshSizes(dextra3);
                      int options = parameters[CbcParam::MIPOPTIONS]->intVal() / 10000;
                      CglStored stored;
                      if (options) {
                        printf("nlp options %d\n", options);
                        /*
                                                  1 - force mini branch and
                           bound 2 - set priorities high on continuous 4 - try
                           adding OA cuts 8 - try doing quadratic linearization
                                                  16 - try expanding knapsacks
                                                              32 - OA cuts
                           strictly concave 64 - no branching at all on bilinear
                           x-x!
                                                */
                        if ((options & 2)) {
                          solver3->setBiLinearPriorities(
                              10, tightenFactor > 0.0 ? tightenFactor : 1.0);
                        } else if (tightenFactor > 0.0) {
                          // set grid size for all continuous bi-linear
                          solver3->setMeshSizes(tightenFactor);
                        }
                        if ((options & 4)) {
                          solver3->setSpecialOptions2(
                              solver3->specialOptions2() | (8 + 4));
                          // say convex
                          solver3->sayConvex((options & 32) == 0);
                        }
                        int extra1 = parameters[CbcParam::EXTRA1]->intVal();
                        if ((options & 1) != 0 && extra1 > 0)
                          solver3->setFixedPriority(extra1);
                        double cutoff = COIN_DBL_MAX;
                        if ((options & 8))
                          cutoff = solver3->linearizedBAB(&stored);
                        if (cutoff < babModel_->getCutoff()) {
                          babModel_->setCutoff(cutoff);
                          // and solution
                          // babModel_->setBestObjectiveValue(solver3->bestObjectiveValue());
                          babModel_->setBestSolution(
                              solver3->bestSolution(), solver3->getNumCols(),
                              solver3->bestObjectiveValue());
                        }
                        if ((options & 64))
                          solver3->setBranchingStrategyOnVariables(16, -1, 4);
                      }
                      solver3->setCbcModel(babModel_);
                      if (stored.sizeRowCuts())
                        babModel_->addCutGenerator(&stored, 1, "Stored");
                      CglTemporary temp;
                      babModel_->addCutGenerator(&temp, 1, "OnceOnly");
                      // choose.setNumberBeforeTrusted(2000);
                      // choose.setNumberStrong(20);
                    }
                    // For temporary testing of heuristics
                    // int testOsiOptions =
                    // parameters[CbcParam::TESTOSI,numberParameters_,parameters)]->intVal();
                    if (testOsiOptions >= 10) {
                      if (testOsiOptions >= 20)
                        testOsiOptions -= 10;
                      printf("*** Temp heuristic with mode %d\n",
                             testOsiOptions - 10);
                      OsiSolverLink *solver3 =
                          dynamic_cast<OsiSolverLink *>(babModel_->solver());
                      assert(solver3);
                      int extra1 = parameters[CbcParam::EXTRA1]->intVal();
                      solver3->setBiLinearPriority(extra1);
                      printf("bilinear priority now %d\n", extra1);
                      int extra2 = parameters[CbcParam::EXTRA2]->intVal();
                      double saveDefault = solver3->defaultBound();
                      solver3->setDefaultBound(static_cast<double>(extra2));
                      double *solution = solver3->heuristicSolution(
                          slpValue > 0 ? slpValue : 40, 1.0e-5,
                          testOsiOptions - 10);
                      solver3->setDefaultBound(saveDefault);
                      if (!solution)
                        printf("Heuristic failed\n");
                    }
#endif
                  } else {
                    // move across
                    babModel_->deleteObjects(false);
                    // babModel_->addObjects(babModel_->solver()->numberObjects(),babModel_->solver()->objects());
                  }
                  CbcBranchDefaultDecision decision;
                  if (babModel_->numberStrong()) {
                    OsiChooseStrong choose(babModel_->solver());
                    choose.setNumberBeforeTrusted(
                        babModel_->numberBeforeTrust());
                    choose.setNumberStrong(babModel_->numberStrong());
                    choose.setShadowPriceMode(testOsiOptions);
                    decision.setChooseMethod(choose);
                  } else {
                    OsiChooseVariable choose(babModel_->solver());
                    decision.setChooseMethod(choose);
                  }
                  babModel_->setBranchingMethod(decision);
                  if (useCosts && testOsiOptions >= 0) {
                    if (newPriorities) {
                      // get rid of
                      delete[] newPriorities;
                      newPriorities = NULL;
                    }
                    int numberColumns = babModel_->getNumCols();
                    int *sort = new int[numberColumns];
                    double *dsort = new double[numberColumns];
                    int *priority = new int[numberColumns];
                    const double *objective = babModel_->getObjCoefficients();
                    const double *lower = babModel_->getColLower();
                    const double *upper = babModel_->getColUpper();
                    const CoinPackedMatrix *matrix =
                        babModel_->solver()->getMatrixByCol();
                    const int *columnLength = matrix->getVectorLengths();
                    int iColumn;
                    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                      sort[iColumn] = iColumn;
                      if (useCosts == 1)
                        dsort[iColumn] = -fabs(objective[iColumn]);
                      else if (useCosts == 2)
                        dsort[iColumn] = iColumn;
                      else if (useCosts == 3)
                        dsort[iColumn] = upper[iColumn] - lower[iColumn];
                      else if (useCosts == 4)
                        dsort[iColumn] = -(upper[iColumn] - lower[iColumn]);
                      else if (useCosts == 5)
                        dsort[iColumn] = -columnLength[iColumn];
                    }
                    CoinSort_2(dsort, dsort + numberColumns, sort);
                    int level = 0;
                    double last = -1.0e100;
                    for (int i = 0; i < numberColumns; i++) {
                      int iPut = sort[i];
                      if (dsort[i] != last) {
                        level++;
                        last = dsort[i];
                      }
                      priority[iPut] = level;
                    }
                    OsiObject **objects = babModel_->objects();
                    int numberObjects = babModel_->numberObjects();
                    for (int iObj = 0; iObj < numberObjects; iObj++) {
                      OsiObject *obj = objects[iObj];
                      int iColumn = obj->columnNumber();
                      if (iColumn >= 0)
                        obj->setPriority(priority[iColumn]);
                    }
                    delete[] priority;
                    delete[] sort;
                    delete[] dsort;
                  }
                }
                checkSOS(babModel_, babModel_->solver());
                if (doSprint > 0) {
                  // Sprint for primal solves
                  ClpSolve::SolveType method = ClpSolve::usePrimalorSprint;
                  ClpSolve::PresolveType presolveType = ClpSolve::presolveOff;
                  int numberPasses = 5;
                  int options[] = {0, 3, 0, 0, 0, 0};
                  int extraInfo[] = {-1, 20, -1, -1, -1, -1};
                  extraInfo[1] = doSprint;
                  int independentOptions[] = {0, 0, 3};
                  ClpSolve clpSolve(method, presolveType, numberPasses, options,
                                    extraInfo, independentOptions);
                  // say use in OsiClp
                  clpSolve.setSpecialOption(6, 1);
                  OsiClpSolverInterface *osiclp =
                      dynamic_cast<OsiClpSolverInterface *>(
                          babModel_->solver());
                  osiclp->setSolveOptions(clpSolve);
                  osiclp->setHintParam(OsiDoDualInResolve, false);
                  // switch off row copy
                  osiclp->getModelPtr()->setSpecialOptions(
                      osiclp->getModelPtr()->specialOptions() | 256);
                  osiclp->getModelPtr()->setInfeasibilityCost(1.0e11);
                }
#ifdef COIN_HAS_LINK
                if (storedAmpl.sizeRowCuts()) {
                  if (preProcess) {
                    const int *originalColumns = process.originalColumns();
                    int numberColumns = babModel_->getNumCols();
                    int *newColumn = new int[numberOriginalColumns];
                    int i;
                    for (i = 0; i < numberOriginalColumns; i++)
                      newColumn[i] = -1;
                    for (i = 0; i < numberColumns; i++) {
                      int iColumn = originalColumns[i];
                      newColumn[iColumn] = i;
                    }
                    int *buildColumn = new int[numberColumns];
                    // Build up valid cuts
                    int nBad = 0;
                    int nCuts = storedAmpl.sizeRowCuts();
                    CglStored newCuts;
                    for (i = 0; i < nCuts; i++) {
                      const OsiRowCut *cut = storedAmpl.rowCutPointer(i);
                      double lb = cut->lb();
                      double ub = cut->ub();
                      int n = cut->row().getNumElements();
                      const int *column = cut->row().getIndices();
                      const double *element = cut->row().getElements();
                      bool bad = false;
                      for (int i = 0; i < n; i++) {
                        int iColumn = column[i];
                        iColumn = newColumn[iColumn];
                        if (iColumn >= 0) {
                          buildColumn[i] = iColumn;
                        } else {
                          bad = true;
                          break;
                        }
                      }
                      if (!bad) {
                        newCuts.addCut(lb, ub, n, buildColumn, element);
                      } else {
                        nBad++;
                      }
                    }
                    storedAmpl = newCuts;
                    if (nBad)
                      printf("%d cuts dropped\n", nBad);
                    delete[] newColumn;
                    delete[] buildColumn;
                  }
                }
#endif
#ifdef CLP_MALLOC_STATISTICS
                malloc_stats();
                malloc_stats2();
#endif
#ifndef CBC_OTHER_SOLVER
                if (outputFormat == 5) {
                  osiclp = dynamic_cast<OsiClpSolverInterface *>(
                      babModel_->solver());
                  lpSolver = osiclp->getModelPtr();
                  lpSolver->setPersistenceFlag(1);
                }
#endif
                // add in lotsizing
                if (info && statusUserFunction_[0] && info->special) {
                  int numberColumns = babModel_->getNumCols();
                  int i;
                  int n = 0;
                  if (preProcess) {
                    const int *originalColumns = process.originalColumns();
                    for (i = 0; i < numberColumns; i++) {
                      int iColumn = originalColumns[i];
                      assert(iColumn >= i);
                      int iType = info->special[iColumn];
                      if (iType) {
                        assert(iType == 1);
                        n++;
                      }
                      info->special[i] = iType;
                    }
                  }
                  if (n) {
                    int numberIntegers = 0;
                    int numberOldObjects = 0;
                    OsiObject **oldObjects = NULL;
                    const double *lower = babModel_->solver()->getColLower();
                    const double *upper = babModel_->solver()->getColUpper();
                    if (testOsiOptions < 0) {
                      // *************************************************************
                      // CbcObjects
                      numberIntegers = babModel_->numberIntegers();
                      /* model may not have created objects
                                            If none then create
                                            */
                      if (!numberIntegers || !babModel_->numberObjects()) {
                        int type = (pseudoUp) ? 1 : 0;
                        babModel_->findIntegers(true, type);
                        numberIntegers = babModel_->numberIntegers();
                      }
                      oldObjects = babModel_->objects();
                      numberOldObjects = babModel_->numberObjects();
                    } else {
                      numberIntegers = testOsiSolver->getNumIntegers();
                      if (!numberIntegers || !testOsiSolver->numberObjects()) {
                        /* model may not have created objects
                                                   If none then create
                                                */
                        testOsiSolver->findIntegers(false);
                        numberIntegers = testOsiSolver->getNumIntegers();
                      }
                      oldObjects = testOsiSolver->objects();
                      numberOldObjects = testOsiSolver->numberObjects();
                    }
                    OsiObject **objects = new OsiObject *[n];
                    n = 0;
                    // set new objects to have one lower priority
                    double ranges[] = {-COIN_DBL_MAX, -1.0, 1.0, COIN_DBL_MAX};
                    for (int iObj = 0; iObj < numberOldObjects; iObj++) {
                      int iColumn = oldObjects[iObj]->columnNumber();
                      if (iColumn >= 0 && info && info->special[iColumn]) {
                        if (lower[iColumn] <= -1.0 && upper[iColumn] >= 0.0) {
                          ranges[0] = lower[iColumn];
                          ranges[3] = upper[iColumn];
                          int priority = oldObjects[iObj]->priority();
                          if (testOsiOptions < 0) {
                            objects[n] = new CbcLotsize(babModel_, iColumn, 2,
                                                        ranges, true);
                          } else {
                            objects[n] = new OsiLotsize(testOsiSolver, iColumn,
                                                        2, ranges, true);
                          }
                          objects[n++]->setPriority(priority - 1);
                        }
                      }
                    }
                    if (testOsiOptions < 0) {
                      babModel_->addObjects(n, objects);
                    } else {
                      testOsiSolver->addObjects(n, objects);
                    }
                    for (i = 0; i < n; i++)
                      delete objects[i];
                    delete[] objects;
                  }
                }
                if (storedAmpl.sizeRowCuts()) {
                  // babModel_->addCutGenerator(&storedAmpl,1,"AmplStored");
                  int numberRowCuts = storedAmpl.sizeRowCuts();
                  for (int i = 0; i < numberRowCuts; i++) {
                    const OsiRowCut *rowCutPointer =
                        storedAmpl.rowCutPointer(i);
                    babModel_->makeGlobalCut(rowCutPointer);
                  }
                }
                // If defaults then increase trust for small models
                if (!strongChanged && false) { // out for now
                  int numberColumns = babModel_->getNumCols();
                  if (numberColumns <= 50)
                    babModel_->setNumberBeforeTrust(1000);
                  else if (numberColumns <= 100)
                    babModel_->setNumberBeforeTrust(100);
                  else if (numberColumns <= 300)
                    babModel_->setNumberBeforeTrust(50);
                }
#ifdef CBC_THREAD
                int numberThreads = parameters[CbcParam::THREADS]->intVal();
                babModel_->setNumberThreads(numberThreads % 100);
		// switch off deterministic if large problem and fastNodeDepth>0
		if (numberThreads/100==2) {
		  numberThreads -= 100;
		  if (babModel_->fastNodeDepth()>0 && babModel_->solver()->getNumRows()>2000) 
		    babModel_->setFastNodeDepth(-999);
		}
                babModel_->setThreadMode((numberThreads%1000) / 100);
#ifdef CBC_USE_OPENMP
		if (numberThreads>100) {
		  int realThreads = numberThreads%100;
		  int pretendThreads;
		  bool ompMode = true;
		  if (numberThreads<1000) {
		    if (numberThreads<800) {
		      omp_set_num_threads(realThreads);
		      pretendThreads = 2*realThreads;
		    } else {
		      babModel_->setThreadMode(1);
		      ompMode = false;
		      pretendThreads = realThreads;
		    }
		  } else {
		    int multiplier = std::max((numberThreads%1000)/100,1);
		    babModel_->setThreadMode(1);
		    pretendThreads = multiplier*realThreads;
		  }
		  char buffer[100];
		  if (ompMode) {
		    omp_set_num_threads(realThreads);
		    sprintf(buffer,"Deterministic parallel using %d OpenMP threads and pretending %d",realThreads,pretendThreads); 
		  } else {
		    sprintf(buffer,"Deterministic parallel using %d threads",realThreads); 
		  }
		  babModel_->messageHandler()->message(CBC_GENERAL, babModel_->messages())
		    << buffer << CoinMessageEol;
		  babModel_->setNumberThreads(pretendThreads);
		}
#endif
#endif
                int returnCode = 0;
                if (callBack != NULL)
                  returnCode = callBack(babModel_, 3);
                if (returnCode) {
                  // exit if user wants
                  delete babModel_;
                  babModel_ = NULL;
                  return returnCode;
                }
#ifndef CBC_OTHER_SOLVER
                osiclp =
                    dynamic_cast<OsiClpSolverInterface *>(babModel_->solver());
                lpSolver = osiclp->getModelPtr();
                int hotits = parameters[CbcParam::MAXHOTITS]->intVal();
                if (hotits > 100) {
                  osiclp->setSpecialOptions(osiclp->specialOptions() & ~32);
                  osiclp->setIntParam(OsiMaxNumIterationHotStart, hotits);
                } else {
                  osiclp->setIntParam(OsiMaxNumIterationHotStart, hotits);
                }
#elif CBC_OTHER_SOLVER == 1
#endif
                if ((experimentFlag >= 1 || strategyFlag >= 1) &&
                    abs(babModel_->fastNodeDepth()) == 1) {
		  int iType = babModel_->fastNodeDepth();
		  int iDepth = iType <0 ? -12 : 5;
		  int iSize = 500; // think harder iType <0 ? 10000 : 500; 
                  if (babModel_->solver()->getNumCols() +
                          babModel_->solver()->getNumRows() < iSize) {
                    babModel_->setFastNodeDepth(iDepth);
		  } else {
 		    if (iDepth==-12)
		      babModel_->setFastNodeDepth(-1);
		    else
		      babModel_->setFastNodeDepth(8);
		  }
                } else if (babModel_->fastNodeDepth() == -999) {
                  babModel_->setFastNodeDepth(-1);
                }
                int heurOptions = parameters[CbcParam::HEUROPTIONS]->intVal();
                if (heurOptions > 100)
                  babModel_->setSpecialOptions(babModel_->specialOptions() |
                                               8192);

#ifndef CBC_OTHER_SOLVER
#ifdef CLP_MULTIPLE_FACTORIZATIONS
                int denseCode = clpParameters[ClpParam::DENSE]->intVal();
                int smallCode = clpParameters[ClpParam::SMALLFACT]->intVal();
                if (bothFlags >= 1) {
                  if (denseCode < 0)
                    denseCode = 40;
                  if (smallCode < 0 &&
                      !lpSolver->factorization()->isDenseOrSmall())
                    smallCode = 40;
                }
                if (denseCode > 0) {
                  lpSolver->factorization()->setGoDenseThreshold(denseCode);
                  assert(osiclp == babModel_->solver());
                  osiclp->setSpecialOptions(osiclp->specialOptions() | 1024);
                }
                if (smallCode > 0 && smallCode > denseCode)
                  lpSolver->factorization()->setGoSmallThreshold(smallCode);
                // if (denseCode>=lpSolver->numberRows()) {
                // lpSolver->factorization()->goDense();
                //}
                if (lpSolver->factorization()->goOslThreshold() > 1000) {
                  // use osl in gomory (may not if CglGomory decides not to)
                  int numberGenerators = babModel_->numberCutGenerators();
                  int nGomory = 0;
                  for (int iGenerator = 0; iGenerator < numberGenerators;
                       iGenerator++) {
                    CbcCutGenerator *generator =
                        babModel_->cutGenerator(iGenerator);
                    CglGomory *gomory =
                        dynamic_cast<CglGomory *>(generator->generator());
                    if (gomory) {
                      if (nGomory < 2) {
                        gomory->useAlternativeFactorization();
                      } else if (gomory->originalSolver()) {
                        OsiClpSolverInterface *clpSolver =
                            dynamic_cast<OsiClpSolverInterface *>(
                                gomory->originalSolver());
                        if (clpSolver) {
                          ClpSimplex *simplex = clpSolver->getModelPtr();
                          simplex->factorization()->setGoOslThreshold(0);
                        }
                      }
                      nGomory++;
                    }
                  }
                }
#endif
#endif
#ifdef CLIQUE_ANALYSIS
                if (!storedAmpl.sizeRowCuts()) {
                  printf("looking at probing\n");
                  babModel_->addCutGenerator(&storedAmpl, 1, "Stored");
                }
#endif
#ifdef SOS_AS_CUTS
                /* SOS as cuts
                                   Could be a bit more sophisticated e.g. only
                   non duplicates Could do something for SOS 2?
                                */
                {
                  OsiClpSolverInterface *clpSolver =
                      dynamic_cast<OsiClpSolverInterface *>(
                          babModel_->solver());
                  if (clpSolver && clpSolver->numberSOS()) {
                    // SOS
                    int numberSOS = clpSolver->numberSOS();
                    const CoinSet *setInfo = clpSolver->setInfo();
                    CglStored sosCuts;
                    const double *lower = clpSolver->getColLower();
                    const double *upper = clpSolver->getColUpper();
                    // Start Cliques
                    // sizes
                    int nEls = 0;
                    for (int i = 0; i < numberSOS; i++)
                      nEls += setInfo[i].numberEntries();
                    double *els = new double[nEls + 2 * numberSOS];
                    for (int i = 0; i < nEls; i++)
                      els[i] = 1.0;
                    double *lo = els + nEls;
                    double *up = lo + numberSOS;
                    // need to get rid of sos
                    ClpSimplex *fakeSimplex =
                        new ClpSimplex(*clpSolver->getModelPtr());
#if 0
				    int numberRows=fakeSimplex->numberRows();
				    int * starts =
				      new int[std::max(numberSOS+1,numberRows)];
				    int * columns = new int[nEls];
				    for (int i=0;i<numberRows;i++)
				      starts[i]=i;
				    fakeSimplex->deleteRows(numberRows,starts);
#else
                    int *starts = new int[numberSOS + 1];
                    int *columns = new int[nEls];
#endif
                    int nAdded = 0;
                    starts[0] = 0;
                    nEls = 0;

                    // just SOS 1 with 0-1
                    for (int i = 0; i < numberSOS; i++) {
                      int type = setInfo[i].setType();
                      if (type == 2)
                        continue;
                      int n = setInfo[i].numberEntries();
                      const int *which = setInfo[i].which();
                      for (int j = 0; j < n; j++) {
                        int iColumn = which[j];
                        if (lower[iColumn] || upper[iColumn] != 1.0) {
                          n = -1;
                          break; // no good
                        }
                      }
                      if (n > 0) {
                        memcpy(columns + nEls, which, n * sizeof(int));
                        lo[nAdded] = -COIN_DBL_MAX;
                        up[nAdded] = 1.0;
                        nAdded++;
                        nEls += n;
                        starts[nAdded] = nEls;
                      }
                    }
                    if (nAdded)
                      fakeSimplex->addRows(nAdded, lo, up, starts, columns,
                                           els);
                    if (nAdded) {
                      OsiClpSolverInterface fakeSolver(fakeSimplex);
                      CglFakeClique fakeGen(&fakeSolver, false);
                      fakeGen.setStarCliqueReport(false);
                      fakeGen.setRowCliqueReport(false);
                      fakeGen.setMinViolation(0.05);
                      babModel_->addCutGenerator(&fakeGen, 1, "SosCuts");
                      // fakeSimplex->writeMps("bad.mps",0,1);
                      sosCuts.setProbingInfo(
                          new CglTreeProbingInfo(&fakeSolver));
                    }
                    delete fakeSimplex;
                    // End Cliques
                    // Start Stored
                    nAdded = 0;
                    for (int i = 0; i < numberSOS; i++) {
                      int type = setInfo[i].setType();
                      int n = setInfo[i].numberEntries();
                      const int *which = setInfo[i].which();
                      double rhs = 0.0;
                      double previous = 0.0;
                      for (int j = 0; j < n; j++) {
                        int iColumn = which[j];
                        if (lower[iColumn]) {
                          n = -1;
                          break; // no good
                        }
                        rhs = std::max(upper[iColumn] + previous, rhs);
                        if (type == 2)
                          previous = upper[iColumn];
                      }
                      if (n > 0) {
                        sosCuts.addCut(0.0, rhs, n, which, els);
                        nAdded++;
                      }
                    }
                    if (nAdded)
                      babModel_->addCutGenerator(&sosCuts, 1, "SosCuts2");
                    // End Stored
                    delete[] els;
                    delete[] columns;
                    delete[] starts;
                  }
                }
#endif
                 if (useSolution > 1) {
                  // use hotstart to try and find solution
                  CbcHeuristicPartial partial(*babModel_, 10000, useSolution);
                  partial.setHeuristicName("Partial solution given");
                  babModel_->addHeuristic(&partial);
                }
                if (logLevel <= 1 && babModel_->solver()->messageHandler() !=
                                         babModel_->messageHandler())
                  babModel_->solver()->setHintParam(OsiDoReducePrint, true,
                                                    OsiHintTry);
#ifdef CBC_TEMP1
                if (osiclp->getModelPtr()->perturbation() == 50)
                  osiclp->getModelPtr()->setPerturbation(52); // try less
#endif
#ifdef JJF_ZERO
                if (osiclp->getNumCols() == 29404) {
                  void restoreSolution(ClpSimplex * lpSolver,
                                       std::string fileName, int mode);
                  restoreSolution(osiclp->getModelPtr(), "debug.file", 0);
                  int numberColumns = osiclp->getNumCols();
                  const double *solution = osiclp->getColSolution();
                  const int *originalColumns = process.originalColumns();
                  for (int i = 0; i < numberColumns; i++) {
                    int iColumn = originalColumns[i];
                    if (saveSolver->isInteger(iColumn)) {
                      double value = solution[i];
                      double value2 = floor(value + 0.5);
                      assert(fabs(value - value2) < 1.0e-3);
                      saveSolver->setColLower(iColumn, value2);
                      saveSolver->setColUpper(iColumn, value2);
                    }
                  }
                  saveSolver->writeMps("fixed");
                  babModel_->setBestSolution(osiclp->getColSolution(),
                                             osiclp->getNumCols(), 1.5325e10);
                } else {
		  // for very fine tuning
		  if (parameters[CbcParam::OPTIONS]->intVal()) {
		    int addOptions = parameters[CbcParam::OPTIONS]->intVal();
		    int options = babModel_->specialOptions();
		    babModel_->setSpecialOptions(options|addOptions);
		  }
                  babModel_->branchAndBound(statistics);
                }
#else
#ifdef ORBITAL
                CbcOrbital orbit(babModel_);
                orbit.morph();
                exit(1);
#endif
                int hOp1 = parameters[CbcParam::HEUROPTIONS]->intVal() / 100000;
                if (hOp1 % 10) {
                  CbcCompareDefault compare;
                  compare.setBreadthDepth(hOp1 % 10);
                  babModel_->setNodeComparison(compare);
                } else if (hOp1 == 10) {
                  CbcCompareObjective compare;
                  babModel_->setNodeComparison(compare);
                }
#if CBC_OTHER_SOLVER == 1
                if (dynamic_cast<OsiCpxSolverInterface *>(babModel_->solver()))
                  babModel_->solver()->messageHandler()->setLogLevel(0);
#endif
                if (parameters[CbcParam::CPX]->modeVal()) {
                  babModel_->setSpecialOptions(babModel_->specialOptions() |
                                               16384);
                  // if (babModel_->fastNodeDepth()==-1)
                  babModel_->setFastNodeDepth(-2); // Use Cplex at root
                }
                int hOp2 = parameters[CbcParam::HEUROPTIONS]->intVal() / 10000;
                if (hOp2 % 10) {
                  babModel_->setSpecialOptions(babModel_->specialOptions() |
                                               16384);
                  if (babModel_->fastNodeDepth() == -1)
                    babModel_->setFastNodeDepth(-2); // Use Cplex at root
                }
                if (experimentFlag >= 5 && experimentFlag < 10000) {
                  CbcModel donor(*babModel_);
                  int options = babModel_->specialOptions();
                  donor.setSpecialOptions(options | 262144);
                  ClpSimplex *lpSolver2;
                  OsiClpSolverInterface *clpSolver2;
                  clpSolver2 =
                      dynamic_cast<OsiClpSolverInterface *>(donor.solver());
                  assert(clpSolver2);
                  lpSolver2 = clpSolver2->getModelPtr();
                  assert(lpSolver2);
                  if (lpSolver->factorization()->isDenseOrSmall()) {
                    lpSolver2->factorization()->forceOtherFactorization(0);
                    lpSolver2->factorization()->setGoOslThreshold(0);
                    lpSolver2->factorization()->setGoDenseThreshold(0);
                    lpSolver2->factorization()->setGoSmallThreshold(0);
                    lpSolver2->allSlackBasis();
                    lpSolver2->initialSolve();
                    int numberGenerators = donor.numberCutGenerators();
                    for (int iGenerator = 0; iGenerator < numberGenerators;
                         iGenerator++) {
                      CbcCutGenerator *generator =
                          donor.cutGenerator(iGenerator);
                      CglGomory *gomory =
                          dynamic_cast<CglGomory *>(generator->generator());
                      if (gomory)
                        gomory->useAlternativeFactorization(false);
                    }
                  } else {
                    printf("code this\n");
                    abort();
                  }
                  babModel_->setSpecialOptions(options | 524288);
                  CglStored *stored = new CglStored(donor.getNumCols());
                  donor.setStoredRowCuts(stored);
                  donor.branchAndBound(0);
                  babModel_->setStoredRowCuts(donor.storedRowCuts());
                  donor.setStoredRowCuts(NULL);
                }
                // We may have priorities from extra variables
                int more2 = parameters[CbcParam::MOREMOREMIPOPTIONS]->intVal();
                if (newPriorities) {
                  if (truncateColumns < babModel_->getNumCols()) {
                    // set new ones as high prority
                    babModel_->passInPriorities(newPriorities, false);
                  }
                  delete[] newPriorities;
                } else if ((more2 & (512 | 1024)) != 0) {
                  babModel_->findIntegers(true);
                  int numberIntegers = babModel_->numberIntegers();
                  int *newPriorities = new int[numberIntegers];
                  int n = numberIntegers -
                          (babModel_->getNumCols() - truncateColumns);
                  for (int i = 0; i < n; i++)
                    newPriorities[i] = babModel_->priority(i);
#if 1
                  int ixxxxxx = parameters[CbcParam::MAXNODES]->intVal();
                  int obj_priority = 1000;
                  int slack_priority = 1000;
                  if (ixxxxxx >= 1000000 && ixxxxxx < 1010000) {
                    ixxxxxx -= 1000000;
                    if (ixxxxxx == 0) {
                      obj_priority = 1000;
                      slack_priority = 1000;
                    } else if (ixxxxxx == 1) {
                      obj_priority = 10000;
                      slack_priority = 10000;
                    } else if (ixxxxxx == 2) {
                      obj_priority = 100;
                      slack_priority = 100;
                    } else if (ixxxxxx == 3) {
                      obj_priority = 100;
                      slack_priority = 10000;
                    } else if (ixxxxxx == 4) {
                      obj_priority = 10000;
                      slack_priority = 100;
                    } else if (ixxxxxx == 5) {
                      obj_priority = 100;
                      slack_priority = 200;
                    } else if (ixxxxxx == 6) {
                      obj_priority = 200;
                      slack_priority = 100;
                    } else {
                      abort();
                    }
                  }
                  if ((more2 & 512) != 0) {
                    newPriorities[n++] = obj_priority;
                  }
                  if ((more2 & 1024) != 0) {
                    for (int i = n; i < numberIntegers; i++)
                      newPriorities[i] = slack_priority;
                  }
#else
#define PRIORITY_TRY 0
#if PRIORITY_TRY == 0
#define OBJ_PRIORITY 1000
#define SLACK_PRIORITY 1000
#elif PRIORITY_TRY == 1
#define OBJ_PRIORITY 10000
#define SLACK_PRIORITY 10000
#elif PRIORITY_TRY == 2
#define OBJ_PRIORITY 100
#define SLACK_PRIORITY 100
#elif PRIORITY_TRY == 3
#define OBJ_PRIORITY 100
#define SLACK_PRIORITY 10000
#elif PRIORITY_TRY == 4
#define OBJ_PRIORITY 10000
#define SLACK_PRIORITY 100
#elif PRIORITY_TRY == 5
#define OBJ_PRIORITY 100
#define SLACK_PRIORITY 200
#elif PRIORITY_TRY == 6
#define OBJ_PRIORITY 200
#define SLACK_PRIORITY 100
#endif
                  if ((more2 & 512) != 0) {
                    newPriorities[n++] = OBJ_PRIORITY;
                  }
                  if ((more2 & 1024) != 0) {
                    for (int i = n; i < numberIntegers; i++)
                      newPriorities[i] = SLACK_PRIORITY;
                  }
#endif
                  babModel_->passInPriorities(newPriorities, false);
                  delete[] newPriorities;
                }
#ifdef JJF_ZERO
                //TODO this doesn't look like it will build at the moment
                int extra5 = parameters[CbcParam::EXTRA5]->intVal();
                if (extra5 > 0) {
                  int numberGenerators = babModel_->numberCutGenerators();
                  for (int iGenerator = 0; iGenerator < numberGenerators;
                       iGenerator++) {
                    CbcCutGenerator *generator =
                        babModel_->cutGenerator(iGenerator);
                    CglGomory *gomory =
                        dynamic_cast<CglGomory *>(generator->generator());
                    if (gomory) {
                      CglGomory gomory2(*gomory);
                      gomory2.useAlternativeFactorization(
                          !gomory->alternativeFactorization());
                      babModel_->addCutGenerator(&gomory2, -99, "Gomory2");
                    }
                  }
                }
#endif
                int specialOptions = parameters[CbcParam::STRONGSTRATEGY]->intVal();
                if (specialOptions >= 0)
                  babModel_->setStrongStrategy(specialOptions);
                if (parameters[CbcParam::CUTOFFCONSTRAINT]->modeVal()) {
                  babModel_->setCutoffAsConstraint(true);
                  int moreOptions = babModel_->moreSpecialOptions();
                  if (parameters[CbcParam::CUTOFFCONSTRAINT]->modeVal() == 4)
                    babModel_->setMoreSpecialOptions(moreOptions | 4194304);
                }
                int multipleRoot = parameters[CbcParam::MULTIPLEROOTS]->intVal();
                if (multipleRoot < 10000) {
                  babModel_->setMultipleRootTries(multipleRoot);
                } else {
                  // will be doing repeated solves and saves
                  int numberGoes = multipleRoot / 10000;
                  multipleRoot -= 10000 * numberGoes;
                  int moreOptions = babModel_->moreSpecialOptions();
                  if (numberGoes < 100) {
                    remove("global.cuts");
                    remove("global.fix");
                    moreOptions |= (67108864 | 134217728);
                  } else {
                    moreOptions |= 67108864 * (numberGoes / 100);
                    numberGoes = numberGoes % 100;
                  }
                  babModel_->setMultipleRootTries(multipleRoot);
                  babModel_->setMoreSpecialOptions(moreOptions);
                  int numberColumns = babModel_->getNumCols();
                  double *bestValues = new double[numberGoes];
                  double **bestSolutions = new double *[numberGoes];
                  int *which = new int[numberGoes];
                  int numberSolutions = 0;
                  buffer.str("");
                  buffer << "Starting " << numberGoes << " passes each with "
                         << multipleRoot % 10 << " solvers";
                  printGeneralMessage(model_, buffer.str());
                  for (int iGo = 0; iGo < numberGoes; iGo++) {
                    buffer.str("");
                    buffer << "Starting pass " << iGo + 1;
                    printGeneralMessage(model_, buffer.str());
                    CbcModel tempModel = *babModel_;
                    tempModel.setMaximumNodes(0);
                    // switch off cuts if none generated
                    int numberGenerators = tempModel.numberCutGenerators();
                    for (int iGenerator = 0; iGenerator < numberGenerators;
                         iGenerator++) {
                      CbcCutGenerator *generator =
                          tempModel.cutGenerator(iGenerator);
                      generator->setSwitchOffIfLessThan(1);
                    }
                    // random
                    tempModel.setRandomSeed(tempModel.getRandomSeed() +
                                            100000000 *
                                                (iGo + 1 + 5 * numberGoes));
                    for (int i = 0; i < tempModel.numberHeuristics(); i++)
                      tempModel.heuristic(i)->setSeed(
                          tempModel.heuristic(i)->getSeed() + 100000000 * iGo);
#ifndef CBC_OTHER_SOLVER
                    OsiClpSolverInterface *solver =
                        dynamic_cast<OsiClpSolverInterface *>(
                            tempModel.solver());
                    ClpSimplex *simplex = solver->getModelPtr();
                    int solverSeed =
                        simplex->randomNumberGenerator()->getSeed();
                    simplex->setRandomSeed(solverSeed + 100000000 * (iGo + 1));
#endif
                    tempModel.branchAndBound();
                    if (tempModel.bestSolution()) {
                      bestSolutions[numberSolutions] = CoinCopyOfArray(
                          tempModel.bestSolution(), numberColumns);
                      bestValues[numberSolutions] =
                          -tempModel.getMinimizationObjValue();
                      which[numberSolutions] = numberSolutions;
                      numberSolutions++;
                    }
                  }
                  // allow solutions
                  double sense = babModel_->solver()->getObjSense();
                  ;
                  CoinSort_2(bestValues, bestValues + numberSolutions, which);
                  babModel_->setMoreSpecialOptions(moreOptions & (~16777216));
                  for (int i = 0; i < numberSolutions; i++) {
                    int k = which[i];
                    if (bestValues[i] < babModel_->getCutoff()) {
                      babModel_->setBestSolution(bestSolutions[k],
                                                 numberColumns,
                                                 -bestValues[i] * sense, true);
                      babModel_->incrementUsed(bestSolutions[k]);
                    }
                    delete[] bestSolutions[k];
                  }
                  babModel_->setMoreSpecialOptions(moreOptions);
                  if (numberSolutions){
                     buffer.str("");
                     buffer << "Ending major passes - best solution "
                            << -bestValues[numberSolutions - 1];
                  } else {
                    buffer.str("");
                    buffer << "Ending major passes - no solution found";
                  }
                  printGeneralMessage(model_, buffer.str());
                  delete[] which;
                  delete[] bestValues;
                  delete[] bestSolutions;
                }
                if (biLinearProblem)
                   babModel_->setSpecialOptions(babModel_->specialOptions() &
                                                (~(512 | 32768)));
                babModel_->setMoreSpecialOptions2(
                    parameters[CbcParam::MOREMOREMIPOPTIONS]->intVal());
#ifdef CBC_HAS_NAUTY
                int nautyAdded = 0;
                {
                  if (parameters[CbcParam::ORBITAL]->modeVal()) {
                    int k = parameters[CbcParam::ORBITAL]->modeVal();
		    if (k < 4) {
                      babModel_->setMoreSpecialOptions2(
                          babModel_->moreSpecialOptions2() | (k * 128));
                    } else if (k == 4) {
#define MAX_NAUTY_PASS 2000
                      nautyAdded =
                          nautiedConstraints(*babModel_, MAX_NAUTY_PASS);
                    } else {
		      assert (k>=5 && k <=9);
                      if (k == 5)
                        babModel_->setMoreSpecialOptions2(
                            babModel_->moreSpecialOptions2() | 128 | 256 |
                            131072);
                      else if (k == 6)
                        babModel_->setMoreSpecialOptions2(
                            babModel_->moreSpecialOptions2() | 128 | 256 |
                            262144);
                      else if (k == 7)
                        babModel_->setMoreSpecialOptions2(
                            babModel_->moreSpecialOptions2() | 128 | 256 |
                            131072 | 262144);
                      else if (k == 8)
                        babModel_->setMoreSpecialOptions2(
                            babModel_->moreSpecialOptions2()|131072|1073741824);
                      else 
                        babModel_->setMoreSpecialOptions2(
                            babModel_->moreSpecialOptions2()|262144|1073741824);
                    }
                  }
                }
#endif
                // Set up pointer to preProcess
                if (preProcessPointer) {
                  babModel_->setPreProcess(preProcessPointer);
                }
#ifndef CBC_OTHER_SOLVER
		{
		  OsiClpSolverInterface *solver = dynamic_cast< OsiClpSolverInterface * >(babModel_->solver());
		  ClpSimplex *simplex = solver->getModelPtr();
		  // if wanted go back to old printing method
		  double value = simplex->getMinIntervalProgressUpdate();
		  if (value<=0.0) {
		    babModel_->setSecsPrintFrequency(-1.0);
		    if (value<-1.0) {
		      babModel_->setPrintFrequency(static_cast<int>(-value));
		    } else {
		      if (babModel_->getNumCols() > 2000)
			babModel_->setPrintFrequency(100);
		      else
			babModel_->setPrintFrequency(1000);
		    }
		  } else {
		    babModel_->setSecsPrintFrequency(value);
		  }
#ifdef CBC_EXPERIMENT_JJF
                  // deal with experiment
                  if (experimentValue > 0 && experimentValue < 5) {
                    if (experimentValue > 1) {
                      // go to end in strong branching
                      solver->setSpecialOptions(solver->specialOptions() & ~32);
                    }
                  }
#endif
                }
#endif
		//if (callBack != NULL)
		//returnCode = callBack(babModel_, 13);
		if (returnCode) {
		  // exit if user wants
		  delete babModel_;
		  babModel_ = NULL;
		  return returnCode;
		}
		// for very fine tuning
		if (parameters[CbcParam::OPTIONS]->intVal()) {
		  int addOptions = parameters[CbcParam::OPTIONS]->intVal();
		  int options = babModel_->specialOptions();
		  babModel_->setSpecialOptions(options|addOptions);
		}
#ifndef CBC_OTHER_SOLVER
		if (doVector) {
		  OsiClpSolverInterface *solver = dynamic_cast< OsiClpSolverInterface * >(babModel_->solver());
		  ClpSimplex *lpSolver = solver->getModelPtr();
                  ClpMatrixBase *matrix = lpSolver->clpMatrix();
                  if (dynamic_cast<ClpPackedMatrix *>(matrix)) {
                    ClpPackedMatrix *clpMatrix =
		      dynamic_cast<ClpPackedMatrix *>(matrix);
                    clpMatrix->makeSpecialColumnCopy();
		    lpSolver->setVectorMode(1);
                  }
                }
#endif
		babModel_->branchAndBound(statistics);
#ifdef CBC_HAS_NAUTY
                if (nautyAdded) {
                  int *which = new int[nautyAdded];
                  int numberOldRows =
                      babModel_->solver()->getNumRows() - nautyAdded;
                  for (int i = 0; i < nautyAdded; i++)
                    which[i] = i + numberOldRows;
                  babModel_->solver()->deleteRows(nautyAdded, which);
                  delete[] which;
                  babModel_->solver()->resolve();
                }
#endif
                if (truncateColumns < babModel_->solver()->getNumCols()) {
                  OsiSolverInterface *solverX = babModel_->solver();
                  int numberColumns = solverX->getNumCols();
                  int numberRows = solverX->getNumRows();
                  int numberDelete = numberColumns - truncateColumns;
                  int *delStuff = new int[numberDelete];
                  for (int i = 0; i < numberDelete; i++)
                    delStuff[i] = i + truncateColumns;
                  solverX->deleteCols(numberDelete, delStuff);
                  numberDelete = numberRows - truncateRows;
                  for (int i = 0; i < numberDelete; i++)
                    delStuff[i] = i + truncateRows;
                  solverX->deleteRows(numberDelete, delStuff);
                  delete[] delStuff;
                  if (truncatedRhsLower) {
                    numberRows = solverX->getNumRows();
                    for (int i = 0; i < numberRows; i++) {
                      solverX->setRowLower(i, truncatedRhsLower[i]);
                      solverX->setRowUpper(i, truncatedRhsUpper[i]);
                    }
                    delete[] truncatedRhsLower;
                    delete[] truncatedRhsUpper;
                  }
                }
                //#define CLP_FACTORIZATION_INSTRUMENT
#ifdef CLP_FACTORIZATION_INSTRUMENT
                extern double factorization_instrument(int type);
                double facTime = factorization_instrument(0);
                printf("Factorization %g seconds\n", facTime);
#endif
#endif
#ifdef COIN_DEVELOP
#ifndef JJF_ONE
                {
                  int numberColumns = babModel_->getNumCols();
                  const double *solution = babModel_->bestSolution();
                  if (solution && numberColumns < 1000) {
                    for (int i = 0; i < numberColumns; i++) {
                      if (solution[i])
                        printf("SOL %d %.18g\n", i, solution[i]);
                    }
                  }
                }
#endif
                void printHistory(const char *file /*,CbcModel * model*/);
                printHistory("branch.log" /*,babModel_*/);
#endif
                returnCode = 0;
                if (callBack != NULL)
                  returnCode = callBack(babModel_, 4);
                if (returnCode) {
                  // exit if user wants
                  model_.moveInfo(*babModel_);
                  delete babModel_;
                  babModel_ = NULL;
                  return returnCode;
                } 
#ifdef CLP_MALLOC_STATISTICS
                malloc_stats();
                malloc_stats2();
#endif
                checkSOS(babModel_, babModel_->solver());
              } else if (cbcParamCode == CbcParam::MIPLIB) {
                int typeOfCuts = babModel_->numberCutGenerators() ? 1 : -1;
                CbcStrategyDefault strategy(typeOfCuts,
                                            babModel_->numberStrong(),
                                            babModel_->numberBeforeTrust());
                // Set up pre-processing
                int translate2[] = {9999, 1, 1, 3, 2, 4, 5, 6, 6};
                if (preProcess) 
                  strategy.setupPreProcessing(translate2[preProcess]);
                babModel_->setStrategy(strategy);
#ifdef CBC_THREAD
                int numberThreads = parameters[CbcParam::THREADS]->intVal();
                babModel_->setNumberThreads(numberThreads % 100);
                babModel_->setThreadMode((numberThreads%1000) / 100);
#endif
#ifndef CBC_OTHER_SOLVER
                if (outputFormat == 5) {
                  osiclp = dynamic_cast<OsiClpSolverInterface *>(
                      babModel_->solver());
                  lpSolver = osiclp->getModelPtr();
                  lpSolver->setPersistenceFlag(1);
                }
#endif
                if (testOsiOptions >= 0) {
                  printf("Testing OsiObject options %d\n", testOsiOptions);
                  CbcBranchDefaultDecision decision;
                  OsiChooseStrong choose(babModel_->solver());
                  choose.setNumberBeforeTrusted(babModel_->numberBeforeTrust());
                  choose.setNumberStrong(babModel_->numberStrong());
                  choose.setShadowPriceMode(testOsiOptions);
                  // babModel_->deleteObjects(false);
                  decision.setChooseMethod(choose);
                  babModel_->setBranchingMethod(decision);
                }
                model_ = *babModel_;
#ifndef CBC_OTHER_SOLVER
                {
                  osiclp =
                      dynamic_cast<OsiClpSolverInterface *>(model_.solver());
                  lpSolver = osiclp->getModelPtr();
                  lpSolver->setSpecialOptions(
                      lpSolver->specialOptions() |
                      IN_BRANCH_AND_BOUND); // say is Cbc (and in branch and
                                            // bound)
                  if (lpSolver->factorization()->goOslThreshold() > 1000) {
                    // use osl in gomory (may not if CglGomory decides not to)
                    int numberGenerators = model_.numberCutGenerators();
                    for (int iGenerator = 0; iGenerator < numberGenerators;
                         iGenerator++) {
                      CbcCutGenerator *generator =
                          model_.cutGenerator(iGenerator);
                      CglGomory *gomory =
                          dynamic_cast<CglGomory *>(generator->generator());
                      if (gomory)
                        gomory->useAlternativeFactorization();
                    }
                  }
                }
#endif
                /* LL: this was done in CoinSolve.cpp: main(argc, argv).
                                   I have moved it here so that the miplib
                   directory location could be passed to CbcClpUnitTest. */
                /* JJF: No need to have 777 flag at all - user
                     says -miplib
                     */
                int extra2 = parameters[CbcParam::EXTRA2]->intVal();
                double stuff[11];
                stuff[0] = parameters[CbcParam::FAKEINCREMENT]->dblVal();
                stuff[1] = parameters[CbcParam::FAKECUTOFF]->dblVal();
                stuff[2] = parameters[CbcParam::DEXTRA3]->dblVal();
                stuff[3] = parameters[CbcParam::DEXTRA4]->dblVal();
                stuff[4] = clpParameters[ClpParam::DENSE]->intVal();
                stuff[5] = parameters[CbcParam::EXTRA1]->intVal();
                stuff[6] = parameters[CbcParam::EXTRA3]->intVal();
                stuff[7] = parameters[CbcParam::DEPTHMINIBAB]->intVal();
                stuff[8] = bothFlags;
                stuff[9] = doVector;
                stuff[10] = clpParameters[ClpParam::SMALLFACT]->intVal();
                if (dominatedCuts){
                  model_.setSpecialOptions(model_.specialOptions() | 64);
                }
                if (parameters[CbcParam::CPX]->modeVal()) {
                  model_.setSpecialOptions(model_.specialOptions() | 16384);
                  // if (model_.fastNodeDepth()==-1)
                  model_.setFastNodeDepth(-2); // Use Cplex at root
                }
                int hOp2 =
                    parameters[CbcParam::HEUROPTIONS]->intVal() / 10000;
                if (hOp2 % 10) {
                  model_.setSpecialOptions(model_.specialOptions() | 16384);
                  if (model_.fastNodeDepth() == -1)
                    model_.setFastNodeDepth(-2); // Use Cplex at root
                }
                int multipleRoot = parameters[CbcParam::MULTIPLEROOTS]->intVal();
                model_.setMultipleRootTries(multipleRoot);
                int specialOptions = parameters[CbcParam::STRONGSTRATEGY]->intVal();
                if (specialOptions >= 0)
                  model_.setStrongStrategy(specialOptions);
                if (!pumpChanged) {
                  // Make more lightweight
                  for (int iHeur = 0; iHeur < model_.numberHeuristics();
                       iHeur++) {
                    CbcHeuristic *heuristic = model_.heuristic(iHeur);
                    CbcHeuristicFPump *pump =
                        dynamic_cast<CbcHeuristicFPump *>(heuristic);
                    if (pump) {
                      CbcHeuristicFPump heuristic4(model_);
                      heuristic4.setFractionSmall(0.5);
                      heuristic4.setMaximumPasses(5);
                      heuristic4.setFeasibilityPumpOptions(30);
                      heuristic4.setWhen(13);
                      heuristic4.setHeuristicName("feasibility pump");
                      // CbcHeuristicFPump & pump2 = pump;
                      *pump = heuristic4;
                    }
                  }
                }
#ifndef CBC_OTHER_SOLVER
                {
                  OsiClpSolverInterface *solver =
                      dynamic_cast<OsiClpSolverInterface *>(model_.solver());
                  ClpSimplex *simplex = solver->getModelPtr();
                  // if wanted go back to old printing method
                  double value = simplex->getMinIntervalProgressUpdate();
                  if (value <= 0.0) {
                    model_.setSecsPrintFrequency(-1.0);
                    if (value < 0.0) {
                      model_.setPrintFrequency(static_cast<int>(-value));
                    }
                  } else {
                    model_.setSecsPrintFrequency(value);
                  }
                }
#endif
                std::string dirMiplib =
                   parameters[CbcParam::DIRMIPLIB]->dirName();
                if (dirMiplib != ""){
                   int returnCode =
                      CbcClpUnitTest(model_, dirMiplib, extra2, stuff,
				     saveInputQueue, callBack, parameters);
                   babModel_ = NULL;
                   return returnCode;
                }
              } else {
                abort(); // can't get here
                // strengthenedModel = babModel_->strengthenedModel();
              }
              currentBranchModel = NULL;
#ifndef CBC_OTHER_SOLVER
              osiclp =
                  dynamic_cast<OsiClpSolverInterface *>(babModel_->solver());
              if (parameters[CbcParam::DEBUGFILE]->fileName() ==
                  "createAfterPre" && babModel_->bestSolution()) {
                lpSolver = osiclp->getModelPtr();
                // move best solution (should be there -- but ..)
                int n = lpSolver->getNumCols();
                memcpy(lpSolver->primalColumnSolution(),
                       babModel_->bestSolution(), n * sizeof(double));
                ClpParamUtils::saveSolution(osiclp->getModelPtr(), "debug.file");
              }
#endif
              statistics_cut_time = 0.0;
	      double direction =
		babModel_->solver()->getObjSense();
                // Print more statistics
                buffer.str("");
                buffer << "Cuts at root node changed objective from "
                       << babModel_->getContinuousObjective()*direction
		       << " to " 
                       << babModel_->rootObjectiveAfterCuts()*direction;
                printGeneralMessage(model_, buffer.str());
                numberGenerators = babModel_->numberCutGenerators();
                // can get here twice!
                if (statistics_number_cuts != NULL)
                  delete[] statistics_number_cuts;
                statistics_number_cuts = new int[numberGenerators];

                if (statistics_name_generators != NULL)
                  delete[] statistics_name_generators;
                statistics_name_generators = new const char *[numberGenerators];

                statistics_number_generators = numberGenerators;

                char timing[30];
                for (iGenerator = 0; iGenerator < numberGenerators;
                     iGenerator++) {
                  CbcCutGenerator *generator =
                      babModel_->cutGenerator(iGenerator);
                  statistics_name_generators[iGenerator] =
                      generator->cutGeneratorName();
                  statistics_number_cuts[iGenerator] =
                     generator->numberCutsInTotal();
                  buffer.str("");
                  buffer << generator->cutGeneratorName() << " was tried "
                         << generator->numberTimesEntered() << " times and created "
                         << generator->numberCutsInTotal() + generator->numberColumnCuts()
                         << " cuts";
                  if (generator->timing()) {
                     buffer << " (" << generator->timeInCutGenerator() <<  " seconds)";
                     statistics_cut_time += generator->timeInCutGenerator();
                  }
                  CglStored *stored = dynamic_cast<CglStored *>(generator->generator());
                  if (stored && !generator->numberCutsInTotal()){
                    continue;
                  }
#ifndef CLP_INVESTIGATE
                  CglImplication *implication =
                      dynamic_cast<CglImplication *>(generator->generator());
                  if (implication && !generator->numberCutsInTotal())
                    continue;
#endif
                    printGeneralMessage(model_, buffer.str());
                }
#ifdef COIN_DEVELOP
                printf("%d solutions found by heuristics\n",
                       babModel_->getNumberHeuristicSolutions());
                // Not really generator but I am feeling lazy
                for (iGenerator = 0; iGenerator < babModel_->numberHeuristics();
                     iGenerator++) {
                  CbcHeuristic *heuristic = babModel_->heuristic(iGenerator);
                  if (heuristic->numRuns()) {
                    // Need to bring others inline
                    buffer.str("");
                    buffer << heuristic->heuristicName() << " was tried "
                           << heuristic->numRuns() << " times out of "
                           << heuristic->numCouldRun() << " and created "
                           << heuristic->numberSolutionsFound()
                           << " solutions" << std::endl;
                    printGeneralMessage(model_, buffer.str());
                  }
                }
#endif
              // adjust time to allow for children on some systems
              time2 = CoinCpuTime() + CoinCpuTimeJustChildren();
              totalTime += time2 - time1;
              // For best solution
              double *bestSolution = NULL;
              // Say in integer
              if (babModel_->status()) {
                // treat as stopped
                integerStatus = 3;
              } else {
                if (babModel_->isProvenOptimal()) {
                  integerStatus = 0;
                } else if (!babModel_->bestSolution()) {
                  // infeasible
                  integerStatus = 6;
                  delete saveSolver;
                  saveSolver = NULL;
                }
              }
              if (babModel_->getMinimizationObjValue() < 1.0e50 &&
                  cbcParamCode == CbcParam::BAB) {
                // post process
                int n;
                if (preProcess) {
                  n = saveSolver->getNumCols();
                  bestSolution = new double[n];
#ifndef CBC_OTHER_SOLVER
                  OsiClpSolverInterface *clpSolver =
                      dynamic_cast<OsiClpSolverInterface *>(
                          babModel_->solver());
#else
                  OsiCpxSolverInterface *clpSolver =
                      dynamic_cast<OsiCpxSolverInterface *>(
                          babModel_->solver());
#endif
                  // Save bounds on processed model
                  const int *originalColumns = process.originalColumns();
                  int numberColumns2 = clpSolver->getNumCols();
                  double *lower2 = new double[n];
                  double *upper2 = new double[n];
                  for (int i = 0; i < n; i++) {
                    lower2[i] = COIN_DBL_MAX;
                    upper2[i] = -COIN_DBL_MAX;
                  }
                  const double *columnLower = clpSolver->getColLower();
                  const double *columnUpper = clpSolver->getColUpper();
                  for (int i = 0; i < numberColumns2; i++) {
                    int jColumn = originalColumns[i];
                    if (jColumn < n) {
                      lower2[jColumn] = columnLower[i];
                      upper2[jColumn] = columnUpper[i];
                    }
                  }
#ifndef CBC_OTHER_SOLVER
                  ClpSimplex *lpSolver = clpSolver->getModelPtr();
                  lpSolver->setSpecialOptions(
                      lpSolver->specialOptions() |
                      IN_BRANCH_AND_BOUND); // say is Cbc (and in branch and
                                            // bound)
#endif
                  // put back any saved solutions
                  putBackOtherSolutions(babModel_, &model_, &process);
                  setPreProcessingMode(babModel_->solver(), 2);
                  process.postProcess(*babModel_->solver());
                  setPreProcessingMode(saveSolver, 0);
#ifdef COIN_DEVELOP
                  if (model_.bestSolution() &&
                      fabs(model_.getMinimizationObjValue() -
                           babModel_->getMinimizationObjValue()) < 1.0e-8) {
                    const double *b1 = model_.bestSolution();
                    const double *b2 = saveSolver->getColSolution();
                    const double *columnLower = saveSolver->getColLower();
                    const double *columnUpper = saveSolver->getColUpper();
                    for (int i = 0; i < n; i++) {
                      if (fabs(b1[i] - b2[i]) > 1.0e-7) {
                        printf("%d %g %g %g %g\n", i, b1[i], b2[i],
                               columnLower[i], columnUpper[i]);
                      }
                    }
                  }
#endif
                  bool tightenB = false;
                  if (testOsiOptions<0) {
                    int n = babModel_->numberObjects();
                    for (int i = 0; i < n; i++) {
                      const OsiObject *obj = babModel_->object(i);
                      if (!dynamic_cast<const CbcSimpleInteger *>(obj)) {
                        tightenB = true;
                        break;
                      }
                    }
                  }
                  // Solution now back in saveSolver
                  // Double check bounds
                  columnLower = saveSolver->getColLower();
                  columnUpper = saveSolver->getColUpper();
		  // skip tighten if not integer if ((process.options()&128)!=0)
		    tightenB = false;
                  int numberChanged = 0;
		  if ((process.options()&256)==0) {
		    for (int i = 0; i < n; i++) {
		      if (!saveSolver->isInteger(i) && !tightenB)
			continue;
		      if (lower2[i] != COIN_DBL_MAX) {
			if (lower2[i] != columnLower[i] || upper2[i] != columnUpper[i]) {
			  if (lower2[i] < columnLower[i] || upper2[i] > columnUpper[i]) {
#ifdef COIN_DEVELOP
			    printf("odd bounds tighter");
			    printf("%d bab bounds %g %g now %g %g\n",
				   i, lower2[i], upper2[i], columnLower[i],
				   columnUpper[i]);
#endif
			  } else {
#ifdef COIN_DEVELOP
			    printf("%d bab bounds %g %g now %g %g\n",
				   i, lower2[i], upper2[i], columnLower[i],
				   columnUpper[i]);
#endif
			    numberChanged++;
			    saveSolver->setColLower(i, lower2[i]);
			    saveSolver->setColUpper(i, upper2[i]);
			  }
			}
		      }
		    }
#ifdef JJF_ZERO
		    // See if sos so we can fix
		    OsiClpSolverInterface *osiclp = dynamic_cast< OsiClpSolverInterface * >(saveSolver);
		    if (osiclp && osiclp->numberSOS()) {
		      // SOS
		      numberSOS = osiclp->numberSOS();
		      const CoinSet *setInfo = osiclp->setInfo();
		      int i;
		      for (i = 0; i < numberSOS; i++) {
			int type = setInfo[i].setType();
			int n = setInfo[i].numberEntries();
			const int *which = setInfo[i].which();
			int first = -1;
			int last = -1;
			for (int j = 0; j < n; j++) {
			  int iColumn = which[j];
			  if (fabs(solution[iColumn]) > 1.0e-7) {
			    last = j;
			    if (first < 0)
			      first = j;
			  }
			}
			assert(last - first < type);
			for (int j = 0; j < n; j++) {
			  if (j < first || j > last) {
			    int iColumn = which[j];
			    saveSolver->setColLower(iColumn, 0.0);
			    saveSolver->setColUpper(iColumn, 0.0);
			  }
			}
		      }
		    }
#endif
		  }
                  delete[] lower2;
                  delete[] upper2;
                  if (numberChanged) {
                    buffer.str("");
                    buffer << numberChanged << " bounds tightened after postprocessing"
                           << std::endl;
                    printGeneralMessage(model_, buffer.str());
                  }
                  // saveSolver->resolve();
                  if (true /*!saveSolver->isProvenOptimal()*/) {
                    // try all slack
                    CoinWarmStartBasis *basis =
                        dynamic_cast<CoinWarmStartBasis *>(
                            babModel_->solver()->getEmptyWarmStart());
                    saveSolver->setWarmStart(basis);
                    delete basis;
                    saveSolver->initialSolve();
#ifdef COIN_DEVELOP
                    saveSolver->writeMps("inf2");
#endif
                    OsiClpSolverInterface *osiclp =
                        dynamic_cast<OsiClpSolverInterface *>(saveSolver);
                    if (osiclp)
                      osiclp->getModelPtr()->checkUnscaledSolution();
                  }
                  // assert(saveSolver->isProvenOptimal());
#ifndef CBC_OTHER_SOLVER
                  // and original solver
                  originalSolver->setDblParam(OsiDualObjectiveLimit,
                                              COIN_DBL_MAX);
                  assert(n >= originalSolver->getNumCols());
                  n = originalSolver->getNumCols();
                  originalSolver->setColLower(saveSolver->getColLower());
                  originalSolver->setColUpper(saveSolver->getColUpper());
                  // basis
                  CoinWarmStartBasis *basis =
                      dynamic_cast<CoinWarmStartBasis *>(
                          babModel_->solver()->getWarmStart());
                  originalSolver->setBasis(*basis);
                  delete basis;
                  originalSolver->resolve();
                  if (!originalSolver->isProvenOptimal()) {
                    // try all slack
                    CoinWarmStartBasis *basis =
                        dynamic_cast<CoinWarmStartBasis *>(
                            babModel_->solver()->getEmptyWarmStart());
                    originalSolver->setBasis(*basis);
                    delete basis;
                    originalSolver->initialSolve();
                    OsiClpSolverInterface *osiclp =
                        dynamic_cast<OsiClpSolverInterface *>(originalSolver);
                    if (osiclp)
                      osiclp->getModelPtr()->checkUnscaledSolution();
                  }
                  // assert(originalSolver->isProvenOptimal());
#endif
                  babModel_->assignSolver(saveSolver);
                  memcpy(bestSolution, originalSolver->getColSolution(),
                         n * sizeof(double));
		  // already set babModel_->setObjValue(babModel_->solver()->getObjValue());
                } else {
                  n = babModel_->solver()->getNumCols();
                  bestSolution = new double[n];
                  memcpy(bestSolution, babModel_->solver()->getColSolution(),
                         n * sizeof(double));
		  babModel_->setObjValue(babModel_->solver()->getObjValue());
#ifndef CBC_OTHER_SOLVER
		  // and put back bounds in very original solver
                  ClpSimplex *original = originalSolver->getModelPtr();
                  double *lower = original->columnLower();
                  double *upper = original->columnUpper();
                  int n2 = original->numberColumns();
		  memcpy(lower,babModel_->solver()->getColLower(),
			 n2*sizeof(double));
		  memcpy(upper,babModel_->solver()->getColUpper(),
			 n2*sizeof(double));
		  originalSolver->resolve();
#endif
                }
                if (returnMode == 1 && model_.numberSavedSolutions() < 2) {
                  model_.deleteSolutions();
                  model_.setBestSolution(bestSolution, n,
                                         babModel_->getMinimizationObjValue());
                }
                babModel_->deleteSolutions();
                babModel_->setBestSolution(
                    bestSolution, n, babModel_->getMinimizationObjValue());
#ifndef CBC_OTHER_SOLVER
                // and put back in very original solver
                {
                  ClpSimplex *original = originalSolver->getModelPtr();
                  double *lower = original->columnLower();
                  double *upper = original->columnUpper();
                  double *solution = original->primalColumnSolution();
                  int n = original->numberColumns();
                  // assert (!n||n==babModel_->solver()->getNumCols());
                  for (int i = 0; i < n; i++) {
                    solution[i] = bestSolution[i];
                    if (originalSolver->isInteger(i)) {
                      lower[i] = solution[i];
                      upper[i] = solution[i];
                    }
                  }
                  // basis
                  CoinWarmStartBasis *basis =
                      dynamic_cast<CoinWarmStartBasis *>(
                          babModel_->solver()->getWarmStart());
                  originalSolver->setBasis(*basis);
                  delete basis;
                  originalSolver->setDblParam(OsiDualObjectiveLimit,
                                              COIN_DBL_MAX);
#ifdef COIN_HAS_LINK
                  if (originalSolver->getMatrixByCol())
                    originalSolver->setHintParam(OsiDoPresolveInResolve, true,
                                                 OsiHintTry);
#else
                  originalSolver->setHintParam(OsiDoPresolveInResolve, true,
                                               OsiHintTry);
#endif
                  originalSolver->resolve();
                  if (!originalSolver->isProvenOptimal()) {
                    // try all slack
                    CoinWarmStartBasis *basis =
                        dynamic_cast<CoinWarmStartBasis *>(
                            babModel_->solver()->getEmptyWarmStart());
                    originalSolver->setBasis(*basis);
                    delete basis;
                    originalSolver->initialSolve();
                    OsiClpSolverInterface *osiclp =
                        dynamic_cast<OsiClpSolverInterface *>(originalSolver);
                    if (osiclp)
                      osiclp->getModelPtr()->checkUnscaledSolution();
#ifdef CLP_INVESTIGATE
                    if (!originalSolver->isProvenOptimal()) {
                      if (saveSolver) {
                        printf(
                            "saveSolver and originalSolver matrices saved\n");
                        saveSolver->writeMps("infA");
                      } else {
                        printf("originalSolver matrix saved\n");
                        originalSolver->writeMps("infB");
                      }
                    }
#endif
                  }
                  // assert(originalSolver->isProvenOptimal());
                }
#endif
                checkSOS(babModel_, babModel_->solver());
              } else if (model_.bestSolution() &&
                         cbcParamCode == CbcParam::BAB &&
                         model_.getMinimizationObjValue() < 1.0e50 &&
                         preProcess) {
                buffer.str("");
                buffer << "Restoring heuristic best solution of "
                       << model_.getMinimizationObjValue();
                printGeneralMessage(model_, buffer.str());
                int n = saveSolver->getNumCols();
                bestSolution = new double[n];
                // Put solution now back in saveSolver
                saveSolver->setColSolution(model_.bestSolution());
                babModel_->assignSolver(saveSolver);
                saveSolver = NULL;
                babModel_->setMinimizationObjValue(
                    model_.getMinimizationObjValue());
                memcpy(bestSolution, babModel_->solver()->getColSolution(),
                       n * sizeof(double));
#ifndef CBC_OTHER_SOLVER
                // and put back in very original solver
                {
                  ClpSimplex *original = originalSolver->getModelPtr();
                  double *lower = original->columnLower();
                  double *upper = original->columnUpper();
                  double *solution = original->primalColumnSolution();
                  int n = original->numberColumns();
                  // assert (!n||n==babModel_->solver()->getNumCols());
                  for (int i = 0; i < n; i++) {
                    solution[i] = bestSolution[i];
                    if (originalSolver->isInteger(i)) {
                      lower[i] = solution[i];
                      upper[i] = solution[i];
                    }
                  }
                  // basis
                  CoinWarmStartBasis *basis =
                      dynamic_cast<CoinWarmStartBasis *>(
                          babModel_->solver()->getWarmStart());
                  originalSolver->setBasis(*basis);
                  delete basis;
                }
#endif
              }
#ifndef CBC_OTHER_SOLVER
              // if (type==STRENGTHEN&&strengthenedModel)
              // clpSolver = dynamic_cast< OsiClpSolverInterface*>
              // (strengthenedModel);
              else if (statusUserFunction_[0])
                clpSolver =
                    dynamic_cast<OsiClpSolverInterface *>(babModel_->solver());
              lpSolver = clpSolver->getModelPtr();
              if (numberChanged) {
                for (int i = 0; i < numberChanged; i++) {
                  int iColumn = changed[i];
                  clpSolver->setContinuous(iColumn);
                }
                delete[] changed;
              }
#endif
              if (cbcParamCode == CbcParam::BAB) {
#ifndef CBC_OTHER_SOLVER
                // move best solution (should be there -- but ..)
                int n = lpSolver->getNumCols();
                if (bestSolution) {
                  memcpy(lpSolver->primalColumnSolution(), bestSolution,
                         n * sizeof(double));
                  // now see what that does to row solution
                  int numberRows = lpSolver->numberRows();
                  double *rowSolution = lpSolver->primalRowSolution();
                  memset(rowSolution, 0, numberRows * sizeof(double));
                  lpSolver->clpMatrix()->times(1.0, bestSolution, rowSolution);
                  lpSolver->setObjectiveValue(babModel_->getObjValue());
                }
                if (parameters[CbcParam::DEBUGFILE]->fileName() ==
                    "create" && bestSolution) {
                   ClpParamUtils::saveSolution(lpSolver, "debug.file");
                }
#else
                if (bestSolution) {
                  model_.solver()->setColSolution(bestSolution);
                }
#endif
                delete saveSolver;
                delete[] bestSolution;
                std::string statusName[] = {
                    "", "Stopped on ", "Run abandoned", "", "", "User ctrl-c"};
                std::string minor[] = {
                    "Optimal solution found",
                    "Linear relaxation infeasible",
                    "Optimal solution found (within gap tolerance)",
                    "node limit",
                    "time limit",
                    "user ctrl-c",
                    "solution limit",
                    "Linear relaxation unbounded",
                    "Problem proven infeasible"};
                int iStat = babModel_->status();
                int iStat2 = babModel_->secondaryStatus();
                if (!iStat && !iStat2 && !bestSolution)
                  iStat2 = 8;
                if (!iStat && iStat2 == 1 && bestSolution)
                  iStat2 = 0; // solution and search completed
                statistics_seconds = time2 - time1;
                statistics_sys_seconds = CoinSysTime();
                statistics_elapsed_seconds = CoinWallclockTime();
                statistics_obj = babModel_->getObjValue();
                statistics_continuous = babModel_->getContinuousObjective();
                statistics_tighter = babModel_->rootObjectiveAfterCuts();
                statistics_nodes = babModel_->getNodeCount();
                statistics_iterations = babModel_->getIterationCount();
                statistics_result = statusName[iStat];
                buffer.str("");
                buffer << std::endl << "Result - "
                       << statusName[iStat].c_str() << minor[iStat2].c_str()
                       << std::endl;
                printGeneralMessage(model_, buffer.str());
		buffer.str("");
                if (babModel_->bestSolution()) {
		    /* This is an important value 
		       - don't lose too much precision but try
		       and make pretty - suggestions welcome -
		       so back to sprintf */
		    char line[20];
		    sprintf(line,"%.12g",babModel_->getObjValue());
                    buffer << "Objective value:                "
                           << line << std::endl;
                } else {
                    buffer << "No feasible solution found\n";
                }
                if (iStat2 >= 2 && iStat2 <= 6) {
                   bool minimizing = babModel_->solver()->getObjSense() > 0.0;
                   std::string tmp = minimizing ? "Lower" : "Upper";
                   buffer << tmp << " bound:                    "
                          << babModel_->getBestPossibleObjValue() << std::endl;
                   if (babModel_->bestSolution()) {
                      buffer << "Gap:                            "
                             << (babModel_->getObjValue() -
                                 babModel_->getBestPossibleObjValue()) /
			fabs(babModel_->getBestPossibleObjValue()) << std::endl;
                   }
                }
                buffer << "Enumerated nodes:               "
                       << babModel_->getNodeCount() << std::endl;
                buffer << "Total iterations:               ";
                buffer << babModel_->getIterationCount() << std::endl;
#if CBC_QUIET == 0
                buffer << "Time (CPU seconds):             "
                       << CoinCpuTime() - time1 << std::endl;
                buffer << "Time (Wallclock seconds):       "
                       << CoinGetTimeOfDay() - time1Elapsed << std::endl;
#endif
                printGeneralMessage(model_, buffer.str());
                int returnCode = 0;
                if (callBack != NULL)
                  returnCode = callBack(babModel_, 5);
                if (returnCode) {
                  // exit if user wants
                  model_.moveInfo(*babModel_);
                  delete babModel_;
                  babModel_ = NULL;
                  return returnCode;
                }
                if (info && statusUserFunction_[0]) {
                  clpSolver = dynamic_cast<OsiClpSolverInterface *>(
                      babModel_->solver());
                  lpSolver = clpSolver->getModelPtr();
                  double value = babModel_->getObjValue();
                  char buf[300];
                  int pos = 0;
                  if (iStat == 0) {
                    if (fabs(value) < 1.0e40) {
                      pos += sprintf(buf + pos, "optimal,");
                    } else {
                      // infeasible
                      iStat = 1;
                      pos += sprintf(buf + pos, "infeasible,");
                    }
                  } else if (iStat == 1) {
                    if (iStat2 != 6)
                      iStat = 3;
                    else
                      iStat = 4;
                    pos += sprintf(buf + pos, "stopped on %s,",
                                   minor[iStat2].c_str());
                  } else if (iStat == 2) {
                    iStat = 7;
                    pos += sprintf(buf + pos, "stopped on difficulties,");
                  } else if (iStat == 5) {
                    iStat = 3;
                    pos += sprintf(buf + pos, "stopped on ctrl-c,");
                  } else {
                    pos += sprintf(buf + pos, "status unknown,");
                    iStat = 6;
                  }
                  info->problemStatus = iStat;
                  info->objValue = value;
                  if (fabs(value) < 1.0e40) {
                    int precision = ampl_obj_prec();
                    if (precision > 0)
                      pos += sprintf(buf + pos, " objective %.*g", precision,
                                     value);
                    else
                      pos += sprintf(buf + pos, " objective %g", value);
                  }
                  sprintf(buf + pos, "\n%d nodes, %d iterations, %g seconds",
                          babModel_->getNodeCount(),
                          babModel_->getIterationCount(), totalTime);
                  if (bestSolution) {
                    free(info->primalSolution);
                    if (!numberKnapsack) {
                      info->primalSolution =
                          (double *)malloc(n * sizeof(double));
                      CoinCopyN(lpSolver->primalColumnSolution(), n,
                                info->primalSolution);
                      int numberRows = lpSolver->numberRows();
                      free(info->dualSolution);
                      info->dualSolution =
                          (double *)malloc(numberRows * sizeof(double));
                      CoinCopyN(lpSolver->dualRowSolution(), numberRows,
                                info->dualSolution);
                    } else {
                      // expanded knapsack
                      info->dualSolution = NULL;
                      int numberColumns = saveCoinModel.numberColumns();
                      info->primalSolution =
                          (double *)malloc(numberColumns * sizeof(double));
                      // Fills in original solution (coinModel length)
                      afterKnapsack(saveTightenedModel, whichColumn,
                                    knapsackStart, knapsackRow, numberKnapsack,
                                    lpSolver->primalColumnSolution(),
                                    info->primalSolution, 1);
                    }
                  } else {
                    info->primalSolution = NULL;
                    info->dualSolution = NULL;
                  }
                  // put buffer into info
                  strcpy(info->buffer, buf);
                }
              } else {
                buffer.str("");
                buffer << "Model strengthened - now has "
                       << clpSolver->getNumRows() << " rows";
                printGeneralMessage(model_, buffer.str());
              }
              time1 = time2;
              if (statusUserFunction_[0]) {
                // keep if going to be destroyed
                OsiSolverInterface *solver = babModel_->solver();
                OsiClpSolverInterface *clpSolver =
                    dynamic_cast<OsiClpSolverInterface *>(solver);
                ClpSimplex *lpSolver2 = clpSolver->getModelPtr();
                if (lpSolver == lpSolver2)
                  babModel_->setModelOwnsSolver(false);
              }
              // delete babModel_;
              // babModel_=NULL;
          } break;
          case CbcParam::IMPORT: {
            if (!statusUserFunction_[0]) {
              free(priorities);
              priorities = NULL;
              free(branchDirection);
              branchDirection = NULL;
              free(pseudoDown);
              pseudoDown = NULL;
              free(pseudoUp);
              pseudoUp = NULL;
              free(solutionIn);
              solutionIn = NULL;
              free(prioritiesIn);
              prioritiesIn = NULL;
              free(sosStart);
              sosStart = NULL;
              free(sosIndices);
              sosIndices = NULL;
              free(sosType);
              sosType = NULL;
              free(sosReference);
              sosReference = NULL;
              free(cut);
              cut = NULL;
              free(sosPriority);
              sosPriority = NULL;
            }
            // delete babModel_;
            // babModel_=NULL;
            // see if we have a file name
            cbcParam->readValue(inputQueue, fileName, &message);
            // TODO Think about how to do this right
            canOpen = false;
            bool absolutePath = true;
            CoinParamUtils::processFile(fileName,
                                 parameters[CbcParam::DIRECTORY]->dirName());
            if (fileName == ""){
               fileName = parameters[CbcParam::IMPORTFILE]->fileName();
            }else{
               parameters[CbcParam::IMPORTFILE]->setFileName(fileName);
            }
            if (fileName[0] != '/' && fileName[0] != '\\' &&
                !strchr(fileName.c_str(), ':')) {
               absolutePath = false;
            }
            // See if gmpl file
            int gmpl = 0;
            std::string gmplData;
            if (fileName == "-") {
              // stdin
              canOpen = true;
            } else if (fileName == "-lp") {
              // stdin
              canOpen = true;
              fileName = "-";
              gmpl = -1; //.lp format
            } else {
              // See if .lp
              const char *c_name = fileName.c_str();
              size_t length = strlen(c_name);
              if ((length > 3 && !strncmp(c_name + length - 3, ".lp", 3)) ||
                  (length > 6 &&
                   !strncmp(c_name + length - 6, ".lp.gz", 6)) ||
                  (length > 7 && !strncmp(c_name + length - 7, ".lp.bz2", 7))){
                 gmpl = -1; // .lp
              }
              length = fileName.size();
              size_t percent = fileName.find('%');
              if (percent < length && percent > 0) {
#ifdef COINUTILS_HAS_GLPK
                 gmpl = 1;
                 fileName = fileName.substr(0, percent);
                 gmplData = fileName.substr(percent + 1);
                 if (!absolutePath){
                    fileName = parameters[CbcParam::DIRECTORY]->dirName() +
                       fileName;
                    gmplData = parameters[CbcParam::DIRECTORY]->dirName() +
                       fileName;
                 }
                 gmpl = (percent < length - 1) ? 2 : 1;
                 printf("GMPL model file %s and data file %s\n",
                        fileName.c_str(), gmplData.c_str());
#else
                 printf("Cbc was not built with GMPL support. Exiting.\n");
                 // This is surely not the right thing to do here. Should we
                 // throw an exceptioon? Exit?
                 abort();
#endif
              }
              if (fileCoinReadable(fileName)) {
                // can open - lets go for it
                canOpen = true;
                if (gmpl == 2) {
                  fp;
                  fp = fopen(gmplData.c_str(), "r");
                  if (fp) {
                    fclose(fp);
                  } else {
                    canOpen = false;
                    buffer.str("");
                    buffer << "Unable to open file " << gmplData.c_str();
                    printGeneralMessage(model_, buffer.str());
                    continue;
                  }
                }
              } else {
                buffer.str("");
                buffer << "Unable to open file " << fileName.c_str();
                printGeneralMessage(model_, buffer.str());
                continue;
              }
            }
            int status;
            numberLotSizing = 0;
            delete[] lotsize;
#ifndef CBC_OTHER_SOLVER
            ClpSimplex *lpSolver = clpSolver->getModelPtr();
            if (!gmpl) {
               status =
                  clpSolver->readMps(fileName.c_str(), keepImportNames != 0,
                                     allowImportErrors != 0);
            } else if (gmpl > 0) {
#if defined(CLP_HAS_GLPK) && defined(COINUTILS_HAS_GLPK)
               status = lpSolver->readGMPL(
                    fileName.c_str(), (gmpl == 2) ? gmplData.c_str() : NULL,
                    keepImportNames != 0, &coin_glp_tran, &coin_glp_prob);
#endif
              } else {
#ifdef KILL_ZERO_READLP
               status = clpSolver->readLp(fileName.c_str(),
                                          lpSolver->getSmallElementValue());
#else
               status = clpSolver->readLp(fileName.c_str(), 1.0e-12);
#endif
            }
#else
            status = clpSolver->readMps(fileName.c_str(), "");
#endif
            if (!status || (status > 0 && allowImportErrors)) {
#ifndef CBC_OTHER_SOLVER
               if (keepImportNames) {
                  lengthName = lpSolver->lengthNames();
                  rowNames = *(lpSolver->rowNames());
                  columnNames = *(lpSolver->columnNames());
               } else {
                  lengthName = 0;
               }
               // really just for testing
               double objScale = clpParameters[ClpParam::OBJSCALE2]->dblVal();
               if (objScale != 1.0) {
                  int iColumn;
                  int numberColumns = lpSolver->numberColumns();
                  double *dualColumnSolution = lpSolver->dualColumnSolution();
                  ClpObjective *obj = lpSolver->objectiveAsObject();
                  assert(dynamic_cast<ClpLinearObjective *>(obj));
                  double offset;
                  double *objective = obj->gradient(NULL, NULL, offset, true);
                  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                     dualColumnSolution[iColumn] *= objScale;
                     objective[iColumn] *= objScale;
                  }
                  int iRow;
                  int numberRows = lpSolver->numberRows();
                  double *dualRowSolution = lpSolver->dualRowSolution();
                  for (iRow = 0; iRow < numberRows; iRow++)
                     dualRowSolution[iRow] *= objScale;
                  lpSolver->setObjectiveOffset(objScale *
                                               lpSolver->objectiveOffset());
               }
               goodModel = true;
               // sets to all slack (not necessary?)
               lpSolver->createStatus();
               // See if sos
               if (clpSolver->numberSOS()) {
                  // SOS
                  numberSOS = clpSolver->numberSOS();
                  const CoinSet *setInfo = clpSolver->setInfo();
                  sosStart = reinterpret_cast<int *>(
                      malloc((numberSOS + 1) * sizeof(int)));
                  sosType = reinterpret_cast<char *>(
                      malloc(numberSOS * sizeof(char)));
                  const double *lower = clpSolver->getColLower();
                  const double *upper = clpSolver->getColUpper();
                  int i;
                  int nTotal = 0;
                  sosStart[0] = 0;
                  for (i = 0; i < numberSOS; i++) {
                     int type = setInfo[i].setType();
                     int n = setInfo[i].numberEntries();
                     sosType[i] = static_cast<char>(type);
                     nTotal += n;
                     sosStart[i + 1] = nTotal;
                  }
                  sosIndices =
                     reinterpret_cast<int *>(malloc(nTotal * sizeof(int)));
                  sosReference = reinterpret_cast<double *>(
                      malloc(nTotal * sizeof(double)));
                  for (i = 0; i < numberSOS; i++) {
                    int n = setInfo[i].numberEntries();
                    const int *which = setInfo[i].which();
                    const double *weights = setInfo[i].weights();
                    int base = sosStart[i];
                    for (int j = 0; j < n; j++) {
                      int k = which[j];
                      // don't allow free
                      if (upper[k] > 1.0e15)
                        clpSolver->setColUpper(k, 1.0e15);
                      if (lower[k] < -1.0e15)
                        clpSolver->setColLower(k, -1.0e15);
                      sosIndices[j + base] = k;
                      sosReference[j + base] =
                          weights ? weights[j] : static_cast<double>(j);
                    }
                  }
               }
               // make sure integer
               // also deal with semi-continuous
               int numberColumns = lpSolver->numberColumns();
               int i;
               for (i = 0; i < numberColumns; i++) {
                  if (clpSolver->integerType(i) > 2)
                     break;
                  if (lpSolver->isInteger(i))
                     clpSolver->setInteger(i);
               }
               if (i < numberColumns) {
                  // semi-continuous
                  clpSolver->setSpecialOptions(clpSolver->specialOptions() |
                                               8388608);
                  int iStart = i;
                  for (i = iStart; i < numberColumns; i++) {
                     if (clpSolver->integerType(i) > 2)
                        numberLotSizing++;
                  }
                  lotsize = new lotStruct[numberLotSizing];
                  numberLotSizing = 0;
                  const double *lower = clpSolver->getColLower();
                  const double *upper = clpSolver->getColUpper();
                  for (i = iStart; i < numberColumns; i++) {
                     if (clpSolver->integerType(i) > 2) {
                        int iType = clpSolver->integerType(i) - 3;
                        if (!iType)
                           clpSolver->setContinuous(i);
                        else
                           clpSolver->setInteger(i);
                        lotsize[numberLotSizing].column = i;
                        lotsize[numberLotSizing].high = upper[i];
                        if (lower[i]) {
                           lotsize[numberLotSizing++].low = lower[i];
                           clpSolver->setColLower(i, 0.0);
                        } else {
                           lotsize[numberLotSizing++].low = 1.0;
                        }
                     }
                  }
               }
#else
               lengthName = 0;
               goodModel = true;
#endif
               time2 = CoinCpuTime();
               totalTime += time2 - time1;
               time1 = time2;
               // Go to canned file if just input file
               if (inputQueue.empty()) {
                  // only if ends .mps
                  char *find =
                     const_cast<char *>(strstr(fileName.c_str(), ".mps"));
                  if (find && find[4] == '\0') {
                     find[1] = 'p';
                     find[2] = 'a';
                     find[3] = 'r';
                     std::ifstream ifs(fileName.c_str());
                     if (ifs.is_open()) {
                        while (!inputQueue.empty()){
                           inputQueue.pop_front();
                        }
                        CoinParamUtils::readFromStream(inputQueue, ifs);
                     }
                  }
               }
            } else {
               // errors
               buffer.str("");
               buffer << "There were " << status << " errors on input";
               printGeneralMessage(model_, buffer.str());
            }
          } break;
          case CbcParam::EXPORT: {
            if (!goodModel) {
              printGeneralMessage(model_, "** Current model not valid\n");
              continue;
            }
            cbcParam->readValue(inputQueue, fileName, &message);
            CoinParamUtils::processFile(fileName,
                                 parameters[CbcParam::DIRECTORY]->dirName());
            if (fileName == ""){
               fileName = parameters[CbcParam::EXPORTFILE]->fileName();
            }else{
               parameters[CbcParam::EXPORTFILE]->setFileName(fileName);
            }
            // If presolve on then save presolved
            bool deleteModel2 = false;
            ClpSimplex *model2 = lpSolver;
            if (dualize && dualize < 3) {
               model2 =
                  static_cast<ClpSimplexOther *>(model2)->dualOfModel();
               buffer.str("");
               buffer << "Dual of model has " << model2->numberRows() << " rows and "
                      <<  model2->numberColumns() << " columns";
               printGeneralMessage(model_, buffer.str());
               model2->setOptimizationDirection(1.0);
            }
#ifndef CBC_OTHER_SOLVER
            if (info && info->numberSos && doSOS && statusUserFunction_[0]) {
               // SOS
               numberSOS = info->numberSos;
               sosStart = info->sosStart;
               sosIndices = info->sosIndices;
               sosReference = info->sosReference;
               clpSolver->setSOSData(numberSOS, info->sosType, sosStart,
                                     sosIndices, sosReference);
            }
            numberSOS = clpSolver->numberSOS();
            if (numberSOS || lotsize)
               preSolve = false;
#endif
            if (preSolve) {
               ClpPresolve pinfo;
               int presolveOptions2 = presolveOptions & ~0x40000000;
               if ((presolveOptions2 & 0xffff) != 0)
                  pinfo.setPresolveActions(presolveOptions2);
               if ((printOptions & 1) != 0)
                  pinfo.statistics();
               double presolveTolerance =
                  clpParameters[ClpParam::PRESOLVETOLERANCE]->dblVal();
               model2 = pinfo.presolvedModel(*lpSolver, presolveTolerance,
                                             true, preSolve);
               if (model2) {
                  printf("Saving presolved model on %s\n", fileName.c_str());
                  deleteModel2 = true;
               } else {
                  printf("Presolved model looks infeasible - saving original "
                         "on %s\n",
                         fileName.c_str());
                  deleteModel2 = false;
                  model2 = lpSolver;
               }
               // see if extension lp
               bool writeLp = false;
               {
                  int lengthName = strlen(fileName.c_str());
                  if (lengthName > 3 &&
                      !strcmp(fileName.c_str() + lengthName - 3, ".lp"))
                     writeLp = true;
               }
               if (!writeLp) {
                  model2->writeMps(fileName.c_str(), (outputFormat - 1) / 2,
                                   1 + ((outputFormat - 1) & 1));
               } else {
                  fp = fopen(fileName.c_str(), "w");
                  if (!fp) {
                     buffer.str("");
                     buffer << "Unable to open file " << fileName.c_str();
                     printGeneralMessage(model_, buffer.str());
                     continue;
                  }
                  OsiClpSolverInterface solver(model2);
                  solver.writeLp(fp, 1.0e-12, 10, 10);
		  fclose(fp);
               }
               if (deleteModel2)
                  delete model2;
            } else {
               printf("Saving model on %s\n", fileName.c_str());
#ifdef COIN_HAS_LINK
               OsiSolverLink *linkSolver =
                  dynamic_cast<OsiSolverLink *>(clpSolver);
               if (!linkSolver || !linkSolver->quadraticModel()) {
#endif
                  // Convert names
                  int iRow;
                  int numberRows = model2->numberRows();
                  int iColumn;
                  int numberColumns = model2->numberColumns();
                  
                  char **rowNames = NULL;
                  char **columnNames = NULL;
                  if (model2->lengthNames()) {
                     rowNames = new char *[numberRows];
                     for (iRow = 0; iRow < numberRows; iRow++) {
                        rowNames[iRow] =
                           CoinStrdup(model2->rowName(iRow).c_str());
                     }
                     
                     columnNames = new char *[numberColumns];
                     for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                        columnNames[iColumn] =
                           CoinStrdup(model2->columnName(iColumn).c_str());
                     }
                  }
                  // see if extension lp
                  bool writeLp = false;
                  {
                     int lengthName = strlen(fileName.c_str());
                     if (lengthName > 3 &&
                         !strcmp(fileName.c_str() + lengthName - 3, ".lp"))
                        writeLp = true;
                  }
                  if (lotsize) {
                     for (int i = 0; i < numberLotSizing; i++) {
                        int iColumn = lotsize[i].column;
                        double low = lotsize[i].low;
                        if (low != 1.0)
                           clpSolver->setColLower(iColumn, low);
                        int type;
                        if (clpSolver->isInteger(iColumn))
                           type = 4;
                        else
                           type = 3;
                        clpSolver->setColumnType(iColumn, type);
                     }
                  }
                  if (!writeLp) {
                     remove(fileName.c_str());
                     // model_.addSOSEtcToSolver();
                     clpSolver->writeMpsNative(
                          fileName.c_str(), const_cast<const char **>(rowNames),
                          const_cast<const char **>(columnNames),
                          (outputFormat - 1) / 2, 1 + ((outputFormat - 1) & 1));
                  } else {
                     fp = fopen(fileName.c_str(), "w");
                     if (!fp) {
                        buffer.str("");
                        buffer << "Unable to open file " << fileName.c_str();
                        printGeneralMessage(model_, buffer.str());
                        continue;
                     }
                     clpSolver->writeLp(fp, 1.0e-12, 10, 10);
                  }
                  if (rowNames) {
                     for (iRow = 0; iRow < numberRows; iRow++) {
                        free(rowNames[iRow]);
                     }
                     delete[] rowNames;
                     for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                        free(columnNames[iColumn]);
                     }
                     delete[] columnNames;
                  }
                  if (lotsize) {
                     for (int i = 0; i < numberLotSizing; i++) {
                        int iColumn = lotsize[i].column;
                        int itype = clpSolver->integerType(iColumn);
                        clpSolver->setColLower(iColumn, 0.0);
                        if (itype == 3)
                           clpSolver->setContinuous(iColumn);
                        else
                           clpSolver->setInteger(iColumn);
                     }
                  }
#ifdef COIN_HAS_LINK
               } else {
                  linkSolver->quadraticModel()->writeMps(
                        fileName.c_str(), (outputFormat - 1) / 2,
                        1 + ((outputFormat - 1) & 1));
               }
#endif
            }
            time2 = CoinCpuTime();
            totalTime += time2 - time1;
            time1 = time2;
          } break;
           case CbcParam::READPRIORITIES:{
            if (!goodModel){
              printGeneralMessage(model_, "** Current model not valid\n");
              continue;
            }
            cbcParam->readValue(inputQueue, fileName, &message);
            CoinParamUtils::processFile(fileName,
                                 parameters[CbcParam::DIRECTORY]->dirName(),
                                        &canOpen);
            if (!canOpen) {
               buffer.str("");
               buffer << "Unable to open file " << fileName.c_str();
               printGeneralMessage(model_, buffer.str());
               continue;
            }
            if (fileName == ""){
               fileName = parameters[CbcParam::PRIORITYFILE]->fileName();
            }else{
               parameters[CbcParam::PRIORITYFILE]->setFileName(fileName);
            }
            fp = fopen(fileName.c_str(), "r");
            std::string headings[] = {"name",     "number", "direction",
                                      "priority", "up",     "down",
                                      "solution", "priin"};
            int got[] = {-1, -1, -1, -1, -1, -1, -1, -1};
            int order[8];
            bool useMasks = false;
            if (strstr(fileName.c_str(), "mask_")) {
               // look more closely
               const char *name = fileName.c_str();
               int length = strlen(name);
               for (int i = length - 1; i >= 0; i--) {
                  if (name[i] == dirsep) {
                     name += i + 1;
                     break;
                  }
               }
               useMasks = !strncmp(name, "mask_", 5);
            }
            assert(sizeof(got) == sizeof(order));
            int nAcross = 0;
            char line[1000];
            int numberColumns = lpSolver->numberColumns();
            if (!fgets(line, 1000, fp)) {
               std::cout << "Odd file " << fileName << std::endl;
            } else {
               char *pos = line;
               char *put = line;
               while (*pos >= ' ' && *pos != '\n') {
                  if (*pos != ' ' && *pos != '\t') {
                     *put = static_cast<char>(tolower(*pos));
                     put++;
                  }
                  pos++;
               }
               *put = '\0';
               pos = line;
               int i;
               bool good = true;
               while (pos) {
                  char *comma = strchr(pos, ',');
                  if (comma)
                     *comma = '\0';
                  for (i = 0; i < static_cast<int>(sizeof(got) / sizeof(int));
                       i++) {
                     if (headings[i] == pos) {
                        if (got[i] < 0) {
                           order[nAcross] = i;
                           got[i] = nAcross++;
                        } else {
                           // duplicate
                           good = false;
                        }
                        break;
                     }
                  }
                  if (i == static_cast<int>(sizeof(got) / sizeof(int)))
                     good = false;
                  if (comma) {
                     *comma = ',';
                     pos = comma + 1;
                  } else {
                     break;
                  }
               }
               if (got[0] < 0 && got[1] < 0)
                  good = false;
               if (got[0] >= 0 && got[1] >= 0)
                  good = false;
               if (got[0] >= 0 && !lpSolver->lengthNames())
                  good = false;
               int numberFields = 99;
               if (good && (strstr(fileName.c_str(), ".mst") ||
                            strstr(fileName.c_str(), ".MST") ||
                            strstr(fileName.c_str(), ".csv"))) {
                  numberFields = 0;
                  for (i = 2; i < static_cast<int>(sizeof(got) / sizeof(int));
                       i++) {
                     if (got[i] >= 0)
                        numberFields++;
                  }
                  if (!numberFields) {
                     // Like Cplex format
                     order[nAcross] = 6;
                     got[6] = nAcross++;
                  }
               }
               if (good) {
                  char **columnNames = new char *[numberColumns];
                  // pseudoDown = NULL;
                  // pseudoUp = NULL;
                  // branchDirection = NULL;
                  // if (got[5]!=-1)
                  pseudoDown = reinterpret_cast<double *>(
                        malloc(numberColumns * sizeof(double)));
                  // if (got[4]!=-1)
                  pseudoUp = reinterpret_cast<double *>(
                        malloc(numberColumns * sizeof(double)));
                  // if (got[2]!=-1)
                  branchDirection = reinterpret_cast<int *>(
                        malloc(numberColumns * sizeof(int)));
                  priorities = reinterpret_cast<int *>(
                        malloc(numberColumns * sizeof(int)));
                  free(solutionIn);
                  solutionIn = NULL;
                  free(prioritiesIn);
                  prioritiesIn = NULL;
                  int iColumn;
                  if (got[6] >= 0) {
                     solutionIn = reinterpret_cast<double *>(
                          malloc(numberColumns * sizeof(double)));
                     for (iColumn = 0; iColumn < numberColumns; iColumn++)
                        solutionIn[iColumn] = -COIN_DBL_MAX;
                  }
                  if (got[7] >= 0 || !numberFields) {
                     prioritiesIn = reinterpret_cast<int *>(
                          malloc(numberColumns * sizeof(int)));
                     for (iColumn = 0; iColumn < numberColumns; iColumn++)
                        prioritiesIn[iColumn] = 10000;
                  }
                  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                     columnNames[iColumn] =
                        CoinStrdup(lpSolver->columnName(iColumn).c_str());
                     // if (got[5]!=-1)
                     pseudoDown[iColumn] = 0.0;
                     // if (got[4]!=-1)
                     pseudoUp[iColumn] = 0.0;
                     // if (got[2]!=-1)
                     branchDirection[iColumn] = 0;
                     priorities[iColumn] = useMasks ? -123456789 : 0;
                  }
                  int nBadPseudo = 0;
                  int nBadDir = 0;
                  int nBadPri = 0;
                  int nBadName = 0;
                  int nBadLine = 0;
                  int nLine = 0;
                  iColumn = -1;
                  int lowestPriority = -COIN_INT_MAX;
                  bool needCard = true;
                  while (!needCard || fgets(line, 1000, fp)) {
                     if (!strncmp(line, "ENDATA", 6) ||
                         !strncmp(line, "endata", 6))
                        break;
                     nLine++;
                     if (!useMasks)
                        iColumn = -1;
                     else
                        needCard = false;
                     double up = 0.0;
                     double down = 0.0;
                     int pri = 0;
                     int dir = 0;
                     double solValue = COIN_DBL_MAX;
                     int priValue = 1000000;
                     char *pos = line;
                     char *put = line;
                     if (!numberFields) {
                        // put in ,
                        for (i = 4; i < 100; i++) {
                           if (line[i] == ' ' || line[i] == '\t') {
                              line[i] = ',';
                              break;
                           }
                        }
                     }
                     while (*pos >= ' ' && *pos != '\n') {
                        if (*pos != ' ' && *pos != '\t') {
                           *put = *pos;
                           put++;
                        }
                        pos++;
                     }
                     *put = '\0';
                     pos = line;
                     for (int i = 0; i < nAcross; i++) {
                        char *comma = strchr(pos, ',');
                        if (comma) {
                           *comma = '\0';
                        } else if (i < nAcross - 1) {
                           nBadLine++;
                           break;
                        }
                        switch (order[i]) {
                           // name
                         case 0:
                           iColumn++;
                           for (; iColumn < numberColumns; iColumn++) {
                              if (priorities[iColumn] != -123456789) {
                                 if (!strcmp(columnNames[iColumn], pos))
                                    break;
                              } else {
                                 // mask (at present ? and trailing *)
                                 const char *name = columnNames[iColumn];
                                 int length = strlen(name);
                                 int lengthMask = strlen(pos);
                                 bool asterisk = pos[lengthMask - 1] == '*';
                                 if (asterisk)
                                    length = lengthMask - 1;
                                 int i;
                                 for (i = 0; i < length; i++) {
                                    if (name[i] != pos[i]) {
                                       if (pos[i] != '?')
                                          break;
                                    }
                                 }
                                 if (i == length)
                                    break;
                              }
                           }
                           if (iColumn == numberColumns) {
                              iColumn = -1;
                              needCard = true;
                           }
                           break;
                           // number
                         case 1:
                           iColumn = atoi(pos);
                           if (iColumn < 0 || iColumn >= numberColumns)
                              iColumn = -1;
                           break;
                           // direction
                         case 2:
                           if (*pos == 'D')
                              dir = -1;
                           else if (*pos == 'U')
                              dir = 1;
                           else if (*pos == 'N')
                              dir = 0;
                           else if (*pos == '1' && *(pos + 1) == '\0')
                              dir = 1;
                           else if (*pos == '0' && *(pos + 1) == '\0')
                              dir = 0;
                           else if (*pos == '1' && *(pos + 1) == '1' &&
                                    *(pos + 2) == '\0')
                              dir = -1;
                           else
                              dir = -2; // bad
                           break;
                           // priority
                         case 3:
                           pri = atoi(pos);
                           lowestPriority = std::max(lowestPriority, pri);
                           break;
                           // up
                         case 4:
                           up = atof(pos);
                           break;
                           // down
                         case 5:
                           down = atof(pos);
                           break;
                           // sol value
                         case 6:
                           solValue = atof(pos);
                           break;
                           // priority in value
                         case 7:
                           priValue = atoi(pos);
                           break;
                        }
                        if (comma) {
                           *comma = ',';
                           pos = comma + 1;
                        }
                     }
                     if (iColumn >= 0) {
                        if (down < 0.0) {
                           nBadPseudo++;
                           down = 0.0;
                        }
                        if (up < 0.0) {
                           nBadPseudo++;
                           up = 0.0;
                        }
                        if (!up)
                           up = down;
                        if (!down)
                           down = up;
                        if (dir < -1 || dir > 1) {
                           nBadDir++;
                           dir = 0;
                        }
                        if (pri < 0) {
                           nBadPri++;
                           pri = 0;
                        }
                        // if (got[5]!=-1)
                        pseudoDown[iColumn] = down;
                        // if (got[4]!=-1)
                        pseudoUp[iColumn] = up;
                        // if (got[2]!=-1)
                        branchDirection[iColumn] = dir;
                        priorities[iColumn] = pri;
                        if (solValue != COIN_DBL_MAX) {
                          assert(solutionIn);
                          solutionIn[iColumn] = solValue;
                        }
                        if (priValue != 1000000) {
                          assert(prioritiesIn);
                          prioritiesIn[iColumn] = priValue;
                        }
                      } else if (!useMasks) {
                        nBadName++;
                      }
                    }
                  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                     if (priorities[iColumn] == -123456789)
                        priorities[iColumn] = lowestPriority + 1;
                  }
                  if (!parameters.noPrinting()) {
                     printf("%d fields and %d records", nAcross, nLine);
                     if (nBadPseudo)
                        printf(" %d bad pseudo costs", nBadPseudo);
                     if (nBadDir)
                        printf(" %d bad directions", nBadDir);
                     if (nBadPri)
                        printf(" %d bad priorities", nBadPri);
                     if (nBadName)
                        printf(" ** %d records did not match on name/sequence",
                               nBadName);
                     printf("\n");
                  }
                  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                     free(columnNames[iColumn]);
                  }
                  delete[] columnNames;
               } else {
                  std::cout << "Duplicate or unknown keyword - or "
                     "name/number fields wrong"
                            << line << std::endl;
               }
            }
            fclose(fp);
           } break;
          case CbcParam::READMIPSTART:
          case CbcParam::READSOL: {
            if (!goodModel){
               printGeneralWarning(model_, "** Current model not valid\n");
               continue;
            }
            cbcParam->readValue(inputQueue, fileName, &message);
            CoinParamUtils::processFile(fileName,
                                 parameters[CbcParam::DIRECTORY]->dirName(),
                                        &canOpen);
            if (fileName == ""){
               fileName = parameters[CbcParam::MIPSTARTFILE]->fileName();
            }else{
               parameters[CbcParam::MIPSTARTFILE]->setFileName(fileName);
            }
            mipStartFile = fileName;
            buffer.str("");
            if (canOpen){
               buffer << "opening mipstart file " << fileName.c_str();
               printGeneralMessage(model_, buffer.str());
            } else {
               buffer << "Unable to open file " << fileName.c_str();
               printGeneralMessage(model_, buffer.str());
               continue;
            }              
            double msObj;
            
            CbcMipStartIO::read(model_.solver(), fileName.c_str(), mipStart,
                                msObj, model_.messageHandler(),
                                model_.messagesPointer());
            // copy to before preprocess if has .before.
            if (strstr(fileName.c_str(), ".before.")) {
               mipStartBefore = mipStart;
               buffer.str("");
               buffer << "file " <<  fileName.c_str() << " will be used before preprocessing.";
               printGeneralMessage(model_, buffer.str());
            }
          } break;
           case CbcParam::DEBUG:{
            if (!goodModel){
               printGeneralWarning(model_, "** Current model not valid\n");
               continue;
            }
            delete[] debugValues;
            debugValues = NULL;
            cbcParam->readValue(inputQueue, fileName, &message);
            if (fileName == "create" || fileName == "createAfterPre") {
               parameters[CbcParam::DEBUGFILE]->setFileName(fileName);
               printGeneralMessage(model_, "Will create a debug file so "
                                   "this run should be a good one");
               break;
            } else if (fileName == "unitTest") {
               printGeneralMessage(model_, "debug will be done using file "
                                   "name of model");
               parameters[CbcParam::DEBUGFILE]->setFileName(fileName);
               break;
            }
            CoinParamUtils::processFile(fileName,
                                 parameters[CbcParam::DIRECTORY]->dirName(),
                                        &canOpen);
            if (!canOpen){
               buffer.str("");
               buffer << "Unable to open file " << fileName.c_str();
               printGeneralMessage(model_, buffer.str());
               continue;
            }
            if (fileName == ""){
               fileName = parameters[CbcParam::DEBUGFILE]->fileName();
            }else{
               parameters[CbcParam::DEBUGFILE]->setFileName(fileName);
            }
            fp = fopen(fileName.c_str(), "rb");
            // can open - lets go for it
            int numRows;
            double obj;
            size_t nRead;
            nRead = fread(&numRows, sizeof(int), 1, fp);
            if (nRead != 1)
               throw("Error in fread");
            nRead = fread(&numberDebugValues, sizeof(int), 1, fp);
            if (nRead != 1)
               throw("Error in fread");
            nRead = fread(&obj, sizeof(double), 1, fp);
            if (nRead != 1)
               throw("Error in fread");
            debugValues = new double[numberDebugValues + numRows];
            nRead = fread(debugValues, sizeof(double), numRows, fp);
            if (nRead != static_cast<size_t>(numRows))
               throw("Error in fread");
            nRead = fread(debugValues, sizeof(double), numRows, fp);
            if (nRead != static_cast<size_t>(numRows))
               throw("Error in fread");
            nRead =
               fread(debugValues, sizeof(double), numberDebugValues, fp);
            if (nRead != static_cast<size_t>(numberDebugValues))
               throw("Error in fread");
            printf("%d doubles read into debugValues\n", numberDebugValues);
#if DEBUG_PREPROCESS > 1
            debugSolution = debugValues;
            debugNumberColumns = numberDebugValues;
#endif
            if (numberDebugValues < 200) {
               for (int i = 0; i < numberDebugValues; i++) {
                  if (clpSolver->isInteger(i) && debugValues[i])
                     printf("%d %g\n", i, debugValues[i]);
               }
            }
           } break;
          case CbcParam::PRINTMASK:
            if ((status = cbcParam->readValue(inputQueue, field, &message))){
               printGeneralMessage(model_, message);
               continue;
            }
            if (cbcParam->setVal(field, &message)){
               printGeneralMessage(model_, message);
               continue;
            }
            break;
            
          case CbcParam::WRITEMODEL:{
            cbcParam->readValue(inputQueue, fileName, &message);
            CoinParamUtils::processFile(fileName,
                                 parameters[CbcParam::DIRECTORY]->dirName());
            if (fileName == ""){
               fileName = parameters[CbcParam::MODELFILE]->fileName();
            }else{
               parameters[CbcParam::MODELFILE]->setFileName(fileName);
            }
            int status;
            // If presolve on then save presolved
            bool deleteModel2 = false;
            ClpSimplex *model2 = lpSolver;
            if (preSolve) {
               ClpPresolve pinfo;
               double presolveTolerance =
                  parameters[ClpParam::PRESOLVETOLERANCE]->dblVal();
               model2 = pinfo.presolvedModel(*lpSolver, presolveTolerance,
                                             false, preSolve);
               if (model2) {
                  printf("Saving presolved model on %s\n",
                         fileName.c_str());
                  deleteModel2 = true;
               } else {
                  printf("Presolved model looks infeasible - saving original on %s\n",
                         fileName.c_str());
                  deleteModel2 = false;
                  model2 = lpSolver;
               }
            } else {
               printf("Saving model on %s\n",
                      fileName.c_str());
            }
            status = model2->saveModel(fileName.c_str());
            if (deleteModel2)
               delete model2;
            if (!status) {
               goodModel = true;
               time2 = CoinCpuTime();
               totalTime += time2 - time1;
               time1 = time2;
            } else {
               // errors
               printGeneralMessage(model_, "There were errors on output");
            }
          } break;
            
          case CbcParam::READMODEL: {

            cbcParam->readValue(inputQueue, fileName, &message);
            CoinParamUtils::processFile(fileName,
                                 parameters[CbcParam::DIRECTORY]->dirName(),
                                        &canOpen);
            if (!canOpen) {
               buffer.str("");
               buffer << "Unable to open file " << fileName.c_str();
               printGeneralMessage(model_, buffer.str());
               continue;
            }
            if (fileName == ""){
               fileName = parameters[CbcParam::MODELFILE]->fileName();
            }else{
               parameters[CbcParam::MODELFILE]->setFileName(fileName);
            }
            int status = lpSolver->restoreModel(fileName.c_str());
            if (!status) {
               goodModel = true;
               time2 = CoinCpuTime();
               totalTime += time2 - time1;
               time1 = time2;
            } else {
               // errors
               buffer.str("");
               buffer << "There were errors on input";
               printGeneralMessage(model_, buffer.str());
            }
          } break;
          case CbcParam::MAXIMIZE:
	    lpSolver->setOptimizationDirection(-1);
            break;
          case CbcParam::MINIMIZE:
            lpSolver->setOptimizationDirection(1);
            break;
           case CbcParam::REVERSE:{
            if (!goodModel){
               continue;
            }
            int iColumn;
            int numberColumns = lpSolver->numberColumns();
            double *dualColumnSolution = lpSolver->dualColumnSolution();
            ClpObjective *obj = lpSolver->objectiveAsObject();
            assert(dynamic_cast<ClpLinearObjective *>(obj));
            double offset;
            double *objective = obj->gradient(NULL, NULL, offset, true);
            for (iColumn = 0; iColumn < numberColumns; iColumn++) {
               dualColumnSolution[iColumn] = dualColumnSolution[iColumn];
               objective[iColumn] = -objective[iColumn];
            }
            int iRow;
            int numberRows = lpSolver->numberRows();
            double *dualRowSolution = lpSolver->dualRowSolution();
            for (iRow = 0; iRow < numberRows; iRow++)
               dualRowSolution[iRow] = dualRowSolution[iRow];
           } break;
          case CbcParam::STDIN:
            interactiveMode = true;
            while (!inputQueue.empty())
               inputQueue.pop_front();
            break;
          case CbcParam::UNITTEST: {
            int returnCode = -1;
            std::string dirMiplib = parameters[CbcParam::DIRMIPLIB]->dirName();
            std::string dirSample = parameters[CbcParam::DIRSAMPLE]->dirName();
            if (dirMiplib != ""){
              returnCode = CbcClpUnitTest(model_, dirMiplib, -3, NULL,
					  saveInputQueue,callBack, parameters);
            } else if (dirSample != ""){
              returnCode = CbcClpUnitTest(model_, dirSample, -2, NULL,
					  saveInputQueue,callBack, parameters);
            } else {
               buffer.str("");
               buffer << "No directory specified for input files"
                      << "Please use dirMiplib or dirSample to set directory name";
               printGeneralMessage(model_, buffer.str());
            }
            babModel_ = NULL;
            return returnCode;
          }
          case CbcParam::USERCBC: {
#ifdef USER_HAS_FAKE_CBC
            // Replace the sample code by whatever you want
            if (!goodModel){
               printGeneralWarning(model_, "** Current model not valid\n");
               continue;
            }
            // Way of using an existing piece of code
            OsiClpSolverInterface *clpSolver =
               dynamic_cast<OsiClpSolverInterface *>(model_.solver());
            ClpSimplex *lpSolver = clpSolver->getModelPtr();
            // set time from integer model
            double timeToGo = model_.getMaximumSeconds();
            lpSolver->setMaximumSeconds(timeToGo);
            fakeMain(*lpSolver, *clpSolver, model);
            // My actual usage has objective only in clpSolver
            double objectiveValue = clpSolver->getObjValue();
            int iStat = lpSolver->status();
            int iStat2 = lpSolver->secondaryStatus();
            // make sure solution back in correct place
            clpSolver =
               dynamic_cast<OsiClpSolverInterface *>(model_.solver());
            lpSolver = clpSolver->getModelPtr();
            if (info && statusUserFunction_[0]) {
               int n = clpSolver->getNumCols();
               double value = objectiveValue * lpSolver->getObjSense();
               char buf[300];
               int pos = 0;
               std::string minor[] = {"",     "", "gap",       "nodes",
                                      "time", "", "solutions", "user ctrl-c"};
               if (iStat == 0) {
                  if (objectiveValue < 1.0e40) {
                     pos += sprintf(buf + pos, "optimal,");
                  } else {
                     // infeasible
                     iStat = 1;
                     pos += sprintf(buf + pos, "infeasible,");
                  }
               } else if (iStat == 1) {
                  if (iStat2 != 6)
                     iStat = 3;
                  else
                     iStat = 4;
                  pos += sprintf(buf + pos, "stopped on %s,",
                                 minor[iStat2].c_str());
               } else if (iStat == 2) {
                  iStat = 7;
                  pos += sprintf(buf + pos, "stopped on difficulties,");
                } else if (iStat == 5) {
                  iStat = 3;
                  pos += sprintf(buf + pos, "stopped on ctrl-c,");
                } else if (iStat == 6) {
                  // bab infeasible
                  pos += sprintf(buf + pos, "integer infeasible,");
                  iStat = 1;
                } else {
                  pos += sprintf(buf + pos, "status unknown,");
                  iStat = 6;
                }
                info->problemStatus = iStat;
                info->objValue = value;
                if (objectiveValue < 1.0e40)
                  pos += sprintf(buf + pos, " objective %.*g", ampl_obj_prec(),
                                 value);
                sprintf(buf + pos, "\n%d nodes, %d iterations",
                        model_.getNodeCount(), model_.getIterationCount());
                if (objectiveValue < 1.0e50) {
                  free(info->primalSolution);
                  info->primalSolution = (double *)malloc(n * sizeof(double));
                  CoinCopyN(lpSolver->primalColumnSolution(), n,
                            info->primalSolution);
                  int numberRows = lpSolver->numberRows();
                  free(info->dualSolution);
                  info->dualSolution =
                      (double *)malloc(numberRows * sizeof(double));
                  CoinCopyN(lpSolver->dualRowSolution(), numberRows,
                            info->dualSolution);
                } else {
                  info->primalSolution = NULL;
                  info->dualSolution = NULL;
                }
                // put buffer into info
                strcpy(info->buffer, buf);
              }
#endif
          } break;
          case CbcParam::HELP:
            std::cout << "Cbc version " << CBC_VERSION << ", build " << __DATE__
                      << std::endl;
#if 0
            std::cout << "Non default values:-" << std::endl;
            std::cout << "Perturbation " << lpSolver->perturbation()
                      << " (default 100)" << std::endl;
	    std::cout << CoinParamUtils::printString("Presolve being done with 5 passes\n\
Dual steepest edge steep/partial on matrix shape and factorization density\n\
Clpnnnn taken out of messages\n\
If Factorization frequency default then done on size of matrix\n\n\
(-)unitTest, (-)netlib or (-)netlibp will do standard tests\n\n\
You can switch to interactive mode at any time so\n\
clp watson.mps -scaling off -primalsimplex\nis the same as\n\
clp watson.mps -\nscaling off\nprimalsimplex");
	    std::endl;
#endif
	    std::cout << "Use ? or ??? to see list of commands" << std::endl;
	    std::cout << "Set allcommands to all to see less common commands" << std::endl;
            break;
          case CbcParam::WRITESTATS: {
            cbcParam->readValue(inputQueue, fileName, &message);
            CoinParamUtils::processFile(fileName,
                                 parameters[CbcParam::DIRECTORY]->dirName());
            if (fileName == ""){
               fileName = parameters[CbcParam::CSVSTATSFILE]->fileName();
            }else{
               parameters[CbcParam::CSVSTATSFILE]->setFileName(fileName);
            }
            int state = 0;
            // FIXME: This needs to fixed up to use modern C++ and to use the inputQueue properly
            char cbuffer[1000];
            fp = fopen(fileName.c_str(), "r");
            if (fp) {
              // file already there
              state = 1;
              char *getBuffer = fgets(cbuffer, 1000, fp);
              if (getBuffer) {
                // assume header there
                state = 2;
              }
              fclose(fp);
            }
            fp = fopen(fileName.c_str(), "a");
            if (!fp) {
               buffer.str("");
               buffer << "Unable to open file " << fileName.c_str();
               printGeneralMessage(model_, buffer.str());
               continue;
            }
            // can open - lets go for it
            // first header if needed
            if (state != 2) {
               fprintf(fp, "Name,result,time,sys,elapsed,objective,continuous,"
                       "tightened,cut_time,nodes,iterations,rows,columns,"
                       "processed_rows,processed_columns");
               for (int i = 0; i < statistics_number_generators; i++)
                  fprintf(fp, ",%s", statistics_name_generators[i]);
               fprintf(fp, ",runtime_options");
               fprintf(fp, "\n");
            }
            strcpy(cbuffer, inputQueue[0].c_str());
            char *slash = cbuffer;
            for (int i = 0; i < static_cast<int>(strlen(cbuffer)); i++) {
               if (cbuffer[i] == '/' || cbuffer[i] == '\\')
                  slash = cbuffer + i + 1;
            }
            fprintf(fp,
                    "%s,%s,%.2f,%.2f,%.2f,%.16g,%g,%g,%.2f,%d,%d,%d,%d,%d,%d",
                    slash, statistics_result.c_str(), statistics_seconds,
                    statistics_sys_seconds, statistics_elapsed_seconds,
                    statistics_obj, statistics_continuous, statistics_tighter,
                    statistics_cut_time, statistics_nodes,
                    statistics_iterations, statistics_nrows, statistics_ncols,
                    statistics_nprocessedrows, statistics_nprocessedcols);
            for (int i = 0; i < statistics_number_generators; i++)
               fprintf(fp, ",%d",
                       statistics_number_cuts != NULL
                       ? statistics_number_cuts[i]
                       : 0);
            fprintf(fp, ",");
            for (int i = 0; i < inputQueue.size(); i++) {
               if (strstr(inputQueue[i].c_str(), ".gz") ||
                   strstr(inputQueue[i].c_str(), ".mps")){
                  continue;
               }
               if (!inputQueue[i].c_str() || !strncmp(inputQueue[i].c_str(), "-csv", 4)){
                  break;
               }
               fprintf(fp, "%s ", inputQueue[i].c_str());
            }
            fprintf(fp, "\n");
            fclose(fp);
          } break;
          case CbcParam::PRINTSOL:
          case CbcParam::WRITESOL:
          case CbcParam::WRITENEXTSOL:
	    
          case CbcParam::WRITEGMPLSOL:{
            if (!goodModel){
               printGeneralWarning(model_, "** Current model not valid\n");
               continue;
            }
              fp = NULL;
              bool append = true;
              ClpSimplex *saveLpSolver = NULL;
              canOpen = false;
              if (cbcParamCode == CbcParam::PRINTSOL) {
                 fp = stdout;
                 fprintf(fp, "\n");
              } else {
                 cbcParam->readValue(inputQueue, fileName, &message);
                 if (fileName == "append"){
                    switch (cbcParamCode){
                     case CbcParam::WRITESOL:
                       fileName =
                          parameters[CbcParam::SOLUTIONFILE]->fileName();
                       break;
                     case CbcParam::WRITENEXTSOL:
                       fileName =
                          parameters[CbcParam::NEXTSOLFILE]->fileName();
                       break;
                     case CbcParam::WRITEGMPLSOL:
                       fileName =
                          parameters[CbcParam::SOLUTIONFILE]->fileName();
                    break;
                    }
                    CoinParamUtils::processFile(fileName,
                                  parameters[CbcParam::DIRECTORY]->dirName(),
                                                &canOpen);
                    if (!canOpen) {
                       buffer.str("");
                       buffer << "Unable to open file " << fileName.c_str();
                       printGeneralMessage(model_, buffer.str());
                       continue;
                    }
                    append = true;
                 } else if (cbcParamCode == CbcParam::WRITESOL) {
                     CoinParamUtils::processFile(fileName,
                                 parameters[CbcParam::DIRECTORY]->dirName());
                    if (fileName == ""){
                       fileName =
                          parameters[CbcParam::SOLUTIONFILE]->fileName();
                    }else{
                       parameters[CbcParam::SOLUTIONFILE]->setFileName(fileName);
                 }
                 } else if (cbcParamCode == CbcParam::WRITENEXTSOL) {
                     CoinParamUtils::processFile(fileName,
                                 parameters[CbcParam::DIRECTORY]->dirName());
                    if (fileName == ""){
                       fileName = parameters[CbcParam::NEXTSOLFILE]->fileName();
                    }else{
                       parameters[CbcParam::NEXTSOLFILE]->setFileName(fileName);
                    }
                 } else if (cbcParamCode == CbcParam::WRITEGMPLSOL) {
                     CoinParamUtils::processFile(fileName,
                                 parameters[CbcParam::DIRECTORY]->dirName());
                    if (fileName == ""){
                       fileName = parameters[CbcParam::GMPLSOLFILE]->fileName();
                    }else{
                       parameters[CbcParam::GMPLSOLFILE]->setFileName(fileName);
                    }
                 }
                 if (!append){
                    fp = fopen(fileName.c_str(), "w");
                 } else {
                    fp = fopen(fileName.c_str(), "a");
                 }
                 if (!fp) {
                    buffer.str("");
                    buffer << "Unable to open file " << fileName.c_str();
                    printGeneralMessage(model_, buffer.str());
                    continue;
                 }
              }


#ifndef CBC_OTHER_SOLVER
              // See if Glpk
                if (cbcParamCode == CbcParam::WRITEGMPLSOL) {
                  int numberRows = lpSolver->getNumRows();
                  int numberColumns = lpSolver->getNumCols();
                  int numberGlpkRows = numberRows + 1;
#ifdef COINUTILS_HAS_GLPK
                  if (coin_glp_prob) {
                    // from gmpl
                    numberGlpkRows = glp_get_num_rows(coin_glp_prob);
                    if (numberGlpkRows != numberRows)
                      printf("Mismatch - cbc %d rows, glpk %d\n", numberRows,
                             numberGlpkRows);
                  }
#endif
                  fprintf(fp, "%d %d\n", numberGlpkRows, numberColumns);
                  int iStat = lpSolver->status();
                  int iStat2 = GLP_UNDEF;
                  bool integerProblem = false;
                  if (integerStatus >= 0) {
                    iStat = integerStatus;
                    integerProblem = true;
                  }
                  if (iStat == 0) {
                    // optimal
                    if (integerProblem)
                      iStat2 = GLP_OPT;
                    else
                      iStat2 = GLP_FEAS;
                  } else if (iStat == 1) {
                    // infeasible
                    iStat2 = GLP_NOFEAS;
                  } else if (iStat == 2) {
                    // unbounded
                    // leave as 1
                  } else if (iStat >= 3 && iStat <= 5) {
                    if (babModel_ && !babModel_->bestSolution())
                      iStat2 = GLP_NOFEAS;
                    else
                      iStat2 = GLP_FEAS;
                  } else if (iStat == 6) {
                    // bab infeasible
                    iStat2 = GLP_NOFEAS;
                  }
                  lpSolver->computeObjectiveValue(false);
                  double objValue = clpSolver->getObjValue();
                  if (integerProblem)
                    fprintf(fp, "%d %g\n", iStat2, objValue);
                  else
                    fprintf(fp, "%d 2 %g\n", iStat2, objValue);
                  if (numberGlpkRows > numberRows) {
                    // objective as row
                    if (integerProblem) {
                      fprintf(fp, "%g\n", objValue);
                    } else {
                      fprintf(fp, "4 %g 1.0\n", objValue);
                    }
                  }
                  int lookup[6] = {4, 1, 3, 2, 4, 5};
                  const double *primalRowSolution =
                      lpSolver->primalRowSolution();
                  const double *dualRowSolution = lpSolver->dualRowSolution();
                  for (int i = 0; i < numberRows; i++) {
                    if (integerProblem) {
                      fprintf(fp, "%g\n", primalRowSolution[i]);
                    } else {
                      fprintf(fp, "%d %g %g\n",
                              lookup[lpSolver->getRowStatus(i)],
                              primalRowSolution[i], dualRowSolution[i]);
                    }
                  }
                  const double *primalColumnSolution =
                      lpSolver->primalColumnSolution();
                  const double *dualColumnSolution =
                      lpSolver->dualColumnSolution();
                  for (int i = 0; i < numberColumns; i++) {
                    if (integerProblem) {
                      fprintf(fp, "%g\n", primalColumnSolution[i]);
                    } else {
                      fprintf(fp, "%d %g %g\n",
                              lookup[lpSolver->getColumnStatus(i)],
                              primalColumnSolution[i], dualColumnSolution[i]);
                    }
                  }
                  fclose(fp);
#ifdef COINUTILS_HAS_GLPK
                  if (coin_glp_prob) {
                    if (integerProblem) {
                      glp_read_mip(coin_glp_prob, fileName.c_str());
                      glp_mpl_postsolve(coin_glp_tran, coin_glp_prob, GLP_MIP);
                    } else {
                      glp_read_sol(coin_glp_prob, fileName.c_str());
                      glp_mpl_postsolve(coin_glp_tran, coin_glp_prob, GLP_SOL);
                    }
                    // free up as much as possible
                    glp_free(coin_glp_prob);
                    glp_mpl_free_wksp(coin_glp_tran);
                    coin_glp_prob = NULL;
                    coin_glp_tran = NULL;
                    // gmp_free_mem();
                    /* check that no memory blocks are still allocated */
                    glp_free_env();
                  }
#endif
                  break;
                }

                bool printingAllAsCsv = false;
                if (printMode == 14) {
                  // when allcsv then set printMode to all and
                  // change the output format to csv
                  printMode = 4;
                  printingAllAsCsv = true;
                }

                if (printMode < 5) {
                  if (cbcParamCode == CbcParam::WRITENEXTSOL) {
                    // save
                    const double *nextBestSolution =
                        model_.savedSolution(currentBestSolution++);
                    if (!nextBestSolution) {
                      buffer.str("");
                      buffer << "All alternative solutions printed";
                      printGeneralMessage(model_, buffer.str());
                      break;
                    } else {
                      buffer.str("");
                      buffer << "Alternative solution - "  
                             << model_.numberSavedSolutions() - 2 << " remaining";
                      printGeneralMessage(model_, buffer.str());
                    }
                    saveLpSolver = lpSolver;
                    assert(clpSolver->getModelPtr() == saveLpSolver);
                    lpSolver = new ClpSimplex(*saveLpSolver);
#ifndef NDEBUG
                    ClpSimplex *oldSimplex = clpSolver->swapModelPtr(lpSolver);
                    assert(oldSimplex == saveLpSolver);
#else
                    clpSolver->swapModelPtr(lpSolver);
#endif
                    double *solution = lpSolver->primalColumnSolution();
                    double *lower = lpSolver->columnLower();
                    double *upper = lpSolver->columnUpper();
                    int numberColumns = lpSolver->numberColumns();
                    memcpy(solution, nextBestSolution,
                           numberColumns * sizeof(double));
                    model_.deleteSavedSolution(1);
                    for (int i = 0; i < numberColumns; i++) {
                      if (clpSolver->isInteger(i)) {
                        double value = floor(solution[i] + 0.5);
                        lower[i] = value;
                        upper[i] = value;
                      }
                    }
                    lpSolver->allSlackBasis();
                    lpSolver->initialSolve();
                  }
                  // Write solution header (suggested by Luigi Poderico)
                  // Refresh solver
                  lpSolver = clpSolver->getModelPtr();
		  ClpQuadraticObjective * quadObj =
		    dynamic_cast<ClpQuadraticObjective *>(lpSolver->objectiveAsObject());
		  double objValue;
		  if (!quadObj) {
		    lpSolver->computeObjectiveValue(false);
		    objValue = lpSolver->getObjValue(); // why not just use model_.getObjValue(); - maybe if no solution
		  } else {
                    double *solution = lpSolver->primalColumnSolution();
		    objValue = quadObj->objectiveValue(lpSolver,solution);
		    // offset
		    objValue -= lpSolver->objectiveOffset();
		  }
                  int iStat = lpSolver->status();
                  int iStat2 = -1;
                  if (integerStatus >= 0) {
                    iStat = integerStatus;
                    iStat2 = babModel_ ? babModel_->secondaryStatus()
		      : model_.secondaryStatus();
                  }
                  if (iStat == 0) {
                    fprintf(fp, "Optimal");
                    if (iStat2 == 2) {
                      fprintf(fp, " (within gap tolerance)");
                    }
                  } else if (iStat == 1) {
                    // infeasible
                    fprintf(fp, "Infeasible");
                  } else if (iStat == 2) {
                    // unbounded
                    fprintf(fp, "Unbounded");
                  } else if (iStat >= 3 && iStat <= 5) {
                    if (iStat == 3) {
                      if (iStat2 == 4) {
                        fprintf(fp, "Stopped on time");
                      } else {
                        fprintf(fp, "Stopped on iterations");
                      }
                    } else if (iStat == 4) {
                      fprintf(fp, "Stopped on difficulties");
                    } else {
                      fprintf(fp, "Stopped on ctrl-c");
                    }
                    if (babModel_ && !babModel_->bestSolution())
                      fprintf(fp, " (no integer solution - continuous used)");
                  } else if (iStat == 6) {
                    // bab infeasible
                    fprintf(fp, "Integer infeasible");
                  } else {
                    fprintf(fp, "Status unknown");
                  }
                  fprintf(fp, " - objective value %.8f\n", objValue);
                }
#endif
                // make fancy later on
                int iRow;
                int numberRows = clpSolver->getNumRows();
                const double *dualRowSolution = clpSolver->getRowPrice();
                const double *primalRowSolution = clpSolver->getRowActivity();
                const double *rowLower = clpSolver->getRowLower();
                const double *rowUpper = clpSolver->getRowUpper();
                double primalTolerance;
                clpSolver->getDblParam(OsiPrimalTolerance, primalTolerance);
                size_t lengthPrint =
                    static_cast<size_t>(std::max(lengthName, 8));
                bool doMask = (parameters[CbcParam::PRINTMASK]->strVal()
                               != "" && lengthName);
                int *maskStarts = NULL;
                int maxMasks = 0;
                char **masks = NULL;
                if (doMask) {
                  int nAst = 0;
                  const char *pMask2 =
                     parameters[CbcParam::PRINTMASK]->strVal().c_str();
                  char pMask[100];
                  size_t iChar;
                  size_t lengthMask = strlen(pMask2);
                  assert(lengthMask < 100);
                  if (*pMask2 == '"') {
                    if (pMask2[lengthMask - 1] != '"') {
                      printf("mismatched \" in mask %s\n", pMask2);
                      break;
                    } else {
                      strcpy(pMask, pMask2 + 1);
                      *strchr(pMask, '"') = '\0';
                    }
                  } else if (*pMask2 == '\'') {
                    if (pMask2[lengthMask - 1] != '\'') {
                      printf("mismatched ' in mask %s\n", pMask2);
                      break;
                    } else {
                      strcpy(pMask, pMask2 + 1);
                      *strchr(pMask, '\'') = '\0';
                    }
                  } else {
                    strcpy(pMask, pMask2);
                  }
                  if (lengthMask > static_cast<size_t>(lengthName)) {
                    printf("mask %s too long - skipping\n", pMask);
                    break;
                  }
                  maxMasks = 1;
                  for (iChar = 0; iChar < lengthMask; iChar++) {
                    if (pMask[iChar] == '*') {
                      nAst++;
                      maxMasks *= (lengthName + 1);
                    }
                  }
                  int nEntries = 1;
                  maskStarts = new int[lengthName + 2];
                  masks = new char *[maxMasks];
                  char **newMasks = new char *[maxMasks];
                  int i;
                  for (i = 0; i < maxMasks; i++) {
                    masks[i] = new char[lengthName + 1];
                    newMasks[i] = new char[lengthName + 1];
                  }
                  strcpy(masks[0], pMask);
                  for (int iAst = 0; iAst < nAst; iAst++) {
                    int nOldEntries = nEntries;
                    nEntries = 0;
                    for (int iEntry = 0; iEntry < nOldEntries; iEntry++) {
                      char *oldMask = masks[iEntry];
                      char *ast = strchr(oldMask, '*');
                      assert(ast);
                      size_t length = strlen(oldMask) - 1;
                      size_t nBefore = ast - oldMask;
                      size_t nAfter = length - nBefore;
                      // and add null
                      nAfter++;
                      for (int i = 0;
                           i <= lengthName - static_cast<int>(length); i++) {
                        char *maskOut = newMasks[nEntries];
                        memcpy(maskOut, oldMask, nBefore);
                        for (int k = 0; k < i; k++)
                          maskOut[k + nBefore] = '?';
                        memcpy(maskOut + nBefore + i, ast + 1, nAfter);
                        nEntries++;
                        assert(nEntries <= maxMasks);
                      }
                    }
                    char **temp = masks;
                    masks = newMasks;
                    newMasks = temp;
                  }
                  // Now extend and sort
                  int *sort = new int[nEntries];
                  for (i = 0; i < nEntries; i++) {
                    char *maskThis = masks[i];
                    size_t length = strlen(maskThis);
                    while (length > 0 && maskThis[length - 1] == ' ')
                      length--;
                    maskThis[length] = '\0';
                    sort[i] = static_cast<int>(length);
                  }
                  CoinSort_2(sort, sort + nEntries, masks);
                  int lastLength = -1;
                  for (i = 0; i < nEntries; i++) {
                    int length = sort[i];
                    while (length > lastLength)
                      maskStarts[++lastLength] = i;
                  }
                  maskStarts[++lastLength] = nEntries;
                  delete[] sort;
                  for (i = 0; i < maxMasks; i++)
                    delete[] newMasks[i];
                  delete[] newMasks;
                }
                if (printMode > 5 && printMode < 12) {
                  ClpSimplex *solver = clpSolver->getModelPtr();
                  int numberColumns = numberPrintingColumns(clpSolver);
                  // int numberColumns = solver->numberColumns();
                  // column length unless rhs ranging
                  int number = numberColumns;
                  if (lpSolver->status()) {
                    fprintf(fp, "**** Results not valid when LP not optimal\n");
                    number = 0;
                  }
                  switch (printMode) {
                    // bound ranging
                  case 6:
                    fprintf(fp, "Bound ranging");
                    break;
                    // rhs ranging
                  case 7:
                    fprintf(fp, "Rhs ranging");
                    if (!lpSolver->status())
                      number = numberRows;
                    break;
                    // objective ranging
                  case 8:
                    fprintf(fp, "Objective ranging");
                    break;
                  }
                  if (lengthName)
                    fprintf(fp, ",name");
                  fprintf(fp, ",increase,variable,decrease,variable\n");
                  int *which = new int[number];
                  if (printMode != 7) {
                    if (!doMask) {
                      for (int i = 0; i < number; i++)
                        which[i] = i;
                    } else {
                      int n = 0;
                      for (int i = 0; i < number; i++) {
                        if (maskMatches(maskStarts, masks, columnNames[i]))
                          which[n++] = i;
                      }
                      if (n) {
                        number = n;
                      } else {
                        printf("No names match - doing all\n");
                        for (int i = 0; i < number; i++)
                          which[i] = i;
                      }
                    }
                  } else {
                    if (!doMask) {
                      for (int i = 0; i < number; i++)
                        which[i] = i + numberColumns;
                    } else {
                      int n = 0;
                      for (int i = 0; i < number; i++) {
                        if (maskMatches(maskStarts, masks, rowNames[i]))
                          which[n++] = i + numberColumns;
                      }
                      if (n) {
                        number = n;
                      } else {
                        printf("No names match - doing all\n");
                        for (int i = 0; i < number; i++)
                          which[i] = i + numberColumns;
                      }
                    }
                  }
                  double *valueIncrease = new double[number];
                  int *sequenceIncrease = new int[number];
                  double *valueDecrease = new double[number];
                  int *sequenceDecrease = new int[number];
                  switch (printMode) {
                    // bound or rhs ranging
                  case 6:
                  case 7:
                    solver->primalRanging(number, which, valueIncrease,
                                          sequenceIncrease, valueDecrease,
                                          sequenceDecrease);
                    break;
                    // objective ranging
                  case 8:
                    solver->dualRanging(number, which, valueIncrease,
                                        sequenceIncrease, valueDecrease,
                                        sequenceDecrease);
                    break;
                  }
                  for (int i = 0; i < number; i++) {
                    int iWhich = which[i];
                    fprintf(fp, "%d,",
                            (iWhich < numberColumns) ? iWhich
                                                     : iWhich - numberColumns);
                    if (lengthName) {
                      const char *name =
                          (printMode == 7)
                              ? rowNames[iWhich - numberColumns].c_str()
                              : columnNames[iWhich].c_str();
                      fprintf(fp, "%s,", name);
                    }
                    if (valueIncrease[i] < 1.0e30) {
                      fprintf(fp, "%.10g,", valueIncrease[i]);
                      int outSequence = sequenceIncrease[i];
                      if (outSequence < numberColumns) {
                        if (lengthName)
                          fprintf(fp, "%s,", columnNames[outSequence].c_str());
                        else
                          fprintf(fp, "C%7.7d,", outSequence);
                      } else {
                        outSequence -= numberColumns;
                        if (lengthName)
                          fprintf(fp, "%s,", rowNames[outSequence].c_str());
                        else
                          fprintf(fp, "R%7.7d,", outSequence);
                      }
                    } else {
                      fprintf(fp, "1.0e100,,");
                    }
                    if (valueDecrease[i] < 1.0e30) {
                      fprintf(fp, "%.10g,", valueDecrease[i]);
                      int outSequence = sequenceDecrease[i];
                      if (outSequence < numberColumns) {
                        if (lengthName)
                          fprintf(fp, "%s", columnNames[outSequence].c_str());
                        else
                          fprintf(fp, "C%7.7d", outSequence);
                      } else {
                        outSequence -= numberColumns;
                        if (lengthName)
                          fprintf(fp, "%s", rowNames[outSequence].c_str());
                        else
                          fprintf(fp, "R%7.7d", outSequence);
                      }
                    } else {
                      fprintf(fp, "1.0e100,");
                    }
                    fprintf(fp, "\n");
                  }
                  if (fp != stdout)
                    fclose(fp);
                  delete[] which;
                  delete[] valueIncrease;
                  delete[] sequenceIncrease;
                  delete[] valueDecrease;
                  delete[] sequenceDecrease;
                  if (masks) {
                    delete[] maskStarts;
                    for (int i = 0; i < maxMasks; i++)
                      delete[] masks[i];
                    delete[] masks;
                  }
                  break;
                }

                char printFormat[50];
                sprintf(printFormat, " %s         %s\n",
                        CLP_QUOTE(CLP_OUTPUT_FORMAT),
                        CLP_QUOTE(CLP_OUTPUT_FORMAT));
                char printIntFormat[50];
                sprintf(printIntFormat, " %s         %s\n",
                        CLP_QUOTE(CLP_INTEGER_OUTPUT_FORMAT),
                        CLP_QUOTE(CLP_OUTPUT_FORMAT));
                if (printMode > 2 && printMode < 5) {
                  for (iRow = 0; iRow < numberRows; iRow++) {
                    int type = printMode - 3;

                    std::string valueHasTolerance;
                    std::string rowName;

                    if (primalRowSolution[iRow] >
                            rowUpper[iRow] + primalTolerance ||
                        primalRowSolution[iRow] <
                            rowLower[iRow] - primalTolerance) {
                      if (printingAllAsCsv) {
                        valueHasTolerance = "**";
                      } else {
                        fprintf(fp, "** ");
                      }
                      type = 2;
                    } else if (fabs(primalRowSolution[iRow]) > 1.0e-8) {
                      type = 1;
                    } else if (numberRows < 50) {
                      type = 3;
                    }
                    if (doMask &&
                        !maskMatches(maskStarts, masks, rowNames[iRow]))
                      type = 0;
                    if (type) {
                      if (!printingAllAsCsv) {
                        fprintf(fp, "%7d ", iRow);
                      }

                      if (lengthName) {
                        if (printingAllAsCsv) {
                          rowName = rowNames[iRow];
                        }
                        else{
                          const char *name = rowNames[iRow].c_str();
                          size_t n = strlen(name);
                          size_t i;
                        
                          for (i = 0; i < n; i++)
                              fprintf(fp, "%c", name[i]);
                          for (; i < lengthPrint; i++)
                              fprintf(fp, " ");
                        }
                      }

                      if (!printingAllAsCsv) {
                        fprintf(fp, printFormat, primalRowSolution[iRow],
                          dualRowSolution[iRow]);
                      }
                      else{
                        fprintf(fp,
                                "%s,%d,%s,%.15g,%.15g\n",
                                valueHasTolerance.c_str(),
                                iRow,
                                rowName.c_str(),
                                primalRowSolution[iRow],
                                dualRowSolution[iRow]);
                        }
                    }
                  }
                }
                int iColumn;
                int numberColumns = numberPrintingColumns(clpSolver);
                const double *dualColumnSolution = clpSolver->getReducedCost();
                const double *primalColumnSolution =
                    clpSolver->getColSolution();
                const double *columnLower = clpSolver->getColLower();
                const double *columnUpper = clpSolver->getColUpper();
                if (printMode != 2 && printMode < 12) {
                  if (printMode == 5) {
                    if (lengthName)
                      fprintf(fp, "name");
                    else
                      fprintf(fp, "number");
                    fprintf(fp, ",solution\n");
                  }
                  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                    int type = (printMode > 3) ? 1 : 0;

                    std::string valueHasTolerance;
                    std::string colName; 

                    if (primalColumnSolution[iColumn] >
                            columnUpper[iColumn] + primalTolerance ||
                        primalColumnSolution[iColumn] <
                            columnLower[iColumn] - primalTolerance) {
                      if (printingAllAsCsv) {
                        valueHasTolerance = "**";
                      } else {
                        fprintf(fp, "** ");
                      }
                      type = 2;
                    } else if (fabs(primalColumnSolution[iColumn]) > 1.0e-8) {
                      type = 1;
                    } else if (numberColumns < 50) {
                      type = 3;
                    }
                    // see if integer
                    if ((!clpSolver->isInteger(iColumn) ||
                         fabs(primalColumnSolution[iColumn]) < 1.0e-8) &&
                        printMode == 1)
                      type = 0;
                    if (doMask &&
                        !maskMatches(maskStarts, masks, columnNames[iColumn]))
                      type = 0;
                    if (type) {
                      if (printMode != 5) {
                        if (!printingAllAsCsv) {
                          fprintf(fp, "%7d ", iColumn);
                        }

                        if (lengthName) {
                          if (printingAllAsCsv) {
                              colName = columnNames[iColumn];
                          } else {
                              const char *name = columnNames[iColumn].c_str();
                              size_t n = strlen(name);
                              size_t i;
                              for (i = 0; i < n; i++)
                                 fprintf(fp, "%c", name[i]);
                              for (; i < lengthPrint; i++)
                                 fprintf(fp, " ");
                          }
                        }
                        if (!printingAllAsCsv) {
			  double value = primalColumnSolution[iColumn];
			  double nearest = floor(value+0.5);
			  if (!clpSolver->isInteger(iColumn)||fabs(value-nearest)>1.0e-8) {
			    fprintf(fp, printFormat, value,
				    dualColumnSolution[iColumn]);
			  } else {
			    // allow for very very large integer values
			    long int iValue = nearest; 
			    fprintf(fp, printIntFormat, iValue,
				    dualColumnSolution[iColumn]);
			  }
                        } else {
                          fprintf(fp,
                            "%s,%d,%s,%.15g,%.15g\n",
                            valueHasTolerance.c_str(),
                            iColumn,
                            colName.c_str(),
                            primalColumnSolution[iColumn],
                            dualColumnSolution[iColumn]);
                        }
                      } else {
                        char temp[100];
                        if (lengthName) {
                          const char *name = columnNames[iColumn].c_str();
                          for (int i = 0; i < lengthName; i++)
                            temp[i] = name[i];
                          temp[lengthName] = '\0';
                        } else {
                          sprintf(temp, "%7d", iColumn);
                        }
                        sprintf(temp + strlen(temp), ", %15.8g",
                                primalColumnSolution[iColumn]);
                        size_t n = strlen(temp);
                        size_t k = 0;
                        for (size_t i = 0; i < n + 1; i++) {
                          if (temp[i] != ' ')
                            temp[k++] = temp[i];
                        }
                        fprintf(fp, "%s\n", temp);
                      }
                    }
                  }

                  if (cbcParamCode == CbcParam::WRITENEXTSOL) {
                    if (saveLpSolver) {
                      clpSolver->swapModelPtr(saveLpSolver);
                      delete lpSolver;
                      lpSolver = saveLpSolver;
                      saveLpSolver = NULL;
                    }
                  }
                } else if (printMode == 2) {
                  // special format suitable for OsiRowCutDebugger
                  int n = 0;
                  bool comma = false;
                  bool newLine = false;
                  fprintf(fp, "\tint intIndicesV[]={\n");
                  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                    if (primalColumnSolution[iColumn] > 0.5 &&
                        model_.solver()->isInteger(iColumn)) {
                      if (comma)
                        fprintf(fp, ",");
                      if (newLine)
                        fprintf(fp, "\n");
                      fprintf(fp, "%d ", iColumn);
                      comma = true;
                      newLine = false;
                      n++;
                      if (n == 10) {
                        n = 0;
                        newLine = true;
                      }
                    }
                  }
                  fprintf(fp, "};\n");
                  n = 0;
                  comma = false;
                  newLine = false;
                  fprintf(fp, "\tdouble intSolnV[]={\n");
                  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                    if (primalColumnSolution[iColumn] > 0.5 &&
                        model_.solver()->isInteger(iColumn)) {
                      if (comma)
                        fprintf(fp, ",");
                      if (newLine)
                        fprintf(fp, "\n");
                      int value =
                          static_cast<int>(primalColumnSolution[iColumn] + 0.5);
                      fprintf(fp, "%d. ", value);
                      comma = true;
                      newLine = false;
                      n++;
                      if (n == 10) {
                        n = 0;
                        newLine = true;
                      }
                    }
                  }
                  fprintf(fp, "};\n");
                } else {
                  // Make up a fake bounds section
                  char outputValue[24];
                  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                    if (printMode == 13 ||
                        model_.solver()->isInteger(iColumn)) {
                      fprintf(fp, " FX BOUND001  ");
                      const char *name = columnNames[iColumn].c_str();
                      size_t n = strlen(name);
                      size_t i;
                      for (i = 0; i < n; i++)
                        fprintf(fp, "%c", name[i]);
                      for (; i < lengthPrint; i++)
                        fprintf(fp, " ");
                      CoinConvertDouble(5, 2, primalColumnSolution[iColumn],
                                        outputValue);
                      fprintf(fp, "  %s\n", outputValue);
                    } else {
                      fprintf(fp, " LO BOUND001  ");
                      const char *name = columnNames[iColumn].c_str();
                      size_t n = strlen(name);
                      size_t i;
                      for (i = 0; i < n; i++)
                        fprintf(fp, "%c", name[i]);
                      for (; i < lengthPrint; i++)
                        fprintf(fp, " ");
                      CoinConvertDouble(5, 2,
                                        std::max(-1.0e30, columnLower[iColumn]),
                                        outputValue);
                      fprintf(fp, "  %s\n", outputValue);
                      fprintf(fp, " UP BOUND001  ");
                      for (i = 0; i < n; i++)
                        fprintf(fp, "%c", name[i]);
                      for (; i < lengthPrint; i++)
                        fprintf(fp, " ");
                      CoinConvertDouble(5, 2,
                                        std::min(1.0e30, columnUpper[iColumn]),
                                        outputValue);
                      fprintf(fp, "  %s\n", outputValue);
                    }
                  }
                  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                    if (primalColumnSolution[iColumn] >
                            columnUpper[iColumn] + primalTolerance ||
                        primalColumnSolution[iColumn] <
                            columnLower[iColumn] - primalTolerance) {
                      fprintf(fp, " FX BOUND002  ");
                      const char *name = columnNames[iColumn].c_str();
                      size_t n = strlen(name);
                      size_t i;
                      for (i = 0; i < n; i++)
                        fprintf(fp, "%c", name[i]);
                      for (; i < lengthPrint; i++)
                        fprintf(fp, " ");
                      CoinConvertDouble(5, 2, primalColumnSolution[iColumn],
                                        outputValue);
                      fprintf(fp, "  %s\n", outputValue);
                    }
                  }
                }
                if (fp != stdout)
                  fclose(fp);
                if (masks) {
                  delete[] maskStarts;
                  for (int i = 0; i < maxMasks; i++)
                    delete[] masks[i];
                  delete[] masks;
                }
          } break;
          case CbcParam::WRITESOLBINARY:
            if (!goodModel){
               printGeneralWarning(model_, "** Current model not valid\n");
               continue;
            }
            cbcParam->readValue(inputQueue, fileName, &message);
            CoinParamUtils::processFile(fileName,
                                 parameters[CbcParam::DIRECTORY]->dirName());
            if (fileName == ""){
               fileName = parameters[CbcParam::SOLUTIONBINARYFILE]->fileName();
            }else{
               parameters[CbcParam::SOLUTIONBINARYFILE]->setFileName(fileName);
            }
            ClpParamUtils::saveSolution(lpSolver, fileName);
            break;
          case CbcParam::ENVIRONMENT: {
#if !defined(_MSC_VER) && !defined(__MSVCRT__)
            // Don't think it will work with Visual Studio
            char *input = getenv("CBC_ENVIRONMENT");
            if (input){
               while (!inputQueue.empty()){
                  inputQueue.pop_front();
               }
               std::istringstream inputStream(input);
               CoinParamUtils::readFromStream(inputQueue, inputStream);
            }
#else
            std::cout << "** Parameter not valid on Windows with Visual Studio"
                      << std::endl;
#endif
            }
            break;
          case CbcParam::DUMMY:
            break;
          default:
            // abort();
            break;
          }
        }
#if 0
      }
#endif
    }
    delete coinModel;
  }
#if CBC_QUIET == 0
  buffer.str("");
  buffer << "Total time (CPU seconds):       " << CoinCpuTime() - time0
         << "   (Wallclock seconds):       "
         << CoinGetTimeOfDay() - time0Elapsed << std::endl;
  printGeneralMessage(model_, buffer.str());
#endif
#ifdef COINUTILS_HAS_GLPK
  if (coin_glp_prob) {
    // free up as much as possible
    glp_free(coin_glp_prob);
    glp_mpl_free_wksp(coin_glp_tran);
    glp_free_env();
    coin_glp_prob = NULL;
    coin_glp_tran = NULL;
  }
#endif
  delete[] lotsize;
  if (statistics_number_cuts != NULL)
    delete[] statistics_number_cuts;

  if (statistics_name_generators != NULL)
    delete[] statistics_name_generators;
    // By now all memory should be freed
#ifdef DMALLOC
    // dmalloc_log_unfreed();
    // dmalloc_shutdown();
#endif
  if (babModel_) {
    model_.moveInfo(*babModel_);
#ifndef CBC_OTHER_SOLVER
    OsiClpSolverInterface *clpSolver0 =
        dynamic_cast<OsiClpSolverInterface *>(babModel_->solver());
    ClpSimplex *lpSolver0 = clpSolver0->getModelPtr();
    OsiClpSolverInterface *clpSolver =
        dynamic_cast<OsiClpSolverInterface *>(model_.solver());
    ClpSimplex *lpSolver = clpSolver->getModelPtr();
    if (lpSolver0 != lpSolver && lpSolver != originalSolver->getModelPtr())
      lpSolver->moveInfo(*lpSolver0);
      // babModel_->setModelOwnsSolver(false);
#endif
  }

  delete babModel_;

  babModel_ = NULL;
  model_.solver()->setWarmStart(NULL);
  buffer.str("");
  buffer << "Total time " << CoinCpuTime() - time0;
  // generalMessageHandler->message(CLP_GENERAL, generalMessages)
  //<< generalPrint
  //<< CoinMessageEol;
  return 0;
}

//###########################################################################
//###########################################################################

int cbcReadAmpl(ampl_info *info, int argc, char **argv, CbcModel &model)
{
   OsiSolverInterface *solver = model.solver();
   memset(info, 0, sizeof(info));
   if (argc > 2 && !strcmp(argv[2], "-AMPL")) {
      // see if log in list
      bool noPrinting = true;
      for (int i = 1; i < argc; i++) {
         if (!strncmp(argv[i], "log", 3)) {
            const char *equals = strchr(argv[i], '=');
            if (equals && atoi(equals + 1) > 0) {
               noPrinting = false;
               info->logLevel = atoi(equals + 1);
               // mark so won't be overWritten
               info->numberRows = -1234567;
               break;
            }
         }
      }
      
      void *voidModel;
      int returnCode = readAmpl(info, argc, const_cast<char **>(argv),
                                &voidModel, "cbc");
      if (returnCode){
         return returnCode;
      }
      // see if log in list (including environment)
      for (int i = 1; i < info->numberArguments; i++) {
         if (!strcmp(info->arguments[i], "log")) {
            if (i < info->numberArguments - 1 && atoi(info->arguments[i + 1]) > 0)
               noPrinting = false;
            break;
         }
      }
      if (noPrinting) {
         model.messageHandler()->setLogLevel(0);
      }
      if (!noPrinting){
         printf("%d rows, %d columns and %d elements\n", info->numberRows,
                info->numberColumns, info->numberElements);
      }
   }
   return 0;
}

//###########################################################################
//  Routines to print statistics.
//###########################################################################

static void breakdown(const char *name, int numberLook, const double *region) {
  double range[] = {
      -COIN_DBL_MAX, -1.0e15,     -1.0e11,  -1.0e8,   -1.0e5,  -1.0e4,  -1.0e3,
      -1.0e2,        -1.0e1,      -1.0,     -1.0e-1,  -1.0e-2, -1.0e-3, -1.0e-4,
      -1.0e-5,       -1.0e-8,     -1.0e-11, -1.0e-15, 0.0,     1.0e-15, 1.0e-11,
      1.0e-8,        1.0e-5,      1.0e-4,   1.0e-3,   1.0e-2,  1.0e-1,  1.0,
      1.0e1,         1.0e2,       1.0e3,    1.0e4,    1.0e5,   1.0e8,   1.0e11,
      1.0e15,        COIN_DBL_MAX};
  int nRanges = static_cast<int>(sizeof(range) / sizeof(double));
  int *number = new int[nRanges];
  memset(number, 0, nRanges * sizeof(int));
  int *numberExact = new int[nRanges];
  memset(numberExact, 0, nRanges * sizeof(int));
  int i;
  for (i = 0; i < numberLook; i++) {
    double value = region[i];
    for (int j = 0; j < nRanges; j++) {
      if (value == range[j]) {
        numberExact[j]++;
        break;
      } else if (value < range[j]) {
        number[j]++;
        break;
      }
    }
  }
  printf("\n%s has %d entries\n", name, numberLook);
  for (i = 0; i < nRanges; i++) {
    if (number[i])
      printf("%d between %g and %g", number[i], range[i - 1], range[i]);
    if (numberExact[i]) {
      if (number[i])
        printf(", ");
      printf("%d exactly at %g", numberExact[i], range[i]);
    }
    if (number[i] + numberExact[i])
      printf("\n");
  }
  delete[] number;
  delete[] numberExact;
}
static void sortOnOther(int *column, const CoinBigIndex *rowStart, int *order,
                        int *other, int nRow, int nInRow, int where) {
  if (nRow < 2 || where >= nInRow)
    return;
  // do initial sort
  int kRow;
  int iRow;
  for (kRow = 0; kRow < nRow; kRow++) {
    iRow = order[kRow];
    other[kRow] = column[rowStart[iRow] + where];
  }
  CoinSort_2(other, other + nRow, order);
  int first = 0;
  iRow = order[0];
  int firstC = column[rowStart[iRow] + where];
  kRow = 1;
  while (kRow < nRow) {
    int lastC = 9999999;
    ;
    for (; kRow < nRow + 1; kRow++) {
      if (kRow < nRow) {
        iRow = order[kRow];
        lastC = column[rowStart[iRow] + where];
      } else {
        lastC = 9999999;
      }
      if (lastC > firstC)
        break;
    }
    // sort
    sortOnOther(column, rowStart, order + first, other, kRow - first, nInRow,
                where + 1);
    firstC = lastC;
    first = kRow;
  }
}

//###########################################################################
//###########################################################################

static void statistics(ClpSimplex *originalModel, ClpSimplex *model) {
  int numberColumns = originalModel->numberColumns();
  const char *integerInformation = originalModel->integerInformation();
  const double *columnLower = originalModel->columnLower();
  const double *columnUpper = originalModel->columnUpper();
  int numberIntegers = 0;
  int numberBinary = 0;
  int iRow, iColumn;
  if (integerInformation) {
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (integerInformation[iColumn]) {
        if (columnUpper[iColumn] > columnLower[iColumn]) {
          numberIntegers++;
          if (columnLower[iColumn] == 0.0 && columnUpper[iColumn] == 1)
            numberBinary++;
        }
      }
    }
    printf("Original problem has %d integers (%d of which binary)\n",
           numberIntegers, numberBinary);
  }
  numberColumns = model->numberColumns();
  int numberRows = model->numberRows();
  columnLower = model->columnLower();
  columnUpper = model->columnUpper();
  const double *rowLower = model->rowLower();
  const double *rowUpper = model->rowUpper();
  const double *objective = model->objective();
  if (model->integerInformation()) {
    const char *integerInformation = model->integerInformation();
    int numberIntegers = 0;
    int numberBinary = 0;
    double *obj = new double[numberColumns];
    int *which = new int[numberColumns];
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      if (columnUpper[iColumn] > columnLower[iColumn]) {
        if (integerInformation[iColumn]) {
          numberIntegers++;
          if (columnLower[iColumn] == 0.0 && columnUpper[iColumn] == 1)
            numberBinary++;
        }
      }
    }
    if (numberColumns != originalModel->numberColumns())
      printf("Presolved problem has %d integers (%d of which binary)\n",
             numberIntegers, numberBinary);
    for (int ifInt = 0; ifInt < 2; ifInt++) {
      for (int ifAbs = 0; ifAbs < 2; ifAbs++) {
        int numberSort = 0;
        int numberZero = 0;
        int numberDifferentObj = 0;
        for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
          if (columnUpper[iColumn] > columnLower[iColumn]) {
            if (!ifInt || integerInformation[iColumn]) {
              obj[numberSort] =
                  (ifAbs) ? fabs(objective[iColumn]) : objective[iColumn];
              which[numberSort++] = iColumn;
              if (!objective[iColumn])
                numberZero++;
            }
          }
        }
        CoinSort_2(obj, obj + numberSort, which);
        double last = obj[0];
        for (int jColumn = 1; jColumn < numberSort; jColumn++) {
          if (fabs(obj[jColumn] - last) > 1.0e-12) {
            numberDifferentObj++;
            last = obj[jColumn];
          }
        }
        numberDifferentObj++;
        printf("==== ");
        if (ifInt)
          printf("for integers ");
        if (!ifAbs)
          printf("%d zero objective ", numberZero);
        else
          printf("absolute objective values ");
        printf("%d different\n", numberDifferentObj);
        bool saveModel = false;
        int target = model->logLevel();
        if (target > 10000) {
          if (ifInt && !ifAbs)
            saveModel = true;
          target -= 10000;
        }

        if (target <= 100)
          target = 12;
        else
          target -= 100;
        if (numberDifferentObj < target) {
          int iLast = 0;
          double last = obj[0];
          for (int jColumn = 1; jColumn < numberSort; jColumn++) {
            if (fabs(obj[jColumn] - last) > 1.0e-12) {
              printf("%d variables have objective of %g\n", jColumn - iLast,
                     last);
              iLast = jColumn;
              last = obj[jColumn];
            }
          }
          printf("%d variables have objective of %g\n", numberSort - iLast,
                 last);
          if (saveModel) {
            int spaceNeeded = numberSort + numberDifferentObj;
            CoinBigIndex *columnAddDummy =
                new CoinBigIndex[numberDifferentObj + 1];
            int *columnAdd = new int[spaceNeeded];
            double *elementAdd = new double[spaceNeeded];
            CoinBigIndex *rowAdd = new CoinBigIndex[2 * numberDifferentObj + 1];
            int *newIsInteger =
                reinterpret_cast<int *>(rowAdd + numberDifferentObj + 1);
            double *objectiveNew = new double[3 * numberDifferentObj];
            double *lowerNew = objectiveNew + numberDifferentObj;
            double *upperNew = lowerNew + numberDifferentObj;
            memset(columnAddDummy, 0,
                   (numberDifferentObj + 1) * sizeof(CoinBigIndex));
            ClpSimplex tempModel = *model;
            int iLast = 0;
            double last = obj[0];
            numberDifferentObj = 0;
            int numberElements = 0;
            rowAdd[0] = 0;
            double *objective = tempModel.objective();
            for (int jColumn = 1; jColumn < numberSort + 1; jColumn++) {
              if (jColumn == numberSort ||
                  fabs(obj[jColumn] - last) > 1.0e-12) {
                // not if just one
                if (jColumn - iLast > 1) {
                  bool allInteger = integerInformation != NULL;
                  int iColumn = which[iLast];
                  objectiveNew[numberDifferentObj] = objective[iColumn];
                  double lower = 0.0;
                  double upper = 0.0;
                  for (int kColumn = iLast; kColumn < jColumn; kColumn++) {
                    iColumn = which[kColumn];
                    objective[iColumn] = 0.0;
                    double lowerValue = columnLower[iColumn];
                    double upperValue = columnUpper[iColumn];
                    double elementValue = -1.0;
                    if (objectiveNew[numberDifferentObj] * objective[iColumn] <
                        0.0) {
                      lowerValue = -columnUpper[iColumn];
                      upperValue = -columnLower[iColumn];
                      elementValue = 1.0;
                    }
                    columnAdd[numberElements] = iColumn;
                    elementAdd[numberElements++] = elementValue;
                    if (integerInformation && !integerInformation[iColumn])
                      allInteger = false;
                    if (lower != -COIN_DBL_MAX) {
                      if (lowerValue != -COIN_DBL_MAX)
                        lower += lowerValue;
                      else
                        lower = -COIN_DBL_MAX;
                    }
                    if (upper != COIN_DBL_MAX) {
                      if (upperValue != COIN_DBL_MAX)
                        upper += upperValue;
                      else
                        upper = COIN_DBL_MAX;
                    }
                  }
                  columnAdd[numberElements] =
                      numberColumns + numberDifferentObj;
                  elementAdd[numberElements++] = 1.0;
                  newIsInteger[numberDifferentObj] = (allInteger) ? 1 : 0;
                  lowerNew[numberDifferentObj] = lower;
                  upperNew[numberDifferentObj] = upper;
                  numberDifferentObj++;
                  rowAdd[numberDifferentObj] = numberElements;
                }
                iLast = jColumn;
                last = obj[jColumn];
              }
            }
            // add columns
            tempModel.addColumns(numberDifferentObj, lowerNew, upperNew,
                                 objectiveNew, columnAddDummy, NULL, NULL);
            // add constraints and make integer if all integer in group
            for (int iObj = 0; iObj < numberDifferentObj; iObj++) {
              lowerNew[iObj] = 0.0;
              upperNew[iObj] = 0.0;
              if (newIsInteger[iObj])
                tempModel.setInteger(numberColumns + iObj);
            }
            tempModel.addRows(numberDifferentObj, lowerNew, upperNew, rowAdd,
                              columnAdd, elementAdd);
            delete[] columnAdd;
            delete[] columnAddDummy;
            delete[] elementAdd;
            delete[] rowAdd;
            delete[] objectiveNew;
            // save
            std::string tempName = model->problemName();
            if (ifInt)
              tempName += "_int";
            if (ifAbs)
              tempName += "_abs";
            tempName += ".mps";
            tempModel.writeMps(tempName.c_str());
          }
        }
      }
    }
    delete[] which;
    delete[] obj;
    printf("===== end objective counts\n");
  }
  CoinPackedMatrix *matrix = model->matrix();
  CoinBigIndex numberElements = matrix->getNumElements();
  const int *columnLength = matrix->getVectorLengths();
  // const CoinBigIndex * columnStart = matrix->getVectorStarts();
  const double *elementByColumn = matrix->getElements();
  int *number = new int[numberRows + 1];
  memset(number, 0, (numberRows + 1) * sizeof(int));
  int numberObjSingletons = 0;
  /* cType
        0 0/inf, 1 0/up, 2 lo/inf, 3 lo/up, 4 free, 5 fix, 6 -inf/0, 7 -inf/up,
        8 0/1
     */
  int cType[9];
  std::string cName[] = {"0.0->inf,",  "0.0->up,",  "lo->inf,",
                         "lo->up,",    "free,",     "fixed,",
                         "-inf->0.0,", "-inf->up,", "0.0->1.0"};
  int nObjective = 0;
  memset(cType, 0, sizeof(cType));
  for (iColumn = 0; iColumn < numberColumns; iColumn++) {
    int length = columnLength[iColumn];
    if (length == 1 && objective[iColumn])
      numberObjSingletons++;
    number[length]++;
    if (objective[iColumn])
      nObjective++;
    if (columnLower[iColumn] > -1.0e20) {
      if (columnLower[iColumn] == 0.0) {
        if (columnUpper[iColumn] > 1.0e20)
          cType[0]++;
        else if (columnUpper[iColumn] == 1.0)
          cType[8]++;
        else if (columnUpper[iColumn] == 0.0)
          cType[5]++;
        else
          cType[1]++;
      } else {
        if (columnUpper[iColumn] > 1.0e20)
          cType[2]++;
        else if (columnUpper[iColumn] == columnLower[iColumn])
          cType[5]++;
        else
          cType[3]++;
      }
    } else {
      if (columnUpper[iColumn] > 1.0e20)
        cType[4]++;
      else if (columnUpper[iColumn] == 0.0)
        cType[6]++;
      else
        cType[7]++;
    }
  }
  /* rType
        0 E 0, 1 E 1, 2 E -1, 3 E other, 4 G 0, 5 G 1, 6 G other,
        7 L 0,  8 L 1, 9 L other, 10 Range 0/1, 11 Range other, 12 free
     */
  int rType[13];
  std::string rName[] = {
      "E 0.0,",          "E 1.0,",       "E -1.0,", "E other,", "G 0.0,",
      "G 1.0,",          "G other,",     "L 0.0,",  "L 1.0,",   "L other,",
      "Range 0.0->1.0,", "Range other,", "Free"};
  memset(rType, 0, sizeof(rType));
  for (iRow = 0; iRow < numberRows; iRow++) {
    if (rowLower[iRow] > -1.0e20) {
      if (rowLower[iRow] == 0.0) {
        if (rowUpper[iRow] > 1.0e20)
          rType[4]++;
        else if (rowUpper[iRow] == 1.0)
          rType[10]++;
        else if (rowUpper[iRow] == 0.0)
          rType[0]++;
        else
          rType[11]++;
      } else if (rowLower[iRow] == 1.0) {
        if (rowUpper[iRow] > 1.0e20)
          rType[5]++;
        else if (rowUpper[iRow] == rowLower[iRow])
          rType[1]++;
        else
          rType[11]++;
      } else if (rowLower[iRow] == -1.0) {
        if (rowUpper[iRow] > 1.0e20)
          rType[6]++;
        else if (rowUpper[iRow] == rowLower[iRow])
          rType[2]++;
        else
          rType[11]++;
      } else {
        if (rowUpper[iRow] > 1.0e20)
          rType[6]++;
        else if (rowUpper[iRow] == rowLower[iRow])
          rType[3]++;
        else
          rType[11]++;
      }
    } else {
      if (rowUpper[iRow] > 1.0e20)
        rType[12]++;
      else if (rowUpper[iRow] == 0.0)
        rType[7]++;
      else if (rowUpper[iRow] == 1.0)
        rType[8]++;
      else
        rType[9]++;
    }
  }
  // Basic statistics
  printf("\n\nProblem has %d rows, %d columns (%d with objective) and %d "
         "elements\n",
         numberRows, numberColumns, nObjective, numberElements);
  if (number[0] + number[1]) {
    printf("There are ");
    if (numberObjSingletons)
      printf("%d singletons with objective ", numberObjSingletons);
    int numberNoObj = number[1] - numberObjSingletons;
    if (numberNoObj)
      printf("%d singletons with no objective ", numberNoObj);
    if (number[0])
      printf("** %d columns have no entries", number[0]);
    printf("\n");
  }
  printf("Column breakdown:\n");
  int k;
  for (k = 0; k < static_cast<int>(sizeof(cType) / sizeof(int)); k++) {
    printf("%d of type %s ", cType[k], cName[k].c_str());
    if (((k + 1) % 3) == 0)
      printf("\n");
  }
  if ((k % 3) != 0)
    printf("\n");
  printf("Row breakdown:\n");
  for (k = 0; k < static_cast<int>(sizeof(rType) / sizeof(int)); k++) {
    printf("%d of type %s ", rType[k], rName[k].c_str());
    if (((k + 1) % 3) == 0)
      printf("\n");
  }
  if ((k % 3) != 0)
    printf("\n");
    //#define SYM
#ifndef SYM
  if (model->logLevel() < 2)
    return;
#endif
  int kMax = model->logLevel() > 3 ? 1000000 : 10;
  k = 0;
  for (iRow = 1; iRow <= numberRows; iRow++) {
    if (number[iRow]) {
      k++;
      printf("%d columns have %d entries\n", number[iRow], iRow);
      if (k == kMax)
        break;
    }
  }
  if (k < numberRows) {
    int kk = k;
    k = 0;
    for (iRow = numberRows; iRow >= 1; iRow--) {
      if (number[iRow]) {
        k++;
        if (k == kMax)
          break;
      }
    }
    if (k > kk) {
      printf("\n    .........\n\n");
      iRow = k;
      k = 0;
      for (; iRow < numberRows; iRow++) {
        if (number[iRow]) {
          k++;
          printf("%d columns have %d entries\n", number[iRow], iRow);
          if (k == kMax)
            break;
        }
      }
    }
  }
  delete[] number;
  printf("\n\n");
  if (model->logLevel() == 63
#ifdef SYM
      || true
#endif
  ) {
    // get column copy
    CoinPackedMatrix columnCopy = *matrix;
    const int *columnLength = columnCopy.getVectorLengths();
    number = new int[numberRows + 1];
    memset(number, 0, (numberRows + 1) * sizeof(int));
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      int length = columnLength[iColumn];
      number[length]++;
    }
    k = 0;
    for (iRow = 1; iRow <= numberRows; iRow++) {
      if (number[iRow]) {
        k++;
      }
    }
    int *row = columnCopy.getMutableIndices();
    const CoinBigIndex *columnStart = columnCopy.getVectorStarts();
    double *element = columnCopy.getMutableElements();
    int *order = new int[numberColumns];
    int *other = new int[numberColumns];
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      int length = columnLength[iColumn];
      order[iColumn] = iColumn;
      other[iColumn] = length;
      CoinBigIndex start = columnStart[iColumn];
      CoinSort_2(row + start, row + start + length, element + start);
    }
    CoinSort_2(other, other + numberColumns, order);
    int jColumn = number[0] + number[1];
    for (iRow = 2; iRow <= numberRows; iRow++) {
      if (number[iRow]) {
        printf("XX %d columns have %d entries\n", number[iRow], iRow);
        int kColumn = jColumn + number[iRow];
        sortOnOther(row, columnStart, order + jColumn, other, number[iRow],
                    iRow, 0);
        // Now print etc
        if (iRow < 500000) {
          for (int lColumn = jColumn; lColumn < kColumn; lColumn++) {
            iColumn = order[lColumn];
            CoinBigIndex start = columnStart[iColumn];
            if (model->logLevel() == 63) {
              printf("column %d %g <= ", iColumn, columnLower[iColumn]);
              for (CoinBigIndex i = start; i < start + iRow; i++)
                printf("( %d, %g) ", row[i], element[i]);
              printf("<= %g\n", columnUpper[iColumn]);
            }
          }
        }
        jColumn = kColumn;
      }
    }
    delete[] order;
    delete[] other;
    delete[] number;
  }
  // get row copy
  CoinPackedMatrix rowCopy = *matrix;
  rowCopy.reverseOrdering();
  const int *rowLength = rowCopy.getVectorLengths();
  number = new int[numberColumns + 1];
  memset(number, 0, (numberColumns + 1) * sizeof(int));
  if (model->logLevel() > 3) {
    // get column copy
    CoinPackedMatrix columnCopy = *matrix;
    const int *columnLength = columnCopy.getVectorLengths();
    const int *row = columnCopy.getIndices();
    const CoinBigIndex *columnStart = columnCopy.getVectorStarts();
    const double *element = columnCopy.getElements();
    const double *elementByRow = rowCopy.getElements();
    const CoinBigIndex *rowStart = rowCopy.getVectorStarts();
    const int *column = rowCopy.getIndices();
    int nPossibleZeroCost = 0;
    int nPossibleNonzeroCost = 0;
    for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
      int length = columnLength[iColumn];
      if (columnLower[iColumn] < -1.0e30 && columnUpper[iColumn] > 1.0e30) {
        if (length == 1) {
          printf("Singleton free %d - cost %g\n", iColumn, objective[iColumn]);
        } else if (length == 2) {
          int iRow0 = row[columnStart[iColumn]];
          int iRow1 = row[columnStart[iColumn] + 1];
          double element0 = element[columnStart[iColumn]];
          double element1 = element[columnStart[iColumn] + 1];
          int n0 = rowLength[iRow0];
          int n1 = rowLength[iRow1];
          printf("Doubleton free %d - cost %g - %g in %srow with %d entries "
                 "and %g in %srow with %d entries\n",
                 iColumn, objective[iColumn], element0,
                 (rowLower[iRow0] == rowUpper[iRow0]) ? "==" : "", n0, element1,
                 (rowLower[iRow1] == rowUpper[iRow1]) ? "==" : "", n1);
        }
      }
      if (length == 1) {
        int iRow = row[columnStart[iColumn]];
        double value = COIN_DBL_MAX;
        for (CoinBigIndex i = rowStart[iRow];
             i < rowStart[iRow] + rowLength[iRow]; i++) {
          int jColumn = column[i];
          if (jColumn != iColumn) {
            if (value != elementByRow[i]) {
              if (value == COIN_DBL_MAX) {
                value = elementByRow[i];
              } else {
                value = -COIN_DBL_MAX;
                break;
              }
            }
          }
        }
        if (!objective[iColumn]) {
          if (model->logLevel() > 4)
            printf("Singleton %d with no objective in row with %d elements - "
                   "rhs %g,%g\n",
                   iColumn, rowLength[iRow], rowLower[iRow], rowUpper[iRow]);
          nPossibleZeroCost++;
        } else if (value != -COIN_DBL_MAX) {
          if (model->logLevel() > 4)
            printf("Singleton %d (%s) with objective in row %d (%s) with %d "
                   "equal elements - rhs %g,%g\n",
                   iColumn, model->getColumnName(iColumn).c_str(), iRow,
                   model->getRowName(iRow).c_str(), rowLength[iRow],
                   rowLower[iRow], rowUpper[iRow]);
          nPossibleNonzeroCost++;
        }
      }
    }
    if (nPossibleZeroCost || nPossibleNonzeroCost)
      printf("%d singletons with zero cost, %d with valid cost\n",
             nPossibleZeroCost, nPossibleNonzeroCost);
    // look for DW
    int *blockStart =
        new int[2 * (numberRows + numberColumns) + 1 + numberRows];
    int *columnBlock = blockStart + numberRows;
    int *nextColumn = columnBlock + numberColumns;
    int *blockCount = nextColumn + numberColumns;
    int *blockEls = blockCount + numberRows + 1;
    int direction[2] = {-1, 1};
    int bestBreak = -1;
    double bestValue = 0.0;
    int iPass = 0;
    int halfway = (numberRows + 1) / 2;
    int firstMaster = -1;
    int lastMaster = -2;
    while (iPass < 2) {
      int increment = direction[iPass];
      int start = increment > 0 ? 0 : numberRows - 1;
      int stop = increment > 0 ? numberRows : -1;
      int numberBlocks = 0;
      int thisBestBreak = -1;
      double thisBestValue = COIN_DBL_MAX;
      int numberRowsDone = 0;
      int numberMarkedColumns = 0;
      int maximumBlockSize = 0;
      for (int i = 0; i < numberRows + 2 * numberColumns; i++)
        blockStart[i] = -1;
      for (int i = 0; i < numberRows + 1; i++)
        blockCount[i] = 0;
      for (int iRow = start; iRow != stop; iRow += increment) {
        int iBlock = -1;
        for (CoinBigIndex j = rowStart[iRow];
             j < rowStart[iRow] + rowLength[iRow]; j++) {
          int iColumn = column[j];
          int whichColumnBlock = columnBlock[iColumn];
          if (whichColumnBlock >= 0) {
            // column marked
            if (iBlock < 0) {
              // put row in that block
              iBlock = whichColumnBlock;
            } else if (iBlock != whichColumnBlock) {
              // merge
              blockCount[iBlock] += blockCount[whichColumnBlock];
              blockCount[whichColumnBlock] = 0;
              int jColumn = blockStart[whichColumnBlock];
              while (jColumn >= 0) {
                columnBlock[jColumn] = iBlock;
                iColumn = jColumn;
                jColumn = nextColumn[jColumn];
              }
              nextColumn[iColumn] = blockStart[iBlock];
              blockStart[iBlock] = blockStart[whichColumnBlock];
              blockStart[whichColumnBlock] = -1;
            }
          }
        }
        int n = numberMarkedColumns;
        if (iBlock < 0) {
          // new block
          if (rowLength[iRow]) {
            numberBlocks++;
            iBlock = numberBlocks;
            int jColumn = column[rowStart[iRow]];
            columnBlock[jColumn] = iBlock;
            blockStart[iBlock] = jColumn;
            numberMarkedColumns++;
            for (CoinBigIndex j = rowStart[iRow] + 1;
                 j < rowStart[iRow] + rowLength[iRow]; j++) {
              int iColumn = column[j];
              columnBlock[iColumn] = iBlock;
              numberMarkedColumns++;
              nextColumn[jColumn] = iColumn;
              jColumn = iColumn;
            }
            blockCount[iBlock] = numberMarkedColumns - n;
          } else {
            // empty
            iBlock = numberRows;
          }
        } else {
          // put in existing block
          int jColumn = blockStart[iBlock];
          for (CoinBigIndex j = rowStart[iRow];
               j < rowStart[iRow] + rowLength[iRow]; j++) {
            int iColumn = column[j];
            assert(columnBlock[iColumn] < 0 || columnBlock[iColumn] == iBlock);
            if (columnBlock[iColumn] < 0) {
              columnBlock[iColumn] = iBlock;
              numberMarkedColumns++;
              nextColumn[iColumn] = jColumn;
              jColumn = iColumn;
            }
          }
          blockStart[iBlock] = jColumn;
          blockCount[iBlock] += numberMarkedColumns - n;
        }
        maximumBlockSize = std::max(maximumBlockSize, blockCount[iBlock]);
        numberRowsDone++;
        if (thisBestValue * numberRowsDone > maximumBlockSize &&
            numberRowsDone > halfway) {
          thisBestBreak = iRow;
          thisBestValue = static_cast<double>(maximumBlockSize) /
                          static_cast<double>(numberRowsDone);
        }
      }
      if (thisBestBreak == stop)
        thisBestValue = COIN_DBL_MAX;
      iPass++;
      if (iPass == 1) {
        bestBreak = thisBestBreak;
        bestValue = thisBestValue;
      } else {
        if (bestValue < thisBestValue) {
          firstMaster = 0;
          lastMaster = bestBreak;
        } else {
          firstMaster = thisBestBreak; // ? +1
          lastMaster = numberRows;
        }
      }
    }
    if (firstMaster < lastMaster) {
      printf("%d master rows %d <= < %d\n", lastMaster - firstMaster,
             firstMaster, lastMaster);
      for (int i = 0; i < numberRows + 2 * numberColumns; i++)
        blockStart[i] = -1;
      for (int i = firstMaster; i < lastMaster; i++)
        blockStart[i] = -2;
      int firstRow = 0;
      int numberBlocks = -1;
      while (true) {
        for (; firstRow < numberRows; firstRow++) {
          if (blockStart[firstRow] == -1)
            break;
        }
        if (firstRow == numberRows)
          break;
        int nRows = 0;
        numberBlocks++;
        int numberStack = 1;
        blockCount[0] = firstRow;
        while (numberStack) {
          int iRow = blockCount[--numberStack];
          for (CoinBigIndex j = rowStart[iRow];
               j < rowStart[iRow] + rowLength[iRow]; j++) {
            int iColumn = column[j];
            int iBlock = columnBlock[iColumn];
            if (iBlock < 0) {
              columnBlock[iColumn] = numberBlocks;
              for (CoinBigIndex k = columnStart[iColumn];
                   k < columnStart[iColumn] + columnLength[iColumn]; k++) {
                int jRow = row[k];
                int rowBlock = blockStart[jRow];
                if (rowBlock == -1) {
                  nRows++;
                  blockStart[jRow] = numberBlocks;
                  blockCount[numberStack++] = jRow;
                }
              }
            }
          }
        }
        if (!nRows) {
          // empty!!
          numberBlocks--;
        }
        firstRow++;
      }
      // adjust
      numberBlocks++;
      for (int i = 0; i < numberBlocks; i++) {
        blockCount[i] = 0;
        nextColumn[i] = 0;
      }
      int numberEmpty = 0;
      int numberMaster = 0;
      memset(blockEls, 0, numberBlocks * sizeof(int));
      for (int iRow = 0; iRow < numberRows; iRow++) {
        int iBlock = blockStart[iRow];
        if (iBlock >= 0) {
          blockCount[iBlock]++;
          blockEls[iBlock] += rowLength[iRow];
        } else {
          if (iBlock == -2)
            numberMaster++;
          else
            numberEmpty++;
        }
      }
      int numberEmptyColumns = 0;
      int numberMasterColumns = 0;
      for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
        int iBlock = columnBlock[iColumn];
        if (iBlock >= 0) {
          nextColumn[iBlock]++;
        } else {
          if (columnLength[iColumn])
            numberMasterColumns++;
          else
            numberEmptyColumns++;
        }
      }
      int largestRows = 0;
      int largestColumns = 0;
      for (int i = 0; i < numberBlocks; i++) {
        if (blockCount[i] + nextColumn[i] > largestRows + largestColumns) {
          largestRows = blockCount[i];
          largestColumns = nextColumn[i];
        }
      }
      bool useful = true;
      if (numberMaster > halfway || largestRows * 3 > numberRows)
        useful = false;
      printf("%s %d blocks (largest %d,%d), %d master rows (%d empty) out of "
             "%d, %d master columns (%d empty) out of %d\n",
             useful ? "**Useful" : "NoGood", numberBlocks, largestRows,
             largestColumns, numberMaster, numberEmpty, numberRows,
             numberMasterColumns, numberEmptyColumns, numberColumns);
      FILE *fp = NULL;
      bool justIntegers = true;
      bool oneFile = true;
      int logLevel = model->logLevel();
      if (logLevel > 19) {
        logLevel -= 2;
        oneFile = true;
        fp = fopen("fake.bnd", "w");
      }
      if (logLevel == 19)
        justIntegers = false;
      for (int i = 0; i < numberBlocks; i++)
        printf("Block %d has %d rows and %d columns (%d elements)\n", i,
               blockCount[i], nextColumn[i], blockEls[i]);
      if (logLevel >= 17 && logLevel <= 21) {
        int *whichRows = new int[numberRows + numberColumns];
        int *whichColumns = whichRows + numberRows;
        char name[20];
        for (int iBlock = 0; iBlock < numberBlocks; iBlock++) {
          sprintf(name, "block%d.mps", iBlock);
          int nRows = 0;
          for (int iRow = 0; iRow < numberRows; iRow++) {
            if (blockStart[iRow] == iBlock)
              whichRows[nRows++] = iRow;
          }
          int nColumns = 0;
          for (int iColumn = 0; iColumn < numberColumns; iColumn++) {
            if (columnBlock[iColumn] == iBlock)
              whichColumns[nColumns++] = iColumn;
          }
          ClpSimplex subset(model, nRows, whichRows, nColumns, whichColumns);
          for (int jRow = 0; jRow < nRows; jRow++) {
            int iRow = whichRows[jRow];
            std::string name = model->getRowName(iRow);
            subset.setRowName(jRow, name);
          }
          int nInteger = 0;
          for (int jColumn = 0; jColumn < nColumns; jColumn++) {
            int iColumn = whichColumns[jColumn];
            if (model->isInteger(iColumn)) {
              subset.setInteger(jColumn);
              nInteger++;
            }
            std::string name = model->getColumnName(iColumn);
            subset.setColumnName(jColumn, name);
          }
          if (logLevel == 17) {
            subset.writeMps(name, 0, 1);
          } else if (nInteger) {
            OsiClpSolverInterface subset2(&subset);
            CbcModel smallModel(subset2);
            smallModel.branchAndBound();
            const double *solution = smallModel.bestSolution();
            if (solution) {
              if (!oneFile) {
                sprintf(name, "block%d.bnd", iBlock);
                fp = fopen(name, "w");
                assert(fp);
              }
              fprintf(fp, "BBB objective %g for block %d\n",
                      smallModel.getObjValue(), iBlock);
              for (int jColumn = 0; jColumn < nColumns; jColumn++) {
                if (subset.isInteger(jColumn) || !justIntegers)
                  fprintf(fp, " FX BOUND1    %.8s  %g\n",
                          subset.getColumnName(jColumn).c_str(),
                          solution[jColumn]);
              }
              if (!oneFile)
                fclose(fp);
            } else {
              printf("***** Problem is infeasible\n");
              abort();
            }
          }
        }
        if (oneFile)
          fclose(fp);
        delete[] whichRows;
      }
    }
    delete[] blockStart;
  }
  for (iRow = 0; iRow < numberRows; iRow++) {
    int length = rowLength[iRow];
    number[length]++;
  }
  if (number[0])
    printf("** %d rows have no entries\n", number[0]);
  k = 0;
  for (iColumn = 1; iColumn <= numberColumns; iColumn++) {
    if (number[iColumn]) {
      k++;
      printf("%d rows have %d entries\n", number[iColumn], iColumn);
      if (k == kMax)
        break;
    }
  }
  if (k < numberColumns) {
    int kk = k;
    k = 0;
    for (iColumn = numberColumns; iColumn >= 1; iColumn--) {
      if (number[iColumn]) {
        k++;
        if (k == kMax)
          break;
      }
    }
    if (k > kk) {
      printf("\n    .........\n\n");
      iColumn = k;
      k = 0;
      for (; iColumn < numberColumns; iColumn++) {
        if (number[iColumn]) {
          k++;
          printf("%d rows have %d entries\n", number[iColumn], iColumn);
          if (k == kMax)
            break;
        }
      }
    }
  }
  if (model->logLevel() == 63
#ifdef SYM
      || true
#endif
  ) {
    int *column = rowCopy.getMutableIndices();
    const CoinBigIndex *rowStart = rowCopy.getVectorStarts();
    double *element = rowCopy.getMutableElements();
    int *order = new int[numberRows];
    int *other = new int[numberRows];
    for (iRow = 0; iRow < numberRows; iRow++) {
      int length = rowLength[iRow];
      order[iRow] = iRow;
      other[iRow] = length;
      CoinBigIndex start = rowStart[iRow];
      CoinSort_2(column + start, column + start + length, element + start);
    }
    CoinSort_2(other, other + numberRows, order);
    int jRow = number[0] + number[1];
    double *weight = new double[numberRows];
    double *randomColumn = new double[numberColumns + 1];
    double *randomRow = new double[numberRows + 1];
    int *sortRow = new int[numberRows];
    int *possibleRow = new int[numberRows];
    int *backRow = new int[numberRows];
    int *stackRow = new int[numberRows];
    int *sortColumn = new int[numberColumns];
    int *possibleColumn = new int[numberColumns];
    int *backColumn = new int[numberColumns];
    int *backColumn2 = new int[numberColumns];
    int *mapRow = new int[numberRows];
    int *mapColumn = new int[numberColumns];
    int *stackColumn = new int[numberColumns];
    double randomLower = CoinDrand48();
    double randomUpper = CoinDrand48();
    double randomInteger = CoinDrand48();
    CoinBigIndex *startAdd = new CoinBigIndex[numberRows + 1];
    int *columnAdd = new int[2 * numberElements];
    double *elementAdd = new double[2 * numberElements];
    int nAddRows = 0;
    startAdd[0] = 0;
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      randomColumn[iColumn] = CoinDrand48();
      backColumn2[iColumn] = -1;
    }
    for (iColumn = 2; iColumn <= numberColumns; iColumn++) {
      if (number[iColumn]) {
        printf("XX %d rows have %d entries\n", number[iColumn], iColumn);
        int kRow = jRow + number[iColumn];
        sortOnOther(column, rowStart, order + jRow, other, number[iColumn],
                    iColumn, 0);
        // Now print etc
        if (iColumn < 500000) {
          int nLook = 0;
          for (int lRow = jRow; lRow < kRow; lRow++) {
            iRow = order[lRow];
            CoinBigIndex start = rowStart[iRow];
            if (model->logLevel() == 63) {
              printf("row %d %g <= ", iRow, rowLower[iRow]);
              for (CoinBigIndex i = start; i < start + iColumn; i++)
                printf("( %d, %g) ", column[i], element[i]);
              printf("<= %g\n", rowUpper[iRow]);
            }
            int first = column[start];
            double sum = 0.0;
            for (CoinBigIndex i = start; i < start + iColumn; i++) {
              int jColumn = column[i];
              double value = element[i];
              jColumn -= first;
              assert(jColumn >= 0);
              sum += value * randomColumn[jColumn];
            }
            if (rowLower[iRow] > -1.0e30 && rowLower[iRow])
              sum += rowLower[iRow] * randomLower;
            else if (!rowLower[iRow])
              sum += 1.234567e-7 * randomLower;
            if (rowUpper[iRow] < 1.0e30 && rowUpper[iRow])
              sum += rowUpper[iRow] * randomUpper;
            else if (!rowUpper[iRow])
              sum += 1.234567e-7 * randomUpper;
            sortRow[nLook] = iRow;
            randomRow[nLook++] = sum;
            // best way is to number unique elements and bounds and use
            if (fabs(sum) > 1.0e4)
              sum *= 1.0e-6;
            weight[iRow] = sum;
          }
          assert(nLook <= numberRows);
          CoinSort_2(randomRow, randomRow + nLook, sortRow);
          randomRow[nLook] = COIN_DBL_MAX;
          double last = -COIN_DBL_MAX;
          int iLast = -1;
          for (int iLook = 0; iLook < nLook + 1; iLook++) {
            if (randomRow[iLook] > last) {
              if (iLast >= 0) {
                int n = iLook - iLast;
                if (n > 1) {
                  // printf("%d rows possible?\n",n);
                }
              }
              iLast = iLook;
              last = randomRow[iLook];
            }
          }
        }
        jRow = kRow;
      }
    }
    CoinPackedMatrix columnCopy = *matrix;
    const int *columnLength = columnCopy.getVectorLengths();
    const int *row = columnCopy.getIndices();
    const CoinBigIndex *columnStart = columnCopy.getVectorStarts();
    const double *elementByColumn = columnCopy.getElements();
    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
      int length = columnLength[iColumn];
      CoinBigIndex start = columnStart[iColumn];
      double sum = objective[iColumn];
      if (columnLower[iColumn] > -1.0e30 && columnLower[iColumn])
        sum += columnLower[iColumn] * randomLower;
      else if (!columnLower[iColumn])
        sum += 1.234567e-7 * randomLower;
      if (columnUpper[iColumn] < 1.0e30 && columnUpper[iColumn])
        sum += columnUpper[iColumn] * randomUpper;
      else if (!columnUpper[iColumn])
        sum += 1.234567e-7 * randomUpper;
      if (model->isInteger(iColumn))
        sum += 9.87654321e-6 * randomInteger;
      for (CoinBigIndex i = start; i < start + length; i++) {
        int iRow = row[i];
        sum += elementByColumn[i] * weight[iRow];
      }
      sortColumn[iColumn] = iColumn;
      randomColumn[iColumn] = sum;
    }
    {
      CoinSort_2(randomColumn, randomColumn + numberColumns, sortColumn);
      for (iColumn = 0; iColumn < numberColumns; iColumn++) {
        int i = sortColumn[iColumn];
        backColumn[i] = iColumn;
      }
      randomColumn[numberColumns] = COIN_DBL_MAX;
      double last = -COIN_DBL_MAX;
      int iLast = -1;
      for (int iLook = 0; iLook < numberColumns + 1; iLook++) {
        if (randomColumn[iLook] > last) {
          if (iLast >= 0) {
            int n = iLook - iLast;
            if (n > 1) {
              // printf("%d columns possible?\n",n);
            }
            for (int i = iLast; i < iLook; i++) {
              possibleColumn[sortColumn[i]] = n;
            }
          }
          iLast = iLook;
          last = randomColumn[iLook];
        }
      }
      for (iRow = 0; iRow < numberRows; iRow++) {
        CoinBigIndex start = rowStart[iRow];
        double sum = 0.0;
        int length = rowLength[iRow];
        for (CoinBigIndex i = start; i < start + length; i++) {
          int jColumn = column[i];
          double value = element[i];
          jColumn = backColumn[jColumn];
          sum += value * randomColumn[jColumn];
          // if (iColumn==23089||iRow==23729)
          // printf("row %d cola %d colb %d value %g rand %g sum %g\n",
          //   iRow,jColumn,column[i],value,randomColumn[jColumn],sum);
        }
        sortRow[iRow] = iRow;
        randomRow[iRow] = weight[iRow];
        randomRow[iRow] = sum;
      }
      CoinSort_2(randomRow, randomRow + numberRows, sortRow);
      for (iRow = 0; iRow < numberRows; iRow++) {
        int i = sortRow[iRow];
        backRow[i] = iRow;
      }
      randomRow[numberRows] = COIN_DBL_MAX;
      last = -COIN_DBL_MAX;
      iLast = -1;
      // Do backward indices from order
      for (iRow = 0; iRow < numberRows; iRow++) {
        other[order[iRow]] = iRow;
      }
      for (int iLook = 0; iLook < numberRows + 1; iLook++) {
        if (randomRow[iLook] > last) {
          if (iLast >= 0) {
            int n = iLook - iLast;
            if (n > 1) {
              // printf("%d rows possible?\n",n);
              // Within group sort as for original "order"
              for (int i = iLast; i < iLook; i++) {
                int jRow = sortRow[i];
                order[i] = other[jRow];
              }
              CoinSort_2(order + iLast, order + iLook, sortRow + iLast);
            }
            for (int i = iLast; i < iLook; i++) {
              possibleRow[sortRow[i]] = n;
            }
          }
          iLast = iLook;
          last = randomRow[iLook];
        }
      }
      // Temp out
      for (int iLook = 0; iLook < numberRows - 1000000; iLook++) {
        iRow = sortRow[iLook];
        CoinBigIndex start = rowStart[iRow];
        int length = rowLength[iRow];
        int numberPossible = possibleRow[iRow];
        for (CoinBigIndex i = start; i < start + length; i++) {
          int jColumn = column[i];
          if (possibleColumn[jColumn] != numberPossible)
            numberPossible = -1;
        }
        int n = numberPossible;
        if (numberPossible > 1) {
          // printf("pppppossible %d\n",numberPossible);
          for (int jLook = iLook + 1; jLook < iLook + numberPossible; jLook++) {
            int jRow = sortRow[jLook];
            CoinBigIndex start2 = rowStart[jRow];
            assert(numberPossible == possibleRow[jRow]);
            assert(length == rowLength[jRow]);
            for (CoinBigIndex i = start2; i < start2 + length; i++) {
              int jColumn = column[i];
              if (possibleColumn[jColumn] != numberPossible)
                numberPossible = -1;
            }
          }
          if (numberPossible < 2) {
            // switch off
            for (int jLook = iLook; jLook < iLook + n; jLook++)
              possibleRow[sortRow[jLook]] = -1;
          }
          // skip rest
          iLook += n - 1;
        } else {
          possibleRow[iRow] = -1;
        }
      }
      for (int iLook = 0; iLook < numberRows; iLook++) {
        iRow = sortRow[iLook];
        int numberPossible = possibleRow[iRow];
        // Only if any integers
        int numberIntegers = 0;
        CoinBigIndex start = rowStart[iRow];
        int length = rowLength[iRow];
        for (CoinBigIndex i = start; i < start + length; i++) {
          int jColumn = column[i];
          if (model->isInteger(jColumn))
            numberIntegers++;
        }
        if (numberPossible > 1 && !numberIntegers) {
          // printf("possible %d - but no integers\n",numberPossible);
        }
        if (numberPossible > 1 && (numberIntegers || false)) {
          //
          printf("possible %d - %d integers\n", numberPossible, numberIntegers);
          int lastLook = iLook;
          int nMapRow = -1;
          for (int jLook = iLook + 1; jLook < iLook + numberPossible; jLook++) {
            // stop if too many failures
            if (jLook > iLook + 10 && nMapRow < 0)
              break;
            // Create identity mapping
            int i;
            for (i = 0; i < numberRows; i++)
              mapRow[i] = i;
            for (i = 0; i < numberColumns; i++)
              mapColumn[i] = i;
            int offset = jLook - iLook;
            int nStackC = 0;
            // build up row and column mapping
            int nStackR = 1;
            stackRow[0] = iLook;
            bool good = true;
            while (nStackR) {
              nStackR--;
              int look1 = stackRow[nStackR];
              int look2 = look1 + offset;
              assert(randomRow[look1] == randomRow[look2]);
              int row1 = sortRow[look1];
              int row2 = sortRow[look2];
              assert(mapRow[row1] == row1);
              assert(mapRow[row2] == row2);
              mapRow[row1] = row2;
              mapRow[row2] = row1;
              CoinBigIndex start1 = rowStart[row1];
              CoinBigIndex offset2 = rowStart[row2] - start1;
              int length = rowLength[row1];
              assert(length == rowLength[row2]);
              for (CoinBigIndex i = start1; i < start1 + length; i++) {
                int jColumn1 = column[i];
                int jColumn2 = column[i + offset2];
                if (randomColumn[backColumn[jColumn1]] !=
                    randomColumn[backColumn[jColumn2]]) {
                  good = false;
                  break;
                }
                if (mapColumn[jColumn1] == jColumn1) {
                  // not touched
                  assert(mapColumn[jColumn2] == jColumn2);
                  if (jColumn1 != jColumn2) {
                    // Put on stack
                    mapColumn[jColumn1] = jColumn2;
                    mapColumn[jColumn2] = jColumn1;
                    stackColumn[nStackC++] = jColumn1;
                  }
                } else {
                  if (mapColumn[jColumn1] != jColumn2 ||
                      mapColumn[jColumn2] != jColumn1) {
                    // bad
                    good = false;
                    printf("bad col\n");
                    break;
                  }
                }
              }
              if (!good)
                break;
              while (nStackC) {
                nStackC--;
                int iColumn = stackColumn[nStackC];
                int iColumn2 = mapColumn[iColumn];
                assert(iColumn != iColumn2);
                int length = columnLength[iColumn];
                assert(length == columnLength[iColumn2]);
                CoinBigIndex start = columnStart[iColumn];
                CoinBigIndex offset2 = columnStart[iColumn2] - start;
                for (CoinBigIndex i = start; i < start + length; i++) {
                  int iRow = row[i];
                  int iRow2 = row[i + offset2];
                  if (mapRow[iRow] == iRow) {
                    // First (but be careful)
                    if (iRow != iRow2) {
                      // mapRow[iRow]=iRow2;
                      // mapRow[iRow2]=iRow;
                      int iBack = backRow[iRow];
                      int iBack2 = backRow[iRow2];
                      if (randomRow[iBack] == randomRow[iBack2] &&
                          iBack2 - iBack == offset) {
                        stackRow[nStackR++] = iBack;
                      } else {
                        // printf("randomRow diff - weights %g %g\n",
                        //     weight[iRow],weight[iRow2]);
                        // bad
                        good = false;
                        break;
                      }
                    }
                  } else {
                    if (mapRow[iRow] != iRow2 || mapRow[iRow2] != iRow) {
                      // bad
                      good = false;
                      printf("bad row\n");
                      break;
                    }
                  }
                }
                if (!good)
                  break;
              }
            }
            // then check OK
            if (good) {
              for (iRow = 0; iRow < numberRows; iRow++) {
                CoinBigIndex start = rowStart[iRow];
                int length = rowLength[iRow];
                if (mapRow[iRow] == iRow) {
                  for (CoinBigIndex i = start; i < start + length; i++) {
                    int jColumn = column[i];
                    backColumn2[jColumn] = static_cast<int>(i - start);
                  }
                  for (CoinBigIndex i = start; i < start + length; i++) {
                    int jColumn = column[i];
                    if (mapColumn[jColumn] != jColumn) {
                      int jColumn2 = mapColumn[jColumn];
                      CoinBigIndex i2 = backColumn2[jColumn2];
                      if (i2 < 0) {
                        good = false;
                      } else if (element[i] != element[i2 + start]) {
                        good = false;
                      }
                    }
                  }
                  for (CoinBigIndex i = start; i < start + length; i++) {
                    int jColumn = column[i];
                    backColumn2[jColumn] = -1;
                  }
                } else {
                  int row2 = mapRow[iRow];
                  assert(iRow = mapRow[row2]);
                  if (rowLower[iRow] != rowLower[row2] ||
                      rowLower[row2] != rowLower[iRow])
                    good = false;
                  CoinBigIndex offset2 = rowStart[row2] - start;
                  for (CoinBigIndex i = start; i < start + length; i++) {
                    int jColumn = column[i];
                    double value = element[i];
                    int jColumn2 = column[i + offset2];
                    double value2 = element[i + offset2];
                    if (value != value2 || mapColumn[jColumn] != jColumn2 ||
                        mapColumn[jColumn2] != jColumn)
                      good = false;
                  }
                }
              }
              if (good) {
                // check rim
                for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                  if (mapColumn[iColumn] != iColumn) {
                    int iColumn2 = mapColumn[iColumn];
                    if (objective[iColumn] != objective[iColumn2])
                      good = false;
                    if (columnLower[iColumn] != columnLower[iColumn2])
                      good = false;
                    if (columnUpper[iColumn] != columnUpper[iColumn2])
                      good = false;
                    if (model->isInteger(iColumn) != model->isInteger(iColumn2))
                      good = false;
                  }
                }
              }
              if (good) {
                // temp
                if (nMapRow < 0) {
                  // const double * solution = model->primalColumnSolution();
                  // find mapped
                  int nMapColumn = 0;
                  for (int i = 0; i < numberColumns; i++) {
                    if (mapColumn[i] > i)
                      nMapColumn++;
                  }
                  nMapRow = 0;
                  int kRow = -1;
                  for (int i = 0; i < numberRows; i++) {
                    if (mapRow[i] > i) {
                      nMapRow++;
                      kRow = i;
                    }
                  }
                  printf("%d columns, %d rows\n", nMapColumn, nMapRow);
                  if (nMapRow == 1) {
                    CoinBigIndex start = rowStart[kRow];
                    int length = rowLength[kRow];
                    printf("%g <= ", rowLower[kRow]);
                    for (CoinBigIndex i = start; i < start + length; i++) {
                      int jColumn = column[i];
                      if (mapColumn[jColumn] != jColumn)
                        printf("* ");
                      printf("%d,%g ", jColumn, element[i]);
                    }
                    printf("<= %g\n", rowUpper[kRow]);
                  }
                }
                // temp
                int row1 = sortRow[lastLook];
                int row2 = sortRow[jLook];
                lastLook = jLook;
                CoinBigIndex start1 = rowStart[row1];
                CoinBigIndex offset2 = rowStart[row2] - start1;
                int length = rowLength[row1];
                assert(length == rowLength[row2]);
                CoinBigIndex put = startAdd[nAddRows];
                double multiplier = length < 11 ? 2.0 : 1.125;
                double value = 1.0;
                for (CoinBigIndex i = start1; i < start1 + length; i++) {
                  int jColumn1 = column[i];
                  int jColumn2 = column[i + offset2];
                  columnAdd[put] = jColumn1;
                  elementAdd[put++] = value;
                  columnAdd[put] = jColumn2;
                  elementAdd[put++] = -value;
                  value *= multiplier;
                }
                nAddRows++;
                startAdd[nAddRows] = put;
              } else {
                printf("ouch - did not check out as good\n");
              }
            }
          }
          // skip rest
          iLook += numberPossible - 1;
        }
      }
    }
    if (nAddRows) {
      double *lower = new double[nAddRows];
      double *upper = new double[nAddRows];
      int i;
      // const double * solution = model->primalColumnSolution();
      for (i = 0; i < nAddRows; i++) {
        lower[i] = 0.0;
        upper[i] = COIN_DBL_MAX;
      }
      printf("Adding %d rows with %d elements\n", nAddRows, startAdd[nAddRows]);
      // ClpSimplex newModel(*model);
      // newModel.addRows(nAddRows,lower,upper,startAdd,columnAdd,elementAdd);
      // newModel.writeMps("modified.mps");
      delete[] lower;
      delete[] upper;
    }
    delete[] startAdd;
    delete[] columnAdd;
    delete[] elementAdd;
    delete[] order;
    delete[] other;
    delete[] randomColumn;
    delete[] weight;
    delete[] randomRow;
    delete[] sortRow;
    delete[] backRow;
    delete[] possibleRow;
    delete[] sortColumn;
    delete[] backColumn;
    delete[] backColumn2;
    delete[] possibleColumn;
    delete[] mapRow;
    delete[] mapColumn;
    delete[] stackRow;
    delete[] stackColumn;
  }
  delete[] number;
  // Now do breakdown of ranges
  breakdown("Elements", static_cast<int>(numberElements), elementByColumn);
  breakdown("RowLower", numberRows, rowLower);
  breakdown("RowUpper", numberRows, rowUpper);
  breakdown("ColumnLower", numberColumns, columnLower);
  breakdown("ColumnUpper", numberColumns, columnUpper);
  breakdown("Objective", numberColumns, objective);
}

//###########################################################################
//###########################################################################

static bool maskMatches(const int *starts, char **masks, std::string &check) {
  // back to char as I am old fashioned
  const char *checkC = check.c_str();
  size_t length = strlen(checkC);
  while (length > 0 && checkC[length - 1] == ' ')
    length--;
  for (int i = starts[length]; i < starts[length + 1]; i++) {
    char *thisMask = masks[i];
    size_t k;
    for (k = 0; k < length; k++) {
      if (thisMask[k] != '?' && thisMask[k] != checkC[k])
        break;
    }
    if (k == length)
      return true;
  }
  return false;
}

//###########################################################################
//###########################################################################

static void clean(char *temp) {
  char *put = temp;
  while (*put >= ' ')
    put++;
  *put = '\0';
}

//###########################################################################
//###########################################################################

static void generateCode(CbcModel * /*model*/, const char *fileName, int type,
                         int preProcess) {
  // options on code generation
  bool sizecode = (type & 4) != 0;
  type &= 3;
  FILE *fp = fopen(fileName, "r");
  assert(fp);
  int numberLines = 0;
#define MAXLINES 5000
#define MAXONELINE 200
  char line[MAXLINES][MAXONELINE];
  strcpy(line[numberLines++], "0#if defined(_MSC_VER)");
  strcpy(line[numberLines++], "0// Turn off compiler warning about long names");
  strcpy(line[numberLines++], "0#  pragma warning(disable:4786)");
  strcpy(line[numberLines++], "0#endif\n");
  strcpy(line[numberLines++], "0#include <cassert>");
  strcpy(line[numberLines++], "0#include <iomanip>");
  strcpy(line[numberLines++], "0#include \"OsiClpSolverInterface.hpp\"");
  strcpy(line[numberLines++], "0#include \"CbcModel.hpp\"");
  strcpy(line[numberLines++], "0#include \"CbcCutGenerator.hpp\"");
  strcpy(line[numberLines++], "0#include \"CbcStrategy.hpp\"");
  strcpy(line[numberLines++], "0#include \"CglPreProcess.hpp\"");
  strcpy(line[numberLines++], "0#include \"CoinTime.hpp\"");
  if (preProcess > 0)
    strcpy(line[numberLines++],
           "0#include \"CglProbing.hpp\""); // possibly redundant
  // To allow generated 5's to be just before branchAndBound - do rest here
  strcpy(line[numberLines++], "5  cbcModel->initialSolve();");
  strcpy(line[numberLines++], "5  if (clpModel->tightenPrimalBounds()!=0) {");
  strcpy(line[numberLines++], "5    std::cout<<\"Problem is infeasible - "
                              "tightenPrimalBounds!\"<<std::endl;");
  strcpy(line[numberLines++], "5    exit(1);");
  strcpy(line[numberLines++], "5  }");
  strcpy(line[numberLines++], "5  clpModel->dual();  // clean up");
  if (sizecode) {
    // override some settings
    strcpy(line[numberLines++], "5  // compute some things using problem size");
    strcpy(line[numberLines++], "5  cbcModel->setMinimumDrop(std::min(5.0e-2,");
    strcpy(
        line[numberLines++],
        "5       fabs(cbcModel->getMinimizationObjValue())*1.0e-3+1.0e-4));");
    strcpy(line[numberLines++], "5  if (cbcModel->getNumCols()<500)");
    strcpy(line[numberLines++],
           "5    cbcModel->setMaximumCutPassesAtRoot(-100); // always do 100 "
           "if possible");
    strcpy(line[numberLines++], "5  else if (cbcModel->getNumCols()<5000)");
    strcpy(
        line[numberLines++],
        "5    cbcModel->setMaximumCutPassesAtRoot(100); // use minimum drop");
    strcpy(line[numberLines++], "5  else");
    strcpy(line[numberLines++],
           "5    cbcModel->setMaximumCutPassesAtRoot(20);");
    strcpy(line[numberLines++], "5  cbcModel->setMaximumCutPasses(1);");
  }
  if (preProcess <= 0) {
    // no preprocessing or strategy
    if (preProcess) {
      strcpy(line[numberLines++], "5  // Preprocessing using CbcStrategy");
      strcpy(line[numberLines++], "5  CbcStrategyDefault strategy(1,5,5);");
      strcpy(line[numberLines++], "5  strategy.setupPreProcessing(1);");
      strcpy(line[numberLines++], "5  cbcModel->setStrategy(strategy);");
    }
  } else {
    int translate[] = {9999, 0, 0, -1, 2, 3, -2};
    strcpy(line[numberLines++], "5  // Hand coded preprocessing");
    strcpy(line[numberLines++], "5  CglPreProcess process;");
    strcpy(line[numberLines++],
           "5  OsiSolverInterface * saveSolver=cbcModel->solver()->clone();");
    strcpy(line[numberLines++], "5  // Tell solver we are in Branch and Cut");
    strcpy(line[numberLines++],
           "5  saveSolver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo) ;");
    strcpy(line[numberLines++], "5  // Default set of cut generators");
    strcpy(line[numberLines++], "5  CglProbing generator1;");
    strcpy(line[numberLines++], "5  generator1.setUsingObjective(1);");
    strcpy(line[numberLines++], "5  generator1.setMaxPass(3);");
    strcpy(line[numberLines++],
           "5  generator1.setMaxProbeRoot(saveSolver->getNumCols());");
    strcpy(line[numberLines++], "5  generator1.setMaxElements(100);");
    strcpy(line[numberLines++], "5  generator1.setMaxLookRoot(50);");
    strcpy(line[numberLines++], "5  generator1.setRowCuts(3);");
    strcpy(line[numberLines++], "5  // Add in generators");
    strcpy(line[numberLines++], "5  process.addCutGenerator(&generator1);");
    strcpy(line[numberLines++],
           "5  process.messageHandler()->setLogLevel(cbcModel->logLevel());");
    strcpy(line[numberLines++], "5  OsiSolverInterface * solver2 = ");
    sprintf(line[numberLines++],
            "5    process.preProcessNonDefault(*saveSolver,%d,10);",
            translate[preProcess]);
    strcpy(line[numberLines++],
           "5  // Tell solver we are not in Branch and Cut");
    strcpy(
        line[numberLines++],
        "5  saveSolver->setHintParam(OsiDoInBranchAndCut,false,OsiHintDo) ;");
    strcpy(line[numberLines++], "5  if (solver2)");
    strcpy(line[numberLines++],
           "5    solver2->setHintParam(OsiDoInBranchAndCut,false,OsiHintDo) ;");
    strcpy(line[numberLines++], "5  if (!solver2) {");
    strcpy(line[numberLines++],
           "5    std::cout<<\"Pre-processing says infeasible!\"<<std::endl;");
    strcpy(line[numberLines++], "5    exit(1);");
    strcpy(line[numberLines++], "5  } else {");
    strcpy(line[numberLines++],
           "5    std::cout<<\"processed model has \"<<solver2->getNumRows()");
    strcpy(line[numberLines++],
           "5	     <<\" rows, \"<<solver2->getNumCols()");
    strcpy(line[numberLines++],
           "5	     <<\" columns and \"<<solver2->getNumElements()");
    strcpy(line[numberLines++],
           "5	     <<\" elements\"<<solver2->getNumElements()<<std::endl;");
    strcpy(line[numberLines++], "5  }");
    strcpy(line[numberLines++], "5  // we have to keep solver2 so pass clone");
    strcpy(line[numberLines++], "5  solver2 = solver2->clone();");
    strcpy(line[numberLines++], "5  cbcModel->assignSolver(solver2);");
    strcpy(line[numberLines++], "5  cbcModel->initialSolve();");
  }
  while (fgets(line[numberLines], MAXONELINE, fp)) {
    assert(numberLines < MAXLINES);
    clean(line[numberLines]);
    numberLines++;
  }
  fclose(fp);
  strcpy(line[numberLines++], "0\nint main (int argc, const char *argv[])\n{");
  strcpy(line[numberLines++], "0  OsiClpSolverInterface solver1;");
  strcpy(line[numberLines++], "0  int status=1;");
  strcpy(line[numberLines++], "0  if (argc<2)");
  strcpy(line[numberLines++],
         "0    std::cout<<\"Please give file name\"<<std::endl;");
  strcpy(line[numberLines++], "0  else");
  strcpy(line[numberLines++], "0    status=solver1.readMps(argv[1],\"\");");
  strcpy(line[numberLines++], "0  if (status) {");
  strcpy(line[numberLines++],
         "0    std::cout<<\"Bad readMps \"<<argv[1]<<std::endl;");
  strcpy(line[numberLines++], "0    exit(1);");
  strcpy(line[numberLines++], "0  }\n");
  strcpy(line[numberLines++], "0  double time1 = CoinCpuTime();");
  strcpy(line[numberLines++], "0  CbcModel model(solver1);");
  strcpy(line[numberLines++], "0  // Now do requested saves and modifications");
  strcpy(line[numberLines++], "0  CbcModel * cbcModel = & model;");
  strcpy(line[numberLines++],
         "0  OsiSolverInterface * osiModel = model.solver();");
  strcpy(line[numberLines++],
         "0  OsiClpSolverInterface * osiclpModel = dynamic_cast< "
         "OsiClpSolverInterface*> (osiModel);");
  strcpy(line[numberLines++],
         "0  ClpSimplex * clpModel = osiclpModel->getModelPtr();");
  // add in comments about messages
  strcpy(line[numberLines++],
         "3  // You can save some time by switching off message building");
  strcpy(line[numberLines++],
         "3  // clpModel->messagesPointer()->setDetailMessages(100,10000,(int "
         "*) NULL);");
  // add in actual solve
  strcpy(line[numberLines++], "5  cbcModel->branchAndBound();");
  strcpy(
      line[numberLines++],
      "8  std::cout<<argv[1]<<\" took \"<<CoinCpuTime()-time1<<\" seconds, \"");
  strcpy(line[numberLines++],
         "8	   <<cbcModel->getNodeCount()<<\" nodes with objective \"");
  strcpy(line[numberLines++], "8	   <<cbcModel->getObjValue()");
  strcpy(line[numberLines++], "8	   <<(!cbcModel->status() ? \" "
                              "Finished\" : \" Not finished\")");
  strcpy(line[numberLines++], "8	   <<std::endl;");
  strcpy(line[numberLines++], "5  // For best solution");
  strcpy(line[numberLines++], "5  int numberColumns = solver1.getNumCols();");
  strcpy(line[numberLines++],
         "5  if (cbcModel->getMinimizationObjValue()<1.0e50) {");
  if (preProcess > 0) {
    strcpy(line[numberLines++], "5    // post process");
    strcpy(line[numberLines++],
           "5    process.postProcess(*cbcModel->solver());");
    strcpy(line[numberLines++], "5    // Solution now back in saveSolver");
    strcpy(line[numberLines++], "5    cbcModel->assignSolver(saveSolver);");
    strcpy(line[numberLines++], "5    "
                                "memcpy(cbcModel->bestSolution(),cbcModel->"
                                "solver()->getColSolution(),");
    strcpy(line[numberLines++], "5	   numberColumns*sizeof(double));");
  }
  strcpy(line[numberLines++], "5    // put back in original solver");
  strcpy(line[numberLines++],
         "5    solver1.setColSolution(cbcModel->bestSolution());");
  strcpy(line[numberLines++],
         "5    const double * solution = solver1.getColSolution();");
  strcpy(line[numberLines++], "8  \n  // Now you would use solution etc etc\n");
  strcpy(line[numberLines++], "5");
  strcpy(line[numberLines++],
         "5    // Get names from solver1 (as OsiSolverInterface may lose)");
  strcpy(line[numberLines++], "5    std::vector<std::string> columnNames = "
                              "*solver1.getModelPtr()->columnNames();");
  strcpy(line[numberLines++], "5    ");
  strcpy(line[numberLines++], "5    int iColumn;");
  strcpy(line[numberLines++], "5    "
                              "std::cout<<std::setiosflags(std::ios::fixed|std:"
                              ":ios::showpoint)<<std::setw(14);");
  strcpy(line[numberLines++], "5    ");
  strcpy(
      line[numberLines++],
      "5    std::cout<<\"--------------------------------------\"<<std::endl;");
  strcpy(line[numberLines++],
         "5    for (iColumn=0;iColumn<numberColumns;iColumn++) {");
  strcpy(line[numberLines++], "5      double value=solution[iColumn];");
  strcpy(line[numberLines++],
         "5      if (fabs(value)>1.0e-7&&solver1.isInteger(iColumn)) ");
  strcpy(line[numberLines++],
         "5	std::cout<<std::setw(6)<<iColumn<<\" \"");
  strcpy(line[numberLines++],
         "5                 <<columnNames[iColumn]<<\" \"");
  strcpy(line[numberLines++], "5                 <<value<<std::endl;");
  strcpy(line[numberLines++], "5    }");
  strcpy(
      line[numberLines++],
      "5    std::cout<<\"--------------------------------------\"<<std::endl;");
  strcpy(line[numberLines++], "5  ");
  strcpy(line[numberLines++], "5    "
                              "std::cout<<std::resetiosflags(std::ios::fixed|"
                              "std::ios::showpoint|std::ios::scientific);");
  strcpy(line[numberLines++], "5  }");
  strcpy(line[numberLines++], "8  return 0;\n}");
  fp = fopen(fileName, "w");
  assert(fp);

  int wanted[9];
  memset(wanted, 0, sizeof(wanted));
  wanted[0] = wanted[3] = wanted[5] = wanted[8] = 1;
  if (type > 0)
    wanted[1] = wanted[6] = 1;
  if (type > 1)
    wanted[2] = wanted[4] = wanted[7] = 1;
  std::string header[9] = {"",
                           "Save values",
                           "Redundant save of default values",
                           "Set changed values",
                           "Redundant set default values",
                           "Solve",
                           "Restore values",
                           "Redundant restore values",
                           "Finish up"};
  for (int iType = 0; iType < 9; iType++) {
    if (!wanted[iType])
      continue;
    int n = 0;
    int iLine;
    for (iLine = 0; iLine < numberLines; iLine++) {
      if (line[iLine][0] == '0' + iType) {
        if (!n && header[iType] != "")
          fprintf(fp, "\n  // %s\n\n", header[iType].c_str());
        n++;
        // skip save and clp as cloned
        if (!strstr(line[iLine], "save") ||
            (!strstr(line[iLine], "clpMo") && !strstr(line[iLine], "_Osi")))
          fprintf(fp, "%s\n", line[iLine] + 1);
      }
    }
  }
  fclose(fp);
  printf("C++ file written to %s\n", fileName);
}
// To allow old way of CbcMain1
int CbcMain1(int argc, const char *argv[],
  CbcModel &model,
	     CbcParameters &parameterData,
             int callBack(CbcModel *currentSolver, int whereFrom),
	     ampl_info *info)
{
  std::deque<std::string> inputQueue;
  CoinParamUtils::formInputQueue(inputQueue, "cbc", argc, const_cast< char ** >(argv));
  return CbcMain1(inputQueue,model,parameterData,callBack,info);
}

int CbcMain1(int argc, const char *argv[],
  CbcModel &model,
             int callBack(CbcModel *currentSolver, int whereFrom),
	     CbcParameters &parameterData)
{
  std::deque<std::string> inputQueue;
  CoinParamUtils::formInputQueue(inputQueue, "cbc", argc, const_cast< char ** >(argv));
  return CbcMain1(inputQueue,model,parameterData,callBack,NULL);
}
//###########################################################################
// Print a general message
//###########################################################################

static void printGeneralMessage(CbcModel &model, const char *message) {
#ifndef DISALLOW_PRINTING
  model.messageHandler()->message(CBC_FPUMP1, model.messages())
      << message << CoinMessageEol;
#endif
}
#ifdef CBC_HAS_NAUTY
#include "CbcSymmetry.hpp"
// returns number of constraints added
static int nautiedConstraints(CbcModel &model, int maxPass) {
  bool changed = true;
  int numberAdded = 0;
  int numberPasses = 0;
  int changeType = 0; //(more2&(128|256))>>7;
  OsiSolverInterface *solverOriginal = model.solver();
  //#define REALLY_CHANGE
#ifdef REALLY_CHANGE
  OsiSolverInterface *solver = solverOriginal;
#else
  int numberOriginalRows = solverOriginal->getNumRows();
  OsiSolverInterface *solver = solverOriginal->clone();
  model.swapSolver(solver);
#endif
  while (changed) {
    changed = false;
    CbcSymmetry symmetryInfo;
    // symmetryInfo.setModel(&model);
    // for now strong is just on counts - use user option
    // int maxN=5000000;
    // OsiSolverInterface * solver = model.solver();
    symmetryInfo.setupSymmetry(&model);
    int numberGenerators;
    if (!numberPasses) {
      numberGenerators = symmetryInfo.statsOrbits(&model, 0);
    } else {
      numberGenerators = symmetryInfo.getNtyInfo()->getNumGenerators();
    }
    if (numberGenerators) {
      const double * solution = solver->getColSolution(); 
      // symmetryInfo.Print_Orbits();
      int numberUsefulOrbits = symmetryInfo.numberUsefulOrbits();
      if (numberUsefulOrbits) {
        symmetryInfo.Compute_Symmetry();
        symmetryInfo.fillOrbits(/*true*/);
        const int *orbits = symmetryInfo.whichOrbit();
        int numberUsefulOrbits = symmetryInfo.numberUsefulOrbits();
        int *counts = new int[2*numberUsefulOrbits];
	int *active = counts+numberUsefulOrbits;
        memset(counts, 0, 2*numberUsefulOrbits * sizeof(int));
        int numberColumns = solver->getNumCols();
        int numberUseful = 0;
        if (changeType == 1) {
          // just 0-1
          for (int i = 0; i < numberColumns; i++) {
            int iOrbit = orbits[i];
            if (iOrbit >= 0) {
              if (solver->isBinary(i)) {
                counts[iOrbit]++;
                numberUseful++;
		double value = solution[i];
		if (fabs(value-floor(value+0.5))>1.0e-5)
		  active[iOrbit]++;
              }
            }
          }
        } else if (changeType == 2) {
          // just integer
          for (int i = 0; i < numberColumns; i++) {
            int iOrbit = orbits[i];
            if (iOrbit >= 0) {
              if (solver->isInteger(i)) {
                counts[iOrbit]++;
                numberUseful++;
		double value = solution[i];
		if (fabs(value-floor(value+0.5))>1.0e-5)
		  active[iOrbit]++;
              }
            }
          }
        } else {
          // all
          for (int i = 0; i < numberColumns; i++) {
            int iOrbit = orbits[i];
            if (iOrbit >= 0) {
              counts[iOrbit]++;
              numberUseful++;
	      double value = solution[i];
	      if (fabs(value-floor(value+0.5))>1.0e-5)
		active[iOrbit]++;
            }
          }
        }
        int iOrbit = -1;
#define LONGEST 1
#if LONGEST
        // choose longest (and most active)
        int maxOrbit = 0;
	int maxActive = 0;
        for (int i = 0; i < numberUsefulOrbits; i++) {
          if (counts[i] > maxOrbit) {
            maxOrbit = counts[i];
	    maxActive = active[i];
            iOrbit = i;
          } else if (counts[i] == maxOrbit && active[i] > maxActive) {
	    maxActive = active[i];
            iOrbit = i;
          }
        }
#else
        // choose closest to 2
        int minOrbit = numberColumns + 1;
        for (int i = 0; i < numberUsefulOrbits; i++) {
          if (counts[i] > 1 && counts[i] < minOrbit) {
            minOrbit = counts[i];
            iOrbit = i;
          }
        }
#endif
        delete[] counts;
        if (!numberUseful)
          break;
        // take largest
        const double *solution = solver->getColSolution();
        double *size = new double[numberColumns];
        int *which = new int[numberColumns];
        int nIn = 0;
        for (int i = 0; i < numberColumns; i++) {
          if (orbits[i] == iOrbit) {
            size[nIn] = -solution[i];
            which[nIn++] = i;
          }
        }
        if (nIn > 1) {
          // printf("Using orbit length %d\n",nIn);
          CoinSort_2(size, size + nIn, which);
          size[0] = 1.0;
          size[1] = -1.0;
#if LONGEST == 0
          solver->addRow(2, which, size, 0.0, COIN_DBL_MAX);
          numberAdded++;
#elif LONGEST == 1
          for (int i = 0; i < nIn - 1; i++) {
            solver->addRow(2, which + i, size, 0.0, COIN_DBL_MAX);
            numberAdded++;
          }
#else
          for (int i = 0; i < nIn - 1; i++) {
            solver->addRow(2, which, size, 0.0, COIN_DBL_MAX);
            which[1] = which[2 + i];
            numberAdded++;
          }
#endif
          numberPasses++;
          if (numberPasses < maxPass)
            changed = true;
        }
        delete[] size;
        delete[] which;
      }
    }
  }
  if (numberAdded) {
    char general[100];
    if (numberPasses < maxPass)
      sprintf(general, "%d constraints added in %d passes", numberAdded,
              numberPasses);
    else
      sprintf(
          general,
          "%d constraints added in %d passes (maximum) - must be better way",
          numberAdded, numberPasses);
    model.messageHandler()->message(CBC_GENERAL, model.messages())
        << general << CoinMessageEol;
#ifdef SAVE_NAUTY
    OsiClpSolverInterface *clpSolver =
        dynamic_cast<OsiClpSolverInterface *>(solver);
    ClpSimplex *lpSolver = clpSolver->getModelPtr();
    char name[100];
    strcpy(name, lpSolver->problemName().c_str());
    strcat(name, "_nauty");
    printf("saving model on %s\n", name);
    solver->writeMps(name);
#endif
  }
#ifndef REALLY_CHANGE
  CbcRowCuts *globalCuts = model.globalCuts();
  int numberRows = solver->getNumRows();
  if (numberRows > numberOriginalRows) {
    const CoinPackedMatrix *rowCopy = solver->getMatrixByRow();
    const int *column = rowCopy->getIndices();
    const int *rowLength = rowCopy->getVectorLengths();
    const CoinBigIndex *rowStart = rowCopy->getVectorStarts();
    const double *elements = rowCopy->getElements();
    const double *rowLower = solver->getRowLower();
    const double *rowUpper = solver->getRowUpper();
    for (int iRow = numberOriginalRows; iRow < numberRows; iRow++) {
      OsiRowCut rc;
      rc.setLb(rowLower[iRow]);
      rc.setUb(rowUpper[iRow]);
      CoinBigIndex start = rowStart[iRow];
      rc.setRow(rowLength[iRow], column + start, elements + start, false);
      globalCuts->addCutIfNotDuplicate(rc);
    }
  }
  model.swapSolver(solverOriginal);
  delete solver;
#endif
  return numberAdded;
}
#endif

/*
  Version 1.00.00 November 16 2005.
  This is to stop me (JJF) messing about too much.
  Tuning changes should be noted here.
  The testing next version may be activated by CBC_NEXT_VERSION
  This applies to OsiClp, Clp etc
  Version 1.00.01 November 24 2005
  Added several classes for advanced users.  This can't affect code (if you
  don't use it) Made some tiny changes (for N way branching) which should not
  change anything. CbcNWay object class - for N way branching this also allows
  use of CbcConsequence class. CbcBranchAllDifferent object class - for
  branching on general integer variables to stop them having same value so
  branches are x >= y+1 and x <= y-1. Added two new Cgl classes -
  CglAllDifferent which does column fixing (too slowly) and CglStored which just
  has a list of cuts which can be activated. Modified preprocess option to SOS
  Version 1.00.02 December 9 2005
  Added use of CbcStrategy to do clean preprocessing
  Added use of referenceSolver for cleaner repetition of Cbc
  Version 1.01.00 February 2 2006
  Added first try at Ampl interface
  Version 1.04 June 2007
  Goes parallel
  Version 2.00 September 2007
  Improvements to feaspump
  Source code changes so up to 2.0
*/
