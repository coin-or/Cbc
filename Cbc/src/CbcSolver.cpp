/* $Id: CbcSolver.cpp 1240 2009-10-02 18:41:44Z forrest $ */
// Copyright (C) 2007, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

/*! \file CbcSolver.cpp
    \brief Second level routines for the cbc stand-alone solver.
*/

#include "CbcConfig.h"
#include "CoinPragma.hpp"

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <iostream>

#include "CoinPragma.hpp"
#include "CoinHelperFunctions.hpp"

#include "CoinMpsIO.hpp"
#include "CoinModel.hpp"

#include "ClpFactorization.hpp"
#include "ClpQuadraticObjective.hpp"
#include "CoinTime.hpp"
#include "ClpSimplex.hpp"
#include "ClpSimplexOther.hpp"
#include "ClpSolve.hpp"
#include "ClpMessage.hpp"
#include "ClpPackedMatrix.hpp"
#include "ClpPlusMinusOneMatrix.hpp"
#include "ClpNetworkMatrix.hpp"
#include "ClpDualRowSteepest.hpp"
#include "ClpDualRowDantzig.hpp"
#include "ClpLinearObjective.hpp"
#include "ClpPrimalColumnSteepest.hpp"
#include "ClpPrimalColumnDantzig.hpp"

#include "ClpPresolve.hpp"
#ifndef COIN_HAS_CBC
#define COIN_HAS_CBC
#endif
#include "CbcOrClpParam.hpp"
#include "OsiRowCutDebugger.hpp"
#include "OsiChooseVariable.hpp"
#include "OsiAuxInfo.hpp"

#include "CbcSolverHeuristics.hpp"
#ifdef COIN_HAS_GLPK
#include "glpk.h"
extern glp_tran* cbc_glp_tran;
extern glp_prob* cbc_glp_prob;
#else
#define GLP_UNDEF 1
#define GLP_FEAS 2
#define GLP_INFEAS 3
#define GLP_NOFEAS 4
#define GLP_OPT 5
#endif

//#define USER_HAS_FAKE_CLP
//#define USER_HAS_FAKE_CBC

//#define CLP_MALLOC_STATISTICS

#ifdef CLP_MALLOC_STATISTICS
#include <malloc.h>
#include <exception>
#include <new>
#include "stolen_from_ekk_malloc.cpp"
static double malloc_times = 0.0;
static double malloc_total = 0.0;
static int malloc_amount[] = {0, 32, 128, 256, 1024, 4096, 16384, 65536, 262144, INT_MAX};
static int malloc_n = 10;
double malloc_counts[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
bool malloc_counts_on = true;
void * operator new (size_t size) throw (std::bad_alloc)
{
    malloc_times ++;
    malloc_total += size;
    int i;
    for (i = 0; i < malloc_n; i++) {
        if ((int) size <= malloc_amount[i]) {
            malloc_counts[i]++;
            break;
        }
    }
# ifdef DEBUG_MALLOC
    void *p;
    if (malloc_counts_on)
        p = stolen_from_ekk_mallocBase(size);
    else
        p = malloc(size);
# else
    void * p = malloc(size);
# endif
    //char * xx = (char *) p;
    //memset(xx,0,size);
    // Initialize random seed
    //CoinSeedRandom(987654321);
    return p;
}
void operator delete (void *p) throw()
{
# ifdef DEBUG_MALLOC
    if (malloc_counts_on)
        stolen_from_ekk_freeBase(p);
    else
        free(p);
# else
    free(p);
# endif
}
static void malloc_stats2()
{
    double average = malloc_total / malloc_times;
    printf("count %g bytes %g - average %g\n", malloc_times, malloc_total, average);
    for (int i = 0; i < malloc_n; i++)
        printf("%g ", malloc_counts[i]);
    printf("\n");
    malloc_times = 0.0;
    malloc_total = 0.0;
    memset(malloc_counts, 0, sizeof(malloc_counts));
    // print results
}
#else	//CLP_MALLOC_STATISTICS
//void stolen_from_ekk_memory(void * dummy,int type)
//{
//}
//bool malloc_counts_on=false;
#endif	//CLP_MALLOC_STATISTICS

//#define DMALLOC
#ifdef DMALLOC
#include "dmalloc.h"
#endif

#ifdef WSSMP_BARRIER
#define FOREIGN_BARRIER
#endif

#ifdef UFL_BARRIER
#define FOREIGN_BARRIER
#endif

#ifdef TAUCS_BARRIER
#define FOREIGN_BARRIER
#endif

static int initialPumpTune = -1;
#include "CoinWarmStartBasis.hpp"

#include "OsiSolverInterface.hpp"
#include "OsiCuts.hpp"
#include "OsiRowCut.hpp"
#include "OsiColCut.hpp"

#ifndef COIN_HAS_LINK
#define COIN_HAS_LINK
#endif
#ifdef COIN_HAS_LINK
#include "CbcLinked.hpp"
#endif

#include "CglPreProcess.hpp"
#include "CglCutGenerator.hpp"
#include "CglGomory.hpp"
#include "CglProbing.hpp"
#include "CglKnapsackCover.hpp"
#include "CglRedSplit.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglTwomir.hpp"
#include "CglDuplicateRow.hpp"
#include "CglStored.hpp"
#include "CglLandP.hpp"
#include "CglResidualCapacity.hpp"

#ifdef ZERO_HALF_CUTS
#include "CglZeroHalf.hpp"
#endif
//#define CGL_WRITEMPS 
#ifdef CGL_WRITEMPS
extern double * debugSolution;
extern int debugNumberColumns;
#endif
#include "CbcModel.hpp"
#include "CbcHeuristic.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcHeuristicPivotAndFix.hpp"
//#include "CbcHeuristicPivotAndComplement.hpp"
#include "CbcHeuristicRandRound.hpp"
#include "CbcHeuristicGreedy.hpp"
#include "CbcHeuristicFPump.hpp"
#include "CbcHeuristicRINS.hpp"
#include "CbcHeuristicDiveCoefficient.hpp"
#include "CbcHeuristicDiveFractional.hpp"
#include "CbcHeuristicDiveGuided.hpp"
#include "CbcHeuristicDiveVectorLength.hpp"
#include "CbcHeuristicDivePseudoCost.hpp"
#include "CbcHeuristicDiveLineSearch.hpp"
#include "CbcTreeLocal.hpp"
#include "CbcCompareActual.hpp"
#include "CbcBranchActual.hpp"
#include "CbcBranchLotsize.hpp"
#include  "CbcOrClpParam.hpp"
#include  "CbcCutGenerator.hpp"
#include  "CbcStrategy.hpp"
#include "CbcBranchCut.hpp"

#include "OsiClpSolverInterface.hpp"

#include "CbcSolverAnalyze.hpp"
#include "CbcSolverExpandKnapsack.hpp"

#include "CbcSolver.hpp"

//#define IN_BRANCH_AND_BOUND (0x01000000|262144)
#define IN_BRANCH_AND_BOUND (0x01000000|262144|128|1024|2048)
//#define IN_BRANCH_AND_BOUND (0x01000000|262144|128)

/*
  CbcStopNow class definitions.
*/

CbcStopNow::CbcStopNow()
{
}
CbcStopNow::~CbcStopNow()
{
}
// Copy constructor
CbcStopNow::CbcStopNow ( const CbcStopNow & )
{
}
// Assignment operator
CbcStopNow &
CbcStopNow::operator=(const CbcStopNow & rhs)
{
    if (this != &rhs) {
    }
    return *this;
}
// Clone
CbcStopNow *
CbcStopNow::clone() const
{
    return new CbcStopNow(*this);
}

/*
  CbcUser class definitions.
*/

// User stuff (base class)
CbcUser::CbcUser()
        : coinModel_(NULL),
        userName_("null")
{
}
CbcUser::~CbcUser()
{
    delete coinModel_;
}
// Copy constructor
CbcUser::CbcUser ( const CbcUser & rhs)
{
    if (rhs.coinModel_)
        coinModel_ = new CoinModel(*rhs.coinModel_);
    else
        coinModel_ = NULL;
    userName_ = rhs.userName_;
}
// Assignment operator
CbcUser &
CbcUser::operator=(const CbcUser & rhs)
{
    if (this != &rhs) {
        if (rhs.coinModel_)
            coinModel_ = new CoinModel(*rhs.coinModel_);
        else
            coinModel_ = NULL;
        userName_ = rhs.userName_;
    }
    return *this;
}

/*
  CbcSolver class definitions
*/

CbcSolver::CbcSolver()
        : babModel_(NULL),
        userFunction_(NULL),
        statusUserFunction_(NULL),
        originalSolver_(NULL),
        originalCoinModel_(NULL),
        cutGenerator_(NULL),
        numberUserFunctions_(0),
        numberCutGenerators_(0),
        startTime_(CoinCpuTime()),
        parameters_(NULL),
        numberParameters_(0),
        doMiplib_(false),
        noPrinting_(false),
        readMode_(1)
{
    callBack_ = new CbcStopNow();
    fillParameters();
}
CbcSolver::CbcSolver(const OsiClpSolverInterface & solver)
        : babModel_(NULL),
        userFunction_(NULL),
        statusUserFunction_(NULL),
        originalSolver_(NULL),
        originalCoinModel_(NULL),
        cutGenerator_(NULL),
        numberUserFunctions_(0),
        numberCutGenerators_(0),
        startTime_(CoinCpuTime()),
        parameters_(NULL),
        numberParameters_(0),
        doMiplib_(false),
        noPrinting_(false),
        readMode_(1)
{
    callBack_ = new CbcStopNow();
    model_ = CbcModel(solver);
    fillParameters();
}
CbcSolver::CbcSolver(const CbcModel & solver)
        : babModel_(NULL),
        userFunction_(NULL),
        statusUserFunction_(NULL),
        originalSolver_(NULL),
        originalCoinModel_(NULL),
        cutGenerator_(NULL),
        numberUserFunctions_(0),
        numberCutGenerators_(0),
        startTime_(CoinCpuTime()),
        parameters_(NULL),
        numberParameters_(0),
        doMiplib_(false),
        noPrinting_(false),
        readMode_(1)
{
    callBack_ = new CbcStopNow();
    model_ = solver;
    fillParameters();
}
CbcSolver::~CbcSolver()
{
    int i;
    for (i = 0; i < numberUserFunctions_; i++)
        delete userFunction_[i];
    delete [] userFunction_;
    for (i = 0; i < numberCutGenerators_; i++)
        delete cutGenerator_[i];
    delete [] cutGenerator_;
    delete [] statusUserFunction_;
    delete originalSolver_;
    delete originalCoinModel_;
    delete babModel_;
    delete [] parameters_;
    delete callBack_;
}
// Copy constructor
CbcSolver::CbcSolver ( const CbcSolver & rhs)
        : model_(rhs.model_),
        babModel_(NULL),
        userFunction_(NULL),
        statusUserFunction_(NULL),
        numberUserFunctions_(rhs.numberUserFunctions_),
        startTime_(CoinCpuTime()),
        parameters_(NULL),
        numberParameters_(rhs.numberParameters_),
        doMiplib_(rhs.doMiplib_),
        noPrinting_(rhs.noPrinting_),
        readMode_(rhs.readMode_)
{
    fillParameters();
    if (rhs.babModel_)
        babModel_ = new CbcModel(*rhs.babModel_);
    userFunction_ = new CbcUser * [numberUserFunctions_];
    int i;
    for (i = 0; i < numberUserFunctions_; i++)
        userFunction_[i] = rhs.userFunction_[i]->clone();
    for (i = 0; i < numberParameters_; i++)
        parameters_[i] = rhs.parameters_[i];
    for (i = 0; i < numberCutGenerators_; i++)
        cutGenerator_[i] = rhs.cutGenerator_[i]->clone();
    callBack_ = rhs.callBack_->clone();
    originalSolver_ = NULL;
    if (rhs.originalSolver_) {
        OsiSolverInterface * temp = rhs.originalSolver_->clone();
        originalSolver_ = dynamic_cast<OsiClpSolverInterface *> (temp);
        assert (originalSolver_);
    }
    originalCoinModel_ = NULL;
    if (rhs.originalCoinModel_)
        originalCoinModel_ = new CoinModel(*rhs.originalCoinModel_);
}
// Assignment operator
CbcSolver &
CbcSolver::operator=(const CbcSolver & rhs)
{
    if (this != &rhs) {
        int i;
        for (i = 0; i < numberUserFunctions_; i++)
            delete userFunction_[i];
        delete [] userFunction_;
        for (i = 0; i < numberCutGenerators_; i++)
            delete cutGenerator_[i];
        delete [] cutGenerator_;
        delete [] statusUserFunction_;
        delete originalSolver_;
        delete originalCoinModel_;
        statusUserFunction_ = NULL;
        delete babModel_;
        delete [] parameters_;
        delete callBack_;
        numberUserFunctions_ = rhs.numberUserFunctions_;
        startTime_ = rhs.startTime_;
        numberParameters_ = rhs.numberParameters_;
        for (i = 0; i < numberParameters_; i++)
            parameters_[i] = rhs.parameters_[i];
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
        userFunction_ = new CbcUser * [numberUserFunctions_];
        for (i = 0; i < numberUserFunctions_; i++)
            userFunction_[i] = rhs.userFunction_[i]->clone();
        callBack_ = rhs.callBack_->clone();
        originalSolver_ = NULL;
        if (rhs.originalSolver_) {
            OsiSolverInterface * temp = rhs.originalSolver_->clone();
            originalSolver_ = dynamic_cast<OsiClpSolverInterface *> (temp);
            assert (originalSolver_);
        }
        originalCoinModel_ = NULL;
        if (rhs.originalCoinModel_)
            originalCoinModel_ = new CoinModel(*rhs.originalCoinModel_);
    }
    return *this;
}
// Get int value
int CbcSolver::intValue(CbcOrClpParameterType type) const
{
    return parameters_[whichParam(type, numberParameters_, parameters_)].intValue();
}
// Set int value
void CbcSolver::setIntValue(CbcOrClpParameterType type, int value)
{
    parameters_[whichParam(type, numberParameters_, parameters_)].setIntValue(value);
}
// Get double value
double CbcSolver::doubleValue(CbcOrClpParameterType type) const
{
    return parameters_[whichParam(type, numberParameters_, parameters_)].doubleValue();
}
// Set double value
void CbcSolver::setDoubleValue(CbcOrClpParameterType type, double value)
{
    parameters_[whichParam(type, numberParameters_, parameters_)].setDoubleValue(value);
}
// User function (NULL if no match)
CbcUser * CbcSolver::userFunction(const char * name) const
{
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
void CbcSolver::fillParameters()
{
    int maxParam = 200;
    CbcOrClpParam * parameters = new CbcOrClpParam [maxParam];
    numberParameters_ = 0 ;
    establishParams(numberParameters_, parameters) ;
    assert (numberParameters_ <= maxParam);
    parameters_ = new CbcOrClpParam [numberParameters_];
    int i;
    for (i = 0; i < numberParameters_; i++)
        parameters_[i] = parameters[i];
    delete [] parameters;
    const char dirsep =  CoinFindDirSeparator();
    std::string directory;
    std::string dirSample;
    std::string dirNetlib;
    std::string dirMiplib;
    if (dirsep == '/') {
        directory = "./";
        dirSample = "../../Data/Sample/";
        dirNetlib = "../../Data/Netlib/";
        dirMiplib = "../../Data/miplib3/";
    } else {
        directory = ".\\";
        dirSample = "..\\..\\..\\..\\Data\\Sample\\";
        dirNetlib = "..\\..\\..\\..\\Data\\Netlib\\";
        dirMiplib = "..\\..\\..\\..\\Data\\miplib3\\";
    }
    std::string defaultDirectory = directory;
    std::string importFile = "";
    std::string exportFile = "default.mps";
    std::string importBasisFile = "";
    std::string importPriorityFile = "";
    std::string debugFile = "";
    std::string printMask = "";
    std::string exportBasisFile = "default.bas";
    std::string saveFile = "default.prob";
    std::string restoreFile = "default.prob";
    std::string solutionFile = "stdout";
    std::string solutionSaveFile = "solution.file";
    int doIdiot = -1;
    int outputFormat = 2;
    int substitution = 3;
    int dualize = 3;
    int preSolve = 5;
    int doSprint = -1;
    int testOsiParameters = -1;
    int createSolver = 0;
    ClpSimplex * lpSolver;
    OsiClpSolverInterface * clpSolver;
    if (model_.solver()) {
        clpSolver = dynamic_cast<OsiClpSolverInterface *> (model_.solver());
        assert (clpSolver);
        lpSolver = clpSolver->getModelPtr();
        assert (lpSolver);
    } else {
        lpSolver = new ClpSimplex();
        clpSolver = new OsiClpSolverInterface(lpSolver, true);
        createSolver = 1 ;
    }
    parameters_[whichParam(CLP_PARAM_ACTION_BASISIN, numberParameters_, parameters_)].setStringValue(importBasisFile);
    parameters_[whichParam(CBC_PARAM_ACTION_PRIORITYIN, numberParameters_, parameters_)].setStringValue(importPriorityFile);
    parameters_[whichParam(CLP_PARAM_ACTION_BASISOUT, numberParameters_, parameters_)].setStringValue(exportBasisFile);
    parameters_[whichParam(CLP_PARAM_ACTION_DEBUG, numberParameters_, parameters_)].setStringValue(debugFile);
    parameters_[whichParam(CLP_PARAM_ACTION_PRINTMASK, numberParameters_, parameters_)].setStringValue(printMask);
    parameters_[whichParam(CLP_PARAM_ACTION_DIRECTORY, numberParameters_, parameters_)].setStringValue(directory);
    parameters_[whichParam(CLP_PARAM_ACTION_DIRSAMPLE, numberParameters_, parameters_)].setStringValue(dirSample);
    parameters_[whichParam(CLP_PARAM_ACTION_DIRNETLIB, numberParameters_, parameters_)].setStringValue(dirNetlib);
    parameters_[whichParam(CBC_PARAM_ACTION_DIRMIPLIB, numberParameters_, parameters_)].setStringValue(dirMiplib);
    parameters_[whichParam(CLP_PARAM_DBL_DUALBOUND, numberParameters_, parameters_)].setDoubleValue(lpSolver->dualBound());
    parameters_[whichParam(CLP_PARAM_DBL_DUALTOLERANCE, numberParameters_, parameters_)].setDoubleValue(lpSolver->dualTolerance());
    parameters_[whichParam(CLP_PARAM_ACTION_EXPORT, numberParameters_, parameters_)].setStringValue(exportFile);
    parameters_[whichParam(CLP_PARAM_INT_IDIOT, numberParameters_, parameters_)].setIntValue(doIdiot);
    parameters_[whichParam(CLP_PARAM_ACTION_IMPORT, numberParameters_, parameters_)].setStringValue(importFile);
    parameters_[whichParam(CLP_PARAM_DBL_PRESOLVETOLERANCE, numberParameters_, parameters_)].setDoubleValue(1.0e-8);
    int iParam = whichParam(CLP_PARAM_INT_SOLVERLOGLEVEL, numberParameters_, parameters_);
    int value = 1;
    clpSolver->messageHandler()->setLogLevel(1) ;
    lpSolver->setLogLevel(1);
    parameters_[iParam].setIntValue(value);
    iParam = whichParam(CLP_PARAM_INT_LOGLEVEL, numberParameters_, parameters_);
    model_.messageHandler()->setLogLevel(value);
    parameters_[iParam].setIntValue(value);
    parameters_[whichParam(CLP_PARAM_INT_MAXFACTOR, numberParameters_, parameters_)].setIntValue(lpSolver->factorizationFrequency());
    parameters_[whichParam(CLP_PARAM_INT_MAXITERATION, numberParameters_, parameters_)].setIntValue(lpSolver->maximumIterations());
    parameters_[whichParam(CLP_PARAM_INT_OUTPUTFORMAT, numberParameters_, parameters_)].setIntValue(outputFormat);
    parameters_[whichParam(CLP_PARAM_INT_PRESOLVEPASS, numberParameters_, parameters_)].setIntValue(preSolve);
    parameters_[whichParam(CLP_PARAM_INT_PERTVALUE, numberParameters_, parameters_)].setIntValue(lpSolver->perturbation());
    parameters_[whichParam(CLP_PARAM_DBL_PRIMALTOLERANCE, numberParameters_, parameters_)].setDoubleValue(lpSolver->primalTolerance());
    parameters_[whichParam(CLP_PARAM_DBL_PRIMALWEIGHT, numberParameters_, parameters_)].setDoubleValue(lpSolver->infeasibilityCost());
    parameters_[whichParam(CLP_PARAM_ACTION_RESTORE, numberParameters_, parameters_)].setStringValue(restoreFile);
    parameters_[whichParam(CLP_PARAM_ACTION_SAVE, numberParameters_, parameters_)].setStringValue(saveFile);
    //parameters_[whichParam(CLP_PARAM_DBL_TIMELIMIT,numberParameters_,parameters_)].setDoubleValue(1.0e8);
    parameters_[whichParam(CBC_PARAM_DBL_TIMELIMIT_BAB, numberParameters_, parameters_)].setDoubleValue(1.0e8);
    parameters_[whichParam(CLP_PARAM_ACTION_SOLUTION, numberParameters_, parameters_)].setStringValue(solutionFile);
    parameters_[whichParam(CLP_PARAM_ACTION_SAVESOL, numberParameters_, parameters_)].setStringValue(solutionSaveFile);
    parameters_[whichParam(CLP_PARAM_INT_SPRINT, numberParameters_, parameters_)].setIntValue(doSprint);
    parameters_[whichParam(CLP_PARAM_INT_SUBSTITUTION, numberParameters_, parameters_)].setIntValue(substitution);
    parameters_[whichParam(CLP_PARAM_INT_DUALIZE, numberParameters_, parameters_)].setIntValue(dualize);
    parameters_[whichParam(CBC_PARAM_INT_NUMBERBEFORE, numberParameters_, parameters_)].setIntValue(model_.numberBeforeTrust());
    parameters_[whichParam(CBC_PARAM_INT_MAXNODES, numberParameters_, parameters_)].setIntValue(model_.getMaximumNodes());
    parameters_[whichParam(CBC_PARAM_INT_STRONGBRANCHING, numberParameters_, parameters_)].setIntValue(model_.numberStrong());
    parameters_[whichParam(CBC_PARAM_DBL_INFEASIBILITYWEIGHT, numberParameters_, parameters_)].setDoubleValue(model_.getDblParam(CbcModel::CbcInfeasibilityWeight));
    parameters_[whichParam(CBC_PARAM_DBL_INTEGERTOLERANCE, numberParameters_, parameters_)].setDoubleValue(model_.getDblParam(CbcModel::CbcIntegerTolerance));
    parameters_[whichParam(CBC_PARAM_DBL_INCREMENT, numberParameters_, parameters_)].setDoubleValue(model_.getDblParam(CbcModel::CbcCutoffIncrement));
    parameters_[whichParam(CBC_PARAM_INT_TESTOSI, numberParameters_, parameters_)].setIntValue(testOsiParameters);
    parameters_[whichParam(CBC_PARAM_INT_FPUMPTUNE, numberParameters_, parameters_)].setIntValue(1003);
    initialPumpTune = 1003;
#ifdef CBC_THREAD
    parameters_[whichParam(CBC_PARAM_INT_THREADS, numberParameters_, parameters_)].setIntValue(0);
#endif
    // Set up likely cut generators and defaults
    parameters_[whichParam(CBC_PARAM_STR_PREPROCESS, numberParameters_, parameters_)].setCurrentOption("sos");
    parameters_[whichParam(CBC_PARAM_INT_MIPOPTIONS, numberParameters_, parameters_)].setIntValue(1057);
    parameters_[whichParam(CBC_PARAM_INT_CUTPASSINTREE, numberParameters_, parameters_)].setIntValue(1);
    parameters_[whichParam(CBC_PARAM_INT_MOREMIPOPTIONS, numberParameters_, parameters_)].setIntValue(-1);
    parameters_[whichParam(CBC_PARAM_INT_MAXHOTITS, numberParameters_, parameters_)].setIntValue(100);
    parameters_[whichParam(CBC_PARAM_STR_CUTSSTRATEGY, numberParameters_, parameters_)].setCurrentOption("on");
    parameters_[whichParam(CBC_PARAM_STR_HEURISTICSTRATEGY, numberParameters_, parameters_)].setCurrentOption("on");
    parameters_[whichParam(CBC_PARAM_STR_NODESTRATEGY, numberParameters_, parameters_)].setCurrentOption("fewest");
    parameters_[whichParam(CBC_PARAM_STR_GOMORYCUTS, numberParameters_, parameters_)].setCurrentOption("ifmove");
    parameters_[whichParam(CBC_PARAM_STR_PROBINGCUTS, numberParameters_, parameters_)].setCurrentOption("ifmove");
    parameters_[whichParam(CBC_PARAM_STR_KNAPSACKCUTS, numberParameters_, parameters_)].setCurrentOption("ifmove");
#ifdef ZERO_HALF_CUTS
    parameters_[whichParam(CBC_PARAM_STR_ZEROHALFCUTS, numberParameters_, parameters_)].setCurrentOption("off");
#endif
    parameters_[whichParam(CBC_PARAM_STR_REDSPLITCUTS, numberParameters_, parameters_)].setCurrentOption("off");
    parameters_[whichParam(CBC_PARAM_STR_CLIQUECUTS, numberParameters_, parameters_)].setCurrentOption("ifmove");
    parameters_[whichParam(CBC_PARAM_STR_MIXEDCUTS, numberParameters_, parameters_)].setCurrentOption("ifmove");
    parameters_[whichParam(CBC_PARAM_STR_FLOWCUTS, numberParameters_, parameters_)].setCurrentOption("ifmove");
    parameters_[whichParam(CBC_PARAM_STR_TWOMIRCUTS, numberParameters_, parameters_)].setCurrentOption("root");
    parameters_[whichParam(CBC_PARAM_STR_LANDPCUTS, numberParameters_, parameters_)].setCurrentOption("off");
    parameters_[whichParam(CBC_PARAM_STR_RESIDCUTS, numberParameters_, parameters_)].setCurrentOption("off");
    parameters_[whichParam(CBC_PARAM_STR_ROUNDING, numberParameters_, parameters_)].setCurrentOption("on");
    parameters_[whichParam(CBC_PARAM_STR_FPUMP, numberParameters_, parameters_)].setCurrentOption("on");
    parameters_[whichParam(CBC_PARAM_STR_GREEDY, numberParameters_, parameters_)].setCurrentOption("on");
    parameters_[whichParam(CBC_PARAM_STR_COMBINE, numberParameters_, parameters_)].setCurrentOption("on");
    parameters_[whichParam(CBC_PARAM_STR_CROSSOVER2, numberParameters_, parameters_)].setCurrentOption("off");
    parameters_[whichParam(CBC_PARAM_STR_PIVOTANDCOMPLEMENT, numberParameters_, parameters_)].setCurrentOption("off");
    parameters_[whichParam(CBC_PARAM_STR_PIVOTANDFIX, numberParameters_, parameters_)].setCurrentOption("off");
    parameters_[whichParam(CBC_PARAM_STR_RANDROUND, numberParameters_, parameters_)].setCurrentOption("off");
    parameters_[whichParam(CBC_PARAM_STR_NAIVE, numberParameters_, parameters_)].setCurrentOption("off");
    parameters_[whichParam(CBC_PARAM_STR_RINS, numberParameters_, parameters_)].setCurrentOption("off");
    parameters_[whichParam(CBC_PARAM_STR_DINS, numberParameters_, parameters_)].setCurrentOption("off");
    parameters_[whichParam(CBC_PARAM_STR_RENS, numberParameters_, parameters_)].setCurrentOption("off");
    parameters_[whichParam(CBC_PARAM_STR_LOCALTREE, numberParameters_, parameters_)].setCurrentOption("off");
    parameters_[whichParam(CBC_PARAM_STR_COSTSTRATEGY, numberParameters_, parameters_)].setCurrentOption("off");
    if (createSolver)
        delete clpSolver;
}

/*
  Initialise a subset of the parameters prior to processing any input from
  the user.

  Why this choice of subset?
*/
/*!
  \todo Guard/replace clp-specific code
*/
void CbcSolver::fillValuesInSolver()
{
    OsiSolverInterface * solver = model_.solver();
    OsiClpSolverInterface * clpSolver =
        dynamic_cast< OsiClpSolverInterface*> (solver);
    assert (clpSolver);
    ClpSimplex * lpSolver = clpSolver->getModelPtr();

    /*
      Why are we reaching into the underlying solver(s) for these settings?
      Shouldn't CbcSolver have its own defaults, which are then imposed on the
      underlying solver?

      Coming at if from the other side, if CbcSolver had the capability to use
      multiple solvers then it definitely makes sense to acquire the defaults from
      the solver (on the assumption that we haven't processed command line
      parameters yet, which can then override the defaults). But then it's more of
      a challenge to avoid solver-specific coding here.
    */
    noPrinting_ = (lpSolver->logLevel() == 0);
    CoinMessageHandler * generalMessageHandler = clpSolver->messageHandler();
    generalMessageHandler->setPrefix(true);

    lpSolver->setPerturbation(50);
    lpSolver->messageHandler()->setPrefix(false);

    parameters_[whichParam(CLP_PARAM_DBL_DUALBOUND, numberParameters_, parameters_)].setDoubleValue(lpSolver->dualBound());
    parameters_[whichParam(CLP_PARAM_DBL_DUALTOLERANCE, numberParameters_, parameters_)].setDoubleValue(lpSolver->dualTolerance());
    /*
      Why are we doing this? We read the log level from parameters_, set it into
      the message handlers for cbc and the underlying solver. Then we read the
      log level back from the handlers and use it to set the values in
      parameters_!
    */
    int iParam = whichParam(CLP_PARAM_INT_SOLVERLOGLEVEL, numberParameters_, parameters_);
    int value = parameters_[iParam].intValue();
    clpSolver->messageHandler()->setLogLevel(value) ;
    lpSolver->setLogLevel(value);
    iParam = whichParam(CLP_PARAM_INT_LOGLEVEL, numberParameters_, parameters_);
    value = parameters_[iParam].intValue();
    model_.messageHandler()->setLogLevel(value);
    parameters_[whichParam(CLP_PARAM_INT_LOGLEVEL, numberParameters_, parameters_)].setIntValue(model_.logLevel());
    parameters_[whichParam(CLP_PARAM_INT_SOLVERLOGLEVEL, numberParameters_, parameters_)].setIntValue(lpSolver->logLevel());
    parameters_[whichParam(CLP_PARAM_INT_MAXFACTOR, numberParameters_, parameters_)].setIntValue(lpSolver->factorizationFrequency());
    parameters_[whichParam(CLP_PARAM_INT_MAXITERATION, numberParameters_, parameters_)].setIntValue(lpSolver->maximumIterations());
    parameters_[whichParam(CLP_PARAM_INT_PERTVALUE, numberParameters_, parameters_)].setIntValue(lpSolver->perturbation());
    parameters_[whichParam(CLP_PARAM_DBL_PRIMALTOLERANCE, numberParameters_, parameters_)].setDoubleValue(lpSolver->primalTolerance());
    parameters_[whichParam(CLP_PARAM_DBL_PRIMALWEIGHT, numberParameters_, parameters_)].setDoubleValue(lpSolver->infeasibilityCost());
    parameters_[whichParam(CBC_PARAM_INT_NUMBERBEFORE, numberParameters_, parameters_)].setIntValue(model_.numberBeforeTrust());
    parameters_[whichParam(CBC_PARAM_INT_MAXNODES, numberParameters_, parameters_)].setIntValue(model_.getMaximumNodes());
    parameters_[whichParam(CBC_PARAM_INT_STRONGBRANCHING, numberParameters_, parameters_)].setIntValue(model_.numberStrong());
    parameters_[whichParam(CBC_PARAM_DBL_INFEASIBILITYWEIGHT, numberParameters_, parameters_)].setDoubleValue(model_.getDblParam(CbcModel::CbcInfeasibilityWeight));
    parameters_[whichParam(CBC_PARAM_DBL_INTEGERTOLERANCE, numberParameters_, parameters_)].setDoubleValue(model_.getDblParam(CbcModel::CbcIntegerTolerance));
    parameters_[whichParam(CBC_PARAM_DBL_INCREMENT, numberParameters_, parameters_)].setDoubleValue(model_.getDblParam(CbcModel::CbcCutoffIncrement));
}
// Add user function
void
CbcSolver::addUserFunction(CbcUser * function)
{
    CbcUser ** temp = new CbcUser * [numberUserFunctions_+1];
    int i;
    for (i = 0; i < numberUserFunctions_; i++)
        temp[i] = userFunction_[i];
    delete [] userFunction_;
    userFunction_ = temp;
    userFunction_[numberUserFunctions_++] = function->clone();
    delete [] statusUserFunction_;
    statusUserFunction_ = NULL;
}
// Set user call back
void
CbcSolver::setUserCallBack(CbcStopNow * function)
{
    delete callBack_;
    callBack_ = function->clone();
}
// Copy of model on initial load (will contain output solutions)
void
CbcSolver::setOriginalSolver(OsiClpSolverInterface * originalSolver)
{
    delete originalSolver_;
    OsiSolverInterface * temp = originalSolver->clone();
    originalSolver_ = dynamic_cast<OsiClpSolverInterface *> (temp);
    assert (originalSolver_);

}
// Copy of model on initial load
void
CbcSolver::setOriginalCoinModel(CoinModel * originalCoinModel)
{
    delete originalCoinModel_;
    originalCoinModel_ = new CoinModel(*originalCoinModel);
}
// Add cut generator
void
CbcSolver::addCutGenerator(CglCutGenerator * generator)
{
    CglCutGenerator ** temp = new CglCutGenerator * [numberCutGenerators_+1];
    int i;
    for (i = 0; i < numberCutGenerators_; i++)
        temp[i] = cutGenerator_[i];
    delete [] cutGenerator_;
    cutGenerator_ = temp;
    cutGenerator_[numberCutGenerators_++] = generator->clone();
}

/*
  The only other solver that's ever been used is cplex, and the use is
  limited -- do the root with clp and all the cbc smarts, then give the
  problem over to cplex to finish. Although the defines can be read in some
  places to allow other options, nothing's been tested and success is
  unlikely.

  CBC_OTHER_SOLVER == 1 is cplex.
*/

#if CBC_OTHER_SOLVER==1
#  ifndef COIN_HAS_CPX
#    error "Configuration did not detect cplex installation."
#  else
#    include "OsiCpxSolverInterface.hpp"
#  endif
#endif

#ifdef COIN_HAS_ASL
#include "Cbc_ampl.h"
#endif

static double totalTime = 0.0;
static void statistics(ClpSimplex * originalModel, ClpSimplex * model);
static bool maskMatches(const int * starts, char ** masks,
                        std::string & check);
static void generateCode(CbcModel * model, const char * fileName, int type, int preProcess);

// dummy fake main programs for UserClp and UserCbc
void fakeMain (ClpSimplex & model, OsiSolverInterface & osiSolver, CbcModel & babSolver);
void fakeMain2 (ClpSimplex & model, OsiClpSolverInterface & osiSolver, int options);

// Allow for interrupts
// But is this threadsafe? (so switched off by option)

#include "CoinSignal.hpp"
static CbcModel * currentBranchModel = NULL;

extern "C" {
    static void signal_handler(int /*whichSignal*/) {
        if (currentBranchModel != NULL) {
            currentBranchModel->setMaximumNodes(0); // stop at next node
            currentBranchModel->setMaximumSeconds(0.0); // stop
        }
        return;
    }
}

//#define CBC_SIG_TRAP
#ifdef CBC_SIG_TRAP
#include <setjmp.h>
static sigjmp_buf cbc_seg_buffer;
extern "C" {
    static void signal_handler_error(int whichSignal) {
        siglongjmp(cbc_seg_buffer, 1);
    }
}
#endif


/*
  Debug checks on special ordered sets.

  This is active only for debugging. The entire body of the routine becomes
  a noop when COIN_DEVELOP is not defined. To avoid compiler warnings, the
  formal parameters also need to go away.
*/
#ifdef COIN_DEVELOP
void checkSOS(CbcModel * babModel, const OsiSolverInterface * solver)
#else
void checkSOS(CbcModel * /*babModel*/, const OsiSolverInterface * /*solver*/)
#endif
{
#ifdef COIN_DEVELOP
    if (!babModel->ownObjects())
        return;
#if COIN_DEVELOP>2
    //const double *objective = solver->getObjCoefficients() ;
    const double *columnLower = solver->getColLower() ;
    const double * columnUpper = solver->getColUpper() ;
    const double * solution = solver->getColSolution();
    //int numberRows = solver->getNumRows();
    //double direction = solver->getObjSense();
    //int iRow,iColumn;
#endif

    // Row copy
    CoinPackedMatrix matrixByRow(*solver->getMatrixByRow());
    //const double * elementByRow = matrixByRow.getElements();
    //const int * column = matrixByRow.getIndices();
    //const CoinBigIndex * rowStart = matrixByRow.getVectorStarts();
    const int * rowLength = matrixByRow.getVectorLengths();

    // Column copy
    CoinPackedMatrix  matrixByCol(*solver->getMatrixByCol());
    const double * element = matrixByCol.getElements();
    const int * row = matrixByCol.getIndices();
    const CoinBigIndex * columnStart = matrixByCol.getVectorStarts();
    const int * columnLength = matrixByCol.getVectorLengths();

    const double * rowLower = solver->getRowLower();
    const double * rowUpper = solver->getRowUpper();
    OsiObject ** objects = babModel->objects();
    int numberObjects = babModel->numberObjects();
    int numberColumns = solver->getNumCols() ;
    for (int iObj = 0; iObj < numberObjects; iObj++) {
        CbcSOS * objSOS =
            dynamic_cast <CbcSOS *>(objects[iObj]) ;
        if (objSOS) {
            int n = objSOS->numberMembers();
            const int * which = objSOS->members();
#if COIN_DEVELOP>2
            const double * weight = objSOS->weights();
#endif
            int type = objSOS->sosType();
            // convexity row?
            int iColumn;
            iColumn = which[0];
            int j;
            int convex = -1;
            for (j = columnStart[iColumn]; j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                int iRow = row[j];
                double value = element[j];
                if (rowLower[iRow] == 1.0 && rowUpper[iRow] == 1.0 &&
                        value == 1.0) {
                    // possible
                    if (rowLength[iRow] == n) {
                        if (convex == -1)
                            convex = iRow;
                        else
                            convex = -2;
                    }
                }
            }
            printf ("set %d of type %d has %d members - possible convexity row %d\n",
                    iObj, type, n, convex);
            for (int i = 0; i < n; i++) {
                iColumn = which[i];
                // Column may have been added
                if (iColumn < numberColumns) {
                    int convex2 = -1;
                    for (j = columnStart[iColumn]; j < columnStart[iColumn] + columnLength[iColumn]; j++) {
                        int iRow = row[j];
                        if (iRow == convex) {
                            double value = element[j];
                            if (value == 1.0) {
                                convex2 = iRow;
                            }
                        }
                    }
                    if (convex2<0 && convex >= 0) {
                        printf("odd convexity row\n");
                        convex = -2;
                    }
#if COIN_DEVELOP>2
                    printf("col %d has weight %g and value %g, bounds %g %g\n",
                           iColumn, weight[i], solution[iColumn], columnLower[iColumn],
                           columnUpper[iColumn]);
#endif
                }
            }
        }
    }
#endif	// COIN_DEVELOP
}

static int dummyCallBack(CbcModel * /*model*/, int /*whereFrom*/)
{
    return 0;
}


/*
  Global parameters for command processing.

  These will need to be moved into an object of some sort in order to make
  this set of calls thread-safe.
*/

int CbcOrClpRead_mode = 1;
FILE * CbcOrClpReadCommand = stdin;
extern int CbcOrClpEnvironmentIndex;
static bool noPrinting = false;



/*
  Wrappers for CbcMain0, CbcMain1. The various forms of callCbc will eventually
  resolve to a call to CbcMain0 followed by a call to callCbc1.
*/
/*
  Simplest calling form: supply just a string with the command options. The
  wrapper creates an OsiClpSolverInterface and calls the next wrapper.
*/
int callCbc(const std::string input2)
{
    char * input3 = CoinStrdup(input2.c_str());
    OsiClpSolverInterface solver1;
    int returnCode = callCbc(input3, solver1);
    free(input3);
    return returnCode;
}

int callCbc(const char * input2)
{
    {
        OsiClpSolverInterface solver1;
        return callCbc(input2, solver1);
    }
}

/*
  Second calling form: supply the command line and an OsiClpSolverInterface.
  the wrapper will create a CbcModel and call the next wrapper.
*/

int callCbc(const std::string input2, OsiClpSolverInterface& solver1)
{
    char * input3 = CoinStrdup(input2.c_str());
    int returnCode = callCbc(input3, solver1);
    free(input3);
    return returnCode;
}

int callCbc(const char * input2, OsiClpSolverInterface& solver1)
{
    CbcModel model(solver1);
    return callCbc(input2, model);
}

/*
  Third calling form: supply the command line and a CbcModel. This wrapper will
  actually call CbcMain0 and then call the next set of wrappers (callCbc1) to
  handle the call to CbcMain1.
*/
int callCbc(const char * input2, CbcModel & babSolver)
{
    CbcMain0(babSolver);
    return callCbc1(input2, babSolver);
}

int callCbc(const std::string input2, CbcModel & babSolver)
{
    char * input3 = CoinStrdup(input2.c_str());
    CbcMain0(babSolver);
    int returnCode = callCbc1(input3, babSolver);
    free(input3);
    return returnCode;
}


/*
  Various overloads of callCbc1. The first pair accepts just a CbcModel and
  supplements it with a dummy callback routine. The second pair allows the
  user to supply a callback. See CbcMain1 for further explanation of the
  callback. The various overloads of callCbc1 resolve to the final version,
  which breaks the string into individual parameter strings (i.e., creates
  something that looks like a standard argv vector).
*/

int callCbc1(const std::string input2, CbcModel & babSolver)
{
    char * input3 = CoinStrdup(input2.c_str());
    int returnCode = callCbc1(input3, babSolver);
    free(input3);
    return returnCode;
}

int callCbc1(const char * input2, CbcModel & model)
{
    return callCbc1(input2, model, dummyCallBack);
}

int callCbc1(const std::string input2, CbcModel & babSolver,
             int callBack(CbcModel * currentSolver, int whereFrom))
{
    char * input3 = CoinStrdup(input2.c_str());
    int returnCode = callCbc1(input3, babSolver, callBack);
    free(input3);
    return returnCode;
}

int callCbc1(const char * input2, CbcModel & model,
             int callBack(CbcModel * currentSolver, int whereFrom))
{
    char * input = CoinStrdup(input2);
    size_t length = strlen(input);
    bool blank = input[0] == '0';
    int n = blank ? 0 : 1;
    for (size_t i = 0; i < length; i++) {
        if (blank) {
            // look for next non blank
            if (input[i] == ' ') {
                continue;
            } else {
                n++;
                blank = false;
            }
        } else {
            // look for next blank
            if (input[i] != ' ') {
                continue;
            } else {
                blank = true;
            }
        }
    }
    char ** argv = new char * [n+2];
    argv[0] = CoinStrdup("cbc");
    size_t i = 0;
    while (input[i] == ' ')
        i++;
    for (int j = 0; j < n; j++) {
        size_t saveI = i;
        for (; i < length; i++) {
            // look for next blank
            if (input[i] != ' ') {
                continue;
            } else {
                break;
            }
        }
        input[i] = '\0';
        argv[j+1] = CoinStrdup(input + saveI);
        while (input[i] == ' ')
            i++;
    }
    argv[n+1] = CoinStrdup("-quit");
    free(input);
    totalTime = 0.0;
    currentBranchModel = NULL;
    CbcOrClpRead_mode = 1;
    CbcOrClpReadCommand = stdin;
    noPrinting = false;
    int returnCode = CbcMain1(n + 2, const_cast<const char **>(argv),
                              model, callBack);
    for (int k = 0; k < n + 2; k++)
        free(argv[k]);
    delete [] argv;
    return returnCode;
}

#define CBCMAXPARAMETERS 200
static CbcOrClpParam parameters[CBCMAXPARAMETERS];
static int numberParameters = 0 ;
CglPreProcess * cbcPreProcessPointer=NULL;

int CbcClpUnitTest (const CbcModel & saveModel,
                    const std::string& dirMiplib, int testSwitch,
                    const double * stuff);

int CbcMain1 (int argc, const char *argv[],
              CbcModel  & model)
{
    return CbcMain1(argc, argv, model, dummyCallBack);
}


/*
  Meaning of whereFrom:
    1 after initial solve by dualsimplex etc
    2 after preprocessing
    3 just before branchAndBound (so user can override)
    4 just after branchAndBound (before postprocessing)
    5 after postprocessing
    6 after a user called heuristic phase
*/

int CbcMain1 (int argc, const char *argv[],
              CbcModel  & model,
              int callBack(CbcModel * currentSolver, int whereFrom))
{
    CbcOrClpParam * parameters_ = parameters;
    int numberParameters_ = numberParameters;
    CbcModel & model_ = model;
#ifdef CBC_USE_INITIAL_TIME
    if (model_.useElapsedTime())
      model_.setDblParam(CbcModel::CbcStartSeconds, CoinGetTimeOfDay());
    else
      model_.setDblParam(CbcModel::CbcStartSeconds, CoinCpuTime());
#endif
    CbcModel * babModel_ = NULL;
    int returnMode = 1;
    CbcOrClpRead_mode = 1;
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
    int * statistics_number_cuts = NULL;
    const char ** statistics_name_generators = NULL;
    int statistics_number_generators = 0;
    memset(statusUserFunction_, 0, numberUserFunctions_*sizeof(int));
    /* Note
       This is meant as a stand-alone executable to do as much of coin as possible.
       It should only have one solver known to it.
    */
    CoinMessageHandler * generalMessageHandler = model_.messageHandler();
    generalMessageHandler->setPrefix(false);
#ifndef CBC_OTHER_SOLVER
    OsiClpSolverInterface * originalSolver = dynamic_cast<OsiClpSolverInterface *> (model_.solver());
    assert (originalSolver);
    // Move handler across if not default
    if (!originalSolver->defaultHandler() && originalSolver->getModelPtr()->defaultHandler())
        originalSolver->getModelPtr()->passInMessageHandler(originalSolver->messageHandler());
    CoinMessages generalMessages = originalSolver->getModelPtr()->messages();
    char generalPrint[10000];
    if (originalSolver->getModelPtr()->logLevel() == 0)
        noPrinting = true;
#elif CBC_OTHER_SOLVER==1
    OsiCpxSolverInterface * originalSolver = dynamic_cast<OsiCpxSolverInterface *> (model_.solver());
    assert (originalSolver);
    OsiClpSolverInterface dummySolver;
    OsiCpxSolverInterface * clpSolver = originalSolver;
    CoinMessages generalMessages = dummySolver.getModelPtr()->messages();
    char generalPrint[10000];
    noPrinting = true;
#endif
    bool noPrinting_ = noPrinting;
    // Say not in integer
    int integerStatus = -1;
    // Say no resolve after cuts
    model_.setResolveAfterTakeOffCuts(false);
    // see if log in list
    for (int i = 1; i < argc; i++) {
        if (!strncmp(argv[i], "log", 3)) {
            const char * equals = strchr(argv[i], '=');
            if (equals && atoi(equals + 1) != 0)
                noPrinting_ = false;
            else
                noPrinting_ = true;
            break;
        } else if (!strncmp(argv[i], "-log", 4) && i < argc - 1) {
            if (atoi(argv[i+1]) != 0)
                noPrinting_ = false;
            else
                noPrinting_ = true;
            break;
        }
    }
    double time0;
    double time0Elapsed = CoinGetTimeOfDay();
    {
        double time1 = CoinCpuTime(), time2;
        time0 = time1;
	double time1Elapsed = time0Elapsed;
        bool goodModel = (originalSolver->getNumCols()) ? true : false;

        // register signal handler
        //CoinSighandler_t saveSignal=signal(SIGINT,signal_handler);
#ifndef CBC_QUIET
        signal(SIGINT, signal_handler);
#endif
        // Set up all non-standard stuff
        int cutPass = -1234567;
        int cutPassInTree = -1234567;
        int tunePreProcess = 0;
        int testOsiParameters = -1;
        // 0 normal, 1 from ampl or MIQP etc (2 allows cuts)
        int complicatedInteger = 0;
        OsiSolverInterface * solver = model_.solver();
        if (noPrinting_)
            setCbcOrClpPrinting(false);
#ifndef CBC_OTHER_SOLVER
        OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
        ClpSimplex * lpSolver = clpSolver->getModelPtr();
        if (noPrinting_) {
            lpSolver->setLogLevel(0);
        }
#else
        ClpSimplex * lpSolver = NULL;
#endif
        // For priorities etc
        int * priorities = NULL;
        int * branchDirection = NULL;
        double * pseudoDown = NULL;
        double * pseudoUp = NULL;
        double * solutionIn = NULL;
        int * prioritiesIn = NULL;
        int numberSOS = 0;
        int * sosStart = NULL;
        int * sosIndices = NULL;
        char * sosType = NULL;
        double * sosReference = NULL;
        int * cut = NULL;
        int * sosPriority = NULL;
        CglStored storedAmpl;
        CoinModel * coinModel = NULL;
        CoinModel saveCoinModel;
        CoinModel saveTightenedModel;
        int * whichColumn = NULL;
        int * knapsackStart = NULL;
        int * knapsackRow = NULL;
        int numberKnapsack = 0;
#ifdef COIN_HAS_ASL
        ampl_info info;
        {
            memset(&info, 0, sizeof(info));
            if (argc > 2 && !strcmp(argv[2], "-AMPL")) {
                statusUserFunction_[0] = 1;
                // see if log in list
                noPrinting_ = true;
                for (int i = 1; i < argc; i++) {
                    if (!strncmp(argv[i], "log", 3)) {
                        const char * equals = strchr(argv[i], '=');
                        if (equals && atoi(equals + 1) > 0) {
                            noPrinting_ = false;
                            info.logLevel = atoi(equals + 1);
                            int log = whichParam(CLP_PARAM_INT_LOGLEVEL, numberParameters_, parameters_);
                            parameters_[log].setIntValue(info.logLevel);
                            // mark so won't be overWritten
                            info.numberRows = -1234567;
                            break;
                        }
                    }
                }

                union {
                    void * voidModel;
                    CoinModel * model;
                } coinModelStart;
                coinModelStart.model = NULL;
                int returnCode = readAmpl(&info, argc, const_cast<char **>(argv), & coinModelStart.voidModel);
                coinModel = coinModelStart.model;
                if (returnCode)
                    return returnCode;
                CbcOrClpRead_mode = 2; // so will start with parameters
                // see if log in list (including environment)
                for (int i = 1; i < info.numberArguments; i++) {
                    if (!strcmp(info.arguments[i], "log")) {
                        if (i < info.numberArguments - 1 && atoi(info.arguments[i+1]) > 0)
                            noPrinting_ = false;
                        break;
                    }
                }
                if (noPrinting_) {
                    model_.messageHandler()->setLogLevel(0);
                    setCbcOrClpPrinting(false);
                }
                if (!noPrinting_)
                    printf("%d rows, %d columns and %d elements\n",
                           info.numberRows, info.numberColumns, info.numberElements);
#ifdef COIN_HAS_LINK
                if (!coinModel) {
#endif
                    solver->loadProblem(info.numberColumns, info.numberRows, info.starts,
                                        info.rows, info.elements,
                                        info.columnLower, info.columnUpper, info.objective,
                                        info.rowLower, info.rowUpper);
                    // take off cuts if ampl wants that
                    if (info.cut && 0) {
                        printf("AMPL CUTS OFF until global cuts fixed\n");
                        info.cut = NULL;
                    }
                    if (info.cut) {
                        int numberRows = info.numberRows;
                        int * whichRow = new int [numberRows];
                        // Row copy
                        const CoinPackedMatrix * matrixByRow = solver->getMatrixByRow();
                        const double * elementByRow = matrixByRow->getElements();
                        const int * column = matrixByRow->getIndices();
                        const CoinBigIndex * rowStart = matrixByRow->getVectorStarts();
                        const int * rowLength = matrixByRow->getVectorLengths();

                        const double * rowLower = solver->getRowLower();
                        const double * rowUpper = solver->getRowUpper();
                        int nDelete = 0;
                        for (int iRow = 0; iRow < numberRows; iRow++) {
                            if (info.cut[iRow]) {
                                whichRow[nDelete++] = iRow;
                                int start = rowStart[iRow];
                                storedAmpl.addCut(rowLower[iRow], rowUpper[iRow],
                                                  rowLength[iRow], column + start, elementByRow + start);
                            }
                        }
                        solver->deleteRows(nDelete, whichRow);
                        delete [] whichRow;
                    }
#ifdef COIN_HAS_LINK
                } else {
#ifndef CBC_OTHER_SOLVER
                    // save
                    saveCoinModel = *coinModel;
                    // load from coin model
                    OsiSolverLink solver1;
                    OsiSolverInterface * solver2 = solver1.clone();
                    model_.assignSolver(solver2, false);
                    OsiSolverLink * si =
                        dynamic_cast<OsiSolverLink *>(model_.solver()) ;
                    assert (si != NULL);
                    si->setDefaultMeshSize(0.001);
                    // need some relative granularity
                    si->setDefaultBound(100.0);
                    double dextra3 = parameters_[whichParam(CBC_PARAM_DBL_DEXTRA3, numberParameters_, parameters_)].doubleValue();
                    if (dextra3)
                        si->setDefaultMeshSize(dextra3);
                    si->setDefaultBound(100000.0);
                    si->setIntegerPriority(1000);
                    si->setBiLinearPriority(10000);
                    CoinModel * model2 = reinterpret_cast<CoinModel *> (coinModel);
                    int logLevel = parameters_[whichParam(CLP_PARAM_INT_LOGLEVEL, numberParameters_, parameters_)].intValue();
                    si->load(*model2, true, logLevel);
                    // redo
                    solver = model_.solver();
                    clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
                    lpSolver = clpSolver->getModelPtr();
                    clpSolver->messageHandler()->setLogLevel(0) ;
                    testOsiParameters = 0;
                    parameters_[whichParam(CBC_PARAM_INT_TESTOSI, numberParameters_, parameters_)].setIntValue(0);
                    complicatedInteger = 1;
                    if (info.cut) {
                        printf("Sorry - can't do cuts with LOS as ruins delicate row order\n");
                        abort();
                        int numberRows = info.numberRows;
                        int * whichRow = new int [numberRows];
                        // Row copy
                        const CoinPackedMatrix * matrixByRow = solver->getMatrixByRow();
                        const double * elementByRow = matrixByRow->getElements();
                        const int * column = matrixByRow->getIndices();
                        const CoinBigIndex * rowStart = matrixByRow->getVectorStarts();
                        const int * rowLength = matrixByRow->getVectorLengths();

                        const double * rowLower = solver->getRowLower();
                        const double * rowUpper = solver->getRowUpper();
                        int nDelete = 0;
                        for (int iRow = 0; iRow < numberRows; iRow++) {
                            if (info.cut[iRow]) {
                                whichRow[nDelete++] = iRow;
                                int start = rowStart[iRow];
                                storedAmpl.addCut(rowLower[iRow], rowUpper[iRow],
                                                  rowLength[iRow], column + start, elementByRow + start);
                            }
                        }
                        solver->deleteRows(nDelete, whichRow);
                        // and special matrix
                        si->cleanMatrix()->deleteRows(nDelete, whichRow);
                        delete [] whichRow;
                    }
#endif
                }
#endif
                // If we had a solution use it
                if (info.primalSolution) {
                    solver->setColSolution(info.primalSolution);
                }
                // status
                if (info.rowStatus) {
                    unsigned char * statusArray = lpSolver->statusArray();
                    int i;
                    for (i = 0; i < info.numberColumns; i++)
                        statusArray[i] = static_cast<unsigned char>(info.columnStatus[i]);
                    statusArray += info.numberColumns;
                    for (i = 0; i < info.numberRows; i++)
                        statusArray[i] = static_cast<unsigned char>(info.rowStatus[i]);
                    CoinWarmStartBasis * basis = lpSolver->getBasis();
                    solver->setWarmStart(basis);
                    delete basis;
                }
                freeArrays1(&info);
                // modify objective if necessary
                solver->setObjSense(info.direction);
                solver->setDblParam(OsiObjOffset, info.offset);
                if (info.offset) {
                    sprintf(generalPrint, "Ampl objective offset is %g",
                            info.offset);
                    generalMessageHandler->message(CLP_GENERAL, generalMessages)
                    << generalPrint
                    << CoinMessageEol;
                }
                // Set integer variables (unless nonlinear when set)
                if (!info.nonLinear) {
                    for (int i = info.numberColumns - info.numberIntegers;
                            i < info.numberColumns; i++)
                        solver->setInteger(i);
                }
                goodModel = true;
                // change argc etc
                argc = info.numberArguments;
                argv = const_cast<const char **>(info.arguments);
            }
        }
#endif
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
        // set reasonable defaults
        int preSolve = 5;
        int preProcess = 4;
        bool useStrategy = false;
        bool preSolveFile = false;
        bool strongChanged = false;
        bool pumpChanged = false;

        double djFix = 1.0e100;
        double tightenFactor = 0.0;
        const char dirsep =  CoinFindDirSeparator();
        std::string directory;
        std::string dirSample;
        std::string dirNetlib;
        std::string dirMiplib;
        if (dirsep == '/') {
            directory = "./";
            dirSample = "../../Data/Sample/";
            dirNetlib = "../../Data/Netlib/";
            dirMiplib = "../../Data/miplib3/";
        } else {
            directory = ".\\";
            dirSample = "..\\..\\..\\..\\Data\\Sample\\";
            dirNetlib = "..\\..\\..\\..\\Data\\Netlib\\";
            dirMiplib = "..\\..\\..\\..\\Data\\miplib3\\";
        }
        std::string defaultDirectory = directory;
        std::string importFile = "";
        std::string exportFile = "default.mps";
        std::string importBasisFile = "";
        std::string importPriorityFile = "";
        std::string debugFile = "";
        std::string printMask = "";
        double * debugValues = NULL;
        int numberDebugValues = -1;
        int basisHasValues = 0;
        std::string exportBasisFile = "default.bas";
        std::string saveFile = "default.prob";
        std::string restoreFile = "default.prob";
        std::string solutionFile = "stdout";
        std::string solutionSaveFile = "solution.file";
        int slog = whichParam(CLP_PARAM_INT_SOLVERLOGLEVEL, numberParameters_, parameters_);
        int log = whichParam(CLP_PARAM_INT_LOGLEVEL, numberParameters_, parameters_);
#ifndef CBC_OTHER_SOLVER
        double normalIncrement = model_.getCutoffIncrement();;
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
        // set default action (0=off,1=on,2=root)
        int gomoryAction = 3;

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
        // set default action (0=off,1=on,2=root)
        int probingAction = 1;

        CglKnapsackCover knapsackGen;
        //knapsackGen.switchOnExpensive();
        //knapsackGen.setMaxInKnapsack(100);
        // set default action (0=off,1=on,2=root)
        int knapsackAction = 3;

        CglRedSplit redsplitGen;
        //redsplitGen.setLimit(100);
        // set default action (0=off,1=on,2=root)
        // Off as seems to give some bad cuts
        int redsplitAction = 0;

        CglFakeClique cliqueGen(NULL, false);
        //CglClique cliqueGen(false,true);
        cliqueGen.setStarCliqueReport(false);
        cliqueGen.setRowCliqueReport(false);
        cliqueGen.setMinViolation(0.1);
        // set default action (0=off,1=on,2=root)
        int cliqueAction = 3;

        // maxaggr,multiply,criterion(1-3)
        CglMixedIntegerRounding2 mixedGen(1, true, 1);
        // set default action (0=off,1=on,2=root)
        int mixedAction = 3;

        CglFlowCover flowGen;
        // set default action (0=off,1=on,2=root)
        int flowAction = 3;

        CglTwomir twomirGen;
        twomirGen.setMaxElements(250);
        // set default action (0=off,1=on,2=root)
        int twomirAction = 2;
#ifndef DEBUG_MALLOC
        CglLandP landpGen;
#endif
        // set default action (0=off,1=on,2=root)
        int landpAction = 0;
        CglResidualCapacity residualCapacityGen;
        // set default action (0=off,1=on,2=root)
        int residualCapacityAction = 0;

#ifdef ZERO_HALF_CUTS
        CglZeroHalf zerohalfGen;
        //zerohalfGen.switchOnExpensive();
        // set default action (0=off,1=on,2=root)
        int zerohalfAction = 0;
#endif

        // Stored cuts
        //bool storedCuts = false;

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
            //twomirGen.setMirScale(1,1);
            //twomirGen.setTwomirScale(1,1);
            //twomirGen.setAMax(2);
#else
            gomoryGen.setAwayAtRoot(0.01);
            twomirGen.setAwayAtRoot(0.01);
            twomirGen.setAway(0.01);
#endif
            int iParam;
            iParam = whichParam(CBC_PARAM_INT_DIVEOPT, numberParameters_, parameters_);
            parameters_[iParam].setIntValue(3);
            iParam = whichParam(CBC_PARAM_INT_FPUMPITS, numberParameters_, parameters_);
            parameters_[iParam].setIntValue(30);
            iParam = whichParam(CBC_PARAM_INT_FPUMPTUNE, numberParameters_, parameters_);
            parameters_[iParam].setIntValue(1005043);
            initialPumpTune = 1005043;
            iParam = whichParam(CLP_PARAM_INT_PROCESSTUNE, numberParameters_, parameters_);
            parameters_[iParam].setIntValue(6);
            tunePreProcess = 6;
            iParam = whichParam(CBC_PARAM_STR_DIVINGC, numberParameters_, parameters_);
            parameters_[iParam].setCurrentOption("on");
            iParam = whichParam(CBC_PARAM_STR_RINS, numberParameters_, parameters_);
            parameters_[iParam].setCurrentOption("on");
            iParam = whichParam(CBC_PARAM_STR_PROBINGCUTS, numberParameters_, parameters_);
            parameters_[iParam].setCurrentOption("on");
            probingAction = 1;
            parameters_[iParam].setCurrentOption("forceOnStrong");
            probingAction = 8;
        }
        std::string field;
#ifndef CBC_QUIET 
        if (!noPrinting_) {
	   sprintf(generalPrint,
		   "Welcome to the CBC MILP Solver \n");
	    if (strcmp(CBC_VERSION, "trunk")){
	       sprintf(generalPrint + strlen(generalPrint),
		       "Version: %s \n", CBC_VERSION);
	    }else{
	       sprintf(generalPrint + strlen(generalPrint),
		       "Version: Trunk (unstable) \n");
	    }
	    sprintf(generalPrint + strlen(generalPrint),
		    "Build Date: %s \n", __DATE__);
#ifdef CBC_SVN_REV
	    sprintf(generalPrint + strlen(generalPrint),
		    "Revision Number: %d \n", CBC_SVN_REV);
#endif
            generalMessageHandler->message(CLP_GENERAL, generalMessages)
            << generalPrint
            << CoinMessageEol;
            // Print command line
            if (argc > 1) {
                bool foundStrategy = false;
                sprintf(generalPrint, "command line - ");
                for (int i = 0; i < argc; i++) {
                    if (!argv[i])
                        break;
                    if (strstr(argv[i], "strat"))
                        foundStrategy = true;
                    sprintf(generalPrint + strlen(generalPrint), "%s ", argv[i]);
                }
                if (!foundStrategy)
                    sprintf(generalPrint + strlen(generalPrint), "(default strategy 1)");
                generalMessageHandler->message(CLP_GENERAL, generalMessages)
                << generalPrint
                << CoinMessageEol;
            }
        }
#endif
        while (1) {
            // next command
            field = CoinReadGetCommand(argc, argv);
            // Reset time
            time1 = CoinCpuTime();
	    time1Elapsed = CoinGetTimeOfDay();
            // adjust field if has odd trailing characters
            char temp [200];
            strcpy(temp, field.c_str());
            int length = static_cast<int>(strlen(temp));
            for (int k = length - 1; k >= 0; k--) {
                if (temp[k] < ' ')
                    length--;
                else
                    break;
            }
            temp[length] = '\0';
            field = temp;
            // exit if null or similar
            if (!field.length()) {
                if (numberGoodCommands == 1 && goodModel) {
                    // we just had file name - do branch and bound
                    field = "branch";
                } else if (!numberGoodCommands) {
                    // let's give the sucker a hint
                    std::cout
                        << "CoinSolver takes input from arguments ( - switches to stdin)"
                        << std::endl
                        << "Enter ? for list of commands or help" << std::endl;
                    field = "-";
                } else {
                    break;
                }
            }

            // see if ? at end
            size_t numberQuery = 0;
            if (field != "?" && field != "???") {
                size_t length = field.length();
                size_t i;
                for (i = length - 1; i > 0; i--) {
                    if (field[i] == '?')
                        numberQuery++;
                    else
                        break;
                }
                field = field.substr(0, length - numberQuery);
            }
            // find out if valid command
            int iParam;
            int numberMatches = 0;
            int firstMatch = -1;
            for ( iParam = 0; iParam < numberParameters_; iParam++ ) {
                int match = parameters_[iParam].matches(field);
                if (match == 1) {
                    numberMatches = 1;
                    firstMatch = iParam;
                    break;
                } else {
                    if (match && firstMatch < 0)
                        firstMatch = iParam;
                    numberMatches += match >> 1;
                }
            }
            if (iParam < numberParameters_ && !numberQuery) {
                // found
                CbcOrClpParam found = parameters_[iParam];
                CbcOrClpParameterType type = found.type();
                int valid;
                numberGoodCommands++;
                if (type == CBC_PARAM_ACTION_BAB && goodModel) {
#ifndef CBC_USE_INITIAL_TIME
		  if (model_.useElapsedTime())
		    model_.setDblParam(CbcModel::CbcStartSeconds, CoinGetTimeOfDay());
		  else
		    model_.setDblParam(CbcModel::CbcStartSeconds, CoinCpuTime());
#endif
                    // check if any integers
#ifndef CBC_OTHER_SOLVER
#ifdef COIN_HAS_ASL
                    if (info.numberSos && doSOS && statusUserFunction_[0]) {
                        // SOS
                        numberSOS = info.numberSos;
                    }
#endif
                    lpSolver = clpSolver->getModelPtr();
                    if (!lpSolver->integerInformation() && !numberSOS &&
                            !clpSolver->numberSOS() && !model_.numberObjects() && !clpSolver->numberObjects())
                        type = CLP_PARAM_ACTION_DUALSIMPLEX;
#endif
                }
                if (type == CBC_PARAM_GENERALQUERY) {
                    bool evenHidden = false;
                    int printLevel =
                        parameters_[whichParam(CLP_PARAM_STR_ALLCOMMANDS,
                                               numberParameters_, parameters_)].currentOptionAsInteger();
                    int convertP[] = {2, 1, 0};
                    printLevel = convertP[printLevel];
                    if ((verbose&8) != 0) {
                        // even hidden
                        evenHidden = true;
                        verbose &= ~8;
                    }
#ifdef COIN_HAS_ASL
                    if (verbose < 4 && statusUserFunction_[0])
                        verbose += 4;
#endif
                    if (verbose < 4) {
                        std::cout << "In argument list keywords have leading - "
                                  ", -stdin or just - switches to stdin" << std::endl;
                        std::cout << "One command per line (and no -)" << std::endl;
                        std::cout << "abcd? gives list of possibilities, if only one + explanation" << std::endl;
                        std::cout << "abcd?? adds explanation, if only one fuller help" << std::endl;
                        std::cout << "abcd without value (where expected) gives current value" << std::endl;
                        std::cout << "abcd value sets value" << std::endl;
                        std::cout << "Commands are:" << std::endl;
                    } else {
                        std::cout << "Cbc options are set within AMPL with commands like:" << std::endl << std::endl;
                        std::cout << "         option cbc_options \"cuts=root log=2 feas=on slog=1\"" << std::endl << std::endl;
                        std::cout << "only maximize, dual, primal, help and quit are recognized without =" << std::endl;
                    }
                    int maxAcross = 10;
                    if ((verbose % 4) != 0)
                        maxAcross = 1;
                    int limits[] = {1, 51, 101, 151, 201, 251, 301, 351, 401};
                    std::vector<std::string> types;
                    types.push_back("Double parameters:");
                    types.push_back("Branch and Cut double parameters:");
                    types.push_back("Integer parameters:");
                    types.push_back("Branch and Cut integer parameters:");
                    types.push_back("Keyword parameters:");
                    types.push_back("Branch and Cut keyword parameters:");
                    types.push_back("Actions or string parameters:");
                    types.push_back("Branch and Cut actions:");
                    int iType;
                    for (iType = 0; iType < 8; iType++) {
                        int across = 0;
                        int lengthLine = 0;
                        if ((verbose % 4) != 0)
                            std::cout << std::endl;
                        std::cout << types[iType] << std::endl;
                        if ((verbose&2) != 0)
                            std::cout << std::endl;
                        for ( iParam = 0; iParam < numberParameters_; iParam++ ) {
                            int type = parameters_[iParam].type();
                            //printf("%d type %d limits %d %d display %d\n",iParam,
                            //     type,limits[iType],limits[iType+1],parameters_[iParam].displayThis());
                            if ((parameters_[iParam].displayThis() >= printLevel || evenHidden) &&
                                    type >= limits[iType]
                                    && type < limits[iType+1]) {
                                // but skip if not useful for ampl (and in ampl mode)
                                if (verbose >= 4 && (parameters_[iParam].whereUsed()&4) == 0)
                                    continue;
                                if (!across) {
                                    if ((verbose&2) != 0)
                                        std::cout << "Command ";
                                }
                                int length = parameters_[iParam].lengthMatchName() + 1;
                                if (lengthLine + length > 80) {
                                    std::cout << std::endl;
                                    across = 0;
                                    lengthLine = 0;
                                }
                                std::cout << " " << parameters_[iParam].matchName();
                                lengthLine += length;
                                across++;
                                if (across == maxAcross) {
                                    across = 0;
                                    if ((verbose % 4) != 0) {
                                        // put out description as well
                                        if ((verbose&1) != 0)
                                            std::cout << " " << parameters_[iParam].shortHelp();
                                        std::cout << std::endl;
                                        if ((verbose&2) != 0) {
                                            std::cout << "---- description" << std::endl;
                                            parameters_[iParam].printLongHelp();
                                            std::cout << "----" << std::endl << std::endl;
                                        }
                                    } else {
                                        std::cout << std::endl;
                                    }
                                }
                            }
                        }
                        if (across)
                            std::cout << std::endl;
                    }
                } else if (type == CBC_PARAM_FULLGENERALQUERY) {
                    std::cout << "Full list of commands is:" << std::endl;
                    int maxAcross = 5;
                    int limits[] = {1, 51, 101, 151, 201, 251, 301, 351, 401};
                    std::vector<std::string> types;
                    types.push_back("Double parameters:");
                    types.push_back("Branch and Cut double parameters:");
                    types.push_back("Integer parameters:");
                    types.push_back("Branch and Cut integer parameters:");
                    types.push_back("Keyword parameters:");
                    types.push_back("Branch and Cut keyword parameters:");
                    types.push_back("Actions or string parameters:");
                    types.push_back("Branch and Cut actions:");
                    int iType;
                    for (iType = 0; iType < 8; iType++) {
                        int across = 0;
                        std::cout << types[iType] << "  ";
                        for ( iParam = 0; iParam < numberParameters_; iParam++ ) {
                            int type = parameters_[iParam].type();
                            if (type >= limits[iType]
                                    && type < limits[iType+1]) {
                                if (!across)
                                    std::cout << "  ";
                                std::cout << parameters_[iParam].matchName() << "  ";
                                across++;
                                if (across == maxAcross) {
                                    std::cout << std::endl;
                                    across = 0;
                                }
                            }
                        }
                        if (across)
                            std::cout << std::endl;
                    }
                } else if (type < 101) {
                    // get next field as double
                    double value = CoinReadGetDoubleField(argc, argv, &valid);
                    if (!valid) {
                        if (type < 51) {
                            int returnCode;
                            const char * message =
                                parameters_[iParam].setDoubleParameterWithMessage(lpSolver, value, returnCode);
                            if (!noPrinting_ && strlen(message)) {
                                generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                << message
                                << CoinMessageEol;
                            }
                        } else if (type < 81) {
                            int returnCode;
                            const char * message =
                                parameters_[iParam].setDoubleParameterWithMessage(model_, value, returnCode);
                            if (!noPrinting_ && strlen(message)) {
                                generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                << message
                                << CoinMessageEol;
                            }
                        } else {
                            int returnCode;
                            const char * message =
                                parameters_[iParam].setDoubleParameterWithMessage(lpSolver, value, returnCode);
                            if (!noPrinting_ && strlen(message)) {
                                generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                << message
                                << CoinMessageEol;
                            }
                            switch (type) {
                            case CBC_PARAM_DBL_DJFIX:
                                djFix = value;
#ifndef CBC_OTHER_SOLVER
                                if (goodModel && djFix < 1.0e20) {
                                    // do some fixing
                                    clpSolver = dynamic_cast< OsiClpSolverInterface*> (model_.solver());
                                    clpSolver->initialSolve();
                                    lpSolver = clpSolver->getModelPtr();
                                    int numberColumns = lpSolver->numberColumns();
                                    int i;
                                    const char * type = lpSolver->integerInformation();
                                    double * lower = lpSolver->columnLower();
                                    double * upper = lpSolver->columnUpper();
                                    double * solution = lpSolver->primalColumnSolution();
                                    double * dj = lpSolver->dualColumnSolution();
                                    int numberFixed = 0;
                                    double dextra4 = parameters_[whichParam(CBC_PARAM_DBL_DEXTRA4, numberParameters_, parameters_)].doubleValue();
                                    if (dextra4)
                                        printf("Multiple for continuous dj fixing is %g\n", dextra4);
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
                                            } else if (value > upper[i] - 1.0e-5 && djValue < -djFix) {
                                                solution[i] = upper[i];
                                                lower[i] = upper[i];
                                                numberFixed++;
                                            }
                                        }
                                    }
                                    sprintf(generalPrint, "%d columns fixed\n", numberFixed);
                                    generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                    << generalPrint
                                    << CoinMessageEol;
                                }
#endif
                                break;
                            case CBC_PARAM_DBL_TIGHTENFACTOR:
                                tightenFactor = value;
                                if (!complicatedInteger)
                                    defaultSettings = false; // user knows what she is doing
                                break;
                            default:
                                break;
                            }
                        }
                    } else if (valid == 1) {
                        std::cout << " is illegal for double parameter " << parameters_[iParam].name() << " value remains " <<
                                  parameters_[iParam].doubleValue() << std::endl;
                    } else {
                        std::cout << parameters_[iParam].name() << " has value " <<
                                  parameters_[iParam].doubleValue() << std::endl;
                    }
                } else if (type < 201) {
                    // get next field as int
                    int value = CoinReadGetIntField(argc, argv, &valid);
                    if (!valid) {
                        if (type < 151) {
                            if (parameters_[iParam].type() == CLP_PARAM_INT_PRESOLVEPASS)
                                preSolve = value;
                            else if (parameters_[iParam].type() == CLP_PARAM_INT_IDIOT)
                                doIdiot = value;
                            else if (parameters_[iParam].type() == CLP_PARAM_INT_SPRINT)
                                doSprint = value;
                            else if (parameters_[iParam].type() == CLP_PARAM_INT_OUTPUTFORMAT)
                                outputFormat = value;
                            else if (parameters_[iParam].type() == CLP_PARAM_INT_SLPVALUE)
                                slpValue = value;
                            else if (parameters_[iParam].type() == CLP_PARAM_INT_CPP)
                                cppValue = value;
                            else if (parameters_[iParam].type() == CLP_PARAM_INT_PRESOLVEOPTIONS)
                                presolveOptions = value;
                            else if (parameters_[iParam].type() == CLP_PARAM_INT_PRINTOPTIONS)
                                printOptions = value;
                            else if (parameters_[iParam].type() == CLP_PARAM_INT_SUBSTITUTION)
                                substitution = value;
                            else if (parameters_[iParam].type() == CLP_PARAM_INT_DUALIZE)
                                dualize = value;
                            else if (parameters_[iParam].type() == CLP_PARAM_INT_PROCESSTUNE)
                                tunePreProcess = value;
                            else if (parameters_[iParam].type() == CLP_PARAM_INT_USESOLUTION)
                                useSolution = value;
                            else if (parameters_[iParam].type() == CLP_PARAM_INT_VERBOSE)
                                verbose = value;
                            int returnCode;
                            const char * message =
                                parameters_[iParam].setIntParameterWithMessage(lpSolver, value, returnCode);
                            if (!noPrinting_ && strlen(message)) {
                                generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                << message
                                << CoinMessageEol;
                            }
                        } else {
                            if (parameters_[iParam].type() == CBC_PARAM_INT_CUTPASS)
                                cutPass = value;
                            else if (parameters_[iParam].type() == CBC_PARAM_INT_CUTPASSINTREE)
                                cutPassInTree = value;
                            else if (parameters_[iParam].type() == CBC_PARAM_INT_STRONGBRANCHING ||
                                     parameters_[iParam].type() == CBC_PARAM_INT_NUMBERBEFORE)
                                strongChanged = true;
                            else if (parameters_[iParam].type() == CBC_PARAM_INT_FPUMPTUNE ||
                                     parameters_[iParam].type() == CBC_PARAM_INT_FPUMPTUNE2 ||
                                     parameters_[iParam].type() == CBC_PARAM_INT_FPUMPITS)
                                pumpChanged = true;
                            else if (parameters_[iParam].type() == CBC_PARAM_INT_EXPERIMENT) {
  			        int addFlags=0;
				if (value>=10) {
				  addFlags = 1048576*(value/10);
				  value = value % 10;
				  parameters[whichParam(CBC_PARAM_INT_EXPERIMENT, numberParameters, parameters)].setIntValue(value);
				}
                                if (value >= 1) {
				    int values[]={24003,280003,792003,24003,24003};
				    if (value>=2&&value<=3) {
				      // swap default diving
				      int iParam = whichParam(CBC_PARAM_STR_DIVINGC, numberParameters_, parameters_);
				      parameters_[iParam].setCurrentOption("off");
				      iParam = whichParam(CBC_PARAM_STR_DIVINGP, numberParameters_, parameters_);
				      parameters_[iParam].setCurrentOption("on");
				    }
				    int extra4 = values[value-1]+addFlags;
                                    parameters[whichParam(CBC_PARAM_INT_EXTRA4, numberParameters, parameters)].setIntValue(extra4);
                                    if (!noPrinting_) {
                                        generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                        << "switching on global root cuts for gomory and knapsack"
                                        << CoinMessageEol;
                                        generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                        << "using OSL factorization"
                                        << CoinMessageEol;
                                        generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                        << "extra options - -rens on -extra4 "
					<<extra4<<" -passc 1000!"
                                        << CoinMessageEol;
                                    }
                                    parameters[whichParam(CBC_PARAM_STR_PROBINGCUTS, numberParameters, parameters)].setCurrentOption("forceOnStrong");
                                    probingAction = 8;
                                    parameters_[whichParam(CBC_PARAM_STR_GOMORYCUTS, numberParameters_, parameters_)].setCurrentOption("onGlobal");
                                    gomoryAction = 5;
                                    parameters_[whichParam(CBC_PARAM_STR_KNAPSACKCUTS, numberParameters_, parameters_)].setCurrentOption("onGlobal");
                                    knapsackAction = 5;
                                    parameters_[whichParam(CLP_PARAM_STR_FACTORIZATION, numberParameters_, parameters_)].setCurrentOption("osl");
                                    lpSolver->factorization()->forceOtherFactorization(3);
                                    parameters_[whichParam(CBC_PARAM_INT_MAXHOTITS, numberParameters_, parameters_)].setIntValue(100);
                                    parameters[whichParam(CBC_PARAM_INT_CUTPASS, numberParameters, parameters)].setIntValue(1000);
                                    cutPass = 1000;
                                    parameters[whichParam(CBC_PARAM_STR_RENS, numberParameters, parameters)].setCurrentOption("on");
                                }
                            } else if (parameters_[iParam].type() == CBC_PARAM_INT_STRATEGY) {
                                if (value == 0) {
                                    gomoryGen.setAwayAtRoot(0.05);
                                    int iParam;
                                    iParam = whichParam(CBC_PARAM_INT_DIVEOPT, numberParameters_, parameters_);
                                    parameters_[iParam].setIntValue(-1);
                                    iParam = whichParam(CBC_PARAM_INT_FPUMPITS, numberParameters_, parameters_);
                                    parameters_[iParam].setIntValue(20);
                                    iParam = whichParam(CBC_PARAM_INT_FPUMPTUNE, numberParameters_, parameters_);
                                    parameters_[iParam].setIntValue(1003);
                                    initialPumpTune = 1003;
                                    iParam = whichParam(CLP_PARAM_INT_PROCESSTUNE, numberParameters_, parameters_);
                                    parameters_[iParam].setIntValue(-1);
                                    tunePreProcess = 0;
                                    iParam = whichParam(CBC_PARAM_STR_DIVINGC, numberParameters_, parameters_);
                                    parameters_[iParam].setCurrentOption("off");
                                    iParam = whichParam(CBC_PARAM_STR_RINS, numberParameters_, parameters_);
                                    parameters_[iParam].setCurrentOption("off");
                                    iParam = whichParam(CBC_PARAM_STR_PROBINGCUTS, numberParameters_, parameters_);
                                    parameters_[iParam].setCurrentOption("on");
                                    probingAction = 1;
                                }
                            }
                            int returnCode;
                            const char * message =
                                parameters_[iParam].setIntParameterWithMessage(model_, value, returnCode);
                            if (!noPrinting_ && strlen(message)) {
                                generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                << message
                                << CoinMessageEol;
                            }
                        }
                    } else if (valid == 1) {
                        std::cout << " is illegal for integer parameter " << parameters_[iParam].name() << " value remains " <<
                                  parameters_[iParam].intValue() << std::endl;
                    } else {
                        std::cout << parameters_[iParam].name() << " has value " <<
                                  parameters_[iParam].intValue() << std::endl;
                    }
                } else if (type < 301) {
                    // one of several strings
                    std::string value = CoinReadGetString(argc, argv);
                    int action = parameters_[iParam].parameterOption(value);
                    if (action < 0) {
                        if (value != "EOL") {
                            // no match
                            parameters_[iParam].printOptions();
                        } else {
                            // print current value
                            std::cout << parameters_[iParam].name() << " has value " <<
                                      parameters_[iParam].currentOption() << std::endl;
                        }
                    } else {
                        const char * message =
                            parameters_[iParam].setCurrentOptionWithMessage(action);
                        if (!noPrinting_ && strlen(message)) {
                            generalMessageHandler->message(CLP_GENERAL, generalMessages)
                            << message
                            << CoinMessageEol;
                        }
                        // for now hard wired
                        switch (type) {
                        case CLP_PARAM_STR_DIRECTION:
                            if (action == 0)
                                lpSolver->setOptimizationDirection(1);
                            else if (action == 1)
                                lpSolver->setOptimizationDirection(-1);
                            else
                                lpSolver->setOptimizationDirection(0);
                            break;
                        case CLP_PARAM_STR_DUALPIVOT:
                            if (action == 0) {
                                ClpDualRowSteepest steep(3);
                                lpSolver->setDualRowPivotAlgorithm(steep);
                            } else if (action == 1) {
                                ClpDualRowDantzig dantzig;
                                //ClpDualRowSteepest dantzig(5);
                                lpSolver->setDualRowPivotAlgorithm(dantzig);
                            } else if (action == 2) {
                                // partial steep
                                ClpDualRowSteepest steep(2);
                                lpSolver->setDualRowPivotAlgorithm(steep);
                            } else {
                                ClpDualRowSteepest steep;
                                lpSolver->setDualRowPivotAlgorithm(steep);
                            }
                            break;
                        case CLP_PARAM_STR_PRIMALPIVOT:
                            if (action == 0) {
                                ClpPrimalColumnSteepest steep(3);
                                lpSolver->setPrimalColumnPivotAlgorithm(steep);
                            } else if (action == 1) {
                                ClpPrimalColumnSteepest steep(0);
                                lpSolver->setPrimalColumnPivotAlgorithm(steep);
                            } else if (action == 2) {
                                ClpPrimalColumnDantzig dantzig;
                                lpSolver->setPrimalColumnPivotAlgorithm(dantzig);
                            } else if (action == 3) {
                                ClpPrimalColumnSteepest steep(4);
                                lpSolver->setPrimalColumnPivotAlgorithm(steep);
                            } else if (action == 4) {
                                ClpPrimalColumnSteepest steep(1);
                                lpSolver->setPrimalColumnPivotAlgorithm(steep);
                            } else if (action == 5) {
                                ClpPrimalColumnSteepest steep(2);
                                lpSolver->setPrimalColumnPivotAlgorithm(steep);
                            } else if (action == 6) {
                                ClpPrimalColumnSteepest steep(10);
                                lpSolver->setPrimalColumnPivotAlgorithm(steep);
                            }
                            break;
                        case CLP_PARAM_STR_SCALING:
                            lpSolver->scaling(action);
                            solver->setHintParam(OsiDoScale, action != 0, OsiHintTry);
                            doScaling = action;
                            break;
                        case CLP_PARAM_STR_AUTOSCALE:
                            lpSolver->setAutomaticScaling(action != 0);
                            break;
                        case CLP_PARAM_STR_SPARSEFACTOR:
                            lpSolver->setSparseFactorization((1 - action) != 0);
                            break;
                        case CLP_PARAM_STR_BIASLU:
                            lpSolver->factorization()->setBiasLU(action);
                            break;
                        case CLP_PARAM_STR_PERTURBATION:
                            if (action == 0)
                                lpSolver->setPerturbation(50);
                            else
                                lpSolver->setPerturbation(100);
                            break;
                        case CLP_PARAM_STR_ERRORSALLOWED:
                            allowImportErrors = action;
                            break;
                        case CLP_PARAM_STR_INTPRINT:
                            printMode = action;
                            break;
                            //case CLP_PARAM_NOTUSED_ALGORITHM:
                            //algorithm  = action;
                            //defaultSettings=false; // user knows what she is doing
                            //abort();
                            //break;
                        case CLP_PARAM_STR_KEEPNAMES:
                            keepImportNames = 1 - action;
                            break;
                        case CLP_PARAM_STR_PRESOLVE:
                            if (action == 0)
                                preSolve = 5;
                            else if (action == 1)
                                preSolve = 0;
                            else if (action == 2)
                                preSolve = 10;
                            else
                                preSolveFile = true;
                            break;
                        case CLP_PARAM_STR_PFI:
                            lpSolver->factorization()->setForrestTomlin(action == 0);
                            break;
                        case CLP_PARAM_STR_FACTORIZATION:
                            lpSolver->factorization()->forceOtherFactorization(action);
                            break;
                        case CLP_PARAM_STR_CRASH:
                            doCrash = action;
                            break;
                        case CLP_PARAM_STR_VECTOR:
                            doVector = action;
                            break;
                        case CLP_PARAM_STR_MESSAGES:
                            lpSolver->messageHandler()->setPrefix(action != 0);
                            break;
                        case CLP_PARAM_STR_CHOLESKY:
                            choleskyType = action;
                            break;
                        case CLP_PARAM_STR_GAMMA:
                            gamma = action;
                            break;
                        case CLP_PARAM_STR_BARRIERSCALE:
                            scaleBarrier = action;
                            break;
                        case CLP_PARAM_STR_KKT:
                            doKKT = action;
                            break;
                        case CLP_PARAM_STR_CROSSOVER:
                            crossover = action;
                            break;
                        case CLP_PARAM_STR_TIME_MODE:
			    model_.setUseElapsedTime(action!=0);
                            break;
                        case CBC_PARAM_STR_SOS:
                            doSOS = action;
                            break;
                        case CBC_PARAM_STR_GOMORYCUTS:
                            defaultSettings = false; // user knows what she is doing
                            gomoryAction = action;
                            break;
                        case CBC_PARAM_STR_PROBINGCUTS:
                            defaultSettings = false; // user knows what she is doing
                            probingAction = action;
                            break;
                        case CBC_PARAM_STR_KNAPSACKCUTS:
                            defaultSettings = false; // user knows what she is doing
                            knapsackAction = action;
                            break;
                        case CBC_PARAM_STR_REDSPLITCUTS:
                            defaultSettings = false; // user knows what she is doing
                            redsplitAction = action;
                            break;
                        case CBC_PARAM_STR_CLIQUECUTS:
                            defaultSettings = false; // user knows what she is doing
                            cliqueAction = action;
                            break;
                        case CBC_PARAM_STR_FLOWCUTS:
                            defaultSettings = false; // user knows what she is doing
                            flowAction = action;
                            break;
                        case CBC_PARAM_STR_MIXEDCUTS:
                            defaultSettings = false; // user knows what she is doing
                            mixedAction = action;
                            break;
                        case CBC_PARAM_STR_TWOMIRCUTS:
                            defaultSettings = false; // user knows what she is doing
                            twomirAction = action;
                            break;
                        case CBC_PARAM_STR_LANDPCUTS:
                            defaultSettings = false; // user knows what she is doing
                            landpAction = action;
                            break;
                        case CBC_PARAM_STR_RESIDCUTS:
                            defaultSettings = false; // user knows what she is doing
                            residualCapacityAction = action;
                            break;
#ifdef ZERO_HALF_CUTS
                        case CBC_PARAM_STR_ZEROHALFCUTS:
                            defaultSettings = false; // user knows what she is doing
                            zerohalfAction = action;
                            break;
#endif
                        case CBC_PARAM_STR_ROUNDING:
                            defaultSettings = false; // user knows what she is doing
                            break;
                        case CBC_PARAM_STR_FPUMP:
                            defaultSettings = false; // user knows what she is doing
                            break;
                        case CBC_PARAM_STR_RINS:
                            break;
                        case CBC_PARAM_STR_DINS:
                            break;
                        case CBC_PARAM_STR_RENS:
                            break;
                        case CBC_PARAM_STR_CUTSSTRATEGY:
                            gomoryAction = action;
                            probingAction = action;
                            knapsackAction = action;
#ifdef ZERO_HALF_CUTS
                            zerohalfAction = action;
#endif
                            cliqueAction = action;
                            flowAction = action;
                            mixedAction = action;
                            twomirAction = action;
                            //landpAction = action;
                            parameters_[whichParam(CBC_PARAM_STR_GOMORYCUTS, numberParameters_, parameters_)].setCurrentOption(action);
                            parameters_[whichParam(CBC_PARAM_STR_PROBINGCUTS, numberParameters_, parameters_)].setCurrentOption(action);
                            parameters_[whichParam(CBC_PARAM_STR_KNAPSACKCUTS, numberParameters_, parameters_)].setCurrentOption(action);
                            parameters_[whichParam(CBC_PARAM_STR_CLIQUECUTS, numberParameters_, parameters_)].setCurrentOption(action);
                            parameters_[whichParam(CBC_PARAM_STR_FLOWCUTS, numberParameters_, parameters_)].setCurrentOption(action);
                            parameters_[whichParam(CBC_PARAM_STR_MIXEDCUTS, numberParameters_, parameters_)].setCurrentOption(action);
                            parameters_[whichParam(CBC_PARAM_STR_TWOMIRCUTS, numberParameters_, parameters_)].setCurrentOption(action);
#ifdef ZERO_HALF_CUTS
                            parameters_[whichParam(CBC_PARAM_STR_ZEROHALFCUTS, numberParameters_, parameters_)].setCurrentOption(action);
#endif
                            if (!action) {
                                redsplitAction = action;
                                parameters_[whichParam(CBC_PARAM_STR_REDSPLITCUTS, numberParameters_, parameters_)].setCurrentOption(action);
                                landpAction = action;
                                parameters_[whichParam(CBC_PARAM_STR_LANDPCUTS, numberParameters_, parameters_)].setCurrentOption(action);
                                residualCapacityAction = action;
                                parameters_[whichParam(CBC_PARAM_STR_RESIDCUTS, numberParameters_, parameters_)].setCurrentOption(action);
                            }
                            break;
                        case CBC_PARAM_STR_HEURISTICSTRATEGY:
                            parameters_[whichParam(CBC_PARAM_STR_ROUNDING, numberParameters_, parameters_)].setCurrentOption(action);
                            parameters_[whichParam(CBC_PARAM_STR_GREEDY, numberParameters_, parameters_)].setCurrentOption(action);
                            parameters_[whichParam(CBC_PARAM_STR_COMBINE, numberParameters_, parameters_)].setCurrentOption(action);
                            //parameters_[whichParam(CBC_PARAM_STR_LOCALTREE,numberParameters_,parameters_)].setCurrentOption(action);
                            parameters_[whichParam(CBC_PARAM_STR_FPUMP, numberParameters_, parameters_)].setCurrentOption(action);
                            parameters_[whichParam(CBC_PARAM_STR_DIVINGC, numberParameters_, parameters_)].setCurrentOption(action);
                            parameters_[whichParam(CBC_PARAM_STR_RINS, numberParameters_, parameters_)].setCurrentOption(action);
                            break;
                        case CBC_PARAM_STR_GREEDY:
                        case CBC_PARAM_STR_DIVINGS:
                        case CBC_PARAM_STR_DIVINGC:
                        case CBC_PARAM_STR_DIVINGF:
                        case CBC_PARAM_STR_DIVINGG:
                        case CBC_PARAM_STR_DIVINGL:
                        case CBC_PARAM_STR_DIVINGP:
                        case CBC_PARAM_STR_DIVINGV:
                        case CBC_PARAM_STR_COMBINE:
                        case CBC_PARAM_STR_PIVOTANDCOMPLEMENT:
                        case CBC_PARAM_STR_PIVOTANDFIX:
                        case CBC_PARAM_STR_RANDROUND:
                        case CBC_PARAM_STR_LOCALTREE:
                        case CBC_PARAM_STR_NAIVE:
                        case CBC_PARAM_STR_CPX:
                            defaultSettings = false; // user knows what she is doing
                            break;
                        case CBC_PARAM_STR_COSTSTRATEGY:
                            useCosts = action;
                            break;
                        case CBC_PARAM_STR_NODESTRATEGY:
                            nodeStrategy = action;
                            break;
                        case CBC_PARAM_STR_PREPROCESS:
                            preProcess = action;
                            break;
                        default:
                            //abort();
                            break;
                        }
                    }
                } else {
                    // action
                    if (type == CLP_PARAM_ACTION_EXIT) {
#ifdef COIN_HAS_ASL
                        if (statusUserFunction_[0]) {
                            if (info.numberIntegers || info.numberBinary) {
                                // integer
                            } else {
                                // linear
                            }
                            writeAmpl(&info);
                            freeArrays2(&info);
                            freeArgs(&info);
                        }
#endif
                        break; // stop all
                    }
                    switch (type) {
                    case CLP_PARAM_ACTION_DUALSIMPLEX:
                    case CLP_PARAM_ACTION_PRIMALSIMPLEX:
                    case CLP_PARAM_ACTION_SOLVECONTINUOUS:
                    case CLP_PARAM_ACTION_BARRIER:
                        if (goodModel) {
                            // Say not in integer
                            integerStatus = -1;
                            double objScale =
                                parameters_[whichParam(CLP_PARAM_DBL_OBJSCALE2, numberParameters_, parameters_)].doubleValue();
                            if (objScale != 1.0) {
                                int iColumn;
                                int numberColumns = lpSolver->numberColumns();
                                double * dualColumnSolution =
                                    lpSolver->dualColumnSolution();
                                ClpObjective * obj = lpSolver->objectiveAsObject();
                                assert(dynamic_cast<ClpLinearObjective *> (obj));
                                double offset;
                                double * objective = obj->gradient(NULL, NULL, offset, true);
                                for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                                    dualColumnSolution[iColumn] *= objScale;
                                    objective[iColumn] *= objScale;;
                                }
                                int iRow;
                                int numberRows = lpSolver->numberRows();
                                double * dualRowSolution =
                                    lpSolver->dualRowSolution();
                                for (iRow = 0; iRow < numberRows; iRow++)
                                    dualRowSolution[iRow] *= objScale;
                                lpSolver->setObjectiveOffset(objScale*lpSolver->objectiveOffset());
                            }
                            ClpSolve::SolveType method;
                            ClpSolve::PresolveType presolveType;
                            ClpSimplex * model2 = lpSolver;
                            if (dualize) {
                                bool tryIt = true;
                                double fractionColumn = 1.0;
                                double fractionRow = 1.0;
                                if (dualize == 3) {
                                    dualize = 1;
                                    int numberColumns = lpSolver->numberColumns();
                                    int numberRows = lpSolver->numberRows();
                                    if (numberRows < 50000 || 5*numberColumns > numberRows) {
                                        tryIt = false;
                                    } else {
                                        fractionColumn = 0.1;
                                        fractionRow = 0.1;
                                    }
                                }
                                if (tryIt) {
                                    model2 = static_cast<ClpSimplexOther *> (model2)->dualOfModel(fractionRow, fractionColumn);
                                    if (model2) {
                                        sprintf(generalPrint, "Dual of model has %d rows and %d columns",
                                                model2->numberRows(), model2->numberColumns());
                                        generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                        << generalPrint
                                        << CoinMessageEol;
                                        model2->setOptimizationDirection(1.0);
                                    } else {
                                        model2 = lpSolver;
                                        dualize = 0;
                                    }
                                } else {
                                    dualize = 0;
                                }
                            }
                            if (noPrinting_)
                                lpSolver->setLogLevel(0);
                            ClpSolve solveOptions;
                            solveOptions.setPresolveActions(presolveOptions);
                            solveOptions.setSubstitution(substitution);
                            if (preSolve != 5 && preSolve) {
                                presolveType = ClpSolve::presolveNumber;
                                if (preSolve < 0) {
                                    preSolve = - preSolve;
                                    if (preSolve <= 100) {
                                        presolveType = ClpSolve::presolveNumber;
                                        sprintf(generalPrint, "Doing %d presolve passes - picking up non-costed slacks",
                                                preSolve);
                                        generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                        << generalPrint
                                        << CoinMessageEol;
                                        solveOptions.setDoSingletonColumn(true);
                                    } else {
                                        preSolve -= 100;
                                        presolveType = ClpSolve::presolveNumberCost;
                                        sprintf(generalPrint, "Doing %d presolve passes - picking up costed slacks",
                                                preSolve);
                                        generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                        << generalPrint
                                        << CoinMessageEol;
                                    }
                                }
                            } else if (preSolve) {
                                presolveType = ClpSolve::presolveOn;
                            } else {
                                presolveType = ClpSolve::presolveOff;
                            }
                            solveOptions.setPresolveType(presolveType, preSolve);
                            if (type == CLP_PARAM_ACTION_DUALSIMPLEX ||
                                    type == CLP_PARAM_ACTION_SOLVECONTINUOUS) {
                                method = ClpSolve::useDual;
                            } else if (type == CLP_PARAM_ACTION_PRIMALSIMPLEX) {
                                method = ClpSolve::usePrimalorSprint;
                            } else {
                                method = ClpSolve::useBarrier;
                                if (crossover == 1) {
                                    method = ClpSolve::useBarrierNoCross;
                                } else if (crossover == 2) {
                                    ClpObjective * obj = lpSolver->objectiveAsObject();
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
                                ClpMatrixBase * matrix = lpSolver->clpMatrix();
                                if (dynamic_cast< ClpPackedMatrix*>(matrix)) {
                                    ClpPackedMatrix * clpMatrix = dynamic_cast< ClpPackedMatrix*>(matrix);
                                    clpMatrix->makeSpecialColumnCopy();
                                }
                            }
                            if (method == ClpSolve::useDual) {
                                // dual
                                if (doCrash)
                                    solveOptions.setSpecialOption(0, 1, doCrash); // crash
                                else if (doIdiot)
                                    solveOptions.setSpecialOption(0, 2, doIdiot); // possible idiot
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
                                            solveOptions.setSpecialOption(1, 9); // all slack or sprint
                                    } else {
                                        if (doSprint == 0)
                                            solveOptions.setSpecialOption(1, 8); // all slack or idiot
                                        else
                                            solveOptions.setSpecialOption(1, 7); // initiative
                                    }
                                }
                                if (basisHasValues == -1)
                                    solveOptions.setSpecialOption(1, 11); // switch off values
                            } else if (method == ClpSolve::useBarrier || method == ClpSolve::useBarrierNoCross) {
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
                            OsiSolverInterface * coinSolver = model_.solver();
                            OsiSolverLink * linkSolver = dynamic_cast< OsiSolverLink*> (coinSolver);
                            if (!linkSolver) {
                                model2->initialSolve(solveOptions);
                            } else {
                                // special solver
                                int testOsiOptions = parameters_[whichParam(CBC_PARAM_INT_TESTOSI, numberParameters_, parameters_)].intValue();
                                double * solution = NULL;
                                if (testOsiOptions < 10) {
                                    solution = linkSolver->nonlinearSLP(slpValue > 0 ? slpValue : 20 , 1.0e-5);
                                } else if (testOsiOptions >= 10) {
                                    CoinModel coinModel = *linkSolver->coinModel();
                                    ClpSimplex * tempModel = approximateSolution(coinModel, slpValue > 0 ? slpValue : 50 , 1.0e-5, 0);
                                    assert (tempModel);
                                    solution = CoinCopyOfArray(tempModel->primalColumnSolution(), coinModel.numberColumns());
                                    model2->setObjectiveValue(tempModel->objectiveValue());
                                    model2->setProblemStatus(tempModel->problemStatus());
                                    model2->setSecondaryStatus(tempModel->secondaryStatus());
                                    delete tempModel;
                                }
                                if (solution) {
                                    memcpy(model2->primalColumnSolution(), solution,
                                           CoinMin(model2->numberColumns(), linkSolver->coinModel()->numberColumns())*sizeof(double));
                                    delete [] solution;
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
                                   -1 - unknown e.g. before solve or if postSolve says not optimal
                                   0 - optimal
                                   1 - primal infeasible
                                   2 - dual infeasible
                                   3 - stopped on iterations or time
                                   4 - stopped due to errors
                                   5 - stopped by event handler (virtual int ClpEventHandler::event()) */
                                /* cbc status
                                   -1 before branchAndBound
                                   0 finished - check isProvenOptimal or isProvenInfeasible to see if solution found
                                   (or check value of best solution)
                                   1 stopped - on maxnodes, maxsols, maxtime
                                   2 difficulties so run was abandoned
                                   (5 event user programmed event occurred) */
                                /* clp secondary status of problem - may get extended
                                   0 - none
                                   1 - primal infeasible because dual limit reached OR probably primal
                                   infeasible but can't prove it (main status 4)
                                   2 - scaled problem optimal - unscaled problem has primal infeasibilities
                                   3 - scaled problem optimal - unscaled problem has dual infeasibilities
                                   4 - scaled problem optimal - unscaled problem has primal and dual infeasibilities
                                   5 - giving up in primal with flagged variables
                                   6 - failed due to empty problem check
                                   7 - postSolve says not optimal
                                   8 - failed due to bad element check
                                   9 - status was 3 and stopped on time
                                   100 up - translation of enum from ClpEventHandler
                                */
                                /* cbc secondary status of problem
                                   -1 unset (status_ will also be -1)
                                   0 search completed with solution
                                   1 linear relaxation not feasible (or worse than cutoff)
                                   2 stopped on gap
                                   3 stopped on nodes
                                   4 stopped on time
                                   5 stopped on user event
                                   6 stopped on solutions
                                   7 linear relaxation unbounded
                                */
                                int iStatus = model2->status();
                                int iStatus2 = model2->secondaryStatus();
                                if (iStatus == 0) {
                                    iStatus2 = 0;
                                    if (found.type() == CBC_PARAM_ACTION_BAB) {
                                        // set best solution in model as no integers
                                        model_.setBestSolution(model2->primalColumnSolution(),
                                                               model2->numberColumns(),
                                                               model2->getObjValue()*
                                                               model2->getObjSense());
                                    }
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
				if ((iStatus == 2 || iStatus2 > 0) &&
				    !noPrinting_) {
				   std::string statusName[] = {"", "Stopped on ", "Run abandoned", "", "", "User ctrl-c"};
				   std::string minor[] = {"Optimal solution found", "Linear relaxation infeasible", "Optimal solution found (within gap tolerance)", "node limit", "time limit", "user ctrl-c", "solution limit", "Linear relaxation unbounded", "Problem proven infeasible"};
				   sprintf(generalPrint, "\nResult - %s%s\n\n", 
					   statusName[iStatus].c_str(), 
					   minor[iStatus2].c_str());
				   sprintf(generalPrint + strlen(generalPrint),
					   "Enumerated nodes:           0\n");
				   sprintf(generalPrint + strlen(generalPrint), 
					   "Total iterations:           0\n");
#ifndef CBC_QUIET 
				   sprintf(generalPrint + strlen(generalPrint),
					   "Time (CPU seconds):         %.2f\n", 
					   CoinCpuTime() - time0);
				   sprintf(generalPrint + strlen(generalPrint),
					   "Time (Wallclock Seconds):   %.2f\n", 
					   CoinGetTimeOfDay()-time0Elapsed);
#endif
				   generalMessageHandler->message(CLP_GENERAL, generalMessages)
				      << generalPrint
				      << CoinMessageEol;
				}
				//assert (lpSolver==clpSolver->getModelPtr());
                                assert (clpSolver == model_.solver());
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
                                int returnCode = static_cast<ClpSimplexOther *> (lpSolver)->restoreFromDual(model2);
                                if (model2->status() == 3)
                                    returnCode = 0;
                                delete model2;
                                if (returnCode && dualize != 2)
                                    lpSolver->primal(1);
                                model2 = lpSolver;
                            }
#ifdef COIN_HAS_ASL
                            if (statusUserFunction_[0]) {
                                double value = model2->getObjValue() * model2->getObjSense();
                                char buf[300];
                                int pos = 0;
                                int iStat = model2->status();
                                if (iStat == 0) {
                                    pos += sprintf(buf + pos, "optimal," );
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
                                info.problemStatus = iStat;
                                info.objValue = value;
                                pos += sprintf(buf + pos, " objective %.*g", ampl_obj_prec(),
                                               value);
                                sprintf(buf + pos, "\n%d iterations",
                                        model2->getIterationCount());
                                free(info.primalSolution);
                                int numberColumns = model2->numberColumns();
                                info.primalSolution = reinterpret_cast<double *> (malloc(numberColumns * sizeof(double)));
                                CoinCopyN(model2->primalColumnSolution(), numberColumns, info.primalSolution);
                                int numberRows = model2->numberRows();
                                free(info.dualSolution);
                                info.dualSolution = reinterpret_cast<double *> (malloc(numberRows * sizeof(double)));
                                CoinCopyN(model2->dualRowSolution(), numberRows, info.dualSolution);
                                CoinWarmStartBasis * basis = model2->getBasis();
                                free(info.rowStatus);
                                info.rowStatus = reinterpret_cast<int *> (malloc(numberRows * sizeof(int)));
                                free(info.columnStatus);
                                info.columnStatus = reinterpret_cast<int *> (malloc(numberColumns * sizeof(int)));
                                // Put basis in
                                int i;
                                // free,basic,ub,lb are 0,1,2,3
                                for (i = 0; i < numberRows; i++) {
                                    CoinWarmStartBasis::Status status = basis->getArtifStatus(i);
                                    info.rowStatus[i] = status;
                                }
                                for (i = 0; i < numberColumns; i++) {
                                    CoinWarmStartBasis::Status status = basis->getStructStatus(i);
                                    info.columnStatus[i] = status;
                                }
                                // put buffer into info
                                strcpy(info.buffer, buf);
                                delete basis;
                            }
#endif
                        } else {
#ifndef DISALLOW_PRINTING
                            std::cout << "** Current model not valid" << std::endl;
#endif
                        }
                        break;
                    case CLP_PARAM_ACTION_STATISTICS:
                        if (goodModel) {
                            // If presolve on look at presolved
                            bool deleteModel2 = false;
                            ClpSimplex * model2 = lpSolver;
                            if (preSolve) {
                                ClpPresolve pinfo;
                                int presolveOptions2 = presolveOptions&~0x40000000;
                                if ((presolveOptions2&0xffff) != 0)
                                    pinfo.setPresolveActions(presolveOptions2);
                                pinfo.setSubstitution(substitution);
                                if ((printOptions&1) != 0)
                                    pinfo.statistics();
                                double presolveTolerance =
                                    parameters_[whichParam(CLP_PARAM_DBL_PRESOLVETOLERANCE, numberParameters_, parameters_)].doubleValue();
                                model2 =
                                    pinfo.presolvedModel(*lpSolver, presolveTolerance,
                                                         true, preSolve);
                                if (model2) {
                                    printf("Statistics for presolved model\n");
                                    deleteModel2 = true;
                                } else {
                                    printf("Presolved model looks infeasible - will use unpresolved\n");
                                    model2 = lpSolver;
                                }
                            } else {
                                printf("Statistics for unpresolved model\n");
                                model2 =  lpSolver;
                            }
                            statistics(lpSolver, model2);
                            if (deleteModel2)
                                delete model2;
                        } else {
#ifndef DISALLOW_PRINTING
                            std::cout << "** Current model not valid" << std::endl;
#endif
                        }
                        break;
                    case CLP_PARAM_ACTION_TIGHTEN:
                        if (goodModel) {
                            int numberInfeasibilities = lpSolver->tightenPrimalBounds();
                            if (numberInfeasibilities)
                                std::cout << "** Analysis indicates model infeasible" << std::endl;
                        } else {
#ifndef DISALLOW_PRINTING
                            std::cout << "** Current model not valid" << std::endl;
#endif
                        }
                        break;
                    case CLP_PARAM_ACTION_PLUSMINUS:
                        if (goodModel) {
                            ClpMatrixBase * saveMatrix = lpSolver->clpMatrix();
                            ClpPackedMatrix* clpMatrix =
                                dynamic_cast< ClpPackedMatrix*>(saveMatrix);
                            if (clpMatrix) {
                                ClpPlusMinusOneMatrix * newMatrix = new ClpPlusMinusOneMatrix(*(clpMatrix->matrix()));
                                if (newMatrix->getIndices()) {
                                    lpSolver->replaceMatrix(newMatrix);
                                    delete saveMatrix;
                                    std::cout << "Matrix converted to +- one matrix" << std::endl;
                                } else {
                                    std::cout << "Matrix can not be converted to +- 1 matrix" << std::endl;
                                }
                            } else {
                                std::cout << "Matrix not a ClpPackedMatrix" << std::endl;
                            }
                        } else {
#ifndef DISALLOW_PRINTING
                            std::cout << "** Current model not valid" << std::endl;
#endif
                        }
                        break;
                    case CLP_PARAM_ACTION_OUTDUPROWS:
                        dominatedCuts = true;
#ifdef JJF_ZERO
                        if (goodModel) {
                            int numberRows = clpSolver->getNumRows();
                            //int nOut = outDupRow(clpSolver);
                            CglDuplicateRow dupcuts(clpSolver);
                            storedCuts = dupcuts.outDuplicates(clpSolver) != 0;
                            int nOut = numberRows - clpSolver->getNumRows();
                            if (nOut && !noPrinting_)
                                sprintf(generalPrint, "%d rows eliminated", nOut);
                            generalMessageHandler->message(CLP_GENERAL, generalMessages)
                            << generalPrint
                            << CoinMessageEol;
                        } else {
#ifndef DISALLOW_PRINTING
                            std::cout << "** Current model not valid" << std::endl;
#endif
                        }
#endif
                        break;
                    case CLP_PARAM_ACTION_NETWORK:
                        if (goodModel) {
                            ClpMatrixBase * saveMatrix = lpSolver->clpMatrix();
                            ClpPackedMatrix* clpMatrix =
                                dynamic_cast< ClpPackedMatrix*>(saveMatrix);
                            if (clpMatrix) {
                                ClpNetworkMatrix * newMatrix = new ClpNetworkMatrix(*(clpMatrix->matrix()));
                                if (newMatrix->getIndices()) {
                                    lpSolver->replaceMatrix(newMatrix);
                                    delete saveMatrix;
                                    std::cout << "Matrix converted to network matrix" << std::endl;
                                } else {
                                    std::cout << "Matrix can not be converted to network matrix" << std::endl;
                                }
                            } else {
                                std::cout << "Matrix not a ClpPackedMatrix" << std::endl;
                            }
                        } else {
#ifndef DISALLOW_PRINTING
                            std::cout << "** Current model not valid" << std::endl;
#endif
                        }
                        break;
                    case CBC_PARAM_ACTION_DOHEURISTIC:
                        if (goodModel) {
                            int vubAction = parameters_[whichParam(CBC_PARAM_INT_VUBTRY, numberParameters_, parameters_)].intValue();
                            if (vubAction != -1) {
                                // look at vubs
                                // extra1 is number of ints to leave free
                                // Just ones which affect >= extra3
                                int extra3 = parameters_[whichParam(CBC_PARAM_INT_EXTRA3, numberParameters_, parameters_)].intValue();
                                /* 2 is cost above which to fix if feasible
                                   3 is fraction of integer variables fixed if relaxing (0.97)
                                   4 is fraction of all variables fixed if relaxing (0.0)
                                */
                                double dextra[6];
                                int extra[5];
                                extra[1] = parameters_[whichParam(CBC_PARAM_INT_EXTRA1, numberParameters_, parameters_)].intValue();
                                int exp1 = parameters_[whichParam(CBC_PARAM_INT_EXPERIMENT, numberParameters_,
                                                                  parameters_)].intValue();
                                if (exp1 == 4 && extra[1] == -1)
                                    extra[1] = 999998;
                                dextra[1] = parameters_[whichParam(CBC_PARAM_DBL_FAKEINCREMENT, numberParameters_, parameters_)].doubleValue();
                                dextra[2] = parameters_[whichParam(CBC_PARAM_DBL_FAKECUTOFF, numberParameters_, parameters_)].doubleValue();
                                dextra[3] = parameters_[whichParam(CBC_PARAM_DBL_DEXTRA3, numberParameters_, parameters_)].doubleValue();
                                dextra[4] = parameters_[whichParam(CBC_PARAM_DBL_DEXTRA4, numberParameters_, parameters_)].doubleValue();
                                dextra[5] = parameters_[whichParam(CBC_PARAM_DBL_DEXTRA5, numberParameters_, parameters_)].doubleValue();
                                if (!dextra[3])
                                    dextra[3] = 0.97;
                                //OsiClpSolverInterface * newSolver =
                                fixVubs(model_, extra3, vubAction, generalMessageHandler,
                                        debugValues, dextra, extra);
                                //assert (!newSolver);
                            }
                            // Actually do heuristics
                            doHeuristics(&model_, 2, parameters_,
                                         numberParameters_, noPrinting_, initialPumpTune);
                            if (model_.bestSolution()) {
                                model_.setProblemStatus(1);
                                model_.setSecondaryStatus(6);
#ifdef COIN_HAS_ASL
                                if (statusUserFunction_[0]) {
                                    double value = model_.getObjValue();
                                    char buf[300];
                                    int pos = 0;
                                    pos += sprintf(buf + pos, "feasible,");
                                    info.problemStatus = 0;
                                    info.objValue = value;
                                    pos += sprintf(buf + pos, " objective %.*g", ampl_obj_prec(),
                                                   value);
                                    sprintf(buf + pos, "\n0 iterations");
                                    free(info.primalSolution);
                                    int numberColumns = lpSolver->numberColumns();
                                    info.primalSolution = reinterpret_cast<double *> (malloc(numberColumns * sizeof(double)));
                                    CoinCopyN(model_.bestSolution(), numberColumns, info.primalSolution);
                                    int numberRows = lpSolver->numberRows();
                                    free(info.dualSolution);
                                    info.dualSolution = reinterpret_cast<double *> (malloc(numberRows * sizeof(double)));
                                    CoinZeroN(info.dualSolution, numberRows);
                                    CoinWarmStartBasis * basis = lpSolver->getBasis();
                                    free(info.rowStatus);
                                    info.rowStatus = reinterpret_cast<int *> (malloc(numberRows * sizeof(int)));
                                    free(info.columnStatus);
                                    info.columnStatus = reinterpret_cast<int *> (malloc(numberColumns * sizeof(int)));
                                    // Put basis in
                                    int i;
                                    // free,basic,ub,lb are 0,1,2,3
                                    for (i = 0; i < numberRows; i++) {
                                        CoinWarmStartBasis::Status status = basis->getArtifStatus(i);
                                        info.rowStatus[i] = status;
                                    }
                                    for (i = 0; i < numberColumns; i++) {
                                        CoinWarmStartBasis::Status status = basis->getStructStatus(i);
                                        info.columnStatus[i] = status;
                                    }
                                    // put buffer into info
                                    strcpy(info.buffer, buf);
                                    delete basis;
                                }
#endif
                            }
                            int returnCode = callBack(&model, 6);
                            if (returnCode) {
                                // exit if user wants
                                delete babModel_;
                                babModel_ = NULL;
                                return returnCode;
                            }
                        }
                        break;
                    case CBC_PARAM_ACTION_MIPLIB:
                        // User can set options - main difference is lack of model and CglPreProcess
                        goodModel = true;
                        /*
                          Run branch-and-cut. First set a few options -- node comparison, scaling.
                          Print elapsed time at the end.
                        */
                    case CBC_PARAM_ACTION_BAB: // branchAndBound
                        // obsolete case STRENGTHEN:
                        if (goodModel) {
                            bool miplib = type == CBC_PARAM_ACTION_MIPLIB;
                            int logLevel = parameters_[slog].intValue();
                            // Reduce printout
                            if (logLevel <= 1) {
                                model_.solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
                            } else {
                                model_.solver()->setHintParam(OsiDoReducePrint, false, OsiHintTry);
                            }
                            {
                                OsiSolverInterface * solver = model_.solver();
#ifndef CBC_OTHER_SOLVER
                                OsiClpSolverInterface * si =
                                    dynamic_cast<OsiClpSolverInterface *>(solver) ;
                                assert (si != NULL);
                                si->getModelPtr()->scaling(doScaling);
                                ClpSimplex * lpSolver = si->getModelPtr();
                                if (doVector) {
                                    ClpMatrixBase * matrix = lpSolver->clpMatrix();
                                    if (dynamic_cast< ClpPackedMatrix*>(matrix)) {
                                        ClpPackedMatrix * clpMatrix = dynamic_cast< ClpPackedMatrix*>(matrix);
                                        clpMatrix->makeSpecialColumnCopy();
                                    }
                                }
#elif CBC_OTHER_SOLVER==1
                                OsiCpxSolverInterface * si =
                                    dynamic_cast<OsiCpxSolverInterface *>(solver) ;
                                assert (si != NULL);
#endif
                                statistics_nrows = si->getNumRows();
                                statistics_ncols = si->getNumCols();
                                statistics_nprocessedrows = si->getNumRows();
                                statistics_nprocessedcols = si->getNumCols();
                                // See if quadratic
#ifndef CBC_OTHER_SOLVER
#ifdef COIN_HAS_LINK
                                if (!complicatedInteger) {
                                    ClpQuadraticObjective * obj = (dynamic_cast< ClpQuadraticObjective*>(lpSolver->objectiveAsObject()));
                                    if (obj) {
                                        preProcess = 0;
                                        int testOsiOptions = parameters_[whichParam(CBC_PARAM_INT_TESTOSI, numberParameters_, parameters_)].intValue();
                                        parameters_[whichParam(CBC_PARAM_INT_TESTOSI, numberParameters_, parameters_)].setIntValue(CoinMax(0, testOsiOptions));
                                        // create coin model
                                        coinModel = lpSolver->createCoinModel();
                                        assert (coinModel);
                                        // load from coin model
                                        OsiSolverLink solver1;
                                        OsiSolverInterface * solver2 = solver1.clone();
                                        model_.assignSolver(solver2, false);
                                        OsiSolverLink * si =
                                            dynamic_cast<OsiSolverLink *>(model_.solver()) ;
                                        assert (si != NULL);
                                        si->setDefaultMeshSize(0.001);
                                        // need some relative granularity
                                        si->setDefaultBound(100.0);
                                        double dextra3 = parameters_[whichParam(CBC_PARAM_DBL_DEXTRA3, numberParameters_, parameters_)].doubleValue();
                                        if (dextra3)
                                            si->setDefaultMeshSize(dextra3);
                                        si->setDefaultBound(1000.0);
                                        si->setIntegerPriority(1000);
                                        si->setBiLinearPriority(10000);
                                        si->setSpecialOptions2(2 + 4 + 8);
                                        CoinModel * model2 = coinModel;
                                        si->load(*model2, true, parameters_[log].intValue());
                                        // redo
                                        solver = model_.solver();
                                        clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
                                        lpSolver = clpSolver->getModelPtr();
                                        clpSolver->messageHandler()->setLogLevel(0) ;
                                        testOsiParameters = 0;
                                        complicatedInteger = 2;  // allow cuts
                                        OsiSolverInterface * coinSolver = model_.solver();
                                        OsiSolverLink * linkSolver = dynamic_cast< OsiSolverLink*> (coinSolver);
                                        if (linkSolver->quadraticModel()) {
                                            ClpSimplex * qp = linkSolver->quadraticModel();
                                            //linkSolver->nonlinearSLP(CoinMax(slpValue,10),1.0e-5);
                                            qp->nonlinearSLP(CoinMax(slpValue, 40), 1.0e-5);
                                            qp->primal(1);
                                            OsiSolverLinearizedQuadratic solver2(qp);
                                            const double * solution = NULL;
                                            // Reduce printout
                                            solver2.setHintParam(OsiDoReducePrint, true, OsiHintTry);
                                            CbcModel model2(solver2);
                                            // Now do requested saves and modifications
                                            CbcModel * cbcModel = & model2;
                                            OsiSolverInterface * osiModel = model2.solver();
                                            OsiClpSolverInterface * osiclpModel = dynamic_cast< OsiClpSolverInterface*> (osiModel);
                                            ClpSimplex * clpModel = osiclpModel->getModelPtr();

                                            // Set changed values

                                            CglProbing probing;
                                            probing.setMaxProbe(10);
                                            probing.setMaxLook(10);
                                            probing.setMaxElements(200);
                                            probing.setMaxProbeRoot(50);
                                            probing.setMaxLookRoot(10);
                                            probing.setRowCuts(3);
                                            probing.setUsingObjective(true);
                                            cbcModel->addCutGenerator(&probing, -1, "Probing", true, false, false, -100, -1, -1);
                                            cbcModel->cutGenerator(0)->setTiming(true);

                                            CglGomory gomory;
                                            gomory.setLimitAtRoot(512);
                                            cbcModel->addCutGenerator(&gomory, -98, "Gomory", true, false, false, -100, -1, -1);
                                            cbcModel->cutGenerator(1)->setTiming(true);

                                            CglKnapsackCover knapsackCover;
                                            cbcModel->addCutGenerator(&knapsackCover, -98, "KnapsackCover", true, false, false, -100, -1, -1);
                                            cbcModel->cutGenerator(2)->setTiming(true);

                                            CglRedSplit redSplit;
                                            cbcModel->addCutGenerator(&redSplit, -99, "RedSplit", true, false, false, -100, -1, -1);
                                            cbcModel->cutGenerator(3)->setTiming(true);

                                            CglClique clique;
                                            clique.setStarCliqueReport(false);
                                            clique.setRowCliqueReport(false);
                                            clique.setMinViolation(0.1);
                                            cbcModel->addCutGenerator(&clique, -98, "Clique", true, false, false, -100, -1, -1);
                                            cbcModel->cutGenerator(4)->setTiming(true);

                                            CglMixedIntegerRounding2 mixedIntegerRounding2;
                                            cbcModel->addCutGenerator(&mixedIntegerRounding2, -98, "MixedIntegerRounding2", true, false, false, -100, -1, -1);
                                            cbcModel->cutGenerator(5)->setTiming(true);

                                            CglFlowCover flowCover;
                                            cbcModel->addCutGenerator(&flowCover, -98, "FlowCover", true, false, false, -100, -1, -1);
                                            cbcModel->cutGenerator(6)->setTiming(true);

                                            CglTwomir twomir;
                                            twomir.setMaxElements(250);
                                            cbcModel->addCutGenerator(&twomir, -99, "Twomir", true, false, false, -100, -1, -1);
                                            cbcModel->cutGenerator(7)->setTiming(true);

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

                                            CbcHeuristicGreedyEquality heuristicGreedyEquality(*cbcModel);
                                            heuristicGreedyEquality.setHeuristicName("greedy equality");
                                            cbcModel->addHeuristic(&heuristicGreedyEquality);

                                            CbcCompareDefault compare;
                                            cbcModel->setNodeComparison(compare);
                                            cbcModel->setNumberBeforeTrust(5);
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
                                            osiclpModel->setHintParam(OsiDoReducePrint, true, OsiHintTry);
                                            // You can save some time by switching off message building
                                            // clpModel->messagesPointer()->setDetailMessages(100,10000,(int *) NULL);

                                            // Solve

                                            cbcModel->initialSolve();
                                            if (clpModel->tightenPrimalBounds() != 0) {
#ifndef DISALLOW_PRINTING
                                                std::cout << "Problem is infeasible - tightenPrimalBounds!" << std::endl;
#endif
                                                break;
                                            }
                                            clpModel->dual();  // clean up
                                            cbcModel->initialSolve();
#ifdef CBC_THREAD
                                            int numberThreads = parameters_[whichParam(CBC_PARAM_INT_THREADS, numberParameters_, parameters_)].intValue();
                                            cbcModel->setNumberThreads(numberThreads % 100);
                                            cbcModel->setThreadMode(CoinMin(numberThreads / 100, 7));
#endif
                                            //setCutAndHeuristicOptions(*cbcModel);
                                            cbcModel->branchAndBound();
                                            OsiSolverLinearizedQuadratic * solver3 = dynamic_cast<OsiSolverLinearizedQuadratic *> (model2.solver());
                                            assert (solver3);
                                            solution = solver3->bestSolution();
                                            double bestObjectiveValue = solver3->bestObjectiveValue();
                                            linkSolver->setBestObjectiveValue(bestObjectiveValue);
                                            linkSolver->setBestSolution(solution, solver3->getNumCols());
                                            CbcHeuristicDynamic3 dynamic(model_);
                                            dynamic.setHeuristicName("dynamic pass thru");
                                            model_.addHeuristic(&dynamic);
                                            // if convex
                                            if ((linkSolver->specialOptions2()&4) != 0) {
                                                int numberColumns = coinModel->numberColumns();
                                                assert (linkSolver->objectiveVariable() == numberColumns);
                                                // add OA cut
                                                double offset;
                                                double * gradient = new double [numberColumns+1];
                                                memcpy(gradient, qp->objectiveAsObject()->gradient(qp, solution, offset, true, 2),
                                                       numberColumns*sizeof(double));
                                                double rhs = 0.0;
                                                int * column = new int[numberColumns+1];
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
                                                storedAmpl.addCut(-COIN_DBL_MAX, offset + 1.0e-7, n, column, gradient);
                                                delete [] gradient;
                                                delete [] column;
                                            }
                                            // could do three way branching round a) continuous b) best solution
                                            printf("obj %g\n", bestObjectiveValue);
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
                                    model_.solver()->setHintParam(OsiDoPresolveInInitial, false, OsiHintTry);
                                    model_.solver()->setHintParam(OsiDoPresolveInResolve, false, OsiHintTry);
                                }
                                double time1a = CoinCpuTime();
                                OsiSolverInterface * solver = model_.solver();
#ifndef CBC_OTHER_SOLVER
                                OsiClpSolverInterface * si =
                                    dynamic_cast<OsiClpSolverInterface *>(solver) ;
                                if (si)
                                    si->setSpecialOptions(si->specialOptions() | 1024);
#endif
                                model_.initialSolve();
#ifndef CBC_OTHER_SOLVER
                                ClpSimplex * clpSolver = si->getModelPtr();
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
                                int returnCode = callBack(&model_, 1);
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
                                    if (!noPrinting_) {
                                        iStatus = clpSolver->status();
                                        const char * msg[] = {"infeasible", "unbounded", "stopped",
                                                              "difficulties", "other"
                                                             };
                                        sprintf(generalPrint, "Problem is %s - %.2f seconds",
                                                msg[iStatus-1], CoinCpuTime() - time1a);
                                        generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                        << generalPrint
                                        << CoinMessageEol;
                                    }
                                    break;
                                }
                                clpSolver->setSpecialOptions(clpSolver->specialOptions() | IN_BRANCH_AND_BOUND); // say is Cbc (and in branch and bound)
#elif CBC_OTHER_SOLVER==1
#endif
                                if (!noPrinting_) {
                                    sprintf(generalPrint, "Continuous objective value is %g - %.2f seconds",
                                            solver->getObjValue(), CoinCpuTime() - time1a);
                                    generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                    << generalPrint
                                    << CoinMessageEol;
                                }
                                if (model_.getMaximumNodes() == -987654321) {
                                    // See if No objective!
                                    int numberColumns = clpSolver->getNumCols();
                                    const double * obj = clpSolver->getObjCoefficients();
                                    const double * lower = clpSolver->getColLower();
                                    const double * upper = clpSolver->getColUpper();
                                    int nObj = 0;
                                    for (int i = 0; i < numberColumns; i++) {
                                        if (upper[i] > lower[i] && obj[i])
                                            nObj++;
                                    }
                                    if (!nObj) {
                                        printf("************No objective!!\n");
                                        model_.setMaximumSolutions(1);
                                        // Column copy
                                        CoinPackedMatrix  matrixByCol(*model_.solver()->getMatrixByCol());
                                        //const double * element = matrixByCol.getElements();
                                        //const int * row = matrixByCol.getIndices();
                                        //const CoinBigIndex * columnStart = matrixByCol.getVectorStarts();
                                        const int * columnLength = matrixByCol.getVectorLengths();
                                        for (int i = 0; i < numberColumns; i++) {
                                            double value = (CoinDrand48() + 0.5) * 10000;
                                            value = 10;
                                            value *= columnLength[i];
                                            int iValue = static_cast<int> (value) / 10;
                                            //iValue=1;
                                            clpSolver->setObjCoeff(i, iValue);
                                        }
                                    }
                                }
#ifndef CBC_OTHER_SOLVER
                                if (!complicatedInteger && preProcess == 0 && clpSolver->tightenPrimalBounds(0.0, 0, true) != 0) {
#ifndef DISALLOW_PRINTING
                                    std::cout << "Problem is infeasible - tightenPrimalBounds!" << std::endl;
#endif
                                    model_.setProblemStatus(0);
                                    model_.setSecondaryStatus(1);
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
                                    const double * rowPrimal = temp.primalRowSolution();
                                    const double * rowLower = temp.rowLower();
                                    const double * rowUpper = temp.rowUpper();
                                    const double * rowScale = temp.rowScale();
                                    int iRow;
                                    for (iRow = 0; iRow < numberRows; iRow++) {
                                        double value = rowPrimal[iRow];
                                        double above = value - rowLower[iRow];
                                        double below = rowUpper[iRow] - value;
                                        if (above < 1.0e12) {
                                            largest = CoinMax(largest, above);
                                        }
                                        if (below < 1.0e12) {
                                            largest = CoinMax(largest, below);
                                        }
                                        if (rowScale) {
                                            double multiplier = rowScale[iRow];
                                            above *= multiplier;
                                            below *= multiplier;
                                        }
                                        if (above < 1.0e12) {
                                            largestScaled = CoinMax(largestScaled, above);
                                        }
                                        if (below < 1.0e12) {
                                            largestScaled = CoinMax(largestScaled, below);
                                        }
                                    }

                                    int numberColumns = temp.numberColumns();
                                    const double * columnPrimal = temp.primalColumnSolution();
                                    const double * columnLower = temp.columnLower();
                                    const double * columnUpper = temp.columnUpper();
                                    const double * columnScale = temp.columnScale();
                                    int iColumn;
                                    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                                        double value = columnPrimal[iColumn];
                                        double above = value - columnLower[iColumn];
                                        double below = columnUpper[iColumn] - value;
                                        if (above < 1.0e12) {
                                            largest = CoinMax(largest, above);
                                        }
                                        if (below < 1.0e12) {
                                            largest = CoinMax(largest, below);
                                        }
                                        if (columnScale) {
                                            double multiplier = 1.0 / columnScale[iColumn];
                                            above *= multiplier;
                                            below *= multiplier;
                                        }
                                        if (above < 1.0e12) {
                                            largestScaled = CoinMax(largestScaled, above);
                                        }
                                        if (below < 1.0e12) {
                                            largestScaled = CoinMax(largestScaled, below);
                                        }
                                    }
#ifdef COIN_DEVELOP
                                    if (!noPrinting_)
                                        std::cout << "Largest (scaled) away from bound " << largestScaled
                                                  << " unscaled " << largest << std::endl;
#endif
                                    clpSolver->setDualBound(CoinMax(1.0001e8, CoinMin(100.0*largest, 1.00001e10)));
                                }
                                si->resolve();  // clean up
#endif
                            }
                            // If user made settings then use them
                            if (!defaultSettings) {
                                OsiSolverInterface * solver = model_.solver();
                                if (!doScaling)
                                    solver->setHintParam(OsiDoScale, false, OsiHintTry);
#ifndef CBC_OTHER_SOLVER
                                OsiClpSolverInterface * si =
                                    dynamic_cast<OsiClpSolverInterface *>(solver) ;
                                assert (si != NULL);
                                // get clp itself
                                ClpSimplex * modelC = si->getModelPtr();
                                //if (modelC->tightenPrimalBounds()!=0) {
                                //std::cout<<"Problem is infeasible!"<<std::endl;
                                //break;
                                //}
                                // bounds based on continuous
                                if (tightenFactor && !complicatedInteger) {
                                    if (modelC->tightenPrimalBounds(tightenFactor) != 0) {
#ifndef DISALLOW_PRINTING
                                        std::cout << "Problem is infeasible!" << std::endl;
#endif
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
                            OsiSolverInterface * saveSolver = NULL;
                            CglPreProcess process;
			    cbcPreProcessPointer = & process;
                            // Say integers in sync
                            bool integersOK = true;
                            delete babModel_;
                            babModel_ = new CbcModel(model_);
#ifndef CBC_OTHER_SOLVER
                            int numberChanged = 0;
                            OsiSolverInterface * solver3 = clpSolver->clone();
                            babModel_->assignSolver(solver3);
                            OsiClpSolverInterface * clpSolver2 = dynamic_cast< OsiClpSolverInterface*> (babModel_->solver());
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
                                    frequency = base + cutoff1 / freq0 + (numberRows - cutoff1) / freq1;
                                else
                                    frequency = base + cutoff1 / freq0 + (cutoff2 - cutoff1) / freq1 + (numberRows - cutoff2) / freq2;
                                lpSolver->setFactorizationFrequency(CoinMin(maximum, frequency));
                            }
#elif CBC_OTHER_SOLVER==1
                            OsiSolverInterface * solver3 = model_.solver()->clone();
                            babModel_->assignSolver(solver3);
#endif
                            time2 = CoinCpuTime();
                            totalTime += time2 - time1;
                            //time1 = time2;
                            double timeLeft = babModel_->getMaximumSeconds();
                            int numberOriginalColumns = babModel_->solver()->getNumCols();
                            if (preProcess == 7) {
                                // use strategy instead
                                preProcess = 0;
                                useStrategy = true;
#ifdef COIN_HAS_LINK
                                // empty out any cuts
                                if (storedAmpl.sizeRowCuts()) {
                                    printf("Emptying ampl stored cuts as internal preprocessing\n");
                                    CglStored temp;
                                    storedAmpl = temp;
                                }
#endif
                            }
                            if (preProcess && type == CBC_PARAM_ACTION_BAB) {
			      // see whether to switch off preprocessing
			      // only allow SOS and integer
			      OsiObject ** objects = babModel_->objects();
			      int numberObjects = babModel_->numberObjects();
			      for (int iObj = 0; iObj < numberObjects; iObj++) {
				CbcSOS * objSOS =
				  dynamic_cast <CbcSOS *>(objects[iObj]) ;
				CbcSimpleInteger * objSimpleInteger =
				  dynamic_cast <CbcSimpleInteger *>(objects[iObj]) ;
				if (!objSimpleInteger&&!objSOS) {
				  preProcess=0;
				  break;
				}
			      }
			    }
                            if (type == CBC_PARAM_ACTION_BAB) {
                                double limit;
                                clpSolver->getDblParam(OsiDualObjectiveLimit, limit);
                                if (clpSolver->getObjValue()*clpSolver->getObjSense() >=
                                        limit*clpSolver->getObjSense())
                                    preProcess = 0;
                            }
                            if (preProcess && type == CBC_PARAM_ACTION_BAB) {
#ifndef CBC_OTHER_SOLVER
                                // See if sos from mps file
                                if (numberSOS == 0 && clpSolver->numberSOS() && doSOS) {
                                    // SOS
                                    numberSOS = clpSolver->numberSOS();
                                    const CoinSet * setInfo = clpSolver->setInfo();
                                    sosStart = new int [numberSOS+1];
                                    sosType = new char [numberSOS];
                                    int i;
                                    int nTotal = 0;
                                    sosStart[0] = 0;
                                    for ( i = 0; i < numberSOS; i++) {
                                        int type = setInfo[i].setType();
                                        int n = setInfo[i].numberEntries();
                                        sosType[i] = static_cast<char>(type);
                                        nTotal += n;
                                        sosStart[i+1] = nTotal;
                                    }
                                    sosIndices = new int[nTotal];
                                    sosReference = new double [nTotal];
                                    for (i = 0; i < numberSOS; i++) {
                                        int n = setInfo[i].numberEntries();
                                        const int * which = setInfo[i].which();
                                        const double * weights = setInfo[i].weights();
                                        int base = sosStart[i];
                                        for (int j = 0; j < n; j++) {
                                            int k = which[j];
                                            sosIndices[j+base] = k;
                                            sosReference[j+base] = weights ? weights[j] : static_cast<double> (j);
                                        }
                                    }
                                }
#endif
                                saveSolver = babModel_->solver()->clone();
                                /* Do not try and produce equality cliques and
                                   do up to 10 passes */
                                OsiSolverInterface * solver2;
                                {
                                    // Tell solver we are in Branch and Cut
                                    saveSolver->setHintParam(OsiDoInBranchAndCut, true, OsiHintDo) ;
                                    // Default set of cut generators
                                    CglProbing generator1;
                                    generator1.setUsingObjective(1);
                                    generator1.setMaxPass(1);
                                    generator1.setMaxPassRoot(1);
                                    generator1.setMaxProbeRoot(CoinMin(3000, saveSolver->getNumCols()));
                                    generator1.setMaxElements(100);
                                    generator1.setMaxElementsRoot(200);
                                    generator1.setMaxLookRoot(50);
                                    if (saveSolver->getNumCols() > 3000)
                                        generator1.setMaxProbeRoot(123);
                                    generator1.setRowCuts(3);
                                    if ((tunePreProcess&1) != 0) {
                                        // heavy probing
                                        generator1.setMaxPassRoot(2);
                                        generator1.setMaxElements(300);
                                        generator1.setMaxProbeRoot(saveSolver->getNumCols());
                                    }
                                    if ((babModel_->specialOptions()&65536) != 0)
                                        process.setOptions(1);
                                    // Add in generators
				    if ((model_.moreSpecialOptions()&65536)==0)
				      process.addCutGenerator(&generator1);
                                    int translate[] = {9999, 0, 0, -3, 2, 3, -2, 9999, 4, 5};
                                    process.passInMessageHandler(babModel_->messageHandler());
                                    //process.messageHandler()->setLogLevel(babModel_->logLevel());
#ifdef COIN_HAS_ASL
                                    if (info.numberSos && doSOS && statusUserFunction_[0]) {
                                        // SOS
                                        numberSOS = info.numberSos;
                                        sosStart = info.sosStart;
                                        sosIndices = info.sosIndices;
                                    }
#endif
                                    if (numberSOS && doSOS) {
                                        // SOS
                                        int numberColumns = saveSolver->getNumCols();
                                        char * prohibited = new char[numberColumns];
                                        memset(prohibited, 0, numberColumns);
                                        int n = sosStart[numberSOS];
                                        for (int i = 0; i < n; i++) {
                                            int iColumn = sosIndices[i];
                                            prohibited[iColumn] = 1;
                                        }
                                        process.passInProhibited(prohibited, numberColumns);
                                        delete [] prohibited;
                                    }
                                    if (!model_.numberObjects() && true) {
                                        /* model may not have created objects
                                           If none then create
                                        */
                                        model_.findIntegers(true);
                                    }
                                    if (model_.numberObjects()) {
                                        OsiObject ** oldObjects = babModel_->objects();
                                        int numberOldObjects = babModel_->numberObjects();
                                        // SOS
                                        int numberColumns = saveSolver->getNumCols();
                                        char * prohibited = new char[numberColumns];
                                        memset(prohibited, 0, numberColumns);
                                        int numberProhibited = 0;
                                        for (int iObj = 0; iObj < numberOldObjects; iObj++) {
                                            CbcSOS * obj =
                                                dynamic_cast <CbcSOS *>(oldObjects[iObj]) ;
                                            if (obj) {
                                                int n = obj->numberMembers();
                                                const int * which = obj->members();
                                                for (int i = 0; i < n; i++) {
                                                    int iColumn = which[i];
                                                    prohibited[iColumn] = 1;
                                                    numberProhibited++;
                                                }
                                            }
                                            CbcLotsize * obj2 =
                                                dynamic_cast <CbcLotsize *>(oldObjects[iObj]) ;
                                            if (obj2) {
                                                int iColumn = obj2->columnNumber();
                                                prohibited[iColumn] = 1;
                                                numberProhibited++;
                                            }
                                        }
                                        if (numberProhibited)
                                            process.passInProhibited(prohibited, numberColumns);
                                        delete [] prohibited;
                                    }
                                    int numberPasses = 10;
                                    if (tunePreProcess >= 1000000) {
                                        numberPasses = (tunePreProcess / 1000000) - 1;
                                        tunePreProcess = tunePreProcess % 1000000;
                                    } else if (tunePreProcess >= 1000) {
                                        numberPasses = (tunePreProcess / 1000) - 1;
                                        tunePreProcess = tunePreProcess % 1000;
                                    }
#ifndef CBC_OTHER_SOLVER
                                    if (doSprint > 0) {
                                        // Sprint for primal solves
                                        ClpSolve::SolveType method = ClpSolve::usePrimalorSprint;
                                        ClpSolve::PresolveType presolveType = ClpSolve::presolveOff;
                                        int numberPasses = 5;
                                        int options[] = {0, 3, 0, 0, 0, 0};
                                        int extraInfo[] = { -1, 20, -1, -1, -1, -1};
                                        extraInfo[1] = doSprint;
                                        int independentOptions[] = {0, 0, 3};
                                        ClpSolve clpSolve(method, presolveType, numberPasses,
                                                          options, extraInfo, independentOptions);
                                        // say use in OsiClp
                                        clpSolve.setSpecialOption(6, 1);
                                        OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (saveSolver);
                                        osiclp->setSolveOptions(clpSolve);
                                        osiclp->setHintParam(OsiDoDualInResolve, false);
                                        // switch off row copy
                                        osiclp->getModelPtr()->setSpecialOptions(osiclp->getModelPtr()->specialOptions() | 256);
                                        osiclp->getModelPtr()->setInfeasibilityCost(1.0e11);
                                    }
#endif
#ifndef CBC_OTHER_SOLVER
                                    {
                                        OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (saveSolver);
                                        osiclp->setSpecialOptions(osiclp->specialOptions() | 1024);
                                        int savePerturbation = osiclp->getModelPtr()->perturbation();
                                        //#define CBC_TEMP1
#ifdef CBC_TEMP1
                                        if (savePerturbation == 50)
                                            osiclp->getModelPtr()->setPerturbation(52); // try less
#endif
					if ((model_.moreSpecialOptions()&65536)!=0)
					  process.setOptions(2+4+8); // no cuts
                                        solver2 = process.preProcessNonDefault(*saveSolver, translate[preProcess], numberPasses,
                                                                               tunePreProcess);
                                        /*solver2->writeMps("after");
                                          saveSolver->writeMps("before");*/
                                        osiclp->getModelPtr()->setPerturbation(savePerturbation);
                                    }
#elif CBC_OTHER_SOLVER==1
                                    solver2 = process.preProcessNonDefault(*saveSolver, translate[preProcess], numberPasses,
                                                                           tunePreProcess);
#endif
                                    integersOK = false; // We need to redo if CbcObjects exist
                                    // Tell solver we are not in Branch and Cut
                                    saveSolver->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo) ;
                                    if (solver2)
                                        solver2->setHintParam(OsiDoInBranchAndCut, false, OsiHintDo) ;
                                }
#ifdef COIN_HAS_ASL
                                if (!solver2 && statusUserFunction_[0]) {
                                    // infeasible
                                    info.problemStatus = 1;
                                    info.objValue = 1.0e100;
                                    sprintf(info.buffer, "infeasible/unbounded by pre-processing");
                                    info.primalSolution = NULL;
                                    info.dualSolution = NULL;
                                    break;
                                }
#endif
                                if (!noPrinting_) {
                                    if (!solver2) {
                                        sprintf(generalPrint, "Pre-processing says infeasible or unbounded");
                                        generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                        << generalPrint
                                        << CoinMessageEol;
                                    } else {
                                        //printf("processed model has %d rows, %d columns and %d elements\n",
                                        //     solver2->getNumRows(),solver2->getNumCols(),solver2->getNumElements());
                                    }
                                }
                                if (!solver2) {
                                    // say infeasible for solution
                                    integerStatus = 6;
                                    model_.setProblemStatus(0);
                                    model_.setSecondaryStatus(1);
                                    babModel_->setProblemStatus(0);
                                    babModel_->setSecondaryStatus(1);
                                } else {
                                    statistics_nprocessedrows = solver2->getNumRows();
                                    statistics_nprocessedcols = solver2->getNumCols();
                                    model_.setProblemStatus(-1);
                                    babModel_->setProblemStatus(-1);
                                }
                                int returnCode = callBack(babModel_, 2);
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
                                    const int * originalColumns = process.originalColumns();
                                    int numberColumns = solver2->getNumCols();
                                    double * bestSolution = babModel_->bestSolution();
                                    const double * oldBestSolution = model_.bestSolution();
                                    for (int i = 0; i < numberColumns; i++) {
                                        int jColumn = originalColumns[i];
                                        bestSolution[i] = oldBestSolution[jColumn];
                                    }
                                }
                                //solver2->resolve();
                                if (preProcess == 2) {
                                    OsiClpSolverInterface * clpSolver2 = dynamic_cast< OsiClpSolverInterface*> (solver2);
                                    ClpSimplex * lpSolver = clpSolver2->getModelPtr();
                                    lpSolver->writeMps("presolved.mps", 0, 1, lpSolver->optimizationDirection());
                                    printf("Preprocessed model (minimization) on presolved.mps\n");
                                }
                                {
                                    // look at new integers
                                    int numberOriginalColumns =
                                        process.originalModel()->getNumCols();
                                    const int * originalColumns = process.originalColumns();
                                    OsiClpSolverInterface * osiclp2 = dynamic_cast< OsiClpSolverInterface*> (solver2);
                                    int numberColumns = osiclp2->getNumCols();
                                    OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (saveSolver);
                                    for (int i = 0; i < numberColumns; i++) {
                                        int iColumn = originalColumns[i];
                                        if (iColumn < numberOriginalColumns) {
                                            if (osiclp2->isInteger(i) && !osiclp->isInteger(iColumn))
                                                osiclp2->setOptionalInteger(i); // say optional
                                        }
                                    }
                                }
                                // we have to keep solver2 so pass clone
                                solver2 = solver2->clone();
                                babModel_->assignSolver(solver2);
                                babModel_->setOriginalColumns(process.originalColumns());
                                babModel_->initialSolve();
                                babModel_->setMaximumSeconds(timeLeft - (CoinCpuTime() - time2));
                            }
                            // now tighten bounds
                            if (!miplib) {
#ifndef CBC_OTHER_SOLVER
                                OsiClpSolverInterface * si =
                                    dynamic_cast<OsiClpSolverInterface *>(babModel_->solver()) ;
                                assert (si != NULL);
                                // get clp itself
                                ClpSimplex * modelC = si->getModelPtr();
                                //if (noPrinting_)
                                //modelC->setLogLevel(0);
                                if (!complicatedInteger && modelC->tightenPrimalBounds() != 0) {
#ifndef DISALLOW_PRINTING
                                    std::cout << "Problem is infeasible!" << std::endl;
#endif
                                    model_.setProblemStatus(0);
                                    model_.setSecondaryStatus(1);
                                    // and in babModel_ if exists
                                    if (babModel_) {
                                        babModel_->setProblemStatus(0);
                                        babModel_->setSecondaryStatus(1);
                                    }
                                    break;
                                }
                                si->resolve();
#elif CBC_OTHER_SOLVER==1
#endif
                            }
                            if (debugValues) {
                                // for debug
                                std::string problemName ;
                                babModel_->solver()->getStrParam(OsiProbName, problemName) ;
                                babModel_->solver()->activateRowCutDebugger(problemName.c_str()) ;
                                twomirGen.probname_ = CoinStrdup(problemName.c_str());
                                // checking seems odd
                                //redsplitGen.set_given_optsol(babModel_->solver()->getRowCutDebuggerAlways()->optimalSolution(),
                                //                         babModel_->getNumCols());
                            }
                            int testOsiOptions = parameters_[whichParam(CBC_PARAM_INT_TESTOSI, numberParameters_, parameters_)].intValue();
                            //#ifdef COIN_HAS_ASL
#ifndef JJF_ONE
                            // If linked then see if expansion wanted
                            {
                                OsiSolverLink * solver3 = dynamic_cast<OsiSolverLink *> (babModel_->solver());
                                int options = parameters_[whichParam(CBC_PARAM_INT_MIPOPTIONS, numberParameters_, parameters_)].intValue() / 10000;
                                if (solver3 || (options&16) != 0) {
                                    if (options) {
                                        /*
                                          1 - force mini branch and bound
                                          2 - set priorities high on continuous
                                          4 - try adding OA cuts
                                          8 - try doing quadratic linearization
                                          16 - try expanding knapsacks
                                        */
                                        if ((options&16)) {
                                            int numberColumns = saveCoinModel.numberColumns();
                                            int numberRows = saveCoinModel.numberRows();
                                            whichColumn = new int[numberColumns];
                                            knapsackStart = new int[numberRows+1];
                                            knapsackRow = new int[numberRows];
                                            numberKnapsack = 10000;
                                            int extra1 = parameters_[whichParam(CBC_PARAM_INT_EXTRA1, numberParameters_, parameters_)].intValue();
                                            int extra2 = parameters_[whichParam(CBC_PARAM_INT_EXTRA2, numberParameters_, parameters_)].intValue();
                                            int logLevel = parameters_[log].intValue();
                                            OsiSolverInterface * solver = expandKnapsack(saveCoinModel, whichColumn, knapsackStart,
                                                                          knapsackRow, numberKnapsack,
                                                                          storedAmpl, logLevel, extra1, extra2,
                                                                          saveTightenedModel);
                                            if (solver) {
#ifndef CBC_OTHER_SOLVER
                                                clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
                                                assert (clpSolver);
                                                lpSolver = clpSolver->getModelPtr();
#endif
                                                babModel_->assignSolver(solver);
                                                testOsiOptions = 0;
                                                // allow gomory
                                                complicatedInteger = 0;
#ifdef COIN_HAS_ASL
                                                // Priorities already done
                                                free(info.priorities);
                                                info.priorities = NULL;
#endif
                                            } else {
                                                numberKnapsack = 0;
                                                delete [] whichColumn;
                                                delete [] knapsackStart;
                                                delete [] knapsackRow;
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
                                int * sort = new int[numberColumns];
                                double * dsort = new double[numberColumns];
                                int * priority = new int [numberColumns];
                                const double * objective = babModel_->getObjCoefficients();
                                const double * lower = babModel_->getColLower() ;
                                const double * upper = babModel_->getColUpper() ;
                                const CoinPackedMatrix * matrix = babModel_->solver()->getMatrixByCol();
                                const int * columnLength = matrix->getVectorLengths();
                                int iColumn;
                                int n = 0;
                                for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                                    if (babModel_->isInteger(iColumn)) {
                                        sort[n] = n;
                                        if (useCosts == 1)
                                            dsort[n++] = -fabs(objective[iColumn]);
                                        else if (useCosts == 2)
                                            dsort[n++] = iColumn;
                                        else if (useCosts == 3)
                                            dsort[n++] = upper[iColumn] - lower[iColumn];
                                        else if (useCosts == 4)
                                            dsort[n++] = -(upper[iColumn] - lower[iColumn]);
                                        else if (useCosts == 5)
                                            dsort[n++] = -columnLength[iColumn];
                                    }
                                }
                                CoinSort_2(dsort, dsort + n, sort);
                                int level = 0;
                                double last = -1.0e100;
                                for (int i = 0; i < n; i++) {
                                    int iPut = sort[i];
                                    if (dsort[i] != last) {
                                        level++;
                                        last = dsort[i];
                                    }
                                    priority[iPut] = level;
                                }
                                babModel_->passInPriorities( priority, false);
                                integersOK = true;
                                delete [] priority;
                                delete [] sort;
                                delete [] dsort;
                            }
                            // Set up heuristics
                            doHeuristics(babModel_, ((!miplib) ? 1 : 10), parameters_,
                                         numberParameters_, noPrinting_, initialPumpTune);
                            if (!miplib) {
                                if (parameters_[whichParam(CBC_PARAM_STR_LOCALTREE, numberParameters_, parameters_)].currentOptionAsInteger()) {
                                    CbcTreeLocal localTree(babModel_, NULL, 10, 0, 0, 10000, 2000);
                                    babModel_->passInTreeHandler(localTree);
                                }
                            }
                            if (type == CBC_PARAM_ACTION_MIPLIB) {
                                if (babModel_->numberStrong() == 5 && babModel_->numberBeforeTrust() == 5)
                                    babModel_->setNumberBeforeTrust(10);
                            }
                            int experimentFlag = parameters_[whichParam(CBC_PARAM_INT_EXPERIMENT, numberParameters_,
                                                             parameters_)].intValue();
                            int strategyFlag = parameters_[whichParam(CBC_PARAM_INT_STRATEGY, numberParameters_,
                                                           parameters_)].intValue();
                            int bothFlags = CoinMax(CoinMin(experimentFlag, 1), strategyFlag);
                            // add cut generators if wanted
                            int switches[30];
                            int accuracyFlag[30];
			    char doAtEnd[30];
			    memset(doAtEnd,0,30);
                            int numberGenerators = 0;
                            int translate[] = { -100, -1, -99, -98, 1, -1098, -999, 1, 1, 1, -1};
                            if (probingAction) {
                                int numberColumns = babModel_->solver()->getNumCols();
                                if (probingAction > 7) {
                                    probingGen.setMaxElements(numberColumns);
                                    probingGen.setMaxElementsRoot(numberColumns);
                                }
                                probingGen.setMaxProbeRoot(CoinMin(2000, numberColumns));
                                probingGen.setMaxProbeRoot(123);
                                probingGen.setMaxProbe(123);
                                probingGen.setMaxLookRoot(20);
                                if (probingAction == 7 || probingAction == 9)
                                    probingGen.setRowCuts(-3); // strengthening etc just at root
                                if (probingAction == 8 || probingAction == 9) {
                                    // Number of unsatisfied variables to look at
                                    probingGen.setMaxProbeRoot(numberColumns);
                                    probingGen.setMaxProbe(numberColumns);
                                    // How far to follow the consequences
                                    probingGen.setMaxLook(50);
                                    probingGen.setMaxLookRoot(50);
                                }
                                if (probingAction == 10) {
                                    probingGen.setMaxPassRoot(2);
                                    probingGen.setMaxProbeRoot(numberColumns);
                                    probingGen.setMaxLookRoot(100);
                                }
                                // If 5 then force on
                                int iAction = translate[probingAction];
                                if (probingAction == 5)
                                    iAction = 1;
                                babModel_->addCutGenerator(&probingGen, iAction, "Probing");
                                accuracyFlag[numberGenerators] = 5;
                                switches[numberGenerators++] = 0;
                            }
                            if (gomoryAction && (complicatedInteger != 1 ||
                                                 (gomoryAction == 1 || gomoryAction >= 4))) {
                                // try larger limit
                                int numberColumns = babModel_->getNumCols();
                                if (gomoryAction == 7) {
                                    gomoryAction = 4;
                                    gomoryGen.setLimitAtRoot(numberColumns);
                                    gomoryGen.setLimit(numberColumns);
                                } else if (gomoryAction == 8) {
                                    gomoryAction = 3;
                                    gomoryGen.setLimitAtRoot(numberColumns);
                                    gomoryGen.setLimit(200);
                                } else if (numberColumns > 5000) {
                                    //#define MORE_CUTS2
#ifdef MORE_CUTS2
                                    // try larger limit
                                    gomoryGen.setLimitAtRoot(numberColumns);
                                    gomoryGen.setLimit(200);
#else
                                    gomoryGen.setLimitAtRoot(2000);
                                    //gomoryGen.setLimit(200);
#endif
                                } else {
#ifdef MORE_CUTS2
                                    // try larger limit
                                    gomoryGen.setLimitAtRoot(numberColumns);
                                    gomoryGen.setLimit(200);
#endif
                                }
                                int cutLength =
                                    parameters_[whichParam(CBC_PARAM_INT_CUTLENGTH, numberParameters_, parameters_)].intValue();
                                if (cutLength != -1) {
                                    gomoryGen.setLimitAtRoot(cutLength);
                                    if (cutLength < 10000000) {
                                        gomoryGen.setLimit(cutLength);
                                    } else {
                                        gomoryGen.setLimit(cutLength % 10000000);
                                    }
                                }
                                int laGomory = parameters_[whichParam(CBC_PARAM_STR_LAGOMORYCUTS, numberParameters_, parameters_)].currentOptionAsInteger();
				int gType = translate[gomoryAction];
				if (!laGomory) {
				  // Normal
				  babModel_->addCutGenerator(&gomoryGen, translate[gomoryAction], "Gomory");
				  accuracyFlag[numberGenerators] = 3;
				  switches[numberGenerators++] = 0;
				} else {
				  laGomory--;
				  int type = (laGomory % 3)+1;
				  int when = laGomory/3;
				  char atEnd = (when<2) ? 1 : 0;
				  int gomoryTypeMajor = 0;
				  if (when<3) {
				    // normal as well
				    babModel_->addCutGenerator(&gomoryGen, gType, "Gomory");
				    accuracyFlag[numberGenerators] = 3;
				    switches[numberGenerators++] = 0;
				    if (when==2)
				      gomoryTypeMajor=10;
				  } else {
				    when--; // so on
				    gomoryTypeMajor=20;
				  }
				  if (!when)
				    gType=-99; // root
				  gomoryGen.passInOriginalSolver(babModel_->solver());
				  if ((type&1) !=0) {
				    // clean
				    gomoryGen.setGomoryType(gomoryTypeMajor+1);
				    babModel_->addCutGenerator(&gomoryGen, gType, "GomoryL1");
				    accuracyFlag[numberGenerators] = 3;
				    doAtEnd[numberGenerators]=atEnd;
				    switches[numberGenerators++] = 0;
				  }
				  if ((type&2) !=0) {
				    // simple
				    gomoryGen.setGomoryType(gomoryTypeMajor+2);
				    babModel_->addCutGenerator(&gomoryGen, gType, "GomoryL2");
				    accuracyFlag[numberGenerators] = 3;
				    doAtEnd[numberGenerators]=atEnd;
				    switches[numberGenerators++] = 0;
				  }
				}
                            }
#ifdef CLIQUE_ANALYSIS
                            if (miplib && !storedAmpl.sizeRowCuts()) {
                                printf("looking at probing\n");
                                babModel_->addCutGenerator(&storedAmpl, 1, "Stored");
                            }
#endif
                            if (knapsackAction) {
                                babModel_->addCutGenerator(&knapsackGen, translate[knapsackAction], "Knapsack");
                                accuracyFlag[numberGenerators] = 1;
                                switches[numberGenerators++] = -2;
                            }
                            if (redsplitAction && !complicatedInteger) {
                                babModel_->addCutGenerator(&redsplitGen, translate[redsplitAction], "Reduce-and-split");
                                accuracyFlag[numberGenerators] = 5;
                                switches[numberGenerators++] = 1;
                            }
                            if (cliqueAction) {
                                babModel_->addCutGenerator(&cliqueGen, translate[cliqueAction], "Clique");
                                accuracyFlag[numberGenerators] = 0;
                                switches[numberGenerators++] = 0;
                            }
                            if (mixedAction) {
                                babModel_->addCutGenerator(&mixedGen, translate[mixedAction], "MixedIntegerRounding2");
                                accuracyFlag[numberGenerators] = 2;
                                switches[numberGenerators++] = 0;
                            }
                            if (flowAction) {
                                babModel_->addCutGenerator(&flowGen, translate[flowAction], "FlowCover");
                                accuracyFlag[numberGenerators] = 2;
                                switches[numberGenerators++] = 1;
                            }
                            if (twomirAction && (complicatedInteger != 1 ||
                                                 (twomirAction == 1 || twomirAction >= 4))) {
                                // try larger limit
                                int numberColumns = babModel_->getNumCols();
                                if (twomirAction == 7) {
                                    twomirAction = 4;
                                    twomirGen.setMaxElements(numberColumns);
                                } else if (numberColumns > 5000 && twomirAction == 4) {
                                    twomirGen.setMaxElements(2000);
                                }
                                babModel_->addCutGenerator(&twomirGen, translate[twomirAction], "TwoMirCuts");
                                accuracyFlag[numberGenerators] = 4;
                                switches[numberGenerators++] = 1;
                            }
#ifndef DEBUG_MALLOC
                            if (landpAction) {
                                babModel_->addCutGenerator(&landpGen, translate[landpAction], "LiftAndProject");
                                accuracyFlag[numberGenerators] = 5;
                                switches[numberGenerators++] = 1;
                            }
#endif
                            if (residualCapacityAction) {
                                babModel_->addCutGenerator(&residualCapacityGen, translate[residualCapacityAction], "ResidualCapacity");
                                accuracyFlag[numberGenerators] = 5;
                                switches[numberGenerators++] = 1;
                            }
#ifdef ZERO_HALF_CUTS
                            if (zerohalfAction) {
                                if (zerohalfAction > 4) {
                                    //zerohalfAction -=4;
                                    zerohalfGen.setFlags(1);
                                }
                                babModel_->addCutGenerator(&zerohalfGen, translate[zerohalfAction], "ZeroHalf");
                                accuracyFlag[numberGenerators] = 5;
                                switches[numberGenerators++] = 2;
                            }
#endif
                            if (dominatedCuts)
                                babModel_->setSpecialOptions(babModel_->specialOptions() | 64);
                            // Say we want timings
                            numberGenerators = babModel_->numberCutGenerators();
                            int iGenerator;
                            int cutDepth =
                                parameters_[whichParam(CBC_PARAM_INT_CUTDEPTH, numberParameters_, parameters_)].intValue();
                            for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
                                CbcCutGenerator * generator = babModel_->cutGenerator(iGenerator);
                                int howOften = generator->howOften();
                                if (howOften == -98 || howOften == -99)
                                    generator->setSwitchOffIfLessThan(switches[iGenerator]);
                                // Use if any at root as more likely later and fairly cheap
                                //if (switches[iGenerator]==-2)
                                //generator->setWhetherToUse(true);
                                generator->setInaccuracy(accuracyFlag[iGenerator]);
				if (doAtEnd[iGenerator]) {
				  generator->setWhetherCallAtEnd(true);
				  generator->setMustCallAgain(true);
				}
                                generator->setTiming(true);
                                if (cutDepth >= 0)
                                    generator->setWhatDepth(cutDepth) ;
                            }
                            // Could tune more
                            if (!miplib) {
                                double minimumDrop =
                                    fabs(babModel_->solver()->getObjValue()) * 1.0e-5 + 1.0e-5;
                                babModel_->setMinimumDrop(CoinMin(5.0e-2, minimumDrop));
                                if (cutPass == -1234567) {
                                    if (babModel_->getNumCols() < 500)
                                        babModel_->setMaximumCutPassesAtRoot(-100); // always do 100 if possible
                                    else if (babModel_->getNumCols() < 5000)
                                        babModel_->setMaximumCutPassesAtRoot(100); // use minimum drop
                                    else
                                        babModel_->setMaximumCutPassesAtRoot(20);
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
                            //if (babModel_->getNumCols()<5000)
                            //babModel_->setNumberStrong(20);
                            // Switch off strong branching if wanted
                            //if (babModel_->getNumCols()>10*babModel_->getNumRows())
                            //babModel_->setNumberStrong(0);
                            if (!noPrinting_) {
                                int iLevel = parameters_[log].intValue();
                                if (iLevel < 0) {
                                    if (iLevel > -10) {
                                        babModel_->setPrintingMode(1);
                                    } else {
                                        babModel_->setPrintingMode(2);
                                        iLevel += 10;
                                        parameters_[log].setIntValue(iLevel);
                                    }
                                    iLevel = -iLevel;
                                }
                                babModel_->messageHandler()->setLogLevel(iLevel);
                                if (babModel_->getNumCols() > 2000 || babModel_->getNumRows() > 1500 ||
                                        babModel_->messageHandler()->logLevel() > 1)
                                    babModel_->setPrintFrequency(100);
                            }

                            babModel_->solver()->setIntParam(OsiMaxNumIterationHotStart,
                                                             parameters_[whichParam(CBC_PARAM_INT_MAXHOTITS, numberParameters_, parameters_)].intValue());
#ifndef CBC_OTHER_SOLVER
                            OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (babModel_->solver());
                            // go faster stripes
                            if ((osiclp->getNumRows() < 300 && osiclp->getNumCols() < 500)) {
                                osiclp->setupForRepeatedUse(2, parameters_[slog].intValue());
                                if (bothFlags >= 1) {
                                    ClpSimplex * lp = osiclp->getModelPtr();
                                    int specialOptions = lp->specialOptions();
                                    lp->setSpecialOptions(specialOptions | (2048 + 4096));
                                }
                            } else {
                                osiclp->setupForRepeatedUse(0, parameters_[slog].intValue());
                            }
                            if (bothFlags >= 2) {
                                ClpSimplex * lp = osiclp->getModelPtr();
                                int specialOptions = lp->specialOptions();
                                lp->setSpecialOptions(specialOptions | (2048 + 4096));
                            }
                            double increment = babModel_->getCutoffIncrement();;
                            int * changed = NULL;
                            if (!miplib && increment == normalIncrement)
                                changed = analyze( osiclp, numberChanged, increment, false, generalMessageHandler, noPrinting);
#elif CBC_OTHER_SOLVER==1
                            double increment = babModel_->getCutoffIncrement();;
#endif
                            if (debugValues) {
                                int numberColumns = babModel_->solver()->getNumCols();
                                if (numberDebugValues == numberColumns) {
                                    // for debug
                                    babModel_->solver()->activateRowCutDebugger(debugValues) ;
                                } else {
                                    int numberOriginalColumns =
                                        process.originalModel()->getNumCols();
                                    if (numberDebugValues <= numberOriginalColumns) {
                                        const int * originalColumns = process.originalColumns();
                                        double * newValues = new double [numberColumns];
                                        // in case preprocess added columns!
                                        // need to find values
                                        OsiSolverInterface * siCopy =
                                            babModel_->solver()->clone();
                                        for (int i = 0; i < numberColumns; i++) {
                                            int jColumn = originalColumns[i];
                                            if (jColumn < numberDebugValues &&
                                                    siCopy->isInteger(i)) {
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
                                                   numberColumns*sizeof(double));
                                        } else {
                                            printf("BAD debug file\n");
                                            siCopy->writeMps("Bad");
                                            exit(22);
                                        }
                                        delete siCopy;
                                        // for debug
                                        babModel_->solver()->activateRowCutDebugger(newValues) ;
                                        delete [] newValues;
                                    } else {
                                        printf("debug file has incorrect number of columns\n");
                                    }
                                }
                            }
                            babModel_->setCutoffIncrement(CoinMax(babModel_->getCutoffIncrement(), increment));
                            // Turn this off if you get problems
                            // Used to be automatically set
                            int mipOptions = parameters_[whichParam(CBC_PARAM_INT_MIPOPTIONS, numberParameters_, parameters_)].intValue() % 10000;
                            if (mipOptions != (1057)) {
                                sprintf(generalPrint, "mip options %d", mipOptions);
                                generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                << generalPrint
                                << CoinMessageEol;
                            }
#ifndef CBC_OTHER_SOLVER
                            osiclp->setSpecialOptions(mipOptions);
#elif CBC_OTHER_SOLVER==1
#endif
                            // probably faster to use a basis to get integer solutions
                            babModel_->setSpecialOptions(babModel_->specialOptions() | 2);
                            currentBranchModel = babModel_;
                            //OsiSolverInterface * strengthenedModel=NULL;
                            if (type == CBC_PARAM_ACTION_BAB ||
                                    type == CBC_PARAM_ACTION_MIPLIB) {
                                if (strategyFlag == 1) {
                                    // try reduced model
                                    babModel_->setSpecialOptions(babModel_->specialOptions() | 512);
                                }
                                if (experimentFlag >= 5 || strategyFlag == 2) {
                                    // try reduced model at root
                                    babModel_->setSpecialOptions(babModel_->specialOptions() | 32768);
                                }
                                {
                                    int depthMiniBab = parameters_[whichParam(CBC_PARAM_INT_DEPTHMINIBAB, numberParameters_, parameters_)].intValue();
                                    if (depthMiniBab != -1)
                                        babModel_->setFastNodeDepth(depthMiniBab);
                                }
                                int extra4 = parameters_[whichParam(CBC_PARAM_INT_EXTRA4, numberParameters_, parameters_)].intValue();
                                if (extra4 >= 0) {
                                    int strategy = extra4 % 10;
                                    extra4 /= 10;
                                    int method = extra4 % 100;
                                    extra4 /= 100;
                                    extra4 = strategy + method * 8 + extra4 * 1024;
                                    babModel_->setMoreSpecialOptions(extra4);
                                }
                                int moreMipOptions = parameters_[whichParam(CBC_PARAM_INT_MOREMIPOPTIONS, numberParameters_, parameters_)].intValue();
                                if (moreMipOptions >= 0) {
                                    sprintf(generalPrint, "more mip options %d", moreMipOptions);
                                    generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                    << generalPrint
                                    << CoinMessageEol;
                                    OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (babModel_->solver());
                                    if (moreMipOptions == 10000) {
                                        // test memory saving
                                        moreMipOptions -= 10000;
                                        ClpSimplex * lpSolver = osiclp->getModelPtr();
                                        lpSolver->setPersistenceFlag(1);
                                        // switch off row copy if few rows
                                        if (lpSolver->numberRows() < 150)
                                            lpSolver->setSpecialOptions(lpSolver->specialOptions() | 256);
                                    }
                                    if (moreMipOptions < 10000 && moreMipOptions) {
                                        if (((moreMipOptions + 1) % 1000000) != 0)
                                            babModel_->setSearchStrategy(moreMipOptions % 1000000);
                                    } else if (moreMipOptions < 100000) {
                                        // try reduced model
                                        babModel_->setSpecialOptions(babModel_->specialOptions() | 512);
                                    }
                                    // go faster stripes
                                    if ( moreMipOptions >= 999999) {
                                        if (osiclp) {
                                            int save = osiclp->specialOptions();
                                            osiclp->setupForRepeatedUse(2, 0);
                                            osiclp->setSpecialOptions(save | osiclp->specialOptions());
                                        }
                                    }
                                }
                            }
                            {
                                int extra1 = parameters_[whichParam(CBC_PARAM_INT_EXTRA1, numberParameters_, parameters_)].intValue();
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
                            if (type == CBC_PARAM_ACTION_BAB) {
#ifdef COIN_HAS_ASL
                                if (statusUserFunction_[0]) {
                                    priorities = info.priorities;
                                    branchDirection = info.branchDirection;
                                    pseudoDown = info.pseudoDown;
                                    pseudoUp = info.pseudoUp;
                                    solutionIn = info.primalSolution;
                                    prioritiesIn = info.priorities;
                                    if (info.numberSos && doSOS) {
                                        // SOS
                                        numberSOS = info.numberSos;
                                        sosStart = info.sosStart;
                                        sosIndices = info.sosIndices;
                                        sosType = info.sosType;
                                        sosReference = info.sosReference;
                                        sosPriority = info.sosPriority;
                                    }
                                }
#endif
                                const int * originalColumns = preProcess ? process.originalColumns() : NULL;
                                if (solutionIn && useSolution >= 0) {
                                    if (!prioritiesIn) {
                                        int n;
                                        if (preProcess) {
                                            int numberColumns = babModel_->getNumCols();
                                            // extend arrays in case SOS
                                            n = originalColumns[numberColumns-1] + 1;
                                        } else {
                                            n = babModel_->getNumCols();
                                        }
                                        prioritiesIn = reinterpret_cast<int *> (malloc(n * sizeof(int)));
                                        for (int i = 0; i < n; i++)
                                            prioritiesIn[i] = 100;
                                    }
                                    if (preProcess) {
                                        int numberColumns = babModel_->getNumCols();
                                        // extend arrays in case SOS
                                        int n = originalColumns[numberColumns-1] + 1;
                                        int nSmaller = CoinMin(n, numberOriginalColumns);
                                        double * solutionIn2 = new double [n];
                                        int * prioritiesIn2 = new int[n];
                                        int i;
                                        for (i = 0; i < nSmaller; i++) {
                                            solutionIn2[i] = solutionIn[i];
                                            prioritiesIn2[i] = prioritiesIn[i];
                                        }
                                        for (; i < n; i++) {
                                            solutionIn2[i] = 0.0;
                                            prioritiesIn2[i] = 1000000;
                                        }
                                        int iLast = -1;
                                        for (i = 0; i < numberColumns; i++) {
                                            int iColumn = originalColumns[i];
                                            assert (iColumn > iLast);
                                            iLast = iColumn;
                                            solutionIn2[i] = solutionIn2[iColumn];
                                            if (prioritiesIn)
                                                prioritiesIn2[i] = prioritiesIn2[iColumn];
                                        }
                                        if (useSolution)
                                            babModel_->setHotstartSolution(solutionIn2, prioritiesIn2);
                                        else
                                            babModel_->setBestSolution(solutionIn2, numberColumns,
                                                                       COIN_DBL_MAX, true);
                                        delete [] solutionIn2;
                                        delete [] prioritiesIn2;
                                    } else {
                                        if (useSolution)
                                            babModel_->setHotstartSolution(solutionIn, prioritiesIn);
                                        else
                                            babModel_->setBestSolution(solutionIn, babModel_->getNumCols(),
                                                                       COIN_DBL_MAX, true);
                                    }
                                }
                                OsiSolverInterface * testOsiSolver = (testOsiOptions >= 0) ? babModel_->solver() : NULL;
                                if (!testOsiSolver) {
                                    // *************************************************************
                                    // CbcObjects
                                    if (preProcess && (process.numberSOS() || babModel_->numberObjects())) {
                                        int numberSOS = process.numberSOS();
                                        int numberIntegers = babModel_->numberIntegers();
                                        /* model may not have created objects
                                           If none then create
                                        */
                                        if (!numberIntegers || !babModel_->numberObjects()) {
                                            int type = (pseudoUp) ? 1 : 0;
                                            babModel_->findIntegers(true, type);
                                            numberIntegers = babModel_->numberIntegers();
                                            integersOK = true;
                                        }
                                        OsiObject ** oldObjects = babModel_->objects();
                                        // Do sets and priorities
                                        OsiObject ** objects = new OsiObject * [numberSOS];
                                        // set old objects to have low priority
                                        int numberOldObjects = babModel_->numberObjects();
                                        int numberColumns = babModel_->getNumCols();
                                        // backward pointer to new variables
                                        // extend arrays in case SOS
                                        assert (originalColumns);
                                        int n = originalColumns[numberColumns-1] + 1;
                                        n = CoinMax(n, CoinMax(numberColumns, numberOriginalColumns));
                                        int * newColumn = new int[n];
                                        int i;
                                        for (i = 0; i < numberOriginalColumns; i++)
                                            newColumn[i] = -1;
                                        for (i = 0; i < numberColumns; i++)
                                            newColumn[originalColumns[i]] = i;
                                        if (!integersOK) {
                                            // Change column numbers etc
                                            int n = 0;
                                            for (int iObj = 0; iObj < numberOldObjects; iObj++) {
                                                int iColumn = oldObjects[iObj]->columnNumber();
                                                if (iColumn < 0 || iColumn >= numberOriginalColumns) {
                                                    oldObjects[n++] = oldObjects[iObj];
                                                } else {
                                                    iColumn = newColumn[iColumn];
                                                    if (iColumn >= 0) {
                                                        CbcSimpleInteger * obj =
                                                            dynamic_cast <CbcSimpleInteger *>(oldObjects[iObj]) ;
                                                        if (obj) {
                                                            obj->setColumnNumber(iColumn);
                                                        } else {
                                                            // only other case allowed is lotsizing
                                                            CbcLotsize * obj2 =
                                                                dynamic_cast <CbcLotsize *>(oldObjects[iObj]) ;
                                                            assert (obj2);
                                                            obj2->setModelSequence(iColumn);
                                                        }
                                                        oldObjects[n++] = oldObjects[iObj];
                                                    } else {
                                                        delete oldObjects[iObj];
                                                    }
                                                }
                                            }
                                            babModel_->setNumberObjects(n);
                                            numberOldObjects = n;
                                            babModel_->zapIntegerInformation();
                                        }
                                        int nMissing = 0;
                                        for (int iObj = 0; iObj < numberOldObjects; iObj++) {
                                            if (process.numberSOS())
                                                oldObjects[iObj]->setPriority(numberColumns + 1);
                                            int iColumn = oldObjects[iObj]->columnNumber();
                                            if (iColumn < 0 || iColumn >= numberOriginalColumns) {
                                                CbcSOS * obj =
                                                    dynamic_cast <CbcSOS *>(oldObjects[iObj]) ;
                                                if (obj) {
                                                    int n = obj->numberMembers();
                                                    int * which = obj->mutableMembers();
                                                    double * weights = obj->mutableWeights();
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
                                                continue;
                                            }
                                            if (originalColumns)
                                                iColumn = originalColumns[iColumn];
                                            if (branchDirection) {
                                                CbcSimpleInteger * obj =
                                                    dynamic_cast <CbcSimpleInteger *>(oldObjects[iObj]) ;
                                                if (obj) {
                                                    obj->setPreferredWay(branchDirection[iColumn]);
                                                } else {
                                                    CbcObject * obj =
                                                        dynamic_cast <CbcObject *>(oldObjects[iObj]) ;
                                                    assert (obj);
                                                    obj->setPreferredWay(branchDirection[iColumn]);
                                                }
                                            }
                                            if (pseudoUp) {
                                                CbcSimpleIntegerPseudoCost * obj1a =
                                                    dynamic_cast <CbcSimpleIntegerPseudoCost *>(oldObjects[iObj]) ;
                                                assert (obj1a);
                                                if (pseudoDown[iColumn] > 0.0)
                                                    obj1a->setDownPseudoCost(pseudoDown[iColumn]);
                                                if (pseudoUp[iColumn] > 0.0)
                                                    obj1a->setUpPseudoCost(pseudoUp[iColumn]);
                                            }
                                        }
                                        if (nMissing) {
                                            sprintf(generalPrint, "%d SOS variables vanished due to pre processing? - check validity?", nMissing);
                                            generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                            << generalPrint
                                            << CoinMessageEol;
                                        }
                                        delete [] newColumn;
                                        const int * starts = process.startSOS();
                                        const int * which = process.whichSOS();
                                        const int * type = process.typeSOS();
                                        const double * weight = process.weightSOS();
                                        int iSOS;
                                        for (iSOS = 0; iSOS < numberSOS; iSOS++) {
                                            int iStart = starts[iSOS];
                                            int n = starts[iSOS+1] - iStart;
                                            objects[iSOS] = new CbcSOS(babModel_, n, which + iStart, weight + iStart,
                                                                       iSOS, type[iSOS]);
                                            // branch on long sets first
                                            objects[iSOS]->setPriority(numberColumns - n);
                                        }
                                        if (numberSOS)
                                            babModel_->addObjects(numberSOS, objects);
                                        for (iSOS = 0; iSOS < numberSOS; iSOS++)
                                            delete objects[iSOS];
                                        delete [] objects;
                                    } else if (priorities || branchDirection || pseudoDown || pseudoUp || numberSOS) {
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
                                            OsiObject ** objects = new OsiObject * [numberSOS];
                                            int iSOS;
                                            if (originalColumns) {
                                                // redo sequence numbers
                                                int numberColumns = babModel_->getNumCols();
                                                int nOld = originalColumns[numberColumns-1] + 1;
                                                int * back = new int[nOld];
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
                                                delete [] back;
                                                if (nMissing) {
                                                    sprintf(generalPrint, "%d SOS variables vanished due to pre processing? - check validity?", nMissing);
                                                    generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                                    << generalPrint
                                                    << CoinMessageEol;
                                                }
                                            }
                                            for (iSOS = 0; iSOS < numberSOS; iSOS++) {
                                                int iStart = sosStart[iSOS];
                                                int n = sosStart[iSOS+1] - iStart;
                                                objects[iSOS] = new CbcSOS(babModel_, n, sosIndices + iStart, sosReference + iStart,
                                                                           iSOS, sosType[iSOS]);
                                                if (sosPriority)
                                                    objects[iSOS]->setPriority(sosPriority[iSOS]);
                                                else if (!prioritiesIn)
                                                    objects[iSOS]->setPriority(10);  // rather than 1000
                                            }
                                            // delete any existing SOS objects
                                            int numberObjects = babModel_->numberObjects();
                                            OsiObject ** oldObjects = babModel_->objects();
                                            int nNew = 0;
                                            for (int i = 0; i < numberObjects; i++) {
                                                OsiObject * objThis = oldObjects[i];
                                                CbcSOS * obj1 =
                                                    dynamic_cast <CbcSOS *>(objThis) ;
                                                OsiSOS * obj2 =
                                                    dynamic_cast <OsiSOS *>(objThis) ;
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
                                            delete [] objects;
                                        }
                                    }
                                    OsiObject ** objects = babModel_->objects();
                                    int numberObjects = babModel_->numberObjects();
                                    for (int iObj = 0; iObj < numberObjects; iObj++) {
                                        // skip sos
                                        CbcSOS * objSOS =
                                            dynamic_cast <CbcSOS *>(objects[iObj]) ;
                                        if (objSOS)
                                            continue;
                                        int iColumn = objects[iObj]->columnNumber();
                                        assert (iColumn >= 0);
                                        if (originalColumns)
                                            iColumn = originalColumns[iColumn];
                                        if (branchDirection) {
                                            CbcSimpleInteger * obj =
                                                dynamic_cast <CbcSimpleInteger *>(objects[iObj]) ;
                                            if (obj) {
                                                obj->setPreferredWay(branchDirection[iColumn]);
                                            } else {
                                                CbcObject * obj =
                                                    dynamic_cast <CbcObject *>(objects[iObj]) ;
                                                assert (obj);
                                                obj->setPreferredWay(branchDirection[iColumn]);
                                            }
                                        }
                                        if (priorities) {
                                            int iPriority = priorities[iColumn];
                                            if (iPriority > 0)
                                                objects[iObj]->setPriority(iPriority);
                                        }
                                        if (pseudoUp && pseudoUp[iColumn]) {
                                            CbcSimpleIntegerPseudoCost * obj1a =
                                                dynamic_cast <CbcSimpleIntegerPseudoCost *>(objects[iObj]) ;
                                            assert (obj1a);
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
                                        //int type = (pseudoUp) ? 1 : 0;
                                        testOsiSolver->findIntegers(false);
                                        numberIntegers = testOsiSolver->getNumIntegers();
                                    }
                                    if (preProcess && process.numberSOS()) {
                                        int numberSOS = process.numberSOS();
                                        OsiObject ** oldObjects = testOsiSolver->objects();
                                        // Do sets and priorities
                                        OsiObject ** objects = new OsiObject * [numberSOS];
                                        // set old objects to have low priority
                                        int numberOldObjects = testOsiSolver->numberObjects();
                                        int numberColumns = testOsiSolver->getNumCols();
                                        for (int iObj = 0; iObj < numberOldObjects; iObj++) {
                                            oldObjects[iObj]->setPriority(numberColumns + 1);
                                            int iColumn = oldObjects[iObj]->columnNumber();
                                            assert (iColumn >= 0);
                                            if (iColumn >= numberOriginalColumns)
                                                continue;
                                            if (originalColumns)
                                                iColumn = originalColumns[iColumn];
                                            if (branchDirection) {
                                                OsiSimpleInteger * obj =
                                                    dynamic_cast <OsiSimpleInteger *>(oldObjects[iObj]) ;
                                                if (obj) {
                                                    obj->setPreferredWay(branchDirection[iColumn]);
                                                } else {
                                                    OsiObject2 * obj =
                                                        dynamic_cast <OsiObject2 *>(oldObjects[iObj]) ;
                                                    if (obj)
                                                        obj->setPreferredWay(branchDirection[iColumn]);
                                                }
                                            }
                                            if (pseudoUp) {
                                                abort();
                                            }
                                        }
                                        const int * starts = process.startSOS();
                                        const int * which = process.whichSOS();
                                        const int * type = process.typeSOS();
                                        const double * weight = process.weightSOS();
                                        int iSOS;
                                        for (iSOS = 0; iSOS < numberSOS; iSOS++) {
                                            int iStart = starts[iSOS];
                                            int n = starts[iSOS+1] - iStart;
                                            objects[iSOS] = new OsiSOS(testOsiSolver, n, which + iStart, weight + iStart,
                                                                       type[iSOS]);
                                            // branch on long sets first
                                            objects[iSOS]->setPriority(numberColumns - n);
                                        }
                                        testOsiSolver->addObjects(numberSOS, objects);
                                        for (iSOS = 0; iSOS < numberSOS; iSOS++)
                                            delete objects[iSOS];
                                        delete [] objects;
                                    } else if (priorities || branchDirection || pseudoDown || pseudoUp || numberSOS) {
                                        if (numberSOS) {
                                            // Do sets and priorities
                                            OsiObject ** objects = new OsiObject * [numberSOS];
                                            int iSOS;
                                            if (originalColumns) {
                                                // redo sequence numbers
                                                int numberColumns = testOsiSolver->getNumCols();
                                                int nOld = originalColumns[numberColumns-1] + 1;
                                                int * back = new int[nOld];
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
                                                delete [] back;
                                                if (nMissing) {
                                                    sprintf(generalPrint, "%d SOS variables vanished due to pre processing? - check validity?", nMissing);
                                                    generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                                    << generalPrint
                                                    << CoinMessageEol;
                                                }
                                            }
                                            for (iSOS = 0; iSOS < numberSOS; iSOS++) {
                                                int iStart = sosStart[iSOS];
                                                int n = sosStart[iSOS+1] - iStart;
                                                objects[iSOS] = new OsiSOS(testOsiSolver, n, sosIndices + iStart, sosReference + iStart,
                                                                           sosType[iSOS]);
                                                if (sosPriority)
                                                    objects[iSOS]->setPriority(sosPriority[iSOS]);
                                                else if (!prioritiesIn)
                                                    objects[iSOS]->setPriority(10);  // rather than 1000
                                            }
                                            // delete any existing SOS objects
                                            int numberObjects = testOsiSolver->numberObjects();
                                            OsiObject ** oldObjects = testOsiSolver->objects();
                                            int nNew = 0;
                                            for (int i = 0; i < numberObjects; i++) {
                                                OsiObject * objThis = oldObjects[i];
                                                OsiSOS * obj1 =
                                                    dynamic_cast <OsiSOS *>(objThis) ;
                                                OsiSOS * obj2 =
                                                    dynamic_cast <OsiSOS *>(objThis) ;
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
                                            delete [] objects;
                                        }
                                    }
                                    OsiObject ** objects = testOsiSolver->objects();
                                    int numberObjects = testOsiSolver->numberObjects();
                                    int logLevel = parameters_[log].intValue();
                                    for (int iObj = 0; iObj < numberObjects; iObj++) {
                                        // skip sos
                                        OsiSOS * objSOS =
                                            dynamic_cast <OsiSOS *>(objects[iObj]) ;
                                        if (objSOS) {
                                            if (logLevel > 2)
                                                printf("Set %d is SOS - priority %d\n", iObj, objSOS->priority());
                                            continue;
                                        }
                                        int iColumn = objects[iObj]->columnNumber();
                                        if (iColumn >= 0) {
                                            if (originalColumns)
                                                iColumn = originalColumns[iColumn];
                                            if (branchDirection) {
                                                OsiSimpleInteger * obj =
                                                    dynamic_cast <OsiSimpleInteger *>(objects[iObj]) ;
                                                if (obj) {
                                                    obj->setPreferredWay(branchDirection[iColumn]);
                                                } else {
                                                    OsiObject2 * obj =
                                                        dynamic_cast <OsiObject2 *>(objects[iObj]) ;
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
                                                printf("Obj %d is int? - priority %d\n", iObj, objects[iObj]->priority());
                                            if (pseudoUp && pseudoUp[iColumn]) {
                                                abort();
                                            }
                                        }
                                    }
                                    // *************************************************************
                                }
                                int statistics = (printOptions > 0) ? printOptions : 0;
#ifdef COIN_HAS_ASL
                                if (!statusUserFunction_[0]) {
#endif
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
#ifdef COIN_HAS_ASL
                                }
#endif
                                if (nodeStrategy) {
                                    // change default
                                    if (nodeStrategy > 2) {
                                        // up or down
                                        int way = (((nodeStrategy - 1) % 1) == 1) ? -1 : +1;
                                        babModel_->setPreferredWay(way);
#ifdef JJF_ZERO
                                        OsiObject ** objects = babModel_->objects();
                                        int numberObjects = babModel_->numberObjects();
                                        for (int iObj = 0; iObj < numberObjects; iObj++) {
                                            CbcObject * obj =
                                                dynamic_cast <CbcObject *>(objects[iObj]) ;
                                            assert (obj);
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
                                        // hybrid was default i.e. mixture of low depth and infeasibility
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
                                    FILE * fp = fopen("user_driver.cpp", "w");
                                    if (fp) {
                                        // generate enough to do BAB
                                        babModel_->generateCpp(fp, 1);
                                        OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (babModel_->solver());
                                        // Make general so do factorization
                                        int factor = osiclp->getModelPtr()->factorizationFrequency();
                                        osiclp->getModelPtr()->setFactorizationFrequency(200);
                                        osiclp->generateCpp(fp);
                                        osiclp->getModelPtr()->setFactorizationFrequency(factor);
                                        //solveOptions.generateCpp(fp);
                                        fclose(fp);
                                        // now call generate code
                                        generateCode(babModel_, "user_driver.cpp", cppValue, prepro);
                                    } else {
                                        std::cout << "Unable to open file user_driver.cpp" << std::endl;
                                    }
                                }
                                if (!babModel_->numberStrong() && babModel_->numberBeforeTrust() > 0)
                                    babModel_->setNumberBeforeTrust(0);
                                if (useStrategy) {
                                    CbcStrategyDefault strategy(1, babModel_->numberStrong(), babModel_->numberBeforeTrust());
                                    strategy.setupPreProcessing(1);
                                    babModel_->setStrategy(strategy);
                                }
                                if (testOsiOptions >= 0) {
                                    sprintf(generalPrint, "Testing OsiObject options %d", testOsiOptions);
                                    generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                    << generalPrint
                                    << CoinMessageEol;
                                    if (!numberSOS) {
                                        babModel_->solver()->findIntegersAndSOS(false);
#ifdef COIN_HAS_LINK
                                        // If linked then pass in model
                                        OsiSolverLink * solver3 = dynamic_cast<OsiSolverLink *> (babModel_->solver());
                                        if (solver3) {
                                            CbcHeuristicDynamic3 serendipity(*babModel_);
                                            serendipity.setHeuristicName("linked");
                                            babModel_->addHeuristic(&serendipity);
                                            double dextra3 = parameters_[whichParam(CBC_PARAM_DBL_DEXTRA3, numberParameters_, parameters_)].doubleValue();
                                            if (dextra3)
                                                solver3->setMeshSizes(dextra3);
                                            int options = parameters_[whichParam(CBC_PARAM_INT_MIPOPTIONS, numberParameters_, parameters_)].intValue() / 10000;
                                            CglStored stored;
                                            if (options) {
                                                printf("nlp options %d\n", options);
                                                /*
                                                  1 - force mini branch and bound
                                                  2 - set priorities high on continuous
                                                  4 - try adding OA cuts
                                                  8 - try doing quadratic linearization
                                                  16 - try expanding knapsacks
                                                              32 - OA cuts strictly concave
                                                  64 - no branching at all on bilinear x-x!
                                                */
                                                if ((options&2)) {
                                                    solver3->setBiLinearPriorities(10, tightenFactor > 0.0 ? tightenFactor : 1.0);
                                                } else if (tightenFactor > 0.0) {
                                                    // set grid size for all continuous bi-linear
                                                    solver3->setMeshSizes(tightenFactor);
                                                }
                                                if ((options&4)) {
                                                    solver3->setSpecialOptions2(solver3->specialOptions2() | (8 + 4));
                                                    // say convex
                                                    solver3->sayConvex((options&32) == 0);
                                                }
                                                int extra1 = parameters_[whichParam(CBC_PARAM_INT_EXTRA1, numberParameters_, parameters_)].intValue();
                                                if ((options&1) != 0 && extra1 > 0)
                                                    solver3->setFixedPriority(extra1);
                                                double cutoff = COIN_DBL_MAX;
                                                if ((options&8))
                                                    cutoff = solver3->linearizedBAB(&stored);
                                                if (cutoff < babModel_->getCutoff()) {
                                                    babModel_->setCutoff(cutoff);
                                                    // and solution
                                                    //babModel_->setBestObjectiveValue(solver3->bestObjectiveValue());
                                                    babModel_->setBestSolution(solver3->bestSolution(), solver3->getNumCols(),
                                                                               solver3->bestObjectiveValue());
                                                }
                                                if ((options&64))
                                                    solver3->setBranchingStrategyOnVariables(16, -1, 4);
                                            }
                                            solver3->setCbcModel(babModel_);
                                            if (stored.sizeRowCuts())
                                                babModel_->addCutGenerator(&stored, 1, "Stored");
                                            CglTemporary temp;
                                            babModel_->addCutGenerator(&temp, 1, "OnceOnly");
                                            //choose.setNumberBeforeTrusted(2000);
                                            //choose.setNumberStrong(20);
                                        }
                                        // For temporary testing of heuristics
                                        //int testOsiOptions = parameters_[whichParam(CBC_PARAM_INT_TESTOSI,numberParameters_,parameters_)].intValue();
                                        if (testOsiOptions >= 10) {
                                            if (testOsiOptions >= 20)
                                                testOsiOptions -= 10;
                                            printf("*** Temp heuristic with mode %d\n", testOsiOptions - 10);
                                            OsiSolverLink * solver3 = dynamic_cast<OsiSolverLink *> (babModel_->solver());
                                            assert (solver3) ;
                                            int extra1 = parameters_[whichParam(CBC_PARAM_INT_EXTRA1, numberParameters_, parameters_)].intValue();
                                            solver3->setBiLinearPriority(extra1);
                                            printf("bilinear priority now %d\n", extra1);
                                            int extra2 = parameters_[whichParam(CBC_PARAM_INT_EXTRA2, numberParameters_, parameters_)].intValue();
                                            double saveDefault = solver3->defaultBound();
                                            solver3->setDefaultBound(static_cast<double> (extra2));
                                            double * solution = solver3->heuristicSolution(slpValue > 0 ? slpValue : 40 , 1.0e-5, testOsiOptions - 10);
                                            solver3->setDefaultBound(saveDefault);
                                            if (!solution)
                                                printf("Heuristic failed\n");
                                        }
#endif
                                    } else {
                                        // move across
                                        babModel_->deleteObjects(false);
                                        //babModel_->addObjects(babModel_->solver()->numberObjects(),babModel_->solver()->objects());
                                    }
                                    CbcBranchDefaultDecision decision;
                                    if (babModel_->numberStrong()) {
                                        OsiChooseStrong choose(babModel_->solver());
                                        choose.setNumberBeforeTrusted(babModel_->numberBeforeTrust());
                                        choose.setNumberStrong(babModel_->numberStrong());
                                        choose.setShadowPriceMode(testOsiOptions);
                                        decision.setChooseMethod(choose);
                                    } else {
                                        OsiChooseVariable choose(babModel_->solver());
                                        decision.setChooseMethod(choose);
                                    }
                                    babModel_->setBranchingMethod(decision);
                                    if (useCosts && testOsiOptions >= 0) {
                                        int numberColumns = babModel_->getNumCols();
                                        int * sort = new int[numberColumns];
                                        double * dsort = new double[numberColumns];
                                        int * priority = new int [numberColumns];
                                        const double * objective = babModel_->getObjCoefficients();
                                        const double * lower = babModel_->getColLower() ;
                                        const double * upper = babModel_->getColUpper() ;
                                        const CoinPackedMatrix * matrix = babModel_->solver()->getMatrixByCol();
                                        const int * columnLength = matrix->getVectorLengths();
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
                                        OsiObject ** objects = babModel_->objects();
                                        int numberObjects = babModel_->numberObjects();
                                        for (int iObj = 0; iObj < numberObjects; iObj++) {
                                            OsiObject * obj = objects[iObj] ;
                                            int iColumn = obj->columnNumber();
                                            if (iColumn >= 0)
                                                obj->setPriority(priority[iColumn]);
                                        }
                                        delete [] priority;
                                        delete [] sort;
                                        delete [] dsort;
                                    }
                                }
                                checkSOS(babModel_, babModel_->solver());
                                if (doSprint > 0) {
                                    // Sprint for primal solves
                                    ClpSolve::SolveType method = ClpSolve::usePrimalorSprint;
                                    ClpSolve::PresolveType presolveType = ClpSolve::presolveOff;
                                    int numberPasses = 5;
                                    int options[] = {0, 3, 0, 0, 0, 0};
                                    int extraInfo[] = { -1, 20, -1, -1, -1, -1};
                                    extraInfo[1] = doSprint;
                                    int independentOptions[] = {0, 0, 3};
                                    ClpSolve clpSolve(method, presolveType, numberPasses,
                                                      options, extraInfo, independentOptions);
                                    // say use in OsiClp
                                    clpSolve.setSpecialOption(6, 1);
                                    OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (babModel_->solver());
                                    osiclp->setSolveOptions(clpSolve);
                                    osiclp->setHintParam(OsiDoDualInResolve, false);
                                    // switch off row copy
                                    osiclp->getModelPtr()->setSpecialOptions(osiclp->getModelPtr()->specialOptions() | 256);
                                    osiclp->getModelPtr()->setInfeasibilityCost(1.0e11);
                                }
#ifdef COIN_HAS_LINK
                                if (storedAmpl.sizeRowCuts()) {
                                    if (preProcess) {
                                        const int * originalColumns = process.originalColumns();
                                        int numberColumns = babModel_->getNumCols();
                                        int * newColumn = new int[numberOriginalColumns];
                                        int i;
                                        for (i = 0; i < numberOriginalColumns; i++)
                                            newColumn[i] = -1;
                                        for (i = 0; i < numberColumns; i++) {
                                            int iColumn = originalColumns[i];
                                            newColumn[iColumn] = i;
                                        }
                                        int * buildColumn = new int[numberColumns];
                                        // Build up valid cuts
                                        int nBad = 0;
                                        int nCuts = storedAmpl.sizeRowCuts();
                                        CglStored newCuts;
                                        for (i = 0; i < nCuts; i++) {
                                            const OsiRowCut * cut = storedAmpl.rowCutPointer(i);
                                            double lb = cut->lb();
                                            double ub = cut->ub();
                                            int n = cut->row().getNumElements();
                                            const int * column = cut->row().getIndices();
                                            const double * element = cut->row().getElements();
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
                                        delete [] newColumn;
                                        delete [] buildColumn;
                                    }
                                }
#endif
#ifdef CLP_MALLOC_STATISTICS
                                malloc_stats();
                                malloc_stats2();
#endif
#ifndef CBC_OTHER_SOLVER
                                if (outputFormat == 5) {
                                    osiclp = dynamic_cast< OsiClpSolverInterface*> (babModel_->solver());
                                    lpSolver = osiclp->getModelPtr();
                                    lpSolver->setPersistenceFlag(1);
                                }
#endif
#ifdef COIN_HAS_ASL
                                // add in lotsizing
                                if (statusUserFunction_[0] && info.special) {
                                    int numberColumns = babModel_->getNumCols();
                                    int i;
                                    int n = 0;
                                    if (preProcess) {
                                        const int * originalColumns = process.originalColumns();
                                        for (i = 0; i < numberColumns; i++) {
                                            int iColumn = originalColumns[i];
                                            assert (iColumn >= i);
                                            int iType = info.special[iColumn];
                                            if (iType) {
                                                assert (iType == 1);
                                                n++;
                                            }
                                            info.special[i] = iType;
                                        }
                                    }
                                    if (n) {
                                        int numberIntegers = 0;
                                        int numberOldObjects = 0;
                                        OsiObject ** oldObjects = NULL;
                                        const double * lower = babModel_->solver()->getColLower();
                                        const double * upper = babModel_->solver()->getColUpper();
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
                                        OsiObject ** objects = new OsiObject * [n];
                                        n = 0;
                                        // set new objects to have one lower priority
                                        double ranges[] = { -COIN_DBL_MAX, -1.0, 1.0, COIN_DBL_MAX};
                                        for (int iObj = 0; iObj < numberOldObjects; iObj++) {
                                            int iColumn = oldObjects[iObj]->columnNumber();
                                            if (iColumn >= 0 && info.special[iColumn]) {
                                                if (lower[iColumn] <= -1.0 && upper[iColumn] >= 0.0) {
                                                    ranges[0] = lower[iColumn];
                                                    ranges[3] = upper[iColumn];
                                                    int priority = oldObjects[iObj]->priority();
                                                    if (testOsiOptions < 0) {
                                                        objects[n] = new CbcLotsize(babModel_, iColumn, 2, ranges, true);
                                                    } else {
                                                        objects[n] = new OsiLotsize(testOsiSolver, iColumn, 2, ranges, true);
                                                    }
                                                    objects[n++]->setPriority (priority - 1);
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
                                        delete [] objects;
                                    }
                                }
#endif
                                if (storedAmpl.sizeRowCuts()) {
                                    //babModel_->addCutGenerator(&storedAmpl,1,"AmplStored");
                                    int numberRowCuts = storedAmpl.sizeRowCuts();
                                    for (int i = 0; i < numberRowCuts; i++) {
                                        const OsiRowCut * rowCutPointer = storedAmpl.rowCutPointer(i);
                                        babModel_->makeGlobalCut(rowCutPointer);
                                    }
                                }
                                // If defaults then increase trust for small models
                                if (!strongChanged) {
                                    int numberColumns = babModel_->getNumCols();
                                    if (numberColumns <= 50)
                                        babModel_->setNumberBeforeTrust(1000);
                                    else if (numberColumns <= 100)
                                        babModel_->setNumberBeforeTrust(100);
                                    else if (numberColumns <= 300)
                                        babModel_->setNumberBeforeTrust(50);
                                }
#ifdef CBC_THREAD
                                int numberThreads = parameters_[whichParam(CBC_PARAM_INT_THREADS, numberParameters_, parameters_)].intValue();
                                babModel_->setNumberThreads(numberThreads % 100);
                                babModel_->setThreadMode(numberThreads / 100);
#endif
                                int returnCode = callBack(babModel_, 3);
                                if (returnCode) {
                                    // exit if user wants
                                    delete babModel_;
                                    babModel_ = NULL;
                                    return returnCode;
                                }
#ifndef CBC_OTHER_SOLVER
                                osiclp = dynamic_cast< OsiClpSolverInterface*> (babModel_->solver());
                                lpSolver = osiclp->getModelPtr();
#elif CBC_OTHER_SOLVER==1
#endif
                                if ((experimentFlag >= 1 || strategyFlag >= 1) && babModel_->fastNodeDepth() == -1) {
                                    if (babModel_->solver()->getNumCols() +
                                            babModel_->solver()->getNumRows() < 500)
                                        babModel_->setFastNodeDepth(-12);
                                } else if (babModel_->fastNodeDepth() == -999) {
                                    babModel_->setFastNodeDepth(-1);
                                }
                                int heurOptions = parameters_[whichParam(CBC_PARAM_INT_HOPTIONS, numberParameters_, parameters_)].intValue();
                                if (heurOptions > 100)
                                    babModel_->setSpecialOptions(babModel_->specialOptions() | 8192);

#ifndef CBC_OTHER_SOLVER
#ifdef CLP_MULTIPLE_FACTORIZATIONS
                                int denseCode = parameters_[whichParam(CBC_PARAM_INT_DENSE, numberParameters_, parameters_)].intValue();
                                int smallCode = parameters_[whichParam(CBC_PARAM_INT_SMALLFACT, numberParameters_, parameters_)].intValue();
                                if (bothFlags >= 1) {
                                    if (denseCode < 0)
                                        denseCode = 40;
                                    if (smallCode < 0 && !lpSolver->factorization()->isDenseOrSmall())
                                        smallCode = 40;
                                }
                                if (denseCode > 0) {
                                    lpSolver->factorization()->setGoDenseThreshold(denseCode);
                                    assert (osiclp == babModel_->solver());
                                    osiclp->setSpecialOptions(osiclp->specialOptions() | 1024);
                                }
                                if (smallCode > 0 && smallCode > denseCode)
                                    lpSolver->factorization()->setGoSmallThreshold(smallCode);
                                //if (denseCode>=lpSolver->numberRows()) {
                                //lpSolver->factorization()->goDense();
                                //}
                                if (lpSolver->factorization()->goOslThreshold() > 1000) {
                                    // use osl in gomory (may not if CglGomory decides not to)
                                    int numberGenerators = babModel_->numberCutGenerators();
				    int nGomory=0;
                                    for (int iGenerator = 0; iGenerator < numberGenerators;
                                            iGenerator++) {
                                        CbcCutGenerator * generator = babModel_->cutGenerator(iGenerator);
                                        CglGomory * gomory = dynamic_cast<CglGomory *>
                                                             (generator->generator());
                                        if (gomory) {
					  if (nGomory<2) {
                                            gomory->useAlternativeFactorization();
					  } else if (gomory->originalSolver()) {
					    OsiClpSolverInterface * clpSolver = dynamic_cast<OsiClpSolverInterface *>(gomory->originalSolver());
					    if (clpSolver) {
					      ClpSimplex * simplex = clpSolver->getModelPtr();
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
                                if (useSolution > 1) {
                                    // use hotstart to try and find solution
                                    CbcHeuristicPartial partial(*babModel_, 10000, useSolution);
                                    partial.setHeuristicName("Partial solution given");
                                    babModel_->addHeuristic(&partial);
                                }
                                if (logLevel <= 1)
                                    babModel_->solver()->setHintParam(OsiDoReducePrint, true, OsiHintTry);
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
                                    const double * solution = osiclp->getColSolution();
                                    const int * originalColumns = process.originalColumns();
                                    for (int i = 0; i < numberColumns; i++) {
                                        int iColumn = originalColumns[i];
                                        if (saveSolver->isInteger(iColumn)) {
                                            double value = solution[i];
                                            double value2 = floor(value + 0.5);
                                            assert (fabs(value - value2) < 1.0e-3);
                                            saveSolver->setColLower(iColumn, value2);
                                            saveSolver->setColUpper(iColumn, value2);
                                        }
                                    }
                                    saveSolver->writeMps("fixed");
                                    babModel_->setBestSolution(osiclp->getColSolution(),
                                                               osiclp->getNumCols(),
                                                               1.5325e10);
                                } else {
                                    babModel_->branchAndBound(statistics);
                                }
#else
#ifdef ORBITAL
                                CbcOrbital orbit(babModel_);
                                orbit.morph();
                                exit(1);
#endif
                                int hOp1 = parameters_[whichParam(CBC_PARAM_INT_HOPTIONS, numberParameters_, parameters_)].intValue() / 100000;
                                if (hOp1 % 10) {
                                    CbcCompareDefault compare;
                                    compare.setBreadthDepth(hOp1 % 10);
                                    babModel_->setNodeComparison(compare);
                                }
#if CBC_OTHER_SOLVER==1
                                if (dynamic_cast<OsiCpxSolverInterface *> (babModel_->solver()))
                                    babModel_->solver()->messageHandler()->setLogLevel(0);
#endif
                                if (parameters_[whichParam(CBC_PARAM_STR_CPX, numberParameters_, parameters_)].currentOptionAsInteger()) {
                                    babModel_->setSpecialOptions(babModel_->specialOptions() | 16384);
                                    //if (babModel_->fastNodeDepth()==-1)
                                    babModel_->setFastNodeDepth(-2); // Use Cplex at root
                                }
                                int hOp2 = parameters_[whichParam(CBC_PARAM_INT_HOPTIONS, numberParameters_, parameters_)].intValue() / 10000;
                                if (hOp2 % 10) {
                                    babModel_->setSpecialOptions(babModel_->specialOptions() | 16384);
                                    if (babModel_->fastNodeDepth() == -1)
                                        babModel_->setFastNodeDepth(-2); // Use Cplex at root
                                }
                                if (experimentFlag >= 5) {
                                    CbcModel donor(*babModel_);
                                    int options = babModel_->specialOptions();
                                    donor.setSpecialOptions(options | 262144);
                                    ClpSimplex * lpSolver2;
                                    OsiClpSolverInterface * clpSolver2;
                                    clpSolver2 =
                                        dynamic_cast<OsiClpSolverInterface *> (donor.solver());
                                    assert (clpSolver2);
                                    lpSolver2 = clpSolver2->getModelPtr();
                                    assert (lpSolver2);
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
                                            CbcCutGenerator * generator = donor.cutGenerator(iGenerator);
                                            CglGomory * gomory = dynamic_cast<CglGomory *>
                                                                 (generator->generator());
                                            if (gomory)
                                                gomory->useAlternativeFactorization(false);
                                        }
                                    } else {
                                        printf("code this\n");
                                        abort();
                                    }
                                    babModel_->setSpecialOptions(options | 524288);
                                    CglStored * stored = new CglStored(donor.getNumCols());
                                    donor.setStoredRowCuts(stored);
                                    donor.branchAndBound(0);
                                    babModel_->setStoredRowCuts(donor.storedRowCuts());
                                    donor.setStoredRowCuts(NULL);
                                }
#ifdef JJF_ZERO
                                int extra5 = parameters_[whichParam(EXTRA5, numberParameters_, parameters_)].intValue();
                                if (extra5 > 0) {
                                    int numberGenerators = babModel_->numberCutGenerators();
                                    for (int iGenerator = 0; iGenerator < numberGenerators;
                                            iGenerator++) {
                                        CbcCutGenerator * generator = babModel_->cutGenerator(iGenerator);
                                        CglGomory * gomory = dynamic_cast<CglGomory *>
                                                             (generator->generator());
                                        if (gomory) {
                                            CglGomory gomory2(*gomory);
                                            gomory2.useAlternativeFactorization(!gomory->alternativeFactorization());
                                            babModel_->addCutGenerator(&gomory2, -99, "Gomory2");
                                        }
                                    }
                                }
#endif
                                babModel_->branchAndBound(statistics);
                                //#define CLP_FACTORIZATION_INSTRUMENT
#ifdef CLP_FACTORIZATION_INSTRUMENT
                                extern double factorization_instrument(int type);
                                double facTime = factorization_instrument(0);
                                printf("Factorization %g seconds\n",
                                       facTime);
#endif
#endif
#ifdef COIN_DEVELOP
#ifndef JJF_ONE
                                {
                                    int numberColumns = babModel_->getNumCols();
                                    const double * solution = babModel_->bestSolution();
                                    if (solution && numberColumns < 1000) {
                                        for (int i = 0; i < numberColumns; i++) {
                                            if (solution[i])
                                                printf("SOL %d %.18g\n", i, solution[i]);
                                        }
                                    }
                                }
#endif
                                void printHistory(const char * file/*,CbcModel * model*/);
                                printHistory("branch.log"/*,babModel_*/);
#endif
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
                            } else if (type == CBC_PARAM_ACTION_MIPLIB) {
                                int typeOfCuts = babModel_->numberCutGenerators() ? 1 : -1;
                                CbcStrategyDefault strategy(typeOfCuts,
                                                            babModel_->numberStrong(),
                                                            babModel_->numberBeforeTrust());
                                // Set up pre-processing
                                int translate2[] = {9999, 1, 1, 3, 2, 4, 5, 6, 6};
                                if (preProcess)
                                    strategy.setupPreProcessing(translate2[ preProcess ]);
                                babModel_->setStrategy(strategy);
#ifdef CBC_THREAD
                                int numberThreads = parameters_[whichParam(CBC_PARAM_INT_THREADS, numberParameters_, parameters_)].intValue();
                                babModel_->setNumberThreads(numberThreads % 100);
                                babModel_->setThreadMode(numberThreads / 100);
#endif
#ifndef CBC_OTHER_SOLVER
                                if (outputFormat == 5) {
                                    osiclp = dynamic_cast< OsiClpSolverInterface*> (babModel_->solver());
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
                                    //babModel_->deleteObjects(false);
                                    decision.setChooseMethod(choose);
                                    babModel_->setBranchingMethod(decision);
                                }
                                model_ = *babModel_;
#ifndef CBC_OTHER_SOLVER
                                {
                                    osiclp = dynamic_cast< OsiClpSolverInterface*> (model_.solver());
                                    lpSolver = osiclp->getModelPtr();
                                    lpSolver->setSpecialOptions(lpSolver->specialOptions() | IN_BRANCH_AND_BOUND); // say is Cbc (and in branch and bound)
                                    if (lpSolver->factorization()->goOslThreshold() > 1000) {
                                        // use osl in gomory (may not if CglGomory decides not to)
                                        int numberGenerators = model_.numberCutGenerators();
                                        for (int iGenerator = 0; iGenerator < numberGenerators;
                                                iGenerator++) {
                                            CbcCutGenerator * generator = model_.cutGenerator(iGenerator);
                                            CglGomory * gomory = dynamic_cast<CglGomory *>
                                                                 (generator->generator());
                                            if (gomory)
                                                gomory->useAlternativeFactorization();
                                        }
                                    }
                                }
#endif
                                /* LL: this was done in CoinSolve.cpp: main(argc, argv).
                                   I have moved it here so that the miplib directory location
                                   could be passed to CbcClpUnitTest. */
                                /* JJF: No need to have 777 flag at all - user
                                   says -miplib
                                */
                                int extra2 = parameters_[whichParam(CBC_PARAM_INT_EXTRA2, numberParameters_, parameters_)].intValue();
                                double stuff[11];
                                stuff[0] = parameters_[whichParam(CBC_PARAM_DBL_FAKEINCREMENT, numberParameters_, parameters_)].doubleValue();
                                stuff[1] = parameters_[whichParam(CBC_PARAM_DBL_FAKECUTOFF, numberParameters_, parameters_)].doubleValue();
                                stuff[2] = parameters_[whichParam(CBC_PARAM_DBL_DEXTRA3, numberParameters_, parameters_)].doubleValue();
                                stuff[3] = parameters_[whichParam(CBC_PARAM_DBL_DEXTRA4, numberParameters_, parameters_)].doubleValue();
                                stuff[4] = parameters_[whichParam(CBC_PARAM_INT_DENSE, numberParameters_, parameters_)].intValue();
                                stuff[5] = parameters_[whichParam(CBC_PARAM_INT_EXTRA1, numberParameters_, parameters_)].intValue();
                                stuff[6] = parameters_[whichParam(CBC_PARAM_INT_EXTRA3, numberParameters_, parameters_)].intValue();
                                stuff[7] = parameters_[whichParam(CBC_PARAM_INT_DEPTHMINIBAB, numberParameters_, parameters_)].intValue();
                                stuff[8] = bothFlags;
                                stuff[9] = doVector;
                                stuff[10] = parameters_[whichParam(CBC_PARAM_INT_SMALLFACT, numberParameters_, parameters_)].intValue();
                                if ( dominatedCuts)
                                    model_.setSpecialOptions(model_.specialOptions() | 64);
                                if (parameters_[whichParam(CBC_PARAM_STR_CPX, numberParameters_, parameters_)].currentOptionAsInteger()) {
                                    model_.setSpecialOptions(model_.specialOptions() | 16384);
                                    //if (model_.fastNodeDepth()==-1)
                                    model_.setFastNodeDepth(-2); // Use Cplex at root
                                }
                                int hOp2 = parameters_[whichParam(CBC_PARAM_INT_HOPTIONS, numberParameters_, parameters_)].intValue() / 10000;
                                if (hOp2 % 10) {
                                    model_.setSpecialOptions(model_.specialOptions() | 16384);
                                    if (model_.fastNodeDepth() == -1)
                                        model_.setFastNodeDepth(-2); // Use Cplex at root
                                }
                                if (!pumpChanged) {
                                    // Make more lightweight
                                    for (int iHeur = 0; iHeur < model_.numberHeuristics(); iHeur++) {
                                        CbcHeuristic * heuristic = model_.heuristic(iHeur);
                                        CbcHeuristicFPump* pump =
                                            dynamic_cast<CbcHeuristicFPump*>(heuristic);
                                        if (pump) {
                                            CbcHeuristicFPump heuristic4(model_);
                                            heuristic4.setFractionSmall(0.5);
                                            heuristic4.setMaximumPasses(5);
                                            heuristic4.setFeasibilityPumpOptions(30);
                                            heuristic4.setWhen(13);
                                            heuristic4.setHeuristicName("feasibility pump");
                                            //CbcHeuristicFPump & pump2 = pump;
                                            *pump = heuristic4;
                                        }
                                    }
                                }
                                int returnCode = CbcClpUnitTest(model_, dirMiplib, extra2, stuff);
                                babModel_ = NULL;
                                return returnCode;
                            } else {
                                abort(); // can't get here
                                //strengthenedModel = babModel_->strengthenedModel();
                            }
                            currentBranchModel = NULL;
#ifndef CBC_OTHER_SOLVER
                            osiclp = dynamic_cast< OsiClpSolverInterface*> (babModel_->solver());
                            if (debugFile == "createAfterPre" && babModel_->bestSolution()) {
                                lpSolver = osiclp->getModelPtr();
                                //move best solution (should be there -- but ..)
                                int n = lpSolver->getNumCols();
                                memcpy(lpSolver->primalColumnSolution(), babModel_->bestSolution(), n*sizeof(double));
                                saveSolution(osiclp->getModelPtr(), "debug.file");
                            }
#endif
                            statistics_cut_time = 0.0;
                            if (!noPrinting_) {
                                // Print more statistics
                                sprintf(generalPrint, "Cuts at root node changed objective from %g to %g",
                                        babModel_->getContinuousObjective(), babModel_->rootObjectiveAfterCuts());
                                generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                << generalPrint
                                << CoinMessageEol;

                                numberGenerators = babModel_->numberCutGenerators();
                                statistics_number_cuts = new int [numberGenerators];;
                                statistics_number_generators = numberGenerators;
                                statistics_name_generators = new const char *[numberGenerators];
                                char timing[30];
                                for (iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
                                    CbcCutGenerator * generator = babModel_->cutGenerator(iGenerator);
                                    statistics_name_generators[iGenerator] =
                                        generator->cutGeneratorName();
                                    statistics_number_cuts[iGenerator] = generator->numberCutsInTotal();
                                    sprintf(generalPrint, "%s was tried %d times and created %d cuts of which %d were active after adding rounds of cuts",
                                            generator->cutGeneratorName(),
                                            generator->numberTimesEntered(),
                                            generator->numberCutsInTotal() +
                                            generator->numberColumnCuts(),
                                            generator->numberCutsActive());
                                    if (generator->timing()) {
                                        sprintf(timing, " (%.3f seconds)", generator->timeInCutGenerator());
                                        strcat(generalPrint, timing);
                                        statistics_cut_time += generator->timeInCutGenerator();
                                    }
                                    CglStored * stored = dynamic_cast<CglStored*>(generator->generator());
                                    if (stored && !generator->numberCutsInTotal())
                                        continue;
#ifndef CLP_INVESTIGATE
                                    CglImplication * implication = dynamic_cast<CglImplication*>(generator->generator());
                                    if (implication)
                                        continue;
#endif
                                    generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                    << generalPrint
                                    << CoinMessageEol;
                                }
#ifdef COIN_DEVELOP
                                printf("%d solutions found by heuristics\n",
                                       babModel_->getNumberHeuristicSolutions());
                                // Not really generator but I am feeling lazy
                                for (iGenerator = 0; iGenerator < babModel_->numberHeuristics(); iGenerator++) {
                                    CbcHeuristic * heuristic = babModel_->heuristic(iGenerator);
                                    if (heuristic->numRuns()) {
                                        // Need to bring others inline
                                        sprintf(generalPrint, "%s was tried %d times out of %d and created %d solutions\n",
                                                heuristic->heuristicName(),
                                                heuristic->numRuns(),
                                                heuristic->numCouldRun(),
                                                heuristic->numberSolutionsFound());
                                        generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                        << generalPrint
                                        << CoinMessageEol;
                                    }
                                }
#endif
                            }
                            // adjust time to allow for children on some systems
                            time2 = CoinCpuTime() + CoinCpuTimeJustChildren();
                            totalTime += time2 - time1;
                            // For best solution
                            double * bestSolution = NULL;
                            // Say in integer
                            if (babModel_->status()) {
                                // treat as stopped
                                integerStatus = 3;
                            } else {
                                if (babModel_->isProvenOptimal()) {
                                    integerStatus = 0;
                                } else {
                                    // infeasible
                                    integerStatus = 6;
                                }
                            }
                            if (babModel_->getMinimizationObjValue() < 1.0e50 && type == CBC_PARAM_ACTION_BAB) {
                                // post process
                                int n;
                                if (preProcess) {
                                    n = saveSolver->getNumCols();
                                    bestSolution = new double [n];
#ifndef CBC_OTHER_SOLVER
                                    OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (babModel_->solver());
#else
                                    OsiCpxSolverInterface * clpSolver = dynamic_cast< OsiCpxSolverInterface*> (babModel_->solver());
#endif
                                    // Save bounds on processed model
                                    const int * originalColumns = process.originalColumns();
                                    int numberColumns2 = clpSolver->getNumCols();
                                    double * solution2 = new double[n];
                                    double * lower2 = new double [n];
                                    double * upper2 = new double [n];
                                    for (int i = 0; i < n; i++) {
                                        solution2[i] = COIN_DBL_MAX;
                                        lower2[i] = COIN_DBL_MAX;
                                        upper2[i] = -COIN_DBL_MAX;
                                    }
                                    const double *columnLower = clpSolver->getColLower() ;
                                    const double * columnUpper = clpSolver->getColUpper() ;
                                    const double * solution = babModel_->bestSolution();
                                    for (int i = 0; i < numberColumns2; i++) {
                                        int jColumn = originalColumns[i];
                                        if (jColumn < n) {
                                            solution2[jColumn] = solution[i];
                                            lower2[jColumn] = columnLower[i];
                                            upper2[jColumn] = columnUpper[i];
                                        }
                                    }
#ifndef CBC_OTHER_SOLVER
                                    ClpSimplex * lpSolver = clpSolver->getModelPtr();
                                    lpSolver->setSpecialOptions(lpSolver->specialOptions() | IN_BRANCH_AND_BOUND); // say is Cbc (and in branch and bound)
#endif
                                    process.postProcess(*babModel_->solver());
#ifdef COIN_DEVELOP
                                    if (model_.bestSolution() && fabs(model_.getMinimizationObjValue() -
                                                                      babModel_->getMinimizationObjValue()) < 1.0e-8) {
                                        const double * b1 = model_.bestSolution();
                                        const double * b2 = saveSolver->getColSolution();
                                        const double * columnLower = saveSolver->getColLower() ;
                                        const double * columnUpper = saveSolver->getColUpper() ;
                                        for (int i = 0; i < n; i++) {
                                            if (fabs(b1[i] - b2[i]) > 1.0e-7) {
                                                printf("%d %g %g %g %g\n", i, b1[i], b2[i],
                                                       columnLower[i], columnUpper[i]);
                                            }
                                        }
                                    }
#endif
                                    bool tightenB = false;
                                    {
                                        int n = babModel_->numberObjects();
                                        for (int i = 0; i < n; i++) {
                                            const OsiObject * obj = babModel_->object(i);
                                            if (!dynamic_cast<const CbcSimpleInteger *>(obj)) {
                                                tightenB = true;
                                                break;
                                            }
                                        }
                                    }
                                    // Solution now back in saveSolver
                                    // Double check bounds
                                    columnLower = saveSolver->getColLower() ;
                                    columnUpper = saveSolver->getColUpper() ;
                                    solution = saveSolver->getColSolution();
                                    int numberChanged = 0;
                                    for (int i = 0; i < n; i++) {
                                        if (!saveSolver->isInteger(i) && !tightenB)
                                            continue;
                                        if (lower2[i] != COIN_DBL_MAX) {
                                            if (lower2[i] != columnLower[i] ||
                                                    upper2[i] != columnUpper[i]) {
                                                if (lower2[i] < columnLower[i] ||
                                                        upper2[i] > columnUpper[i]) {
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
                                    OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (saveSolver);
                                    if (osiclp && osiclp->numberSOS()) {
                                        // SOS
                                        numberSOS = osiclp->numberSOS();
                                        const CoinSet * setInfo = osiclp->setInfo();
                                        int i;
                                        for ( i = 0; i < numberSOS; i++) {
                                            int type = setInfo[i].setType();
                                            int n = setInfo[i].numberEntries();
                                            const int * which = setInfo[i].which();
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
                                            assert (last - first < type);
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
                                    delete [] solution2;
                                    delete [] lower2;
                                    delete [] upper2;
                                    if (numberChanged) {
                                        sprintf(generalPrint, "%d bounds tightened after postprocessing\n",
                                                numberChanged);
                                        generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                        << generalPrint
                                        << CoinMessageEol;
                                    }
                                    saveSolver->resolve();
                                    if (!saveSolver->isProvenOptimal()) {
                                        // try all slack
                                        CoinWarmStartBasis * basis = dynamic_cast<CoinWarmStartBasis *> (babModel_->solver()->getEmptyWarmStart());
                                        saveSolver->setWarmStart(basis);
                                        delete basis;
                                        saveSolver->initialSolve();
#ifdef COIN_DEVELOP
                                        saveSolver->writeMps("inf2");
#endif
					OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (saveSolver);
					if (osiclp)
					  osiclp->getModelPtr()->checkUnscaledSolution();
                                    }
                                    assert (saveSolver->isProvenOptimal());
#ifndef CBC_OTHER_SOLVER
                                    // and original solver
                                    originalSolver->setDblParam(OsiDualObjectiveLimit, COIN_DBL_MAX);
                                    assert (n >= originalSolver->getNumCols());
                                    n = originalSolver->getNumCols();
                                    originalSolver->setColLower(saveSolver->getColLower());
                                    originalSolver->setColUpper(saveSolver->getColUpper());
                                    // basis
                                    CoinWarmStartBasis * basis = dynamic_cast<CoinWarmStartBasis *> (babModel_->solver()->getWarmStart());
                                    originalSolver->setBasis(*basis);
                                    delete basis;
                                    originalSolver->resolve();
                                    if (!originalSolver->isProvenOptimal()) {
                                        // try all slack
                                        CoinWarmStartBasis * basis = dynamic_cast<CoinWarmStartBasis *> (babModel_->solver()->getEmptyWarmStart());
                                        originalSolver->setBasis(*basis);
                                        delete basis;
                                        originalSolver->initialSolve();
					OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (originalSolver);
					if (osiclp)
					  osiclp->getModelPtr()->checkUnscaledSolution();
                                    }
                                    assert (originalSolver->isProvenOptimal());
#endif
                                    babModel_->assignSolver(saveSolver);
                                    memcpy(bestSolution, babModel_->solver()->getColSolution(), n*sizeof(double));
                                } else {
                                    n = babModel_->solver()->getNumCols();
                                    bestSolution = new double [n];
                                    memcpy(bestSolution, babModel_->solver()->getColSolution(), n*sizeof(double));
                                }
                                if (returnMode == 1) {
                                    model_.deleteSolutions();
                                    model_.setBestSolution(bestSolution, n, babModel_->getMinimizationObjValue());
                                }
                                babModel_->deleteSolutions();
                                babModel_->setBestSolution(bestSolution, n, babModel_->getMinimizationObjValue());
#ifndef CBC_OTHER_SOLVER
                                // and put back in very original solver
                                {
                                    ClpSimplex * original = originalSolver->getModelPtr();
                                    double * lower = original->columnLower();
                                    double * upper = original->columnUpper();
                                    double * solution = original->primalColumnSolution();
                                    int n = original->numberColumns();
                                    //assert (!n||n==babModel_->solver()->getNumCols());
                                    for (int i = 0; i < n; i++) {
                                        solution[i] = bestSolution[i];
                                        if (originalSolver->isInteger(i)) {
                                            lower[i] = solution[i];
                                            upper[i] = solution[i];
                                        }
                                    }
                                    // basis
                                    CoinWarmStartBasis * basis = dynamic_cast<CoinWarmStartBasis *> (babModel_->solver()->getWarmStart());
                                    originalSolver->setBasis(*basis);
                                    delete basis;
                                    originalSolver->setDblParam(OsiDualObjectiveLimit, COIN_DBL_MAX);
                                    originalSolver->resolve();
                                    if (!originalSolver->isProvenOptimal()) {
                                        // try all slack
                                        CoinWarmStartBasis * basis = dynamic_cast<CoinWarmStartBasis *> (babModel_->solver()->getEmptyWarmStart());
                                        originalSolver->setBasis(*basis);
                                        delete basis;
                                        originalSolver->initialSolve();
					OsiClpSolverInterface * osiclp = dynamic_cast< OsiClpSolverInterface*> (originalSolver);
					if (osiclp)
					  osiclp->getModelPtr()->checkUnscaledSolution();
#ifdef CLP_INVESTIGATE
                                        if (!originalSolver->isProvenOptimal()) {
                                            if (saveSolver) {
                                                printf("saveSolver and originalSolver matrices saved\n");
                                                saveSolver->writeMps("infA");
                                            } else {
                                                printf("originalSolver matrix saved\n");
                                                originalSolver->writeMps("infB");
                                            }
                                        }
#endif
                                    }
                                    assert (originalSolver->isProvenOptimal());
                                }
#endif
                                checkSOS(babModel_, babModel_->solver());
                            } else if (model_.bestSolution() && type == CBC_PARAM_ACTION_BAB && model_.getMinimizationObjValue() < 1.0e50 && preProcess) {
                                sprintf(generalPrint, "Restoring heuristic best solution of %g", model_.getMinimizationObjValue());
                                generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                << generalPrint
                                << CoinMessageEol;
                                int n = saveSolver->getNumCols();
                                bestSolution = new double [n];
                                // Put solution now back in saveSolver
                                saveSolver->setColSolution(model_.bestSolution());
                                babModel_->assignSolver(saveSolver);
				saveSolver=NULL;
                                babModel_->setMinimizationObjValue(model_.getMinimizationObjValue());
                                memcpy(bestSolution, babModel_->solver()->getColSolution(), n*sizeof(double));
#ifndef CBC_OTHER_SOLVER
                                // and put back in very original solver
                                {
                                    ClpSimplex * original = originalSolver->getModelPtr();
                                    double * lower = original->columnLower();
                                    double * upper = original->columnUpper();
                                    double * solution = original->primalColumnSolution();
                                    int n = original->numberColumns();
                                    //assert (!n||n==babModel_->solver()->getNumCols());
                                    for (int i = 0; i < n; i++) {
                                        solution[i] = bestSolution[i];
                                        if (originalSolver->isInteger(i)) {
                                            lower[i] = solution[i];
                                            upper[i] = solution[i];
                                        }
                                    }
                                    // basis
                                    CoinWarmStartBasis * basis = dynamic_cast<CoinWarmStartBasis *> (babModel_->solver()->getWarmStart());
                                    originalSolver->setBasis(*basis);
                                    delete basis;
                                }
#endif
                            }
#ifndef CBC_OTHER_SOLVER
                            //if (type==CBC_PARAM_ACTION_STRENGTHEN&&strengthenedModel)
                            //clpSolver = dynamic_cast< OsiClpSolverInterface*> (strengthenedModel);
#ifdef COIN_HAS_ASL
                            else if (statusUserFunction_[0])
                                clpSolver = dynamic_cast< OsiClpSolverInterface*> (babModel_->solver());
#endif
                            lpSolver = clpSolver->getModelPtr();
                            if (numberChanged) {
                                for (int i = 0; i < numberChanged; i++) {
                                    int iColumn = changed[i];
                                    clpSolver->setContinuous(iColumn);
                                }
                                delete [] changed;
                            }
#endif
                            if (type == CBC_PARAM_ACTION_BAB) {
#ifndef CBC_OTHER_SOLVER
                                //move best solution (should be there -- but ..)
                                int n = lpSolver->getNumCols();
                                if (bestSolution) {
                                    memcpy(lpSolver->primalColumnSolution(), bestSolution, n*sizeof(double));
                                    // now see what that does to row solution
                                    int numberRows = lpSolver->numberRows();
                                    double * rowSolution = lpSolver->primalRowSolution();
                                    memset (rowSolution, 0, numberRows*sizeof(double));
                                    lpSolver->clpMatrix()->times(1.0, bestSolution, rowSolution);
                                    lpSolver->setObjectiveValue(babModel_->getObjValue());
                                }
                                if (debugFile == "create" && bestSolution) {
                                    saveSolution(lpSolver, "debug.file");
                                }
#else
                                if (bestSolution) {
                                    model_.solver()->setColSolution(bestSolution);
                                }
#endif
				delete saveSolver;
                                delete [] bestSolution;
				std::string statusName[] = {"", "Stopped on ", "Run abandoned", "", "", "User ctrl-c"};
				std::string minor[] = {"Optimal solution found", "Linear relaxation infeasible", "Optimal solution found (within gap tolerance)", "node limit", "time limit", "user ctrl-c", "solution limit", "Linear relaxation unbounded", "Problem proven infeasible"};
                                int iStat = babModel_->status();
                                int iStat2 = babModel_->secondaryStatus();
                                if (!iStat && !iStat2 && !bestSolution)
                                    iStat2 = 8;
                                statistics_seconds = time2 - time1;
                                statistics_sys_seconds = CoinSysTime();
                                statistics_elapsed_seconds = CoinWallclockTime();
                                statistics_obj = babModel_->getObjValue();
                                statistics_continuous = babModel_->getContinuousObjective();
                                statistics_tighter = babModel_->rootObjectiveAfterCuts();
                                statistics_nodes = babModel_->getNodeCount();
                                statistics_iterations = babModel_->getIterationCount();;
                                statistics_result = statusName[iStat];;
                                if (!noPrinting_) {
				    sprintf(generalPrint, "\nResult - %s%s\n", 
					    statusName[iStat].c_str(), 
					    minor[iStat2].c_str());
				    generalMessageHandler->message(CLP_GENERAL, generalMessages)
				       << generalPrint
				       << CoinMessageEol;
				    if (babModel_->bestSolution()){
				      sprintf(generalPrint, 
					      "Objective value:                %.8f\n", 
					      babModel_->getObjValue());
				    }else{
				      sprintf(generalPrint,
					      "No feasible solution found\n");
				    }
				    if (iStat2 >= 2 && iStat2 <=6){
				       sprintf(generalPrint + strlen(generalPrint), 
					       "Lower bound:                    %.3f\n", 
					       babModel_->getBestPossibleObjValue());
				       if (babModel_->bestSolution()){
					  sprintf(generalPrint + strlen(generalPrint), 
						  "Gap:                            %.2f\n", 
						  (babModel_->getObjValue()-babModel_->getBestPossibleObjValue())/babModel_->getBestPossibleObjValue());
				       }
				    }
				    sprintf(generalPrint + strlen(generalPrint), 
					    "Enumerated nodes:               %d\n", 
					    babModel_->getNodeCount());
				    sprintf(generalPrint + strlen(generalPrint), 
					    "Total iterations:               %d\n", 
					    babModel_->getIterationCount());
#ifndef CBC_QUIET 
				    sprintf(generalPrint + strlen(generalPrint), 
					    "Time (CPU seconds):             %.2f\n",
					    CoinCpuTime() - time1);
				    sprintf(generalPrint + strlen(generalPrint),
					    "Time (Wallclock seconds):       %.2f\n", 
					    CoinGetTimeOfDay() - time1Elapsed);
#endif
				    generalMessageHandler->message(CLP_GENERAL, generalMessages)
				       << generalPrint
				       << CoinMessageEol;
                                }
                                int returnCode = callBack(babModel_, 5);
                                if (returnCode) {
                                    // exit if user wants
                                    model_.moveInfo(*babModel_);
                                    delete babModel_;
                                    babModel_ = NULL;
                                    return returnCode;
                                }
                            } else {
                                std::cout << "Model strengthened - now has " << clpSolver->getNumRows()
                                          << " rows" << std::endl;
                            }
                            time1 = time2;
#ifdef COIN_HAS_ASL
                            if (statusUserFunction_[0]) {
                                // keep if going to be destroyed
                                OsiSolverInterface * solver = babModel_->solver();
                                OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
                                ClpSimplex * lpSolver2 = clpSolver->getModelPtr();
                                if (lpSolver == lpSolver2)
                                    babModel_->setModelOwnsSolver(false);
                            }
#endif
                            //delete babModel_;
                            //babModel_=NULL;
                        } else {
#ifndef DISALLOW_PRINTING
                            std::cout << "** Current model not valid" << std::endl ;
#endif
                        }
                        break ;
                    case CLP_PARAM_ACTION_IMPORT: {
#ifdef COIN_HAS_ASL
                        if (!statusUserFunction_[0]) {
#endif
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
#ifdef COIN_HAS_ASL
                        }
#endif
                        //delete babModel_;
                        //babModel_=NULL;
                        // get next field
                        field = CoinReadGetString(argc, argv);
                        if (field == "$") {
                            field = parameters_[iParam].stringValue();
                        } else if (field == "EOL") {
                            parameters_[iParam].printString();
                            break;
                        } else {
                            parameters_[iParam].setStringValue(field);
                        }
                        std::string fileName;
                        bool canOpen = false;
                        // See if gmpl file
                        int gmpl = 0;
                        std::string gmplData;
                        if (field == "-") {
                            // stdin
                            canOpen = true;
                            fileName = "-";
                        } else {
                            // See if .lp
                            {
                                const char * c_name = field.c_str();
                                size_t length = strlen(c_name);
                                if (length > 3 && !strncmp(c_name + length - 3, ".lp", 3))
                                    gmpl = -1; // .lp
                            }
                            bool absolutePath;
                            if (dirsep == '/') {
                                // non Windows (or cygwin)
                                absolutePath = (field[0] == '/');
                            } else {
                                //Windows (non cycgwin)
                                absolutePath = (field[0] == '\\');
                                // but allow for :
                                if (strchr(field.c_str(), ':'))
                                    absolutePath = true;
                            }
                            if (absolutePath) {
                                fileName = field;
				size_t length = field.size();
				size_t percent = field.find('%');
				if (percent < length && percent > 0) {
				  gmpl = 1;
				  fileName = field.substr(0, percent);
				  gmplData = field.substr(percent + 1);
				  if (percent < length - 1)
				    gmpl = 2; // two files
				  printf("GMPL model file %s and data file %s\n",
					 fileName.c_str(), gmplData.c_str());
				}
                            } else if (field[0] == '~') {
                                char * environVar = getenv("HOME");
                                if (environVar) {
                                    std::string home(environVar);
                                    field = field.erase(0, 1);
                                    fileName = home + field;
                                } else {
                                    fileName = field;
                                }
                            } else {
                                fileName = directory + field;
                                // See if gmpl (model & data) - or even lp file
                                size_t length = field.size();
                                size_t percent = field.find('%');
                                if (percent<length && percent>0) {
                                    gmpl = 1;
                                    fileName = directory + field.substr(0, percent);
                                    gmplData = directory + field.substr(percent + 1);
                                    if (percent < length - 1)
                                        gmpl = 2; // two files
                                    printf("GMPL model file %s and data file %s\n",
                                           fileName.c_str(), gmplData.c_str());
                                }
                            }
                            std::string name = fileName;
                            if (fileCoinReadable(name)) {
                                // can open - lets go for it
                                canOpen = true;
                                if (gmpl == 2) {
                                    FILE *fp;
                                    fp = fopen(gmplData.c_str(), "r");
                                    if (fp) {
                                        fclose(fp);
                                    } else {
                                        canOpen = false;
                                        std::cout << "Unable to open file " << gmplData << std::endl;
                                    }
                                }
                            } else {
                                std::cout << "Unable to open file " << fileName << std::endl;
                            }
                        }
                        if (canOpen) {
                            int status;
#ifndef CBC_OTHER_SOLVER
                            ClpSimplex * lpSolver = clpSolver->getModelPtr();
                            if (!gmpl) {
                                status = clpSolver->readMps(fileName.c_str(),
                                                            keepImportNames != 0,
                                                            allowImportErrors != 0);
                            } else if (gmpl > 0) {
                                status = lpSolver->readGMPL(fileName.c_str(),
                                                            (gmpl == 2) ? gmplData.c_str() : NULL,
                                                            keepImportNames != 0);
                            } else {
#ifdef KILL_ZERO_READLP
                                status = lpSolver->readLp(fileName.c_str(), lpSolver->getSmallElementValue());
#else
                                status = lpSolver->readLp(fileName.c_str(), 1.0e-12);
#endif
                            }
#else
                            status = clpSolver->readMps(fileName.c_str(), "");
#endif
                            if (!status || (status > 0 && allowImportErrors)) {
#ifndef CBC_OTHER_SOLVER
                                if (keepImportNames && gmpl <= 0) {
                                    lengthName = lpSolver->lengthNames();
                                    rowNames = *(lpSolver->rowNames());
                                    columnNames = *(lpSolver->columnNames());
                                } else {
                                    lengthName = 0;
                                }
                                goodModel = true;
                                // sets to all slack (not necessary?)
                                lpSolver->createStatus();
                                // make sure integer
                                int numberColumns = lpSolver->numberColumns();
                                for (int i = 0; i < numberColumns; i++) {
                                    if (lpSolver->isInteger(i))
                                        clpSolver->setInteger(i);
                                }
#else
                                lengthName = 0;
                                goodModel = true;
#endif
                                time2 = CoinCpuTime();
                                totalTime += time2 - time1;
                                time1 = time2;
                                // Go to canned file if just input file
                                if (CbcOrClpRead_mode == 2 && argc == 2) {
                                    // only if ends .mps
                                    char * find = const_cast<char *>(strstr(fileName.c_str(), ".mps"));
                                    if (find && find[4] == '\0') {
                                        find[1] = 'p';
                                        find[2] = 'a';
                                        find[3] = 'r';
                                        FILE *fp = fopen(fileName.c_str(), "r");
                                        if (fp) {
                                            CbcOrClpReadCommand = fp; // Read from that file
                                            CbcOrClpRead_mode = -1;
                                        }
                                    }
                                }
                            } else {
                                // errors
                                std::cout << "There were " << status <<
                                          " errors on input" << std::endl;
                            }
                        }
                    }
                    break;
                    case CLP_PARAM_ACTION_MODELIN:
#ifndef CBC_OTHER_SOLVER
#ifdef COIN_HAS_LINK
                        {
                            // get next field
                            field = CoinReadGetString(argc, argv);
                            if (field == "$") {
                                field = parameters_[iParam].stringValue();
                            } else if (field == "EOL") {
                                parameters_[iParam].printString();
                                break;
                            } else {
                                parameters_[iParam].setStringValue(field);
                            }
                            std::string fileName;
                            bool canOpen = false;
                            if (field == "-") {
                                // stdin
                                canOpen = true;
                                fileName = "-";
                            } else {
                                bool absolutePath;
                                if (dirsep == '/') {
                                    // non Windows (or cygwin)
                                    absolutePath = (field[0] == '/');
                                } else {
                                    //Windows (non cycgwin)
                                    absolutePath = (field[0] == '\\');
                                    // but allow for :
                                    if (strchr(field.c_str(), ':'))
                                        absolutePath = true;
                                }
                                if (absolutePath) {
                                    fileName = field;
                                } else if (field[0] == '~') {
                                    char * environVar = getenv("HOME");
                                    if (environVar) {
                                        std::string home(environVar);
                                        field = field.erase(0, 1);
                                        fileName = home + field;
                                    } else {
                                        fileName = field;
                                    }
                                } else {
                                    fileName = directory + field;
                                }
                                FILE *fp = fopen(fileName.c_str(), "r");
                                if (fp) {
                                    // can open - lets go for it
                                    fclose(fp);
                                    canOpen = true;
                                } else {
                                    std::cout << "Unable to open file " << fileName << std::endl;
                                }
                            }
                            if (canOpen) {
                                CoinModel coinModel(fileName.c_str(), 2);
                                // load from coin model
                                OsiSolverLink solver1;
                                OsiSolverInterface * solver2 = solver1.clone();
                                model_.assignSolver(solver2, false);
                                OsiSolverLink * si =
                                    dynamic_cast<OsiSolverLink *>(model_.solver()) ;
                                assert (si != NULL);
                                si->setDefaultMeshSize(0.001);
                                // need some relative granularity
                                si->setDefaultBound(100.0);
                                double dextra3 = parameters_[whichParam(CBC_PARAM_DBL_DEXTRA3, numberParameters_, parameters_)].doubleValue();
                                if (dextra3)
                                    si->setDefaultMeshSize(dextra3);
                                si->setDefaultBound(100.0);
                                si->setIntegerPriority(1000);
                                si->setBiLinearPriority(10000);
                                CoinModel * model2 = &coinModel;
                                si->load(*model2);
                                // redo
                                solver = model_.solver();
                                clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
                                lpSolver = clpSolver->getModelPtr();
                                clpSolver->messageHandler()->setLogLevel(0) ;
                                testOsiParameters = 0;
                                complicatedInteger = 2;
                            }
                        }
#endif
#endif
                        break;
                    case CLP_PARAM_ACTION_EXPORT:
                        if (goodModel) {
                            // get next field
                            field = CoinReadGetString(argc, argv);
                            if (field == "$") {
                                field = parameters_[iParam].stringValue();
                            } else if (field == "EOL") {
                                parameters_[iParam].printString();
                                break;
                            } else {
                                parameters_[iParam].setStringValue(field);
                            }
                            std::string fileName;
                            bool canOpen = false;
                            if (field[0] == '/' || field[0] == '\\') {
                                fileName = field;
                            } else if (field[0] == '~') {
                                char * environVar = getenv("HOME");
                                if (environVar) {
                                    std::string home(environVar);
                                    field = field.erase(0, 1);
                                    fileName = home + field;
                                } else {
                                    fileName = field;
                                }
                            } else {
                                fileName = directory + field;
                            }
                            FILE *fp = fopen(fileName.c_str(), "w");
                            if (fp) {
                                // can open - lets go for it
                                fclose(fp);
                                canOpen = true;
                            } else {
                                std::cout << "Unable to open file " << fileName << std::endl;
                            }
                            if (canOpen) {
                                // If presolve on then save presolved
                                bool deleteModel2 = false;
                                ClpSimplex * model2 = lpSolver;
                                if (dualize && dualize < 3) {
                                    model2 = static_cast<ClpSimplexOther *> (model2)->dualOfModel();
                                    sprintf(generalPrint, "Dual of model has %d rows and %d columns",
                                            model2->numberRows(), model2->numberColumns());
                                    generalMessageHandler->message(CLP_GENERAL, generalMessages)
                                    << generalPrint
                                    << CoinMessageEol;
                                    model2->setOptimizationDirection(1.0);
                                }
#ifndef CBC_OTHER_SOLVER
#ifdef COIN_HAS_ASL
                                if (info.numberSos && doSOS && statusUserFunction_[0]) {
                                    // SOS
                                    numberSOS = info.numberSos;
                                    sosStart = info.sosStart;
                                    sosIndices = info.sosIndices;
                                    sosReference = info.sosReference;
                                    preSolve = false;
                                    clpSolver->setSOSData(numberSOS, info.sosType, sosStart, sosIndices, sosReference);
                                }
#endif
#endif
                                if (preSolve) {
                                    ClpPresolve pinfo;
                                    int presolveOptions2 = presolveOptions&~0x40000000;
                                    if ((presolveOptions2&0xffff) != 0)
                                        pinfo.setPresolveActions(presolveOptions2);
                                    if ((printOptions&1) != 0)
                                        pinfo.statistics();
                                    double presolveTolerance =
                                        parameters_[whichParam(CLP_PARAM_DBL_PRESOLVETOLERANCE, numberParameters_, parameters_)].doubleValue();
                                    model2 =
                                        pinfo.presolvedModel(*lpSolver, presolveTolerance,
                                                             true, preSolve);
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
                                    model2->writeMps(fileName.c_str(), (outputFormat - 1) / 2, 1 + ((outputFormat - 1)&1));
                                    if (deleteModel2)
                                        delete model2;
                                } else {
                                    printf("Saving model on %s\n",
                                           fileName.c_str());
                                    if (numberSOS) {
                                        // Convert names
                                        int iRow;
                                        int numberRows = model2->numberRows();
                                        int iColumn;
                                        int numberColumns = model2->numberColumns();

                                        char ** rowNames = NULL;
                                        char ** columnNames = NULL;
                                        if (model2->lengthNames()) {
                                            rowNames = new char * [numberRows];
                                            for (iRow = 0; iRow < numberRows; iRow++) {
                                                rowNames[iRow] =
                                                    CoinStrdup(model2->rowName(iRow).c_str());
                                            }

                                            columnNames = new char * [numberColumns];
                                            for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                                                columnNames[iColumn] =
                                                    CoinStrdup(model2->columnName(iColumn).c_str());
                                            }
                                        }
                                        clpSolver->writeMpsNative(fileName.c_str(), const_cast<const char **> (rowNames), const_cast<const char **> (columnNames),
                                                                  (outputFormat - 1) / 2, 1 + ((outputFormat - 1)&1));
                                        if (rowNames) {
                                            for (iRow = 0; iRow < numberRows; iRow++) {
                                                free(rowNames[iRow]);
                                            }
                                            delete [] rowNames;
                                            for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                                                free(columnNames[iColumn]);
                                            }
                                            delete [] columnNames;
                                        }
                                    } else {
#ifdef COIN_HAS_LINK
                                        OsiSolverLink * linkSolver = dynamic_cast< OsiSolverLink*> (clpSolver);
                                        if (!linkSolver || !linkSolver->quadraticModel())
                                            model2->writeMps(fileName.c_str(), (outputFormat - 1) / 2, 1 + ((outputFormat - 1)&1));
                                        else
                                            linkSolver->quadraticModel()->writeMps(fileName.c_str(), (outputFormat - 1) / 2, 1 + ((outputFormat - 1)&1));
#endif
                                    }
                                }
                                time2 = CoinCpuTime();
                                totalTime += time2 - time1;
                                time1 = time2;
                            }
                        } else {
#ifndef DISALLOW_PRINTING
                            std::cout << "** Current model not valid" << std::endl;
#endif
                        }
                        break;
                    case CLP_PARAM_ACTION_BASISIN:
                        if (goodModel) {
                            // get next field
                            field = CoinReadGetString(argc, argv);
                            if (field == "$") {
                                field = parameters_[iParam].stringValue();
                            } else if (field == "EOL") {
                                parameters_[iParam].printString();
                                break;
                            } else {
                                parameters_[iParam].setStringValue(field);
                            }
                            std::string fileName;
                            bool canOpen = false;
                            if (field == "-") {
                                // stdin
                                canOpen = true;
                                fileName = "-";
                            } else {
                                if (field[0] == '/' || field[0] == '\\') {
                                    fileName = field;
                                } else if (field[0] == '~') {
                                    char * environVar = getenv("HOME");
                                    if (environVar) {
                                        std::string home(environVar);
                                        field = field.erase(0, 1);
                                        fileName = home + field;
                                    } else {
                                        fileName = field;
                                    }
                                } else {
                                    fileName = directory + field;
                                }
                                FILE *fp = fopen(fileName.c_str(), "r");
                                if (fp) {
                                    // can open - lets go for it
                                    fclose(fp);
                                    canOpen = true;
                                } else {
                                    std::cout << "Unable to open file " << fileName << std::endl;
                                }
                            }
                            if (canOpen) {
#ifndef CBC_OTHER_SOLVER
                                int values = lpSolver->readBasis(fileName.c_str());
                                if (values == 0)
                                    basisHasValues = -1;
                                else
                                    basisHasValues = 1;
                                assert (lpSolver == clpSolver->getModelPtr());
                                clpSolver->setWarmStart(NULL);
#endif
                            }
                        } else {
#ifndef DISALLOW_PRINTING
                            std::cout << "** Current model not valid" << std::endl;
#endif
                        }
                        break;
                    case CBC_PARAM_ACTION_PRIORITYIN:
                        if (goodModel) {
                            // get next field
                            field = CoinReadGetString(argc, argv);
                            if (field == "$") {
                                field = parameters_[iParam].stringValue();
                            } else if (field == "EOL") {
                                parameters_[iParam].printString();
                                break;
                            } else {
                                parameters_[iParam].setStringValue(field);
                            }
                            std::string fileName;
                            if (field[0] == '/' || field[0] == '\\') {
                                fileName = field;
                            } else if (field[0] == '~') {
                                char * environVar = getenv("HOME");
                                if (environVar) {
                                    std::string home(environVar);
                                    field = field.erase(0, 1);
                                    fileName = home + field;
                                } else {
                                    fileName = field;
                                }
                            } else {
                                fileName = directory + field;
                            }
                            FILE *fp = fopen(fileName.c_str(), "r");
                            if (fp) {
                                // can open - lets go for it
                                std::string headings[] = {"name", "number", "direction", "priority", "up", "down",
                                                          "solution", "priin"
                                                         };
                                int got[] = { -1, -1, -1, -1, -1, -1, -1, -1};
                                int order[8];
                                assert(sizeof(got) == sizeof(order));
                                int nAcross = 0;
                                char line[1000];
                                int numberColumns = lpSolver->numberColumns();
                                if (!fgets(line, 1000, fp)) {
                                    std::cout << "Odd file " << fileName << std::endl;
                                } else {
                                    char * pos = line;
                                    char * put = line;
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
                                        char * comma = strchr(pos, ',');
                                        if (comma)
                                            *comma = '\0';
                                        for (i = 0; i < static_cast<int> (sizeof(got) / sizeof(int)); i++) {
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
                                        if (i == static_cast<int> (sizeof(got) / sizeof(int)))
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
                                    if (good && (strstr(fileName.c_str(), ".mst") || strstr(fileName.c_str(), ".MST") || strstr(fileName.c_str(), ".csv"))) {
                                        numberFields = 0;
                                        for (i = 2; i < static_cast<int> (sizeof(got) / sizeof(int)); i++) {
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
                                        char ** columnNames = new char * [numberColumns];
                                        pseudoDown = reinterpret_cast<double *> (malloc(numberColumns * sizeof(double)));
                                        pseudoUp = reinterpret_cast<double *> (malloc(numberColumns * sizeof(double)));
                                        branchDirection = reinterpret_cast<int *> (malloc(numberColumns * sizeof(int)));
                                        priorities = reinterpret_cast<int *> (malloc(numberColumns * sizeof(int)));
                                        free(solutionIn);
                                        solutionIn = NULL;
                                        free(prioritiesIn);
                                        prioritiesIn = NULL;
                                        int iColumn;
                                        if (got[6] >= 0) {
                                            solutionIn = reinterpret_cast<double *> (malloc(numberColumns * sizeof(double)));
                                            for (iColumn = 0; iColumn < numberColumns; iColumn++)
                                                solutionIn[iColumn] = -COIN_DBL_MAX;
                                        }
                                        if (got[7] >= 0 || !numberFields) {
                                            prioritiesIn = reinterpret_cast<int *> (malloc(numberColumns * sizeof(int)));
                                            for (iColumn = 0; iColumn < numberColumns; iColumn++)
                                                prioritiesIn[iColumn] = 10000;
                                        }
                                        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                                            columnNames[iColumn] =
                                                CoinStrdup(lpSolver->columnName(iColumn).c_str());
                                            pseudoDown[iColumn] = 0.0;
                                            pseudoUp[iColumn] = 0.0;
                                            branchDirection[iColumn] = 0;
                                            priorities[iColumn] = 0;
                                        }
                                        int nBadPseudo = 0;
                                        int nBadDir = 0;
                                        int nBadPri = 0;
                                        int nBadName = 0;
                                        int nBadLine = 0;
                                        int nLine = 0;
                                        while (fgets(line, 1000, fp)) {
                                            if (!strncmp(line, "ENDATA", 6))
                                                break;
                                            nLine++;
                                            iColumn = -1;
                                            double up = 0.0;
                                            double down = 0.0;
                                            int pri = 0;
                                            int dir = 0;
                                            double solValue = COIN_DBL_MAX;
                                            int priValue = 1000000;
                                            char * pos = line;
                                            char * put = line;
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
                                                char * comma = strchr(pos, ',');
                                                if (comma) {
                                                    *comma = '\0';
                                                } else if (i < nAcross - 1) {
                                                    nBadLine++;
                                                    break;
                                                }
                                                switch (order[i]) {
                                                    // name
                                                case 0:
                                                    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                                                        if (!strcmp(columnNames[iColumn], pos))
                                                            break;
                                                    }
                                                    if (iColumn == numberColumns)
                                                        iColumn = -1;
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
                                                    else if (*pos == '1' && *(pos + 1) == '1' && *(pos + 2) == '\0')
                                                        dir = -1;
                                                    else
                                                        dir = -2; // bad
                                                    break;
                                                    // priority
                                                case 3:
                                                    pri = atoi(pos);
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
                                                pseudoDown[iColumn] = down;
                                                pseudoUp[iColumn] = up;
                                                branchDirection[iColumn] = dir;
                                                priorities[iColumn] = pri;
                                                if (solValue != COIN_DBL_MAX) {
                                                    assert (solutionIn);
                                                    solutionIn[iColumn] = solValue;
                                                }
                                                if (priValue != 1000000) {
                                                    assert (prioritiesIn);
                                                    prioritiesIn[iColumn] = priValue;
                                                }
                                            } else {
                                                nBadName++;
                                            }
                                        }
                                        if (!noPrinting_) {
                                            printf("%d fields and %d records", nAcross, nLine);
                                            if (nBadPseudo)
                                                printf(" %d bad pseudo costs", nBadPseudo);
                                            if (nBadDir)
                                                printf(" %d bad directions", nBadDir);
                                            if (nBadPri)
                                                printf(" %d bad priorities", nBadPri);
                                            if (nBadName)
                                                printf(" ** %d records did not match on name/sequence", nBadName);
                                            printf("\n");
                                        }
                                        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                                            free(columnNames[iColumn]);
                                        }
                                        delete [] columnNames;
                                    } else {
                                        std::cout << "Duplicate or unknown keyword - or name/number fields wrong" << line << std::endl;
                                    }
                                }
                                fclose(fp);
                            } else {
                                std::cout << "Unable to open file " << fileName << std::endl;
                            }
                        } else {
#ifndef DISALLOW_PRINTING
                            std::cout << "** Current model not valid" << std::endl;
#endif
                        }
                        break;
                    case CLP_PARAM_ACTION_DEBUG:
                        if (goodModel) {
                            delete [] debugValues;
                            debugValues = NULL;
                            // get next field
                            field = CoinReadGetString(argc, argv);
                            if (field == "$") {
                                field = parameters_[iParam].stringValue();
                            } else if (field == "EOL") {
                                parameters_[iParam].printString();
                                break;
                            } else {
                                parameters_[iParam].setStringValue(field);
                                debugFile = field;
                                if (debugFile == "create" ||
                                        debugFile == "createAfterPre") {
                                    printf("Will create a debug file so this run should be a good one\n");
                                    break;
                                }
                            }
                            std::string fileName;
                            if (field[0] == '/' || field[0] == '\\') {
                                fileName = field;
                            } else if (field[0] == '~') {
                                char * environVar = getenv("HOME");
                                if (environVar) {
                                    std::string home(environVar);
                                    field = field.erase(0, 1);
                                    fileName = home + field;
                                } else {
                                    fileName = field;
                                }
                            } else {
                                fileName = directory + field;
                            }
                            FILE *fp = fopen(fileName.c_str(), "rb");
                            if (fp) {
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
                                debugValues = new double[numberDebugValues+numRows];
                                nRead = fread(debugValues, sizeof(double), numRows, fp);
                                if (nRead != static_cast<size_t>(numRows))
                                    throw("Error in fread");
                                nRead = fread(debugValues, sizeof(double), numRows, fp);
                                if (nRead != static_cast<size_t>(numRows))
                                    throw("Error in fread");
                                nRead = fread(debugValues, sizeof(double), numberDebugValues, fp);
                                if (nRead != static_cast<size_t>(numberDebugValues))
                                    throw("Error in fread");
                                printf("%d doubles read into debugValues\n", numberDebugValues);
#ifdef CGL_WRITEMPS
				debugSolution = debugValues;
				debugNumberColumns = numberDebugValues;
#endif
                                if (numberDebugValues < 200) {
                                    for (int i = 0; i < numberDebugValues; i++) {
                                        if (clpSolver->isInteger(i) && debugValues[i])
                                            printf("%d %g\n", i, debugValues[i]);
                                    }
                                }
                                fclose(fp);
                            } else {
                                std::cout << "Unable to open file " << fileName << std::endl;
                            }
                        } else {
#ifndef DISALLOW_PRINTING
                            std::cout << "** Current model not valid" << std::endl;
#endif
                        }
                        break;
                    case CLP_PARAM_ACTION_PRINTMASK:
                        // get next field
                    {
                        std::string name = CoinReadGetString(argc, argv);
                        if (name != "EOL") {
                            parameters_[iParam].setStringValue(name);
                            printMask = name;
                        } else {
                            parameters_[iParam].printString();
                        }
                    }
                    break;
                    case CLP_PARAM_ACTION_BASISOUT:
                        if (goodModel) {
                            // get next field
                            field = CoinReadGetString(argc, argv);
                            if (field == "$") {
                                field = parameters_[iParam].stringValue();
                            } else if (field == "EOL") {
                                parameters_[iParam].printString();
                                break;
                            } else {
                                parameters_[iParam].setStringValue(field);
                            }
                            std::string fileName;
                            bool canOpen = false;
                            if (field[0] == '/' || field[0] == '\\') {
                                fileName = field;
                            } else if (field[0] == '~') {
                                char * environVar = getenv("HOME");
                                if (environVar) {
                                    std::string home(environVar);
                                    field = field.erase(0, 1);
                                    fileName = home + field;
                                } else {
                                    fileName = field;
                                }
                            } else {
                                fileName = directory + field;
                            }
                            FILE *fp = fopen(fileName.c_str(), "w");
                            if (fp) {
                                // can open - lets go for it
                                fclose(fp);
                                canOpen = true;
                            } else {
                                std::cout << "Unable to open file " << fileName << std::endl;
                            }
                            if (canOpen) {
                                ClpSimplex * model2 = lpSolver;
                                model2->writeBasis(fileName.c_str(), outputFormat > 1, outputFormat - 2);
                                time2 = CoinCpuTime();
                                totalTime += time2 - time1;
                                time1 = time2;
                            }
                        } else {
#ifndef DISALLOW_PRINTING
                            std::cout << "** Current model not valid" << std::endl;
#endif
                        }
                        break;
                    case CLP_PARAM_ACTION_SAVE: {
                        // get next field
                        field = CoinReadGetString(argc, argv);
                        if (field == "$") {
                            field = parameters_[iParam].stringValue();
                        } else if (field == "EOL") {
                            parameters_[iParam].printString();
                            break;
                        } else {
                            parameters_[iParam].setStringValue(field);
                        }
                        std::string fileName;
                        bool canOpen = false;
                        if (field[0] == '/' || field[0] == '\\') {
                            fileName = field;
                        } else if (field[0] == '~') {
                            char * environVar = getenv("HOME");
                            if (environVar) {
                                std::string home(environVar);
                                field = field.erase(0, 1);
                                fileName = home + field;
                            } else {
                                fileName = field;
                            }
                        } else {
                            fileName = directory + field;
                        }
                        FILE *fp = fopen(fileName.c_str(), "wb");
                        if (fp) {
                            // can open - lets go for it
                            fclose(fp);
                            canOpen = true;
                        } else {
                            std::cout << "Unable to open file " << fileName << std::endl;
                        }
                        if (canOpen) {
                            int status;
                            // If presolve on then save presolved
                            bool deleteModel2 = false;
                            ClpSimplex * model2 = lpSolver;
                            if (preSolve) {
                                ClpPresolve pinfo;
                                double presolveTolerance =
                                    parameters_[whichParam(CLP_PARAM_DBL_PRESOLVETOLERANCE, numberParameters_, parameters_)].doubleValue();
                                model2 =
                                    pinfo.presolvedModel(*lpSolver, presolveTolerance,
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
                                std::cout << "There were errors on output" << std::endl;
                            }
                        }
                    }
                    break;
                    case CLP_PARAM_ACTION_RESTORE: {
                        // get next field
                        field = CoinReadGetString(argc, argv);
                        if (field == "$") {
                            field = parameters_[iParam].stringValue();
                        } else if (field == "EOL") {
                            parameters_[iParam].printString();
                            break;
                        } else {
                            parameters_[iParam].setStringValue(field);
                        }
                        std::string fileName;
                        bool canOpen = false;
                        if (field[0] == '/' || field[0] == '\\') {
                            fileName = field;
                        } else if (field[0] == '~') {
                            char * environVar = getenv("HOME");
                            if (environVar) {
                                std::string home(environVar);
                                field = field.erase(0, 1);
                                fileName = home + field;
                            } else {
                                fileName = field;
                            }
                        } else {
                            fileName = directory + field;
                        }
                        FILE *fp = fopen(fileName.c_str(), "rb");
                        if (fp) {
                            // can open - lets go for it
                            fclose(fp);
                            canOpen = true;
                        } else {
                            std::cout << "Unable to open file " << fileName << std::endl;
                        }
                        if (canOpen) {
                            int status = lpSolver->restoreModel(fileName.c_str());
                            if (!status) {
                                goodModel = true;
                                time2 = CoinCpuTime();
                                totalTime += time2 - time1;
                                time1 = time2;
                            } else {
                                // errors
                                std::cout << "There were errors on input" << std::endl;
                            }
                        }
                    }
                    break;
                    case CLP_PARAM_ACTION_MAXIMIZE:
                        lpSolver->setOptimizationDirection(-1);
                        break;
                    case CLP_PARAM_ACTION_MINIMIZE:
                        lpSolver->setOptimizationDirection(1);
                        break;
                    case CLP_PARAM_ACTION_ALLSLACK:
                        lpSolver->allSlackBasis(true);
                        break;
                    case CLP_PARAM_ACTION_REVERSE:
                        if (goodModel) {
                            int iColumn;
                            int numberColumns = lpSolver->numberColumns();
                            double * dualColumnSolution =
                                lpSolver->dualColumnSolution();
                            ClpObjective * obj = lpSolver->objectiveAsObject();
                            assert(dynamic_cast<ClpLinearObjective *> (obj));
                            double offset;
                            double * objective = obj->gradient(NULL, NULL, offset, true);
                            for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                                dualColumnSolution[iColumn] = dualColumnSolution[iColumn];
                                objective[iColumn] = -objective[iColumn];
                            }
                            int iRow;
                            int numberRows = lpSolver->numberRows();
                            double * dualRowSolution =
                                lpSolver->dualRowSolution();
                            for (iRow = 0; iRow < numberRows; iRow++)
                                dualRowSolution[iRow] = dualRowSolution[iRow];
                        }
                        break;
                    case CLP_PARAM_ACTION_DIRECTORY: {
                        std::string name = CoinReadGetString(argc, argv);
                        if (name != "EOL") {
                            size_t length = name.length();
                            if (length > 0 && name[length-1] == dirsep) {
                                directory = name;
                            } else {
                                directory = name + dirsep;
                            }
                            parameters_[iParam].setStringValue(directory);
                        } else {
                            parameters_[iParam].printString();
                        }
                    }
                    break;
                    case CLP_PARAM_ACTION_DIRSAMPLE: {
                        std::string name = CoinReadGetString(argc, argv);
                        if (name != "EOL") {
                            size_t length = name.length();
                            if (length > 0 && name[length-1] == dirsep) {
                                dirSample = name;
                            } else {
                                dirSample = name + dirsep;
                            }
                            parameters_[iParam].setStringValue(dirSample);
                        } else {
                            parameters_[iParam].printString();
                        }
                    }
                    break;
                    case CLP_PARAM_ACTION_DIRNETLIB: {
                        std::string name = CoinReadGetString(argc, argv);
                        if (name != "EOL") {
                            size_t length = name.length();
                            if (length > 0 && name[length-1] == dirsep) {
                                dirNetlib = name;
                            } else {
                                dirNetlib = name + dirsep;
                            }
                            parameters_[iParam].setStringValue(dirNetlib);
                        } else {
                            parameters_[iParam].printString();
                        }
                    }
                    break;
                    case CBC_PARAM_ACTION_DIRMIPLIB: {
                        std::string name = CoinReadGetString(argc, argv);
                        if (name != "EOL") {
                            size_t length = name.length();
                            if (length > 0 && name[length-1] == dirsep) {
                                dirMiplib = name;
                            } else {
                                dirMiplib = name + dirsep;
                            }
                            parameters_[iParam].setStringValue(dirMiplib);
                        } else {
                            parameters_[iParam].printString();
                        }
                    }
                    break;
                    case CLP_PARAM_ACTION_STDIN:
                        CbcOrClpRead_mode = -1;
                        break;
                    case CLP_PARAM_ACTION_NETLIB_DUAL:
                    case CLP_PARAM_ACTION_NETLIB_EITHER:
                    case CLP_PARAM_ACTION_NETLIB_BARRIER:
                    case CLP_PARAM_ACTION_NETLIB_PRIMAL:
                    case CLP_PARAM_ACTION_NETLIB_TUNE: {
                        printf("unit test is now only from clp - does same thing\n");
                        //return(22);
                    }
                    break;
                    case CLP_PARAM_ACTION_UNITTEST: {
                        CbcClpUnitTest(model_, dirSample, -2, NULL);
                    }
                    break;
                    case CLP_PARAM_ACTION_FAKEBOUND:
                        if (goodModel) {
                            // get bound
                            double value = CoinReadGetDoubleField(argc, argv, &valid);
                            if (!valid) {
                                std::cout << "Setting " << parameters_[iParam].name() <<
                                          " to DEBUG " << value << std::endl;
                                int iRow;
                                int numberRows = lpSolver->numberRows();
                                double * rowLower = lpSolver->rowLower();
                                double * rowUpper = lpSolver->rowUpper();
                                for (iRow = 0; iRow < numberRows; iRow++) {
                                    // leave free ones for now
                                    if (rowLower[iRow] > -1.0e20 || rowUpper[iRow] < 1.0e20) {
                                        rowLower[iRow] = CoinMax(rowLower[iRow], -value);
                                        rowUpper[iRow] = CoinMin(rowUpper[iRow], value);
                                    }
                                }
                                int iColumn;
                                int numberColumns = lpSolver->numberColumns();
                                double * columnLower = lpSolver->columnLower();
                                double * columnUpper = lpSolver->columnUpper();
                                for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                                    // leave free ones for now
                                    if (columnLower[iColumn] > -1.0e20 ||
                                            columnUpper[iColumn] < 1.0e20) {
                                        columnLower[iColumn] = CoinMax(columnLower[iColumn], -value);
                                        columnUpper[iColumn] = CoinMin(columnUpper[iColumn], value);
                                    }
                                }
                            } else if (valid == 1) {
                                abort();
                            } else {
                                std::cout << "enter value for " << parameters_[iParam].name() <<
                                          std::endl;
                            }
                        }
                        break;
                    case CLP_PARAM_ACTION_REALLY_SCALE:
                        if (goodModel) {
                            ClpSimplex newModel(*lpSolver,
                                                lpSolver->scalingFlag());
                            printf("model really really scaled\n");
                            *lpSolver = newModel;
                        }
                        break;
                    case CLP_PARAM_ACTION_USERCLP:
#ifdef USER_HAS_FAKE_CLP
                        // Replace the sample code by whatever you want
                        if (goodModel) {
                            // Way of using an existing piece of code
                            OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (model_.solver());
                            ClpSimplex * lpSolver = clpSolver->getModelPtr();
                            // set time from integer model
                            double timeToGo = model_.getMaximumSeconds();
                            lpSolver->setMaximumSeconds(timeToGo);
                            int extra1 = parameters_[whichParam(CBC_PARAM_INT_EXTRA1, numberParameters_, parameters_)].intValue();
                            fakeMain2(*lpSolver, *clpSolver, extra1);
                            lpSolver = clpSolver->getModelPtr();
#ifdef COIN_HAS_ASL
                            // My actual usage has objective only in clpSolver
                            //double objectiveValue=clpSolver->getObjValue();
                            //int iStat = lpSolver->status();
                            //int iStat2 = lpSolver->secondaryStatus();
#endif
                        }
#endif
                        break;
                    case CBC_PARAM_ACTION_USERCBC:
#ifdef USER_HAS_FAKE_CBC
                        // Replace the sample code by whatever you want
                        if (goodModel) {
                            // Way of using an existing piece of code
                            OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (model_.solver());
                            ClpSimplex * lpSolver = clpSolver->getModelPtr();
                            // set time from integer model
                            double timeToGo = model_.getMaximumSeconds();
                            lpSolver->setMaximumSeconds(timeToGo);
                            fakeMain(*lpSolver, *clpSolver, model);
#ifdef COIN_HAS_ASL
                            // My actual usage has objective only in clpSolver
                            double objectiveValue = clpSolver->getObjValue();
                            int iStat = lpSolver->status();
                            int iStat2 = lpSolver->secondaryStatus();
#endif
                            // make sure solution back in correct place
                            clpSolver = dynamic_cast< OsiClpSolverInterface*> (model_.solver());
                            lpSolver = clpSolver->getModelPtr();
#ifdef COIN_HAS_ASL
                            if (statusUserFunction_[0]) {
                                int n = clpSolver->getNumCols();
                                double value = objectiveValue * lpSolver->getObjSense();
                                char buf[300];
                                int pos = 0;
                                std::string minor[] = {"", "", "gap", "nodes", "time", "", "solutions", "user ctrl-c"};
                                if (iStat == 0) {
                                    if (objectiveValue < 1.0e40) {
                                        pos += sprintf(buf + pos, "optimal," );
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
                                    pos += sprintf(buf + pos, "stopped on %s,", minor[iStat2].c_str());
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
                                info.problemStatus = iStat;
                                info.objValue = value;
                                if (objectiveValue < 1.0e40)
                                    pos += sprintf(buf + pos, " objective %.*g", ampl_obj_prec(),
                                                   value);
                                sprintf(buf + pos, "\n%d nodes, %d iterations",
                                        model_.getNodeCount(),
                                        model_.getIterationCount());
                                if (objectiveValue < 1.0e50) {
                                    free(info.primalSolution);
                                    info.primalSolution = (double *) malloc(n * sizeof(double));
                                    CoinCopyN(lpSolver->primalColumnSolution(), n, info.primalSolution);
                                    int numberRows = lpSolver->numberRows();
                                    free(info.dualSolution);
                                    info.dualSolution = (double *) malloc(numberRows * sizeof(double));
                                    CoinCopyN(lpSolver->dualRowSolution(), numberRows, info.dualSolution);
                                } else {
                                    info.primalSolution = NULL;
                                    info.dualSolution = NULL;
                                }
                                // put buffer into info
                                strcpy(info.buffer, buf);
                            }
#endif
                        }
#endif
                        break;
                    case CLP_PARAM_ACTION_HELP:
                        std::cout << "Cbc version " << CBC_VERSION
                                  << ", build " << __DATE__ << std::endl;
                        std::cout << "Non default values:-" << std::endl;
                        std::cout << "Perturbation " << lpSolver->perturbation() << " (default 100)"
                                  << std::endl;
                        CoinReadPrintit(
                            "Presolve being done with 5 passes\n\
Dual steepest edge steep/partial on matrix shape and factorization density\n\
Clpnnnn taken out of messages\n\
If Factorization frequency default then done on size of matrix\n\n\
(-)unitTest, (-)netlib or (-)netlibp will do standard tests\n\n\
You can switch to interactive mode at any time so\n\
clp watson.mps -scaling off -primalsimplex\nis the same as\n\
clp watson.mps -\nscaling off\nprimalsimplex"
                        );
                        break;
                    case CLP_PARAM_ACTION_CSVSTATISTICS: {
                        // get next field
                        field = CoinReadGetString(argc, argv);
                        if (field == "$") {
                            field = parameters_[iParam].stringValue();
                        } else if (field == "EOL") {
                            parameters_[iParam].printString();
                            break;
                        } else {
                            parameters_[iParam].setStringValue(field);
                        }
                        std::string fileName;
                        if (field[0] == '/' || field[0] == '\\') {
                            fileName = field;
                        } else if (field[0] == '~') {
                            char * environVar = getenv("HOME");
                            if (environVar) {
                                std::string home(environVar);
                                field = field.erase(0, 1);
                                fileName = home + field;
                            } else {
                                fileName = field;
                            }
                        } else {
                            fileName = directory + field;
                        }
                        int state = 0;
                        char buffer[1000];
                        FILE *fp = fopen(fileName.c_str(), "r");
                        if (fp) {
                            // file already there
                            state = 1;
                            char * getBuffer = fgets(buffer, 1000, fp);
                            if (getBuffer) {
                                // assume header there
                                state = 2;
                            }
                            fclose(fp);
                        }
                        fp = fopen(fileName.c_str(), "a");
                        if (fp) {
                            // can open - lets go for it
                            // first header if needed
                            if (state != 2) {
                                fprintf(fp, "Name,result,time,sys,elapsed,objective,continuous,tightened,cut_time,nodes,iterations,rows,columns,processed_rows,processed_columns");
                                for (int i = 0; i < statistics_number_generators; i++)
                                    fprintf(fp, ",%s", statistics_name_generators[i]);
                                fprintf(fp, ",runtime_options");
                                fprintf(fp, "\n");
                            }
                            strcpy(buffer, argv[1]);
                            char * slash = buffer;
                            for (int i = 0; i < static_cast<int>(strlen(buffer)); i++) {
                                if (buffer[i] == '/' || buffer[i] == '\\')
                                    slash = buffer + i + 1;
                            }
                            fprintf(fp, "%s,%s,%.2f,%.2f,%.2f,%.16g,%g,%g,%.2f,%d,%d,%d,%d,%d,%d",
                                    slash, statistics_result.c_str(), statistics_seconds,
                                    statistics_sys_seconds, statistics_elapsed_seconds,
                                    statistics_obj,
                                    statistics_continuous, statistics_tighter, statistics_cut_time, statistics_nodes,
                                    statistics_iterations, statistics_nrows, statistics_ncols,
                                    statistics_nprocessedrows, statistics_nprocessedcols);
                            for (int i = 0; i < statistics_number_generators; i++)
                                fprintf(fp, ",%d", statistics_number_cuts[i]);
                            fprintf(fp, ",");
                            for (int i = 1; i < argc; i++) {
                                if (strstr(argv[i], ".gz") || strstr(argv[i], ".mps"))
                                    continue;
                                if (!argv[i] || !strncmp(argv[i], "-csv", 4))
                                    break;
                                fprintf(fp, "%s ", argv[i]);
                            }
                            fprintf(fp, "\n");
                            fclose(fp);
                        } else {
                            std::cout << "Unable to open file " << fileName << std::endl;
                        }
                    }
                    break;
                    case CLP_PARAM_ACTION_SOLUTION:
                    case CLP_PARAM_ACTION_GMPL_SOLUTION:
                        if (goodModel) {
                            // get next field
                            field = CoinReadGetString(argc, argv);
			    bool append = false;
			    if (field == "append$") {
			      field = "$";
			      append = true;
			    }
                            if (field == "$") {
                                field = parameters_[iParam].stringValue();
                            } else if (field == "EOL") {
                                parameters_[iParam].printString();
                                break;
                            } else {
                                parameters_[iParam].setStringValue(field);
                            }
                            std::string fileName;
                            FILE *fp = NULL;
                            if (field == "-" || field == "EOL" || field == "stdout") {
                                // stdout
                                fp = stdout;
                            } else if (field == "stderr") {
                                // stderr
                                fp = stderr;
                            } else {
                                bool absolutePath;
                                if (dirsep == '/') {
                                    // non Windows (or cygwin)
                                    absolutePath = (field[0] == '/');
                                } else {
                                    //Windows (non cycgwin)
                                    absolutePath = (field[0] == '\\');
                                    // but allow for :
                                    if (strchr(field.c_str(), ':'))
                                        absolutePath = true;
                                }
                                if (absolutePath) {
                                    fileName = field;
                                } else if (field[0] == '~') {
                                    char * environVar = getenv("HOME");
                                    if (environVar) {
                                        std::string home(environVar);
                                        field = field.erase(0, 1);
                                        fileName = home + field;
                                    } else {
                                        fileName = field;
                                    }
                                } else {
                                    fileName = directory + field;
                                }
				if (!append)
				  fp = fopen(fileName.c_str(), "w");
				else
				  fp = fopen(fileName.c_str(), "a");
                            }
                            if (fp) {
#ifndef CBC_OTHER_SOLVER
			      // See if Glpk 
			      if (type == CLP_PARAM_ACTION_GMPL_SOLUTION) {
				int numberRows = lpSolver->getNumRows();
				int numberColumns = lpSolver->getNumCols();
				int numberGlpkRows=numberRows+1;
#ifdef COIN_HAS_GLPK
				if (cbc_glp_prob) {
				  // from gmpl
				  numberGlpkRows=glp_get_num_rows(cbc_glp_prob);
				  if (numberGlpkRows!=numberRows)
				    printf("Mismatch - cbc %d rows, glpk %d\n",
					   numberRows,numberGlpkRows);
				}
#endif
				fprintf(fp,"%d %d\n",numberGlpkRows,
					numberColumns);
				int iStat = lpSolver->status();
				int iStat2 = GLP_UNDEF;
				bool integerProblem = false;
				if (integerStatus >= 0){
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
				  fprintf(fp,"%d %g\n",iStat2,objValue);
				else
				  fprintf(fp,"%d 2 %g\n",iStat2,objValue);
				if (numberGlpkRows > numberRows) {
				  // objective as row
				  if (integerProblem) {
				    fprintf(fp,"%g\n",objValue);
				  } else {
				    fprintf(fp,"4 %g 1.0\n",objValue);
				  }
				}
				int lookup[6]=
				  {4,1,3,2,4,5};
				const double * primalRowSolution =
				  lpSolver->primalRowSolution();
				const double * dualRowSolution =
				  lpSolver->dualRowSolution();
				for (int i=0;i<numberRows;i++) {
				  if (integerProblem) {
				    fprintf(fp,"%g\n",primalRowSolution[i]);
				  } else {
				    fprintf(fp,"%d %g %g\n",lookup[lpSolver->getRowStatus(i)],
								   primalRowSolution[i],dualRowSolution[i]);
				  }
				}
				const double * primalColumnSolution =
				  lpSolver->primalColumnSolution();
				const double * dualColumnSolution =
				  lpSolver->dualColumnSolution();
				for (int i=0;i<numberColumns;i++) {
				  if (integerProblem) {
				    fprintf(fp,"%g\n",primalColumnSolution[i]);
				  } else {
				    fprintf(fp,"%d %g %g\n",lookup[lpSolver->getColumnStatus(i)],
								   primalColumnSolution[i],dualColumnSolution[i]);
				  }
				}
				fclose(fp);
#ifdef COIN_HAS_GLPK
				if (cbc_glp_prob) {
				  if (integerProblem) {
				    glp_read_mip(cbc_glp_prob,fileName.c_str());
				    glp_mpl_postsolve(cbc_glp_tran,
						      cbc_glp_prob,
						      GLP_MIP);
				  } else {
				    glp_read_sol(cbc_glp_prob,fileName.c_str());
				    glp_mpl_postsolve(cbc_glp_tran,
						      cbc_glp_prob,
						      GLP_SOL);
				  }
				  // free up as much as possible
				  glp_free(cbc_glp_prob);
				  glp_mpl_free_wksp(cbc_glp_tran);
				  cbc_glp_prob = NULL;
				  cbc_glp_tran = NULL;
				  //gmp_free_mem();
				  /* check that no memory blocks are still allocated */
				  glp_free_env();
				}
#endif
				break;
			      }
                                if (printMode < 5) {
                                    // Write solution header (suggested by Luigi Poderico)
                                    lpSolver->computeObjectiveValue(false);
                                    double objValue = lpSolver->getObjValue();
                                    int iStat = lpSolver->status();
				    int iStat2 = -1;
                                    if (integerStatus >= 0){
                                        iStat = integerStatus;
					iStat2 = babModel_->secondaryStatus();
				    }
                                    if (iStat == 0) {
                                        fprintf(fp, "Optimal" );
					if (iStat2 == 2){
					   fprintf(fp, " (within gap tolerance)" );
					}
                                    } else if (iStat == 1) {
                                        // infeasible
                                        fprintf(fp, "Infeasible" );
                                    } else if (iStat == 2) {
                                        // unbounded
                                        fprintf(fp, "Unbounded" );
                                    } else if (iStat >= 3 && iStat <= 5) {
				        if (iStat == 3) {
					    if (iStat2 == 4){
					        fprintf(fp, "Stopped on time" );
					    }else{
					        fprintf(fp, "Stopped on iterations" );
					    }
					} else if (iStat == 4){
					    fprintf(fp, "Stopped on difficulties" );
					} else {
					    fprintf(fp, "Stopped on ctrl-c" );
					}
                                        if (babModel_ && !babModel_->bestSolution())
                                            fprintf(fp, " (no integer solution - continuous used)");
                                    } else if (iStat == 6) {
				         // bab infeasible
				         fprintf(fp, "Integer infeasible" );
                                    } else {
				         fprintf(fp, "Status unknown" );
                                    }
                                    fprintf(fp, " - objective value %.8f\n", objValue);
                                }
#endif
                                // make fancy later on
                                int iRow;
                                int numberRows = clpSolver->getNumRows();
                                const double * dualRowSolution = clpSolver->getRowPrice();
                                const double * primalRowSolution =
                                    clpSolver->getRowActivity();
                                const double * rowLower = clpSolver->getRowLower();
                                const double * rowUpper = clpSolver->getRowUpper();
                                double primalTolerance ;
                                clpSolver->getDblParam(OsiPrimalTolerance, primalTolerance);
                                size_t lengthPrint = static_cast<size_t>(CoinMax(lengthName, 8));
                                bool doMask = (printMask != "" && lengthName);
                                int * maskStarts = NULL;
                                int maxMasks = 0;
                                char ** masks = NULL;
                                if (doMask) {
                                    int nAst = 0;
                                    const char * pMask2 = printMask.c_str();
                                    char pMask[100];
                                    size_t iChar;
                                    size_t lengthMask = strlen(pMask2);
                                    assert (lengthMask < 100);
                                    if (*pMask2 == '"') {
                                        if (pMask2[lengthMask-1] != '"') {
                                            printf("mismatched \" in mask %s\n", pMask2);
                                            break;
                                        } else {
                                            strcpy(pMask, pMask2 + 1);
                                            *strchr(pMask, '"') = '\0';
                                        }
                                    } else if (*pMask2 == '\'') {
                                        if (pMask2[lengthMask-1] != '\'') {
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
                                    maskStarts = new int[lengthName+2];
                                    masks = new char * [maxMasks];
                                    char ** newMasks = new char * [maxMasks];
                                    int i;
                                    for (i = 0; i < maxMasks; i++) {
                                        masks[i] = new char[lengthName+1];
                                        newMasks[i] = new char[lengthName+1];
                                    }
                                    strcpy(masks[0], pMask);
                                    for (int iAst = 0; iAst < nAst; iAst++) {
                                        int nOldEntries = nEntries;
                                        nEntries = 0;
                                        for (int iEntry = 0; iEntry < nOldEntries; iEntry++) {
                                            char * oldMask = masks[iEntry];
                                            char * ast = strchr(oldMask, '*');
                                            assert (ast);
                                            size_t length = strlen(oldMask) - 1;
                                            size_t nBefore = ast - oldMask;
                                            size_t nAfter = length - nBefore;
                                            // and add null
                                            nAfter++;
                                            for (int i = 0; i <= lengthName - static_cast<int>(length); i++) {
                                                char * maskOut = newMasks[nEntries];
                                                memcpy(maskOut, oldMask, nBefore);
                                                for (int k = 0; k < i; k++)
                                                    maskOut[k+nBefore] = '?';
                                                memcpy(maskOut + nBefore + i, ast + 1, nAfter);
                                                nEntries++;
                                                assert (nEntries <= maxMasks);
                                            }
                                        }
                                        char ** temp = masks;
                                        masks = newMasks;
                                        newMasks = temp;
                                    }
                                    // Now extend and sort
                                    int * sort = new int[nEntries];
                                    for (i = 0; i < nEntries; i++) {
                                        char * maskThis = masks[i];
                                        size_t length = strlen(maskThis);
                                        while (length > 0 && maskThis[length-1] == ' ')
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
                                    delete [] sort;
                                    for (i = 0; i < maxMasks; i++)
                                        delete [] newMasks[i];
                                    delete [] newMasks;
                                }
				if (printMode > 5) {
				  ClpSimplex * solver = clpSolver->getModelPtr();
				  int numberColumns = solver->numberColumns();
				  // column length unless rhs ranging
				  int number = numberColumns;
				  switch (printMode) {
				    // bound ranging
				  case 6:
				    fprintf(fp,"Bound ranging");
				    break;
				    // rhs ranging
				  case 7:
				    fprintf(fp,"Rhs ranging");
				    number = numberRows;
				    break;
				    // objective ranging
				  case 8:
				    fprintf(fp,"Objective ranging");
				    break;
				  }
				  if (lengthName)
				    fprintf(fp,",name");
				  fprintf(fp,",increase,variable,decrease,variable\n");
				  int * which = new int [ number];
				  if (printMode != 7) {
				    if (!doMask) {
				      for (int i = 0; i < number;i ++)
					which[i]=i;
				    } else {
				      int n = 0;
				      for (int i = 0; i < number;i ++) {
					if (maskMatches(maskStarts,masks,columnNames[i]))
					  which[n++]=i;
				      }
				      if (n) {
					number=n;
				      } else {
					printf("No names match - doing all\n");
					for (int i = 0; i < number;i ++)
					  which[i]=i;
				      }
				    }
				  } else {
				    if (!doMask) {
				    for (int i = 0; i < number;i ++)
				      which[i]=i+numberColumns;
				    } else {
				      int n = 0;
				      for (int i = 0; i < number;i ++) {
					if (maskMatches(maskStarts,masks,rowNames[i]))
					  which[n++]=i+numberColumns;
				      }
				      if (n) {
					number=n;
				      } else {
					printf("No names match - doing all\n");
					for (int i = 0; i < number;i ++)
					  which[i]=i+numberColumns;
				      }
				    }
				  }
				  double * valueIncrease = new double [ number];
				  int * sequenceIncrease = new int [ number];
				  double * valueDecrease = new double [ number];
				  int * sequenceDecrease = new int [ number];
				  switch (printMode) {
				    // bound or rhs ranging
				  case 6:
				  case 7:
				    solver->primalRanging(numberRows,
								 which, valueIncrease, sequenceIncrease,
								 valueDecrease, sequenceDecrease);
				    break;
				    // objective ranging
				  case 8:
				    solver->dualRanging(number,
							       which, valueIncrease, sequenceIncrease,
							       valueDecrease, sequenceDecrease);
				    break;
				  }
				  for (int i = 0; i < number; i++) {
				    int iWhich = which[i];
				    fprintf(fp, "%d,", (iWhich<numberColumns) ? iWhich : iWhich-numberColumns);
				    if (lengthName) {
				      const char * name = (printMode==7) ? rowNames[iWhich-numberColumns].c_str() : columnNames[iWhich].c_str();
				      fprintf(fp,"%s,",name);
				    }
				    if (valueIncrease[i]<1.0e30) {
				      fprintf(fp, "%.10g,", valueIncrease[i]);
				      int outSequence = sequenceIncrease[i];
				      if (outSequence<numberColumns) {
					if (lengthName)
					  fprintf(fp,"%s,",columnNames[outSequence].c_str());
					else
					  fprintf(fp,"C%7.7d,",outSequence);
				      } else {
					outSequence -= numberColumns;
					if (lengthName)
					  fprintf(fp,"%s,",rowNames[outSequence].c_str());
					else
					  fprintf(fp,"R%7.7d,",outSequence);
				      }
				    } else {
				      fprintf(fp,"1.0e100,,");
				    }
				    if (valueDecrease[i]<1.0e30) {
				      fprintf(fp, "%.10g,", valueDecrease[i]);
				      int outSequence = sequenceDecrease[i];
				      if (outSequence<numberColumns) {
					if (lengthName)
					  fprintf(fp,"%s",columnNames[outSequence].c_str());
					else
					  fprintf(fp,"C%7.7d",outSequence);
				      } else {
					outSequence -= numberColumns;
					if (lengthName)
					  fprintf(fp,"%s",rowNames[outSequence].c_str());
					else
					  fprintf(fp,"R%7.7d",outSequence);
				      }
				    } else {
				      fprintf(fp,"1.0e100,");
				    }
				    fprintf(fp,"\n");
				  }
				  if (fp != stdout)
				    fclose(fp);
				  delete [] which;
				  delete [] valueIncrease;
				  delete [] sequenceIncrease;
				  delete [] valueDecrease;
				  delete [] sequenceDecrease;
				  if (masks) {
				    delete [] maskStarts;
				    for (int i = 0; i < maxMasks; i++)
				      delete [] masks[i];
				    delete [] masks;
				  }
				  break;
				}
                                if (printMode > 2 && printMode < 5) {
                                    for (iRow = 0; iRow < numberRows; iRow++) {
                                        int type = printMode - 3;
                                        if (primalRowSolution[iRow] > rowUpper[iRow] + primalTolerance ||
                                                primalRowSolution[iRow] < rowLower[iRow] - primalTolerance) {
                                            fprintf(fp, "** ");
                                            type = 2;
                                        } else if (fabs(primalRowSolution[iRow]) > 1.0e-8) {
                                            type = 1;
                                        } else if (numberRows < 50) {
                                            type = 3;
                                        }
                                        if (doMask && !maskMatches(maskStarts, masks, rowNames[iRow]))
                                            type = 0;
                                        if (type) {
                                            fprintf(fp, "%7d ", iRow);
                                            if (lengthName) {
                                                const char * name = rowNames[iRow].c_str();
                                                size_t n = strlen(name);
                                                size_t i;
                                                for (i = 0; i < n; i++)
                                                    fprintf(fp, "%c", name[i]);
                                                for (; i < lengthPrint; i++)
                                                    fprintf(fp, " ");
                                            }
                                            fprintf(fp, " %15.8g        %15.8g\n", primalRowSolution[iRow],
                                                    dualRowSolution[iRow]);
                                        }
                                    }
                                }
                                int iColumn;
                                int numberColumns = clpSolver->getNumCols();
                                const double * dualColumnSolution =
                                    clpSolver->getReducedCost();
                                const double * primalColumnSolution =
                                    clpSolver->getColSolution();
                                const double * columnLower = clpSolver->getColLower();
                                const double * columnUpper = clpSolver->getColUpper();
                                if (printMode != 2) {
                                    if (printMode == 5) {
                                        if (lengthName)
                                            fprintf(fp, "name");
                                        else
                                            fprintf(fp, "number");
                                        fprintf(fp, ",solution\n");
                                    }
                                    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                                        int type = (printMode > 3) ? 1 : 0;
                                        if (primalColumnSolution[iColumn] > columnUpper[iColumn] + primalTolerance ||
                                                primalColumnSolution[iColumn] < columnLower[iColumn] - primalTolerance) {
                                            fprintf(fp, "** ");
                                            type = 2;
                                        } else if (fabs(primalColumnSolution[iColumn]) > 1.0e-8) {
                                            type = 1;
                                        } else if (numberColumns < 50) {
                                            type = 3;
                                        }
                                        // see if integer
                                        if ((!clpSolver->isInteger(iColumn) || fabs(primalColumnSolution[iColumn]) < 1.0e-8)
                                                && printMode == 1)
                                            type = 0;
                                        if (doMask && !maskMatches(maskStarts, masks,
                                                                   columnNames[iColumn]))
                                            type = 0;
                                        if (type) {
                                            if (printMode != 5) {
                                                fprintf(fp, "%7d ", iColumn);
                                                if (lengthName) {
                                                    const char * name = columnNames[iColumn].c_str();
                                                    size_t n = strlen(name);
                                                    size_t i;
                                                    for (i = 0; i < n; i++)
                                                        fprintf(fp, "%c", name[i]);
                                                    for (; i < lengthPrint; i++)
                                                        fprintf(fp, " ");
                                                }
                                                fprintf(fp, " %15.8g        %15.8g\n",
                                                        primalColumnSolution[iColumn],
                                                        dualColumnSolution[iColumn]);
                                            } else {
                                                char temp[100];
                                                if (lengthName) {
                                                    const char * name = columnNames[iColumn].c_str();
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
                                } else {
                                    // special format suitable for OsiRowCutDebugger
                                    int n = 0;
                                    bool comma = false;
                                    bool newLine = false;
                                    fprintf(fp, "\tint intIndicesV[]={\n");
                                    for (iColumn = 0; iColumn < numberColumns; iColumn++) {
                                        if (primalColumnSolution[iColumn] > 0.5 && model_.solver()->isInteger(iColumn)) {
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
                                    for ( iColumn = 0; iColumn < numberColumns; iColumn++) {
                                        if (primalColumnSolution[iColumn] > 0.5 && model_.solver()->isInteger(iColumn)) {
                                            if (comma)
                                                fprintf(fp, ",");
                                            if (newLine)
                                                fprintf(fp, "\n");
                                            int value = static_cast<int> (primalColumnSolution[iColumn] + 0.5);
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
                                }
                                if (fp != stdout)
                                    fclose(fp);
                                if (masks) {
                                    delete [] maskStarts;
                                    for (int i = 0; i < maxMasks; i++)
                                        delete [] masks[i];
                                    delete [] masks;
                                }
                            } else {
                                std::cout << "Unable to open file " << fileName << std::endl;
                            }
                        } else {
#ifndef DISALLOW_PRINTING
                            std::cout << "** Current model not valid" << std::endl;
#endif
                        }
                        break;
                    case CLP_PARAM_ACTION_SAVESOL:
                        if (goodModel) {
                            // get next field
                            field = CoinReadGetString(argc, argv);
                            if (field == "$") {
                                field = parameters_[iParam].stringValue();
                            } else if (field == "EOL") {
                                parameters_[iParam].printString();
                                break;
                            } else {
                                parameters_[iParam].setStringValue(field);
                            }
                            std::string fileName;
                            if (field[0] == '/' || field[0] == '\\') {
                                fileName = field;
                            } else if (field[0] == '~') {
                                char * environVar = getenv("HOME");
                                if (environVar) {
                                    std::string home(environVar);
                                    field = field.erase(0, 1);
                                    fileName = home + field;
                                } else {
                                    fileName = field;
                                }
                            } else {
                                fileName = directory + field;
                            }
                            saveSolution(lpSolver, fileName);
                        } else {
#ifndef DISALLOW_PRINTING
                            std::cout << "** Current model not valid" << std::endl;
#endif
                        }
                        break;
                    case CLP_PARAM_ACTION_DUMMY:
                        break;
                    case CLP_PARAM_ACTION_ENVIRONMENT:
                        CbcOrClpEnvironmentIndex = 0;
                        break;
                    case CLP_PARAM_ACTION_PARAMETRICS:
                          if (goodModel) {
                               // get next field
                               field = CoinReadGetString(argc, argv);
                               if (field == "$") {
                                    field = parameters[iParam].stringValue();
                               } else if (field == "EOL") {
                                    parameters[iParam].printString();
                                    break;
                               } else {
                                    parameters[iParam].setStringValue(field);
                               }
                               std::string fileName;
                               //bool canOpen = false;
                               if (field[0] == '/' || field[0] == '\\') {
                                    fileName = field;
                               } else if (field[0] == '~') {
                                    char * environVar = getenv("HOME");
                                    if (environVar) {
                                         std::string home(environVar);
                                         field = field.erase(0, 1);
                                         fileName = home + field;
                                    } else {
                                         fileName = field;
                                    }
                               } else {
                                    fileName = directory + field;
                               }
			       static_cast<ClpSimplexOther *> (lpSolver)->parametrics(fileName.c_str());
			       time2 = CoinCpuTime();
			       totalTime += time2 - time1;
			       time1 = time2;
                          } else {
                               std::cout << "** Current model not valid" << std::endl;
                          }
                          break;
                    default:
                        abort();
                    }
                }
            } else if (!numberMatches) {
                std::cout << "No match for " << field << " - ? for list of commands"
                          << std::endl;
            } else if (numberMatches == 1) {
                if (!numberQuery) {
                    std::cout << "Short match for " << field << " - completion: ";
                    std::cout << parameters_[firstMatch].matchName() << std::endl;
                } else if (numberQuery) {
                    std::cout << parameters_[firstMatch].matchName() << " : ";
                    std::cout << parameters_[firstMatch].shortHelp() << std::endl;
                    if (numberQuery >= 2)
                        parameters_[firstMatch].printLongHelp();
                }
            } else {
                if (!numberQuery)
                    std::cout << "Multiple matches for " << field << " - possible completions:"
                              << std::endl;
                else
                    std::cout << "Completions of " << field << ":" << std::endl;
                for ( iParam = 0; iParam < numberParameters_; iParam++ ) {
                    int match = parameters_[iParam].matches(field);
                    if (match && parameters_[iParam].displayThis()) {
                        std::cout << parameters_[iParam].matchName();
                        if (numberQuery >= 2)
                            std::cout << " : " << parameters_[iParam].shortHelp();
                        std::cout << std::endl;
                    }
                }
            }
        }
    }
#ifndef CBC_QUIET 
    sprintf(generalPrint ,
	    "Total time (CPU seconds):       %.2f   (Wallclock seconds):       %.2f\n", 
	    CoinCpuTime() - time0,
	    CoinGetTimeOfDay() - time0Elapsed);
    generalMessageHandler->message(CLP_GENERAL, generalMessages)
      << generalPrint
      << CoinMessageEol;
#endif
#ifdef COIN_HAS_GLPK
    if (cbc_glp_prob) {
      // free up as much as possible
      glp_free(cbc_glp_prob);
      glp_mpl_free_wksp(cbc_glp_tran);
      glp_free_env(); 
      cbc_glp_prob = NULL;
      cbc_glp_tran = NULL;
    }
#endif
    delete [] statistics_number_cuts;
    delete [] statistics_name_generators;
    // By now all memory should be freed
#ifdef DMALLOC
    //dmalloc_log_unfreed();
    //dmalloc_shutdown();
#endif
    if (babModel_) {
        model_.moveInfo(*babModel_);
#ifndef CBC_OTHER_SOLVER
        OsiClpSolverInterface * clpSolver0 = dynamic_cast< OsiClpSolverInterface*> (babModel_->solver());
        ClpSimplex * lpSolver0 = clpSolver0->getModelPtr();
        OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (model_.solver());
        ClpSimplex * lpSolver = clpSolver->getModelPtr();
        if (lpSolver0 != lpSolver && lpSolver != originalSolver->getModelPtr())
            lpSolver->moveInfo(*lpSolver0);
        //babModel_->setModelOwnsSolver(false);
#endif
    }
#ifdef CBC_SIG_TRAP
    // On Sun sometimes seems to be error - try and get round it
    CoinSighandler_t saveSignal = SIG_DFL;
    // register signal handler
    saveSignal = signal(SIGSEGV, signal_handler_error);
    // to force failure!babModel_->setNumberObjects(20000);
    if (!sigsetjmp(cbc_seg_buffer, 1)) {
#endif
        delete babModel_;
#ifdef CBC_SIG_TRAP
    } else {
        std::cerr << "delete babModel_ failed" << std::endl;
    }
#endif
    babModel_ = NULL;
    model_.solver()->setWarmStart(NULL);
    //sprintf(generalPrint, "Total time %.2f", CoinCpuTime() - time0);
    //generalMessageHandler->message(CLP_GENERAL, generalMessages)
    //<< generalPrint
    //<< CoinMessageEol;
    return 0;
}



int CbcMain (int argc, const char *argv[],
             CbcModel  & model)
{
    CbcMain0(model);
    return CbcMain1(argc, argv, model);
}


void CbcMain0 (CbcModel  & model)
{
#ifndef CBC_OTHER_SOLVER
    OsiClpSolverInterface * originalSolver = dynamic_cast<OsiClpSolverInterface *> (model.solver());
#elif CBC_OTHER_SOLVER==1
    OsiCpxSolverInterface * originalSolver = dynamic_cast<OsiCpxSolverInterface *> (model.solver());
    // Dummy solvers
    OsiClpSolverInterface dummySolver;
    ClpSimplex * lpSolver = dummySolver.getModelPtr();
    OsiCpxSolverInterface * clpSolver = originalSolver;
#endif
    assert (originalSolver);
    CoinMessageHandler * generalMessageHandler = originalSolver->messageHandler();
    generalMessageHandler->setPrefix(true);
#ifndef CBC_OTHER_SOLVER
    OsiSolverInterface * solver = model.solver();
    OsiClpSolverInterface * clpSolver = dynamic_cast< OsiClpSolverInterface*> (solver);
    ClpSimplex * lpSolver = clpSolver->getModelPtr();
    lpSolver->setPerturbation(50);
    lpSolver->messageHandler()->setPrefix(false);
#endif
    establishParams(numberParameters, parameters) ;
    const char dirsep =  CoinFindDirSeparator();
    std::string directory;
    std::string dirSample;
    std::string dirNetlib;
    std::string dirMiplib;
    if (dirsep == '/') {
        directory = "./";
        dirSample = "../../Data/Sample/";
        dirNetlib = "../../Data/Netlib/";
        dirMiplib = "../../Data/miplib3/";
    } else {
        directory = ".\\";
        dirSample = "..\\..\\..\\..\\Data\\Sample\\";
        dirNetlib = "..\\..\\..\\..\\Data\\Netlib\\";
        dirMiplib = "..\\..\\..\\..\\Data\\miplib3\\";
    }
    std::string defaultDirectory = directory;
    std::string importFile = "";
    std::string exportFile = "default.mps";
    std::string importBasisFile = "";
    std::string importPriorityFile = "";
    std::string debugFile = "";
    std::string printMask = "";
    std::string exportBasisFile = "default.bas";
    std::string saveFile = "default.prob";
    std::string restoreFile = "default.prob";
    std::string solutionFile = "stdout";
    std::string solutionSaveFile = "solution.file";
    int doIdiot = -1;
    int outputFormat = 2;
    int substitution = 3;
    int dualize = 3;
    int preSolve = 5;
    int doSprint = -1;
    int testOsiParameters = -1;
    parameters[whichParam(CLP_PARAM_ACTION_BASISIN, numberParameters, parameters)].setStringValue(importBasisFile);
    parameters[whichParam(CBC_PARAM_ACTION_PRIORITYIN, numberParameters, parameters)].setStringValue(importPriorityFile);
    parameters[whichParam(CLP_PARAM_ACTION_BASISOUT, numberParameters, parameters)].setStringValue(exportBasisFile);
    parameters[whichParam(CLP_PARAM_ACTION_DEBUG, numberParameters, parameters)].setStringValue(debugFile);
    parameters[whichParam(CLP_PARAM_ACTION_PRINTMASK, numberParameters, parameters)].setStringValue(printMask);
    parameters[whichParam(CLP_PARAM_ACTION_DIRECTORY, numberParameters, parameters)].setStringValue(directory);
    parameters[whichParam(CLP_PARAM_ACTION_DIRSAMPLE, numberParameters, parameters)].setStringValue(dirSample);
    parameters[whichParam(CLP_PARAM_ACTION_DIRNETLIB, numberParameters, parameters)].setStringValue(dirNetlib);
    parameters[whichParam(CBC_PARAM_ACTION_DIRMIPLIB, numberParameters, parameters)].setStringValue(dirMiplib);
    parameters[whichParam(CLP_PARAM_DBL_DUALBOUND, numberParameters, parameters)].setDoubleValue(lpSolver->dualBound());
    parameters[whichParam(CLP_PARAM_DBL_DUALTOLERANCE, numberParameters, parameters)].setDoubleValue(lpSolver->dualTolerance());
    parameters[whichParam(CLP_PARAM_ACTION_EXPORT, numberParameters, parameters)].setStringValue(exportFile);
    parameters[whichParam(CLP_PARAM_INT_IDIOT, numberParameters, parameters)].setIntValue(doIdiot);
    parameters[whichParam(CLP_PARAM_ACTION_IMPORT, numberParameters, parameters)].setStringValue(importFile);
    parameters[whichParam(CLP_PARAM_DBL_PRESOLVETOLERANCE, numberParameters, parameters)].setDoubleValue(1.0e-8);
    int slog = whichParam(CLP_PARAM_INT_SOLVERLOGLEVEL, numberParameters, parameters);
    int log = whichParam(CLP_PARAM_INT_LOGLEVEL, numberParameters, parameters);
    parameters[slog].setIntValue(1);
    clpSolver->messageHandler()->setLogLevel(1) ;
    model.messageHandler()->setLogLevel(1);
    lpSolver->setLogLevel(1);
    parameters[log].setIntValue(1);
    parameters[whichParam(CLP_PARAM_INT_MAXFACTOR, numberParameters, parameters)].setIntValue(lpSolver->factorizationFrequency());
    parameters[whichParam(CLP_PARAM_INT_MAXITERATION, numberParameters, parameters)].setIntValue(lpSolver->maximumIterations());
    parameters[whichParam(CLP_PARAM_INT_OUTPUTFORMAT, numberParameters, parameters)].setIntValue(outputFormat);
    parameters[whichParam(CLP_PARAM_INT_PRESOLVEPASS, numberParameters, parameters)].setIntValue(preSolve);
    parameters[whichParam(CLP_PARAM_INT_PERTVALUE, numberParameters, parameters)].setIntValue(lpSolver->perturbation());
    parameters[whichParam(CLP_PARAM_DBL_PRIMALTOLERANCE, numberParameters, parameters)].setDoubleValue(lpSolver->primalTolerance());
    parameters[whichParam(CLP_PARAM_DBL_PRIMALWEIGHT, numberParameters, parameters)].setDoubleValue(lpSolver->infeasibilityCost());
    parameters[whichParam(CLP_PARAM_ACTION_RESTORE, numberParameters, parameters)].setStringValue(restoreFile);
    parameters[whichParam(CLP_PARAM_ACTION_SAVE, numberParameters, parameters)].setStringValue(saveFile);
    //parameters[whichParam(CLP_PARAM_DBL_TIMELIMIT,numberParameters,parameters)].setDoubleValue(1.0e8);
    parameters[whichParam(CBC_PARAM_DBL_TIMELIMIT_BAB, numberParameters, parameters)].setDoubleValue(1.0e8);
    parameters[whichParam(CLP_PARAM_ACTION_SOLUTION, numberParameters, parameters)].setStringValue(solutionFile);
    parameters[whichParam(CLP_PARAM_ACTION_SAVESOL, numberParameters, parameters)].setStringValue(solutionSaveFile);
    parameters[whichParam(CLP_PARAM_INT_SPRINT, numberParameters, parameters)].setIntValue(doSprint);
    parameters[whichParam(CLP_PARAM_INT_SUBSTITUTION, numberParameters, parameters)].setIntValue(substitution);
    parameters[whichParam(CLP_PARAM_INT_DUALIZE, numberParameters, parameters)].setIntValue(dualize);
    model.setNumberBeforeTrust(10);
    parameters[whichParam(CBC_PARAM_INT_NUMBERBEFORE, numberParameters, parameters)].setIntValue(5);
    parameters[whichParam(CBC_PARAM_INT_MAXNODES, numberParameters, parameters)].setIntValue(model.getMaximumNodes());
    model.setNumberStrong(5);
    parameters[whichParam(CBC_PARAM_INT_STRONGBRANCHING, numberParameters, parameters)].setIntValue(model.numberStrong());
    parameters[whichParam(CBC_PARAM_DBL_INFEASIBILITYWEIGHT, numberParameters, parameters)].setDoubleValue(model.getDblParam(CbcModel::CbcInfeasibilityWeight));
    parameters[whichParam(CBC_PARAM_DBL_INTEGERTOLERANCE, numberParameters, parameters)].setDoubleValue(model.getDblParam(CbcModel::CbcIntegerTolerance));
    parameters[whichParam(CBC_PARAM_DBL_INCREMENT, numberParameters, parameters)].setDoubleValue(model.getDblParam(CbcModel::CbcCutoffIncrement));
    parameters[whichParam(CBC_PARAM_INT_TESTOSI, numberParameters, parameters)].setIntValue(testOsiParameters);
    parameters[whichParam(CBC_PARAM_INT_FPUMPTUNE, numberParameters, parameters)].setIntValue(1003);
    initialPumpTune = 1003;
#ifdef CBC_THREAD
    parameters[whichParam(CBC_PARAM_INT_THREADS, numberParameters, parameters)].setIntValue(0);
#endif
    // Set up likely cut generators and defaults
    parameters[whichParam(CBC_PARAM_STR_PREPROCESS, numberParameters, parameters)].setCurrentOption("sos");
    parameters[whichParam(CBC_PARAM_INT_MIPOPTIONS, numberParameters, parameters)].setIntValue(1057);
    parameters[whichParam(CBC_PARAM_INT_CUTPASSINTREE, numberParameters, parameters)].setIntValue(1);
    parameters[whichParam(CBC_PARAM_INT_MOREMIPOPTIONS, numberParameters, parameters)].setIntValue(-1);
    parameters[whichParam(CBC_PARAM_INT_MAXHOTITS, numberParameters, parameters)].setIntValue(100);
    parameters[whichParam(CBC_PARAM_STR_CUTSSTRATEGY, numberParameters, parameters)].setCurrentOption("on");
    parameters[whichParam(CBC_PARAM_STR_HEURISTICSTRATEGY, numberParameters, parameters)].setCurrentOption("on");
    parameters[whichParam(CBC_PARAM_STR_NODESTRATEGY, numberParameters, parameters)].setCurrentOption("fewest");
    parameters[whichParam(CBC_PARAM_STR_GOMORYCUTS, numberParameters, parameters)].setCurrentOption("ifmove");
    parameters[whichParam(CBC_PARAM_STR_PROBINGCUTS, numberParameters, parameters)].setCurrentOption("ifmove");
    parameters[whichParam(CBC_PARAM_STR_KNAPSACKCUTS, numberParameters, parameters)].setCurrentOption("ifmove");
#ifdef ZERO_HALF_CUTS
    parameters[whichParam(CBC_PARAM_STR_ZEROHALFCUTS, numberParameters, parameters)].setCurrentOption("off");
#endif
    parameters[whichParam(CBC_PARAM_STR_REDSPLITCUTS, numberParameters, parameters)].setCurrentOption("off");
    parameters[whichParam(CBC_PARAM_STR_CLIQUECUTS, numberParameters, parameters)].setCurrentOption("ifmove");
    parameters[whichParam(CBC_PARAM_STR_MIXEDCUTS, numberParameters, parameters)].setCurrentOption("ifmove");
    parameters[whichParam(CBC_PARAM_STR_FLOWCUTS, numberParameters, parameters)].setCurrentOption("ifmove");
    parameters[whichParam(CBC_PARAM_STR_TWOMIRCUTS, numberParameters, parameters)].setCurrentOption("root");
    parameters[whichParam(CBC_PARAM_STR_LANDPCUTS, numberParameters, parameters)].setCurrentOption("off");
    parameters[whichParam(CBC_PARAM_STR_RESIDCUTS, numberParameters, parameters)].setCurrentOption("off");
    parameters[whichParam(CBC_PARAM_STR_ROUNDING, numberParameters, parameters)].setCurrentOption("on");
    parameters[whichParam(CBC_PARAM_STR_FPUMP, numberParameters, parameters)].setCurrentOption("on");
    parameters[whichParam(CBC_PARAM_STR_GREEDY, numberParameters, parameters)].setCurrentOption("on");
    parameters[whichParam(CBC_PARAM_STR_COMBINE, numberParameters, parameters)].setCurrentOption("on");
    parameters[whichParam(CBC_PARAM_STR_CROSSOVER2, numberParameters, parameters)].setCurrentOption("off");
    parameters[whichParam(CBC_PARAM_STR_PIVOTANDCOMPLEMENT, numberParameters, parameters)].setCurrentOption("off");
    parameters[whichParam(CBC_PARAM_STR_PIVOTANDFIX, numberParameters, parameters)].setCurrentOption("off");
    parameters[whichParam(CBC_PARAM_STR_RANDROUND, numberParameters, parameters)].setCurrentOption("off");
    parameters[whichParam(CBC_PARAM_STR_NAIVE, numberParameters, parameters)].setCurrentOption("off");
    parameters[whichParam(CBC_PARAM_STR_RINS, numberParameters, parameters)].setCurrentOption("off");
    parameters[whichParam(CBC_PARAM_STR_DINS, numberParameters, parameters)].setCurrentOption("off");
    parameters[whichParam(CBC_PARAM_STR_RENS, numberParameters, parameters)].setCurrentOption("off");
    parameters[whichParam(CBC_PARAM_STR_LOCALTREE, numberParameters, parameters)].setCurrentOption("off");
    parameters[whichParam(CBC_PARAM_STR_COSTSTRATEGY, numberParameters, parameters)].setCurrentOption("off");
}

/*
  Routines to print statistics.
*/
static void breakdown(const char * name, int numberLook, const double * region)
{
    double range[] = {
        -COIN_DBL_MAX,
        -1.0e15, -1.0e11, -1.0e8, -1.0e5, -1.0e4, -1.0e3, -1.0e2, -1.0e1,
        -1.0,
        -1.0e-1, -1.0e-2, -1.0e-3, -1.0e-4, -1.0e-5, -1.0e-8, -1.0e-11, -1.0e-15,
        0.0,
        1.0e-15, 1.0e-11, 1.0e-8, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1,
        1.0,
        1.0e1, 1.0e2, 1.0e3, 1.0e4, 1.0e5, 1.0e8, 1.0e11, 1.0e15,
        COIN_DBL_MAX
    };
    int nRanges = static_cast<int> (sizeof(range) / sizeof(double));
    int * number = new int[nRanges];
    memset(number, 0, nRanges*sizeof(int));
    int * numberExact = new int[nRanges];
    memset(numberExact, 0, nRanges*sizeof(int));
    int i;
    for ( i = 0; i < numberLook; i++) {
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
            printf("%d between %g and %g", number[i], range[i-1], range[i]);
        if (numberExact[i]) {
            if (number[i])
                printf(", ");
            printf("%d exactly at %g", numberExact[i], range[i]);
        }
        if (number[i] + numberExact[i])
            printf("\n");
    }
    delete [] number;
    delete [] numberExact;
}
static void statistics(ClpSimplex * originalModel, ClpSimplex * model)
{
    int numberColumns = originalModel->numberColumns();
    const char * integerInformation  = originalModel->integerInformation();
    const double * columnLower = originalModel->columnLower();
    const double * columnUpper = originalModel->columnUpper();
    int numberIntegers = 0;
    int numberBinary = 0;
    int iRow, iColumn;
    if (integerInformation) {
        for (iColumn = 0; iColumn < numberColumns; iColumn++) {
            if (integerInformation[iColumn]) {
                if (columnUpper[iColumn] > columnLower[iColumn]) {
                    numberIntegers++;
                    if (columnUpper[iColumn] == 0.0 && columnLower[iColumn] == 1)
                        numberBinary++;
                }
            }
        }
    }
    numberColumns = model->numberColumns();
    int numberRows = model->numberRows();
    columnLower = model->columnLower();
    columnUpper = model->columnUpper();
    const double * rowLower = model->rowLower();
    const double * rowUpper = model->rowUpper();
    const double * objective = model->objective();
    CoinPackedMatrix * matrix = model->matrix();
    CoinBigIndex numberElements = matrix->getNumElements();
    const int * columnLength = matrix->getVectorLengths();
    //const CoinBigIndex * columnStart = matrix->getVectorStarts();
    const double * elementByColumn = matrix->getElements();
    int * number = new int[numberRows+1];
    memset(number, 0, (numberRows + 1)*sizeof(int));
    int numberObjSingletons = 0;
    /* cType
       0 0/inf, 1 0/up, 2 lo/inf, 3 lo/up, 4 free, 5 fix, 6 -inf/0, 7 -inf/up,
       8 0/1
    */
    int cType[9];
    std::string cName[] = {"0.0->inf,", "0.0->up,", "lo->inf,", "lo->up,", "free,", "fixed,", "-inf->0.0,",
                           "-inf->up,", "0.0->1.0"
                          };
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
    std::string rName[] = {"E 0.0,", "E 1.0,", "E -1.0,", "E other,", "G 0.0,", "G 1.0,", "G other,",
                           "L 0.0,", "L 1.0,", "L other,", "Range 0.0->1.0,", "Range other,", "Free"
                          };
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
    printf("\n\nProblem has %d rows, %d columns (%d with objective) and %d elements\n",
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
    for (k = 0; k < static_cast<int> (sizeof(cType) / sizeof(int)); k++) {
        printf("%d of type %s ", cType[k], cName[k].c_str());
        if (((k + 1) % 3) == 0)
            printf("\n");
    }
    if ((k % 3) != 0)
        printf("\n");
    printf("Row breakdown:\n");
    for (k = 0; k < static_cast<int> (sizeof(rType) / sizeof(int)); k++) {
        printf("%d of type %s ", rType[k], rName[k].c_str());
        if (((k + 1) % 3) == 0)
            printf("\n");
    }
    if ((k % 3) != 0)
        printf("\n");
    if (model->logLevel() < 2)
        return ;
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
    delete [] number;
    printf("\n\n");
    // get row copy
    CoinPackedMatrix rowCopy = *matrix;
    rowCopy.reverseOrdering();
    //const int * column = rowCopy.getIndices();
    const int * rowLength = rowCopy.getVectorLengths();
    //const CoinBigIndex * rowStart = rowCopy.getVectorStarts();
    //const double * element = rowCopy.getElements();
    number = new int[numberColumns+1];
    memset(number, 0, (numberColumns + 1)*sizeof(int));
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
    delete [] number;
    // Now do breakdown of ranges
    breakdown("Elements", numberElements, elementByColumn);
    breakdown("RowLower", numberRows, rowLower);
    breakdown("RowUpper", numberRows, rowUpper);
    breakdown("ColumnLower", numberColumns, columnLower);
    breakdown("ColumnUpper", numberColumns, columnUpper);
    breakdown("Objective", numberColumns, objective);
}


static bool maskMatches(const int * starts, char ** masks,
                        std::string & check)
{
    // back to char as I am old fashioned
    const char * checkC = check.c_str();
    size_t length = strlen(checkC);
    while (length > 0 && checkC[length-1] == ' ')
        length--;
    for (int i = starts[length]; i < starts[length+1]; i++) {
        char * thisMask = masks[i];
        size_t k;
        for ( k = 0; k < length; k++) {
            if (thisMask[k] != '?' && thisMask[k] != checkC[k])
                break;
        }
        if (k == length)
            return true;
    }
    return false;
}

static void clean(char * temp)
{
    char * put = temp;
    while (*put >= ' ')
        put++;
    *put = '\0';
}

static void generateCode(CbcModel * /*model*/, const char * fileName, int type, int preProcess)
{
    // options on code generation
    bool sizecode = (type & 4) != 0;
    type &= 3;
    FILE * fp = fopen(fileName, "r");
    assert (fp);
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
        strcpy(line[numberLines++], "0#include \"CglProbing.hpp\""); // possibly redundant
    // To allow generated 5's to be just before branchAndBound - do rest here
    strcpy(line[numberLines++], "5  cbcModel->initialSolve();");
    strcpy(line[numberLines++], "5  if (clpModel->tightenPrimalBounds()!=0) {");
    strcpy(line[numberLines++], "5    std::cout<<\"Problem is infeasible - tightenPrimalBounds!\"<<std::endl;");
    strcpy(line[numberLines++], "5    exit(1);");
    strcpy(line[numberLines++], "5  }");
    strcpy(line[numberLines++], "5  clpModel->dual();  // clean up");
    if (sizecode) {
        // override some settings
        strcpy(line[numberLines++], "5  // compute some things using problem size");
        strcpy(line[numberLines++], "5  cbcModel->setMinimumDrop(CoinMin(5.0e-2,");
        strcpy(line[numberLines++], "5       fabs(cbcModel->getMinimizationObjValue())*1.0e-3+1.0e-4));");
        strcpy(line[numberLines++], "5  if (cbcModel->getNumCols()<500)");
        strcpy(line[numberLines++], "5    cbcModel->setMaximumCutPassesAtRoot(-100); // always do 100 if possible");
        strcpy(line[numberLines++], "5  else if (cbcModel->getNumCols()<5000)");
        strcpy(line[numberLines++], "5    cbcModel->setMaximumCutPassesAtRoot(100); // use minimum drop");
        strcpy(line[numberLines++], "5  else");
        strcpy(line[numberLines++], "5    cbcModel->setMaximumCutPassesAtRoot(20);");
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
        strcpy(line[numberLines++], "5  OsiSolverInterface * saveSolver=cbcModel->solver()->clone();");
        strcpy(line[numberLines++], "5  // Tell solver we are in Branch and Cut");
        strcpy(line[numberLines++], "5  saveSolver->setHintParam(OsiDoInBranchAndCut,true,OsiHintDo) ;");
        strcpy(line[numberLines++], "5  // Default set of cut generators");
        strcpy(line[numberLines++], "5  CglProbing generator1;");
        strcpy(line[numberLines++], "5  generator1.setUsingObjective(1);");
        strcpy(line[numberLines++], "5  generator1.setMaxPass(3);");
        strcpy(line[numberLines++], "5  generator1.setMaxProbeRoot(saveSolver->getNumCols());");
        strcpy(line[numberLines++], "5  generator1.setMaxElements(100);");
        strcpy(line[numberLines++], "5  generator1.setMaxLookRoot(50);");
        strcpy(line[numberLines++], "5  generator1.setRowCuts(3);");
        strcpy(line[numberLines++], "5  // Add in generators");
        strcpy(line[numberLines++], "5  process.addCutGenerator(&generator1);");
        strcpy(line[numberLines++], "5  process.messageHandler()->setLogLevel(cbcModel->logLevel());");
        strcpy(line[numberLines++], "5  OsiSolverInterface * solver2 = ");
        sprintf(line[numberLines++], "5    process.preProcessNonDefault(*saveSolver,%d,10);", translate[preProcess]);
        strcpy(line[numberLines++], "5  // Tell solver we are not in Branch and Cut");
        strcpy(line[numberLines++], "5  saveSolver->setHintParam(OsiDoInBranchAndCut,false,OsiHintDo) ;");
        strcpy(line[numberLines++], "5  if (solver2)");
        strcpy(line[numberLines++], "5    solver2->setHintParam(OsiDoInBranchAndCut,false,OsiHintDo) ;");
        strcpy(line[numberLines++], "5  if (!solver2) {");
        strcpy(line[numberLines++], "5    std::cout<<\"Pre-processing says infeasible!\"<<std::endl;");
        strcpy(line[numberLines++], "5    exit(1);");
        strcpy(line[numberLines++], "5  } else {");
        strcpy(line[numberLines++], "5    std::cout<<\"processed model has \"<<solver2->getNumRows()");
        strcpy(line[numberLines++], "5	     <<\" rows, \"<<solver2->getNumCols()");
        strcpy(line[numberLines++], "5	     <<\" columns and \"<<solver2->getNumElements()");
        strcpy(line[numberLines++], "5	     <<\" elements\"<<solver2->getNumElements()<<std::endl;");
        strcpy(line[numberLines++], "5  }");
        strcpy(line[numberLines++], "5  // we have to keep solver2 so pass clone");
        strcpy(line[numberLines++], "5  solver2 = solver2->clone();");
        strcpy(line[numberLines++], "5  cbcModel->assignSolver(solver2);");
        strcpy(line[numberLines++], "5  cbcModel->initialSolve();");
    }
    while (fgets(line[numberLines], MAXONELINE, fp)) {
        assert (numberLines < MAXLINES);
        clean(line[numberLines]);
        numberLines++;
    }
    fclose(fp);
    strcpy(line[numberLines++], "0\nint main (int argc, const char *argv[])\n{");
    strcpy(line[numberLines++], "0  OsiClpSolverInterface solver1;");
    strcpy(line[numberLines++], "0  int status=1;");
    strcpy(line[numberLines++], "0  if (argc<2)");
    strcpy(line[numberLines++], "0    std::cout<<\"Please give file name\"<<std::endl;");
    strcpy(line[numberLines++], "0  else");
    strcpy(line[numberLines++], "0    status=solver1.readMps(argv[1],\"\");");
    strcpy(line[numberLines++], "0  if (status) {");
    strcpy(line[numberLines++], "0    std::cout<<\"Bad readMps \"<<argv[1]<<std::endl;");
    strcpy(line[numberLines++], "0    exit(1);");
    strcpy(line[numberLines++], "0  }\n");
    strcpy(line[numberLines++], "0  double time1 = CoinCpuTime();");
    strcpy(line[numberLines++], "0  CbcModel model(solver1);");
    strcpy(line[numberLines++], "0  // Now do requested saves and modifications");
    strcpy(line[numberLines++], "0  CbcModel * cbcModel = & model;");
    strcpy(line[numberLines++], "0  OsiSolverInterface * osiModel = model.solver();");
    strcpy(line[numberLines++], "0  OsiClpSolverInterface * osiclpModel = dynamic_cast< OsiClpSolverInterface*> (osiModel);");
    strcpy(line[numberLines++], "0  ClpSimplex * clpModel = osiclpModel->getModelPtr();");
    // add in comments about messages
    strcpy(line[numberLines++], "3  // You can save some time by switching off message building");
    strcpy(line[numberLines++], "3  // clpModel->messagesPointer()->setDetailMessages(100,10000,(int *) NULL);");
    // add in actual solve
    strcpy(line[numberLines++], "5  cbcModel->branchAndBound();");
    strcpy(line[numberLines++], "8  std::cout<<argv[1]<<\" took \"<<CoinCpuTime()-time1<<\" seconds, \"");
    strcpy(line[numberLines++], "8	   <<cbcModel->getNodeCount()<<\" nodes with objective \"");
    strcpy(line[numberLines++], "8	   <<cbcModel->getObjValue()");
    strcpy(line[numberLines++], "8	   <<(!cbcModel->status() ? \" Finished\" : \" Not finished\")");
    strcpy(line[numberLines++], "8	   <<std::endl;");
    strcpy(line[numberLines++], "5  // For best solution");
    strcpy(line[numberLines++], "5  int numberColumns = solver1.getNumCols();");
    strcpy(line[numberLines++], "5  if (cbcModel->getMinimizationObjValue()<1.0e50) {");
    if (preProcess > 0) {
        strcpy(line[numberLines++], "5    // post process");
        strcpy(line[numberLines++], "5    process.postProcess(*cbcModel->solver());");
        strcpy(line[numberLines++], "5    // Solution now back in saveSolver");
        strcpy(line[numberLines++], "5    cbcModel->assignSolver(saveSolver);");
        strcpy(line[numberLines++], "5    memcpy(cbcModel->bestSolution(),cbcModel->solver()->getColSolution(),");
        strcpy(line[numberLines++], "5	   numberColumns*sizeof(double));");
    }
    strcpy(line[numberLines++], "5    // put back in original solver");
    strcpy(line[numberLines++], "5    solver1.setColSolution(cbcModel->bestSolution());");
    strcpy(line[numberLines++], "5    const double * solution = solver1.getColSolution();");
    strcpy(line[numberLines++], "8  \n  // Now you would use solution etc etc\n");
    strcpy(line[numberLines++], "5");
    strcpy(line[numberLines++], "5    // Get names from solver1 (as OsiSolverInterface may lose)");
    strcpy(line[numberLines++], "5    std::vector<std::string> columnNames = *solver1.getModelPtr()->columnNames();");
    strcpy(line[numberLines++], "5    ");
    strcpy(line[numberLines++], "5    int iColumn;");
    strcpy(line[numberLines++], "5    std::cout<<std::setiosflags(std::ios::fixed|std::ios::showpoint)<<std::setw(14);");
    strcpy(line[numberLines++], "5    ");
    strcpy(line[numberLines++], "5    std::cout<<\"--------------------------------------\"<<std::endl;");
    strcpy(line[numberLines++], "5    for (iColumn=0;iColumn<numberColumns;iColumn++) {");
    strcpy(line[numberLines++], "5      double value=solution[iColumn];");
    strcpy(line[numberLines++], "5      if (fabs(value)>1.0e-7&&solver1.isInteger(iColumn)) ");
    strcpy(line[numberLines++], "5	std::cout<<std::setw(6)<<iColumn<<\" \"");
    strcpy(line[numberLines++], "5                 <<columnNames[iColumn]<<\" \"");
    strcpy(line[numberLines++], "5                 <<value<<std::endl;");
    strcpy(line[numberLines++], "5    }");
    strcpy(line[numberLines++], "5    std::cout<<\"--------------------------------------\"<<std::endl;");
    strcpy(line[numberLines++], "5  ");
    strcpy(line[numberLines++], "5    std::cout<<std::resetiosflags(std::ios::fixed|std::ios::showpoint|std::ios::scientific);");
    strcpy(line[numberLines++], "5  }");
    strcpy(line[numberLines++], "8  return 0;\n}");
    fp = fopen(fileName, "w");
    assert (fp);

    int wanted[9];
    memset(wanted, 0, sizeof(wanted));
    wanted[0] = wanted[3] = wanted[5] = wanted[8] = 1;
    if (type > 0)
        wanted[1] = wanted[6] = 1;
    if (type > 1)
        wanted[2] = wanted[4] = wanted[7] = 1;
    std::string header[9] = { "", "Save values", "Redundant save of default values", "Set changed values",
                              "Redundant set default values", "Solve", "Restore values", "Redundant restore values", "Finish up"
                            };
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
                if (!strstr(line[iLine], "save") || (!strstr(line[iLine], "clpMo") &&
                                                     !strstr(line[iLine], "_Osi")))
                    fprintf(fp, "%s\n", line[iLine] + 1);
            }
        }
    }
    fclose(fp);
    printf("C++ file written to %s\n", fileName);
}
/*
  Version 1.00.00 November 16 2005.
  This is to stop me (JJF) messing about too much.
  Tuning changes should be noted here.
  The testing next version may be activated by CBC_NEXT_VERSION
  This applies to OsiClp, Clp etc
  Version 1.00.01 November 24 2005
  Added several classes for advanced users.  This can't affect code (if you don't use it)
  Made some tiny changes (for N way branching) which should not change anything.
  CbcNWay object class - for N way branching this also allows use of CbcConsequence class.
  CbcBranchAllDifferent object class - for branching on general integer variables
  to stop them having same value so branches are x >= y+1 and x <= y-1.
  Added two new Cgl classes - CglAllDifferent which does column fixing (too slowly)
  and CglStored which just has a list of cuts which can be activated.
  Modified preprocess option to SOS
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
