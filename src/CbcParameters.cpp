/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/

#include "CoinPragma.hpp"

#include <cassert>

#include "CbcParameters.hpp"
#include "CbcParamUtils.hpp"

/*
  Constructor for parameters class.

  Set up defaults. Note that prototypes for
  cut generators and heuristics will be created on demand; see the access
  functions.

  Once this structure settles down, simple intialisation should move up to
  the standard `:' block. In the meantime, this avoids complaints about
  ordering.
*/

CbcParameters::CbcParameters() : parameters_(CbcParam::LASTPARAM), model_(0)
{

   init(DefaultStrategy);
      
}

//###########################################################################
//###########################################################################

CbcParameters::CbcParameters(int strategy) :
   parameters_(CbcParam::LASTPARAM), model_(0)
{

   init(strategy);
}

//###########################################################################
//###########################################################################

void CbcParameters::init(int strategy){
   
  for (int i = 0; i < parameters_.size(); i++){
     parameters_[i] = new CbcParam();
  }

  dfltDirectory_ = "";
  lastMpsIn_ = "";
  allowImportErrors_ = false;
  lastSolnOut_ = "stdout";
  printMode_ = 0;
  printMask_ = "";
  noPrinting_ = false;
  printWelcome_ = true;
  useSignalHandler_ = false;

  verbose_ = 0;
  paramsProcessed_ = 0;
  defaultSettings_ = true;

  debugCreate_ = "";
  debugFile_ = "";
  debugSol_.numCols_ = -1;
  debugSol_.values_ = 0;

  printOpt_ = 0;

  /*
      Assigning us_en to cur_lang_ is entirely bogus, but CoinMessages::Language
      does not provide an `unspecified' code.
    */
  msgHandler_ = new CoinMessageHandler();
  ourMsgHandler_ = true;
  cur_lang_ = CoinMessages::us_en;
  msgs_ = 0;
  logLvl_ = 0;

  totalTime_ = 0.0;

  model_ = 0;
  dfltSolver_ = 0;
  goodModel_ = false;
  bab_.majorStatus_ = CbcParameters::BACInvalid;
  bab_.minorStatus_ = CbcParameters::BACmInvalid;
  bab_.where_ = CbcParameters::BACwInvalid;
  bab_.haveAnswer_ = false;
  bab_.answerSolver_ = 0;

  preProcess_ = CbcParameters::IPPSOS;
  cutDepth_ = -1;

  probing_.mode_ = CbcParameters::CGIfMove;
  probing_.proto_ = 0;
  probing_.usingObjective_ = true;
  probing_.maxPass_ = 3;
  probing_.maxPassRoot_ = 3;
  probing_.maxProbe_ = 10;
  probing_.maxProbeRoot_ = 50;
  probing_.maxLook_ = 10;
  probing_.maxLookRoot_ = 50;
  probing_.maxElements_ = 200;
  probing_.rowCuts_ = 3;

  clique_.mode_ = CbcParameters::CGIfMove;
  clique_.proto_ = 0;
  clique_.starCliqueReport_ = false;
  clique_.rowCliqueReport_ = false;
  clique_.minViolation_ = 0.1;

  flow_.mode_ = CbcParameters::CGIfMove;
  flow_.proto_ = 0;

  gomory_.mode_ = CbcParameters::CGIfMove;
  gomory_.proto_ = 0;
  gomory_.limit_ = 50;
  gomory_.limitAtRoot_ = 512;

  knapsack_.mode_ = CbcParameters::CGIfMove;
  knapsack_.proto_ = 0;

  // landp_mode_ = CbcParameters::CGOff ;
  // landp_.proto_ = 0 ;

  mir_.mode_ = CbcParameters::CGIfMove;
  mir_.proto_ = 0;

#if 0
  oddHole_.mode_ = CbcParameters::CGOff;
  oddHole_.proto_ = 0;
#endif

  redSplit_.mode_ = CbcParameters::CGRoot;
  redSplit_.proto_ = 0;

  twomir_.mode_ = CbcParameters::CGRoot;
  twomir_.proto_ = 0;
  twomir_.maxElements_ = 250;

  fpump_.mode_ = CbcParameters::HeurOn;
  fpump_.proto_ = 0;
  fpump_.initialTune_ = -1;

  combine_.mode_ = CbcParameters::HeurOn;
  combine_.proto_ = 0;
  combine_.trySwap_ = 1;

  greedyCover_.mode_ = CbcParameters::HeurOn;
  greedyCover_.proto_ = 0;
  greedyEquality_.mode_ = CbcParameters::HeurOn;
  greedyEquality_.proto_ = 0;

  localTree_.mode_ = CbcParameters::HeurOff;
  localTree_.proto_ = 0;
  localTree_.soln_ = 0;
  localTree_.range_ = 10;
  localTree_.typeCuts_ = 0;
  localTree_.maxDiverge_ = 0;
  localTree_.timeLimit_ = 10000;
  localTree_.nodeLimit_ = 2000;
  localTree_.refine_ = true;

  rounding_.mode_ = CbcParameters::HeurOn;
  rounding_.proto_ = 0;

  djFix_.mode_ = CbcParameters::ParamOff;
  djFix_.threshold_ = 1.0e100;

  artVar_.mode_ = CbcParameters::ParamOff;
  artVar_.threshold_ = 0.0;

  priorityMode_ = CbcParameters::BPOff;
  /*
      The value for numBeforeTrust is as recommended by Achterberg. Cbc's
      implementation doesn't really have a parameter equivalent to Achterberg's
      dynamic limit on number of strong branching evaluations, so go with a
     fairly large default. As of 06.12.16, the magic number for shadow price
     mode meant `use shadow prices (penalties, I think) if there's no strong
     branching info'.
    */
  chooseStrong_.numBeforeTrust_ = 8;
  chooseStrong_.numStrong_ = 100;
  chooseStrong_.shadowPriceMode_ = 1;

  addCbcParams();
  addCbcModelParams();
  setDefaults(strategy);

  return;
}

//###########################################################################
//###########################################################################

/*
  Note that we don't want to delete dfltSolver_ here because it's just a
  copy of the pointer held in the solvers map over in CbcGenSolvers.cpp.
*/
CbcParameters::~CbcParameters() {
   // TODO Do we own pointer here?
   for (int i = 0; i < parameters_.size(); i++){
      delete parameters_[i];
   }

#ifndef CBC_CLUMSY_CODING
   if (model_)
    delete model_;
#endif
  if (bab_.answerSolver_)
    delete bab_.answerSolver_;

  if (probing_.proto_)
    delete probing_.proto_;
  if (clique_.proto_)
    delete clique_.proto_;
  if (flow_.proto_)
    delete flow_.proto_;
  if (gomory_.proto_)
    delete gomory_.proto_;
  if (knapsack_.proto_)
    delete knapsack_.proto_;
  if (mir_.proto_)
    delete mir_.proto_;
#if 0
  if (oddHole_.proto_)
    delete oddHole_.proto_;
#endif
  if (redSplit_.proto_)
    delete redSplit_.proto_;
  if (twomir_.proto_)
    delete twomir_.proto_;

  if (fpump_.proto_)
    delete fpump_.proto_;
  if (combine_.proto_)
    delete combine_.proto_;
  if (greedyCover_.proto_)
    delete greedyCover_.proto_;
  if (greedyEquality_.proto_)
    delete greedyEquality_.proto_;
  if (rounding_.proto_)
    delete rounding_.proto_;

  if (msgHandler_ && ourMsgHandler_)
    delete msgHandler_;
  if (msgs_)
    delete msgs_;

  return;
}

//###########################################################################
//###########################################################################

int CbcParameters::matches(std::string field, int &numberMatches){
   int firstMatch = -1;
   for (int iParam = 0; iParam < (int)parameters_.size(); iParam++) {
      int match = parameters_[iParam]->matches(field);
      if (match == 1) {
         numberMatches = 1;
         return iParam;
      } else {
         if (match){
            if (firstMatch < 0){
               firstMatch = iParam;
            }
            numberMatches++;
         }
      }
   }
   return firstMatch < 0 ? CbcParam::INVALID : firstMatch;
}   

//###########################################################################
//###########################################################################

// some help strings that repeat for many options
#define CUTS_LONGHELP                                                          \
  "Value 'on' enables the cut generator and CBC will try it in the branch "    \
  "and cut tree (see cutDepth on how to fine tune the behavior). Value "       \
  "'root' lets CBC run the cut generator generate only at the root node. "     \
  "Value 'ifmove' lets CBC use the cut generator in the tree if it looks as "  \
  "if it is doing some good and moves the objective value. Value 'forceon' "   \
  "turns on the cut generator and forces CBC to use it at every node."

#define HEURISTICS_LONGHELP                                                    \
  "Value 'on' means to use the heuristic in each node of the tree, i.e. "      \
  "after preprocessing. Value 'before' means use the heuristic only if "       \
  "option doHeuristics is used. Value 'both' means to use the heuristic if "   \
  "option doHeuristics is used and during solve."

//###########################################################################
//###########################################################################

/*
  Function to add Cbc parameters to the Cbc parameter
  vector. Where needed, defaults are drawn from CbcParameters.
  This function is a friend of CbcParameters.
*/

void CbcParameters::addCbcParams() {

  addCbcSolverStrParams();
  addCbcSolverHelpParams();
  addCbcSolverActionParams();
  addCbcSolverKwdParams();
  addCbcSolverDblParams();
  addCbcSolverIntParams();
  addCbcSolverBoolParams();
  addCbcSolverCutParams();
  addCbcSolverHeurParams();
  addCbcModelParams();

  for (int code = CbcParam::FIRSTPARAM + 1; code < CbcParam::LASTPARAM;
       code++) {
    getParam(code)->setParameters(this);
    getParam(code)->setModel(model_);
    getParam(code)->setParamCode(code);
  }

  return;
}

//###########################################################################
//###########################################################################

void CbcParameters::setDefaults(int strategy) {

   for (int code = CbcParam::FIRSTSTRINGPARAM + 1;
       code < CbcParam::LASTSTRINGPARAM; code++) {
      getParam(code)->setDefault(dfltDirectory_);
  }

  // Now set up  parameters according to overall strategies
  switch (strategy) {
   case CbcParameters::DefaultStrategy:
     parameters_[CbcParam::DEBUG]->setDefault("");
     parameters_[CbcParam::GMPLSOLFILE]->setDefault(std::string("gmpl.sol"));
     parameters_[CbcParam::MODELFILE]->setDefault(std::string("prob.mod"));
     parameters_[CbcParam::NEXTSOLFILE]->setDefault(std::string("next.sol"));
     parameters_[CbcParam::PRINTMASK]->setDefault("");
     parameters_[CbcParam::SOLUTIONFILE]->setDefault(std::string("opt.sol"));
     parameters_[CbcParam::PRIORITYIN]->setDefault(std::string("priorities.txt"));
     parameters_[CbcParam::COMMANDPRINTLEVEL]->setDefault("high");
     parameters_[CbcParam::CLQSTRENGTHENING]->setDefault("after");
     parameters_[CbcParam::BRANCHPRIORITY]->setDefault("off");
     parameters_[CbcParam::CUTOFFCONSTRAINT]->setDefault("off");
     parameters_[CbcParam::INTPRINT]->setDefault("normal");
     parameters_[CbcParam::NODESTRATEGY]->setDefault("hybrid");
     parameters_[CbcParam::ORBITAL]->setDefault("off");
     parameters_[CbcParam::PREPROCESS]->setDefault("off");
     parameters_[CbcParam::SOSPRIORITIZE]->setDefault("off");
     parameters_[CbcParam::STRATEGY]->setDefault("default");
     parameters_[CbcParam::TIMEMODE]->setDefault("cpu");
     parameters_[CbcParam::USECGRAPH]->setDefault("on");
     parameters_[CbcParam::ARTIFICIALCOST]->setDefault(getArtVarThreshold());
     parameters_[CbcParam::DEXTRA3]->setDefault(0.0);
     parameters_[CbcParam::DEXTRA4]->setDefault(0.0);
     parameters_[CbcParam::DEXTRA5]->setDefault(0.0);
     parameters_[CbcParam::DJFIX]->setDefault(getDjFixThreshold());
     parameters_[CbcParam::FAKECUTOFF]->setDefault(0.0);
     parameters_[CbcParam::FAKEINCREMENT]->setDefault(0.0);
     parameters_[CbcParam::SMALLBAB]->setDefault(0.5);
     parameters_[CbcParam::TIGHTENFACTOR]->setDefault(0.0);
     parameters_[CbcParam::BKPIVOTINGSTRATEGY]->setDefault(3);
     parameters_[CbcParam::BKMAXCALLS]->setDefault(1000);
     parameters_[CbcParam::BKCLQEXTMETHOD]->setDefault(4);
     parameters_[CbcParam::CPP]->setDefault(0);
     parameters_[CbcParam::CUTDEPTH]->setDefault(getCutDepth());
     parameters_[CbcParam::CUTLENGTH]->setDefault(-1);
     parameters_[CbcParam::CUTPASSINTREE]->setDefault(10);
     parameters_[CbcParam::DEPTHMINIBAB]->setDefault(-1);
     parameters_[CbcParam::DIVEOPT]->setDefault(-1);
     parameters_[CbcParam::DIVEOPTSOLVES]->setDefault(100);
     parameters_[CbcParam::DUMMY]->setDefault(0);
     parameters_[CbcParam::EXPERIMENT]->setDefault(0);
     parameters_[CbcParam::EXTRA1]->setDefault(-1);
     parameters_[CbcParam::EXTRA2]->setDefault(-1);
     parameters_[CbcParam::EXTRA3]->setDefault(-1);
     parameters_[CbcParam::EXTRA4]->setDefault(-1);
     parameters_[CbcParam::EXTRAVARIABLES]->setDefault(0);
     parameters_[CbcParam::FPUMPITS]->setDefault(getFeasPumpIters());
     parameters_[CbcParam::FPUMPTUNE]->setDefault(0);
     parameters_[CbcParam::FPUMPTUNE2]->setDefault(0);
     parameters_[CbcParam::HEUROPTIONS]->setDefault(0);
     parameters_[CbcParam::LOGLEVEL]->setDefault(getLogLevel());
     parameters_[CbcParam::LPLOGLEVEL]->setDefault(getLpLogLevel());
     parameters_[CbcParam::MAXHOTITS]->setDefault(0);
     parameters_[CbcParam::MAXSAVEDSOLS]->setDefault(1);
     parameters_[CbcParam::MAXSLOWCUTS]->setDefault(10);
     parameters_[CbcParam::MOREMOREMIPOPTIONS]->setDefault(0);
     parameters_[CbcParam::MULTIPLEROOTS]->setDefault(0);
     parameters_[CbcParam::ODDWEXTMETHOD]->setDefault(2);
     parameters_[CbcParam::OUTPUTFORMAT]->setDefault(2);
     parameters_[CbcParam::PRINTOPTIONS]->setDefault(0);
     parameters_[CbcParam::PROCESSTUNE]->setDefault(0);
     parameters_[CbcParam::RANDOMSEED]->setDefault(-1);
     parameters_[CbcParam::STRONGSTRATEGY]->setDefault(0);
     parameters_[CbcParam::TESTOSI]->setDefault(-1);
#ifdef CBC_THREAD
     parameters_[CbcParam::THREADS]->setDefault(0);
#endif
     parameters_[CbcParam::USERCBC]->setDefault(0);
     parameters_[CbcParam::VERBOSE]->setDefault(verbose_);
     parameters_[CbcParam::VUBTRY]->setDefault(-1);
     parameters_[CbcParam::CPX]->setDefault("off");
     parameters_[CbcParam::DOHEURISTIC]->setDefault("off");
     parameters_[CbcParam::ERRORSALLOWED]->setDefault("off");
     parameters_[CbcParam::MESSAGES]->setDefault("off");
     parameters_[CbcParam::PREPROCNAMES]->setDefault("off");
     parameters_[CbcParam::SOS]->setDefault("off");
     parameters_[CbcParam::USESOLUTION]->setDefault("off");
     parameters_[CbcParam::CUTSTRATEGY]->setDefault("off");
     parameters_[CbcParam::FLOWCUTS]->setDefault("off");
     parameters_[CbcParam::GMICUTS]->setDefault("off");
     parameters_[CbcParam::GOMORYCUTS]->setDefault("off");
     parameters_[CbcParam::KNAPSACKCUTS]->setDefault("off");
     parameters_[CbcParam::LAGOMORYCUTS]->setDefault("off");
     parameters_[CbcParam::LANDPCUTS]->setDefault("off");
     parameters_[CbcParam::LATWOMIRCUTS]->setDefault("off");
     parameters_[CbcParam::MIRCUTS]->setDefault("off");
     parameters_[CbcParam::ODDWHEELCUTS]->setDefault("off");
     parameters_[CbcParam::PROBINGCUTS]->setDefault("off");
     parameters_[CbcParam::REDSPLITCUTS]->setDefault("off");
     parameters_[CbcParam::REDSPLIT2CUTS]->setDefault("off");
     parameters_[CbcParam::RESIDCAPCUTS]->setDefault("off");
     parameters_[CbcParam::TWOMIRCUTS]->setDefault("off");
     parameters_[CbcParam::ZEROHALFCUTS]->setDefault("off");
     parameters_[CbcParam::COMBINE]->setDefault("off");
     parameters_[CbcParam::CROSSOVER]->setDefault("off");
     parameters_[CbcParam::DINS]->setDefault("off");
     parameters_[CbcParam::DIVINGC]->setDefault("off");
     parameters_[CbcParam::DIVINGF]->setDefault("off");
     parameters_[CbcParam::DIVINGG]->setDefault("off");
     parameters_[CbcParam::DIVINGL]->setDefault("off");
     parameters_[CbcParam::DIVINGP]->setDefault("off");
     parameters_[CbcParam::DIVINGS]->setDefault("off");
     parameters_[CbcParam::DIVINGV]->setDefault("off");
     parameters_[CbcParam::DW]->setDefault("off");
     parameters_[CbcParam::FPUMP]->setDefault("off");
     parameters_[CbcParam::GREEDY]->setDefault("off");
     parameters_[CbcParam::HEURISTICSTRATEGY]->setDefault("off");
     parameters_[CbcParam::LOCALTREE]->setDefault("off");
     parameters_[CbcParam::NAIVE]->setDefault("off");
     parameters_[CbcParam::PIVOTANDFIX]->setDefault("off");
#if 0
     parameters_[CbcParam::PIVOTANDCOMPLEMENT]->setDefault("off");
#endif
     parameters_[CbcParam::PROXIMITY]->setDefault("off");
     parameters_[CbcParam::RANDROUND]->setDefault("off");
     parameters_[CbcParam::RENS]->setDefault("off");
     parameters_[CbcParam::RINS]->setDefault("off");
     parameters_[CbcParam::ROUNDING]->setDefault("off");
     parameters_[CbcParam::VND]->setDefault("off");
     parameters_[CbcParam::ALLOWABLEGAP]->setDefault(1.0e-12);
     parameters_[CbcParam::CUTOFF]->setDefault(1.0e50);
     parameters_[CbcParam::DIRECTION]->setDefault("minimize");
     parameters_[CbcParam::INCREMENT]->setDefault(1.0e-4);
     parameters_[CbcParam::INFEASIBILITYWEIGHT]->setDefault(0.0);
     parameters_[CbcParam::INTEGERTOLERANCE]->setDefault(1.0e-6);
     parameters_[CbcParam::LOGLEVEL]->setDefault(1);
     parameters_[CbcParam::MAXIMIZE]->setType(CoinParam::paramAct);
     parameters_[CbcParam::MAXNODES]->setDefault(COIN_INT_MAX);
     parameters_[CbcParam::MAXNODESNOTIMPROVING]->setDefault(COIN_INT_MAX);
     parameters_[CbcParam::MAXSECONDSNOTIMPROVING]->setDefault(COIN_DBL_MAX);
     parameters_[CbcParam::MAXSOLS]->setDefault(COIN_INT_MAX);
     parameters_[CbcParam::MINIMIZE]->setType(CoinParam::paramAct);
     parameters_[CbcParam::MIPOPTIONS]->setDefault(0);
     parameters_[CbcParam::MOREMIPOPTIONS]->setDefault(0);
#if 0
     parameters_[CbcParam::NUMBERMINI]->setDefault(0);
#endif
     parameters_[CbcParam::NUMBERANALYZE]->setDefault(0);
     parameters_[CbcParam::REVERSE]->setType(CoinParam::paramAct);
     parameters_[CbcParam::CUTPASS]->setDefault(100);
     parameters_[CbcParam::GAPRATIO]->setDefault(0.0);
     parameters_[CbcParam::TIMELIMIT]->setDefault( 1.0e11);
     parameters_[CbcParam::STRONGBRANCHING]->setDefault(0);
     parameters_[CbcParam::NUMBERBEFORE]->setDefault(0);
     break;
   default:
     std::cout << "Unknown strategy!" << std::endl;
     break;
  }

  // Set all parameter values to their defaults to begin with
   
  for (int code = CbcParam::FIRSTPARAM + 1;
       code < CbcParam::LASTPARAM; code++) {
     if (getParam(code)->type() != CoinParam::paramInvalid &&
         getParam(code)->type() != CoinParam::paramAct){ 
        getParam(code)->restoreDefault();
     }
  }
  
}
   
//###########################################################################
//###########################################################################

#ifdef CBC_CLUMSY_CODING
// Synchronize Cbc (and Clp) model - Int and Dbl 
void CbcParameters::synchronizeModel() {
  if (goodModel_) {
    assert (model_);
    // Integer parameters
    int intValue;
    int modelIntValue;
    parameters_[CbcParam::MAXNODES]->getVal(intValue);
    //#define PRINT_CBC_CHANGES
#ifdef PRINT_CBC_CHANGES
    modelIntValue = model_->getIntParam(CbcModel::CbcMaxNumNode);
    if (intValue!=modelIntValue)
      printf("changing MAXNODES from %d to %d at line %d\n",modelIntValue,intValue,__LINE__+1);
#endif
    model_->setIntParam(CbcModel::CbcMaxNumNode, intValue);
    parameters_[CbcParam::MAXNODESNOTIMPROVING]->getVal(intValue);
#ifdef PRINT_CBC_CHANGES
    modelIntValue = model_->getIntParam(CbcModel::CbcMaxNodesNotImproving);
    if (intValue!=modelIntValue)
      printf("changing MAXNODESNOTIMPROVING from %d to %d at line %d\n",modelIntValue,intValue,__LINE__+1);
#endif
    model_->setIntParam(CbcModel::CbcMaxNodesNotImproving, intValue);
    parameters_[CbcParam::MAXSOLS]->getVal(intValue);
#ifdef PRINT_CBC_CHANGES
    modelIntValue = model_->getIntParam(CbcModel::CbcMaxNumSol);
    if (intValue!=modelIntValue)
      printf("changing MAXSOLS from %d to %d at line %d\n",modelIntValue,intValue,__LINE__+1);
#endif
    model_->setIntParam(CbcModel::CbcMaxNumSol, intValue);
    // ?case CBC_PARAM_INT_MAXSAVEDSOLS:
    parameters_[CbcParam::STRONGBRANCHING]->getVal(intValue);
#ifdef PRINT_CBC_CHANGES
    modelIntValue = model_->numberStrong();
    if (intValue!=modelIntValue)
      printf("changing STRONGBRANCHING from %d to %d at line %d\n",modelIntValue,intValue,__LINE__+1);
#endif
    model_->setNumberStrong(intValue);
    parameters_[CbcParam::NUMBERBEFORE]->getVal(intValue);
#ifdef PRINT_CBC_CHANGES
    modelIntValue = model_->numberBeforeTrust();
    if (intValue!=modelIntValue)
      printf("changing NUMBERBEFORE from %d to %d at line %d\n",modelIntValue,intValue,__LINE__+1);
#endif
    model_->setNumberBeforeTrust(intValue);
    parameters_[CbcParam::NUMBERANALYZE]->getVal(intValue);
    modelIntValue = model_->numberAnalyzeIterations();
#ifdef PRINT_CBC_CHANGES
    if (intValue!=modelIntValue)
      printf("changing NUMBERANALYZE from %d to %d at line %d\n",modelIntValue,intValue,__LINE__+1);
#endif
    model_->setNumberAnalyzeIterations(intValue);
    parameters_[CbcParam::CUTPASSINTREE]->getVal(intValue);
#ifdef PRINT_CBC_CHANGES
    modelIntValue = model_->getMaximumCutPasses();
    if (intValue!=modelIntValue)
      printf("changing CUTPASSINTREE from %d to %d at line %d\n",modelIntValue,intValue,__LINE__+1);
#endif
    model_->setMaximumCutPasses(intValue);
    parameters_[CbcParam::CUTPASS]->getVal(intValue);
#ifdef PRINT_CBC_CHANGES
    modelIntValue = model_->getMaximumCutPassesAtRoot();
    if (intValue!=modelIntValue)
      printf("changing CUTPASS from %d to %d at line %d\n",modelIntValue,intValue,__LINE__+1);
#endif
    model_->setMaximumCutPassesAtRoot(intValue);
#ifdef CBC_THREAD
    parameters_[CbcParam::THREADS]->getVal(intValue);
#ifdef PRINT_CBC_CHANGES
    modelIntValue = model_->getNumberThreads();
    if (intValue!=modelIntValue)
      printf("changing THREADS from %d to %d at line %d\n",modelIntValue,intValue,__LINE__+1);
#endif
    model_->setNumberThreads(intValue);
#endif
    parameters_[CbcParam::RANDOMSEED]->getVal(intValue);
#ifdef PRINT_CBC_CHANGES
    modelIntValue = model_->getRandomSeed();
    if (intValue!=modelIntValue)
      printf("changing RANDOMSEED from %d to %d at line %d\n",modelIntValue,intValue,__LINE__+1);
#endif
    model_->setRandomSeed(intValue);
    // Double parameters
    double doubleValue;
    double modelDoubleValue;
    parameters_[CbcParam::INFEASIBILITYWEIGHT]->getVal(doubleValue);
#ifdef PRINT_CBC_CHANGES
    modelDoubleValue = model_->getDblParam(CbcModel::CbcInfeasibilityWeight);
    if (doubleValue!=modelDoubleValue)
      printf("changing INFEASIBILITYWEIGHT from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
    model_->setDblParam(CbcModel::CbcInfeasibilityWeight, intValue);
    parameters_[CbcParam::INTEGERTOLERANCE]->getVal(doubleValue);
#ifdef PRINT_CBC_CHANGES
    modelDoubleValue = model_->getDblParam(CbcModel::CbcIntegerTolerance);
    if (doubleValue!=modelDoubleValue)
      printf("changing INTEGERTOLERANCE from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
    model_->setDblParam(CbcModel::CbcIntegerTolerance, doubleValue);
    parameters_[CbcParam::INCREMENT]->getVal(doubleValue);
#ifdef PRINT_CBC_CHANGES
    modelDoubleValue = model_->getDblParam(CbcModel::CbcCutoffIncrement);
    if (doubleValue!=modelDoubleValue)
      printf("changing INCREMENT from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
    model_->setDblParam(CbcModel::CbcCutoffIncrement, doubleValue);
    parameters_[CbcParam::ALLOWABLEGAP]->getVal(doubleValue);
#ifdef PRINT_CBC_CHANGES
    modelDoubleValue = model_->getDblParam(CbcModel::CbcAllowableGap);
    if (doubleValue!=modelDoubleValue)
      printf("changing ALLOWABLEGAP from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
    model_->setDblParam(CbcModel::CbcAllowableGap, doubleValue);
    parameters_[CbcParam::GAPRATIO]->getVal(doubleValue);
#ifdef PRINT_CBC_CHANGES
    if (doubleValue!=modelDoubleValue)
    modelDoubleValue = model_->getDblParam(CbcModel::CbcAllowableFractionGap);
      printf("changing GAPRATIO from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
    model_->setDblParam(CbcModel::CbcAllowableFractionGap, doubleValue);
    parameters_[CbcParam::CUTOFF]->getVal(doubleValue);
#ifdef PRINT_CBC_CHANGES
    modelDoubleValue = model_->getCutoff();
    if (doubleValue!=modelDoubleValue)
      printf("changing CUTOFF from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
    model_->setCutoff(doubleValue);
    parameters_[CbcParam::TIMELIMIT]->getVal(doubleValue);
#ifdef PRINT_CBC_CHANGES
    modelDoubleValue = model_->getDblParam(CbcModel::CbcMaximumSeconds);
    if (doubleValue!=modelDoubleValue)
      printf("changing TIMELIMIT from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
    model_->setDblParam(CbcModel::CbcMaximumSeconds, doubleValue);
    parameters_[CbcParam::MAXSECONDSNOTIMPROVING]->getVal(doubleValue);
#ifdef PRINT_CBC_CHANGES
    modelDoubleValue = model_->getDblParam(CbcModel::CbcMaxSecondsNotImproving);
    if (doubleValue!=modelDoubleValue)
      printf("changing MAXSECONDSNOTIMPROVING from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
    model_->setDblParam(CbcModel::CbcMaxSecondsNotImproving, doubleValue);
    if (clpParameters_.getModel())
      clpParameters_.synchronizeModel();
  }
}
#endif

//###########################################################################
//###########################################################################

void CbcParameters::addCbcSolverHelpParams() {
  for (int code = CbcParam::FIRSTHELPPARAM + 1;
       code < CbcParam::LASTHELPPARAM; code++) {
    getParam(code)->setPushFunc(CbcParamUtils::doHelpParam);
    getParam(code)->setType(CoinParam::paramAct);
  }
  parameters_[CbcParam::GENERALQUERY]->setup(
      "?", "Print a list of commands", CoinParam::displayPriorityNone);

  parameters_[CbcParam::FULLGENERALQUERY]->setup(
      "???", "Print a list with *all* commands, even those hidden with `?'",
      CoinParam::displayPriorityNone);

  // Need display parameter to resolve ambiguity
  parameters_[CbcParam::HELP]->setup(
      "help", "Print out version, non-standard options and some help",
      "This prints out some help to get a user started. If you're seeing this "
      "message, you should be past that stage.",
      CoinParam::displayPriorityHigh);
}

//###########################################################################
//###########################################################################

void CbcParameters::addCbcSolverActionParams() {

  for (int code = CbcParam::FIRSTACTIONPARAM + 1;
       code < CbcParam::LASTACTIONPARAM; code++) {
    getParam(code)->setType(CoinParam::paramAct);
  }

  parameters_[CbcParam::BAB]->setup(
      "solv!e", "invoke branch and cut to solve the current problem",
      "This does branch and cut. There are many parameters which can affect "
      "the performance.  First just try with default cbcSettings and look "
      "carefully at the log file.  Did cuts help?  Did they take too long?  "
      "Look at output to see which cuts were effective and then do some "
      "tuning.  You will see that the options for cuts are off, on, root and "
      "ifmove.  Off is obvious, on means that this cut generator will be tried "
      "in the branch and cut tree (you can fine tune using 'depth').  Root "
      "means just at the root node while 'ifmove' means that cuts will be used "
      "in the tree if they look as if they are doing some good and moving the "
      "objective value.  If pre-processing reduced the size of the problem or "
      "strengthened many coefficients then it is probably wise to leave it on. "
      " Switch off heuristics which did not provide solutions.  The other "
      "major area to look at is the search.  Hopefully good solutions were "
      "obtained fairly early in the search so the important point is to select "
      "the best variable to branch on.  See whether strong branching did a "
      "good job - or did it just take a lot of iterations.  Adjust the "
      "strongBranching and trustPseudoCosts parameters.",
      CoinParam::displayPriorityHigh);
  parameters_[CbcParam::BAB]->setPushFunc(CbcParamUtils::doBaCParam);

  parameters_[CbcParam::ENVIRONMENT]->setup(
      "environ!ment", "Read commands from environment",
      "This starts reading from environment variable COIN_ENVIRONMENT.",
      CoinParam::displayPriorityNone);
  parameters_[CbcParam::ENVIRONMENT]->setPushFunc(CbcParamUtils::doNothingParam);

  parameters_[CbcParam::EXIT]->setup(
      "end", "Stops execution",
      "This stops execution; end, exit, quit and stop are synonyms.",
      CoinParam::displayPriorityHigh);
  parameters_[CbcParam::EXIT]->setPushFunc(CbcParamUtils::doExitParam);

  parameters_[CbcParam::EXPORT]->setup(
      "export", "Export model as mps file",
      "This will write an MPS format file to the given file name.  It will use "
      "the default directory given by 'directory'.  A name of '$' will use the "
      "previous value for the name.  This is initialized to 'default.mps'. It "
      "can be useful to get rid of the original names and go over to using "
      "Rnnnnnnn and Cnnnnnnn.  This can be done by setting 'keepnames' off "
      "before importing mps file.",
      CoinParam::displayPriorityHigh);

  parameters_[CbcParam::IMPORT]->setup(
      "import", "Import model from file", 
      "This will read an MPS format file from the given file name.  It will "
      "use the default directory given by 'directory'.  A name of '$' will use "
      "the previous value for the name.  This is initialized to '', i.e., it "
      "must be set.  If you have libgz then it can read compressed files "
      "'xxxxxxxx.gz'.",
      CoinParam::displayPriorityHigh);
  parameters_[CbcParam::IMPORT]->setPushFunc(CbcParamUtils::doImportParam);

  parameters_[CbcParam::MIPLIB]->setup("miplib", "Do some of miplib test set", "",
                            CoinParam::displayPriorityHigh);

  parameters_[CbcParam::OUTDUPROWS]->setup(
      "outDup!licates", "Takes duplicate rows, etc., out of the integer model",
      "", CoinParam::displayPriorityNone);

  parameters_[CbcParam::PRINTVERSION]->setup(
      "version", "Print version", "", CoinParam::displayPriorityHigh);
  parameters_[CbcParam::PRINTVERSION]->setPushFunc(CbcParamUtils::doVersionParam);

  parameters_[CbcParam::PRINTSOL]->setup(
      "writeS!olution", "writes solution to file (or stdout)",
      "This will write a binary solution file to the file set by solutionFile.",
      CoinParam::displayPriorityHigh);

  parameters_[CbcParam::READMIPSTART]->setup(
      "mipS!tart", "reads an initial feasible solution from file",
      "The MIPStart allows one to enter an initial integer feasible solution "
      "to CBC. Values of the main decision variables which are active (have "
      "non-zero values) in this solution are specified in a text  file. The "
      "text file format used is the same of the solutions saved by CBC, but "
      "not all fields are required to be filled. First line may contain the "
      "solution status and will be ignored, remaining lines contain column "
      "indexes, names and values as in this example:\n\n Stopped on iterations "
      "- objective value 57597.00000000\n      0  x(1,1,2,2)               1 "
      "\n      1  x(3,1,3,2)               1 \n      5  v(5,1)                 "
      "  2 \n      33 x(8,1,5,2)               1 \n      ...\n\n Column "
      "indexes are also ignored since pre-processing can change them. There is "
      "no need to include values for continuous or integer auxiliary "
      "variables, since they can be computed based on main decision variables. "
      "Starting CBC with an integer feasible solution can dramatically improve "
      "its performance: several MIP heuristics (e.g. RINS) rely on having at "
      "least one feasible solution available and can start immediately if the "
      "user provides one. Feasibility Pump (FP) is a heuristic which tries to "
      "overcome the problem of taking too long to find feasible solution (or "
      "not finding at all), but it not always succeeds. If you provide one "
      "starting solution you will probably save some time by disabling FP. "
      "\n\n Knowledge specific to your problem can be considered to write an "
      "external module to quickly produce an initial feasible solution - some "
      "alternatives are the implementation of simple greedy heuristics or the "
      "solution (by CBC for example) of a simpler model created just to find a "
      "feasible solution. \n\n Silly options added.  If filename ends .low "
      "then integers not mentioned are set low - also .high, .lowcheap, "
      ".highcheap, .lowexpensive, .highexpensive where .lowexpensive sets "
      "costed ones to make expensive others low. Also if filename starts "
      "empty. then no file is read at all - just actions done. \n\n Question "
      "and suggestions regarding MIPStart can be directed to\n "
      "haroldo.santos@gmail.com. ");

  parameters_[CbcParam::READSOL]->setup(
      "readS!olution", "reads solution from file (synonym for mipStart)",
      "This will read a binary solution file from the file set by solutionFile."
      "This is a synonym for mipStart. See its help for more information",
      CoinParam::displayPriorityHigh);

  parameters_[CbcParam::READMODEL]->setup(
      "readM!odel", "writes problem to file", 
      "This will save the problem to the file name set by modelFile "
      "for future use by readModel.", CoinParam::displayPriorityHigh);

  parameters_[CbcParam::SHOWUNIMP]->setup(
      "unimp!lemented", "Report unimplemented commands.", "",
      CoinParam::displayPriorityNone);
  parameters_[CbcParam::SHOWUNIMP]->setPushFunc(CbcParamUtils::doUnimplementedParam);

  parameters_[CbcParam::SOLVECONTINUOUS]->setup(
      "initialS!olve", "Solve to continuous optimum",
      "This just solves the problem to the continuous optimum, without adding "
      "any cuts.",
      CoinParam::displayPriorityHigh);

  parameters_[CbcParam::STATISTICS]->setup(
      "stat!istics", "Print some statistics",
      "This command prints some statistics for the current model. If log level "
      ">1 then more is printed. These are for presolved model if presolve on "
      "(and unscaled).",
      CoinParam::displayPriorityHigh);

#if 0
  // Need to figure out what to do here. Same parameter can't have two names...
  parameters_[CbcParam::STDIN]->setup( "-", "Switch to interactive command line mode", ""
                            CoinParam::displayPriorityNone);
  parameters_[CbcParam::STDIN]->setPushFunc(CbcParamUtils::doNothingParam);
#endif

  parameters_[CbcParam::STDIN]->setup(
      "stdin", "Switch to interactive command line mode", "",
      CoinParam::displayPriorityNone);

  parameters_[CbcParam::STRENGTHEN]->setup(
      "strengthen", "Create strengthened problem",
      "This creates a new problem by applying the root node cuts. All tight "
      "constraints will be in resulting problem.",
      CoinParam::displayPriorityHigh);

  parameters_[CbcParam::UNITTEST]->setup(
      "unitTest", "Do unit test", "This exercises the unit test.",
      CoinParam::displayPriorityHigh);

  parameters_[CbcParam::WRITEGMPLSOL]->setup(
      "writeGSolu!tion", "Puts glpk solution to file",
      "Will write a glpk solution file to the given file name.  It will use "
      "the default directory given by 'directory'.  A name of '$' will use the "
      "previous value for the name.  This is initialized to 'stdout' (this "
      "defaults to ordinary solution if stdout). If problem created from gmpl "
      "model - will do any reports.",
      CoinParam::displayPriorityHigh);

  parameters_[CbcParam::WRITEMODEL]->setup(
      "writeM!odel", "writes solution to file (or stdout)", 
      "This will write the problem to the file name set by modelFile"
      "for future use by readModel.", CoinParam::displayPriorityHigh);

  parameters_[CbcParam::WRITENEXTSOL]->setup(
      "nextB!estSolution", "Prints next best saved solution to file",
      "To write best solution, just use writeSolution.  This prints next best (if "
      "exists) and then deletes it. This will write a primitive solution file "
      "to the given file name.  It will use the directory set by "
      "nextBestSolutionFile. The amount of output can be varied "
      "using printi!ngOptions or printMask.");

  parameters_[CbcParam::WRITESOL]->setup(
      "writeS!olution", "writes solution to file (or stdout)",
      "This will write a binary solution file to the file set by solutionFile.",
      CoinParam::displayPriorityHigh);

}

//###########################################################################
//###########################################################################

void CbcParameters::addCbcSolverStrParams() {

  for (int code = CbcParam::FIRSTSTRINGPARAM + 1;
       code < CbcParam::LASTSTRINGPARAM; code++) {
    getParam(code)->setType(CoinParam::paramStr);
  }

  parameters_[CbcParam::CSVSTATISTICS]->setup(
      "csv!Statistics", "Create one line of statistics",
      "This appends statistics to given file name.  It will use the default "
      "directory given by 'directory'.  A name of '$' will use the previous "
      "value for the name.  This is initialized to '', i.e. it must be set.  "
      "Adds header if file empty or does not exist.",
      CoinParam::displayPriorityLow);
  parameters_[CbcParam::CSVSTATISTICS]->setPushFunc(CbcParamUtils::pushCbcSolverStrParam);

  parameters_[CbcParam::DEBUG]->setup(
      "debug!In", "Read/write valid solution from/to file", 
      "This will read a solution file from the given file name.  It will use "
      "the default directory given by 'directory'.  A name of '$' will use the "
      "previous value for the name.  This is initialized to '', i.e. it must "
      "be set.\n\nIf set to create it will create a file called debug.file "
      "after B&C search; if set to createAfterPre it will create the file "
      "before undoing preprocessing.\n\nThe idea is that if you suspect a bad "
      "cut generator and you did not use preprocessing you can do a good run "
      "with debug set to 'create' and then switch on the cuts you suspect and "
      "re-run with debug set to 'debug.file'  Similarly if you do use "
      "preprocessing, but use createAfterPre.  The create case has the same "
      "effect as saveSolution.",
      CoinParam::displayPriorityNone);
  parameters_[CbcParam::DEBUG]->setPushFunc(CbcParamUtils::doDebugParam);

  parameters_[CbcParam::DIRECTORY]->setup(
      "directory", "Set Default directory for import etc.",
      "This sets the directory which import, export, saveModel, restoreModel "
      "etc. will use. It is initialized to the current directory.");
  parameters_[CbcParam::DIRECTORY]->setPushFunc(CbcParamUtils::pushCbcSolverStrParam);

  parameters_[CbcParam::DIRSAMPLE]->setup(
      "dirSample", "Set directory where the COIN-OR sample problems are.",
      "This sets the directory where the COIN-OR sample problems reside. It is "
      "used only when -unitTest is passed to clp. clp will pick up the test "
      "problems from this directory. It is initialized to '../../Data/Sample'",
      CoinParam::displayPriorityLow);
  parameters_[CbcParam::DIRSAMPLE]->setPushFunc(CbcParamUtils::pushCbcSolverStrParam);

  parameters_[CbcParam::DIRNETLIB]->setup(
      "dirNetlib", "Set directory where the netlib problems are.",
      "This sets the directory where the netlib problems reside. One can get "
      "the netlib problems from COIN-OR or from the main netlib site. This "
      "parameter is used only when -netlib is passed to cbc. cbc will pick up "
      "the netlib problems from this directory. If clp is built without zlib "
      "support then the problems must be uncompressed. It is initialized to "
      "'../../Data/Netlib'",
      CoinParam::displayPriorityLow);
  parameters_[CbcParam::DIRNETLIB]->setPushFunc(CbcParamUtils::pushCbcSolverStrParam);

  parameters_[CbcParam::DIRMIPLIB]->setup(
      "dirMiplib", "Set directory where the miplib 2003 problems are.",
      "This sets the directory where the miplib 2003 problems reside. One can "
      "get the miplib problems from COIN-OR or from the main miplib site. This "
      "parameter is used only when -miplib is passed to cbc. cbc will pick up "
      "the miplib problems from this directory. If cbc is built without zlib "
      "support then the problems must be uncompressed. It is initialized to "
      "'../../Data/miplib3'",
      CoinParam::displayPriorityLow);
  parameters_[CbcParam::DIRMIPLIB]->setPushFunc(CbcParamUtils::pushCbcSolverStrParam);

  parameters_[CbcParam::GMPLSOLFILE]->setup(
      "gmplSolutionF!ile", "sets name for file to store GMPL solution in",
      "This will set the name the GMPL solution will be written to and read from. "
      "This is initialized to 'gmpl.sol'. ",
      CoinParam::displayPriorityHigh);

  parameters_[CbcParam::MODELFILE]->setup(
      "modelF!ile", "sets name for file to store model in",
      "This will set the name the model will be written to and read from. "
      "This is initialized to 'prob.mod'. ",
      CoinParam::displayPriorityHigh);

  parameters_[CbcParam::NEXTSOLFILE]->setup(
      "nextSolutionF!ile", "sets name for file to store suboptimal solutions in",
      "This will set the name solutions will be written to and read from "
      "This is initialized to 'next.sol'. ",
      CoinParam::displayPriorityHigh);

  parameters_[CbcParam::PRINTMASK]->setup(
      "printM!ask", "Control printing of solution with a regular expression",
      "If set then only those names which match mask are printed in a "
      "solution. '?' matches any character and '*' matches any set of "
      "characters.  The default is '' (unset) so all variables are printed. "
      "This is only active if model has names.");
  parameters_[CbcParam::PRINTMASK]->setPushFunc(CbcParamUtils::doPrintMaskParam);

  parameters_[CbcParam::SOLUTIONFILE]->setup(
      "solutionF!ile", "sets name for file to store solution in",
      "This will set the name the solution will be saved to and read from. "
      "By default, solutions are written to 'opt.sol'. To print to stdout, "
      "use printSolution.", CoinParam::displayPriorityHigh);

  parameters_[CbcParam::PRIORITYIN]->setup(
      "prio!rityIn", "Import priorities etc from file",
      "This will read a file with priorities from the given file name.  It "
      "will use the default directory given by 'directory'.  A name of '$' "
      "will use the previous value for the name.  This is initialized to '', "
      "i.e. it must be set.  This can not read from compressed files. File is "
      "in csv format with allowed headings - name, number, priority, "
      "direction, up, down, solution.  Exactly one of name and number must be "
      "given.");
  parameters_[CbcParam::PRIORITYIN]->setPushFunc(CbcParamUtils::pushCbcSolverStrParam);
}

//###########################################################################
//###########################################################################

void CbcParameters::addCbcSolverKwdParams() {
  for (int code = CbcParam::FIRSTKWDPARAM + 1;
       code < CbcParam::LASTKWDPARAM; code++) {
    getParam(code)->setPushFunc(CbcParamUtils::pushCbcSolverKwdParam);
    getParam(code)->setType(CoinParam::paramKwd);
  }

  parameters_[CbcParam::COMMANDPRINTLEVEL]->setup(
      "allC!ommands", "What priority level of commands to print", 
      "For the sake of your sanity, only the more useful and simple commands "
      "are printed out on ?.",
      CoinParam::displayPriorityNone);
  parameters_[CbcParam::COMMANDPRINTLEVEL]->appendKwd("high", CbcParameters::displayHigh);
  parameters_[CbcParam::COMMANDPRINTLEVEL]->appendKwd("all", CbcParameters::displayAll);
  parameters_[CbcParam::COMMANDPRINTLEVEL]->appendKwd("highlow", CbcParameters::displayLowHigh);

  parameters_[CbcParam::CLQSTRENGTHENING]->setup(
      "clqstr!engthen",
      "Whether and when to perform Clique Strengthening preprocessing routine", "");
  parameters_[CbcParam::CLQSTRENGTHENING]->appendKwd("after", CbcParameters::ClqStrAfter);
  parameters_[CbcParam::CLQSTRENGTHENING]->appendKwd("off", CbcParameters::ClqStrOff);
  parameters_[CbcParam::CLQSTRENGTHENING]->appendKwd("before", CbcParameters::ClqStrBefore);

  parameters_[CbcParam::BRANCHPRIORITY]->setup(
      "cost!Strategy", "Whether to use costs or column order as priorities",
      "This orders the variables in order of their absolute costs - with "
      "largest cost ones being branched on first.  This primitive strategy can "
      "be surprisingly effective.  The column order option is obviously not on "
      "costs but it's easy to implement.");
  parameters_[CbcParam::BRANCHPRIORITY]->appendKwd("off", CbcParameters::BPOff);
  parameters_[CbcParam::BRANCHPRIORITY]->appendKwd("pri!orities", CbcParameters::BPCost);
  parameters_[CbcParam::BRANCHPRIORITY]->appendKwd("column!Order", CbcParameters::BPOrder);

  parameters_[CbcParam::CUTOFFCONSTRAINT]->setup(
      "constraint!fromCutoff", "Whether to use cutoff as constraint", 
      "For some problems, cut generators and general branching work better if "
      "the problem would be infeasible if the cost is too high. "
      "If this option is enabled, the objective function is added as a "
      "constraint which right hand side is set to the current cutoff value "
      "(objective value of best known solution)");
  parameters_[CbcParam::CUTOFFCONSTRAINT]->appendKwd("off", CbcParameters::COOff);
  parameters_[CbcParam::CUTOFFCONSTRAINT]->appendKwd("on", CbcParameters::COOn);
  parameters_[CbcParam::CUTOFFCONSTRAINT]->appendKwd("variable", CbcParameters::COVariable);
  parameters_[CbcParam::CUTOFFCONSTRAINT]->appendKwd("forcevariable", CbcParameters::COForceVariable);
  parameters_[CbcParam::CUTOFFCONSTRAINT]->appendKwd("conflict", CbcParameters::COConflict);

  parameters_[CbcParam::INTPRINT]->setup(
      "printi!ngOptions", "Print options", 
      "This changes the amount and format of printing a solution:\n normal - "
      "nonzero column variables \ninteger - nonzero integer column variables\n "
      "special - in format suitable for OsiRowCutDebugger\n rows - nonzero "
      "column variables and row activities\n all - all column variables and "
      "row activities.\n\n For non-integer problems 'integer' and 'special' "
      "act like 'normal'.  Also see printMask for controlling output.");
  parameters_[CbcParam::INTPRINT]->appendKwd("normal", CbcParameters::PMNormal);
  parameters_[CbcParam::INTPRINT]->appendKwd("integer", CbcParameters::PMInteger);
  parameters_[CbcParam::INTPRINT]->appendKwd("special", CbcParameters::PMSpecial);
  parameters_[CbcParam::INTPRINT]->appendKwd("rows", CbcParameters::PMRows);
  parameters_[CbcParam::INTPRINT]->appendKwd("all", CbcParameters::PMAll);

  parameters_[CbcParam::NODESTRATEGY]->setup(
      "node!Strategy",
      "What strategy to use to select the next node from the branch and cut "
      "tree",
      "Normally before a feasible solution is found, CBC will choose a node "
      "with fewest infeasibilities. Alternatively, one may choose tree-depth "
      "as the criterion. This requires the minimal amount of memory, but may "
      "take a long time to find the best solution. Additionally, one may "
      "specify whether up or down branches must be selected first (the up-down "
      "choice will carry on after a first solution has been bound). The choice "
      "'hybrid' does breadth first on small depth nodes and then switches to "
      "'fewest'.");
  parameters_[CbcParam::NODESTRATEGY]->appendKwd("hybrid", CbcParameters::NSHybrid);
  parameters_[CbcParam::NODESTRATEGY]->appendKwd("fewest", CbcParameters::NSFewest);
  parameters_[CbcParam::NODESTRATEGY]->appendKwd("depth", CbcParameters::NSDepth);
  parameters_[CbcParam::NODESTRATEGY]->appendKwd("upfewest", CbcParameters::NSUpFewest);
  parameters_[CbcParam::NODESTRATEGY]->appendKwd("downfewest", CbcParameters::NSDownFewest);
  parameters_[CbcParam::NODESTRATEGY]->appendKwd("updepth", CbcParameters::NSUpDepth);
  parameters_[CbcParam::NODESTRATEGY]->appendKwd("downdepth", CbcParameters::NSDownDepth);

  parameters_[CbcParam::ORBITAL]->setup(
      "Orbit!alBranching", "Whether to try orbital branching", 
      "This switches on Orbital branching. Value 'on' just adds orbital, "
      "'strong' tries extra fixing in strong branching.");
  parameters_[CbcParam::ORBITAL]->appendKwd("off", CbcParameters::OBOff);
  parameters_[CbcParam::ORBITAL]->appendKwd("slowish", CbcParameters::OBSlowish);
  parameters_[CbcParam::ORBITAL]->appendKwd("strong", CbcParameters::OBStrong);
  parameters_[CbcParam::ORBITAL]->appendKwd("force", CbcParameters::OBForce);
  parameters_[CbcParam::ORBITAL]->appendKwd("simple", CbcParameters::OBSimple);
  parameters_[CbcParam::ORBITAL]->appendKwd("on", CbcParameters::OBOn);
  parameters_[CbcParam::ORBITAL]->appendKwd("moreprinting", CbcParameters::OBMorePrinting);  

  parameters_[CbcParam::PREPROCESS]->setup(
      "preprocess", "Whether to use integer preprocessing", 
      "This tries to reduce size of the model in a similar way to presolve and "
      "it also tries to strengthen the model. This can be very useful and is "
      "worth trying.  save option saves on file presolved.mps.  equal will "
      "turn <= cliques into ==.  sos will create sos sets if all 0-1 in sets "
      "(well one extra is allowed) and no overlaps.  trysos is same but allows "
      "any number extra. equalall will turn all valid inequalities into "
      "equalities with integer slacks. strategy is as on but uses "
      "CbcStrategy.");
  parameters_[CbcParam::PREPROCESS]->appendKwd("off", CbcParameters::IPPOff);
  parameters_[CbcParam::PREPROCESS]->appendKwd("on", CbcParameters::IPPOn);
  parameters_[CbcParam::PREPROCESS]->appendKwd("save", CbcParameters::IPPSave);
  parameters_[CbcParam::PREPROCESS]->appendKwd("equal", CbcParameters::IPPEqual);
  parameters_[CbcParam::PREPROCESS]->appendKwd("sos", CbcParameters::IPPSOS);
  parameters_[CbcParam::PREPROCESS]->appendKwd("trysos", CbcParameters::IPPTrySOS);
  parameters_[CbcParam::PREPROCESS]->appendKwd("equalall", CbcParameters::IPPEqualAll);
  parameters_[CbcParam::PREPROCESS]->appendKwd("strategy", CbcParameters::IPPStrategy);
  parameters_[CbcParam::PREPROCESS]->appendKwd("aggregate", CbcParameters::IPPAggregate);
  parameters_[CbcParam::PREPROCESS]->appendKwd("forcesos", CbcParameters::IPPForceSOS);
  parameters_[CbcParam::PREPROCESS]->appendKwd("stop!aftersaving", CbcParameters::IPPStopAfterSaving);
  parameters_[CbcParam::PREPROCESS]->appendKwd("equalallstop", CbcParameters::IPPEqualAllStop);

  parameters_[CbcParam::SOSPRIORITIZE]->setup(
      "sosP!rioritize", "How to deal with SOS priorities", 
      "This sets priorities for SOS.  Values 'high' and 'low' just set a "
      "priority relative to the for integer variables.  Value 'orderhigh' "
      "gives first highest priority to the first SOS and integer variables a "
      "low priority.  Value 'orderlow' gives integer variables a high priority "
      "then SOS in order.");
  parameters_[CbcParam::SOSPRIORITIZE]->appendKwd("off", CbcParameters::SOSOff);
  parameters_[CbcParam::SOSPRIORITIZE]->appendKwd("high", CbcParameters::SOSHigh);
  parameters_[CbcParam::SOSPRIORITIZE]->appendKwd("low", CbcParameters::SOSLow);
  parameters_[CbcParam::SOSPRIORITIZE]->appendKwd("orderhigh", CbcParameters::SOSOrderHigh);
  parameters_[CbcParam::SOSPRIORITIZE]->appendKwd("orderlow", CbcParameters::SOSOrderLow);

  parameters_[CbcParam::STRATEGY]->setup(
      "strat!egy", "Switches on groups of features", 
      "This turns on newer features. "
      "Default uses Gomory cuts with a tolerance of 0.01 at the root "
      "node, does a possible restart after 100 nodes if many variables could "
      "be fixed, activates a diving and RINS heuristic, and makes the "
      "feasibility pump more aggressive."); // This does not apply to unit tests
                                            // (where 'experiment' may have
                                            // similar effects)
  parameters_[CbcParam::STRATEGY]->appendKwd("default", CbcParameters::StrategyDefault);
  parameters_[CbcParam::STRATEGY]->appendKwd("easy", CbcParameters::StrategyEasy);
  parameters_[CbcParam::STRATEGY]->appendKwd("aggressive", CbcParameters::StrategyAggressive);

  parameters_[CbcParam::TIMEMODE]->setup(
      "timeM!ode", "Whether to use CPU or elapsed time",
      "cpu uses CPU time for stopping, while elapsed uses elapsed time. (On "
      "Windows, elapsed time is always used).");
  parameters_[CbcParam::TIMEMODE]->appendKwd( "cpu", CbcParameters::ClockCpu);
  parameters_[CbcParam::TIMEMODE]->appendKwd("elapsed", CbcParameters::ClockElapsed);

  parameters_[CbcParam::USECGRAPH]->setup(
      "cgraph",
      "Whether to use the conflict graph-based preprocessing and cut "
      "separation routines.",
      "This switches the conflict graph-based preprocessing and cut separation "
      "routines (CglBKClique, CglOddWheel and CliqueStrengthening) on or off. "
      "Values: \n\t off: turns these routines off;\n\t on: turns these "
      "routines on; \n\t clq: turns these routines off and enables the cut "
      "separator of CglClique.");
  parameters_[CbcParam::USECGRAPH]->appendKwd("on", CbcParameters::CGraphOn);
  parameters_[CbcParam::USECGRAPH]->appendKwd("off", CbcParameters::CGraphOff);
  parameters_[CbcParam::USECGRAPH]->appendKwd("clq", CbcParameters::CGraphClique);
}

//###########################################################################
//###########################################################################

void CbcParameters::addCbcSolverDblParams() {
  for (int code = CbcParam::FIRSTDBLPARAM + 1;
       code < CbcParam::LASTDBLPARAM; code++) {
    getParam(code)->setPushFunc(CbcParamUtils::pushCbcSolverDblParam);
  }

  parameters_[CbcParam::ARTIFICIALCOST]->setup(
      "artif!icialCost",
      "Costs >= this treated as artificials in feasibility pump", 0.0,
      COIN_DBL_MAX, "", CoinParam::displayPriorityLow);

  parameters_[CbcParam::DEXTRA3]->setup(
      "dextra3", "Extra double parameter 3", -COIN_DBL_MAX, COIN_DBL_MAX,
      "", CoinParam::displayPriorityNone);

  parameters_[CbcParam::DEXTRA4]->setup(
      "dextra4", "Extra double parameter 4", -COIN_DBL_MAX, COIN_DBL_MAX,
      "", CoinParam::displayPriorityNone);

  parameters_[CbcParam::DEXTRA5]->setup(
      "dextra5", "Extra double parameter 5", -COIN_DBL_MAX, COIN_DBL_MAX,
      "", CoinParam::displayPriorityNone);

  parameters_[CbcParam::DJFIX]->setup(
      "fix!OnDj", "Try heuristic that fixes variables based on reduced costs",
      -1.0e20, 1.0e20,
      "If set, integer variables with reduced costs greater than the specified "
      "value will be fixed before branch and bound - use with extreme "
      "caution!");

  parameters_[CbcParam::FAKECUTOFF]->setup(
      "pumpC!utoff", "Fake cutoff for use in feasibility pump", -COIN_DBL_MAX,
      COIN_DBL_MAX,
      "A value of 0.0 means off. Otherwise, add a constraint forcing objective "
      "below this value in feasibility pump",
      CoinParam::displayPriorityLow);

  parameters_[CbcParam::FAKEINCREMENT]->setup(
      "pumpI!ncrement", "Fake increment for use in feasibility pump",
      -COIN_DBL_MAX, COIN_DBL_MAX, 
      "A value of 0.0 means off. Otherwise, add a constraint forcing objective "
      "below this value in feasibility pump",
      CoinParam::displayPriorityLow);

  parameters_[CbcParam::SMALLBAB]->setup(
      "fraction!forBAB", "Fraction in feasibility pump", 1.0e-5, 1.1,
      "After a pass in the feasibility pump, variables which have not moved "
      "about are fixed and if the preprocessed model is smaller than this "
      "fraction of the original problem, a few nodes of branch and bound are "
      "done on the reduced problem.",
      CoinParam::displayPriorityLow);

  parameters_[CbcParam::TIGHTENFACTOR]->setup(
      "tighten!Factor",
      "Tighten bounds using value times largest activity at continuous "
      "solution",
      0.0, COIN_DBL_MAX, "This sleazy trick can help on some problems.");
}

//###########################################################################
//###########################################################################

void CbcParameters::addCbcSolverIntParams() {

  for (int code = CbcParam::FIRSTINTPARAM + 1;
       code < CbcParam::LASTINTPARAM; code++) {
    getParam(code)->setPushFunc(CbcParamUtils::pushCbcSolverIntParam);
  }

  parameters_[CbcParam::BKPIVOTINGSTRATEGY]->setup(
      "bkpivot!ing", "Pivoting strategy used in Bron-Kerbosch algorithm", 0, 6);

  parameters_[CbcParam::BKMAXCALLS]->setup(
      "bkmaxcalls",
      "Maximum number of recursive calls made by Bron-Kerbosch algorithm", 1,
      COIN_INT_MAX);

  parameters_[CbcParam::BKCLQEXTMETHOD]->setup(
      "bkclqext!method",
      "Strategy used to extend violated cliques found by BK Clique Cut "
      "Separation routine",
      0, 5,
      "Sets the method used in the extension module of BK Clique Cut "
      "Separation routine: 0=no extension; 1=random; 2=degree; 3=modified "
      "degree; 4=reduced cost(inversely proportional); 5=reduced "
      "cost(inversely proportional) + modified degree");

  parameters_[CbcParam::CPP]->setup(
      "cpp!Generate", "Generates C++ code", 0, 4,
      "Once you like what the stand-alone solver does then this allows you to "
      "generate user_driver.cpp which approximates the code.  0 gives simplest "
      "driver, 1 generates saves and restores, 2 generates saves and restores "
      "even for variables at default value. 4 bit in cbc generates size "
      "dependent code rather than computed values.");

  parameters_[CbcParam::CUTDEPTH]->setup(
      "cutD!epth", "Depth in tree at which to do cuts", -1, 999999,
      "Cut generators may be off, on only at the root, on if they look useful, "
      "and on at some interval.  If they are done every node then that is "
      "that, but it may be worth doing them every so often.  The original "
      "method was every so many nodes but it is more logical to do it whenever "
      "depth in tree is a multiple of K.  This option does that and defaults "
      "to -1 (off).");

  parameters_[CbcParam::CUTLENGTH]->setup(
      "cutL!ength", "Length of a cut", -1, COIN_INT_MAX,
      "At present this only applies to Gomory cuts. -1 (default) leaves as is. "
      "Any value >0 says that all cuts <= this length can be generated both at "
      "root node and in tree. 0 says to use some dynamic lengths.  If value "
      ">=10,000,000 then the length in tree is value%10000000 - so 10000100 "
      "means unlimited length at root and 100 in tree.");

  parameters_[CbcParam::CUTPASSINTREE]->setup(
      "passT!reeCuts",
      "Number of rounds that cut generators are applied in the tree",
      -COIN_INT_MAX, COIN_INT_MAX);

  parameters_[CbcParam::DEPTHMINIBAB]->setup(
      "depth!MiniBab", "Depth at which to try mini branch-and-bound",
      -COIN_INT_MAX, COIN_INT_MAX,
      "Rather a complicated parameter but can be useful. -1 means off for "
      "large problems but on as if -12 for problems where rows+columns<500, -2 "
      "means use Cplex if it is linked in.  Otherwise if negative then go into "
      "depth first complete search fast branch and bound when depth>= -value-2 "
      "(so -3 will use this at depth>=1).  This mode is only switched on after "
      "500 nodes.  If you really want to switch it off for small problems then "
      "set this to -999.  If >=0 the value doesn't matter very much.  The code "
      "will do approximately 100 nodes of fast branch and bound every now and "
      "then at depth>=5. The actual logic is too twisted to describe here.");

  parameters_[CbcParam::DIVEOPT]->setup(
      "diveO!pt", "Diving options", -1, 20,
      "If >2 && <20 then modify diving options -	 \n\t3 only at root "
      "and if no solution,	 \n\t4 only at root and if this heuristic has "
      "not got solution,	 \n\t5 decay only if no solution,	 \n\t6 "
      "if depth <3 or decay,	 \n\t7 run up to 2 times if solution found 4 "
      "otherwise,	 \n\t>10 All only at root (DivingC normal as "
      "value-10),	 \n\t>20 All with value-20).",
      CoinParam::displayPriorityLow);

  parameters_[CbcParam::DIVEOPTSOLVES]->setup(
                                              "diveS!olves", "Diving solve option", -1, 200000,
      "If >0 then do up to this many solves. However, the last digit is "
      "ignored and used for extra options: 1-3 enables fixing of satisfied "
      "integer variables (but not at bound), where 1 switches this off for "
      "that dive if the dive goes infeasible, and 2 switches it off "
      "permanently if the dive goes infeasible.",
      CoinParam::displayPriorityLow);

  parameters_[CbcParam::DUMMY]->setup(
      "sleep", "for debug", 0, 9999,
      "If passed to solver from ampl, then ampl will wait so that you can copy "
      ".nl file for debug.",
      CoinParam::displayPriorityNone);

  parameters_[CbcParam::EXPERIMENT]->setup(
      "exper!iment", "Whether to use testing features", -1, 200000,
      "Defines how adventurous you want to be in using new ideas. 0 then no "
      "new ideas, 1 fairly sensible, 2 a bit dubious, 3 you are on your own!",
      CoinParam::displayPriorityLow);

  parameters_[CbcParam::EXTRA1]->setup(
      "extra1", "Extra integer parameter 1", -COIN_INT_MAX, COIN_INT_MAX,
      "", CoinParam::displayPriorityLow);

  parameters_[CbcParam::EXTRA2]->setup(
      "extra2", "Extra integer parameter 2", -COIN_INT_MAX, COIN_INT_MAX,
      "", CoinParam::displayPriorityLow);

  parameters_[CbcParam::EXTRA3]->setup(
      "extra3", "Extra integer parameter 3", -COIN_INT_MAX, COIN_INT_MAX, 
      "", CoinParam::displayPriorityLow);

  parameters_[CbcParam::EXTRA4]->setup(
      "extra4", "Extra integer parameter 4", -COIN_INT_MAX, COIN_INT_MAX, 
      "", CoinParam::displayPriorityLow);

  parameters_[CbcParam::EXTRAVARIABLES]->setup(
      "extraV!ariables", "Allow creation of extra integer variables",
      -COIN_INT_MAX, COIN_INT_MAX,
      "Switches on a trivial re-formulation that introduces extra integer "
      "variables to group together variables with same cost.",
      CoinParam::displayPriorityLow);

  parameters_[CbcParam::FPUMPITS]->setup(
      "passF!easibilityPump", "How many passes in feasibility pump", 0, 10000,
      "This fine tunes the Feasibility Pump heuristic by doing more or fewer "
      "passes.");

  parameters_[CbcParam::FPUMPTUNE]->setup(
      "pumpT!une", "Dubious ideas for feasibility pump", 0, 100000000,
      "This fine tunes Feasibility Pump     \n\t>=10000000 use as objective "
      "weight switch     \n\t>=1000000 use as accumulate switch     \n\t>=1000 "
      "use index+1 as number of large loops     \n\t==100 use objvalue "
      "+0.05*fabs(objvalue) as cutoff OR fakeCutoff if set     \n\t%100 == "
      "10,20 affects how each solve is done     \n\t1 == fix ints at bounds, 2 "
      "fix all integral ints, 3 and continuous at bounds. If accumulate is on "
      "then after a major pass, variables which have not moved are fixed and a "
      "small branch and bound is tried.");

  parameters_[CbcParam::FPUMPTUNE2]->setup(
      "moreT!une", "Yet more dubious ideas for feasibility pump", 0, 100000000,
      "Yet more ideas for Feasibility Pump     \n\t/100000 == 1 use box "
      "constraints and original obj in cleanup     \n\t/1000 == 1 Pump will "
      "run twice if no solution found     \n\t/1000 == 2 Pump will only run "
      "after root cuts if no solution found     \n\t/1000 >10 as above but "
      "even if solution found     \n\t/100 == 1,3.. exact 1.0 for objective "
      "values     \n\t/100 == 2,3.. allow more iterations per pass     \n\t n "
      "fix if value of variable same for last n iterations.",
      CoinParam::displayPriorityNone);

  parameters_[CbcParam::HEUROPTIONS]->setup(
      "hOp!tions", "Heuristic options", -COIN_INT_MAX, COIN_INT_MAX,
      "Value 1 stops heuristics immediately if the allowable gap has been "
      "reached. Other values are for the feasibility pump - 2 says do exact "
      "number of passes given, 4 only applies if an initial cutoff has been "
      "given and says relax after 50 passes, while 8 will adapt the cutoff rhs "
      "after the first solution if it looks as if the code is stalling.",
      CoinParam::displayPriorityLow);

  parameters_[CbcParam::LOGLEVEL]->setup(
      "log!Level", "Level of detail in CBC output.", -1, 999999,
      "If set to 0 then there should be no output in normal circumstances. A "
      "value of 1 is probably the best value for most uses, while 2 and 3 give "
      "more information.");

  parameters_[CbcParam::LPLOGLEVEL]->setup(
      "log!Level", "Level of detail in LP solver output.", -1, 999999,
      "If set to 0 then there should be no output in normal circumstances. A "
      "value of 1 is probably the best value for most uses, while 2 and 3 give "
      "more information.");

  parameters_[CbcParam::MAXHOTITS]->setup(
      "hot!StartMaxIts", "Maximum iterations on hot start",
      0, COIN_INT_MAX);

  parameters_[CbcParam::MAXSAVEDSOLS]->setup(
      "maxSaved!Solutions", "Maximum number of solutions to save", 0,
      COIN_INT_MAX, "Number of solutions to save.");

  parameters_[CbcParam::MAXSLOWCUTS]->setup(
      "slow!cutpasses", "Maximum number of rounds for slower cut generators",
      -1, COIN_INT_MAX, 
      "Some cut generators are fairly slow - this limits the number of times "
      "they are tried. The cut generators identified as 'may be slow' at "
      "present are Lift and project cuts and both versions of Reduce and Split "
      "cuts.");

  parameters_[CbcParam::MOREMOREMIPOPTIONS]->setup(
      "more2!MipOptions", "More more dubious options for mip", -1, COIN_INT_MAX,
      "", CoinParam::displayPriorityNone);

  parameters_[CbcParam::MULTIPLEROOTS]->setup(
      "multiple!RootPasses",
      "Do multiple root passes to collect cuts and solutions", 0, COIN_INT_MAX,
      "Solve (in parallel, if enabled) the root phase this number of times, "
      "each with its own different seed, and collect all solutions and cuts "
      "generated. The actual format is aabbcc where aa is the number of extra "
      "passes; if bb is non zero, then it is number of threads to use "
      "(otherwise uses threads setting); and cc is the number of times to do "
      "root phase. The solvers do not interact with each other.  However if "
      "extra passes are specified then cuts are collected and used in later "
      "passes - so there is interaction there. Some parts of this "
      "implementation have their origin in idea of Andrea Lodi, Matteo "
      "Fischetti, Michele Monaci, Domenico Salvagnin, and Andrea Tramontani.",
      CoinParam::displayPriorityNone);

  parameters_[CbcParam::ODDWEXTMETHOD]->setup(
      "oddwext!method",
      "Strategy used to search for wheel centers for the cuts found by Odd "
      "Wheel Cut Separation routine",
      0, 2,
      "Sets the method used in the extension module of Odd Wheel Cut "
      "Separation routine: 0=no extension; 1=one variable; 2=clique");

  parameters_[CbcParam::OUTPUTFORMAT]->setup(
      "output!Format", "Which output format to use", 1, 6,
      "Normally export will be done using normal representation for numbers "
      "and two values per line.  You may want to do just one per line (for "
      "grep or suchlike) and you may wish to save with absolute accuracy using "
      "a coded version of the IEEE value. A value of 2 is normal. Otherwise, "
      "odd values give one value per line, even values two.  Values of 1 and 2 "
      "give normal format, 3 and 4 give greater precision, 5 and 6 give IEEE "
      "values.  When exporting a basis, 1 does not save values, 2 saves "
      "values, 3 saves with greater accuracy and 4 saves in IEEE format.");

  parameters_[CbcParam::PRINTOPTIONS]->setup(
      "pO!ptions", "Dubious print options", 0, COIN_INT_MAX,
      "If this is greater than 0 then presolve will give more information and "
      "branch and cut will give statistics");

  parameters_[CbcParam::PROCESSTUNE]->setup(
      "tune!PreProcess", "Dubious tuning parameters for preprocessing", 0,
      COIN_INT_MAX,
      "Format aabbcccc - \n If aa then this is number of major passes (i.e. "
      "with presolve) \n If bb and bb>0 then this is number of minor passes "
      "(if unset or 0 then 10) \n cccc is bit set \n 0 - 1 Heavy probing \n 1 "
      "- 2 Make variables integer if possible (if obj value)\n 2 - 4 As above "
      "but even if zero objective value\n 7 - 128 Try and create cliques\n 8 - "
      "256 If all +1 try hard for dominated rows\n 9 - 512 Even heavier "
      "probing \n 10 - 1024 Use a larger feasibility tolerance in presolve\n "
      "11 - 2048 Try probing before creating cliques\n 12 - 4096 Switch off "
      "duplicate column checking for integers \n 13 - 8192 Allow scaled "
      "duplicate column checking \n \n     Now aa 99 has special meaning i.e. "
      "just one simple presolve.",
      CoinParam::displayPriorityLow);

  parameters_[CbcParam::RANDOMSEED]->setup(
      "randomC!bcSeed", "Random seed for Cbc", -1, COIN_INT_MAX,
      "Allows initialization of the random seed for pseudo-random numbers used "
      "in heuristics such as the Feasibility Pump to decide whether to round "
      "up or down. The special value of 0 lets Cbc use the time of the day for "
      "the initial seed.");

  parameters_[CbcParam::STRONGSTRATEGY]->setup(
      "expensive!Strong", "Whether to do even more strong branching", 0,
      COIN_INT_MAX,
      "Strategy for extra strong branching. 0 is normal strong branching. 1, "
      "2, 4, and 6 does strong branching on all fractional variables if at the "
      "root node (1), at depth less than modifier (2), objective equals best "
      "possible (4), or at depth less than modifier and objective equals best "
      "possible (6). 11, 12, 14, and 16 are like 1, 2, 4, and 6, "
      "respecitively, but do strong branching on all integer (incl. "
      "non-fractional) variables. Values >= 100 are used to specify a depth "
      "limit (value/100), otherwise 5 is used. If the values >= 100, then "
      "above rules are applied to value%100.",
      CoinParam::displayPriorityNone);

  parameters_[CbcParam::TESTOSI]->setup("testO!si", "Test OsiObject stuff",
                                            -1, COIN_INT_MAX, "",
                                            CoinParam::displayPriorityNone);

#ifdef CBC_THREAD
  parameters_[CbcParam::THREADS]->setup(
      "thread!s", "Number of threads to try and use", -100, 100000,
      "To use multiple threads, set threads to number wanted.  It may be "
      "better to use one or two more than number of cpus available.  If 100+n "
      "then n threads and search is repeatable (maybe be somewhat slower), if "
      "200+n use threads for root cuts, 400+n threads used in sub-trees.",
      CoinParam::displayPriorityLow);
#endif

  parameters_[CbcParam::USERCBC]->setup(
      "userCbc", "Hand coded Cbc stuff", 0, COIN_INT_MAX,
      "There are times (e.g., when using AMPL interface) when you may wish to "
      "do something unusual.  Look for USERCBC in main driver and modify "
      "sample code.",
      CoinParam::displayPriorityNone);

  parameters_[CbcParam::VERBOSE]->setup(
      "verbose", "Switches on longer help on single ?", 0, 15,
      "Set to 1 to get short help with ? list, 2 to get long help.",
      CoinParam::displayPriorityNone);

  parameters_[CbcParam::VUBTRY]->setup(
      "vub!heuristic", "Type of VUB heuristic", -2, 20,
      "This heuristic tries to fix some integer variables.",
      CoinParam::displayPriorityNone);
}

//###########################################################################
//###########################################################################

void CbcParameters::addCbcSolverBoolParams() {
  for (int code = CbcParam::FIRSTBOOLPARAM + 1;
       code < CbcParam::LASTBOOLPARAM; code++) {
    getParam(code)->setType(CoinParam::paramKwd);
    getParam(code)->appendKwd("off", CbcParameters::ParamOff);
    getParam(code)->appendKwd("on", CbcParameters::ParamOn);
    getParam(code)->setPushFunc(CbcParamUtils::pushCbcSolverKwdParam);
  }

  parameters_[CbcParam::CPX]->setup(
      "cplex!Use", "Whether to use Cplex!",
      "If the user has Cplex, but wants to use some of Cbc's heuristics then "
      "you can!  If this is on, then Cbc will get to the root node and then "
      "hand over to Cplex.  If heuristics find a solution this can be "
      "significantly quicker.  You will probably want to switch off Cbc's cuts "
      "as Cplex thinks they are genuine constraints.  It is also probable that "
      "you want to switch off preprocessing, although for difficult problems "
      "it is worth trying both.");

  parameters_[CbcParam::DOHEURISTIC]->setup(
      "doH!euristic", "Do heuristics before any preprocessing",
      "Normally heuristics are done in branch and bound.  It may be useful to "
      "do them outside. Only those heuristics with 'both' or 'before' set will "
      "run. Doing this may also set cutoff, which can help with "
      "preprocessing.");

  parameters_[CbcParam::ERRORSALLOWED]->setup(
      "error!sAllowed", "Whether to allow import errors",
      "The default is not to use any model which had errors when reading the "
      "mps file.  Setting this to 'on' will allow all errors from which the "
      "code can recover simply by ignoring the error.  There are some errors "
      "from which the code can not recover, e.g., no ENDATA.  This has to be "
      "set before import, i.e., -errorsAllowed on -import xxxxxx.mps.");

  parameters_[CbcParam::MESSAGES]->setup(
      "mess!ages", "Controls whether standardised message prefix is printed",
      "By default, messages have a standard prefix, such as:\n   Clp0005 2261  "
      "Objective 109.024 Primal infeas 944413 (758)\nbut this program turns "
      "this off to make it look more friendly.  It can be useful to turn them "
      "back on if you want to be able to 'grep' for particular messages or if "
      "you intend to override the behavior of a particular message.");

  parameters_[CbcParam::PREPROCNAMES]->setup(
      "PrepN!ames", "If column names will be kept in pre-processed model",
      "Normally the preprocessed model has column names replaced by new names "
      "C0000... Setting this option to on keeps original names in variables "
      "which still exist in the preprocessed problem");

  parameters_[CbcParam::SOS]->setup(
      "sos!Options", "Whether to use SOS from AMPL",
      "Normally if AMPL says there are SOS variables they should be used, but "
      "sometimes they should be turned off - this does so.");

  parameters_[CbcParam::USESOLUTION]->setup(
      "force!Solution", "Whether to use given solution as crash for BAB",
      "If on then tries to branch to solution given by AMPL or priorities "
      "file.");
}

//###########################################################################
//###########################################################################

void CbcParameters::addCbcSolverCutParams() {
  for (int code = CbcParam::FIRSTCUTPARAM + 1;
       code < CbcParam::LASTCUTPARAM; code++) {
    getParam(code)->setType(CoinParam::paramKwd);
    getParam(code)->setPushFunc(CbcParamUtils::pushCbcSolverKwdParam);
  }

  parameters_[CbcParam::CLIQUECUTS]->setup(
      "clique!Cuts", "Whether to use clique cuts",
      "This switches on clique cuts (either at root or in entire tree). See "
      "branchAndCut for information on options.");

  parameters_[CbcParam::CUTSTRATEGY]->setup(
      "cuts!OnOff", "Switches all cuts on or off", 
      "This can be used to switch on or off all cuts (apart from Reduce and "
      "Split).  Then you can set individual ones off or on.  See branchAndCut "
      "for information on options.");

  parameters_[CbcParam::FLOWCUTS]->setup(
      "flow!CoverCuts", "Whether to use Flow Cover cuts",
      "This switches on flow cover cuts (either at root or in entire "
      "tree)." CUTS_LONGHELP);

  parameters_[CbcParam::GMICUTS]->setup(
      "GMI!Cuts", "Whether to use alternative Gomory cuts", 
      CUTS_LONGHELP " This version is by Giacomo Nannicini and may be more "
                    "robust than gomoryCuts.");

  parameters_[CbcParam::GOMORYCUTS]->setup(
      "gomory!Cuts", "Whether to use Gomory cuts",
      "The original cuts - beware of imitations!  Having gone out of favor, "
      "they are now more fashionable as LP solvers are more robust and they "
      "interact well with other cuts.  They will almost always give cuts "
      "(although in this executable they are limited as to number of variables "
      "in cut).  However the cuts may be dense so it is worth experimenting "
      "(Long allows any length). " CUTS_LONGHELP
      " Reference: https://github.com/coin-or/Cgl/wiki/CglGomory");

  parameters_[CbcParam::KNAPSACKCUTS]->setup(
      "knapsack!Cuts", "Whether to use Knapsack cuts",
      "This switches on knapsack cuts (either at root or in entire tree). See "
      "branchAndCut for information on options.");

  parameters_[CbcParam::LAGOMORYCUTS]->setup(
      "lagomory!Cuts", "Whether to use Lagrangean Gomory cuts", 
      "This is a gross simplification of 'A Relax-and-Cut Framework for "
      "Gomory's Mixed-Integer Cuts' by Matteo Fischetti & Domenico Salvagnin.  "
      "This simplification just uses original constraints while modifying "
      "objective using other cuts. So you don't use messy constraints "
      "generated by Gomory etc. A variant is to allow non messy cuts e.g. "
      "clique cuts. So 'only' does this while 'clean' also allows integral "
      "valued cuts.  'End' is recommended and waits until other cuts have "
      "finished before it does a few passes. The length options for gomory "
      "cuts are used.");

  parameters_[CbcParam::LANDPCUTS]->setup(
      "lift!AndProjectCuts", "Whether to use lift-and-project cuts",
      "This switches on lift-and-project cuts (either at root or in entire "
      "tree). See branchAndCut for information on options.");

  parameters_[CbcParam::LATWOMIRCUTS]->setup(
      "latwomir!Cuts", "Whether to use Lagrangean Twomir cuts",
      "This is a Lagrangean relaxation for Twomir cuts.  See lagomoryCuts for "
      "description of options.");

  parameters_[CbcParam::MIRCUTS]->setup(
      "mixed!IntegerRoundingCuts", "Whether to use Mixed Integer Rounding cuts",
      "This switches on mixed integer rounding cuts (either at root or in "
      "entire tree).  See branchAndCut for information on options.");

  parameters_[CbcParam::ODDWHEELCUTS]->setup(
      "oddwheel!Cuts", "Whether to use odd wheel cuts", 
      "This switches on odd-wheel inequalities (either at root or in entire "
      "tree).");

  parameters_[CbcParam::PROBINGCUTS]->setup(
      "probing!Cuts", "Whether to use Probing cuts", 
      "This switches on probing cuts (either at root or in entire tree). See "
      "branchAndCut for information on options.");

  parameters_[CbcParam::REDSPLITCUTS]->setup(
      "reduce!AndSplitCuts", "Whether to use Reduce-and-Split cuts",
      "This switches on reduce and split cuts (either at root or in entire "
      "tree). See branchAndCut for information on options.");

  parameters_[CbcParam::REDSPLIT2CUTS]->setup(
      "reduce2!AndSplitCuts", "Whether to use Reduce-and-Split cuts - style 2",
      "This switches on reduce and split cuts (either at root or in entire "
      "tree). See branchAndCut for information on options.");

  parameters_[CbcParam::RESIDCAPCUTS]->setup(
      "residual!CapacityCuts", "Whether to use Residual Capacity cuts",
      CUTS_LONGHELP
      " Reference: https://github.com/coin-or/Cgl/wiki/CglResidualCapacity");

  parameters_[CbcParam::TWOMIRCUTS]->setup(
      "two!MirCuts", "Whether to use Two phase Mixed Integer Rounding cuts",
      "This switches on two phase mixed integer rounding cuts (either at root "
      "or in entire tree). See branchAndCut for information on options.");

  parameters_[CbcParam::ZEROHALFCUTS]->setup(
      "zero!HalfCuts", "Whether to use zero half cuts", 
      CUTS_LONGHELP " This implementation was written by Alberto Caprara.");

  // Populate the keyword lists
  for (int code = CbcParam::FIRSTCUTPARAM + 1;
       code < CbcParam::LASTCUTPARAM; code++) {
    // First the common keywords
    switch (code) {
     case CbcParam::CUTSTRATEGY:
     case CbcParam::CLIQUECUTS:
     case CbcParam::FLOWCUTS:
     case CbcParam::GMICUTS:
     case CbcParam::GOMORYCUTS:
     case CbcParam::KNAPSACKCUTS:
     case CbcParam::LANDPCUTS:
     case CbcParam::MIRCUTS:
     case CbcParam::ODDWHEELCUTS:
     case CbcParam::PROBINGCUTS:
     case CbcParam::RESIDCAPCUTS:
     case CbcParam::TWOMIRCUTS:
     case CbcParam::ZEROHALFCUTS:
        parameters_[code]->appendKwd("off", CbcParameters::CGOff);
        parameters_[code]->appendKwd("on", CbcParameters::CGOn);
        parameters_[code]->appendKwd("root", CbcParameters::CGRoot);
        parameters_[code]->appendKwd("ifmove", CbcParameters::CGIfMove);
        parameters_[code]->appendKwd("forceon", CbcParameters::CGForceOn);
        break;
     default: 
        break;
    }
    
    // Now, add some additional keywords for different classes
    switch (code) {
     case CbcParam::GOMORYCUTS: 
        parameters_[code]->appendKwd("onglobal", CbcParameters::CGOnGlobal);
        parameters_[code]->appendKwd("forceandglobal", CbcParameters::CGForceAndGlobal);
        parameters_[code]->appendKwd("forcelongon", CbcParameters::CGForceLongOn);
        parameters_[code]->appendKwd("longer", CbcParameters::CGLonger);
        parameters_[code]->appendKwd("shorter", CbcParameters::CGShorter);
        break;
     case CbcParam::GMICUTS: 
        parameters_[code]->appendKwd("long", CbcParameters::CGLong);
        parameters_[code]->appendKwd("longroot", CbcParameters::CGLongRoot);
        parameters_[code]->appendKwd("longifmove", CbcParameters::CGLongIfMove);
        parameters_[code]->appendKwd("forcelongon", CbcParameters::CGForceLongOn);
        parameters_[code]->appendKwd("longendonly", CbcParameters::CGLongEndOnly);
        break;
     case CbcParam::LAGOMORYCUTS: 
        parameters_[code]->appendKwd("root", CbcParameters::CGRoot);
        parameters_[code]->appendKwd("onlyaswellroot", CbcParameters::CGOnlyAsWellRoot);
        parameters_[code]->appendKwd("cleanaswellroot", CbcParameters::CGCleanAsWellRoot);
        parameters_[code]->appendKwd("bothaswellroot", CbcParameters::CGCleanBothAsWellRoot);
        // Here, we intentionally drop through to the next set
     case CbcParam::LATWOMIRCUTS: 
        parameters_[code]->appendKwd("endonlyroot", CbcParameters::CGEndOnlyRoot);
        parameters_[code]->appendKwd("endcleanroot", CbcParameters::CGEndCleanRoot);
        parameters_[code]->appendKwd("endonly", CbcParameters::CGEndOnly);
        parameters_[code]->appendKwd("endclean", CbcParameters::CGEndClean);
        parameters_[code]->appendKwd("endboth", CbcParameters::CGEndBoth);
        parameters_[code]->appendKwd("onlyaswell", CbcParameters::CGOnlyAsWell);
        parameters_[code]->appendKwd("cleanaswell", CbcParameters::CGCleanAsWell);
        parameters_[code]->appendKwd("bothaswell", CbcParameters::CGBothAsWell);
        parameters_[code]->appendKwd("onlyinstead", CbcParameters::CGOnlyInstead);
        parameters_[code]->appendKwd("cleaninstead", CbcParameters::CGCleanInstead);
        parameters_[code]->appendKwd("bothinstead", CbcParameters::CGBothInstead);
        break;
     case CbcParam::ODDWHEELCUTS: 
        parameters_[code]->appendKwd("onglobal", CbcParameters::CGOnGlobal);
        break;
     case CbcParam::PROBINGCUTS: 
        parameters_[code]->appendKwd("forceonbut", CbcParameters::CGForceOnBut);
        break;
     case CbcParam::REDSPLITCUTS:
     case CbcParam::REDSPLIT2CUTS: 
        parameters_[code]->appendKwd("on", CbcParameters::CGOn);
        parameters_[code]->appendKwd("root", CbcParameters::CGRoot);
        parameters_[code]->appendKwd("longon", CbcParameters::CGLongOn);
        parameters_[code]->appendKwd("longroot", CbcParameters::CGLongRoot);
        break;
     default:
       break;
    }
  }
}

//###########################################################################
//###########################################################################

void CbcParameters::addCbcSolverHeurParams() {
  for (int code = CbcParam::FIRSTHEURPARAM + 1;
       code < CbcParam::LASTHEURPARAM; code++) {
    getParam(code)->setType(CoinParam::paramKwd);
    getParam(code)->setPushFunc(CbcParamUtils::pushCbcSolverKwdParam);
  }

  parameters_[CbcParam::COMBINE]->setup(
      "combine!Solutions", "Whether to use combine solution heuristic", 
      "This switches on a heuristic which does branch and cut on the problem "
      "given by just using variables which have appeared in one or more "
      "solutions. It is obviously only tried after two or more solutions.");

  parameters_[CbcParam::CROSSOVER]->setup(
      "combine2!Solutions", "Whether to use crossover solution heuristic",
      HEURISTICS_LONGHELP);

  parameters_[CbcParam::DINS]->setup(
      "Dins", "Whether to try Distance Induced Neighborhood Search", HEURISTICS_LONGHELP);

  parameters_[CbcParam::DIVINGC]->setup(
      "DivingC!oefficient", "Whether to try Coefficient diving heuristic",
      HEURISTICS_LONGHELP);

  parameters_[CbcParam::DIVINGF]->setup(
      "DivingF!ractional", "Whether to try Fractional diving heuristic", HEURISTICS_LONGHELP);

  parameters_[CbcParam::DIVINGG]->setup(
      "DivingG!uided", "Whether to try Guided diving heuristic", HEURISTICS_LONGHELP);

  parameters_[CbcParam::DIVINGL]->setup(
      "DivingL!ineSearch", "Whether to try Linesearch diving heuristic", HEURISTICS_LONGHELP);

  parameters_[CbcParam::DIVINGP]->setup(
      "DivingP!seudocost", "Whether to try Pseudocost diving heuristic", HEURISTICS_LONGHELP);

  parameters_[CbcParam::DIVINGS]->setup(
      "DivingS!ome", "Whether to try Diving heuristics", 
      "This switches on a random diving heuristic at various times. One may "
      "prefer to individually turn diving heuristics on or off. ");

  parameters_[CbcParam::DIVINGV]->setup(
      "DivingV!ectorLength", "Whether to try Vectorlength diving heuristic",
      HEURISTICS_LONGHELP);

  parameters_[CbcParam::DW]->setup(
      "dw!Heuristic", "Whether to try Dantzig Wolfe heuristic", 
      "This heuristic is very very compute intensive. It tries to find a "
      "Dantzig Wolfe structure and use that. " HEURISTICS_LONGHELP);

  parameters_[CbcParam::FPUMP]->setup(
      "feas!ibilityPump", "Whether to try Feasibility Pump", 
      "This switches on feasibility pump heuristic at root. This is due to "
      "Fischetti and Lodi and uses a sequence of LPs to try and get an integer "
      "feasible solution.  Some fine tuning is available by "
      "passFeasibilityPump.");

  parameters_[CbcParam::GREEDY]->setup(
      "greedy!Heuristic", "Whether to use a greedy heuristic",
      "Switches on a pair of greedy heuristic which will try and obtain a "
      "solution.  It may just fix a percentage of variables and then try a "
      "small branch and cut run.");

  parameters_[CbcParam::HEURISTICSTRATEGY]->setup(
      "heur!isticsOnOff", "Switches most heuristics on or off", 
      "This can be used to switch on or off all heuristics.  Then you can set "
      "individual ones off or on.  CbcTreeLocal is not included as it "
      "dramatically alters search.");

  parameters_[CbcParam::LOCALTREE]->setup(
      "local!TreeSearch", "Whether to use local tree search", 
      "This switches on a local search algorithm when a solution is found.  "
      "This is from Fischetti and Lodi and is not really a heuristic although "
      "it can be used as one. When used from this program it has limited "
      "functionality.  It is not controlled by heuristicsOnOff.");

  parameters_[CbcParam::NAIVE]->setup(
      "naive!Heuristics", "Whether to try some stupid heuristic", 
      "This is naive heuristics which, e.g., fix all integers with costs to "
      "zero!. " HEURISTICS_LONGHELP,
      CoinParam::displayPriorityLow);

  parameters_[CbcParam::PIVOTANDFIX]->setup(
      "pivotAndF!ix", "Whether to try Pivot and Fix heuristic", HEURISTICS_LONGHELP);

#if 0
  parameters_[CbcParam::PIVOTANDCOMPLEMENT]->setup(
      "pivotAndC!omplement", "Whether to try Pivot and Complement heuristic",
      HEURISTICS_LONGHELP);
#endif

  parameters_[CbcParam::PROXIMITY]->setup(
      "proximity!Search", "Whether to do proximity search heuristic", 
      "This heuristic looks for a solution close to the incumbent solution "
      "(Fischetti and Monaci, 2012). The idea is to define a sub-MIP without "
      "additional constraints but with a modified objective function intended "
      "to attract the search in the proximity of the incumbent. The approach "
      "works well for 0-1 MIPs whose solution landscape is not too irregular "
      "(meaning the there is reasonable probability of finding an improved "
      "solution by flipping a small number of binary variables), in particular "
      "when it is applied to the first heuristic solutions found at the root "
      "node. " HEURISTICS_LONGHELP); // Can also set different maxNode
                                     // cbcSettings by plusnnnn (and are
                                     // 'on'(on==30)).

  parameters_[CbcParam::RANDROUND]->setup(
      "randomi!zedRounding", "Whether to try randomized rounding heuristic",
      HEURISTICS_LONGHELP);

  parameters_[CbcParam::RENS]->setup(
      "Rens", "Whether to try Relaxation Enforced Neighborhood Search", 
      HEURISTICS_LONGHELP " Value 'on' just does 50 nodes. 200, 1000, and "
                          "10000 does that many nodes.");

  parameters_[CbcParam::RINS]->setup(
      "Rins", "Whether to try Relaxed Induced Neighborhood Search",
      HEURISTICS_LONGHELP);

  parameters_[CbcParam::ROUNDING]->setup(
      "round!ingHeuristic", "Whether to use Rounding heuristic", 
      "This switches on a simple (but effective) rounding heuristic at each "
      "node of tree.");

  parameters_[CbcParam::VND]->setup(
      "Vnd!VariableNeighborhoodSearch",
      "Whether to try Variable Neighborhood Search", HEURISTICS_LONGHELP);

  // Populate the keyword lists
  for (int code = CbcParam::FIRSTHEURPARAM + 1;
       code < CbcParam::LASTHEURPARAM; code++) {
    // First the common keywords
    switch (code) {
     case CbcParam::HEURISTICSTRATEGY:
     case CbcParam::CROSSOVER:
     case CbcParam::DINS:
     case CbcParam::DIVINGC:
     case CbcParam::DIVINGF:
     case CbcParam::DIVINGG:
     case CbcParam::DIVINGL:
     case CbcParam::DIVINGP:
     case CbcParam::DIVINGS:
     case CbcParam::DIVINGV:
     case CbcParam::DW:
     case CbcParam::NAIVE:
     case CbcParam::PIVOTANDFIX:
     case CbcParam::PIVOTANDCOMPLEMENT:
     case CbcParam::PROXIMITY:
     case CbcParam::RANDROUND:
     case CbcParam::RENS:
     case CbcParam::RINS:
     case CbcParam::ROUNDING:
     case CbcParam::VND: 
       parameters_[code]->appendKwd("off", CbcParameters::CGOff);
       parameters_[code]->appendKwd("on", CbcParameters::CGOn);
       parameters_[code]->appendKwd("both", CbcParameters::CGRoot);
       parameters_[code]->appendKwd("before", CbcParameters::CGIfMove);
       break;
     default:
       break;
    }
    // Check the unique keywords
    switch (code) {
     case CbcParam::COMBINE:
        parameters_[code]->appendKwd("off", CbcParameters::HeurOff);
        parameters_[code]->appendKwd("on", CbcParameters::HeurOn);
        break;
     case CbcParam::DINS:
        parameters_[code]->appendKwd("often", CbcParameters::HeurOften);
        break;
     case CbcParam::FPUMP:
        parameters_[code]->appendKwd("off", CbcParameters::HeurOff);
        parameters_[code]->appendKwd("on", CbcParameters::HeurOn);
        break;
     case CbcParam::GREEDY: 
        parameters_[code]->appendKwd("off", CbcParameters::HeurOff);
        parameters_[code]->appendKwd("on", CbcParameters::HeurOn);
        parameters_[code]->appendKwd("root", CbcParameters::HeurRoot);
        break;
     case CbcParam::LOCALTREE: 
        parameters_[code]->appendKwd("off", CbcParameters::HeurOff);
        parameters_[code]->appendKwd("on", CbcParameters::HeurOn);
     case CbcParam::PROXIMITY: 
        parameters_[code]->appendKwd("10", CbcParameters::HeurTen);
        parameters_[code]->appendKwd("100", CbcParameters::HeurOneHundred);
        parameters_[code]->appendKwd("300", CbcParameters::HeurThreeHundred);
        break;
     case CbcParam::RENS: 
        parameters_[code]->appendKwd("200", CbcParameters::HeurTwoHundred);
        parameters_[code]->appendKwd("1000", CbcParameters::HeurOneThousand);
        parameters_[code]->appendKwd("10000", CbcParameters::HeurTenThousand);
        parameters_[code]->appendKwd("dj", CbcParameters::HeurDj);
        parameters_[code]->appendKwd("djbefore", CbcParameters::HeurDjBefore);
        parameters_[code]->appendKwd("usesolution", CbcParameters::HeurUseSolution);
        break;
     case CbcParam::RINS: 
        parameters_[code]->appendKwd("often", CbcParameters::HeurOften);
        break;
     case CbcParam::VND: 
        parameters_[code]->appendKwd("intree", CbcParameters::HeurInTree);
        break;
     default:
       break;
    }
  }
}

//###########################################################################
//###########################################################################

/* Function to set up cbc (CbcModel) parameters.  */

void CbcParameters::addCbcModelParams()
{

  parameters_[CbcParam::ALLOWABLEGAP]->setup(
      "allow!ableGap",
      "Stop when gap between best possible and incumbent is less than this",
      0.0, 1.0e20, 
      "If the gap between best solution and best possible solution is less "
      "than this then the search will be terminated. Also see ratioGap.");

  parameters_[CbcParam::CUTOFF]->setup(
      "cuto!ff", "All solutions must be better than this", -1.0e60, 1.0e60,
      "All solutions must be better than this value (in a minimization sense). "
      " This is also set by cbc whenever it obtains a solution and is set to "
      "the value of the objective for the solution minus the cutoff "
      "increment.");

  parameters_[CbcParam::DIRECTION]->setup(
      "direction", "Minimize or maximize", 
      "The default is minimize - use 'direction maximize' for "
      "maximization.\nYou can also use the parameters_ 'maximize' or "
      "'minimize'.");
  parameters_[CbcParam::DIRECTION]->setType(CoinParam::paramKwd);
  parameters_[CbcParam::DIRECTION]->appendKwd("min!imize", CbcParameters::OptDirMinimize);
  parameters_[CbcParam::DIRECTION]->appendKwd("max!imize", CbcParameters::OptDirMaximize);
  parameters_[CbcParam::DIRECTION]->appendKwd("zero", CbcParameters::OptDirZero);
  
  parameters_[CbcParam::INCREMENT]->setup(
      "inc!rement",
      "A new solution must be at least this much better than the incumbent",
      -1.0e20, 1.0e20, 
      "Whenever a solution is found the bound on future solutions is set to "
      "the objective of the solution (in a minimization sense) plus the "
      "specified increment.  If this option is not specified, the code will "
      "try and work out an increment.  E.g., if all objective coefficients are "
      "multiples of 0.01 and only integer variables have entries in objective "
      "then the increment can be set to 0.01.  Be careful if you set this "
      "negative!");

  parameters_[CbcParam::INFEASIBILITYWEIGHT]->setup(
      "inf!easibilityWeight",
      "Each integer infeasibility is expected to cost this much", 0.0, 1.0e20,
      "A primitive way of deciding which node to explore next.  Satisfying "
      "each integer infeasibility is expected to cost this much.");

  parameters_[CbcParam::INTEGERTOLERANCE]->setup(
      "integerT!olerance",
      "For an optimal solution, no integer variable may be farther than this "
      "from an integer value",
      1.0e-20, 0.5, 
      "When checking a solution for feasibility, if the difference between the "
      "value of a variable and the nearest integer is less than the integer "
      "tolerance, the value is considered to be integral. Beware of setting "
      "this smaller than the primal tolerance.");

  parameters_[CbcParam::LOGLEVEL]->setup(
      "bclog!Level", "Level of detail in Coin branch and Cut output", -1, 63,
      "If set to 0 then there should be no output in normal circumstances. A "
      "value of 1 is probably the best value for most uses, while 2 and 3 give "
      "more information.");

  parameters_[CbcParam::MAXIMIZE]->setup(
      "max!imize", "Set optimization direction to maximize",
      "The default is minimize - use 'maximize' for maximization.\n A synonym "
      "for 'direction maximize'.",
      CoinParam::displayPriorityHigh);
  parameters_[CbcParam::MAXIMIZE]->setType(CoinParam::paramAct);

  parameters_[CbcParam::MAXNODES]->setup(
      "maxN!odes", "Maximum number of nodes to evaluate", 1, COIN_INT_MAX, 
      "This is a repeatable way to limit search.  Normally using time is "
      "easier but then the results may not be repeatable.");

  parameters_[CbcParam::MAXNODESNOTIMPROVING]->setup(
      "maxNNI!FS",
      "Maximum number of nodes to be processed without improving the incumbent "
      "solution.",
      -1, COIN_INT_MAX, 
      "This criterion specifies that when a feasible solution is available, "
      "the search should continue only if better feasible solutions were "
      "produced in the last nodes.");

  parameters_[CbcParam::MAXSECONDSNOTIMPROVING]->setup(
      "secni!fs", "maximum seconds without improving the incumbent solution",
      -1.0, COIN_DBL_MAX, 
      "With this stopping criterion, after a feasible solution is found, the "
      "search should continue only if the incumbent solution was updated "
      "recently, the tolerance is specified here. A discussion on why this "
      "criterion can be useful is included here: "
      "https://yetanothermathprogrammingconsultant.blogspot.com/2019/11/"
      "mip-solver-stopping-criteria.html .");

  parameters_[CbcParam::MAXSOLS]->setup(
      "maxSo!lutions", "Maximum number of feasible solutions to get", 1,
      COIN_INT_MAX, 
      "You may want to stop after (say) two solutions or an hour. This is "
      "checked every node in tree, so it is possible to get more solutions "
      "from heuristics.");

  parameters_[CbcParam::MINIMIZE]->setup(
      "min!imize", "Set optimization direction to minimize",
      "The default is minimize - use 'maximize' for maximization.\nThis should "
      "only be necessary if you have previously set maximization. A synonym "
      "for 'direction minimize'.");
  parameters_[CbcParam::MINIMIZE]->setType(CoinParam::paramAct);

  parameters_[CbcParam::MIPOPTIONS]->setup(
      "mipO!ptions", "Dubious options for mip", 0, COIN_INT_MAX, "",
      CoinParam::displayPriorityNone);

  parameters_[CbcParam::MOREMIPOPTIONS]->setup(
      "more!MipOptions", "More dubious options for mip", -1, COIN_INT_MAX, 
      "", CoinParam::displayPriorityNone);

#if 0
  parameters_[CbcParam::NUMBERMINI]->setup("miniT!ree",
      "Size of fast mini tree", 0, COIN_INT_MAX,
      "The idea is that I can do a small tree fast. This is a first try and will" 
      "hopefully become more sophisticated.", CoinParam::displayPriorityNone);
#endif

  parameters_[CbcParam::NUMBERANALYZE]->setup(
      "numberA!nalyze", "Number of analysis iterations", -COIN_INT_MAX,
      COIN_INT_MAX, 
      "This says how many iterations to spend at the root node analyzing the "
      "problem.  This is a first try and will hopefully become more "
      "sophisticated.",
      CoinParam::displayPriorityNone);

  parameters_[CbcParam::REVERSE]->setup(
      "reverse", "Reverses sign of objective",
      "Useful for testing if maximization works correctly",
      CoinParam::displayPriorityNone);
  parameters_[CbcParam::REVERSE]->setType(CoinParam::paramAct);

  parameters_[CbcParam::CUTPASS]->setup(
      "passC!uts", "Number of cut passes at root node", -999999, 999999,
      "The default is 100 passes if less than 500 columns, 100 passes (but "
      "stop if the drop is small) if less than 5000 columns, 20 otherwise.");

  parameters_[CbcParam::GAPRATIO]->setup(
      "ratio!Gap",
      "Stop when the gap between the best possible solution and the incumbent "
      "is less than this fraction of the larger of the two",
      0.0, 1.0e20,
      "If the gap between the best solution and the best possible solution is "
      "less than this fraction of the objective value at the root node then "
      "the search will terminate.  See 'allowableGap' for a way of using "
      "absolute value rather than fraction.");

  parameters_[CbcParam::TIMELIMIT]->setup(
      "sec!onds", "Maximum seconds for branch and cut", -1.0, 1.0e12,
      "After this many seconds the program will act as if maximum nodes had "
      "been reached.");

  parameters_[CbcParam::STRONGBRANCHING]->setup(
      "strong!Branching", "Number of variables to look at in strong branching",
      0, 999999,
      "In order to decide which variable to branch on, the code will choose up "
      "to this number of unsatisfied variables and try mini up and down "
      "branches.  The most effective one is chosen. If a variable is branched "
      "on many times then the previous average up and down costs may be used - "
      "see number before trust.");

  parameters_[CbcParam::NUMBERBEFORE]->setup(
      "trust!PseudoCosts", "Number of branches before we trust pseudocosts", -1,
      2000000,
      "Using strong branching computes pseudo-costs.  After this many times "
      "for a variable we just trust the pseudo costs and do not do any more "
      "strong branching.");

  for (int code = CbcParam::FIRSTMODELPARAM + 1;
       code < CbcParam::LASTMODELPARAM; code++) {
     if (getParam(code)->type() == CoinParam::paramInt){
        getParam(code)->setPushFunc(CbcParamUtils::pushCbcModelIntParam);
     }else{
        getParam(code)->setPushFunc(CbcParamUtils::pushCbcModelDblParam);
     }
  }
}

//###########################################################################
//###########################################################################

/*
  Access functions for cut generators. These support lazy
  creation --- if mode_ is other than CGOff, an object is created if
  necessary and a pointer is stored in proto_. The pointer is returned as
  a generic CglCutGenerator or CbcHeuristic. The return value of the function
  is the value of mode_.

  Because the model may have changed, the default for heuristics is to delete
  any existing object and create a new one. This can be suppressed if desired.
*/

CbcParameters::CGMode CbcParameters::getClique(CglCutGenerator *&gen) {
  if (clique_.mode_ != CbcParameters::CGOff && clique_.proto_ == 0) {
    clique_.proto_ = new CglClique();
    clique_.proto_->setStarCliqueReport(clique_.starCliqueReport_);
    clique_.proto_->setRowCliqueReport(clique_.rowCliqueReport_);
    clique_.proto_->setMinViolation(clique_.minViolation_);
  }
  gen = dynamic_cast<CglCutGenerator *>(clique_.proto_);

  return (clique_.mode_);
}

CbcParameters::CGMode CbcParameters::getFlow(CglCutGenerator *&gen)

{
  if (flow_.mode_ != CbcParameters::CGOff && flow_.proto_ == 0) {
    flow_.proto_ = new CglFlowCover();
  }
  gen = dynamic_cast<CglCutGenerator *>(flow_.proto_);

  return (flow_.mode_);
}

CbcParameters::CGMode CbcParameters::getGomory(CglCutGenerator *&gen) {
  if (gomory_.mode_ != CbcParameters::CGOff && gomory_.proto_ == 0) {
    gomory_.proto_ = new CglGomory();
    gomory_.proto_->setLimitAtRoot(gomory_.limitAtRoot_);
    gomory_.proto_->setLimit(gomory_.limit_);
  }
  gen = dynamic_cast<CglCutGenerator *>(gomory_.proto_);

  return (gomory_.mode_);
}

CbcParameters::CGMode CbcParameters::getKnapsack(CglCutGenerator *&gen) {
  if (knapsack_.mode_ != CbcParameters::CGOff && knapsack_.proto_ == 0) {
    knapsack_.proto_ = new CglKnapsackCover();
  }
  gen = dynamic_cast<CglCutGenerator *>(knapsack_.proto_);

  return (knapsack_.mode_);
}

CbcParameters::CGMode CbcParameters::getMir(CglCutGenerator *&gen) {
  if (mir_.mode_ != CbcParameters::CGOff && mir_.proto_ == 0) {
    mir_.proto_ = new CglMixedIntegerRounding2();
  }
  gen = dynamic_cast<CglCutGenerator *>(mir_.proto_);

  return (mir_.mode_);
}

CbcParameters::CGMode CbcParameters::getProbing(CglCutGenerator *&gen) {
  if (probing_.mode_ != CbcParameters::CGOff && probing_.proto_ == 0) {
    probing_.proto_ = new CglProbing();
    probing_.proto_->setUsingObjective(probing_.usingObjective_);
    probing_.proto_->setMaxPass(probing_.maxPass_);
    probing_.proto_->setMaxPassRoot(probing_.maxPassRoot_);
    probing_.proto_->setMaxProbe(probing_.maxProbe_);
    probing_.proto_->setMaxProbeRoot(probing_.maxProbeRoot_);
    probing_.proto_->setMaxLook(probing_.maxLook_);
    probing_.proto_->setMaxLookRoot(probing_.maxLookRoot_);
    probing_.proto_->setMaxElements(probing_.maxElements_);
    probing_.proto_->setRowCuts(probing_.rowCuts_);
  }
  gen = dynamic_cast<CglCutGenerator *>(probing_.proto_);

  return (probing_.mode_);
}

CbcParameters::CGMode CbcParameters::getRedSplit(CglCutGenerator *&gen)

{
  if (redSplit_.mode_ != CbcParameters::CGOff && redSplit_.proto_ == 0) {
    redSplit_.proto_ = new CglRedSplit();
  }
  gen = dynamic_cast<CglCutGenerator *>(redSplit_.proto_);

  return (redSplit_.mode_);
}

CbcParameters::CGMode CbcParameters::getTwomir(CglCutGenerator *&gen) {
  if (twomir_.mode_ != CbcParameters::CGOff && twomir_.proto_ == 0) {
    twomir_.proto_ = new CglTwomir();
    twomir_.proto_->setMaxElements(twomir_.maxElements_);
  }
  gen = dynamic_cast<CglCutGenerator *>(twomir_.proto_);

  return (twomir_.mode_);
}

CbcParameters::HeurMode CbcParameters::getFeasPump(CbcHeuristic *&gen,
                                                        bool alwaysCreate)

{
  if (fpump_.mode_ != CbcParameters::HeurOff &&
      (fpump_.proto_ == 0 || alwaysCreate)) {
    if (fpump_.proto_) {
      delete fpump_.proto_;
    }
    fpump_.proto_ = new CbcHeuristicFPump(*model_);
    fpump_.proto_->setMaximumPasses(fpump_.iters_);
  }
  gen = dynamic_cast<CbcHeuristic *>(fpump_.proto_);

  return (fpump_.mode_);
}

/*
  Access functions for heuristics. These support lazy
  creation --- if mode_ is other than HeurOff, an object is created if
  necessary and a pointer is stored in proto_. The pointer is returned as
  a generic CbcHeuristic. The return value of the function
  is the value of mode_.

  Because the model may have changed, the default for heuristics is to delete
  any existing object and create a new one. This can be suppressed if desired.
*/

CbcParameters::HeurMode CbcParameters::getCombine(CbcHeuristic *&gen,
                                                       bool alwaysCreate)

{
  if (combine_.mode_ != CbcParameters::HeurOff &&
      (combine_.proto_ == 0 || alwaysCreate)) {
    if (combine_.proto_) {
      delete combine_.proto_;
    }
    //TODO Should we be passing a pointer here? Otherwise, making a copy
    combine_.proto_ = new CbcHeuristicLocal(*model_);
    combine_.proto_->setSearchType(combine_.trySwap_);
  }
  gen = dynamic_cast<CbcHeuristic *>(combine_.proto_);

  return (combine_.mode_);
}

CbcParameters::HeurMode CbcParameters::getGreedyCover(CbcHeuristic *&gen,
                                                           bool alwaysCreate)

{
  if (greedyCover_.mode_ != CbcParameters::HeurOff &&
      (greedyCover_.proto_ == 0 || alwaysCreate)) {
    if (greedyCover_.proto_) {
      delete greedyCover_.proto_;
    }
    greedyCover_.proto_ = new CbcHeuristicGreedyCover(*model_);
  }
  gen = dynamic_cast<CbcHeuristic *>(greedyCover_.proto_);

  return (greedyCover_.mode_);
}

CbcParameters::HeurMode
CbcParameters::getGreedyEquality(CbcHeuristic *&gen,
                                     bool alwaysCreate)

{
  if (greedyEquality_.mode_ != CbcParameters::HeurOff &&
      (greedyEquality_.proto_ == 0 || alwaysCreate)) {
    if (greedyEquality_.proto_) {
      delete greedyEquality_.proto_;
    }
    greedyEquality_.proto_ = new CbcHeuristicGreedyEquality(*model_);
  }
  gen = dynamic_cast<CbcHeuristic *>(greedyEquality_.proto_);

  return (greedyEquality_.mode_);
}

CbcParameters::HeurMode CbcParameters::getRounding(CbcHeuristic *&gen,
                                                        bool alwaysCreate)

{
  if (rounding_.mode_ != CbcParameters::HeurOff &&
      (rounding_.proto_ == 0 || alwaysCreate)) {
    if (rounding_.proto_) {
      delete rounding_.proto_;
    }
    rounding_.proto_ = new CbcRounding(*model_);
  }
  gen = dynamic_cast<CbcHeuristic *>(rounding_.proto_);

  return (rounding_.mode_);
}

CbcParameters::HeurMode
CbcParameters::getLocalTree(CbcTreeLocal *&localTree, 
                                bool alwaysCreate)

{
  if (localTree_.mode_ != CbcParameters::HeurOff &&
      (localTree_.proto_ == 0 || alwaysCreate)) {
    if (localTree_.proto_) {
      delete localTree_.proto_;
    }
    localTree_.proto_ = new CbcTreeLocal(
        model_, localTree_.soln_, localTree_.range_, localTree_.typeCuts_,
        localTree_.maxDiverge_, localTree_.timeLimit_, localTree_.nodeLimit_,
        localTree_.refine_);
  }
  localTree = localTree_.proto_;

  return (localTree_.mode_);
}

/*
  A bunch of little translation helper routines leading up to a version of
  setBaBStatus that figures it all out given a CbcModel and BACWhere code.
  This translation needs to be centralised to avoid sprinkling magic numbers
  all through the code.

  Be a bit careful with the translation routines --- they aren't sensitive to
  where the search stopped.
*/

CbcParameters::BACMajorStatus CbcParameters::translateMajor(int status)

{
  switch (status) {
  case -1: {
    return (CbcParameters::BACNotRun);
  }
  case 0: {
    return (CbcParameters::BACFinish);
  }
  case 1: {
    return (CbcParameters::BACStop);
  }
  case 2: {
    return (CbcParameters::BACAbandon);
  }
  case 5: {
    return (CbcParameters::BACUser);
  }
  default: { return (CbcParameters::BACInvalid); }
  }
}

CbcParameters::BACMinorStatus CbcParameters::translateMinor(int status)

{
  switch (status) {
  case -1: {
    return (CbcParameters::BACmInvalid);
  }
  case 0: {
    return (CbcParameters::BACmFinish);
  }
  case 1: {
    return (CbcParameters::BACmInfeas);
  }
  case 2: {
    return (CbcParameters::BACmGap);
  }
  case 3: {
    return (CbcParameters::BACmNodeLimit);
  }
  case 4: {
    return (CbcParameters::BACmTimeLimit);
  }
  case 5: {
    return (CbcParameters::BACmUser);
  }
  case 6: {
    return (CbcParameters::BACmSolnLimit);
  }
  case 7: {
    return (CbcParameters::BACmUbnd);
  }
  default: { return (CbcParameters::BACmOther); }
  }
}

/*
  A bit different --- given an OSI, use its interrogation functions to choose
  an appropriate BACMinorStatus code. Not everything matches up, eh?
*/
CbcParameters::BACMinorStatus
CbcParameters::translateMinor(const OsiSolverInterface *osi)

{
  if (osi->isProvenOptimal()) {
    return (CbcParameters::BACmFinish);
  } else if (osi->isProvenPrimalInfeasible()) {
    return (CbcParameters::BACmInfeas);
  } else if (osi->isProvenDualInfeasible()) {
    return (CbcParameters::BACmUbnd);
  } else {
    return (CbcParameters::BACmOther);
  }
}

/*
  A routine to set the bab_ status block given a CbcModel and an indication
  of where we're at in the search. Really, this is just a big mapping from
  CbcModel codes to CbcGeneric codes.
*/

void CbcParameters::setBaBStatus(CbcParameters::BACWhere where,
                                     bool haveAnswer,
                                     OsiSolverInterface *answerSolver)

{
  CbcParameters::BACMajorStatus major;
  CbcParameters::BACMinorStatus minor;

  major = translateMajor(model_->status());

  if (where == CbcParameters::BACwBareRoot ||
      where == CbcParameters::BACwIPPRelax) {
    minor = translateMinor(model_->solver());
  } else {
    minor = translateMinor(model_->secondaryStatus());
  }

  setBaBStatus(major, minor, where, haveAnswer, answerSolver);

  return;
}

/*
  Last, but not least, a routine to print the result.
*/

void CbcParameters::printBaBStatus()

{
  std::cout << "BAC result: stopped ";

  switch (bab_.where_) {
  case CbcParameters::BACwNotStarted: {
    std::cout << "before root relaxation";
    break;
  }
  case CbcParameters::BACwBareRoot: {
    std::cout << "after root relaxation";
    break;
  }
  case CbcParameters::BACwIPP: {
    std::cout << "after integer preprocessing";
    break;
  }
  case CbcParameters::BACwIPPRelax: {
    std::cout << "after solving preprocessed relaxation";
    break;
  }
  case CbcParameters::BACwBAC: {
    std::cout << "after branch-and-cut";
    break;
  }
  default: {
    std::cout << "!!invalid phase code!!";
    break;
  }
  }

  std::cout << std::endl << "    Branch-and-cut ";

  switch (bab_.majorStatus_) {
  case CbcParameters::BACNotRun: {
    std::cout << "never got started";
    break;
  }
  case CbcParameters::BACFinish: {
    std::cout << "finished";
    break;
  }
  case CbcParameters::BACStop: {
    std::cout << "stopped on a limit";
    break;
  }
  case CbcParameters::BACAbandon: {
    std::cout << "was abandoned";
    break;
  }
  case CbcParameters::BACUser: {
    std::cout << "stopped due to a user event";
    break;
  }
  default: {
    std::cout << "!!invalid major status code!!";
    break;
  }
  }

  std::cout << "; minor status is ";

  switch (bab_.minorStatus_) {
  case CbcParameters::BACmFinish: {
    std::cout << "optimal";
    break;
  }
  case CbcParameters::BACmInfeas: {
    std::cout << "infeasible";
    break;
  }
  case CbcParameters::BACmUbnd: {
    std::cout << "unbounded";
    break;
  }
  case CbcParameters::BACmGap: {
    std::cout << "reached specified integrality gap.";
    break;
  }
  case CbcParameters::BACmNodeLimit: {
    std::cout << "reached node limit";
    break;
  }
  case CbcParameters::BACmTimeLimit: {
    std::cout << "reached time limit";
    break;
  }
  case CbcParameters::BACmSolnLimit: {
    std::cout << "reached limit on number of solutions";
    break;
  }
  case CbcParameters::BACmUser: {
    std::cout << "stopped due to a user event";
    break;
  }
  case CbcParameters::BACmOther: {
    std::cout << "other";
    break;
  }
  default: {
    std::cout << "!!invalid minor status code!!";
    break;
  }
  }

  std::cout << "." << std::endl;
}

/* vi: softtabstop=2 shiftwidth=2 expandtab tabstop=2
 */
