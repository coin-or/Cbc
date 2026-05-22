/*
  Copyright (C) 2007, Lou Hafer, International Business Machines Corporation
  and others.  All Rights Reserved.

  This code is licensed under the terms of the Eclipse Public License (EPL).

*/

#include <cassert>

#include "ClpParamUtils.hpp"
#include "ClpParameters.hpp"
#include "ClpFactorization.hpp"

//###########################################################################
//###########################################################################

/*
  Constructor for settings class.
*/

ClpParameters::ClpParameters(bool cbcMode)
  : parameters_(ClpParam::LASTPARAM), model_(0)
{
   cbcMode_ = cbcMode;
   init(DefaultStrategy);
}

//###########################################################################
//###########################################################################

ClpParameters::ClpParameters(int strategy, bool cbcMode) :
  parameters_(ClpParam::LASTPARAM), model_(0) {

   cbcMode_ = cbcMode;
   init(strategy);

}

//###########################################################################
//###########################################################################

void ClpParameters::init(int strategy){

  for (int i = 0; i < parameters_.size(); i++){
     parameters_[i] = new ClpParam();
  }
  char dirsep = CoinFindDirSeparator();
  dfltDirectory_ = (dirsep == '/' ? "./" : ".\\");

  addClpParams();
  setDefaults(strategy);

}

//###########################################################################
//###########################################################################

ClpParameters::~ClpParameters() {
  for (int i = 0; i < parameters_.size(); i++){
     delete parameters_[i];
  }
}

ClpParameters::ClpParameters(const ClpParameters &rhs)
  : parameters_(rhs.parameters_.size()), model_(rhs.model_)
{
  cbcMode_ = rhs.cbcMode_;
  for (int i = 0; i < (int)rhs.parameters_.size(); i++)
    parameters_[i] = rhs.parameters_[i] ? rhs.parameters_[i]->clone() : nullptr;
  dfltDirectory_ = rhs.dfltDirectory_;
}

ClpParameters &ClpParameters::operator=(const ClpParameters &rhs)
{
  if (this == &rhs)
    return *this;
  for (int i = 0; i < (int)parameters_.size(); i++)
    delete parameters_[i];
  parameters_.resize(rhs.parameters_.size());
  for (int i = 0; i < (int)rhs.parameters_.size(); i++)
    parameters_[i] = rhs.parameters_[i] ? rhs.parameters_[i]->clone() : nullptr;
  dfltDirectory_ = rhs.dfltDirectory_;
  model_ = rhs.model_;
  cbcMode_ = rhs.cbcMode_;
  return *this;
}

//###########################################################################
//###########################################################################

int ClpParameters::matches(std::string field, int &numberMatches) {
  int firstMatch = -1;
  for (int iParam = 0; iParam < (int)parameters_.size(); iParam++) {
    int match = parameters_[iParam]->matches(field);
    if (match == 1) {
      numberMatches = 1;
      return iParam;
    } else {
      if (match) {
        if (firstMatch < 0) {
          firstMatch = iParam;
        }
        numberMatches++;
      }
    }
  }
  return firstMatch < 0 ? ClpParam::INVALID : firstMatch;
}

//###########################################################################
//###########################################################################

void ClpParameters::addClpParams() {

  addClpStrParams();
  addClpDirParams();
  addClpFileParams();
  addClpHelpParams();
  addClpActionParams();
  addClpKwdParams();
  addClpDblParams();
  addClpIntParams();
  addClpBoolParams();

  for (int code = ClpParam::FIRSTPARAM + 1; code < ClpParam::LASTPARAM;
       code++) {
     getParam(code)->setParameters(this);
     getParam(code)->setParamCode(code);
  }

  // --- Assign semantic topics to every Clp parameter ---

  // Bulk-set: file params → I/O
  for (int i = ClpParam::FIRSTFILEPARAM + 1; i < ClpParam::LASTFILEPARAM; i++)
    parameters_[i]->setTopic("I/O");

  // Bulk-set: directory params → I/O
  for (int i = ClpParam::FIRSTDIRECTORYPARAM + 1; i < ClpParam::LASTDIRECTORYPARAM; i++)
    parameters_[i]->setTopic("I/O");

  // Action params — Solving
  for (int code : {ClpParam::DUALSIMPLEX, ClpParam::PRIMALSIMPLEX,
                    ClpParam::EITHERSIMPLEX, ClpParam::BARRIER,
                    ClpParam::SOLVE, ClpParam::PARAMETRICS,
                    ClpParam::NETWORK, ClpParam::PLUSMINUS,
                    ClpParam::ALLSLACK, ClpParam::TIGHTEN,
                    ClpParam::REALLY_SCALE})
    parameters_[code]->setTopic("Solving");

  // Action params — I/O
  for (int code : {ClpParam::IMPORT, ClpParam::EXPORT,
                    ClpParam::BASISIN, ClpParam::BASISOUT,
                    ClpParam::READMODEL, ClpParam::READMODEL_OLD,
                    ClpParam::WRITEMODEL, ClpParam::WRITEMODEL_OLD,
                    ClpParam::PRINTSOL, ClpParam::WRITESOL,
                    ClpParam::WRITESOL_OLD, ClpParam::WRITESOLBINARY,
                    ClpParam::WRITESOLBINARY_OLD,
                    ClpParam::WRITEGMPLSOL, ClpParam::WRITEGMPLSOL_OLD})
    parameters_[code]->setTopic("I/O");

  // Action params — other
  parameters_[ClpParam::STATISTICS]->setTopic("Output");
  parameters_[ClpParam::REVERSE]->setTopic("Solving");
  parameters_[ClpParam::MINIMIZE]->setTopic("Solving");
  parameters_[ClpParam::MAXIMIZE]->setTopic("Solving");
  parameters_[ClpParam::OUTDUPROWS]->setTopic("LP Presolve");

  // Keyword params — Simplex
  for (int code : {ClpParam::DUALPIVOT, ClpParam::PRIMALPIVOT,
                    ClpParam::CRASH, ClpParam::FACTORIZATION})
    parameters_[code]->setTopic("Simplex");

  // Keyword params — Barrier
  for (int code : {ClpParam::CHOLESKY, ClpParam::GAMMA,
                    ClpParam::BARRIERSCALE, ClpParam::CROSSOVER})
    parameters_[code]->setTopic("Barrier");

  // Keyword params — other
  parameters_[ClpParam::SCALING]->setTopic("Scaling");
  parameters_[ClpParam::PRESOLVE]->setTopic("LP Presolve");
  parameters_[ClpParam::BIASLU]->setTopic("Simplex");
  parameters_[ClpParam::DIRECTION]->setTopic("Solving");
  parameters_[ClpParam::INTPRINT]->setTopic("Output");
  parameters_[ClpParam::COMMANDPRINTLEVEL]->setTopic("Output");
  parameters_[ClpParam::VECTOR]->setTopic("Solving");

  // Double params — Tolerances
  for (int code : {ClpParam::DUALTOLERANCE, ClpParam::PRIMALTOLERANCE,
                    ClpParam::ZEROTOLERANCE})
    parameters_[code]->setTopic("Tolerances");
  parameters_[ClpParam::PRESOLVETOLERANCE]->setTopic("LP Presolve");

  // Double params — Simplex
  for (int code : {ClpParam::DUALBOUND, ClpParam::PRIMALWEIGHT,
                    ClpParam::PSI})
    parameters_[code]->setTopic("Simplex");

  // Double params — Scaling
  for (int code : {ClpParam::OBJSCALE, ClpParam::OBJSCALE2,
                    ClpParam::RHSSCALE})
    parameters_[code]->setTopic("Scaling");

  // Double params — other
  parameters_[ClpParam::TIMELIMIT]->setTopic("Stopping");
  parameters_[ClpParam::PROGRESS]->setTopic("Output");
  parameters_[ClpParam::FAKEBOUND]->setTopic("Strategy");

  // Integer params — Simplex
  for (int code : {ClpParam::MAXFACTOR, ClpParam::MAXITERATION,
                    ClpParam::SPRINT, ClpParam::IDIOT,
                    ClpParam::SLPVALUE, ClpParam::PERTVALUE,
                    ClpParam::SPECIALOPTIONS, ClpParam::MORESPECIALOPTIONS,
                    ClpParam::DENSE,
                    ClpParam::SMALLFACT})
    parameters_[code]->setTopic("Simplex");
  parameters_[ClpParam::SUBSTITUTION]->setTopic("LP Presolve");

  // Integer params — LP Presolve
  parameters_[ClpParam::PRESOLVEPASS]->setTopic("LP Presolve");

  // Integer params — Solving (model transformations)
  for (int code : {ClpParam::DUALIZE, ClpParam::DECOMPOSE_BLOCKS})
    parameters_[code]->setTopic("Solving");
  parameters_[ClpParam::CPP]->setTopic("Output");

  // Integer params — Output
  for (int code : {ClpParam::LOGLEVEL, ClpParam::OUTPUTFORMAT,
                    ClpParam::PRINTOPTIONS, ClpParam::VERBOSE})
    parameters_[code]->setTopic("Output");

  // Integer params — other
  parameters_[ClpParam::RANDOMSEED]->setTopic("Solving");
  parameters_[ClpParam::THREADS]->setTopic("Parallelism");

  // Bool params
  for (int code : {ClpParam::AUTOSCALE, ClpParam::SPARSEFACTOR,
                    ClpParam::PFI, ClpParam::KKT})
    parameters_[code]->setTopic("Simplex");
  parameters_[ClpParam::PERTURBATION]->setTopic("Simplex");
  parameters_[ClpParam::KEEPNAMES]->setTopic("I/O");
  parameters_[ClpParam::ERRORSALLOWED]->setTopic("I/O");
  parameters_[ClpParam::MESSAGES]->setTopic("Output");
  parameters_[ClpParam::BUFFER_MODE]->setTopic("Output");

  // String params
  parameters_[ClpParam::PRINTMASK]->setTopic("Output");

  return;
}

//###########################################################################
//###########################################################################

void ClpParameters::setDefaults(int strategy) {

   if (!cbcMode_){
      // First, set up the psrameters that have no strategy implications
      for (int code = ClpParam::FIRSTDIRECTORYPARAM + 1;
           code < ClpParam::LASTDIRECTORYPARAM; code++) {
         getParam(code)->setDefault(dfltDirectory_);
      }
   } else {
     // change name of time limit from seconds to lpseconds
     parameters_[ClpParam::TIMELIMIT]->setName("lpsec!onds");
   }

   parameters_[ClpParam::BASISFILE]->setDefault(std::string("default.bas"));
   parameters_[ClpParam::PARAMETRICSFILE]->setDefault(std::string("default.mps"));

   if (!cbcMode_){
      parameters_[ClpParam::EXPORTFILE]->setDefault(std::string("default.mps"));
      parameters_[ClpParam::GMPLSOLFILE]->setDefault(std::string("gmpl.sol"));
      parameters_[ClpParam::IMPORTFILE]->setDefault(std::string("default.mps"));
      parameters_[ClpParam::PRINTMASK]->setDefault("");
      parameters_[ClpParam::MODELFILE]->setDefault(std::string("default.prob"));
      parameters_[ClpParam::SOLUTIONFILE]->setDefault(std::string("solution.sln"));
      parameters_[ClpParam::SOLUTIONBINARYFILE]->setDefault(std::string("solution.file"));
   }

   // Now set up  parameters according to overall strategies
   switch (strategy) {
    case ClpParameters::DefaultStrategy:
      parameters_[ClpParam::COMMANDPRINTLEVEL]->setDefault("more");
      parameters_[ClpParam::BARRIERSCALE]->setDefault("off");
      parameters_[ClpParam::BIASLU]->setDefault("LX");
      parameters_[ClpParam::CHOLESKY]->setDefault("native");
      parameters_[ClpParam::CRASH]->setDefault("off");
      parameters_[ClpParam::CROSSOVER]->setDefault("on");
      parameters_[ClpParam::DIRECTION]->setDefault("min!imize");
      parameters_[ClpParam::DUALPIVOT]->setDefault("auto!matic");
      parameters_[ClpParam::FACTORIZATION]->setDefault("normal");
      parameters_[ClpParam::GAMMA]->setDefault("off");
      parameters_[ClpParam::PRESOLVE]->setDefault("on");
      parameters_[ClpParam::PRIMALPIVOT]->setDefault("auto!matic");
      parameters_[ClpParam::INTPRINT]->setDefault("normal");
      parameters_[ClpParam::SCALING]->setDefault("auto!matic");
#ifndef COIN_AVX2
      parameters_[ClpParam::VECTOR]->setDefault("off");
#else
      parameters_[ClpParam::VECTOR]->setDefault("off");
#endif
      parameters_[ClpParam::DUALBOUND]->setDefault(0.0);
      parameters_[ClpParam::FAKEBOUND]->setDefault(0.0);
      parameters_[ClpParam::FAKEBOUND]->setDefault(0.0);
      parameters_[ClpParam::OBJSCALE]->setDefault(1.0);
      parameters_[ClpParam::PRESOLVETOLERANCE]->setDefault(0.09);
      parameters_[ClpParam::PRESOLVETOLERANCE]->setDefault(0.09);
      parameters_[ClpParam::PRIMALTOLERANCE]->setDefault(0.0);
      parameters_[ClpParam::PRIMALWEIGHT]->setDefault(0.0);
      parameters_[ClpParam::PSI]->setDefault(-0.5);
      parameters_[ClpParam::PROGRESS]->setDefault(0.7);
      parameters_[ClpParam::OBJSCALE2]->setDefault(1.0);
      parameters_[ClpParam::RHSSCALE]->setDefault(1.0);
      parameters_[ClpParam::TIMELIMIT]->setDefault(-1.0);
      parameters_[ClpParam::ZEROTOLERANCE]->setDefault(1.0e-20);
      parameters_[ClpParam::CPP]->setDefault(0);
      parameters_[ClpParam::DECOMPOSE_BLOCKS]->setDefault(0);
      parameters_[ClpParam::DENSE]->setDefault(-1);
      parameters_[ClpParam::DUALIZE]->setDefault(0);
      parameters_[ClpParam::IDIOT]->setDefault(0);
      parameters_[ClpParam::MAXFACTOR]->setDefault(0);
      parameters_[ClpParam::MAXITERATION]->setDefault(0);
      parameters_[ClpParam::MORESPECIALOPTIONS]->setDefault(0);
      parameters_[ClpParam::PRESOLVEPASS]->setDefault(0);
      parameters_[ClpParam::PERTVALUE]->setDefault(0);
      parameters_[ClpParam::RANDOMSEED]->setDefault(1234567);
      parameters_[ClpParam::SLPVALUE]->setDefault(0);
      parameters_[ClpParam::SMALLFACT]->setDefault(-1);
      parameters_[ClpParam::SPECIALOPTIONS]->setDefault(0);
      parameters_[ClpParam::SPRINT]->setDefault(0);
      parameters_[ClpParam::SUBSTITUTION]->setDefault(3);
#ifdef CLP_THREAD
      parameters_[ClpParam::SUBSTITUTION]->setDefault(3);
#endif
      parameters_[ClpParam::AUTOSCALE]->setDefault("off");
      parameters_[ClpParam::BUFFER_MODE]->setDefault("on");
      parameters_[ClpParam::ERRORSALLOWED]->setDefault("off");
      parameters_[ClpParam::KEEPNAMES]->setDefault("on");
      parameters_[ClpParam::KKT]->setDefault("off");
      parameters_[ClpParam::MESSAGES]->setDefault("off");
      parameters_[ClpParam::PERTURBATION]->setDefault("on");
      parameters_[ClpParam::PFI]->setDefault("off");
      parameters_[ClpParam::SPARSEFACTOR]->setDefault("on");
      if (!cbcMode_){
         parameters_[ClpParam::LOGLEVEL]->setDefault(1);
         parameters_[ClpParam::OUTPUTFORMAT]->setDefault(0);
         parameters_[ClpParam::PRINTOPTIONS]->setDefault(0);
         parameters_[ClpParam::PROGRESSITER]->setDefault(0);
         parameters_[ClpParam::VERBOSE]->setDefault(0);
      }
      break;

    default:
      std::cout << "Unknown strategy!" << std::endl;
      break;
   }

}

//###########################################################################
//###########################################################################

#ifdef CBC_CLUMSY_CODING
//#define PRINT_CLP_CHANGES
// Synchronize Clp model - Int and Dbl
void ClpParameters::synchronizeModel() {
  if (model_) {
    // Integer parameters
    int intValue;
    int modelIntValue;
    parameters_[ClpParam::MAXFACTOR]->getVal(intValue);
#ifdef PRINT_CLP_CHANGES
    modelIntValue = model_->factorization()->maximumPivots();
    if (intValue!=modelIntValue)
      printf("changing ? from %d to %d at line %d\n",modelIntValue,intValue,__LINE__+1);
#endif
    model_->factorization()->maximumPivots(intValue);
    parameters_[ClpParam::PERTVALUE]->getVal(intValue);
#ifdef PRINT_CLP_CHANGES
    modelIntValue = model_->perturbation();
    if (intValue!=modelIntValue)
      printf("changing ? from %d to %d at line %d\n",modelIntValue,intValue,__LINE__+1);
#endif
    model_->setPerturbation(intValue);
    parameters_[ClpParam::MAXITERATION]->getVal(intValue);
#ifdef PRINT_CLP_CHANGES
    modelIntValue = model_->maximumIterations();
    if (intValue!=modelIntValue)
      printf("changing ? from %d to %d at line %d\n",modelIntValue,intValue,__LINE__+1);
#endif
    model_->setMaximumIterations(intValue);
    parameters_[ClpParam::SPECIALOPTIONS]->getVal(intValue);
#ifdef PRINT_CLP_CHANGES
    modelIntValue = model_->specialOptions();
    if (intValue!=modelIntValue)
      printf("changing ? from %d to %d at line %d\n",modelIntValue,intValue,__LINE__+1);
#endif
    model_->setSpecialOptions(intValue);
    parameters_[ClpParam::RANDOMSEED]->getVal(intValue);
    model_->setRandomSeed(intValue);
    parameters_[ClpParam::MORESPECIALOPTIONS]->getVal(intValue);
#ifdef PRINT_CLP_CHANGES
    modelIntValue = model_->moreSpecialOptions();
    if (intValue!=modelIntValue)
      printf("changing ? from %d to %d at line %d\n",modelIntValue,intValue,__LINE__+1);
#endif
    model_->setMoreSpecialOptions(intValue);
    parameters_[ClpParam::LOGLEVEL]->getVal(intValue);
#ifdef PRINT_CLP_CHANGES
    modelIntValue = model_->logLevel();
    if (intValue!=modelIntValue)
      printf("changing ? from %d to %d at line %d\n",modelIntValue,intValue,__LINE__+1);
#endif
    model_->setLogLevel(intValue);
    if (intValue > 2)
      model_->factorization()->messageLevel(8);
    else
      model_->factorization()->messageLevel(0);
    // Double parameters
    double doubleValue;
    double modelDoubleValue;
    parameters_[ClpParam::DUALTOLERANCE]->getVal(doubleValue);
#ifdef PRINT_CLP_CHANGES
    modelDoubleValue = model_->dualTolerance();
    if (doubleValue!=modelDoubleValue)
      printf("changing ? from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
    model_->setDualTolerance(doubleValue);
    parameters_[ClpParam::PRIMALTOLERANCE]->getVal(doubleValue);
#ifdef PRINT_CLP_CHANGES
    modelDoubleValue = model_->primalTolerance();
    if (doubleValue!=modelDoubleValue)
      printf("changing ? from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
    model_->setPrimalTolerance(doubleValue);
    parameters_[ClpParam::ZEROTOLERANCE]->getVal(doubleValue);
#ifdef PRINT_CLP_CHANGES
    modelDoubleValue = model_->getSmallElementValue();
    if (doubleValue!=modelDoubleValue)
      printf("changing ? from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
    model_->setSmallElementValue(doubleValue);
    parameters_[ClpParam::DUALBOUND]->getVal(doubleValue);
#ifdef PRINT_CLP_CHANGES
    modelDoubleValue = model_->dualBound();
    if (doubleValue!=modelDoubleValue)
      printf("changing ? from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
    model_->setDualBound(doubleValue);
    parameters_[ClpParam::PRIMALWEIGHT]->getVal(doubleValue);
#ifdef PRINT_CLP_CHANGES
    modelDoubleValue = model_->infeasibilityCost();
    if (doubleValue!=modelDoubleValue)
      printf("changing ? from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
    model_->setInfeasibilityCost(doubleValue);
    parameters_[ClpParam::TIMELIMIT]->getVal(doubleValue);
#ifdef PRINT_CLP_CHANGES
    modelDoubleValue = model_->maximumSeconds();
    if (doubleValue!=modelDoubleValue)
      printf("changing ? from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
    model_->setMaximumSeconds(doubleValue);
    parameters_[ClpParam::OBJSCALE]->getVal(doubleValue);
#ifdef PRINT_CLP_CHANGES
    modelDoubleValue = model_->objectiveScale();
    if (doubleValue!=modelDoubleValue)
      printf("changing ? from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
    model_->setObjectiveScale(doubleValue);
    parameters_[ClpParam::RHSSCALE]->getVal(doubleValue);
#ifdef PRINT_CLP_CHANGES
    modelDoubleValue = model_->rhsScale();
    if (doubleValue!=modelDoubleValue)
      printf("changing ? from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
    model_->setRhsScale(doubleValue);
    parameters_[ClpParam::PRESOLVETOLERANCE]->getVal(doubleValue);
#ifdef PRINT_CLP_CHANGES
    model_->getDblParam(ClpPresolveTolerance, doubleValue);
    if (doubleValue!=modelDoubleValue)
      printf("changing ? from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
    model_->setDblParam(ClpPresolveTolerance, doubleValue);
    parameters_[ClpParam::PROGRESS]->getVal(doubleValue);
#ifdef PRINT_CLP_CHANGES
    modelDoubleValue = model_->getMinIntervalProgressUpdate();
    if (doubleValue!=modelDoubleValue)
      printf("changing ? from %g to %g at line %d\n",modelDoubleValue,doubleValue,__LINE__+1);
#endif
#if CLP_OLD_PROGRESS
    // I dislike new way
    if (doubleValue < 0.0)
#endif
    model_->setMinIntervalProgressUpdate(doubleValue);
  }
}
#endif

//###########################################################################
//###########################################################################

void ClpParameters::addClpHelpParams() {

  // All these parameters duplicate ones in Cbc
  if (cbcMode_){
     return;
  }

  for (int code = ClpParam::FIRSTHELPPARAM + 1; code < ClpParam::LASTHELPPARAM;
       code++) {
    getParam(code)->setPushFunc(ClpParamUtils::doHelpParam);
    getParam(code)->setType(CoinParam::paramAct);
  }

  parameters_[ClpParam::GENERALQUERY]->setup(
      "?", "Print a list of commands", "", CoinParam::displayPriorityNone);

  parameters_[ClpParam::FULLGENERALQUERY]->setup(
      "???", "Print a list with *all* commands, even those hidden with `?'", "",
      CoinParam::displayPriorityNone);
}

//###########################################################################
//###########################################################################

void ClpParameters::addClpActionParams() {

  // First, we add the parameters that are unique to Cbc
  for (int code = ClpParam::FIRSTACTIONPARAM + 1;
       code < ClpParam::LASTCLPACTIONPARAM; code++) {
    getParam(code)->setType(CoinParam::paramAct);
  }

  parameters_[ClpParam::ALLSLACK]->setup(
      "allS!lack", "Set basis back to all slack and reset solution",
      "Mainly useful for tuning purposes.  Normally the first dual or primal "
      "will be using an all slack basis anyway.",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::BARRIER]->setup(
      "barr!ier", "Solve using primal dual predictor corrector algorithm",
      "This command solves the current model using the  primal dual predictor "
      "corrector algorithm. You may want to link in an alternative ordering "
      "and factorization. It will also solve models with quadratic objectives.",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::BASISIN]->setup(
      "basisI!n", "Import basis from bas file",
      "This will read an MPS format basis file from the given file name.  It "
      "will use the default directory given by 'directory'.  A name of '$' "
      "will use the previous value for the name. This is initialized to '', "
      "i.e. it must be set.  If you have libz then it can read compressed "
      "files 'xxxxxxxx.gz' or xxxxxxxx.bz2.",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::BASISOUT]->setup(
      "basisO!ut", "Export basis as bas file",
      "This will write an MPS format basis file to the given file name.  It "
      "will use the default directory given by 'directory'.  A name of '$' "
      "will use the previous value for the name.  This is initialized to "
      "'default.bas'.",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::DUALSIMPLEX]->setup(
      "dualS!implex", "Do dual simplex algorithm",
      "This command solves the continuous relaxation of the current model "
      "using the dual steepest edge algorithm. The time and iterations may be "
      "affected by settings such as presolve, scaling, crash and also by dual "
      "pivot method, fake bound on variables and dual and primal tolerances.",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::EITHERSIMPLEX]->setup(
      "either!Simplex", "Do dual or primal simplex algorithm",
      "This command solves the continuous relaxation of the current model "
      "using the dual or primal algorithm, based on a dubious analysis of "
      "model.",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::GUESS]->setup(
      "guess", "[DEPRECATED] Guesses at good parameters",
      "Deprecated: use '-lpMethod=auto' instead, which uses a Random Forest "
      "model trained on 207 instance features to automatically select the "
      "best LP algorithm and settings.\n"
      "This action looks at basic model statistics and heuristically sets a "
      "few LP parameters, but is superseded by the ML-based auto selection.",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::NETLIB_EITHER]->setup(
      "netlib", "Solve entire netlib test set",
      "This exercises the unit test for clp and then solves the netlib test "
      "set using dual or primal. The user can set options before e.g. clp "
      "-presolve off -netlib",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::NETLIB_BARRIER]->setup(
      "netlibB!arrier", "Solve entire netlib test set with barrier",
      "This exercises the unit test for clp and then solves the netlib test "
      "set using barrier. The user can set options before e.g. clp -kkt on "
      "-netlib",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::NETLIB_DUAL]->setup(
      "netlibD!ual", "Solve entire netlib test set (dual)",
      "This exercises the unit test for clp and then solves the netlib test "
      "set using dual. The user can set options before e.g. clp -presolve off "
      "-netlib",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::NETLIB_PRIMAL]->setup(
      "netlibP!rimal", "Solve entire netlib test set (primal)",
      "This exercises the unit test for clp and then solves the netlib test "
      "set using primal. The user can set options before e.g. clp -presolve "
      "off -netlibp",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::NETLIB_TUNE]->setup(
      "netlibT!une", "Solve entire netlib test set with 'best' algorithm",
      "This exercises the unit test for clp and then solves the netlib test "
      "set using whatever works best. I know this is cheating but it also "
      "stresses the code better by doing a mixture of stuff. The best "
      "algorithm was chosen on a Linux ThinkPad using native cholesky with "
      "University of Florida ordering.",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::NETWORK]->setup(
      "network", "Tries to make network matrix",
      "Clp will go faster if the matrix can be converted to a network.  The "
      "matrix operations may be a bit faster with more efficient storage, but "
      "the main advantage comes from using a network factorization. It will "
      "probably not be as fast as a specialized network code.",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::PARAMETRICS]->setup(
      "para!metrics", "Import data from file and do parametrics",
      "This will read a file with parametric data from the given file name and "
      "then do parametrics. It will use the default directory given by "
      "'directory'. A name of '$' will use the previous value for the name. "
      "This is initialized to '', i.e. it must be set.  This can not read from "
      "compressed files. File is in modified csv format - a line ROWS will be "
      "followed by rows data while a line COLUMNS will be followed by column "
      "data.  The last line should be ENDATA. The ROWS line must exist and is "
      "in the format ROWS, inital theta, final theta, interval theta, n where "
      "n is 0 to get CLPI0062 message at interval or at each change of theta "
      "and 1 to get CLPI0063 message at each iteration.  If interval theta is "
      "0.0 or >= final theta then no interval reporting.  n may be missed out "
      "when it is taken as 0.  If there is Row data then there is a headings "
      "line with allowed headings - name, number, lower(rhs change), upper(rhs "
      "change), rhs(change).  Either the lower and upper fields should be "
      "given or the rhs field. The optional COLUMNS line is followed by a "
      "headings line with allowed headings - name, number, objective(change), "
      "lower(change), upper(change). Exactly one of name and number must be "
      "given for either section and missing ones have value 0.0.",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::PLUSMINUS]->setup(
      "plus!Minus", "Tries to make +- 1 matrix",
      "Clp will go slightly faster if the matrix can be converted so that the "
      "elements are not stored and are known to be unit.  The main advantage "
      "is memory use.  Clp may automatically see if it can convert the problem "
      "so you should not need to use this.",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::PRIMALSIMPLEX]->setup(
      "primalS!implex", "Do primal simplex algorithm",
      "This command solves the continuous relaxation of the current model "
      "using the primal algorithm. The default is to use exact devex. The time "
      "and iterations may be affected by settings such as presolve, scaling, "
      "crash and also by column selection  method, infeasibility weight and "
      "dual and primal tolerances.",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::REALLY_SCALE]->setup("reallyS!cale",
                                            "Scales model in place", "",
                                            CoinParam::displayPriorityLow);

  parameters_[ClpParam::REVERSE]->setup("reverse",
      "Reverses sign of objective",
      "Useful for testing if maximization works correctly",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::TIGHTEN]->setup("tightLP",
                                       "Poor person's preSolve for now", "",
                                       CoinParam::displayPriorityNone);

  parameters_[ClpParam::USERCLP]->setup(
      "userClp", "Hand coded Clp stuff",
      "There are times e.g. when using AMPL interface when you may wish to do "
      "something unusual. Look for USERCLP in main driver and modify sample "
      "code.",
      CoinParam::displayPriorityNone);

  if (cbcMode_){
     // Remaining parameters are duplicates and don't get added in Cbc mode.
     return;
  }

  // Add remaining parameters now
  for (int code = ClpParam::LASTCLPACTIONPARAM + 1;
       code < ClpParam::LASTACTIONPARAM; code++) {
    getParam(code)->setType(CoinParam::paramAct);
  }

  parameters_[ClpParam::DUMMY]->setup(
      "sleep", "for debug",
      "If passed to solver fom ampl, then ampl will wait so that you can copy "
      ".nl file for debug.",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::END]->setup(
      "end", "Stops clp execution",
      "This stops execution ; end, exit, quit and stop are synonyms",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::QUIT]->setup("quit", "Stops clp execution",
      "This stops the execution of Clp, end, exit, quit and stop are synonyms",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::EXIT]->setup("exit", "Stops clp execution",
      "This stops the execution of Clp, end, exit, quit and stop are synonyms",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::STOP]->setup("stop", "Stops clp execution",
      "This stops the execution of Clp, end, exit, quit and stop are synonyms",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::ENVIRONMENT]->setup(
      "environ!ment", "Read commands from environment",
      "This starts reading from environment variable CLP_ENVIRONMENT.",
      CoinParam::displayPriorityNone);

  parameters_[ClpParam::EXPORT]->setup(
      "export", "Export model as mps file",
      "This will write an MPS format file to the given file name.  It will use "
      "the default directory given by 'directory'.  A name of '$' will use the "
      "previous value for the name.  This is initialized to 'default.mps'. It "
      "can be useful to get rid of the original names and go over to using "
      "Rnnnnnnn and Cnnnnnnn.  This can be done by setting 'keepnames' off "
      "before importing mps file.",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::IMPORT]->setup(
      "import", "Import model from mps file",
      "This will read an MPS format file from the given file name.  It will "
      "use the default directory given by 'directory'.  A name of '$' will use "
      "the previous value for the name. This is initialized to '', i.e. it "
      "must be set.  If you have libgz then it can read compressed files "
      "'xxxxxxxx.gz' or 'xxxxxxxx.bz2'. If 'keepnames' is off, then names are "
      "dropped -> Rnnnnnnn and Cnnnnnnn.",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::MINIMIZE]->setup(
      "min!imize", "Set optimization direction to minimize",
      "The default is minimize - use 'maximize' for maximization.\n This "
      "should only be necessary if you have previously set maximization You "
      "can also use the parameters 'direction minimize'.",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::MAXIMIZE]->setup(
      "max!imize", "Set optimization direction to maximize",
      "The default is minimize - use 'maximize' for maximization.\n"
      "You can also use the parameters 'direction maximize'.",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::OUTDUPROWS]->setup(
      "outDup!licates", "takes duplicate rows etc out of integer model", "",
      CoinParam::displayPriorityNone);

  parameters_[ClpParam::PRINTSOL]->setup(
      "printS!olution", "prints solution to stdout",
      "This will write a binary solution file to the file set by solFile.",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::PRINTVERSION]->setup(
      "version", "Print version", "", CoinParam::displayPriorityHigh);

  parameters_[ClpParam::READMODEL]->setup(
      "readM!odel", "Reads problem from a binary save file",
      "This will read the problem saved by 'writeModel' from the file name "
      "set by 'modelFile'.",
      CoinParam::displayPriorityHigh);

  // For backward compatibility
  parameters_[ClpParam::READMODEL_OLD]->setup(
      "restoreM!odel", "Reads problem from a binary save file (synonym for "
      "readModel)",
      "This will read the problem saved by 'writeModel' from the file name "
      "set by 'modelFile'.",
      CoinParam::displayPriorityNone);

  parameters_[ClpParam::SOLVE]->setup(
      "solv!e", "Do dual or primal simplex algorithm",
      "This command solves the continuous relaxation of the current model "
      "using the dual or primal algorithm, based on a dubious analysis of "
      "model.",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::STATISTICS]->setup(
      "stat!istics", "Print some statistics",
      "This command prints some statistics for the current model. If log level "
      ">1 then more is printed. These are for presolved model if presolve on "
      "(and unscaled).",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::UNITTEST]->setup("unitTest", "Do unit test",
                                         "This exercises the unit test for clp",
                                         CoinParam::displayPriorityLow);

  parameters_[ClpParam::WRITEGMPLSOL]->setup(
      "writeGSolu!tion", "Puts glpk solution to file",
      "Will write a glpk solution file to the given file name.  It will use "
      "the default directory given by 'directory'.  A name of '$' will use the "
      "previous value for the name.  This is initialized to 'stdout' (this "
      "defaults to ordinary solution if stdout). If problem created from gmpl "
      "model - will do any reports.",
      CoinParam::displayPriorityHigh);

  // For backward compatibility
  parameters_[ClpParam::WRITEGMPLSOL_OLD]->setup(
      "gsolu!tion", "Puts glpk solution to file (synonym for writeGSolu!tion)",
      "Will write a glpk solution file to the given file name.  It will use "
      "the default directory given by 'directory'.  A name of '$' will use the "
      "previous value for the name.  This is initialized to 'stdout' (this "
      "defaults to ordinary solution if stdout). If problem created from gmpl "
      "model - will do any reports.",
      CoinParam::displayPriorityNone);

  parameters_[ClpParam::WRITEMODEL]->setup(
      "writeM!odel", "save model to binary file",
      "This will write the problem in binary foramt to the file name set by "
      "'modelFile' for future use by readModel.",
      CoinParam::displayPriorityHigh);

  // For backward compatibility
  parameters_[ClpParam::WRITEMODEL_OLD]->setup(
      "saveM!odel", "save model to binary file (synonym for writeModel)",
      "This will write the problem in binary foramt to the file name set by "
      "'modelFile' for future use by readModel.",
      CoinParam::displayPriorityNone);

  parameters_[ClpParam::WRITESOL]->setup(
      "writeS!olution", "writes solution to file (or stdout)",
      "This will write a primitive solution file to the file set by "
      "'solFile'. The amount of output can be varied using "
      "'printingOptions' or 'printMask'.",
      CoinParam::displayPriorityHigh);

  // For backward compatibility
  parameters_[ClpParam::WRITESOL_OLD]->setup(
      "solu!tion", "writes solution to file (or stdout) (synonym for "
      "writeSolution).",
      "This will write a primitive solution file to the file set by "
      "'solFile'. The amount of output can be varied using "
      "'printingOptions' or 'printMask'.",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::WRITESOLBINARY]->setup(
      "writeSolB!inary", "writes solution to file in binary format",
      "This will write a binary solution file to the file set by "
      "'solBinaryFile'. To read the file use fread(int) twice to pick up "
      "number of rows and columns, then fread(double) to pick up objective "
      "value, then pick up row activities, row duals, column activities and "
      "reduced costs - see bottom of ClpParamUtils.cpp for code that reads or "
      "writes file. If name contains '_fix_read_', then does not write but "
      "reads and will fix all variables",
      CoinParam::displayPriorityHigh);

  // For backward compatibility
  parameters_[ClpParam::WRITESOLBINARY_OLD]->setup(
      "saveS!olution", "writes solution to file in binary format (synonym for "
      "writeSolBinary",
      "This will write a binary solution file to the file set by "
      "'solBinaryFile'. To read the file use fread(int) twice to pick up "
      "number of rows and columns, then fread(double) to pick up objective "
      "value, then pick up row activities, row duals, column activities and "
      "reduced costs - see bottom of ClpParamUtils.cpp for code that reads or "
      "writes file. If name contains '_fix_read_', then does not write but "
      "reads and will fix all variables",
      CoinParam::displayPriorityNone);
}

//###########################################################################
//###########################################################################

void ClpParameters::addClpFileParams() {

  for (int code = ClpParam::FIRSTFILEPARAM + 1;
       code < ClpParam::LASTCLPFILEPARAM; code++) {
    getParam(code)->setType(CoinParam::paramFile);
  }

  parameters_[ClpParam::BASISFILE]->setup(
      "basisF!ile", "sets the name for file for reading/writing the basis",
      "This will read an MPS format basis file from the given file name.  It "
      "will use the default directory given by 'directory'.  If no name is "
      "specified, the previous value will be used. This is initialized to '', "
      "i.e. it must be set.  If you have libz then it can read compressed "
      "files 'xxxxxxxx.gz' or xxxxxxxx.bz2.",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::PARAMETRICSFILE]->setup(
      "paramF!ile", "set name of file to import parametrics data from",
      "This will read a file with parametric data from the given file name and "
      "then do parametrics. It will use the default directory given by "
      "'directory'. A name of '$' will use the previous value for the name. "
      "This is initialized to '', i.e. it must be set.  This can not read from "
      "compressed files. File is in modified csv format - a line ROWS will be "
      "followed by rows data while a line COLUMNS will be followed by column "
      "data.  The last line should be ENDATA. The ROWS line must exist and is "
      "in the format ROWS, inital theta, final theta, interval theta, n where "
      "n is 0 to get CLPI0062 message at interval or at each change of theta "
      "and 1 to get CLPI0063 message at each iteration.  If interval theta is "
      "0.0 or >= final theta then no interval reporting.  n may be missed out "
      "when it is taken as 0.  If there is Row data then there is a headings "
      "line with allowed headings - name, number, lower(rhs change), upper(rhs "
      "change), rhs(change).  Either the lower and upper fields should be "
      "given or the rhs field. The optional COLUMNS line is followed by a "
      "headings line with allowed headings - name, number, objective(change), "
      "lower(change), upper(change). Exactly one of name and number must be "
      "given for either section and missing ones have value 0.0.",
      CoinParam::displayPriorityHigh);

  // Remaining parameters are duplicates
  if (cbcMode_){
     return;
  }

  for (int code = ClpParam::LASTCLPFILEPARAM + 1;
       code < ClpParam::LASTFILEPARAM; code++) {
    getParam(code)->setType(CoinParam::paramFile);
  }

  parameters_[ClpParam::EXPORTFILE]->setup(
      "exportF!ile", "Sets name for file to export model to",
      "This will set the name of the model will be written to and read from. "
      "This is initialized to 'export.mps'. ",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::GMPLSOLFILE]->setup(
      "gmplSolutionF!ile", "Sets name for file to store GMPL solution in",
      "This will set the name the GMPL solution will be written to and read "
      "from. This is initialized to 'gmpl.sol'. ",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::IMPORTFILE]->setup(
      "importF!ile", "Sets name for file to import model from",
      "This will set the name of the model to be read in with the import "
      "command. This is initialized to 'import.mps'",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::MODELFILE]->setup(
      "modelF!ile", "Sets name for file to store model in",
      "This will set the name the model will be written to and read from. "
      "This is initialized to 'prob.mod'. ",
      CoinParam::displayPriorityHigh);

  parameters_[ClpParam::SOLUTIONBINARYFILE]->setup(
      "solBinaryF!ile",
      "Sets name for file to store solution in binary format",
      "This will set the name the solution will be saved to and read from. "
      "By default, binary solutions are written to 'solution.file'."
      "use printSolution.", CoinParam::displayPriorityHigh);

  parameters_[ClpParam::SOLUTIONFILE]->setup(
      "solF!ile", "Sets name for file to store solution in",
      "This will set the name the solution will be saved to and read from. "
      "By default, solutions are written to 'opt.sol'. To print to stdout, "
      "use printSolution.", CoinParam::displayPriorityHigh);
}

//###########################################################################
//###########################################################################

void ClpParameters::addClpDirParams() {

  if (cbcMode_){
     return;
  }

  for (int code = ClpParam::FIRSTDIRECTORYPARAM + 1;
       code < ClpParam::LASTDIRECTORYPARAM; code++) {
    getParam(code)->setType(CoinParam::paramDir);
  }

  parameters_[ClpParam::DIRECTORY]->setup(
      "directory", "Set Default directory for import etc.",
      "This sets the directory which import, export, saveModel, restoreModel "
      "etc. will use. It is initialized to the current directory.");

  parameters_[ClpParam::DIRSAMPLE]->setup(
      "dirSample", "Set directory where the COIN-OR sample problems are.",
      "This sets the directory where the COIN-OR sample problems reside. It is "
      "used only when -unitTest is passed to clp. clp will pick up the test "
      "problems from this directory. It is initialized to '../../Data/Sample'",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::DIRNETLIB]->setup(
      "dirNetlib", "Set directory where the netlib problems are.",
      "This sets the directory where the netlib problems reside. One can get "
      "the netlib problems from COIN-OR or from the main netlib site. This "
      "parameter is used only when -netlib is passed to cbc. cbc will pick up "
      "the netlib problems from this directory. If clp is built without zlib "
      "support then the problems must be uncompressed. It is initialized to "
      "'../../Data/Netlib'",
      CoinParam::displayPriorityLow);

}

//###########################################################################
//###########################################################################

void ClpParameters::addClpStrParams() {

  if (cbcMode_){
     return;
  }

  for (int code = ClpParam::FIRSTSTRINGPARAM + 1;
       code < ClpParam::LASTSTRINGPARAM; code++) {
    getParam(code)->setType(CoinParam::paramStr);
  }

  parameters_[ClpParam::PRINTMASK]->setup(
      "printM!ask", "Control printing of solution with a regular expression",
      "If set then only those names which match mask are printed in a "
      "solution. '?' matches any character and '*' matches any set of "
      "characters.  The default is '' (unset) so all variables are printed. "
      "This is only active if model has names.");

}

//###########################################################################
//###########################################################################

void ClpParameters::addClpKwdParams() {
  for (int code = ClpParam::FIRSTKWDPARAM + 1; code < ClpParam::LASTKWDPARAM;
       code++) {
    getParam(code)->setPushFunc(ClpParamUtils::pushClpKwdParam);
    getParam(code)->setType(CoinParam::paramKwd);
  }


  parameters_[ClpParam::COMMANDPRINTLEVEL]->setup(
      "allC!ommands", "What priority level of commands to print",
      "For the sake of your sanity, only the more useful and simple commands "
      "are printed out on ?.");
  parameters_[ClpParam::COMMANDPRINTLEVEL]->appendKwd("all");
  parameters_[ClpParam::COMMANDPRINTLEVEL]->appendKwd("more");
  parameters_[ClpParam::COMMANDPRINTLEVEL]->appendKwd("important");

  parameters_[ClpParam::BIASLU]->setup(
      "biasLU", "Whether factorization biased towards U", "",
      CoinParam::displayPriorityNone);
  parameters_[ClpParam::BIASLU]->appendKwd("UU");
  parameters_[ClpParam::BIASLU]->appendKwd("UX");
  parameters_[ClpParam::BIASLU]->appendKwd("LX");
  parameters_[ClpParam::BIASLU]->appendKwd("LL");

  parameters_[ClpParam::BARRIERSCALE]->setup(
      "bscale", "Whether to scale in barrier (and ordering speed)",
      "", CoinParam::displayPriorityNone);
  parameters_[ClpParam::BARRIERSCALE]->appendKwd("off");
  parameters_[ClpParam::BARRIERSCALE]->appendKwd("on");
  parameters_[ClpParam::BARRIERSCALE]->appendKwd("off1");
  parameters_[ClpParam::BARRIERSCALE]->appendKwd("on1");
  parameters_[ClpParam::BARRIERSCALE]->appendKwd("off2");
  parameters_[ClpParam::BARRIERSCALE]->appendKwd("on2");

  parameters_[ClpParam::CHOLESKY]->setup(
      "chol!esky", "Which cholesky algorithm",
      "For a barrier code to be effective it needs a good Cholesky ordering "
      "and factorization. The native ordering and factorization is not state "
      "of the art, although acceptable. You may want to link in one from "
      "another source.  See Makefile.locations for some possibilities.");
  parameters_[ClpParam::CHOLESKY]->appendKwd("native");
  parameters_[ClpParam::CHOLESKY]->appendKwd("dense");
  parameters_[ClpParam::CHOLESKY]->appendKwd("fudge!Long_dummy");
  parameters_[ClpParam::CHOLESKY]->appendKwd("wssmp_dummy");
#if defined(CLP_HAS_AMD) || defined(CLP_HAS_CHOLMOD)
  parameters_[ClpParam::CHOLESKY]->appendKwd("Uni!versityOfFlorida");
#else
  parameters_[ClpParam::CHOLESKY]->appendKwd("Uni!versityOfFlorida_dummy");
#endif
  parameters_[ClpParam::CHOLESKY]->appendKwd("Taucs_dummy");
  parameters_[ClpParam::CHOLESKY]->appendKwd("Mumps_dummy");
#ifdef PARDISO_BARRIER
  parameters_[ClpParam::CHOLESKY]->appendKwd("Pardiso");
#else
  parameters_[ClpParam::CHOLESKY]->appendKwd("Pardiso_dummy");
#endif

  parameters_[ClpParam::CRASH]->setup(
      "crash", "Whether to create basis for problem",
      "If crash is set to 'on' and there is an all slack basis then Clp will "
      "flip or put structural variables into the basis with the aim of getting "
      "dual feasible.  On average, dual simplex seems to perform better "
      "without it and there are alternative types of 'crash' for primal "
      "simplex, e.g. 'idiot' or 'sprint'. A variant due to Solow and Halim "
      "which is as 'on' but just flips is also available.");
  parameters_[ClpParam::CRASH]->appendKwd("off");
  parameters_[ClpParam::CRASH]->appendKwd("on");
  parameters_[ClpParam::CRASH]->appendKwd("so!low_halim");
  parameters_[ClpParam::CRASH]->appendKwd("lots");
  parameters_[ClpParam::CRASH]->appendKwd("free");
  parameters_[ClpParam::CRASH]->appendKwd("zero");
  parameters_[ClpParam::CRASH]->appendKwd("single!ton");
#ifdef CLP_INHERIT_MODE
  parameters_[ClpParam::CRASH]->appendKwd("dual");
  parameters_[ClpParam::CRASH]->appendKwd("dw");
  parameters_[ClpParam::CRASH]->appendKwd("idiot");
#else
  parameters_[ClpParam::CRASH]->appendKwd("idiot1");
  parameters_[ClpParam::CRASH]->appendKwd("idiot2");
  parameters_[ClpParam::CRASH]->appendKwd("idiot3");
  parameters_[ClpParam::CRASH]->appendKwd("idiot4");
  parameters_[ClpParam::CRASH]->appendKwd("idiot5");
  parameters_[ClpParam::CRASH]->appendKwd("idiot6");
  parameters_[ClpParam::CRASH]->appendKwd("idiot7");
#endif

  parameters_[ClpParam::CROSSOVER]->setup(
      "cross!over",
      "Whether to get a basic solution with the simplex algorithm after the "
      "barrier algorithm finished",
      "Interior point algorithms do not obtain a basic solution. This option "
      "will crossover to a basic solution suitable for ranging or branch and "
      "cut. With the current state of the solver for quadratic programs it may "
      "be a good idea to switch off crossover for this case (and maybe "
      "presolve as well) - the option 'maybe' does this.");
  parameters_[ClpParam::CROSSOVER]->appendKwd("off");
  parameters_[ClpParam::CROSSOVER]->appendKwd("on");
  parameters_[ClpParam::CROSSOVER]->appendKwd("maybe");
  parameters_[ClpParam::CROSSOVER]->appendKwd("presolve");

  parameters_[ClpParam::DIRECTION]->setup(
      "direction", "Minimize or Maximize",
      "The default is minimize - use 'direction maximize' for maximization.\n "
      "You can also use the parameters 'maximize' or 'minimize'.");
  parameters_[ClpParam::DIRECTION]->appendKwd("min!imize");
  parameters_[ClpParam::DIRECTION]->appendKwd("max!imize");
  parameters_[ClpParam::DIRECTION]->appendKwd("zero");

  parameters_[ClpParam::DUALPIVOT]->setup(
      "dualP!ivot", "Dual pivot choice algorithm",
      "The Dantzig method is simple but its use is deprecated.  Steepest is "
      "the method of choice and there are two variants which keep all weights "
      "updated but only scan a subset each iteration. Partial switches this on "
      "while automatic decides at each iteration based on information about "
      "the factorization. The PE variants add the Positive Edge criterion. "
      "This selects incoming variables to try to avoid degenerate moves. See "
      "also option psi.",
      CoinParam::displayPriorityLow);
  parameters_[ClpParam::DUALPIVOT]->appendKwd("auto!matic");
  parameters_[ClpParam::DUALPIVOT]->appendKwd("dant!zig");
  parameters_[ClpParam::DUALPIVOT]->appendKwd("partial");
  parameters_[ClpParam::DUALPIVOT]->appendKwd("steep!est");
  parameters_[ClpParam::DUALPIVOT]->appendKwd("PEsteep!est");
  parameters_[ClpParam::DUALPIVOT]->appendKwd("PEdantzig");

  parameters_[ClpParam::FACTORIZATION]->setup(
      "fact!orization", "Which factorization to use",
      "The default is to use the normal CoinFactorization, but other choices "
      "are a dense one, OSL's, or one designed for small problems."
  );
  parameters_[ClpParam::FACTORIZATION]->appendKwd("normal");
  parameters_[ClpParam::FACTORIZATION]->appendKwd("dense");
  parameters_[ClpParam::FACTORIZATION]->appendKwd("simple");
  parameters_[ClpParam::FACTORIZATION]->appendKwd("osl");

  parameters_[ClpParam::GAMMA]->setup("gamma!(Delta)",
                                     "Whether to regularize barrier",
                                     "", CoinParam::displayPriorityLow);
  parameters_[ClpParam::GAMMA]->appendKwd("off");
  parameters_[ClpParam::GAMMA]->appendKwd("on");
  parameters_[ClpParam::GAMMA]->appendKwd("gamma");
  parameters_[ClpParam::GAMMA]->appendKwd("delta");
  parameters_[ClpParam::GAMMA]->appendKwd("onstrong");
  parameters_[ClpParam::GAMMA]->appendKwd("gammastrong");
  parameters_[ClpParam::GAMMA]->appendKwd("deltastrong");

  parameters_[ClpParam::PRESOLVE]->setup(
      "presolve", "Whether to presolve problem",
      "Presolve analyzes the model to find such things as redundant equations, "
      "equations which fix some variables, equations which can be transformed "
      "into bounds, etc. For the initial solve of any problem this is worth "
      "doing unless one knows that it will have no effect. Option 'on' will "
      "normally do 5 passes, while using 'more' will do 10.  If the problem is "
      "very large one can let CLP write the original problem to file by using "
      "'file'.");
  parameters_[ClpParam::PRESOLVE]->appendKwd("on");
  parameters_[ClpParam::PRESOLVE]->appendKwd("off");
  parameters_[ClpParam::PRESOLVE]->appendKwd("more");
  parameters_[ClpParam::PRESOLVE]->appendKwd("file");

  parameters_[ClpParam::PRIMALPIVOT]->setup(
      "primalP!ivot", "Primal pivot choice algorithm",
      "The Dantzig method is simple but its use is deprecated.  Exact devex is "
      "the method of choice and there are two variants which keep all weights "
      "updated but only scan a subset each iteration. Partial switches this on "
      "while 'change' initially does 'dantzig' until the factorization becomes "
      "denser. This is still a work in progress. The PE variants add the "
      "Positive Edge criterion. This selects incoming variables to try to "
      "avoid degenerate moves. See also Towhidi, M., Desrosiers, J., Soumis, "
      "F., The positive edge criterion within COIN-OR's CLP; Omer, J., "
      "Towhidi, M., Soumis, F., The positive edge pricing rule for the dual "
      "simplex.",
      CoinParam::displayPriorityLow);
  parameters_[ClpParam::PRIMALPIVOT]->appendKwd("auto!matic");
  parameters_[ClpParam::PRIMALPIVOT]->appendKwd("exa!ct");
  parameters_[ClpParam::PRIMALPIVOT]->appendKwd("dant!zig");
  parameters_[ClpParam::PRIMALPIVOT]->appendKwd("part!ial");
  parameters_[ClpParam::PRIMALPIVOT]->appendKwd("steep!est");
  parameters_[ClpParam::PRIMALPIVOT]->appendKwd("change");
  parameters_[ClpParam::PRIMALPIVOT]->appendKwd("sprint");
  parameters_[ClpParam::PRIMALPIVOT]->appendKwd("PEsteep!est");
  parameters_[ClpParam::PRIMALPIVOT]->appendKwd("PEdantzig");

  parameters_[ClpParam::INTPRINT]->setup(
      "printi!ngOptions", "Print options",
      "This changes the amount and format of printing a solution:\nnormal - "
      "nonzero column variables \n integer - nonzero integer column "
      "variables\n special - in format suitable for OsiRowCutDebugger\n rows - "
      "nonzero column variables and row activities\n all - all column "
      "variables and row activities.\n\n For non-integer problems 'integer' "
      "and 'special' act like 'normal'. Also see printMask for controlling "
      "output.");
  parameters_[ClpParam::INTPRINT]->appendKwd("normal");
  parameters_[ClpParam::INTPRINT]->appendKwd("integer");
  parameters_[ClpParam::INTPRINT]->appendKwd("special");
  parameters_[ClpParam::INTPRINT]->appendKwd("rows");
  parameters_[ClpParam::INTPRINT]->appendKwd("all");
  parameters_[ClpParam::INTPRINT]->appendKwd("csv");
  parameters_[ClpParam::INTPRINT]->appendKwd("bound!ranging");
  parameters_[ClpParam::INTPRINT]->appendKwd("rhs!ranging");
  parameters_[ClpParam::INTPRINT]->appendKwd("objective!ranging");
  parameters_[ClpParam::INTPRINT]->appendKwd("stats");
  parameters_[ClpParam::INTPRINT]->appendKwd("boundsint");
  parameters_[ClpParam::INTPRINT]->appendKwd("boundsall");
  parameters_[ClpParam::INTPRINT]->appendKwd("fixint");
  parameters_[ClpParam::INTPRINT]->appendKwd("fixall");

  parameters_[ClpParam::SCALING]->setup(
      "scal!ing", "Whether to scale problem",
      "Scaling can help in solving problems which might otherwise fail because of "
      "lack of accuracy.  It can also reduce the number of iterations. "
      "It is not applied if the range of elements is small.  When the solution "
      "is evaluated in the unscaled problem, it is possible that small primal and/or "
      "dual infeasibilities occur. \n"
      " - 'equilibrium' uses the largest element for scaling. \n"
      " - 'geometric' uses the squareroot of the product of largest and smallest element.\n"
      " - 'auto' lets CLP choose a method that gives the best ratio of the largest element "
      "to the smallest one.", CoinParam::displayPriorityLow);
  parameters_[ClpParam::SCALING]->appendKwd("off");
  parameters_[ClpParam::SCALING]->appendKwd("equi!librium");
  parameters_[ClpParam::SCALING]->appendKwd("geo!metric");
  parameters_[ClpParam::SCALING]->appendKwd("auto!matic");
  parameters_[ClpParam::SCALING]->appendKwd("dynamic");
  parameters_[ClpParam::SCALING]->appendKwd("rows!only");

#ifndef COIN_AVX2
  parameters_[ClpParam::VECTOR]->setup(
      "vector", "Whether to use vector? Form of matrix in simplex",
      "If this is on ClpPackedMatrix uses extra column copy in odd format.",
      CoinParam::displayPriorityLow);
  parameters_[ClpParam::VECTOR]->appendKwd("off");
  parameters_[ClpParam::VECTOR]->appendKwd("on");
#else
  parameters_[ClpParam::VECTOR]->setup(
      "vector", "Try and use vector instructions in simplex",
      "At present only for Intel architectures - but could be extended. Uses "
      "avx2 or avx512 instructions. Uses different storage for matrix - can be "
      "of benefit without instruction set on some problems. I may add pool to "
      "switch on a pool matrix",
      CoinParam::displayPriorityLow);
  parameters_[ClpParam::VECTOR]->appendKwd("off");
  parameters_[ClpParam::VECTOR]->appendKwd("on");
  parameters_[ClpParam::VECTOR]->appendKwd("ones");
#endif
}

//###########################################################################
//###########################################################################

void ClpParameters::addClpDblParams() {
  for (int code = ClpParam::FIRSTDBLPARAM + 1; code < ClpParam::LASTDBLPARAM;
       code++) {
    getParam(code)->setPushFunc(ClpParamUtils::pushClpDblParam);
    getParam(code)->setType(CoinParam::paramDbl);
  }

  parameters_[ClpParam::DUALBOUND]->setup(
      "dualB!ound",
      "Initially algorithm acts as if no gap between bounds exceeds this value",
      1.0e-20, 1.0e20,
      "The dual algorithm in Clp is a single phase algorithm as opposed to a "
      "two phase algorithm where you first get feasible then optimal.  If a "
      "problem has both upper and lower bounds then it is trivial to get dual "
      "feasible by setting non basic variables to correct bound.  If the gap "
      "between the upper and lower bounds of a variable is more than the value "
      "of dualBound Clp introduces fake bounds so that it can make the problem "
      "dual feasible.  This has the same effect as a composite objective "
      "function in the primal algorithm.  Too high a value may mean more "
      "iterations, while too low a bound means the code may go all the way and "
      "then have to increase the bounds.  OSL had a heuristic to adjust "
      "bounds, maybe we need that here.");

  parameters_[ClpParam::DUALTOLERANCE]->setup("dualT!olerance", "For an optimal solution no dual infeasibility may exceed this value", 1.0e-20, COIN_DBL_MAX, "Normally the default tolerance is fine, but one may want to increase it a bit if the dual simplex algorithm seems to be having a hard time. One method which can be faster is to use a large tolerance e.g. 1.0e-4 and the dual simplex algorithm and then to clean up the problem using the primal simplex algorithm with the correct tolerance (remembering to switch off presolve for this final short clean up phase).");

  parameters_[ClpParam::FAKEBOUND]->setup(
	  "fakeB!ound", "All bounds <= this value - DEBUG", 1.0, 1.0e15,
	   "",CoinParam::displayPriorityLow);

  parameters_[ClpParam::OBJSCALE]->setup(
      "objective!Scale", "Scale factor to apply to objective", -COIN_DBL_MAX,
      COIN_DBL_MAX,
      "If the objective function has some very large values, you may wish to "
      "scale them internally by this amount.  It can also be set by autoscale. "
      "It is applied after scaling.  You are unlikely to need this.",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::PRESOLVETOLERANCE]->setup(
      "preT!olerance", "Tolerance to use in presolve", 1.0e-20, COIN_DBL_MAX,
      "One may want to increase this tolerance if presolve says the problem is "
      "infeasible and one has awkward numbers and is sure that the problem is "
      "really feasible.");

  parameters_[ClpParam::PRIMALTOLERANCE]->setup(
      "primalT!olerance",
      "For a feasible solution no primal infeasibility, i.e., constraint "
      "violation, may exceed this value",
      1.0e-20, COIN_DBL_MAX,
      "Normally the default tolerance is fine, but one may want to increase it "
      "a bit if the primal simplex algorithm seems to be having a hard time.");

  parameters_[ClpParam::PRIMALWEIGHT]->setup(
      "primalW!eight",
      "Initially algorithm acts as if it costs this much to be infeasible",
      1.0e-20, COIN_DBL_MAX,
      "The primal algorithm in Clp is a single phase algorithm as opposed to a "
      "two phase algorithm where you first get feasible then optimal.  So Clp "
      "is minimizing this weight times the sum of primal infeasibilities plus "
      "the true objective function (in minimization sense). Too high a value "
      "may mean more iterations, while too low a value means the algorithm may "
      "iterate into the wrong directory for long and then has to increase the "
      "weight in order to get feasible."); // OSL had a heuristic to adjust
                                           // bounds, maybe we need that here.

  parameters_[ClpParam::PSI]->setup(
      "psi", "Two-dimension pricing factor for Positive Edge criterion", -1.1,
      1.1,
      "The Positive Edge criterion has been added to select incoming variables "
      "to try and avoid degenerate moves. Variables not in the promising set "
      "have their infeasibility weight multiplied by psi, so 0.01 would mean "
      "that if there were any promising variables, then they would always be "
      "chosen, while 1.0 effectively switches the algorithm off. There are "
      "two ways of switching this feature on. One way is to set psi to a "
      "positive value and then the Positive Edge criterion will be used for "
      "both primal and dual simplex. The other way is to select PEsteepest in "
      "dualpivot choice (for example), then the absolute value of psi is used. "
      "Code donated by Jeremy Omer. See Towhidi, M., Desrosiers, J., Soumis, F., "
      "The positive edge criterion within COIN-OR's CLP; Omer, J., Towhidi, M., "
      "Soumis, F., The positive edge pricing rule for the dual simplex."); // Until this settles down it is only implemented in CLP.

  parameters_[ClpParam::PROGRESS]->setup(
      "progress!Interval", "Time interval for printing progress",
      -COIN_DBL_MAX, COIN_DBL_MAX,
      "This sets a minimum interval for some printing - elapsed seconds");

  parameters_[ClpParam::OBJSCALE2]->setup(
      "reallyO!bjectiveScale", "Scale factor to apply to objective in place",
      -COIN_DBL_MAX, COIN_DBL_MAX,
      "You can set this to -1.0 to test maximization or other to stress code",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::RHSSCALE]->setup(
      "rhs!Scale", "Scale factor to apply to rhs and bounds", -COIN_DBL_MAX,
      COIN_DBL_MAX,
      "If the rhs or bounds have some very large meaningful values, you may "
      "wish to scale them internally by this amount.  It can also be set by "
      "autoscale.  This should not be needed.",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::TIMELIMIT]->setup(
      "sec!onds", "Maximum seconds", -1.0, COIN_DBL_MAX,
      "After this many seconds clp will act as if maximum iterations had been "
      "reached (if value >=0).");

  parameters_[ClpParam::ZEROTOLERANCE]->setup(
      "zeroT!olerance",
      "Kill all coefficients whose absolute value is less than this value",
      1.0e-100, 1.0e-5,
      "This applies to reading mps files (and also lp files if "
      "KILL_ZERO_READLP defined)");
}

//###########################################################################
//###########################################################################

void ClpParameters::addClpIntParams() {
  for (int code = ClpParam::FIRSTINTPARAM + 1; code < ClpParam::LASTCLPINTPARAM;
       code++) {
    getParam(code)->setPushFunc(ClpParamUtils::pushClpIntParam);
    getParam(code)->setType(CoinParam::paramInt);
  }

  parameters_[ClpParam::CPP]->setup(
      "cpp!Generate", "Generates C++ code", -1, 50000,
      "Once you like what the stand-alone solver does then this allows you to "
      "generate user_driver.cpp which approximates the code. 0 gives simplest "
      "driver, 1 generates saves and restores, 2 generates saves and restores "
      "even for variables at default value. 4 bit in cbc generates size "
      "dependent code rather than computed values. This is now deprecated as "
      "you can call stand-alone solver - see Cbc/examples/driver4.cpp.",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::DECOMPOSE_BLOCKS]->setup(
      "decomp!ose", "Whether to try decomposition", -COIN_INT_MAX, COIN_INT_MAX,
      "0 - off, 1 choose blocks >1 use as blocks Dantzig Wolfe if primal, "
      "Benders if dual - uses sprint pass for number of passes",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::DENSE]->setup(
      "dense!Threshold", "Threshold for using dense factorization", -1, 10000,
      "If processed problem <= this use dense factorization",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::DUALIZE]->setup("dualize", "Solves dual reformulation",
                                       0, 4, "Don't even think about it.",
                                       CoinParam::displayPriorityLow);

  parameters_[ClpParam::IDIOT]->setup(
      "idiot!Crash", "Whether to try idiot crash", -1, COIN_INT_MAX,
      "This is a type of 'crash' which works well on some homogeneous "
      "problems. It works best on problems with unit elements and rhs but will "
      "do something to any model.  It should only be used before the primal "
      "simplex algorithm.  It can be set to -1 when the code decides for "
      "itself whether to use it, 0 to switch off, or n > 0 to do n passes.");

  parameters_[ClpParam::MAXFACTOR]->setup(
      "maxF!actor", "Maximum number of iterations between refactorizations", 1,
      COIN_INT_MAX,
      "If this is left at its default value of 200 then CLP will guess a  "
      "value to use.  CLP may decide to re-factorize earlier for accuracy.");

  parameters_[ClpParam::MAXITERATION]->setup(
      "maxIt!erations", "Maximum number of iterations before stopping", 0,
      COIN_INT_MAX,
      "This can be used for testing purposes.  The corresponding library "
      "call\n \tsetMaximumIterations(value)\n can be useful.  If the code "
      "stops on seconds or by an interrupt this will be treated as stopping on "
      "maximum iterations. This is ignored in branchAndCut - use maxN!odes.");

  parameters_[ClpParam::MORESPECIALOPTIONS]->setup(
      "moreS!pecialOptions", "Yet more dubious options for Simplex", 0,
      COIN_INT_MAX,
      "See ClpSimplex.hpp.");
  parameters_[ClpParam::MORESPECIALOPTIONS]->appendKwd("keep!DualOrPrimal#If you ask for dual you will always get dual (and for primal)",8192);
  parameters_[ClpParam::MORESPECIALOPTIONS]->appendKwd("clean!Scaled#Make sure unscaled problem is feasible if scaled problem is feasible",134217728);

  parameters_[ClpParam::PRESOLVEPASS]->setup(
      "passP!resolve", "How many passes in presolve", -200, 100,
      "Normally Presolve does 10 passes but you may want to do less to make it "
      "more lightweight or do more if improvements are still being made.  As "
      "Presolve will return if nothing is being taken out, you should not "
      "normally need to use this fine tuning.",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::PERTVALUE]->setup("pertV!alue", "Method of perturbation",
                                         -5000, 102, "",
                                         CoinParam::displayPriorityLow);

  parameters_[ClpParam::RANDOMSEED]->setup(
      "randomS!eed", "Random seed for Clp", 0, COIN_INT_MAX,
      "Initialization of the random seed for pseudo-random numbers used to "
      "break ties in degenerate problems. This may yield a different "
      "continuous optimum and, in the context of Cbc, different cuts and "
      "heuristic solutions. The special value of 0 lets CLP use the time of "
      "the day for the initial seed.");

  parameters_[ClpParam::SLPVALUE]->setup(
      "slp!Value", "Number of slp passes before primal", -50000, 50000,
      "If you are solving a quadratic problem using primal then it may be "
      "helpful to do some sequential Lps to get a good approximate solution.",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::SMALLFACT]->setup(
      "small!Factorization", "Threshold for using small factorization", -1,
      10000, "If processed problem <= this use small factorization",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::SPECIALOPTIONS]->setup(
      "special!Options", "Dubious options for Simplex - see ClpSimplex.hpp", 0,
      COIN_INT_MAX, "", CoinParam::displayPriorityLow);

  parameters_[ClpParam::SPRINT]->setup(
      "sprint!Crash", "Whether to try sprint crash", -1, COIN_INT_MAX,
      "For long and thin problems this method may solve a series of small "
      "problems created by taking a subset of the columns.  The idea as "
      "'Sprint' was introduced by J. Forrest after an LP code of that name of "
      "the 60's which tried the same tactic (not totally successfully). CPLEX "
      "calls it 'sifting'.  -1 lets CLP automatically choose the number of "
      "passes, 0 is off, n is number of passes");

  parameters_[ClpParam::SUBSTITUTION]->setup(
      "subs!titution", "How long a column to substitute for in presolve", 0,
      10000,
      "Normally Presolve gets rid of 'free' variables when there are no more "
      "than 3 coefficients in a row.  If you increase this, the number of rows "
      "may decrease but the number of coefficients may increase.",
      CoinParam::displayPriorityNone);

#ifdef CLP_THREAD
  parameters_[ClpParam::THREADS]->setup(
      "thread!s", "Number of threads to try and use", -100, 100000,
      "To use multiple threads, set threads to number wanted.  It may be "
      "better to use one or two more than number of cpus available.  If 100+n "
      "then n threads and search is repeatable (maybe be somewhat slower), if "
      "200+n use threads for root cuts, 400+n threads used in sub-trees.");
#endif

  if (cbcMode_){
     return;
  }

  for (int code = ClpParam::LASTCLPINTPARAM + 1;
       code < ClpParam::LASTINTPARAM; code++) {
    getParam(code)->setPushFunc(ClpParamUtils::pushClpIntParam);
    getParam(code)->setType(CoinParam::paramInt);
  }

  parameters_[ClpParam::LOGLEVEL]->setup(
      "log!Level", "Level of detailin Clp output", -63, 63,
      "If 0 then there should be no output in normal circumstances.  1 is "
      "probably the best value for most uses, while 2 and 3 give more "
      "information.");

  parameters_[ClpParam::OUTPUTFORMAT]->setup(
      "output!Format", "Which output format to use", 1, 6,
      "Normally export will be done using normal representation for numbers "
      "and two values per line.  You may want to do just one per line (for "
      "grep or suchlike) and you may wish to save with absolute accuracy using "
      "a coded version of the IEEE value. A value of 2 is normal. otherwise "
      "odd values gives one value per line, even two.  Values 1,2 give normal "
      "format, 3,4 gives greater precision, while 5,6 give IEEE values.  When "
      "used for exporting a basis 1 does not save values, 2 saves values, 3 "
      "with greater accuracy and 4 in IEEE.");

  parameters_[ClpParam::PRINTOPTIONS]->setup(
      "pO!ptions", "Dubious print options", 0, COIN_INT_MAX,
      "If this is > 0 then presolve will give more information and branch and "
      "cut will give statistics",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::PROGRESSITER]->setup(
      "progressIter!ations",
      "Print progress every N iterations (0 = time-based only)", 0, COIN_INT_MAX,
      "When set to a positive value, prints a progress row every N iterations "
      "instead of (or in addition to) time-based printing. "
      "Useful for deterministic output that does not depend on machine speed.");

  parameters_[ClpParam::VERBOSE]->setup(
      "verbose", "Switches on longer help on single ?", 0, 31,
      "Set to 1 to get short help with ? list, 2 to get long help, 3 for both. "
      " (add 4 to just get ampl ones).",
      CoinParam::displayPriorityNone);
}

//###########################################################################
//###########################################################################

void ClpParameters::addClpBoolParams() {
  for (int code = ClpParam::FIRSTBOOLPARAM + 1; code < ClpParam::LASTBOOLPARAM;
       code++) {
    getParam(code)->setType(CoinParam::paramKwd);
    getParam(code)->appendKwd("off", ClpParameters::ParamOff);
    getParam(code)->appendKwd("on", ClpParameters::ParamOn);
    getParam(code)->setPushFunc(ClpParamUtils::pushClpKwdParam);
  }

  parameters_[ClpParam::AUTOSCALE]->setup(
      "auto!Scale",
      "Whether to scale objective, rhs and bounds of problem if they look odd",
      "If you think you may get odd objective values or large equality rows "
      "etc then it may be worth setting this true.  It is still experimental "
      "and you may prefer to use objective!Scale and rhs!Scale.",
      CoinParam::displayPriorityNone);

  parameters_[ClpParam::BUFFER_MODE]->setup(
      "buff!eredMode", "Whether to flush print buffer",
      "Default is on, off switches on unbuffered output");

  parameters_[ClpParam::ERRORSALLOWED]->setup(
      "error!sAllowed", "Whether to allow import errors",
      "The default is not to use any model which had errors when reading the "
      "mps file. Setting this to 'on' will allow all errors from which the "
      "code can recover simply by ignoring the error.  There are some errors "
      "from which the code can not recover e.g. no ENDATA.  This has to be set "
      "before import i.e. -errorsAllowed on -import xxxxxx.mps.");

  parameters_[ClpParam::KEEPNAMES]->setup(
      "keepN!ames", "Whether to keep names from import",
      "It saves space to get rid of names so if you need to you can set this "
      "to off. This needs to be set before the import of model - so -keepnames "
      "off -import xxxxx.mps.");

  parameters_[ClpParam::KKT]->setup(
      "KKT", "Whether to use KKT factorization in barrier", "",
      CoinParam::displayPriorityLow);

  parameters_[ClpParam::MESSAGES]->setup(
      "mess!ages", "Controls if Clpnnnn is printed",
      "The default behavior is to put out messages such as:\n Clp0005 2261  "
      "Objective 109.024 Primal infeas 944413 (758)\n but this program turns "
      "this off to make it look more friendly.  It can be useful to turn them "
      "back on if you want to be able to 'grep' for particular messages or if "
      "you intend to override the behavior of a particular message.");

  parameters_[ClpParam::PERTURBATION]->setup(
      "perturb!ation", "Whether to perturb the problem",
      "Perturbation helps to stop cycling, but CLP uses other measures for "
      "this. However, large problems and especially ones with unit elements "
      "and unit right hand sides or costs benefit from perturbation.  Normally "
      "CLP tries to be intelligent, but one can switch this off.");

  parameters_[ClpParam::PFI]->setup(
      "PFI", "Whether to use Product Form of Inverse in simplex",
      "By default clp uses Forrest-Tomlin L-U update.  If you are masochistic "
      "you can switch it off.",
      CoinParam::displayPriorityNone);

  parameters_[ClpParam::SPARSEFACTOR]->setup(
      "spars!eFactor", "Whether factorization treated as sparse", "",
      CoinParam::displayPriorityNone);
}

//###########################################################################
//###########################################################################
